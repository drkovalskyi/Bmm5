"""
L1 Trigger Efficiency Study

This script measures the L1 trigger turn-on curves as a function of
offline muon pT in bins of eta.  The information is extracted from the
matched L1 trigger objects using a binned ML fit.

Sample selection requirements:
- Use primary datasets (PDs) built using non-muon triggers.
- Explicitly specify triggers to be used for selecting unbiased events
  in complex PDs like EGamma, which use L1_Mu6_DoubleEG12er2p5 as a
  seed for diphoton triggers.
- Use certified data if possible to avoid surprises. ZeroBias samples
  may be not certified.


Sample format:
- Bmm NanoAOD or its skim containing the following branches
  - MuonId_.*
  - Muon_.*
  - HLT_.*

Offline muon selection:
- Loose muon ID
- isGlobalMuon - to match Outside-In HLT algorithm
- isTrackerMuon - to clean up and match Inside-Out HLT algorithm

It is recommended to repeat the measurement using tighter muon
selection requirements, such as the Soft MVA ID, to verify the
stability of the results. Ideally, the actual analysis muon selection
should be used.

"""
import os
import sys
import json
import re
import ROOT
import tdrstyle
from array import array

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)
ROOT.gPad.SetGrid()
# ROOT.gPad.SetLogx()

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/"
output_path = "/eos/home-d/dmytro/www/plots/Run3/l1_turn-on_curves/";

samples = {
    'Bsmm-22EE': {
        'samples': [
            data_path + "/BsToMuMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root",
        ],
        'triggers': [], 
    },
    'ZeroBias-22E': {
        'samples': [
            data_path + "/ZeroBias+Run2022E-PromptReco-v1+MINIAOD/*.root",
        ],
        'triggers': [], 
    },
}

data = dict()

def split_range(start, end, step):
    """
    Split era range into list of bins
    Example: split_range(-2.4, 2.4, 0.2)
    [(-2.4, -2.2), (-2.2, -2.0), ..., (2.0, 2.2), (2.2, 2.4)]
    """
    result = []
    current = start
    while current < end:
        result.append((round(current, 2), round(current + step, 2)))
        current += step
    return result


def get_data(name, tree="Events"):
    if name in data:
        return data[name]
    chain = ROOT.TChain(tree)
    n_files = 0 
    for sample in samples[name]['samples']:
        n_files += chain.Add(sample)
    if n_files == 0:
        print("No files found for " + name)
        return None
    data[name] = chain
    return data[name]
        

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

hists = dict()

def measure_trigger_object_efficiency(sample, suffix, preselection):
    chain  = get_data(sample)
    if not chain: return
    name = sample + suffix

    eta_bins = split_range(-2.4, 2.4, 0.2)

    fout = ROOT.TFile.Open("results/l1_turn-on_%s.root" % name, "recreate")

    for (eta_min, eta_max) in eta_bins:
        h_name = f"h_all_eta{eta_min:+.1f}to{eta_max:+.1f}_{name}"
        h_all = ROOT.TH1F(h_name, "", 36, 2, 20)
        h_all.Sumw2()
        selection = f"Muon_eta > {eta_min} && Muon_eta < {eta_max}"
        if preselection != "":
            selection += "&&" + preselection
        chain.Draw(f"Muon_pt>>{h_name}", selection)

        print_canvas(h_name, output_path)
        h_all.SetDirectory(0)
        hists[h_name] = h_all

        h_name = f"h_trig_eta{eta_min:+.1f}to{eta_max:+.1f}_{name}"
        h_trig = h_all.Clone(h_name)
        chain.Draw(f"Muon_pt>>{h_name}", selection + "&& MuonId_l1_quality >= 12" ) # single muon quality
        print_canvas("%s" % (h_name), output_path)
        hists[h_name] = h_trig

        h_name = f"h_eff_eta{eta_min:+.1f}to{eta_max:+.1f}_{name}"
        h_eff = h_trig.Clone(h_name)
        h_eff.Divide(h_trig, h_all, 1, 1, "B")
        h_eff.SetMinimum(0)
        h_eff.SetMaximum(1.1)
        h_eff.Draw("hist")

        print_canvas("%s" % (h_name), output_path)
        
        hists[h_name] = h_eff
        
    fout.Close()

# # def _measure_trigger_object_efficiency_mm(chain, bin_name, name, eta_min, eta_max, preselection):
# #     h_name = "h_all_%s_%s"    % (bin_name, name)
# #     h_name2 = "h_all_%s_%s_2" % (bin_name, name)
# #     h_all = ROOT.TH1F(h_name, "", 32, 4, 20)
# #     h_all.Sumw2()
# #     h_all2 = ROOT.TH1F(h_name2, "", 32, 4, 20)
# #     h_all2.Sumw2()
# #     selection = preselection
# #     if selection != "":
# #         selection += "&&" 
# #     chain.Draw("%s>>%s" % ("m1pt", h_name), selection + "m1eta > %s && m1eta < %s" % (eta_min, eta_max))
# #     chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + "m2eta > %s && m2eta < %s" % (eta_min, eta_max))
# #     h_all.Add(h_all2)
# #     h_all.Draw("")

# #     print_canvas("%s" % (h_name), output_path)
# #     h_all.SetDirectory(0)
# #     hists[h_name] = h_all

# #     h_name  = "h_trig_%s_%s" % (bin_name, name)
# #     h_name2 = "h_trig_%s_%s_2" % (bin_name, name)
# #     h_trig  = h_all.Clone(h_name)
# #     h_trig2 = h_all.Clone(h_name2)
# #     chain.Draw("%s>>%s" % ("m1pt", h_name), selection + "m1_hlt_pt>0 && m1eta > %s && m1eta < %s" % (eta_min, eta_max))
# #     chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + "m2_hlt_pt>0 && m2eta > %s && m2eta < %s" % (eta_min, eta_max))
# #     h_trig.Add(h_trig2)

# #     hists[h_name] = h_trig

# #     h_name = "h_eff_%s_%s" % (bin_name, name)
# #     h_eff = h_trig.Clone(h_name)
# #     h_eff.Divide(h_trig, h_all, 1, 1, "B")
# #     h_eff.SetMinimum(0)
# #     h_eff.SetMaximum(1.1)
# #     h_eff.Draw("hist")

# #     print_canvas("%s" % (h_name), output_path)
# #     # h_eff.SetDirectory(fout)
# #     # fout.cd()
# #     # h_eff.Write()

# #     hists[h_name] = h_eff

# def _measure_trigger_object_efficiency_mm(chain, eta_bin_name, name, eta_cut, preselection, tag_selection, trigger, fout):
#     h_name = "h_all_%s_%s"    % (eta_bin_name, name)
#     h_name2 = "h_all_%s_%s_2" % (eta_bin_name, name)

#     pt_bins = [4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.25, 5.5, 5.75, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0, 20.0]

#     h_all = ROOT.TH1F(h_name, ";p_{T}", 32, 4, 20)
#     # h_all = ROOT.TH1F(h_name, ";p_{T}", len(pt_bins)-1, array('d', pt_bins))
#     h_all.Sumw2()
#     h_all2 = h_all.Clone(h_name2)
#     selection = preselection
#     if selection != "":
#         selection += "&&"
#     tag1_selection = tag_selection.format(index=1)
#     tag2_selection = tag_selection.format(index=2)

#     chain.Draw("%s>>%s" % ("m1pt", h_name),  selection + tag2_selection + "&&" + eta_cut.format(index=1))
#     chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + tag1_selection + "&&" + eta_cut.format(index=2))
#     h_all.Add(h_all2)
#     h_all.Draw("")

#     print_canvas("%s" % (h_name), output_path)
#     h_all.SetDirectory(0)
#     hists[h_name] = h_all

#     h_name  = "h_trig_%s_%s"   % (eta_bin_name, name)
#     h_name2 = "h_trig_%s_%s_2" % (eta_bin_name, name)
#     h_trig  = h_all.Clone(h_name)
#     h_trig2 = h_all.Clone(h_name2)
#     chain.Draw("%s>>%s" % ("m1pt", h_name),  selection + tag2_selection + "&& %s && %s" % (eta_cut.format(index=1), trigger))
#     chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + tag1_selection + "&& %s && %s" % (eta_cut.format(index=2), trigger))
#     h_trig.Add(h_trig2)

#     hists[h_name] = h_trig

#     h_name = "h_eff_%s_%s" % (eta_bin_name, name)
#     h_eff = h_trig.Clone(h_name)
#     h_eff.Divide(h_trig, h_all, 1, 1, "B")
#     h_eff.SetMinimum(0)
#     h_eff.SetMaximum(1.1)
#     h_eff.Draw("hist e")

#     print_canvas("%s" % (h_name), output_path)
#     h_eff.SetDirectory(fout)
#     fout.cd()
#     h_eff.Write()
#     h_eff.SetDirectory(0)
    
#     hists[h_name] = h_eff
    
    
    
# def measure_trigger_object_efficiency_mm(sample, suffix, preselection, tag_selection, trigger):
#     chain  = get_data(sample, "mm")
#     if not chain: return
#     name = sample + suffix

#     # eta_bin_step = 1
#     eta_bin_step = 5
#     eta_bins = range(-15, 15, eta_bin_step)
#     eta_bin_size = 0.1 * eta_bin_step

#     fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

#     # _measure_trigger_object_efficiency_mm(chain, "", name, -1.5, 1.5, preselection)
#     _measure_trigger_object_efficiency_mm(chain, "",      name, "abs(m{index}eta) < 1.5", preselection, tag_selection, trigger, fout)
#     _measure_trigger_object_efficiency_mm(chain, "chan0", name, "abs(m{index}eta) < 0.7", preselection, tag_selection, trigger, fout)
#     _measure_trigger_object_efficiency_mm(chain, "chan1", name, "abs(m{index}eta) > 0.7", preselection, tag_selection, trigger, fout)
    
#     # for ieta in eta_bins:
#     #     eta = ieta * 0.1
#     #     _measure_trigger_object_efficiency_mm(chain, "eta%s" % ieta, name, "m{index}eta > %s && m{index}eta < %s" % (eta, eta + eta_bin_size), preselection, tag_selection, trigger)
        
#     fout.Close()
    

# def measure_trigger_object_efficiency_2D(sample, suffix, binning, draw_options, selection):
#     chain  = get_data(sample)
#     if not chain: return
#     name = sample + suffix

#     fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

#     h_name = "h2_all_%s" % name
#     # h_all = ROOT.TH2F(h_name, "", 32, 4, 20, 30, -1.5, 1.5)
#     h_all = ROOT.TH2F(h_name, "", *binning)
#     h_all.Sumw2()
#     chain.Draw("probe_eta:probe_pt>>%s" % h_name, selection)

#     h_all.SetDirectory(0)
#     hists[h_name] = h_all

#     h_name = "h2_trig_%s" % name
#     h_trig = h_all.Clone(h_name)
#     if selection != "":
#         selection += "&&"
#     selection += "probe_hlt_pt>0"
#     chain.Draw("probe_eta:probe_pt>>%s" % h_name, selection )
#     # print_canvas("%s" % (h_name), output_path)
#     hists[h_name] = h_trig

#     h_name = "h2_eff_%s" % name
#     h_eff = h_trig.Clone(h_name)
#     h_eff.Divide(h_trig, h_all, 1, 1, "B")
#     h_eff.SetMinimum(0)
#     h_eff.SetMaximum(1.0)
#     h_eff.SetMarkerSize(0.5)
#     h_eff.Draw(draw_options)

#     print_canvas("%s" % (h_name), output_path)
#     # h_eff.SetDirectory(fout)
#     fout.cd()
#     h_eff.Write()
        
#     hists[h_name] = h_eff
        
#     fout.Close()

# def measure_trigger_object_efficiency_2D_generic(sample, suffix, var_x, var_y, binning, draw_options, selection):
#     chain  = get_data(sample)
#     if not chain: return
#     name = sample + suffix

#     h_name = "h2_all_%s" % name
#     h_all = ROOT.TH2F(h_name, "", *binning)
#     h_all.Sumw2()
#     h_all.SetMarkerSize(0.5)
#     chain.Draw("%s:%s>>%s" % (var_y, var_x, h_name), selection, draw_options)
#     print_canvas("%s" % (h_name), output_path)

#     h_all.SetDirectory(0)
#     hists[h_name] = h_all

#     h_name = "h2_trig_%s" % name
#     h_trig = h_all.Clone(h_name)
#     if selection != "":
#         selection += "&&"
#     selection += "probe_hlt_pt>0"
#     chain.Draw("%s:%s>>%s" % (var_y, var_x, h_name), selection )
#     # print_canvas("%s" % (h_name), output_path)
#     hists[h_name] = h_trig

#     h_name = "h2_eff_%s" % name
#     h_eff = h_trig.Clone(h_name)
#     h_eff.Divide(h_trig, h_all, 1, 1, "B")
#     h_eff.SetMinimum(0)
#     h_eff.SetMaximum(1)
#     h_eff.SetMarkerSize(0.5)
#     h_eff.Draw(draw_options)

#     print_canvas("%s" % (h_name), output_path)
#     # h_eff.SetDirectory(fout)
        
#     hists[h_name] = h_eff

# # def measure_trigger_object_efficiency_2D_mm(sample, suffix, binning, draw_options, selection):
# #     chain  = get_data(sample, "mm")
# #     if not chain: return
# #     name = sample + suffix
    
# #     fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

# #     h_name  = "h2_all_%s" % name
# #     h_name2 = "h2_all_%s_2" % name
# #     # h_all = ROOT.TH2F(h_name, "", 32, 4, 20, 30, -1.5, 1.5)
# #     h_all = ROOT.TH2F(h_name, "", *binning)
# #     h_all.Sumw2()
# #     h_all2 = ROOT.TH2F(h_name2, "", *binning)
# #     h_all2.Sumw2()
# #     chain.Draw("m1eta:m1pt>>%s" % h_name, selection)
# #     chain.Draw("m2eta:m2pt>>%s" % h_name2, selection)
# #     h_all.Add(h_all2)
# #     h_all.Draw("colz")
# #     print_canvas("%s" % (h_name), output_path)
    
# #     h_all.SetDirectory(0)
# #     hists[h_name] = h_all

# #     h_name  = "h2_trig_%s" % name
# #     h_name2 = "h2_trig_%s_2" % name
# #     h_trig  = h_all.Clone(h_name)
# #     h_trig2 = h_all.Clone(h_name2)
# #     if selection != "":
# #         selection += "&&"
# #     chain.Draw("m1eta:m1pt>>%s" % h_name,  selection + "m1_hlt_pt>0")
# #     chain.Draw("m2eta:m2pt>>%s" % h_name2, selection + "m2_hlt_pt>0")
# #     h_trig.Add(h_trig2)
# #     hists[h_name] = h_trig

# #     h_name = "h2_eff_%s" % name
# #     h_eff = h_trig.Clone(h_name)
# #     h_eff.Divide(h_trig, h_all, 1, 1, "B")
# #     h_eff.SetMinimum(0)
# #     h_eff.SetMaximum(1.0)
# #     h_eff.SetMarkerSize(0.5)
# #     h_eff.Draw(draw_options)

# #     print_canvas("%s" % (h_name), output_path)
# #     # h_eff.SetDirectory(fout)
# #     fout.cd()
# #     h_eff.Write()
        
# #     hists[h_name] = h_eff
        
# #     fout.Close()

# def measure_trigger_object_efficiency_2D_mm(sample, suffix, var_x, var_y, binning,
#                                             draw_options, preselection, tag_selection, trigger):
#     ## Data source
#     chain  = get_data(sample, "mm")
#     if not chain: return
    
#     name = sample + suffix

#     ## Selection
#     selection = preselection
#     if selection != "":
#         selection += "&&"
#     tag1_selection = tag_selection.format(index=1)
#     tag2_selection = tag_selection.format(index=2)

#     ## Histograms
#     h_name  = "h2_all_%s"   % name
#     h_name2 = "h2_all_%s_2" % name
#     h_all = ROOT.TH2F(h_name, "", *binning)
#     h_all.Sumw2()
#     h_all.SetMarkerSize(0.5)
#     h_all2 = ROOT.TH2F(h_name2, "", *binning)
#     h_all2.Sumw2()

#     ## Denominator
#     chain.Draw("%s:%s>>%s" % (var_y.format(index=1), var_x.format(index=1), h_name),
#                selection + tag2_selection, draw_options)
#     chain.Draw("%s:%s>>%s" % (var_y.format(index=2), var_x.format(index=2), h_name2),
#                selection + tag1_selection, draw_options)
#     h_all.Add(h_all2)
#     h_all.Draw(draw_options)
#     print_canvas("%s" % (h_name), output_path)

#     h_all.SetDirectory(0)
#     hists[h_name] = h_all

#     ## Numerator
#     h_name  = "h2_trig_%s"   % name
#     h_name2 = "h2_trig_%s_2" % name
#     h_trig  = h_all.Clone(h_name)
#     h_trig2 = h_all.Clone(h_name2)
#     chain.Draw("%s:%s>>%s" % (var_y.format(index=1), var_x.format(index=1), h_name),
#                selection + tag2_selection + "&&" + trigger, draw_options)
#     chain.Draw("%s:%s>>%s" % (var_y.format(index=2), var_x.format(index=2), h_name2),
#                selection + tag1_selection + "&&" + trigger, draw_options)
#     h_trig.Add(h_trig2)
#     # print_canvas("%s" % (h_name), output_path)
#     hists[h_name] = h_trig

#     ## Efficiency
#     h_name = "h2_eff_%s" % name
#     h_eff = h_trig.Clone(h_name)
#     h_eff.Divide(h_trig, h_all, 1, 1, "B")
#     h_eff.SetMinimum(0)
#     h_eff.SetMaximum(1)
#     h_eff.SetMarkerSize(0.5)
#     h_eff.Draw(draw_options)

#     print_canvas("%s" % (h_name), output_path)
#     # h_eff.SetDirectory(fout)
        
#     hists[h_name] = h_eff
    
    
# # def fit_efficiency(sample):
# #     name = sample
# #     ROOT.gPad.SetGrid()

# #     eta_bin_size = 0.5
# #     eta_bins = np.arange(-1.5, 1.5, eta_bin_size)
# #     fits = dict()
    
# #     fin = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name)

# #     m_erf = ROOT.TF1("m_erf", "[0]*([1]*(TMath::Erf((x - [2])*[3])+1) + (1-[1])*(TMath::Erf((x - [4])*([5]))+1))", 0 , 30);
# #     m_erf.SetParLimits(1,0,1)
# #     m_erf.SetParLimits(2,1,10)
# #     m_erf.SetParLimits(3,0,2)
# #     m_erf.SetParLimits(4,1,10)
# #     m_erf.SetParLimits(5,0,2)
    
# #     for eta in eta_bins:
# #         h_name = "h_eff_eta%s_%s" % (eta, name)
# #         h_eff = fin.Get(h_name)

# #         h_eff.Fit("m_erf")
# #         h_eff.Fit("m_erf","M")
# #         h_eff.SetMaximum(1.3)
# #         h_eff.Draw()
# #         print_canvas("%s" % (h_name), output_path)

# #         fit_name = '%s' % (eta)
# #         fits[fit_name] = []
# #         for i in range(6):
# #             fits[fit_name].append((m_erf.GetParameter(i),m_erf.GetParError(i)))

# #     json.dump(fits, open("results/trigger_object_efficiency_%s.json" % name, "w"))
# #     fin.Close()
        
ROOT.gStyle.SetPaintTextFormat(".3f");
    
for sample in samples:

    measure_trigger_object_efficiency(sample, "_loose", "Muon_looseId")

# def make_overlay_plot(name, h1_name, name1, h2_name, name2):
#     legend = ROOT.TLegend(0.70,0.65,0.85,0.77)
#     legend.SetShadowColor(ROOT.kWhite)
#     legend.SetLineColor(ROOT.kWhite)
#     # legend.SetFillColor(10)
#     # colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kOrange+5, ROOT.kGreen+3]
#     # scale = 1.2

#     if h1_name not in hists:
#         print "%s is no available. Skip" % h1_name
#         return 
#     if h2_name not in hists:
#         print "%s is no available. Skip" % h2_name
#         return 
#     h1 = hists[h1_name]
#     h1.SetLineColor(ROOT.kRed)
#     h1.SetLineWidth(2)
#     h1.SetMarkerStyle(20)
#     h1.SetMarkerColor(ROOT.kRed)
#     # h.GetXaxis().SetTitle("nPV")
#     h1.Draw("hist")
#     legend.AddEntry(h1, name1)

#     h2 = hists[h2_name]
#     h2.SetLineColor(ROOT.kBlack)
#     h2.SetLineWidth(2)
#     h2.SetMarkerStyle(20)
#     h2.SetMarkerColor(ROOT.kBlack)
#     # h.GetXaxis().SetTitle("nPV")
#     h2.Draw("same e")
#     legend.AddEntry(h2, name2)
    
#     legend.Draw()
#     print_canvas(name, output_path)
    
#     ratio_plot = ROOT.TRatioPlot(h2, h1)
#     # ratio_plot.SetH1DrawOpt("hist e")
#     ratio_plot.SetH1DrawOpt("e")
#     ratio_plot.SetH2DrawOpt("hist")
#     ratio_plot.Draw()
#     ratio_plot.SetSeparationMargin(0.03)
#     ratio_plot.GetLowerRefGraph().SetMinimum(0.4)
#     ratio_plot.GetLowerRefGraph().SetMaximum(1.6)
#     # ratio_plot.GetXaxis().SetTitleSize()
#     # SetBottomMargin(2.0)
#     # c1.SetBottomMargin(2.0)
#     # rp->GetLowerRefYaxis()->SetRange(...)
#     # rp->SetH1DrawOpt("E");
#     ratio_plot.GetLowerRefYaxis().SetTitle("Data/MC")
        
#     ratio_plot.GetLowerRefYaxis().SetTitleSize()
#     ratio_plot.GetLowerRefYaxis().SetTitleOffset(1.1)
#     ratio_plot.GetLowerRefYaxis().SetLabelSize(0.035)
#     ratio_plot.GetLowYaxis().SetNdivisions(503)
        
#     ratio_plot.GetLowerRefXaxis().SetTitleSize()
#     ratio_plot.GetLowerRefXaxis().SetTitleOffset()
#     ratio_plot.GetLowerRefXaxis().SetLabelSize(0.035)
    
#     ratio_plot.GetUpperRefYaxis().SetTitle("")
#     ratio_plot.GetUpperRefYaxis().SetTitleSize()
#     ratio_plot.GetUpperRefYaxis().SetTitleOffset()
#     ratio_plot.GetUpperRefYaxis().SetLabelSize(0.035)
        
#     ratio_plot.GetUpperRefXaxis().SetTitleSize()
#     ratio_plot.GetUpperRefXaxis().SetTitleOffset(0)
#     ratio_plot.GetUpperRefXaxis().SetLabelSize(0.035)

#     c1.Update()
#     legend.Draw()
#     print_canvas(name + "_ratio", output_path)

# make_overlay_plot("h_eff_chan0_2018-SM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2018-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2018-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2018-SM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2018-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2018-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2017-SM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2017-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2017-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2017-SM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2017-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2017-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2016GH-SM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2016GH-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2016GH-SM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2016GH-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2016BF-SM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2016BF-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2016BF-SM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2016BF-SM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan0_2018-DM-charm_Bs",
# #                   "h_eff_chan0_Bsmm-2018-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan0_Run2018-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan1_2018-DM-charm_Bs",
# #                   "h_eff_chan1_Bsmm-2018-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan1_Run2018-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan0_2017-DM-charm_Bs",
# #                   "h_eff_chan0_Bsmm-2017-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan0_Run2017-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan1_2017-DM-charm_Bs",
# #                   "h_eff_chan1_Bsmm-2017-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan1_Run2017-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan0_2016GH-DM-charm_Bs",
# #                   "h_eff_chan0_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan0_Run2016GH-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan1_2016GH-DM-charm_Bs",
# #                   "h_eff_chan1_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan1_Run2016GH-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan0_2016BF-DM-charm_Bs",
# #                   "h_eff_chan0_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan0_Run2016BF-DM-charm_Bs", "Data")

# # make_overlay_plot("h_eff_chan1_2016BF-DM-charm_Bs",
# #                   "h_eff_chan1_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
# #                   "h_eff_chan1_Run2016BF-DM-charm_Bs", "Data")

