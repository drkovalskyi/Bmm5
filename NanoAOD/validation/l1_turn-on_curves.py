import os
import sys
import json
import re
import ROOT
import tdrstyle
from array import array
import subprocess

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)
ROOT.gPad.SetGrid()
# ROOT.gPad.SetLogx()

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/"
output_path = "/eos/home-d/dmytro/www/plots/tmp/Run3/l1_turn-on_curves/";

datasets =  [
    data_path + '/ParkingDoubleMuonLowMass0+Run2024C-PromptReco-v1+MINIAOD/11*.root',
    # data_path3 + '/ParkingDoubleMuonLowMass0+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass0+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass0+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass0+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass0+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass1+Run2024C-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass1+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass1+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass1+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass1+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass1+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass2+Run2024C-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass2+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass2+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass2+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass2+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass2+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass3+Run2024C-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass3+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass3+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass3+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass3+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass3+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass4+Run2024C-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass4+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass4+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass4+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass4+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass4+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass5+Run2024C-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass5+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass5+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass5+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass5+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass5+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass6+Run2024C-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass6+Run2024D-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass6+Run2024E-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass6+Run2024E-PromptReco-v2+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass6+Run2024F-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass6+Run2024G-PromptReco-v1+MINIAOD/',
    # data_path3 + '/ParkingDoubleMuonLowMass7+Run2024C-PromptReco-v1+MINIAOD/',
]

def get_data(tree="Events"):
    chain = ROOT.TChain(tree)
    n_files = 0 
    for entry in datasets:
        n_files += chain.Add(entry)
    if n_files == 0:
        print("No files found for " + name)
        return None
    return chain

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))


def trigger_efficiency():
    chain  = get_data(sample)
    if not chain: return

    h_name = "h_all"
    h_all = ROOT.TH1F(h_name, "", 36, 2, 20)
    h_all.Sumw2()
    selection = "nMuon==3&&nmm>0 && abs(mm_kin_mass[0]-3.1)<0.1 && mm_mu1_index[0]==0 " \
        "&& mm_mu2_index[0]==1 && abs(mm_mu1_eta[0])>1.0 && abs(mm_mu2_eta[0])>1.0 && abs(Muon_eta[2])<0.8"

    chain.Draw(f"Muon_pt[2]>>{h_name}", selection)

    print_canvas(h_name, output_path)
    h_all.SetDirectory(0)

    h_name = "h_trig"
    h_trig = h_all.Clone(h_name)
    chain.Draw(f"Muon_pt[2]>>{h_name}", selection + "&& MuonId_l1_quality >= 12" ) # single muon quality
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


h_name = "h_all"
h_all = rdf2.Define(rdf_name + "_pt", f"Muon_pt[{rdf_name}]").Histo1D((h_name, "", 32, 2, 10), f"{rdf_name}_pt")
        # h_all.Sumw2()
        histos[h_name] = h_all
        print(f"booked {h_name}")

        h_name = f"h_trig_eta{eta_min:+.1f}to{eta_max:+.1f}_{suffix}_{sample_name}"
        selection += f"&& {test_selection}"
        rdf_name += "_trig"
        rdf3 = rdf.Define(rdf_name, selection)
        h_trig = rdf3.Define(rdf_name + "_pt", f"Muon_pt[{rdf_name}]").Histo1D((h_name, "", 32, 2, 10), rdf_name + "_pt")
        # h_trig.Sumw2()
        histos[h_name] = h_trig
        print(f"booked {h_name}")

        category_name = f"eta{eta_min:+.1f}to{eta_max:+.1f}_{suffix}"
        h_name = f"h_eff_{category_name}_{sample_name}"
        eff_histos[h_name] = (h_trig, h_all)
        if category_name not in histos_by_category:
            histos_by_category[category_name] = dict()
        histos_by_category[category_name][sample_name] = h_name


def make_trigger_efficiency_plots(h_name, h_trig, h_all):
    
    h_eff = h_trig.Clone(h_name)
    h_eff.Divide(h_trig.GetPtr(), h_all.GetPtr(), 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.SetMarkerStyle(20)
    h_eff.SetMarkerSize(2)
    h_eff.SetMaximum(1.1)
    h_eff.Draw("e0")
    histos[h_name] = h_eff
    print_canvas("%s" % (h_name), output_path)


def make_overlay_plots(sample_name1, legend_name1, sample_name2, legend_name2, ratio_label = "Data/MC"):
    for category_name, histograms in histos_by_category.items():
        if sample_name1 not in histograms:
            continue
        else:
            h1_name = histograms[sample_name1]

        if sample_name2 not in histograms:
            continue
        else:
            h2_name = histograms[sample_name2]
        
        legend = ROOT.TLegend(0.50,0.65,0.85,0.77)
        legend.SetShadowColor(ROOT.kWhite)
        legend.SetLineColor(ROOT.kWhite)
        # legend.SetFillColor(10)
        # colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kOrange+5, ROOT.kGreen+3]
        # scale = 1.2

        h1 = histos[h1_name]
        h1.SetLineColor(ROOT.kRed)
        h1.SetLineWidth(2)
        h1.SetMarkerStyle(20)
        h1.SetMarkerColor(ROOT.kRed)
        # h.GetXaxis().SetTitle("nPV")
        h1.Draw("hist")
        legend.AddEntry(h1, legend_name1)

        h2 = histos[h2_name]
        h2.SetLineColor(ROOT.kBlack)
        h2.SetLineWidth(2)
        h2.SetMarkerStyle(20)
        h2.SetMarkerColor(ROOT.kBlack)
        # h.GetXaxis().SetTitle("nPV")
        h2.Draw("same e")
        legend.AddEntry(h2, legend_name2)

        legend.Draw()
        print_canvas(f"overlay_{category_name}_{sample_name1}_{sample_name2}", output_path)

        ratio_plot = ROOT.TRatioPlot(h2, h1)
        # ratio_plot.SetH1DrawOpt("hist e")
        ratio_plot.SetH1DrawOpt("e")
        ratio_plot.SetH2DrawOpt("hist")
        ratio_plot.Draw()
        ratio_plot.SetSeparationMargin(0.03)
        ratio_plot.GetLowerRefGraph().SetMinimum(0.7)
        ratio_plot.GetLowerRefGraph().SetMaximum(1.3)
        # ratio_plot.GetXaxis().SetTitleSize()
        # SetBottomMargin(2.0)
        # c1.SetBottomMargin(2.0)
        # rp->GetLowerRefYaxis()->SetRange(...)
        # rp->SetH1DrawOpt("E");
        ratio_plot.GetLowerRefYaxis().SetTitle(f"{ratio_label}")

        ratio_plot.GetLowerRefYaxis().SetTitleSize()
        ratio_plot.GetLowerRefYaxis().SetTitleOffset(1.1)
        ratio_plot.GetLowerRefYaxis().SetLabelSize(0.035)
        ratio_plot.GetLowYaxis().SetNdivisions(503)

        ratio_plot.GetLowerRefXaxis().SetTitleSize()
        ratio_plot.GetLowerRefXaxis().SetTitleOffset()
        ratio_plot.GetLowerRefXaxis().SetLabelSize(0.035)

        ratio_plot.GetUpperRefYaxis().SetTitle("")
        ratio_plot.GetUpperRefYaxis().SetTitleSize()
        ratio_plot.GetUpperRefYaxis().SetTitleOffset()
        ratio_plot.GetUpperRefYaxis().SetLabelSize(0.035)

        ratio_plot.GetUpperRefXaxis().SetTitleSize()
        ratio_plot.GetUpperRefXaxis().SetTitleOffset(0)
        ratio_plot.GetUpperRefXaxis().SetLabelSize(0.035)

        c1.Update()
        legend.Draw()
        print_canvas(f"overlay_{category_name}_{sample_name1}_{sample_name2}_ratio", output_path)
    
        
ROOT.gStyle.SetPaintTextFormat(".3f");

eta_step = 0.4
for sample_name in samples:
    print(f"processing {sample_name}")
    
    single_muon_rdf = get_data_single_muon(sample_name)
    jpsi_l1_rdf = get_data_jpsi_l1(sample_name)
    jpsi_hlt_rdf = get_data_jpsi_hlt(sample_name)
    
    print(f"RDFs are defined")

    # if single_muon_rdf:
    #     book_histograms(single_muon_rdf, sample_name, "base",
    #                     "Muon_looseId&&Muon_isTracker&&Muon_isGlobal&&Muon_nStations>=2&&abs(Muon_dxybs)<0.1",
    #                     "MuonId_l1_quality >= 12", eta_step)
    #     book_histograms(single_muon_rdf, sample_name, "base_medium",
    #                     "Muon_mediumId&&Muon_isTracker&&Muon_isGlobal&&Muon_nStations>=2&&abs(Muon_dxybs)<0.1",
    #                     "MuonId_l1_quality >= 12", eta_step)
    #     book_histograms(single_muon_rdf, sample_name, "loose",
    #                     "Muon_looseId",
    #                     "MuonId_l1_quality >= 12", eta_step)
    #     book_histograms(single_muon_rdf, sample_name, "medium",
    #                     "Muon_mediumId",
    #                     "MuonId_l1_quality >= 12", eta_step)
    #     book_histograms(single_muon_rdf, sample_name, "tight",
    #                     "Muon_tightId",
    #                     "MuonId_l1_quality >= 12", eta_step)

    # if jpsi_l1_rdf:
    #     book_histograms(jpsi_l1_rdf, sample_name, "jpsi_loose",
    #                     "Muon_looseId",
    #                     "MuonId_l1_quality >= 12")
        
    if jpsi_hlt_rdf:
        book_histograms(jpsi_hlt_rdf, sample_name, "jpsi_probe_hlt_eff_looseid",
                        "L1_jpsi_probe",
                        "HLT_fired", eta_step)

ROOT.RDF.RunGraphs(histos.values())
print("done with processing")
# print(histos.keys())

# primary histograms
for h_name, histo in histos.items():
    histo.Draw()
    print_canvas(h_name, output_path)

# derived histograms
for h_name, (h_trig, h_all) in eff_histos.items():
    make_trigger_efficiency_plots(h_name, h_trig, h_all)

make_overlay_plots("B2JpsiK", "MC", "ParkingDoubleEl", "Data")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2022C", "Data-Run2022C")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2022D", "Data-Run2022D")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2022E", "Data-Run2022E")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2022F", "Data-Run2022F")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2022", "Data-Run2022")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2023", "Data-Run2023")
make_overlay_plots("BdToJpsiKstar", "MC-Summer22EE", "ParkingDoubleMu_2024", "Data-Run2024")

