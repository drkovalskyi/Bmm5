#!/bin/env python
import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv8-516/mc_bjpsik_vs_bsmm"
quick_check = False # quick check uses only the first file for each sample

mm_cuts = [
    # skimming requirements
    "mm_mu1_index>=0", "mm_mu2_index>=0",
    "abs(Muon_eta[mm_mu1_index])<1.4", "Muon_pt[mm_mu1_index]>4",
    "abs(Muon_eta[mm_mu2_index])<1.4", "Muon_pt[mm_mu2_index]>4",
    "abs(mm_kin_mass-5.4)<0.5", "mm_kin_sl3d>4", "mm_kin_vtx_chi2dof<5",
    # muon id
    "Muon_softMva[mm_mu1_index]>0.45 ", " Muon_softMva[mm_mu2_index]>0.45",
    # gen matching
    "mm_gen_pdgId!=0",
    # tuning
    # "mm_kin_pt>10.*5.3/3.1", "mm_kin_pt<20*5.3/3.1",
    "mm_kin_alpha<0.2",
    
]
mm_selection = "&&".join(mm_cuts)

bkmm_cuts = [
    # skimming requirements
    "mm_mu1_index[bkmm_mm_index]>=0", "mm_mu2_index[bkmm_mm_index]>=0",
    "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4",
    "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4",
    "mm_kin_sl3d[bkmm_mm_index]>4", "mm_kin_vtx_chi2dof[bkmm_mm_index]<5",
    "abs(bkmm_jpsimc_mass-5.4)<0.5", 
    "bkmm_jpsimc_vtx_chi2dof<5", # need to get rid of it
    "bkmm_jpsimc_alpha<0.2",
    # muon id
    "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 ", " Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45",
    # gen matching
    "bkmm_gen_pdgId!=0",
    # kaon
    "bkmm_kaon_pt<1.5",
    # tuning
    # "mm_kin_pt[bkmm_mm_index]>10", "mm_kin_pt[bkmm_mm_index]<20",
]
bkmm_selection = "&&".join(bkmm_cuts)


# "bkmm_kaon_sdxy_bs>5&&bkmm_kaon_pt>1&&abs(bkmm_kaon_eta)<1.4" 

samples = {
    'BsToMuMu':{
        'files':[
        ],
        'color':ROOT.kBlue,
        'type':'bmm',
        'legend':'B #rightarrow #mu#mu',
    },
    'BuToJpsiK':{
        'files':[
        ],
        'color':ROOT.kRed,
        'type':'bjpsik',
        'legend':'B #rightarrow J/#psiK'
    }
}

n_max = 100
bmm_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/"
for i,f in enumerate(subprocess.check_output("find %s/ -type f -name '*.root'" % (bmm_path), shell=True).split("\n")):
    if i >= n_max: break
    if f: samples['BsToMuMu']['files'].append(f)
    
bjpsik_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/"
for i,f in enumerate(subprocess.check_output("find %s/ -type f -name '*.root'" % (bjpsik_path), shell=True).split("\n")):
    if i >= n_max: break
    if f: samples['BuToJpsiK']['files'].append(f)


# read list of files from a file instead

# with open("Charmonium+Run2018D-PromptReco-v2+MINIAOD.list") as f:
#     samples['Charmonium_Run2018D_PromptReco']['files'] = []
#     for file in f:
#         samples['Charmonium_Run2018D_PromptReco']['files'].append(file.strip())
        
# Load data

mask = '\w\d\_\.'
for name,sample in samples.items():
    if re.search("[^%s]"%mask,name):
        raise Exception("Illigal symbol used for sample name. Allowed %s" % mask) 
    chain = ROOT.TChain("Events")
    for entry in sample['files']:
        chain.Add(entry)
        if quick_check: break
    sample['events'] = chain
    
print "Number of events:"
for name,sample in samples.items():
    sample['nAll'] = sample['events'].GetEntries()
    print "\t%s: \t%u" % (name,sample['nAll'])

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def plot_generic_1D(selections,hist_title,file_name,vars,nbins=100,xmin=0,xmax=100):
    c1 = TCanvas("c1", "c1", 800, 800)
    max_value = 0
    for name,sample in samples.items():
        selection = selections[sample['type']]
        var = vars[sample['type']]
        hist = ROOT.TH1D("hist",hist_title,nbins,xmin,xmax)
        hist.SetLineColor(sample['color'])
        hist.SetLineWidth(2)
        sample['nSelected'] = sample['events'].Draw("%s>>hist"%var,selection)
        # print_canvas("%s_%s"%(file_name,name), output_path)
        if hist.GetEntries()>0:
            hist.Scale(1/hist.GetEntries()) # normalize
        if max_value < hist.GetMaximum(): max_value = hist.GetMaximum()
        hist.SetDirectory(0)
        sample['hist'] = hist
    
    legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
    legend.SetFillStyle(0)
    legend.SetLineWidth(0)

    first_plot = True
    for name,sample in samples.items():
        sample['hist'].SetMinimum(0)
        sample['hist'].SetMaximum(max_value*1.2)
        legend.AddEntry(sample['hist'], sample['legend'])
        if first_plot:
            sample['hist'].Draw("hist")
            first_plot = False
        else:
            sample['hist'].Draw("hist same")
    legend.Draw()
    print_canvas(file_name, output_path)
    print "Number of selected events:"
    for name,sample in samples.items():
        print "\t%s: \t%u out of %u" % (name,sample['nSelected'],sample['nAll'])

def integrate(source_hist, destination_hist, left_to_right = True):
    assert(source_hist.GetNbinsX() == destination_hist.GetNbinsX())
    assert(source_hist.Integral(0,-1) > 0)
    for bin in range(source_hist.GetNbinsX() + 1):
        if left_to_right:
            destination_hist.SetBinContent(bin, source_hist.Integral(bin,-1)/source_hist.Integral(0,-1))
        else:
            destination_hist.SetBinContent(bin, 1 - source_hist.Integral(bin,-1)/source_hist.Integral(0,-1))

def plot_integral_1D(selections, hist_title, file_name, vars, 
                     nbins=100, xmin=0, xmax=100, 
                     left_to_right = True):
    """Make disitributions of the target variables and make intergral
       plots. By default the integral is from x till xmax.
    """
    c1 = TCanvas("c1", "c1", 800, 800)
    max_value = 0
    for name,sample in samples.items():
        selection = selections[sample['type']]
        var = vars[sample['type']]
        hist = ROOT.TH1D("hist", hist_title, nbins, xmin, xmax)
        hist_integrated = ROOT.TH1D("hist_integrated", hist_title, nbins, xmin, xmax)
        hist_integrated.SetLineColor(sample['color'])
        hist_integrated.SetLineWidth(2)
        sample['nSelected'] = sample['events'].Draw("%s>>hist" % var, selection)
        integrate(hist, hist_integrated, left_to_right)
        # print_canvas("%s_%s"%(file_name,name), output_path)
        if max_value < hist_integrated.GetMaximum(): max_value = hist_integrated.GetMaximum()
        hist_integrated.SetDirectory(0)
        hist.SetDirectory(0)
        sample['hist'] = hist_integrated
    
    legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
    legend.SetFillStyle(0)
    legend.SetLineWidth(0)

    first_plot = True
    for name,sample in samples.items():
        sample['hist'].SetMinimum(0)
        sample['hist'].SetMaximum(max_value * 1.2)
        legend.AddEntry(sample['hist'], sample['legend'])
        if first_plot:
            sample['hist'].Draw("hist")
            first_plot = False
        else:
            sample['hist'].Draw("hist same")
    legend.Draw()
    print_canvas(file_name, output_path)
    print "Number of selected events:"
    for name,sample in samples.items():
        print "\t%s: \t%u out of %u" % (name,sample['nSelected'],sample['nAll'])

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

selections = {
    'bmm':mm_selection,
    'bjpsik':bkmm_selection,
}

# Variables used in MVA
# * bkmm_jpsimc_alpha
# * bkmm_jpsimc_cosAlphaXY
# * bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr
# * bkmm_jpsimc_pvip
# * bkmm_bmm_m1iso
# * bkmm_bmm_m2iso
# * bkmm_bmm_iso
# * bkmm_bmm_nBMTrks
# * bkmm_bmm_otherVtxMaxProb1
# * bkmm_bmm_otherVtxMaxProb2
# * mm_kin_vtx_chi2dof
# * mm_kin_sl3d
#
# Should also consider
# * bkmm_jpsimc_vtx_chi2dof
# * bkmm_jpsimc_sl3d

# plot_generic_1D(selections, "#mu#mu;P_{T}, [GeV]", "01_mm_pt",
#                 {'bmm':'mm_kin_pt', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'}, 100, 0, 100)

plot_generic_1D(selections, "#mu#mu vertex displacement significance;#sigma", "02_mm_sl3d",
                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]'}, 100, 0, 100)
plot_generic_1D(selections, "#mu#muK vertex displacement significance;#sigma", "02_kmm_sl3d",
                {'bmm':'mm_kin_sl3d', 'bjpsik':'bkmm_jpsimc_sl3d'}, 100, 0, 100)
# plot_generic_1D(selections, "#mu#mu vertex displacement significance (#mu#muK is scaled by 1.5);#sigma", "02_mm_sl3d_scaled",
#                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]*1.5'}, 100, 0, 100)
scale = 1.60
plot_generic_1D({'bmm':mm_selection + "&&mm_kin_sl3d>4*%s" % scale, 'bjpsik':bkmm_selection},
                "#mu#mu vertex displacement significance (#mu#muK is scaled by %s);#sigma" % scale, 
                "02_mm_sl3d_scaled_and_matched_selection",
                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]*%s' % scale}, 100, 0, 100)

plot_generic_1D(selections, "Pointing angle 3D;#alpha_{3D}", "03_alpha",
                {'bmm':'mm_kin_alpha', 'bjpsik':'bkmm_jpsimc_alpha'}, 100, 0, 0.2)
plot_generic_1D(selections, "Pointing angle BS;#alpha_{BS}", "03_alphaBS",
                {'bmm':'mm_kin_alphaBS', 'bjpsik':'bkmm_jpsimc_alphaBS'}, 110, 0, 0.2)

plot_generic_1D(selections, "Impact parameter significance", "04_spvip",
                {'bmm':'mm_kin_pvip/mm_kin_pvipErr', 'bjpsik':'bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr'},
                100, 0, 5)
plot_generic_1D(selections, "Impact parameter 3D", "04_pvip",
                {'bmm':'mm_kin_pvip', 'bjpsik':'bkmm_jpsimc_pvip'}, 100, 0, 0.02)

plot_generic_1D(selections, "#mu#mu isolation", "05_iso",
                {'bmm':'mm_iso', 'bjpsik':'bkmm_bmm_iso'}, 120, 0, 1.2)
plot_generic_1D(selections, "#mu1 isolation", "05_m1iso",
                {'bmm':'mm_m1iso', 'bjpsik':'bkmm_bmm_m1iso'}, 120, 0, 1.2)
plot_generic_1D(selections, "#mu2 isolation", "05_m2iso",
                {'bmm':'mm_m2iso', 'bjpsik':'bkmm_bmm_m2iso'}, 120, 0, 1.2)

plot_generic_1D(selections, "#chi/nDof for #mu#mu vertex", "06_mm_chi2dof",
                {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'mm_kin_vtx_chi2dof[bkmm_mm_index]'}, 100, 0, 5)
plot_generic_1D(selections, "#chi/nDof for #mu#muK vertex", "06_mmK_chi2dof",
                {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'bkmm_jpsimc_vtx_chi2dof'}, 100, 0, 5)

plot_generic_1D(selections, "nBMTrks", "07_nBMTrks",
                {'bmm':'mm_nBMTrks','bjpsik':'min(bkmm_bmm_nBMTrks,9)'}, 10, 0, 10)
plot_generic_1D(selections, "otherVtxMaxProb1", "07_otherVtxMaxProb1",
                {'bmm':'mm_otherVtxMaxProb1', 'bjpsik':'bkmm_bmm_otherVtxMaxProb1'}, 120, 0, 1.2)
plot_generic_1D(selections, "otherVtxMaxProb2", "07_otherVtxMaxProb2",
                {'bmm':'mm_otherVtxMaxProb2', 'bjpsik':'bkmm_bmm_otherVtxMaxProb2'}, 120, 0, 1.2)

# plot_generic_1D(selections, "BDT Matched", "09_bdt_matched",
#                 {'bmm':'mm_bdt', 'bjpsik':'bkmm_bmm_bdt'}, 100, -1.5, 1.5)
# plot_generic_1D(selections, "BDT Raw", "09_bdt_raw",
#                 {'bmm':'mm_bdt', 'bjpsik':'mm_bdt[bkmm_mm_index]'}, 100, -1.5, 1.5)

plot_generic_1D(selections,"MVA Matched", "09_mva_matched",
                {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)
plot_generic_1D(selections,"MVA Matched", "09_mva_matched_zoomed",
                {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 101, 0.9, 1.01)
plot_generic_1D(selections,"MVA Raw", "09_mva_raw",
                {'bmm':'mm_mva', 'bjpsik':'mm_mva[bkmm_mm_index]'}, 110, 0, 1.1)

# plot_integral_1D(selections,"MVA Matched Efficiency", "09_mva_matched_eff",
#                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)

