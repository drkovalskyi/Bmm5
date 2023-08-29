#!/bin/env python
import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

# import tdrstyle
# Set the TDR style
# tdrstyle.setTDRStyle()
# ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/dmm/mc_dzpipi_production"
quick_check = False # quick check uses only the first file for each sample

mc_selection = (
    ## basedline
    'dstar_hh_index>=0 && hh_had1_pdgId*hh_had2_pdgId==-211*211 && '
    'abs(dstar_gen_pdgId)==413 && abs(hh_gen_had1_pdgId)==211 && abs(hh_gen_had2_pdgId)==211 && '
    'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
    ## mass cuts
    '&& dstar_dm_pv>0.140 && dstar_dm_pv<0.155'
    '&& hh_kin_mass[dstar_hh_index]>1.81 && hh_kin_mass[dstar_hh_index]<1.94'
    ## other
    # 'hh_kin_vtx_prob[dstar_hh_index]>0.01'
    # 'dstar_pv_with_pion_prob>0.1'
    # 'hh_kin_sl3d[dstar_hh_index]>3'
    # 'hh_kin_alpha[dstar_hh_index]<0.1'
)

bkg_selection = (
    ## basedline
    'dstar_hh_index>=0 && hh_had1_pdgId*hh_had2_pdgId==-211*211 && '
    'hh_had1_pt[dstar_hh_index]>4 && hh_had2_pt[dstar_hh_index]>4'
    ## mass cuts
    '&& dstar_dm_pv>0.140 && dstar_dm_pv<0.155'
    '&& hh_kin_mass[dstar_hh_index]>1.90'
)

bkg_dmm_selection = (
    ## basedline
    'dstar_mm_index>=0 && '
    'mm_mu1_pt[dstar_mm_index]>4 && mm_mu2_pt[dstar_mm_index]>4'
    ## mass cuts
    '&& dstar_dm_pv>0.140 && dstar_dm_pv<0.155'
    '&& mm_kin_mass[dstar_mm_index]>1.90'
)


samples = {
    'dzpipi_b':{
        'files':[
        ],
        'color':ROOT.kBlack,
        'fill':ROOT.kYellow,
        'type':'dzpipi_b',
        'legend':'D^{0} #rightarrow #pi#pi (b-quark)',
    },
    'dzpipi_c':{
        'files':[
        ],
        'color':ROOT.kBlue,
        'type':'dzpipi_c',
        'legend':'D^{0} #rightarrow #pi#pi (c-quark)',
    },
    'data_bkg':{
        'files':[
        ],
        'color':ROOT.kBlack,
        'type':'data_bkg',
        'legend':'Data hh background',
    },
    'data_dzmm':{
        'files':[
        ],
        'color':ROOT.kRed,
        # 'style':ROOT.kDashed,
        'type':'data_dzmm',
        'legend':'Data #mu#mu background',
    },
}

n_max = 100
path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5+MINIAODSIM/"
for i,f in enumerate(subprocess.check_output("find %s/ -type f -name '*.root'" % (path), shell=True, encoding='utf8').split("\n")):
    # if i >= n_max: break
    if f: samples['dzpipi_b']['files'].append(f)

directories = [
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022C-PromptReco-v1+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022D-PromptReco-v1+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022D-PromptReco-v2+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022D-PromptReco-v3+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022E-PromptReco-v1+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022F-PromptReco-v1+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzpipi/ZeroBias+Run2022G-PromptReco-v1+MINIAOD/",
]
for path in directories:
    for i,f in enumerate(subprocess.check_output("find %s/ -type f -name '*.root'" % (path), shell=True, encoding='utf8').split("\n")):
        # if i >= n_max: break
        if f: samples['data_bkg']['files'].append(f)

directories = [
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzmm/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/",
    "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/dzmm/ParkingDoubleMuonLowMass1+Run2022C-PromptReco-v1+MINIAOD/",
]
for path in directories:
    for i,f in enumerate(subprocess.check_output("find %s/ -type f -name '*.root'" % (path), shell=True, encoding='utf8').split("\n")):
        # if i >= n_max: break
        if f: samples['data_dzmm']['files'].append(f)
        
# bjpsik_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/"
# for i,f in enumerate(subprocess.check_output("find %s/ -type f -name '*.root'" % (bjpsik_path), shell=True).split("\n")):
#     if i >= n_max: break
#     if f: samples['BuToJpsiK']['files'].append(f)
samples['dzpipi_c']['files'].append('/eos/cms/store/group/phys_muon/dmytro/tmp/DstarToD0Pi_D0To2Pi_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen.root')

# read list of files from a file instead

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
    
print("Number of events:")
for name,sample in samples.items():
    sample['nAll'] = sample['events'].GetEntries()
    print("\t%s: \t%u" % (name,sample['nAll']))

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def plot_generic_1D(selections,hist_title,file_name,vars,nbins=100,xmin=0,xmax=100,show_overflow=False):
    c1 = TCanvas("c1", "c1", 800, 800)
    max_value = 0
    for name,sample in samples.items():
        selection = selections[sample['type']]
        if sample['type'] not in vars:
            continue
        var = vars[sample['type']]
        hist = ROOT.TH1D("hist",hist_title,nbins,xmin,xmax)
        if show_overflow:
            hist.GetXaxis().SetRangeUser(xmin - 1e-5, xmax + 1e-5)
        hist.SetLineColor(sample['color'])
        if 'fill' in sample:
            hist.SetFillColor(sample['fill'])
            hist.SetLineWidth(1)
        else:
            hist.SetLineWidth(3)
            if 'style' in sample:
                hist.SetLineStyle(sample['style'])
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
        if sample['type'] not in vars:
            continue
        sample['hist'].SetMinimum(0)
        sample['hist'].SetMaximum(max_value * 1.3)
        legend.AddEntry(sample['hist'], sample['legend'])
        if first_plot:
            sample['hist'].Draw("hist")
            first_plot = False
        else:
            sample['hist'].Draw("hist same")
    legend.Draw()
    print_canvas(file_name, output_path)
    print("Number of selected events:")
    for name,sample in samples.items():
        if sample['type'] not in vars:
            continue
        print("\t%s: \t%u out of %u" % (name,sample['nSelected'],sample['nAll']))

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
        sample['hist'].SetMaximum(max_value * 1.3)
        legend.AddEntry(sample['hist'], sample['legend'])
        if first_plot:
            sample['hist'].Draw("hist")
            first_plot = False
        else:
            sample['hist'].Draw("hist same")
    legend.Draw()
    print_canvas(file_name, output_path)
    print("Number of selected events:")
    for name,sample in samples.items():
        print("\t%s: \t%u out of %u" % (name,sample['nSelected'],sample['nAll']))

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

selections = {
    'dzpipi_c':mc_selection,
    'dzpipi_b':(
        mc_selection +
        " && (abs(dstar_gen_mpdgId)==511 || abs(dstar_gen_mpdgId)==521 || abs(dstar_gen_mpdgId)==531)"
    ),
    'data_bkg':bkg_selection,
    'data_dzmm':bkg_dmm_selection
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

plot_generic_1D(selections, "#Delta m with IP constrained soft pion;#Delta m, [GeV]", "01_dm",
                {'dzpipi_c':'dstar_dm_pv', 'dzpipi_b':'dstar_dm_pv', 'data_bkg':'dstar_dm_pv'},
                30, 0.140, 0.155)

plot_generic_1D(selections, "#Delta m with free soft pion;#Delta m, [GeV]", "01_dm_free",
                {'dzpipi_c':'dstar_dm_free', 'dzpipi_b':'dstar_dm_free', 'data_bkg':'dstar_dm_free'},
                30, 0.140, 0.155)

plot_generic_1D(selections, "Vertex contrained D^{0} mass;m, [GeV]", "01_d0mass",
                {'dzpipi_c':'hh_kin_mass[dstar_hh_index]',
                 'dzpipi_b':'hh_kin_mass[dstar_hh_index]',
                 # 'data_bkg':'hh_kin_mass[dstar_hh_index]'
                 },
                30, 1.80, 1.95)

plot_generic_1D(selections, "Pointing angle 3D;#alpha_{3D}", "02_alpha",
                {'dzpipi_c':'hh_kin_alpha[dstar_hh_index]',
                 'dzpipi_b':'hh_kin_alpha[dstar_hh_index]',
                 'data_bkg':'hh_kin_alpha[dstar_hh_index]',
                 'data_dzmm':'mm_kin_alpha[dstar_mm_index]'},
                64, 0., 3.2)

plot_generic_1D(selections, "#pi#pi vertex displacement significance;#sigma", "02_hh_sl3d",
                {'dzpipi_c':'hh_kin_sl3d[dstar_hh_index]',
                 'dzpipi_b':'hh_kin_sl3d[dstar_hh_index]',
                 'data_bkg':'hh_kin_sl3d[dstar_hh_index]',
                 'data_dzmm':'mm_kin_sl3d[dstar_mm_index]'},
                50, 0, 50, True)

plot_generic_1D(selections, "Leading pion ;p_{T}, [GeV]", "03_pion1_pt",
                {'dzpipi_c':'hh_had1_pt[dstar_hh_index]',
                 'dzpipi_b':'hh_had1_pt[dstar_hh_index]',
                 'data_bkg':'hh_had1_pt[dstar_hh_index]'},
                40, 0, 20, True)

plot_generic_1D(selections, "Sub-leading pion ;p_{T}, [GeV]", "03_pion2_pt",
                {'dzpipi_c':'hh_had2_pt[dstar_hh_index]',
                 'dzpipi_b':'hh_had2_pt[dstar_hh_index]',
                 'data_bkg':'hh_had2_pt[dstar_hh_index]'},
                40, 0, 20, True)

plot_generic_1D(selections, "Dstar PV with soft pion vertex probability", "04_dstar_pv_with_pion_prob",
                {'dzpipi_c':'dstar_pv_with_pion_prob',
                 'dzpipi_b':'dstar_pv_with_pion_prob',
                 'data_bkg':'dstar_pv_with_pion_prob',
                 'data_dzmm':'dstar_pv_with_pion_prob'},
                50, 0, 1.0)

plot_generic_1D(selections, "D^{0} vertex probability", "04_d0_vtx_prob",
                {'dzpipi_c':'hh_kin_vtx_prob[dstar_hh_index]',
                 'dzpipi_b':'hh_kin_vtx_prob[dstar_hh_index]',
                 'data_bkg':'hh_kin_vtx_prob[dstar_hh_index]',
                 'data_dzmm':'mm_kin_vtx_prob[dstar_mm_index]'},
                50, 0, 1.0)

# plot_generic_1D(selections, "Pointing angle 3D;#alpha_{3D}", "03_alpha",
#                 {'bmm':'mm_kin_alpha', 'bjpsik':'bkmm_jpsimc_alpha'}, 100, 0, 0.2)
# plot_generic_1D(selections, "Pointing angle BS;#alpha_{BS}", "03_alphaBS",
#                 {'bmm':'mm_kin_alphaBS', 'bjpsik':'bkmm_jpsimc_alphaBS'}, 110, 0, 0.2)

# plot_generic_1D(selections, "Impact parameter significance", "04_spvip",
#                 {'bmm':'mm_kin_pvip/mm_kin_pvipErr', 'bjpsik':'bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr'},
#                 100, 0, 5)
# plot_generic_1D(selections, "Impact parameter 3D", "04_pvip",
#                 {'bmm':'mm_kin_pvip', 'bjpsik':'bkmm_jpsimc_pvip'}, 100, 0, 0.02)

# plot_generic_1D(selections, "#mu#mu isolation", "05_iso",
#                 {'bmm':'mm_iso', 'bjpsik':'bkmm_bmm_iso'}, 120, 0, 1.2)
# plot_generic_1D(selections, "#mu1 isolation", "05_m1iso",
#                 {'bmm':'mm_m1iso', 'bjpsik':'bkmm_bmm_m1iso'}, 120, 0, 1.2)
# plot_generic_1D(selections, "#mu2 isolation", "05_m2iso",
#                 {'bmm':'mm_m2iso', 'bjpsik':'bkmm_bmm_m2iso'}, 120, 0, 1.2)

# plot_generic_1D(selections, "#chi/nDof for #mu#mu vertex", "06_mm_chi2dof",
#                 {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'mm_kin_vtx_chi2dof[bkmm_mm_index]'}, 100, 0, 5)
# plot_generic_1D(selections, "#chi/nDof for #mu#muK vertex", "06_mmK_chi2dof",
#                 {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'bkmm_jpsimc_vtx_chi2dof'}, 100, 0, 5)

# plot_generic_1D(selections, "nBMTrks", "07_nBMTrks",
#                 {'bmm':'mm_nBMTrks','bjpsik':'min(bkmm_bmm_nBMTrks,9)'}, 10, 0, 10)
# plot_generic_1D(selections, "otherVtxMaxProb1", "07_otherVtxMaxProb1",
#                 {'bmm':'mm_otherVtxMaxProb1', 'bjpsik':'bkmm_bmm_otherVtxMaxProb1'}, 120, 0, 1.2)
# plot_generic_1D(selections, "otherVtxMaxProb2", "07_otherVtxMaxProb2",
#                 {'bmm':'mm_otherVtxMaxProb2', 'bjpsik':'bkmm_bmm_otherVtxMaxProb2'}, 120, 0, 1.2)

# # plot_generic_1D(selections, "BDT Matched", "09_bdt_matched",
# #                 {'bmm':'mm_bdt', 'bjpsik':'bkmm_bmm_bdt'}, 100, -1.5, 1.5)
# # plot_generic_1D(selections, "BDT Raw", "09_bdt_raw",
# #                 {'bmm':'mm_bdt', 'bjpsik':'mm_bdt[bkmm_mm_index]'}, 100, -1.5, 1.5)

# plot_generic_1D(selections,"MVA Matched", "09_mva_matched",
#                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)
# plot_generic_1D(selections,"MVA Matched", "09_mva_matched_zoomed",
#                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 101, 0.9, 1.01)
# plot_generic_1D(selections,"MVA Raw", "09_mva_raw",
#                 {'bmm':'mm_mva', 'bjpsik':'mm_mva[bkmm_mm_index]'}, 110, 0, 1.1)

# plot_integral_1D(selections,"MVA Matched Efficiency", "09_mva_matched_eff",
#                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)

