"""Trigger efficiency study

Trigger efficiency plots for Ks->mumu events.

The input data should be in Bmm5 NanoAOD format.
"""

import sys
import ROOT
if not hasattr(ROOT.RDataFrame, "DefaultValueFor"):
    sys.exit(f"Error: ROOT.RDataFrame.DefaultValueFor is not available in ROOT {ROOT.gROOT.GetVersion()}. Please use ROOT version 6.34 or later.")
from collections import defaultdict
from pprint import pprint
import os
import re
import glob
import json
import math
import tdrstyle
import numpy as np

ROOT.ROOT.EnableImplicitMT()
# max_files = 10
max_files = 999999
path   = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/"
path2  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/532/mm/"
output_path = "/eos/home-d/dmytro/www/plots/tmp/kmm_mc_trigger_efficiency"


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))


def load_data(file_patterns):
    chain = ROOT.TChain("Events")
    files = []
    for file_pattern in file_patterns:
        files.extend(glob.glob(file_pattern))
    print("Number of files:", len(files))
    if len(files) == 0:
        return None
    for file in files[:max_files]:
        chain.Add(file)
    n_files = chain.GetListOfFiles().GetEntries()
    print("Number of files in the chain:", n_files)
    if n_files == 0:
        return None
    return chain


def process_sample(file_patterns, selection=None, nbins=32):
    
    chain = load_data(file_patterns)
    if chain == None:
        return dict()
            
    rdf = ROOT.RDataFrame(chain)

    # Make sure all triggers have a default value
    triggers = [
        # HLT
        "HLT_DoubleMu4_3_LowMass", "HLT_Mu4_L1DoubleMu", "HLT_Mu0_L1DoubleMu",
        "HLT_Mu3_PFJet40", "HLT_Mu8",
        # Run2022
        "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", "L1_DoubleMu4_SQ_OS_dR_Max1p2",
        # Run2023
        "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6",
        "L1_DoubleMu4_SQ_OS_dR_Max1p2",
        # Run2024
        "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
        "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"
    ]
    for trigger in triggers:
        rdf = rdf.DefaultValueFor(trigger, False)

    # Extract information from other branches
    rdf = rdf.Define("mm_mu1_mediumId", "Take(Muon_mediumId,             mm_mu1_index)")
    rdf = rdf.Define("mm_mu1_npixels",  "Take(MuonId_nPixels,            mm_mu1_index)")
    rdf = rdf.Define("mm_mu2_mediumId", "Take(Muon_mediumId,             mm_mu2_index)")
    rdf = rdf.Define("mm_mu2_npixels",  "Take(MuonId_nPixels,            mm_mu2_index)")
        
    # Offline selection
    rdf = rdf.Define("cands", "mm_mu1_index>=0 && mm_mu2_index>=0 && "
                     "mm_mu1_pt > 4 && mm_mu2_pt > 3 && "
                     "(mm_mu1_npixels > 0 && mm_mu2_npixels > 0) &&" # trigger requirement
                     # " mm_kin_lxy > 1 &&"
                     # "mm_mu1_pt > 4 && mm_mu2_pt > 4 &&"
                     # "abs(mm_mu1_eta)<1.4 && abs(mm_mu2_eta)<1.4 &&"
                     # "(abs(mm_mu1_eta)>1.4 || abs(mm_mu2_eta)>1.4) &&"
                     "mm_mu1_mediumId && mm_mu2_mediumId && "
                     "mm_mu1_pdgId * mm_mu2_pdgId == -169 && "
                     "abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01")

    if selection:
        rdf = rdf.Filter(selection)
    rdf = rdf.Filter("Sum(cands)>0")
    
    histos = dict()

    rdf_lxy = rdf.Define("lxy", "mm_kin_lxy[cands]")
    histos["h_probe_lxy"] = rdf_lxy.Histo1D(("h_probe_lxy","", nbins, 0, 16), "lxy")
    histos["h_test_lxy"]  = rdf_lxy.Filter("HLT_DoubleMu4_3_LowMass").Histo1D(("h_test_lxy","", nbins, 0, 16), "lxy")
    rdf_l1 = rdf_lxy.Filter("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2")
    histos["h_test_l1_lxy"]  = rdf_l1.Histo1D(("h_test_l1_lxy","", nbins, 0, 16), "lxy")
    histos["h_test_hlt_lxy"] = rdf_l1.Filter("HLT_DoubleMu4_3_LowMass").Histo1D(("h_test_hlt_lxy","", nbins, 0, 16), "lxy")
    ROOT.RDF.RunGraphs(histos.values())
    return histos

        
def make_trigger_efficiency_plots(name, h_trig, h_all):
    h_eff = h_trig.Clone(name)
    h_eff.Divide(h_trig.GetPtr(), h_all.GetPtr(), 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.SetMaximum(1.0)
    h_eff.SetMarkerStyle(20)
    h_eff.SetMarkerColor(ROOT.kBlue)
    h_eff.SetMarkerSize(0.5)
    h_eff.Draw("e0")
    return h_eff

def plot_results(histos, prefix):
    histos["h_probe_lxy"].GetXaxis().SetTitle("L_{xy}, cm")
    histos["h_probe_lxy"].SetMinimum(0)
    histos["h_probe_lxy"].SetMarkerStyle(20)
    histos["h_probe_lxy"].SetMarkerSize(0.5)
    histos["h_probe_lxy"].Draw("e0")
    print_canvas(f"{prefix}_probe_lxy", output_path)
    
    h_eff_lxy = make_trigger_efficiency_plots("_eff_lxy", histos["h_test_lxy"], histos["h_probe_lxy"])
    h_eff_lxy.GetXaxis().SetTitle("L_{xy}, cm")
    h_eff_lxy.Draw("e0")
    print_canvas(f"{prefix}_eff_lxy", output_path)
    
    h_eff_l1_lxy = make_trigger_efficiency_plots("h_eff_l1_lxy", histos["h_test_l1_lxy"], histos["h_probe_lxy"])
    h_eff_l1_lxy.GetXaxis().SetTitle("L_{xy}, cm")
    h_eff_l1_lxy.Draw("e0")
    print_canvas(f"{prefix}_eff_l1_lxy", output_path)

    histos["h_eff_hlt_lxy"] = make_trigger_efficiency_plots("h_eff_hlt_lxy", histos["h_test_hlt_lxy"], histos["h_test_l1_lxy"])
    histos["h_eff_hlt_lxy"].GetXaxis().SetTitle("L_{xy}, cm")
    histos["h_eff_hlt_lxy"].Draw("e0")
    print_canvas(f"{prefix}_eff_hlt_vs_l1_lxy", output_path)
    

####################################################################

patterns_K0sToMuMu = [
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*.root",
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v1+MINIAODSIM/*.root",
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4+MINIAODSIM/*.root",
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4+MINIAODSIM/*.root",
    # path + "/K0sToMuMu_Fil-K0s_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v1+MINIAODSIM/*.root",
]

patterns_HLTPhysics = [
    path2 + "HLTPhysics+Run2022C-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2022D-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2022D-PromptReco-v2+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2022E-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2022F-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2022G-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2023C-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2023C-PromptReco-v2+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2023C-PromptReco-v3+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2023C-PromptReco-v4+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2023D-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2023D-PromptReco-v2+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024C-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024D-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024E-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024E-PromptReco-v2+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024F-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024G-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024H-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024I-PromptReco-v1+MINIAOD/*.root",
    path2 + "HLTPhysics+Run2024I-PromptReco-v2+MINIAOD/*.root",
]

ROOT.gROOT.SetBatch(True)
# ROOT.gStyle.SetOptStat(0)
tdrstyle.setTDRStyle()
c1 = ROOT.TCanvas("c1","", 500, 500)

mc_histos = process_sample(patterns_K0sToMuMu)

if mc_histos:
    plot_results(mc_histos, "mc")

mc_histos2 = process_sample(patterns_K0sToMuMu, nbins=16)
if mc_histos2:
    plot_results(mc_histos2, "mc2")

data_histos = process_sample(patterns_HLTPhysics, selection="L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2", nbins=16)
if data_histos:
    plot_results(data_histos, "data")

if mc_histos2 and data_histos:
    mc_histos2["h_eff_hlt_lxy"].SetMarkerColor(ROOT.kRed)
    mc_histos2["h_eff_hlt_lxy"].Draw("e0")
    data_histos["h_eff_hlt_lxy"].Draw("e0 same")
    
    print_canvas("data_n_mc_eff_hlt_vs_l1_lxy", output_path)

    h_ratio = mc_histos2["h_eff_hlt_lxy"].Clone("h_eff_ratio_hlt_lxy")
    h_ratio.Divide(data_histos["h_eff_hlt_lxy"], mc_histos2["h_eff_hlt_lxy"])
    h_ratio.SetMinimum(0)
    h_ratio.SetMaximum(1.3)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerColor(ROOT.kBlue)
    h_ratio.SetMarkerSize(0.5)
    h_ratio.GetXaxis().SetTitle("L_{xy}, cm")
    fit_result = h_ratio.Fit("pol3", "S")
    corr_function = h_ratio.GetFunction("pol3")
    h_ratio.Draw("e0")
    print_canvas("h_eff_ratio_hlt_lxy", output_path)

    # Compute average efficiency correction
    numerator = 0.0
    denominator = 0.0
    h = mc_histos2["h_probe_lxy"]
    n_bins = h.GetNbinsX()

    # Accumulate weighted powers of x (for error propagation)
    max_order = 3
    x_moments = np.zeros(max_order + 1)

    for bin in range(1, n_bins + 1):
        x = h.GetBinCenter(bin)
        w = h.GetBinContent(bin)
        f = corr_function.Eval(x)

        # Weighted sum for the average
        numerator += f * w
        denominator += w

        # Store x^i * weight
        for i in range(max_order + 1):
            x_moments[i] += w * x**i

    if denominator > 0:
        avg_corr = numerator / denominator
        norm_moments = x_moments / denominator
    else:
        avg_corr = 0.0
        norm_moments = np.zeros_like(x_moments)

    # Error propagation using the covariance matrix
    cov = fit_result.GetCovarianceMatrix()
    variance = 0.0
    for i in range(max_order + 1):
        for j in range(max_order + 1):
            variance += norm_moments[i] * cov[i][j] * norm_moments[j]

    avg_corr_err = np.sqrt(variance)

    print(f"Average efficiency correction: {avg_corr:.4f} ± {avg_corr_err:.4f}")
    
    # for bin in range(1, mc_histos2["h_probe_lxy"].GetNbinsX() + 1):
    #     x = h.GetBinCenter(bin)
    #     weight = h.GetBinContent(bin)
    #     corr = corr_function.Eval(x)
    #     numerator += corr * weight
    #     denominator += weight

    # if denominator > 0:
    #     avg_corr = numerator / denominator
    # else:
    #     avg_corr = 0.0
        
    # print(f"Average efficiency correction: {avg_corr:.4f}")
