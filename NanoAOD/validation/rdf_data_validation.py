import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob
from collections import defaultdict
import pprint
from math import sqrt
import numpy as np
import json
import pickle

ROOT.ROOT.EnableImplicitMT()

create_used_lumi_json = False
restrict_to_processed_lumis = True
process_data = False
recompute_results = False
aggregate_eras = True

# if eras list is not empty, only listed eras will be used to make plots
# eras = ["Run2023C", "Run2023D", "Run2024C", "Run2024D", "Run2025G"]
# eras = ["Run2024C", "Run2025G"]
pds = ["ParkingDoubleMuonLowMass"]
# pds = []
eras = []
eras.extend([f"Run2023{chr(c)}" for c in range(ord("C"), ord("D") + 1)])
eras.extend([f"Run2024{chr(c)}" for c in range(ord("C"), ord("I") + 1)])
eras.extend([f"Run2025{chr(c)}" for c in range(ord("C"), ord("G") + 1)])

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/535/"
# output_path = "/eos/home-d/dmytro/www/plots/2025/run3_data_validation"
output_path = "/eos/home-d/dmytro/www/plots/tmp/2025/run3_data_validation"

# histos_file = "rdf_data_validation_new.root"
# histos_file = "rdf_data_validation.root"
histos_file = "rdf_data_validation.root"
file_fit_results = f"{output_path}/results.json"
json_folder = "json/"


min_jpsi_mass = 2.8
max_jpsi_mass = 3.3
nbins_jpsi = 50
min_jpsik_mass = 5.15
max_jpsik_mass = 5.55
nbins_jpsik = 80
min_jpsiphi_mass = 5.2
max_jpsiphi_mass = 5.5
nbins_jpsiphi = 60
rrv_jpsi_mass  = ROOT.RooRealVar("rrv_jpsi_mass", "", (max_jpsi_mass + min_jpsi_mass) / 2, min_jpsi_mass, max_jpsi_mass)
rrv_jpsik_mass = ROOT.RooRealVar("rrv_jpsik_mass", "", (max_jpsik_mass + min_jpsik_mass) / 2, min_jpsik_mass, max_jpsik_mass)
rrv_jpsiphi_mass = ROOT.RooRealVar("rrv_jpsiphi_mass", "", (max_jpsiphi_mass + min_jpsiphi_mass) / 2, min_jpsiphi_mass, max_jpsiphi_mass)
rrv_jpsi_mass_binning  = ROOT.RooFit.Binning(nbins_jpsi, min_jpsi_mass, max_jpsi_mass)
rrv_jpsik_mass_binning = ROOT.RooFit.Binning(nbins_jpsik, min_jpsik_mass, max_jpsik_mass)
rrv_jpsiphi_mass_binning = ROOT.RooFit.Binning(nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass)
print_level = 0
# pattern = "11*.root"
# pattern = "^11.*.root"
pattern = "^[^\/]+.root$"

histos = defaultdict(lambda: defaultdict(dict))
histograms = defaultdict(dict)

results = defaultdict(lambda: defaultdict(dict))
certified_run_lumi_set = dict()

# Integrated lumi for data with muon certification
# lumi : ( integrated lumi in fb^-1, relative uncertainty in %)
# https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3
# Muon
lumi = {
    "Run2022C": ( 6.25, 1.3),
    "Run2022D": ( 3.34, 1.3),
    "Run2022E": ( 6.08, 1.3),
    "Run2022F": (18.30, 1.3),
    "Run2022G": ( 3.11, 1.3),
    "Run2023C": (18.39, 1.4),
    "Run2023D": ( 9.82, 1.4),
    "Run2024C": ( 7.70, 4.0),
    "Run2024D": ( 8.37, 4.0),
    "Run2024E": (11.61, 4.0),
    "Run2024F": (28.00, 4.0),
    "Run2024G": (38.47, 4.0),
    "Run2024H": ( 5.48, 4.0),
    "Run2024I": (11.49, 4.0),
    "Run2025C": (20.03, 4.0),
}    
# # Golden
# lumi = {
#     "Run2022C": ( 5.01, 1.3),
#     "Run2022D": ( 2.97, 1.3),
#     "Run2022E": ( 5.81, 1.3),
#     "Run2022F": (17.78, 1.3),
#     "Run2022G": ( 3.08, 1.3),
# }

def unpack_histo_key(key):
    parts = key.split("__")
    if len(parts) != 3:
        return (None,) * 3
    return tuple(parts)

def histo_key(pd, era, selection):
    return f"{pd}__{era}__{selection}";

def add_hist(hist_full_name, hist):
    print(hist_full_name)
    (pd, era, selection_name) = unpack_histo_key(hist_full_name)
    if selection_name == None:
        return
    histos[pd][era][selection_name] = hist
    
def load_histos():
    histos.clear() 
    f = ROOT.TFile(histos_file)
    try:
        keys = f.GetListOfKeys()
        if not keys:
            print(f"Cannot find keys in the {histos_file}")
            return False

        for k in keys:
            hist = k.ReadObj()
            hist_types = ["TH1", "TH2"]
            supported = False
            for hist_type in hist_types:
                if hist.InheritsFrom(hist_type):
                    supported = True
                    break
            if supported:
                f.Remove(hist)
                # if not re.search("(jpsi|jpsiphi)_pix1_mass", k.GetName()):
                #     continue
                # if re.search("Run2022", k.GetName()):
                #     continue
                # if not re.search("jpsi_vloose_trig_HLT_Mu0", k.GetName()):
                #     continue
                # if not re.search("jpsiphi_loose_mass", k.GetName()):
                #     continue
                add_hist(k.GetName(), hist)
    finally:
        f.Close()

    return True

ROOT.gInterpreter.Declare(r"""
#include <ROOT/RVec.hxx>
using ROOT::VecOps::RVec;

int count_bits_7(unsigned int mask) {
    unsigned int m = mask & 0x7F;      // 0b1111111
    return __builtin_popcount(m);      // works with Cling/GCC
}

RVec<int> get_nlayers(const RVec<unsigned int>& algoMask) {
    RVec<int> out(algoMask.size());
    for (size_t i = 0; i < algoMask.size(); ++i) {
        out[i] = count_bits_7(algoMask[i]);
    }
    return out;
}
""")


def book_histos(rdf, pd, era, trigger=None):
    # Add defaults
    rdf = rdf.DefaultValueFor("HLT_DoubleMu4_3_LowMass", False)
    rdf = rdf.DefaultValueFor("HLT_Mu0_L1DoubleMu", False)
    rdf = rdf.DefaultValueFor("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6", False)
    
    rdf = rdf.Define("certified", "passed_lumi_mask(run, luminosityBlock)")
    rdf = rdf.Filter("certified == 1", "passed data certification")
    rdf = rdf.Define("bkkmm_mu1_pixelPattern", "Take(mm_mu1_pixelPattern,    bkkmm_mm_index)")
    rdf = rdf.Define("bkkmm_mu2_pixelPattern", "Take(mm_mu2_pixelPattern,    bkkmm_mm_index)")
    rdf = rdf.Define("mm_mu1_mediumId",        "Take(Muon_mediumId,             mm_mu1_index)")
    rdf = rdf.Define("mm_mu2_mediumId",        "Take(Muon_mediumId,             mm_mu2_index)")
    rdf = rdf.Define("mm_mu1_algoMask",        "Take(MuonId_algoMask,           mm_mu1_index)")
    rdf = rdf.Define("mm_mu2_algoMask",        "Take(MuonId_algoMask,           mm_mu2_index)")
    # rdf = rdf.Define("mm_mu1_nPixelLayers",    "Take(get_nlayers(MuonId_algoMask), mm_mu1_index)")
    # rdf = rdf.Define("mm_mu2_nPixelLayers",    "Take(get_nlayers(MuonId_algoMask), mm_mu2_index)")
    rdf = rdf.Define("mm_mu1_nPixelLayers",    "get_nlayers(mm_mu1_pixelPattern)")
    rdf = rdf.Define("mm_mu2_nPixelLayers",    "get_nlayers(mm_mu2_pixelPattern)")
        
    if trigger != None:
        rdf = rdf.Filter(trigger)

    # if BPHNano:
    #     rdf = rdf.Define("jpsi_cands", f"MuMu_svprob>0.1 && MuMu_fit_mass>{min_jpsi_mass} && MuMu_fit_mass<{max_jpsi_mass}")
    #     rdf = rdf.Define("jpsi_mass","MuMu_fit_mass[jpsi_cands]")

    selections = dict()
    
    #### Jpsi

    ### Very loose - no kinematic fit
    
    selections["jpsi_vvloose"] = f"mm_mass>{min_jpsi_mass} && mm_mass<{max_jpsi_mass}"
    selections["jpsi_vloose"] = selections["jpsi_vvloose"] + \
        " && mm_mu1_mediumId && mm_mu2_mediumId"
    selections["jpsi_vloose_pix1"] = selections["jpsi_vloose"] + \
        "&& (mm_mu1_pixelPattern&1)==1 && (mm_mu2_pixelPattern&1)==1"
    selections["jpsi_vloose_tight_vtx"] = selections["jpsi_vloose"] + \
        "&& mm_kin_vtx_prob>0.1"
    selections["jpsi_vloose_not_tight_vtx"] = selections["jpsi_vloose"] + \
        "&& mm_kin_vtx_prob<0.1"
    selections["jpsi_vloose_trig"] = selections["jpsi_vloose"] + \
        "&& mm_mu1_pt>4 && mm_mu2_pt>3"

    for selection_name, selection in selections.items():
        # define selection
        rdf = rdf.Define(selection_name, selection)
        # define variables
        rdf = rdf.Define(f"{selection_name}_mass", f"mm_mass[{selection_name}]")
        # book histograms
        histos[pd][era][f"{selection_name}_mass"] = \
            rdf.Histo1D((f"{selection_name}_mass",";Mass, GeV", \
                         nbins_jpsi, min_jpsi_mass, max_jpsi_mass),
                        f"{selection_name}_mass")
    ## Trigger studies
    if pd == "EGamma":
        histos[pd][era]["jpsi_vloose_trig_HLT_Mu0_L1DoubleMu_mass"] = \
            rdf.Filter("HLT_Mu0_L1DoubleMu").Histo1D(
                (f"jpsi_vloose_trig_HLT_Mu0_L1DoubleMu_mass",";Mass, GeV", \
                 nbins_jpsi, min_jpsi_mass, max_jpsi_mass),
                f"jpsi_vloose_trig_mass")
        histos[pd][era]["jpsi_vloose_trig_HLT_Mu0_L1DoubleMu_not_HLT_DoubleMu4_3_LowMass_mass"] = \
            rdf.Filter("HLT_Mu0_L1DoubleMu&&!HLT_DoubleMu4_3_LowMass").Histo1D(
                (f"jpsi_vloose_trig_HLT_Mu0_L1DoubleMu_notHLT_DoubleMu4_3_LowMass_mass",
                 ";Mass, GeV", \
                 nbins_jpsi, min_jpsi_mass, max_jpsi_mass),
                f"jpsi_vloose_trig_mass")

    if pd == "ParkingDoubleMuonLowMass":
        histos[pd][era]["jpsi_vvloose_fixed_trig_mass"] = \
            rdf.Filter("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6&&HLT_DoubleMu4_3_LowMass").Histo1D(
                (f"jpsi_vvloose_fixed_trig_mass",";Mass, GeV", \
                 nbins_jpsi, min_jpsi_mass, max_jpsi_mass),
                f"jpsi_vvloose_mass")

    ### Nominal - kinematic fits
    
    selections.clear()

    # 31 bit mask exlcuding
    # * muonSeededStepInOut = 13,
    # * muonSeededStepOutIn = 14,
    # '0b1111111111111111001111111111111'
    trk_algo_mask = 2147459071
    
    selections["jpsi"] = f"mm_kin_vtx_prob>0.1" + \
        f"&& mm_kin_mass>{min_jpsi_mass} && mm_kin_mass<{max_jpsi_mass}" + \
        "&& mm_mu1_mediumId && mm_mu2_mediumId"
    selections["jpsi_lt0p8"]       = selections["jpsi"] + "&& abs(mm_kin_eta)<0.8"
    selections["jpsi_gt0p8_lt1p4"] = selections["jpsi"] + "&& abs(mm_kin_eta)>0.8 && abs(mm_kin_eta)<1.4"
    selections["jpsi_gt1p4"]       = selections["jpsi"] + "&& abs(mm_kin_eta)>1.4"
    selections["jpsi_displaced"]   = selections["jpsi"] + "&& mm_kin_sl3d>5"
    selections["jpsi_pix1"]        = selections["jpsi"] + \
        "&& (mm_mu1_pixelPattern&1)==1 && (mm_mu2_pixelPattern&1)==1"
    selections["jpsi_pix1_lt0p8"]       = selections["jpsi_lt0p8"] + \
        "&& (mm_mu1_pixelPattern&1)==1 && (mm_mu2_pixelPattern&1)==1"
    selections["jpsi_pix1_gt0p8_lt1p4"] = selections["jpsi_gt0p8_lt1p4"] + \
        "&& (mm_mu1_pixelPattern&1)==1 && (mm_mu2_pixelPattern&1)==1"
    selections["jpsi_pix1_gt1p4"]       = selections["jpsi_gt1p4"] + \
        "&& (mm_mu1_pixelPattern&1)==1 && (mm_mu2_pixelPattern&1)==1"
    selections["jpsi_pix1_displaced"]   = selections["jpsi_displaced"] + \
        "&& (mm_mu1_pixelPattern&1)==1 && (mm_mu2_pixelPattern&1)==1"

    for i in range(2):
        selections[f"jpsi_mu{i+1}_pix1_mualgo"] = selections["jpsi"] + \
            f"&& (mm_mu{i+1}_pixelPattern&1)==1 && (mm_mu{i+1}_algoMask&8192)>0"
        selections[f"jpsi_mu{i+1}_pix1_mualgo_trkalgo"] = selections[f"jpsi_mu{i+1}_pix1_mualgo"] + \
            f"&& (mm_mu{i+1}_algoMask&{trk_algo_mask})>0"
        selections[f"jpsi_mu{i+1}_nopix1_mualgo"] = selections["jpsi"] + \
            f"&& (mm_mu{i+1}_pixelPattern&1)==0 && (mm_mu{i+1}_algoMask&8192)>0"
        selections[f"jpsi_mu{i+1}_nopix1_mualgo_trkalgo"] = selections[f"jpsi_mu{i+1}_nopix1_mualgo"] + \
            f"&& (mm_mu{i+1}_algoMask&{trk_algo_mask})>0"

        for j in range(4):
            selections[f"jpsi_mu{i+1}_nopix1_bpix{j}_mualgo"] = selections["jpsi"] + \
                f"&& (mm_mu{i+1}_pixelPattern&1)==0 && (mm_mu{i+1}_algoMask&8192)>0" + \
                f"&& mm_mu{i+1}_nPixelLayers=={j}"
            selections[f"jpsi_mu{i+1}_nopix1_bpix{j}_mualgo_trkalgo"] = \
                selections[f"jpsi_mu{i+1}_nopix1_bpix{j}_mualgo"] + \
                f"&& (mm_mu{i+1}_algoMask&{trk_algo_mask})>0"

    for selection_name, selection in selections.items():
        # define selection
        rdf = rdf.Define(selection_name, selection)
        # define variables
        rdf = rdf.Define(f"{selection_name}_mass", f"mm_kin_mass[{selection_name}]")
        # book histograms
        histos[pd][era][f"{selection_name}_mass"] = rdf.Histo1D((f"{selection_name}_mass",";Mass, GeV", \
                                                       nbins_jpsi, min_jpsi_mass, max_jpsi_mass),
                                                      f"{selection_name}_mass")
        match = re.search("^jpsi_mu(\d)_nopix1_mualgo$", selection_name)
        if match:
            rdf = rdf.Define(f"{selection_name}_nPixelLayers",
                             f"mm_mu{match.group(1)}_nPixelLayers[{selection_name}]")
            histos[pd][era][f"{selection_name}_nPixelLayers"] = \
                rdf.Histo1D((f"{selection_name}_nPixelLayers",";N", 10,0,10),
                            f"{selection_name}_nPixelLayers")
            

    if pd == "ParkingDoubleMuonLowMass":
        histos[pd][era]["jpsi_fixed_trig_mass"] = \
            rdf.Filter("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6&&HLT_DoubleMu4_3_LowMass").Histo1D(
                (f"jpsi_fixed_trig_mass",";Mass, GeV", \
                 nbins_jpsi, min_jpsi_mass, max_jpsi_mass),
                f"jpsi_mass")
        
        # rdf_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6 = rdf.Filter("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6")
        # histos["jpsi_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6"][era]  = \
        #     rdf_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6.Histo1D(("h_jpsi_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
        #                                                    "Vertex constrained dimuon mass;Mass, GeV",
        #                                                    nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "jpsi_mass")
        # rdf_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = rdf.Filter("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4")
        # histos["jpsi_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"][era]  = \
        #     rdf_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4.Histo1D(("h_jpsi_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        #                                                    "Vertex constrained dimuon mass;Mass, GeV",
        #                                                    nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "jpsi_mass")
    
        
    ### BsToJpsiPhi
    
    # rdf = rdf.Define("jpsiphi_cands", f"bkkmm_jpsikk_vtx_prob>0.1 && bkkmm_jpsikk_sl3d>5 && "
    #                  f"abs(bkkmm_jpsikk_alpha)<0.01 && bkkmm_jpsikk_mass>{min_jpsiphi_mass} && "
    #                  f"bkkmm_jpsikk_mass<{max_jpsiphi_mass} && abs(bkkmm_kk_mass-1.02)<0.01")
    
    # rdf = rdf.Define("jpsiphi_mass","bkkmm_jpsikk_mass[jpsiphi_cands]")
    selections.clear()
    selections["jpsiphi_loose"] = f"bkkmm_jpsikk_mass>{min_jpsiphi_mass} && " + \
        f"bkkmm_jpsikk_mass<{max_jpsiphi_mass} && abs(bkkmm_kk_mass-1.02)<0.01"
    selections["jpsiphi_loose_vtx"] = selections["jpsiphi_loose"] + \
        "&& bkkmm_jpsikk_vtx_prob>0.1"
    selections["jpsiphi_loose_vtx_displaced"] = selections["jpsiphi_loose_vtx"] + \
        "&& bkkmm_jpsikk_sl3d>5"
    selections["jpsiphi"] = selections["jpsiphi_loose_vtx_displaced"] + \
        "&& abs(bkkmm_jpsikk_alpha)<0.01"
    selections["jpsiphi_lt0p8"]       = selections["jpsiphi"] + "&& abs(bkkmm_jpsikk_eta)<0.8"
    selections["jpsiphi_gt0p8_lt1p4"] = selections["jpsiphi"] + \
        "&& abs(bkkmm_jpsikk_eta)>0.8 && abs(bkkmm_jpsikk_eta)<1.4"
    selections["jpsiphi_gt1p4"]       = selections["jpsiphi"] + "&& abs(bkkmm_jpsikk_eta)>1.4"
    selections["jpsiphi_pix1"] = selections["jpsiphi"] + \
        "&& (bkkmm_mu1_pixelPattern&1)==1 && (bkkmm_mu2_pixelPattern&1)==1"
    selections["jpsiphi_kaon2m"] = selections["jpsiphi"] + \
        "&& bkkmm_kaon1_pt<2 && bkkmm_kaon1_pt<2"
    selections["jpsiphi_kaon2"] = selections["jpsiphi"] + \
        "&& bkkmm_kaon1_pt>2 && bkkmm_kaon1_pt>2"
    selections["jpsiphi_kaon3"] = selections["jpsiphi"] + \
        "&& bkkmm_kaon1_pt>3 && bkkmm_kaon1_pt>3"
    selections["jpsiphi_kaon4"] = selections["jpsiphi"] + \
        "&& bkkmm_kaon1_pt>4 && bkkmm_kaon1_pt>4"
    selections["jpsiphi_kaon5"] = selections["jpsiphi"] + \
        "&& bkkmm_kaon1_pt>5 && bkkmm_kaon1_pt>5"

    for selection_name, selection in selections.items():
        # define selection
        rdf = rdf.Define(selection_name, selection)
        # define variables
        rdf = rdf.Define(f"{selection_name}_mass", f"bkkmm_jpsikk_mass[{selection_name}]")
        # book histograms
        histos[pd][era][f"{selection_name}_mass"] = \
            rdf.Histo1D((f"{selection_name}_mass",";Mass, GeV", \
                         nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass),
                        f"{selection_name}_mass")

    if pd == "ParkingDoubleMuonLowMass":
        histos[pd][era]["jpsiphi_loose_fixed_trig_mass"] = \
            rdf.Filter("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6&&HLT_DoubleMu4_3_LowMass").Histo1D(
                (f"jpsiphi_loose_fixed_trig_mass",";Mass, GeV", \
                 nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass),
                f"jpsiphi_loose_mass")
        histos[pd][era]["jpsiphi_fixed_trig_mass"] = \
            rdf.Filter("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6&&HLT_DoubleMu4_3_LowMass").Histo1D(
                (f"jpsiphi_fixed_trig_mass",";Mass, GeV", \
                 nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass),
                f"jpsiphi_mass")

    # histos["jpsiphi"][era] = rdf.Histo1D(("h_jpsiphi","BtoJpsiPhi;Mass, GeV", \
    #                                     nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass), "jpsiphi_mass")
    # histos["jpsiphi_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6"][era]  = \
    #     rdf_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6.Histo1D(("h_jpsiphi_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
    #                                                    "BtoJpsiPhi;Mass, GeV",
    #                                                    nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass), "jpsiphi_mass")


        # rdf = rdf.Define("jpsik_cands", f"bkmm_jpsimc_vtx_prob>0.1 && bkmm_jpsimc_sl3d>5 && "
    #                  f"abs(bkmm_jpsimc_alpha)<0.01 && bkmm_jpsimc_mass>{min_jpsik_mass} && bkmm_jpsimc_mass<{max_jpsik_mass} && bkmm_kaon_pt>3")
    # rdf = rdf.Define("jpsik_mass","bkmm_jpsimc_mass[jpsik_cands]")

    # histos["jpsik"][era] = rdf.Histo1D(("h_jpsik","BtoJpsiK;Mass, GeV", \
    #                                     nbins_jpsik, min_jpsik_mass, max_jpsik_mass), "jpsik_mass")
    # histos["jpsik_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6"][era]  = \
    #     rdf_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6.Histo1D(("h_jpsik_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
    #                                                    "BtoJpsiK;Mass, GeV",
    #                                                    nbins_jpsik, min_jpsik_mass, max_jpsik_mass), "jpsik_mass")
    # histos["jpsi_lt0p8"][era]  = rdf.Histo1D(("h_jpsi_lt0p8","Vertex constrained dimuon mass;Mass, GeV", \
    #                                     nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "jpsi_mass_lt0p8")
    # histos["jpsi_gt1p4"][era]  = rdf.Histo1D(("h_jpsi_gt1p4","Vertex constrained dimuon mass;Mass, GeV", \
    #                                     nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "jpsi_mass_gt1p4")
    # histos["jpsi_gt1p4_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6"][era]  = \
    #     rdf_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6.Histo1D(("h_jpsi_gt1p4_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
    #                                                    "Vertex constrained dimuon mass;Mass, GeV",
    #                                                    nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "jpsi_mass_gt1p4")

    h_list = []
    for pd in histos:
        for era in histos:
            for selection, h in histos[pd][era].items():
                h.SetLineColor(ROOT.kBlack)
                h.SetLineWidth(2)
                h_list.append(h)
    ROOT.RDF.RunGraphs(h_list)
    n_events = rdf.Count().GetValue()
    print(n_events)
    return n_events

    
def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print(f"{path}/{output_name_without_extention}.png")
    canvas.Print(f"{path}/{output_name_without_extention}.pdf")
    canvas.Print(f"{path}/{output_name_without_extention}.root")


def build_model(workspace_name, mass_var, peak=3.09, search_width=0.04, max_exp_c=0, min_sigma=0.04):
    """Build fit model and save it in a workspace"""

    # JpsiPhi signal pdf
    # bias     = ROOT.RooRealVar("bias", "bias", 0, -0.1, 0.1)
    # sigma    = ROOT.RooRealVar("sigma", "sigma", 0.0001, 0., 0.01)
    # gaussM   = ROOT.RooGaussModel("gaussM", "signal pdf", mass_var, bias, sigma)
    # ref_data = ROOT.RooDataHist("ref_data", "", ROOT.RooArgList(mass_var), ref_jpsik_hist)
    # ref_pdf  = ROOT.RooHistPdf("ref_pdf", "theoretical lineshape", ROOT.RooArgSet(mass_var), ref_data, 2)
    # mass_var.setBins(10000, "fft");
    # sig      = ROOT.RooFFTConvPdf("sig", "smeared distribution", mass_var, ref_pdf, gaussM)

    # # JpsiPi background pdf
    # jpsipi_data = ROOT.RooDataHist("jpsipi_data", "", ROOT.RooArgList(mass_var), ref_jpsipi_hist)
    # jpsipi_pdf  = ROOT.RooHistPdf("jpsipi_pdf", "theoretical lineshape", ROOT.RooArgSet(mass_var), jpsipi_data, 2)
    ## BsToJpsiPhi signal
    
    #sigmaG_John = ROOT.RooRealVar(n_+"sigmaG_John"," sigma ",0.01, 0.001, 1.)
    #gaus_John = ROOT.RooGaussian(n_+"gaus_John","", m, sig_mu, sigmaG_John)
    #JohnG_frac = ROOT.RooRealVar(n_+"JohnG_frac","",0.3,0.,1.0)

    # if workspace_name in ["jpsi", "jpsik"]:
    sig_mu     = ROOT.RooRealVar("sig_mu", "mu", peak, peak - search_width, peak + search_width)
    sig_sigma  = ROOT.RooRealVar("sig_sigma", "sigma", 0.05, min_sigma, 1.0)
    sig_gamma  = ROOT.RooRealVar("sig_gamma", "gamma", 0.2, 0, 1.0)
    sig_delta  = ROOT.RooRealVar("sig_delta", "delta", 1, 0.1, 10)
    sig = ROOT.RooJohnson("sig", "signal", mass_var, sig_mu, sig_sigma, sig_gamma, sig_delta)

    # if workspace_name == "jpsiphi":
    #     # multi-gaussian
    #     G1_mean   = ROOT.RooRealVar("sig_G1_mean",   "", peak, peak - search_width, peak + search_width)
    #     G1_sigma  = ROOT.RooRealVar("sig_G1_sigma",  "", 0.03, 0.005, 0.10)
    #     G2_sigmaL = ROOT.RooRealVar("sig_G2_sigmaL", "", 0.03, 0.01, 0.10)
    #     G2_sigmaR = ROOT.RooRealVar("sig_G2_sigmaR", "", 0.03, 0.01, 0.10)
    #     # sig = ROOT.RooGaussian(  "sig", "", mass_var, G1_mean, G1_sigma)
    #     G1 = ROOT.RooGaussian(  "sig_G1", "", mass_var, G1_mean, G1_sigma)
    #     G2 = ROOT.RooBifurGauss("sig_G2", "", mass_var, G1_mean, G2_sigmaL, G2_sigmaR)
    #     G2_fract = ROOT.RooRealVar("sig_G2_fract","",0.2,0.0,1.0)
    #     sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_fract))
        
    # # multi-gaussian
    # G1_mean  = ROOT.RooRealVar("sig_G1_mean",  "", peak, peak - search_width, peak + search_width)
    # G1_sigma = ROOT.RooRealVar("sig_G1_sigma", "", 0.03, 0.001, 0.10)
    # G2_scale = ROOT.RooRealVar("sig_G2_scale", "", 2.5, 0.2, 7.5)
    # G3_scale = ROOT.RooRealVar("sig_G3_scale", "", 3.0, 0.5, 6.7)
    # G2_sigma = ROOT.RooProduct("sig_G2_sigma", "", ROOT.RooArgList(G1_sigma,G2_scale))
    # G3_sigma = ROOT.RooProduct("sig_G3_sigma", "", ROOT.RooArgList(G1_sigma,G3_scale))
    # G1 = ROOT.RooGaussian("sig_G1", "", mass_var, G1_mean, G1_sigma)
    # G2 = ROOT.RooGaussian("sig_G2", "", mass_var, G1_mean, G2_sigma)
    # G3 = ROOT.RooGaussian("sig_G3", "", mass_var, G1_mean, G3_sigma)
    
    # G2_fract = ROOT.RooRealVar("sig_G2_fract","",0.3,0.0,1.0)
    # G3_fract = ROOT.RooRealVar("sig_G3_fract","",0.2,0.0,1.0)
    # # sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_fract))
    # sig  = ROOT.RooAddPdf("sig"," ",ROOT.RooArgList(G3,G2,G1),ROOT.RooArgList(G2_fract,G3_fract))

    # # CB
    # sig_mean  = ROOT.RooRealVar("sig_mean",  "", peak, peak - search_width, peak + search_width)
    # sig_sigma = ROOT.RooRealVar("sig_sigma", "sigma", 0.03, 0.001, 0.10)
    # sig_tail  = ROOT.RooRealVar("sig_tail",  "tail", 2.8, 0.1, 10.0)
    # sig_pow   = ROOT.RooRealVar("sig_pow0",  "pow", 3, 0, 50)
    # # sig = ROOT.RooCBShape("sig", "signal", mass_var, sig_mean, sig_sigma, sig_tail, sig_pow)
    # sig_cb = ROOT.RooCBShape("sig_cb", "signal", mass_var, sig_mean, sig_sigma, sig_tail, sig_pow)
    
    # G2_scale = ROOT.RooRealVar("sig_G2_scale", "", 2.5, 0.2, 7.5)
    # G3_scale = ROOT.RooRealVar("sig_G3_scale", "", 3.0, 0.5, 6.7)
    # G2_sigma = ROOT.RooProduct("sig_G2_sigma", "", ROOT.RooArgList(sig_sigma,G2_scale))
    # G3_sigma = ROOT.RooProduct("sig_G3_sigma", "", ROOT.RooArgList(sig_sigma,G3_scale))
    # G2 = ROOT.RooGaussian("sig_G2", "", mass_var, sig_mean, G2_sigma)
    # G3 = ROOT.RooGaussian("sig_G3", "", mass_var, sig_mean, G3_sigma)
    
    # G2_fract = ROOT.RooRealVar("sig_G2_fract","",0.3,0.0,1.0)
    # G3_fract = ROOT.RooRealVar("sig_G3_fract","",0.2,0.0,1.0)
    # # sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_fract))
    # # sig  = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G3,G2,sig_cb), ROOT.RooArgList(G2_fract,G3_fract))
    # sig  = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,sig_cb), ROOT.RooArgList(G2_fract))

    
    ## Combinatorial background
    
    # a0    = ROOT.RooRealVar("a0", "a0", -0.8, -1,  1.0)
    # a1    = ROOT.RooRealVar("a1", "a1", 0.0, -0.3, 0.3)
    # bkg   = ROOT.RooChebychev("bkg", "Background", mass_var, ROOT.RooArgList(a0, a1))

    # b0 = ROOT.RooRealVar("b0","b0", 0.4, 1e-5, 1.)
    # # b1 = ROOT.RooRealVar("b1","b1", 0.5, 0, 1.)
    # bkg = ROOT.RooBernstein("bkg","Background", mass_var, ROOT.RooArgList(b0))
    
    exp_c = ROOT.RooRealVar("exp_c","exp_c", max_exp_c, -5., max_exp_c)
    bkg   = ROOT.RooExponential("bkg", "Background", mass_var, exp_c)
    
    Nsig  = ROOT.RooRealVar("Nsig", "Nsig", 1000, 0, 1e9)
    Nbkg  = ROOT.RooRealVar("Nbkg", "Nbkg", 0, 0, 1e9)
    # Njpsipi = ROOT.RooFormulaVar("Njpsipi", "Njpsipi", "@0*%s" % jpsipi_fraction, ROOT.RooArgList(Nsig))
    
    # model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg,jpsipi_pdf), ROOT.RooArgList(Nsig,Nbkg,Njpsipi))
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg), ROOT.RooArgList(Nsig,Nbkg))

    ws = ROOT.RooWorkspace(workspace_name, "")
    getattr(ws,'import')(model)
    return ws

def make_plots():
    for pd in histos:
        if len(pds) > 0 and pd not in pds:
            continue
        for era in histos[pd]:
            if len(eras) > 0 and era not in eras:
                continue
            for name, h in histos[pd][era].items():
                # print(name)
                final_state = "jpsi"
                if re.search('jpsik', name):
                    final_state = "jpsik"
                elif re.search('jpsiphi', name):
                    final_state = "jpsiphi"

                if pd == "EGamma" and final_state != "jpsi":
                    continue
                
                h.Draw("hist")
                print_canvas(f"{pd}_{era}_{name}", f"{output_path}/histograms/{pd}")

                # fit only mass distributions
                if not re.search("_mass$", name):
                    continue
                
                if final_state == "jpsik":
                    ws = ws_jpsik
                    mass = rrv_jpsik_mass
                    mass_binning = rrv_jpsik_mass_binning
                elif final_state == "jpsiphi":
                    ws = ws_jpsiphi
                    mass = rrv_jpsiphi_mass
                    mass_binning = rrv_jpsiphi_mass_binning
                else:
                    ws = ws_jpsi
                    mass = rrv_jpsi_mass
                    mass_binning = rrv_jpsi_mass_binning

                model = ws.pdf("model")

                data = ROOT.RooDataHist("data", "", ROOT.RooArgList(mass), h)
            
                # prefit
                print("======================= prefit =============================")
                ws.var("Nsig").setVal(h.GetEntries() * 0.9)
                ws.var("Nbkg").setVal(h.GetEntries() * 0.1)
                if final_state in ['jpsi', 'jpsik', 'jpsiphi']:
                    sig_sigma_min, sig_sigma_max = \
                        ws.var("sig_sigma").getMin(), ws.var("sig_sigma").getMax()
                    ws.var("sig_sigma").setMin(0.01)
                    ws.var("sig_sigma").setMax(0.1)
                    ws.var("sig_gamma").setVal(0.0)
                    ws.var("sig_gamma").setConstant(True)
                    ws.var("sig_delta").setVal(1.0)
                    ws.var("sig_delta").setConstant(True)
                    
                else:
                    ws.var("sig_G2_fract").setVal(0.0)
                    ws.var("sig_G2_fract").setConstant(True)
                    ws.var("sig_G2_sigmaL").setConstant(True)
                    ws.var("sig_G2_sigmaR").setConstant(True)
                    
                model.fitTo(data, ROOT.RooFit.NumCPU(8),
                            ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
                            ROOT.RooFit.PrintLevel(print_level))

                # final fit
                print("======================= final fit ==============================")
                if final_state in ['jpsi', 'jpsik', 'jpsiphi']:
                    ws.var("sig_sigma").setMin(sig_sigma_min)
                    ws.var("sig_sigma").setMax(sig_sigma_max)
                    ws.var("sig_gamma").setConstant(False)
                    ws.var("sig_delta").setConstant(False)
                else:
                    ws.var("sig_G2_fract").setConstant(False)
                    ws.var("sig_G2_sigmaL").setConstant(False)
                    ws.var("sig_G2_sigmaR").setConstant(False)
                model.fitTo(data,  ROOT.RooFit.NumCPU(8),
                            ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
                            ROOT.RooFit.PrintLevel(print_level))

                ## Plot results

                frame = mass.frame()
                data.plotOn(frame, mass_binning)
                frame.SetMaximum(frame.GetMaximum() * 1.2)
                # model.plotOn(frame, ROOT.RooFit.Components("sig"), ROOT.RooFit.LineColor(ROOT.kRed))
                model.plotOn(frame, ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))
                model.plotOn(frame)
                # print("chiSquare: ", frame.chiSquare(6))
                # print("chiSquare: ", frame.chiSquare("model","data", 6))

                model.paramOn(frame, ROOT.RooFit.Layout(0.7, 0.95, 0.92))
                frame.getAttText().SetTextSize(0.02)
                frame.Draw()
                print_canvas(f"{pd}_{era}_{name}_fit", f"{output_path}/fits/{pd}")

                results[pd][era][name.removesuffix("_mass")] = \
                    (ws.var("Nsig").getVal(), ws.var("Nsig").getError())

def add_files(pd, era, path):
    if len(pds) > 0 and pd not in pds:
        return
    if len(eras) > 0 and era not in eras:
        return
    nfiles = 0
    for root_dir, _, files in os.walk(path):
        for file in files:
            if re.search(f"{pattern}", file):
                full_path = os.path.join(root_dir, file)
                # print(full_path)
                nfiles += 1
                data_files[pd][era].append(full_path)
    print(f"Found {nfiles} files matching /{pattern}/ pattern for {path}")


# Function to convert JSON to custom string format
def json_to_custom_format(filename):
    try:
        with open(filename, 'r') as file:
            json_data = json.load(file)
        
        parts = []
        for run, lumi_ranges in json_data.items():
            lumi_parts = [f"{start}-{end}" for start, end in lumi_ranges]
            parts.append(f"{run}:{','.join(lumi_parts)}")
        return ';'.join(parts)
    except json.JSONDecodeError as e:
        print(f"JSON parsing error: {e}")
        return None

def add_certification_information(filename):
    with open(filename, 'r') as file:
        json_data = json.load(file)

        """
        Convert {"379029": [[1,25],[58,67]], ...}
        into {379029: {1,2,...,25,58,...,67}, ...}.
        """
        for run_str, ranges in json_data.items():
            run = int(run_str)
            certified_run_lumi_set[run] = set()
            for pair in ranges:
                if not (isinstance(pair, (list, tuple)) and len(pair) == 2):
                    raise ValueError(f"Bad range for run {run}: {pair}")
                a, b = int(pair[0]), int(pair[1])
                certified_run_lumi_set[run].update(range(a, b + 1))

    
def load_lumi_masks():
    # Load lumi masks
    if len(certified_run_lumi_set) != 0:
        print("Certification information is loaded already. Skip it.")
        return
    lumi_mask_string = ""

    if restrict_to_processed_lumis:
        # load processing mask
        certification_files = []
        for file in glob.glob(f"{json_folder}/*.json"):
            certification_files.append(file)
            print(file)
            
    else:
        certification_files = [
            # "/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
            "/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Muon.json",
            "/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Muon.json",
            "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Muon.json",
            "/eos/user/c/cmsdqm/www/CAF/certification/Collisions25/Cert_Collisions2025_391658_398595_Muon.json",
        ]
    for json_filename in certification_files:
        if lumi_mask_string != "":
            lumi_mask_string += ';'
        lumi_mask_string += json_to_custom_format(json_filename)
        add_certification_information(json_filename)

    # Check if conversion succeeded
    if lumi_mask_string:
        # Declare the LumiMask class in ROOT using the header file
        with open('LumiMask.h', 'r') as file:
            lumi_mask_code = file.read()
        ROOT.gInterpreter.Declare(lumi_mask_code)

        # Create the LumiMask object using the custom string format
        lumi_mask = ROOT.LumiMask.fromCustomString(lumi_mask_string, 0, 0)
    else:
        raise Exception("Failed to convert JSON to custom string format.")

    ROOT.gInterpreter.Declare(f'''
    std::string lumi_mask_string;

    bool passed_lumi_mask(unsigned int run, unsigned int lumi) {{
        static LumiMask lumi_mask = LumiMask::fromCustomString(lumi_mask_string);
        return lumi_mask.accept(run, lumi);
    }}
    ''')

    ROOT.gInterpreter.ProcessLine(f'lumi_mask_string = "{lumi_mask_string}";')

def save_histograms():
    f = ROOT.TFile(histos_file, "RECREATE")
    for pd in histos:
        for era in histos[pd]:
            for selection, hist in histos[pd][era].items():
                hist_name = histo_key(pd, era, hist.GetName())
                hist.Write(hist_name)
    f.Close()

def make_evolution_plot(data, pd, name, y_axis_title,
                        format = ".1f", max_scale = 1.4, include_error=True):

    labels = sorted(data.keys())
    n_bins = len(labels)

    # Create histogram with labeled bins
    hist = ROOT.TH1F(f"{name}", "", n_bins, 0.5, n_bins + 0.5)
    hist.SetDirectory(0)
    
    for i, era in enumerate(labels):
        value, error = data[era]
        bin_idx = i + 1
        hist.SetBinContent(bin_idx, value)
        hist.SetBinError(bin_idx, error)
        hist.GetXaxis().SetBinLabel(bin_idx, era)

    # Style and draw
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1)
    hist.SetLineWidth(2)
    hist.SetStats(0)
    hist.SetMinimum(0)
    hist.GetYaxis().SetTitle(y_axis_title)

    hist.Draw("E1P")  # Error bars and point markers

    # set limits to make room for text
    hist.SetMaximum(hist.GetMaximum() * max_scale)

    # Add labels with values ± errors
    latex = ROOT.TLatex()
    latex.SetTextAlign(22)
    latex.SetTextSize(0.03)

    for i in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(i)
        y = hist.GetBinContent(i)
        err = hist.GetBinError(i)
        if include_error:
            label = f"{y:{format}} #pm {err:{format}}"
        else:
            label = f"{y:{format}}"
        # Draw text well above the marker and error bar
        y_pos = y + err + 0.05 * hist.GetMaximum()  # adjust offset here
        latex.DrawLatex(x, y_pos, label)


    c.Update()
    print_canvas(f"{pd}_{name}_{labels[0]}-{labels[-1]}", f"{output_path}/evolution")
    return hist

def make_chains(chain_name):
    chains = defaultdict(lambda: defaultdict(lambda: ROOT.TChain(chain_name)))
    for pd in data_files:
        for era, file_patterns in data_files[pd].items():
            for pattern in file_patterns:
                chains[pd][era].Add(pattern)
    return chains

def merge_into_ranges(sorted_vals):
    """Given a sorted list of unique integers, return [[start, end], ...]."""
    if not sorted_vals:
        return []
    ranges = []
    start = prev = sorted_vals[0]
    for v in sorted_vals[1:]:
        if v == prev + 1:
            prev = v
            continue
        ranges.append([start, prev])
        start = prev = v
    ranges.append([start, prev])
    return ranges

def get_run_to_lumi_ranges(chain):
    """Extract set of lumis present for each run."""
    run_to_lumi_set = {}
    nlumis = 0
    for entry in chain: 
        run = int(entry.run)
        lumi = int(entry.luminosityBlock)
        # Skip not certified data - we won't use it
        if run not in certified_run_lumi_set or lumi not in certified_run_lumi_set[run]:
            continue
        run_to_lumi_set.setdefault(run, set()).add(lumi)
        nlumis += 1
    print(f"Total number of lumi sections: {nlumis}")

    # Convert sets to sorted ranges
    run_to_lumi_ranges = {}
    for run, lumis in run_to_lumi_set.items():
        sorted_lumis = sorted(lumis)
        run_to_lumi_ranges[str(run)] = merge_into_ranges(sorted_lumis)
    return run_to_lumi_ranges


    
####################################################################################

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.045, "X")
ROOT.gStyle.SetLabelSize(0.045, "Y")
ROOT.gStyle.SetTitleSize(0.045, "X")
ROOT.gStyle.SetTitleSize(0.045, "Y")
ROOT.gStyle.SetTitleOffset(1.2, "X")
ROOT.gStyle.SetTitleOffset(0.8, "Y")
# ROOT.gStyle.SetPadLeftMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.10)
# ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadBottomMargin(0.10)

c = ROOT.TCanvas("c", "", 800, 800)

ws_jpsi = build_model("jpsi", rrv_jpsi_mass)
ws_jpsik = build_model("jpsik", rrv_jpsik_mass, 5.3, 0.05, -3)
ws_jpsiphi = build_model("jpsiphi", rrv_jpsiphi_mass, 5.36, 0.05, 0, 0.005)

# define data
data_files = defaultdict(lambda: defaultdict(list))

if restrict_to_processed_lumis:
    # load luminosity information
    with open(f"{json_folder}/luminosities.pkl", "rb") as f:
        lumi_map = pickle.load(f)
        # print(lumi_map)

        lumi.clear()
        for file, integrated_lumi in lumi_map.items():
            match = re.search("ParkingDoubleMuonLowMass-(Run\w+)\.json", file)
            if match:
                # WARNING: assume that only 1 out 8 PDs are used
                lumi[match.group(1)] = (integrated_lumi / 8, 4.0)
    print(lumi)


if process_data or create_used_lumi_json:
    
    ## Run2022
    add_files("ParkingDoubleMuonLowMass", "Run2022C", f"{input_path}/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2022D", f"{input_path}/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2022D", f"{input_path}/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v2+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2022E", f"{input_path}/ParkingDoubleMuonLowMass0+Run2022E-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2022F", f"{input_path}/ParkingDoubleMuonLowMass0+Run2022F-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2022G", f"{input_path}/ParkingDoubleMuonLowMass0+Run2022G-PromptReco-v1+MINIAOD/")


    ## Run2023
    for i in range(4):
        add_files("ParkingDoubleMuonLowMass", "Run2023C", f"{input_path}/ParkingDoubleMuonLowMass0+Run2023C-PromptReco-v{i+1}+MINIAOD/")
    for i in range(2):
        add_files("ParkingDoubleMuonLowMass", "Run2023D", f"{input_path}/ParkingDoubleMuonLowMass0+Run2023D-PromptReco-v{i+1}+MINIAOD/")

    ## Run2024
    add_files("ParkingDoubleMuonLowMass", "Run2024C", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024C-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024D", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024D-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024E", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024E-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024E", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024E-PromptReco-v2+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024F", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024F-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024G", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024G-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024H", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024H-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024I", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024I-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2024I", f"{input_path}/ParkingDoubleMuonLowMass0+Run2024I-PromptReco-v2+MINIAOD/")

    for i in range(2):
        add_files("EGamma", "Run2024C", f"{input_path}/EGamma{i}+Run2024C-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024D", f"{input_path}/EGamma{i}+Run2024D-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024E", f"{input_path}/EGamma{i}+Run2024E-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024E", f"{input_path}/EGamma{i}+Run2024E-PromptReco-v2+MINIAOD/")
        add_files("EGamma", "Run2024F", f"{input_path}/EGamma{i}+Run2024F-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024G", f"{input_path}/EGamma{i}+Run2024G-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024H", f"{input_path}/EGamma{i}+Run2024H-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024I", f"{input_path}/EGamma{i}+Run2024I-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2024I", f"{input_path}/EGamma{i}+Run2024I-PromptReco-v2+MINIAOD/")
    
    ## 2025
    add_files("ParkingDoubleMuonLowMass", "Run2025C", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025C-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2025C", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025C-PromptReco-v2+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2025D", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025D-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2025E", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025E-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2025F", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025F-PromptReco-v1+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2025F", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025F-PromptReco-v2+MINIAOD/")
    add_files("ParkingDoubleMuonLowMass", "Run2025G", f"{input_path}/ParkingDoubleMuonLowMass0+Run2025G-PromptReco-v1+MINIAOD/")
    
    for i in range(3):
        add_files("EGamma", "Run2025C", f"{input_path}/EGamma{i}+Run2025C-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2025C", f"{input_path}/EGamma{i}+Run2025C-PromptReco-v2+MINIAOD/")
        add_files("EGamma", "Run2025D", f"{input_path}/EGamma{i}+Run2025D-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2025E", f"{input_path}/EGamma{i}+Run2025E-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2025F", f"{input_path}/EGamma{i}+Run2025F-PromptReco-v1+MINIAOD/")
        add_files("EGamma", "Run2025F", f"{input_path}/EGamma{i}+Run2025F-PromptReco-v2+MINIAOD/")
        add_files("EGamma", "Run2025G", f"{input_path}/EGamma{i}+Run2025G-PromptReco-v1+MINIAOD/")

if process_data:
    print("Processing data")
    load_lumi_masks()

    chains = make_chains("Events")

    for pd in chains:
        if len(pds) > 0 and pd not in pds:
            continue
        for era, chain in chains[pd].items():
            if len(eras) > 0 and era not in eras:
                continue
            print(f"Processing {era}")
            if chain.GetListOfFiles().GetEntries() == 0:
                print("No files selected. Skip it")
                continue
            rdf = ROOT.RDataFrame(chain)
            if pd == "ParkingDoubleMuonLowMass":
                n1 = book_histos(rdf, pd, era, "HLT_DoubleMu4_3_LowMass")
            else:
                n1 = book_histos(rdf, pd, era)
            print(f"n1: {n1}")

    save_histograms()

if create_used_lumi_json:
    print("Extracting processed lumi information")
    load_lumi_masks()
    
    chains = make_chains("LuminosityBlocks")

    for pd in chains:
        if len(pds) > 0 and pd not in pds:
            continue
        for era, chain in chains[pd].items():
            if len(eras) > 0 and era not in eras:
                continue
            print(f"Processing {era}")
            run_to_lumi_ranges = get_run_to_lumi_ranges(chain)
            with open(f"{json_folder}/{pd}-{era}.json", "w") as f:
                json.dump(run_to_lumi_ranges, f, indent=4, sort_keys=True)

load_histos()

if recompute_results or not os.path.exists(file_fit_results):

    make_plots()
    # pprint(results)
    with open(file_fit_results, "w") as f:
        json.dump(results, f, indent=4)
else:
    with open(file_fit_results) as f:
        results = json.load(f)

# prepare aggregated histograms by era
if aggregate_eras:
    for pd in list(results):
        for era in list(results[pd]):
            if len(eras) >0 and era not in eras:
                del results[pd][era]
                continue

            # merge Run2023 results
            if re.search("Run2023", era):
                if "Run2023" not in results[pd]:
                    results[pd]["Run2023"] = dict()
                for selection in results[pd][era]:
                    if selection not in results[pd]["Run2023"]:
                        results[pd]["Run2023"][selection] = [0,0]
                    # print(f"pd: {pd}, era: {era}, selection: {selection}")
                    results[pd]["Run2023"][selection][0] += results[pd][era][selection][0]
                    results[pd]["Run2023"][selection][1] += (results[pd][era][selection][1])**2
                del results[pd][era]
        # postprocess Run2023
        if "Run2023" in results[pd]:
            for selection in results[pd]["Run2023"]:
                results[pd]["Run2023"][selection] = \
                    (results[pd]["Run2023"][selection][0], sqrt(results[pd]["Run2023"][selection][1]))
        
        # sort eras
        results[pd] = {k: results[pd][k] for k in sorted(results[pd])}
    # print(results)

    # lumi
    for era in list(lumi):
        # merge Run2023 results
        if re.search("Run2023", era):
            if "Run2023" not in lumi:
                lumi["Run2023"] = [0,0]
            lumi["Run2023"][0] += lumi[era][0]
            del lumi[era]
    # postprocess Run2023
    if "Run2023" in lumi:
        lumi["Run2023"] = (lumi["Run2023"][0], 4)
        
    # sort eras
    lumi = {k: lumi[k] for k in sorted(lumi)}
    # print(lumi)
    
# Use different proportions for evolution plots

c = ROOT.TCanvas("c1", "", 1800, 600)


### Yield per lumi

reference_era = "Run2023"
relative_evolution_plots = dict()
for selection in ["jpsi", "jpsi_lt0p8", "jpsi_gt0p8_lt1p4", "jpsi_gt1p4",
                  "jpsi_vvloose_fixed_trig", "jpsiphi_loose_fixed_trig",
                  "jpsi_fixed_trig", "jpsiphi_fixed_trig",
                  "jpsiphi", "jpsiphi_lt0p8", "jpsiphi_gt0p8_lt1p4", "jpsiphi_gt1p4"]:
    for pd in results:
        # print(results)
        data = dict()
        if reference_era in results[pd]:
            data_relative = dict()
            reference_values = dict()
            if pd not in relative_evolution_plots:
                relative_evolution_plots[pd] = dict()
            for sel, value in results[pd][reference_era].items():
                reference_values[sel] = value[0]/lumi[reference_era][0]
            
        for era in results[pd]:
            if era not in lumi:
                print(f"WARNING: no lumi info for {era}. Skip it")
                continue
            if selection not in results[pd][era]:
                continue
            result = results[pd][era][selection]
            event_yield = result[0] / lumi[era][0]
            event_yield_err = 0
            if result[0] > 0:
                event_yield_err = sqrt((result[1]/result[0])**2 + (lumi[era][1]/100)**2) * event_yield
            data[era] = (event_yield, event_yield_err)
            if reference_era in results[pd]:
                if era != reference_era:
                    data_relative[era] = (event_yield / reference_values[selection],
                                          event_yield_err / reference_values[selection])
        if len(data) > 0:
            make_evolution_plot(data, pd, selection, "Events*fb", ",.0f", 1.4, False)
            if reference_era in results[pd]:
                relative_evolution_plots[pd][selection] = \
                    make_evolution_plot(data_relative, pd, f"{selection}_relative", "Efficiency", ".4f", 1.3, False)

### Stacked relative evolution plots
colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kOrange+5, ROOT.kGreen+3]
if reference_era in results[pd]:
    selections = {
        "jpsi": {
            "jpsi_lt0p8" : "|#eta| < 0.8",
            "jpsi_gt0p8_lt1p4" : "0.8 < |#eta| < 1.4",
            "jpsi_gt1p4" : "|#eta| > 1.4"
        },
        "jpsiphi": {
            "jpsiphi_lt0p8" : "|#eta| < 0.8",
            "jpsiphi_gt0p8_lt1p4" : "0.8 < |#eta| < 1.4",
            "jpsiphi_gt1p4" : "|#eta| > 1.4"
        },
    }
    for name, selection_list in selections.items():
        for pd, plots in relative_evolution_plots.items():
            labels = None
            legend = ROOT.TLegend(0.70,0.75,0.85,0.87)
            legend.SetShadowColor(ROOT.kWhite)
            legend.SetLineColor(ROOT.kWhite)
            for i, (selection, title) in enumerate(selection_list.items()):
                h = relative_evolution_plots[pd][selection]
                if labels == None:
                    labels = (h.GetXaxis().GetBinLabel(1),
                              h.GetXaxis().GetBinLabel(h.GetNbinsX()))
                h.SetMarkerColor(colors[i])
                h.SetLineColor(colors[i])
                legend.AddEntry(h, title)
                if i == 0:
                    h.Draw("PE")
                else:
                    h.Draw("PE same")
            legend.Draw()
            print_canvas(f"{pd}_{name}_{labels[0]}-{labels[1]}", f"{output_path}/evolution")

    
        
### Efficiencies and fractions

ratio_plots = [
    ("jpsi", "jpsi_vvloose", "jpsi_tight_vtx_efficiency", 1.0),
    ("jpsiphi_loose_vtx", "jpsiphi_loose", "jpsiphi_vtx_efficiency", 1.),
    ("jpsiphi_loose_vtx_displaced", "jpsiphi_loose_vtx", "jpsiphi_displacement_efficiency", 1.),
    ("jpsiphi", "jpsiphi_loose_vtx_displaced", "jpsiphi_alpha_efficiency", 1.),
    # ("jpsi_pix1", "jpsi", "jpsi_pix_layer1_efficiency", 1.0),
    # ("jpsi_pix1", "jpsi", "jpsi_pix_layer1_efficiency", 1.0),
    # ("jpsi_displaced", "jpsi", "jpsi_displacement_efficiency", 1.0),
    # ("jpsi_pix1_displaced", "jpsi_pix1", "jpsi_pix1_displacement_efficiency", 1.0),
    # ("jpsi_vloose_pix1", "jpsi_vloose", "jpsi_vloose_pix_layer1_efficiency"),
    # ("jpsi_vloose_tight_vtx", "jpsi_vloose", "jpsi_vloose_tight_vtx_efficiency"),
    # ("jpsi_vloose_not_tight_vtx", "jpsi_vloose", "jpsi_vloose_tight_vtx_inefficiency"),
    # ("jpsi_mu1_pix1_mualgo_trkalgo", "jpsi_mu1_pix1_mualgo", "jpsi_mu1_pix1_trkmualgo_ratio"),
    # ("jpsi_mu2_pix1_mualgo_trkalgo", "jpsi_mu2_pix1_mualgo", "jpsi_mu2_pix1_trkmualgo_ratio"),
    # ("jpsi_mu1_nopix1_mualgo_trkalgo", "jpsi_mu1_nopix1_mualgo", "jpsi_mu1_nopix1_trkmualgo_ratio"),
    # ("jpsi_mu2_nopix1_mualgo_trkalgo", "jpsi_mu2_nopix1_mualgo", "jpsi_mu2_nopix1_trkmualgo_ratio"),
    # ("jpsi_vloose_trig_HLT_Mu0_L1DoubleMu_notHLT_DoubleMu4_3_LowMass",
    #  "jpsi_vloose_trig_HLT_Mu0_L1DoubleMu",
    #  "jpsi_vloose_trig_dimuon_inefficiency")
    # ("jpsiphi", "jpsi", "jpsiphi_to_jpsi_ratio", 1000.)
    # ("jpsiphi_pix1", "jpsi_pix1", "jpsiphi_pix1_to_jpsi_pix1_ratio", 1000.)
    # ("jpsiphi_loose", "jpsi", "jpsiphi_loose_to_jpsi_ratio", 1000.)
    # ("jpsiphi_kaon2m", "jpsi", "jpsiphi_kaon2m_to_jpsi_ratio", 1000.),
    # ("jpsiphi_kaon2",  "jpsi", "jpsiphi_kaon2_to_jpsi_ratio", 1000.),
    # ("jpsiphi_kaon3",  "jpsi", "jpsiphi_kaon3_to_jpsi_ratio", 1000.),
    # ("jpsiphi_kaon4",  "jpsi", "jpsiphi_kaon4_to_jpsi_ratio", 1000.),
]

for num, den, name, scale in ratio_plots:
    for pd in results:
        # print(results)
        data = dict()
        for era, result in results[pd].items():
            if num not in result:
                print(f"WARNING: {num} results are not available. Skip it")
                print(result)
                continue
            if den not in result:
                print(f"WARNING: {den} results are not available. Skip it")
                print(result)
                continue
            if result[den][0] <= 0:
                continue
            value = result[num][0] / result[den][0]
            error = sqrt((result[num][1] / result[num][0])**2 +
                         (result[den][1] / result[den][0])**2) * value
            data[era] = (value * scale, error * scale)
        if len(data.keys()) > 0:
            make_evolution_plot(data, pd, name, "Efficiency", ".3f")
