import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob
if not hasattr(ROOT.RDataFrame, "DefaultValueFor"):
    sys.exit(f"Error: ROOT.RDataFrame.DefaultValueFor is not available in ROOT {ROOT.gROOT.GetVersion()}. Please use ROOT version 6.34 or later.")

ROOT.ROOT.EnableImplicitMT()

input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/'
# input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/533/'
# input_path = '/data/dmytro/bmm6/NanoAOD/529/'
output_path = "/eos/home-d/dmytro/www/plots/mc_trigger_impact_bjpsik_vs_bsjpsiphi"
# output_path = "/eos/home-d/dmytro/www/plots/tmp/mc_trigger_impact_bjpsik_vs_bsjpsiphi5"
# xrootd_server = "eoscms.cern.ch:1094"
xrootd_server = None
histos = dict()
chains = []
kaon_pt_values = [1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 3.0]
# muon_selection = "mu1_pt > 2 && mu2_pt > 2 && (mu1_pt < 4 || mu2_pt < 3)"
# muon_selection_name = "22_m1pt.lt.4.or.m2pt.lt.3"
muon_selection = "mu1_pt > 4 && mu2_pt > 3"
muon_selection_name = "43"
nbins = 22
min_pt = 8
max_pt = 30
# nfiles = 1
nfiles = 999999
trigger = "HLT_DoubleMu4_3_LowMass"
# trigger = "HLT_DoubleMu4_3_LowMass&&(L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4||L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6||L1_DoubleMu4p5_SQ_OS_dR_Max1p2)"
# trigger = "HLT_Dimuon0_Jpsi"
# trigger = "HLT_DoubleMu2_Jpsi_LowPt"
trigger_name = trigger
# trigger_name = "HLT_DoubleMu4_3_LowMass_CorrectL1Seeds"

def xrd_ls(path):
    if xrootd_server:
        command = f"xrdfs {xrootd_server} ls {path}"
    else:
        command = f"find {path} -type f -name \"*.root\""
    output = subprocess.check_output(command, text=True, shell=True)
    files = []
    for line in output.splitlines():
        if re.search('\.root$', line):
            files.append(line)
    return files[:nfiles]


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))


def get_data(datasets):
    chain = ROOT.TChain("Events")
    chains.append(chain)
    n_files = 0
    for dataset in datasets:
        for file in xrd_ls(f"{input_path}/{dataset}"):
            if xrootd_server:
                n_files += chain.Add(f"root://{xrootd_server}/{file}")
            else:
                n_files += chain.Add(f"{file}")
    print(f"Number of files: {n_files}")
    # print(f"Number of events: {chain.GetEntries()}")
    return ROOT.RDataFrame(chain)


def process_BsToJPsiPhi_data():
    datasets = [
        # 529
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
        "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",

        # # 2024
        # "BsToJPsiPhi-JPsiToMuMu-PhiToKK_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/",
        
        
    ]
    rdf = get_data(datasets)
    rdf = rdf.Define("mu1_pt",       "Take(mm_mu1_pt, bkkmm_mm_index)")
    rdf = rdf.Define("mu2_pt",       "Take(mm_mu2_pt, bkkmm_mm_index)")
    rdf = rdf.Define("mu1_index",    "Take(mm_mu1_index, bkkmm_mm_index)")
    rdf = rdf.Define("mu2_index",    "Take(mm_mu2_index, bkkmm_mm_index)")
    rdf = rdf.Define("mu1_mediumId", "Take(Muon_mediumId, mu1_index)")
    rdf = rdf.Define("mu2_mediumId", "Take(Muon_mediumId, mu2_index)")
    rdf = rdf.Define("vtx_prob",     "Take(mm_kin_vtx_prob, bkkmm_mm_index)")

    rdf = rdf.Define("BsToJPsiPhi", "bkkmm_jpsikk_vtx_prob>0.025 && bkkmm_jpsikk_sl3d>3" \
                     "&& bkkmm_jpsikk_alpha<0.1 && abs(bkkmm_jpsikk_mass-5.4)<0.5" \
                     "&& abs(bkkmm_kk_mass-1.02)<0.03" \
                     "&& bkkmm_gen_pdgId!=0" \
                     "&& mu1_mediumId && mu2_mediumId" \
                     f"&& {muon_selection}" \
                     "&& vtx_prob>0.01")
    h_reco = rdf.Define("b_pt", "bkkmm_jpsikk_pt[BsToJPsiPhi]").\
        Histo1D(("h_reco_BsToJPsiPhi", "", nbins, min_pt, max_pt), "b_pt")
    histos['h_reco_BsToJPsiPhi'] = h_reco
    
    h_trig = rdf.Filter(trigger).Define("b_pt", "bkkmm_jpsikk_pt[BsToJPsiPhi]").\
        Histo1D(("h_trig_BsToJPsiPhi", "", nbins, min_pt, max_pt), "b_pt")
    histos['h_trig_BsToJPsiPhi'] = h_trig
    ROOT.RDF.RunGraphs([h_reco, h_trig])
    n_events = rdf.Count().GetValue()
    return n_events

def process_BuToJpsiK_data():
    datasets = [
        # 529
        "BuToJpsiK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/",
        "BuToJpsiK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",
        "BuToJpsiK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
        "BuToJpsiK_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",
        "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",
        "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
        "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        "ButoJpsiK_Jpsito2Mu_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",

        # # 2024
        # "BuToJpsiK-JpsiToMuMu_Fil-MuPt2_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/",
    ]
    rdf = get_data(datasets)
    rdf = rdf.Define("mu1_pt",          "Take(mm_mu1_pt, bkmm_mm_index)")
    rdf = rdf.Define("mu2_pt",          "Take(mm_mu2_pt, bkmm_mm_index)")
    rdf = rdf.Define("mu1_index",       "Take(mm_mu1_index, bkmm_mm_index)")
    rdf = rdf.Define("mu2_index",       "Take(mm_mu2_index, bkmm_mm_index)")
    rdf = rdf.Define("mu1_mediumId",    "Take(Muon_mediumId, mu1_index)")
    rdf = rdf.Define("mu2_mediumId",    "Take(Muon_mediumId, mu2_index)")
    rdf = rdf.Define("vtx_prob", "Take(mm_kin_vtx_prob, bkmm_mm_index)")

    rdfs = []
    histos_to_process = []
    for i, pt in enumerate(kaon_pt_values):
        rdfs.append(rdf.Define(f"BuToJpsiK_{i}", "bkmm_jpsimc_vtx_prob>0.025 && bkmm_jpsimc_sl3d>3" \
                               "&& bkmm_jpsimc_alpha<0.1 && abs(bkmm_jpsimc_mass-5.4)<0.5" \
                               f"&& bkmm_gen_pdgId!=0 && bkmm_kaon_pt>{pt}" \
                               "&& mu1_mediumId && mu2_mediumId" \
                               f"&& {muon_selection}" \
                               "&& vtx_prob>0.01"))
        h_reco = rdfs[i].Define("b_pt", f"bkmm_jpsimc_pt[BuToJpsiK_{i}]").\
            Histo1D((f"h_reco_BuToJpsiK_{pt}", "", nbins, min_pt, max_pt), "b_pt")
        histos[f"h_reco_BuToJpsiK_{pt}"] = h_reco
        histos_to_process.append(h_reco)
        
        h_trig = rdfs[i].Filter(trigger).Define("b_pt", f"bkmm_jpsimc_pt[BuToJpsiK_{i}]").\
            Histo1D((f"h_trig_BuToJpsiK_{pt}", "", nbins, min_pt, max_pt), "b_pt")
        histos[f"h_trig_BuToJpsiK_{pt}"] = h_trig
        histos_to_process.append(h_trig)
        
    ROOT.RDF.RunGraphs(histos_to_process)
    n_events = rdf.Count().GetValue()
    return n_events


def plot_ratio(name, trig1, ref1, trig2, ref2):
    h1 = trig1.Clone(f"h1_{name}")
    h1.Divide(trig1, ref1, 1, 1, "B")
    h2 = trig2.Clone(f"h2_{name}")
    h2.Divide(trig2, ref2, 1, 1, "B")

    hist = trig1.Clone(name)
    hist.Divide(h1, h2)
    hist.SetMinimum(0.5)
    hist.SetMaximum(1.5)
    hist.SetMarkerColor(ROOT.kBlue)
    hist.SetMarkerStyle(20)
    hist.GetXaxis().SetTitle("p^{B}_{T}, [GeV]")
    # hist.GetXaxis().SetRangeUser(10, 30)
    hist.Draw("e0")
    print_canvas(name, output_path)

n_BsToJPsiPhi = process_BsToJPsiPhi_data()
print(f"Number of BsToJPsiPhi events: {n_BsToJPsiPhi}")
n_BuToJpsiK = process_BuToJpsiK_data()
print(f"Number of BuToJpsiK events: {n_BuToJpsiK}")

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.045, "X")
ROOT.gStyle.SetLabelSize(0.045, "Y")
ROOT.gStyle.SetTitleSize(0.045, "X")
ROOT.gStyle.SetTitleSize(0.045, "Y")
ROOT.gStyle.SetTitleOffset(1.2, "X")
ROOT.gStyle.SetTitleOffset(1.2, "Y")
ROOT.gStyle.SetPadLeftMargin(0.15)
ROOT.gStyle.SetPadBottomMargin(0.15)

c1 = TCanvas("c1", "c1", 600, 600)
ROOT.gPad.SetGrid()

histos[f'h_reco_BsToJPsiPhi'].Draw()
print_canvas(f'h{muon_selection_name}_reco_BsToJPsiPhi', output_path)
histos[f'h_trig_BsToJPsiPhi'].Draw()
print_canvas(f'h{muon_selection_name}_trig_{trigger_name}_BsToJPsiPhi', output_path)

for i, pt in enumerate(kaon_pt_values):
    histos[f'h_reco_BuToJpsiK_{pt}'].Draw()
    print_canvas(f'h{muon_selection_name}_reco_BuToJpsiK_{pt}', output_path)
    histos[f'h_trig_BuToJpsiK_{pt}'].Draw()
    print_canvas(f'h{muon_selection_name}_trig_{trigger_name}_BuToJpsiK_{pt}', output_path)
    
    plot_ratio(f'h{muon_selection_name}_trig_{trigger_name}_ratio_k{pt}',
               histos['h_trig_BsToJPsiPhi'].GetPtr(), histos['h_reco_BsToJPsiPhi'].GetPtr(),
               histos[f'h_trig_BuToJpsiK_{pt}'].GetPtr(), histos[f'h_reco_BuToJpsiK_{pt}'].GetPtr())
