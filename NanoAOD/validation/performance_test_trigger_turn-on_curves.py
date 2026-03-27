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
import subprocess

ROOT.ROOT.EnableImplicitMT()

use_xrootd = True
use_caching = True
xrootd_server = "eoscms.cern.ch:1094"
xcache_server = "localhost:1094"

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)
ROOT.gPad.SetGrid()
# ROOT.gPad.SetLogx()

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/"
# data_path = "/data/dmytro/bmm6/NanoAOD/529/"
data_path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/"
# data_path2 = "/data/dmytro/bmm6/NanoAOD/531/"
output_path = "/eos/home-d/dmytro/www/plots/tmp/Run3/l1_turn-on_curves/";

samples = {
    # 'Bsmm-22EE': {
    #     'datasets': [
    #         data_path + "/BsToMuMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root",
    #     ],
    #     'triggers': [], 
    # },
    # 'ZeroBias-22E': {
    #     'datasets': [
    #         data_path + "/ZeroBias+Run2022E-PromptReco-v1+MINIAOD/*.root",
    #     ],
    #     'triggers': [], 
    # },
    'ParkingDoubleEl': {
        'datasets': [
            data_path2 + '/ParkingDoubleElectronLowMass0+Run2022C-PromptReco-v1+MINIAOD/',
            # data_path2 + '/ParkingDoubleElectronLowMass1+Run2022C-PromptReco-v1+MINIAOD/',
            # data_path2 + '/ParkingDoubleElectronLowMass2+Run2022C-PromptReco-v1+MINIAOD/',
            # data_path2 + '/ParkingDoubleElectronLowMass3+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass4+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass5+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass0+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass1+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass2+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass3+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass4+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass5+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass0+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass1+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass2+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass3+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass4+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass5+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass0+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass1+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass2+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass3+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass4+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleElectronLowMass5+Run2022F-PromptReco-v1+MINIAOD/*.root',
        ],
        'save_triggers': [],
        'all_triggers_are_safe': True,
    },
    'ParkingDoubleMu': {
        'datasets': [
            data_path2 + '/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/',
            # data_path2 + '/ParkingDoubleMuonLowMass1+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass2+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass3+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass4+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass5+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass6+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass7+Run2022C-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass1+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass2+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass3+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass4+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass5+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass6+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass7+Run2022D-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass0+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass1+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass2+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass3+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass4+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass5+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass6+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass7+Run2022E-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass0+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass1+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass2+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass3+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass4+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass5+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass6+Run2022F-PromptReco-v1+MINIAOD/*.root',
            # data_path2 + '/ParkingDoubleMuonLowMass7+Run2022F-PromptReco-v1+MINIAOD/*.root',
        ],
        'save_triggers': ["HLT_Mu4_L1DoubleMu"],
        'all_triggers_are_safe': False,
    },
    # 'Bsmm_test': {
    #     'datasets': [
    #         data_path + "/BsToMuMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/1*.root",
    #     ],
    #     'triggers': [], 
    #     # 'eta_step': 1.6
    #     'eta_step': 9999
    # },
    'B2JpsiK': {
        'datasets': [
            data_path + "/ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        ],
        'save_triggers': [],
        'all_triggers_are_safe': True,
    },
    # 'Ksmm': {
    #     'datasets': [
    #         data_path2 + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15_ext1-v2+MINIAODSIM/*.root",
    #     ],
    #     'save_triggers': [],
    #     'all_triggers_are_safe': True,
    # },
        
}

# Valida sample names
mask = '\w\d\_\.'
for sample_name in samples:
    if re.search(f"[^{mask}]", sample_name):
        raise Exception(f"Illigal symbol used in sample name {sample_name}. Allowed regexp: '{mask}'") 


data = dict()
rdfs = dict()
single_muon_rdfs = dict()
jpsi_l1_rdfs = dict()
jpsi_hlt_rdfs = dict()
histos = dict()
eff_histos = dict()
histos_by_category = dict()

def xrd_ls(path):
    command = ["xrdfs", xrootd_server, "ls", path]
    output = subprocess.check_output(command, text=True)
    files = []
    for line in output.splitlines():
        if re.search('\.root$', line):
            files.append(line)
    return files

# print(xrd_ls("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD"))
# exit()

def split_range(start, end, step):
    """
    Split era range into list of bins
    Example: split_range(-2.4, 2.4, 0.2)
    [(-2.4, -2.2), (-2.2, -2.0), ..., (2.0, 2.2), (2.2, 2.4)]
    """
    result = []
    current = start
    while current < end:
        result.append((round(current, 2), min(round(current + step, 2),end)))
        current += step
    return result


def get_data(name, tree="Events"):
    if name in rdfs:
        return rdfs[name]
    chain = ROOT.TChain(tree)
    n_files = 0 
    for dataset in samples[name]['datasets']:
        if use_xrootd:
            for file in xrd_ls(dataset):
                if use_caching:
                    n_files += chain.Add(f"root://{xcache_server}/{file}")
                else:
                    n_files += chain.Add(f"root://{xrootd_server}/{file}")
        else:
            n_files += chain.Add(f"{dataset}/*.root")
    if n_files == 0:
        print("No files found for " + name)
        return None
    data[name] = chain
    rdfs[name] = ROOT.RDataFrame(chain)
    return rdfs[name]


def get_data_single_muon(name, tree="Events"):
    if name in rdfs:
        return rdfs[name]
    if not samples[name]['all_triggers_are_safe']:
        return None
    rdf = get_data(name)
    single_muon_rdfs[name] = rdf
    return single_muon_rdfs[name]


def get_data_jpsi_l1(name, tree="Events"):
    if name in jpsi_l1_rdfs:
        return jpsi_l1_rdfs[name]
    if not samples[name]['all_triggers_are_safe']:
        return None
    rdf = get_data(name)
    rdf = rdf.Define("jpsis", "abs(mm_kin_mass-3.09)<0.06 && mm_kin_vtx_prob>0.1")
    rdf = rdf.Redefine("Muon_pt", "Concatenate(Take(Muon_pt, mm_mu1_index[jpsis]), Take(Muon_pt, mm_mu2_index[jpsis]))")
    rdf = rdf.Redefine("Muon_eta", "Concatenate(Take(Muon_eta, mm_mu1_index[jpsis]), Take(Muon_eta, mm_mu2_index[jpsis]))")
    rdf = rdf.Redefine("MuonId_l1_quality", "Concatenate(Take(MuonId_l1_quality, mm_mu1_index[jpsis]), Take(MuonId_l1_quality, mm_mu2_index[jpsis]))")
    rdf = rdf.Redefine("Muon_looseId", "Concatenate(Take(Muon_looseId, mm_mu1_index[jpsis]), Take(Muon_looseId, mm_mu2_index[jpsis]))")
    rdf = rdf.Redefine("Muon_mediumId", "Concatenate(Take(Muon_mediumId, mm_mu1_index[jpsis]), Take(Muon_mediumId, mm_mu2_index[jpsis]))")
    rdf = rdf.Redefine("Muon_tightId", "Concatenate(Take(Muon_tightId, mm_mu1_index[jpsis]), Take(Muon_tightId, mm_mu2_index[jpsis]))")
    jpsi_l1_rdfs[name] = rdf
    return jpsi_l1_rdfs[name]


def get_data_jpsi_hlt(name, tree="Events"):

    reference_trigger = "HLT_Mu4_L1DoubleMu"
    
    if name in jpsi_hlt_rdfs:
        return jpsi_hlt_rdfs[name]
    
    if not samples[name]['all_triggers_are_safe'] and \
       reference_trigger not in samples[name]['save_triggers']:
        return None
    
    rdf = get_data(name)
    rdf = rdf.Filter(reference_trigger)
    rdf = rdf.Define("jpsis_hlt", "abs(mm_kin_mass-3.09)<0.06 && mm_kin_vtx_prob>0.1")
    
    # Define jpsi based information
    rdf = rdf.Define("jpsi_mu1_pt",       "Take(Muon_pt,                   mm_mu1_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu1_eta",      "Take(Muon_eta,                  mm_mu1_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu1_looseId",  "Take(Muon_looseId,              mm_mu1_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu1_mediumId", "Take(Muon_mediumId,             mm_mu1_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu2_pt",       "Take(Muon_pt,                   mm_mu2_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu2_eta",      "Take(Muon_eta,                  mm_mu2_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu2_looseId",  "Take(Muon_looseId,              mm_mu2_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_mu2_mediumId", "Take(Muon_mediumId,             mm_mu2_index[jpsis_hlt])")

    # probe is defined as the other muon firing HLT_Mu4_L1DoubleMu
    rdf = rdf.Define("jpsi_probe1",       "Take(MuonId_HLT_Mu4_L1DoubleMu, mm_mu2_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_test1",        "Take(MuonId_HLT_Mu0_L1DoubleMu, mm_mu1_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_probe2",       "Take(MuonId_HLT_Mu4_L1DoubleMu, mm_mu1_index[jpsis_hlt])")
    rdf = rdf.Define("jpsi_test2",        "Take(MuonId_HLT_Mu0_L1DoubleMu, mm_mu2_index[jpsis_hlt])")
    rdf = rdf.Define("L1_jpsi_probe",     "Concatenate(jpsi_probe1, jpsi_probe2)")
    rdf = rdf.Define("HLT_fired",         "Concatenate(jpsi_test1,  jpsi_test2)")
    rdf = rdf.Redefine("Muon_pt",         "Concatenate(jpsi_mu1_pt, jpsi_mu2_pt)")
    rdf = rdf.Redefine("Muon_eta",        "Concatenate(jpsi_mu1_eta, jpsi_mu2_eta)")

    jpsi_hlt_rdfs[name] = rdf
    return jpsi_hlt_rdfs[name]


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))


def book_histograms(rdf, sample_name, suffix, preselection,
                    test_selection, eta_step=9999.):
    eta_bins = split_range(-2.4, 2.4, eta_step)

    count = 0
    for counter, (eta_min, eta_max) in enumerate(eta_bins):
        h_name = f"h_all_eta{eta_min:+.1f}to{eta_max:+.1f}_{suffix}_{sample_name}"
        rdf_name = f"rdf_all_eta{counter}_{suffix}_{sample_name}"
        selection = f"Muon_eta > {eta_min} && Muon_eta < {eta_max}"
        if preselection != "":
            selection += "&&" + preselection
        rdf2 = rdf.Define(rdf_name, selection)
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


def make_overlay_plots(sample_name1, legend_name1, sample_name2, legend_name2):
    for category_name, histograms in histos_by_category.items():
        if sample_name1 not in histograms:
            continue
        else:
            h1_name = histograms[sample_name1]

        if sample_name2 not in histograms:
            continue
        else:
            h2_name = histograms[sample_name2]
        
        legend = ROOT.TLegend(0.70,0.65,0.85,0.77)
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
        print_canvas(f"overlay_{category_name}", output_path)

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
        ratio_plot.GetLowerRefYaxis().SetTitle(f"{legend_name2}/{legend_name1}")

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
        print_canvas(f"overlay_{category_name}_ratio", output_path)
    
        
ROOT.gStyle.SetPaintTextFormat(".3f");

eta_step = 0.4
for sample_name in samples:
    print(f"processing {sample_name}")
    
    single_muon_rdf = get_data_single_muon(sample_name)
    jpsi_l1_rdf = get_data_jpsi_l1(sample_name)
    # jpsi_hlt_rdf = get_data_jpsi_hlt(sample_name)
    
    print(f"RDFs are defined")

    if single_muon_rdf:
        book_histograms(single_muon_rdf, sample_name, "base",
                        "Muon_looseId&&Muon_isTracker&&Muon_isGlobal&&Muon_nStations>=2&&abs(Muon_dxybs)<0.1",
                        "MuonId_l1_quality >= 12", eta_step)
        book_histograms(single_muon_rdf, sample_name, "base_medium",
                        "Muon_mediumId&&Muon_isTracker&&Muon_isGlobal&&Muon_nStations>=2&&abs(Muon_dxybs)<0.1",
                        "MuonId_l1_quality >= 12", eta_step)
        book_histograms(single_muon_rdf, sample_name, "loose",
                        "Muon_looseId",
                        "MuonId_l1_quality >= 12", eta_step)
        book_histograms(single_muon_rdf, sample_name, "medium",
                        "Muon_mediumId",
                        "MuonId_l1_quality >= 12", eta_step)
        book_histograms(single_muon_rdf, sample_name, "tight",
                        "Muon_tightId",
                        "MuonId_l1_quality >= 12", eta_step)

    # if jpsi_l1_rdf:
    #     book_histograms(jpsi_l1_rdf, sample_name, "jpsi_loose",
    #                     "Muon_looseId",
    #                     "MuonId_l1_quality >= 12")
        
    # if jpsi_hlt_rdf:
    #     book_histograms(jpsi_hlt_rdf, sample_name, "jpsi_probe_hlt_eff_looseid",
    #                     "L1_jpsi_probe",
    #                     "HLT_fired")

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
make_overlay_plots("Ksmm", "MC", "ParkingDoubleMu", "Data")
