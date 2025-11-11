import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob
from collections import defaultdict
import pprint
from math import sqrt
import numpy as np
import json
import tdrstyle

# Set the TDR style
tdrstyle.setTDRStyle()

ROOT.ROOT.EnableImplicitMT()

#
# In this study we compare PromptReco and ReReco for Run2025C
# - only common luminosity blocks are considered
# - BuToJpsiK reconstruction with tighter than normal cuts to extrapolate to Bmm selection
#

study = "HLT_DoubleMu4_3_LowMass_Muon"
input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/535/"
output_path = "/eos/home-d/dmytro/www/plots/2025/Run2025C_re-reco_validation"
# output_path = "/eos/home-d/dmytro/www/plots/tmp/2025/Run2025C_re-reco_validation"

histos_file = "rdf_re-reco_validation.root"
create_lumi_mask = False
process_data = False

min_jpsi_mass = 2.8
max_jpsi_mass = 3.3
nbins_jpsi = 50
min_jpsik_mass = 5.17
max_jpsik_mass = 5.55
nbins_jpsik = 76
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

histos = defaultdict(                      # pd
    lambda: defaultdict(                   # version
        lambda: defaultdict(               # process
            lambda: defaultdict(dict)      # selection -> {hist_name: hist}
        )
    )
)
histograms = defaultdict(dict)

def unpack_histo_key(key):
    parts = key.split("__")
    if len(parts) != 5:
        return (None,) * 5
    return tuple(parts)

def histo_key(pd, version, process, selection, hist_name):
    return f"{pd}__{version}__{process}__{selection}__{hist_name}";

def add_hist(hist_full_name, hist):
    (pd, version, process, selection, hist_name) = unpack_histo_key(hist_full_name)
    if process == None:
        return
    histos[pd][version][process][selection][hist_name] = hist

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
                # if not re.search(".*__jpsi__loose__", hist.GetName()):
                #     continue
                # if re.search("vtx_prob", hist.GetName()):
                #     continue
                # if re.search("(gt|lt)\d", hist.GetName()):
                #     continue
                f.Remove(hist)
                add_hist(hist.GetName(), hist)
    finally:
        f.Close()

    return True

def book_histo1D(rdf, pd, version, process, selection, hist_name,
                 hist_title, nbins, xmin, xmax, var):
    hist_full_name = histo_key(pd, version, process, selection, hist_name)
    rdf_hist = rdf.Histo1D((hist_full_name, hist_title, nbins, xmin, xmax), var)
    add_hist(hist_full_name, rdf_hist)
    
def book_histo2D(rdf, pd, version, process, selection, hist_name, hist_title,
                 nbinsx, xmin, xmax, nbinsy, ymin, ymax, varx, vary):
    hist_full_name = histo_key(pd, version, process, selection, hist_name)
    rdf_hist = rdf.Histo2D((hist_full_name, hist_title,
                            nbinsx, xmin, xmax,
                            nbinsy, ymin, ymax),
                           varx, vary)
    add_hist(hist_full_name, rdf_hist)

def book_histos(rdf, pd, version, trigger=None):
    # Add defaults
    rdf = rdf.DefaultValueFor("HLT_DoubleMu4_3_LowMass", False)
    rdf = rdf.DefaultValueFor("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6", False)
    
    rdf = rdf.Define("certified", "passed_lumi_mask(run, luminosityBlock)")
    rdf = rdf.Filter("certified == 1", "passed data certification")
    rdf = rdf.Define("mm_mu1_mediumId", "Take(Muon_mediumId,             mm_mu1_index)")
    rdf = rdf.Define("mm_mu2_mediumId", "Take(Muon_mediumId,             mm_mu2_index)")
        
    if trigger != None:
        rdf = rdf.Filter(trigger)

    #### Jpsi
    jpsi_loose_selection =   f"mm_kin_mass>{min_jpsi_mass} && mm_kin_mass<{max_jpsi_mass} && mm_mu1_mediumId && mm_mu2_mediumId && mm_kin_vtx_prob>0.001"
    jpsi_nominal_selection = f"{jpsi_loose_selection} && mm_kin_vtx_prob>0.1"

    ## Loose
    # Selection
    rdf = rdf.Define("jpsi_loose",               jpsi_loose_selection)
    rdf = rdf.Define("jpsi_loose_lt0p8",         f"{jpsi_loose_selection} && abs(mm_kin_eta)<0.8")
    rdf = rdf.Define("jpsi_loose_gt0p8_lt1p4",   f"{jpsi_loose_selection} && abs(mm_kin_eta)>0.8 && abs(mm_kin_eta)<1.4")
    rdf = rdf.Define("jpsi_loose_gt1p4",         f"{jpsi_loose_selection} && abs(mm_kin_eta)>1.4")
    
    for sel in ['loose', 'loose_lt0p8', 'loose_gt0p8_lt1p4', 'loose_gt1p4']:
        # Variables
        rdf = rdf.Define(f"jpsi_{sel}_mass",     f"mm_kin_mass[jpsi_{sel}]")
        rdf = rdf.Define(f"jpsi_{sel}_vtx_prob", f"mm_kin_vtx_prob[jpsi_{sel}]")
        # Book histograms
        book_histo1D(rdf, pd, version, "jpsi", sel, "mass", "Vertex constrained dimuon mass;Mass, GeV",
                     nbins_jpsi, min_jpsi_mass, max_jpsi_mass, f"jpsi_{sel}_mass")
        book_histo2D(rdf, pd, version, "jpsi", sel, "vtx_prob_vs_mass", "Jpsi;Mass, GeV;Probability",
                     nbins_jpsi, min_jpsi_mass, max_jpsi_mass, 50, 0.0, 1.0,
                     f"jpsi_{sel}_mass", f"jpsi_{sel}_vtx_prob")
    ## Nominal
    # Selection
    rdf = rdf.Define("jpsi_nominal",             jpsi_nominal_selection)
    rdf = rdf.Define("jpsi_nominal_lt0p8",       f"{jpsi_nominal_selection} && abs(mm_kin_eta)<0.8")
    rdf = rdf.Define("jpsi_nominal_gt0p8_lt1p4", f"{jpsi_nominal_selection} && abs(mm_kin_eta)>0.8 && abs(mm_kin_eta)<1.4")
    rdf = rdf.Define("jpsi_nominal_gt1p4",       f"{jpsi_nominal_selection} && abs(mm_kin_eta)>1.4")
    
    for sel in ['nominal', 'nominal_lt0p8', 'nominal_gt0p8_lt1p4', 'nominal_gt1p4']:
        # Variables
        rdf = rdf.Define(f"jpsi_{sel}_mass",     f"mm_kin_mass[jpsi_{sel}]")
        rdf = rdf.Define(f"jpsi_{sel}_vtx_prob", f"mm_kin_vtx_prob[jpsi_{sel}]")
        # Variables
        book_histo1D(rdf, pd, version, "jpsi", sel, "mass", "Vertex constrained dimuon mass;Mass, GeV",
                     nbins_jpsi, min_jpsi_mass, max_jpsi_mass, f"jpsi_{sel}_mass")
        book_histo2D(rdf, pd, version, "jpsi", sel, "vtx_prob_vs_mass", "Jpsi;Mass, GeV;Probability",
                     nbins_jpsi, min_jpsi_mass, max_jpsi_mass, 18, 0.1, 1.0,
                     f"jpsi_{sel}_mass", f"jpsi_{sel}_vtx_prob")

    #### BuToJpsiK
    if pd == "ParkingDoubleMuonLowMass0":
        jpsik_loose_selection = f"bkmm_jpsimc_vtx_prob>0.01 && bkmm_jpsimc_sl3d>3 && " \
            f"abs(bkmm_jpsimc_alpha)<0.1 && bkmm_jpsimc_mass>{min_jpsik_mass} && bkmm_jpsimc_mass<{max_jpsik_mass} && bkmm_kaon_pt>3"
        jpsik_tight_selection = f"bkmm_jpsimc_vtx_prob>0.1 && bkmm_jpsimc_sl3d>5 && abs(bkmm_jpsimc_alpha)<0.01 && {jpsik_loose_selection}"

        ## Loose
        # Selection
        rdf = rdf.Define("jpsik_loose",               jpsik_loose_selection)
        rdf = rdf.Define("jpsik_loose_lt0p8",         f"{jpsik_loose_selection} && abs(bkmm_jpsimc_eta)<0.8")
        rdf = rdf.Define("jpsik_loose_gt0p8_lt1p4",   f"{jpsik_loose_selection} && abs(bkmm_jpsimc_eta)>0.8 && abs(bkmm_jpsimc_eta)<1.4")
        rdf = rdf.Define("jpsik_loose_gt1p4",         f"{jpsik_loose_selection} && abs(bkmm_jpsimc_eta)>1.4")
    
        for sel in ['loose', 'loose_lt0p8', 'loose_gt0p8_lt1p4', 'loose_gt1p4']:
            # Variables
            rdf = rdf.Define(f"jpsik_{sel}_mass",     f"bkmm_jpsimc_mass[jpsik_{sel}]")
            rdf = rdf.Define(f"jpsik_{sel}_pvip",     f"bkmm_jpsimc_pvip[jpsik_{sel}]")
            rdf = rdf.Define(f"jpsik_{sel}_pvlip",    f"bkmm_jpsimc_pvlip[jpsik_{sel}]")
            rdf = rdf.Define(f"jpsik_{sel}_alpha",    f"bkmm_jpsimc_alpha[jpsik_{sel}]")
            rdf = rdf.Define(f"jpsik_{sel}_alphaBS",  f"bkmm_jpsimc_alphaBS[jpsik_{sel}]")
            rdf = rdf.Define(f"jpsik_{sel}_vtx_prob", f"bkmm_jpsimc_vtx_prob[jpsik_{sel}]")
            # Book histograms
            book_histo1D(rdf, pd, version, "jpsik", sel, "mass", "BtoJpsiK;Mass, GeV",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, f"jpsik_{sel}_mass")
            book_histo2D(rdf, pd, version, "jpsik", sel, "pvip_vs_mass", "BtoJpsiK;Mass, GeV;IP_{PV}",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, 50, 0, 0.02,
                         f"jpsik_{sel}_mass", f"jpsik_{sel}_pvip")
            book_histo2D(rdf, pd, version, "jpsik", sel, "pvlip_vs_mass", "BtoJpsiK;Mass, GeV;IP_{PV,long}",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, 50, 0, 0.02,
                         f"jpsik_{sel}_mass", f"jpsik_{sel}_pvlip")
            book_histo2D(rdf, pd, version, "jpsik", sel, "alpha_vs_mass", "BtoJpsiK;Mass, GeV;#alpha_{3D}",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, 50, 0, 0.1,
                         f"jpsik_{sel}_mass", f"jpsik_{sel}_alpha")
            book_histo2D(rdf, pd, version, "jpsik", sel, "alphaBS_vs_mass", "BtoJpsiK;Mass, GeV;#alpha_{BS}",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, 50, 0, 0.1,
                         f"jpsik_{sel}_mass", f"jpsik_{sel}_alphaBS")
            book_histo2D(rdf, pd, version, "jpsik", sel, "vtx_prob_vs_mass", "BtoJpsiK;Mass, GeV;Probability",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, 50, 0.0, 1.0,
                         f"jpsik_{sel}_mass", f"jpsik_{sel}_vtx_prob")
        ## Tight
        # Selection
        rdf = rdf.Define("jpsik_tight",               jpsik_tight_selection)
        rdf = rdf.Define("jpsik_tight_lt0p8",         f"{jpsik_tight_selection} && abs(bkmm_jpsimc_eta)<0.8")
        rdf = rdf.Define("jpsik_tight_gt0p8_lt1p4",   f"{jpsik_tight_selection} && abs(bkmm_jpsimc_eta)>0.8 && abs(bkmm_jpsimc_eta)<1.4")
        rdf = rdf.Define("jpsik_tight_gt1p4",         f"{jpsik_tight_selection} && abs(bkmm_jpsimc_eta)>1.4")
        for sel in ['tight', 'tight_lt0p8', 'tight_gt0p8_lt1p4', 'tight_gt1p4']:
            # Variables
            rdf = rdf.Define(f"jpsik_{sel}_mass",     f"bkmm_jpsimc_mass[jpsik_{sel}]")
            book_histo1D(rdf, pd, version, "jpsik", sel, "mass", "BtoJpsiK;Mass, GeV",
                         nbins_jpsik, min_jpsik_mass, max_jpsik_mass, f"jpsik_{sel}_mass")

        ### BuToJpsiPhi
        rdf = rdf.Define("jpsiphi_cands", f"bkkmm_jpsikk_vtx_prob>0.1 && bkkmm_jpsikk_sl3d>5 && "
                         f"abs(bkkmm_jpsikk_alpha)<0.01 && bkkmm_jpsikk_mass>{min_jpsiphi_mass} && "
                         f"bkkmm_jpsikk_mass<{max_jpsiphi_mass} && abs(bkkmm_kk_mass-1.02)<0.01")
        rdf = rdf.Define("jpsiphi_mass","bkkmm_jpsikk_mass[jpsiphi_cands]")

        book_histo1D(rdf, pd, version, "jpsiphi", "loose", "mass", "BtoJpsiPhi;Mass, GeV",
                     nbins_jpsiphi, min_jpsiphi_mass, max_jpsiphi_mass, "jpsiphi_mass")
    
    h_list = []
    for process in histos[pd][version]:
        for selection in histos[pd][version][process]:
            for h_name, h in histos[pd][version][process][selection].items():
                # h.SetLineColor(ROOT.kBlack)
                # h.SetLineWidth(2)
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


def build_model(workspace_name, mass_var, peak=3.09, search_width=0.04, max_exp_c=0):
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
    
    # sigmaG_John = ROOT.RooRealVar(n_+"sigmaG_John"," sigma ",0.01, 0.001, 1.)
    # gaus_John = ROOT.RooGaussian(n_+"gaus_John","", m, sig_mu, sigmaG_John)
    # JohnG_frac = ROOT.RooRealVar(n_+"JohnG_frac","",0.3,0.,1.0)
    
    # sig_mu     = ROOT.RooRealVar("sig_mu", "mu", peak, peak - search_width, peak + search_width)
    # sig_lambda = ROOT.RooRealVar("sig_lambda", "lambda", 0.03, 0.001, 0.10)
    # sig_gamma  = ROOT.RooRealVar("sig_gamma", "gamma", 1, 0, 5)
    # sig_delta  = ROOT.RooRealVar("sig_delta", "delta", 1, 0.1, 10)
    # sig = ROOT.RooJohnson("sig", "signal", mass_var, sig_mu, sig_lambda, sig_gamma, sig_delta)

    # multi-gaussian
    G1_mean  = ROOT.RooRealVar("sig_G1_mean",  "", peak, peak - search_width, peak + search_width)
    G1_sigma = ROOT.RooRealVar("sig_G1_sigma", "", 0.03, 0.001, 0.10)
    G2_scale = ROOT.RooRealVar("sig_G2_scale", "", 1.5, 1.1, 3.5)
    G3_scale = ROOT.RooRealVar("sig_G3_scale", "", 3.0, 1.0, 6.7)
    G2_sigma = ROOT.RooProduct("sig_G2_sigma", "", ROOT.RooArgList(G1_sigma,G2_scale))
    G3_sigma = ROOT.RooProduct("sig_G3_sigma", "", ROOT.RooArgList(G1_sigma,G3_scale))
    G1 = ROOT.RooGaussian("sig_G1", "", mass_var, G1_mean, G1_sigma)
    G2 = ROOT.RooGaussian("sig_G2", "", mass_var, G1_mean, G2_sigma)
    G3 = ROOT.RooGaussian("sig_G3", "", mass_var, G1_mean, G3_sigma)
    
    G2_fract = ROOT.RooRealVar("sig_G2_fract","",0.2,0.0,0.5)
    G3_fract = ROOT.RooRealVar("sig_G3_fract","",0.2,0.0,1.0)
    sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_fract))
    # sig  = ROOT.RooAddPdf("sig"," ",ROOT.RooArgList(G3,G2,G1),ROOT.RooArgList(G2_fract,G3_fract))

    # # CB
    # sig_mean  = ROOT.RooRealVar("sig_mean",  "", peak, peak - search_width, peak + search_width)
    # sig_sigma = ROOT.RooRealVar("sig_sigma", "sigma", 0.03, 0.001, 0.10)
    # sig_tail  = ROOT.RooRealVar("sig_tail",  "tail", 2.8, 0.1, 10.0)
    # sig_pow   = ROOT.RooRealVar("sig_pow0",  "pow", 3, 1, 50)
    # # sig = ROOT.RooCBShape("sig", "signal", mass_var, sig_mean, sig_sigma, sig_tail, sig_pow)
    # sig_cb = ROOT.RooCBShape("sig_cb", "signal", mass_var, sig_mean, sig_sigma, sig_tail, sig_pow)
    
    # G2_scale = ROOT.RooRealVar("sig_G2_scale", "", 1.5, 1.0, 2.5)
    # G3_scale = ROOT.RooRealVar("sig_G3_scale", "", 3.0, 0.5, 6.7)
    # G2_sigma = ROOT.RooProduct("sig_G2_sigma", "", ROOT.RooArgList(sig_sigma,G2_scale))
    # G3_sigma = ROOT.RooProduct("sig_G3_sigma", "", ROOT.RooArgList(sig_sigma,G3_scale))
    # G2 = ROOT.RooGaussian("sig_G2", "", mass_var, sig_mean, G2_sigma)
    # G3 = ROOT.RooGaussian("sig_G3", "", mass_var, sig_mean, G3_sigma)
    
    # G2_fract = ROOT.RooRealVar("sig_G2_fract","",0.1,0.0,0.3)
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
    
    exp_c = ROOT.RooRealVar("exp_c","exp_c", max_exp_c, -1000., max_exp_c)
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
        results = defaultdict(                 # pd
            lambda: defaultdict(               # process
                lambda: defaultdict(           # selection
                    lambda: defaultdict(dict)  # hist_name -> {version: hist}
                )
            )
        )
        for version in histos[pd]:
            for process in histos[pd][version]:
                for selection, histograms in histos[pd][version][process].items():
                    if not "mass" in histograms:
                        # skip it for now
                        continue
            
                    ### Fit mass distribution to get model parameters

                    h_ref = histograms["mass"]
                    results[pd][process][selection]["mass"][version] = h_ref

                    if process == "jpsik":
                        ws = ws_jpsik
                        mass = rrv_jpsik_mass
                        mass_binning = rrv_jpsik_mass_binning
                        nbins = nbins_jpsik
                    elif process == "jpsiphi":
                        ws = ws_jpsiphi
                        mass = rrv_jpsiphi_mass
                        mass_binning = rrv_jpsiphi_mass_binning
                        nbins = nbins_jpsiphi
                    else:
                        ws = ws_jpsi
                        mass = rrv_jpsi_mass
                        mass_binning = rrv_jpsi_mass_binning
                        nbins = nbins_jpsi

                    model = ws.pdf("model")

                    data = ROOT.RooDataHist("data", "", ROOT.RooArgList(mass), h_ref)

                    # prefit
                    ws.var("Nsig").setVal(h_ref.GetEntries() * 0.99)
                    ws.var("Nbkg").setVal(h_ref.GetEntries() * 0.01)
                    # ws.var("sig_G2_fract").setVal(0.2)
                    # ws.var("sig_G2_fract").setConstant(True)
                    # ws.var("sig_G2_scale").setConstant(True)
                    # ws.var("sig_G3_fract").setVal(0.0)
                    # ws.var("sig_G3_fract").setConstant(True)
                    # ws.var("sig_G3_scale").setConstant(True)
                    # model.fitTo(data,  ROOT.RooFit.NumCPU(8),
                    #             ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
                    #             ROOT.RooFit.PrintLevel(print_level))

                    # final fit
                    # ws.var("sig_G2_fract").setConstant(False)
                    # ws.var("sig_G2_scale").setConstant(False)
                    # ws.var("sig_G3_fract").setConstant(False)
                    # ws.var("sig_G3_scale").setConstant(False)
                    fr = model.fitTo(
                        data,
                        ROOT.RooFit.NumCPU(8),
                        ROOT.RooFit.Extended(ROOT.kTRUE),
                        # ROOT.RooFit.Minos(ROOT.kFALSE),
                        ROOT.RooFit.Minos(ROOT.kTRUE),
                        ROOT.RooFit.PrintLevel(print_level),
                        ROOT.RooFit.Save(True)
                    )
                    fr.Print("V")
                    #     h.Draw("hist")
                    #     print_canvas(f"{sample_name}_{name}", f"{output_path}/{study}/histograms")


                    ### Plot fit results

                    frame = mass.frame()
                    data.plotOn(frame, mass_binning)
                    frame.SetMaximum(frame.GetMaximum() * 1.2)
                    # model.plotOn(frame, ROOT.RooFit.Components("sig"), ROOT.RooFit.LineColor(ROOT.kRed))
                    model.plotOn(frame, ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))
                    model.plotOn(frame)
                    print("chiSquare: ", frame.chiSquare(6))
                    print("chiSquare: ", frame.chiSquare("model","data", 6))

                    model.paramOn(frame, ROOT.RooFit.Layout(0.7, 0.95, 0.92))
                    frame.getAttText().SetTextSize(0.02)
                    frame.Draw()
                    print_canvas(f"{pd}_{version}_{process}_{selection}", f"{output_path}/fits")

                    # results[name][sample_name] = (ws.var("Nsig").getVal(), ws.var("Nsig").getError())

                    ### Compute sWeights for sPlots

                    # Covariance matrix for the yields only
                    cov = fr.reducedCovarianceMatrix(ROOT.RooArgList(ws.var("Nsig"), ws.var("Nbkg")))
                    # print(f"cov[0][0]: {cov[0][0]}")
                    # print(f"cov[0][1]: {cov[0][1]}")

                    # Normally sWeights are computed per event as a function
                    # of mass. In our case instead of events we use 2D
                    # histograms and sWeights are computed as a function of
                    # the the mass bin instead.

                    # Compute binned PDFs
                    Fsig = [0.0] * nbins
                    Fbkg = [0.0] * nbins
                    sig_pdf = ws.pdf("sig")
                    bkg_pdf = ws.pdf("bkg")

                    for ib in range(1, nbins + 1):
                        xl = h_ref.GetXaxis().GetBinLowEdge(ib)
                        xh = h_ref.GetXaxis().GetBinUpEdge(ib)
                        bin_name = f"bin_{ib}"
                        mass.setRange(bin_name, xl, xh)

                        # Integral of pdf over this bin (normalization handled by NormSet)
                        isig = sig_pdf.createIntegral(
                            ROOT.RooArgSet(mass),
                            ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)),
                            ROOT.RooFit.Range(bin_name)
                        ).getVal()

                        ibkg = bkg_pdf.createIntegral(
                            ROOT.RooArgSet(mass),
                            ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)),
                            ROOT.RooFit.Range(bin_name)
                        ).getVal()

                        # You can keep them as integrals (no need to divide by bin width)
                        Fsig[ib - 1] = isig
                        Fbkg[ib - 1] = ibkg

                    # print(Fsig)
                    # print(Fbkg)
                    
                    # Compute binned sWeights for signal
                    wSig = [0.0] * nbins

                    NS = ws.var("Nsig").getVal()
                    NB = ws.var("Nbkg").getVal()
                    # print(f"NS: {NS}")
                    # print(f"NB: {NB}")

                    for ib in range(nbins):
                        denom = NS * Fsig[ib] + NB * Fbkg[ib]
                        if denom <= 0:
                            wSig[ib] = 0.0
                            continue
                        numSig = cov[0][0] * Fsig[ib] + cov[0][1] * Fbkg[ib]
                        wSig[ib] = numSig / denom

                    # print(wSig)
                    
                    ### Make plots

                    for histo_name, hist in histograms.items():
                        if hist.InheritsFrom("TH2"):
                            
                            ### Make projections using 1-sigma signal window
                            
                            mean  = ws.var("sig_G1_mean").getVal()
                            sigma = ws.var("sig_G1_sigma").getVal()
                            xaxis = hist.GetXaxis()
                            bin_x_min = xaxis.FindBin(mean - sigma)
                            bin_x_max = xaxis.FindBin(mean + sigma)
                            hProjY = hist.ProjectionY(f"{hist.GetName()}_hProjY", bin_x_min, bin_x_max)
                            hProjY.SetMinimum(0)
                            hProjY.SetLineWidth(2)
                            hProjY.Draw()
                            
                            print_canvas(f"{pd}_{version}_{process}_{selection}_{histo_name}", f"{output_path}/projections")
                            results[pd][process][selection][f"{histo_name}_projection"][version] = hProjY
                            
                            ### Make sPlots
                            
                            # skip vtx_prob - doesn't look right
                            if re.search("vtx_prob", hist.GetName()):
                                continue
                            hSplotY = hist.ProjectionY(f"{hist.GetName()}_hSplotY")
                            hSplotY.SetDirectory(0)
                            hSplotY.Reset("ICE")
                            hSplotY.Sumw2()

                            # Build sWeighted y histogram
                            for ix in range(1, nbins + 1):
                                w = wSig[ix - 1]           # sWeight for this x-bin
                                if w == 0:
                                    continue

                                for iy in range(1, hist.GetNbinsY() + 1):
                                    n = hist.GetBinContent(ix, iy)
                                    if n <= 0:
                                        continue

                                    # Add weighted content
                                    old = hSplotY.GetBinContent(iy)
                                    hSplotY.SetBinContent(iy, old + w * n)

                                    # Poisson variance propagation: var += w^2 * n
                                    old_err2 = hSplotY.GetBinError(iy)**2
                                    new_err2 = old_err2 + (w * w) * n
                                    hSplotY.SetBinError(iy, sqrt(new_err2))
                                    
                            hSplotY.Draw()
                            hSplotY.Print("V")
                            print_canvas(f"{pd}_{version}_{process}_{selection}_{histo_name}", f"{output_path}/splots")
                            results[pd][process][selection][histo_name][version] = hSplotY

        ### Ratio plots
        ratio_plots = list()
        for process in results[pd]:
            for selection in results[pd][process]:
                for histo_name, histograms in results[pd][process][selection].items():
                    if "rereco" not in histograms or "prompt" not in histograms:
                        print("ERROR: missing histograms")
                        for n, h in histograms.items():
                            print(f"{n}: {h.GetName()}")
                        continue
                    histograms["rereco"].SetMarkerStyle(20)
                    # histograms["rereco"].SetMinimum(0)
                    histograms["prompt"].SetLineColor(ROOT.kRed)
                    histograms["prompt"].SetLineWidth(2)
                    # histograms["prompt"].SetMinimum(0)

                    ratio_plot = ROOT.TRatioPlot(histograms["rereco"], histograms["prompt"])
                    ratio_plots.append(ratio_plot) # keep it memory to avoid seg fault
                    ratio_plot.SetH1DrawOpt("e")
                    ratio_plot.SetH2DrawOpt("hist")
                    ratio_plot.Draw()
                    ratio_plot.SetSeparationMargin(0.03)
                    if pd == "ParkingDoubleMuonLowMass0" and process == "jpsi":
                        ratio_plot.GetLowerRefGraph().SetMinimum(0.8)
                        ratio_plot.GetLowerRefGraph().SetMaximum(1.2)
                    else:
                        ratio_plot.GetLowerRefGraph().SetMinimum(0.5)
                        ratio_plot.GetLowerRefGraph().SetMaximum(2.0)
                    # ratio_plot.GetXaxis().SetTitleSize()
                    # SetBottomMargin(2.0)           
                    # c1.SetBottomMargin(2.0)
                    # rp->GetLowerRefYaxis()->SetRange(...)
                    # rp->SetH1DrawOpt("E");

                    ratio_plot.GetLowerRefYaxis().SetTitle("rereco/prompt")

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

                    c.Update()
                    c.cd()
                    # right position
                    legend = ROOT.TLegend(0.70,0.75,0.85,0.85)
                    # left position
                    # legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
                    legend.SetFillStyle(0)
                    legend.SetLineWidth(0)
                    legend.SetBorderSize(1)

                    legend.AddEntry(histograms["rereco"], "Rereco", "p")
                    legend.AddEntry(histograms["prompt"], "Prompt", "l")
                    legend.Draw()

                    if histo_name == "mass":
                        print_canvas(f"{pd}_{process}_{selection}_{histo_name}", f"{output_path}/simple_ratios")
                    elif re.search("hSplotY", histograms["prompt"].GetName()):
                        print_canvas(f"{pd}_{process}_{selection}_{histo_name}", f"{output_path}/splot_ratios")
                    elif re.search("hProjY", histograms["prompt"].GetName()):
                        print_canvas(f"{pd}_{process}_{selection}_{histo_name}", f"{output_path}/proj_ratios")
                    else:
                        print_canvas(f"{pd}_{process}_{selection}_{histo_name}", f"{output_path}")
            

def add_files(chain, path):
    nfiles = 0
    for root_dir, _, files in os.walk(path):
        for file in files:
            if re.search(f"{pattern}", file):
                full_path = os.path.join(root_dir, file)
                # print(full_path)
                nfiles += 1
                chain.Add(full_path)
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

def intersect_by_run(a_by_run, b_by_run):
    out = {}
    common_runs = set(a_by_run) & set(b_by_run)
    for r in common_runs:
        lumis = a_by_run[r] & b_by_run[r]
        if lumis:
            out[r] = lumis
    return out
    
def get_run_to_lumi_set(chain):
    """Extract set of lumis present for each run."""
    run_to_lumi_set = {}
    nlumis = 0
    for entry in chain: 
        run = int(entry.run)
        lumi = int(entry.luminosityBlock)
        run_to_lumi_set.setdefault(run, set()).add(lumi)
        nlumis += 1
    print(f"Total number of lumi sections: {nlumis}")
    return run_to_lumi_set

def get_common_lumis(lumi_set1, lumi_set2):
    """Intersect lumis per run where run exists in both and return merged run-lumi ranges"""
    
    common_run_to_lumi_set = {}
    common_runs = set(lumi_set1) & set(lumi_set2)
    ncommon_lumis = 0
    for run in common_runs:
        lumis = lumi_set1[run] & lumi_set2[run]
        if lumis:
            common_run_to_lumi_set[run] = lumis
            ncommon_lumis += len(lumis)
    print(f"Total number of common lumi sections: {ncommon_lumis}")
    
    # Convert sets to sorted ranges
    run_to_lumi_ranges = {}
    for run, lumis in common_run_to_lumi_set.items():
        sorted_lumis = sorted(lumis)
        run_to_lumi_ranges[str(run)] = merge_into_ranges(sorted_lumis)
    return run_to_lumi_ranges

def load_lumi_mask(file):
    # Load lumi masks
    lumi_mask_string = json_to_custom_format(file)

    # Check if conversion succeeded
    if lumi_mask_string:
        # Create the LumiMask object using the custom string format
        lumi_mask = ROOT.LumiMask.fromCustomString(lumi_mask_string, 0, 0)
    else:
        raise Exception("Failed to convert JSON to custom string format.")

    # ROOT.gInterpreter.ProcessLine(f'lumi_mask_string = "{lumi_mask_string}";')
    ROOT.gInterpreter.ProcessLine(f'set_lumi_mask("{lumi_mask_string}");')

def save_histograms():
    f = ROOT.TFile(histos_file, "RECREATE")
    for pd in histos:
        for version in histos[pd]:
            for process in histos[pd][version]:
                for selection in histos[pd][version][process]:
                    for histo_name, hist in histos[pd][version][process][selection].items():
                        print(f"Saving {hist.GetName()}")
                        hist.Write(hist.GetName())
    f.Close()

    
####################################################################################

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

c = ROOT.TCanvas("c", "", 900, 900)

ws_jpsi = build_model("jpsi", rrv_jpsi_mass)
ws_jpsik = build_model("jpsik", rrv_jpsik_mass, 5.3, 0.05, -3)
ws_jpsiphi = build_model("jpsiphi", rrv_jpsiphi_mass, 5.4, 0.05)

data = {
    'ParkingDoubleMuonLowMass0':{
        'prompt': [
            f"{input_path}/ParkingDoubleMuonLowMass0+Run2025C-PromptReco-v1+MINIAOD/",
            f"{input_path}/ParkingDoubleMuonLowMass0+Run2025C-PromptReco-v2+MINIAOD/",
        ],
        'rereco': [
            f"{input_path}/ParkingDoubleMuonLowMass0+Run2025C-TrkRadDamage-v1+MINIAOD/"
        ]            
    },
    'EGamma0':{
        'prompt': [
            f"{input_path}/EGamma0+Run2025C-PromptReco-v1+MINIAOD/",
            f"{input_path}/EGamma0+Run2025C-PromptReco-v2+MINIAOD/",
        ],
        'rereco': [
            f"{input_path}/EGamma0+Run2025C-TrkRadDamage-v2+MINIAOD/"
        ]            
    },
}

# Declare the LumiMask class in ROOT using the header file
with open('LumiMask.h', 'r') as file:
    lumi_mask_code = file.read()
ROOT.gInterpreter.Declare(lumi_mask_code)

ROOT.gInterpreter.Declare(f'''
// std::string lumi_mask_string;
static LumiMask lumi_mask(std::vector<LumiMask::LumiBlockRange>{{}});

void set_lumi_mask(std::string lumi_mask_string) {{
    lumi_mask = LumiMask::fromCustomString(lumi_mask_string);
}}

bool passed_lumi_mask(unsigned int run, unsigned int lumi) {{
    return lumi_mask.accept(run, lumi);
}}
''')

if create_lumi_mask:
    print("Creating lumi masks")

    for pd in data:
    
        # Find common lumi blocks
        chain_lumi_prompt = ROOT.TChain("LuminosityBlocks")
        for path in data[pd]["prompt"]:
            add_files(chain_lumi_prompt, path)
        lumi_set1 = get_run_to_lumi_set(chain_lumi_prompt)
    
        chain_lumi_rereco = ROOT.TChain("LuminosityBlocks")
        for path in data[pd]["rereco"]:
            add_files(chain_lumi_rereco, path)
        lumi_set2 = get_run_to_lumi_set(chain_lumi_rereco)
    
        lumi_mask = get_common_lumis(lumi_set1, lumi_set2)
        lumi_mask_file = f"rdf_re-reco_validation-{pd}.json"

        with open(lumi_mask_file, "w") as f:
            json.dump(lumi_mask, f, indent=4, sort_keys=True)

if process_data:

    for pd in data:
        trigger = None
        if pd == "ParkingDoubleMuonLowMass0":
            trigger = "HLT_DoubleMu4_3_LowMass"
        print("Loading lumi mask")
        load_lumi_mask(f"rdf_re-reco_validation-{pd}.json")
    
        print("Processing data")
        chain_prompt = ROOT.TChain("Events")
        for path in data[pd]["prompt"]:
            add_files(chain_prompt, path)
        rdf = ROOT.RDataFrame(chain_prompt)
        n = book_histos(rdf, pd, "prompt", trigger)
        print(f"n: {n}")
    
        chain_rereco = ROOT.TChain("Events")
        for path in data[pd]["rereco"]:
            add_files(chain_rereco, path)
        rdf = ROOT.RDataFrame(chain_rereco)
        n = book_histos(rdf, pd, "rereco", trigger)
        print(f"n: {n}")

    save_histograms()

load_histos()
make_plots()
