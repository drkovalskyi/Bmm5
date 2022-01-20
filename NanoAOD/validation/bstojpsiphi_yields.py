"""Bmm5 BstoJpsiPhi normalization fits

Unbinned maximum likelihood fit using fit-bkkmm flat ntuples. The
signal shape is a non-parameteric PDF build from the simulated JpsiPhi
mass distribution convoluted with a Gaussian. The JpsiK* pdf is also a
non-parameteric pdf, but without convolution. The JpsiK* yield is
fixed to the signal yields and the coefficient is automatically
computed based on the branching fractions of the two processes and
their relative efficiency.

The fits are heavy and memory demanding due to a large number of
events in each Era for each category. If the script fails running out
of memomory you can re-run it for one Era at a time.

"""
import os, re, copy, json
import ROOT
import math 
import tdrstyle
from ROOT.RooFit import Binning
from ROOT import RooRealVar

# Set the TDR style
tdrstyle.setTDRStyle()

version = 518

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/bsjpsiphi-yields/";

recompute_results = True
result_file = output_path + "/results.json"
results = dict()
if os.path.exists(result_file):
    results = json.load(open(result_file))
    
min_mass = 5.25
max_mass = 5.60
mass = ROOT.RooRealVar("m", "",  min_mass, max_mass)
nbins = int((max_mass-min_mass)/0.005)
mass_ref_binning  = Binning(2 * nbins, min_mass, max_mass)
mass_plot_binning = Binning(nbins, min_mass, max_mass)
chan = RooRealVar("chan","", 0)

jpsikstar_bf = 1.27e-3 * 2 / 3
jpsiphi_bf   = 1.08e-3 * 0.492
fsfd_mc = 0.28
fsfd = 0.25

# Silence RooFit info messages
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

# RooMinimizer::PrintLevel
# None =-1 , Reduced =0 , Normal =1 , ExtraForProblem =2 ,  Maximum =3
# Default 1
print_level = 0

# mc_match = ROOT.RooCategory("mc_match", "mc_match")
# mc_match.defineType("B+", 521)
# mc_match.defineType("B-", -521)

# variables = [
#     ( RooRealVar("HLT_DoubleMu4_3_Jpsi","", 0),                                None,                   None   ), # UInt_t is not supported as input for RooCategory
#     ( RooRealVar("HLT_DoubleMu4_3_Jpsi_Displaced","", 0),                      None,                   None   ), # UInt_t is not supported as input for RooCategory
#     ( RooRealVar("chan","", 0),                                                None,                   None   ), # UInt_t is not supported as input for RooCategory
#     ( RooRealVar("HLT_DoubleMu4_3_Jpsi_ps","", 0),                             None,                   None   ), 
#     ( RooRealVar("HLT_DoubleMu4_3_Jpsi_Displaced_ps","", 0),                   None,                   None   ), 
#     ( m1q,                                                                     None,                   None   ), 
#     ( m2q,                                                                     None,                   None   ), 
#     ( RooRealVar("m1bdt", "", 0),                                              None,                   None   ), 
#     ( RooRealVar("m2bdt", "", 0),                                              None,                   None   ), 
# ]

# roovars = []
# for var in variables:
#     roovars.append(var[0])
# roovars_mc = copy.deepcopy(roovars)
# roovars_mc.append(mc_match)

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/%u/fit-bkkmm/" % version
mc_jpsiphi_2018   = data_path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*root"
mc_jpsiphi_2017   = data_path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*root"
mc_jpsiphi_2016BF = data_path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*root"
mc_jpsiphi_2016GH = data_path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"

mc_jpsikstar_2018   = data_path + "/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*root"
mc_jpsikstar_2017   = data_path + "/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*root"
mc_jpsikstar_2016BF = data_path + "/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*root"
mc_jpsikstar_2016GH = data_path + "/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*root"
    
datasets = {

    "Run2018":{
        "mc_jpsiphi": mc_jpsiphi_2018,
        "mc_jpsiphi_tree": "bspsiphiMc",
        "mc_jpsikstar": mc_jpsikstar_2018,
        "mc_jpsikstar_tree": "bdpsikstarMc",
        "data_tree": "bspsiphiData",
        "trigger": "HLT_DoubleMu4_3_Jpsi",
        "fix_resolution": True,
        "data": [
            data_path + "/Charmonium+Run2018A-12Nov2019_UL2018_rsb-v1+MINIAOD/*.root",
            data_path + "/Charmonium+Run2018B-12Nov2019_UL2018-v1+MINIAOD/*.root",
            data_path + "/Charmonium+Run2018C-12Nov2019_UL2018_rsb_v3-v1+MINIAOD/*.root",
            data_path + "/Charmonium+Run2018D-12Nov2019_UL2018-v1+MINIAOD/*.root",
        ],
    },
    # "Run2017":{
    #     "mc_jpsiphi": mc_jpsiphi_2017,
    #     "mc_jpsiphi_tree": "bspsiphiMc",
    #     "mc_jpsikstar": mc_jpsikstar_2017,
    #     "mc_jpsikstar_tree": "bdpsikstarMc",
    #     "data_tree": "bspsiphiData",
    #     "trigger": "HLT_DoubleMu4_3_Jpsi_Displaced",
    #     "fix_resolution": True,
    #     "data": [
    #         data_path + "/Charmonium+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
    #         data_path + "/Charmonium+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
    #         data_path + "/Charmonium+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
    #         data_path + "/Charmonium+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
    #         data_path + "/Charmonium+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root"
    #     ]
    # },
    # "Run2016BF":{
    #     "mc_jpsiphi": mc_jpsiphi_2016BF,
    #     "mc_jpsiphi_tree": "bspsiphiMc",
    #     "mc_jpsikstar": mc_jpsikstar_2016BF,
    #     "mc_jpsikstar_tree": "bdpsikstarMc",
    #     "data_tree": "bspsiphiData",
    #     "trigger": "HLT_DoubleMu4_3_Jpsi_Displaced",
    #     "fix_resolution": True,
    #     "data": [
    #         data_path + "Charmonium+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016F-21Feb2020_UL2016-v1+MINIAOD/*.root"
    #     ]
    # },
    # "Run2016GH":{
    #     "mc_jpsiphi": mc_jpsiphi_2016GH,
    #     "mc_jpsiphi_tree": "bspsiphiMc",
    #     "mc_jpsikstar": mc_jpsikstar_2016GH,
    #     "mc_jpsikstar_tree": "bdpsikstarMc",
    #     "data_tree": "bspsiphiData",
    #     "trigger": "HLT_DoubleMu4_3_Jpsi_Displaced",
    #     "fix_resolution": True,
    #     "data": [
    #         data_path + "Charmonium+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    #     ]
    # },
}

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    """Print canvas in different formats"""
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

def make_dataset(input_tree, name, cuts="", weight=None):
    """Build RooDataSet from a TTree or TChain"""
    
    # pre-filter tree to memory
    # ROOT.gROOT.cd()
    ftmp = ROOT.TFile.Open("/tmp/dmytro/tmp.root","RECREATE")
    tree = input_tree.CopyTree(cuts)
    print "Nubmer of events after filter:", tree.GetEntries()
    
    var_set = ROOT.RooArgSet()
    var_set.add(mass)
    var_set.add(chan)

    if weight == None:
        data = ROOT.RooDataSet(name, "", var_set, ROOT.RooFit.Import(tree))
    else:
        var_set.add(weight)
        data = ROOT.RooDataSet(name, "", var_set, ROOT.RooFit.Import(tree), ROOT.RooFit.WeightVar(weight))
        
    # if cuts != "":
    #    data = data.reduce(cuts)

    # data.Print("V")
    print "Input tree has ", tree.GetEntries(), "entries. The derived dataset has ", data.sumEntries()

    return data

def build_model(mass_var, ref_jpsiphi_hist, ref_jpsikstar_hist, jpsikstar_fraction):
    """Build fit model and save it in a workspace"""

    # JpsiPhi signal pdf
    bias     = ROOT.RooRealVar("bias", "bias", 0, -0.01, 0.01)
    sigma    = ROOT.RooRealVar("sigma", "sigma", 0.0001, 0., 0.01)
    gaussM   = ROOT.RooGaussModel("gaussM", "signal pdf", mass_var, bias, sigma)
    ref_data = ROOT.RooDataHist("ref_data", "", ROOT.RooArgList(mass_var), ref_jpsiphi_hist)
    
    ref_pdf  = ROOT.RooHistPdf("ref_pdf", "theoretical lineshape", ROOT.RooArgSet(mass_var), ref_data, 2)
    mass_var.setBins(10000, "fft");
    sig      = ROOT.RooFFTConvPdf("sig", "smeared distribution", mass_var, ref_pdf, gaussM)

    # JpsiKstar background pdf
    jpsikstar_data = ROOT.RooDataHist("jpsikstar_data", "", ROOT.RooArgList(mass_var), ref_jpsikstar_hist)
    jpsikstar_pdf  = ROOT.RooHistPdf("jpsikstar_pdf", "theoretical lineshape", ROOT.RooArgSet(mass_var), jpsikstar_data, 2)
    
    a0    = ROOT.RooRealVar("a0", "a0", 0.0, -1.0, 1.0)
    a1    = ROOT.RooRealVar("a1", "a1", 0.0, -0.5, 0.5)
    bkg   = ROOT.RooChebychev("bkg", "Background", mass_var, ROOT.RooArgList(a0, a1))

    # b0 = ROOT.RooRealVar("b0","b0", 0.5, 0, 1.)
    # # b1 = ROOT.RooRealVar("b1","b1", 0.5, 0, 1.)
    # bkg = ROOT.RooBernstein("bkg","Background", mass_var, ROOT.RooArgList(b0))
    
    # exp_c = ROOT.RooRealVar("exp_c","exp_c", -1, -1000., 0.)
    # bkg   = ROOT.RooExponential("bkg", "Background", mass_var, exp_c)
    
    Nsig  = ROOT.RooRealVar("Nsig", "Nsig", 1000, 0, 1e9)
    Nbkg  = ROOT.RooRealVar("Nbkg", "Nbkg", 0, 0, 1e9)
    Njpsikstar = ROOT.RooFormulaVar("Njpsikstar", "Njpsikstar", "@0*%s" % jpsikstar_fraction, ROOT.RooArgList(Nsig))
    
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg,jpsikstar_pdf), ROOT.RooArgList(Nsig,Nbkg,Njpsikstar))

    ws = ROOT.RooWorkspace("ws","")
    getattr(ws,'import')(model)

    return ws

def get_generated_number_of_events(file_pattern):
    n = 0
    chain = ROOT.TChain("info")
    chain.Add(file_pattern)
    for entry in chain:
        n += entry.n_gen_all
    return n

############################################################################

ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 800, 800)

for dataset, info in datasets.items():
    name = dataset
    if name in results and not recompute_results:
        print "Results are already available for %s. Skip" % name
        continue
    
    results[name] = dict()
    selection = "(run == 1 || certified_muon) && abs(kk_mass - 1.02)<0.01 && m>%s && m<%s" % (min_mass, max_mass)
    
    print "Processing", name

    ## Get MC datasets

    # B to Jpsi Phi
    chain_mc_jpsiphi = ROOT.TChain(info['mc_jpsiphi_tree'])
    chain_mc_jpsiphi.Add(info['mc_jpsiphi'])
    print "Number of Jpsiphi MC events:", chain_mc_jpsiphi.GetEntries()
    if chain_mc_jpsiphi.GetEntries() == 0:
        raise Exception("No Jpsiphi MC events found for for " + name)
    ds_mc_jpsiphi = make_dataset(chain_mc_jpsiphi, "jpsiphi", selection + "&& mc_match!=0")

    n_gen_jpsiphi = get_generated_number_of_events(info['mc_jpsiphi'])
    eff_jpsiphi   = ds_mc_jpsiphi.sumEntries() / n_gen_jpsiphi
    print "Efficiency Jpsiphistar: ", eff_jpsiphi

    # B to Jpsi Pi
    chain_mc_jpsikstar = ROOT.TChain(info['mc_jpsikstar_tree'])
    chain_mc_jpsikstar.Add(info['mc_jpsikstar'])
    print "Number of Jpsikstar MC events:", chain_mc_jpsikstar.GetEntries()
    if chain_mc_jpsikstar.GetEntries() == 0:
        raise Exception("No Jpsikstar MC events found for for " + name)
    ds_mc_jpsikstar = make_dataset(chain_mc_jpsikstar, "jpsikstar", selection + "&& mc_match!=0")

    n_gen_jpsikstar = get_generated_number_of_events(info['mc_jpsikstar'])
    eff_jpsikstar =  ds_mc_jpsikstar.sumEntries() / n_gen_jpsikstar
    print "Efficiency Jpsikstar: ", eff_jpsikstar

    jpsikstar_fraction = jpsikstar_bf / jpsiphi_bf * eff_jpsikstar / eff_jpsiphi * fsfd_mc / fsfd
    print "Jpsikstar yield is fixed to JpsiPhi at", jpsikstar_fraction
    
    ## Get Data datasets
    
    chain_data = ROOT.TChain(info['data_tree'])
    for pattern in info['data']:
        chain_data.Add(pattern)
    print "Number of Data events:", chain_data.GetEntries()
    if chain_data.GetEntries() == 0:
        raise Exception("No Data events found for for " + name)
    weight = RooRealVar("%s_ps" % info['trigger'] ,"", 0)
    ds_data = make_dataset(chain_data, "data", selection, weight)
    
    for ich in range(2):
        if ich == 0:
            ds_mc_jpsiphi_ch   = ds_mc_jpsiphi.reduce("chan<0.5")
            ds_mc_jpsikstar_ch  = ds_mc_jpsikstar.reduce("chan<0.5")
            ds_data_ch = ds_data.reduce("chan<0.5")
        else:
            ds_mc_jpsiphi_ch   = ds_mc_jpsiphi.reduce("chan>0.5")
            ds_mc_jpsikstar_ch  = ds_mc_jpsikstar.reduce("chan>0.5")
            ds_data_ch = ds_data.reduce("chan>0.5")

        print "Number of Jpsiphi MC events in chan%u: %d"   % (ich, ds_mc_jpsiphi_ch.sumEntries())
        print "Number of Jpsikstar MC events in chan%u: %d"  % (ich, ds_mc_jpsikstar_ch.sumEntries())
        print "Number of Data events in chan%u: %d" % (ich, ds_data_ch.sumEntries())

        ## Get reference signal mass distribution
        
        jpsiphi_hist = ROOT.RooAbsData.createHistogram(ds_mc_jpsiphi_ch, 'jpsiphi', mass, mass_ref_binning)
        jpsiphi_hist.SetDirectory(0)
        jpsiphi_hist.Draw()
        print_canvas(name + "_mc_jpsiphi_chan%u_mass_ref" % ich, output_path)
        
        jpsikstar_hist = ROOT.RooAbsData.createHistogram(ds_mc_jpsikstar_ch, 'jpsikstar', mass, mass_ref_binning)
        jpsikstar_hist.SetDirectory(0)
        jpsikstar_hist.Draw()
        print_canvas(name + "_mc_jpsikstar_chan%u_mass_ref" % ich, output_path)

        ws = build_model(mass, jpsiphi_hist, jpsikstar_hist, jpsikstar_fraction)

        ## Fit Data
    
        model = ws.pdf("model")
        model.fitTo(ds_data_ch,  ROOT.RooFit.NumCPU(8), ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE), ROOT.RooFit.PrintLevel(print_level))

        if info['fix_resolution']:
            ws.var("sigma").setConstant(True)
            model.fitTo(ds_data_ch,  ROOT.RooFit.NumCPU(8), ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE), ROOT.RooFit.PrintLevel(print_level))

        ## Plot results
    
        frame = mass.frame()
        ds_data_ch.plotOn(frame, mass_plot_binning)
        frame.SetMaximum(frame.GetMaximum() * 1.2)
        bkgs = ROOT.RooArgSet()
        bkgs.add(ws.pdf("bkg"))
        bkgs.add(ws.pdf("jpsikstar_pdf"))
        model.plotOn(frame, ROOT.RooFit.Components(bkgs),  ROOT.RooFit.LineColor(ROOT.kRed))
        model.plotOn(frame, ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        model.plotOn(frame)
        if info['fix_resolution']:
            print "chiSquare/nDOF: ", frame.chiSquare(5)
        else:
            print "chiSquare/nDOF: ", frame.chiSquare(6)
            
        # print "chiSquare: ", frame.chiSquare("model","data", 6)

        model.paramOn(frame, ROOT.RooFit.Layout(0.7, 0.95, 0.92))
        frame.getAttText().SetTextSize(0.02)
        frame.Draw()
        print_canvas(name + "_chan%u_mass_fit" % ich, output_path)

        ## Store results
        # {"2016BF": {"chan0": {"val": 305019.3, "err": 1167.6},
        results[name]["chan%u" % ich] = {"val": ws.var("Nsig").getVal(), "err": ws.var("Nsig").getError()}
        
        ws.Delete() # Cleanup

    ## Save results
    json.dump(results, open(output_path + "/results.json", "w"))
    
print results


