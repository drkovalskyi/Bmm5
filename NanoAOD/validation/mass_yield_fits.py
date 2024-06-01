import os, re, copy, json
import ROOT
import math 
import tdrstyle
from ROOT.RooFit import Binning
from ROOT import RooRealVar
from fwhm_calculator import compute_fwhm

"""
Bmm5 BstoJpsiPhi fit
Unbinned maximum likelihood fit using bkkmm flat ntuples. 
"""

# Set the TDR style
tdrstyle.setTDRStyle()

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/scouting/bsjpsiphi/";

recompute_results = True
result_file = output_path + "/results.json"
results = dict()
if os.path.exists(result_file):
    results = json.load(open(result_file))
    
# Silence RooFit info messages
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
print_level = 0

variables = [
    ( RooRealVar("HLT_DoubleMu4_3_LowMass","", 0),   None, None), # UInt_t is not supported as input for RooCategory
    ( RooRealVar("DST_PFScouting_DoubleMuon","", 0), None, None), # UInt_t is not supported as input for RooCategory
]

roovars = []
for var in variables:
    roovars.append(var[0])

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/"
    
datasets = {

    # "Run2024D-Scouting":{
    #     "lumi": "8.156",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     # "trigger": "HLT_DoubleMu4_3_Jpsi",
    #     # "fix_resolution": False,
    #     "data": [
    #         # data_path + "ScoutingPFRun3+Run2024C-v1+HLTSCOUT/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ScoutingPFRun3+Run2024D-v1+HLTSCOUT/*.root",
    #     ],
    # },

    # "Run2024D-Scouting-DST":{
    #     "lumi": "8.156",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     # "trigger": "HLT_DoubleMu4_3_Jpsi",
    #     # "fix_resolution": False,
    #     "data": [
    #         # data_path + "ScoutingPFRun3+Run2024C-v1+HLTSCOUT/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ScoutingPFRun3+Run2024D-v1+HLTSCOUT/*.root",
    #     ],
    # },

    # "Run2024E-Scouting-Vtx-ZB-Jpsi":{
    #     "lumi": "6.946",
    #     "unit": "pb",
    #     "data_tree": "jpsiData",
    #     "min_mass": 2.8,
    #     "max_mass": 3.4,
    #     "nbins": 80,
    #     "peak": 3.1,
    #     "search_window": 0.1,
    #     # "trigger": "HLT_DoubleMu4_3_Jpsi",
    #     # "fix_resolution": False,
    #     "data": [
    #         data_path + "crab-140x-mm-Vtx-ZB/jpsi/ScoutingPFRun3+Run2024E-v1+HLTSCOUT/*.root",
    #     ],
    # },
    
    # "Run2024E-Scouting-Vtx-ZB-Jpsi_HLT_DoubleMu4_3_LowMass":{
    #     "lumi": "6.946",
    #     "unit": "pb",
    #     "data_tree": "jpsiData",
    #     "min_mass": 2.8,
    #     "max_mass": 3.4,
    #     "nbins": 80,
    #     "peak": 3.1,
    #     "search_window": 0.1,
    #     "trigger": "HLT_DoubleMu4_3_LowMass",
    #     # "fix_resolution": False,
    #     "data": [
    #         data_path + "crab-140x-mm-Vtx-ZB/jpsi/ScoutingPFRun3+Run2024E-v1+HLTSCOUT/*.root",
    #     ],
    # },

    # "Run2024E-Scouting-Vtx-ZB-Jpsi_DST_PFScouting_DoubleMuon":{
    #     "lumi": "6.946",
    #     "unit": "pb",
    #     "data_tree": "jpsiData",
    #     "min_mass": 2.8,
    #     "max_mass": 3.4,
    #     "nbins": 80,
    #     "peak": 3.1,
    #     "search_window": 0.1,
    #     "trigger": "DST_PFScouting_DoubleMuon",
    #     # "fix_resolution": False,
    #     "data": [
    #         data_path + "crab-140x-mm-Vtx-ZB/jpsi/ScoutingPFRun3+Run2024E-v1+HLTSCOUT/*.root",
    #     ],
    # },
    
    # "Run2024D-Scouting-Vtx":{
    #     "lumi": "8.574",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     # "trigger": "HLT_DoubleMu4_3_Jpsi",
    #     # "fix_resolution": False,
    #     "data": [
    #         data_path + "crab-140x-mm-Vtx/bkkmm/ScoutingPFRun3+Run2024D-v1+HLTSCOUT/*.root",
    #     ],
    # },
    
    # "Run2024D-Scouting-NoVtx-DST":{
    #     "lumi": "8.327",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     "min_mass": 5.20,
    #     "max_mass": 5.56,
    #     "nbins": 72,
    #     "peak": 5.36,
    #     "search_window": 0.04,
    #     "trigger": "DST_PFScouting_DoubleMuon",
    #     # "fix_resolution": False,
    #     "data": [
    #         data_path + "crab-140x-mm/bkkmm/ScoutingPFRun3+Run2024D-v1+HLTSCOUT/*.root",
    #     ],
    # },
    
    # "Run2024-Scouting":{
    #     "lumi": "16.0",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     "data": [
    #         data_path + "crab-140x/bkkmm/ScoutingPFRun3+Run2024C-v1+HLTSCOUT/*.root",
    #         data_path + "crab-140x/bkkmm/ScoutingPFRun3+Run2024D-v1+HLTSCOUT/*.root",
    #     ],
    # },

    # "Run2024-ScoutingMon":{
    #     "lumi": "13.8",
    #     "unit": "pb",
    #     "data_tree": "bspsiphiData",
    #     "data": [
    #         data_path + "crab-140x-mm/bkkmm/ScoutingPFMonitor+Run2024C-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-mm/bkkmm/ScoutingPFMonitor+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #     ],
    # },
    
    # "Run2024D-Parking":{
    #     "lumi": "7.8",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     "data": [
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass0+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass1+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass2+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass3+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass4+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass5+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass6+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm/ParkingDoubleMuonLowMass7+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #     ],
    # },
    # "Run2024D-Parking-DST":{
    #     "lumi": "7.8",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     "data": [
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass0+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass1+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass2+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass3+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass4+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass5+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass6+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #         data_path + "crab-140x-bkkmm/bkkmm-dst/ParkingDoubleMuonLowMass7+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #     ],
    # },
    
    # "Run2024D-Parking-DST":{
    #     "lumi": "0.992",
    #     "unit": "fb",
    #     "data_tree": "bspsiphiData",
    #     "data": [
    #         data_path + "crab-140x-mm/bkkmm-dst/ParkingDoubleMuonLowMass0+Run2024D-PromptReco-v1+MINIAOD/*.root",
    #     ],
    # },
    
    "Run2024D-Parking-HLT":{
        "lumi": "0.992",
        "unit": "fb",
        "data_tree": "bspsiphiData",
        "min_mass": 5.20,
        "max_mass": 5.56,
        "nbins": 72,
        "peak": 5.36,
        "search_window": 0.04,
        "trigger": "HLT_DoubleMu4_3_LowMass",
        "data": [
            data_path + "crab-140x-mm/bkkmm/ParkingDoubleMuonLowMass0+Run2024D-PromptReco-v1+MINIAOD/*.root",
        ],
    },
}


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    """Print canvas in different formats"""
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

def make_dataset(input_tree, name, mass_var, other_vars, cuts="", weight=""):
    """Build RooDataSet from a TTree or TChain"""
    
    # pre-filter tree to memory
    # ROOT.gROOT.cd()
    ftmp = ROOT.TFile.Open("/tmp/dmytro/tmp.root","RECREATE")
    tree = input_tree.CopyTree(cuts)
    print("Nubmer of events after filter:", tree.GetEntries())
    
    var_set = ROOT.RooArgSet()
    # for var in other_vars:
    #     var_set.add(var)
    var_set.add(mass_var)

    if weight == "":
        data = ROOT.RooDataSet(name, "", var_set, ROOT.RooFit.Import(tree))
    else:
        data = ROOT.RooDataSet(name, "", var_set, ROOT.RooFit.Import(tree), ROOT.RooFit.WeightVar(weight))
        
    # if cuts != "":
    #    data = data.reduce(cuts)

    # data.Print("V")
    print("Input tree has ", tree.GetEntries(), "entries. The derived dataset has ", data.sumEntries())

    return data

def build_model(mass_var, peak=5.36, search_width=0.04):
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
    
    # sig_mu     = ROOT.RooRealVar("sig_mu", "mu", 5.36, 5.3, 5.4)
    # sig_lambda = ROOT.RooRealVar("sig_lambda", "lambda", 0.003, 0.001, 0.1)
    # sig_gamma  = ROOT.RooRealVar("sig_gamma", "gamma", 1, 0, 5)
    # sig_delta  = ROOT.RooRealVar("sig_delta", "delta", 1, 0.1, 10)
    # sig = ROOT.RooJohnson("sig", "signal", m, sig_mu, sig_lambda, sig_gamma, sig_delta)

    # multi-gaussian
    # peak = 5.36
    G1_mean  = ROOT.RooRealVar("sig_G1_mean",  "", peak, peak - search_width, peak + search_width)
    G1_sigma = ROOT.RooRealVar("sig_G1_sigma", "", 0.03, 0.001, 0.10)
    G2_scale = ROOT.RooRealVar("sig_G2_scale", "", 2.5, 0.2, 7.5)
    G3_scale = ROOT.RooRealVar("sig_G3_scale", "", 3.0, 0.5, 6.7)
    G2_sigma = ROOT.RooProduct("sig_G2_sigma", "", ROOT.RooArgList(G1_sigma,G2_scale))
    G3_sigma = ROOT.RooProduct("sig_G3_sigma", "", ROOT.RooArgList(G1_sigma,G3_scale))
    G1 = ROOT.RooGaussian("sig_G1", "", mass_var, G1_mean, G1_sigma)
    G2 = ROOT.RooGaussian("sig_G2", "", mass_var, G1_mean, G2_sigma)
    G3 = ROOT.RooGaussian("sig_G3", "", mass_var, G1_mean, G3_sigma)
    
    G2_fract = ROOT.RooRealVar("sig_G2_fract","",0.3,0.0,1.0)
    G3_fract = ROOT.RooRealVar("sig_G3_fract","",0.2,0.0,1.0)
    # sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_fract))
    sig  = ROOT.RooAddPdf("sig"," ",ROOT.RooArgList(G3,G2,G1),ROOT.RooArgList(G2_fract,G3_fract))

    # # CB
    # peak = 5.36
    # sig_mean  = ROOT.RooRealVar("sig_mean",  "", peak, peak-0.04, peak+0.04)
    # sig_sigma = ROOT.RooRealVar("sig_sigma", "sigma", 0.03, 0.001, 0.10)
    # sig_tail  = ROOT.RooRealVar("sig_tail",  "tail", 2.8, 0.1, 10.0)
    # sig_pow   = ROOT.RooRealVar("sig_pow0",  "pow", 3, 0, 50)
    # sig = ROOT.RooCBShape("sig", "signal", mass_var, sig_mean, sig_sigma, sig_tail, sig_pow)
    
    ## Combinatorial background
    
    a0    = ROOT.RooRealVar("a0", "a0", -0.8, -1,  1.0)
    a1    = ROOT.RooRealVar("a1", "a1", 0.0, -0.3, 0.3)
    bkg   = ROOT.RooChebychev("bkg", "Background", mass_var, ROOT.RooArgList(a0, a1))

    # b0 = ROOT.RooRealVar("b0","b0", 0.4, 1e-5, 1.)
    # # b1 = ROOT.RooRealVar("b1","b1", 0.5, 0, 1.)
    # bkg = ROOT.RooBernstein("bkg","Background", mass_var, ROOT.RooArgList(b0))
    
    # exp_c = ROOT.RooRealVar("exp_c","exp_c", -1, -1000., 0.)
    # bkg   = ROOT.RooExponential("bkg", "Background", mass_var, exp_c)
    
    Nsig  = ROOT.RooRealVar("Nsig", "Nsig", 1000, 0, 1e9)
    Nbkg  = ROOT.RooRealVar("Nbkg", "Nbkg", 0, 0, 1e9)
    # Njpsipi = ROOT.RooFormulaVar("Njpsipi", "Njpsipi", "@0*%s" % jpsipi_fraction, ROOT.RooArgList(Nsig))
    
    # model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg,jpsipi_pdf), ROOT.RooArgList(Nsig,Nbkg,Njpsipi))
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg), ROOT.RooArgList(Nsig,Nbkg))

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
# c1 = ROOT.TCanvas("c1","c1", 800, 800)
import cmsstyle as CMS
CMS.SetExtraText("Preliminary")
CMS.SetEnergy("13.6")

for dataset, info in list(datasets.items()):
    name = dataset
    if name in results and not recompute_results:
        print("Results are already available for %s. Skip" % name)
        continue
    
    # min_mass = 5.20
    # max_mass = 5.56
    min_mass = info["min_mass"]
    max_mass = info["max_mass"]
    selection = "m>%s && m<%s" % (min_mass, max_mass)
    if "trigger" in info:
        selection += " && %s>0" % (info['trigger'])
    
    print("Processing", name)

    ## Get Data datasets
    
    chain_data = ROOT.TChain(info['data_tree'])
    for pattern in info['data']:
        print(pattern)
        chain_data.Add(pattern)

    print("Number of Data events:", chain_data.GetEntries())
    if chain_data.GetEntries() == 0:
        raise Exception("No Data events found for for " + name)

    nbins = info["nbins"]

    mass = ROOT.RooRealVar("m", "",  min_mass, max_mass)
    mass_plot_binning = Binning(nbins, min_mass, max_mass)

    ds_data = make_dataset(chain_data, "data", mass, roovars, selection)
    
    ws = build_model(mass, info["peak"], info["search_window"])

    ## Fit Data
    
    model = ws.pdf("model")

    # prefit
    ws.var("Nsig").setVal(chain_data.GetEntries() * 0.9)
    ws.var("Nbkg").setVal(chain_data.GetEntries() * 0.1)
    ws.var("sig_G2_fract").setVal(0.0)
    ws.var("sig_G2_fract").setConstant(True)
    ws.var("sig_G2_scale").setConstant(True)
    ws.var("sig_G3_fract").setVal(0.0)
    ws.var("sig_G3_fract").setConstant(True)
    ws.var("sig_G3_scale").setConstant(True)
    model.fitTo(ds_data,  ROOT.RooFit.NumCPU(8),
                ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
                ROOT.RooFit.PrintLevel(print_level))

    # final fit
    ws.var("sig_G2_fract").setConstant(False)
    ws.var("sig_G2_scale").setConstant(False)
    ws.var("sig_G3_fract").setConstant(False)
    ws.var("sig_G3_scale").setConstant(False)
    model.fitTo(ds_data,  ROOT.RooFit.NumCPU(8),
                ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
                ROOT.RooFit.PrintLevel(print_level))
    
    fwhm = compute_fwhm(ws.pdf("sig"), ws.var("m"), n_steps = 10000)
    print(f"Full Width at Half Maximum (FWHM): {fwhm}")
    
    ## Plot results
    
    frame = mass.frame()
    ds_data.plotOn(frame, mass_plot_binning)
    # frame.SetMaximum(frame.GetMaximum() * 1.2)
    bkgs = ROOT.RooArgSet()
    bkgs.add(ws.pdf("bkg"))
    model.plotOn(frame, ROOT.RooFit.Components(bkgs),  ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame, ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    model.plotOn(frame)
    print("chiSquare: ", frame.chiSquare(6))

    CMS.SetLumi(info["lumi"], info["unit"])
    c1 = CMS.cmsCanvas('c1', min_mass, max_mass, 0, frame.GetMaximum(),
                       'm (GeV)', '', square = CMS.kSquare, extraSpace=0.01, iPos=0)
    
    # model.paramOn(frame, ROOT.RooFit.Layout(0.7, 0.95, 0.92), ROOT.RooFit.Parameters(ws.var("Nsig"),ws.var("Nbkg")))
    # model.paramOn(frame, 
    #               ROOT.RooFit.Layout(0.65, 0.9, 0.9),
    #               ROOT.RooFit.Parameters(ROOT.RooArgSet(ws.var("Nsig"),
    #                                                     ws.var("Nbkg")))
    # )
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.DrawLatexNDC(0.65, 0.80, "N_{sig} = %0.0f \pm %0.0f" % (ws.var("Nsig").getVal(), ws.var("Nsig").getError()))
    latex.DrawLatexNDC(0.65, 0.75, "N_{bkg} = %0.0f \pm %0.0f" % (ws.var("Nbkg").getVal(), ws.var("Nbkg").getError()))
    latex.DrawLatexNDC(0.65, 0.70, "FWHM = %0.1f MeV" % (fwhm * 1e3) )    

    # model.paramOn(frame, 
    #               ROOT.RooArgSet(ws.var("Nsig"),
    #                              ws.var("Nbkg"),
    #                              rf_fwhm),
    #               ROOT.RooFit.Layout(0.65, 0.9, 0.9),
    #               )
    # frame.getAttText().SetTextSize(0.03)
    frame.Draw("same")
    
    print_canvas(name + "_mass_fit", output_path)
    # CMS.SaveCanvas(c1, os.path.join(output_path, name + "_mass_fit_cms.png"))

    ## Store results
    # {"2016BF": {"val": 305019.3, "err": 1167.6},
    results[name] = {"val": ws.var("Nsig").getVal(), "err": ws.var("Nsig").getError()}

    # ws.Delete() # Cleanup

## Save results
json.dump(results, open(output_path + "/results.json", "w"))
    
print(results)


