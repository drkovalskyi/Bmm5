import ROOT, subprocess

min_jpsi_mass = 2.85
max_jpsi_mass = 3.35
nbins_jpsi = 100
mass  = ROOT.RooRealVar("mass", "", (max_jpsi_mass + min_jpsi_mass) / 2, min_jpsi_mass, max_jpsi_mass)
mass_binning  = ROOT.RooFit.Binning(nbins_jpsi, min_jpsi_mass, max_jpsi_mass)

ROOT.ROOT.EnableImplicitMT()

ROOT.gInterpreter.Declare(
'''
using namespace ROOT::VecOps;
RVec<double> pair_mass(const RVec<double> &pt, 
                  const RVec<double> &eta, 
                  const RVec<double> &phi,
                  const RVec<int>    &charge,
                  const double       &mass) {
  RVec<double> res;
  res.reserve(pt.size()*pt.size());

  if (pt.size()>1) {
    for (std::size_t i=0; i<pt.size()-1; ++i) {
       if (pt[i] < 4) continue;
       ROOT::Math::PtEtaPhiMVector m1(pt[i], eta[i], phi[i], mass);
       for (std::size_t j=i+1; j<pt.size(); ++j) {
         if (pt[j] < 4) continue;
         if (charge[i] * charge[j] > 0) continue;
         ROOT::Math::PtEtaPhiMVector m2(pt[j], eta[j], phi[j], mass);
         res.emplace_back((m1+m2).M());
       }
    }
  }
  return res;
}

// https://github.com/cms-sw/cmssw/blob/7cbe615827ebfc80e770cd875b8e68c2e2b7a084/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L114
RVec<double> pair_mass_trig(const RVec<double> &pt, 
                  const RVec<double> &eta, 
                  const RVec<double> &phi,
                  const RVec<int>    &ids,
                  const double       &mass) {
  RVec<double> res;
  res.reserve(pt.size()*pt.size());

  if (pt.size()>1) {
    for (std::size_t i=0; i<pt.size()-1; ++i) {
       if (ids[i] != 13) continue;
       if (pt[i] < 5) continue;
       ROOT::Math::PtEtaPhiMVector m1(pt[i], eta[i], phi[i], mass);
       for (std::size_t j=i+1; j<pt.size(); ++j) {
         if (ids[j] != 13) continue;
         if (pt[j] < 5) continue;
         ROOT::Math::PtEtaPhiMVector m2(pt[j], eta[j], phi[j], mass);
         res.emplace_back((m1+m2).M());
       }
    }
  }
  return res;
}
'''
)

def build_model(workspace_name, mass_var, peak=3.09, search_width=0.04):
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
    sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_fract))
    # sig  = ROOT.RooAddPdf("sig"," ",ROOT.RooArgList(G3,G2,G1),ROOT.RooArgList(G2_fract,G3_fract))

    # # CB
    # sig_mean  = ROOT.RooRealVar("sig_mean",  "", peak, peak - search_width, peak + search_width)
    # sig_sigma = ROOT.RooRealVar("sig_sigma", "sigma", 0.03, 0.001, 0.10)
    # sig_tail  = ROOT.RooRealVar("sig_tail",  "tail", 1, 0.1, 2.0)
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
    # sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,sig_cb), ROOT.RooArgList(G2_fract))
    # sig  = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G3,G2,sig_cb), ROOT.RooArgList(G2_fract,G3_fract))
    
    ## Combinatorial background
    
    # a0    = ROOT.RooRealVar("a0", "a0", -0.8, -1,  1.0)
    # a1    = ROOT.RooRealVar("a1", "a1", 0.0, -0.3, 0.3)
    # bkg   = ROOT.RooChebychev("bkg", "Background", mass_var, ROOT.RooArgList(a0, a1))

    # b0 = ROOT.RooRealVar("b0","b0", 0.4, 1e-5, 1.)
    # # b1 = ROOT.RooRealVar("b1","b1", 0.5, 0, 1.)
    # bkg = ROOT.RooBernstein("bkg","Background", mass_var, ROOT.RooArgList(b0))
    
    exp_c = ROOT.RooRealVar("exp_c","exp_c", -1, -1000., 0.)
    bkg   = ROOT.RooExponential("bkg", "Background", mass_var, exp_c)
    
    Nsig  = ROOT.RooRealVar("Nsig", "Nsig", 1000, 0, 1e9)
    Nbkg  = ROOT.RooRealVar("Nbkg", "Nbkg", 0, 0, 1e9)
    # Njpsipi = ROOT.RooFormulaVar("Njpsipi", "Njpsipi", "@0*%s" % jpsipi_fraction, ROOT.RooArgList(Nsig))
    
    # model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg,jpsipi_pdf), ROOT.RooArgList(Nsig,Nbkg,Njpsipi))
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg), ROOT.RooArgList(Nsig,Nbkg))

    ws = ROOT.RooWorkspace(workspace_name, "")
    getattr(ws,'import')(model)
    return ws

tree = ROOT.TChain("Events")


directories = [
    # "/eos/cms/store/data/Run2025C/ZeroBias/NANOAOD/PromptReco-v1/000/393/"
    # "/eos/cms/store/data/Run2025C/ZeroBias/NANOAOD/PromptReco-v1/",
    # "/eos/cms/store/data/Run2025C/ZeroBias/NANOAOD/PromptReco-v2/"
    # "/eos/cms/store/data/Run2025C/HLTPhysics/NANOAOD/PromptReco-v1/",
    # "/eos/cms/store/data/Run2025C/HLTPhysics/NANOAOD/PromptReco-v2/",
    # "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass0/NANOAOD/PromptReco-v1/000/393/"
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass0/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass0/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass1/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass1/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass2/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass2/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass3/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass3/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass4/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass4/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass5/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass5/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass6/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass6/NANOAOD/PromptReco-v2/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass7/NANOAOD/PromptReco-v1/",
    "/eos/cms/store/data/Run2025C/ParkingDoubleMuonLowMass7/NANOAOD/PromptReco-v2/",
]

nfiles = 0
for dir in directories:
    command = f'find {dir} -type f -path "*NANOAOD*.root"'
    result = subprocess.check_output(command, shell=True, text=True)
    for file in result.splitlines():
        nfiles += tree.Add(file)

print("Number of files:", nfiles)

rdf = ROOT.RDataFrame(tree)
rdf = rdf.DefaultValueFor("HLT_DoubleMu4_3_LowMass", False)
rdf = rdf.DefaultValueFor("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6", False)
rdf = rdf.DefaultValueFor("HLT_Mu0_L1DoubleMu", False)
# rdf = rdf.Filter("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6")
# rdf = rdf.Filter("HLT_Mu0_L1DoubleMu")
rdf = rdf.Filter("HLT_Mu4_L1DoubleMu")
# rdf2 = rdf.Filter("!HLT_DoubleMu4_3_LowMass")
rdf2 = rdf.Filter("!HLT_DoubleMu4_3_LowMass")

# rdf = rdf.Define("candMass", "pair_mass(Muon_pt, Muon_eta, Muon_phi, Muon_charge, 0.1057)")
rdf = rdf.Define("candMass", "pair_mass_trig(TrigObj_pt, TrigObj_eta, TrigObj_phi, TrigObj_id, 0.1057)")
rdf2 = rdf2.Define("candMass2", "pair_mass_trig(TrigObj_pt, TrigObj_eta, TrigObj_phi, TrigObj_id, 0.1057)")

h = rdf.Histo1D(("h", "Mass;Mass (GeV/c^{2});Events", nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "candMass")
h2 = rdf2.Histo1D(("h", "Mass;Mass (GeV/c^{2});Events", nbins_jpsi, min_jpsi_mass, max_jpsi_mass), "candMass2")

c = ROOT.TCanvas("c", "", 1600, 800)
c.Divide(2,1)
c.cd(1)
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

h.SetFillColor(ROOT.kMagenta)
h.Draw()

# c.SetLogy()
# c.Draw()

ws = build_model("jpsi", mass)

model = ws.pdf("model")
data = ROOT.RooDataHist("data", "", ROOT.RooArgList(mass), h.GetPtr())
data2 = ROOT.RooDataHist("data2", "", ROOT.RooArgList(mass), h2.GetPtr())

# preselection fit
ws.var("Nsig").setVal(h.GetEntries() * 0.9)
ws.var("Nbkg").setVal(h.GetEntries() * 0.1)
# ws.var("sig_G2_fract").setVal(0.0)
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
model.fitTo(data,  ROOT.RooFit.NumCPU(8),
            ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
            ROOT.RooFit.PrintLevel(0))

## Plot results
    
frame = mass.frame()
data.plotOn(frame)
# data.plotOn(frame, mass_binning)
frame.SetMaximum(frame.GetMaximum() * 1.2)
# model.plotOn(frame, ROOT.RooFit.Components("sig"), ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(frame, ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))
model.plotOn(frame)
print("chiSquare: ", frame.chiSquare(6))
print("chiSquare: ", frame.chiSquare("model","data", 6))

model.paramOn(frame, ROOT.RooFit.Layout(0.7, 0.95, 0.92))
frame.getAttText().SetTextSize(0.02)
frame.Draw()

c.cd(2)
# failed events fit
# ws.var("sig_G1_mean").setConstant(True)
# ws.var("sig_G1_sigma").setConstant(True)
# ws.var("sig_G2_scale").setConstant(True)
# ws.var("sig_G2_fract").setConstant(True)
model.fitTo(data2,  ROOT.RooFit.NumCPU(8),
            ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE),
            ROOT.RooFit.PrintLevel(0))

frame2 = mass.frame()
data2.plotOn(frame2)
frame2.SetMaximum(frame2.GetMaximum() * 1.2)
model.plotOn(frame2, ROOT.RooFit.Components("bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed))
model.plotOn(frame2)
print("chiSquare: ", frame2.chiSquare(6))
print("chiSquare: ", frame2.chiSquare("model","data2", 6))

model.paramOn(frame2, ROOT.RooFit.Layout(0.7, 0.95, 0.92))
frame2.getAttText().SetTextSize(0.02)
frame2.Draw()
