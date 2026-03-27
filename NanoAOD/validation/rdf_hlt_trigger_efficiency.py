import ROOT
from ROOT import RDataFrame

ROOT.ROOT.EnableImplicitMT()

def simple_dump(rdf, name):
    data = rdf.AsNumpy([name])
    for i, rvec in enumerate(data[name]):
        print(f"Entry {i}: ", list(rvec))

def dump(rdf, column_names):
    data = rdf.AsNumpy(column_names)
    n_entries = len(next(iter(data.values())))
    for i in range(n_entries):
        entry_values = {col: list(data[col][i]) for col in column_names}
        entry_str = ", ".join(f"{col} = {entry_values[col]}" for col in column_names)
        print(f"Entry {i}: {entry_str}")

def make_trigger_efficiency_plots(name, h_trig, h_all):
    h_eff = h_trig.Clone(name)
    h_eff.Divide(h_trig.GetPtr(), h_all.GetPtr(), 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.Draw("e0")
    return h_eff

chain = ROOT.TChain("Events")
chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/11*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass1+Run2022C-PromptReco-v1+MINIAOD/*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass2+Run2022C-PromptReco-v1+MINIAOD/*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass3+Run2022C-PromptReco-v1+MINIAOD/*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass4+Run2022C-PromptReco-v1+MINIAOD/*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass5+Run2022C-PromptReco-v1+MINIAOD/*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass6+Run2022C-PromptReco-v1+MINIAOD/*.root")
# chain.Add("/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ParkingDoubleMuonLowMass7+Run2022C-PromptReco-v1+MINIAOD/*.root")

rdf = RDataFrame(chain)

rdf = rdf.Filter("HLT_Mu4_L1DoubleMu")
# print(f"Number of events triggered by HLT_Mu4_L1DoubleMu: {rdf.Count().GetValue()}")

rdf = rdf.Define("jpsis_hlt", "abs(mm_kin_mass-3.09)<0.06 && mm_kin_vtx_prob>0.1")
# simple_dump(rdf, "jpsis_hlt")

# rdf = rdf.Filter("Sum(jpsis_hlt)>0")
# print(f"Number of events with reconstructed Jpsis: {rdf.Count().GetValue()}")

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
rdf = rdf.Define("jpsi_probes",       "Concatenate(jpsi_probe1, jpsi_probe2)")
rdf = rdf.Define("jpsi_tests",        "Concatenate(jpsi_test1,  jpsi_test2)")
rdf = rdf.Define("jpsi_mu_pt",        "Concatenate(jpsi_mu1_pt, jpsi_mu2_pt)")
rdf = rdf.Define("jpsi_mu_eta",       "Concatenate(jpsi_mu1_eta, jpsi_mu2_eta)")

# dump(rdf, ["jpsi_probe1", "jpsi_tag1", "jpsi_probe2", "jpsi_tag2", "jpsi_tags", "jpsi_probes"])

h_probe_eta0p8 = rdf.Define("probe_mask", "jpsi_probes>0 && abs(jpsi_mu_eta)<0.8").Define("pt", "jpsi_mu_pt[probe_mask]").Histo1D(("h_probe_eta0p8","", 36, 2, 10), "pt")
h_test_eta0p8  = rdf.Define("test_mask", "jpsi_probes>0 && abs(jpsi_mu_eta)<0.8 && jpsi_tests>0").Define("pt", "jpsi_mu_pt[test_mask]").Histo1D(("h_test_eta0p8","", 36, 2, 10), "pt")

h_probe_eta1p2 = rdf.Define("probe_mask", "jpsi_probes>0 && abs(jpsi_mu_eta)>0.8 && abs(jpsi_mu_eta)<1.2").Define("pt", "jpsi_mu_pt[probe_mask]").Histo1D(("h_probe_eta1p2","", 36, 2, 10), "pt")
h_test_eta1p2  = rdf.Define("test_mask", "jpsi_probes>0 && abs(jpsi_mu_eta)>0.8 && abs(jpsi_mu_eta)<1.2 && jpsi_tests>0").Define("pt", "jpsi_mu_pt[test_mask]").Histo1D(("h_test_eta1p2","", 36, 2, 10), "pt")

h_probe_eta2p1 = rdf.Define("probe_mask", "jpsi_probes>0 && abs(jpsi_mu_eta)>1.2 && abs(jpsi_mu_eta)<2.1").Define("pt", "jpsi_mu_pt[probe_mask]").Histo1D(("h_probe_eta2p1","", 36, 2, 10), "pt")
h_test_eta2p1  = rdf.Define("test_mask", "jpsi_probes>0 && abs(jpsi_mu_eta)>1.2 && abs(jpsi_mu_eta)<2.1 && jpsi_tests>0").Define("pt", "jpsi_mu_pt[test_mask]").Histo1D(("h_test_eta2p1","", 36, 2, 10), "pt")


ROOT.RDF.RunGraphs([h_probe_eta0p8, h_test_eta0p8, h_probe_eta1p2, h_test_eta1p2,
                    h_probe_eta2p1, h_test_eta2p1])

h_eff_eta0p8 = make_trigger_efficiency_plots("h_eff_eta0p8", h_test_eta0p8, h_probe_eta0p8)
h_eff_eta1p2 = make_trigger_efficiency_plots("h_eff_eta1p2", h_test_eta1p2, h_probe_eta1p2)
h_eff_eta2p1 = make_trigger_efficiency_plots("h_eff_eta2p1", h_test_eta2p1, h_probe_eta2p1)

c1 = ROOT.TCanvas("c1","", 1500, 500)
c1.Divide(3,1)
c1.cd(1)
h_eff_eta0p8.Draw("e0")
c1.cd(2)
h_eff_eta1p2.Draw("e0")
c1.cd(3)
h_eff_eta2p1.Draw("e0")
