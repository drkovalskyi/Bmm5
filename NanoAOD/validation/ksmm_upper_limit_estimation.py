import ROOT, glob, os, re, pprint
from math import sqrt

selection = "abs(m-0.45)<0.1 && mm_kin_slxy>10 && mm_kin_lxy>1 && mm_kin_alpha<0.001"
signal_box_3sigma = "abs(m-0.498) < 3 * 6.22730e-03"
signal_box_1sigma = "abs(m-0.498) < 6.22730e-03"

def LoadChain(name, sample):
    chain = ROOT.TChain(name)
    for entry in sample:
        chain.Add(entry)
    return chain

def GetInfo(name, sample, selections):
    results = dict()

    info = LoadChain("info", sample)
    results["n_gen_all"] = None
    results["n_processed"] = 0
    for entry in info:
        results["n_processed"] += entry.n_processed
        if hasattr(entry, "n_gen_all"):
            if results["n_gen_all"] == None:
                results["n_gen_all"] = entry.n_gen_all
            else:
                results["n_gen_all"] = entry.n_gen_all

    events = LoadChain(name, sample)
    results["n_preselected"] = events.GetEntries()

    for name, selection in selections.items():
        results[name] = events.GetEntries(selection)

    return results


path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/529/ksmm/"
mc = [
    path + "K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v1+MINIAODSIM/*root",
    path + "K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v1+MINIAODSIM/*.root",
]

mc_results = GetInfo("ksmmMc", mc, {
    '3sigma': selection + "&&" + signal_box_3sigma,
    '1sigma': selection + "&&" + signal_box_1sigma
})
pprint.pprint(mc_results)

data = [ ]
for f in glob.glob("%s/*" % path):
    if not re.search("ParkingDoubleMuonLowMass", f): continue
    data.append("%s/*.root" % f)

data_results = GetInfo("ksmmData", data, {
    'sideband': selection + "&&" + signal_box_3sigma + "&& !(" + signal_box_1sigma + ")",
})

pprint.pprint(data_results)

print("Expect number of background events: %0.1f" % (data_results["sideband"]/(3-1)))
print("Expect upper limit on the signal event count in the Gaussian limit: %0.1f" % (sqrt(data_results["sideband"]/(3-1))*1.64))
