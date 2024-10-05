import ROOT
import glob
import pprint

prefix = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/"
input_paths = [
    prefix + "DstarToD0Pi_D0To2Pi_PiFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v1+MINIAODSIM/",
]

types = {
    999: 'Unknown',
    0: 'NotMatched',
    1: 'MatchedPunchthrough',
    11: 'MatchedElectron',
    4: 'MatchedPrimaryMuon',
    3: 'MatchedMuonFromHeavyFlavour',
    2: 'MatchedMuonFromLightFlavour',
    -1: 'GhostPunchthrough',
    -11: 'GhostElectron',
    -4: 'GhostPrimaryMuon',
    -3: 'GhostMuonFromHeavyFlavour',
    -2: 'GhostMuonFromLightFlavour'
    }

chain = ROOT.TChain("Events")
for input_path in input_paths:
    for f in glob.glob("%s/*.root" % input_path):
        chain.Add(f)
print("Number of events:", chain.GetEntries())

results = dict()

def save_entry(name, sim_type):
    if name not in results:
        results[name] = dict()
    if sim_type not in results[name]:
        results[name][sim_type] = 0
    results[name][sim_type] += 1

for event in chain:
    for i_mu in range(event.nMuon):
        # Keep only muons of interest
        if event.Muon_pt[i_mu] < 4.0: continue

        # Ignore not identified muon types due to lack of pileup matching
        sim_type = event.MuonId_simType[i_mu]
        if sim_type == 0: continue
        
        if event.Muon_looseId[i_mu] and event.Muon_isGlobal[i_mu] and event.Muon_isTracker[i_mu]:
            save_entry("D0mm selection", sim_type)
            
        if event.Muon_looseId[i_mu]:
            save_entry("Loose ID", sim_type)

        if event.Muon_mediumId[i_mu]:
            save_entry("Medium ID", sim_type)

for id in sorted(results):
    print(id)
    info = sorted(results[id].items(), key=lambda item: item[1], reverse=True)
    total_sum = sum(v1 for k1, v1 in info) 
    for k1, v1 in info:
        print("\t%s: %u (%0.2f%%)" % (types[k1], v1, 100.*v1/total_sum))
