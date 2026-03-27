import ROOT
import glob
from collections import defaultdict
import os

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/535/fit-bkmm/" + \
    "InclusiveDileptonMinBias_Fil-DoubleMuOS43_TuneCP5Plus_13p6TeV_pythia8+RunIII2024Summer24MiniAOD-Filter_DoubleMuOS43_140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/"

pdgIds = {
    511:'B0', 521:'B+', 531:'Bs', 541:'Bc',
    411:'D0', 421:'D+', 431:'Ds', 443:'Jpsi',
    10313:'K1(1270)', 10323:'K1(1270)+',
    321:'K+', 313:'K*', 333:'phi', 323:'K*+',
    315:'K*2(1430)', 331:"eta'", 221:'eta',
    2212:'p', 211:'pi+', 111:'pi0', 213:'rho+',
    22:'gamma', 310:'K_S'
}

def decipher(signature):
    result = ""
    reminder = signature
    if signature > 0:
        for id, name in pdgIds.items():
            while reminder % id == 0:
                if result != "":
                    result += " "
                result += name
                reminder = reminder / id
    return (result, reminder)

def name(pdgId):
    id = abs(pdgId)
    if id not in pdgIds:
        if id not in unknown_pdgIds:
            unknown_pdgIds.append(id)
        return "other"
    return pdgIds[id]

def report(signature, count, fraction):
    parent, parent_reminder = decipher(signature[0])
    decay, decay_reminder = decipher(signature[1])
    if decay_reminder not in (0, 1):
        decay += f" {decay_reminder}"
    if len(parent) > 0 and len(decay) > 0:
        print(f"{parent} -> {decay:30s} {count} ({100.0 * fraction:0.1f}%)")


chain = ROOT.TChain("inclusiveMc")
for f in glob.glob("%s/*.root" % input_path):
# for f in glob.glob("%s/10*.root" % input_path):
    chain.Add(f)

all_processes = defaultdict(int)
all_processes_bmass = defaultdict(int)
all_processes_left_sideband = defaultdict(int)
all_processes_transition = defaultdict(int)

nevents = chain.GetEntries()
nevents_bmass = 0
nevents_left_sideband = 0
nevents_transition = 0

print(f"N selected candidates: {nevents}")

for event in chain:
    signature = (abs(event.decay_parent), abs(event.decay_signature))
    all_processes[signature] += 1

    if abs(event.m-5.28) < 0.1:
        nevents_bmass += 1
        all_processes_bmass[signature] += 1
    elif event.m < 5.15:
        nevents_left_sideband += 1
        all_processes_left_sideband[signature] += 1
    elif event.m < 5.18:
        nevents_transition += 1
        all_processes_transition[signature] += 1
        
        
print("All processes")
for signature, count in sorted(all_processes.items(), key=lambda x: x[1], reverse=True):
    if count > 10:
        report(signature, count, float(count) / nevents)

print(f"\nMass peak abs(m-5.28)<0.1: {nevents_bmass}")
shown_fraction = 0
for signature, count in sorted(all_processes_bmass.items(), key=lambda x: x[1], reverse=True):
    fraction = float(count) / nevents_bmass
    if fraction > 0.001:
        shown_fraction += fraction
        report(signature, count, fraction)
print(f"Other {100 * (1 - shown_fraction):0.1f}%")

print(f"\nLeft sideband: {nevents_left_sideband}")
shown_fraction = 0
for signature, count in sorted(all_processes_left_sideband.items(), key=lambda x: x[1], reverse=True):
    fraction = float(count) / nevents_left_sideband
    if fraction > 0.001:
        shown_fraction += fraction
        report(signature, count, fraction)
print(f"Other {100 * (1 - shown_fraction):0.1f}%")

print(f"\nTransition [5.15, 5.18]: {nevents_transition}")
shown_fraction = 0
for signature, count in sorted(all_processes_transition.items(), key=lambda x: x[1], reverse=True):
    fraction = float(count) / nevents_transition
    if fraction > 0.001:
        shown_fraction += fraction
        report(signature, count, fraction)
print(f"Other {100 * (1 - shown_fraction):0.1f}%")
