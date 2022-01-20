#!/bin/env python
from ROOT import TFile,TTree

f = TFile("/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/0C239543-D4D7-C948-8F0B-8DFB8ADF369B.root")
t = f.Get("Events")

print "Number of events:", t.GetEntries()
n = 0
for event in t:
    matched = False
    for pdgId in event.genbmm_pdgId:
        if abs(pdgId) == 531:
            matched = True
            break
    if matched: continue
    
    print "Found unmatched event: %u \tlumi: %u" % (event.event, event.luminosityBlock)
    n += 1
    for i, pt in enumerate(event.GenPart_pt):
        mother = 0
        if event.GenPart_genPartIdxMother[i] >= 0:
            mother = event.GenPart_pdgId[event.GenPart_genPartIdxMother[i]]
        if event.GenPart_status[i] == 1:
            print "i: %2u \tpt: %5.1f \tpdgId: %d \tmother: %d \tstatus: %d" % (i, pt, event.GenPart_pdgId[i], mother, event.GenPart_status[i])

    if n > 10: break
