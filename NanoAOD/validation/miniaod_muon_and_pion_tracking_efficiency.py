import ROOT
import sys
from DataFormats.FWLite import Events, Handle
from math import *
import math

from ROOT import reco
ROOT.FWLiteEnabler.enable()

def isTrackerSeeded(trk):
    # reco.TrackBase enums are exposed in Python as reco.TrackBase.<algoName>
    # print(trk.algoMaskUL() & 0x7FFFFFFF)
    mask = trk.algoMaskUL()
    return (
        (mask & reco.TrackBase.ctf) or
        (mask & reco.TrackBase.initialStep) or
        (mask & reco.TrackBase.lowPtTripletStep) or
        (mask & reco.TrackBase.pixelPairStep) or
        (mask & reco.TrackBase.detachedTripletStep) or
        (mask & reco.TrackBase.mixedTripletStep) or
        (mask & reco.TrackBase.pixelLessStep) or
        (mask & reco.TrackBase.tobTecStep) or
        (mask & reco.TrackBase.jetCoreRegionalStep) or
        (mask & reco.TrackBase.lowPtQuadStep) or
        (mask & reco.TrackBase.highPtTripletStep) or
        (mask & reco.TrackBase.detachedQuadStep)
    )

def deltaPhi(phi1,phi2):
    return acos(cos(phi2-phi1))

def deltaR(p1,p2):
    return sqrt(pow(deltaPhi(p1.phi(),p2.phi()),2)+pow(p2.eta()-p1.eta(),2))


events = Events ([
    # '/eos/cms/store/user/dmytro/tmp/store+relval+CMSSW_14_0_18+RelValSinglePiFlatPt0p7To10+MINIAODSIM+140X_mcRun3_2022_realistic_v12_STD_noPU_2022_reMC-v1+2590000+2b7aafd3-09f5-479a-b47a-2e2d61dd5d86.root'
    # '/eos/cms/store/user/dmytro/tmp/store+relval+CMSSW_14_0_18+RelValSingleMuPt10+MINIAODSIM+140X_mcRun3_2023_realistic_v9_STD_noPU_2023_reMC-v1+2590000+ba833f2d-5e6a-45b6-a94e-f2a60856a57a.root'
    '/eos/cms/store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv3+BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+124X_mcRun3_2022_realistic_postEE_v1-v2+2550000+0f450e50-38cc-4dd3-a874-26b4a956b959.root'
])

handlePruned  = Handle ("std::vector<reco::GenParticle>")
labelPruned = ("prunedGenParticles")
handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
labelPacked = ("packedGenParticles")

muonHandle, muonLabel = Handle("std::vector<pat::Muon>"),"slimmedMuons"
pf_handle = Handle("std::vector<pat::PackedCandidate>")
pf_label  = ("packedPFCandidates")

# loop over events
n_gen_pions = [0, 0, 0]
n_gen_muons = [0, 0, 0]
n_pf_pions = [0, 0, 0]
n_pf_muons = [0, 0, 0]
n_muons_pfcands = [0, 0, 0]

min_pt = 5.0

def match(gen_particles, pf_cand):
    for gen in gen_particles:
        if deltaR(pf, gen) > 0.1: continue
        if abs(pf.pt() / gen.pt() - 1) > 0.02: continue
        return gen
    return None

def get_eta_bin(eta):
    eta_bin = 0
    if abs(eta) > 0.8 and abs(eta) < 1.6:
        eta_bin = 1
    elif abs(eta) > 1.6 and abs(eta) < 2.4:
        eta_bin = 2
    elif abs(eta) > 2.4:
        eta_bin = 3
    return eta_bin

for i, event in enumerate(events):
    # if i > 1000: break
    ### Get collections
	
    # Muons
    event.getByLabel(muonLabel, muonHandle)
    muons = muonHandle.product()

    # PFCandidates
    event.getByLabel(pf_label, pf_handle)
    pfs = pf_handle.product()

    # Gen particles
    event.getByLabel (labelPacked, handlePacked)
    packed = handlePacked.product()
    event.getByLabel (labelPruned, handlePruned)
    pruned = handlePruned.product()
    
    gen_mus = []
    gen_pis = []
    for pa in packed:
        if pa.pt() < min_pt:
            continue
        eta_bin = get_eta_bin(pa.eta())
        if eta_bin >= len(n_gen_muons):
            continue
        if abs(pa.pdgId()) == 13:
            n_gen_muons[eta_bin] += 1
            gen_mus.append(pa)
        if abs(pa.pdgId()) == 211:
            n_gen_pions[eta_bin] += 1
            gen_pis.append(pa)

    for pf in pfs:
        if pf.charge() == 0:
            continue
        if not pf.hasTrackDetails():
            continue
        if not isTrackerSeeded(pf.bestTrack()):
            continue
        
        gen_mu = match(gen_mus, pf)
        if gen_mu != None:
            eta_bin = get_eta_bin(gen_mu.eta())
            n_pf_muons[eta_bin] += 1
            if abs(pf.pdgId()) == 13:
                n_muons_pfcands[eta_bin] += 1
        gen_pi = match(gen_pis, pf)
        if gen_pi != None:
            eta_bin = get_eta_bin(gen_pi.eta())
            n_pf_pions[eta_bin] += 1

for eta_bin in range(len(n_pf_muons)):
    if eta_bin == 0:
        print("|eta| in [0.0, 0.8]")
    elif eta_bin == 1:
        print("\n|eta| in [0.8, 1.6]")
    else:
        print("\n|eta| in [1.6, 2.4]")
    print(f"n_gen_muons: {n_gen_muons[eta_bin]}")
    print(f"n_pf_muons: {n_pf_muons[eta_bin]}")
    if n_gen_muons[eta_bin] > 0:
        p = float(n_pf_muons[eta_bin]) / n_gen_muons[eta_bin]
        print(f"eff: ${p * 100:0.2f} \pm {math.sqrt(p * (1-p)/n_gen_muons[eta_bin]) * 100:0.2f}$")
    print(f"n_muons_pfcands: {n_muons_pfcands[eta_bin]}\n")
    
    print(f"n_gen_pions: {n_gen_pions[eta_bin]}")
    print(f"n_pf_pions: {n_pf_pions[eta_bin]}")
    if n_gen_pions[eta_bin] > 0:
        p = float(n_pf_pions[eta_bin]) / n_gen_pions[eta_bin]
        print(f"eff: ${p * 100:0.2f} \pm {math.sqrt(p * (1-p)/n_gen_pions[eta_bin]) * 100:0.2f}$")


