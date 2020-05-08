from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def nanoAOD_customizeBxToMuMu(process):
    process.load('Bmm5.NanoAOD.BxToMuMu_cff')
    process.load('Bmm5.NanoAOD.V0ForMuonFake_cff')
    # Data 
    process.nanoSequence   = cms.Sequence( process.nanoSequence + process.BxToMuMuSequence + process.BxToMuMuTables + process.V0ForMuonFakeSequence + process.V0ForMuonFakeTables)
    # MC
    process.nanoSequenceMC = cms.Sequence( process.nanoSequenceMC + process.BxToMuMuMcSequence + process.BxToMuMuMcTables + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    return process
