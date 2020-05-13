from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def nanoAOD_customizeBxToMuMu(process):
    process.load('Bmm5.NanoAOD.BxToMuMu_cff')
    # Data 
    process.nanoSequence   = cms.Sequence( process.nanoSequence + process.BxToMuMuSequence + process.BxToMuMuTables )
    # MC
    process.nanoSequenceMC = cms.Sequence( process.nanoSequenceMC + process.BxToMuMuMcSequence + process.BxToMuMuMcTables )
    return process

def nanoAOD_customizeV0ForMuonFake(process):
    process.load('Bmm5.NanoAOD.V0ForMuonFake_cff')
    # Data 
    process.nanoSequence   = cms.Sequence( process.nanoSequence + process.V0ForMuonFakeSequence + process.V0ForMuonFakeTables)
    # MC
    process.nanoSequenceMC = cms.Sequence( process.nanoSequenceMC + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    return process
