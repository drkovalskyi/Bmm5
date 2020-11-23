from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def nanoAOD_customizeBxToMuMu(process):
    process.load('Bmm5.NanoAOD.BxToMuMu_cff')
    process.load('Bmm5.NanoAOD.UpdateSlimmedMuons_cff')
    process.load('PhysicsTools.NanoAOD.muons_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.BxToMuMuSequence + process.BxToMuMuTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.BxToMuMuMcSequence + process.BxToMuMuMcTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
    return process

def nanoAOD_customizeV0ForMuonFake(process):
    process.load('Bmm5.NanoAOD.BmmV0ForMuonFake_cff')
    process.load('Bmm5.NanoAOD.UpdateSlimmedMuons_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.V0ForMuonFakeSequence + process.V0ForMuonFakeTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
    return process
