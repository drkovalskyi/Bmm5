from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def nanoAOD_customizeDileptonPlusX(process):
    process.load('Bmm5.NanoAOD.DileptonPlusX_cff')
    process.load('Bmm5.NanoAOD.UpdateSlimmedMuons_cff')
    process.load('PhysicsTools.NanoAOD.muons_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.DileptonPlusXSequence + process.DileptonPlusXTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.DileptonPlusXMcSequence + process.DileptonPlusXMcTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)

    # enforce process name
    # process.load('PhysicsTools.NanoAOD.globals_cff')
    process.genFilterTable.src = cms.InputTag("genFilterEfficiencyProducer")

    # keep all genparticles
    process.finalGenParticles.select = cms.vstring("keep *")
    
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

def nanoAOD_customizeBmmMuonId(process):
    process.load('Bmm5.NanoAOD.BmmMuonId_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.nanoSequence + process.BmmMuonIdSequence + process.BmmMuonIdTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.nanoSequenceMC + process.BmmMuonIdMcSequence + process.BmmMuonIdMcTables)

    process.load('PhysicsTools.NanoAOD.muons_cff')
    # process.muonTable.variables.mvaMuID = Var("mvaIDValue()",float,doc="MVA-based ID score ",precision=6)
    return process

def nanoAOD_keepLowPtMuons(process):
    process.muonTable.doc = cms.string("slimmedMuons after basic selection (pt > 2 || (pt > 2 && (passed(\'CutBasedIdLoose\') || passed(\'SoftCutBasedId\') || passed(\'SoftMvaId\') || passed(\'CutBasedIdGlobalHighPt\') || passed(\'CutBasedIdTrkHighPt\'))))")

    process.finalMuons.cut = cms.string("pt > 2 || (pt > 2 && (passed(\'CutBasedIdLoose\') || passed(\'SoftCutBasedId\') || passed(\'SoftMvaId\') || passed(\'CutBasedIdGlobalHighPt\') || passed(\'CutBasedIdTrkHighPt\')))")

    return process
