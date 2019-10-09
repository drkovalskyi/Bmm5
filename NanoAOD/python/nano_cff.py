from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *
def nanoAOD_customizeBxToMuMu(process):
    # process = nanoAOD_customizeCommon(process)
    # Data 
    process.nanoSequence   = cms.Sequence( process.nanoSequence + BxToMuMuSequence + BxToMuMuTables)
    # MC
    process.nanoSequenceMC = cms.Sequence( process.nanoSequenceMC + BxToMuMuMcSequence + BxToMuMuMcTables)
    return process
