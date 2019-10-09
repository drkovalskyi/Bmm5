from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def nanoAOD_customizeBxToMuMu(process):
    process.load('Bmm5.NanoAOD.BxToMuMu_cff')
    # Data 
    process.nanoSequence   = cms.Sequence( process.nanoSequence + process.BxToMuMuSequence + process.BxToMuMuTables)
    # MC
    process.nanoSequenceMC = cms.Sequence( process.nanoSequenceMC + process.BxToMuMuMcSequence + process.BxToMuMuMcTables)
    return process
