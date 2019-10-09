import FWCore.ParameterSet.Config as cms
BxToMuMuFilter = cms.EDFilter("BxFilter",
    src = cms.InputTag("BxToMuMu","DiMuon"),
    minNumber = cms.uint32( 1 )
)

BxToMuMuFilterMc = cms.EDFilter("BxFilter",
    src = cms.InputTag("BxToMuMuMc","DiMuon"),
    minNumber = cms.uint32( 1 )
)

BxToMuMuFilterSequence   = cms.Sequence(BxToMuMuFilter)
BxToMuMuFilterSequenceMC = cms.Sequence(BxToMuMuFilterMc)
