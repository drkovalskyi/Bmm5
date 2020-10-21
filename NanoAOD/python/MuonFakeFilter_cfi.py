import FWCore.ParameterSet.Config as cms
muonFakeFilter = cms.EDFilter(
    "MuonFakeFilter",
    muonCollection = cms.InputTag("slimmedMuons"),
    prunedGenParticleCollection = cms.InputTag("prunedGenParticles"),
    minPt = cms.double( 4.0 ),
    maxEta = cms.double( 1.4 ),
    maxDR = cms.double( 0.1 ),
    maxDPt = cms.double( 0.02 )
)
