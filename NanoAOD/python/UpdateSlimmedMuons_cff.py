import FWCore.ParameterSet.Config as cms
slimmedMuons = cms.EDProducer(
    "MuonWithSoftMvaProducer",
    input = cms.InputTag("slimmedMuons", processName=cms.InputTag.skipCurrentProcess()),
    softMvaTrainingFile = cms.FileInPath("RecoMuon/MuonIdentification/data/TMVA-muonid-bmm4-B-25.weights.xml"),
    computeMuonIDMVA = cms.bool(True),
    mvaIDTrainingFile = cms.FileInPath("Bmm5/NanoAOD/data/mvaID.onnx"),
)
