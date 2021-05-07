from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms
import re

# NOTE: 
#    All instances of FlatTableProducers must end with Table in their
#    names so that their product match the keep patterns in the default
#    event content. Otherwise you need to modify outputCommands in
#    NanoAODEDMEventContent or provide a custom event content to the
#    output module

# can use cms.PSet.clone() method instead
def merge_psets(*argv):
    result = cms.PSet()
    for pset in argv:
        if isinstance(pset, cms._Parameterizable):
            for name in pset.parameters_().keys():
                value = getattr(pset,name)
                type = value.pythonTypeName()
                setattr(result,name,value)
    return result

BmmMuonId = cms.EDProducer("BmmMuonIdProducer",
    muonCollection = cms.InputTag("linkedObjects","muons"),
    prunedGenParticleCollection = cms.InputTag("prunedGenParticles"),
    packedGenParticleCollection = cms.InputTag("packedGenParticles"),
    softMuonMva = cms.FileInPath('Bmm5/NanoAOD/data/muon_mva/Run2018-20210430-2004-Event0.model'),
    isMC = cms.bool(False)
)

BmmMuonIdMc = BmmMuonId.clone( isMC = cms.bool(True) ) 

BmmMuonIdVariables = cms.PSet(
    trkKink             = Var("userFloat('trkKink')",             float, doc = "Inner track kink chi2"),
    glbTrackProbability = Var("userFloat('glbTrackProbability')", float, doc = "Log probability of the global fit"),
    chi2LocalPosition   = Var("userFloat('chi2LocalPosition')",   float, doc = "chi2 for STA-TK matching by local position"),
    glbNormChi2         = Var("userFloat('glbNormChi2')",         float, doc = "Normalized chi2 of the global fit"),
    trkValidFrac        = Var("userFloat('trkValidFrac')",        float, doc = "Fraction of valid hits for inner track"),
    
    match1_dX           = Var("userFloat('match1_dX')",           float, doc = "Station 1 local segment-track dX"),
    match1_pullX        = Var("userFloat('match1_pullX')",        float, doc = "Station 1 local segment-track dX/dErr"),
    match1_pullDxDz     = Var("userFloat('match1_pullDxDz')",     float, doc = "Station 1 local segment-track direction matching in x"),
    match1_dY           = Var("userFloat('match1_dY')",           float, doc = "Station 1 local segment-track dY"),
    match1_pullY        = Var("userFloat('match1_pullY')",        float, doc = "Station 1 local segment-track dY/dErr"),
    match1_pullDyDz     = Var("userFloat('match1_pullDyDz')",     float, doc = "Station 1 local segment-track direction matching in y"),
    match2_dX           = Var("userFloat('match2_dX')",           float, doc = "Station 2 local segment-track dX"),
    match2_pullX        = Var("userFloat('match2_pullX')",        float, doc = "Station 2 local segment-track dX/dErr"),
    match2_pullDxDz     = Var("userFloat('match2_pullDxDz')",     float, doc = "Station 2 local segment-track direction matching in x"),
    match2_dY           = Var("userFloat('match2_dY')",           float, doc = "Station 2 local segment-track dY"),
    match2_pullY        = Var("userFloat('match2_pullY')",        float, doc = "Station 2 local segment-track dY/dErr"),
    match2_pullDyDz     = Var("userFloat('match2_pullDyDz')",     float, doc = "Station 2 local segment-track direction matching in y"),

    nPixels             = Var("userInt('nPixels')",                 int, doc = "Number of valid pixel hits"),
    nValidHits          = Var("userInt('nValidHits')",              int, doc = "Number of valid hits"),
    nLostHitsInner      = Var("userInt('nLostHitsInner')",          int, doc = "Number of lost hits before tracker track"),
    nLostHitsOn         = Var("userInt('nLostHitsOn')",             int, doc = "Number of lost hits on tracker track"),
    nLostHitsOuter      = Var("userInt('nLostHitsOuter')",          int, doc = "Number of lost hits after tracker track"),
    
    trkLayers           = Var("userInt('trkLayers')",               int, doc = "Number of layers with measurements in tracker track"),
    trkLostLayersInner  = Var("userInt('trkLostLayersInner')",      int, doc = "Number of lost layers befor tracker track"),
    trkLostLayersOn     = Var("userInt('trkLostLayersOn')",         int, doc = "Number of lost layers on tracker track"),
    trkLostLayersOuter  = Var("userInt('trkLostLayersOuter')",      int, doc = "Number of lost layers after tracker track"),

    highPurity          = Var("userInt('highPurity')",              int, doc = "High purity inner track"),
    newSoftMuonMva      = Var("userFloat('newSoftMuonMva')",      float, doc = "New softMuonMva"),

)

BmmMuonIdMcVariables = merge_psets(
    BmmMuonIdVariables,
    cms.PSet(
        simType             = Var("userInt('simType')",        int, doc = "reco::MuonSimType"),
        simExtType          = Var("userInt('simExtType')",     int, doc = "reco::ExtendedMuonSimType"),
        simPdgId            = Var("userInt('simPdgId')",       int, doc = "SIM particle pdgId"),
        simMotherPdgId      = Var("userInt('simMotherPdgId')", int, doc = "SIM particle mother pdgId"),
        simProdRho          = Var("userFloat('simProdRho')", float, doc = "SIM particle production vertex"),
        simProdZ            = Var("userFloat('simProdZ')",   float, doc = "SIM particle production vertex"),
        ),
)

BmmMuonIdTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BmmMuonId","muons"),
    cut=cms.string(""),
    name=cms.string("MuonId"),
    doc=cms.string("Muon Id Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BmmMuonIdVariables
)
BmmMuonIdMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BmmMuonIdMc","muons"),
    cut=cms.string(""),
    name=cms.string("MuonId"),
    doc=cms.string("Muon Id Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BmmMuonIdMcVariables
)

BmmMuonIdSequence   = cms.Sequence(BmmMuonId)
BmmMuonIdMcSequence = cms.Sequence(BmmMuonIdMc)
BmmMuonIdTables     = cms.Sequence(BmmMuonIdTable)
BmmMuonIdMcTables   = cms.Sequence(BmmMuonIdMcTable)
