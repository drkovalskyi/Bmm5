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

triggers = [
    'HLT_Dimuon0_LowMass_L1_0er1p5',
    'HLT_DoubleMu4_3_Bs',
    'HLT_DoubleMu4_3_Jpsi',
    'HLT_DoubleMu4_3_Jpsi_Displaced',
    'HLT_DoubleMu4_Jpsi_Displaced',
    'HLT_DoubleMu4_Jpsi_NoVertexing',
    'HLT_Dimuon6_Jpsi_NoVertexing',
    'HLT_Dimuon0_LowMass',
    'HLT_DoubleMu0',
    'HLT_IsoMu24',
    'HLT_IsoMu27',
    'HLT_Mu8',
    'HLT_Mu17'
]

BmmMuonId = cms.EDProducer(
    "BmmMuonIdProducer",
    muonCollection = cms.InputTag("linkedObjects","muons"),
    prunedGenParticleCollection = cms.InputTag("prunedGenParticles"),
    packedGenParticleCollection = cms.InputTag("packedGenParticles"),
    trigger = cms.InputTag("slimmedPatTrigger"),
    softMuonMva = cms.FileInPath('Bmm5/NanoAOD/data/muon_mva/Run2018-20210430-2004-Event0.model'),
    isMC = cms.bool(False),
    triggers = cms.vstring(triggers),
    triggerCollection = cms.string("hltIterL3MuonCandidates"),
    l1Src = cms.InputTag("gmtStage2Digis:Muon"),
)

from Configuration.Eras.Modifier_run2_muon_2016_cff import run2_muon_2016
run2_muon_2016.toModify(BmmMuonId,
    triggerCollection = cms.string("hltL3MuonCandidates")
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
    hlt_pt              = Var("userFloat('hlt_pt')",              float, doc = "HLT pt"),
    hlt_dr              = Var("userFloat('hlt_dr')",              float, doc = "HLT dR"),
    l1_pt               = Var("userFloat('l1_pt')",               float, doc = "L1 pt"),
    l1_mpt              = Var("userFloat('l1_mpt')",              float, doc = "L1 pt (cross-check)"),
    l1_etaAtVtx         = Var("userFloat('l1_etaAtVtx')",         float, doc = "L1 eta at vertex"),
    l1_phiAtVtx         = Var("userFloat('l1_phiAtVtx')",         float, doc = "L1 phi at vertex"),
    l1_eta              = Var("userFloat('l1_eta')",              float, doc = "L1 eta"),
    l1_phi              = Var("userFloat('l1_phi')",              float, doc = "L1 phi"),
    l1_quality          = Var("userInt('l1_quality')",              int, doc = "L1 quality"),
)

for trigger in triggers:
    setattr(BmmMuonIdVariables, trigger, Var("userInt('%s')" % trigger, int, doc = "Used in trigger"))


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
