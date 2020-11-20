from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms

# NOTE: 
#    All instances of FlatTableProducers must end with Table in their
#    names so that their product match keep patterns in the default
#    event content. Otherwise you need to modify outputCommands in
#    NanoAODEDMEventContent or provide a custom event content to the
#    output module

def merge_psets(*argv):
    result = cms.PSet()
    for pset in argv:
        if isinstance(pset, cms._Parameterizable):
            for name in pset.parameters_().keys():
                value = getattr(pset,name)
                type = value.pythonTypeName()
                setattr(result,name,value)
    return result

V0ForMuonFake = cms.EDProducer(
    "V0Producer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection = cms.InputTag("packedPFCandidates"),
    packedGenParticleCollection = cms.InputTag("packedGenParticles"),
    minMuonPt  = cms.double(3.5),
    maxMuonEta = cms.double(1.4),
    minPionPt  = cms.double(1.0),
    maxPionEta = cms.double(1.4),
    minKsMass  = cms.double(0.45),
    maxKsMass  = cms.double(0.55),
    minKsPreselectMass = cms.double(0.4),
    maxKsPreselectMass = cms.double(0.6),
    maxTwoTrackDOCA = cms.double(0.1),
    maxLxy = cms.double(999),
    minSigLxy = cms.double(5),
    minVtxProb = cms.double(0.001),
    minCosAlpha = cms.double(0.9),
    minDisplaceTrackSignificance = cms.double(1),
    isMC = cms.bool(False)
)

V0ForMuonFakeMC = V0ForMuonFake.clone( isMC = cms.bool(True) ) 

V0ForMuonFakeVariables = cms.PSet(
    mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks"),
    trk1_pt      = Var("userFloat('trk1_pt')",         float, doc = "Track 1 pt"),
    trk1_eta     = Var("userFloat('trk1_eta')",        float, doc = "Track 1 eta"),
    trk1_phi     = Var("userFloat('trk1_phi')",        float, doc = "Track 1 phi"),
    trk1_mu_index = Var("userInt('trk1_mu_index')",      int, doc = "Matched muon index for track 1"),
    trk2_pt      = Var("userFloat('trk2_pt')",         float, doc = "Track 2 pt"),
    trk2_eta     = Var("userFloat('trk2_eta')",        float, doc = "Track 2 eta"),
    trk2_phi     = Var("userFloat('trk2_phi')",        float, doc = "Track 2 phi"),
    trk2_mu_index = Var("userInt('trk2_mu_index')",      int, doc = "Matched muon index for track 2"),
    trk1_sip     = Var("userFloat('trk1_sip')",        float, doc = "Track 1 2D impact parameter significance wrt Beam Spot"),
    trk2_sip     = Var("userFloat('trk2_sip')",        float, doc = "Track 2 2D impact parameter significance wrt Beam Spot"),
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability"),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2"),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass"),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt"),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta"),
    kin_phi      = Var("userFloat('kin_phi')",         float, doc = "Kinematic fit: vertex refitted phi"),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error"),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot"),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot"),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS"),
)

V0ForMuonFakeVariablesMC = merge_psets(
    V0ForMuonFakeVariables,
    cms.PSet(
        gen_trk1_pdgId  = Var("userInt(  'gen_trk1_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_trk1_mpdgId = Var("userInt(  'gen_trk1_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_trk1_pt     = Var("userFloat('gen_trk1_pt')",     float,   doc = "Gen match: first track pt"),
        gen_trk2_pdgId  = Var("userInt(  'gen_trk2_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_trk2_mpdgId = Var("userInt(  'gen_trk2_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_trk2_pt     = Var("userFloat('gen_trk2_pt')",     float,   doc = "Gen match: second track pt"),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt"),
        ),
)

V0ForMuonFakeTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFake","Ks"),
    cut=cms.string(""),
    name=cms.string("ks"),
    doc=cms.string("Ks Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = V0ForMuonFakeVariables
)

V0ForMuonFakeMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFakeMC","Ks"),
    cut=cms.string(""),
    name=cms.string("ks"),
    doc=cms.string("Ks Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = V0ForMuonFakeVariablesMC
)

V0ForMuonFakeSequence   = cms.Sequence(V0ForMuonFake)
V0ForMuonFakeMcSequence = cms.Sequence(V0ForMuonFakeMC)
V0ForMuonFakeTables     = cms.Sequence(V0ForMuonFakeTable)
V0ForMuonFakeMcTables   = cms.Sequence(V0ForMuonFakeMcTable)
