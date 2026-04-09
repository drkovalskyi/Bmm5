from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms

# Precision settings (mantissa bits, -1 = full float32)
full_precision   = -1   # mass, pt, eta, phi, vtx_prob
medium_precision = 12   # vertex pos/err, angles, DCA, isolation
low_precision    =  6   # not critical variables

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
    "BmmV0Producer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection = cms.InputTag("packedPFCandidates"),
    packedGenParticleCollection = cms.InputTag("packedGenParticles"),
    minMuonPt  = cms.double(2.0),
    maxMuonEta = cms.double(2.4),
    minPionPt  = cms.double(1.0),
    maxPionEta = cms.double(2.4),
    minKsMass  = cms.double(0.45),
    maxKsMass  = cms.double(0.55),
    minKsPreselectMass = cms.double(0.4),
    maxKsPreselectMass = cms.double(0.6),
    minPhiMass  = cms.double(1.00),
    maxPhiMass  = cms.double(1.04),
    minPhiPreselectMass = cms.double(0.9),
    maxPhiPreselectMass = cms.double(1.1),
    minDsMass  = cms.double(1.91),
    maxDsMass  = cms.double(2.03),
    minDstarPreselectMass = cms.double(1.8),
    maxDstarPreselectMass = cms.double(2.2),
    minDstarMass  = cms.double(1.9),
    maxDstarMass  = cms.double(2.1),
    minDsPreselectMass = cms.double(1.8),
    maxDsPreselectMass = cms.double(2.1),
    minD0Mass  = cms.double(1.8),
    maxD0Mass  = cms.double(1.9),
    minD0PreselectMass = cms.double(1.6),
    maxD0PreselectMass = cms.double(2.0),
    minLambdaMass  = cms.double(1.05),
    maxLambdaMass  = cms.double(1.15),
    minLambdaPreselectMass = cms.double(1.0),
    maxLambdaPreselectMass = cms.double(1.2),
    maxTwoTrackDOCA = cms.double(0.1),
    maxLxy = cms.double(999),
    minSigLxy = cms.double(5),
    minVtxProb = cms.double(0.001),
    minCosAlpha = cms.double(0.9),
    minDisplaceTrackSignificance = cms.double(1),
    isMC = cms.bool(False)
)

V0ForMuonFakeMC = V0ForMuonFake.clone( isMC = cms.bool(True) ) 

# KsToPiPi

KsForMuonFakeVariables = cms.PSet(
    mass         = Var("mass",                         float, doc = "Unfit invariant mass", precision = full_precision),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks", precision = medium_precision),
    trk1_pt      = Var("userFloat('trk1_pt')",         float, doc = "Track 1 pt", precision = full_precision),
    trk1_eta     = Var("userFloat('trk1_eta')",        float, doc = "Track 1 eta", precision = full_precision),
    trk1_phi     = Var("userFloat('trk1_phi')",        float, doc = "Track 1 phi", precision = full_precision),
    trk1_mu_index = Var("userInt('trk1_mu_index')",      int, doc = "Matched muon index for track 1"),
    trk2_pt      = Var("userFloat('trk2_pt')",         float, doc = "Track 2 pt", precision = full_precision),
    trk2_eta     = Var("userFloat('trk2_eta')",        float, doc = "Track 2 eta", precision = full_precision),
    trk2_phi     = Var("userFloat('trk2_phi')",        float, doc = "Track 2 phi", precision = full_precision),
    trk2_mu_index = Var("userInt('trk2_mu_index')",      int, doc = "Matched muon index for track 2"),
    trk1_sip     = Var("userFloat('trk1_sip')",        float, doc = "Track 1 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    trk2_sip     = Var("userFloat('trk2_sip')",        float, doc = "Track 2 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability", precision = full_precision),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2", precision = medium_precision),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass", precision = full_precision),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt", precision = full_precision),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta", precision = full_precision),
    kin_phi      = Var("userFloat('kin_phi')",         float, doc = "Kinematic fit: vertex refitted phi", precision = full_precision),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error", precision = medium_precision),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot", precision = full_precision),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot", precision = full_precision),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS", precision = full_precision),
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS", precision = medium_precision),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV", precision = medium_precision),
)

KsForMuonFakeVariablesMC = merge_psets(
    KsForMuonFakeVariables,
    cms.PSet(
        gen_trk1_pdgId  = Var("userInt(  'gen_trk1_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_trk1_mpdgId = Var("userInt(  'gen_trk1_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_trk1_pt     = Var("userFloat('gen_trk1_pt')",     float,   doc = "Gen match: first track pt", precision = full_precision),
        gen_trk2_pdgId  = Var("userInt(  'gen_trk2_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_trk2_mpdgId = Var("userInt(  'gen_trk2_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_trk2_pt     = Var("userFloat('gen_trk2_pt')",     float,   doc = "Gen match: second track pt", precision = full_precision),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass", precision = full_precision),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt", precision = full_precision),
        ),
)

KsForMuonFakeTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFake","Ks"),
    cut=cms.string(""),
    name=cms.string("ks"),
    doc=cms.string("Ks Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = KsForMuonFakeVariables
)

KsForMuonFakeMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFakeMC","Ks"),
    cut=cms.string(""),
    name=cms.string("ks"),
    doc=cms.string("Ks Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = KsForMuonFakeVariablesMC
)

# D0ToKPi

D0ForMuonFakeVariables = cms.PSet(
    mass         = Var("mass",                         float, doc = "Unfit invariant mass", precision = full_precision),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks", precision = medium_precision),
    kaon_pt      = Var("userFloat('kaon_pt')",         float, doc = "Kaon pt", precision = full_precision),
    kaon_eta     = Var("userFloat('kaon_eta')",        float, doc = "Kaon eta", precision = full_precision),
    kaon_phi     = Var("userFloat('kaon_phi')",        float, doc = "Kaon phi", precision = full_precision),
    kaon_charge  = Var("userInt('kaon_charge')",         int, doc = "Kaon charge"),
    kaon_mu_index = Var("userInt('kaon_mu_index')",      int, doc = "Matched muon index for track 1"),
    pion_pt      = Var("userFloat('pion_pt')",         float, doc = "Pion pt", precision = full_precision),
    pion_eta     = Var("userFloat('pion_eta')",        float, doc = "Pion eta", precision = full_precision),
    pion_phi     = Var("userFloat('pion_phi')",        float, doc = "Pion phi", precision = full_precision),
    pion_charge  = Var("userInt('pion_charge')",         int, doc = "Pion charge"),
    pion_mu_index = Var("userInt('pion_mu_index')",      int, doc = "Matched muon index for track 2"),
    kaon_sip     = Var("userFloat('kaon_sip')",        float, doc = "Kaon 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    pion_sip     = Var("userFloat('pion_sip')",        float, doc = "Pion 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability", precision = full_precision),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2", precision = medium_precision),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass", precision = full_precision),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt", precision = full_precision),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta", precision = full_precision),
    kin_phi      = Var("userFloat('kin_phi')",         float, doc = "Kinematic fit: vertex refitted phi", precision = full_precision),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error", precision = medium_precision),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot", precision = full_precision),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot", precision = full_precision),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS", precision = full_precision),
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS", precision = medium_precision),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV", precision = medium_precision),
    dstar_pion_pt   = Var("userFloat('dstar_pion_pt')",      float, doc = "DstarToD0Pi: slow pion pt", precision = full_precision),
    dstar_pion_eta  = Var("userFloat('dstar_pion_eta')",     float, doc = "DstarToD0Pi: slow pion eta", precision = full_precision),
    dstar_pion_phi  = Var("userFloat('dstar_pion_phi')",     float, doc = "DstarToD0Pi: slow pion phi", precision = full_precision),
    dstar_pion_charge = Var("userInt('dstar_pion_charge')",    int, doc = "DstarToD0Pi: slow pion charge"),
    dstar_mass      = Var("userFloat('dstar_mass')",         float, doc = "DstarToD0Pi: 3-body mass with vertex constraint", precision = full_precision),
    dstar_vtx_prob  = Var("userFloat('dstar_vtx_prob')",     float, doc = "DstarToD0Pi: vertex probability", precision = full_precision),
    dstar_vtx_chi2dof = Var("userFloat('dstar_vtx_chi2dof')", float, doc = "DstarToD0Pi: vertex normalized Chi^2", precision = medium_precision),
    dstar_pt        = Var("userFloat('dstar_pt')",           float, doc = "DstarToD0Pi: vertex refitted pt", precision = full_precision),
    dstar_eta       = Var("userFloat('dstar_eta')",          float, doc = "DstarToD0Pi: vertex refitted eta", precision = full_precision),
    dstar_phi       = Var("userFloat('dstar_phi')",          float, doc = "DstarToD0Pi: vertex refitted phi", precision = full_precision),
    dstar_massErr   = Var("userFloat('dstar_massErr')",      float, doc = "DstarToD0Pi: vertex refitted mass error", precision = medium_precision),
    dstar_lxy       = Var("userFloat('dstar_lxy')",          float, doc = "DstarToD0Pi: vertex displacement in XY plane wrt Beam Spot", precision = medium_precision),
    dstar_slxy      = Var("userFloat('dstar_sigLxy')",       float, doc = "DstarToD0Pi: vertex displacement significance in XY plane wrt Beam Spot", precision = medium_precision),
    dstar_cosAlphaXY = Var("userFloat('dstar_cosAlphaXY')",    float, doc = "DstarToD0Pi: cosine of pointing angle in XY wrt BS", precision = full_precision),
    dstar_sipBS     = Var("userFloat('dstar_sipBS')",        float, doc = "DstarToD0Pi: impact parameter significance of the candidate trajectory in XY wrt BS", precision = medium_precision),
    dstar_sipPV     = Var("userFloat('dstar_sipPV')",        float, doc = "DstarToD0Pi: impact parameter significance of the candidate trajectory in 3D wrt PV", precision = medium_precision),
)

D0ForMuonFakeVariablesMC = merge_psets(
    D0ForMuonFakeVariables,
    cms.PSet(
        gen_kaon_pdgId  = Var("userInt(  'gen_kaon_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_kaon_mpdgId = Var("userInt(  'gen_kaon_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_kaon_pt     = Var("userFloat('gen_kaon_pt')",     float,   doc = "Gen match: first track pt", precision = full_precision),
        gen_pion_pdgId  = Var("userInt(  'gen_pion_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_pion_mpdgId = Var("userInt(  'gen_pion_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_pion_pt     = Var("userFloat('gen_pion_pt')",     float,   doc = "Gen match: second track pt", precision = full_precision),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass", precision = full_precision),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt", precision = full_precision),
        dstar_pion_gen_pt = Var("userFloat('dstar_pion_gen_pt')", float, doc = "DstarToD0Pi: slow pion gen pt", precision = full_precision),
        dstar_gen_pdgId = Var("userInt('dstar_gen_pdgId')",     int,   doc = "DstarToD0Pi: gen pdg Id"),
        ),
)

D0ForMuonFakeTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFake","D0"),
    cut=cms.string(""),
    name=cms.string("d0"),
    doc=cms.string("D0s Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = D0ForMuonFakeVariables
)

D0ForMuonFakeMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFakeMC","D0"),
    cut=cms.string(""),
    name=cms.string("d0"),
    doc=cms.string("D0 Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = D0ForMuonFakeVariablesMC
)

# PhiToKK and DsToPhiPi

PhiForMuonFakeVariables = cms.PSet(
    mass         = Var("mass",                         float, doc = "Unfit invariant mass", precision = full_precision),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks", precision = medium_precision),
    trk1_pt      = Var("userFloat('trk1_pt')",         float, doc = "Track 1 pt", precision = full_precision),
    trk1_eta     = Var("userFloat('trk1_eta')",        float, doc = "Track 1 eta", precision = full_precision),
    trk1_phi     = Var("userFloat('trk1_phi')",        float, doc = "Track 1 phi", precision = full_precision),
    trk1_mu_index = Var("userInt('trk1_mu_index')",      int, doc = "Matched muon index for track 1"),
    trk2_pt      = Var("userFloat('trk2_pt')",         float, doc = "Track 2 pt", precision = full_precision),
    trk2_eta     = Var("userFloat('trk2_eta')",        float, doc = "Track 2 eta", precision = full_precision),
    trk2_phi     = Var("userFloat('trk2_phi')",        float, doc = "Track 2 phi", precision = full_precision),
    trk2_mu_index = Var("userInt('trk2_mu_index')",      int, doc = "Matched muon index for track 2"),
    trk1_sip     = Var("userFloat('trk1_sip')",        float, doc = "Track 1 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    trk2_sip     = Var("userFloat('trk2_sip')",        float, doc = "Track 2 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability", precision = full_precision),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2", precision = medium_precision),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass", precision = full_precision),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt", precision = full_precision),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta", precision = full_precision),
    kin_phi      = Var("userFloat('kin_phi')",         float, doc = "Kinematic fit: vertex refitted phi", precision = full_precision),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error", precision = medium_precision),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot", precision = full_precision),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot", precision = full_precision),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS", precision = full_precision),
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS", precision = medium_precision),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV", precision = medium_precision),
    ds_pion_pt   = Var("userFloat('ds_pion_pt')",      float, doc = "DsToPhiPi: pion pt", precision = full_precision),
    ds_pion_eta  = Var("userFloat('ds_pion_eta')",     float, doc = "DsToPhiPi: pion eta", precision = full_precision),
    ds_pion_phi  = Var("userFloat('ds_pion_phi')",     float, doc = "DsToPhiPi: pion phi", precision = full_precision),
    ds_pion_mu_index = Var("userInt('ds_pion_mu_index')",     float, doc = "DsToPhiPi: pion muon index", precision = medium_precision),
    ds_mass      = Var("userFloat('ds_mass')",         float, doc = "DsToPhiPi: 3-body mass with vertex constraint", precision = full_precision),
    ds_vtx_prob  = Var("userFloat('ds_vtx_prob')",     float, doc = "DsToPhiPi: vertex probability", precision = full_precision),
    ds_vtx_chi2dof = Var("userFloat('ds_vtx_chi2dof')", float, doc = "DsToPhiPi: vertex normalized Chi^2", precision = medium_precision),
    ds_pt        = Var("userFloat('ds_pt')",           float, doc = "DsToPhiPi: vertex refitted pt", precision = full_precision),
    ds_eta       = Var("userFloat('ds_eta')",          float, doc = "DsToPhiPi: vertex refitted eta", precision = full_precision),
    ds_phi       = Var("userFloat('ds_phi')",          float, doc = "DsToPhiPi: vertex refitted phi", precision = full_precision),
    ds_massErr   = Var("userFloat('ds_massErr')",      float, doc = "DsToPhiPi: vertex refitted mass error", precision = medium_precision),
    ds_lxy       = Var("userFloat('ds_lxy')",          float, doc = "DsToPhiPi: vertex displacement in XY plane wrt Beam Spot", precision = medium_precision),
    ds_slxy      = Var("userFloat('ds_sigLxy')",       float, doc = "DsToPhiPi: vertex displacement significance in XY plane wrt Beam Spot", precision = medium_precision),
    ds_cosAlphaXY = Var("userFloat('ds_cosAlphaXY')",    float, doc = "DsToPhiPi: cosine of pointing angle in XY wrt BS", precision = full_precision),
    ds_sipBS     = Var("userFloat('ds_sipBS')",        float, doc = "DsToPhiPi: impact parameter significance of the candidate trajectory in XY wrt BS", precision = medium_precision),
    ds_sipPV     = Var("userFloat('ds_sipPV')",        float, doc = "DsToPhiPi: impact parameter significance of the candidate trajectory in 3D wrt PV", precision = medium_precision),
)

PhiForMuonFakeVariablesMC = merge_psets(
    PhiForMuonFakeVariables,
    cms.PSet(
        gen_trk1_pdgId  = Var("userInt(  'gen_trk1_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_trk1_mpdgId = Var("userInt(  'gen_trk1_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_trk1_pt     = Var("userFloat('gen_trk1_pt')",     float,   doc = "Gen match: first track pt", precision = full_precision),
        gen_trk2_pdgId  = Var("userInt(  'gen_trk2_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_trk2_mpdgId = Var("userInt(  'gen_trk2_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_trk2_pt     = Var("userFloat('gen_trk2_pt')",     float,   doc = "Gen match: second track pt", precision = full_precision),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass", precision = full_precision),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt", precision = full_precision),
        ),
)

PhiForMuonFakeTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFake","Phi"),
    cut=cms.string(""),
    name=cms.string("phi"),
    doc=cms.string("Phi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = PhiForMuonFakeVariables
)

PhiForMuonFakeMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFakeMC","Phi"),
    cut=cms.string(""),
    name=cms.string("phi"),
    doc=cms.string("Phi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = PhiForMuonFakeVariablesMC
)

# LambdaToPPi

LambdaForMuonFakeVariables = cms.PSet(
    mass         = Var("mass",                         float, doc = "Unfit invariant mass", precision = full_precision),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks", precision = medium_precision),
    proton_pt      = Var("userFloat('proton_pt')",         float, doc = "Proton pt", precision = full_precision),
    proton_eta     = Var("userFloat('proton_eta')",        float, doc = "Proton eta", precision = full_precision),
    proton_phi     = Var("userFloat('proton_phi')",        float, doc = "Proton phi", precision = full_precision),
    proton_mu_index = Var("userInt('proton_mu_index')",      int, doc = "Matched muon index for track 1"),
    pion_pt      = Var("userFloat('pion_pt')",         float, doc = "Pion pt", precision = full_precision),
    pion_eta     = Var("userFloat('pion_eta')",        float, doc = "Pion eta", precision = full_precision),
    pion_phi     = Var("userFloat('pion_phi')",        float, doc = "Pion phi", precision = full_precision),
    pion_mu_index = Var("userInt('pion_mu_index')",      int, doc = "Matched muon index for track 2"),
    proton_sip     = Var("userFloat('proton_sip')",        float, doc = "Proton 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    pion_sip     = Var("userFloat('pion_sip')",        float, doc = "Pion 2D impact parameter significance wrt Beam Spot", precision = medium_precision),
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability", precision = full_precision),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2", precision = medium_precision),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass", precision = full_precision),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt", precision = full_precision),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta", precision = full_precision),
    kin_phi      = Var("userFloat('kin_phi')",         float, doc = "Kinematic fit: vertex refitted phi", precision = full_precision),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error", precision = medium_precision),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot", precision = full_precision),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot", precision = full_precision),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS", precision = full_precision),
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS", precision = medium_precision),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV", precision = medium_precision),
)

LambdaForMuonFakeVariablesMC = merge_psets(
    LambdaForMuonFakeVariables,
    cms.PSet(
        gen_proton_pdgId  = Var("userInt(  'gen_proton_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_proton_mpdgId = Var("userInt(  'gen_proton_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_proton_pt     = Var("userFloat('gen_proton_pt')",     float,   doc = "Gen match: first track pt", precision = full_precision),
        gen_pion_pdgId  = Var("userInt(  'gen_pion_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_pion_mpdgId = Var("userInt(  'gen_pion_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_pion_pt     = Var("userFloat('gen_pion_pt')",     float,   doc = "Gen match: second track pt", precision = full_precision),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass", precision = full_precision),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt", precision = full_precision),
        ),
)

LambdaForMuonFakeTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFake","Lambda"),
    cut=cms.string(""),
    name=cms.string("lambda"),
    doc=cms.string("Lambdas Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = LambdaForMuonFakeVariables
)

LambdaForMuonFakeMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("V0ForMuonFakeMC","Lambda"),
    cut=cms.string(""),
    name=cms.string("lambda"),
    doc=cms.string("Lambda Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = LambdaForMuonFakeVariablesMC
)


V0ForMuonFakeSequence   = cms.Sequence(V0ForMuonFake)
V0ForMuonFakeMcSequence = cms.Sequence(V0ForMuonFakeMC)
V0ForMuonFakeTables     = cms.Sequence(KsForMuonFakeTable + D0ForMuonFakeTable + PhiForMuonFakeTable + LambdaForMuonFakeTable)
V0ForMuonFakeMcTables   = cms.Sequence(KsForMuonFakeMcTable + D0ForMuonFakeMcTable + PhiForMuonFakeMcTable + LambdaForMuonFakeMcTable)
