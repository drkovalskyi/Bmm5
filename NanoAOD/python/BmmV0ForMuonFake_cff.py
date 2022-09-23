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
    "BmmV0Producer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
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
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS"),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV"),
)

KsForMuonFakeVariablesMC = merge_psets(
    KsForMuonFakeVariables,
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
    mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks"),
    kaon_pt      = Var("userFloat('kaon_pt')",         float, doc = "Kaon pt"),
    kaon_eta     = Var("userFloat('kaon_eta')",        float, doc = "Kaon eta"),
    kaon_phi     = Var("userFloat('kaon_phi')",        float, doc = "Kaon phi"),
    kaon_charge  = Var("userInt('kaon_charge')",         int, doc = "Kaon charge"),
    kaon_mu_index = Var("userInt('kaon_mu_index')",      int, doc = "Matched muon index for track 1"),
    pion_pt      = Var("userFloat('pion_pt')",         float, doc = "Pion pt"),
    pion_eta     = Var("userFloat('pion_eta')",        float, doc = "Pion eta"),
    pion_phi     = Var("userFloat('pion_phi')",        float, doc = "Pion phi"),
    pion_charge  = Var("userInt('pion_charge')",         int, doc = "Pion charge"),
    pion_mu_index = Var("userInt('pion_mu_index')",      int, doc = "Matched muon index for track 2"),
    kaon_sip     = Var("userFloat('kaon_sip')",        float, doc = "Kaon 2D impact parameter significance wrt Beam Spot"),
    pion_sip     = Var("userFloat('pion_sip')",        float, doc = "Pion 2D impact parameter significance wrt Beam Spot"),
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
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS"),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV"),
    dstar_pion_pt   = Var("userFloat('dstar_pion_pt')",      float, doc = "DstarToD0Pi: slow pion pt"),
    dstar_pion_eta  = Var("userFloat('dstar_pion_eta')",     float, doc = "DstarToD0Pi: slow pion eta"),
    dstar_pion_phi  = Var("userFloat('dstar_pion_phi')",     float, doc = "DstarToD0Pi: slow pion phi"),
    dstar_pion_charge = Var("userInt('dstar_pion_charge')",    int, doc = "DstarToD0Pi: slow pion charge"),
    dstar_mass      = Var("userFloat('dstar_mass')",         float, doc = "DstarToD0Pi: 3-body mass with vertex constraint"),
    dstar_vtx_prob  = Var("userFloat('dstar_vtx_prob')",     float, doc = "DstarToD0Pi: vertex probability"),
    dstar_vtx_chi2dof = Var("userFloat('dstar_vtx_chi2dof')", float, doc = "DstarToD0Pi: vertex normalized Chi^2"),
    dstar_pt        = Var("userFloat('dstar_pt')",           float, doc = "DstarToD0Pi: vertex refitted pt"),
    dstar_eta       = Var("userFloat('dstar_eta')",          float, doc = "DstarToD0Pi: vertex refitted eta"),
    dstar_phi       = Var("userFloat('dstar_phi')",          float, doc = "DstarToD0Pi: vertex refitted phi"),
    dstar_massErr   = Var("userFloat('dstar_massErr')",      float, doc = "DstarToD0Pi: vertex refitted mass error"),
    dstar_lxy       = Var("userFloat('dstar_lxy')",          float, doc = "DstarToD0Pi: vertex displacement in XY plane wrt Beam Spot"),
    dstar_slxy      = Var("userFloat('dstar_sigLxy')",       float, doc = "DstarToD0Pi: vertex displacement significance in XY plane wrt Beam Spot"),
    dstar_cosAlphaXY = Var("userFloat('dstar_cosAlphaXY')",    float, doc = "DstarToD0Pi: cosine of pointing angle in XY wrt BS"),
    dstar_sipBS     = Var("userFloat('dstar_sipBS')",        float, doc = "DstarToD0Pi: impact parameter significance of the candidate trajectory in XY wrt BS"),
    dstar_sipPV     = Var("userFloat('dstar_sipPV')",        float, doc = "DstarToD0Pi: impact parameter significance of the candidate trajectory in 3D wrt PV"),
)

D0ForMuonFakeVariablesMC = merge_psets(
    D0ForMuonFakeVariables,
    cms.PSet(
        gen_kaon_pdgId  = Var("userInt(  'gen_kaon_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_kaon_mpdgId = Var("userInt(  'gen_kaon_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_kaon_pt     = Var("userFloat('gen_kaon_pt')",     float,   doc = "Gen match: first track pt"),
        gen_pion_pdgId  = Var("userInt(  'gen_pion_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_pion_mpdgId = Var("userInt(  'gen_pion_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_pion_pt     = Var("userFloat('gen_pion_pt')",     float,   doc = "Gen match: second track pt"),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt"),
        dstar_pion_gen_pt = Var("userFloat('dstar_pion_gen_pt')", float, doc = "DstarToD0Pi: slow pion gen pt"),
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
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS"),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV"),
    ds_pion_pt   = Var("userFloat('ds_pion_pt')",      float, doc = "DsToPhiPi: pion pt"),
    ds_pion_eta  = Var("userFloat('ds_pion_eta')",     float, doc = "DsToPhiPi: pion eta"),
    ds_pion_phi  = Var("userFloat('ds_pion_phi')",     float, doc = "DsToPhiPi: pion phi"),
    ds_pion_mu_index = Var("userInt('ds_pion_mu_index')",     float, doc = "DsToPhiPi: pion muon index"),
    ds_mass      = Var("userFloat('ds_mass')",         float, doc = "DsToPhiPi: 3-body mass with vertex constraint"),
    ds_vtx_prob  = Var("userFloat('ds_vtx_prob')",     float, doc = "DsToPhiPi: vertex probability"),
    ds_vtx_chi2dof = Var("userFloat('ds_vtx_chi2dof')", float, doc = "DsToPhiPi: vertex normalized Chi^2"),
    ds_pt        = Var("userFloat('ds_pt')",           float, doc = "DsToPhiPi: vertex refitted pt"),
    ds_eta       = Var("userFloat('ds_eta')",          float, doc = "DsToPhiPi: vertex refitted eta"),
    ds_phi       = Var("userFloat('ds_phi')",          float, doc = "DsToPhiPi: vertex refitted phi"),
    ds_massErr   = Var("userFloat('ds_massErr')",      float, doc = "DsToPhiPi: vertex refitted mass error"),
    ds_lxy       = Var("userFloat('ds_lxy')",          float, doc = "DsToPhiPi: vertex displacement in XY plane wrt Beam Spot"),
    ds_slxy      = Var("userFloat('ds_sigLxy')",       float, doc = "DsToPhiPi: vertex displacement significance in XY plane wrt Beam Spot"),
    ds_cosAlphaXY = Var("userFloat('ds_cosAlphaXY')",    float, doc = "DsToPhiPi: cosine of pointing angle in XY wrt BS"),
    ds_sipBS     = Var("userFloat('ds_sipBS')",        float, doc = "DsToPhiPi: impact parameter significance of the candidate trajectory in XY wrt BS"),
    ds_sipPV     = Var("userFloat('ds_sipPV')",        float, doc = "DsToPhiPi: impact parameter significance of the candidate trajectory in 3D wrt PV"),
)

PhiForMuonFakeVariablesMC = merge_psets(
    PhiForMuonFakeVariables,
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
    mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks"),
    proton_pt      = Var("userFloat('proton_pt')",         float, doc = "Proton pt"),
    proton_eta     = Var("userFloat('proton_eta')",        float, doc = "Proton eta"),
    proton_phi     = Var("userFloat('proton_phi')",        float, doc = "Proton phi"),
    proton_mu_index = Var("userInt('proton_mu_index')",      int, doc = "Matched muon index for track 1"),
    pion_pt      = Var("userFloat('pion_pt')",         float, doc = "Pion pt"),
    pion_eta     = Var("userFloat('pion_eta')",        float, doc = "Pion eta"),
    pion_phi     = Var("userFloat('pion_phi')",        float, doc = "Pion phi"),
    pion_mu_index = Var("userInt('pion_mu_index')",      int, doc = "Matched muon index for track 2"),
    proton_sip     = Var("userFloat('proton_sip')",        float, doc = "Proton 2D impact parameter significance wrt Beam Spot"),
    pion_sip     = Var("userFloat('pion_sip')",        float, doc = "Pion 2D impact parameter significance wrt Beam Spot"),
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
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS"),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV"),
)

LambdaForMuonFakeVariablesMC = merge_psets(
    LambdaForMuonFakeVariables,
    cms.PSet(
        gen_proton_pdgId  = Var("userInt(  'gen_proton_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_proton_mpdgId = Var("userInt(  'gen_proton_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_proton_pt     = Var("userFloat('gen_proton_pt')",     float,   doc = "Gen match: first track pt"),
        gen_pion_pdgId  = Var("userInt(  'gen_pion_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_pion_mpdgId = Var("userInt(  'gen_pion_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_pion_pt     = Var("userFloat('gen_pion_pt')",     float,   doc = "Gen match: second track pt"),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt"),
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
