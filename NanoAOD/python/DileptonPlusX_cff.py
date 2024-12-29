from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms
import re
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

# NOTE: 
#    All instances of FlatTableProducers must end with Table in their
#    names so that their product match keep patterns in the default
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

def copy_pset(pset, replace_dict=dict(), drop_patterns=[]):
    result = cms.PSet()
    for name in pset.parameters_().keys():
        drop_param = False
        for pattern in drop_patterns:
            if re.search(pattern, name):
                drop_param = True
                break
        if drop_param:
            continue
        new_name = name
        for pattern,repl in replace_dict.items():
            new_name = re.sub(pattern,repl,new_name)
        if getattr(pset, name).pythonTypeName() in ('cms.PSet', 'cms.untracked.PSet'):
            value = getattr(pset, name).clone()
            for sname in value.parameters_().keys():
                svalue = getattr(value,sname)
                stype = svalue.pythonTypeName()
                if stype in ('cms.string', 'cms.untracked.string'):
                    new_svalue = svalue.value()
                    for pattern,repl in replace_dict.items():
                        new_svalue = re.sub(pattern,repl,new_svalue)
                    setattr(value,sname,new_svalue)
            setattr(result,new_name,value)
    return result

PrimaryVertexInfo = cms.EDProducer(
    "PrimaryVertexInformation",
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vertexScores = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandCollection = cms.InputTag("packedPFCandidates"),
    isMC = cms.bool(False),
    storeTracks = cms.bool(False),
)

PrimaryVertexInfoMc = PrimaryVertexInfo.clone( isMC = cms.bool(True) ) 

PrimaryVertexInfoTableVariables = cms.PSet(
    ntrks        = Var("userInt('ntrks')",         int, doc = "Number of tracks used in PV fit with tight association"),
    ndof         = Var("userFloat('ndof')",      float, doc = "number of degree of freedom for the vertex fit"),
    score        = Var("userFloat('score')",     float, doc = "Score of the primary vertex"),
    x            = Var("userFloat('x')",         float, doc = "Vertex position: x"),
    xErr         = Var("userFloat('xErr')",      float, doc = "Vertex position: xErr"),
    y            = Var("userFloat('y')",         float, doc = "Vertex position: y"),
    yErr         = Var("userFloat('yErr')",      float, doc = "Vertex position: yErr"),
    z            = Var("userFloat('z')",         float, doc = "Vertex position: z"),
    zErr         = Var("userFloat('zErr')",      float, doc = "Vertex position: zErr"),
    sumpt        = Var("userFloat('sumpt')",     float, doc = "Sum pt of tracks used in the vertex fit"),
    sumpt2       = Var("userFloat('sumpt2')",    float, doc = "Sum pt^2 of tracks used in the vertex fit"),
)

PrimaryVertexInfoMcTableVariables = merge_psets(
    PrimaryVertexInfoTableVariables,
)

PrimaryVertexInfoTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("PrimaryVertexInfo","pvs"),
    cut=cms.string(""),
    name=cms.string("pvs"),
    doc=cms.string("PrimaryVertexInfo"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = PrimaryVertexInfoTableVariables
)

PrimaryVertexInfoMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("PrimaryVertexInfoMc","pvs"),
    cut=cms.string(""),
    name=cms.string("pvs"),
    doc=cms.string("PrimaryVertexInfo"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = PrimaryVertexInfoMcTableVariables
)

###########################################################################################

Dileptons = cms.EDProducer(
    "DileptonPlusXProducer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
    # linkedObjects represent subset of original collections that are
    # kept in NanoAOD. Use them in order to link with NanoAOD objects.
    muonCollection = cms.InputTag("linkedObjects","muons"),
    electronCollection = cms.InputTag("linkedObjects","electrons"),
    photonCollection = cms.InputTag("linkedObjects","photons"),
    conversionCollection = cms.InputTag("oniaPhotonCandidates","conversions"),
    PFCandCollection = cms.InputTag("packedPFCandidates"),
    prunedGenParticleCollection = cms.InputTag("prunedGenParticles"),
    nanoGenParticleCollection = cms.InputTag("finalGenParticles"),
    MuonMinPt = cms.double(1.),
    MuonMaxEta = cms.double(2.4),
    KaonMinPt = cms.double(1.0),
    # KaonMinPt=cms.double(-1),
    KaonMaxEta = cms.double(2.4),
    ElectronMinPt = cms.double(2.),
    ElectronMaxEta = cms.double(2.4),
    KaonMinDCASig=cms.double(-1.),
    DiLeptonChargeCheck=cms.bool(False),
    minBKllMass = cms.double(4.0),
    maxBKllMass = cms.double(6.0),
    minLLGammaMass = cms.double(3.0),
    maxLLGammaMass = cms.double(6.0),
    minGammaPt = cms.double(1.0),
    minBKKllMass = cms.double(4.5),
    maxBKKllMass = cms.double(6.0),
    maxTwoTrackDOCA = cms.double(0.1),
    bdtEvent0 = cms.FileInPath('Bmm5/NanoAOD/data/TMVA-100-Events0_BDT.weights.xml'),
    bdtEvent1 = cms.FileInPath('Bmm5/NanoAOD/data/TMVA-100-Events1_BDT.weights.xml'),
    bdtEvent2 = cms.FileInPath('Bmm5/NanoAOD/data/TMVA-100-Events2_BDT.weights.xml'),
    xgbEvent0 = cms.FileInPath('Bmm5/NanoAOD/data/Run2017-2018-20200515-1144-Event0.model'),
    xgbEvent1 = cms.FileInPath('Bmm5/NanoAOD/data/Run2017-2018-20200515-1143-Event1.model'),
    xgbEvent2 = cms.FileInPath('Bmm5/NanoAOD/data/Run2017-2018-20200515-1143-Event2.model'),
    isMC = cms.bool(False),
    # injectMatchedBtohh = cms.bool(True),
    injectMatchedBtohh = cms.bool(False),
    injectBtohh = cms.bool(False),
    injectJpsiTracks = cms.bool(False),
    recoMuMuGamma = cms.bool(True),
    recoMuMuGammaConv = cms.bool(True),
    recoElElX = cms.bool(False),
    recoElMu = cms.bool(True),
    recoDstar = cms.bool(True),
    recoD0pipi = cms.bool(True),
    recoD0Kpi = cms.bool(True),
    recoKspipi = cms.bool(True),
    recoKstar = cms.bool(True),
    minBhhHadronPt = cms.double(4.0),
    maxBhhHadronEta = cms.double(1.4),
    minDhhHadronPt = cms.double(3.0),
    maxDhhHadronEta = cms.double(2.4),
    minKsHadronPt = cms.double(3.0),
    maxKsHadronEta = cms.double(2.4),
    minJpsiHadronPt = cms.double(2.0),
    minBhhMass   = cms.double(4.9),
    maxBhhMass   = cms.double(5.9),
    minBhhSigLxy = cms.double(4.0),
    minBhhVtxProb  = cms.double(0.01),
    recoMuMuPi = cms.bool(True),
    minD0Mass = cms.double(1.75),
    maxD0Mass = cms.double(1.95),
    minDmmMass = cms.double(1.25),
    maxDmmMass = cms.double(2.45),
    minKsMass = cms.double(0.45),
    maxKsMass = cms.double(0.55),
    minKstarMass  = cms.double(0.7),
    maxKstarMass  = cms.double(1.1),
    minDm = cms.double(0.135),
    maxDm = cms.double(0.160),
)

DileptonsMc = Dileptons.clone( isMC = cms.bool(True) ) 

##################################################################################
###
###                              Common PSets
###
##################################################################################

kinematic_common_pset = cms.PSet(
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
    kin_alphaBS  = Var("userFloat('kin_alphaBS')",     float, doc = "Kinematic fit: pointing angle in XY wrt BS"),
    kin_alphaBSErr = Var("userFloat('kin_alphaBSErr')", float, doc = "Kinematic fit: pointing angle in XY wrt BS"),
    kin_vtx_x    = Var("userFloat('kin_vtx_x')",       float, doc = "Kinematic fit: vertex x"),
    kin_vtx_xErr = Var("userFloat('kin_vtx_xErr')",    float, doc = "Kinematic fit: vertex x-error"),
    kin_vtx_y    = Var("userFloat('kin_vtx_y')",       float, doc = "Kinematic fit: vertex y"),
    kin_vtx_yErr = Var("userFloat('kin_vtx_yErr')",    float, doc = "Kinematic fit: vertex y-error"),
    kin_vtx_z    = Var("userFloat('kin_vtx_z')",       float, doc = "Kinematic fit: vertex y"),
    kin_vtx_zErr = Var("userFloat('kin_vtx_zErr')",    float, doc = "Kinematic fit: vertex y-error"),
)

kinematic_displacement_pset = cms.PSet(
    kin_alpha    = Var("userFloat('kin_alpha')",       float, doc = "Kinematic fit: pointing angle in 3D wrt PV"),
    kin_alphaErr = Var("userFloat('kin_alphaErr')",    float, doc = "Kinematic fit: pointing angle uncertainty in 3D wrt PV"),
    kin_l3d      = Var("userFloat('kin_l3d')",         float, doc = "Kinematic fit: decay length wrt Primary Vertex in 3D"),
    kin_sl3d     = Var("userFloat('kin_sl3d')",        float, doc = "Kinematic fit: decay length significance wrt Primary Vertex in 3D"),
    kin_pv_z     = Var("userFloat('kin_pv_z')",        float, doc = "Kinematic fit: primary vertex z"),
    kin_pv_zErr  = Var("userFloat('kin_pv_zErr')",     float, doc = "Kinematic fit: primary vertex z-error"),
    kin_pvip     = Var("userFloat('kin_pvip')",        float, doc = "Kinematic fit: impact parameter wrt Primary Vertex in 3D"),
    kin_spvip    = Var("userFloat('kin_spvip')",       float, doc = "Kinematic fit: impact parameter significance wrt Primary Vertex in 3D"),
    kin_pvipErr  = Var("userFloat('kin_pvipErr')",     float, doc = "Kinematic fit: impact parameter uncertainty wrt Primary Vertex in 3D"),
    kin_pvlip    = Var("userFloat('kin_pvlip')",       float, doc = "Kinematic fit: longitudinal impact parameter wrt Primary Vertex"),
    kin_pvlipSig = Var("userFloat('kin_pvlipSig')",    float, doc = "Kinematic fit: longitudinal impact parameter significance wrt Primary Vertex"),
    kin_pvlipErr = Var("userFloat('kin_pvlipErr')",    float, doc = "Kinematic fit: longitudinal impact parameter uncertainty wrt Primary Vertex"),
    kin_pvIndex  = Var("userInt('kin_pvIndex')",         int, doc = "Kinematic fit: primary vertex index"),
    kin_tau      = Var("userFloat('kin_tau')",         float, doc = "Kinematic fit: decay time wrt Primary Vertex in 3D"),
    kin_taue     = Var("userFloat('kin_taue')",        float, doc = "Kinematic fit: decay time error wrt Primary Vertex in 3D"),
    kin_tauxy    = Var("userFloat('kin_tauxy')",       float, doc = "Kinematic fit: decay time wrt Primary Vertex in XY"),
    kin_tauxye   = Var("userFloat('kin_tauxye')",      float, doc = "Kinematic fit: decay time error wrt Primary Vertex in XY"),
)

kinematic_pset = merge_psets(
    kinematic_common_pset,
    kinematic_displacement_pset,
    copy_pset(kinematic_displacement_pset, {"kin_":"kin_pv2_"}),
    copy_pset(kinematic_displacement_pset, {"kin_":"kin_refit_"}),
)

isolation_pset = cms.PSet(
    nTrks            = Var("userInt('nTrks')",             int,   doc = "Number of tracks compatible with the vertex by vertex probability"),
    nBMTrks          = Var("userInt('nBMTrks')",           int,   doc = "Number of tracks more compatible with the mm vertex than with PV by doca significance"),
    nDisTrks         = Var("userInt('nDisTrks')",          int,   doc = "Number of displaced tracks compatible with the vertex by vertex probability"),
    closetrk         = Var("userInt('closetrk')",          int,   doc = "Number of tracks compatible with the vertex by doca"),
    closetrks1       = Var("userInt('closetrks1')",        int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 1"),
    closetrks2       = Var("userInt('closetrks2')",        int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 2"),
    closetrks3       = Var("userInt('closetrks3')",        int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 3"),
    docatrk          = Var("userFloat('docatrk')",         float, doc = "Distance of closest approach of a track to the vertex"),
    m1iso            = Var("userFloat('m1iso')",           float, doc = "Muon isolation the way it's done in Bmm4"),
    m2iso            = Var("userFloat('m2iso')",           float, doc = "Muon isolation the way it's done in Bmm4"),
    iso              = Var("userFloat('iso')",             float, doc = "B isolation the way it's done in Bmm4"),
    otherVtxMaxProb  = Var("userFloat('otherVtxMaxProb')", float, doc = "Max vertexing probability of one of the muons with a random track with minPt=0.5GeV"),
    otherVtxMaxProb1 = Var("userFloat('otherVtxMaxProb1')", float, doc = "Max vertexing probability of one of the muons with a random track with minPt=1.0GeV"),
    otherVtxMaxProb2 = Var("userFloat('otherVtxMaxProb2')", float, doc = "Max vertexing probability of one of the muons with a random track with minPt=2.0GeV"),
)

def make_track_info_pset(prefix):
    return cms.PSet(
        **{
            f"{prefix}trkValidFrac":       Var(f"userFloat('{prefix}trkValidFrac')", float,
                                               doc = "Fraction of valid hits for inner track"),
            f"{prefix}trkNormChi2":        Var(f"userFloat('{prefix}trkNormChi2')", float,
                                               doc = "Normalized chi2 of the inner fit"),
            f"{prefix}pixelPattern":       Var(f"userInt('{prefix}pixelPattern')", int,
                                               doc = "Masks: barrel 0b1111, endcap 0b1110000"),
            f"{prefix}nPixels":            Var(f"userInt('{prefix}nPixels')", int,
                                               doc = "Number of valid pixel hits"),
            f"{prefix}nValidHits":         Var(f"userInt('{prefix}nValidHits')", int,
                                               doc = "Number of valid hits"),
            f"{prefix}nLostHitsInner":     Var(f"userInt('{prefix}nLostHitsInner')", int,
                                               doc = "Number of lost hits before tracker track"),
            f"{prefix}nLostHitsOn":        Var(f"userInt('{prefix}nLostHitsOn')", int,
                                               doc = "Number of lost hits on tracker track"),
            f"{prefix}nLostHitsOuter":     Var(f"userInt('{prefix}nLostHitsOuter')", int,
                                               doc = "Number of lost hits after tracker track"),
            f"{prefix}trkLayers":          Var(f"userInt('{prefix}trkLayers')", int,
                                               doc = "Number of layers with measurements in tracker track"),
            f"{prefix}trkLostLayersInner": Var(f"userInt('{prefix}trkLostLayersInner')", int,
                                               doc = "Number of lost layers befor tracker track"),
            f"{prefix}trkLostLayersOn":    Var(f"userInt('{prefix}trkLostLayersOn')", int,
                                               doc = "Number of lost layers on tracker track"),
            f"{prefix}trkLostLayersOuter": Var(f"userInt('{prefix}trkLostLayersOuter')", int,
                                               doc = "Number of lost layers after tracker track"),
            }
    )

##################################################################################
###
###                              Dilepton Info
###
##################################################################################

DileptonsDiMuonTableVariables = merge_psets(
    cms.PSet(
        mu1_index    = Var("userInt('mu1_index')",         int,   doc = "Index of corresponding leading muon"),
        mu1_pdgId    = Var("userInt('mu1_pdgId')",         int,   doc = "Leading muon candidate pdgId. Used in hadron fake studies"),
        mu1_pt       = Var("userFloat('mu1_pt')",          float, doc = "Leading muon pt"),
        mu1_eta      = Var("userFloat('mu1_eta')",         float, doc = "Leading muon eta"),
        mu1_phi      = Var("userFloat('mu1_phi')",         float, doc = "Leading muon phi"),
        mu2_index    = Var("userInt('mu2_index')",         int,   doc = "Index of corresponding subleading muon"),
        mu2_pdgId    = Var("userInt('mu2_pdgId')",         int,   doc = "Trailing muon candidate pdgId. Used in hadron fake studies"),
        mu2_pt       = Var("userFloat('mu2_pt')",          float, doc = "Trailing muon pt"),
        mu2_eta      = Var("userFloat('mu2_eta')",         float, doc = "Trailing muon eta"),
        mu2_phi      = Var("userFloat('mu2_phi')",         float, doc = "Trailing muon phi"),
        mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
        doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of muons"),
        bdt          = Var("userFloat('bdt')",             float, doc = "Bmm4 BDT"),
        mva          = Var("userFloat('mva')",             float, doc = "XGBoost based Bmm5 MVA"),
        # Kinematic Fit daugter info
        kin_mu1_pt    = Var("userFloat('kin_mu1pt')",       float, doc = "Kinematic fit: refitted muon 1 pt"),
        kin_mu1_eta   = Var("userFloat('kin_mu1eta')",      float, doc = "Kinematic fit: refitted muon 1 eta"),
        kin_mu1_phi   = Var("userFloat('kin_mu1phi')",      float, doc = "Kinematic fit: refitted muon 1 phi"),
        kin_mu2_pt    = Var("userFloat('kin_mu2pt')",       float, doc = "Kinematic fit: refitted muon 2 pt"),
        kin_mu2_eta   = Var("userFloat('kin_mu2eta')",      float, doc = "Kinematic fit: refitted muon 2 eta"),
        kin_mu2_phi   = Var("userFloat('kin_mu2phi')",      float, doc = "Kinematic fit: refitted muon 2 phi"),
        ),
    kinematic_pset,
    copy_pset(kinematic_common_pset, {"kin_":"kinpc_"}),
    isolation_pset,
    make_track_info_pset("mu1_"),
    make_track_info_pset("mu2_"),
)

DileptonsDiMuonMcTableVariables = merge_psets(
    DileptonsDiMuonTableVariables,
    cms.PSet(
        gen_mu1_pdgId      = Var("userInt(  'gen_mu1_pdgId')",        int,   doc = "Gen match: first daughter pdg Id"),
        gen_mu1_index      = Var("userInt(  'gen_mu1_index')",        int,   doc = "Gen match: first daughter index int GenPart"),
        gen_mu1_mpdgId     = Var("userInt(  'gen_mu1_mpdgId')",       int,   doc = "Gen match: first daughter mother pdg Id"),
        gen_mu1_pt         = Var("userFloat('gen_mu1_pt')",         float,   doc = "Gen match: first daughter pt"),
        gen_mu2_pdgId      = Var("userInt(  'gen_mu2_pdgId')",        int,   doc = "Gen match: second daughter pdg Id"),
        gen_mu2_index      = Var("userInt(  'gen_mu2_index')",        int,   doc = "Gen match: second daughter index in GenPart"),
        gen_mu2_mpdgId     = Var("userInt(  'gen_mu2_mpdgId')",       int,   doc = "Gen match: second daughter mother pdg Id"),
        gen_mu2_pt         = Var("userFloat('gen_mu2_pt')",         float,   doc = "Gen match: second daughter pt"),
        gen_pdgId         = Var("userInt(  'gen_pdgId')",           int,   doc = "Gen match: dilepton pdg Id"),
        gen_index         = Var("userInt(  'gen_index')",           int,   doc = "Gen match: dilepton index in GenPart"),
        gen_cpdgId        = Var("userInt(  'gen_cpdgId')",          int,   doc = "Gen match: common mother pdg Id"),
        gen_cindex        = Var("userInt(  'gen_cindex')",          int,   doc = "Gen match: common mother index in GenPart"),
        gen_mpdgId        = Var("userInt(  'gen_mpdgId')",          int,   doc = "Gen match: dilepton mother pdg Id"),
        gen_mass          = Var("userFloat('gen_mass')",          float,   doc = "Gen match: dilepton mass"),
        gen_pt            = Var("userFloat('gen_pt')",            float,   doc = "Gen match: dilepton pt"),
        gen_prod_x        = Var("userFloat('gen_prod_x')",        float,   doc = "Gen match: dilepton mother production vertex x"),
        gen_prod_y        = Var("userFloat('gen_prod_y')",        float,   doc = "Gen match: dilepton mother production vertex y"),
        gen_prod_z        = Var("userFloat('gen_prod_z')",        float,   doc = "Gen match: dilepton mother production vertex z"),
        gen_vtx_x         = Var("userFloat('gen_vtx_x')",         float,   doc = "Gen match: dilepton vertex x"),
        gen_vtx_y         = Var("userFloat('gen_vtx_y')",         float,   doc = "Gen match: dilepton vertex y"),
        gen_vtx_z         = Var("userFloat('gen_vtx_z')",         float,   doc = "Gen match: dilepton vertex z"),
        gen_l3d           = Var("userFloat('gen_l3d')",           float,   doc = "Gen match: dilepton decay legnth 3D"),
        gen_lxy           = Var("userFloat('gen_lxy')",           float,   doc = "Gen match: dilepton decay legnth XY"),
        gen_tau           = Var("userFloat('gen_tau')",           float,   doc = "Gen match: dilepton decay time 3D"),
        gen_doca          = Var("userFloat('gen_doca')",          float,   doc = "Gen match: dilepton distance of closest approach"),
        gen_alpha_p_phi   = Var("userFloat('gen_alpha_p_phi')",   float,   doc = "Difference in direction between reconstructed and generated B in phi"),
        gen_alpha_p_theta = Var("userFloat('gen_alpha_p_theta')", float,   doc = "Difference in direction between reconstructed and generated B in theta"),
        gen_alpha_ip      = Var("userFloat('gen_alpha_ip')",      float,   doc = "Pointing angle due to reco IP for gen vertex and momentum"),
        gen_alpha_vtx     = Var("userFloat('gen_alpha_vtx')",     float,   doc = "Pointing angle due to reco vtx for gen IP and momentum"),
    ),
)

DileptonsDiMuonTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","MuMu"),
    cut=cms.string(""),
    name=cms.string("mm"),
    doc=cms.string("Dimuon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsDiMuonTableVariables
)

DileptonsDiMuonMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","MuMu"),
    cut=cms.string(""),
    name=cms.string("mm"),
    doc=cms.string("Dimuon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsDiMuonMcTableVariables
)

### Relevant Track Information

TrackTableVariables = cms.PSet(
    pt       = Var("userFloat('pt')",          float, doc = "pt"),
    eta      = Var("userFloat('eta')",         float, doc = "eta"),
    phi      = Var("userFloat('phi')",         float, doc = "phi"),
    charge   = Var("userInt('charge')",          int, doc = "charge"),
    vtx_x    = Var("userFloat('vtx_x')",       float, doc = "reference point"),
    vtx_y    = Var("userFloat('vtx_y')",       float, doc = "reference point"),
    vtx_z    = Var("userFloat('vtx_z')",       float, doc = "reference point"),
    sip      = Var("userFloat('bs_sip')",      float, doc = "Track impact parameter wrt BS: significance"),
    trkValidFrac   = Var("userFloat('trkValidFrac')",   float, doc = "N_valid/(N_valid + N_lost + N_lostIn + N_lostOut)"),
    trkNormChi2    = Var("userFloat('trkNormChi2')",    float, doc = "Normalized chi2"),
    pixelPattern   = Var("userInt('pixelPattern')",       int, doc = "Masks: barrel 0b1111, endcap 0b1110000"),
    nPixels        = Var("userInt('nPixels')",            int, doc = "Number Of Valid Pixel Hits"),
    nValidHits     = Var("userInt('nValidHits')",         int, doc = "Number Of Valid Tracker Hits"),
    trkLayers           = Var("userInt('trkLayers')",          int, doc = "Tracker Layers With Measurement"),
    trkLostLayersInner  = Var("userInt('trkLostLayersInner')", int, doc = "Missing inner layers"),
    trkLostLayersOn     = Var("userInt('trkLostLayersOn')",    int, doc = "Missing layers on the track"),
    trkLostLayersOuter  = Var("userInt('trkLostLayersOuter')", int, doc = "Missing outter layers"),
)

TrackMcTableVariables = merge_psets(
    TrackTableVariables,
    cms.PSet(
        gen_pdgId       = Var("userInt('gen_pdgId')",          int, doc = "Gen match pdgId"),
        gen_mpdgId      = Var("userInt('gen_mpdgId')",         int, doc = "Gen match mother pdgId"),
    ),
)

DileptonTrackTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","InterestingTracks"),
    cut=cms.string(""),
    name=cms.string("trk"),
    doc=cms.string("Track Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = TrackTableVariables
)
DileptonTrackMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","InterestingTracks"),
    cut=cms.string(""),
    name=cms.string("trk"),
    doc=cms.string("Track Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = TrackMcTableVariables
)

IsoTableVariables = cms.PSet(
    trk_index       = Var("userInt('trk_index')",           int, doc = "Track index"),
    trk_pv_index    = Var("userInt('trk_pv_index')",        int, doc = "Track PV index"),
    ref_pv_index    = Var("userInt('ref_pv_index')",        int, doc = "Reference candidate PV index"),
    mm_index        = Var("userInt('mm_index')",            int, doc = "mm index"),
    hh_index        = Var("userInt('hh_index')",            int, doc = "hh index"),
    dr              = Var("userFloat('dr')",              float, doc = "deltaR(track,dimuon)"),
    pv_ip3d_status  = Var("userInt('pv_ip3d_status')",      int, doc = "Track impact parameter wrt PV: status"),
    pv_ip3d         = Var("userFloat('pv_ip3d')",         float, doc = "Track impact parameter wrt PV: value"),
    pv_sip3d        = Var("userFloat('pv_sip3d')",        float, doc = "Track impact parameter wrt PV: significance"),
    d1_doca         = Var("userFloat('d1_doca')",         float, doc = "Track distance of closest approach wrt daughter 1"),
    d1_vtx_prob     = Var("userFloat('d1_vtx_prob')",     float, doc = "Track + daughter1 vertex: probability"),
    d1_vtx_x        = Var("userFloat('d1_vtx_x')",        float, doc = "Track + daughter1 vertex: x"),
    d1_vtx_y        = Var("userFloat('d1_vtx_y')",        float, doc = "Track + daughter1 vertex: y"),
    d1_vtx_z        = Var("userFloat('d1_vtx_z')",        float, doc = "Track + daughter1 vertex: z"),
    d1_vtx_xErr     = Var("userFloat('d1_vtx_xErr')",     float, doc = "Track + daughter1 vertex: xErr"),
    d1_vtx_yErr     = Var("userFloat('d1_vtx_yErr')",     float, doc = "Track + daughter1 vertex: yErr"),
    d1_vtx_zErr     = Var("userFloat('d1_vtx_zErr')",     float, doc = "Track + daughter1 vertex: zErr"),
    d1_vtx_alphaBS  = Var("userFloat('d1_vtx_alphaBS')",  float, doc = "Track + daughter1 vertex: pointing angle 2D"),
    d1_vtx_alpha    = Var("userFloat('d1_vtx_alpha')",    float, doc = "Track + daughter1 vertex: pointing angle 3D"),
    d1_vtx_pv_l3d   = Var("userFloat('d1_vtx_pv_l3d')",   float, doc = "Track + daughter1 vertex: decay length 3D wrt PV"),
    d1_vtx_pv_sl3d  = Var("userFloat('d1_vtx_pv_sl3d')",  float, doc = "Track + daughter1 vertex: decay length significance 3D wrt PV"),
    d1_vtx_ref_l3d  = Var("userFloat('d1_vtx_ref_l3d')",  float, doc = "Track + daughter1 vertex: decay length 3D wrt ref vertex"),
    d1_vtx_ref_sl3d = Var("userFloat('d1_vtx_ref_sl3d')", float, doc = "Track + daughter1 vertex: decay length significance 3D wrt ref vertex"),
    d2_doca         = Var("userFloat('d2_doca')",         float, doc = "Track distance of closest approach wrt daughter 2"),
    d2_vtx_prob     = Var("userFloat('d2_vtx_prob')",     float, doc = "Track + daughter2 vertex: probability"),
    d2_vtx_x        = Var("userFloat('d2_vtx_x')",        float, doc = "Track + daughter2 vertex: x"),
    d2_vtx_y        = Var("userFloat('d2_vtx_y')",        float, doc = "Track + daughter2 vertex: y"),
    d2_vtx_z        = Var("userFloat('d2_vtx_z')",        float, doc = "Track + daughter2 vertex: z"),
    d2_vtx_xErr     = Var("userFloat('d2_vtx_xErr')",     float, doc = "Track + daughter2 vertex: xErr"),
    d2_vtx_yErr     = Var("userFloat('d2_vtx_yErr')",     float, doc = "Track + daughter2 vertex: yErr"),
    d2_vtx_zErr     = Var("userFloat('d2_vtx_zErr')",     float, doc = "Track + daughter2 vertex: zErr"),
    d2_vtx_alpha    = Var("userFloat('d2_vtx_alpha')",    float, doc = "Track + daughter2 vertex: pointing angle 3D"),
    d2_vtx_alphaBS  = Var("userFloat('d2_vtx_alphaBS')",  float, doc = "Track + daughter2 vertex: pointing angle 2D"),
    d2_vtx_pv_l3d   = Var("userFloat('d2_vtx_pv_l3d')",   float, doc = "Track + daughter2 vertex: decay length 3D wrt PV"),
    d2_vtx_pv_sl3d  = Var("userFloat('d2_vtx_pv_sl3d')",  float, doc = "Track + daughter2 vertex: decay length significance 3D wrt PV"),
    d2_vtx_ref_l3d  = Var("userFloat('d2_vtx_ref_l3d')",  float, doc = "Track + daughter2 vertex: decay length 3D wrt ref vertex"),
    d2_vtx_ref_sl3d = Var("userFloat('d2_vtx_ref_sl3d')", float, doc = "Track + daughter2 vertex: decay length significance 3D wrt ref vertex"),
    vtx_prob        = Var("userFloat('vtx_prob')",        float, doc = "Track + reference candidate vertex: probability"),
)

IsoMcTableVariables = merge_psets(
    IsoTableVariables,
    cms.PSet(),
)

DileptonIsoTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","Iso"),
    cut=cms.string(""),
    name=cms.string("iso"),
    doc=cms.string("Isolation information"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = IsoTableVariables
)
DileptonIsoMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","Iso"),
    cut=cms.string(""),
    name=cms.string("iso"),
    doc=cms.string("Isolation information"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = IsoMcTableVariables
)


### hh

DileptonsHHTableVariables = copy_pset(DileptonsDiMuonTableVariables, {"mu1_":"had1_", "mu2_":"had2_"})
DileptonsHHMcTableVariables = copy_pset(DileptonsDiMuonMcTableVariables, {"mu1_":"had1_", "mu2_":"had2_"})

DileptonsHHTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","HH"),
    cut=cms.string(""),
    name=cms.string("hh"),
    doc=cms.string("hh Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsHHTableVariables
)

DileptonsHHMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","HH"),
    cut=cms.string(""),
    name=cms.string("hh"),
    doc=cms.string("hh Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsHHMcTableVariables
)

DileptonsElElTableVariables = copy_pset(DileptonsDiMuonTableVariables, {"mu1_":"el1_", "mu2_":"el2_"})
DileptonsElElMcTableVariables = copy_pset(DileptonsDiMuonMcTableVariables, {"mu1_":"el1_", "mu2_":"el2_"})

DileptonsElElTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","ElEl"),
    cut=cms.string(""),
    name=cms.string("ee"),
    doc=cms.string("ee Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsElElTableVariables
)

DileptonsElElMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","ElEl"),
    cut=cms.string(""),
    name=cms.string("ee"),
    doc=cms.string("ee Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsElElMcTableVariables
)

def fix_parameter_names(pset):
    new_pset = cms.PSet()
    for p in pset.parameters_():
        new_name = p
        if "el1_" in p:
            new_name = p.replace("el1_", "el_")
        elif "mu2_" in p:
            new_name = p.replace("mu2_", "mu_")
        setattr(new_pset, new_name, getattr(pset, p))
    return new_pset

DileptonsElMuTableVariables = fix_parameter_names(copy_pset(DileptonsDiMuonTableVariables, {"mu1_":"el1_"}))
DileptonsElMuMcTableVariables = fix_parameter_names(copy_pset(DileptonsDiMuonMcTableVariables, {"mu1_":"el1_"}))
# print(DileptonsElMuMcTableVariables)


DileptonsElMuTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","ElMu"),
    cut=cms.string(""),
    name=cms.string("em"),
    doc=cms.string("em Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsElMuTableVariables
)

DileptonsElMuMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","ElMu"),
    cut=cms.string(""),
    name=cms.string("em"),
    doc=cms.string("em Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsElMuMcTableVariables
)


##################################################################################
###
###                              B to K l l
###
##################################################################################

DileptonsKmumuTableVariables =  merge_psets(
    copy_pset(kinematic_pset,{"kin_":"nomc_"}),
    copy_pset(kinematic_pset,{"kin_":"jpsimc_"}),
    cms.PSet(
        mm_index       = Var("userInt('mm_index')",          int,   doc = "Index of dimuon pair"),
        kaon_charge    = Var("userInt('kaon_charge')",       int,   doc = "Kaon charge"),
        kaon_pt        = Var("userFloat('kaon_pt')",         float, doc = "Kaon pt"),
        kaon_eta       = Var("userFloat('kaon_eta')",        float, doc = "Kaon eta"),
        kaon_phi       = Var("userFloat('kaon_phi')",        float, doc = "Kaon phi"),
        kaon_dxy_bs    = Var("userFloat('kaon_dxy_bs')",     float, doc = "Kaon impact parameter wrt the beam spot"),
        kaon_sdxy_bs   = Var("userFloat('kaon_sdxy_bs')",    float, doc = "Kaon impact parameter significance wrt the beam spot"),
        kaon_l1_doca   = Var("userFloat('kaon_l1_doca')",    float, doc = "Kaon distance of closest approach to lepton1"),
        kaon_l2_doca   = Var("userFloat('kaon_l2_doca')",    float, doc = "Kaon distance of closest approach to lepton2"),
        bmm_nTrks      = Var("userInt('bmm_nTrks')",         int,   doc = "Number of tracks compatible with the vertex by vertex probability (BtoJpsiK as Bmm)"),
        bmm_nBMTrks    = Var("userInt('bmm_nBMTrks')",       int,   doc = "Number of tracks more compatible with the mm vertex than with PV by doca significance (BtoJpsiK as Bmm)"),
        bmm_nDisTrks   = Var("userInt('bmm_nDisTrks')",      int,   doc = "Number of displaced tracks compatible with the vertex by vertex probability (BtoJpsiK as Bmm)"),
        bmm_closetrk   = Var("userInt('bmm_closetrk')",      int,   doc = "Number of tracks compatible with the vertex by doca (BtoJpsiK as Bmm)"),
        bmm_closetrks1 = Var("userInt('bmm_closetrks1')",    int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 1 (BtoJpsiK as Bmm)"),
        bmm_closetrks2 = Var("userInt('bmm_closetrks2')",    int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 2 (BtoJpsiK as Bmm)"),
        bmm_closetrks3 = Var("userInt('bmm_closetrks3')",    int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 3 (BtoJpsiK as Bmm)"),
        bmm_docatrk    = Var("userFloat('bmm_docatrk')",     float, doc = "Distance of closest approach of a track to the vertex (BtoJpsiK as Bmm)"),
        bmm_m1iso      = Var("userFloat('bmm_m1iso')",       float, doc = "Muon isolation the way it's done in Bmm4 (BtoJpsiK as Bmm)"),
        bmm_m2iso      = Var("userFloat('bmm_m2iso')",       float, doc = "Muon isolation the way it's done in Bmm4 (BtoJpsiK as Bmm)"),
        bmm_iso        = Var("userFloat('bmm_iso')",         float, doc = "B isolation the way it's done in Bmm4 (BtoJpsiK as Bmm)"),
        bmm_bdt        = Var("userFloat('bmm_bdt')",         float, doc = "BDT (BtoJpsiK as Bmm)"),
        bmm_otherVtxMaxProb = Var("userFloat('bmm_otherVtxMaxProb')", float, doc = "Max vertexing probability of one of the leptons with a random track with minPt=0.5GeV (BtoJpsiK as Bmm)"),
        bmm_otherVtxMaxProb1 = Var("userFloat('bmm_otherVtxMaxProb1')", float, doc = "Max vertexing probability of one of the leptons with a random track with minPt=1.0GeV (BtoJpsiK as Bmm)"),
        bmm_otherVtxMaxProb2 = Var("userFloat('bmm_otherVtxMaxProb2')", float, doc = "Max vertexing probability of one of the leptons with a random track with minPt=2.0GeV (BtoJpsiK as Bmm)"),
        bmm_mva        = Var("userFloat('bmm_mva')",         float, doc = "MVA (BtoJpsiK as Bmm)"),
        # Kinematic Fit daugter info
        nomc_kaon_pt    = Var("userFloat('nomc_kaon1_pt')",       float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 pt"),
        nomc_kaon_eta   = Var("userFloat('nomc_kaon1_eta')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 eta"),
        nomc_kaon_phi   = Var("userFloat('nomc_kaon1_phi')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 phi"),
        jpsimc_kaon_pt    = Var("userFloat('jpsimc_kaon1_pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 pt"),
        jpsimc_kaon_eta   = Var("userFloat('jpsimc_kaon1_eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 eta"),
        jpsimc_kaon_phi   = Var("userFloat('jpsimc_kaon1_phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 phi"),
    )
)

DileptonsKmumuMcTableVariables = merge_psets(
    DileptonsKmumuTableVariables,
    cms.PSet(
        gen_kaon_pdgId  = Var("userInt('gen_kaon_pdgId')",    int,   doc = "Gen match: kaon pdg Id"),
        gen_kaon_mpdgId = Var("userInt('gen_kaon_mpdgId')",   int,   doc = "Gen match: kaon mother pdg Id"),
        gen_kaon_pt     = Var("userFloat('gen_kaon_pt')",     float, doc = "Gen match: kaon pt"),
        gen_pdgId       = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: kmm pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float, doc = "Gen match: kmm mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float, doc = "Gen match: kmm pt"),
        gen_prod_x      = Var("userFloat('gen_prod_x')",      float, doc = "Gen match: kmm mother production vertex x"),
        gen_prod_y      = Var("userFloat('gen_prod_y')",      float, doc = "Gen match: kmm mother production vertex y"),
        gen_prod_z      = Var("userFloat('gen_prod_z')",      float, doc = "Gen match: kmm mother production vertex z"),
        gen_l3d         = Var("userFloat('gen_l3d')",         float, doc = "Gen match: kmm decay legnth 3D"),
        gen_lxy         = Var("userFloat('gen_lxy')",         float, doc = "Gen match: kmm decay legnth XY"),
        gen_tau         = Var("userFloat('gen_tau')",         float, doc = "Gen match: kmm decay time 3D"),
    )
)
        

DileptonsKmumuTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","BToKmumu"),
    cut=cms.string(""),
    name=cms.string("bkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKmumuTableVariables
)

DileptonsKmumuMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","BToKmumu"),
    cut=cms.string(""),
    name=cms.string("bkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKmumuMcTableVariables
)

DileptonsKeeTableVariables = copy_pset(DileptonsKmumuTableVariables, {"mm_index":"ee_index"})
DileptonsKeeMcTableVariables = copy_pset(DileptonsKmumuMcTableVariables, {"mm_index":"ee_index"})

DileptonsKeeTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","BToKee"),
    cut=cms.string(""),
    name=cms.string("bkee"),
    doc=cms.string("BToKee Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKeeTableVariables
)

DileptonsKeeMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","BToKee"),
    cut=cms.string(""),
    name=cms.string("bkee"),
    doc=cms.string("BToKee Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKeeMcTableVariables
)

##################################################################################
###
###                              Dstar to D0 pi
###                               '-> D0 to Mu Mu
###                               '-> D0 to pi pi
###                               '-> D0 to K pi
###
##################################################################################

DileptonsDstarTableVariables =  merge_psets(
    cms.PSet(
        mm_index       = Var("userInt('mm_index')",          int,   doc = "Index of muon pair"),
        hh_index       = Var("userInt('hh_index')",          int,   doc = "Index of hadron pair"),
        pion_charge    = Var("userInt('pion_charge')",       int,   doc = "Pion charge"),
        mass           = Var("userFloat('mass')",            float, doc = "Raw dstar mass"),
        dm_raw         = Var("userFloat('dm_raw')",          float, doc = "Raw dm"),
        dm_free        = Var("userFloat('dm_free')",         float, doc = "dm with vertexed d0 and raw soft pion"),
        dm_pv          = Var("userFloat('dm_prompt')",       float, doc = "dm with vertexed d0 and refitted with primary vertex constraint soft pion"),
        pion_pt        = Var("userFloat('pion_pt')",         float, doc = "Pion pt"),
        pion_eta       = Var("userFloat('pion_eta')",        float, doc = "Pion eta"),
        pion_phi       = Var("userFloat('pion_phi')",        float, doc = "Pion phi"),
        pion_dxy_bs    = Var("userFloat('pion_dxy_bs')",     float, doc = "Pion impact parameter wrt the beam spot"),
        pion_sdxy_bs   = Var("userFloat('pion_sdxy_bs')",    float, doc = "Pion impact parameter significance wrt the beam spot"),
        pv_prob        = Var("userFloat('pv_prob')",         float, doc = "PV refit probability"),
        pv_sum_pt      = Var("userFloat('pv_sum_pt')",       float, doc = "PV sum pt"),
        pv_sum_pt2     = Var("userFloat('pv_sum_pt2')",      float, doc = "PV sum pt^2"),
        pv_ntrks       = Var("userInt('pv_ntrks')",          int,   doc = "PV number of tracks"),
        pv_with_pion_prob = Var("userFloat('pv_with_pion_prob')", float, doc = "PV refit probability with soft pion"),
        # pion_l1_doca   = Var("userFloat('pion_l1_doca')",    float, doc = "Pion distance of closest approach to lepton1"),
        # pion_l2_doca   = Var("userFloat('pion_l2_doca')",    float, doc = "Pion distance of closest approach to lepton2"),
        # Kinematic Fit daugter info
        # nomc_kaon1pt    = Var("userFloat('nomc_kaon1pt')",       float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 pt"),
        # nomc_kaon1eta   = Var("userFloat('nomc_kaon1eta')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 eta"),
        # nomc_kaon1phi   = Var("userFloat('nomc_kaon1phi')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 phi"),
        # jpsimc_kaon1pt    = Var("userFloat('jpsimc_kaon1pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 pt"),
        # jpsimc_kaon1eta   = Var("userFloat('jpsimc_kaon1eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 eta"),
        # jpsimc_kaon1phi   = Var("userFloat('jpsimc_kaon1phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 phi"),
    )
)

DileptonsDstarMcTableVariables = merge_psets(
    DileptonsDstarTableVariables,
    cms.PSet(
        gen_pion_pdgId  = Var("userInt('gen_pion_pdgId')",    int,   doc = "Gen match: pion pdg Id"),
        gen_pion_mpdgId = Var("userInt('gen_pion_mpdgId')",   int,   doc = "Gen match: pion mother pdg Id"),
        gen_pion_pt     = Var("userFloat('gen_pion_pt')",     float, doc = "Gen match: pion pt"),
        gen_pdgId       = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: dstar pdg Id"),
        gen_mpdgId      = Var("userInt('gen_mpdgId')",        int,   doc = "Gen match: dstar mother pdg Id"),
        gen_cpdgId      = Var("userInt('gen_cpdgId')",        int,   doc = "Gen match: dstar commont ancestor pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float, doc = "Gen match: kmm mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float, doc = "Gen match: kmm pt"),
        gen_prod_x      = Var("userFloat('gen_prod_x')",      float, doc = "Gen match: kmm mother production vertex x"),
        gen_prod_y      = Var("userFloat('gen_prod_y')",      float, doc = "Gen match: kmm mother production vertex y"),
        gen_prod_z      = Var("userFloat('gen_prod_z')",      float, doc = "Gen match: kmm mother production vertex z"),
        gen_l3d         = Var("userFloat('gen_l3d')",         float, doc = "Gen match: kmm decay legnth 3D"),
        gen_lxy         = Var("userFloat('gen_lxy')",         float, doc = "Gen match: kmm decay legnth XY"),
        gen_tau         = Var("userFloat('gen_tau')",         float, doc = "Gen match: kmm decay time 3D"),
    )
)
        

DileptonsDstarTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","Dstar"),
    cut=cms.string(""),
    name=cms.string("dstar"),
    doc=cms.string("Dstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsDstarTableVariables
)

DileptonsDstarMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","Dstar"),
    cut=cms.string(""),
    name=cms.string("dstar"),
    doc=cms.string("Dstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsDstarMcTableVariables
)
##################################################################################
###
###                              Kstar to Ks pi
###                               '-> Ks to Mu Mu
###                               '-> Ks to pi pi
###
##################################################################################

DileptonsKstarTableVariables =  merge_psets(
    kinematic_pset,
    cms.PSet(
        mm_index       = Var("userInt('mm_index')",          int,   doc = "Index of muon pair"),
        hh_index       = Var("userInt('hh_index')",          int,   doc = "Index of hadron pair"),
        pion_charge    = Var("userInt('pion_charge')",       int,   doc = "Pion charge"),
        mass           = Var("userFloat('mass')",            float, doc = "Raw dstar mass"),
        ks_dist         = Var("userFloat('ks_dist')" ,       float,   doc = "Absolute distance between Kstar and Ks vertices"),
        ks_distErr      = Var("userFloat('ks_distErr')" ,    float,   doc = "Uncertainty on distance between Kstar and Ks vertices"),
        pion_pt        = Var("userFloat('pion_pt')",         float, doc = "Pion pt"),
        pion_eta       = Var("userFloat('pion_eta')",        float, doc = "Pion eta"),
        pion_phi       = Var("userFloat('pion_phi')",        float, doc = "Pion phi"),
        pion_dxy_bs    = Var("userFloat('pion_dxy_bs')",     float, doc = "Pion impact parameter wrt the beam spot"),
        pion_sdxy_bs   = Var("userFloat('pion_sdxy_bs')",    float, doc = "Pion impact parameter significance wrt the beam spot"),
        pv_prob        = Var("userFloat('pv_prob')",         float, doc = "PV refit probability"),
        pv_sum_pt      = Var("userFloat('pv_sum_pt')",       float, doc = "PV sum pt"),
        pv_sum_pt2     = Var("userFloat('pv_sum_pt2')",      float, doc = "PV sum pt^2"),
        pv_ntrks       = Var("userInt('pv_ntrks')",          int,   doc = "PV number of tracks"),
        pv_with_pion_prob = Var("userFloat('pv_with_pion_prob')", float, doc = "PV refit probability with soft pion"),
        # pion_l1_doca   = Var("userFloat('pion_l1_doca')",    float, doc = "Pion distance of closest approach to lepton1"),
        # pion_l2_doca   = Var("userFloat('pion_l2_doca')",    float, doc = "Pion distance of closest approach to lepton2"),
        # Kinematic Fit daugter info
        # nomc_kaon1pt    = Var("userFloat('nomc_kaon1pt')",       float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 pt"),
        # nomc_kaon1eta   = Var("userFloat('nomc_kaon1eta')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 eta"),
        # nomc_kaon1phi   = Var("userFloat('nomc_kaon1phi')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 phi"),
        # jpsimc_kaon1pt    = Var("userFloat('jpsimc_kaon1pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 pt"),
        # jpsimc_kaon1eta   = Var("userFloat('jpsimc_kaon1eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 eta"),
        # jpsimc_kaon1phi   = Var("userFloat('jpsimc_kaon1phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 phi"),
    )
)

DileptonsKstarMcTableVariables = merge_psets(
    DileptonsKstarTableVariables,
    cms.PSet(
        gen_pion_pdgId  = Var("userInt('gen_pion_pdgId')",    int,   doc = "Gen match: pion pdg Id"),
        gen_pion_mpdgId = Var("userInt('gen_pion_mpdgId')",   int,   doc = "Gen match: pion mother pdg Id"),
        gen_pion_pt     = Var("userFloat('gen_pion_pt')",     float, doc = "Gen match: pion pt"),
        gen_pdgId       = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: kstar pdg Id"),
        gen_mpdgId      = Var("userInt('gen_mpdgId')",        int,   doc = "Gen match: kstar mother pdg Id"),
        gen_cpdgId      = Var("userInt('gen_cpdgId')",        int,   doc = "Gen match: kstar common ancestor pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float, doc = "Gen match: kstar mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float, doc = "Gen match: kstar pt"),
        gen_prod_x      = Var("userFloat('gen_prod_x')",      float, doc = "Gen match: kstar mother production vertex x"),
        gen_prod_y      = Var("userFloat('gen_prod_y')",      float, doc = "Gen match: kstar mother production vertex y"),
        gen_prod_z      = Var("userFloat('gen_prod_z')",      float, doc = "Gen match: kstar mother production vertex z"),
        gen_l3d         = Var("userFloat('gen_l3d')",         float, doc = "Gen match: kstar decay legnth 3D"),
        gen_lxy         = Var("userFloat('gen_lxy')",         float, doc = "Gen match: kstar decay legnth XY"),
        gen_tau         = Var("userFloat('gen_tau')",         float, doc = "Gen match: kstar decay time 3D"),
    )
)
        

DileptonsKstarTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","Kstar"),
    cut=cms.string(""),
    name=cms.string("kstar"),
    doc=cms.string("Kstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKstarTableVariables
)

DileptonsKstarMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","Kstar"),
    cut=cms.string(""),
    name=cms.string("kstar"),
    doc=cms.string("Kstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKstarMcTableVariables
)

##################################################################################
###
###                              Bs to l l K K
###
##################################################################################

DileptonsKKmumuTableVariables =  merge_psets(
    kinematic_pset,
    copy_pset(kinematic_pset, {"kin_":"jpsikk_"}),
    copy_pset(kinematic_pset, {"kin_":"phill_"}),
    cms.PSet(
        mm_index        = Var("userInt('mm_index')",           int,   doc = "Index of dimuon pair"),
        kaon1_charge    = Var("userInt('kaon1_charge')",       int,   doc = "Kaon1 charge"),
        kaon1_pt        = Var("userFloat('kaon1_pt')",         float, doc = "Kaon1 pt"),
        kaon1_eta       = Var("userFloat('kaon1_eta')",        float, doc = "Kaon1 eta"),
        kaon1_phi       = Var("userFloat('kaon1_phi')",        float, doc = "Kaon1 phi"),
        kaon1_dxy_bs    = Var("userFloat('kaon1_dxy_bs')",     float, doc = "Kaon1 impact parameter wrt the beam spot"),
        kaon1_sdxy_bs   = Var("userFloat('kaon1_sdxy_bs')",    float, doc = "Kaon1 impact parameter significance wrt the beam spot"),
        # kaon1_mu1_doca  = Var("userFloat('kaon1_mu1_doca')",   float, doc = "Kaon1 distance of closest approach to muon1"),
        # kaon1_mu2_doca  = Var("userFloat('kaon1_mu2_doca')",   float, doc = "Kaon1 distance of closest approach to muon2"),
        kaon2_charge    = Var("userInt('kaon2_charge')",       int,   doc = "Kaon2 charge"),
        kaon2_pt        = Var("userFloat('kaon2_pt')",         float, doc = "Kaon2 pt"),
        kaon2_eta       = Var("userFloat('kaon2_eta')",        float, doc = "Kaon2 eta"),
        kaon2_phi       = Var("userFloat('kaon2_phi')",        float, doc = "Kaon2 phi"),
        kaon2_dxy_bs    = Var("userFloat('kaon2_dxy_bs')",     float, doc = "Kaon2 impact parameter wrt the beam spot"),
        kaon2_sdxy_bs   = Var("userFloat('kaon2_sdxy_bs')",    float, doc = "Kaon2 impact parameter significance wrt the beam spot"),
        kk_mass         = Var("userFloat('kk_mass')",          float, doc = "Mass of two kaons"),
        jpsikk_kk_mass = Var("userFloat('jpsikk_kk_mass')",    float, doc = "Mass of two kaons refitted"),
        # Kinematic Fit daugter info
        jpsikk_kaon1_pt    = Var("userFloat('jpsikk_kaon1_pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 pt"),
        jpsikk_kaon1_eta   = Var("userFloat('jpsikk_kaon1_eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 eta"),
        jpsikk_kaon1_phi   = Var("userFloat('jpsikk_kaon1_phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 phi"),
        jpsikk_kaon2_pt    = Var("userFloat('jpsikk_kaon2_pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 2 pt"),
        jpsikk_kaon2_eta   = Var("userFloat('jpsikk_kaon2_eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 2 eta"),
        jpsikk_kaon2_phi   = Var("userFloat('jpsikk_kaon2_phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 2 phi"),

        # kaon2_mu1_doca  = Var("userFloat('kaon2_mu1_doca')",   float, doc = "Kaon2 distance of closest approach to muon1"),
        # kaon2_mu2_doca  = Var("userFloat('kaon2_mu2_doca')",   float, doc = "Kaon2 distance of closest approach to muon2"),
    #     bmm_nTrks      = Var("userInt('bmm_nTrks')",         int,   doc = "Number of tracks compatible with the vertex by vertex probability (BtoJpsiK as Bmm)"),
    #     bmm_nBMTrks    = Var("userInt('bmm_nBMTrks')",       int,   doc = "Number of tracks more compatible with the mm vertex than with PV by doca significance (BtoJpsiK as Bmm)"),
    #     bmm_nDisTrks   = Var("userInt('bmm_nDisTrks')",      int,   doc = "Number of displaced tracks compatible with the vertex by vertex probability (BtoJpsiK as Bmm)"),
    #     bmm_closetrk   = Var("userInt('bmm_closetrk')",      int,   doc = "Number of tracks compatible with the vertex by doca (BtoJpsiK as Bmm)"),
    #     bmm_closetrks1 = Var("userInt('bmm_closetrks1')",    int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 1 (BtoJpsiK as Bmm)"),
    #     bmm_closetrks2 = Var("userInt('bmm_closetrks2')",    int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 2 (BtoJpsiK as Bmm)"),
    #     bmm_closetrks3 = Var("userInt('bmm_closetrks3')",    int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 3 (BtoJpsiK as Bmm)"),
    #     bmm_docatrk    = Var("userFloat('bmm_docatrk')",     float, doc = "Distance of closest approach of a track to the vertex (BtoJpsiK as Bmm)"),
    #     bmm_m1iso      = Var("userFloat('bmm_m1iso')",       float, doc = "Muon isolation the way it's done in Bmm4 (BtoJpsiK as Bmm)"),
    #     bmm_m2iso      = Var("userFloat('bmm_m2iso')",       float, doc = "Muon isolation the way it's done in Bmm4 (BtoJpsiK as Bmm)"),
    #     bmm_iso        = Var("userFloat('bmm_iso')",         float, doc = "B isolation the way it's done in Bmm4 (BtoJpsiK as Bmm)"),
    #     bmm_bdt        = Var("userFloat('bmm_bdt')",         float, doc = "BDT (BtoJpsiK as Bmm)"),
    )
)

DileptonsKKmumuMcTableVariables = merge_psets(
    DileptonsKKmumuTableVariables,
    cms.PSet(
        gen_kaon1_pdgId  = Var("userInt('gen_kaon1_pdgId')",    int,   doc = "Gen match: kaon1 pdg Id"),
        gen_kaon1_mpdgId = Var("userInt('gen_kaon1_mpdgId')",   int,   doc = "Gen match: kaon1 mother pdg Id"),
        gen_kaon1_pt     = Var("userFloat('gen_kaon1_pt')",     float, doc = "Gen match: kaon1 pt"),
        gen_kaon2_pdgId  = Var("userInt('gen_kaon2_pdgId')",    int,   doc = "Gen match: kaon2 pdg Id"),
        gen_kaon2_mpdgId = Var("userInt('gen_kaon2_mpdgId')",   int,   doc = "Gen match: kaon2 mother pdg Id"),
        gen_kaon2_pt     = Var("userFloat('gen_kaon2_pt')",     float, doc = "Gen match: kaon2 pt"),
        gen_pdgId       = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: kkmm pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float, doc = "Gen match: kkmm mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float, doc = "Gen match: kkmm pt"),
        gen_prod_x      = Var("userFloat('gen_prod_x')",      float, doc = "Gen match: kkmm mother production vertex x"),
        gen_prod_y      = Var("userFloat('gen_prod_y')",      float, doc = "Gen match: kkmm mother production vertex y"),
        gen_prod_z      = Var("userFloat('gen_prod_z')",      float, doc = "Gen match: kkmm mother production vertex z"),
        gen_l3d         = Var("userFloat('gen_l3d')",         float, doc = "Gen match: kkmm decay legnth 3D"),
        gen_lxy         = Var("userFloat('gen_lxy')",         float, doc = "Gen match: kkmm decay legnth XY"),
        gen_tau         = Var("userFloat('gen_tau')",         float, doc = "Gen match: kkmm decay time 3D"),
    )
)

DileptonsKKmumuTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","BToKKmumu"),
    cut=cms.string(""),
    name=cms.string("bkkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKKmumuTableVariables
)

DileptonsKKmumuMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","BToKKmumu"),
    cut=cms.string(""),
    name=cms.string("bkkmm"),
    doc=cms.string("BToKKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKKmumuMcTableVariables
)

DileptonsKKeeTableVariables = copy_pset(DileptonsKKmumuTableVariables, {"mm_index":"ee_index"})
DileptonsKKeeMcTableVariables = copy_pset(DileptonsKKmumuMcTableVariables, {"mm_index":"ee_index"})

DileptonsKKeeTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","BToKKee"),
    cut=cms.string(""),
    name=cms.string("bkkee"),
    doc=cms.string("BToKKee Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKKeeTableVariables
)

DileptonsKKeeMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","BToKKee"),
    cut=cms.string(""),
    name=cms.string("bkkee"),
    doc=cms.string("BToKKee Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsKKeeMcTableVariables
)


##################################################################################
###
###                              B to Mu Mu Gamma
###
##################################################################################

DileptonsMuMuGammaTableVariables =  merge_psets(
    copy_pset(kinematic_pset,{"kin_":"nomc_"}),
    copy_pset(kinematic_pset,{"kin_":"jpsimc_"}),
    cms.PSet(
        mm_index        = Var("userInt('mm_index')",           int,   doc = "Index of dimuon pair"),
        ph_index        = Var("userInt('ph_index')",           int,   doc = "Index of photon"),
        ph_pt           = Var("userFloat('ph_pt')",          float,   doc = "Photon pt"),
        ph_eta          = Var("userFloat('ph_eta')",         float,   doc = "Photon eta"),
        ph_phi          = Var("userFloat('ph_phi')",         float,   doc = "Photon phi"),
        mass            = Var("userFloat('mass')" ,          float,   doc = "Mass - no fit"),
    )
)

DileptonsMuMuGammaMcTableVariables = merge_psets(
    DileptonsMuMuGammaTableVariables,
    cms.PSet(
        gen_ph_pdgId  = Var("userInt('gen_ph_pdgId')",      int,   doc = "Gen match: photon pdg Id"),
        gen_ph_mpdgId = Var("userInt('gen_ph_mpdgId')",     int,   doc = "Gen match: photon mother pdg Id"),
        gen_ph_pt     = Var("userFloat('gen_ph_pt')",       float, doc = "Gen match: photon pt"),
        gen_pdgId     = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: mmg pdg Id"),
        gen_mass      = Var("userFloat('gen_mass')",        float, doc = "Gen match: mmg mass"),
        gen_pt        = Var("userFloat('gen_pt')",          float, doc = "Gen match: mmg pt"),
        gen_prod_x    = Var("userFloat('gen_prod_x')",      float, doc = "Gen match: mmg mother production vertex x"),
        gen_prod_y    = Var("userFloat('gen_prod_y')",      float, doc = "Gen match: mmg mother production vertex y"),
        gen_prod_z    = Var("userFloat('gen_prod_z')",      float, doc = "Gen match: mmg mother production vertex z"),
        gen_l3d       = Var("userFloat('gen_l3d')",         float, doc = "Gen match: mmg decay legnth 3D"),
        gen_lxy       = Var("userFloat('gen_lxy')",         float, doc = "Gen match: mmg decay legnth XY"),
        gen_tau       = Var("userFloat('gen_tau')",         float, doc = "Gen match: mmg decay time 3D"),
    )
)

DileptonsMuMuGammaTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","BToMuMuGamma"),
    cut=cms.string(""),
    name=cms.string("mmg"),
    doc=cms.string("BToMuMuGamma Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsMuMuGammaTableVariables
)

DileptonsMuMuGammaMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","BToMuMuGamma"),
    cut=cms.string(""),
    name=cms.string("mmg"),
    doc=cms.string("BToMuMuGamma Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DileptonsMuMuGammaMcTableVariables
)
##################################################################################
###
###                              Mu Mu Mu
###
##################################################################################

Dileptons3MuTableVariables =  merge_psets(
    copy_pset(kinematic_pset),
    cms.PSet(
        mm_index        = Var("userInt('mm_index')",           int,   doc = "Index of dimuon pair"),
        mu3_index       = Var("userInt('mu3_index')",          int,   doc = "Index of 3d muon"),
        mu3_pt          = Var("userFloat('mu3_pt')",         float,   doc = "3d muon pt"),
        mu3_eta         = Var("userFloat('mu3_eta')",        float,   doc = "3d muon eta"),
        mu3_phi         = Var("userFloat('mu3_phi')",        float,   doc = "3d muon phi"),
        mass            = Var("userFloat('mass')" ,          float,   doc = "Mass - no fit"),
        mm_dist         = Var("userFloat('mm_dist')" ,       float,   doc = "Absolute distance between 2mu and 3mu vertices"),
        mm_distErr      = Var("userFloat('mm_distErr')" ,    float,   doc = "Uncertainty on distance between 2mu and 3mu vertices"),
    )
)

Dileptons3MuMcTableVariables = merge_psets(
    Dileptons3MuTableVariables,
    cms.PSet(
        gen_mu3_pdgId   = Var("userInt(  'gen_mu3_pdgId')",        int,   doc = "Gen match: 3d muon pdg Id"),
        gen_mu3_index   = Var("userInt(  'gen_mu3_index')",        int,   doc = "Gen match: 3d muon index int GenPart"),
        gen_mu3_mpdgId  = Var("userInt(  'gen_mu3_mpdgId')",       int,   doc = "Gen match: 3d muon mother pdg Id"),
        gen_mu3_pt      = Var("userFloat('gen_mu3_pt')",         float,   doc = "Gen match: 3d muon pt"),
        gen_pdgId       = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: mmg pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float, doc = "Gen match: mmg mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float, doc = "Gen match: mmg pt"),
        gen_prod_x      = Var("userFloat('gen_prod_x')",      float, doc = "Gen match: mmg mother production vertex x"),
        gen_prod_y      = Var("userFloat('gen_prod_y')",      float, doc = "Gen match: mmg mother production vertex y"),
        gen_prod_z      = Var("userFloat('gen_prod_z')",      float, doc = "Gen match: mmg mother production vertex z"),
        gen_l3d         = Var("userFloat('gen_l3d')",         float, doc = "Gen match: mmg decay legnth 3D"),
        gen_lxy         = Var("userFloat('gen_lxy')",         float, doc = "Gen match: mmg decay legnth XY"),
        gen_tau         = Var("userFloat('gen_tau')",         float, doc = "Gen match: mmg decay time 3D"),
    )
)

Dileptons3MuTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("Dileptons","MuMuMu"),
    cut=cms.string(""),
    name=cms.string("mmm"),
    doc=cms.string("Triple Muon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = Dileptons3MuTableVariables
)

Dileptons3MuMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DileptonsMc","MuMuMu"),
    cut=cms.string(""),
    name=cms.string("mmm"),
    doc=cms.string("Triple Muon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = Dileptons3MuMcTableVariables
)

##################################################################################
###
###                              B Gen Info
###
##################################################################################

BxToMuMuGen = cms.EDProducer(
    "GenBmmProducer",
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection=cms.InputTag("packedPFCandidates"),
)

BxToMuMuGenVars = cms.PSet(
    signature    = Var("userInt('signature')",     int, doc = "Product of PDG ids of all daughters"),

    # bhadron
    pdgId        = Var("userInt('pdgId')",         int, doc = "PDG id of the b-hadron"),
    pt           = Var("userFloat('pt')",        float, doc = "Pt of the b-hadron"),
    eta          = Var("userFloat('eta')",       float, doc = "Eta of the b-hadron"),
    phi          = Var("userFloat('phi')",       float, doc = "Phi of the b-hadron"),
    mass         = Var("userFloat('mass')",      float, doc = "Mass of the b-hadron"),

    # dimuon
    mu1_pdgId        = Var("userInt('mu1_pdgId')",         int, doc = "PDG id of mu1"),
    mu1_index        = Var("userInt('mu1_index')",         int, doc = "Reco muon index"),
    mu1_good         = Var("userInt('mu1_good')",          int, doc = "Reco muon preselection"),
    mu1_pt           = Var("userFloat('mu1_pt')",        float, doc = "Pt of mu1"),
    mu1_eta          = Var("userFloat('mu1_eta')",       float, doc = "Eta of mu1"),
    mu1_phi          = Var("userFloat('mu1_phi')",       float, doc = "Phi of mu1"),
    mu2_pdgId        = Var("userInt('mu2_pdgId')",         int, doc = "PDG id of mu2"),
    mu2_index        = Var("userInt('mu2_index')",         int, doc = "Reco muon index"),
    mu2_good         = Var("userInt('mu2_good')",          int, doc = "Reco muon preselection"),
    mu2_pt           = Var("userFloat('mu2_pt')",        float, doc = "Pt of mu2"),
    mu2_eta          = Var("userFloat('mu2_eta')",       float, doc = "Eta of mu2"),
    mu2_phi          = Var("userFloat('mu2_phi')",       float, doc = "Phi of mu2"),
    dimuon_mass      = Var("userFloat('dimuon_mass')",   float, doc = "Mass of dimuon"),

    # other daughters
    dau3_pdgId        = Var("userInt('dau3_pdgId')",         int, doc = "PDG id of dau3"),
    dau3_pt           = Var("userFloat('dau3_pt')",        float, doc = "Pt of dau3"),
    dau3_eta          = Var("userFloat('dau3_eta')",       float, doc = "Eta of dau3"),
    dau3_phi          = Var("userFloat('dau3_phi')",       float, doc = "Phi of dau3"),
    dau3_reco_pt      = Var("userFloat('dau3_reco_pt')",   float, doc = "Pt of dau3 PFCandidate match"),
    dau3_reco_eta     = Var("userFloat('dau3_reco_eta')",  float, doc = "Eta of dau3 PFCandidate match"),
    dau3_reco_phi     = Var("userFloat('dau3_reco_phi')",  float, doc = "Phi of dau3 PFCandidate match"),
    dau4_pdgId        = Var("userInt('dau4_pdgId')",         int, doc = "PDG id of dau4"),
    dau4_pt           = Var("userFloat('dau4_pt')",        float, doc = "Pt of dau4"),
    dau4_eta          = Var("userFloat('dau4_eta')",       float, doc = "Eta of dau4"),
    dau4_phi          = Var("userFloat('dau4_phi')",       float, doc = "Phi of dau4"),
    dau4_reco_pt      = Var("userFloat('dau4_reco_pt')",   float, doc = "Pt of dau4 PFCandidate match"),
    dau4_reco_eta     = Var("userFloat('dau4_reco_eta')",  float, doc = "Eta of dau4 PFCandidate match"),
    dau4_reco_phi     = Var("userFloat('dau4_reco_phi')",  float, doc = "Phi of dau4 PFCandidate match"),
    
    # radiation
    rad_p            = Var("userFloat('rad_p')",         float, doc = "P of radiation sum"),
    rad_pt           = Var("userFloat('rad_pt')",        float, doc = "Pt of radiation sum"),
    rad_eta          = Var("userFloat('rad_eta')",       float, doc = "Eta of radiation sum"),
    rad_phi          = Var("userFloat('rad_phi')",       float, doc = "Phi of radiation sum"),

)

BxToMuMuGenTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuGen","genbmm"),
    cut=cms.string(""),
    name=cms.string("genbmm"),
    doc=cms.string("GenInfo Variables for Bmm5"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuGenVars
)

BxToMuMuGenSummaryVars = cms.PSet(
    n_b          = Var("userInt('n_b')",          int, doc = "Number of b-quarks in hard scatter"),
    n_anti_b     = Var("userInt('n_anti_b')",     int, doc = "Number of anti b-quarks in hard scatter"),
    process_type = Var("userInt('process_type')", int, doc = "Process type for Bmm analysis")
)

BxToMuMuGenSummaryTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuGen","gensummary"),
    cut=cms.string(""),
    name=cms.string("gensummary"),
    doc=cms.string("GenInfo Variables for Bmm5"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuGenSummaryVars
)

##################################################################################
###
###                              D Gen Info
###
##################################################################################

DstarGen = cms.EDProducer(
    "GenDstarProducer",
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection=cms.InputTag("packedPFCandidates"),
)

DstarGenVars = cms.PSet(
    signature    = Var("userInt('signature')",     int, doc = "Product of PDG ids of all daughters"),

    # chadron
    pdgId        = Var("userInt('pdgId')",         int, doc = "PDG id of the c-hardon"),
    mpdgId       = Var("userInt('mpdgId')",        int, doc = "PDG id of the mother"),
    pt           = Var("userFloat('pt')",        float, doc = "Pt of the c-hardon"),
    eta          = Var("userFloat('eta')",       float, doc = "Eta of the c-hardon"),
    phi          = Var("userFloat('phi')",       float, doc = "Phi of the c-hardon"),
    mass         = Var("userFloat('mass')",      float, doc = "Mass of the c-hardon"),

    dau1_pdgId        = Var("userInt('dau1_pdgId')",         int, doc = "PDG id of dau1"),
    dau1_mpdgId       = Var("userInt('dau1_mpdgId')",        int, doc = "PDG id of mother"),
    dau1_pt           = Var("userFloat('dau1_pt')",        float, doc = "Pt of dau1"),
    dau1_eta          = Var("userFloat('dau1_eta')",       float, doc = "Eta of dau1"),
    dau1_phi          = Var("userFloat('dau1_phi')",       float, doc = "Phi of dau1"),
    dau1_reco_pt      = Var("userFloat('dau1_reco_pt')",   float, doc = "Pt of dau1 PFCandidate match"),
    dau1_reco_eta     = Var("userFloat('dau1_reco_eta')",  float, doc = "Eta of dau1 PFCandidate match"),
    dau1_reco_phi     = Var("userFloat('dau1_reco_phi')",  float, doc = "Phi of dau1 PFCandidate match"),
    dau1_reco_id      = Var("userInt('dau1_reco_id')",       int, doc = "Track quality of dau1 PFCandidate match"),

    dau2_pdgId        = Var("userInt('dau2_pdgId')",         int, doc = "PDG id of dau2"),
    dau2_mpdgId       = Var("userInt('dau2_mpdgId')",        int, doc = "PDG id of mother"),
    dau2_pt           = Var("userFloat('dau2_pt')",        float, doc = "Pt of dau2"),
    dau2_eta          = Var("userFloat('dau2_eta')",       float, doc = "Eta of dau2"),
    dau2_phi          = Var("userFloat('dau2_phi')",       float, doc = "Phi of dau2"),
    dau2_reco_pt      = Var("userFloat('dau2_reco_pt')",   float, doc = "Pt of dau2 PFCandidate match"),
    dau2_reco_eta     = Var("userFloat('dau2_reco_eta')",  float, doc = "Eta of dau2 PFCandidate match"),
    dau2_reco_phi     = Var("userFloat('dau2_reco_phi')",  float, doc = "Phi of dau2 PFCandidate match"),
    dau2_reco_id      = Var("userInt('dau2_reco_id')",       int, doc = "Track quality of dau2 PFCandidate match"),

    dau3_pdgId        = Var("userInt('dau3_pdgId')",         int, doc = "PDG id of dau3"),
    dau3_mpdgId       = Var("userInt('dau3_mpdgId')",        int, doc = "PDG id of mother"),
    dau3_pt           = Var("userFloat('dau3_pt')",        float, doc = "Pt of dau3"),
    dau3_eta          = Var("userFloat('dau3_eta')",       float, doc = "Eta of dau3"),
    dau3_phi          = Var("userFloat('dau3_phi')",       float, doc = "Phi of dau3"),
    dau3_reco_pt      = Var("userFloat('dau3_reco_pt')",   float, doc = "Pt of dau3 PFCandidate match"),
    dau3_reco_eta     = Var("userFloat('dau3_reco_eta')",  float, doc = "Eta of dau3 PFCandidate match"),
    dau3_reco_phi     = Var("userFloat('dau3_reco_phi')",  float, doc = "Phi of dau3 PFCandidate match"),
    dau3_reco_id      = Var("userInt('dau3_reco_id')",       int, doc = "Track quality of dau3 PFCandidate match"),
    
    # radiation
    rad_p            = Var("userFloat('rad_p')",         float, doc = "P of radiation sum"),
    rad_pt           = Var("userFloat('rad_pt')",        float, doc = "Pt of radiation sum"),
    rad_eta          = Var("userFloat('rad_eta')",       float, doc = "Eta of radiation sum"),
    rad_phi          = Var("userFloat('rad_phi')",       float, doc = "Phi of radiation sum"),

)

DstarGenTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DstarGen","gendstar"),
    cut=cms.string(""),
    name=cms.string("gendstar"),
    doc=cms.string("GenInfo Variables for Dmm"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = DstarGenVars
)

##################################################################################
###
###                             Prescale Info
###
##################################################################################

prescaleTable = cms.EDProducer("TriggerPrescaleProducer",
      trigger       = cms.InputTag("TriggerResults","","HLT"),
      prescales     = cms.InputTag('patTrigger'),
      triggerNames  = cms.vstring('HLT_DoubleMu4_3_Bs', 'HLT_DoubleMu4_3_Jpsi',
                                  'HLT_DoubleMu4_Jpsi_Displaced', 'HLT_DoubleMu4_Jpsi_NoVertexing',
                                  'HLT_DoubleMu4_3_Jpsi_Displaced',
                                  'HLT_DoubleMu4_3_LowMass', 'HLT_Mu4_L1DoubleMu', 'HLT_Mu0_L1DoubleMu',
                                  'HLT_DoubleMu4_LowMass_Displaced')
)

DileptonPlusXSequence   = cms.Sequence(Dileptons * PrimaryVertexInfo)
DileptonPlusXMcSequence = cms.Sequence(DileptonsMc * PrimaryVertexInfoMc * BxToMuMuGen * DstarGen )
DileptonPlusXTables     = cms.Sequence(DileptonsDiMuonTable   * DileptonsHHTable    * DileptonsElElTable     *
                                       DileptonsElMuTable     * DileptonsKmumuTable * DileptonsKeeTable      *
                                       DileptonsKKmumuTable   * DileptonsKKeeTable  * DileptonsDstarTable    *
                                       Dileptons3MuTable      * DileptonsKstarTable * DileptonTrackTable     *
                                       DileptonIsoTable       * 
                                       DileptonsMuMuGammaTable * PrimaryVertexInfoTable * prescaleTable)

DileptonPlusXMcTables   = cms.Sequence(DileptonsDiMuonMcTable * DileptonsHHMcTable     * DileptonsElElMcTable *
                                       DileptonsElMuMcTable   * DileptonsKmumuMcTable  * DileptonsKeeMcTable  *
                                       DileptonsKKmumuMcTable * DileptonsKKeeMcTable   * DileptonsDstarMcTable *
                                       PrimaryVertexInfoMcTable * DileptonsMuMuGammaMcTable * BxToMuMuGenTable *
                                       Dileptons3MuMcTable    * DileptonsKstarMcTable  * DileptonTrackMcTable *
                                       DileptonIsoMcTable     *
                                       BxToMuMuGenSummaryTable * DstarGenTable * prescaleTable)
