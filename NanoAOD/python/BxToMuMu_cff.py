from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms
import re

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

def copy_pset(pset,replace_dict):
    result = cms.PSet()
    for name in pset.parameters_().keys():
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

##################################################################################
###
###                              Bx to mu mu
###
##################################################################################

BxToMuMu = cms.EDProducer("BxToMuMuProducer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection=cms.InputTag("packedPFCandidates"),
    prunedGenParticleCollection = cms.InputTag("prunedGenParticles"),
    MuonMinPt=cms.double(1.),
    MuonMaxEta=cms.double(2.4),
    KaonMinPt=cms.double(1.),
    KaonMaxEta=cms.double(2.4),
    KaonMinDCASig=cms.double(-1.),
    DiMuonChargeCheck=cms.bool(False),
    minBKmmMass = cms.double(4.5),
    maxBKmmMass = cms.double(6.0),
    minBKKmmMass = cms.double(4.5),
    maxBKKmmMass = cms.double(6.0),
    maxTwoTrackDOCA = cms.double(0.1),
    bdtEvent0 = cms.FileInPath('Bmm5/NanoAOD/data/TMVA-100-Events0_BDT.weights.xml'),
    bdtEvent1 = cms.FileInPath('Bmm5/NanoAOD/data/TMVA-100-Events1_BDT.weights.xml'),
    bdtEvent2 = cms.FileInPath('Bmm5/NanoAOD/data/TMVA-100-Events2_BDT.weights.xml'),
    isMC = cms.bool(False)
)

BxToMuMuMc = BxToMuMu.clone( isMC = cms.bool(True) ) 

kinematic_pset = cms.PSet(
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability"),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2"),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass"),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt"),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta"),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error"),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot"),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot"),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS"),
    kin_alpha = Var("userFloat('kin_alpha')",          float, doc = "Kinematic fit: pointing angle in 3D wrt PV"),
    kin_vtx_x    = Var("userFloat('kin_vtx_x')",       float, doc = "Kinematic fit: vertex x"),
    kin_vtx_xErr = Var("userFloat('kin_vtx_xErr')",    float, doc = "Kinematic fit: vertex x-error"),
    kin_vtx_y    = Var("userFloat('kin_vtx_y')",       float, doc = "Kinematic fit: vertex y"),
    kin_vtx_yErr = Var("userFloat('kin_vtx_yErr')",    float, doc = "Kinematic fit: vertex y-error"),
    kin_vtx_z    = Var("userFloat('kin_vtx_z')",       float, doc = "Kinematic fit: vertex y"),
    kin_vtx_zErr = Var("userFloat('kin_vtx_zErr')",    float, doc = "Kinematic fit: vertex y-error"),
    kin_pv_z     = Var("userFloat('kin_pv_z')",        float, doc = "Kinematic fit: primary vertex z"),
    kin_pv_zErr  = Var("userFloat('kin_pv_zErr')",     float, doc = "Kinematic fit: primary vertex z-error"),
    kin_l3d      = Var("userFloat('kin_l3d')",         float, doc = "Kinematic fit: decay length wrt Primary Vertex in 3D"),
    kin_sl3d     = Var("userFloat('kin_sl3d')",        float, doc = "Kinematic fit: decay length significance wrt Primary Vertex in 3D"),
    kin_tau      = Var("userFloat('kin_tau')",         float, doc = "Kinematic fit: decay time wrt Primary Vertex in 3D"),
    kin_taue     = Var("userFloat('kin_taue')",        float, doc = "Kinematic fit: decay time error wrt Primary Vertex in 3D"),
    kin_tauxy    = Var("userFloat('kin_tauxy')",       float, doc = "Kinematic fit: decay time wrt Primary Vertex in XY"),
    kin_tauxye   = Var("userFloat('kin_tauxye')",      float, doc = "Kinematic fit: decay time error wrt Primary Vertex in XY"),
    kin_pvip     = Var("userFloat('kin_pvip')",        float, doc = "Kinematic fit: impact parameter wrt Primary Vertex in 3D"),
    kin_pvipErr  = Var("userFloat('kin_pvipErr')",  float, doc = "Kinematic fit: impact parameter uncertainty wrt Primary Vertex in 3D"),
    kin_pvlip    = Var("userFloat('kin_pvlip')",       float, doc = "Kinematic fit: longitudinal impact parameter wrt Primary Vertex"),
    kin_pvlipErr = Var("userFloat('kin_pvlipErr')",    float, doc = "Kinematic fit: longitudinal impact parameter uncertainty wrt Primary Vertex"),
    kin_pv2lip   = Var("userFloat('kin_pv2lip')",       float, doc = "Kinematic fit: longitudinal impact parameter wrt Second best Primary Vertex"),
    kin_pv2lipErr = Var("userFloat('kin_pv2lipErr')",    float, doc = "Kinematic fit: longitudinal impact parameter uncertainty wrt Second best Primary Vertex"),
)

BxToMuMuDiMuonTableVariables = merge_psets(
    cms.PSet(
        mu1_index    = Var("userInt('mu1_index')",         int,   doc = "Index of corresponding leading muon"),
        mu2_index    = Var("userInt('mu2_index')",         int,   doc = "Index of corresponding subleading muon"),
        mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
        doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of muons"),
        nTrks        = Var("userInt('nTrks')",             int,   doc = "Number of tracks compatible with the vertex by vertex probability"),
        nBMTrks      = Var("userInt('nBMTrks')",           int,   doc = "Number of tracks more compatible with the mm vertex than with PV by doca significance"),
        nDisTrks     = Var("userInt('nDisTrks')",          int,   doc = "Number of displaced tracks compatible with the vertex by vertex probability"),
        closetrk     = Var("userInt('closetrk')",          int,   doc = "Number of tracks compatible with the vertex by doca"),
        closetrks1   = Var("userInt('closetrks1')",        int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 1"),
        closetrks2   = Var("userInt('closetrks2')",        int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 2"),
        closetrks3   = Var("userInt('closetrks3')",        int,   doc = "Number of tracks compatible with the vertex with doca signifance less than 3"),
        docatrk      = Var("userFloat('docatrk')",         float, doc = "Distance of closest approach of a track to the vertex"),
        m1iso        = Var("userFloat('m1iso')",           float, doc = "Muon isolation the way it's done in Bmm4"),
        m2iso        = Var("userFloat('m2iso')",           float, doc = "Muon isolation the way it's done in Bmm4"),
        iso          = Var("userFloat('iso')",             float, doc = "B isolation the way it's done in Bmm4"),
        bdt          = Var("userFloat('bdt')",             float, doc = "BDT"),
        otherVtxMaxProb = Var("userFloat('otherVtxMaxProb')", float, doc = "Max vertexing probability of one of the muons with a random track with minPt=0.5GeV"),
        otherVtxMaxProb1 = Var("userFloat('otherVtxMaxProb1')", float, doc = "Max vertexing probability of one of the muons with a random track with minPt=1.0GeV"),
        otherVtxMaxProb2 = Var("userFloat('otherVtxMaxProb2')", float, doc = "Max vertexing probability of one of the muons with a random track with minPt=2.0GeV"),
        # Kalman Fit
        kal_valid    = Var("userInt('kalman_valid')",      int,   doc = "Kalman vertex fit validity"),
        kal_vtx_prob = Var("userFloat('kalman_vtx_prob')", float, doc = "Kalman fit vertex probability"),
        kal_mass     = Var("userFloat('kalman_mass')",     float, doc = "Kalman vertex refitted mass"),
        kal_lxy      = Var("userFloat('kalman_lxy')",      float, doc = "Kalman fit vertex displacement in XY plane"),
        kal_slxy     = Var("userFloat('kalman_sigLxy')",   float, doc = "Kalman fit vertex displacement significance in XY plane"),
        # Kinematic Fit daugter info
        kin_mu1pt    = Var("userFloat('kin_mu1pt')",       float, doc = "Kinematic fit: refitted muon 1 pt"),
        kin_mu1eta   = Var("userFloat('kin_mu1eta')",      float, doc = "Kinematic fit: refitted muon 1 eta"),
        kin_mu1phi   = Var("userFloat('kin_mu1phi')",      float, doc = "Kinematic fit: refitted muon 1 phi"),
        kin_mu2pt    = Var("userFloat('kin_mu2pt')",       float, doc = "Kinematic fit: refitted muon 2 pt"),
        kin_mu2eta   = Var("userFloat('kin_mu2eta')",      float, doc = "Kinematic fit: refitted muon 2 eta"),
        kin_mu2phi   = Var("userFloat('kin_mu2phi')",      float, doc = "Kinematic fit: refitted muon 2 phi"),
        ),
    kinematic_pset
)

BxToMuMuDiMuonMcTableVariables = merge_psets(
    BxToMuMuDiMuonTableVariables,
    cms.PSet(
        gen_mu1_pdgId  = Var("userInt(  'gen_mu1_pdgId')",    int,   doc = "Gen match: first muon pdg Id"),
        gen_mu1_mpdgId = Var("userInt(  'gen_mu1_mpdgId')",   int,   doc = "Gen match: first muon mother pdg Id"),
        gen_mu1_pt     = Var("userFloat('gen_mu1_pt')",     float,   doc = "Gen match: first muon pt"),
        gen_mu2_pdgId  = Var("userInt(  'gen_mu2_pdgId')",    int,   doc = "Gen match: second muon pdg Id"),
        gen_mu2_mpdgId = Var("userInt(  'gen_mu2_mpdgId')",   int,   doc = "Gen match: second muon mother pdg Id"),
        gen_mu2_pt     = Var("userFloat('gen_mu2_pt')",     float,   doc = "Gen match: second muon pt"),
        gen_pdgId      = Var("userInt(  'gen_pdgId')",        int,   doc = "Gen match: dimuon pdg Id"),
        gen_mpdgId     = Var("userInt(  'gen_mpdgId')",       int,   doc = "Gen match: dimuon mother pdg Id"),
        gen_mass       = Var("userFloat('gen_mass')",       float,   doc = "Gen match: dimuon mass"),
        gen_pt         = Var("userFloat('gen_pt')",         float,   doc = "Gen match: dimuon pt"),
        gen_prod_z     = Var("userFloat('gen_prod_z')",     float,   doc = "Gen match: dimuon mother production vertex z"),
        gen_vtx_x      = Var("userFloat('gen_vtx_x')",      float,   doc = "Gen match: dimuon vertex x"),
        gen_vtx_y      = Var("userFloat('gen_vtx_y')",      float,   doc = "Gen match: dimuon vertex y"),
        gen_vtx_z      = Var("userFloat('gen_vtx_z')",      float,   doc = "Gen match: dimuon vertex z"),
        gen_l3d        = Var("userFloat('gen_l3d')",        float,   doc = "Gen match: dimuon decay legnth 3D"),
        gen_lxy        = Var("userFloat('gen_lxy')",        float,   doc = "Gen match: dimuon decay legnth XY"),
        ),
)

BxToMuMuDiMuonTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMu","DiMuon"),
    cut=cms.string(""),
    name=cms.string("mm"),
    doc=cms.string("Dimuon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuDiMuonTableVariables
)

BxToMuMuDiMuonMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuMc","DiMuon"),
    cut=cms.string(""),
    name=cms.string("mm"),
    doc=cms.string("Dimuon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuDiMuonMcTableVariables
)

##################################################################################
###
###                              B to K mu mu
###
##################################################################################

BxToMuMuBToKmumuTableVariables =  merge_psets(
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
        kaon_mu1_doca  = Var("userFloat('kaon_mu1_doca')",   float, doc = "Kaon distance of closest approach to muon1"),
        kaon_mu2_doca  = Var("userFloat('kaon_mu2_doca')",   float, doc = "Kaon distance of closest approach to muon2"),
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
        # Kinematic Fit daugter info
        nomc_kaon1pt    = Var("userFloat('nomc_kaon1pt')",       float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 pt"),
        nomc_kaon1eta   = Var("userFloat('nomc_kaon1eta')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 eta"),
        nomc_kaon1phi   = Var("userFloat('nomc_kaon1phi')",      float, doc = "Kinematic fit (no Jpsi mass constraint): refitted kaon 1 phi"),
        jpsimc_kaon1pt    = Var("userFloat('jpsimc_kaon1pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 pt"),
        jpsimc_kaon1eta   = Var("userFloat('jpsimc_kaon1eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 eta"),
        jpsimc_kaon1phi   = Var("userFloat('jpsimc_kaon1phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 phi"),
    )
)

BxToMuMuBToKmumuMcTableVariables = merge_psets(
    BxToMuMuBToKmumuTableVariables,
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
    )
)
        

BxToMuMuBToKmumuTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMu","BToKmumu"),
    cut=cms.string(""),
    name=cms.string("bkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuBToKmumuTableVariables
)

BxToMuMuBToKmumuMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuMc","BToKmumu"),
    cut=cms.string(""),
    name=cms.string("bkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuBToKmumuMcTableVariables
)

##################################################################################
###
###                              Bs to Jpsi K K
###
##################################################################################

BxToMuMuBToKKmumuTableVariables =  merge_psets(
    copy_pset(kinematic_pset,{"kin_":"jpsikk_"}),
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
        jpsikk_kaon1pt    = Var("userFloat('jpsikk_kaon1pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 pt"),
        jpsikk_kaon1eta   = Var("userFloat('jpsikk_kaon1eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 eta"),
        jpsikk_kaon1phi   = Var("userFloat('jpsikk_kaon1phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 1 phi"),
        jpsikk_kaon2pt    = Var("userFloat('jpsikk_kaon2pt')",       float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 2 pt"),
        jpsikk_kaon2eta   = Var("userFloat('jpsikk_kaon2eta')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 2 eta"),
        jpsikk_kaon2phi   = Var("userFloat('jpsikk_kaon2phi')",      float, doc = "Kinematic fit (with Jpsi mass constraint): refitted kaon 2 phi"),

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

BxToMuMuBToKKmumuMcTableVariables = merge_psets(
    BxToMuMuBToKKmumuTableVariables,
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
    )
)

BxToMuMuBToKKmumuTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMu","BToKKmumu"),
    cut=cms.string(""),
    name=cms.string("bkkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuBToKKmumuTableVariables
)

BxToMuMuBToKKmumuMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuMc","BToKKmumu"),
    cut=cms.string(""),
    name=cms.string("bkkmm"),
    doc=cms.string("BToKKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuBToKKmumuMcTableVariables
)

##################################################################################
###
###                              Gen Summary Info
###
##################################################################################

BxToMuMuGen = cms.EDProducer("GenBmmProducer")

BxToMuMuGenVars = cms.PSet(
    # bhadron
    pdgId        = Var("userInt('pdgId')",         int, doc = "PDG id of the b-hadron"),
    pt           = Var("userFloat('pt')",        float, doc = "Pt of the b-hadron"),
    eta          = Var("userFloat('eta')",       float, doc = "Eta of the b-hadron"),
    phi          = Var("userFloat('phi')",       float, doc = "Phi of the b-hadron"),
    mass         = Var("userFloat('mass')",      float, doc = "Mass of the b-hadron"),

    # dimuon
    mu1_pdgId        = Var("userInt('mu1_pdgId')",         int, doc = "PDG id of mu1"),
    mu1_pt           = Var("userFloat('mu1_pt')",        float, doc = "Pt of mu1"),
    mu1_eta          = Var("userFloat('mu1_eta')",       float, doc = "Eta of mu1"),
    mu1_phi          = Var("userFloat('mu1_phi')",       float, doc = "Phi of mu1"),
    mu2_pdgId        = Var("userInt('mu2_pdgId')",         int, doc = "PDG id of mu2"),
    mu2_pt           = Var("userFloat('mu2_pt')",        float, doc = "Pt of mu2"),
    mu2_eta          = Var("userFloat('mu2_eta')",       float, doc = "Eta of mu2"),
    mu2_phi          = Var("userFloat('mu2_phi')",       float, doc = "Phi of mu2"),
    dimuon_mass      = Var("userFloat('dimuon_mass')",   float, doc = "Mass of dimuon"),

    # kaons
    kaon1_pdgId        = Var("userInt('kaon1_pdgId')",         int, doc = "PDG id of kaon1"),
    kaon1_pt           = Var("userFloat('kaon1_pt')",        float, doc = "Pt of kaon1"),
    kaon1_eta          = Var("userFloat('kaon1_eta')",       float, doc = "Eta of kaon1"),
    kaon1_phi          = Var("userFloat('kaon1_phi')",       float, doc = "Phi of kaon1"),
    kaon2_pdgId        = Var("userInt('kaon2_pdgId')",         int, doc = "PDG id of kaon2"),
    kaon2_pt           = Var("userFloat('kaon2_pt')",        float, doc = "Pt of kaon2"),
    kaon2_eta          = Var("userFloat('kaon2_eta')",       float, doc = "Eta of kaon2"),
    kaon2_phi          = Var("userFloat('kaon2_phi')",       float, doc = "Phi of kaon2"),
    kk_mass            = Var("userFloat('kk_mass')",         float, doc = "Mass of kk"),
    
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

BxToMuMuSequence   = cms.Sequence(BxToMuMu)
BxToMuMuMcSequence = cms.Sequence(BxToMuMuMc * BxToMuMuGen )
BxToMuMuTables     = cms.Sequence(BxToMuMuDiMuonTable   * BxToMuMuBToKmumuTable * BxToMuMuBToKKmumuTable)
BxToMuMuMcTables   = cms.Sequence(BxToMuMuDiMuonMcTable * BxToMuMuBToKmumuMcTable * BxToMuMuBToKKmumuMcTable * BxToMuMuGenTable * BxToMuMuGenSummaryTable)
