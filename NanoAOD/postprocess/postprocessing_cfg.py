# Config file for the postprocessor
from resources_cfg import resources

workdir = "/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv8/src/Bmm5/NanoAOD/postprocess/"
version = 523
input_location = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD"
# output_location = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing-NEW"
output_location = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing"
xrootd_prefix = "root://eoscms.cern.ch:/"
web_report_path = "/afs/cern.ch/user/d/dmytro/www/public_html/bmm5/postprocessing/"

debug = False
tmp_prefix = "tmpPPNA"
require_log_for_success = True

# Test RegEx at https://www.regextester.com/

active_tasks = {
    # Skims with NanoAOD event content
    'NanoAOD-skims':[
        ## 'bkmm', 'trig', 'ks', 'lambda', 'phi', 'ds'
        # 'mm', 'ks', 'phi', 'lambda', 'bkmm',
        # 'trig', # 'mm'
        # 'mm-vtx'
        # 'bmm'
        # 'ks', 'lambda', 'phi'
    ],

    # Skims with task specific content
    'Skims':[
        # 'ks','phi', 'dstar', 'dstar2'
        # 'ks', 'phi', 'dstar'
        # 'ksmm'
    ],

    # Flat ntuples with task specific content
    'FlatNtuples':[
        ## bmm_mva_jpsik, muon_mva, bmm_mva
        # 'fit', 'fit-bkmm', 'bmm_mva_jpsik', 'bmm_mva', 'muon_mva'
        # 'fit', 'fit-bkmm', 'bmm_mva', 'bmm_mva_jpsik', 'dimuon'
        # 'bmm_mva_jpsik'
        # 'bmm_mva_jpsik_loose_vtx'
        # 'trig-eff',
        # 'trig-info'
        # 'fit-bkmm',
        # 'fit'
        # 'fit-bkkmm'
        'muon_mva'

    ]
}

cuts = {
    "fit":
    "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    "Muon_charge[mm_mu1_index] * Muon_charge[mm_mu2_index] < 0 and "\
    "Muon_softMva[mm_mu1_index] > 0.45 and "\
    "abs(mm_kin_mu1eta)<1.4 and "\
    "mm_kin_mu1pt>4 and "\
    "Muon_softMva[mm_mu2_index] > 0.45 and "\
    "abs(mm_kin_mu2eta)<1.4 and "\
    "mm_kin_mu2pt>4 and "\
    "abs(mm_kin_mass-5.4)<0.5 and "\
    "mm_kin_sl3d>6 and "\
    "mm_kin_pt>5.0 and mm_kin_vtx_prob>0.025 and "\
    "HLT_DoubleMu4_3_Bs",

    "fit-bkmm" :
    "mm_mu1_index[bkmm_mm_index]>=0 and "\
    "mm_mu2_index[bkmm_mm_index]>=0 and "\
    "Muon_charge[mm_mu1_index[bkmm_mm_index]] * Muon_charge[mm_mu2_index[bkmm_mm_index]] < 0 and "\
    "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 and "\
    "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 and "\
    "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 and "\
    "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 and "\
    "Muon_softMva[mm_mu1_index[bkmm_mm_index]] > 0.45 and "\
    "Muon_softMva[mm_mu2_index[bkmm_mm_index]] > 0.45 and "\
    "mm_kin_pt[bkmm_mm_index]>7.0 and "\
    "mm_kin_alphaBS[bkmm_mm_index]<0.4 and "\
    "mm_kin_vtx_prob[bkmm_mm_index]>0.1 and "\
    "bkmm_jpsimc_vtx_prob>0.025 and "\
    "mm_kin_sl3d[bkmm_mm_index]>4 and "\
    "abs(bkmm_jpsimc_mass-5.4)<0.5",

    "fit-bkkmm" :
    "mm_mu1_index[bkkmm_mm_index]>=0 and "\
    "mm_mu2_index[bkkmm_mm_index]>=0 and "\
    "Muon_charge[mm_mu1_index[bkkmm_mm_index]] * Muon_charge[mm_mu2_index[bkkmm_mm_index]] < 0 and "\
    "abs(Muon_eta[mm_mu1_index[bkkmm_mm_index]])<1.4 and "\
    "Muon_pt[mm_mu1_index[bkkmm_mm_index]]>4 and "\
    "abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4 and "\
    "Muon_pt[mm_mu2_index[bkkmm_mm_index]]>4 and "\
    "Muon_softMva[mm_mu1_index[bkkmm_mm_index]] > 0.45 and "\
    "Muon_softMva[mm_mu2_index[bkkmm_mm_index]] > 0.45 and "\
    "mm_kin_alphaBS[bkkmm_mm_index]<0.4 and "\
    "mm_kin_vtx_prob[bkkmm_mm_index]>0.1 and "\
    "bkkmm_jpsikk_sl3d>4 and "\
    "bkkmm_jpsikk_vtx_prob>0.025 and "\
    "abs(bkkmm_jpsikk_mass-5.4)<0.5 and "\
    "abs(bkkmm_jpsikk_kk_mass-1.02)<0.03"
}

common_branches = 'PV_npvs|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock'

tasks = [

    ############################
    # Skims
    ############################
    
    {
        'input_pattern':'BuToJpsiK',
        'processor':'Skimmer',
        # trigger based
        "cut" : "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && mm_kin_pt[bkmm_mm_index]>7.0 && mm_kin_alphaXY[bkmm_mm_index]<0.4 && mm_kin_vtx_prob[bkmm_mm_index]>0.1 && bkmm_jpsimc_vtx_prob>0.025 && mm_kin_sl3d[bkmm_mm_index]>4",
        'name':'bkmm',
        'type':'NanoAOD-skims',
        'files_per_job':10
    },
    {
        'input_pattern':'Charmonium',
        'processor':'Skimmer',
        # trigger based
        "cut" : "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && mm_kin_pt[bkmm_mm_index]>7.0 && mm_kin_alphaXY[bkmm_mm_index]<0.4 && mm_kin_vtx_prob[bkmm_mm_index]>0.1 && bkmm_jpsimc_vtx_prob>0.025 && mm_kin_sl3d[bkmm_mm_index]>4",
        'name':'bkmm',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    # {
    #     'input_pattern':'DoubleMuon',
    #     'processor':'Skimmer',
    #     'cut':'HLT_Mu8 || HLT_Mu17',
    #     'name':'trig',
    #     'type':'NanoAOD-skims',
    #     'files_per_job':50
    # },
    {
        'input_pattern':'Charmonium.Run2016',
        'processor':'Skimmer',
        'cut':'HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing',
        'name':'trig',
        'type':'NanoAOD-skims',
        'files_per_job':25
    },
    {
        'input_pattern':'DoubleMuon.Run2016H',
        'processor':'Skimmer',
        'cut':'HLT_DoubleMu0',
        'name':'trig',
        'type':'NanoAOD-skims',
        'files_per_job':25
    },
    {
        'input_pattern':'Charmonium.Run2018|Charmonium.Run2017',
        'processor':'Skimmer',
        'cut':'HLT_Dimuon0_LowMass_L1_0er1p5',
        'name':'trig',
        'type':'NanoAOD-skims',
        'files_per_job':25
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999',
        'name':'ks',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'DoubleElectron|EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'SimpleSkimmer',
        'cut':'ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999',
        'name':'ks',
        "keep": "^(ks_.*|nks|Muon_.*|nMuon|run|event|luminosityBlock|HLT_Ele30_WPTight_Gsf)$",
        'type':'Skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'lambda_kin_slxy>3 && lambda_kin_sipPV<3 && lambda_proton_sip>2 && lambda_pion_sip>2',
        'name':'lambda',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'((phi_trk1_pt > 4 && phi_trk2_pt > 3.0) || ( phi_trk2_pt > 4 && phi_trk1_pt > 3.0)) &&  phi_kin_vtx_prob>0.3 && phi_kin_sipPV<1  && phi_doca<0.004 && phi_kin_lxy<4',
        'name':'phi',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'DoubleElectron|EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        'processor':'SimpleSkimmer',
        'cut':'((phi_trk1_pt > 4 && phi_trk2_pt > 3.0) || ( phi_trk2_pt > 4 && phi_trk1_pt > 3.0)) &&  phi_kin_vtx_prob>0.3 && phi_kin_sipPV<1  && phi_doca<0.004 && phi_kin_lxy<4',
        'name':'phi',
        'type':'Skims',
        "keep": "^(phi_.*|nphi|Muon_.*|nMuon|" + common_branches + ")$",
        'files_per_job':200
    },
    # {
    #     'input_pattern':'DoubleMuon|MINIAODSIM',
    #     'processor':'SimpleSkimmer',
    #     'cut':'d0_dstar_pion_pt > 0',
    #     'name':'dstar',
    #     'type':'Skims',
    #     "keep": "^(d0_.*|nd0|mm_.*|nmm|Muon_.*|nMuon|run|event|luminosityBlock)$",
    #     'files_per_job':100
    # },
    {
        'input_pattern':'DoubleMuon|ZeroBias|MINIAODSIM',
        'processor':'SimpleSkimmer',
        'cut':'dstar_dm_pv > 0',
        'name':'dstar',
        'type':'Skims',
        "keep": "^(dstar_.*|ndstar|d0_.*|nd0|mm_.*|nmm|Muon_.*|nMuon|hh_.*|nhh|run|event|luminosityBlock)$",
        'files_per_job':50
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'phi_ds_pion_pt>0',
        'name':'ds',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'nmm>0',
        'name':'mm',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'InclusiveDileptonMinBias',
        'processor':'SimpleSkimmer',
        'cut':'mm_mass > 0',
        'name':'mm',
        'type':'Skims',
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + common_branches + ")$",
        'files_per_job':100
    },
    {
        'input_pattern':'ParkingDoubleMuonLowMass',
        'processor':'SimpleSkimmer',
        'cut':'mm_mass > 0',
        'name':'mm',
        'type':'Skims',
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + common_branches + ")$",
        'files_per_job':20
    },
    {
        'input_pattern':'ParkingDoubleMuonLowMass',
        'processor':'SimpleSkimmer',
        'cut':'mm_kin_slxy>10&&mm_kin_lxy>1&&mm_kin_alpha<0.01&&mm_mu1_pdgId==-mm_mu2_pdgId&&abs(mm_kin_mass-0.5)<0.2',
        'name':'ksmm',
        'type':'Skims',
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + common_branches + ")$",
        'files_per_job':200
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'mm_mu1_index>=0&&mm_mu2_index>=0&&mm_kin_vtx_prob>0.025&&(Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index])==-1',
        'name':'mm-vtx',
        'type':'NanoAOD-skims',
        'files_per_job':20
    },
    {
        'input_pattern':'LambdaB|B.ToPiPi|B.ToKPi|B.ToKK',
        'processor':'Skimmer',
        'cut':"abs(mm_kin_mu1eta)<1.4 && "\
        "mm_kin_mu1pt>4 && "\
        "abs(mm_kin_mu2eta)<1.4 && "\
        "mm_kin_mu2pt>4 && "\
        "abs(mm_kin_mass-5.4)<0.5 && "\
        "mm_kin_sl3d>6 && "\
        "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025 &&"\
        "mm_mu1_index>=0 && mm_mu2_index>=0 && HLT_DoubleMu4_3_Bs",
        'name':'bmm',
        'type':'NanoAOD-skims',
        'files_per_job':20
    },
    

    ############################
    # Flat Ntuples
    ############################
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|SingleMuon|DoubleMuon',
        'processor':'FlatNtupleForTrigInfo',
        'tree_name':'mm',
        'name':'trig-info',
        'type':'FlatNtuples',
        'files_per_job':25
    },
    {
        'input_pattern':'BsToMuMu_BMuonFilter|BuToJpsiK_BMuonFilter',
        'processor':'FlatNtupleForTrigInfo',
        'tree_name':'mm',
        'name':'trig-info',
        'type':'FlatNtuples',
        'files_per_job':10
    },
    {
        'input_pattern':'SingleMuon',
        'processor':'FlatNtupleForTrigEfficiency',
        'tree_name': 'muons',
        'type':'FlatNtuples',
        "require_muon_tag" : True, 
        'tag_triggers' : ["HLT_IsoMu24", "HLT_IsoMu27"],
        'cut' : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        'name':'trig-eff',
        'files_per_job':25
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'FlatNtupleForTrigEfficiency',
        'tree_name': 'muons',
        'type':'FlatNtuples',
        "require_muon_tag" : False, 
        'tag_triggers' : [
            "L1_Mu6_DoubleEG10er2p5", "L1_Mu6_DoubleEG17er2p5",
            "L1_Mu6_DoubleEG12er2p5", "L1_Mu6_DoubleEG15er2p5",
            "L1_SingleMu3", "L1_SingleMu5", "L1_SingleMu7",
            "L1_Mu6_DoubleEG10", "L1_Mu6_DoubleEG17"
        ],
        'cut' : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        'name':'trig-eff',
        'files_per_job':25
    },
    {
        'input_pattern':'DoubleMuon',
        'processor':'FlatNtupleForTrigEfficiency',
        'tree_name': 'muons',
        'type':'FlatNtuples',
        "require_muon_tag" : True, 
        'tag_triggers' : ["HLT_Mu8", "HLT_Mu17"],
        'cut' : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        'name':'trig-eff',
        'files_per_job':25
    },
    {
        'input_pattern':'BsToMuMu_BMuonFilter',
        'processor':'FlatNtupleForTrigEfficiency',
        'tree_name': 'muons',
        'type':'FlatNtuples',
        "require_muon_tag" : True, 
        'tag_triggers' : ["HLT_Mu8"],
        'cut' : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        'name':'trig-eff',
        'files_per_job':20
    },

    ############### bmm_mva_jpsik ###############
    {
        'input_pattern':'BuToJpsiK',
        'processor':'FlatNtupleForBmmMvaJpsiK',
        'signal_only': True,
        'tree_name': "mva",
        'name':'bmm_mva_jpsik',
        'type':'FlatNtuples',
        'files_per_job':20
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForBmmMvaJpsiK',
        'signal_only': False,
        'tree_name': "mva",
        'name':'bmm_mva_jpsik',
        'type':'FlatNtuples',
        'files_per_job':20
    },
    
    ############### bmm_mva_jpsik_loose_vtx ###############
    {
        'input_pattern':'BuToJpsiK',
        'processor':'FlatNtupleForBmmMvaJpsiK',
        'signal_only': True,
        'loose_vtx': True,
        'tree_name': "mva",
        'name':'bmm_mva_jpsik_loose_vtx',
        'type':'FlatNtuples',
        'files_per_job':20
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForBmmMvaJpsiK',
        'signal_only': False,
        'loose_vtx': True,
        'tree_name': "mva",
        'name':'bmm_mva_jpsik_loose_vtx',
        'type':'FlatNtuples',
        'files_per_job':20
    },

    ################ fit ################
    {
        'input_pattern':'BsToMuMu',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bsmmMc",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
        
    },
    {
        'input_pattern':'BdToMuMu',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bmmMc",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
        
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bmmData",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ### Btohh exclusive

    {
        'input_pattern':'BsToKK_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bskkMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BsToKPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bskpiMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BsToPiPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bspipiMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BdToKK_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bdkkMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BdToKPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bdkpiMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BdToPiPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bdpipiMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ### Btohh inclusive
    
    {
        'input_pattern':'BTohh',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "btohhMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ###
    
    {
        'input_pattern':'LambdaBToPPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "lbppiMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'LambdaBToPK_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "lbpkMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BsToKMuNu_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bskmunuMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BdToPiMuNu_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bdpimunuMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'LambdaBToPMuNu_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "lbpmunuMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BdToMuMuPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bdpimumuMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BuToMuMuPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bupimumuMcBg",
        "blind" : False,
        "cut" : cuts['fit'],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ################ fit-bkmm ################
    {
        'input_pattern':'BuToJpsiK',
        'processor':'FlatNtupleForMLFit',
        'name':'fit-bkmm',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bupsikMc",
        "blind" : False,
        "cut" : cuts["fit-bkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },
    {
        'input_pattern':'BuToJpsiPi',
        'processor':'FlatNtupleForMLFit',
        'name':'fit-bkmm',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bupsipiMc",
        "blind" : False,
        "cut" : cuts["fit-bkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForMLFit',
        'name':'fit-bkmm',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bupsikData",
        "blind" : False,
        "cut" : cuts["fit-bkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },

    ################ fit-bkkmm ################
    {
        'input_pattern':'BsToJPsiPhi',
        'processor':'FlatNtupleForMLFit',
        'name':'fit-bkkmm',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bspsiphiMc",
        "blind" : False,
        "cut" : cuts["fit-bkkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "best_candidate": "",
    },
    {
        'input_pattern':'BdToJpsiKstar',
        'processor':'FlatNtupleForMLFit',
        'name':'fit-bkkmm',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bdpsikstarMc",
        "blind" : False,
        "cut" : cuts["fit-bkkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "best_candidate": "",
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForMLFit',
        'name':'fit-bkkmm',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bspsiphiData",
        "blind" : False,
        "cut" : cuts["fit-bkkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "best_candidate": "",
    },
    
    ################ dimuon ################
    # {
    #     'input_pattern':'BsToMuMu',
    #     'processor':'FlatNtupleForMLFit',
    #     'name':'fit',
    #     'type':'FlatNtuples',
    #     'files_per_job':20,
    #     "tree_name" : "bsmmMc",
    #     "blind" : False,
    #     "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    #             # "Muon_softMvaId[mm_mu1_index] and "\
    #             "abs(mm_kin_mu1eta)<1.4 and "\
    #             "mm_kin_mu1pt>4 and "\
    #             # "Muon_softMvaId[mm_mu2_index] and "\
    #             "abs(mm_kin_mu2eta)<1.4 and "\
    #             "mm_kin_mu2pt>4 and "\
    #             "abs(mm_kin_mass-5.4)<0.5 and "\
    #             "mm_kin_sl3d>6 and "\
    #             "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025",
    #     "final_state" : "mm",
    #     # "best_candidate": "mm_kin_pt",
    #     "best_candidate": "",
        
    # },
    # {
    #     'input_pattern':'BdToMuMu',
    #     'processor':'FlatNtupleForMLFit',
    #     'name':'fit',
    #     'type':'FlatNtuples',
    #     'files_per_job':20,
    #     "tree_name" : "bmmMc",
    #     "blind" : False,
    #     "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    #             # "Muon_softMvaId[mm_mu1_index] and "\
    #             "abs(mm_kin_mu1eta)<1.4 and "\
    #             "mm_kin_mu1pt>4 and "\
    #             # "Muon_softMvaId[mm_mu2_index] and "\
    #             "abs(mm_kin_mu2eta)<1.4 and "\
    #             "mm_kin_mu2pt>4 and "\
    #             "abs(mm_kin_mass-5.4)<0.5 and "\
    #             "mm_kin_sl3d>6 and "\
    #             "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025",
    #     "final_state" : "mm",
    #     # "best_candidate": "mm_kin_pt",
    #     "best_candidate": "",
        
    # },
    {
        'input_pattern':'Jpsi',
        'processor':'FlatNtupleForMLFit',
        'name':'dimuon',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "jpsiMC",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMva[mm_mu1_index]>0.45 and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMva[mm_mu2_index]>0.45 and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mass>2.9 and "\
                "mm_kin_mass<4.0 and "\
                "mm_kin_vtx_prob>0.1",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'YnSToMuMu',
        'processor':'FlatNtupleForMLFit',
        'name':'dimuon',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "ynsMC",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMva[mm_mu1_index]>0.45 and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMva[mm_mu2_index]>0.45 and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mass>8.5 and "\
                "mm_kin_mass<11.5 and "\
                "mm_kin_vtx_prob>0.1",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForMLFit',
        'name':'dimuon',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "mmData",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMva[mm_mu1_index]>0.45 and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMva[mm_mu2_index]>0.45 and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mass>2.9 and "\
                "mm_kin_mass<4.0 and "\
                "mm_kin_vtx_prob>0.1",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'MuOnia',
        'processor':'FlatNtupleForMLFit',
        'name':'dimuon',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "mmData",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMva[mm_mu1_index]>0.45 and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMva[mm_mu2_index]>0.45 and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mass>8.0 and "\
                "mm_kin_mass<12.0 and "\
                "mm_kin_vtx_prob>0.1",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
        
    ################ muon_mva ################
    {
        'input_pattern':'InclusiveDileptonMinBias', # 'QCD_Pt.*?_MuEnriched|Charmonium|BuToJpsiK',
        'processor':'FlatNtupleForMuonMVA',
        'name':'muon_mva',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "muons",
    },

    ################ bmm_mva ################
    {
        'input_pattern':'Charmonium',
        'processor':'FlatNtupleForBmmMva',
        'signal_only': False,
        'tree_name': "mva",
        'blind': True,
        'name':'bmm_mva',
        'type':'FlatNtuples',
        'files_per_job':20
    },
    {
        'input_pattern':'BsToMuMu_BMuonFilter',
        'processor':'FlatNtupleForBmmMva',
        'signal_only': True,
        'tree_name': "mva",
        'blind': False,
        'name':'bmm_mva',
        'type':'FlatNtuples',
        'files_per_job':20
    },
]
