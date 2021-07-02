# Config file for the postprocessor
from resources_cfg import resources

workdir = "/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/postprocess/"
version = 514
input_location = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD"
output_location = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing"
xrootd_prefix = "root://eoscms.cern.ch:/"
web_report_path = "/afs/cern.ch/user/d/dmytro/www/public_html/bmm5/postprocessing/"

debug = False
tmp_prefix = "tmpPPNA"
require_log_for_success = True

# Test RegEx at https://www.regextester.com/

active_tasks = {
    'NanoAOD-skims':[
        # 'bkmm', 'trig', 'ks', 'lambda', 'phi', 'ds'
        'mm'
    ],
    'FlatNtuples':[
        # bmm_mva_jpsik, muon_mva, bmm_mva
        'fit', 'fit-bkmm', 'bmm_mva_jpsik'
    ]
}

tasks = [

    ############################
    # Skims
    ############################
    
    {
        'input_pattern':'BuToJpsiK',
        'processor':'Skimmer',
        # trigger based
        "cut" : "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && mm_kin_pt[bkmm_mm_index]>7.0 && mm_kin_alphaXY[bkmm_mm_index]<0.4 && mm_kin_vtx_prob[bkmm_mm_index]>0.1 && bkmm_jpsimc_vtx_prob>0.1 && mm_kin_sl3d[bkmm_mm_index]>4",
        'name':'bkmm',
        'type':'NanoAOD-skims',
        'files_per_job':10
    },
    {
        'input_pattern':'Charmonium',
        'processor':'Skimmer',
        # trigger based
        "cut" : "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && mm_kin_pt[bkmm_mm_index]>7.0 && mm_kin_alphaXY[bkmm_mm_index]<0.4 && mm_kin_vtx_prob[bkmm_mm_index]>0.1 && bkmm_jpsimc_vtx_prob>0.1 && mm_kin_sl3d[bkmm_mm_index]>4",
        'name':'bkmm',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'SingleMuon',
        'processor':'Skimmer',
        'cut':'HLT_Mu8>0 or HLT_Mu3_PFJet40>0',
        'name':'trig',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'nks>0',
        'name':'ks',
        'type':'NanoAOD-skims',
        'files_per_job':100
    },
    {
        'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM',
        # 'input_pattern':'EGamma|DoubleEG|SingleElectron|SinglePhoton',
        'processor':'Skimmer',
        'cut':'nlambda>0',
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


    ############################
    # Flat Ntuples
    ############################

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

    ################ fit ################
    {
        'input_pattern':'BsToMuMu',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bsmmMc",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "blind" : True,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5 and "\
                "HLT_DoubleMu4_3_Bs",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BsToKK_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bskkMcBg",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'LambdaBToPPi_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "lbppiMcBg",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        'input_pattern':'BsToKPiNu_',
        'processor':'FlatNtupleForMLFit',
        'name':'fit',
        'type':'FlatNtuples',
        'files_per_job':20,
        "tree_name" : "bspimunuMcBg",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
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
        "cut" :
            "mm_mu1_index[bkmm_mm_index]>=0 and "\
            "mm_mu2_index[bkmm_mm_index]>=0 and "\
            "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 and "\
            "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 and "\
            "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 and "\
            "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 and "\
            "mm_kin_pt[bkmm_mm_index]>7.0 and "\
            "mm_kin_alphaBS[bkmm_mm_index]<0.4 and "\
            "mm_kin_vtx_prob[bkmm_mm_index]>0.1 and "\
            "bkmm_jpsimc_vtx_prob>0.1 and "\
            "mm_kin_sl3d[bkmm_mm_index]>4 and "\
            "abs(bkmm_jpsimc_mass-5.4)<0.5",
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
        "cut" :
            "mm_mu1_index[bkmm_mm_index]>=0 and "\
            "mm_mu2_index[bkmm_mm_index]>=0 and "\
            "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 and "\
            "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 and "\
            "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 and "\
            "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 and "\
            "mm_kin_pt[bkmm_mm_index]>7.0 and "\
            "mm_kin_alphaBS[bkmm_mm_index]<0.4 and "\
            "mm_kin_vtx_prob[bkmm_mm_index]>0.1 and "\
            "bkmm_jpsimc_vtx_prob>0.1 and "\
            "mm_kin_sl3d[bkmm_mm_index]>4 and "\
            "abs(bkmm_jpsimc_mass-5.4)<0.5",
        "final_state" : "bkmm",
        "best_candidate": "",
    },

    ################ muon_mva ################
    {
        'input_pattern':'QCD_Pt.*?_MuEnriched|Charmonium|BuToJpsiK',
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
