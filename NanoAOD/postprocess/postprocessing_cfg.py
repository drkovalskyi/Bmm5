# Config file for the postprocessor
from resources_cfg import resources

workdir = "/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv8/src/Bmm5/NanoAOD/postprocess/"
version = 'crab-140x-mm'

input_location = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/" + str(version)
# input_location = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/523/bkmm"

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
    "NanoAOD-skims":[
        ## "bkmm", "trig", "ks", "lambda", "phi", "ds"
        # "mm", "ks", "phi", "lambda", "bkmm",
        # "trig", # "mm"
        # "mm-vtx"
        # "bmm"
        # "ks", "lambda", "phi"
    ],

    # Skims with task specific content
    "Skims":[
        # "ks","phi", "dstar", "dstar2"
        # "ks", "phi", "dstar"
        # "ks", "phi"
        # "ksmm"
        #"em", "dzpipi"
        # "dzmm"
        # "bkmm"
        # "dstar_mm"
        # "trig"
        # "mm_mva0p9"
        # "tau3mu"
        # "mm_vtx"
        # "dzkpimm"
    ],

    # Flat ntuples with task specific content
    "FlatNtuples":[
        ## bmm_mva_jpsik, muon_mva, bmm_mva
        # "fit", "fit-bkmm", "bmm_mva_jpsik", "bmm_mva", "muon_mva"
        # "fit", "fit-bkmm", "bmm_mva", "bmm_mva_jpsik", "dimuon"
        # "bmm_mva_jpsik"
        # "bmm_mva_jpsik_loose_vtx"
        # "trig-eff",
        # "trig-info"
        # "fit"
        # "fit-bkkmm"
        # "muon_mva"
        # "fit-em"
        # "dstar", "dzpipi", "dzkpi"
        # "dzpipi"
        # "dzpipi_otherZB"
        # "dzkpi"
        # "ksmm",
        # "kspipi"
        # "bkkmm"
        "bkmm"
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

    "fit-em":
    "em_mu_index>=0 and em_el_index>=0 and "\
    "em_mu_pt>5 and em_el_pt>5 and "\
    "Muon_mediumId[em_mu_index] and Electron_mvaNoIso_WPL[em_el_index] and "\
    "abs(em_kin_mass-5.4)<1.0 and "\
    "em_kin_vtx_prob>0.025",

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

    "bkmm_skim" : (
        "bkmm_jpsimc_vtx_prob>0.025 and bkmm_jpsimc_sl3d>3"
        " and abs(bkmm_jpsimc_mass-5.4)<0.5 && bkmm_jpsimc_alpha<0.2"
    ),
    
    # BuToJpsiK validation
    "bkmm-validation" : (
        "mm_mu1_index[bkmm_mm_index]>=0 and "
        "mm_mu2_index[bkmm_mm_index]>=0 and "
        "mm_mu1_pdgId[bkmm_mm_index] * mm_mu2_pdgId[bkmm_mm_index] < 0 and "
        "bkmm_jpsimc_vtx_prob>0.1 and "
        "bkmm_jpsimc_sl3d>5 and "
        "abs(bkmm_jpsimc_alpha)<0.01 and "
        "abs(bkmm_jpsimc_mass-5.4)<0.5"
    ),
    
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
    "abs(bkkmm_jpsikk_kk_mass-1.02)<0.03",

    # BsToJpsiPhi validation
    "bkkmm" : (
        "mm_mu1_index[bkkmm_mm_index]>=0 and "
        "mm_mu2_index[bkkmm_mm_index]>=0 and "
        "mm_mu1_pdgId[bkkmm_mm_index] * mm_mu2_pdgId[bkkmm_mm_index] < 0 and "
        "abs(bkkmm_jpsikk_alpha)<0.01 and "
        "bkkmm_jpsikk_sl3d>5 and "
        "bkkmm_jpsikk_vtx_prob>0.1 and "
        "abs(bkkmm_jpsikk_mass-5.4)<0.5 and "
        "abs(bkkmm_kk_mass-1.02)<0.01"
    ),

    "dstar_dzpipi" : (
        "dstar_hh_index>=0 and "
        "max(hh_had1_pt[dstar_hh_index], hh_had2_pt[dstar_hh_index])>4 and "
        "min(hh_had1_pt[dstar_hh_index], hh_had2_pt[dstar_hh_index])>3 and "
        "hh_kin_alpha[dstar_hh_index]<0.1 and hh_kin_sl3d[dstar_hh_index]>3 and "
        "hh_kin_vtx_prob[dstar_hh_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
        "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
        "hh_kin_mass[dstar_hh_index]>1.81 and hh_kin_mass[dstar_hh_index]<1.94 and "
        # "hh_kin_mass[dstar_mm_index]>1.65 and hh_kin_mass[dstar_mm_index]<2.4 and "
        "hh_had1_pdgId[dstar_hh_index] * hh_had2_pdgId[dstar_hh_index] == - 211 * 211"
    ),
    "dstar_dzkpi" : (
        "dstar_hh_index>=0 and "
        "max(hh_had1_pt[dstar_hh_index], hh_had2_pt[dstar_hh_index])>4 and "
        "min(hh_had1_pt[dstar_hh_index], hh_had2_pt[dstar_hh_index])>3 and "
        # "hh_had1_pt[dstar_hh_index]>4 and hh_had2_pt[dstar_hh_index]>3 and "
        "hh_kin_alpha[dstar_hh_index]<0.1 and hh_kin_sl3d[dstar_hh_index]>3 and "
        "hh_kin_vtx_prob[dstar_hh_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
        "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
        "hh_kin_mass[dstar_hh_index]>1.81 and hh_kin_mass[dstar_hh_index]<1.94 and "
        # "hh_kin_mass[dstar_mm_index]>1.65 and hh_kin_mass[dstar_mm_index]<2.4 and "
        "hh_had1_pdgId[dstar_hh_index] * hh_had2_pdgId[dstar_hh_index] == - 321 * 211"
    ),
    
    # "dstar_dzpipi_loose" : (
    #     "dstar_hh_index>=0 and "
    #     "hh_had1_pt[dstar_hh_index]>4 and hh_had2_pt[dstar_hh_index]>4 and "
    #     "hh_kin_alpha[dstar_hh_index]<0.1 and hh_kin_sl3d[dstar_hh_index]>3 and "
    #     "hh_kin_vtx_prob[dstar_hh_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
    #     "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
    #     "hh_kin_mass[dstar_hh_index]>1.81 and hh_kin_mass[dstar_hh_index]<1.94 and "
    #     "hh_had1_pdgId[dstar_hh_index] * hh_had2_pdgId[dstar_hh_index] == - 211 * 211"
    # ),
    "dstar_dzmm" : (
        "dstar_mm_index>=0 and "
        "mm_mu1_pt[dstar_mm_index]>4 and mm_mu2_pt[dstar_mm_index]>3 and "
        "mm_kin_alpha[dstar_mm_index]<0.1 and mm_kin_sl3d[dstar_mm_index]>3 and "
        "mm_kin_vtx_prob[dstar_mm_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
        # "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
        # "mm_kin_mass[dstar_mm_index]>1.76 and mm_kin_mass[dstar_mm_index]<1.94 and "
        "mm_kin_mass[dstar_mm_index]>1.5 and mm_kin_mass[dstar_mm_index]<2.4 and "
        "Muon_isGlobal[mm_mu1_index[dstar_mm_index]] and "
        "Muon_isGlobal[mm_mu2_index[dstar_mm_index]]"
        # "Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 and "
        # "Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45"
    ),
    # "dstar_dzmm_loose" : (
    #     "dstar_mm_index>=0 and "
    #     "mm_mu1_pt[dstar_mm_index]>4 and mm_mu2_pt[dstar_mm_index]>3 and "
    #     "mm_kin_alpha[dstar_mm_index]<0.1 and mm_kin_sl3d[dstar_mm_index]>3 and "
    #     "mm_kin_vtx_prob[dstar_mm_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
    #     "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
    #     "mm_kin_mass[dstar_mm_index]>1.65 and mm_kin_mass[dstar_mm_index]<2.4 and "
    #     "Muon_isGlobal[mm_mu1_index[dstar_mm_index]] and "
    #     "Muon_isGlobal[mm_mu2_index[dstar_mm_index]]"
    # ),
    
    # selection for cut-based analysis
    "ksmm": (
        "mm_mu1_index>=0 and mm_mu2_index>=0 and "
        "Muon_charge[mm_mu1_index] * Muon_charge[mm_mu2_index] < 0 and "
        "Muon_softMva[mm_mu1_index] > 0.45 and mm_mu1_pt>4 and "
        "Muon_softMva[mm_mu2_index] > 0.45 and mm_mu2_pt>4 and "
        # "abs(mm_kin_mass-0.45)<0.1 and "
        "abs(mm_kin_mass-0.50)<0.15 and "
        # "mm_kin_slxy>10 and mm_kin_lxy>1 and mm_kin_alpha<0.01 and "
        "mm_kin_slxy>3 and mm_kin_alpha<0.1 and "
        "HLT_DoubleMu4_3_LowMass"
    ),
    "kspipi": (
        "hh_had1_pdgId * hh_had2_pdgId == -211*211 and "
        "hh_had1_pt>4 and hh_had2_pt>4 and "
        "abs(hh_kin_mass-0.45) < 0.1 and "
        "hh_kin_slxy>10 and hh_kin_lxy>1 and hh_kin_alpha<0.001 and "
        "hh_kin_vtx_prob>0.01"
    ),
}

cuts["bkmm"] = (
    "mm_mu1_index[bkmm_mm_index]>=0 and mm_mu2_index[bkmm_mm_index]>=0 and "
    "Muon_charge[mm_mu1_index[bkmm_mm_index]] * Muon_charge[mm_mu2_index[bkmm_mm_index]] < 0 and "
    "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 and Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 and "
    "Muon_softMva[mm_mu1_index[bkmm_mm_index]] > 0.45 and Muon_softMva[mm_mu2_index[bkmm_mm_index]] > 0.45 and "
    "mm_kin_vtx_prob[bkmm_mm_index]>0.1 and " + cuts["bkmm_skim"]
)    


common_branches = "PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock"

tasks = [

    ############################
    # Skims
    ############################
    
    {
        "input_pattern":"BuToJpsiK",
        "processor":"Skimmer",
        # trigger based
        "cut" : "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && mm_kin_pt[bkmm_mm_index]>7.0 && mm_kin_alphaXY[bkmm_mm_index]<0.4 && mm_kin_vtx_prob[bkmm_mm_index]>0.1 && bkmm_jpsimc_vtx_prob>0.025 && mm_kin_sl3d[bkmm_mm_index]>4",
        "name":"bkmm",
        "type":"NanoAOD-skims",
        "files_per_job":10
    },
    {
        "input_pattern":"Charmonium",
        "processor":"Skimmer",
        # trigger based
        "cut" : "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && mm_kin_pt[bkmm_mm_index]>7.0 && mm_kin_alphaXY[bkmm_mm_index]<0.4 && mm_kin_vtx_prob[bkmm_mm_index]>0.1 && bkmm_jpsimc_vtx_prob>0.025 && mm_kin_sl3d[bkmm_mm_index]>4",
        "name":"bkmm",
        "type":"NanoAOD-skims",
        "files_per_job":100
    },

    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"SimpleSkimmer",
        "cut": cuts["bkmm_skim"],
        "name":"bkmm",
        "type":"Skims",
        "keep": "^(mm_.*|HLT_DoubleMu4_3_LowMass|nmm|bkmm_.*|nbkmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" +
        common_branches + ")$",
        "files_per_job":200
    },


    # {
    #     "input_pattern":"DoubleMuon",
    #     "processor":"Skimmer",
    #     "cut":"HLT_Mu8 || HLT_Mu17",
    #     "name":"trig",
    #     "type":"NanoAOD-skims",
    #     "files_per_job":50
    # },
    {
        "input_pattern":"Charmonium.Run2016",
        "processor":"Skimmer",
        "cut":"HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
        "name":"trig",
        "type":"NanoAOD-skims",
        "files_per_job":25
    },
    {
        "input_pattern":"DoubleMuon.Run2016H",
        "processor":"Skimmer",
        "cut":"HLT_DoubleMu0",
        "name":"trig",
        "type":"NanoAOD-skims",
        "files_per_job":25
    },
    {
        "input_pattern":"Charmonium.Run2018|Charmonium.Run2017",
        "processor":"Skimmer",
        "cut":"HLT_Dimuon0_LowMass_L1_0er1p5",
        "name":"trig",
        "type":"NanoAOD-skims",
        "files_per_job":25
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM",
        # "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"Skimmer",
        "cut":"ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
        "name":"ks",
        "type":"NanoAOD-skims",
        "files_per_job":100
    },
    {
        "input_pattern":"DoubleElectron|EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM",
        # "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"SimpleSkimmer",
        "cut":"ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
        "name":"ks",
        "keep": "^(ks_.*|nks|Muon_.*|nMuon|MuonId_.*|nMuonId|run|event|luminosityBlock|HLT_Ele30_WPTight_Gsf)$",
        "type":"Skims",
        "files_per_job":100,
        "comments":"Dmm and Muon paper setup (2024-01-01)"
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM",
        # "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"Skimmer",
        "cut":"lambda_kin_slxy>3 && lambda_kin_sipPV<3 && lambda_proton_sip>2 && lambda_pion_sip>2",
        "name":"lambda",
        "type":"NanoAOD-skims",
        "files_per_job":100
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM",
        # "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"Skimmer",
        "cut":"((phi_trk1_pt > 4 && phi_trk2_pt > 3.0) || ( phi_trk2_pt > 4 && phi_trk1_pt > 3.0)) &&  phi_kin_vtx_prob>0.3 && phi_kin_sipPV<1  && phi_doca<0.004 && phi_kin_lxy<4",
        "name":"phi",
        "type":"NanoAOD-skims",
        "files_per_job":100
    },
    {
        "input_pattern":"DoubleElectron|EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM",
        "processor":"SimpleSkimmer",
        "cut":"((phi_trk1_pt > 4 && phi_trk2_pt > 3.0) || ( phi_trk2_pt > 4 && phi_trk1_pt > 3.0)) &&  phi_kin_vtx_prob>0.3 && phi_kin_sipPV<1  && phi_doca<0.004 && phi_kin_lxy<4",
        "name":"phi",
        "type":"Skims",
        "keep": "^(phi_.*|nphi|Muon_.*|nMuon|MuonId_.*|nMuonId|" + common_branches + ")$",
        "files_per_job":200
    },
    # {
    #     "input_pattern":"DoubleMuon|MINIAODSIM",
    #     "processor":"SimpleSkimmer",
    #     "cut":"d0_dstar_pion_pt > 0",
    #     "name":"dstar",
    #     "type":"Skims",
    #     "keep": "^(d0_.*|nd0|mm_.*|nmm|Muon_.*|nMuon|run|event|luminosityBlock)$",
    #     "files_per_job":100
    # },
    {
        "input_pattern":"DoubleMuon|ZeroBias|MINIAODSIM",
        "processor":"SimpleSkimmer",
        "cut":"dstar_dm_pv > 0",
        "name":"dstar",
        "type":"Skims",
        "keep": "^(dstar_.*|ndstar|d0_.*|nd0|mm_.*|nmm|Muon_.*|nMuon|hh_.*|nhh|run|event|luminosityBlock)$",
        "files_per_job":50
    },
    {
        "input_pattern":"MINIAODSIM",
        "processor":"SimpleSkimmer",
        "cut":"dstar_hh_index >= 0 && dstar_gen_pdgId!=0",
        "name":"dstar_hh",
        "type":"Skims",
        "keep": "^(dstar_.*|ndstar|d0_.*|nd0|mm_.*|nmm|Muon_.*|nMuon|hh_.*|nhh|" + common_branches + ")$",
        "files_per_job":200
    },
    {
        "input_pattern":"^(DstarToD0Pi_D0To2Pi_SoftQCD|InclusiveDileptonMinBias|DY|B)",
        "processor":"SimpleSkimmer",
        "cut":"abs(hh_gen_had1_pdgId)==211&&abs(hh_gen_had2_pdgId)==211&&abs(hh_gen_pdgId)==421&&hh_had1_pt>4&&hh_had2_pt>4",
        "name":"dzpipi",
        "type":"Skims",
        "keep": "^(dstar_.*|ndstar|hh_.*|nhh|" + common_branches + ")$",
        "files_per_job":1000
    },
    {
        "input_pattern":"ZeroBias",
        "processor":"SimpleSkimmer",
        "cut":"hh_had1_pdgId*hh_had2_pdgId==-211*211&&hh_had1_pt>4&&hh_had2_pt>4",
        "name":"dzpipi",
        "type":"Skims",
        "keep": "^(dstar_.*|ndstar|hh_.*|nhh|" + common_branches + ")$",
        "files_per_job":50
    },
    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_mass>1.81 && mm_kin_mass<1.94 && mm_mu1_pt>4 && mm_mu2_pt>4",
        "name":"dzmm",
        "type":"Skims",
        "keep": "^(dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|HLT_.*" + common_branches + ")$",
        "files_per_job":50
    },
    {
        "input_pattern":"Inclusive",
        "processor":"SimpleSkimmer",
        "cut":"abs(mm_gen_pdgId)==421",
        "name":"dzkpimm",
        "type":"Skims",
        "keep": "^(GenPart_.*|nGenPart|dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|HLT_.*" + common_branches + ")$",
        "files_per_job":200
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton|MINIAODSIM",
        # "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"Skimmer",
        "cut":"phi_ds_pion_pt>0",
        "name":"ds",
        "type":"NanoAOD-skims",
        "files_per_job":100
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"Skimmer",
        "cut":"nmm>0",
        "name":"mm",
        "type":"NanoAOD-skims",
        "files_per_job":100
    },
    {
        "input_pattern":"InclusiveDileptonMinBias",
        "processor":"SimpleSkimmer",
        "cut":"mm_mass > 0",
        "name":"mm",
        "type":"Skims",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + common_branches + ")$",
        "files_per_job":100
    },
    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"SimpleSkimmer",
        "cut":"mm_mass > 0",
        "name":"mm",
        "type":"Skims",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + common_branches + ")$",
        "files_per_job":20
    },
    {
        "input_pattern":"ParkingBPH|InclusiveDileptonMinBias",
        "processor":"SimpleSkimmer",
        "cut":"abs(em_kin_mass-5.4)<1.0",
        "name":"em",
        "type":"Skims",
        "keep": "^(em_.*|nem|HLT_Mu.*_IP.*|Muon_.*|nMuon|MuonId_.*|nMuonId|Electron_.*|nElectron|npvs|pvs_.*|" + common_branches + ")$",
        "files_per_job":250
    },
    {
        "input_pattern":"InclusiveDileptonMinBias",
        "processor":"SimpleSkimmer",
        "cut":"Muon_pt > 0",
        "name":"muons",
        "type":"Skims",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + common_branches + ")$",
        "files_per_job":100
    },
    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_vtx_prob>0.01&&mm_kin_slxy>10&&mm_kin_lxy>1&&mm_kin_alpha<0.01&&mm_mu1_pdgId==-mm_mu2_pdgId&&abs(mm_kin_mass-0.5)<0.2",
        "name":"ksmm",
        "type":"Skims",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|" + common_branches + ")$",
        "files_per_job":100
    },
    {
        "input_pattern":"InclusiveDileptonMinBias",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_vtx_prob>0.01&&mm_kin_slxy>10&&mm_kin_lxy>1&&mm_kin_alpha<0.01&&mm_mu1_pdgId==-mm_mu2_pdgId&&abs(mm_kin_mass-0.5)<0.2",
        "name":"ksmm",
        "type":"Skims",
        "keep": "^(GenPart_.*|nGenPart|mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|" + common_branches + ")$",
        "files_per_job":100
    },
    {
        "input_pattern":"K0sToMuMu",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_vtx_prob>0.01&&mm_kin_slxy>10&&mm_kin_lxy>1&&mm_kin_alpha<0.01&&mm_mu1_pdgId==-mm_mu2_pdgId&&abs(mm_kin_mass-0.5)<0.2",
        "name":"ksmm",
        "type":"Skims",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|" + common_branches + ")$",
        "files_per_job":10
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"Skimmer",
        "cut":"mm_mu1_index>=0&&mm_mu2_index>=0&&mm_kin_vtx_prob>0.025&&(Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index])==-1",
        "name":"mm-vtx",
        "type":"NanoAOD-skims",
        "files_per_job":20
    },
    {
        "input_pattern":"LambdaB|B.ToPiPi|B.ToKPi|B.ToKK",
        "processor":"Skimmer",
        "cut":"abs(mm_kin_mu1eta)<1.4 && "\
        "mm_kin_mu1pt>4 && "\
        "abs(mm_kin_mu2eta)<1.4 && "\
        "mm_kin_mu2pt>4 && "\
        "abs(mm_kin_mass-5.4)<0.5 && "\
        "mm_kin_sl3d>6 && "\
        "mm_kin_pt>5.0 && mm_kin_vtx_prob>0.025 &&"\
        "mm_mu1_index>=0 && mm_mu2_index>=0 && HLT_DoubleMu4_3_Bs",
        "name":"bmm",
        "type":"NanoAOD-skims",
        "files_per_job":20
    },
    

    ############################
    # Flat Ntuples
    ############################
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton|SingleMuon|DoubleMuon",
        "processor":"FlatNtupleForTrigInfo",
        "tree_name":"mm",
        "name":"trig-info",
        "type":"FlatNtuples",
        "files_per_job":25
    },
    {
        "input_pattern":"BsToMuMu_BMuonFilter|BuToJpsiK_BMuonFilter",
        "processor":"FlatNtupleForTrigInfo",
        "tree_name":"mm",
        "name":"trig-info",
        "type":"FlatNtuples",
        "files_per_job":10
    },
    {
        "input_pattern":"SingleMuon",
        "processor":"FlatNtupleForTrigEfficiency",
        "tree_name": "muons",
        "type":"FlatNtuples",
        "require_muon_tag" : True, 
        "tag_triggers" : ["HLT_IsoMu24", "HLT_IsoMu27"],
        "cut" : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        "name":"trig-eff",
        "files_per_job":25
    },
    {
        "input_pattern":"EGamma|DoubleEG|SingleElectron|SinglePhoton",
        "processor":"FlatNtupleForTrigEfficiency",
        "tree_name": "muons",
        "type":"FlatNtuples",
        "require_muon_tag" : False, 
        "tag_triggers" : [
            "L1_Mu6_DoubleEG10er2p5", "L1_Mu6_DoubleEG17er2p5",
            "L1_Mu6_DoubleEG12er2p5", "L1_Mu6_DoubleEG15er2p5",
            "L1_SingleMu3", "L1_SingleMu5", "L1_SingleMu7",
            "L1_Mu6_DoubleEG10", "L1_Mu6_DoubleEG17"
        ],
        "cut" : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        "name":"trig-eff",
        "files_per_job":25
    },
    {
        "input_pattern":"DoubleMuon",
        "processor":"FlatNtupleForTrigEfficiency",
        "tree_name": "muons",
        "type":"FlatNtuples",
        "require_muon_tag" : True, 
        "tag_triggers" : ["HLT_Mu8", "HLT_Mu17"],
        "cut" : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        "name":"trig-eff",
        "files_per_job":25
    },
    {
        "input_pattern":"BsToMuMu_BMuonFilter",
        "processor":"FlatNtupleForTrigEfficiency",
        "tree_name": "muons",
        "type":"FlatNtuples",
        "require_muon_tag" : True, 
        "tag_triggers" : ["HLT_Mu8"],
        "cut" : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20",
        "name":"trig-eff",
        "files_per_job":20
    },

    ############### bmm_mva_jpsik ###############
    # {
    #     "input_pattern":"BuToJpsiK",
    #     "processor":"FlatNtupleForBmmMvaJpsiK",
    #     "signal_only": True,
    #     "tree_name": "mva",
    #     "name":"bmm_mva_jpsik",
    #     "type":"FlatNtuples",
    #     "files_per_job":20
    # },
    # {
    #     "input_pattern":"Charmonium",
    #     "processor":"FlatNtupleForBmmMvaJpsiK",
    #     "signal_only": False,
    #     "tree_name": "mva",
    #     "name":"bmm_mva_jpsik",
    #     "type":"FlatNtuples",
    #     "files_per_job":20
    # },
    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForBmmMvaJpsiK",
        "signal_only": False,
        "tree_name": "mva",
        "name":"bmm_mva_jpsik",
        "type":"FlatNtuples",
        "files_per_job":10
    },
    
    ############### bmm_mva_jpsik_loose_vtx ###############
    {
        "input_pattern":"BuToJpsiK",
        "processor":"FlatNtupleForBmmMvaJpsiK",
        "signal_only": True,
        "loose_vtx": True,
        "tree_name": "mva",
        "name":"bmm_mva_jpsik_loose_vtx",
        "type":"FlatNtuples",
        "files_per_job":20
    },
    {
        "input_pattern":"Charmonium",
        "processor":"FlatNtupleForBmmMvaJpsiK",
        "signal_only": False,
        "loose_vtx": True,
        "tree_name": "mva",
        "name":"bmm_mva_jpsik_loose_vtx",
        "type":"FlatNtuples",
        "files_per_job":20
    },

    ################ fit ################
    {
        "input_pattern":"BsToMuMu",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bsmmMc",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
        
    },
    {
        "input_pattern":"BdToMuMu",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bmmMc",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
        
    },
    {
        "input_pattern":"Charmonium",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bmmData",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ### Btohh exclusive

    {
        "input_pattern":"BsToKK_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bskkMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BsToKPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bskpiMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BsToPiPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bspipiMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BdToKK_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bdkkMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BdToKPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bdkpiMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BdToPiPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bdpipiMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ### Btohh inclusive
    
    {
        "input_pattern":"BTohh",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "btohhMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ###
    
    {
        "input_pattern":"LambdaBToPPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "lbppiMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"LambdaBToPK_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "lbpkMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BsToKMuNu_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bskmunuMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BdToPiMuNu_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bdpimunuMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"LambdaBToPMuNu_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "lbpmunuMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BdToMuMuPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bdpimumuMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },
    {
        "input_pattern":"BuToMuMuPi_",
        "processor":"FlatNtupleForMLFit",
        "name":"fit",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bupimumuMcBg",
        "blind" : False,
        "cut" : cuts["fit"],
        "final_state" : "mm",
        # "best_candidate": "mm_kin_pt",
        "best_candidate": "",
    },

    ################ fit-dmm ################
    # {
    #     "input_pattern":"BsToMuMu",
    #     "processor":"FlatNtupleForMLFit",
    #     "name":"fit",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
    #     "tree_name" : "bsmmMc",
    #     "blind" : False,
    #     "cut" : cuts["fit"],
    #     "final_state" : "mm",
    #     # "best_candidate": "mm_kin_pt",
    #     "best_candidate": "",
        
    # },
    # {
    #     "input_pattern":"BdToMuMu",
    #     "processor":"FlatNtupleForMLFit",
    #     "name":"fit",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
    #     "tree_name" : "bmmMc",
    #     "blind" : False,
    #     "cut" : cuts["fit"],
    #     "final_state" : "mm",
    #     # "best_candidate": "mm_kin_pt",
    #     "best_candidate": "",
        
    # },
    # {
    #     "input_pattern":"Charmonium",
    #     "processor":"FlatNtupleForMLFit",
    #     "name":"fit",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
    #     "tree_name" : "bmmData",
    #     "blind" : False,
    #     "cut" : cuts["fit"],
    #     "final_state" : "mm",
    #     # "best_candidate": "mm_kin_pt",
    #     "best_candidate": "",
    # },
    
    ################ fit-em ################
    {
        "input_pattern":"BsToEMu",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-em",
        "type":"FlatNtuples",
        "files_per_job":6,
        "tree_name" : "bsemMc",
        "blind" : False,
        "cut" : cuts["fit-em"],
        "final_state" : "em",
        "best_candidate": "",
        
    },
    {
        "input_pattern":"BdToEMu",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-em",
        "type":"FlatNtuples",
        "files_per_job":6,
        "tree_name" : "bdemMc",
        "blind" : False,
        "cut" : cuts["fit-em"],
        "final_state" : "em",
        "best_candidate": "",
        
    },
    {
        "input_pattern":"Skim.*ParkingBPH",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-em",
        "type":"FlatNtuples",
        "files_per_job":50,
        "tree_name" : "bemData",
        "blind" : False,
        "cut" : cuts["fit-em"],
        "final_state" : "em",
        "best_candidate": "",
    },
    
    ################ fit-bkmm ################
    {
        "input_pattern":"BuToJpsiK",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-bkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bupsikMc",
        "blind" : False,
        "cut" : cuts["fit-bkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },
    {
        "input_pattern":"BuToJpsiPi",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-bkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bupsipiMc",
        "blind" : False,
        "cut" : cuts["fit-bkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },
    {
        "input_pattern":"Charmonium",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-bkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bupsikData",
        "blind" : False,
        "cut" : cuts["fit-bkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },
    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForMLFit",
        "name":"bkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bupsikData",
        "blind" : False,
        "cut" : cuts["bkmm"],
        "triggers": ["HLT_DoubleMu4_3_LowMass"],
        "final_state" : "bkmm",
        "best_candidate": "",
    },

    ################ fit-bkkmm ################
    {
        "input_pattern":"BsToJPsiPhi",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-bkkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bspsiphiMc",
        "blind" : False,
        "cut" : cuts["fit-bkkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "best_candidate": "",
    },
    {
        "input_pattern":"BdToJpsiKstar",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-bkkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bdpsikstarMc",
        "blind" : False,
        "cut" : cuts["fit-bkkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "best_candidate": "",
    },
    {
        "input_pattern":"Charmonium",
        "processor":"FlatNtupleForMLFit",
        "name":"fit-bkkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bspsiphiData",
        "blind" : False,
        "cut" : cuts["fit-bkkmm"], "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "best_candidate": "",
    },

    {
        "input_pattern":"ScoutingPFRun3|ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForMLFit",
        "name":"bkkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bspsiphiData",
        "blind" : False,
        "cut" : cuts["bkkmm"],
        # "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkkmm",
        "pre-selection":"abs(bkkmm_kk_mass-1.02)<0.01&&bkkmm_jpsikk_sl3d>5",
        "pre-selection-keep":"^(mm_.*|nmm|bkkmm_.*|nbkkmm|HLT_*|" + common_branches + ")$",
    },

    {
        "input_pattern":"ScoutingPFRun3|ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForMLFit",
        "name":"bkmm",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "bupsikData",
        "blind" : False,
        "cut" : cuts["bkmm-validation"],
        # "triggers": ["HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_3_Jpsi_Displaced"],
        "final_state" : "bkmm",
        "pre-selection":"bkmm_jpsimc_sl3d>5 && abs(bkmm_jpsimc_mass-5.4)<0.5 && abs(bkmm_jpsimc_alpha)<0.01",
        "pre-selection-keep":"^(mm_.*|nmm|bkmm_.*|nbkmm|HLT_*|" + common_branches + ")$",
    },
    
    ################ Dimuon ################
    # {
    #     "input_pattern":"BsToMuMu",
    #     "processor":"FlatNtupleForMLFit",
    #     "name":"fit",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
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
    #     "input_pattern":"BdToMuMu",
    #     "processor":"FlatNtupleForMLFit",
    #     "name":"fit",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
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
        "input_pattern":"Jpsi",
        "processor":"FlatNtupleForMLFit",
        "name":"dimuon",
        "type":"FlatNtuples",
        "files_per_job":20,
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
        "input_pattern":"YnSToMuMu",
        "processor":"FlatNtupleForMLFit",
        "name":"dimuon",
        "type":"FlatNtuples",
        "files_per_job":20,
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
        "input_pattern":"Charmonium",
        "processor":"FlatNtupleForMLFit",
        "name":"dimuon",
        "type":"FlatNtuples",
        "files_per_job":20,
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
        "input_pattern":"MuOnia",
        "processor":"FlatNtupleForMLFit",
        "name":"dimuon",
        "type":"FlatNtuples",
        "files_per_job":20,
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
        "input_pattern":"InclusiveDileptonMinBias", # "QCD_Pt.*?_MuEnriched|Charmonium|BuToJpsiK",
        "processor":"FlatNtupleForMuonMVA",
        "name":"muon_mva",
        "type":"FlatNtuples",
        "files_per_job":10,
        "tree_name" : "muons",
    },

    ################ bmm_mva ################
    {
        "input_pattern":"Charmonium",
        "processor":"FlatNtupleForBmmMva",
        "signal_only": False,
        "tree_name": "mva",
        "blind": True,
        "name":"bmm_mva",
        "type":"FlatNtuples",
        "files_per_job":20
    },
    {
        "input_pattern":"BsToMuMu_BMuonFilter",
        "processor":"FlatNtupleForBmmMva",
        "signal_only": True,
        "tree_name": "mva",
        "blind": False,
        "name":"bmm_mva",
        "type":"FlatNtuples",
        "files_per_job":20
    },
    ################ dstar ##################
    {
        "input_pattern":"DstarToD0Pi_D0To2Pi_PiFilter",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzpipi",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "dzpipiMC",
        "blind" : False,
        "cut" : cuts["dstar_dzpipi"],
        "final_state" : "dzpipi",
        "best_candidate": "",
    },
    {
        "input_pattern":"DstarToD0Pi_D0To2Mu_MuFilter",
        "processor":"FlatNtupleForDstarFit",
        "name":"dstar",
        "type":"FlatNtuples",
        "files_per_job":10,
        "tree_name" : "dzmmMC",
        "blind" : False,
        "cut" : cuts["dstar_dzmm"],
        "final_state" : "dzmm",
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        "best_candidate": "",
    },
    {
        "input_pattern":"DstarToD0Pi_D0To2Pi_PiToMu",
        "processor":"FlatNtupleForDstarFit",
        "name":"dstar",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "dzpipimmMC",
        "blind" : False,
        "cut" : cuts["dstar_dzmm"],
        "final_state" : "dzmm",
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        "best_candidate": "",
    },
    {
        "input_pattern":"DstarToD0Pi_D0ToKMuNu",
        "processor":"FlatNtupleForDstarFit",
        "name":"dstar",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "dzkmunuMC",
        "blind" : False,
        "cut" : cuts["dstar_dzmm"],
        "final_state" : "dzmm",
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        "best_candidate": "",
    },
    {
        "input_pattern":"DstarToD0Pi_D0ToPiMuNu",
        "processor":"FlatNtupleForDstarFit",
        "name":"dstar",
        "type":"FlatNtuples",
        "files_per_job":20,
        "tree_name" : "dzpimunuMC",
        "blind" : False,
        "cut" : cuts["dstar_dzmm"],
        "final_state" : "dzmm",
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        "best_candidate": "",
    },
    {
        "input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForDstarFit",
        "name":"dstar",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzmmData",
        "blind" : False,
        "cut" : cuts["dstar_dzmm"],
        "final_state" : "dzmm",
        "best_candidate": "",
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        # "pre-selection":"dstar_mm_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection":"dstar_mm_index>=0",
        "pre-selection-keep":"^(dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
    },
    {
        "input_pattern":"ZeroBias",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzpipi",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzpipiData",
        "blind" : False,
        "cut" : cuts["dstar_dzpipi"],
        "final_state" : "dzpipi",
        "best_candidate": "",
        "triggers":["HLT_ZeroBias"],
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|HLT_ZeroBias|" + common_branches + ")$",
    },
    {
        "input_pattern":"ZeroBias",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzpipi_otherZB",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzpipiData",
        "blind" : False,
        "cut" : cuts["dstar_dzpipi"],
        "final_state" : "dzpipi",
        "best_candidate": "",
        "triggers":["HLT_ZeroBias_FirstBXAfterTrain", "HLT_ZeroBias_FirstCollisionAfterAbortGap",
                    "HLT_ZeroBias_FirstCollisionInTrain", "HLT_ZeroBias_IsolatedBunches",
                    "HLT_ZeroBias_LastCollisionInTrain"],
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|HLT_ZeroBias.*|" + common_branches + ")$",
    },
    {
        "input_pattern":"ParkingDoubleElectronLowMass|ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzpipi",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzpipiData",
        "blind" : False,
        "cut" : cuts["dstar_dzpipi"],
        "final_state" : "dzpipi",
        "best_candidate": "",
        "triggers":[],
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|" + common_branches + ")$",
    },
    {
        "input_pattern":"ZeroBias",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzkpi",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzkpiData",
        "blind" : False,
        "cut" : cuts["dstar_dzkpi"],
        "final_state" : "dzkpi",
        "best_candidate": "",
        "triggers":["HLT_ZeroBias"],
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|HLT_ZeroBias|" + common_branches + ")$",
    },
    {
        "input_pattern":"ZeroBias",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzkpi_otherZB",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzkpiData",
        "blind" : False,
        "cut" : cuts["dstar_dzkpi"],
        "final_state" : "dzkpi",
        "best_candidate": "",
        "triggers":["HLT_ZeroBias_FirstBXAfterTrain", "HLT_ZeroBias_FirstCollisionAfterAbortGap",
                    "HLT_ZeroBias_FirstCollisionInTrain", "HLT_ZeroBias_IsolatedBunches",
                    "HLT_ZeroBias_LastCollisionInTrain"],
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|HLT_ZeroBias.*|" + common_branches + ")$",
    },
    {
        "input_pattern":"ParkingDoubleElectronLowMass|ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzkpi",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzkpiData",
        "blind" : False,
        "cut" : cuts["dstar_dzkpi"],
        "final_state" : "dzkpi",
        "best_candidate": "",
        "triggers":[],
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|" + common_branches + ")$",
    },
    {
        "input_pattern":"Inclusive|DStartoD0Pi_D0toKPi|DstarToD0Pi_D0ToKPi",
        "processor":"FlatNtupleForDstarFit",
        "name":"dzkpi",
        "type":"FlatNtuples",
        "files_per_job":200,
        "tree_name" : "dzkpiMC",
        "blind" : False,
        "cut" : cuts["dstar_dzkpi"] + " and abs(dstar_gen_pdgId)==413 and hh_gen_had1_pdgId[dstar_hh_index]*hh_gen_had2_pdgId[dstar_hh_index]==-321*211 and abs(hh_gen_pdgId[dstar_hh_index])==421",
        "final_state" : "dzkpi",
        "best_candidate": "",
        "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|" + common_branches + ")$",
    },
    {
        "input_pattern":"Inclusive",
        "processor":"FlatNtupleForDstarFit",
        "name":"dstar",
        "type":"FlatNtuples",
        "files_per_job":100,
        "tree_name" : "dzmmMcBkg",
        "blind" : False,
        "cut" : cuts["dstar_dzmm"],
        "final_state" : "dzmm",
        "best_candidate": "",
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        # "pre-selection":"dstar_mm_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection":"dstar_mm_index>=0",
        "pre-selection-keep":"^(GenPart_.*|nGenPart|dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
    },
    ################ dstar loose ##################
    # {
    #     "input_pattern":"DstarToD0Pi_D0To2Mu_MuFilter",
    #     "processor":"FlatNtupleForDstarFit",
    #     "name":"dstar_loose",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
    #     "tree_name" : "dzmmMC",
    #     "blind" : False,
    #     "cut" : cuts["dstar_dzmm_loose"],
    #     "final_state" : "dzmm",
    #     "triggers":["HLT_DoubleMu4_3_LowMass"],
    #     "best_candidate": "",
    # },
    # {
    #     "input_pattern":"DstarToD0Pi_D0To2Pi_PiToMu",
    #     "processor":"FlatNtupleForDstarFit",
    #     "name":"dstar_loose",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
    #     "tree_name" : "dzpipimmMC",
    #     "blind" : False,
    #     "cut" : cuts["dstar_dzmm_loose"],
    #     "final_state" : "dzmm",
    #     "triggers":["HLT_DoubleMu4_3_LowMass"],
    #     "best_candidate": "",
    # },
    # {
    #     "input_pattern":"ParkingDoubleMuonLowMass",
    #     "processor":"FlatNtupleForDstarFit",
    #     "name":"dstar_loose",
    #     "type":"FlatNtuples",
    #     "files_per_job":100,
    #     "tree_name" : "dzmmData",
    #     "blind" : False,
    #     "cut" : cuts["dstar_dzmm_loose"],
    #     "final_state" : "dzmm",
    #     "best_candidate": "",
    #     "triggers":["HLT_DoubleMu4_3_LowMass"],
    #     "pre-selection":"dstar_mm_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
    #     "pre-selection-keep":"^(dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
    # },
    # {
    #     "input_pattern":"DstarToD0Pi_D0To2Pi_PiFilter",
    #     "processor":"FlatNtupleForDstarFit",
    #     "name":"dstar_loose",
    #     "type":"FlatNtuples",
    #     "files_per_job":20,
    #     "tree_name" : "dzpipiMC",
    #     "blind" : False,
    #     "cut" : cuts["dstar_dzpipi_loose"],
    #     "final_state" : "dzpipi",
    #     "best_candidate": "",
    # },
    # {
    #     "input_pattern":"ZeroBias",
    #     "processor":"FlatNtupleForDstarFit",
    #     "name":"dstar_loose",
    #     "type":"FlatNtuples",
    #     "files_per_job":100,
    #     "tree_name" : "dzpipiData",
    #     "blind" : False,
    #     "cut" : cuts["dstar_dzpipi_loose"],
    #     "final_state" : "dzpipi",
    #     "best_candidate": "",
    #     "triggers":["HLT_ZeroBias"],
    #     "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
    #     "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|HLT_ZeroBias|" + common_branches + ")$",
    # },
    # {
    #     "input_pattern":"ParkingDoubleElectronLowMass|ParkingDoubleMuonLowMass",
    #     "processor":"FlatNtupleForDstarFit",
    #     "name":"dzpipi_loose",
    #     "type":"FlatNtuples",
    #     "files_per_job":100,
    #     "tree_name" : "dzpipiData",
    #     "blind" : False,
    #     "cut" : cuts["dstar_dzpipi_loose"],
    #     "final_state" : "dzpipi",
    #     "best_candidate": "",
    #     "triggers":[],
    #     "pre-selection":"dstar_hh_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
    #     "pre-selection-keep":"^(dstar_.*|ndstar|hh_.*|nhh|" + common_branches + ")$",
    # },
    
    ###################################################################
    #      Other
    ###################################################################
    # {
    #     "input_pattern":"InclusiveDileptonMinBias",
    #     "processor":"SimpleSkimmer",
    #     "cut":"dstar_mm_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
    #     "name":"dstar_mm",
    #     "type":"Skims",
    #     "files_per_job":50,
    #     "keep": "^(dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|" + common_branches + ")$",
    # },
    {
        "input_pattern":"Muon.*MINIAOD$",
        "processor":"SimpleSkimmer",
        "cut":"HLT_Mu4_L1DoubleMu || HLT_Mu3_PFJet40 || HLT_Mu8",
        "name":"trig",
        "candidate_loop":False,
        "type":"Skims",
        "files_per_job":100,
        "keep": "^(dstar_.*|ndstar|mm_.*|nmm|Muon_.*|nMuon|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|HLT_Mu3_PFJet40|HLT_Mu8|nJet|Jet_pt|Jet_eta|Jet_phi|" + common_branches + ")$",
    },
    {
        "input_pattern":"InclusiveDileptonMinBias",
        "processor":"SimpleSkimmer",
        "cut":"mm_mva>0.9",
        "name":"mm_mva0p9",
        "candidate_loop":True,
        "type":"Skims",
        "files_per_job":500,
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|HLT_Mu3_PFJet40|HLT_Mu8|HLT_Mu12_IP6.*|" + common_branches + ")$",
    },
    {
        "input_pattern":"ZeroBias|DoubleElectron",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_vtx_prob>0.1",
        "name":"mm_vtx",
        "candidate_loop":True,
        "type":"Skims",
        "files_per_job":100,
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|bkmm_.*|nbkmm|bkkmm_.*|nbkkmm|" + common_branches + ")$",
    },
    # {
    #     "input_pattern":"InclusiveDileptonMinBias",
    #     "processor":"SimpleSkimmer",
    #     "cut":"mm_kin_mass<2",
    #     "name":"tau3mu",
    #     "candidate_loop":True,
    #     "type":"Skims",
    #     "files_per_job":100,
    #     "keep": "^(mm_.*|nmm|MuonId_.*|nMuonId|Muon_.*|nMuon|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|HLT_Mu3_PFJet40|HLT_Mu8|HLT_Mu12_IP6.*|" + common_branches + ")$",
    # },
    {
        "input_pattern":"DstoTau_Tauto3Mu_3MuFilter",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_mass<2",
        "name":"tau3mu",
        "candidate_loop":True,
        "type":"Skims",
        "files_per_job":1,
        "keep": "^(mm_.*|nmm|MuonId_.*|nMuonId|Muon_.*|nMuon|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|HLT_Mu3_PFJet40|HLT_Mu8|HLT_Mu12_IP6.*|" + common_branches + ")$",
    },
    {
        "input_pattern":"ParkingDoubleMuonLowMass0.*Run2022C",
        "processor":"SimpleSkimmer",
        "cut":"mm_kin_mass<2",
        "name":"tau3mu",
        "candidate_loop":True,
        "type":"Skims",
        "files_per_job":20,
        "keep": "^(mm_.*|nmm|MuonId_.*|nMuonId|Muon_.*|nMuon|HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|HLT_Mu3_PFJet40|HLT_Mu8|" + common_branches + ")$",
    },
    ##############################################
    #                Ks to mu mu
    ##############################################
    {
     	"input_pattern":"K0sToMuMu",
        "processor":"FlatNtupleForMLFit",
        "name":"ksmm",
	"type":"FlatNtuples",
	"files_per_job":50,
        "tree_name" : "ksmmMc",
        "blind" : False,
        "cut" : cuts["ksmm"],
        "final_state" : "mm",
        "mm_extra_info": True,
    },
    {
     	"input_pattern":"ParkingDoubleMuonLowMass",
        "processor":"FlatNtupleForMLFit",
        "name":"ksmm",
	"type":"FlatNtuples",
	"files_per_job":100,
        "tree_name" : "ksmmData",
        "blind" : False,
        "cut" : cuts["ksmm"],
        "final_state" : "mm",
        "pre-selection":"abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01 && mm_kin_slxy>3 && mm_kin_alpha<0.1", 
        "pre-selection-keep":"^(mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
        "mm_extra_info": True,
    },
    {
     	"input_pattern":"InclusiveDileptonMinBias",
        "processor":"FlatNtupleForMLFit",
        "name":"ksmm",
	"type":"FlatNtuples",
	"files_per_job":100,
        "tree_name" : "ksmmMcBkg",
        "blind" : False,
        "cut" : cuts["ksmm"],
        "final_state" : "mm",
        "pre-selection":"abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01 && mm_kin_slxy>3 && mm_kin_alpha<0.1", 
        "pre-selection-keep":"^(mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
        "mm_extra_info": True,
    },
    {
     	"input_pattern":"K0sTo2PiTo2Mu_K0sFilter_MuFilter",
        "processor":"FlatNtupleForMLFit",
        "name":"ksmm",
	"type":"FlatNtuples",
	"files_per_job":10,
        "tree_name" : "kspipimmMc",
        "blind" : False,
        "cut" : cuts["ksmm"],
        "final_state" : "mm",
        "pre-selection":"abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01 && mm_kin_slxy>3 && mm_kin_alpha<0.1", 
        "pre-selection-keep":"^(mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
        "mm_extra_info": True,
    },
    {
     	"input_pattern":"InclusiveDileptonMinBias",
        "processor":"FlatNtupleForMLFit",
        "name":"kspipi",
	"type":"FlatNtuples",
	"files_per_job":200,
        "tree_name" : "kspipiMc",
        "blind" : False,
        "cut" : cuts["kspipi"],
        "final_state" : "hh",
        "pre-selection":"hh_gen_pdgId==310", 
        "pre-selection-keep":"^(hh_.*|nhh|" + common_branches + ")$",
    },
    {
     	"input_pattern":"ZeroBias",
        "processor":"FlatNtupleForMLFit",
        "name":"kspipi",
	"type":"FlatNtuples",
	"files_per_job":200,
        "tree_name" : "kspipiData",
        "blind" : False,
        "cut" : cuts["kspipi"],
        "triggers": ["HLT_ZeroBias"],
        "final_state" : "hh",
        "pre-selection":"abs(hh_kin_mass-0.45)<0.1 && hh_kin_slxy>10 && hh_kin_vtx_prob>0.01", 
        "pre-selection-keep":"^(hh_.*|nhh|HLT_ZeroBias|" + common_branches + ")$",
    },

]
