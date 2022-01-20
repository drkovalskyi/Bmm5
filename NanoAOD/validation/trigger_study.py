"""Trigger Properties

The script performs a study of trigger properties using MC
simulations. The main goal is to find the primary factors affecting
the trigger efficiency.

"""
import os, sys, json
import ROOT
import tdrstyle
import numpy as np

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/"
# skim_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/NanoAOD-skims/518/mm/"
skim_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing-NEW/NanoAOD-skims/518/mm/"
# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/trigger_info_tmp/";
output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/trigger_info/";

selections = {
    'bsmm_noid':
    'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
    'mm_kin_pt>5 &&' +
    'mm_kin_vtx_prob>0.025 &&' +
    'abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 &&' +
    'Muon_pt[mm_mu1_index]>4.0 && Muon_pt[mm_mu2_index]>4.0 &&' +
    'mm_kin_mass>4.6 && mm_kin_mass<5.9 &&' +
    'Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1',
    
    'bjpsik_noid':
    "mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && "\
    "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && "\
    "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && "\
    "mm_kin_pt[bkmm_mm_index]>7.0 && "\
    "mm_kin_alphaBS[bkmm_mm_index]<0.4 && "\
    "mm_kin_vtx_prob[bkmm_mm_index]>0.1 && "\
    "bkmm_jpsimc_vtx_prob>0.025 && "\
    "mm_kin_sl3d[bkmm_mm_index]>4 && "\
    "abs(bkmm_jpsimc_mass-5.4)<0.5"
}

selections['bsmm']            = selections['bsmm_noid'] + \
    "&& Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45"
selections['bsmm_mc']         = selections['bsmm'] + "&& mm_gen_pdgId!=0"
selections['bsmm_mc_L1_DoubleMu0er1p5_SQ_OS'] = selections['bsmm_mc'] + "&& L1_DoubleMu0er1p5_SQ_OS"
selections['bsmm_mc_HLT_Dimuon0_LowMass_L1_0er1p5'] = selections['bsmm_mc'] + "&& HLT_Dimuon0_LowMass_L1_0er1p5"
selections['bsmm_mc_L1_and_HLT_Objects'] = selections['bsmm_mc'] + "&& MuonId_l1_quality[mm_mu1_index]==12 && MuonId_l1_quality[mm_mu2_index]==12 && MuonId_hlt_pt[mm_mu1_index]>0 && MuonId_hlt_pt[mm_mu2_index]>0"
selections['bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects'] = selections['bsmm_mc_L1_and_HLT_Objects'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"
selections['bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4'] = selections['bsmm_mc'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"

selections['bsmm_mc_ch0']     = selections['bsmm_mc'] + "&& abs(Muon_eta[mm_mu1_index])<0.7 && abs(Muon_eta[mm_mu2_index])<0.7"
selections['bsmm_mc_ch1']     = selections['bsmm_mc'] + "&&(abs(Muon_eta[mm_mu1_index])>0.7 || abs(Muon_eta[mm_mu2_index])>0.7)"
selections['bsmm_mc_chp0']    = selections['bsmm_mc'] + "&& (abs(Muon_eta[mm_mu1_index])<0.7 || abs(Muon_eta[mm_mu2_index])<0.7)"
selections['bsmm_mc_chp1']    = selections['bsmm_mc'] + "&& (abs(Muon_eta[mm_mu1_index])>0.7 && abs(Muon_eta[mm_mu2_index])>0.7)"
selections['bsmm_mc_mva0.9']  = selections['bsmm_mc'] + "&& mm_mva>0.9"
selections['bsmm_mc_mva0.99'] = selections['bsmm_mc'] + "&& mm_mva>0.99"

selections['bjpsik']          = selections['bjpsik_noid'] + "&& Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 " + \
    "&& Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45"
selections['bjpsik_mc']       = selections['bjpsik'] + "&& bkmm_gen_pdgId!=0"
selections['bjpsik_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4'] = selections['bjpsik_mc'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"
selections['bjpsik_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_Objects'] = selections['bjpsik_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4'] + "&& MuonId_hlt_pt[mm_mu1_index[bkmm_mm_index]]>0 && MuonId_hlt_pt[mm_mu2_index[bkmm_mm_index]]>0"
selections['bjpsik_mc_L1_and_HLT_Objects'] = selections['bjpsik_mc'] + "&& MuonId_l1_quality[mm_mu1_index[bkmm_mm_index]]==12 && MuonId_l1_quality[mm_mu2_index[bkmm_mm_index]]==12 && MuonId_hlt_pt[mm_mu1_index[bkmm_mm_index]]>0 && MuonId_hlt_pt[mm_mu2_index[bkmm_mm_index]]>0"


selections['bjpsik_mc_L1_DoubleMu0er1p5_SQ_OS'] = selections['bjpsik_mc'] + "&& L1_DoubleMu0er1p5_SQ_OS"
selections['bjpsik_mc_HLT_Dimuon0_LowMass_L1_0er1p5'] = selections['bjpsik_mc'] + "&& HLT_Dimuon0_LowMass_L1_0er1p5"

selections['bjpsik_mc_ch0']   = selections['bjpsik_mc'] + "&& abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<0.7 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<0.7"
selections['bjpsik_mc_ch1']   = selections['bjpsik_mc'] + "&&(abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])>0.7 || abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])>0.7)"
selections['bjpsik_mc_chp0']  = selections['bjpsik_mc'] + "&& (abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<0.7 || abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<0.7)"
selections['bjpsik_mc_chp1']  = selections['bjpsik_mc'] + "&& (abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])>0.7 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])>0.7)"

variables = {
    'mm_pt':{
        'var': 'mm_kin_pt',
        'title': 'p_{T}(#mu#mu)',
        'min': 0,
        'max': 30,
        'nbins': 60
    },
    'mm_m1_pt':{
        'var': 'Muon_pt[mm_mu1_index]',
        'title': 'p_{T}(#mu_{1})',
        'min': 0,
        'max': 20,
        'nbins': 40
    },
    'mm_m2_pt':{
        'var': 'Muon_pt[mm_mu2_index]',
        'title': 'p_{T}(#mu_{2})',
        'min': 0,
        'max': 20,
        'nbins': 40
    },
    'mm_m1_eta':{
        'var': 'Muon_eta[mm_mu1_index]',
        'title': '#eta(#mu_{1})',
        'min': -1.5,
        'max':  1.5,
        'nbins': 30
    },
    'mm_m2_eta':{
        'var': 'Muon_eta[mm_mu2_index]',
        'title': '#eta(#mu_{2})',
        'min': -1.5,
        'max':  1.5,
        'nbins': 30
    },
    'mm_mva':{
        'var': 'mm_mva',
        'title': 'MVA_{B}',
        'min': 0.9,
        'max': 1,
        'nbins': 50
    },
    'mmk_pt':{
        'var': 'mm_kin_pt[bkmm_mm_index]',
        'title': 'p_{T}(#mu#mu)',
        'min': 0,
        'max': 30,
        'nbins': 60
    },
    'mmk_pt_lowstat':{
        'var': 'mm_kin_pt[bkmm_mm_index]',
        'title': 'p_{T}(#mu#mu)',
        'min': 0,
        'max': 30,
        'nbins': 12
    },
    'mmk_m1_pt':{
        'var': 'Muon_pt[mm_mu1_index[bkmm_mm_index]]',
        'title': 'p_{T}(#mu_{1})',
        'min': 0,
        'max': 20,
        'nbins': 40
    },
    'mmk_m1_pt_lowstat':{
        'var': 'Muon_pt[mm_mu1_index[bkmm_mm_index]]',
        'title': 'p_{T}(#mu_{1})',
        'min': 0,
        'max': 20,
        'nbins': 8
    },
    'mmk_m2_pt':{
        'var': 'Muon_pt[mm_mu2_index[bkmm_mm_index]]',
        'title': 'p_{T}(#mu_{2})',
        'min': 0,
        'max': 20,
        'nbins': 40
    },
    'mmk_m2_pt_lowstat':{
        'var': 'Muon_pt[mm_mu2_index[bkmm_mm_index]]',
        'title': 'p_{T}(#mu_{2})',
        'min': 0,
        'max': 20,
        'nbins': 8
    },
    'mmk_m1_eta':{
        'var': 'Muon_eta[mm_mu1_index[bkmm_mm_index]]',
        'title': '#eta(#mu_{1})',
        'min': -1.5,
        'max':  1.5,
        'nbins': 30
    },
    'mmk_m1_eta_lowstat':{
        'var': 'Muon_eta[mm_mu1_index[bkmm_mm_index]]',
        'title': '#eta(#mu_{1})',
        'min': -1.5,
        'max':  1.5,
        'nbins': 6
    },
    'mmk_m2_eta':{
        'var': 'Muon_eta[mm_mu2_index[bkmm_mm_index]]',
        'title': '#eta(#mu_{2})',
        'min': -1.5,
        'max':  1.5,
        'nbins': 30
    },
    'mmk_m2_eta_lowstat':{
        'var': 'Muon_eta[mm_mu2_index[bkmm_mm_index]]',
        'title': '#eta(#mu_{2})',
        'min': -1.5,
        'max':  1.5,
        'nbins': 6
    },
    'mmk_mva':{
        'var': 'bkmm_bmm_mva',
        'title': 'MVA_{B}',
        'min': 0.9,
        'max': 1,
        'nbins': 50
    },
    'mm_mva_mu': {
        'var': 'min(Muon_softMva[mm_mu1_index], Muon_softMva[mm_mu2_index])',
        'title': 'MVA_{#mu}',
        'min': 0,
        'max': 0.8,
        'nbins': 40
    },
    'mmk_mva_mu': {
        'var': 'min(Muon_softMva[mm_mu1_index[bkmm_mm_index]], Muon_softMva[mm_mu2_index[bkmm_mm_index]])',
        'title': 'MVA_{#mu}',
        'min': 0,
        'max': 0.8,
        'nbins': 40
    },
}


samples = {
    'Bsmm-2018': [
        data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2017':[	
	data_path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2016BF':[
	data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2016GH':[
	data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
    ],
    'Bjpsik-2018': [
        data_path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root",
    ],
    'Bjpsik-2017':[
        data_path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root",
    ],
    'Bjpsik-2016BF':[
    data_path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
    ],
    'Bjpsik-2016GH':[
        data_path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root",
    ],
    
    'Run2018': [
	skim_path + "/EGamma+Run2018A-12Nov2019_UL2018-v2+MINIAOD/*.root",
	skim_path + "/EGamma+Run2018B-12Nov2019_UL2018-v2+MINIAOD/*.root",
	skim_path + "/EGamma+Run2018C-12Nov2019_UL2018-v2+MINIAOD/*.root",
        skim_path + "/EGamma+Run2018D-12Nov2019_UL2018-v4+MINIAOD/*.root"
    ],
}

studies = [
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': 'bjpsik_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': 'bjpsik_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_Objects',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': 'bjpsik_mc_L1_and_HLT_Objects',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },

    {
        'trigger': 'HLT_DoubleMu4_3_Bs',
        'type':'plot',
        'selection': 'bsmm_mc',
        'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
        'sample':'Bsmm-2018'
    },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'plot',
    #     'selection': 'bsmm_mc_L1_and_HLT_Objects',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'type':'plot',
    #     'selection': 'bsmm_mc_L1_and_HLT_Objects',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    
    {
        'trigger': 'HLT_DoubleMu4_3_Bs',
        'type':'plot',
        'selection': 'bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects',
        'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
        'sample':'Bsmm-2018'
    },
    {
        'trigger': 'HLT_DoubleMu4_3_Bs',
        'type':'plot',
        'selection': 'bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
        'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
        'sample':'Bsmm-2018'
    },
    
    
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'efficiency',
    #     'selection': 'bsmm',
    #     'sample':'Bsmm-2018'
    # },

    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'fit',
    #     'selection': 'bsmm',
    #     'sample':'Bsmm-2018'
    # },

    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'efficiency',
    #     'selection': 'bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
    #     'sample':'Bsmm-2018'
    # },

    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'fit',
    #     'selection': 'bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
    #     'sample':'Bsmm-2018'
    # },


    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'efficiency',
    #     'selection': 'bsmm_mc',
    #     'sample':'Bsmm-2018'
    # },

    # {
    #     'fits':'results/Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs.json'
    #     'type':'predict',
    #     'selection': 'bsmm_mc',
    #     'sample':'Bsmm-2018'
    # },
    
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'type':'efficiency',
    #     'selection': 'bsmm',
    #     'sample':'Run2018'
    # },
    

    
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bsmm_mc',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bsmm_mc_ch0',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bsmm_mc_ch1',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bjpsik_mc_ch0',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bjpsik_mc_ch1',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bsmm_mc_L1_DoubleMu0er1p5_SQ_OS',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass_L1_0er1p5',
    #     'selection': 'bjpsik_mc_L1_DoubleMu0er1p5_SQ_OS',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'L1_DoubleMu0er1p5_SQ_OS',
    #     'selection': 'bsmm_mc',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'L1_DoubleMu0er1p5_SQ_OS',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm_mc',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm_mc',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2017'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi_Displaced',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2017'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm_mc',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2016BF'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi_Displaced',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2016BF'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm_mc',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Bsmm-2016GH'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi_Displaced',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_m1_pt','mmk_m2_pt', 'mmk_pt', 'mmk_m1_eta', 'mmk_m2_eta'],
    #     'sample':'Bjpsik-2016GH'
    # },
    
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt', 'mm_mva', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Run2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm_mc_mva0.9',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': 'bsmm_mc_mva0.99',
    #     'variables':['mm_pt', 'mm_m1_pt', 'mm_m2_pt'],
    #     'sample':'Bsmm-2018'
    # },
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': 'bjpsik',
    #     'prescale': 'prescale_HLT_DoubleMu4_3_Jpsi',
    #     'variables':['mmk_pt_lowstat', 'mmk_m1_pt_lowstat', 'mmk_m2_pt_lowstat',
    #                  'mmk_m1_eta_lowstat', 'mmk_m2_eta_lowstat'],
    #     'sample':'Run2018'
    # },


    # 'Bjpsik-2018_DoubleMu4_Jpsi_NoVertexing':
    # {
    #     'trigger': 'HLT_DoubleMu4_Jpsi_NoVertexing',
    #     'selection': selections['bjpsik_mc'],
    #     'variables':['mmk_pt', 'mmk_m1_pt', 'mmk_m2_pt', 'mmk_mva'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass',
    #     'selection': 'bjpsik_mc',
    #     'variables':['mmk_pt', 'mmk_m1_pt', 'mmk_m2_pt', 'mmk_mva'],
    #     'sample':'Bjpsik-2018'
    # },
    # {
    #     'trigger': 'HLT_Dimuon0_LowMass',
    #     'selection': 'bjpsik',
    #     'variables':['mmk_pt', 'mmk_m1_pt', 'mmk_m2_pt', 'mmk_mva'],
    #     'sample':'Run2018'
    # },
    # 'Bsmm-2018-noid':
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Bs',
    #     'selection': selections['bsmm_noid'],
    #     'variables':['mm_mva_mu'],
    #     'sample':'Bsmm-2018'
    # },
    # 'Bjpsik-2018-noid':
    # {
    #     'trigger': 'HLT_DoubleMu4_3_Jpsi',
    #     'selection': selections['bjpsik_noid'],
    #     'variables':['mmk_mva_mu'],
    #     'sample':'Bjpsik-2018'
    # },
]

data = dict()

def get_data(name):
    if name in data:
        return data[name]
    chain = ROOT.TChain("Events")
    for sample in samples[name]:
        chain.Add(sample)
    data[name] = chain
    return data[name]
        

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

hists = dict()

def draw_shades(h_eff, h_ref):
    colors = [ ROOT.kGray, ROOT.kGray + 1, ROOT.kGray + 2, ROOT.kGray + 3 ]

    h_eff.Draw("hist e0")

    nbins = h_ref.GetNbinsX()
    for i_color, color in enumerate(colors):

        h_name = "h_shade_%u" % i_color
        h_shade = h_eff.Clone(h_name)
        h_shade.SetDirectory(0)
        hists[h_name] = h_shade
        h_shade.SetFillColor(color)
        
        bin_low = None
        bin_high = None
        integral = h_ref.Integral(0, nbins + 1)

        shade_fraction = float(i_color + 1) / (len(colors) + 1) / 2
    
        for i in range(1, nbins + 1):
            if h_ref.Integral(1, i) / integral > shade_fraction and bin_low == None:
                bin_low = i
            if h_ref.Integral(1, i) / integral > 1 - shade_fraction and bin_high == None:
                bin_high = i

        if bin_low == None:
            bin_low = 1
        if bin_high == None:
            bin_high = nbins + 1

        h_shade.GetXaxis().SetRange(bin_low, bin_high - 1)
        
        h_shade.Draw("hist e0 same")

def plot(study):
    chain  = get_data(study['sample'])
    name = '%s_%s_%s' % (study['sample'], study['selection'], study['trigger'])
    
    for var in study['variables']:
        var_name = "%s-%s" % (name, var)
        v = variables[var]
        
        h_all = ROOT.TH1F("h_all_%s" % var_name, "", v['nbins'], v['min'], v['max'])
        h_all.GetXaxis().SetTitle(v['title'])
        h_all.Sumw2()
        h_all.SetLineWidth(2)
        chain.Draw("%s>>h_all_%s" % (v['var'], var_name), selections[study['selection']])

        print_canvas("%s" % (var_name), output_path)

        h_trig = h_all.Clone("h_trig_%s" % var_name)
        weight = 1
        if 'prescale' in study:
            weight = study['prescale']
        chain.Draw("%s>>h_trig_%s" % (v['var'], var_name), "(" + selections[study['selection']] + "&&" + study['trigger'] + ")*%s" % weight)

        h_eff = h_trig.Clone("h_eff_%s" % var_name)
        h_eff.Divide(h_trig, h_all, 1, 1, "B")
        h_eff.SetMinimum(0)
        h_eff.SetMaximum(1.1)
        # h_eff.Draw("hist")

        draw_shades(h_eff, h_all)
                    
        print_canvas("%s_eff" % (var_name), output_path)

        h_all.SetDirectory(0)
        hists["h_all_%s" % var_name] = h_all
        hists["h_trig_%s" % var_name] = h_trig
        hists["h_eff_%s" % var_name] = h_eff

def measure_trigger_object_efficiency(study):
    chain  = get_data(study['sample'])
    name = '%s_%s_%s' % (study['sample'], study['selection'], study['trigger'])

    eta_bin_size = 0.5
    eta_bins = np.arange(-1.5, 1.5, eta_bin_size)

    fout = ROOT.TFile.Open("results/%s.root" % name, "recreate")

    for eta in eta_bins:
        for mu in ['mu1', 'mu2']:
            h_name = "h_all_%s_eta%s_%s" % (mu, eta, name)
            h_all = ROOT.TH1F(h_name, "", 32, 4, 20)
            h_all.Sumw2()
            selection = selections[study['selection']] + "&& Muon_eta[mm_%s_index] > %s && Muon_eta[mm_%s_index] < %s" % (mu, eta, mu, eta + eta_bin_size) 
            chain.Draw("%s>>%s" % ("Muon_pt[mm_%s_index]" % mu, h_name), selection)

            print_canvas("%s" % (h_name), output_path)
            h_all.SetDirectory(0)
            hists[h_name] = h_all

            h_name = "h_trig_%s_eta%s_%s" % (mu, eta, name)
            h_trig = h_all.Clone(h_name)
            weight = 1
            if 'prescale' in study:
                weight = study['prescale']
            chain.Draw("%s>>%s" % ("Muon_pt[mm_%s_index]" % mu, h_name), selection + "&& MuonId_hlt_pt[mm_%s_index]>0" % mu)
            # print_canvas("%s" % (h_name), output_path)
            hists[h_name] = h_trig

            h_name = "h_eff_%s_eta%s_%s" % (mu, eta, name)
            h_eff = h_trig.Clone(h_name)
            h_eff.Divide(h_trig, h_all, 1, 1, "B")
            h_eff.SetMinimum(0)
            h_eff.SetMaximum(1.1)
            h_eff.Draw("hist")

            print_canvas("%s" % (h_name), output_path)
            # h_eff.SetDirectory(fout)
            fout.cd()
            h_eff.Write()
        
        
        hists[h_name] = h_eff
        
    fout.Close()

def fit_efficiency(study):
    name = '%s_%s_%s' % (study['sample'], study['selection'], study['trigger'])
    ROOT.gPad.SetGrid()

    eta_bin_size = 0.5
    eta_bins = np.arange(-1.5, 1.5, eta_bin_size)
    fits = dict()
    
    fin = ROOT.TFile.Open("results/%s.root" % name)

    m_erf = ROOT.TF1("m_erf", "[0]*([1]*(TMath::Erf((x - [2])*[3])+1) + (1-[1])*(TMath::Erf((x - [4])*([5]))+1))", 0 , 30);
    m_erf.SetParLimits(1,0,1)
    m_erf.SetParLimits(2,1,10)
    m_erf.SetParLimits(3,0,2)
    m_erf.SetParLimits(4,1,10)
    m_erf.SetParLimits(5,0,2)
    
    for eta in eta_bins:
        for mu in ['mu1', 'mu2']:
            h_name = "h_eff_%s_eta%s_%s" % (mu, eta, name)
            h_eff = fin.Get(h_name)

            h_eff.Fit("m_erf")
            h_eff.Fit("m_erf","M")
            h_eff.SetMaximum(1.3)
            h_eff.Draw()
            print_canvas("%s" % (h_name), output_path)

            fit_name = '%s_%s' % (mu, eta)
            fits[fit_name] = []
            for i in range(6):
                fits[fit_name].append((m_erf.GetParameter(i),m_erf.GetParError(i)))

    json.dump(fits, open("results/%s.json" % name, "w"))
    fin.Close()
        
    
for study in studies:
    if study['type'] == 'plot':
        plot(study)
        continue
    
    if study['type'] == 'efficiency':
        measure_trigger_object_efficiency(study)
        continue

    if study['type'] == 'fit':
        fit_efficiency(study)
        continue
    
        
# Extra plots

def plot_pdfs(id1, name1, id2, name2, xtitle, output_name, left = False):

    try:
        h1 = hists[id1]
        h2 = hists[id2]
        h1.Scale(1/h1.Integral())
        h2.Scale(1/h2.Integral())
        
        m = max(h1.GetMaximum(), h2.GetMaximum()) * 1.2
        h1.SetMaximum(m)
        h2.SetMaximum(m)

        h1.Draw("hist")
        h2.SetMarkerStyle(20)
        h2.SetMarkerColor(ROOT.kBlue)
        h2.Draw("same e")

        if left:
            legend = ROOT.TLegend(0.2,0.80,0.4,0.92)
        else:
            legend = ROOT.TLegend(0.78,0.75,0.85,0.87)
        legend.SetTextSize(0.03)
        legend.SetShadowColor(ROOT.kWhite)
        legend.SetLineColor(ROOT.kWhite)
        legend.AddEntry(h1, name1)
        legend.AddEntry(h2, name2)
        legend.Draw()
    
        print_canvas(output_name, output_path)
    except:
        pass

def overlay_eff_plots(input, output_name, left = False):
    try:
        if left:
            legend = ROOT.TLegend(0.2,0.80,0.4,0.92)
        else:
            legend = ROOT.TLegend(0.78,0.75,0.85,0.87)
        legend.SetTextSize(0.03)
        legend.SetShadowColor(ROOT.kWhite)
        legend.SetLineColor(ROOT.kWhite)

        ROOT.gPad.SetGrid()
        for i, hist in enumerate(input):
            h = hists[hist[0]]
            h.SetMaximum(1.3)
            h.SetLineWidth(2)
            h.SetLineColor(hist[2])
            h.SetMarkerColor(hist[2])
            if i == 0:
                h.Draw("hist")
            else:
                h.Draw("hist same")

            legend.AddEntry(h, hist[1])
            legend.Draw()

        print_canvas(output_name, output_path)
        ROOT.gPad.SetGrid(0,0)
    except:
        pass

# 2018    
plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m1_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{1})", "2018_m1_pt")

plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m2_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{2})", "2018_m2_pt")

plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu#mu})", "2018_mm_pt")

plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m1_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{1})", "2018_m1_eta")

plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m2_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{2})", "2018_m2_eta")

# 2017
plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{1})", "2017_m1_pt")

plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{2})", "2017_m2_pt")

plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu#mu})", "2017_mm_pt")

plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{1})", "2017_m1_eta")

plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{2})", "2017_m2_eta")

# 2016BF
plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{1})", "2016BF_m1_pt")

plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{2})", "2016BF_m2_pt")

plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu#mu})", "2016BF_mm_pt")

plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{1})", "2016BF_m1_eta")

plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{2})", "2016BF_m2_eta")

# 2016GH
plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{1})", "2016GH_m1_pt")

plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu_{2})", "2016GH_m2_pt")

plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_pt", "B #rightarrow J/#psiK",
          "p_{T}(#mu#mu})", "2016GH_mm_pt")

plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{1})", "2016GH_m1_eta")

plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
          "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_eta", "B #rightarrow J/#psiK",
          "#eta(#mu_{2})", "2016GH_m2_eta")

# other
plot_pdfs("h_eff_Bsmm-2018_bsmm_mc_HLT_Dimuon0_LowMass_L1_0er1p5-mm_m1_pt", "Dimuon0_LowMass_L1_0er1p5",
          "h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "HLT_DoubleMu4_3_Bs",
          "#eta(#mu_{2})", "2018_hlt_mm_m1_pt", True)
plot_pdfs("h_eff_Bsmm-2018_bsmm_mc_HLT_Dimuon0_LowMass_L1_0er1p5-mm_m2_pt", "Dimuon0_LowMass_L1_0er1p5",
          "h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "HLT_DoubleMu4_3_Bs",
          "#eta(#mu_{2})", "2018_hlt_mm_m2_pt", True)

overlay_eff_plots(
    [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt",
      "Offline", ROOT.kBlack),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m1_pt",
      "Offline + L1", ROOT.kRed),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m1_pt",
      "Offline + L1 + HLT objects", ROOT.kBlue),
     ],
    "2018_bsmm_trigger_efficiency_mm_m1_pt", True
    )

overlay_eff_plots(
    [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt",
      "Offline", ROOT.kBlack),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m2_pt",
      "Offline + L1", ROOT.kRed),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m2_pt",
      "Offline + L1 + HLT objects", ROOT.kBlue),
     ],
    "2018_bsmm_trigger_efficiency_mm_m2_pt", True
    )

overlay_eff_plots(
    [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta",
      "Offline", ROOT.kBlack),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m1_eta",
      "Offline + L1", ROOT.kRed),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m1_eta",
      "Offline + L1 + HLT objects", ROOT.kBlue),
     ],
    "2018_bsmm_trigger_efficiency_mm_m1_eta", True
    )

overlay_eff_plots(
    [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta",
      "Offline", ROOT.kBlack),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m2_eta",
      "Offline + L1", ROOT.kRed),
     ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m2_eta",
      "Offline + L1 + HLT objects", ROOT.kBlue),
     ],
    "2018_bsmm_trigger_efficiency_mm_m2_eta", True
    )

