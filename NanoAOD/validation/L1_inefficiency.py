"""Trigger Properties

The script performs a study of the L1 inefficiency in data for two
offline event selections:
- two good quality muons from a common vertex
- Bmm like events with MVA_B>0.99

"""
import os, sys, json
import ROOT
import tdrstyle
import numpy as np

debug = False

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/"
mc_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/519/"
# data_path = "/tmp/dmytro/"
output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Run3/L1_trigger/";

selections = {
    'mm':
    'HLT_DoubleMu4_3_LowMass && L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 &&' +
    'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
    'mm_mu1_pt>4 && mm_mu1_pt>3 &&' +
    'Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45 &&' +
    'mm_kin_vtx_prob>0.1 &&' +
    'Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1',

    'bkkmm':
    'mm_mu1_index[bkkmm_mm_index]>=0 && mm_mu2_index[bkkmm_mm_index]>=0 &&' +
    'Muon_softMva[mm_mu1_index[bkkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkkmm_mm_index]]>0.45 &&' +
    'mm_kin_vtx_prob[bkkmm_mm_index]>0.1 &&' +
    'bkkmm_gen_pdgId!=0'
}

selections['mm44']         = selections['mm'] + "&& mm_mu1_pt>4 && mm_mu1_pt>4"
selections['mm_high_mass'] = selections['mm'] + "&& mm_kin_mass>5.5 && nMuon==2 && MuonId_l1_pt[mm_mu1_index]>0 && MuonId_l1_pt[mm_mu2_index]>0"
selections['mm_bmm_mass']  = selections['mm'] + "&& mm_kin_mass>5.2 && mm_kin_mass<5.5"
selections['mm_jpsi_mass'] = selections['mm'] + "&& mm_kin_mass>3.0 && mm_kin_mass<3.2"
selections['bmm']          = selections['mm'] + "&& mm_mva>0.9"
selections['single_mm']    = "nMuon==2 &&" + selections['mm']
selections['bkkmm_er_1.4'] = selections['bkkmm'] + "&& abs(mm_mu1_eta[bkkmm_mm_index])<1.4 && abs(mm_mu2_eta[bkkmm_mm_index])<1.4"
selections['mm_jpsi_mass_er_1.4']    = selections['mm_jpsi_mass'] + "&& abs(mm_mu1_eta)<1.4 && abs(mm_mu2_eta)<1.4"

variables = {
    'mm_mass':{
        'var': 'mm_kin_mass',
        'title': 'm(#mu#mu)',
        'min': 0,
        'max': 8.5,
        'nbins': 85
    },
    'mm_pt':{
        'var': 'mm_kin_pt',
        'title': 'p_{T}(#mu#mu)',
        'min': 0,
        'max': 30,
        'nbins': 60
    },
    'mmkk_pt':{
        'var': 'bkkmm_jpsikk_pt',
        'title': 'p_{T}(B_{S})',
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
        'min': -2.4,
        'max':  2.4,
        'nbins': 28
    },
    'mm_m2_eta':{
        'var': 'Muon_eta[mm_mu2_index]',
        'title': '#eta(#mu_{2})',
        'min': -2.4,
        'max':  2.4,
        'nbins': 28
    },
    'mm_dEta':{
        'var': 'abs(Muon_eta[mm_mu2_index]-Muon_eta[mm_mu1_index])',
        'title': '#Delta#eta(#mu#mu)',
        'min': 0.0,
        'max': 2.0,
        'nbins': 40
    },
    'mm_L1_dEta':{
        'var': 'abs(MuonId_l1_eta[mm_mu2_index]-MuonId_l1_eta[mm_mu1_index])',
        'title': 'L1 #Delta#eta(#mu#mu)',
        'min': 0.0,
        'max': 5.0,
        'nbins': 50
    },
    'mm_dPhi':{
        'var': 'acos(cos(Muon_phi[mm_mu2_index]-Muon_phi[mm_mu1_index]))',
        'title': '#Delta#phi(#mu#mu)',
        'min': 0.0,
        'max': 5.0,
        'nbins': 50
    },
    'mm_dR':{
        'var': 'sqrt(pow(acos(cos(Muon_phi[mm_mu2_index]-Muon_phi[mm_mu1_index])),2)+pow(Muon_eta[mm_mu2_index]-Muon_eta[mm_mu1_index],2))',
        'title': 'dR(#mu#mu)',
        'min': 0.0,
        'max': 5.0,
        'nbins': 50
    },
}


samples = {
    'Run2022C':[
	# data_path + "/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/fb98660d-a04a-4c8d-be5a-4e071f051468.root"
	# data_path + "/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/*.root"
        data_path + "/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/1*.root"
    ],
    'BsToJPsiPhi':[
        # mc_path + '/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/1091A216-A8D3-344E-9F51-0971D7793618.root'
        mc_path + '/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/1*root'
    ]
}

vars = ['mm_pt']
# vars = ['mm_dEta', 'mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta', 'mm_dPhi', 'mm_dR']
 
studies = [
    # {
    #     'test_name': 'L1_dEta_Max1p2',
    #     'test': 'abs(MuonId_l1_eta[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_eta[mm_mu1_index[bkkmm_mm_index]]) < 1.2',
    #     'type':'plot',
    #     'selection': 'bkkmm',
    #     'variables': vars,
    #     'sample':'BsToJPsiPhi'
    # },
    # {
    #     'test_name': 'L1_dR_Max1p6',
    #     'test': 'sqrt(pow(acos(cos(MuonId_l1_phi[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_phi[mm_mu1_index[bkkmm_mm_index]])),2)+pow(MuonId_l1_eta[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_eta[mm_mu1_index[bkkmm_mm_index]],2)) < 1.6',
    #     'type':'plot',
    #     'selection': 'bkkmm',
    #     'variables': vars,
    #     'sample':'BsToJPsiPhi'
    # },
    # {
    #     'test_name': 'L1_dR_Max1p4',
    #     'test': 'sqrt(pow(acos(cos(MuonId_l1_phi[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_phi[mm_mu1_index[bkkmm_mm_index]])),2)+pow(MuonId_l1_eta[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_eta[mm_mu1_index[bkkmm_mm_index]],2)) < 1.4',
    #     'type':'plot',
    #     'selection': 'bkkmm',
    #     'variables': vars,
    #     'sample':'BsToJPsiPhi'
    # },
    # {
    #     'test_name': 'L1_dEta_Max1p2',
    #     'test': 'abs(MuonId_l1_eta[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_eta[mm_mu1_index[bkkmm_mm_index]]) < 1.2',
    #     'type':'plot',
    #     'selection': 'bkkmm_er_1.4',
    #     'variables': vars,
    #     'sample':'BsToJPsiPhi'
    # },
    # {
    #     'test_name': 'L1_dR_Max1p6',
    #     'test': 'sqrt(pow(acos(cos(MuonId_l1_phi[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_phi[mm_mu1_index[bkkmm_mm_index]])),2)+pow(MuonId_l1_eta[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_eta[mm_mu1_index[bkkmm_mm_index]],2)) < 1.6',
    #     'type':'plot',
    #     'selection': 'bkkmm_er_1.4',
    #     'variables': vars,
    #     'sample':'BsToJPsiPhi'
    # },
    # {
    #     'test_name': 'L1_dR_Max1p4',
    #     'test': 'sqrt(pow(acos(cos(MuonId_l1_phi[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_phi[mm_mu1_index[bkkmm_mm_index]])),2)+pow(MuonId_l1_eta[mm_mu2_index[bkkmm_mm_index]]-MuonId_l1_eta[mm_mu1_index[bkkmm_mm_index]],2)) < 1.4',
    #     'type':'plot',
    #     'selection': 'bkkmm_er_1.4',
    #     'variables': vars,
    #     'sample':'BsToJPsiPhi'
    # },
    {
        'test_name': 'L1_dEta_Max1p2',
        'test': 'abs(MuonId_l1_eta[mm_mu2_index]-MuonId_l1_eta[mm_mu1_index]) < 1.2',
        'type':'plot',
        'selection': 'mm_jpsi_mass_er_1.4',
        'variables': vars,
        'sample':'Run2022C'
    },
    {
        'test_name': 'L1_dR_Max1p6',
        'test': 'sqrt(pow(acos(cos(MuonId_l1_phi[mm_mu2_index]-MuonId_l1_phi[mm_mu1_index])),2)+pow(MuonId_l1_eta[mm_mu2_index]-MuonId_l1_eta[mm_mu1_index],2)) < 1.6',
        'type':'plot',
        'selection': 'mm_jpsi_mass_er_1.4',
        'variables': vars,
        'sample':'Run2022C'
    },
    # {
    #     'test_name': 'L1_dEta_Max1p2',
    #     'test': 'abs(MuonId_l1_eta[mm_mu2_index]-MuonId_l1_eta[mm_mu1_index]) < 1.2',
    #     'type':'plot',
    #     'selection': 'mm_jpsi_mass',
    #     'variables': vars,
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_dR_Max1p6',
    #     'test': 'sqrt(pow(acos(cos(MuonId_l1_phi[mm_mu2_index]-MuonId_l1_phi[mm_mu1_index])),2)+pow(MuonId_l1_eta[mm_mu2_index]-MuonId_l1_eta[mm_mu1_index],2)) < 1.6',
    #     'type':'plot',
    #     'selection': 'mm_jpsi_mass',
    #     'variables': vars,
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm',
    #     'variables': vars,
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm_bmm_mass',
    #     'variables': vars,
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
    #     'test': 'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
    #     'type':'plot',
    #     'selection': 'mm',
    #     'variables': vars,
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
    #     'test': 'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
    #     'type':'plot',
    #     'selection': 'mm_bmm_mass',
    #     'variables': vars,
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm',
    #     'variables':['mm_dEta', 'mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta', 'mm_dPhi', 'mm_dR'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'single_mm',
    #     'variables':['mm_dEta', 'mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta', 'mm_dPhi', 'mm_dR'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm_high_mass',
    #     'variables':['mm_L1_dEta', 'mm_dEta', 'mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta', 'mm_dPhi', 'mm_dR'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm',
    #     'variables':['mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm44',
    #     'variables':['mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm_high_mass',
    #     'variables':['mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'mm_bmm_mass',
    #     'variables':['mm_dPhi', 'mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta', 'mm_dEta', 'mm_dR'],
    #     'sample':'Run2022C'
    # },
    # {
    #     'test_name': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'test': 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4_SQ_OS_dR_Max1p2',
    #     'type':'plot',
    #     'selection': 'bmm',
    #     'variables':['mm_mass', 'mm_m1_pt', 'mm_m2_pt', 'mm_m1_eta', 'mm_m2_eta'],
    #     'sample':'Run2022C'
    # },

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
    name = '%s_%s_%s' % (study['sample'], study['selection'], study['test_name'])
    
    for var in study['variables']:
        var_name = "%s-%s" % (name, var)
        v = variables[var]
        
        h_all = ROOT.TH1F("h_all_%s" % var_name, "", v['nbins'], v['min'], v['max'])
        h_all.GetXaxis().SetTitle(v['title'])
        h_all.Sumw2()
        h_all.SetLineWidth(2)
        ndraw = chain.Draw("%s>>h_all_%s" % (v['var'], var_name), selections[study['selection']])
        if debug:
            print("Selection: " + selections[study['selection']])
            print("Selected: %u" % ndraw) 

        print_canvas("%s" % (var_name), output_path)

        h_trig = h_all.Clone("h_trig_%s" % var_name)
        weight = 1
        if 'prescale' in study:
            weight = study['prescale']
        chain.Draw("%s>>h_trig_%s" % (v['var'], var_name), "(" + selections[study['selection']] + "&& (" + study['test'] + "))*%s" % weight)

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
    
        
# # Extra plots

def plot_pdfs(id1, name1, id2, name2, xtitle, output_name, left = False):

    try:
        h1 = hists[id1]
        h2 = hists[id2]
        h1.Scale(1/h1.Integral())
        h2.Scale(1/h2.Integral())
        
        m = max(h1.GetMaximum(), h2.GetMaximum()) * 1.2
        h1.SetMaximum(m)
        h1.SetLineWidth(3)
        h2.SetMaximum(m)
        h2.SetLineWidth(3)

        h1.Draw("hist")
        h2.SetMarkerStyle(20)
        h2.SetMarkerColor(ROOT.kBlue)
        h2.SetLineColor(ROOT.kBlue)
        h2.Draw("same hist")

        if left:
            legend = ROOT.TLegend(0.2,0.80,0.4,0.92)
        else:
            legend = ROOT.TLegend(0.58,0.75,0.85,0.87)
        legend.SetTextSize(0.03)
        legend.SetShadowColor(ROOT.kWhite)
        legend.SetLineColor(ROOT.kWhite)
        legend.AddEntry(h1, name1)
        legend.AddEntry(h2, name2)
        legend.Draw()
    
        print_canvas(output_name, output_path)
    except:
        pass

# def overlay_eff_plots(input, output_name, left = False):
#     try:
#         if left:
#             legend = ROOT.TLegend(0.2,0.80,0.4,0.92)
#         else:
#             legend = ROOT.TLegend(0.78,0.75,0.85,0.87)
#         legend.SetTextSize(0.03)
#         legend.SetShadowColor(ROOT.kWhite)
#         legend.SetLineColor(ROOT.kWhite)

#         ROOT.gPad.SetGrid()
#         for i, hist in enumerate(input):
#             h = hists[hist[0]]
#             h.SetMaximum(1.3)
#             h.SetLineWidth(2)
#             h.SetLineColor(hist[2])
#             h.SetMarkerColor(hist[2])
#             if i == 0:
#                 h.Draw("hist")
#             else:
#                 h.Draw("hist same")

#             legend.AddEntry(h, hist[1])
#             legend.Draw()

#         print_canvas(output_name, output_path)
#         ROOT.gPad.SetGrid(0,0)
#     except:
#         pass

plot_pdfs("h_all_Run2022C_mm_bmm_mass_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2-mm_dEta", 
          "m(#mu#mu) in [5.2, 5.5]",
          "h_all_Run2022C_mm_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2-mm_dEta",
          "m(#mu#mu) in [0.2, 8.5]",
          "#Delta#eta(#mu#mu)", "Run2022C_mm_dEta")

plot_pdfs("h_all_Run2022C_mm_bmm_mass_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2-mm_dR", 
          "m(#mu#mu) in [5.2, 5.5]",
          "h_all_Run2022C_mm_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_OR_L1_DoubleMu4_SQ_OS_dR_Max1p2-mm_dR",
          "m(#mu#mu) in [0.2, 8.5]",
          "dR(#mu#mu)", "Run2022C_mm_dR")



# # 2018    
# plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m1_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{1})", "2018_m1_pt")

# plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m2_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{2})", "2018_m2_pt")

# plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu#mu})", "2018_mm_pt")

# plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m1_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{1})", "2018_m1_eta")

# plot_pdfs("h_all_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2018_bjpsik_mc_HLT_DoubleMu4_3_Jpsi-mmk_m2_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{2})", "2018_m2_eta")

# # 2017
# plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{1})", "2017_m1_pt")

# plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{2})", "2017_m2_pt")

# plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu#mu})", "2017_mm_pt")

# plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{1})", "2017_m1_eta")

# plot_pdfs("h_all_Bsmm-2017_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2017_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{2})", "2017_m2_eta")

# # 2016BF
# plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{1})", "2016BF_m1_pt")

# plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{2})", "2016BF_m2_pt")

# plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu#mu})", "2016BF_mm_pt")

# plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{1})", "2016BF_m1_eta")

# plot_pdfs("h_all_Bsmm-2016BF_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016BF_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{2})", "2016BF_m2_eta")

# # 2016GH
# plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{1})", "2016GH_m1_pt")

# plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu_{2})", "2016GH_m2_pt")

# plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_pt", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_pt", "B #rightarrow J/#psiK",
#           "p_{T}(#mu#mu})", "2016GH_mm_pt")

# plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m1_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{1})", "2016GH_m1_eta")

# plot_pdfs("h_all_Bsmm-2016GH_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta", "B_{s} #rightarrow #mu#mu",
#           "h_all_Bjpsik-2016GH_bjpsik_mc_HLT_DoubleMu4_3_Jpsi_Displaced-mmk_m2_eta", "B #rightarrow J/#psiK",
#           "#eta(#mu_{2})", "2016GH_m2_eta")

# # other
# plot_pdfs("h_eff_Bsmm-2018_bsmm_mc_HLT_Dimuon0_LowMass_L1_0er1p5-mm_m1_pt", "Dimuon0_LowMass_L1_0er1p5",
#           "h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt", "HLT_DoubleMu4_3_Bs",
#           "#eta(#mu_{2})", "2018_hlt_mm_m1_pt", True)
# plot_pdfs("h_eff_Bsmm-2018_bsmm_mc_HLT_Dimuon0_LowMass_L1_0er1p5-mm_m2_pt", "Dimuon0_LowMass_L1_0er1p5",
#           "h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt", "HLT_DoubleMu4_3_Bs",
#           "#eta(#mu_{2})", "2018_hlt_mm_m2_pt", True)

# overlay_eff_plots(
#     [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_pt",
#       "Offline", ROOT.kBlack),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m1_pt",
#       "Offline + L1", ROOT.kRed),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m1_pt",
#       "Offline + L1 + HLT objects", ROOT.kBlue),
#      ],
#     "2018_bsmm_trigger_efficiency_mm_m1_pt", True
#     )

# overlay_eff_plots(
#     [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_pt",
#       "Offline", ROOT.kBlack),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m2_pt",
#       "Offline + L1", ROOT.kRed),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m2_pt",
#       "Offline + L1 + HLT objects", ROOT.kBlue),
#      ],
#     "2018_bsmm_trigger_efficiency_mm_m2_pt", True
#     )

# overlay_eff_plots(
#     [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m1_eta",
#       "Offline", ROOT.kBlack),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m1_eta",
#       "Offline + L1", ROOT.kRed),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m1_eta",
#       "Offline + L1 + HLT objects", ROOT.kBlue),
#      ],
#     "2018_bsmm_trigger_efficiency_mm_m1_eta", True
#     )

# overlay_eff_plots(
#     [("h_eff_Bsmm-2018_bsmm_mc_HLT_DoubleMu4_3_Bs-mm_m2_eta",
#       "Offline", ROOT.kBlack),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_HLT_DoubleMu4_3_Bs-mm_m2_eta",
#       "Offline + L1", ROOT.kRed),
#      ("h_eff_Bsmm-2018_bsmm_mc_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_L1_and_HLT_Objects_HLT_DoubleMu4_3_Bs-mm_m2_eta",
#       "Offline + L1 + HLT objects", ROOT.kBlue),
#      ],
#     "2018_bsmm_trigger_efficiency_mm_m2_eta", True
#     )

