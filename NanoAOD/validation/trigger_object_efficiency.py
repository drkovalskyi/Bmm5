"""Trigger Properties

The script performs a study of trigger properties using MC
simulations. The main goal is to find the primary factors affecting
the trigger efficiency.

"""
import os, sys, json, re
import ROOT
import tdrstyle
from array import array

tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 1200, 1200)
ROOT.gPad.SetGrid()
# ROOT.gPad.SetLogx()

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing-NEW/FlatNtuples/518/trig-eff"
data_path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing-NEW/FlatNtuples/518/trig-info"
# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/trigger_object_efficiency_test10/";
output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/Bmm/AN/trigger_object_efficiency/";
# output_path = "/eos/user/d/dmytro/www/plots/Bmm/AN/trigger_object_efficiency_test/"

samples = {
    'Bsmm-2018': [
        data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2018-MC-charm': [
        data_path2 + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2017':[	
        data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2017-MC-charm':[	
        data_path2 + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2016BF':[
        data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2016GH':[
        data_path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
    ],
    'Bsmm-2016BF-MC-charm':[
        data_path2 + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
    ],
    'Bsmm-2016GH-MC-charm':[
        data_path2 + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
    ],
    'Run2018': [
        data_path + "/SingleMuon+Run2018A-12Nov2019_UL2018-v5+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2018B-12Nov2019_UL2018-v3+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2018C-12Nov2019_UL2018-v3+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2018D-12Nov2019_UL2018-v8+MINIAOD/*.root",
    ],
    'Run2018-DM': [
        data_path + "/DoubleMuon+Run2018A-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2018B-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2018C-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2018D-12Nov2019_UL2018-v3+MINIAOD/*.root",
    ],
    'Run2018-SM-charm': [
        data_path2 + "/SingleMuon+Run2018A-12Nov2019_UL2018-v5+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2018B-12Nov2019_UL2018-v3+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2018C-12Nov2019_UL2018-v3+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2018D-12Nov2019_UL2018-v8+MINIAOD/*.root",
    ],
    'Run2018-DM-charm': [
        data_path2 + "/DoubleMuon+Run2018A-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2018B-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2018C-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2018D-12Nov2019_UL2018-v3+MINIAOD/*.root",
    ],
    'Run2018-EG-charm': [
        data_path2 + "/EGamma+Run2018A-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path2 + "/EGamma+Run2018B-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path2 + "/EGamma+Run2018C-12Nov2019_UL2018-v2+MINIAOD/*.root",
        data_path2 + "/EGamma+Run2018D-12Nov2019_UL2018-v4+MINIAOD/*.root",
    ],
    'Run2017': [
        # data_path + "/SingleMuon+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
        # data_path + "/DoubleMuon+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        # data_path + "/DoubleMuon+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        # data_path + "/DoubleMuon+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        # data_path + "/DoubleMuon+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        # data_path + "/DoubleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root", 
    ],
    'Run2017-DM': [
        # data_path + "/DoubleMuon+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path + "/DoubleMuon+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path + "/DoubleMuon+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path + "/DoubleMuon+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path + "/DoubleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root", 
    ],
    'Run2017-SM-charm': [
        # data_path + "/SingleMuon+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
    ],
    'Run2017-DM-charm': [
        # data_path + "/DoubleMuon+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path2 + "/DoubleMuon+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path2 + "/DoubleMuon+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path2 + "/DoubleMuon+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root", 
        data_path2 + "/DoubleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root", 
    ],
    'Run2017-EG': [
        # data_path + "/DoubleEG+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
        # data_path + "/SingleElectron+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2017F-09Aug2019_UL2017_rsb-v2+MINIAOD/*.root",
        # data_path + "/SinglePhoton+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
    ],
    'Run2017-EG-charm': [
        # data_path2 + "/DoubleEG+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
        # data_path2 + "/SingleElectron+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2017F-09Aug2019_UL2017_rsb-v2+MINIAOD/*.root",
        # data_path2 + "/SinglePhoton+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
    ],
    'Run2016BF': [
        data_path + "/SingleMuon+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    ],
    'Run2016GH': [
        data_path + "/SingleMuon+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/SingleMuon+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    ],
    'Run2016BF-SM-charm': [
        data_path2 + "/SingleMuon+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    ],
    'Run2016GH-SM-charm': [
        data_path2 + "/SingleMuon+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/SingleMuon+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    ],
    'Run2016BF-DM': [
        data_path + "/DoubleMuon+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    ],
    'Run2016GH-DM': [
        data_path + "/DoubleMuon+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/DoubleMuon+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    ],
    'Run2016BF-DM-charm': [
        data_path2 + "/DoubleMuon+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    ],
    'Run2016GH-DM-charm': [
        data_path2 + "/DoubleMuon+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/DoubleMuon+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    ],
    'Run2016BF-EG': [
        # data_path + "/DoubleEG+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016F-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        # data_path + "/SinglePhoton+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    ],
    'Run2016GH-EG': [
        data_path + "/DoubleEG+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/DoubleEG+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/SingleElectron+Run2016H-21Feb2020_UL2016-v2+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path + "/SinglePhoton+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    ],
    'Run2016BF-EG-charm': [
        # data_path2 + "/DoubleEG+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016F-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        # data_path2 + "/SinglePhoton+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
    ],
    'Run2016GH-EG-charm': [
        data_path2 + "/DoubleEG+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/DoubleEG+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/SingleElectron+Run2016H-21Feb2020_UL2016-v2+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
        data_path2 + "/SinglePhoton+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
    ],
    # 'Run2018A': [ data_path + "/SingleMuon+Run2018A-12Nov2019_UL2018-v5+MINIAOD/*.root" ],
    # 'Run2018B': [ data_path + "/SingleMuon+Run2018B-12Nov2019_UL2018-v3+MINIAOD/*.root" ],
    # 'Run2018C': [ data_path + "/SingleMuon+Run2018C-12Nov2019_UL2018-v3+MINIAOD/*.root" ],
    # 'Run2018D': [ data_path + "/SingleMuon+Run2018D-12Nov2019_UL2018-v8+MINIAOD/*.root" ],
    
    # 'Run2017B': [ data_path + "/SingleMuon+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root" ],
    # 'Run2017C': [ data_path + "/SingleMuon+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root" ],
    # 'Run2017D': [ data_path + "/SingleMuon+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root" ],
    # 'Run2017E': [ data_path + "/SingleMuon+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root" ],
    # 'Run2017F': [ data_path + "/SingleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root" ],
    # 'Run2016B': [ data_path + "/SingleMuon+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root" ],
    # 'Run2016C': [ data_path + "/SingleMuon+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root" ],
    # 'Run2016D': [ data_path + "/SingleMuon+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root" ],
    # 'Run2016E': [ data_path + "/SingleMuon+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root" ],
    # 'Run2016F': [ data_path + "/SingleMuon+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root" ],
    # 'Run2016G': [ data_path + "/SingleMuon+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root" ],
    # 'Run2016H': [ data_path + "/SingleMuon+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root" ]
}

data = dict()

def get_data(name, tree="muons"):
    if name in data:
        return data[name]
    chain = ROOT.TChain(tree)
    n_files = 0 
    for sample in samples[name]:
        n_files += chain.Add(sample)
    if n_files == 0:
        print "No files found for " + name
        return None
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

def measure_trigger_object_efficiency(sample, suffix, preselection):
    chain  = get_data(sample)
    if not chain: return
    name = sample + suffix

    # eta_bin_step = 1
    eta_bin_step = 5
    eta_bins = range(-15, 15, eta_bin_step)
    eta_bin_size = 0.1 * eta_bin_step

    fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

    for ieta in eta_bins:
        eta = ieta * 0.1
        h_name = "h_all_eta%s_%s" % (ieta, name)
        h_all = ROOT.TH1F(h_name, "", 32, 4, 20)
        h_all.Sumw2()
        selection = "probe_eta > %s && probe_eta < %s" % (eta, eta + eta_bin_size)
        if preselection != "":
            selection += "&&" + preselection
        chain.Draw("%s>>%s" % ("probe_pt", h_name), selection)

        print_canvas("%s" % (h_name), output_path)
        h_all.SetDirectory(0)
        hists[h_name] = h_all

        h_name = "h_trig_eta%s_%s" % (ieta, name)
        h_trig = h_all.Clone(h_name)
        chain.Draw("%s>>%s" % ("probe_pt", h_name), selection + "&& probe_hlt_pt>0" )
        print_canvas("%s" % (h_name), output_path)
        hists[h_name] = h_trig

        h_name = "h_eff_eta%s_%s" % (ieta, name)
        h_eff = h_trig.Clone(h_name)
        h_eff.Divide(h_trig, h_all, 1, 1, "B")
        h_eff.SetMinimum(0)
        h_eff.SetMaximum(1.1)
        h_eff.Draw("hist")

        print_canvas("%s" % (h_name), output_path)
        # h_eff.SetDirectory(fout)
        # fout.cd()
        # h_eff.Write()
        
        hists[h_name] = h_eff
        
    fout.Close()

# def _measure_trigger_object_efficiency_mm(chain, bin_name, name, eta_min, eta_max, preselection):
#     h_name = "h_all_%s_%s"    % (bin_name, name)
#     h_name2 = "h_all_%s_%s_2" % (bin_name, name)
#     h_all = ROOT.TH1F(h_name, "", 32, 4, 20)
#     h_all.Sumw2()
#     h_all2 = ROOT.TH1F(h_name2, "", 32, 4, 20)
#     h_all2.Sumw2()
#     selection = preselection
#     if selection != "":
#         selection += "&&" 
#     chain.Draw("%s>>%s" % ("m1pt", h_name), selection + "m1eta > %s && m1eta < %s" % (eta_min, eta_max))
#     chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + "m2eta > %s && m2eta < %s" % (eta_min, eta_max))
#     h_all.Add(h_all2)
#     h_all.Draw("")

#     print_canvas("%s" % (h_name), output_path)
#     h_all.SetDirectory(0)
#     hists[h_name] = h_all

#     h_name  = "h_trig_%s_%s" % (bin_name, name)
#     h_name2 = "h_trig_%s_%s_2" % (bin_name, name)
#     h_trig  = h_all.Clone(h_name)
#     h_trig2 = h_all.Clone(h_name2)
#     chain.Draw("%s>>%s" % ("m1pt", h_name), selection + "m1_hlt_pt>0 && m1eta > %s && m1eta < %s" % (eta_min, eta_max))
#     chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + "m2_hlt_pt>0 && m2eta > %s && m2eta < %s" % (eta_min, eta_max))
#     h_trig.Add(h_trig2)

#     hists[h_name] = h_trig

#     h_name = "h_eff_%s_%s" % (bin_name, name)
#     h_eff = h_trig.Clone(h_name)
#     h_eff.Divide(h_trig, h_all, 1, 1, "B")
#     h_eff.SetMinimum(0)
#     h_eff.SetMaximum(1.1)
#     h_eff.Draw("hist")

#     print_canvas("%s" % (h_name), output_path)
#     # h_eff.SetDirectory(fout)
#     # fout.cd()
#     # h_eff.Write()

#     hists[h_name] = h_eff

def _measure_trigger_object_efficiency_mm(chain, eta_bin_name, name, eta_cut, preselection, tag_selection, trigger, fout):
    h_name = "h_all_%s_%s"    % (eta_bin_name, name)
    h_name2 = "h_all_%s_%s_2" % (eta_bin_name, name)

    pt_bins = [4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.25, 5.5, 5.75, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0, 20.0]

    h_all = ROOT.TH1F(h_name, ";p_{T}", 32, 4, 20)
    # h_all = ROOT.TH1F(h_name, ";p_{T}", len(pt_bins)-1, array('d', pt_bins))
    h_all.Sumw2()
    h_all2 = h_all.Clone(h_name2)
    selection = preselection
    if selection != "":
        selection += "&&"
    tag1_selection = tag_selection.format(index=1)
    tag2_selection = tag_selection.format(index=2)

    chain.Draw("%s>>%s" % ("m1pt", h_name),  selection + tag2_selection + "&&" + eta_cut.format(index=1))
    chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + tag1_selection + "&&" + eta_cut.format(index=2))
    h_all.Add(h_all2)
    h_all.Draw("")

    print_canvas("%s" % (h_name), output_path)
    h_all.SetDirectory(0)
    hists[h_name] = h_all

    h_name  = "h_trig_%s_%s"   % (eta_bin_name, name)
    h_name2 = "h_trig_%s_%s_2" % (eta_bin_name, name)
    h_trig  = h_all.Clone(h_name)
    h_trig2 = h_all.Clone(h_name2)
    chain.Draw("%s>>%s" % ("m1pt", h_name),  selection + tag2_selection + "&& %s && %s" % (eta_cut.format(index=1), trigger))
    chain.Draw("%s>>%s" % ("m2pt", h_name2), selection + tag1_selection + "&& %s && %s" % (eta_cut.format(index=2), trigger))
    h_trig.Add(h_trig2)

    hists[h_name] = h_trig

    h_name = "h_eff_%s_%s" % (eta_bin_name, name)
    h_eff = h_trig.Clone(h_name)
    h_eff.Divide(h_trig, h_all, 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.SetMaximum(1.1)
    h_eff.Draw("hist e")

    print_canvas("%s" % (h_name), output_path)
    h_eff.SetDirectory(fout)
    fout.cd()
    h_eff.Write()
    h_eff.SetDirectory(0)
    
    hists[h_name] = h_eff
    
    
    
def measure_trigger_object_efficiency_mm(sample, suffix, preselection, tag_selection, trigger):
    chain  = get_data(sample, "mm")
    if not chain: return
    name = sample + suffix

    # eta_bin_step = 1
    eta_bin_step = 5
    eta_bins = range(-15, 15, eta_bin_step)
    eta_bin_size = 0.1 * eta_bin_step

    fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

    # _measure_trigger_object_efficiency_mm(chain, "", name, -1.5, 1.5, preselection)
    _measure_trigger_object_efficiency_mm(chain, "",      name, "abs(m{index}eta) < 1.5", preselection, tag_selection, trigger, fout)
    _measure_trigger_object_efficiency_mm(chain, "chan0", name, "abs(m{index}eta) < 0.7", preselection, tag_selection, trigger, fout)
    _measure_trigger_object_efficiency_mm(chain, "chan1", name, "abs(m{index}eta) > 0.7", preselection, tag_selection, trigger, fout)
    
    # for ieta in eta_bins:
    #     eta = ieta * 0.1
    #     _measure_trigger_object_efficiency_mm(chain, "eta%s" % ieta, name, "m{index}eta > %s && m{index}eta < %s" % (eta, eta + eta_bin_size), preselection, tag_selection, trigger)
        
    fout.Close()
    

def measure_trigger_object_efficiency_2D(sample, suffix, binning, draw_options, selection):
    chain  = get_data(sample)
    if not chain: return
    name = sample + suffix

    fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

    h_name = "h2_all_%s" % name
    # h_all = ROOT.TH2F(h_name, "", 32, 4, 20, 30, -1.5, 1.5)
    h_all = ROOT.TH2F(h_name, "", *binning)
    h_all.Sumw2()
    chain.Draw("probe_eta:probe_pt>>%s" % h_name, selection)

    h_all.SetDirectory(0)
    hists[h_name] = h_all

    h_name = "h2_trig_%s" % name
    h_trig = h_all.Clone(h_name)
    if selection != "":
        selection += "&&"
    selection += "probe_hlt_pt>0"
    chain.Draw("probe_eta:probe_pt>>%s" % h_name, selection )
    # print_canvas("%s" % (h_name), output_path)
    hists[h_name] = h_trig

    h_name = "h2_eff_%s" % name
    h_eff = h_trig.Clone(h_name)
    h_eff.Divide(h_trig, h_all, 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.SetMaximum(1.0)
    h_eff.SetMarkerSize(0.5)
    h_eff.Draw(draw_options)

    print_canvas("%s" % (h_name), output_path)
    # h_eff.SetDirectory(fout)
    fout.cd()
    h_eff.Write()
        
    hists[h_name] = h_eff
        
    fout.Close()

def measure_trigger_object_efficiency_2D_generic(sample, suffix, var_x, var_y, binning, draw_options, selection):
    chain  = get_data(sample)
    if not chain: return
    name = sample + suffix

    h_name = "h2_all_%s" % name
    h_all = ROOT.TH2F(h_name, "", *binning)
    h_all.Sumw2()
    h_all.SetMarkerSize(0.5)
    chain.Draw("%s:%s>>%s" % (var_y, var_x, h_name), selection, draw_options)
    print_canvas("%s" % (h_name), output_path)

    h_all.SetDirectory(0)
    hists[h_name] = h_all

    h_name = "h2_trig_%s" % name
    h_trig = h_all.Clone(h_name)
    if selection != "":
        selection += "&&"
    selection += "probe_hlt_pt>0"
    chain.Draw("%s:%s>>%s" % (var_y, var_x, h_name), selection )
    # print_canvas("%s" % (h_name), output_path)
    hists[h_name] = h_trig

    h_name = "h2_eff_%s" % name
    h_eff = h_trig.Clone(h_name)
    h_eff.Divide(h_trig, h_all, 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.SetMaximum(1)
    h_eff.SetMarkerSize(0.5)
    h_eff.Draw(draw_options)

    print_canvas("%s" % (h_name), output_path)
    # h_eff.SetDirectory(fout)
        
    hists[h_name] = h_eff

# def measure_trigger_object_efficiency_2D_mm(sample, suffix, binning, draw_options, selection):
#     chain  = get_data(sample, "mm")
#     if not chain: return
#     name = sample + suffix
    
#     fout = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name, "recreate")

#     h_name  = "h2_all_%s" % name
#     h_name2 = "h2_all_%s_2" % name
#     # h_all = ROOT.TH2F(h_name, "", 32, 4, 20, 30, -1.5, 1.5)
#     h_all = ROOT.TH2F(h_name, "", *binning)
#     h_all.Sumw2()
#     h_all2 = ROOT.TH2F(h_name2, "", *binning)
#     h_all2.Sumw2()
#     chain.Draw("m1eta:m1pt>>%s" % h_name, selection)
#     chain.Draw("m2eta:m2pt>>%s" % h_name2, selection)
#     h_all.Add(h_all2)
#     h_all.Draw("colz")
#     print_canvas("%s" % (h_name), output_path)
    
#     h_all.SetDirectory(0)
#     hists[h_name] = h_all

#     h_name  = "h2_trig_%s" % name
#     h_name2 = "h2_trig_%s_2" % name
#     h_trig  = h_all.Clone(h_name)
#     h_trig2 = h_all.Clone(h_name2)
#     if selection != "":
#         selection += "&&"
#     chain.Draw("m1eta:m1pt>>%s" % h_name,  selection + "m1_hlt_pt>0")
#     chain.Draw("m2eta:m2pt>>%s" % h_name2, selection + "m2_hlt_pt>0")
#     h_trig.Add(h_trig2)
#     hists[h_name] = h_trig

#     h_name = "h2_eff_%s" % name
#     h_eff = h_trig.Clone(h_name)
#     h_eff.Divide(h_trig, h_all, 1, 1, "B")
#     h_eff.SetMinimum(0)
#     h_eff.SetMaximum(1.0)
#     h_eff.SetMarkerSize(0.5)
#     h_eff.Draw(draw_options)

#     print_canvas("%s" % (h_name), output_path)
#     # h_eff.SetDirectory(fout)
#     fout.cd()
#     h_eff.Write()
        
#     hists[h_name] = h_eff
        
#     fout.Close()

def measure_trigger_object_efficiency_2D_mm(sample, suffix, var_x, var_y, binning,
                                            draw_options, preselection, tag_selection, trigger):
    ## Data source
    chain  = get_data(sample, "mm")
    if not chain: return
    
    name = sample + suffix

    ## Selection
    selection = preselection
    if selection != "":
        selection += "&&"
    tag1_selection = tag_selection.format(index=1)
    tag2_selection = tag_selection.format(index=2)

    ## Histograms
    h_name  = "h2_all_%s"   % name
    h_name2 = "h2_all_%s_2" % name
    h_all = ROOT.TH2F(h_name, "", *binning)
    h_all.Sumw2()
    h_all.SetMarkerSize(0.5)
    h_all2 = ROOT.TH2F(h_name2, "", *binning)
    h_all2.Sumw2()

    ## Denominator
    chain.Draw("%s:%s>>%s" % (var_y.format(index=1), var_x.format(index=1), h_name),
               selection + tag2_selection, draw_options)
    chain.Draw("%s:%s>>%s" % (var_y.format(index=2), var_x.format(index=2), h_name2),
               selection + tag1_selection, draw_options)
    h_all.Add(h_all2)
    h_all.Draw(draw_options)
    print_canvas("%s" % (h_name), output_path)

    h_all.SetDirectory(0)
    hists[h_name] = h_all

    ## Numerator
    h_name  = "h2_trig_%s"   % name
    h_name2 = "h2_trig_%s_2" % name
    h_trig  = h_all.Clone(h_name)
    h_trig2 = h_all.Clone(h_name2)
    chain.Draw("%s:%s>>%s" % (var_y.format(index=1), var_x.format(index=1), h_name),
               selection + tag2_selection + "&&" + trigger, draw_options)
    chain.Draw("%s:%s>>%s" % (var_y.format(index=2), var_x.format(index=2), h_name2),
               selection + tag1_selection + "&&" + trigger, draw_options)
    h_trig.Add(h_trig2)
    # print_canvas("%s" % (h_name), output_path)
    hists[h_name] = h_trig

    ## Efficiency
    h_name = "h2_eff_%s" % name
    h_eff = h_trig.Clone(h_name)
    h_eff.Divide(h_trig, h_all, 1, 1, "B")
    h_eff.SetMinimum(0)
    h_eff.SetMaximum(1)
    h_eff.SetMarkerSize(0.5)
    h_eff.Draw(draw_options)

    print_canvas("%s" % (h_name), output_path)
    # h_eff.SetDirectory(fout)
        
    hists[h_name] = h_eff
    
    
# def fit_efficiency(sample):
#     name = sample
#     ROOT.gPad.SetGrid()

#     eta_bin_size = 0.5
#     eta_bins = np.arange(-1.5, 1.5, eta_bin_size)
#     fits = dict()
    
#     fin = ROOT.TFile.Open("results/trigger_object_efficiency_%s.root" % name)

#     m_erf = ROOT.TF1("m_erf", "[0]*([1]*(TMath::Erf((x - [2])*[3])+1) + (1-[1])*(TMath::Erf((x - [4])*([5]))+1))", 0 , 30);
#     m_erf.SetParLimits(1,0,1)
#     m_erf.SetParLimits(2,1,10)
#     m_erf.SetParLimits(3,0,2)
#     m_erf.SetParLimits(4,1,10)
#     m_erf.SetParLimits(5,0,2)
    
#     for eta in eta_bins:
#         h_name = "h_eff_eta%s_%s" % (eta, name)
#         h_eff = fin.Get(h_name)

#         h_eff.Fit("m_erf")
#         h_eff.Fit("m_erf","M")
#         h_eff.SetMaximum(1.3)
#         h_eff.Draw()
#         print_canvas("%s" % (h_name), output_path)

#         fit_name = '%s' % (eta)
#         fits[fit_name] = []
#         for i in range(6):
#             fits[fit_name].append((m_erf.GetParameter(i),m_erf.GetParError(i)))

#     json.dump(fits, open("results/trigger_object_efficiency_%s.json" % name, "w"))
#     fin.Close()
        
ROOT.gStyle.SetPaintTextFormat(".3f");
    
for sample in samples:

    # if not re.search("2018", sample): continue
        
    l1_trigger = "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"
    if re.search("2016", sample):
        l1_trigger = "L1_DoubleMu0er1p6_dEta_Max1p8_OS"

    # selection = l1_trigger + " && tag_q*probe_q==-1 && tag_probe_l1_dr<1.4 && tag_probe_l1_dr>0.1 &&" +\
    #     "abs(tag_eta)<1.4 && abs(probe_eta)<1.4 && tag_l1_quality==12 && probe_l1_quality==12"
    
    ## measure_trigger_object_efficiency_2D(sample, "_zoom",   (8, 4, 8, 30, -1.5, 1.5), "colz texte")
    
    ## Dimuon
    # if re.search('Run201.[^\-]*$', sample):
    #     selection = l1_trigger + " && tag_q*probe_q==-1 && tag_probe_dr<2.0 && tag_probe_dr>0.01 &&" +\
    #         "abs(tag_eta)<1.4 && abs(probe_eta)<1.4 && tag_l1_quality==12 && probe_l1_quality==12 && tag_pt>17 && nl1==2"
    #     measure_trigger_object_efficiency(sample, selection)
    #     measure_trigger_object_efficiency_2D(sample,      "", (32, 4, 20, 30, -1.5, 1.5), "colz", selection)

    # # measure_trigger_object_efficiency_2D_generic(sample, "_dr_pt", "probe_pt", "tag_probe_dr",
    # #                                              (16, 4, 20, 20, 0, 2.0), "colz text", selection)
    # # measure_trigger_object_efficiency_2D_generic(sample, "_l1dr_pt", "probe_pt", "tag_probe_l1_dr",
    # #                                              (16, 4, 20, 20, 0, 2.0), "colz text", selection)

    # ## Single muon
    # if re.search('\-EG$', sample):
    #     selection = "abs(probe_eta)<1.4 && probe_l1_quality==12 && nl1==1"
    #     measure_trigger_object_efficiency(sample, "_single_EG", selection)
    #     measure_trigger_object_efficiency_2D(sample, "_single_EG", (32, 4, 20, 30, -1.5, 1.5), "colz", selection)

    # ====================================================================================
    
    ## Dimuon Charmonium
    if re.search('EG\-charm$', sample):
        selection = l1_trigger + "&& m1_l1_pt>0 && m2_l1_pt>0"
        # measure_trigger_object_efficiency_mm(sample, "_charm", selection)
        # measure_trigger_object_efficiency_2D_mm(sample, "_charm", (32, 4, 20, 30, -1.5, 1.5), "colz", selection)
        measure_trigger_object_efficiency_mm(sample, "_Bs", selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                            "m{index}_hlt_pt>0", "HLT_DoubleMu4_3_Bs")
    
    if re.search('SM\-charm$', sample):
        selection = l1_trigger + "&& m1_l1_pt>0 && m2_l1_pt>0"
        measure_trigger_object_efficiency_mm(sample, "_Bs", selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                             "(m{index}_HLT_IsoMu24||m{index}_HLT_IsoMu27)", "HLT_DoubleMu4_3_Bs")
        measure_trigger_object_efficiency_2D_mm(sample, "_eta_pt_Bs", "m{index}pt", "m{index}eta",
                                                (32, 4, 20, 30, -1.5, 1.5), "colz",
                                                selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                                "(m{index}_HLT_IsoMu24||m{index}_HLT_IsoMu27)", "HLT_DoubleMu4_3_Bs")
        measure_trigger_object_efficiency_2D_mm(sample, "_dr_pt_Bs", "m{index}pt", "sqrt(pow(acos(cos(m2phi-m1phi)),2)+pow(m2eta-m1eta,2))",
                                                (16, 4, 20, 20, 0, 2.0), "colz",
                                                selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                                "(m{index}_HLT_IsoMu24||m{index}_HLT_IsoMu27)", "HLT_DoubleMu4_3_Bs")
        
    if re.search('DM\-charm$', sample):
        selection = l1_trigger + "&& m1_l1_pt>0 && m2_l1_pt>0"
        measure_trigger_object_efficiency_mm(sample, "_Bs", selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                             "(m{index}_HLT_Mu8||m{index}_HLT_Mu17)", "HLT_DoubleMu4_3_Bs")
    if re.search('MC\-charm$', sample):
        selection = l1_trigger + "&& m1_l1_pt>0 && m2_l1_pt>0"
        measure_trigger_object_efficiency_mm(sample, "_SM_Bs", selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                             "(m{index}_HLT_Mu8||m{index}_HLT_Mu17)", "HLT_DoubleMu4_3_Bs")
        measure_trigger_object_efficiency_mm(sample, "_DM_Bs", selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                             "(m{index}_HLT_IsoMu24||m{index}_HLT_IsoMu27)", "HLT_DoubleMu4_3_Bs")
        measure_trigger_object_efficiency_mm(sample, "_EG_Bs", selection + "&& good_for_HLT_DoubleMu4_3_Bs",
                                             "m{index}_hlt_pt>0", "HLT_DoubleMu4_3_Bs")
    
    # ====================================================================================
        
    # measure_trigger_object_efficiency_2D_generic(sample, "_dr_pt_sg", "probe_pt", "tag_probe_dr",
    #                                              (16, 4, 20, 20, 0, 2.0), "colz text", selection + "&& (tag_q-probe_q)*(tag_phi-probe_phi)<0")
    # measure_trigger_object_efficiency_2D_generic(sample, "_dr_pt_cw", "probe_pt", "tag_probe_dr",
    #                                             (16, 4, 20, 20, 0, 2.0), "colz text", selection + "&& (tag_q-probe_q)*(tag_phi-probe_phi)>0")

    
    ## measure_trigger_object_efficiency(sample)

    ## fit_efficiency(sample)

def make_overlay_plot(name, h1_name, name1, h2_name, name2):
    legend = ROOT.TLegend(0.70,0.65,0.85,0.77)
    legend.SetShadowColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    # legend.SetFillColor(10)
    # colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kOrange+5, ROOT.kGreen+3]
    # scale = 1.2

    if h1_name not in hists:
        print "%s is no available. Skip" % h1_name
        return 
    if h2_name not in hists:
        print "%s is no available. Skip" % h2_name
        return 
    h1 = hists[h1_name]
    h1.SetLineColor(ROOT.kRed)
    h1.SetLineWidth(2)
    h1.SetMarkerStyle(20)
    h1.SetMarkerColor(ROOT.kRed)
    # h.GetXaxis().SetTitle("nPV")
    h1.Draw("hist")
    legend.AddEntry(h1, name1)

    h2 = hists[h2_name]
    h2.SetLineColor(ROOT.kBlack)
    h2.SetLineWidth(2)
    h2.SetMarkerStyle(20)
    h2.SetMarkerColor(ROOT.kBlack)
    # h.GetXaxis().SetTitle("nPV")
    h2.Draw("same e")
    legend.AddEntry(h2, name2)
    
    legend.Draw()
    print_canvas(name, output_path)
    
    ratio_plot = ROOT.TRatioPlot(h2, h1)
    # ratio_plot.SetH1DrawOpt("hist e")
    ratio_plot.SetH1DrawOpt("e")
    ratio_plot.SetH2DrawOpt("hist")
    ratio_plot.Draw()
    ratio_plot.SetSeparationMargin(0.03)
    ratio_plot.GetLowerRefGraph().SetMinimum(0.4)
    ratio_plot.GetLowerRefGraph().SetMaximum(1.6)
    # ratio_plot.GetXaxis().SetTitleSize()
    # SetBottomMargin(2.0)
    # c1.SetBottomMargin(2.0)
    # rp->GetLowerRefYaxis()->SetRange(...)
    # rp->SetH1DrawOpt("E");
    ratio_plot.GetLowerRefYaxis().SetTitle("Data/MC")
        
    ratio_plot.GetLowerRefYaxis().SetTitleSize()
    ratio_plot.GetLowerRefYaxis().SetTitleOffset(1.1)
    ratio_plot.GetLowerRefYaxis().SetLabelSize(0.035)
    ratio_plot.GetLowYaxis().SetNdivisions(503)
        
    ratio_plot.GetLowerRefXaxis().SetTitleSize()
    ratio_plot.GetLowerRefXaxis().SetTitleOffset()
    ratio_plot.GetLowerRefXaxis().SetLabelSize(0.035)
    
    ratio_plot.GetUpperRefYaxis().SetTitle("")
    ratio_plot.GetUpperRefYaxis().SetTitleSize()
    ratio_plot.GetUpperRefYaxis().SetTitleOffset()
    ratio_plot.GetUpperRefYaxis().SetLabelSize(0.035)
        
    ratio_plot.GetUpperRefXaxis().SetTitleSize()
    ratio_plot.GetUpperRefXaxis().SetTitleOffset(0)
    ratio_plot.GetUpperRefXaxis().SetLabelSize(0.035)

    c1.Update()
    legend.Draw()
    print_canvas(name + "_ratio", output_path)

make_overlay_plot("h_eff_chan0_2018-SM-charm_Bs",
                  "h_eff_chan0_Bsmm-2018-MC-charm_SM_Bs", "MC",
                  "h_eff_chan0_Run2018-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan1_2018-SM-charm_Bs",
                  "h_eff_chan1_Bsmm-2018-MC-charm_SM_Bs", "MC",
                  "h_eff_chan1_Run2018-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan0_2017-SM-charm_Bs",
                  "h_eff_chan0_Bsmm-2017-MC-charm_SM_Bs", "MC",
                  "h_eff_chan0_Run2017-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan1_2017-SM-charm_Bs",
                  "h_eff_chan1_Bsmm-2017-MC-charm_SM_Bs", "MC",
                  "h_eff_chan1_Run2017-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan0_2016GH-SM-charm_Bs",
                  "h_eff_chan0_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
                  "h_eff_chan0_Run2016GH-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan1_2016GH-SM-charm_Bs",
                  "h_eff_chan1_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
                  "h_eff_chan1_Run2016GH-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan0_2016BF-SM-charm_Bs",
                  "h_eff_chan0_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
                  "h_eff_chan0_Run2016BF-SM-charm_Bs", "Data")

make_overlay_plot("h_eff_chan1_2016BF-SM-charm_Bs",
                  "h_eff_chan1_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
                  "h_eff_chan1_Run2016BF-SM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2018-DM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2018-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2018-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2018-DM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2018-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2018-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2017-DM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2017-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2017-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2017-DM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2017-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2017-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2016GH-DM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2016GH-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2016GH-DM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2016GH-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2016GH-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan0_2016BF-DM-charm_Bs",
#                   "h_eff_chan0_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan0_Run2016BF-DM-charm_Bs", "Data")

# make_overlay_plot("h_eff_chan1_2016BF-DM-charm_Bs",
#                   "h_eff_chan1_Bsmm-2016BF-MC-charm_SM_Bs", "MC",
#                   "h_eff_chan1_Run2016BF-DM-charm_Bs", "Data")

