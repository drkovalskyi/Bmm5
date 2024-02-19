"""Trigger efficiency study

The script measures the absolute trigger efficiency of dimuon
B-physics triggers in MC and Data. For Data only datasets selected by
non-muon triggers can be used.

The input data should be in Bmm5 NanoAOD format.

The script is fairly slow especially for large MC samples. For testing
you may want to run on just a few files.

The results are stored in two two ways:
- json file with efficiency measurements
- histograms of PV distributions to track conditions of triggers with 
  dynamic prescales

The json file is used in an update mode, i.e. new results override old
results for specific studies keeping results from other studies
unaffected.

"""

import sys, os, subprocess, re, json
import ROOT
from math import *
from tdrstyle import *
from Bmm5.NanoAOD.selection import *
from collections import Counter
import numpy

# Modes
# - flat - use flat ntuples made for the study
# - nano - use NanoAOD

# mode = "flat"
mode = "nano"
recompute_results = False

path_skim1 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/526/trig/"
path_skim2 = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing-NEW/NanoAOD-skims/518/trig/"
path_flat  = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing-NEW/FlatNtuples/518/trig-info/"
path_nano  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/526/"

# split_channels = True
split_channels = False
output = "results/summary.json"
nbins = 60

if mode == "flat":
	path1 = path_flat
	path2 = path_flat
	path3 = path_flat
else:
	path1 = path_nano
	path2 = path_skim1
	path3 = path_skim2

results = dict()
if not os.path.exists('results'):
	os.mkdir('results')
if os.path.exists(output):
	results = json.load(open(output))
	
samples = {
	'Data ParkingDoubleMuonLowMass - Run2022C':[
		path2 + "/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass1+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass2+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass3+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass4+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass5+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass6+Run2022C-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass7+Run2022C-PromptReco-v1+MINIAOD/*root",
	],
	'Data ParkingDoubleMuonLowMass - Run2022D':[
		path2 + "/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass1+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass1+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass2+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass2+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass3+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass3+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass4+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass4+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass5+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass5+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass6+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass6+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass7+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass7+Run2022D-PromptReco-v2+MINIAOD/*root",
	],
	'Data ParkingDoubleMuonLowMass - Run2022E':[
		path2 + "/ParkingDoubleMuonLowMass0+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass1+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass2+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass3+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass4+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass5+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass6+Run2022E-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass7+Run2022E-PromptReco-v1+MINIAOD/*root",
	],
	'Data ParkingDoubleMuonLowMass - Run2022F':[
		path2 + "/ParkingDoubleMuonLowMass0+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass1+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass2+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass3+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass4+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass5+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass6+Run2022F-PromptReco-v1+MINIAOD/*root",
		path2 + "/ParkingDoubleMuonLowMass7+Run2022F-PromptReco-v1+MINIAOD/*root",
	],
	'Data ParkingDoubleMuonLowMass - Run2023C':[
        path2 + "ParkingDoubleMuonLowMass0+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass0+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass0+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass0+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass1+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass1+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass1+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass1+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass2+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass2+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass2+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass2+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass3+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass3+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass3+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass3+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass4+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass4+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass4+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass4+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass5+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass5+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass5+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass5+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass6+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass6+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass6+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass6+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass7+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass7+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass7+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass7+Run2023C-PromptReco-v4+MINIAOD/*root",
	],
	'Data ParkingDoubleMuonLowMass - Run2023D':[
        path2 + "ParkingDoubleMuonLowMass0+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass0+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass1+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass1+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass2+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass2+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass3+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass3+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass4+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass4+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass5+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass5+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass6+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass6+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass7+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "ParkingDoubleMuonLowMass7+Run2023D-PromptReco-v2+MINIAOD/*root",
	],
	'Data SingleMuon - Run2022C':[
		path2 + "/SingleMuon+Run2022C-PromptReco-v1+MINIAOD/*root",
	],
	'Data DoubleMuon - Run2022C':[
		path2 + "/DoubleMuon+Run2022C-PromptReco-v1+MINIAOD/*root",
	],
	'Data Muon - Run2022D':[
		path2 + "/Muon+Run2022D-PromptReco-v1+MINIAOD/*root",
		path2 + "/Muon+Run2022D-PromptReco-v2+MINIAOD/*root",
		path2 + "/Muon+Run2022D-PromptReco-v3+MINIAOD/*root",
	],
	'Data Muon - Run2022E':[
		path2 + "/Muon+Run2022E-PromptReco-v1+MINIAOD/*root",
	],
	'Data Muon - Run2022F':[
		path2 + "/Muon+Run2022F-PromptReco-v1+MINIAOD/*root",
	],
	'D0ToMuMu - 2022':[	
		path1 + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/*.root",
	],
	'Data Muon - Run2023C':[
        path2 + "Muon0+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "Muon0+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "Muon0+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "Muon0+Run2023C-PromptReco-v4+MINIAOD/*root",
        path2 + "Muon1+Run2023C-PromptReco-v1+MINIAOD/*root",
        path2 + "Muon1+Run2023C-PromptReco-v2+MINIAOD/*root",
        path2 + "Muon1+Run2023C-PromptReco-v3+MINIAOD/*root",
        path2 + "Muon1+Run2023C-PromptReco-v4+MINIAOD/*root",
	],
	'Data Muon - Run2023D':[
        path2 + "Muon0+Run2023D-PromptReco-v1+MINIAOD/*root",
        path2 + "Muon0+Run2023D-PromptReco-v2+MINIAOD/*root",
        path2 + "Muon1+Run2023D-PromptReco-v1+MINIAOD/*root",
	],
}
samples['Data ParkingDoubleMuonLowMass - Run2022'] = []
samples['Data ParkingDoubleMuonLowMass - Run2022'].extend(samples['Data ParkingDoubleMuonLowMass - Run2022C'])
samples['Data ParkingDoubleMuonLowMass - Run2022'].extend(samples['Data ParkingDoubleMuonLowMass - Run2022D'])
samples['Data ParkingDoubleMuonLowMass - Run2022'].extend(samples['Data ParkingDoubleMuonLowMass - Run2022E'])
samples['Data ParkingDoubleMuonLowMass - Run2022'].extend(samples['Data ParkingDoubleMuonLowMass - Run2022F'])
samples['Data ParkingDoubleMuonLowMass - Run2023'] = []
samples['Data ParkingDoubleMuonLowMass - Run2023'].extend(samples['Data ParkingDoubleMuonLowMass - Run2023C'])
samples['Data ParkingDoubleMuonLowMass - Run2023'].extend(samples['Data ParkingDoubleMuonLowMass - Run2023D'])
samples['Data Muon - Run2022'] = []
samples['Data Muon - Run2022'].extend(samples['Data Muon - Run2022D'])
samples['Data Muon - Run2022'].extend(samples['Data Muon - Run2022E'])
samples['Data Muon - Run2022'].extend(samples['Data Muon - Run2022F'])
samples['Data Muon - Run2023'] = []
samples['Data Muon - Run2023'].extend(samples['Data Muon - Run2023C'])
samples['Data Muon - Run2023'].extend(samples['Data Muon - Run2023D'])

cuts = dict()

if mode == "nano":
	cuts['dimuon'] = 'mm_mu1_index>=0 && mm_mu2_index>=0' + \
		'&& Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45' + \
		'&& Muon_pt[mm_mu1_index]>4.0 && Muon_pt[mm_mu2_index]>4.0' + \
		'&& Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1' + \
		'&& mm_kin_vtx_prob>0.025'
	
	cuts['d0mm_loose'] = cuts['dimuon'] + '&& mm_kin_mass>1.50 && mm_kin_mass<2.50'
	cuts['d0mm']       = cuts['dimuon'] + '&& mm_kin_mass>1.76 && mm_kin_mass<1.94'
	
else:
	cuts['HLT_DoubleMu4_3_Bs'] = "good_for_HLT_DoubleMu4_3_Bs"
	cuts['HLT_DoubleMu4_3_Jpsi'] = "good_for_HLT_DoubleMu4_3_Jpsi"
	cuts['HLT_DoubleMu4_3_Jpsi_Displaced'] = "good_for_HLT_DoubleMu4_3_Jpsi_Displaced"

studies = {
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2022C':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2022C',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2022D':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2022D',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2022E':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2022E',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2022F':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2022F',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2022':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2022',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - 2022 MC':{
	#  	'samples':'D0ToMuMu - 2022',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	'HLT_DoubleMu4_3_LowMass - 2022 MC':{
	 	'samples':'D0ToMuMu - 2022',
	 	'trigger':'HLT_DoubleMu4_3_LowMass',
	 	'cut':cuts['d0mm'],
	},


	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2023C':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2023C',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2023D':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2023D',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },
	# 'HLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu - Run2023':{
	#  	'samples':'Data ParkingDoubleMuonLowMass - Run2023',
	# 	'trigger':'HLT_DoubleMu4_3_LowMass',
	# 	'cut':cuts['d0mm'] + ' && HLT_Mu4_L1DoubleMu',
	# },


	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2022C':{
	#  	'samples':'Data SingleMuon - Run2022C',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2022D':{
	#  	'samples':'Data Muon - Run2022D',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2022E':{
	#  	'samples':'Data Muon - Run2022E',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2022F':{
	#  	'samples':'Data Muon - Run2022F',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2022':{
	#  	'samples':'Data Muon - Run2022',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - 2022 MC':{
	#  	'samples':'D0ToMuMu - 2022',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	'L1 - 2022 MC':{
	 	'samples':'D0ToMuMu - 2022',
		'trigger':'(L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4 || L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 || L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 || L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 || L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 || L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4p5_SQ_OS_dR_Max1p2 || L1_DoubleMu4_SQ_OS_dR_Max1p2)',
		'cut':cuts['d0mm_loose']
	},
	'HLT_Mu0_L1DoubleMu wrt L1 - 2022 MC':{
	 	'samples':'D0ToMuMu - 2022',
		'trigger':'HLT_Mu0_L1DoubleMu',
		'cut':cuts['d0mm_loose'] + ' && (L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4 || L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 || L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 || L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 || L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 || L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 || L1_DoubleMu4p5_SQ_OS_dR_Max1p2 || L1_DoubleMu4_SQ_OS_dR_Max1p2)'
	},
	
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2023C':{
	#  	'samples':'Data Muon - Run2023C',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2023D':{
	#  	'samples':'Data Muon - Run2023D',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu3_PFJet40 - Run2023':{
	#  	'samples':'Data Muon - Run2023',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu3_PFJet40',
	# },
	
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2022C':{
	#  	'samples':'Data DoubleMuon - Run2022C',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2022D':{
	#  	'samples':'Data Muon - Run2022D',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2022E':{
	#  	'samples':'Data Muon - Run2022E',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2022F':{
	#  	'samples':'Data Muon - Run2022F',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2022':{
	#  	'samples':'Data Muon - Run2022',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - 2022 MC':{
	#  	'samples':'D0ToMuMu - 2022',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu8 - 2022 MC':{
	#  	'samples':'D0ToMuMu - 2022',
	# 	'trigger':'HLT_Mu8',
	# 	'cut':cuts['d0mm_loose'],
	# },

	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2023C':{
	#  	'samples':'Data Muon - Run2023C',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2023D':{
	#  	'samples':'Data Muon - Run2023D',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },
	# 'HLT_Mu0_L1DoubleMu wrt HLT_Mu8 - Run2023':{
	#  	'samples':'Data Muon - Run2023',
	# 	'trigger':'HLT_Mu0_L1DoubleMu',
	# 	'cut':cuts['d0mm_loose'] + ' && HLT_Mu8',
	# },

	# # ======= HLT_DoubleMu4_3_Bs ==========
	
	# 'HLT_DoubleMu4_3_Bs - 2016BF':{
	# 	'samples':'Data - 2016BF',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'HLT_DoubleMu4_3_Bs - 2016GH':{
	# 	'samples':'Data - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'HLT_DoubleMu4_3_Bs - 2017':{
	# 	'samples':'Data - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'HLT_DoubleMu4_3_Bs - 2018':{
	# 	'samples':'Data - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'HLT_DoubleMu4_3_Bs - 2018 MC':{
	# 	'samples':'BsToMuMu - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'HLT_DoubleMu4_3_Bs - 2017 MC':{
	# 	'samples':'BsToMuMu - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'HLT_DoubleMu4_3_Bs - 2016BF MC':{
	# 	'samples':'BsToMuMu - 2016BF',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'HLT_DoubleMu4_3_Bs - 2016GH MC':{
	# 	'samples':'BsToMuMu - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# # ======= HLT_DoubleMu4_3_Jpsi ==========

	# 'HLT_DoubleMu4_3_Jpsi - 2018 MC':{
	# 	'samples':'BuToJpsiK - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'efficiency': ['trigger_object_efficiency_Bsmm-2018-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2018-SM-charm_Bs'],
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi - 2018':{
	# 	'samples':'Data - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi'],
	# },

	# # ======= HLT_DoubleMu4_3_Jpsi_Displaced ==========

	# 'HLT_DoubleMu4_3_Jpsi_Displaced - 2017':{
	# 	'samples':'Data - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced - 2016BF':{
	# 	'samples':'Data - 2016BF',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced - 2016GH':{
	# 	'samples':'Data - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced - 2017 MC':{
	# 	'samples':'BuToJpsiK - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'efficiency': ['trigger_object_efficiency_Bsmm-2017-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2017-SM-charm_Bs'],
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced - 2016BF MC':{
	# 	'samples':'BuToJpsiK - 2016BF',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'efficiency': ['trigger_object_efficiency_Bsmm-2016BF-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2016BF-SM-charm_Bs'],
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'] + "&& L1_DoubleMu0er1p6_dEta_Max1p8_OS && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced - 2016GH MC':{
	# 	'samples':'BuToJpsiK - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'efficiency': ['trigger_object_efficiency_Bsmm-2016GH-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2016GH-SM-charm_Bs'],
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'] + "&& L1_DoubleMu0er1p6_dEta_Max1p8_OS && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# },

	# ======= HLT_DoubleMu4_Jpsi_NoVertexing ==========

	# 'HLT_DoubleMu4_Jpsi_NoVertexing - 2018':{
	# 	'samples':'Charmonium - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_DoubleMu4_Jpsi_NoVertexing && ' + cuts['HLT_DoubleMu4_3_Jpsi'],
	# },

	# 'HLT_DoubleMu4_Jpsi_NoVertexing - 2018 MC':{
	# 	'samples':'BuToJpsiK - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_DoubleMu4_Jpsi_NoVertexing && ' + cuts['HLT_DoubleMu4_3_Jpsi'],
	# },
	
	# 'HLT_Dimuon0_LowMass - Jpsi - 2018':{
	# 	'samples':'Charmonium - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_Dimuon0_LowMass && ' + cuts['HLT_DoubleMu4_3_Jpsi'],
	# },

	# 'HLT_Dimuon0_LowMass - Bs - 2018':{
	# 	'samples':'Charmonium - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':'HLT_Dimuon0_LowMass && ' + cuts['HLT_DoubleMu4_3_Bs'],
	# },
	
	# 'HLT_DoubleMu4_Jpsi_NoVertexing - 2017':{
	# 	'samples':'Data - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_DoubleMu4_Jpsi_NoVertexing && ' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },

	# 'HLT_DoubleMu4_3_Jpsi wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Jpsi - 2018':{
	# 	'samples':'Charmonium - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Jpsi'],
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Jpsi - 2018 MC':{
	# 	'samples':'BuToJpsiK - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Jpsi'],
	# },

	# 'HLT_DoubleMu4_3_Jpsi wrt HLT_Dimuon0_LowMass_L1_0er1p5 displaced trig eff for Jpsi - 2018':{
	# 	'samples':'Charmonium - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi wrt HLT_Dimuon0_LowMass_L1_0er1p5 displaced trig eff for Jpsi - 2018 MC':{
	# 	'samples':'BuToJpsiK - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },

	
	# 'HLT_DoubleMu4_3_Bs wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Bs - 2018':{
	# 	'samples':'Charmonium - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Bs'],
	# },
	
	# 'HLT_DoubleMu4_3_Bs wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Bs - 2018 MC':{
	# 	'samples':'BsToMuMu - 2018',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Jpsi - 2017':{
	# 	'samples':'Charmonium - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Jpsi - 2017 MC':{
	# 	'samples':'BuToJpsiK - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	
	# 'HLT_DoubleMu4_3_Bs wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Bs - 2017':{
	# 	'samples':'Charmonium - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Bs'],
	# },
	
	# 'HLT_DoubleMu4_3_Bs wrt HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Bs - 2017 MC':{
	# 	'samples':'BsToMuMu - 2017',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 &&' + cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing trig eff for Jpsi - 2016BF':{
	# 	'samples':'Charmonium - 2016BF',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing trig eff for Jpsi - 2016BF MC':{
	# 	'samples':'BuToJpsiK - 2016BF',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing trig eff for Jpsi - 2016GH':{
	# 	'samples':'Charmonium - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing trig eff for Jpsi - 2016GH MC':{
	# 	'samples':'BuToJpsiK - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_DoubleMu0 trig eff for Jpsi - 2016H':{
	# 	'samples':'DoubleMuon - 2016H',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':' HLT_DoubleMu0 && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	
	# 'HLT_DoubleMu4_3_Jpsi_Displaced wrt HLT_DoubleMu0 trig eff for Jpsi - 2016H MC':{
	# 	'samples':'BuToJpsiK - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	# 	'cut':'HLT_DoubleMu0 && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'HLT_DoubleMu4_3_Bs wrt HLT_DoubleMu0 trig eff for Bs - 2016H':{
	# 	'samples':'DoubleMuon - 2016H',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':' HLT_DoubleMu0 && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Bs'],
	# },
	
	# 'HLT_DoubleMu4_3_Bs wrt HLT_DoubleMu0 trig eff for Bs - 2016H MC':{
	# 	'samples':'BsToMuMu - 2016GH',
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'cut':'HLT_DoubleMu0 && L1_DoubleMu0er1p6_dEta_Max1p8 &&' + cuts['HLT_DoubleMu4_3_Bs'],
	# },
	
	# # 'HLT_Dimuon0_LowMass_L1_0er1p5 trig eff for Bs - 2018 MC':{
	# # 	'samples':'BsToMuMu - 2018',
	# # 	'trigger':'HLT_Dimuon0_LowMass_L1_0er1p5',
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# # },
	
	# # 'L1_DoubleMu0er1p5_SQ trig eff for Jpsi - 2018 MC 2':{
	# # 	'samples':'BuToJpsiK - 2018',
	# # 	'trigger':'L1_DoubleMu0er1p5_SQ',
	# # 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Match'],
	# # },

	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Bs - 2018 MC':{
	# 	'samples':'BsToMuMu - 2018',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Bs - 2018':{
	# 	'samples':'Data - 2018',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Jpsi - 2018 MC':{
	# 	'samples':'BuToJpsiK - 2018',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi'],
	# },
	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Jpsi - 2018':{
	# 	'samples':'Data - 2018',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi'],
	# },

	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Bs - 2017 MC':{
	# 	'samples':'BsToMuMu - 2017',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Bs - 2017':{
	# 	'samples':'Data - 2017',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Jpsi - 2017 MC':{
	# 	'samples':'BuToJpsiK - 2017',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 trig eff for Jpsi - 2017':{
	# 	'samples':'Data - 2017',
	# 	'trigger':'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },


	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Bs - 2016BF MC':{
	# 	'samples':'BsToMuMu - 2016BF',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Bs - 2016BF':{
	# 	'samples':'Data - 2016BF',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Jpsi - 2016BF MC':{
	# 	'samples':'BuToJpsiK - 2016BF',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Jpsi - 2016BF':{
	# 	'samples':'Data - 2016BF',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },

	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Bs - 2016GH MC':{
	# 	'samples':'BsToMuMu - 2016GH',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },
	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Bs - 2016GH':{
	# 	'samples':'Data - 2016GH',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Bs'],
	# },

	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Jpsi - 2016GH MC':{
	# 	'samples':'BuToJpsiK - 2016GH',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	# 'L1_DoubleMu0er1p6_dEta_Max1p8_OS trig eff for Jpsi - 2016GH':{
	# 	'samples':'Data - 2016GH',
	# 	'trigger':'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
	# 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Displaced'],
	# },
	

	# # 'L1_DoubleMu0er1p5_SQ_OS trig eff for Jpsi - 2018 MC 2':{
	# # 	'samples':'BuToJpsiK - 2018',
	# # 	'trigger':'L1_DoubleMu0er1p5_SQ_OS',
	# # 	'cut':cuts['HLT_DoubleMu4_3_Jpsi_Match'],
	# # },

	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2016BF MC':{
	# # 	'samples':'BsToMuMu - 2016BF',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2016BF-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2016BF-SM-charm_Bs'],
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p6_dEta_Max1p8_OS && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },
	
	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2016GH MC':{
	# # 	'samples':'BsToMuMu - 2016GH',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2016GH-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2016GH-SM-charm_Bs'],
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p6_dEta_Max1p8_OS && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },
	
	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2017 MC':{
	# # 	'samples':'BsToMuMu - 2017',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2017-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2017-SM-charm_Bs'],
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },
	
	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2018 MC':{
	# # 	'samples':'BsToMuMu - 2018',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2018-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2018-SM-charm_Bs'],
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },

	
	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2018':{
	# # 	'samples':'Data - 2018',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2018-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2018-SM-charm_Bs'],
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },

	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2017':{
	# # 	'samples':'Data - 2017',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2017-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2017-SM-charm_Bs'],
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },
	
	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2016GH':{
	# # 	'samples':'Data - 2016GH',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2016GH-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2016GH-SM-charm_Bs'],
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p6_dEta_Max1p8_OS && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },
	
	# # 'HLT_DoubleMu4_3_Bs trig eff wrt L1 - 2016BF':{
	# # 	'samples':'Data - 2016BF',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'efficiency': ['trigger_object_efficiency_Bsmm-2016BF-MC-charm_SM_Bs', 'trigger_object_efficiency_Run2016BF-SM-charm_Bs'],
	# # 	'cut':cuts['HLT_DoubleMu4_3_Bs'] + "&& L1_DoubleMu0er1p6_dEta_Max1p8_OS && nl1==2 && m1_l1_pt>0 && m2_l1_pt>0 && m1pt<20 && m2pt<20",
	# # },
	
	
	# # 'HLT_Dimuon0_LowMass_L1_0er1p5 - Bs - 2018':{
	# # 	'samples':'Charmonium - 2018',
	# # 	'trigger':'HLT_DoubleMu4_3_Bs',
	# # 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && ' + cuts['HLT_DoubleMu4_3_Bs'],
	# # },

	# # 'HLT_Dimuon0_LowMass_L1_0er1p5 - Jpsi Loose Vtx - 2018':{
	# # 	'samples':'Charmonium - 2018',
	# # 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# # 	'cut':'HLT_Dimuon0_LowMass_L1_0er1p5 && ' + cuts['HLT_DoubleMu4_3_Jpsi_LooseVtx'],
	# # },
	
}



output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/trigger_efficiency_run3"

ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 800,600)


def load_data(name):
	if mode == "nano":
		chain = ROOT.TChain("Events")
	else:
		chain = ROOT.TChain("mm")
	for sample in samples[name]:
		chain.Add(sample)

	print("Total number of events:", chain.GetEntries())
	return chain

def apply_trigger_object_efficiency2D(chain, selection, filename):
	if mode != "flat": return
	
	f_eff = ROOT.TFile("results/%s.root" % filename)

	cut = convert(chain, selection)

	# print cut

	print("using", filename)
	h_name = re.sub("trigger_object_efficiency", "h2_eff", filename)
	h2 = f_eff.Get(h_name)

	n = 0
	n_weighted = 0
	events = dict()
	for event in chain:
		cut = cut.format(tree='event')
		if not eval(cut): continue
		
		evt = "%u-%u" % (event.run, event.evt)
		if evt in events: continue
		n += 1
		
		events[evt] = 1
			
		mu1_pt = event.m1pt
		if mu1_pt >= 20:
			mu1_pt = 19.99
		mu2_pt = event.m2pt
		if mu2_pt >= 20:
			mu2_pt = 19.99
		
		mu1_eff = h2.GetBinContent(h2.FindBin(mu1_pt, event.m1eta))
		mu2_eff = h2.GetBinContent(h2.FindBin(mu2_pt, event.m2eta))

		n_weighted += mu1_eff * mu2_eff
		# if filename == 'trigger_object_efficiency_Run2017':
		# 	print "mu1: %4.1f %4.1f eff: %4.2f mu2: %4.1f %4.1f eff: %4.2f" % (mu1_pt, event.m1eta, mu1_eff, mu2_pt, event.m2eta, mu2_eff)
	
	print("Selected candidates:", n)
	print("Efficiency: %0.1f" % (100. * n_weighted/ n))
	# print "Candidate multiplicity: %0.2f" % (float(n)/len(events))

def apply_trigger_object_efficiency(chain, selection, filename):
	if mode != "flat": return
	
	f_eff = ROOT.TFile("results/%s.root" % filename)

	cut = convert(chain, selection)

	# print cut

	print("using", filename)
	name = re.sub("trigger_object_efficiency_", "", filename)
	h_chan0 = f_eff.Get("h_eff_chan0_" + name)
	h_chan1 = f_eff.Get("h_eff_chan1_" + name)

	n = 0
	n_weighted = 0
	nbins = h_chan0.GetNbinsX()
	event_counts = numpy.zeros((2*nbins, 2*nbins))
	events = dict()
	for event in chain:
		cut = cut.format(tree='event')
		if not eval(cut): continue
		
		evt = "%u-%u" % (event.run, event.evt)
		if evt in events: continue
		n += 1
		
		events[evt] = 1
			
		mu1_pt = event.m1pt
		if mu1_pt >= 20:
			mu1_pt = 19.99
		mu2_pt = event.m2pt
		if mu2_pt >= 20:
			mu2_pt = 19.99

		if abs(event.m1eta) < 0.7:
			mu1_ch  = 0
			mu1_bin = h_chan0.FindBin(mu1_pt)
			mu1_eff = h_chan0.GetBinContent(mu1_bin)
		else:
			mu1_ch  = 1
			mu1_bin = h_chan1.FindBin(mu1_pt)
			mu1_eff = h_chan1.GetBinContent(mu1_bin)
			
		if abs(event.m2eta) < 0.7:
			mu2_ch  = 0
			mu2_bin = h_chan0.FindBin(mu2_pt)
			mu2_eff = h_chan0.GetBinContent(mu2_bin)
		else:
			mu2_ch  = 1
			mu2_bin = h_chan1.FindBin(mu2_pt)
			mu2_eff = h_chan1.GetBinContent(mu2_bin)

		event_counts[mu1_ch * nbins + mu1_bin - 1][mu2_ch * nbins + mu2_bin - 1] += 1

		n_weighted += mu1_eff * mu2_eff
		
		# if filename == 'trigger_object_efficiency_Run2017':
		# 	print "mu1: %4.1f %4.1f eff: %4.2f mu2: %4.1f %4.1f eff: %4.2f" % (mu1_pt, event.m1eta, mu1_eff, mu2_pt, event.m2eta, mu2_eff)
	
	print("Selected candidates:", n)
	
	sum_rel_err2 = 0
	h1 = None
	h2 = None
	for mu1 in range(2 * nbins):
		if mu1 < nbins:
			h1 = h_chan0
			m1 = mu1
		else:
			h1 = h_chan1
			m1 = mu1 - nbins
			
		for mu2 in range(2 * nbins):
			# print "%2u %2u: %u" % (mu1, mu2, event_counts[mu1][mu2])
			if event_counts[mu1][mu2] == 0: continue
			
			if mu2 < nbins:
				h2 = h_chan0
				m2 = mu2
			else:
				h2 = h_chan1
				m2 = mu2 - nbins
				
			if mu1 == mu2:
				sum_rel_err2 += pow(2 * h1.GetBinError(m1 + 1) / h1.GetBinContent(m1 + 1), 2)
			else:
				sum_rel_err2 += pow(h1.GetBinError(m1 + 1) / h1.GetBinContent(m1 + 1), 2)
				sum_rel_err2 += pow(h2.GetBinError(m2 + 1) / h2.GetBinContent(m2 + 1), 2)
			
	print("Efficiency: %0.1f \pm %0.1f" % (100. * n_weighted/ n,
										   100. * n_weighted/ n * sqrt(sum_rel_err2 / n))) 
	# print "Candidate multiplicity: %0.2f" % (float(n)/len(events))


	
for name, info in sorted(studies.items()):
	print("\nProcessing", name)
	
	chain = None
	trigger = info['trigger']

	for ch in range(2):
		file_name = re.sub('\s+', '', name)
		study_name = name
		
		cut = info['cut']

		if split_channels:
			print("Channel:", ch)
			file_name += "_%u" % ch
			study_name += " (channel: %u)" % ch
			if mode == "nano":
				if ch == 0:
					cut += "&& abs(Muon_eta[mm_mu1_index])<0.7 && abs(Muon_eta[mm_mu2_index])<0.7"
				else:
					cut += "&& (abs(Muon_eta[mm_mu1_index])>0.7 || abs(Muon_eta[mm_mu2_index])>0.7)"
			else:
				cut += "&& chan==%u" % ch

		if study_name in results and not recompute_results:
			print("Results are already available. Skip the study")
			print("%s efficiency: %0.2f \pm %0.2f %%" % (trigger, 100. * results[study_name]['eff'], 100. * results[study_name]['eff_err']))
			if not split_channels:
				break
			else:
				continue

		if chain == None:
			chain = load_data(info['samples'])

		prescale = ""
		if hasattr(chain, 'prescale_%s' % trigger):
			prescale = "*prescale_%s" % trigger

		if mode == "nano":
			f = ROOT.TFile.Open('results/' + file_name + ".root", "recreate")
			h_off = ROOT.TH1F("h_off","h_off", nbins, 0, nbins)
			h_off.Sumw2()
			chain.Draw("PV_npvsGood>>h_off", cut)
			h_off.Write()
			h_off.SetDirectory(0)

			h_off_trig = ROOT.TH1F("h_off_trig","h_off_trig", nbins, 0, nbins)
			h_off_trig.Sumw2()

			chain.Draw("PV_npvsGood>>h_off_trig", "(%s)%s" % (cut + "&&" + trigger, prescale))
			h_off_trig.Write()
			h_off_trig.SetDirectory(0)
			f.Close()
		else:
			h_off = ROOT.TH1F("h_off","h_off", nbins, 0, 100)
			h_off.Sumw2()
			chain.Draw("pt>>h_off", cut)

			h_off_trig = h_off.Clone("h_off_trig")

			chain.Draw("pt>>h_off_trig", "(%s)%s" % (cut + "&&" + trigger, prescale))

			h_off_trig2 = h_off.Clone("h_off_trig2")

			chain.Draw("pt>>h_off_trig2", cut + "&& m1_hlt_pt>0 && m2_hlt_pt>0")
			
			h_off_trig3 = h_off.Clone("h_off_trig3")

			chain.Draw("pt>>h_off_trig3", "(%s)%s" % (cut + "&& m1_hlt_pt>0 && m2_hlt_pt>0 &&" + trigger, prescale))

		if h_off_trig.Integral(0,nbins + 1) > 2 * h_off_trig.GetEntries():
			# prescaled case - use standard error estimation
			
			err_off = ROOT.Double(0)
			n_off = h_off.IntegralAndError(0, nbins + 1, err_off)

			err_off_trig = ROOT.Double(0)
			n_off_trig = h_off_trig.IntegralAndError(0, nbins + 1, err_off_trig)

			eff = float(n_off_trig) / n_off
			eff_err = eff * sqrt((err_off/n_off)**2 + (err_off_trig/n_off_trig)**2)
		else:
			# unprescaled case - use binomial error estimation

			n_off = h_off.Integral(0, nbins + 1)
			n_off_trig = h_off_trig.Integral(0, nbins + 1)

			eff = float(n_off_trig) / n_off
			eff_err = sqrt(eff * (1 - eff) / n_off)

		print("n_off:", n_off)
		print("%s efficiency: %0.2f \pm %0.2f %%" % (trigger, 100. * eff, 100. * eff_err))

		# err_off_trig2 = ROOT.Double(0)
		# n_off_trig2 = h_off_trig2.IntegralAndError(0, nbins + 1, err_off_trig2)

		# err_off_trig3 = ROOT.Double(0)
		# n_off_trig3 = h_off_trig3.IntegralAndError(0, nbins + 1, err_off_trig3)

		# eff2 = float(n_off_trig2) / n_off
		# eff_err2 = eff2 * sqrt((err_off/n_off)**2 + (err_off_trig2/n_off_trig)**2)
		# print "%s efficiency: %0.1f +/- %0.1f %%" % ("m1_hlt_pt>0 && m2_hlt_pt>0",
		# 											 100. * eff2,
		# 											 100. * eff_err2)

		# eff3 = float(n_off_trig3) / n_off
		# eff_err3 = eff3 * sqrt((err_off/n_off)**2 + (err_off_trig3/n_off_trig)**2)
		# print "%s efficiency: %0.1f +/- %0.1f %%" % ("m1_hlt_pt>0 && m2_hlt_pt>0 && " + trigger,
		# 											 100. * eff3,
		# 											 100. * eff_err3)

		if 'efficiency' in info:
			for to_eff in info['efficiency']:
				apply_trigger_object_efficiency(chain, cut, to_eff)

		results[study_name] = {
			'eff':eff,
			'eff_err':eff_err
		}
			
		if not split_channels:
			break

## Save results
json.dump(results, open(output, 'w'))

sys.exit()

## Produce report
print("\nSummary")
for study, info in sorted(results.items()):
	print("%s \t: %0.1f +/- %0.1f %%" % (study, 100. * info['eff'], 100. * info['eff_err']))

## Ratios
def compute_ratio(bmm, jpsik):
	bmm_rel_err = results[bmm]['eff_err']/results[bmm]['eff']
	jpsik_rel_err = results[jpsik]['eff_err']/results[jpsik]['eff']
	ratio = results[jpsik]['eff']/results[bmm]['eff']
	ratio_err = sqrt(bmm_rel_err * bmm_rel_err + jpsik_rel_err * jpsik_rel_err) * ratio
	print("%s / %s: \t%0.2f\pm%0.2f" % (jpsik, bmm, ratio, ratio_err))

print("\nRatios")
for bmm in sorted(results):
	match = re.search("HLT_DoubleMu4_3_Bs(.*)$", bmm)
	if not match: continue
	suffix = match.group(1)
	for jpsik in results:
		if jpsik == bmm:
			continue
		if not re.search("^\S+%s$" % re.escape(suffix), jpsik):
			continue
		compute_ratio(bmm, jpsik)
	
	

	
## Make plots
sys.exit()

setTDRStyle()
c1 = ROOT.TCanvas("c1", "c1", 800, 800)
def make_plot(name, filename_mc, filename_data, plots, rebin=1):
	f_mc = ROOT.TFile.Open(filename_mc)
	f_data = ROOT.TFile.Open(filename_data)
	normalization = None
	legend = ROOT.TLegend(0.70,0.75,0.85,0.87)
	legend.SetShadowColor(ROOT.kWhite)
	legend.SetLineColor(ROOT.kWhite)
	# legend.SetFillColor(10)
	colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kOrange+5, ROOT.kGreen+3]
	scale = 1.2
	for idx, plot in enumerate(plots):
		if plot['type'] == "mc":
			f = f_mc
		else:
			f = f_data
		h = f.Get(plot['name'])
		if rebin > 1:
			h.Rebin(rebin)
		h.SetLineColor(colors[idx])
		h.SetLineWidth(2)
		h.SetMarkerStyle(20)
		h.SetMarkerColor(colors[idx])
		h.GetXaxis().SetTitle("nPV")
		if normalization:
			h.Scale(normalization/h.Integral())
			h.Draw("same")
		else:
			normalization = h.Integral()
			h.SetMaximum(h.GetMaximum()*scale)
			h.Draw()
		legend.AddEntry(h, plot['label'])
	legend.Draw()
	c1.Print(name)
	

make_plot("HLT_DoubleMu4_3_Bs-2018.pdf", "HLT_DoubleMu4_3_Bs-2018MC.root", "HLT_DoubleMu4_3_Bs-2018.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Bs-2017.pdf", "HLT_DoubleMu4_3_Bs-2017MC.root", "HLT_DoubleMu4_3_Bs-2017.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Bs-2016BF.pdf", "HLT_DoubleMu4_3_Bs-2016BFMC.root", "HLT_DoubleMu4_3_Bs-2016BF.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Bs-2016GH.pdf", "HLT_DoubleMu4_3_Bs-2016GHMC.root", "HLT_DoubleMu4_3_Bs-2016GH.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi-2018.pdf", "HLT_DoubleMu4_3_Jpsi-2018MC.root", "HLT_DoubleMu4_3_Jpsi-2018.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2017.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2017MC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2017.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016BF.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BFMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BF.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016GH.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GHMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GH.root", 
		  plots = [
			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
			  {'type':'data','name':'h_off_trig','label':'Data'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2017-Data.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2017MC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2017.root", 
		  plots = [
			  {'type':'Data', 'name':'h_off', 'label':'All'}, 
			  {'type':'Data','name':'h_off_trig','label':'Triggered'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016BF-Data.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BFMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BF.root", 
		  plots = [
			  {'type':'Data', 'name':'h_off', 'label':'All'}, 
			  {'type':'Data','name':'h_off_trig','label':'Triggered'}
		  ],
		  rebin = 3
)
make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016GH-Data.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GHMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GH.root", 
		  plots = [
			  {'type':'Data', 'name':'h_off', 'label':'All'}, 
			  {'type':'Data','name':'h_off_trig','label':'Triggered'}
		  ],
		  rebin = 3
)
		  
		
# Local Variables:
# indent-tabs-mode: 1
# tab-width: 4
# python-indent: 4
# End:
