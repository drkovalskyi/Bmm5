import sys, os, subprocess
from ROOT import TChain, TFile, TTree, TH1, TROOT, TDirectory, TPad, TCanvas, TColor
from math import *

path_skim = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/NanoAOD-skims/513/mm/"
path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/513/"
studies = {
	'HLT_DoubleMu4_3_Bs - 2016 APV problem':{
		'samples':[
            path_skim + "/SingleElectron+Run2016B-17Jul2018_ver1-v1+MINIAOD/",
            path_skim + "/SingleElectron+Run2016B-17Jul2018_ver2-v1+MINIAOD/",
            path_skim + "/SingleElectron+Run2016C-17Jul2018-v1+MINIAOD/",
            path_skim + "/SingleElectron+Run2016D-17Jul2018-v1+MINIAOD/",
            path_skim + "/SingleElectron+Run2016E-17Jul2018-v1+MINIAOD/",
            path_skim + "/SingleElectron+Run2016F-17Jul2018-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016B-17Jul2018_ver1-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016B-17Jul2018_ver2-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016C-17Jul2018-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016D-17Jul2018-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016E-17Jul2018-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016F-17Jul2018-v1+MINIAOD/",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
		'prescale':1.0
	},
	'HLT_DoubleMu4_3_Bs - 2016 fixed':{
		'samples':[
            path_skim + "/SingleElectron+Run2016G-17Jul2018-v1+MINIAOD/",
            path_skim + "/SingleElectron+Run2016H-17Jul2018-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016G-17Jul2018-v1+MINIAOD/",
            path_skim + "/SinglePhoton+Run2016H-17Jul2018-v1+MINIAOD/",			
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
		'prescale':1.0
	},

	# 'HLT_DoubleMu4_3_Jpsi - 2018 MC':{
	# 	'samples':[
	# 		path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs - 2018 MC':{
	# 	'samples':[
	# 		path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs - 2017 MC':{
	# 	'samples':[
	# 		path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs - 2016 MC':{
	# 	'samples':[
	# 		path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1+MINIAODSIM/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs - 2018':{
	# 	'samples':[
	# 		path_skim + "/EGamma+Run2018A-17Sep2018-v2+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018B-17Sep2018-v1+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018C-17Sep2018-v1+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018D-22Jan2019-v2+MINIAOD/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs - 2017':{
	# 	'samples':[
    #         path_skim + "/DoubleEG+Run2017B-31Mar2018-v1+MINIAOD",
    #         path_skim + "/DoubleEG+Run2017C-31Mar2018-v1+MINIAOD",
    #         path_skim + "/DoubleEG+Run2017D-31Mar2018-v1+MINIAOD",
    #         path_skim + "/DoubleEG+Run2017E-31Mar2018-v1+MINIAOD",
    #         path_skim + "/DoubleEG+Run2017F-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SingleElectron+Run2017B-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SingleElectron+Run2017C-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SingleElectron+Run2017D-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SingleElectron+Run2017E-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SingleElectron+Run2017F-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SinglePhoton+Run2017B-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SinglePhoton+Run2017C-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SinglePhoton+Run2017D-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SinglePhoton+Run2017E-31Mar2018-v1+MINIAOD",
    #         path_skim + "/SinglePhoton+Run2017F-31Mar2018-v1+MINIAOD",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs - 2016':{
	# 	'samples':[
    #         path_skim + "/SingleElectron+Run2016B-17Jul2018_ver1-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016B-17Jul2018_ver2-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016C-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016D-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016E-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016F-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016G-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SingleElectron+Run2016H-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016B-17Jul2018_ver1-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016B-17Jul2018_ver2-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016C-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016D-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016E-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016F-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016G-17Jul2018-v1+MINIAOD/",
    #         path_skim + "/SinglePhoton+Run2016H-17Jul2018-v1+MINIAOD/",			
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Bs_v14':{
	# 	'samples':[
	# 		path_skim + "/EGamma+Run2018B-17Sep2018-v1+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018C-17Sep2018-v1+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018D-22Jan2019-v2+MINIAOD/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Bs',
	# 	'prescale':1.0
	# },
	# 'HLT_DoubleMu4_3_Jpsi - 2018':{
	# 	'samples':[
	#  		path_skim + "/EGamma+Run2018A-17Sep2018-v2+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018B-17Sep2018-v1+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018C-17Sep2018-v1+MINIAOD/",
	# 		path_skim + "/EGamma+Run2018D-22Jan2019-v2+MINIAOD/",
	# 	],
	# 	'trigger':'HLT_DoubleMu4_3_Jpsi',
	# 	'prescale':10.0
	# },
}

triggers = {
	'HLT_DoubleMu4_3_Bs':{
		'cuts':
		# 'sqrt(pow(abs(mm_mu2_eta-mm_mu1_eta),2) + pow(acos(cos(mm_mu2_phi-mm_mu1_phi)),2)) &&' +
		'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
		'mm_kin_pt>5 && mm_kin_vtx_prob>0.1 && Muon_softMvaId[mm_mu1_index] && Muon_softMvaId[mm_mu2_index] &&' +
		'abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 && Muon_pt[mm_mu1_index]>4.0 &&' +
		'Muon_pt[mm_mu2_index]>4.0 && mm_kin_mass>4.6 && mm_kin_mass<5.9 && Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1'
	},
	'HLT_DoubleMu4_3_Jpsi':{
		'cuts':
		'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
		'cos(mm_kin_alphaBS)>0.9 && mm_kin_pt>7 && mm_kin_vtx_prob>0.1 && Muon_softMvaId[mm_mu1_index] && ' +
		'Muon_softMvaId[mm_mu2_index] && abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 &&' +
		'Muon_pt[mm_mu1_index]>4.0 && Muon_pt[mm_mu2_index]>4.0 && abs(mm_kin_mass-3.1)<0.1 &&' +
		'Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1'
	}
}

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-513/trigger_efficiency"

# ROOT.gROOT.SetBatch(True)

for name, info in studies.items():
	print "Processing", name
	
	chain = TChain("Events")
	for sample in info['samples']:
		files = subprocess.check_output("find %s/ -type f -name '*.root'" % (sample), shell=True).split("\n")
		for f in files:
			if f != "":
				chain.Add(f)
				# break # use just the first file

	print "Total number of events:", chain.GetEntries()

	trigger = info['trigger']

	n_off = chain.GetEntries(triggers[trigger]['cuts'])
	n_off_trig = chain.GetEntries(triggers[trigger]['cuts'] + "&&" + trigger)
	eff = float(n_off_trig) / n_off
	eff_err = sqrt(eff * (1 - eff) / n_off)
	print "%s efficiency: %0.1f +/- %0.1f %%" % (trigger,
												 100. * info['prescale'] * eff,
												 100. * info['prescale'] * eff_err)
		
# Local Variables:
# indent-tabs-mode: 1
# tab-width: 4
# python-indent: 4
# End:
