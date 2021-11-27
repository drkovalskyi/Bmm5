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

path_skim = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/NanoAOD-skims/516/mm/"
path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/"
split_channels = True
output = "results/summary.json"
recompute_results = False

results = dict()
if not os.path.exists('results'):
	os.mkdir('results')
if os.path.exists(output):
	results = json.load(open(output))
	
studies = {
	# ======= HLT_DoubleMu4_3_Bs ==========
	
	'HLT_DoubleMu4_3_Bs - 2016BF':{
		'samples':[
            path_skim + "/DoubleEG+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016F-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},
	'HLT_DoubleMu4_3_Bs - 2016GH':{
		'samples':[
            path_skim + "/DoubleEG+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016H-21Feb2020_UL2016-v2+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},
	'HLT_DoubleMu4_3_Bs - 2017':{
		'samples':[
            path_skim + "/DoubleEG+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017F-09Aug2019_UL2017_rsb-v2+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},
	'HLT_DoubleMu4_3_Bs - 2018':{
		'samples':[
	 		path_skim + "/EGamma+Run2018A-12Nov2019_UL2018-v2+MINIAOD/*.root",
			path_skim + "/EGamma+Run2018B-12Nov2019_UL2018-v2+MINIAOD/*.root",
			path_skim + "/EGamma+Run2018C-12Nov2019_UL2018-v2+MINIAOD/*.root",
			path_skim + "/EGamma+Run2018D-12Nov2019_UL2018-v4+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},

	'HLT_DoubleMu4_3_Bs - 2018 MC':{
		'samples':[
			path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},
	'HLT_DoubleMu4_3_Bs - 2017 MC':{
		'samples':[
			path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},
	'HLT_DoubleMu4_3_Bs - 2016BF MC':{
		'samples':[
			path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},
	'HLT_DoubleMu4_3_Bs - 2016GH MC':{
		'samples':[
			path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
		],
		'trigger':'HLT_DoubleMu4_3_Bs',
	},

	# ======= HLT_DoubleMu4_3_Jpsi ==========

	'HLT_DoubleMu4_3_Jpsi - 2018 MC':{
		'samples':[
			# path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/E9A429AF-972C-9B4D-98EF-CD3AD16A5062.root"
			path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi',
	},
	
	'HLT_DoubleMu4_3_Jpsi - 2018':{
		'samples':[
	 		path_skim + "/EGamma+Run2018A-12Nov2019_UL2018-v2+MINIAOD/*.root",
			path_skim + "/EGamma+Run2018B-12Nov2019_UL2018-v2+MINIAOD/*.root",
			path_skim + "/EGamma+Run2018C-12Nov2019_UL2018-v2+MINIAOD/*.root",
			path_skim + "/EGamma+Run2018D-12Nov2019_UL2018-v4+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi',
	},

	# ======= HLT_DoubleMu4_3_Jpsi_Displaced ==========

	'HLT_DoubleMu4_3_Jpsi_Displaced - 2017':{
		'samples':[
            path_skim + "/DoubleEG+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2017F-09Aug2019_UL2017_rsb-v2+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	},
	'HLT_DoubleMu4_3_Jpsi_Displaced - 2016BF':{
		'samples':[
            path_skim + "/DoubleEG+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016F-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016B-21Feb2020_ver1_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	},
	'HLT_DoubleMu4_3_Jpsi_Displaced - 2016GH':{
		'samples':[
            path_skim + "/DoubleEG+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/DoubleEG+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/SingleElectron+Run2016H-21Feb2020_UL2016-v2+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            path_skim + "/SinglePhoton+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	},

	'HLT_DoubleMu4_3_Jpsi_Displaced - 2017 MC':{
		'samples':[
            path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	},
	'HLT_DoubleMu4_3_Jpsi_Displaced - 2016BF MC':{
		'samples':[
            path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	},
	'HLT_DoubleMu4_3_Jpsi_Displaced - 2016GH MC':{
		'samples':[
            path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root",
		],
		'trigger':'HLT_DoubleMu4_3_Jpsi_Displaced',
	},
}

triggers = {
	'HLT_DoubleMu4_3_Bs':{
		'cuts':
		# 'sqrt(pow(abs(mm_mu2_eta-mm_mu1_eta),2) + pow(acos(cos(mm_mu2_phi-mm_mu1_phi)),2)) &&' +
		'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
		'mm_kin_pt>5 && mm_kin_vtx_prob>0.025 && Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45 &&' +
		'abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 && Muon_pt[mm_mu1_index]>4.0 &&' +
		'Muon_pt[mm_mu2_index]>4.0 && mm_kin_mass>4.6 && mm_kin_mass<5.9 && Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1'
	},
	'HLT_DoubleMu4_3_Jpsi':{
		'cuts':
		'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
		'mm_kin_alphaBS<0.4 && mm_kin_pt>7 && mm_kin_vtx_prob>0.1 && Muon_softMva[mm_mu1_index]>0.45 && ' +
		'Muon_softMva[mm_mu2_index]>0.45 && abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 &&' +
		'Muon_pt[mm_mu1_index]>4.0 && Muon_pt[mm_mu2_index]>4.0 && abs(mm_kin_mass-3.1)<0.1 &&' +
		'Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1'
	},
	'HLT_DoubleMu4_3_Jpsi_Displaced':{
		'cuts':
		'mm_mu1_index>=0 && mm_mu2_index>=0 &&' +
		'mm_kin_alphaBS<0.4 && mm_kin_pt>7 && mm_kin_vtx_prob>0.1 && Muon_softMva[mm_mu1_index]>0.45 && ' +
		'Muon_softMva[mm_mu2_index]>0.45 && abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 &&' +
		'Muon_pt[mm_mu1_index]>4.0 && Muon_pt[mm_mu2_index]>4.0 && abs(mm_kin_mass-3.1)<0.1 &&' +
		'Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1 &&' +
		'mm_kin_slxy > 4.0'
	}
}

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv8-516/trigger_efficiency"

ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 800,600)


def load_data():
	chain = ROOT.TChain("Events")
	for sample in info['samples']:
		chain.Add(sample)

	print "Total number of events:", chain.GetEntries()
	return chain

for name, info in sorted(studies.items()):
	print "Processing", name
	
	chain = None
	trigger = info['trigger']

	for ch in range(2):
		file_name = re.sub('\s+', '', name)
		study_name = name
		extra_cut = ""
		if split_channels:
			print "Channel:", ch
			file_name += "_%u" % ch
			study_name += " (channel: %u)" % ch
			if ch == 0:
				extra_cut = "&& abs(Muon_eta[mm_mu1_index])<0.7 && abs(Muon_eta[mm_mu2_index])<0.7"
			else:
				extra_cut = "&& (abs(Muon_eta[mm_mu1_index])>0.7 || abs(Muon_eta[mm_mu2_index])>0.7)"

		if study_name in results and not recompute_results:
			print "Results are already available. Skip the study"
			continue

		if chain == None:
			chain = load_data()
		f = ROOT.TFile.Open('results/' + file_name + ".root", "recreate")
		nbins = 60
		h_off = ROOT.TH1F("h_off","h_off", nbins, 0, nbins)
		h_off.Sumw2()
		chain.Draw("PV_npvsGood>>h_off", "%s" % (triggers[trigger]['cuts'] + extra_cut))
		h_off.Write()
		
		h_off_trig = ROOT.TH1F("h_off_trig","h_off_trig", nbins, 0, nbins)
		h_off_trig.Sumw2()
		chain.Draw("PV_npvsGood>>h_off_trig",
				   "(%s)*prescale_%s" % (triggers[trigger]['cuts'] + extra_cut + "&&" + trigger, trigger))
		h_off_trig.Write()

		# n_off_trig = chain.GetEntries(triggers[trigger]['cuts'] + "&&" + trigger)
		err_off = ROOT.Double(0)
		n_off = h_off.IntegralAndError(1, nbins, err_off)
	
		err_off_trig = ROOT.Double(0)
		n_off_trig = h_off_trig.IntegralAndError(1, nbins, err_off_trig)
	
		eff = float(n_off_trig) / n_off
		eff_err = eff * sqrt((err_off/n_off)**2 + (err_off_trig/n_off_trig)**2)
		print "%s efficiency: %0.1f +/- %0.1f %%" % (trigger,
													 100. * eff,
													 100. * eff_err)
		results[study_name] = {
			'eff':eff,
			'eff_err':eff_err
		}
		f.Close()
		if not split_channels:
			break

## Save results
json.dump(results, open(output, 'w'))

## Produce report
print "\nSummary"
for study, info in sorted(results.items()):
	print "%s \t: %0.1f +/- %0.1f %%" % (study, 100. * info['eff'], 100. * info['eff_err'])

## Ratios
def compute_ratio(bmm, jpsik):
	bmm_rel_err = results[bmm]['eff_err']/results[bmm]['eff']
	jpsik_rel_err = results[jpsik]['eff_err']/results[jpsik]['eff']
	ratio = results[jpsik]['eff']/results[bmm]['eff']
	ratio_err = sqrt(bmm_rel_err * bmm_rel_err + jpsik_rel_err * jpsik_rel_err) * ratio
	print "%s / %s: \t%0.2f \pm %0.2f" % (jpsik, bmm, ratio, ratio_err)

print "\nRatios"
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
