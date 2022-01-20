"""PileUp efficiency study

The script extracts MC efficiency as a function of nPV. It can also
compute weighted efficiency for any nPV distribution.

The script is fairly slow (~3h), so you may want to use a subset of
data for quicker checks.

"""

import sys, os, subprocess, re, json
import ROOT
from math import *
from tdrstyle import *

path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/"
nbins = 60

hlt_studies = json.load(open("results/summary.json"))


samples = {
	'Bsmm 2018':[
		path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/*.root",
	],
	'Bsmm 2017':[
		path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v1+MINIAODSIM/*.root",
	],
	'Bsmm 2016BF':[
		path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
	],
	'Bsmm 2016GH':[
		path + "/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
	],
	'Bjpsik 2018':[
		# path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/E9A429AF-972C-9B4D-98EF-CD3AD16A5062.root"
		path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
	],
	'Bjpsik 2017':[
        path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root",
	],
	'Bjpsik 2016BF':[
        path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root",
	],
	'Bjpsik 2016GH':[
        path + "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root",
	],
}

bsmm_preselection = 'mm_gen_pdgId != 0 && ' + \
	'mm_mu1_index>=0 && mm_mu2_index>=0 && ' + \
	'Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45 && ' + \
	'Muon_charge[mm_mu1_index]*Muon_charge[mm_mu2_index]==-1 && ' + \
	'abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4 && ' + \
	'Muon_pt[mm_mu1_index]>4.0 && Muon_pt[mm_mu2_index]>4.0 && ' + \
	'mm_kin_pt>5 && mm_kin_vtx_prob>0.025 &&' + \
    'abs(mm_kin_mass-5.4)<0.5 && mm_kin_sl3d>6'

bjpsik_preselection = "bkmm_gen_pdgId != 0 && " + \
	"mm_mu1_index[bkmm_mm_index]>=0 && "\
    "mm_mu2_index[bkmm_mm_index]>=0 && "\
	'Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45 && ' + \
	'Muon_charge[mm_mu1_index[bkmm_mm_index]]*Muon_charge[mm_mu2_index[bkmm_mm_index]]==-1 && ' + \
    "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && "\
    "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && "\
    "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && "\
    "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 && "\
    "mm_kin_pt[bkmm_mm_index]>7.0 && "\
    "mm_kin_alphaBS[bkmm_mm_index]<0.4 && "\
    "mm_kin_vtx_prob[bkmm_mm_index]>0.1 && "\
    "bkmm_jpsimc_vtx_prob>0.025 && "\
    "mm_kin_sl3d[bkmm_mm_index]>4 && "\
    "abs(bkmm_jpsimc_mass-5.4)<0.5"

selections = {
	'Bsmm':bsmm_preselection,
	'Bsmm MVA>0.90':bsmm_preselection + "&& mm_mva>0.90",
	'Bsmm MVA>0.99':bsmm_preselection + "&& mm_mva>0.99",
	'Bjpsik':bjpsik_preselection,
	'Bjpsik MVA>0.90':bjpsik_preselection + "&& bkmm_bmm_mva>0.90",
	'Bjpsik MVA>0.99':bjpsik_preselection + "&& bkmm_bmm_mva>0.99",
}

jpsi_triggers = {
	'2016BF':'HLT_DoubleMu4_3_Jpsi_Displaced',
	'2016GH':'HLT_DoubleMu4_3_Jpsi_Displaced',
	'2017':'HLT_DoubleMu4_3_Jpsi_Displaced',
	'2018':'HLT_DoubleMu4_3_Jpsi',
}

studies = dict()

for era in ['2016BF', '2016GH', '2017', '2018']:
# for era in ['2018']:
	studies['%s channel 0 Bsmm preselection' % era] = {
		'sample':'Bsmm %s' % era,
		'selection':'Bsmm',
		'hist':'results/HLT_DoubleMu4_3_Bs-%s_0.root' % era
	}

	studies['%s channel 1 Bsmm preselection' % era] = {
		'sample':'Bsmm %s' % era,
		'selection':'Bsmm',
		'hist':'results/HLT_DoubleMu4_3_Bs-%s_1.root' % era
	}
	
	studies['%s channel 0 Bjpsik preselection' % era] = {
		'sample':'Bjpsik %s' % era,
		'selection':'Bjpsik',
		'hist':'results/%s-%s_0.root' % (jpsi_triggers[era], era)
	}
	
	studies['%s channel 1 Bjpsik preselection' % era] = {
		'sample':'Bjpsik %s' % era,
		'selection':'Bjpsik',
		'hist':'results/%s-%s_1.root' % (jpsi_triggers[era], era)
	}

	#  MVA>0.90
	studies['%s channel 0 Bsmm MVA>0.90' % era] = {
		'sample':'Bsmm %s' % era,
		'selection':'Bsmm MVA>0.90',
		'hist':'results/HLT_DoubleMu4_3_Bs-%s_0.root' % era
	}

	studies['%s channel 1 Bsmm MVA>0.90' % era] = {
		'sample':'Bsmm %s' % era,
		'selection':'Bsmm MVA>0.90',
		'hist':'results/HLT_DoubleMu4_3_Bs-%s_1.root' % era
	}
	
	studies['%s channel 0 Bjpsik MVA>0.90' % era] = {
		'sample':'Bjpsik %s' % era,
		'selection':'Bjpsik MVA>0.90',
		'hist':'results/%s-%s_0.root' % (jpsi_triggers[era], era)
	}
	
	studies['%s channel 1 Bjpsik MVA>0.90' % era] = {
		'sample':'Bjpsik %s' % era,
		'selection':'Bjpsik MVA>0.90',
		'hist':'results/%s-%s_1.root' % (jpsi_triggers[era], era)
	}

	#  MVA>0.99
	studies['%s channel 0 Bsmm MVA>0.99' % era] = {
		'sample':'Bsmm %s' % era,
		'selection':'Bsmm MVA>0.99',
		'hist':'results/HLT_DoubleMu4_3_Bs-%s_0.root' % era
	}

	studies['%s channel 1 Bsmm MVA>0.99' % era] = {
		'sample':'Bsmm %s' % era,
		'selection':'Bsmm MVA>0.99',
		'hist':'results/HLT_DoubleMu4_3_Bs-%s_1.root' % era
	}
	
	studies['%s channel 0 Bjpsik MVA>0.99' % era] = {
		'sample':'Bjpsik %s' % era,
		'selection':'Bjpsik MVA>0.99',
		'hist':'results/%s-%s_0.root' % (jpsi_triggers[era], era)
	}
	
	studies['%s channel 1 Bjpsik MVA>0.99' % era] = {
		'sample':'Bjpsik %s' % era,
		'selection':'Bjpsik MVA>0.99',
		'hist':'results/%s-%s_1.root' % (jpsi_triggers[era], era)
	}
	

ROOT.gROOT.SetBatch(True)

data = dict()

def load_data(name):
	if name in data:
		return data[name]
	
	chain = ROOT.TChain("Events")
	for sample in samples[name]:
		chain.Add(sample)

	data[name] = chain

	return data[name]

for study, info in sorted(studies.items()):
	print study
	
	chain = load_data(info['sample'])
	
	f = ROOT.TFile(info['hist'])

	print "Total number of events:", chain.GetEntries()

	h_all = ROOT.TH1F("h_all","h_all", nbins, 0, nbins)
	h_all.Sumw2()
	chain.Draw("PV_npvsGood>>h_all", "", "goff")

	h_pre = ROOT.TH1F("h_pre","h_pre", nbins, 0, nbins)
	h_pre.Sumw2()
	chain.Draw("PV_npvsGood>>h_pre", selections[info['selection']], "goff")

	h_eff = h_pre.Clone("h_eff")
	h_tmp = h_pre.Clone("h_tmp")
	h_eff.Divide(h_pre, h_all, 1, 1, "B")

	h_tmp.Multiply(h_eff, h_all, 1, 1/h_all.Integral())
	eff_mc_0 = h_tmp.Integral()

	h_tmp.Multiply(h_eff, h_pre, 1, 1/h_pre.Integral())
	eff_mc_1 = h_tmp.Integral()

	# h_off = f.Get("h_off")
	# h_tmp.Multiply(h_eff, h_off, 1, 1/h_off.Integral())
	# eff_data_all = h_tmp.Integral()

	# print "Relative efficiency deviation (data/mc-1): %0.1f %%" % (100.*(eff_data_all/eff_mc_0-1))

	h_off_trig = f.Get("h_off_trig")
	h_tmp.Multiply(h_eff, h_off_trig, 1, 1/h_off_trig.Integral())
	eff_data_trigger = h_tmp.Integral()

	print "Relative efficiency deviation for trigger (data/mc-1): %0.1f %%" % (100.*(eff_data_trigger/eff_mc_0-1))

	print


# trigger = info['trigger']

# 	for ch in range(2):
# 		file_name = re.sub('\s+', '', name)
# 		study_name = name
# 		extra_cut = ""
# 		if split_channels:
# 			print "Channel:", ch
# 			file_name += "_%u" % ch
# 			study_name += " (channel: %u)" % ch
# 			if ch == 0:
# 				extra_cut = "&& abs(Muon_eta[mm_mu1_index])<0.7 && abs(Muon_eta[mm_mu2_index])<0.7"
# 			else:
# 				extra_cut = "&& (abs(Muon_eta[mm_mu1_index])>0.7 || abs(Muon_eta[mm_mu2_index])>0.7)"

# 		if study_name in results and not recompute_results:
# 			print "Results are already available. Skip the study"
# 			continue

# 		if chain == None:
# 			chain = load_data()
# 		f = ROOT.TFile.Open('results/' + file_name + ".root", "recreate")
# 		nbins = 60
# 		h_off = ROOT.TH1F("h_off","h_off", nbins, 0, nbins)
# 		h_off.Sumw2()
# 		chain.Draw("PV_npvsGood>>h_off", "%s" % (triggers[trigger]['cuts'] + extra_cut))
# 		h_off.Write()
		
# 		h_off_trig = ROOT.TH1F("h_off_trig","h_off_trig", nbins, 0, nbins)
# 		h_off_trig.Sumw2()
# 		chain.Draw("PV_npvsGood>>h_off_trig",
# 				   "(%s)*prescale_%s" % (triggers[trigger]['cuts'] + extra_cut + "&&" + trigger, trigger))
# 		h_off_trig.Write()

# 		# n_off_trig = chain.GetEntries(triggers[trigger]['cuts'] + "&&" + trigger)
# 		err_off = ROOT.Double(0)
# 		n_off = h_off.IntegralAndError(1, nbins, err_off)
	
# 		err_off_trig = ROOT.Double(0)
# 		n_off_trig = h_off_trig.IntegralAndError(1, nbins, err_off_trig)
	
# 		eff = float(n_off_trig) / n_off
# 		eff_err = eff * sqrt((err_off/n_off)**2 + (err_off_trig/n_off_trig)**2)
# 		print "%s efficiency: %0.1f +/- %0.1f %%" % (trigger,
# 													 100. * eff,
# 													 100. * eff_err)
# 		results[study_name] = {
# 			'eff':eff,
# 			'eff_err':eff_err
# 		}
# 		f.Close()
# 		if not split_channels:
# 			break

# ## Save results
# json.dump(results, open(output, 'w'))

# ## Produce report
# print "\nSummary"
# for study, info in sorted(results.items()):
# 	print "%s \t: %0.1f +/- %0.1f %%" % (study, 100. * info['eff'], 100. * info['eff_err'])


# ## Make plots
# sys.exit()

# setTDRStyle()
# c1 = ROOT.TCanvas("c1", "c1", 800, 800)
# def make_plot(name, filename_mc, filename_data, plots, rebin=1):
# 	f_mc = ROOT.TFile.Open(filename_mc)
# 	f_data = ROOT.TFile.Open(filename_data)
# 	normalization = None
# 	legend = ROOT.TLegend(0.70,0.75,0.85,0.87)
# 	legend.SetShadowColor(ROOT.kWhite)
# 	legend.SetLineColor(ROOT.kWhite)
# 	# legend.SetFillColor(10)
# 	colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kOrange+5, ROOT.kGreen+3]
# 	scale = 1.2
# 	for idx, plot in enumerate(plots):
# 		if plot['type'] == "mc":
# 			f = f_mc
# 		else:
# 			f = f_data
# 		h = f.Get(plot['name'])
# 		if rebin > 1:
# 			h.Rebin(rebin)
# 		h.SetLineColor(colors[idx])
# 		h.SetLineWidth(2)
# 		h.SetMarkerStyle(20)
# 		h.SetMarkerColor(colors[idx])
# 		h.GetXaxis().SetTitle("nPV")
# 		if normalization:
# 			h.Scale(normalization/h.Integral())
# 			h.Draw("same")
# 		else:
# 			normalization = h.Integral()
# 			h.SetMaximum(h.GetMaximum()*scale)
# 			h.Draw()
# 		legend.AddEntry(h, plot['label'])
# 	legend.Draw()
# 	c1.Print(name)
	

# make_plot("HLT_DoubleMu4_3_Bs-2018.pdf", "HLT_DoubleMu4_3_Bs-2018MC.root", "HLT_DoubleMu4_3_Bs-2018.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Bs-2017.pdf", "HLT_DoubleMu4_3_Bs-2017MC.root", "HLT_DoubleMu4_3_Bs-2017.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Bs-2016BF.pdf", "HLT_DoubleMu4_3_Bs-2016BFMC.root", "HLT_DoubleMu4_3_Bs-2016BF.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Bs-2016GH.pdf", "HLT_DoubleMu4_3_Bs-2016GHMC.root", "HLT_DoubleMu4_3_Bs-2016GH.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi-2018.pdf", "HLT_DoubleMu4_3_Jpsi-2018MC.root", "HLT_DoubleMu4_3_Jpsi-2018.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2017.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2017MC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2017.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016BF.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BFMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BF.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016GH.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GHMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GH.root", 
# 		  plots = [
# 			  {'type':'mc', 'name':'h_off_trig', 'label':'MC'}, 
# 			  {'type':'data','name':'h_off_trig','label':'Data'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2017-Data.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2017MC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2017.root", 
# 		  plots = [
# 			  {'type':'Data', 'name':'h_off', 'label':'All'}, 
# 			  {'type':'Data','name':'h_off_trig','label':'Triggered'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016BF-Data.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BFMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016BF.root", 
# 		  plots = [
# 			  {'type':'Data', 'name':'h_off', 'label':'All'}, 
# 			  {'type':'Data','name':'h_off_trig','label':'Triggered'}
# 		  ],
# 		  rebin = 3
# )
# make_plot("HLT_DoubleMu4_3_Jpsi_Displaced-2016GH-Data.pdf", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GHMC.root", "HLT_DoubleMu4_3_Jpsi_Displaced-2016GH.root", 
# 		  plots = [
# 			  {'type':'Data', 'name':'h_off', 'label':'All'}, 
# 			  {'type':'Data','name':'h_off_trig','label':'Triggered'}
# 		  ],
# 		  rebin = 3
# )
		  
		
# Local Variables:
# indent-tabs-mode: 1
# tab-width: 4
# python-indent: 4
# End:
