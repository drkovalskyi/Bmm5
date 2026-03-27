import ROOT
import sys, os, subprocess
from DataFormats.FWLite import Events, Handle
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from math import *

# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/muon_fake_sources_bhh_medium"
# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/muon_fake_sources_bhh_loose"
# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/muon_fake_sources_bhh_mva"
# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/muon_fake_sources_mu_enriched _loose"
# output_path = "/eos/home-d/dmytro/www/plots/tmp/muon_fake_sources_ttbar_semi_Run3Summer22EE"
# output_path = "/eos/home-d/dmytro/www/plots/tmp/muon_fake_sources_ttbar_semi_RunIISummer20UL18"
# output_path = "/eos/home-d/dmytro/www/plots/tmp/muon_fake_sources_minbias_ee_Run3Summer22"
# output_path = "/eos/home-d/dmytro/www/plots/muon_fake_sources/run2-bhh_softmvaid"
output_path = "/eos/home-d/dmytro/www/plots/muon_fake_sources/run2-bhh_loose"
dump_info = False
# dump_info = True
min_pt = 4
nmax_files = 1000

muon_id = None
# muon_id = "MediumId"
# muon_id = "SoftMvaId"

use_gen_filter = False

files = []

# find files
path = "/eos/cms/store/user/dmytro/"
# path = "/eos/cms/store/mc/"
pds = [
	# 'JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8'

	# MuEnriched - lots of data
	# 'QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_MuonFakeSkim_QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8_1603447990',
	# 'QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_MuonFakeSkim_QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8_1603446796',
	# 'QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_MuonFakeSkim_QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8_1603448253',
	# 'QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_MuonFakeSkim_QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8_1603448316',
	# 'QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_MuonFakeSkim_QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_1603448633'
	
	# Bhh
	'BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603710814',
	'BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603711796',
	'BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603711442',
	'BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603711503',
	'BsToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_BsToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603731003',
	'LambdaBToPK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_LambdaBToPK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603711266',
	'LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_MuonFakeSkim_LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_1603710918'

	# "Run3Summer22EEMiniAODv4/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/Poisson60ForMUOVal_130X_mcRun3_2022_realistic_postEE_v6-v2/"
	# "RunIISummer20UL18MiniAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PUForMUOVal_106X_upgrade2018_realistic_v16_L1v1_ext1-v2/"
	# "Run3Summer22MiniAODv4/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/validDigi_130X_mcRun3_2022_realistic_v5-v4/"
]

for pd in pds:
    # output = subprocess.check_output("find %s/%s -type f -name '*.root'|grep crab_MuonFakeSkim" % (path, pd), shell=True)
    output = subprocess.check_output("find %s/%s -type f -name '*.root'" % (path, pd), shell=True)
    for f in output.decode("utf-8").split("\n"):
        if f.strip():
            files.append(f)

print("Number of files found: %u" % len(files))
print("Using %u files" % min(nmax_files, len(files)))
events = Events(files[:nmax_files])

# events = Events (
# 	[
# 		# '/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/muon_fake_skim.root'
# 		# '/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_RunIIAutumn18MiniAOD_muon_fake_skim.root',
# 		#'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM.root'
# 		# '/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/test.root'
# 		'/eos/cms/store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv4+TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8+MINIAODSIM+Poisson60ForMUOVal_130X_mcRun3_2022_realistic_postEE_v6-v2+2520000+00ac69e4-f54a-4193-9e4b-4c9236b6f3da.root'
# 	]
# )

def isAncestor(a,p) :
	if a == p : 
		return True
	for i in range(0,p.numberOfMothers()) :
		if isAncestor(a,p.mother(i)) :
			return True
	return False

def find_parent(cand):
    # look for a parent with a different PDG id
    parent = cand
    # do a fixed depth loop to avoid endless loops
    for i in range(100):
        parent = parent.mother()
        if not parent: break
        if parent.pdgId()!=cand.pdgId():
            return parent
    return None

def isGoodMuon(muon):
	if not muon.isTrackerMuon(): return False
	if not muon.isGlobalMuon(): return False
	if not muon.isLooseMuon(): return False
	if not muon.innerTrack().quality(ROOT.reco.Track.highPurity): return False 
	if muon.pt() < min_pt or abs(muon.eta()) > 1.4: return False
	if muon_id != None:
		if muon_id == "SoftMvaId":
			# if not muon.passed(ROOT.reco.Muon.SoftMvaId): return False
			if muon.softMvaValue() < 0.58: return False
		elif muon_id == "MediumId":
			if not muon.isMediumMuon(): return False
		else:
			raise Exception("Uknown muon_id: %s" % muon_id)
	return True

def deltaPhi(phi1,phi2):
    return acos(cos(phi2-phi1))

def deltaR(p1,p2):
    return sqrt(pow(deltaPhi(p1.phi(),p2.phi()),2)+pow(p2.eta()-p1.eta(),2))

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))
    canvas.Print("%s/%s.C"%(path,output_name_without_extention))

def make_sim_type_hist(name):
	h = ROOT.TH1D(name, "Muon matching based on simulated hits", 9, 0, 9)
	h.GetXaxis().SetBinLabel(1, "No Match")
	h.GetXaxis().SetBinLabel(2, "Punch through")
	h.GetXaxis().SetBinLabel(3, "#mu from pion")
	h.GetXaxis().SetBinLabel(4, "#mu from kaon")
	h.GetXaxis().SetBinLabel(5, "#mu from c-hadron")
	h.GetXaxis().SetBinLabel(6, "#mu from b-hadron")
	h.GetXaxis().SetBinLabel(7, "#mu from W/Z")
	h.GetXaxis().SetBinLabel(8, "#mu from #tau")
	h.GetXaxis().SetBinLabel(9, "Other #mu")
	h.SetFillColor(ROOT.kMagenta)
	
	return h
	
fake_types = {
	211:'pion',
	321:'kaon',
	2212:'proton'
}

def get_sim_type(muon):
	if muon.simPdgId() == 0: return 'not_matched'
	if abs(muon.simPdgId()) == 13:
		if abs(muon.simMotherPdgId()) == 211: return 'pion'
		if abs(muon.simMotherPdgId()) == 321: return 'kaon'
		if abs(muon.simMotherPdgId()) == 2212: return 'proton'
		if abs(muon.simFlavour()) == 4: return 'c-hadron'
		if abs(muon.simFlavour()) == 5: return 'b-hadron'
		if abs(muon.simFlavour()) == 13: return 'prompt muon'
		if abs(muon.simFlavour()) == 15: return 'tau'
		return 'other'
	return 'punch_through'

sim_type_to_value = {
	'not_matched':0.5,
	'punch_through':1.5,
	'pion':2.5,
	'kaon':3.5,
	'c-hadron':4.5,
	'b-hadron':5.5,
	'prompt muon':6.5,
	'tau':7.5,
	'other':8.5
}


ROOT.gROOT.SetBatch(True)

handlePruned  = Handle ("std::vector<reco::GenParticle>")
handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
labelPruned = ("prunedGenParticles")
labelPacked = ("packedGenParticles")

muonHandle, muonLabel = Handle("std::vector<pat::Muon>"),"slimmedMuons"

h_sim_match_all = make_sim_type_hist("h_sim_match_all")

h_sim_match = dict()
h_genid_not_matched = dict()
h_sim_relative_pt = dict()
h_sim_relative_pt_type_matched = dict()
h_sim_relative_pt_wrong_type_matched = dict()
h_sim_relative_pt_matched_to_other = dict()
h_sim_relative_pt_matched_to_punchthrough = dict()
h_sim_decay_rho = dict()
h_sim_decay_rho_type_matched = dict()
h_sim_decay_rho_wrong_type_matched = dict()
h_sim_decay_rho_matched_to_other = dict()
h_sim_decay_rho_matched_to_punchthrough = dict()
h_sim_decay_rho_matched_to_other_same_pt = dict()
for id,name in list(fake_types.items()):
	h_sim_match[name] = make_sim_type_hist("h_sim_match_%s" % name)

	h_genid_not_matched[name] = ROOT.TH1D("h_genid_not_matched_%s" % name,
									"Gen |pdgId| for muons not matched by simulated hits", 350, 0, 350)
	h_genid_not_matched[name].SetLineColor(ROOT.kBlue)
	h_genid_not_matched[name].SetLineWidth(2)
	h_genid_not_matched[name].GetXaxis().SetTitle("|pdgId|")
	
	h_sim_relative_pt[name] = ROOT.TH1D("h_sim_relative_pt_%s" % name,
										"Relative Pt of sim-matched particle", 100, 0, 2)
	h_sim_relative_pt[name].SetLineColor(ROOT.kBlue)
	h_sim_relative_pt[name].SetLineWidth(2)
	h_sim_relative_pt[name].GetXaxis().SetTitle("Pt_{sim}/Pt_{reco}")
	h_sim_relative_pt_type_matched[name] = ROOT.TH1D("h_sim_relative_pt_type_matched_%s" % name,
										"Relative Pt of sim particle (type matched)", 100, 0, 2)
	h_sim_relative_pt_type_matched[name].SetLineColor(ROOT.kBlue)
	h_sim_relative_pt_type_matched[name].SetLineWidth(2)
	h_sim_relative_pt_type_matched[name].GetXaxis().SetTitle("Pt_{sim}/Pt_{reco}")
	h_sim_relative_pt_wrong_type_matched[name] = ROOT.TH1D("h_sim_relative_pt_wrong_type_matched_%s" % name,
														   "Relative Pt of sim particle (wrong match)", 100, 0, 2)
	h_sim_relative_pt_wrong_type_matched[name].SetLineColor(ROOT.kBlue)
	h_sim_relative_pt_wrong_type_matched[name].SetLineWidth(2)
	h_sim_relative_pt_wrong_type_matched[name].GetXaxis().SetTitle("Pt_{sim}/Pt_{reco}")
	h_sim_relative_pt_matched_to_other[name] = ROOT.TH1D("h_sim_relative_pt_matched_to_other_%s" % name,
														   "Relative Pt of sim particle (other match)", 100, 0, 2)
	h_sim_relative_pt_matched_to_other[name].SetLineColor(ROOT.kBlue)
	h_sim_relative_pt_matched_to_other[name].SetLineWidth(2)
	h_sim_relative_pt_matched_to_other[name].GetXaxis().SetTitle("Pt_{sim}/Pt_{reco}")
	h_sim_relative_pt_matched_to_punchthrough[name] = ROOT.TH1D("h_sim_relative_pt_matched_to_punchthrough_%s" % name,
														   "Relative Pt of sim particle (punchthrough match)", 100, 0, 2)
	h_sim_relative_pt_matched_to_punchthrough[name].SetLineColor(ROOT.kBlue)
	h_sim_relative_pt_matched_to_punchthrough[name].SetLineWidth(2)
	h_sim_relative_pt_matched_to_punchthrough[name].GetXaxis().SetTitle("Pt_{sim}/Pt_{reco}")

	h_sim_decay_rho[name] = ROOT.TH1D("h_sim_decay_rho_%s" % name,
									  "Decay radius", 100, 0, 400)
	h_sim_decay_rho[name].SetLineColor(ROOT.kBlue)
	h_sim_decay_rho[name].SetLineWidth(2)
	h_sim_decay_rho[name].GetXaxis().SetTitle("#rho")
	h_sim_decay_rho_type_matched[name] = ROOT.TH1D("h_sim_decay_rho_type_matched_%s" % name,
												   "Decay radius for type matched", 100, 0, 400)
	h_sim_decay_rho_type_matched[name].SetLineColor(ROOT.kBlue)
	h_sim_decay_rho_type_matched[name].SetLineWidth(2)
	h_sim_decay_rho_type_matched[name].GetXaxis().SetTitle("#rho")
	h_sim_decay_rho_wrong_type_matched[name] = ROOT.TH1D("h_sim_decay_rho_wrong_type_matched_%s" % name,
												   "Decay radius for matched to wrong type", 100, 0, 400)
	h_sim_decay_rho_wrong_type_matched[name].SetLineColor(ROOT.kBlue)
	h_sim_decay_rho_wrong_type_matched[name].SetLineWidth(2)
	h_sim_decay_rho_wrong_type_matched[name].GetXaxis().SetTitle("#rho")
	h_sim_decay_rho_matched_to_other[name] = ROOT.TH1D("h_sim_decay_rho_matched_to_other_%s" % name,
												   "Decay radius for matched to other type", 100, 0, 400)
	h_sim_decay_rho_matched_to_other[name].SetLineColor(ROOT.kBlue)
	h_sim_decay_rho_matched_to_other[name].SetLineWidth(2)
	h_sim_decay_rho_matched_to_other[name].GetXaxis().SetTitle("#rho")
	h_sim_decay_rho_matched_to_other_same_pt[name] = ROOT.TH1D("h_sim_decay_rho_matched_to_other_same_pt_%s" % name,
												   "Decay radius for matched to 'other' type with |Pt_sim/Pt_reco-1|<0.05", 100, 0, 400)
	h_sim_decay_rho_matched_to_other_same_pt[name].SetLineColor(ROOT.kBlue)
	h_sim_decay_rho_matched_to_other_same_pt[name].SetLineWidth(2)
	h_sim_decay_rho_matched_to_other_same_pt[name].GetXaxis().SetTitle("#rho")
	h_sim_decay_rho_matched_to_punchthrough[name] = ROOT.TH1D("h_sim_decay_rho_matched_to_punchthrough_%s" % name,
												   "Decay radius for matched to punchthrough type", 100, 0, 400)
	h_sim_decay_rho_matched_to_punchthrough[name].SetLineColor(ROOT.kBlue)
	h_sim_decay_rho_matched_to_punchthrough[name].SetLineWidth(2)
	h_sim_decay_rho_matched_to_punchthrough[name].GetXaxis().SetTitle("#rho")
	
# loop over events
count= 0
for event in events:
    # muons
	event.getByLabel(muonLabel, muonHandle)
	muons = muonHandle.product()
	has_interesting_muon = False
	for muon in muons:
		if not isGoodMuon(muon): continue
		has_interesting_muon = True

	if not has_interesting_muon: continue

	# all muon candidates
	h_sim_match_all.Fill(sim_type_to_value[get_sim_type(muon)])
	
	# check if we have a match
	event.getByLabel (labelPacked, handlePacked)
	event.getByLabel (labelPruned, handlePruned)
	# get the product
	packed = handlePacked.product()
	pruned = handlePruned.product()

	if use_gen_filter:
		# select di-electron events
		n_electrons = 0
		for p in pruned:
			if p.status() != 1: continue
			if abs(p.pdgId()) != 11: continue
			if abs(p.eta()) > 2.5: continue
			if p.pt() < 4.0: continue
			n_electrons += 1

		if n_electrons < 2:
			continue
	
	interesting_event = False

	for p in pruned:
		if abs(p.pdgId()) not in (321, 211, 2212): continue
		fake_type = fake_types[abs(p.pdgId())]
		# mother = find_parent(p)
		# if not mother or abs(mother.pdgId()) != 511: continue
		muon = None
		for m in muons:
			if deltaR(m,p)>0.1: continue
			if abs(m.pt()/p.pt()-1)>0.02: continue
			if not isGoodMuon(m): continue
			muon = m
			break
		if not muon: continue
		h_sim_match[fake_type].Fill(sim_type_to_value[get_sim_type(muon)])
		if muon.simPdgId()!=0:
			h_sim_relative_pt[fake_type].Fill(min(muon.simPt()/p.pt(), 1.999))
			h_sim_decay_rho[fake_type].Fill(min(muon.simProdRho(), 399.999))
			if fake_type == get_sim_type(muon):
				h_sim_relative_pt_type_matched[fake_type].Fill(min(muon.simPt()/p.pt(), 1.999))
				h_sim_decay_rho_type_matched[fake_type].Fill(min(muon.simProdRho(), 399.999))
			else:
				h_sim_relative_pt_wrong_type_matched[fake_type].Fill(min(muon.simPt()/p.pt(), 1.999))
				h_sim_decay_rho_wrong_type_matched[fake_type].Fill(min(muon.simProdRho(), 399.999))
				if get_sim_type(muon) == 'other':
					h_sim_relative_pt_matched_to_other[fake_type].Fill(min(muon.simPt()/p.pt(), 1.999))
					h_sim_decay_rho_matched_to_other[fake_type].Fill(min(muon.simProdRho(), 399.999))
					if abs(muon.simPt()/p.pt() - 1) < 0.05:
						h_sim_decay_rho_matched_to_other_same_pt[fake_type].Fill(min(muon.simProdRho(), 399.999))
			if get_sim_type(muon) == "punch_through":
				h_sim_relative_pt_matched_to_punchthrough[fake_type].Fill(min(muon.simPt()/p.pt(), 1.999))
				h_sim_decay_rho_matched_to_punchthrough[fake_type].Fill(min(muon.simProdRho(), 399.999))
				
		else:
			if muon.genParticle():
				h_genid_not_matched[fake_type].Fill(min(abs(muon.genParticle().pdgId()),399.999))
			else:
				h_genid_not_matched[fake_type].Fill(0)
			
		interesting_event = True
		
		status = p.statusFlags().fromHardProcessBeforeFSR()
		if dump_info:
			print("Reco Muon pt : %s  eta : %s   phi: %s" % (muon.pt(), muon.eta(), muon.phi()))
			print("\tsim_pt: %s sim_flavour: %s sim_pid: %s sim_mother_pid: %s sim_prod_rho: %f sim_prod_z: %f" % (muon.simPt(), muon.simFlavour(), muon.simPdgId(), muon.simMotherPdgId(),
																									muon.simProdRho(), muon.simProdZ()))
			print("Gen pdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(p.pdgId(), p.pt(), p.eta(), p.phi(), status))    
			mother = find_parent(p)
			if mother:
				status = mother.statusFlags().fromHardProcessBeforeFSR()
				print("\tMother PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(mother.pdgId(),mother.pt(),mother.eta(),mother.phi(),status))    
				mother = find_parent(mother)
				if mother:
					status = mother.statusFlags().fromHardProcessBeforeFSR()
					print("\tGrand Mother PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(mother.pdgId(),mother.pt(),mother.eta(),mother.phi(),status))    

	if not interesting_event: continue

	count += 1

	# https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
	# if dump_info:
	# 	print("Event dump")
	# 	print("Run: %u, Event: %u" % (event.eventAuxiliary().run(),event.eventAuxiliary().event()))
	# 	for p in pruned :
	# 		if not abs(p.pdgId()) in [511,521,531]: continue

	# 		final_b = True
	# 		for dau in p.daughterRefVector():
	# 			if dau.pdgId() == -p.pdgId():
	# 				final_b = False
	# 				break
	# 		if not final_b: continue

	# 		signature = 1
	# 		for pa in packed:
	# 			mother = pa.mother(0)
	# 			if mother and isAncestor(p,mother) :
	# 				if pa.pdgId()!=22: signature *= pa.pdgId()

	# 		# if not signature in [13*13*321*321, 13*13, 13*13*321, -13*13*321]: continue

	# 		# d_p4 = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D(ROOT.Double))()
	# 		d_p4 = ROOT.reco.Candidate.LorentzVector()
	# 		rad_p4 = ROOT.reco.Candidate.LorentzVector()
	# 		status = p.statusFlags().fromHardProcessBeforeFSR()
	# 		print("     PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(p.pdgId(),p.pt(),p.eta(),p.phi(),status))    
	# 		for dau in p.daughterRefVector():
	# 			status = dau.statusFlags().fromHardProcessBeforeFSR()
	# 			print("     dau PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(dau.pdgId(),dau.pt(),dau.eta(),dau.phi(),status))    
	# 		if p.mother():
	# 			status = p.mother().statusFlags().fromHardProcessBeforeFSR()
	# 			print("     mother PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(p.mother().pdgId(),p.mother().pt(),p.mother().eta(),p.mother().phi(),status))    
	# 		for pa in packed:
	# 			mother = pa.mother(0)
	# 			if mother and isAncestor(p,mother) :
	# 				print("          PdgId : %s   pt : %s  eta : %s   phi : %s" %(pa.pdgId(),pa.pt(),pa.eta(),pa.phi()))
	# 				d_p4 += pa.p4()
	# 				if pa.pdgId()==22: rad_p4 += pa.p4()
	# 		print("     delta: %0.5f%%" % (100.*(p.p4()-d_p4).P()/p.p4().P()))
	# 		print("     radiation: %0.2f%%" % (100.*rad_p4.P()/p.p4().P()))
	if dump_info and count >= 100: 
		sys.exit()

c1 = TCanvas("c1", "c1", 800, 800)
h_sim_match_all.Draw()
print_canvas("sim_match_all", output_path)
for id,name in list(fake_types.items()):
	h_sim_match[name].Draw()
	print_canvas("sim_match_%s" % name, output_path)
	h_genid_not_matched[name].Draw()
	print_canvas("h_genid_not_matched_%s" % name, output_path)
	h_sim_relative_pt[name].Draw()
	print_canvas("sim_relative_pt_%s" % name, output_path)
	h_sim_relative_pt_type_matched[name].Draw()
	print_canvas("sim_relative_pt_type_matched_%s" % name, output_path)
	h_sim_relative_pt_wrong_type_matched[name].Draw()
	print_canvas("sim_relative_pt_wrong_type_matched_%s" % name, output_path)
	h_sim_relative_pt_matched_to_other[name].Draw()
	print_canvas("sim_relative_pt_matched_to_other_%s" % name, output_path)
	h_sim_relative_pt_matched_to_punchthrough[name].Draw()
	print_canvas("sim_relative_pt_matched_to_punchthrough_%s" % name, output_path)
	h_sim_decay_rho[name].Draw()
	print_canvas("sim_decay_rho_%s" % name, output_path)
	h_sim_decay_rho_type_matched[name].Draw()
	print_canvas("h_sim_decay_rho_type_matched_%s" % name, output_path)
	h_sim_decay_rho_wrong_type_matched[name].Draw()
	print_canvas("h_sim_decay_rho_wrong_type_matched_%s" % name, output_path)
	h_sim_decay_rho_matched_to_other[name].Draw()
	print_canvas("h_sim_decay_rho_matched_to_other_%s" % name, output_path)
	h_sim_decay_rho_matched_to_other_same_pt[name].Draw()
	print_canvas("h_sim_decay_rho_matched_to_other_same_pt_%s" % name, output_path)
	h_sim_decay_rho_matched_to_other[name].Draw()
	print_canvas("h_sim_decay_rho_matched_to_punchthrough_%s" % name, output_path)

	
# Local Variables:
# indent-tabs-mode: 1
# tab-width: 4
# python-indent: 4
# End:
