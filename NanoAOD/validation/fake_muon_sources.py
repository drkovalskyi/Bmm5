import ROOT
import sys, os
from DataFormats.FWLite import Events, Handle
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from math import *

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/muon_fake_sources"
dump_info = False
min_pt = 10

events = Events (
	[
		# '/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/muon_fake_skim.root'
		'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_RunIIAutumn18MiniAOD_muon_fake_skim.root',
		'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/test/LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM.root'
	]
)

def isAncestor(a,p) :
	if a == p : 
		return True
	for i in xrange(0,p.numberOfMothers()) :
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
	if not muon.isLooseMuon(): return False;
	if not muon.innerTrack().quality(ROOT.reco.Track.highPurity): return False 
	if muon.pt() < min_pt or abs(muon.eta()) > 1.4: return False
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
	return 'other'

sim_type_to_value = {
	'not_matched':0.5,
	'pion':1.5,
	'kaon':2.5,
	'proton':3.5,
	'other':4.5
}


ROOT.gROOT.SetBatch(True)

handlePruned  = Handle ("std::vector<reco::GenParticle>")
handlePacked  = Handle ("std::vector<pat::PackedGenParticle>")
labelPruned = ("prunedGenParticles")
labelPacked = ("packedGenParticles")

muonHandle, muonLabel = Handle("std::vector<pat::Muon>"),"slimmedMuons"

h_sim_match = dict()
h_sim_relative_pt = dict()
h_sim_relative_pt_type_matched = dict()
h_sim_decay_rho = dict()
h_sim_decay_rho_type_matched = dict()
for id,name in fake_types.items():
	h_sim_match[name] = ROOT.TH1D("h_sim_match_%s" % name,
								  "Muon matching based on simulated hits", 5, 0, 5)
	h_sim_match[name].GetXaxis().SetBinLabel(1, "No Match")
	h_sim_match[name].GetXaxis().SetBinLabel(2, "Muon from Pion")
	h_sim_match[name].GetXaxis().SetBinLabel(3, "Muon from Kaon")
	h_sim_match[name].GetXaxis().SetBinLabel(4, "Muon from Proton")
	h_sim_match[name].GetXaxis().SetBinLabel(5, "Other")

	h_sim_relative_pt[name] = ROOT.TH1D("h_sim_relative_pt_%s" % name,
										"Relative Pt of sim-matched particle", 100, 0, 2)
	h_sim_relative_pt_type_matched[name] = ROOT.TH1D("h_sim_relative_pt_type_matched_%s" % name,
										"Relative Pt of sim particle (type matched)", 100, 0, 2)
	h_sim_decay_rho[name] = ROOT.TH1D("h_sim_decay_rho_%s" % name,
									  "Decay radius", 100, 0, 400)
	h_sim_decay_rho_type_matched[name] = ROOT.TH1D("h_sim_decay_rho_type_matched_%s" % name,
												   "Decay radius for matched types", 100, 0, 400)
	
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

	# check if we have a match

	event.getByLabel (labelPacked, handlePacked)
	event.getByLabel (labelPruned, handlePruned)
	# get the product
	packed = handlePacked.product()
	pruned = handlePruned.product()

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
			
		interesting_event = True
		
		status = p.statusFlags().fromHardProcessBeforeFSR()
		if dump_info:
			print "Reco Muon pt : %s  eta : %s   phi: %s" % (muon.pt(), muon.eta(), muon.phi())
			print "\tsim_pt: %s  sim_pid: %s sim_mother_pid: %s sim_prod_rho: %f sim_prod_z: %f" % (muon.simPt(), muon.simPdgId(), muon.simMotherPdgId(),
																									muon.simProdRho(), muon.simProdZ())
			print "Gen pdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(p.pdgId(), p.pt(), p.eta(), p.phi(), status)    
			mother = find_parent(p)
			if mother:
				status = mother.statusFlags().fromHardProcessBeforeFSR()
				print "\tMother PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(mother.pdgId(),mother.pt(),mother.eta(),mother.phi(),status)    
				mother = find_parent(mother)
				if mother:
					status = mother.statusFlags().fromHardProcessBeforeFSR()
					print "\tGrand Mother PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(mother.pdgId(),mother.pt(),mother.eta(),mother.phi(),status)    

	if not interesting_event: continue

	count += 1

	# https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
	if dump_info:
		print "Event dump"
		print "Run: %u, Event: %u" % (event.eventAuxiliary().run(),event.eventAuxiliary().event())
		for p in pruned :
			if not abs(p.pdgId()) in [511,521,531]: continue

			final_b = True
			for dau in p.daughterRefVector():
				if dau.pdgId() == -p.pdgId():
					final_b = False
					break
			if not final_b: continue

			signature = 1
			for pa in packed:
				mother = pa.mother(0)
				if mother and isAncestor(p,mother) :
					if pa.pdgId()!=22: signature *= pa.pdgId()

			# if not signature in [13*13*321*321, 13*13, 13*13*321, -13*13*321]: continue

			# d_p4 = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D(ROOT.Double))()
			d_p4 = ROOT.reco.Candidate.LorentzVector()
			rad_p4 = ROOT.reco.Candidate.LorentzVector()
			status = p.statusFlags().fromHardProcessBeforeFSR()
			print "     PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(p.pdgId(),p.pt(),p.eta(),p.phi(),status)    
			for dau in p.daughterRefVector():
				status = dau.statusFlags().fromHardProcessBeforeFSR()
				print "     dau PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(dau.pdgId(),dau.pt(),dau.eta(),dau.phi(),status)    
			if p.mother():
				status = p.mother().statusFlags().fromHardProcessBeforeFSR()
				print "     mother PdgId : %s   pt : %s  eta : %s   phi : %s status: %d" %(p.mother().pdgId(),p.mother().pt(),p.mother().eta(),p.mother().phi(),status)    
			for pa in packed:
				mother = pa.mother(0)
				if mother and isAncestor(p,mother) :
					print "          PdgId : %s   pt : %s  eta : %s   phi : %s" %(pa.pdgId(),pa.pt(),pa.eta(),pa.phi())
					d_p4 += pa.p4()
					if pa.pdgId()==22: rad_p4 += pa.p4()
			print "     delta: %0.5f%%" % (100.*(p.p4()-d_p4).P()/p.p4().P())
			print "     radiation: %0.2f%%" % (100.*rad_p4.P()/p.p4().P())
	if dump_info and count >= 100: 
		sys.exit()

c1 = TCanvas("c1", "c1", 800, 800)
for id,name in fake_types.items():
	h_sim_match[name].Draw()
	print_canvas("sim_match_%s" % name, output_path)
	h_sim_relative_pt[name].Draw()
	print_canvas("sim_relative_pt_%s" % name, output_path)
	h_sim_relative_pt_type_matched[name].Draw()
	print_canvas("sim_relative_pt_type_matched_%s" % name, output_path)
	h_sim_decay_rho[name].Draw()
	print_canvas("sim_decay_rho_%s" % name, output_path)
	h_sim_decay_rho_type_matched[name].Draw()
	print_canvas("h_sim_decay_rho_type_matched_%s" % name, output_path)

	
# Local Variables:
# indent-tabs-mode: 1
# tab-width: 4
# python-indent: 4
# End:
