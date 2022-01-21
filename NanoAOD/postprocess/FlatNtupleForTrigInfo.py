from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
import multiprocessing
from datetime import datetime
import hashlib

class FlatNtupleForTrigInfo(FlatNtupleBase):
    """Flat ROOT ntuple producer for dimuon trigger studies"""

    triggers = [
        'HLT_DoubleMu4_3_Bs',
        'HLT_DoubleMu4_3_Jpsi',
        'HLT_DoubleMu4_Jpsi_Displaced',
        'HLT_DoubleMu4_Jpsi_NoVertexing',
        'HLT_DoubleMu4_3_Jpsi_Displaced',
        'HLT_Dimuon6_Jpsi_NoVertexing',
        'HLT_DoubleMu4_JpsiTrk_Displaced',       # Jpsi
        'HLT_DoubleMu4_PsiPrimeTrk_Displaced',   # Psi2S
        'HLT_Dimuon0_Upsilon_NoVertexing',       # Upsilon 2017-2018
        'HLT_Dimuon8_Upsilon_Barrel',            # Upsilon 2016
        'L1_DoubleMu0er1p5_SQ_OS',
        'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
        'L1_DoubleMu0er1p6_dEta_Max1p8_OS',
        'HLT_Mu17', 'HLT_Mu8',
        'HLT_IsoMu24', 'HLT_IsoMu27'
    ]
    
    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing informatione
        for parameter in ['input', 'tree_name']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

    def _is_good_dimuon(self, i):
        if self.event.mm_mu1_index[i] < 0 or self.event.mm_mu2_index[i] < 0: return False
	if self.event.Muon_pt[self.event.mm_mu1_index[i]] < 4.0: return False
	if self.event.Muon_pt[self.event.mm_mu2_index[i]] < 4.0: return False
        if abs(self.event.Muon_eta[self.event.mm_mu1_index[i]]) > 1.4: return False
        if abs(self.event.Muon_eta[self.event.mm_mu2_index[i]]) > 1.4: return False
        if self.event.Muon_charge[self.event.mm_mu1_index[i]] * self.event.Muon_charge[self.event.mm_mu2_index[i]] != -1: return False
        if self.event.Muon_softMva[self.event.mm_mu1_index[i]] < 0.45: return False
        if self.event.Muon_softMva[self.event.mm_mu2_index[i]] < 0.45: return False
        return True
    
    def _is_good_for_HLT_DoubleMu4_3_Bs(self, i):
        if self.event.mm_kin_mass[i] < 4.6 or self.event.mm_kin_mass[i] > 5.9: return False 
	if self.event.mm_kin_pt[i] < 5 or self.event.mm_kin_vtx_prob[i] < 0.025: return False
        return True

    def _is_good_for_HLT_DoubleMu4_3_Jpsi(self, i):
        if abs(self.event.mm_kin_mass[i] - 3.1) > 0.1: return False 
	if self.event.mm_kin_pt[i] < 7 or self.event.mm_kin_vtx_prob[i] < 0.1: return False
	if self.event.mm_kin_alphaBS[i] > 0.4: return False
        return True
    
    def _is_good_for_HLT_DoubleMu4_3_Jpsi_Displaced(self, i):
        if not self._is_good_for_HLT_DoubleMu4_3_Jpsi(i): return False
	if self.event.mm_kin_slxy[i] < 4: return False
        return True
    
    def __select_candidates(self, candidates):
        """Select candidates to be stored"""
        
        if len(candidates) < 2:
            return candidates

        # Find the best candidate
        best_candidate = None
        max_value = None
        for i in candidates:
            if not best_candidate or max_value < self.event.mm_kin_pt[i]:
                best_candidate = i
                max_value = self.event.mm_kin_pt[i]

        return [best_candidate]
    
    def _process_events(self):
        """Event loop"""

        for event in self.input_tree:
            self.event = event           
            candidates = []

            for cand in range(self.event.nmm):
                if not self._is_good_dimuon(cand): continue

                if self._is_good_for_HLT_DoubleMu4_3_Bs(cand) or \
                   self._is_good_for_HLT_DoubleMu4_3_Jpsi(cand) or \
                   self._is_good_for_HLT_DoubleMu4_3_Jpsi_Displaced(cand):
                    candidates.append(cand)

            # Find canidates to be stored
            for cand in self.__select_candidates(candidates):
                self._fill_tree(cand)

    def _configure_output_tree(self):

        self.tree.addBranch('run',         'UInt_t', 0)
        self.tree.addBranch('evt',      'ULong64_t', 0)
        self.tree.addBranch('ls',          'UInt_t', 0)
        self.tree.addBranch('PV_npvsGood',  'Int_t', 0)
        self.tree.addBranch('nl1',                'Int_t', 0, "Number of L1 objects with any quality")
        self.tree.addBranch('nl1_SQ',             'Int_t', 0, "Number of L1 objects with single quality")
        
        self.tree.addBranch('pt',         'Float_t', 0, "dimuon pt")
        self.tree.addBranch('eta',        'Float_t', 0, "dimuon eta")
        self.tree.addBranch('m',          'Float_t', 0, "dimuon mass")
        self.tree.addBranch('me',         'Float_t', 0, "mass error of the B candidate")
        self.tree.addBranch('chan',        'UInt_t', 0, "0 if |eta|<0.7 for both muons. 1 otherwise")

        self.tree.addBranch('m1pt',       'Float_t', 0)
        self.tree.addBranch('m1eta',      'Float_t', 0)
        self.tree.addBranch('m1phi',      'Float_t', 0)
        self.tree.addBranch('m2pt',       'Float_t', 0)
        self.tree.addBranch('m2eta',      'Float_t', 0)
        self.tree.addBranch('m2phi',      'Float_t', 0)

        self.tree.addBranch('m1_hlt_pt',  'Float_t', 0)
        self.tree.addBranch('m1_l1_pt',   'Float_t', 0)
        self.tree.addBranch('m1_l1_qual',   'Int_t', 0)
        self.tree.addBranch('m2_hlt_pt',  'Float_t', 0)
        self.tree.addBranch('m2_l1_pt',   'Float_t', 0)
        self.tree.addBranch('m2_l1_qual',   'Int_t', 0)

        self.tree.addBranch('m1_HLT_Mu8',            'Int_t', 0)
        self.tree.addBranch('m1_HLT_Mu17',           'Int_t', 0)
        self.tree.addBranch('m1_HLT_IsoMu24',        'Int_t', 0)
        self.tree.addBranch('m1_HLT_IsoMu27',        'Int_t', 0)
        self.tree.addBranch('m1_HLT_DoubleMu4_3_Bs', 'Int_t', 0)
        self.tree.addBranch('m2_HLT_Mu8',            'Int_t', 0)
        self.tree.addBranch('m2_HLT_Mu17',           'Int_t', 0)
        self.tree.addBranch('m2_HLT_IsoMu24',        'Int_t', 0)
        self.tree.addBranch('m2_HLT_IsoMu27',        'Int_t', 0)
        self.tree.addBranch('m2_HLT_DoubleMu4_3_Bs', 'Int_t', 0)
        
        self.tree.addBranch('good_for_HLT_DoubleMu4_3_Bs',  'UInt_t', 0, "offline selection fits the trigger")
        self.tree.addBranch('good_for_HLT_DoubleMu4_3_Jpsi',  'UInt_t', 0, "offline selection fits the trigger")
        self.tree.addBranch('good_for_HLT_DoubleMu4_3_Jpsi_Displaced',  'UInt_t', 0, "offline selection fits the trigger")
        for trigger in self.triggers:
            self.tree.addBranch(trigger, 'UInt_t', 0)
            self.tree.addBranch("%s_ps" % trigger, 'Int_t', -1, "Prescale. 0 - Off, -1 - no information")

    def _fill_tree(self, cand):
        self.tree.reset()

        ## event info
        self.tree['run'] = self.event.run
        self.tree['ls']  = self.event.luminosityBlock
        self.tree['evt'] = self.event.event
        self.tree['PV_npvsGood'] = self.event.PV_npvsGood

        nl1 = 0
        nl1_SQ = 0
        for index, l1_pt in enumerate(self.event.MuonId_l1_pt):
            if l1_pt >= 0:
                nl1 += 1
                if self.event.MuonId_l1_quality[index] == 12:
                    nl1_SQ += 1
        self.tree['nl1']              = nl1
        self.tree['nl1_SQ']           = nl1_SQ
                    
        self.tree['pt']     = self.event.mm_kin_pt[cand]
        self.tree['eta']    = self.event.mm_kin_eta[cand]
        self.tree['m']      = self.event.mm_kin_mass[cand]
        self.tree['me']     = self.event.mm_kin_massErr[cand]
        
        mu1 = self.event.mm_mu1_index[cand]
        mu2 = self.event.mm_mu2_index[cand]
        
        self.tree['m1pt']  = self.event.Muon_pt[mu1]
        self.tree['m1eta'] = self.event.Muon_eta[mu1]
        self.tree['m1phi'] = self.event.Muon_phi[mu1]
        self.tree['m2pt']  = self.event.Muon_pt[mu2]
        self.tree['m2eta'] = self.event.Muon_eta[mu2]
        self.tree['m2phi'] = self.event.Muon_phi[mu2]
        
        self.tree['m1_hlt_pt']  = self.event.MuonId_hlt_pt[mu1]
        self.tree['m1_l1_pt']   = self.event.MuonId_l1_pt[mu1]
        self.tree['m1_l1_qual'] = self.event.MuonId_l1_quality[mu1]
        self.tree['m2_hlt_pt']  = self.event.MuonId_hlt_pt[mu2]
        self.tree['m2_l1_pt']   = self.event.MuonId_l1_pt[mu2]
        self.tree['m2_l1_qual'] = self.event.MuonId_l1_quality[mu2]

        self.tree['m1_HLT_Mu8']            = self.event.MuonId_HLT_Mu8[mu1]
        self.tree['m1_HLT_Mu17']           = self.event.MuonId_HLT_Mu17[mu1]
        self.tree['m1_HLT_IsoMu24']        = self.event.MuonId_HLT_IsoMu24[mu1]
        self.tree['m1_HLT_IsoMu27']        = self.event.MuonId_HLT_IsoMu27[mu1]
        self.tree['m1_HLT_DoubleMu4_3_Bs'] = self.event.MuonId_HLT_DoubleMu4_3_Bs[mu1]
        self.tree['m2_HLT_Mu8']            = self.event.MuonId_HLT_Mu8[mu2]
        self.tree['m2_HLT_Mu17']           = self.event.MuonId_HLT_Mu17[mu2]
        self.tree['m2_HLT_IsoMu24']        = self.event.MuonId_HLT_IsoMu24[mu2]
        self.tree['m2_HLT_IsoMu27']        = self.event.MuonId_HLT_IsoMu27[mu2]
        self.tree['m2_HLT_DoubleMu4_3_Bs'] = self.event.MuonId_HLT_DoubleMu4_3_Bs[mu2]
        
        if abs(self.tree['m1eta']) < 0.7 and abs(self.tree['m2eta']) < 0.7:
            self.tree['chan'] = 0
        else:
            self.tree['chan'] = 1
            
        self.tree['good_for_HLT_DoubleMu4_3_Bs'] = self._is_good_for_HLT_DoubleMu4_3_Bs(cand)
        self.tree['good_for_HLT_DoubleMu4_3_Jpsi'] = self._is_good_for_HLT_DoubleMu4_3_Jpsi(cand)
        self.tree['good_for_HLT_DoubleMu4_3_Jpsi_Displaced'] = self._is_good_for_HLT_DoubleMu4_3_Jpsi_Displaced(cand)

        for trigger in self.triggers:
            if hasattr(self.event, trigger):
                self.tree[trigger] = getattr(self.event, trigger)
            if hasattr(self.event, "prescale_" + trigger):
                self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
        
        self.tree.fill()

if __name__ == "__main__":

    ### create a test job

    job = {
        "input": [
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/006D33F8-E3CB-C74D-A623-7BF9D2695526.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/EGamma+Run2018D-12Nov2019_UL2018-v4+MINIAOD/8012D764-1747-4D44-9212-1D68EDBED646.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/SingleMuon+Run2018D-12Nov2019_UL2018-v8+MINIAOD/7B0C6A4B-3B98-2342-9925-85DD641A1D6F.root"
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/DoubleMuon+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/8B8F2D5F-097D-4B49-BC79-12B3A1523B64.root"
            ],
        "tree_name" : "mm",
      }  
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = FlatNtupleForMLFit("/tmp/dmytro/9081fd38604b24f4c6035628a89e18ed.job")
    p = FlatNtupleForTrigInfo(file_name)
    print p.__dict__
        
    p.process()
