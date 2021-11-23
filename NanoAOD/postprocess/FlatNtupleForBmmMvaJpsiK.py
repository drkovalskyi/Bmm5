from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
from Bmm5.MVA.mtree import MTree
import multiprocessing
from datetime import datetime
import hashlib
from ROOT import TFile, TTree

class FlatNtupleForBmmMvaJpsiK(FlatNtupleBase):
    """Produce flat ROOT ntuples for Bmm5 MVA training using JpsiK events reconstructed like Bmm"""

    def _validate_inputs(self):
        """Task specific input validation"""
        # check for missing information
        for parameter in ['input', 'signal_only', 'tree_name']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

    def _process_events(self):
        """Event loop"""

        for event_index, event in enumerate(self.input_tree):
            if self.limit > 0 and event_index >= self.limit: break

            self.event = event
            candidates = self._select_candidates()
            # _analyze_selection(info,statisitics)
            for cand in self._good_candidates(candidates):
                if hasattr(self.event, 'bkmm_gen_pdgId'):
                    if self.job_info['signal_only'] and self.event.bkmm_gen_pdgId[cand]==0: continue
                self._fill_tree(cand)

    def _select_candidates(self):
        candidates = []
        if not self.event.nmm>0: return candidates
        for i in range(self.event.nbkmm):
            selection = []
            mu1_index = self.event.mm_mu1_index[self.event.bkmm_mm_index[i]]
            mu2_index = self.event.mm_mu2_index[self.event.bkmm_mm_index[i]]
            selection.append( ('good_muon1', self._good_muon(mu1_index)) )
            selection.append( ('good_muon2', self._good_muon(mu2_index)) )
            selection.append( ('mass_cut',  abs(self.event.bkmm_jpsimc_mass[i] - 5.4) < 0.5) )
            selection.append( ('valid_fit', self.event.bkmm_jpsimc_valid[i]) )
            # selection.append( ('kaon_impact_parameter_significance', self.event.bkmm_kaon_sdxy_bs[i] > 15) )
            selection.append( ('mm_vtx_prob', self.event.mm_kin_vtx_prob[self.event.bkmm_mm_index[i]] > 0.1) )
            selection.append( ('mmk_vtx_prob', self.event.bkmm_jpsimc_vtx_prob[i] > 0.025) )
            # trigger driven, but tighter for purity
            selection.append( ('decay_length_significance', self.event.bkmm_jpsimc_sl3d[i] > 4.0) )
            selection.append( ('alpha', self.event.bkmm_jpsimc_alpha[i] < 0.2) )
            # phase space tunning
            selection.append( ('kaon_pt',  self.event.bkmm_kaon_pt[i] < 1.5))
            q1 = 0
            if mu1_index >= 0: q1 = self.event.Muon_charge[mu1_index]
            q2 = 0
            if mu2_index >= 0: q2 = self.event.Muon_charge[mu2_index]
            selection.append( ('charge', q1*q2 == -1) )

            candidates.append(selection)
        return candidates

    def _good_candidates(self, candidates):
        good_candidates = []
        for i in range(len(candidates)):
            good_candidate = True
            for selection in candidates[i]:
                if not selection[1]:
                    good_candidate = False
                    break
            if good_candidate: good_candidates.append(i)
        return good_candidates

    def _configure_output_tree(self):
        ## event info
        self.tree.addBranch( 'evt_run',           'UInt_t', 0)
        self.tree.addBranch( 'evt_lumi',          'UInt_t', 0)
        self.tree.addBranch( 'evt_event',         'ULong64_t', 0)
        self.tree.addBranch( 'evt_nvtx',          'UInt_t',  0, "Number of good reconstructed primary vertices")

        ## mm info
        self.tree.addBranch('mm_kin_mass',        'Float_t', 0, "Vertex constrained dimuon mass")
        self.tree.addBranch('mm_kin_eta',         'Float_t', 0, "Vertex constrained dimuon eta")
        self.tree.addBranch('mm_kin_pt',          'Float_t', 0, "Vertex constrained dimuon pt")
        self.tree.addBranch('mm_mu1_pt',          'Float_t', 0)
        self.tree.addBranch('mm_mu2_pt',          'Float_t', 0)
        self.tree.addBranch('mm_mu1_eta',         'Float_t', 0)
        self.tree.addBranch('mm_mu2_eta',         'Float_t', 0)
        self.tree.addBranch('mm_mu1_phi',         'Float_t', 0)
        self.tree.addBranch('mm_mu2_phi',         'Float_t', 0)
        self.tree.addBranch('mm_kin_l3d',         'Float_t', 0, "Decay length wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_kin_sl3d',        'Float_t', 0, "Decay length significance wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_kin_slxy',        'Float_t', 0, "Decay length significance wrt Beam Spot in XY plain")

        self.tree.addBranch('mm_kin_alpha',       'Float_t', 0, "Pointing angle in 3D wrt PV")
        self.tree.addBranch('mm_kin_alphaErr',    'Float_t', 0, "Pointing angle in 3D wrt PV uncertainty")
        self.tree.addBranch('mm_kin_alphaSig',    'Float_t', 999, "Pointing angle in 3D wrt PV significance")
        self.tree.addBranch('mm_kin_alphaBS',     'Float_t', 0, "Cosine of pointing angle in XY wrt BS")
        self.tree.addBranch('mm_kin_alphaBSErr',  'Float_t', 0, "Cosine of pointing angle in XY wrt BS uncertainty")
        self.tree.addBranch('mm_kin_alphaBSSig',  'Float_t', 999, "Cosine of pointing angle in XY wrt BS significance")
        self.tree.addBranch('mm_kin_spvip',       'Float_t', 0, "Significance of impact parameter wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_kin_pvip',        'Float_t', 0, "Impact parameter wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_iso',             'Float_t', 0, "B isolation the way it's done in Bmm4")
        self.tree.addBranch('mm_kin_vtx_chi2dof', 'Float_t', 0, "Normalized chi2 of the dimuon vertex kinematic fit")
        self.tree.addBranch('mm_otherVtxMaxProb1','Float_t', 0, "Max vertexing probability of one of the muons with a random track with minPt=1.0GeV")
        self.tree.addBranch('mm_otherVtxMaxProb2','Float_t', 0, "Max vertexing probability of one of the muons with a random track with minPt=2.0GeV")
        self.tree.addBranch('mm_m1iso',           'Float_t', 0, "Muon isolation the way it's done in Bmm4")
        self.tree.addBranch('mm_m2iso',           'Float_t', 0, "Muon isolation the way it's done in Bmm4")
        self.tree.addBranch('mm_nBMTrks',         'UInt_t',  0, "Number of tracks more compatible with the mm vertex than with PV by doca significance")

        self.tree.addBranch('trigger',            'UInt_t',  0, "Main analysis trigger")
        self.tree.addBranch('mm_mva',             'Float_t', 0, "MVA")
        self.tree.addBranch('bhad_pt',            'Float_t', 0, "B hadron pt for MC matched decays")
  
    def _fill_tree(self, cand):
        self.tree.reset()

        ## event info
        self.tree['evt_run']   = self.event.run
        self.tree['evt_lumi']  = self.event.luminosityBlock
        self.tree['evt_event'] = self.event.event
        self.tree['evt_nvtx']  = self.event.PV_npvsGood

        ## mm info
        self.tree['mm_kin_mass']    = self.event.bkmm_jpsimc_mass[cand]
        self.tree['mm_kin_eta']     = self.event.bkmm_jpsimc_eta[cand]
        self.tree['mm_kin_pt']      = self.event.bkmm_jpsimc_pt[cand]
        self.tree['mm_mu1_pt']      = self.event.Muon_pt[self.event.mm_mu1_index[self.event.bkmm_mm_index[cand]]]
        self.tree['mm_mu2_pt']      = self.event.Muon_pt[self.event.mm_mu2_index[self.event.bkmm_mm_index[cand]]]
        self.tree['mm_mu1_eta']     = self.event.Muon_eta[self.event.mm_mu1_index[self.event.bkmm_mm_index[cand]]]
        self.tree['mm_mu2_eta']     = self.event.Muon_eta[self.event.mm_mu2_index[self.event.bkmm_mm_index[cand]]]
        self.tree['mm_mu1_phi']     = self.event.Muon_phi[self.event.mm_mu1_index[self.event.bkmm_mm_index[cand]]]
        self.tree['mm_mu2_phi']     = self.event.Muon_phi[self.event.mm_mu2_index[self.event.bkmm_mm_index[cand]]]
        self.tree['mm_kin_sl3d']    = self.event.bkmm_jpsimc_sl3d[cand]*1.6

        self.tree['mm_kin_alpha']   = self.event.bkmm_jpsimc_alpha[cand]
        self.tree['mm_kin_alphaErr'] = self.event.bkmm_jpsimc_alphaErr[cand]
        if self.event.bkmm_jpsimc_alphaErr[cand] > 0:
            self.tree['mm_kin_alphaSig'] = self.event.bkmm_jpsimc_alpha[cand] / self.event.bkmm_jpsimc_alphaErr[cand]
        self.tree['mm_kin_alphaBS'] = self.event.bkmm_jpsimc_alphaBS[cand]
        self.tree['mm_kin_alphaBSErr'] = self.event.bkmm_jpsimc_alphaBSErr[cand]
        if self.event.bkmm_jpsimc_alphaBSErr[cand] > 0:
            self.tree['mm_kin_alphaBSSig'] = self.event.bkmm_jpsimc_alphaBS[cand] / self.event.bkmm_jpsimc_alphaBSErr[cand]
        if self.event.bkmm_jpsimc_pvipErr[cand] > 0:
            self.tree['mm_kin_spvip'] = self.event.bkmm_jpsimc_pvip[cand]/self.event.bkmm_jpsimc_pvipErr[cand]
        self.tree['mm_kin_pvip']    = self.event.bkmm_jpsimc_pvip[cand]
        self.tree['mm_iso']         = self.event.bkmm_bmm_iso[cand]
        self.tree['mm_kin_vtx_chi2dof'] = self.event.mm_kin_vtx_chi2dof[self.event.bkmm_mm_index[cand]]
        self.tree['mm_otherVtxMaxProb1'] = self.event.bkmm_bmm_otherVtxMaxProb1[cand]
        self.tree['mm_otherVtxMaxProb2'] = self.event.bkmm_bmm_otherVtxMaxProb2[cand]
        self.tree['mm_m1iso']       = self.event.bkmm_bmm_m1iso[cand]
        self.tree['mm_m2iso']       = self.event.bkmm_bmm_m2iso[cand]
        self.tree['mm_nBMTrks']     = self.event.bkmm_bmm_nBMTrks[cand]
        self.tree['mm_mva']         = self.event.bkmm_bmm_mva[cand]
        

        trigger_2018 = 0
        if hasattr(self.event, 'HLT_DoubleMu4_3_Jpsi'):
            trigger_2018 = self.event.HLT_DoubleMu4_3_Jpsi

        trigger_2017 = 0
        trigger_2016 = 0
        if hasattr(self.event, 'HLT_DoubleMu4_3_Jpsi_Displaced'):
            trigger_2017 = self.event.HLT_DoubleMu4_3_Jpsi_Displaced
            trigger_2016 = self.event.HLT_DoubleMu4_3_Jpsi_Displaced

        self.tree['trigger'] = trigger_2018 | trigger_2017 | trigger_2016

        if hasattr(self.event, 'bkmm_gen_pdgId'):
            if self.job_info['signal_only'] and self.event.bkmm_gen_pdgId[cand]!=0:
                self.tree['bhad_pt'] = self.event.bkmm_gen_pt[cand]
        
        self.tree.fill()

    def _good_muon(self, index):
        if index < 0: return False
        if not (abs(self.event.Muon_eta[index]) < 1.4): return False
        if not (self.event.Muon_pt[index] > 4): return False
        if not self.event.Muon_softId[index]: return False
        return True

def unit_test():
    path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/"
    job = {
        "input": [
            path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/3161083A-28F2-244D-AEE2-0DD084D47E23.root"
        ],
        "signal_only": True,
        "tree_name": "mva",
    }

    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    p = FlatNtupleForBmmMvaJpsiK(file_name)
    p.process()

if __name__ == "__main__":
    unit_test()
