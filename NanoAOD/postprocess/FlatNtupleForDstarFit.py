from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
import multiprocessing
from datetime import datetime
import hashlib

class FlatNtupleForDstarFit(FlatNtupleBase):
    """Flat ROOT ntuple producer for Bmm5 UML fit"""

    triggers_to_store = [
        # Run 3
        'HLT_ZeroBias',
        'HLT_DoubleMu4_3_LowMass',
    ]

    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing information
        for parameter in ['input', 'blind', 'cut', 'final_state', 'best_candidate']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

    def __select_candidates(self, candidates):
        """Select candidates to be stored"""
        
        if len(candidates) == 0 or self.job_info['best_candidate'] == "":
            return candidates

        # Find the best candidate
        best_candidate = None
        max_value = None
        values = getattr(self.event, self.job_info['best_candidate'])
        for i in candidates:
            if not best_candidate or max_value < values[i]:
                best_candidate = i
                max_value = values[i]

        return [best_candidate]
        
            
    def _process_events(self):
        """Event loop"""

        parsed_cut = self.get_cut()
                    
        for event_index, event in enumerate(self.input_tree):
            self.event = event           
            candidates = []

            # Trigger requirements
            if 'triggers' in self.job_info and len(self.job_info['triggers']) > 0:
                passed_trigger = False
                for trigger in  self.job_info['triggers']:
                    if hasattr(self.event, trigger):
                        if getattr(self.event, trigger):
                            passed_trigger = True
                            break
                if not passed_trigger:
                    continue
            
            # Find candidates the satisfy the selection requirements
            n = self.event.ndstar
            for cand in range(n):
                # if self.job_info['blind']:
                #     if self.job_info['final_state'] == 'mm':
                #         if self.event.mm_kin_mass[cand] < 5.50 and \
                #            self.event.mm_kin_mass[cand] > 5.15:
                #             continue
                #     if self.job_info['final_state'] == 'em':
                #         if self.event.mm_kin_mass[cand] < 5.70 and \
                #            self.event.mm_kin_mass[cand] > 5.10:
                #             continue
                
                # cut = self.job_info['cut'].format(cand=cand, tree="self.event")
                format_dict = { 'ndstar' : cand,
                                'tree' : 'self.event'}
                cut = parsed_cut.format(**format_dict)
                if not eval(cut):
                    continue

                candidates.append(cand)

            # Find canidates to be stored
            cands = self.__select_candidates(candidates)
            for cand in cands:
                self._fill_tree(cand, len(cands))
            

    def _configure_output_tree(self):
        ## event info
        self.tree.addBranch('run',         'UInt_t', 0, "Run number")
        self.tree.addBranch('ls',          'UInt_t', 0, "Luminosity section")
        self.tree.addBranch('evt',      'ULong64_t', 0, "Event number")
        self.tree.addBranch('npv',         'UInt_t', 0, "Number of good reconstructed primary vertices")
        self.tree.addBranch('npu',         'UInt_t', 0, "number of pileup interactions that have been added to the event in the current bunch crossing")
        self.tree.addBranch('npu_mean',   'Float_t', 0, "tru mean number of pileup interactions")
        self.tree.addBranch('certified_muon',    'Int_t', 0, "Event passed Muon Certification")
        self.tree.addBranch('certified_golden',  'Int_t', 0, "Event passed Golden Certification")
        self.tree.addBranch('n',           'UInt_t', 0, "Number of candidates")
        
        self.tree.addBranch('dm',         'Float_t', 0, "m(D*) - m(D0)")
        
        self.tree.addBranch('chan',        'UInt_t', 0, "0: Kpi, 1: pipi, 2: mumu")
        self.tree.addBranch('mc_match',     'Int_t', 0, "PdgId of the MC matched D meson")

        self.tree.addBranch('dstar_vtx_prob', 'Float_t', 0, "D* PV with soft pion probability")
        # self.tree.addBranch('dstar_pt',       'Float_t', 0, "D* pt")
        # self.tree.addBranch('dstar_m',        'Float_t', 0, "D* mass")
        # self.tree.addBranch('dstar_me',       'Float_t', 0, "D* mass error")
        self.tree.addBranch('dstar_pi_pt',    'Float_t', 0, "Soft pion pt")
        self.tree.addBranch('dstar_pi_eta',   'Float_t', 0, "Soft pion eta")
        self.tree.addBranch('dstar_pi_phi',   'Float_t', 0, "Soft pion phi")

        self.tree.addBranch('d0_vtx_prob', 'Float_t', 0, "D0 vtx probability")
        self.tree.addBranch('d0_pt',       'Float_t', 0, "D0 pt")
        self.tree.addBranch('d0_eta',      'Float_t', 0, "D0 eta")
        self.tree.addBranch('d0_phi',      'Float_t', 0, "D0 phi")
        self.tree.addBranch('d0_m',        'Float_t', 0, "D0 mass")
        self.tree.addBranch('d0_me',       'Float_t', 0, "D0 mass error")
        self.tree.addBranch('d0_d1_pt',    'Float_t', 0, "D0 daughter1 pt")
        self.tree.addBranch('d0_d1_eta',   'Float_t', 0, "D0 daughter1 eta")
        self.tree.addBranch('d0_d1_phi',   'Float_t', 0, "D0 daughter1 phi")
        self.tree.addBranch('d0_d2_pt',    'Float_t', 0, "D0 daughter2 pt")
        self.tree.addBranch('d0_d2_eta',   'Float_t', 0, "D0 daughter2 eta")
        self.tree.addBranch('d0_d2_phi',   'Float_t', 0, "D0 daughter2 phi")

        self.tree.addBranch('d0_alpha',    'Float_t', 0, "D0 pointing angle")
        self.tree.addBranch('d0_sl3d',     'Float_t', 0, "D0 significance of flight length 3D")
        
        for trigger in self.triggers_to_store:
            self.tree.addBranch(trigger, 'Int_t', -1, "Trigger decision: 1 - fired, 0 - didn't fire, -1 - no information")
            self.tree.addBranch("%s_ps" % trigger, 'UInt_t', 999999, "Prescale. 0 - Off, 999999 - no information")
            # self.tree.addBranch("%s_matched" % trigger, 'Int_t', 0,  "matched to the trigger objets")

    def _fill_tree(self, cand, ncands):
        self.tree.reset()

        ## event info
        self.tree['run'] = self.event.run
        self.tree['ls']  = self.event.luminosityBlock
        self.tree['evt'] = self.event.event
        self.tree['npv'] = self.event.PV_npvsGood
        if hasattr(self.event, 'Pileup_nTrueInt'):
            self.tree['npu']      = self.event.Pileup_nPU
            self.tree['npu_mean'] = self.event.Pileup_nTrueInt
            self.tree['mc_match'] = self.event.dstar_gen_pdgId[cand]

        self.tree['certified_muon']   = self._is_certified(self.event, "muon")
        self.tree['certified_golden'] = self._is_certified(self.event, "golden")
        self.tree['n']   = ncands

        if self.job_info['final_state'] == 'dzpipi':
            # pipi final state
            if self.event.dstar_hh_index[cand] >= 0:
                hh_index = self.event.dstar_hh_index[cand]
                if self.event.hh_had1_pdgId[hh_index] * self.event.hh_had2_pdgId[hh_index] == - 211 * 211:
                
                    self.tree['chan'] = 1
                    self.tree['dm'] = self.event.dstar_dm_pv[cand]

                    self.tree['dstar_vtx_prob'] = self.event.dstar_pv_with_pion_prob[cand]
                    self.tree['dstar_pi_pt']    = self.event.dstar_pion_pt[cand]
                    self.tree['dstar_pi_eta']   = self.event.dstar_pion_eta[cand]
                    self.tree['dstar_pi_phi']   = self.event.dstar_pion_phi[cand]

                    self.tree['d0_vtx_prob']    = self.event.hh_kin_vtx_prob[hh_index]
                    self.tree['d0_pt']          = self.event.hh_kin_pt[hh_index]
                    self.tree['d0_eta']         = self.event.hh_kin_eta[hh_index]
                    self.tree['d0_phi']         = self.event.hh_kin_phi[hh_index]
                    self.tree['d0_m']           = self.event.hh_kin_mass[hh_index]
                    self.tree['d0_me']          = self.event.hh_kin_massErr[hh_index]
                    self.tree['d0_d1_pt']       = self.event.hh_kin_had1_pt[hh_index]
                    self.tree['d0_d1_eta']      = self.event.hh_kin_had1_eta[hh_index]
                    self.tree['d0_d1_phi']      = self.event.hh_kin_had1_phi[hh_index]
                    self.tree['d0_d2_pt']       = self.event.hh_kin_had2_pt[hh_index]
                    self.tree['d0_d2_eta']      = self.event.hh_kin_had2_eta[hh_index]
                    self.tree['d0_d2_phi']      = self.event.hh_kin_had2_phi[hh_index]

                    self.tree['d0_alpha']       = self.event.hh_kin_alpha[hh_index]
                    self.tree['d0_sl3d']        = self.event.hh_kin_sl3d[hh_index]
        elif self.job_info['final_state'] == 'dzmm':
            # mm final state
            if self.event.dstar_mm_index[cand] >= 0:
        
                self.tree['chan'] = 2
                mm_index = self.event.dstar_mm_index[cand]

                self.tree['dm'] = self.event.dstar_dm_pv[cand]

                self.tree['dstar_vtx_prob'] = self.event.dstar_pv_with_pion_prob[cand]
                self.tree['dstar_pi_pt']    = self.event.dstar_pion_pt[cand]
                self.tree['dstar_pi_eta']   = self.event.dstar_pion_eta[cand]
                self.tree['dstar_pi_phi']   = self.event.dstar_pion_phi[cand]

                self.tree['d0_vtx_prob']    = self.event.mm_kin_vtx_prob[mm_index]
                self.tree['d0_pt']          = self.event.mm_kin_pt[mm_index]
                self.tree['d0_eta']         = self.event.mm_kin_eta[mm_index]
                self.tree['d0_phi']         = self.event.mm_kin_phi[mm_index]
                self.tree['d0_m']           = self.event.mm_kin_mass[mm_index]
                self.tree['d0_me']          = self.event.mm_kin_massErr[mm_index]
                self.tree['d0_d1_pt']       = self.event.mm_kin_mu1_pt[mm_index]
                self.tree['d0_d1_eta']      = self.event.mm_kin_mu1_eta[mm_index]
                self.tree['d0_d1_phi']      = self.event.mm_kin_mu1_phi[mm_index]
                self.tree['d0_d2_pt']       = self.event.mm_kin_mu2_pt[mm_index]
                self.tree['d0_d2_eta']      = self.event.mm_kin_mu2_eta[mm_index]
                self.tree['d0_d2_phi']      = self.event.mm_kin_mu2_phi[mm_index]

                self.tree['d0_alpha']       = self.event.mm_kin_alpha[mm_index]
                self.tree['d0_sl3d']        = self.event.mm_kin_sl3d[mm_index]
        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        for trigger in self.triggers_to_store:
            if hasattr(self.event, trigger):
                self.tree[trigger] = getattr(self.event, trigger)
            if hasattr(self.event, "prescale_" + trigger):
                self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
        
        self.tree.fill()

if __name__ == "__main__":

    ### create a test job
    
    # job = {
    #     "input": [
    #         # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/089fcee0-0260-470b-8f08-a458129f2c4a.root"
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "dzmmMC",
    #     "blind" : False,
    #     "triggers":['HLT_DoubleMu4_3_LowMass'],
    #     "cut" :
    #         "dstar_mm_index>=0 and Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 and "\
    #         "Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45 and "\
    #         "mm_mu1_pt[dstar_mm_index]>4 and mm_mu2_pt[dstar_mm_index]>4 and "\
    #         "mm_kin_alpha[dstar_mm_index]<0.1 and mm_kin_sl3d[dstar_mm_index]>3 and "\
    #         "mm_kin_vtx_prob[dstar_mm_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "\
    #         "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "\
    #         "mm_kin_mass[dstar_mm_index]>1.81 and mm_kin_mass[dstar_mm_index]<1.94",
    #     "final_state" : "dstar",
    #     "best_candidate": "",
    #   }

    # job = {
    #     "input": [
    #         'root://eoscms.cern.ch://eos/cms/store/group/phys_muon/dmytro/tmp/DstarToD0Pi_D0To2Pi_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen.root'
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "dzpipiMC",
    #     "blind" : False,
    #     "triggers":[],
    #     "cut" : (
    #         "dstar_hh_index>=0 and "
    #         "hh_had1_pt[dstar_hh_index]>4 and hh_had2_pt[dstar_hh_index]>4 and "
    #         "hh_kin_alpha[dstar_hh_index]<0.1 and hh_kin_sl3d[dstar_hh_index]>3 and "
    #         "hh_kin_vtx_prob[dstar_hh_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
    #         "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
    #         "hh_kin_mass[dstar_hh_index]>1.81 and hh_kin_mass[dstar_hh_index]<1.94 and "
    #         "hh_had1_pdgId[dstar_hh_index] * hh_had2_pdgId[dstar_hh_index] == - 211 * 211"
    #     ),
    #     "final_state" : "dzpipi",
    #     "best_candidate": "",
    #   }

    # job = {
    #     "input": [
    #         'root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2+MINIAODSIM/1211fb7f-8f51-41fb-87f4-266a5db414e5.root',
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "dzmmMC",
    #     "blind" : False,
    #     "triggers":[],
    #     "cut" : (
    #         "dstar_mm_index>=0 and "
    #         "mm_mu1_pt[dstar_mm_index]>4 and mm_mu2_pt[dstar_mm_index]>4 and "
    #         "mm_kin_alpha[dstar_mm_index]<0.1 and mm_kin_sl3d[dstar_mm_index]>3 and "
    #         "mm_kin_vtx_prob[dstar_mm_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
    #         "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
    #         "mm_kin_mass[dstar_mm_index]>1.81 and mm_kin_mass[dstar_mm_index]<1.94"
    #     ),
    #     "final_state" : "dzmm",
    #     "best_candidate": "",
    #   }

    common_branches = 'PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock'
    
    job = {
        "input": [
            'root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523/ParkingDoubleMuonLowMass5+Run2022C-PromptReco-v1+MINIAOD/3237cbb0-122b-4d28-bebf-dea5da307147.root',
        ],
        "signal_only" : False,
        "tree_name" : "dzmmData",
        "blind" : False,
        "triggers":["HLT_DoubleMu4_3_LowMass"],
        "pre-selection":"dstar_mm_index>=0 && dstar_dm_pv>0.140 && dstar_dm_pv<0.155",
        "pre-selection-keep":"^(dstar_.*|ndstar|mm_.*|nmm|HLT_DoubleMu4_3_LowMass|" + common_branches + ")$",
        "cut" : (
            "dstar_mm_index>=0 and "
            "mm_mu1_pt[dstar_mm_index]>4 and mm_mu2_pt[dstar_mm_index]>4 and "
            "mm_kin_alpha[dstar_mm_index]<0.1 and mm_kin_sl3d[dstar_mm_index]>3 and "
            "mm_kin_vtx_prob[dstar_mm_index]>0.01 and dstar_pv_with_pion_prob>0.1 and "
            "dstar_dm_pv>0.140 and dstar_dm_pv<0.155 and "
            "mm_kin_mass[dstar_mm_index]>1.81 and mm_kin_mass[dstar_mm_index]<1.94"
        ),
        "final_state" : "dzmm",
        "best_candidate": "",
      }
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    p = FlatNtupleForDstarFit("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/523/dstar/ParkingDoubleMuonLowMass3+Run2022F-PromptReco-v1+MINIAOD/058a22ecf4cf8182f9ce50ef6afa4971.job")
    # p = FlatNtupleForDstarFit(file_name)

    print(p.__dict__)
        
    p.process()
