from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
import multiprocessing
from datetime import datetime
import hashlib

class FlatNtupleForMLFit(FlatNtupleBase):
    """Flat ROOT ntuple producer for Bmm5 UML fit"""

    leaf_counts = { 
        'mm':    'nmm',
        'bkmm':  'nbkmm',
        'bkkmm': 'nbkkmm'
    }

    triggers = [
        'HLT_DoubleMu4_3_Bs',
        'HLT_DoubleMu4_3_Jpsi',
        'HLT_DoubleMu4_Jpsi_Displaced',
        'HLT_DoubleMu4_Jpsi_NoVertexing',
        'HLT_DoubleMu4_3_Jpsi_Displaced',
        'HLT_Dimuon6_Jpsi_NoVertexing'
    ]

    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing information
        for parameter in ['input', 'blind', 'cut', 'final_state', 'best_candidate']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

    # def __get(self, name, index=None):
    #     """Return branch value"""

    #     var = self.branch_map[self.job_info['final_state']][name]

    #     return eval(var.format(tree="self.event", cand=index))
            
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

            # Find candidates the satisfy the selection requirements
            n = getattr(self.event, self.leaf_counts[self.job_info['final_state']])
            for cand in range(n):
                if self.job_info['final_state'] == 'mm':
                    if self.job_info['blind']:
                        if self.event.mm_kin_mass[cand] < 5.50 and \
                           self.event.mm_kin_mass[cand] > 5.15:
                            continue
                
                # cut = self.job_info['cut'].format(cand=cand, tree="self.event")
                format_dict = { self.leaf_counts[self.job_info['final_state']] : cand,
                                'tree' : 'self.event'}
                cut = parsed_cut.format(**format_dict)
                if not eval(cut):
                    continue

                candidates.append(cand)

            # Find canidates to be stored
            for cand in self.__select_candidates(candidates):
                self._fill_tree(cand)
            

    def _configure_output_tree(self):
        ## event info
        self.tree.addBranch('run',         'UInt_t', 0)
        self.tree.addBranch('evt',      'ULong64_t', 0)
        self.tree.addBranch('ls',          'UInt_t', 0)
        # ps
        # cw8
        self.tree.addBranch('bdt',        'Float_t', 0, "BDT score")
        self.tree.addBranch('pt',         'Float_t', 0, "pt of the B candidate")
        self.tree.addBranch('eta',        'Float_t', 0, "eta of the B candidate")
        self.tree.addBranch('phi',        'Float_t', 0, "phi of the B candidate")
        self.tree.addBranch('m',          'Float_t', 0, "mass of the B candidate")
        self.tree.addBranch('me',         'Float_t', 0, "mass error of the B candidate")
        self.tree.addBranch('tau',        'Float_t', 0, "decay time of the B candidate")
        self.tree.addBranch('taue',       'Float_t', 0, "decay time error of the B candidate")
        self.tree.addBranch('tauxy',      'Float_t', 0, "2D decay time")
        self.tree.addBranch('tauxye',     'Float_t', 0, "2D decay time error")
        self.tree.addBranch('gtau',       'Float_t', 0, "generated decay time")

        self.tree.addBranch('muid',        'UInt_t', 0, "MVA muon selection for both muons")
        self.tree.addBranch('chan',        'UInt_t', 0, "0 if |eta|<0.7 for both muons. 1 otherwise")

        self.tree.addBranch('m1pt',       'Float_t', 0)
        self.tree.addBranch('m1eta',      'Float_t', 0)
        self.tree.addBranch('m1phi',      'Float_t', 0)
        self.tree.addBranch('m1q',          'Int_t', 0, "Muon charge. 0 for hadrons")
        self.tree.addBranch('m1bdt',      'Float_t', 0, "Muon BDT. 1 for hadrons")

        self.tree.addBranch('m2pt',       'Float_t', 0)
        self.tree.addBranch('m2eta',      'Float_t', 0)
        self.tree.addBranch('m2phi',      'Float_t', 0)
        self.tree.addBranch('m2q',          'Int_t', 0, "Muon charge. 0 for hadrons")
        self.tree.addBranch('m2bdt',      'Float_t', 0, "Muon BDT. 1 for hadrons")

        self.tree.addBranch('mc_match',     'Int_t', 0, "PdgId of the MC matched B meson")

        for trigger in self.triggers:
            self.tree.addBranch(trigger, 'UInt_t', 0)
            self.tree.addBranch("%s_ps" % trigger, 'UInt_t', 999999, "Prescale. 0 - Off, 999999 - no information")
        
    def _fill_tree(self, cand):
        self.tree.reset()

        ## event info
        self.tree['run'] = self.event.run
        self.tree['ls']  = self.event.luminosityBlock
        self.tree['evt'] = self.event.event

        if self.job_info['final_state'] == 'mm':
            # B to mm
            self.tree['bdt']    = self.event.mm_mva[cand]
            self.tree['pt']     = self.event.mm_kin_pt[cand]
            self.tree['eta']    = self.event.mm_kin_eta[cand]
            self.tree['phi']    = self.event.mm_kin_phi[cand]
            self.tree['m']      = self.event.mm_kin_mass[cand]
            self.tree['me']     = self.event.mm_kin_massErr[cand]
            self.tree['tau']    = self.event.mm_kin_tau[cand]
            self.tree['taue']   = self.event.mm_kin_taue[cand]
            self.tree['tauxy']  = self.event.mm_kin_tauxy[cand]
            self.tree['tauxye'] = self.event.mm_kin_tauxye[cand]
            if hasattr(self.event, 'mm_gen_tau'):
                self.tree['gtau'] = self.event.mm_gen_tau[cand]
                self.tree['mc_match'] = self.event.mm_gen_pdgId[cand]

            mu1 = self.event.mm_mu1_index[cand]
            mu2 = self.event.mm_mu2_index[cand]
        
            if mu1 >= 0:
                self.tree['m1pt']  = self.event.Muon_pt[mu1]
                self.tree['m1eta'] = self.event.Muon_eta[mu1]
                self.tree['m1phi'] = self.event.Muon_phi[mu1]
                self.tree['m1q']   = self.event.Muon_charge[mu1]
                self.tree['m1bdt'] = self.event.Muon_softMva[mu1]
            else: 
                self.tree['m1pt']  = self.event.mm_mu1_pt[cand]
                self.tree['m1eta'] = self.event.mm_mu1_eta[cand]
                self.tree['m1phi'] = self.event.mm_mu1_phi[cand]
                self.tree['m1bdt'] = 1

            if mu2 >= 0:
                self.tree['m2pt']  = self.event.Muon_pt[mu2]
                self.tree['m2eta'] = self.event.Muon_eta[mu2]
                self.tree['m2phi'] = self.event.Muon_phi[mu2]
                self.tree['m2q']   = self.event.Muon_charge[mu2]
                self.tree['m2bdt'] = self.event.Muon_softMva[mu2]
            else: 
                self.tree['m2pt']  = self.event.mm_mu2_pt[cand]
                self.tree['m2eta'] = self.event.mm_mu2_eta[cand]
                self.tree['m2phi'] = self.event.mm_mu2_phi[cand]
                self.tree['m2bdt'] = 1

            if mu1 >= 0 and mu2 >= 0:
                self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
                
        elif self.job_info['final_state'] == 'bkmm':
            # B to Jpsi K
            self.tree['pt']     = self.event.bkmm_jpsimc_pt[cand]
            self.tree['eta']    = self.event.bkmm_jpsimc_eta[cand]
            self.tree['phi']    = self.event.bkmm_jpsimc_phi[cand]
            self.tree['m']      = self.event.bkmm_jpsimc_mass[cand]
            self.tree['me']     = self.event.bkmm_jpsimc_massErr[cand]
            self.tree['tau']    = self.event.bkmm_jpsimc_tau[cand]
            self.tree['taue']   = self.event.bkmm_jpsimc_taue[cand]
            self.tree['tauxy']  = self.event.bkmm_jpsimc_tauxy[cand]
            self.tree['tauxye'] = self.event.bkmm_jpsimc_tauxye[cand]
            if hasattr(self.event, 'bkmm_gen_tau'):
                self.tree['gtau'] = self.event.bkmm_gen_tau[cand]
                self.tree['mc_match'] = self.event.bkmm_gen_pdgId[cand]

            mu1 = self.event.mm_mu1_index[self.event.bkmm_mm_index[cand]]
            mu2 = self.event.mm_mu2_index[self.event.bkmm_mm_index[cand]]
        
            if mu1 >= 0:
                self.tree['m1pt']  = self.event.Muon_pt[mu1]
                self.tree['m1eta'] = self.event.Muon_eta[mu1]
                self.tree['m1phi'] = self.event.Muon_phi[mu1]
                self.tree['m1q']   = self.event.Muon_charge[mu1]
                self.tree['m1bdt'] = self.event.Muon_softMva[mu1]

            if mu2 >= 0:
                self.tree['m2pt']  = self.event.Muon_pt[mu2]
                self.tree['m2eta'] = self.event.Muon_eta[mu2]
                self.tree['m2phi'] = self.event.Muon_phi[mu2]
                self.tree['m2q']   = self.event.Muon_charge[mu2]
                self.tree['m2bdt'] = self.event.Muon_softMva[mu2]

            if mu1 >= 0 and mu2 >= 0:
                self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
                
        elif self.job_info['final_state'] == 'bkkmm':
            # B to Jpsi Phi 
            self.tree['pt']     = self.event.bkkmm_jpsikk_pt[cand]
            self.tree['eta']    = self.event.bkkmm_jpsikk_eta[cand]
            self.tree['phi']    = self.event.bkkmm_jpsikk_phi[cand]
            self.tree['m']      = self.event.bkkmm_jpsikk_mass[cand]
            self.tree['me']     = self.event.bkkmm_jpsikk_massErr[cand]

            self.tree['tau']    = self.event.bkkmm_jpsikk_tau[cand]
            self.tree['taue']   = self.event.bkkmm_jpsikk_taue[cand]
            self.tree['tauxy']  = self.event.bkkmm_jpsikk_tauxy[cand]
            self.tree['tauxye'] = self.event.bkkmm_jpsikk_tauxye[cand]
            if hasattr(self.event, 'bkmm_gen_tau'):
                self.tree['gtau'] = self.event.bkkmm_gen_tau[cand]
                self.tree['mc_match'] = self.event.bkkmm_gen_pdgId[cand]

            mu1 = self.event.mm_mu1_index[self.event.bkkmm_mm_index[cand]]
            mu2 = self.event.mm_mu2_index[self.event.bkkmm_mm_index[cand]]
        
            if mu1 >= 0:
                self.tree['m1pt']  = self.event.Muon_pt[mu1]
                self.tree['m1eta'] = self.event.Muon_eta[mu1]
                self.tree['m1phi'] = self.event.Muon_phi[mu1]
                self.tree['m1q']   = self.event.Muon_charge[mu1]
                self.tree['m1bdt'] = self.event.Muon_softMva[mu1]

            if mu2 >= 0:
                self.tree['m2pt']  = self.event.Muon_pt[mu2]
                self.tree['m2eta'] = self.event.Muon_eta[mu2]
                self.tree['m2phi'] = self.event.Muon_phi[mu2]
                self.tree['m2q']   = self.event.Muon_charge[mu2]
                self.tree['m2bdt'] = self.event.Muon_softMva[mu2]

            if mu1 >= 0 and mu2 >= 0:
                self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]

        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])
            
        if abs(self.tree['m1eta']) < 0.7 and abs(self.tree['m2eta']) < 0.7:
            self.tree['chan'] = 0
        else:
            self.tree['chan'] = 1

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
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/515/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/02F7319D-D4D6-8340-B94F-2D882775B406.root",
            ],
        "signal_only" : True,
        "tree_name" : "bsmmMc",
        "blind" : False,
        "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
                "Muon_softMvaId[mm_mu1_index] and "\
                "abs(mm_kin_mu1eta)<1.4 and "\
                "mm_kin_mu1pt>4 and "\
                "Muon_softMvaId[mm_mu2_index] and "\
                "abs(mm_kin_mu2eta)<1.4 and "\
                "mm_kin_mu2pt>4 and "\
                "abs(mm_kin_mass-5.4)<0.5 and "\
                "mm_kin_sl3d>4 and "\
                "mm_kin_vtx_chi2dof<5",
        "final_state" : "mm",
        "best_candidate": "mm_kin_pt",
        # "best_candidate": "",
      }  

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/515/Charmonium+Run2018D-PromptReco-v2+MINIAOD/6E7D2677-4622-F047-86A2-91E4D847CD88.root",
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "bupsikData",
    #     "blind" : False,
    #     "cut" :
    #         "mm_mu1_index[bkmm_mm_index]>=0 and "\
    #         "mm_mu2_index[bkmm_mm_index]>=0 and "\
    #         "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 and "\
    #         "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 and "\
    #         "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 and "\
    #         "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 and "\
    #         "mm_kin_pt[bkmm_mm_index]>7.0 and "\
    #         "mm_kin_alphaBS[bkmm_mm_index]<0.4 and "\
    #         "mm_kin_vtx_prob[bkmm_mm_index]>0.1 and "\
    #         "bkmm_jpsimc_vtx_prob>0.1 and "\
    #         "mm_kin_sl3d[bkmm_mm_index]>4 and "\
    #         "abs(bkmm_jpsimc_mass-5.4)<0.5",
    #     "final_state" : "bkmm",
    #     # "best_candidate": "bkmm_jpsimc_pt",
    #     "best_candidate": "",
    #   }  

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/512/BsToJpsiPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/3BF074E4-960E-4E4D-80FA-D0B67D2EFC7F.root"
    #     ],
    #     "signal_only" : True,
    #     "tree_name" : "bspsiphiMc",
    #     "blind" : False,
    #     "cut" :
    #         "mm_mu1_index[bkkmm_mm_index]>=0 and "\
    #         "mm_mu2_index[bkkmm_mm_index]>=0 and "\
    #         "Muon_softMvaId[mm_mu1_index[bkkmm_mm_index]] and "\
    #         "abs(Muon_eta[mm_mu1_index[bkkmm_mm_index]])<1.4 and "\
    #         "Muon_pt[mm_mu1_index[bkkmm_mm_index]]>4 and "\
    #         "Muon_softMvaId[mm_mu2_index[bkkmm_mm_index]] and "\
    #         "abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4 and "\
    #         "Muon_pt[mm_mu2_index[bkkmm_mm_index]]>4 and "\
    #         "abs(bkkmm_jpsikk_mass-5.3)<0.5 and "\
    #         "bkkmm_jpsikk_sl3d>4 and "\
    #         "bkkmm_jpsikk_vtx_chi2dof<5",
    #     "final_state" : "bkkmm",
    #     # "best_candidate": "bkkmm_jpsikk_pt",
    #     "best_candidate": "",
    #   }  
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    p = FlatNtupleForMLFit(file_name)
    print p.__dict__
        
    p.process()
