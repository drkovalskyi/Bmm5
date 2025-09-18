from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
import multiprocessing
from datetime import datetime
import hashlib
import ROOT

LorentzVector = ROOT.ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')

def is_b_hadron(pdgId):
    # b meson
    if int(abs(pdgId)/100) % 10 == 5:
        return True
    # b meson
    if int(abs(pdgId)/1000) % 10 == 5:
        return True
    return False

class FlatNtupleForMLFit(FlatNtupleBase):
    """Flat ROOT ntuple producer for Bmm5 UML fit"""

    final_states = [
        'mm', 'bkmm', 'bkkmm', 'em', 'hh', 'jpsimm', 'bkstarmm',
    ]
    
    leaf_counts = { 
        'mm':    'nmm',
        'hh':    'nhh',
        'bkmm':  'nbkmm',
        'bkkmm': 'nbkkmm',
        'bkstarmm': 'nbkkmm',
        'em':    'nem',
    }

    triggers_to_store = [
        # Run 2
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
        'HLT_Mu9_IP5_part0',
        'HLT_Mu9_IP5_part1',
        'HLT_Mu9_IP5_part2',
        'HLT_Mu9_IP5_part3',
        'HLT_Mu9_IP5_part4',
        'HLT_Mu12_IP6_part0',
        'HLT_Mu12_IP6_part1',
        'HLT_Mu12_IP6_part2',
        'HLT_Mu12_IP6_part3',
        'HLT_Mu12_IP6_part4',

        # Run 3
        'HLT_Mu12_IP6',
        'HLT_DoubleMu4_3_LowMass',
        'HLT_DoubleMu2_Jpsi_LowPt',
        'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
        'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
        'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
        'L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2'
    ]
    
    mm_extra_floats = ["mm_kin_alpha", "mm_kin_alphaBS", "mm_kin_spvip", "mm_kin_pvip", 
		       "mm_iso", "mm_m1iso", "mm_m2iso", "mm_kin_sl3d", "mm_kin_vtx_chi2dof", 
		       "mm_otherVtxMaxProb1", "mm_otherVtxMaxProb2", "mm_kin_slxy", "mm_kin_lxy"]
    mm_extra_ints = ["mm_nBMTrks",]

    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing information
        if 'best_candidate' not in self.job_info:
            self.job_info['best_candidate'] = ""

        for parameter in ['input', 'blind', 'cut', 'final_state']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)
            
            if self.job_info['final_state'] not in self.final_states:
                raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        # set branches to be kept if pre-skimmed
        if "pre-selection-keep" not in self.job_info:
            if self.job_info['final_state'] == "bkstarmm":
                self.job_info["pre-selection-keep"] = "^(" + \
                    "GenPart_.*|nGenPart|mm_.*|nmm|bkkmm_.*|nbkkmm|" + \
                    "Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + \
                    "HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_DoubleMu2_Jpsi_LowPt|HLT_DoubleMu4_3_LowMass_SS|HLT_Mu0_L1DoubleMu|HLT_ZeroBias|" + \
                    "PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock" + \
                    ")$"
            elif 'pre-selection' in self.job_info:
                raise Exception("'pre-selection-keep' needs to be provided if 'pre-selection' is used to gain any benefits from preselection")

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
            n = getattr(self.event, self.leaf_counts[self.job_info['final_state']])
            for cand in range(n):
                if self.job_info['blind']:
                    if self.job_info['final_state'] == 'mm':
                        if self.event.mm_kin_mass[cand] < 5.50 and \
                           self.event.mm_kin_mass[cand] > 5.15:
                            continue
                    if self.job_info['final_state'] == 'em':
                        if self.event.mm_kin_mass[cand] < 5.70 and \
                           self.event.mm_kin_mass[cand] > 5.10:
                            continue

                # cut = self.job_info['cut'].format(cand=cand, tree="self.event")
                format_dict = { self.leaf_counts[self.job_info['final_state']] : cand,
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
        self.tree.addBranch('run',         'UInt_t', 0)
        self.tree.addBranch('evt',      'ULong64_t', 0)
        self.tree.addBranch('npv',         'UInt_t', 0, "Number of good reconstructed primary vertices")
        self.tree.addBranch('ls',          'UInt_t', 0)
        self.tree.addBranch('n',           'UInt_t', 0, "Number of candidates")
        self.tree.addBranch('npu',         'UInt_t', 0, "number of pileup interactions that have been added to the event in the current bunch crossing")
        self.tree.addBranch('npu_mean',   'Float_t', 0, "tru mean number of pileup interactions")
        self.tree.addBranch('certified_muon',    'Int_t', 0, "Event passed Muon Certification")
        self.tree.addBranch('certified_golden',  'Int_t', 0, "Event passed Golden Certification")

        # ps
        # cw8
        self.tree.addBranch('bdt',        'Float_t', 0, "BDT score")
        self.tree.addBranch('pt',         'Float_t', 0, "pt of the B candidate")
        self.tree.addBranch('eta',        'Float_t', 0, "eta of the B candidate")
        self.tree.addBranch('phi',        'Float_t', 0, "phi of the B candidate")
        self.tree.addBranch('m',          'Float_t', 0, "mass of the B candidate")
        self.tree.addBranch('m2',         'Float_t', 0, "mass of the B candidate (alternative hypothesis)")
        self.tree.addBranch('m_raw',      'Float_t', 0, "raw mass of the B candidate")
        self.tree.addBranch('me',         'Float_t', 0, "mass error of the B candidate")
        self.tree.addBranch('tau',        'Float_t', 0, "decay time of the B candidate")
        self.tree.addBranch('taue',       'Float_t', 0, "decay time error of the B candidate")
        self.tree.addBranch('tauxy',      'Float_t', 0, "2D decay time")
        self.tree.addBranch('tauxye',     'Float_t', 0, "2D decay time error")
        self.tree.addBranch('gtau',       'Float_t', 0, "generated decay time")

        self.tree.addBranch('chan',        'UInt_t', 0, "0 if |eta|<0.7 for both muons. 1 otherwise")

        if self.job_info['final_state'] in ['mm', 'bkmm', 'bkkmm', 'bkstarmm']:
            self.tree.addBranch('muid',        'UInt_t', 0, "MVA muon selection for both muons")

            self.tree.addBranch('m1pt',       'Float_t', 0)
            self.tree.addBranch('m1eta',      'Float_t', 0)
            self.tree.addBranch('m1phi',      'Float_t', 0)
            self.tree.addBranch('m1q',          'Int_t', 0, "Muon charge. 0 for hadrons")
            self.tree.addBranch('m1mc',         'Int_t', 0, "PDG id of the mother of the MC matched muon")
            self.tree.addBranch('m1bdt',      'Float_t', 0, "Muon BDT. 1 for hadrons")

            self.tree.addBranch('m2pt',       'Float_t', 0)
            self.tree.addBranch('m2eta',      'Float_t', 0)
            self.tree.addBranch('m2phi',      'Float_t', 0)
            self.tree.addBranch('m2q',          'Int_t', 0, "Muon charge. 0 for hadrons")
            self.tree.addBranch('m2mc',         'Int_t', 0, "PDG id of the mother of the MC matched muon")
            self.tree.addBranch('m2bdt',      'Float_t', 0, "Muon BDT. 1 for hadrons")
            
            self.tree.addBranch('kaon_pt',       'Float_t', 0)
            self.tree.addBranch('kaon_eta',      'Float_t', 0)
            self.tree.addBranch('kaon_phi',      'Float_t', 0)
            self.tree.addBranch('kaon_sdxy_bs',  'Float_t', -1)
            self.tree.addBranch('kaon_mc',         'Int_t', 0, "PDG id of the MC gen particle matched to the kaon candidate")
            self.tree.addBranch('decay_signature', 'Long64_t', 0, "Decay signature: product of pdgId of parent daughters")
            self.tree.addBranch('decay_parent',    'Int_t', 0, "Decay parent: mother for mm and grandma for mmX")
            self.tree.addBranch('decay_ndau',      'Int_t', 0, "Number of daughters")
            
            if self.job_info['final_state'] == "mm" and \
               "mm_extra_info" in self.job_info and \
               self.job_info["mm_extra_info"] == True:
                for var in self.mm_extra_floats:
                    self.tree.addBranch(var, "Float_t", 0)
                for var in self.mm_extra_ints:
                    self.tree.addBranch(var, "Int_t", 0)
            
            if self.job_info['final_state'] == "bkstarmm":
                self.tree.addBranch('pion_pt',       'Float_t', 0)
                self.tree.addBranch('pion_eta',      'Float_t', 0)
                self.tree.addBranch('pion_phi',      'Float_t', 0)
                self.tree.addBranch('pion_sdxy_bs',  'Float_t', -1)
                self.tree.addBranch('pion_mc',         'Int_t', 0, "PDG id of the MC gen particle matched to the pion candidate")
                self.tree.addBranch('kstar_mass',    'Float_t', 0, "kpi pair with the closest mass to Kstar")
                self.tree.addBranch('kstar2_mass',   'Float_t', 0, "Alternative kpi hypothese assignment")
                self.tree.addBranch('phi_veto',       'UInt_t', 0)

            if self.job_info['final_state'] == "bkkmm":
                self.tree.addBranch('kaon2_pt',       'Float_t', 0)
                self.tree.addBranch('kaon2_eta',      'Float_t', 0)
                self.tree.addBranch('kaon2_phi',      'Float_t', 0)
                self.tree.addBranch('kaon2_sdxy_bs',  'Float_t', -1)
                self.tree.addBranch('kaon2_mc',         'Int_t', 0, "PDG id of the MC gen particle matched to the kaon candidate")
                self.tree.addBranch('kstar_veto',      'UInt_t', 0)

            if self.job_info['final_state'] in ["bkstarmm", "bkkmm"]:
                self.tree.addBranch('kk_mass',    'Float_t', 0, "Mass of two tracks with kaon-kaon hyposthesis")
                self.tree.addBranch('kpi_mass',   'Float_t', 0, "Mass of two tracks with kaon-pion hyposthesis")
                self.tree.addBranch('pik_mass',   'Float_t', 0, "Mass of two tracks with pion-kaon hyposthesis")
                self.tree.addBranch('pipi_mass',  'Float_t', 0, "Mass of two tracks with pion-pion hyposthesis")

        elif self.job_info['final_state'] == 'em':
            self.tree.addBranch('id',          'UInt_t', 0, "Tight electron and muon selections")

            self.tree.addBranch('elpt',       'Float_t', 0)
            self.tree.addBranch('eleta',      'Float_t', 0)
            self.tree.addBranch('elphi',      'Float_t', 0)
            self.tree.addBranch('elq',          'Int_t', 0, "Electron charge. 0 for hadrons")
            self.tree.addBranch('elmc',         'Int_t', 0, "PDG id of the mother of the MC matched muon")

            self.tree.addBranch('mupt',       'Float_t', 0)
            self.tree.addBranch('mueta',      'Float_t', 0)
            self.tree.addBranch('muphi',      'Float_t', 0)
            self.tree.addBranch('muq',          'Int_t', 0, "Muon charge. 0 for hadrons")
            self.tree.addBranch('mumc',         'Int_t', 0, "PDG id of the mother of the MC matched muon")
        elif self.job_info['final_state'] == 'hh':
            self.tree.addBranch('h1pt',       'Float_t', 0)
            self.tree.addBranch('h1eta',      'Float_t', 0)
            self.tree.addBranch('h1phi',      'Float_t', 0)
            self.tree.addBranch('h1q',          'Int_t', 0, "Hadron charge")

            self.tree.addBranch('h2pt',       'Float_t', 0)
            self.tree.addBranch('h2eta',      'Float_t', 0)
            self.tree.addBranch('h2phi',      'Float_t', 0)
            self.tree.addBranch('h2q',          'Int_t', 0, "Hadron charge")
        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        self.tree.addBranch('mc_match',     'Int_t', 0, "PdgId of the MC matched B meson")
        self.tree.addBranch('mc_bhh',       'Int_t', 0, "PdgId of B meson in Btohh MC events")

        for trigger in self.triggers_to_store:
            self.tree.addBranch(trigger, 'Int_t', -1, "Trigger decision: 1 - fired, 0 - didn't fire, -1 - no information")
            self.tree.addBranch("%s_ps" % trigger, 'Float_t', 999999, "Prescale. 0 - Off, 999999 - no information")
            self.tree.addBranch("%s_matched" % trigger, 'Int_t', 0,  "matched to the trigger objets")

    def _tag_bhh(self):
        """Determine B meson pdgId for Btohh mixed samples"""
        if not hasattr(self.event, "GenPart_pdgId"): return 0
        for i, id in enumerate(self.event.GenPart_pdgId):
            if abs(id) != 13: continue

            mother_idx = self.event.GenPart_genPartIdxMother[i]
            if mother_idx < 0: continue
            if abs(self.event.GenPart_pdgId[mother_idx]) not in (211, 321): continue

            grandmother_idx = self.event.GenPart_genPartIdxMother[mother_idx]
            if grandmother_idx < 0: continue

            grandmother_pdgId = self.event.GenPart_pdgId[grandmother_idx]
            if abs(grandmother_pdgId) not in (511, 531): continue
            return grandmother_pdgId
        return 0

    def _compute_mass_hypotheses_from_kk(self, cand):
        kaon1_p4 = LorentzVector(self.event.bkkmm_kaon1_pt[cand],
                                 self.event.bkkmm_kaon1_eta[cand],
                                 self.event.bkkmm_kaon1_phi[cand], 0.497648)
        pion1_p4 = LorentzVector(self.event.bkkmm_kaon1_pt[cand],
                                 self.event.bkkmm_kaon1_eta[cand],
                                 self.event.bkkmm_kaon1_phi[cand], 0.139570)
        
        kaon2_p4 = LorentzVector(self.event.bkkmm_kaon2_pt[cand],
                                 self.event.bkkmm_kaon2_eta[cand],
                                 self.event.bkkmm_kaon2_phi[cand], 0.497648)
        pion2_p4 = LorentzVector(self.event.bkkmm_kaon2_pt[cand],
                                 self.event.bkkmm_kaon2_eta[cand],
                                 self.event.bkkmm_kaon2_phi[cand], 0.139570)
        
        b_p4     = LorentzVector(self.event.bkkmm_jpsikk_pt[cand],
                                 self.event.bkkmm_jpsikk_eta[cand],
                                 self.event.bkkmm_jpsikk_phi[cand],
                                 self.event.bkkmm_jpsikk_mass[cand])
            
        kpi_mass  = (kaon1_p4 + pion2_p4).mass()
        pik_mass  = (kaon2_p4 + pion1_p4).mass()
        pipi_mass = (pion1_p4 + pion2_p4).mass()

        b_kpi_mass = (b_p4 - kaon2_p4 + pion2_p4).mass()
        b_pik_mass = (b_p4 - kaon1_p4 + pion1_p4).mass()
        b_pipi_mass = (b_p4 - kaon1_p4 - kaon2_p4 + pion1_p4 + pion2_p4).mass()

        return (kpi_mass, pik_mass, pipi_mass, b_kpi_mass, b_pik_mass, b_pipi_mass)

    
    def _fill_tree(self, cand, ncands):
        self.tree.reset()
        save_cand = True
        kstar_pdg_mass = 0.89167
        
        ## event info
        self.tree['run'] = self.event.run
        self.tree['ls']  = self.event.luminosityBlock
        self.tree['evt'] = self.event.event
        self.tree['npv'] = ord(self.event.PV_npvsGood) if isinstance(self.event.PV_npvsGood, str) else self.event.PV_npvsGood

        self.tree['certified_muon']   = self._is_certified_event(self.event, "muon")
        self.tree['certified_golden'] = self._is_certified_event(self.event, "golden")
        self.tree['n']   = ncands

        ## MC info
        mc = False
        if hasattr(self.event, 'Pileup_nTrueInt'):
            mc = True
            self.tree['npu']      = self.event.Pileup_nPU
            self.tree['npu_mean'] = self.event.Pileup_nTrueInt
            self.tree['mc_bhh']   = self._tag_bhh()

            # b-hadron decay signature
            if self.job_info['final_state'] in ['mm', 'bkmm', 'bkkmm', 'bkstarmm'] and hasattr(self.event, 'nGenPart'):
                mm_index = None
                if self.job_info['final_state'] == 'mm':
                    mm_index = cand
                elif self.job_info['final_state'] == 'bkmm':
                    mm_index = self.event.bkmm_mm_index[cand]
                elif self.job_info['final_state'] in ['bkkmm', 'bkstarmm']:
                    mm_index = self.event.bkkmm_mm_index[cand]
                    
                # check if we have a common ancestor for the dimuon 
                if mm_index != None :
                    icgen = self.event.mm_gen_cindex[mm_index]
                    # print(f"icgen: {icgen}")
                    if icgen >= 0:
                        decay_index = None
                        if self.job_info['final_state'] == "mm":
                            decay_index = icgen
                        else:
                            # find grand mother
                            # print(f"GenPart_pdgId[icgen]: {self.event.GenPart_pdgId[icgen]}")
                            # print(f"GenPart_genPartIdxMother[icgen]: {self.event.GenPart_genPartIdxMother[icgen]}")
                            if self.event.GenPart_genPartIdxMother[icgen] >= 0:
                                decay_index = self.event.GenPart_genPartIdxMother[icgen]
                                # print(f"GenPart_pdgId[decay_index]: {self.event.GenPart_pdgId[decay_index]}")
                        if decay_index != None:
                            decay_parent = self.event.GenPart_pdgId[decay_index]
                            decay_signature = 1
                            decay_ndau = 0
                            for igen in range(self.event.nGenPart):
                                if self.event.GenPart_genPartIdxMother[igen] == decay_index:
                                    # print(f"\tdaughter: {self.event.GenPart_pdgId[igen]}")
                                    decay_signature *= self.event.GenPart_pdgId[igen]
                                    decay_ndau += 1 
                            self.tree['decay_signature'] = decay_signature 
                            self.tree['decay_parent'] = decay_parent
                            self.tree['decay_ndau'] = decay_ndau

                    
        if self.job_info['final_state'] == 'mm':
            # B to mm
            self.tree['bdt']    = self.event.mm_mva[cand]
            self.tree['pt']     = self.event.mm_kin_pt[cand]
            self.tree['eta']    = self.event.mm_kin_eta[cand]
            self.tree['phi']    = self.event.mm_kin_phi[cand]
            self.tree['m']      = self.event.mm_kin_mass[cand]
            self.tree['m_raw']  = self.event.mm_mass[cand]
            self.tree['me']     = self.event.mm_kin_massErr[cand]
            self.tree['tau']    = self.event.mm_kin_tau[cand]
            self.tree['taue']   = self.event.mm_kin_taue[cand]
            self.tree['tauxy']  = self.event.mm_kin_tauxy[cand]
            self.tree['tauxye'] = self.event.mm_kin_tauxye[cand]
            if mc:
                self.tree['gtau'] = self.event.mm_gen_tau[cand]
                self.tree['mc_match'] = self.event.mm_gen_pdgId[cand]

            if "mm_extra_info" in self.job_info and \
               self.job_info["mm_extra_info"] == True:
                for var in self.mm_extra_floats:
                    self.tree[var] = getattr(self.event, var)[cand]
                for var in self.mm_extra_ints:
                    self.tree[var] = getattr(self.event, var)[cand]

            mu1 = self.event.mm_mu1_index[cand]
            mu2 = self.event.mm_mu2_index[cand]

            self.tree['m1pt']  = self.event.mm_mu1_pt[cand]
            self.tree['m1eta'] = self.event.mm_mu1_eta[cand]
            self.tree['m1phi'] = self.event.mm_mu1_phi[cand]
            self.tree['m2pt']  = self.event.mm_mu2_pt[cand]
            self.tree['m2eta'] = self.event.mm_mu2_eta[cand]
            self.tree['m2phi'] = self.event.mm_mu2_phi[cand]

            if mu1 >= 0:
                self.tree['m1q']   = self.event.Muon_charge[mu1]
                if mc:
                    self.tree['m1mc']  = self.event.mm_gen_mu1_mpdgId[cand]
                self.tree['m1bdt'] = self.event.Muon_softMva[mu1]
            else: 
                self.tree['m1bdt'] = 1

            if mu2 >= 0:
                self.tree['m2q']   = self.event.Muon_charge[mu2]
                if mc:
                    self.tree['m2mc']  = self.event.mm_gen_mu2_mpdgId[cand]
                self.tree['m2bdt'] = self.event.Muon_softMva[mu2]
            else: 
                self.tree['m2bdt'] = 1

            if mu1 >= 0 and mu2 >= 0:
                self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]

        elif self.job_info['final_state'] == 'hh':
            # B to hh
            self.tree['bdt']    = self.event.hh_mva[cand]
            self.tree['pt']     = self.event.hh_kin_pt[cand]
            self.tree['eta']    = self.event.hh_kin_eta[cand]
            self.tree['phi']    = self.event.hh_kin_phi[cand]
            self.tree['m']      = self.event.hh_kin_mass[cand]
            self.tree['me']     = self.event.hh_kin_massErr[cand]
            self.tree['tau']    = self.event.hh_kin_tau[cand]
            self.tree['taue']   = self.event.hh_kin_taue[cand]
            self.tree['tauxy']  = self.event.hh_kin_tauxy[cand]
            self.tree['tauxye'] = self.event.hh_kin_tauxye[cand]
            if mc:
                self.tree['gtau'] = self.event.hh_gen_tau[cand]
                self.tree['mc_match'] = self.event.hh_gen_pdgId[cand]

            self.tree['h1pt']  = self.event.hh_had1_pt[cand]
            self.tree['h1eta'] = self.event.hh_had1_eta[cand]
            self.tree['h1phi'] = self.event.hh_had1_phi[cand]
            self.tree['h1q']   = 1 if self.event.hh_had1_pdgId[cand] > 0 else 0

            self.tree['h2pt']  = self.event.hh_had2_pt[cand]
            self.tree['h2eta'] = self.event.hh_had2_eta[cand]
            self.tree['h2phi'] = self.event.hh_had2_phi[cand]
            self.tree['h2q']   = 1 if self.event.hh_had2_pdgId[cand] > 0 else 0

        elif self.job_info['final_state'] == 'bkmm':
            # B to Jpsi K
            self.tree['pt']     = self.event.bkmm_jpsimc_pt[cand]
            self.tree['eta']    = self.event.bkmm_jpsimc_eta[cand]
            self.tree['phi']    = self.event.bkmm_jpsimc_phi[cand]
            self.tree['m']      = self.event.bkmm_jpsimc_mass[cand]
            self.tree['m_raw']  = self.event.bkmm_nomc_mass[cand]
            self.tree['me']     = self.event.bkmm_jpsimc_massErr[cand]
            self.tree['tau']    = self.event.bkmm_jpsimc_tau[cand]
            self.tree['taue']   = self.event.bkmm_jpsimc_taue[cand]
            self.tree['tauxy']  = self.event.bkmm_jpsimc_tauxy[cand]
            self.tree['tauxye'] = self.event.bkmm_jpsimc_tauxye[cand]
            self.tree['kaon_pt']     = self.event.bkmm_kaon_pt[cand]
            self.tree['kaon_eta']    = self.event.bkmm_kaon_eta[cand]
            self.tree['kaon_phi']    = self.event.bkmm_kaon_phi[cand]
            self.tree['kaon_sdxy_bs'] = self.event.bkmm_kaon1_sdxy_bs[cand]
            if mc:
                self.tree['gtau'] = self.event.bkmm_gen_tau[cand]
                self.tree['mc_match'] = self.event.bkmm_gen_pdgId[cand]

            mm_index = self.event.bkmm_mm_index[cand]
            mu1 = self.event.mm_mu1_index[mm_index]
            mu2 = self.event.mm_mu2_index[mm_index]

            if mu1 >= 0:
                self.tree['m1pt']  = self.event.mm_mu1_pt[mm_index]
                self.tree['m1eta'] = self.event.mm_mu1_eta[mm_index]
                self.tree['m1phi'] = self.event.mm_mu1_phi[mm_index]
                try:
                    self.tree['m1q']   = self.event.Muon_charge[mu1]
                    self.tree['m1bdt'] = self.event.Muon_softMva[mu1]
                except:
                    pass
                if mc:
                    self.tree['m1mc']  = self.event.mm_gen_mu1_mpdgId[mm_index]

            if mu2 >= 0:
                self.tree['m2pt']  = self.event.mm_mu2_pt[mm_index]
                self.tree['m2eta'] = self.event.mm_mu2_eta[mm_index]
                self.tree['m2phi'] = self.event.mm_mu2_phi[mm_index]
                try:
                    self.tree['m2q']   = self.event.Muon_charge[mu2]
                    self.tree['m2bdt'] = self.event.Muon_softMva[mu2]
                except:
                    pass
                if mc:
                    self.tree['m2mc']  = self.event.mm_gen_mu2_mpdgId[mm_index]

            if mu1 >= 0 and mu2 >= 0:
                try:
                    self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
                except:
                    pass

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
            if mc:
                self.tree['gtau'] = self.event.bkkmm_gen_tau[cand]
                self.tree['mc_match'] = self.event.bkkmm_gen_pdgId[cand]

            mm_index = self.event.bkkmm_mm_index[cand]
            mu1 = self.event.mm_mu1_index[mm_index]
            mu2 = self.event.mm_mu2_index[mm_index]
            
            self.tree['m1pt']  = self.event.mm_mu1_pt[mm_index]
            self.tree['m1eta'] = self.event.mm_mu1_eta[mm_index]
            self.tree['m1phi'] = self.event.mm_mu1_phi[mm_index]
            self.tree['m2pt']  = self.event.mm_mu2_pt[mm_index]
            self.tree['m2eta'] = self.event.mm_mu2_eta[mm_index]
            self.tree['m2phi'] = self.event.mm_mu2_phi[mm_index]

            (kpi_mass, pik_mass, pipi_mass, b_kpi_mass, b_pik_mass, b_pipi_mass) = self._compute_mass_hypotheses_from_kk(cand)
            self.tree['kk_mass'] = self.event.bkkmm_kk_mass[cand]
            self.tree['kpi_mass'] = kpi_mass 
            self.tree['pik_mass'] = pik_mass 
            self.tree['pipi_mass'] = pipi_mass 
            
            if abs(kpi_mass - kstar_pdg_mass) < abs(pik_mass - kstar_pdg_mass):
                kstar_mass  = kpi_mass
            else:
                kstar_mass  = pik_mass
            if kstar_mass > 0.846 and kstar_mass < 0.946:
                self.tree['kstar_veto'] = 1
                
            if 'veto_kstar' in self.job_info and self.job_info['veto_kstar'] == True:
                if self.tree['kstar_veto']:
                    save_cand = False
            
            self.tree['kaon_sdxy_bs'] = self.event.bkkmm_kaon1_sdxy_bs[cand]
            self.tree['kaon2_sdxy_bs'] = self.event.bkkmm_kaon2_sdxy_bs[cand]
            
            try:
                if mu1 >= 0:
                    self.tree['m1q']   = self.event.Muon_charge[mu1]
                    if mc:
                        self.tree['m1mc']  = self.event.mm_gen_mu1_mpdgId[mm_index]
                    self.tree['m1bdt'] = self.event.Muon_softMva[mu1]

                if mu2 >= 0:
                    self.tree['m2q']   = self.event.Muon_charge[mu2]
                    if mc:
                        self.tree['m2mc']  = self.event.mm_gen_mu2_mpdgId[mm_index]
                    self.tree['m2bdt'] = self.event.Muon_softMva[mu2]

                if mu1 >= 0 and mu2 >= 0:
                    self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
            except:
                pass

        elif self.job_info['final_state'] == 'bkstarmm':
            # B to Jpsi Kstar 
            self.tree['pt']     = self.event.bkkmm_jpsikk_pt[cand]
            self.tree['eta']    = self.event.bkkmm_jpsikk_eta[cand]
            self.tree['phi']    = self.event.bkkmm_jpsikk_phi[cand]

            (kpi_mass, pik_mass, pipi_mass, b_kpi_mass, b_pik_mass, b_pipi_mass) = self._compute_mass_hypotheses_from_kk(cand)
            self.tree['kk_mass'] = self.event.bkkmm_kk_mass[cand]
            self.tree['pipi_mass'] = pipi_mass
            
            if abs(kpi_mass - kstar_pdg_mass) < abs(pik_mass - kstar_pdg_mass):
                self.tree['m']  = b_kpi_mass
                self.tree['m2'] = b_pik_mass
                self.tree['kstar_mass']  = kpi_mass
                self.tree['kstar2_mass'] = pik_mass
                self.tree['kaon_pt'] = self.event.bkkmm_kaon1_pt[cand]
                self.tree['kaon_eta'] = self.event.bkkmm_kaon1_pt[cand]
                self.tree['kaon_phi'] = self.event.bkkmm_kaon1_pt[cand]
                self.tree['kaon_sdxy_bs'] = self.event.bkkmm_kaon1_sdxy_bs[cand]
                self.tree['pion_pt'] = self.event.bkkmm_kaon2_pt[cand]
                self.tree['pion_eta'] = self.event.bkkmm_kaon2_pt[cand]
                self.tree['pion_phi'] = self.event.bkkmm_kaon2_pt[cand]
                self.tree['pion_sdxy_bs'] = self.event.bkkmm_kaon2_sdxy_bs[cand]
                if mc:
                    self.tree['kaon_mc'] = self.event.bkkmm_gen_kaon1_pdgId[cand]
                    self.tree['pion_mc'] = self.event.bkkmm_gen_kaon2_pdgId[cand]
            else:
                self.tree['m']  = b_pik_mass
                self.tree['m2'] = b_kpi_mass
                self.tree['kstar_mass']  = pik_mass
                self.tree['kstar2_mass'] = kpi_mass
                self.tree['kaon_pt'] = self.event.bkkmm_kaon2_pt[cand]
                self.tree['kaon_eta'] = self.event.bkkmm_kaon2_pt[cand]
                self.tree['kaon_phi'] = self.event.bkkmm_kaon2_pt[cand]
                self.tree['kaon_sdxy_bs'] = self.event.bkkmm_kaon2_sdxy_bs[cand]
                self.tree['pion_pt'] = self.event.bkkmm_kaon1_pt[cand]
                self.tree['pion_eta'] = self.event.bkkmm_kaon1_pt[cand]
                self.tree['pion_phi'] = self.event.bkkmm_kaon1_pt[cand]
                self.tree['pion_sdxy_bs'] = self.event.bkkmm_kaon1_sdxy_bs[cand]
                if mc:
                    self.tree['kaon_mc'] = self.event.bkkmm_gen_kaon2_pdgId[cand]
                    self.tree['pion_mc'] = self.event.bkkmm_gen_kaon1_pdgId[cand]
            
            self.tree['me']     = self.event.bkkmm_jpsikk_massErr[cand]

            # apply cuts
            if (self.tree['kstar_mass'] < 0.846 or self.tree['kstar_mass'] > 0.946) or \
               (self.tree['m'] < 4.9 or self.tree['m'] > 5.9):
                save_cand = False
                
            if abs(self.event.bkkmm_kk_mass[cand] - 1.02) < 0.01:
                self.tree['phi_veto'] = 1
                
            if 'veto_phi' in self.job_info and self.job_info['veto_phi'] == True:
                if self.tree['phi_veto']:
                    save_cand = False
            
            mass_scale = self.tree['m'] / self.event.bkkmm_jpsikk_mass[cand]
            
            self.tree['tau']    = self.event.bkkmm_jpsikk_tau[cand] * mass_scale
            self.tree['taue']   = self.event.bkkmm_jpsikk_taue[cand] * mass_scale
            self.tree['tauxy']  = self.event.bkkmm_jpsikk_tauxy[cand] * mass_scale
            self.tree['tauxye'] = self.event.bkkmm_jpsikk_tauxye[cand] * mass_scale
            if mc:
                self.tree['gtau'] = self.event.bkkmm_gen_tau[cand]
                self.tree['mc_match'] = self.event.bkkmm_gen_pdgId[cand]

            mm_index = self.event.bkkmm_mm_index[cand]
            mu1 = self.event.mm_mu1_index[mm_index]
            mu2 = self.event.mm_mu2_index[mm_index]
            
            self.tree['m1pt']  = self.event.mm_mu1_pt[mm_index]
            self.tree['m1eta'] = self.event.mm_mu1_eta[mm_index]
            self.tree['m1phi'] = self.event.mm_mu1_phi[mm_index]
            self.tree['m2pt']  = self.event.mm_mu2_pt[mm_index]
            self.tree['m2eta'] = self.event.mm_mu2_eta[mm_index]
            self.tree['m2phi'] = self.event.mm_mu2_phi[mm_index]

            try:
                if mu1 >= 0:
                    self.tree['m1q']   = self.event.Muon_charge[mu1]
                    if mc:
                        self.tree['m1mc']  = self.event.mm_gen_mu1_mpdgId[mm_index]
                    self.tree['m1bdt'] = self.event.Muon_softMva[mu1]

                if mu2 >= 0:
                    self.tree['m2q']   = self.event.Muon_charge[mu2]
                    if mc:
                        self.tree['m2mc']  = self.event.mm_gen_mu2_mpdgId[mm_index]
                    self.tree['m2bdt'] = self.event.Muon_softMva[mu2]

                if mu1 >= 0 and mu2 >= 0:
                    self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
            except:
                pass

        elif self.job_info['final_state'] == 'em':
            # B to em
            self.tree['bdt']    = self.event.em_mva[cand]
            self.tree['pt']     = self.event.em_kin_pt[cand]
            self.tree['eta']    = self.event.em_kin_eta[cand]
            self.tree['phi']    = self.event.em_kin_phi[cand]
            self.tree['m']      = self.event.em_kin_mass[cand]
            self.tree['me']     = self.event.em_kin_massErr[cand]
            self.tree['tau']    = self.event.em_kin_tau[cand]
            self.tree['taue']   = self.event.em_kin_taue[cand]
            self.tree['tauxy']  = self.event.em_kin_tauxy[cand]
            self.tree['tauxye'] = self.event.em_kin_tauxye[cand]
            if mc:
                self.tree['gtau'] = self.event.em_gen_tau[cand]
                self.tree['mc_match'] = self.event.em_gen_pdgId[cand]

            el = self.event.em_el_index[cand]
            mu = self.event.em_mu_index[cand]

            if el >= 0:
                self.tree['elpt']  = self.event.Electron_pt[el]
                self.tree['eleta'] = self.event.Electron_eta[el]
                self.tree['elphi'] = self.event.Electron_phi[el]
                self.tree['elq']   = self.event.Electron_charge[el]
                if mc:
                    self.tree['elmc']  = self.event.em_gen_el_mpdgId[cand]
            else: 
                self.tree['elpt']  = self.event.em_el_pt[cand]
                self.tree['eleta'] = self.event.em_el_eta[cand]
                self.tree['elphi'] = self.event.em_el_phi[cand]

            if mu >= 0:
                self.tree['mupt']  = self.event.Muon_pt[mu]
                self.tree['mueta'] = self.event.Muon_eta[mu]
                self.tree['muphi'] = self.event.Muon_phi[mu]
                self.tree['muq']   = self.event.Muon_charge[mu]
                if mc:
                    self.tree['mumc']  = self.event.em_gen_mu_mpdgId[cand]
            else: 
                self.tree['mupt']  = self.event.em_mu_pt[cand]
                self.tree['mueta'] = self.event.em_mu_eta[cand]
                self.tree['muphi'] = self.event.em_mu_phi[cand]

            if el >= 0 and mu >= 0:
                self.tree['id']   = self.event.Electron_mvaNoIso_WP90[el] and self.event.Muon_softMvaId[mu]
        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        if self.job_info['final_state'] in ['mm', 'bkmm', 'bkkmm', 'bkstarmm']:
            if abs(self.tree['m1eta']) < 0.7 and abs(self.tree['m2eta']) < 0.7:
                self.tree['chan'] = 0
            else:
                self.tree['chan'] = 1

            for trigger in self.triggers_to_store:
                if hasattr(self.event, trigger):
                    self.tree[trigger] = getattr(self.event, trigger)
                if hasattr(self.event, "prescale_" + trigger):
                    self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
                if hasattr(self.event, "MuonId_" + trigger):
                    if mu1 >= 0 and mu2 >= 0:
                        self.tree[trigger + "_matched"] = getattr(self.event, "MuonId_" + trigger)[mu1] and getattr(self.event, "MuonId_" + trigger)[mu2] 
        elif self.job_info['final_state'] == 'em':
            if abs(self.tree['eleta']) < 0.7 and abs(self.tree['mueta']) < 0.7:
                self.tree['chan'] = 0
            else:
                self.tree['chan'] = 1

            for trigger in self.triggers_to_store:
                if hasattr(self.event, trigger):
                    self.tree[trigger] = getattr(self.event, trigger)
                if hasattr(self.event, "prescale_" + trigger):
                    self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
        elif self.job_info['final_state'] == 'hh':
            for trigger in self.triggers_to_store:
                if hasattr(self.event, trigger):
                    self.tree[trigger] = getattr(self.event, trigger)
                if hasattr(self.event, "prescale_" + trigger):
                    self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        if save_cand:
            self.tree.fill()

if __name__ == "__main__":

    ### create a test job

    common_branches = "PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock"

    cuts = dict()
    
    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/517/BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/06387049-7142-D749-A447-C97E700210BF.root",
    #         ],
    #     "signal_only" : True,
    #     "tree_name" : "bsmmMc",
    #     "blind" : False,
    #     "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    #             "Muon_softMvaId[mm_mu1_index] and "\
    #             "abs(mm_kin_mu1eta)<1.4 and "\
    #             "mm_kin_mu1pt>4 and "\
    #             "Muon_softMvaId[mm_mu2_index] and "\
    #             "abs(mm_kin_mu2eta)<1.4 and "\
    #             "mm_kin_mu2pt>4 and "\
    #             "abs(mm_kin_mass-5.4)<0.5 and "\
    #             "mm_kin_sl3d>4 and "\
    #             "mm_kin_vtx_chi2dof<5",
    #     "final_state" : "mm",
    #     "best_candidate": "mm_kin_pt",
    #     # "best_candidate": "",
    #   }

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/517/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v3+MINIAODSIM/04C35C04-7851-6A4C-A80A-6FC132E25025.root",
    #         ],
    #     "signal_only" : False,
    #     "tree_name" : "btohhMcBg",
    #     "blind" : False,
    #     "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    #             "Muon_softMvaId[mm_mu1_index] and "\
    #             "abs(mm_kin_mu1eta)<1.4 and "\
    #             "mm_kin_mu1pt>4 and "\
    #             "Muon_softMvaId[mm_mu2_index] and "\
    #             "abs(mm_kin_mu2eta)<1.4 and "\
    #             "mm_kin_mu2pt>4 and "\
    #             "abs(mm_kin_mass-5.4)<0.5 and "\
    #             "mm_kin_sl3d>4 and "\
    #             "mm_kin_vtx_chi2dof<5",
    #     "final_state" : "mm",
    #     "best_candidate": "mm_kin_pt",
    #     # "best_candidate": "",
    #   }  

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/954c0f64-a35a-4eff-9676-15cc8f5d0bb1.root"
    #         ],
    #     "signal_only" : False,
    #     "tree_name" : "bmmData",
    #     "blind" : False,
    #     "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    #             "Muon_charge[mm_mu1_index] * Muon_charge[mm_mu2_index] < 0 and "\
    #             "Muon_softMva[mm_mu1_index] > 0.45 and "\
    #             "Muon_softMva[mm_mu2_index] > 0.45 and "\
    #             "abs(mm_kin_mass-5.4)<0.5 and "\
    #             "mm_kin_sl3d>6 and "\
    #             "mm_kin_vtx_prob>0.025 and "\
    #             "HLT_DoubleMu4_3_Bs"
    #             ,
    #     "final_state" : "mm",
    #     # "best_candidate": "mm_kin_pt",
    #     "best_candidate": "",
    #   }  

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/516/BTohh_hToMuNu_BsBdMixture_modHadLifetime_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v3+MINIAODSIM/04C35C04-7851-6A4C-A80A-6FC132E25025.root",
    #         ],
    #     "signal_only" : False,
    #     "tree_name" : "btohhMcBg",
    #     "blind" : False,
    #     "cut" : "mm_mu1_index>=0 and mm_mu2_index>=0 and "\
    #             "abs(mm_kin_mu1eta)<1.4 and "\
    #             "mm_kin_mu1pt>4 and "\
    #             "abs(mm_kin_mu2eta)<1.4 and "\
    #             "mm_kin_mu2pt>4 and "\
    #             "abs(mm_kin_mass-5.4)<0.5 and "\
    #             "mm_kin_sl3d>4 and "\
    #             "mm_kin_vtx_prob>0.025",
    #     "final_state" : "mm",
    #     "best_candidate": "mm_kin_pt",
    #     # "best_candidate": "",
    #   }  

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/4e83dd1a-a5db-42a9-a30a-2f618d520065.root",
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "bupsikData",
    #     "blind" : False,
    #     "cut" : (
    #         "mm_mu1_index[bkmm_mm_index]>=0 and "
    #         "mm_mu2_index[bkmm_mm_index]>=0 and "
    #         "Muon_charge[mm_mu1_index[bkmm_mm_index]] * Muon_charge[mm_mu2_index[bkmm_mm_index]] < 0 and "
    #         "Muon_mediumId[mm_mu1_index[bkmm_mm_index]] and "
    #         "Muon_mediumId[mm_mu2_index[bkmm_mm_index]] and "
    #         "Muon_isGlobal[mm_mu1_index[bkmm_mm_index]] and "
    #         "Muon_isGlobal[mm_mu2_index[bkmm_mm_index]] and "
    #         "mm_kin_vtx_prob[bkmm_mm_index]>0.01 and "
    #         "bkmm_jpsimc_vtx_prob>0.025 and "
    #         "bkmm_jpsimc_sl3d>3 and "
    #         "abs(bkmm_jpsimc_alpha) < 0.1 and "
    #         "abs(bkmm_jpsimc_mass-5.4)<0.5"
    #     ),
    #     "final_state" : "bkmm",
    #     "best_candidate": "",
    #     "pre-selection":"bkmm_jpsimc_sl3d>3",
    #     "pre-selection-keep":"^(bkmm_.*|nbkmm|mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|HLT_DoubleMu2_Jpsi_LowPt|L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5|L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4|L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4|L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2|HLT_DoubleMu4_3_LowMass_SS|" + common_branches + ")$",
    #   }  

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/515/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0440EA45-411F-E04F-92A6-B3DDA0AA9DDD.root",
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "bupsikMC",
    #     "blind" : False,
    #     "cut" :
        
    #     "mm_mu1_index[bkmm_mm_index]>=0 and "\
    #     "mm_mu2_index[bkmm_mm_index]>=0 and "\
    #     "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 and "\
    #     "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 and "\
    #     "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 and "\
    #     "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4 and "\
    #     "mm_kin_pt[bkmm_mm_index]>7.0 and "\
    #     "mm_kin_alphaBS[bkmm_mm_index]<0.4 and "\
    #     "mm_kin_vtx_prob[bkmm_mm_index]>0.1 and "\
    #     "bkmm_jpsimc_vtx_prob>0.1 and "\
    #     "mm_kin_sl3d[bkmm_mm_index]>4 and "\
    #     "abs(bkmm_jpsimc_mass-5.4)<0.5",

    #     "final_state" : "bkmm",
    #     "best_candidate": "bkmm_jpsimc_pt",
    #     # "best_candidate": "",
    #   }  
    
    cuts["fit-bkkmm-all"] = \
        "mm_mu1_index[bkkmm_mm_index]>=0 and " \
        "mm_mu2_index[bkkmm_mm_index]>=0 and " \
        "Muon_charge[mm_mu1_index[bkkmm_mm_index]] * Muon_charge[mm_mu2_index[bkkmm_mm_index]] < 0 and " \
        "Muon_mediumId[mm_mu1_index[bkkmm_mm_index]] and " \
        "Muon_mediumId[mm_mu2_index[bkkmm_mm_index]] and " \
        "Muon_isGlobal[mm_mu1_index[bkkmm_mm_index]] and " \
        "Muon_isGlobal[mm_mu2_index[bkkmm_mm_index]] and " \
        "mm_kin_vtx_prob[bkkmm_mm_index]>0.01 and " \
        "bkkmm_jpsikk_vtx_prob>0.025 and " \
        "bkkmm_jpsikk_sl3d>3 and " \
        "abs(bkkmm_jpsikk_alpha) < 0.1 and " \
        "abs(bkkmm_jpsikk_mass-5.4)<0.5 and " \
        "abs(bkkmm_kk_mass-1.02)<0.01"

    cuts["fit-bkkmm"] = cuts["fit-bkkmm-all"] + " and " + \
        "Muon_pt[mm_mu1_index[bkkmm_mm_index]] > 4 and " + \
        "Muon_pt[mm_mu2_index[bkkmm_mm_index]] > 3"

    cuts["fit-bkkmm2"] = cuts["fit-bkkmm-all"] + " and " + \
        "(Muon_pt[mm_mu1_index[bkkmm_mm_index]] < 4 or " + \
        "Muon_pt[mm_mu2_index[bkkmm_mm_index]] < 3)"

    cuts["fit-bkstarmm-all"] = \
        "mm_mu1_index[bkkmm_mm_index]>=0 and " \
        "mm_mu2_index[bkkmm_mm_index]>=0 and " \
        "Muon_charge[mm_mu1_index[bkkmm_mm_index]] * Muon_charge[mm_mu2_index[bkkmm_mm_index]] < 0 and " \
        "Muon_mediumId[mm_mu1_index[bkkmm_mm_index]] and " \
        "Muon_mediumId[mm_mu2_index[bkkmm_mm_index]] and " \
        "Muon_isGlobal[mm_mu1_index[bkkmm_mm_index]] and " \
        "Muon_isGlobal[mm_mu2_index[bkkmm_mm_index]] and " \
        "mm_kin_vtx_prob[bkkmm_mm_index]>0.01 and " \
        "bkkmm_jpsikk_vtx_prob>0.025 and " \
        "bkkmm_jpsikk_sl3d>3 and " \
        "abs(bkkmm_jpsikk_alpha) < 0.1"

    cuts["fit-bkstarmm"] = cuts["fit-bkstarmm-all"] + " and " + \
        "Muon_pt[mm_mu1_index[bkkmm_mm_index]] > 4 and " + \
        "Muon_pt[mm_mu2_index[bkkmm_mm_index]] > 3"

    cuts["fit-bkstarmm2"] = cuts["fit-bkstarmm-all"] + " and " + \
        "(Muon_pt[mm_mu1_index[bkkmm_mm_index]] < 4 or " + \
        "Muon_pt[mm_mu2_index[bkkmm_mm_index]] < 3)"

    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/BsToJPsiPhi-JpsiToMuMu-PhiToKK_Fil-MuPt2_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/67309e83-705c-44bc-9044-8bcbadf57451.root",
    #     ],
    #     "tree_name" : "bspsiphiMc",
    #     "blind" : False,
    #     "cut" : cuts["fit-bkkmm"],
    #     "final_state" : "bkkmm",
    #     "best_candidate": "",
    #   }  
    
    # job = {
    #     "input": [
    #         "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/4e83dd1a-a5db-42a9-a30a-2f618d520065.root",
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "bspsiphiData",
    #     "blind" : False,
    #     "cut" : (
    #         "mm_mu1_index[bkkmm_mm_index]>=0 and "
    #         "mm_mu2_index[bkkmm_mm_index]>=0 and "
    #         "Muon_charge[mm_mu1_index[bkkmm_mm_index]] * Muon_charge[mm_mu2_index[bkkmm_mm_index]] < 0 and "
    #         "Muon_mediumId[mm_mu1_index[bkkmm_mm_index]] and "
    #         "Muon_mediumId[mm_mu2_index[bkkmm_mm_index]] and "
    #         "Muon_isGlobal[mm_mu1_index[bkkmm_mm_index]] and "
    #         "Muon_isGlobal[mm_mu2_index[bkkmm_mm_index]] and "
    #         "mm_kin_vtx_prob[bkkmm_mm_index]>0.01 and "
    #         "bkkmm_jpsikk_vtx_prob>0.025 and "
    #         "bkkmm_jpsikk_sl3d>3 and "
    #         "abs(bkkmm_jpsikk_alpha) < 0.1 and "        
    #         "abs(bkkmm_jpsikk_mass-5.4)<0.5 and "
    #         "abs(bkkmm_kk_mass-1.02)<0.02"
    #     ),
    #     "final_state" : "bkkmm",
    #     "best_candidate": "",
    #     "pre-selection":"abs(bkkmm_kk_mass-1.02)<0.02&&bkkmm_jpsikk_sl3d>3",
    #     "pre-selection-keep":"^(bkkmm_.*|nbkkmm|mm_.*|nmm|Muon_.*|nMuon|HLT_DoubleMu4_3_LowMass|HLT_DoubleMu2_Jpsi_LowPt|L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5|L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4|L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4|L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2|" + common_branches + ")$",
    #   }  

    job = {
        "input": [
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/BdToJpsiKstar_JpsiToMuMu_MuFilter_Pt-2_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/49ba731a-6beb-4e9f-bda1-6061c0dabc61.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/BsToJPsiPhi-JpsiToMuMu-PhiToKK_Fil-MuPt2_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/67309e83-705c-44bc-9044-8bcbadf57451.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/BsToJPsiPhi-JPsiToMuMu-PhiToKK_Fil-EtaPt_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/3d4290e4-defc-41df-9cce-a3fcee45ee60.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv4-validDigi_130X_mcRun3_2022_realistic_v5-v4+MINIAODSIM/34014200-dc80-44c5-800b-d498f51bdd6f.root"
            # "/tmp/dmytro/49ba731a-6beb-4e9f-bda1-6061c0dabc61.root",
            # "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/ParkingDoubleMuonLowMass6+Run2024F-PromptReco-v1+MINIAOD/001c6ac4-2f19-4931-900d-4b4ac4b0049a.root",
        ],
        "signal_only" : False,
        "tree_name" : "bdpsikstarMc",
        "blind" : False,
        "cut" : cuts["fit-bkstarmm"],
        "triggers": ["HLT_DoubleMu4_3_LowMass"],
        "final_state" : "bkstarmm",
        "best_candidate": "",
        "veto_phi": True,
        "pre-selection":"bkkmm_jpsikk_vtx_prob>0.025 && bkkmm_jpsikk_sl3d>3 && abs(bkkmm_jpsikk_alpha) < 0.1",
      }
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = FlatNtupleForMLFit("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing-NEW/FlatNtuples/532/fit-bkstarmm/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/44c57b0f05d4ca5a0579216260518dc5.job")
    p = FlatNtupleForMLFit(file_name)

    print(p.__dict__)
        
    p.process()
