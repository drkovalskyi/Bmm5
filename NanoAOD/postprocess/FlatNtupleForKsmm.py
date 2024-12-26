from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
import multiprocessing
from datetime import datetime
import hashlib
import numpy as np

def compute_dip_angle(x, y, z):
    pt = np.sqrt(x**2 + y**2)
    dip_angle = np.arctan2(z, pt)
    
    return dip_angle

def compute_pseudo_rapidity(x, y, z):
    pt = np.sqrt(x**2 + y**2)
    theta = np.arctan2(pt, z)
    eta = -np.log(np.tan(theta / 2.0))
    
    return eta

def expected_missing_hits(vtx_x, vtx_y, vtx_z, pixel_2017=True):
    """
    Compute max expected missing layer measurements in front of vertex
    Mask definition: barrel 0b1111, endcap 0b1110000
    least significant bits represent inner layers/disks

    WARNING: very rough pixel geometry approximation
    """
    mask = 0b0000000
    rho = math.sqrt(pow(vtx_x, 2) + pow(vtx_y, 2))
    dip_angle = np.degrees(compute_dip_angle(vtx_x, vtx_y, vtx_z))
    tolerance = 0.5
    
    if pixel_2017:
        # ignore beam line area
        if rho < 3.0 + tolerance:
            return 0
        elif rho < 6.8 + tolerance:
            if abs(vtx_z) < 29.1 + tolerance:
                return 0b0000001
            if abs(vtx_z) < 39.6 + tolerance:
                return 0b0010001
            if abs(vtx_z) < 51.6 + tolerance:
                return 0b0110001
            else:
                return 0b1111111
        elif rho < 10.9 + tolerance:
            if abs(vtx_z) < 29.1 + tolerance:
                return 0b0000011
            if abs(vtx_z) < 39.6 + tolerance:
                return 0b0010011
            if abs(vtx_z) < 51.6 + tolerance:
                return 0b0110011
            else:
                return 0b1111111
        elif rho < 16.0 + tolerance:
            if abs(vtx_z) < 29.1 + tolerance:
                return 0b0000111
            if abs(vtx_z) < 39.6 + tolerance:
                return 0b0010111
            if abs(vtx_z) < 51.6 + tolerance:
                return 0b0110111
            else:
                return 0b1111111
        else:
            return 0b1111111
    return 0
            

class FlatNtupleForKsmm(FlatNtupleBase):
    """Flat ROOT ntuple producer for Bmm5 UML fit"""

    final_states = [
        'mm', 'hh'
    ]
    
    leaf_counts = { 
        'mm':    'nmm',
        'hh':    'nhh',
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
		       "mm_otherVtxMaxProb1", "mm_otherVtxMaxProb2", "mm_kin_slxy", "mm_kin_lxy",
                       "mm_kin_vtx_z",
                       # "mm_gen_vtx_x", "mm_gen_vtx_y", "mm_gen_vtx_z"
                       ]
    mm_extra_ints = ["mm_nBMTrks",]

    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing information
        if 'best_candidate' not in self.job_info:
            self.job_info['best_candidate'] = ""

        for parameter in ['input', 'cut', 'final_state']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)
            
            if self.job_info['final_state'] not in self.final_states:
                raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        # set branches to be kept if pre-skimmed
        if "pre-selection-keep" not in self.job_info:
            self.job_info["pre-selection-keep"] = "^(" + \
                "GenPart_.*|nGenPart|mm_.*|nmm|hh_.*|nhh|trk_.*|ntrk|iso_.*|niso|" + \
                "Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + \
                "HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|HLT_ZeroBias|" + \
                "PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock" + \
                ")$"

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
        self.tree.addBranch('pt',         'Float_t', 0, "pt of the Ks candidate")
        self.tree.addBranch('eta',        'Float_t', 0, "eta of the Ks candidate")
        self.tree.addBranch('phi',        'Float_t', 0, "phi of the Ks candidate")
        self.tree.addBranch('m',          'Float_t', 0, "mass of the Ks candidate")
        self.tree.addBranch('m_raw',      'Float_t', 0, "raw mass of the Ks candidate")
        self.tree.addBranch('me',         'Float_t', 0, "mass error of the Ks candidate")
        self.tree.addBranch('tau',        'Float_t', 0, "decay time of the Ks candidate")
        self.tree.addBranch('taue',       'Float_t', 0, "decay time error of the Ks candidate")
        self.tree.addBranch('tauxy',      'Float_t', 0, "2D decay time")
        self.tree.addBranch('tauxye',     'Float_t', 0, "2D decay time error")
        self.tree.addBranch('gtau',       'Float_t', 0, "generated decay time")

        self.tree.addBranch('chan',        'UInt_t', 0, "0 if |eta|<0.7 for both muons. 1 otherwise")

        self.tree.addBranch('muid',        'UInt_t', 0, "MVA muon selection for both muons")

        self.tree.addBranch('d1pt',       'Float_t', 0)
        self.tree.addBranch('d1eta',      'Float_t', 0)
        self.tree.addBranch('d1phi',      'Float_t', 0)
        self.tree.addBranch('d1q',          'Int_t', 0, "First daughter charge")
        self.tree.addBranch('d1mc_pdgId',   'Int_t', 0, "PDG id of the MC match")
        self.tree.addBranch('d1mc_mpdgId',  'Int_t', 0, "PDG id of the mother of the MC match")
        self.tree.addBranch('d1muId',     'Float_t', 0, "Muon MVA value. 0 for hadrons")
        
        self.tree.addBranch('d2pt',       'Float_t', 0)
        self.tree.addBranch('d2eta',      'Float_t', 0)
        self.tree.addBranch('d2phi',      'Float_t', 0)
        self.tree.addBranch('d2q',          'Int_t', 0, "Second daughter charge")
        self.tree.addBranch('d2mc_pdgId',   'Int_t', 0, "PDG id of the MC match")
        self.tree.addBranch('d2mc_mpdgId',  'Int_t', 0, "PDG id of the mother of the MC match")
        self.tree.addBranch('d2muId',     'Float_t', 0, "Muon MVA value. 0 for hadrons")
            
        if self.job_info['final_state'] == "mm" and \
           "mm_extra_info" in self.job_info and \
           self.job_info["mm_extra_info"] == True:
            for var in self.mm_extra_floats:
                self.tree.addBranch(var, "Float_t", 0)
            for var in self.mm_extra_ints:
                self.tree.addBranch(var, "Int_t", 0)


        # Vertex quality block
        self.tree.addBranch('vtx_qual_prob',             'Float_t', 0, "dimuon vertex probability")
        
        self.tree.addBranch('vtx_qual_d1_nLostHitsInner',  'Int_t', 0, "Innner lost hits for daughter 1")
        self.tree.addBranch('vtx_qual_d2_nLostHitsInner',  'Int_t', 0, "Innner lost hits for daughter 2")
        self.tree.addBranch('vtx_qual_d1_veto',            'Int_t', 0, "Number of unexpected pixel layers wits measurements for daughter 1")
        self.tree.addBranch('vtx_qual_d2_veto',            'Int_t', 0, "Number of unexpected pixel layers wits measurements for daughter 2")
        self.tree.addBranch('vtx_qual_d1_pixel_pattern',   'Int_t', 0, "Pixel layer measurements for daughter 1")
        self.tree.addBranch('vtx_qual_d2_pixel_pattern',   'Int_t', 0, "Pixel layer measurements for daughter 2")
        
        self.tree.addBranch('vtx_qual_doca',             'Float_t', 0, "distance between two muons at the point of closest approach")
        self.tree.addBranch('vtx_qual_flightdist_err',   'Float_t', 0, "uncertainty on the flight distance wrt PV in 3D")
        self.tree.addBranch('vtx_qual_flightdist2D_err', 'Float_t', 0, "uncertainty on the flight distance wrt beam spot in 2D")
        self.tree.addBranch('vtx_qual_veto_mask',          'Int_t', 0, "Veto mask based on dimuon vertex position")
            
        # Flight length
        self.tree.addBranch('displacement_l3d',              'Float_t', 0, "Distance from dimuon vertex to PV in 3D")
        self.tree.addBranch('displacement_sl3d',             'Float_t', 0, "Significance of distance from dimuon vertex to PV in 3D")
        self.tree.addBranch('displacement_lxy',              'Float_t', 0, "Distance from dimuon vertex to beam spot in 2D")
        self.tree.addBranch('displacement_slxy',             'Float_t', 0, "Significance of distance from dimuon vertex to beam spot in 2D")
        self.tree.addBranch('displacement_pseudo_decaytime', 'Float_t', 0, "Proper decay time caclulated with the PDG Ks mass value (ps)")

        # Isolation information
        self.tree.addBranch('iso_vtx_tight_ntracks',           'Int_t', 0, "Number of tracks compatible with dimuon vertex (tight)")
        self.tree.addBranch('iso_vtx_tight_sum_pt',          'Float_t', 0, "Sum pt of tracks compatible with dimuon vertex (tight)")
        self.tree.addBranch('iso_vtx_loose_ntracks',           'Int_t', 0, "Number of tracks compatible with dimuon vertex (loose)")
        self.tree.addBranch('iso_vtx_loose_sum_pt',          'Float_t', 0, "Sum pt of tracks compatible with dimuon vertex (loose)")
        self.tree.addBranch('iso_tight_ntracks',            'Int_t', 0, "Number of tracks compatible with either muons (tight)")
        self.tree.addBranch('iso_tight_sum_pt',           'Float_t', 0, "Sum pt of tracks compatible with either muons (tight)")
        self.tree.addBranch('iso_loose_ntracks',            'Int_t', 0, "Number of tracks compatible with either muons (loose)")
        self.tree.addBranch('iso_loose_sum_pt',           'Float_t', 0, "Sum pt of tracks compatible with either muons (loose)")

        # Production information
        self.tree.addBranch('prod_alpha',                    'Float_t', 0, "Pointing angle wrt PV in 3D")
        self.tree.addBranch('prod_alpha_significance',       'Float_t', 0, "Pointing angle significance wrt PV in 3D")
        self.tree.addBranch('prod_alphaBS',                  'Float_t', 0, "Pointing angle wrt beam spot in 2D")
        self.tree.addBranch('prod_alphaBS_significance',     'Float_t', 0, "Pointing angle significance wrt beam spot in 2D")
        self.tree.addBranch('prod_pvip',                     'Float_t', 0, "Impact parameter wrt PV")
        self.tree.addBranch('prod_spvip',                    'Float_t', 0, "Significance of impact parameter wrt PV")
            
        ## heavy flavor tagging
        # FIXME
            
        if self.job_info['final_state'] not in ['mm', 'hh']:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        self.tree.addBranch('mc_match',      'Int_t', 0, "PdgId of the MC match")
        self.tree.addBranch('mc_vtx_x',    'Float_t', 0, "MC matched gen vertex position")
        self.tree.addBranch('mc_vtx_y',    'Float_t', 0, "MC matched gen vertex position")
        self.tree.addBranch('mc_vtx_z',    'Float_t', 0, "MC matched gen vertex position")

        for trigger in self.triggers_to_store:
            self.tree.addBranch(trigger, 'Int_t', -1, "Trigger decision: 1 - fired, 0 - didn't fire, -1 - no information")
            self.tree.addBranch("%s_ps" % trigger, 'Float_t', 999999, "Prescale. 0 - Off, 999999 - no information")
            self.tree.addBranch("%s_matched" % trigger, 'Int_t', 0,  "matched to the trigger objets")

    def _get(self, name):
        final_state = self.job_info['final_state']
        return getattr(self.event, f"{final_state}_{name}")
    
    def _fill_tree(self, cand, ncands):
        self.tree.reset()

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

        final_state = self.job_info['final_state']
        if final_state in ['mm', 'hh']:
            self.tree['bdt']    = self._get("mva")[cand]
            self.tree['pt']     = self._get("kin_pt")[cand]
            self.tree['eta']    = self._get("kin_eta")[cand]
            self.tree['phi']    = self._get("kin_phi")[cand]
            self.tree['m']      = self._get("kin_mass")[cand]
            self.tree['m_raw']  = self._get("mass")[cand]
            self.tree['me']     = self._get("kin_massErr")[cand]
            self.tree['tau']    = self._get("kin_tau")[cand]
            self.tree['taue']   = self._get("kin_taue")[cand]
            self.tree['tauxy']  = self._get("kin_tauxy")[cand]
            self.tree['tauxye'] = self._get("kin_tauxye")[cand]
            
            if mc:
                self.tree['gtau'] = self._get("gen_tau")[cand]
                self.tree['mc_match'] = self._get("gen_pdgId")[cand]
                self.tree['mc_vtx_x'] = self._get("gen_vtx_x")[cand]
                self.tree['mc_vtx_y'] = self._get("gen_vtx_y")[cand]
                self.tree['mc_vtx_z'] = self._get("gen_vtx_z")[cand]

            mu1 = -1
            mu2 = -1
            
            if final_state == 'mm':
                if "mm_extra_info" in self.job_info and \
                   self.job_info["mm_extra_info"] == True:
                    for var in self.mm_extra_floats:
                        self.tree[var] = getattr(self.event, var)[cand]
                    for var in self.mm_extra_ints:
                        self.tree[var] = getattr(self.event, var)[cand]

                mu1 = self.event.mm_mu1_index[cand]
                mu2 = self.event.mm_mu2_index[cand]

                self.tree['d1pt']  = self.event.mm_mu1_pt[cand]
                self.tree['d1eta'] = self.event.mm_mu1_eta[cand]
                self.tree['d1phi'] = self.event.mm_mu1_phi[cand]
                self.tree['d2pt']  = self.event.mm_mu2_pt[cand]
                self.tree['d2eta'] = self.event.mm_mu2_eta[cand]
                self.tree['d2phi'] = self.event.mm_mu2_phi[cand]
            
                if mu1 >= 0:
                    self.tree['d1q']   = self.event.Muon_charge[mu1]
                    if mc:
                        self.tree['d1mc_pdgId']   = self.event.mm_gen_mu1_pdgId[cand]
                        self.tree['d1mc_mpdgId']  = self.event.mm_gen_mu1_mpdgId[cand]
                    self.tree['d1muId'] = self.event.Muon_softMva[mu1]
                    self.tree['vtx_qual_d1_nLostHitsInner'] = self.event.MuonId_nLostHitsInner[mu1]
                else: 
                    self.tree['d1muId'] = 0

                if mu2 >= 0:
                    self.tree['d2q']   = self.event.Muon_charge[mu2]
                    if mc:
                        self.tree['d2mc_pdgId']   = self.event.mm_gen_mu2_pdgId[cand]
                        self.tree['d2mc_mpdgId']  = self.event.mm_gen_mu2_mpdgId[cand]
                    self.tree['d2muId'] = self.event.Muon_softMva[mu2]
                    self.tree['vtx_qual_d2_nLostHitsInner'] = self.event.MuonId_nLostHitsInner[mu2]
                else: 
                    self.tree['d2muId'] = 0
                
                if mu1 >= 0 and mu2 >= 0:
                    self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
                
                self.tree['vtx_qual_d1_pixel_pattern'] = self.event.MuonId_pixelPattern[mu1]
                self.tree['vtx_qual_d2_pixel_pattern'] = self.event.MuonId_pixelPattern[mu2]
            else:
                ## hh
                self.tree['d1pt']  = self.event.hh_had1_pt[cand]
                self.tree['d1eta'] = self.event.hh_had1_eta[cand]
                self.tree['d1phi'] = self.event.hh_had1_phi[cand]
                self.tree['d2pt']  = self.event.hh_had2_pt[cand]
                self.tree['d2eta'] = self.event.hh_had2_eta[cand]
                self.tree['d2phi'] = self.event.hh_had2_phi[cand]
            
                self.tree['d1q']   = self.event.hh_had1_pdgId[cand] > 0
                self.tree['d2q']   = self.event.hh_had2_pdgId[cand] > 0
                if mc:
                    self.tree['d1mc_pdgId']   = self.event.hh_gen_had1_pdgId[cand]
                    self.tree['d1mc_mpdgId']  = self.event.hh_gen_had1_mpdgId[cand]
                    self.tree['d2mc_pdgId']   = self.event.hh_gen_had2_pdgId[cand]
                    self.tree['d2mc_mpdgId']  = self.event.hh_gen_had2_mpdgId[cand]
                self.tree['d1muId'] = 0
                self.tree['d2muId'] = 0
                self.tree['vtx_qual_d1_nLostHitsInner'] = 0 # FIXME
                self.tree['vtx_qual_d2_nLostHitsInner'] = 0 # FIXME
                self.tree['vtx_qual_d1_pixel_pattern'] = 0  # FIXME
                self.tree['vtx_qual_d2_pixel_pattern'] = 0  # FIXME
                
                self.tree['muid']   = 0


            # Vertex quality block
            self.tree['vtx_qual_prob'] = self._get("kin_vtx_prob")[cand]
            veto_mask = expected_missing_hits(self._get("kin_vtx_x")[cand],
                                              self._get("kin_vtx_y")[cand],
                                              self._get("kin_vtx_z")[cand])
            self.tree['vtx_qual_veto_mask'] = veto_mask
            

            n_veto_d1 = 0
            n_veto_d2 = 0
            for i in range(7):
                index = ( 1 << i )
                if ( veto_mask & self.tree['vtx_qual_d1_pixel_pattern'] & index ) > 0:
                    n_veto_d1 += 1
                if ( veto_mask & self.tree['vtx_qual_d2_pixel_pattern'] & index ) > 0:
                    n_veto_d2 += 1
            self.tree['vtx_qual_d1_veto'] = n_veto_d1
            self.tree['vtx_qual_d2_veto'] = n_veto_d2

            self.tree['vtx_qual_doca'] = self._get("doca")[cand]
            if self._get("kin_sl3d")[cand] > 0:
                self.tree['vtx_qual_flightdist_err'] = self._get("kin_l3d")[cand] / self._get("kin_sl3d")[cand]
            if self._get("kin_slxy")[cand] > 0:
                self.tree['vtx_qual_flightdist2D_err'] = self._get("kin_lxy")[cand] / self._get("kin_slxy")[cand]
            
            # Flight length
            self.tree['displacement_l3d']  = self._get("kin_l3d")[cand]
            self.tree['displacement_sl3d'] = self._get("kin_sl3d")[cand]
            self.tree['displacement_lxy']  = self._get("kin_lxy")[cand]
            self.tree['displacement_slxy'] = self._get("kin_slxy")[cand]
            massOverC = 0.4976 / 29.979e-3 # GeV ps/cm
            self.tree['displacement_pseudo_decaytime'] = self._get("kin_lxy")[cand] / self._get("kin_pt")[cand] * massOverC
            
            # Isolation information
            self.tree['iso_vtx_tight_ntracks'] = 0
            self.tree['iso_vtx_tight_sum_pt'] = 0
            self.tree['iso_vtx_loose_ntracks'] = 0
            self.tree['iso_vtx_loose_sum_pt'] = 0
            self.tree['iso_tight_ntracks'] = 0
            self.tree['iso_tight_sum_pt'] = 0
            self.tree['iso_loose_ntracks'] = 0
            self.tree['iso_loose_sum_pt'] = 0
            for i in range(self.event.niso):
                if final_state == 'mm' and self.event.iso_mm_index[i] != cand:
                    continue
                if final_state == 'hh' and self.event.iso_hh_index[i] != cand:
                    continue
                # Tracks compatible with mm vertex (tight)
                if self.event.iso_vtx_prob[i] > 0.1:
                    self.tree['iso_vtx_tight_ntracks'] += 1
                    self.tree['iso_vtx_tight_sum_pt'] += self.event.trk_pt[self.event.iso_trk_index[i]]

                # Tracks compatible with mm vertex (loose)
                if self.event.iso_vtx_prob[i] > 0.001:
                    self.tree['iso_vtx_loose_ntracks'] += 1
                    self.tree['iso_vtx_loose_sum_pt'] += self.event.trk_pt[self.event.iso_trk_index[i]]

                if self.event.iso_dr[i] < 1.5 and \
                   (self.event.iso_d1_doca[i] < 0.02 or self.event.iso_d2_doca[i] < 0.02) and \
                   self.event.trk_sip[self.event.iso_trk_index[i]] > 2.0 and \
                   self.event.trk_pt[self.event.iso_trk_index[i]] > 1.0:

                    # Tracks compatible with either muon (tight)
                    if self.event.iso_d1_vtx_prob[i] > 0.1 or self.event.iso_d2_vtx_prob[i] > 0.1:
                        self.tree['iso_tight_ntracks'] += 1
                        self.tree['iso_tight_sum_pt'] += self.event.trk_pt[self.event.iso_trk_index[i]]

                    # Tracks compatible with either muon (loose)
                    if self.event.iso_d1_vtx_prob[i] > 0.001 or self.event.iso_d2_vtx_prob[i] > 0.001:
                        self.tree['iso_loose_ntracks'] += 1
                        self.tree['iso_loose_sum_pt'] += self.event.trk_pt[self.event.iso_trk_index[i]]

            # Production information
            self.tree['prod_alpha'] = self._get("kin_alpha")[cand]
            if self._get("kin_alphaErr")[cand] > 0:
                self.tree['prod_alpha_significance'] = self._get("kin_alpha")[cand] / self._get("kin_alphaErr")[cand]
            self.tree['prod_alphaBS'] = self._get("kin_alphaBS")[cand]
            if self._get("kin_alphaBSErr")[cand] > 0:
                self.tree['prod_alphaBS_significance'] = self._get("kin_alphaBS")[cand] / self._get("kin_alphaBSErr")[cand]
            self.tree['prod_pvip']  = self._get("kin_pvip")[cand]
            self.tree['prod_spvip'] = self._get("kin_spvip")[cand]

            # Other information
            if abs(self.tree['d1eta']) < 0.7 and abs(self.tree['d2eta']) < 0.7:
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
                        self.tree[trigger + "_matched"] = getattr(self.event, "MuonId_" + trigger)[mu1] and \
                            getattr(self.event, "MuonId_" + trigger)[mu2] 
        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        self.tree.fill()

if __name__ == "__main__":

    ### create a test job

    common_branches = "PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock"

    # job = {
    #     "input": [
    #         "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4+MINIAODSIM/affd8ce1-091d-41cb-af18-2f1d2e9c7e3d.root"
    #     ],
    #     "signal_only" : False,
    #     "tree_name" : "ksmm",
    #     "blind" : False,
    #     "cut" :
    #         "mm_mu1_index>=0 and mm_mu2_index>=0 and "
    #         "Muon_charge[mm_mu1_index] * Muon_charge[mm_mu2_index] < 0 and "
    #         "abs(mm_kin_mass-0.45)<0.2 and "
    #         "mm_kin_slxy>10 and mm_kin_lxy>1 and mm_kin_alpha<0.1 and "
    #         "mm_kin_vtx_prob>0.01 and "
    #         "HLT_DoubleMu4_3_LowMass",
    #     "final_state" : "mm",
    #     "best_candidate": "",
    #     "mm_extra_info": True,
    #     "pre-selection":"abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01 && mm_kin_slxy>3 && mm_kin_alpha<0.1", 
    #   }

    job = {
        "input": [
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/531/ZeroBias+Run2024C-PromptReco-v1+MINIAOD/e4d100e2-3a51-4097-a3d4-4ec637e57e76.root",
        ],
        "signal_only" : False,
        "tree_name" : "kspipiData",
        "blind" : False,
        "cut" :
           "hh_had1_pdgId * hh_had2_pdgId == -211*211 and hh_had1_pt>4 and "
           "hh_had2_pt>3 and abs(hh_kin_mass-0.50)<0.15 and hh_kin_slxy>3 and "
           "hh_kin_alpha<0.1 and hh_kin_vtx_prob>0.01",
        "final_state" : "hh",
        "best_candidate": "",
        "triggers": ["HLT_ZeroBias"],
        "mm_extra_info": False,
        "pre-selection":"abs(hh_kin_mass-0.50)<0.15 && hh_kin_slxy>3 && hh_kin_alpha<0.1",
      }
    

    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = FlatNtupleForKsmm("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/531/ksmm/ParkingDoubleMuonLowMass6+Run2024F-PromptReco-v1+MINIAOD/7d7d0aeac7a4f0c892f453f1d71c72b9.job")
    # p = FlatNtupleForKsmm("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing-NEW/FlatNtuples/531/ksmm/ParkingDoubleMuonLowMass7+Run2022F-PromptReco-v1+MINIAOD/792208be3d8481d404f870f8f7058f2e.job")
    # p = FlatNtupleForKsmm("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/531/kspipi/ZeroBias+Run2024C-PromptReco-v1+MINIAOD/f95441959fea119e92b59b9258f646e4.job")
    p = FlatNtupleForKsmm("/tmp/dmytro/7d7d0aeac7a4f0c892f453f1d71c72b9.job")
    
    # p = FlatNtupleForKsmm(file_name)

    print(p.__dict__)
        
    p.process()
