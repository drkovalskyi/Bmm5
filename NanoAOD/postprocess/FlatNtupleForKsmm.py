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
        self.job_info["pre-selection-keep"] = "^(" + \
            "GenPart_.*|nGenPart|mm_.*|nmm|trk_.*|ntrk|mmiso_.*|nmmiso|" + \
            "Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + \
            "HLT_Mu4_L1DoubleMu|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|" + \
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

        if self.job_info['final_state'] in ['mm']:
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
            
            if self.job_info['final_state'] == "mm" and \
               "mm_extra_info" in self.job_info and \
               self.job_info["mm_extra_info"] == True:
                for var in self.mm_extra_floats:
                    self.tree.addBranch(var, "Float_t", 0)
                for var in self.mm_extra_ints:
                    self.tree.addBranch(var, "Int_t", 0)

            # Dimuon vertex quality block
            self.tree.addBranch('vtx_qual_prob',             'Float_t', 0, "dimuon vertex probability")
            self.tree.addBranch('vtx_qual_mu1_nLostHitsInner', 'Int_t', 0, "Innner lost hits for mu1")
            self.tree.addBranch('vtx_qual_mu2_nLostHitsInner', 'Int_t', 0, "Innner lost hits for mu2")
            self.tree.addBranch('vtx_qual_mu1_veto',           'Int_t', 0, "Number of unexpected pixel layers wits measurements for muon 1")
            self.tree.addBranch('vtx_qual_mu2_veto',           'Int_t', 0, "Number of unexpected pixel layers wits measurements for muon 2")
            self.tree.addBranch('vtx_qual_doca',             'Float_t', 0, "distance between two muons at the point of closest approach")
            self.tree.addBranch('vtx_qual_flightdist_err',   'Float_t', 0, "uncertainty on the flight distance wrt PV in 3D")
            self.tree.addBranch('vtx_qual_flightdist2D_err', 'Float_t', 0, "uncertainty on the flight distance wrt beam spot in 2D")
            self.tree.addBranch('vtx_qual_veto_mask',          'Int_t', 0, "Veto mask based on dimuon vertex position")
            self.tree.addBranch('vtx_qual_mu1_pixel_pattern',  'Int_t', 0, "Pixel layer measurements for mu1")
            self.tree.addBranch('vtx_qual_mu2_pixel_pattern',  'Int_t', 0, "Pixel layer measurements for mu2")
            
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
            self.tree.addBranch('iso_mu_tight_ntracks',            'Int_t', 0, "Number of tracks compatible with either muons (tight)")
            self.tree.addBranch('iso_mu_tight_sum_pt',           'Float_t', 0, "Sum pt of tracks compatible with either muons (tight)")
            self.tree.addBranch('iso_mu_loose_ntracks',            'Int_t', 0, "Number of tracks compatible with either muons (loose)")
            self.tree.addBranch('iso_mu_loose_sum_pt',           'Float_t', 0, "Sum pt of tracks compatible with either muons (loose)")

            # Production information
            self.tree.addBranch('prod_alpha',                    'Float_t', 0, "Pointing angle wrt PV in 3D")
            self.tree.addBranch('prod_alpha_significance',       'Float_t', 0, "Pointing angle significance wrt PV in 3D")
            self.tree.addBranch('prod_alphaBS',                  'Float_t', 0, "Pointing angle wrt beam spot in 2D")
            self.tree.addBranch('prod_alphaBS_significance',     'Float_t', 0, "Pointing angle significance wrt beam spot in 2D")
            self.tree.addBranch('prod_pvip',                     'Float_t', 0, "Impact parameter wrt PV")
            self.tree.addBranch('prod_spvip',                    'Float_t', 0, "Significance of impact parameter wrt PV")
            
            ## heavy flavor tagging
            
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

        self.tree.addBranch('mc_match',      'Int_t', 0, "PdgId of the MC match")
        self.tree.addBranch('mc_vtx_x',    'Float_t', 0, "MC matched gen vertex position")
        self.tree.addBranch('mc_vtx_y',    'Float_t', 0, "MC matched gen vertex position")
        self.tree.addBranch('mc_vtx_z',    'Float_t', 0, "MC matched gen vertex position")

        for trigger in self.triggers_to_store:
            self.tree.addBranch(trigger, 'Int_t', -1, "Trigger decision: 1 - fired, 0 - didn't fire, -1 - no information")
            self.tree.addBranch("%s_ps" % trigger, 'Float_t', 999999, "Prescale. 0 - Off, 999999 - no information")
            self.tree.addBranch("%s_matched" % trigger, 'Int_t', 0,  "matched to the trigger objets")

    def _fill_tree(self, cand, ncands):
        self.tree.reset()

        try:
            ## event info
            self.tree['run'] = self.event.run
            self.tree['ls']  = self.event.luminosityBlock
            self.tree['evt'] = self.event.event
            self.tree['npv'] = ord(self.event.PV_npvsGood) if isinstance(self.event.PV_npvsGood, str) else self.event.PV_npvsGood

            self.tree['certified_muon']   = self._is_certified(self.event, "muon")
            self.tree['certified_golden'] = self._is_certified(self.event, "golden")
            self.tree['n']   = ncands

            ## MC info
            if hasattr(self.event, 'Pileup_nTrueInt'):
                self.tree['npu']      = self.event.Pileup_nPU
                self.tree['npu_mean'] = self.event.Pileup_nTrueInt
        except:
            pass

        if self.job_info['final_state'] == 'mm':
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
            if hasattr(self.event, 'mm_gen_tau'):
                self.tree['gtau'] = self.event.mm_gen_tau[cand]
                self.tree['mc_match'] = self.event.mm_gen_pdgId[cand]
                self.tree['mc_vtx_x'] = self.event.mm_gen_vtx_x[cand]
                self.tree['mc_vtx_y'] = self.event.mm_gen_vtx_y[cand]
                self.tree['mc_vtx_z'] = self.event.mm_gen_vtx_z[cand]

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
            
            try:
                if mu1 >= 0:
                    self.tree['m1q']   = self.event.Muon_charge[mu1]
                    if hasattr(self.event, 'mm_gen_mu1_mpdgId'):
                        self.tree['m1mc']  = self.event.mm_gen_mu1_mpdgId[cand]
                    self.tree['m1bdt'] = self.event.Muon_softMva[mu1]
                    self.tree['vtx_qual_mu1_nLostHitsInner'] = self.event.MuonId_nLostHitsInner[mu1]
                else: 
                    self.tree['m1bdt'] = 1

                if mu2 >= 0:
                    self.tree['m2q']   = self.event.Muon_charge[mu2]
                    if hasattr(self.event, 'mm_gen_mu2_mpdgId'):
                        self.tree['m2mc']  = self.event.mm_gen_mu2_mpdgId[cand]
                    self.tree['m2bdt'] = self.event.Muon_softMva[mu2]
                    self.tree['vtx_qual_mu2_nLostHitsInner'] = self.event.MuonId_nLostHitsInner[mu2]
                else: 
                    self.tree['m2bdt'] = 1

                if mu1 >= 0 and mu2 >= 0:
                    self.tree['muid']   = self.event.Muon_softMvaId[mu1] and self.event.Muon_softMvaId[mu2]
            except:
                pass

            # Dimuon vertex quality block
            self.tree['vtx_qual_prob'] = self.event.mm_kin_vtx_prob[cand]
            veto_mask = expected_missing_hits(self.event.mm_kin_vtx_x[cand],
                                              self.event.mm_kin_vtx_y[cand],
                                              self.event.mm_kin_vtx_z[cand])
            self.tree['vtx_qual_veto_mask'] = veto_mask
            self.tree['vtx_qual_mu1_pixel_pattern'] = self.event.MuonId_pixelPattern[mu1]
            self.tree['vtx_qual_mu2_pixel_pattern'] = self.event.MuonId_pixelPattern[mu2]

            n_veto_mu1 = 0
            n_veto_mu2 = 0
            for i in range(7):
                index = ( 1 << i )
                if ( veto_mask & self.event.MuonId_pixelPattern[mu1] & index ) > 0:
                    n_veto_mu1 += 1
                if ( veto_mask & self.event.MuonId_pixelPattern[mu2] & index ) > 0:
                    n_veto_mu2 += 1
            self.tree['vtx_qual_mu1_veto'] = n_veto_mu1
            self.tree['vtx_qual_mu2_veto'] = n_veto_mu2

            self.tree['vtx_qual_doca'] = self.event.mm_doca[cand]
            if self.event.mm_kin_sl3d[cand] > 0:
                self.tree['vtx_qual_flightdist_err'] = self.event.mm_kin_l3d[cand] / self.event.mm_kin_sl3d[cand]
            if self.event.mm_kin_slxy[cand] > 0:
                self.tree['vtx_qual_flightdist2D_err'] = self.event.mm_kin_lxy[cand] / self.event.mm_kin_slxy[cand]
            
            # Flight length
            self.tree['displacement_l3d']  = self.event.mm_kin_l3d[cand]
            self.tree['displacement_sl3d'] = self.event.mm_kin_sl3d[cand]
            self.tree['displacement_lxy']  = self.event.mm_kin_lxy[cand]
            self.tree['displacement_slxy'] = self.event.mm_kin_slxy[cand]
            massOverC = 0.4976 / 29.979e-3 # GeV ps/cm
            self.tree['displacement_pseudo_decaytime'] = self.event.mm_kin_lxy[cand] / self.event.mm_kin_pt[cand] * massOverC
            
            # Isolation information
            self.tree['iso_vtx_tight_ntracks'] = 0
            self.tree['iso_vtx_tight_sum_pt'] = 0
            self.tree['iso_vtx_loose_ntracks'] = 0
            self.tree['iso_vtx_loose_sum_pt'] = 0
            self.tree['iso_mu_tight_ntracks'] = 0
            self.tree['iso_mu_tight_sum_pt'] = 0
            self.tree['iso_mu_loose_ntracks'] = 0
            self.tree['iso_mu_loose_sum_pt'] = 0
            for i in range(self.event.nmmiso):
                if self.event.mmiso_mm_index[i] != cand:
                    continue
                # Tracks compatible with mm vertex (tight)
                if self.event.mmiso_mm_vtx_prob[i] > 0.1:
                    self.tree['iso_vtx_tight_ntracks'] += 1
                    self.tree['iso_vtx_tight_sum_pt'] += self.event.trk_pt[self.event.mmiso_trk_index[i]]

                # Tracks compatible with mm vertex (loose)
                if self.event.mmiso_mm_vtx_prob[i] > 0.001:
                    self.tree['iso_vtx_loose_ntracks'] += 1
                    self.tree['iso_vtx_loose_sum_pt'] += self.event.trk_pt[self.event.mmiso_trk_index[i]]

                # Tracks compatible with either muon (tight)
                if self.event.mmiso_mu1_vtx_prob[i] > 0.1 or self.event.mmiso_mu2_vtx_prob[i] > 0.1:
                    self.tree['iso_mu_tight_ntracks'] += 1
                    self.tree['iso_mu_tight_sum_pt'] += self.event.trk_pt[self.event.mmiso_trk_index[i]]

                # Tracks compatible with either muon (loose)
                if self.event.mmiso_mu1_vtx_prob[i] > 0.001 or self.event.mmiso_mu2_vtx_prob[i] > 0.001:
                    self.tree['iso_mu_loose_ntracks'] += 1
                    self.tree['iso_mu_loose_sum_pt'] += self.event.trk_pt[self.event.mmiso_trk_index[i]]

            # Production information
            self.tree['prod_alpha'] = self.event.mm_kin_alpha[cand]
            if self.event.mm_kin_alphaErr[cand] > 0:
                self.tree['prod_alpha_significance'] = self.event.mm_kin_alpha[cand] / self.event.mm_kin_alphaErr[cand]
            self.tree['prod_alphaBS'] = self.event.mm_kin_alphaBS[cand]
            if self.event.mm_kin_alphaBSErr[cand] > 0:
                self.tree['prod_alphaBS_significance'] = self.event.mm_kin_alphaBS[cand] / self.event.mm_kin_alphaBSErr[cand]
            self.tree['prod_pvip'] = self.event.mm_kin_pvip[cand]
            self.tree['prod_spvip'] = self.event.mm_kin_spvip[cand]
                    
        elif self.job_info['final_state'] == 'hh':
            # hh
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
            if hasattr(self.event, 'hh_gen_tau'):
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

        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        if self.job_info['final_state'] in ['mm']:
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
        elif self.job_info['final_state'] == 'hh':
            for trigger in self.triggers_to_store:
                if hasattr(self.event, trigger):
                    self.tree[trigger] = getattr(self.event, trigger)
                if hasattr(self.event, "prescale_" + trigger):
                    self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
        else:
            raise Exception("Unsupported final state: %s" % self.job_info['final_state'])

        self.tree.fill()

if __name__ == "__main__":

    ### create a test job

    common_branches = "PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock"

    input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/530/ParkingDoubleMuonLowMass0+Run2022D-PromptReco-v1+MINIAOD/"
    input_path2 = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/530/ksmm/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv4-validDigi_130X_mcRun3_2022_realistic_v5-v4+MINIAODSIM/"
    
    job = {
        "input": [
            # "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/530/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/28491152-7072-4c6d-ac38-6e20c57ae346.root",
            # input_path + "01340b25-f1c3-40bd-8c09-225da18b30e3.root",
            # input_path + "2fa3ea36-eb49-41aa-8496-71de4ed04620.root",
            # input_path + "632b1e27-5f2b-49f3-acb9-c47e545e7547.root",
            # input_path + "99f4a5dc-cb1b-41a7-92a8-959acb241375.root",
            # input_path + "c84cc13c-3994-4b22-a606-64ff50e167da.root",
            # input_path + "02773df7-b597-49bb-9512-3b34f4fde72c.root",
            # input_path + "2fbdab90-acb6-4b55-83a9-631ad71a1611.root",
            # input_path + "659b27a1-4b2e-4aa9-9746-8ac41530e413.root",
            # input_path + "9a54368e-e153-4ec9-861a-5c1238306cf9.root",
            # input_path + "c868ac55-5c22-4b75-9fae-cbd3592d18b3.root",
            # input_path + "0434f813-58bf-4939-a1de-7a9316938df7.root",
            # input_path + "30aa4ecc-ba74-4241-8cc1-f04e95294670.root",
            # input_path + "66783f3e-8d0d-43ec-97fe-c10dfa92e095.root",
            # input_path + "9bd3f34b-fa42-447b-b6fc-b046060bf6ec.root",
            # input_path + "c9ddcac3-5d7c-43bc-91e1-a0b23b46462c.root"
            
            input_path2 + "a066bd55dfc8b89df455746c28853bf6.root",
            input_path2 + "a11192928591ecc63d2c983ad3bafd46.root",
            input_path2 + "a1493042e87f67fb1c77ba66dd66bfbc.root",
            input_path2 + "a24546b594e05f1109d204eaf8960b88.root",
            input_path2 + "a362633465f127f5465f3a16789f64af.root",
            input_path2 + "a3765ef6de92d86864f1d4f0dc739a15.root",
            input_path2 + "a4d22291b7b95cdcab895eefd1d29395.root",
            input_path2 + "a944734ffccf72a9572ba5698052608c.root",
            input_path2 + "a99a7511d9a696dd11ea0c6545ffdb07.root",
            input_path2 + "a9e0c20dd3d627a296cef163ad57f2ea.root",
            input_path2 + "aa4624f2d3640df442dc87c48ba25728.root",
            input_path2 + "ab13185a8760930af5fe25a052dc9144.root",
            input_path2 + "ad357ad845d1324792d177631feca87f.root",
            input_path2 + "ad6fd14f5247ad19e83b220b223f07e7.root",
            input_path2 + "adc8d7b00a80f92c9376513ccc85b737.root",
            input_path2 + "ae93fb1b17d79d80f3f4aa2faf3b4fa0.root",
            input_path2 + "aecebfc80dc8c2dfddc0f42c7ea4fa9c.root",
            input_path2 + "af1d28e789e75504c38b9f9fa865829a.root"
        ],
        "signal_only" : False,
        "tree_name" : "ksmm",
        "blind" : False,
        "cut" :
            "mm_mu1_index>=0 and mm_mu2_index>=0 and "
            "Muon_charge[mm_mu1_index] * Muon_charge[mm_mu2_index] < 0 and "
            "abs(mm_kin_mass-0.45)<0.2 and "
            "mm_kin_slxy>10 and mm_kin_lxy>1 and mm_kin_alpha<0.1 and "
            "mm_kin_vtx_prob>0.01 and "
            "HLT_DoubleMu4_3_LowMass",
        "final_state" : "mm",
        "best_candidate": "",
        "mm_extra_info": True,
        "pre-selection":"abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01 && mm_kin_slxy>3 && mm_kin_alpha<0.1", 
      }


    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = FlatNtupleForKsmm("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/526/ksmm/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v1+MINIAODSIM/73c7e720b2098d39f840aaabb26a4d6b.job")
    p = FlatNtupleForKsmm(file_name)

    print(p.__dict__)
        
    p.process()
