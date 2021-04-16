from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, math, json
import multiprocessing
from datetime import datetime
import hashlib

class FlatNtupleForMuonMVA(FlatNtupleBase):
    """Flat ROOT ntuple producer for Bmm5 Muon MVA"""

    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing information
        for parameter in ['input', 'tree_name']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

    def _get_jpsis(self):
        """Get clean jpsi to be used as a control region"""

            
    def _process_events(self):
        """Event loop"""

        for event_index, event in enumerate(self.input_tree):
            if self.limit > 0 and event_index >= self.limit: break
            
            self.event = event
            jpsi_candidates = self._get_jpsis()
            
            candidates = []

            for i in range(self.event.nMuon):
                if not self.event.Muon_isTracker[i]: continue
                if not self.event.Muon_isGlobal[i]: continue

                keep_muon = False

                data = dict()
                # MC specific requirements
                if hasattr(self.event, 'MuonId_simPdgId'):
                    if abs(self.event.MuonId_simPdgId[i]) == 13 and \
                       abs(self.event.MuonId_simMotherPdgId[i]) in [211, 321, 443, 531]:
                        keep_muon = True
                        data['sim_pdgId']  = self.event.MuonId_simPdgId[i]
                        data['sim_mpdgId'] = self.event.MuonId_simMotherPdgId[i]
                        if abs(self.event.MuonId_simMotherPdgId[i]) in [211, 321]:
                            data['sim_type'] = 1
                        elif abs(self.event.MuonId_simMotherPdgId[i]) in [531]:
                            data['sim_type'] = 2
                        elif abs(self.event.MuonId_simMotherPdgId[i]) in [443]:
                            data['sim_type'] = 3
                            
                # Add Jpsi muons

                
                if keep_muon:
                    self._fill_tree(i, data)

    def _configure_output_tree(self):
        ## event info
        self.tree.addBranch('run',         'UInt_t', 0)
        self.tree.addBranch('evt',      'ULong64_t', 0)
        self.tree.addBranch('ls',          'UInt_t', 0)

        self.tree.addBranch('pt',         'Float_t', 0)
        self.tree.addBranch('eta',        'Float_t', 0)
        self.tree.addBranch('phi',        'Float_t', 0)

        self.tree.addBranch('bmmBaseId',   'UInt_t', 0)
        self.tree.addBranch('looseId',     'UInt_t', 0)
        self.tree.addBranch('mediumId',    'UInt_t', 0)
        self.tree.addBranch('tightId',     'UInt_t', 0)
        self.tree.addBranch('softId',      'UInt_t', 0)
        self.tree.addBranch('softMvaId',   'UInt_t', 0)
        self.tree.addBranch('softMva',    'Float_t', 0)

        self.tree.addBranch('trkKink',             'Float_t', 0, "Inner track kink chi2")
        self.tree.addBranch('glbTrackProbability', 'Float_t', 0, "Log probability of the global fit")
        self.tree.addBranch('chi2LocalPosition',   'Float_t', 0, "chi2 for STA-TK matching by local position")
        self.tree.addBranch('glbNormChi2',         'Float_t', 0, "Normalized chi2 of the global fit")
        self.tree.addBranch('trkValidFrac',        'Float_t', 0, "Fraction of valid hits for inner track")
        self.tree.addBranch('match1_dX',           'Float_t', 0, "Station 1 local segment-track dX")
        self.tree.addBranch('match1_pullX',        'Float_t', 0, "Station 1 local segment-track dX/dErr")
        self.tree.addBranch('match1_pullDxDz',     'Float_t', 0, "Station 1 local segment-track direction matching in x")
        self.tree.addBranch('match1_dY',           'Float_t', 0, "Station 1 local segment-track dY")
        self.tree.addBranch('match1_pullY',        'Float_t', 0, "Station 1 local segment-track dY/dErr")
        self.tree.addBranch('match1_pullDyDz',     'Float_t', 0, "Station 1 local segment-track direction matching in y")
        self.tree.addBranch('match2_dX',           'Float_t', 0, "Station 2 local segment-track dX")
        self.tree.addBranch('match2_pullX',        'Float_t', 0, "Station 2 local segment-track dX/dErr")
        self.tree.addBranch('match2_pullDxDz',     'Float_t', 0, "Station 2 local segment-track direction matching in x")
        self.tree.addBranch('match2_dY',           'Float_t', 0, "Station 2 local segment-track dY")
        self.tree.addBranch('match2_pullY',        'Float_t', 0, "Station 2 local segment-track dY/dErr")
        self.tree.addBranch('match2_pullDyDz',     'Float_t', 0, "Station 2 local segment-track direction matching in y")

        self.tree.addBranch('nPixels',              'UInt_t', 0, "Number of valid pixel hits")
        self.tree.addBranch('nValidHits',           'UInt_t', 0, "Number of valid hits")
        self.tree.addBranch('nLostHitsInner',       'UInt_t', 0, "Number of lost hits before tracker track")
        self.tree.addBranch('nLostHitsOn',          'UInt_t', 0, "Number of lost hits on tracker track")
        self.tree.addBranch('nLostHitsOuter',       'UInt_t', 0, "Number of lost hits after tracker track")

        self.tree.addBranch('trkLayers',            'UInt_t', 0, "Number of layers with measurements in tracker track")
        self.tree.addBranch('trkLostLayersInner',   'UInt_t', 0, "Number of lost layers befor tracker track")
        self.tree.addBranch('trkLostLayersOn',      'UInt_t', 0, "Number of lost layers on tracker track")
        self.tree.addBranch('trkLostLayersOuter',   'UInt_t', 0, "Number of lost layers after tracker track")

        self.tree.addBranch('highPurity',           'UInt_t', 0, "High purity inner track")
        
        self.tree.addBranch('sim_type',             'UInt_t', 0, "1 - decay in flight, 2 - heavy flavor decays, 3 - jpsi")
        self.tree.addBranch('sim_pdgId',             'Int_t', 0, "PDG id of SIM match")
        self.tree.addBranch('sim_mpdgId',            'Int_t', 0, "PDG id of SIM mother")

    def _fill_tree(self, i, data):
        self.tree.reset()

        ## event info
        self.tree['run'] = self.event.run
        self.tree['ls']  = self.event.luminosityBlock
        self.tree['evt'] = self.event.event

        self.tree['pt']     = self.event.Muon_pt[i]
        self.tree['eta']    = self.event.Muon_eta[i]
        self.tree['phi']    = self.event.Muon_phi[i]

        self.tree['bmmBaseId']  = 0
        if self.event.Muon_isTracker[i] and self.event.Muon_isTracker[i] and \
           self.event.MuonId_highPurity[i] and self.event.Muon_pt[i] > 4.0 and \
           abs(self.event.Muon_eta[i]) < 1.4:
            self.tree['bmmBaseId']  = 1
        
        self.tree['looseId']    = self.event.Muon_looseId[i]
        self.tree['mediumId']   = self.event.Muon_mediumId[i]
        self.tree['tightId']    = self.event.Muon_tightId[i]
        self.tree['softId']     = self.event.Muon_softId[i]
        self.tree['softMvaId']  = self.event.Muon_softMvaId[i]

        self.tree['softMva']  = self.event.Muon_softMva[i]
                  
        self.tree['trkKink']             = self.event.MuonId_trkKink[i]
        self.tree['glbTrackProbability'] = self.event.MuonId_glbTrackProbability[i]
        self.tree['chi2LocalPosition']   = self.event.MuonId_chi2LocalPosition[i]
        self.tree['glbNormChi2']         = self.event.MuonId_glbNormChi2[i]
        self.tree['trkValidFrac']        = self.event.MuonId_trkValidFrac[i]
        self.tree['match1_dX']           = self.event.MuonId_match1_dX[i]
        self.tree['match1_pullX']        = self.event.MuonId_match1_pullX[i]
        self.tree['match1_pullDxDz']     = self.event.MuonId_match1_pullDxDz[i]
        self.tree['match1_dY']           = self.event.MuonId_match1_dY[i]
        self.tree['match1_pullY']        = self.event.MuonId_match1_pullY[i]
        self.tree['match1_pullDyDz']     = self.event.MuonId_match1_pullDyDz[i]
        self.tree['match2_dX']           = self.event.MuonId_match2_dX[i]
        self.tree['match2_pullX']        = self.event.MuonId_match2_pullX[i]
        self.tree['match2_pullDxDz']     = self.event.MuonId_match2_pullDxDz[i]
        self.tree['match2_dY']           = self.event.MuonId_match2_dY[i]
        self.tree['match2_pullY']        = self.event.MuonId_match2_pullY[i]
        self.tree['match2_pullDyDz']     = self.event.MuonId_match2_pullDyDz[i]

        self.tree['nPixels']             = self.event.MuonId_nPixels[i]
        self.tree['nValidHits']          = self.event.MuonId_nValidHits[i]
        self.tree['nLostHitsInner']      = self.event.MuonId_nLostHitsInner[i]
        self.tree['nLostHitsOn']         = self.event.MuonId_nLostHitsOn[i]
        self.tree['nLostHitsOuter']      = self.event.MuonId_nLostHitsOuter[i]

        self.tree['trkLayers']           = self.event.MuonId_trkLayers[i]
        self.tree['trkLostLayersInner']  = self.event.MuonId_trkLostLayersInner[i]
        self.tree['trkLostLayersOn']     = self.event.MuonId_trkLostLayersOn[i]
        self.tree['trkLostLayersOuter']  = self.event.MuonId_trkLostLayersOuter[i]

        self.tree['highPurity']          = self.event.MuonId_highPurity[i]

        for key, value in data.items():
            self.tree[key] = value

        self.tree.fill()

if __name__ == "__main__":

    job = {
        "input": [
            "/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/BsToMuMu_bmm_fakes_and_ids.root"
        ],
        "tree_name" : "muons",
      }  

    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    p = FlatNtupleForMuonMVA(file_name)
    p.limit = 1000
        
    p.process()
