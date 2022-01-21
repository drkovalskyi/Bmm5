from PostProcessingBase import FlatNtupleBase

import os, re, sys, time, subprocess, json
from math import *
import multiprocessing
from datetime import datetime
import hashlib

def deltaPhi(phi1, phi2):
    return acos(cos(phi2 - phi1))

def deltaR(eta1, phi1, eta2, phi2):
    return sqrt(pow(deltaPhi(phi1, phi2), 2) + pow(eta2 - eta1, 2))

class FlatNtupleForTrigEfficiency(FlatNtupleBase):
    """Flat ROOT ntuple producer for muon trigger object reconstruction efficiency study"""

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
        'L1_Mu6_DoubleEG15er2p5',
        'HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60',
        'HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5',
        'HLT_Mu15_IsoVVVL_PFHT450_PFMET50',
        'HLT_Mu15_IsoVVVL_PFHT450',
        'HLT_Mu15_IsoVVVL_PFHT600',
        'HLT_IsoMu24',
        'HLT_IsoMu27',
        'HLT_IsoMu20',
        'HLT_Mu8',
        'HLT_Mu15',
        'HLT_Mu17',
    ]
    
    def _validate_inputs(self):
        """Task specific input validation"""

        # check for missing information
        for parameter in ['input', 'tag_triggers', 'cut', 'tree_name', 'require_muon_tag']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

    def _is_tag(self, i):
        if self.event.Muon_pt < 5: return False
        if not self.event.Muon_isTracker[i]: return False
        if not self.event.Muon_isGlobal[i]: return False

        for trigger in self.job_info['tag_triggers']:
            if hasattr(self.event, "MuonId_%s" % trigger):
                trig_info = getattr(self.event, "MuonId_%s" % trigger)
                if trig_info[i]: return True
        return False
        
            
    def _process_events(self):
        """Event loop"""

        parsed_cut = self.get_cut()
                    
        for event in self.input_tree:
            self.event = event           

            # check for firing triggers
            if len(self.job_info['tag_triggers']) > 0:
                fired = False
                for trigger in self.job_info['tag_triggers']:
                    if hasattr(event, trigger):
                        trig_info = getattr(event, trigger)
                        if trig_info:
                            fired = True
                            break
                if not fired: continue

            # find probe muons
            probes = []
            tag_prob_pairs = []
            itag = -1
            while itag < self.event.nMuon:
                if self.job_info['require_muon_tag']:
                    if itag < 0:
                        itag = 0
                    if not self._is_tag(itag):
                        continue
                for iprobe in range(self.event.nMuon):
                    if iprobe == itag:
                        continue
                    if iprobe in probes:
                        continue
                
                    # cut = self.job_info['cut'].format(cand=cand, tree="self.event")
                    format_dict = { 'nMuon' : iprobe, 'tree' : 'self.event'}
                    cut = parsed_cut.format(**format_dict)
                    if not eval(cut):
                        continue

                    probes.append(iprobe)
                    tag_prob_pairs.append((itag, iprobe))
                if itag < 0:
                    break
                else:
                    itag += 1

            # Store probes
            for tag, probe in tag_prob_pairs:
                self._fill_tree(tag, probe)

    def _configure_output_tree(self):

        self.tree.addBranch('run',               'UInt_t', 0)
        self.tree.addBranch('evt',            'ULong64_t', 0)
        self.tree.addBranch('ls',                'UInt_t', 0)
        self.tree.addBranch('PV_npvsGood',        'Int_t', 0)
        self.tree.addBranch('nl1',                'Int_t', 0, "Number of L1 objects with any quality")
        self.tree.addBranch('nl1_SQ',             'Int_t', 0, "Number of L1 objects with single quality")
                
        self.tree.addBranch('probe_pt',         'Float_t', 0, "probe pt")
        self.tree.addBranch('probe_eta',        'Float_t', 0, "probe eta")
        self.tree.addBranch('probe_phi',        'Float_t', 0, "probe phi")
        self.tree.addBranch('probe_index',       'UInt_t', 0, "probe index")
        self.tree.addBranch('probe_q',            'Int_t', 0, "probe charge")

        self.tree.addBranch('tag_pt',           'Float_t', 0, "tag pt")
        self.tree.addBranch('tag_eta',          'Float_t', 0, "tag eta")
        self.tree.addBranch('tag_phi',          'Float_t', 0, "tag phi")
        self.tree.addBranch('tag_index',         'UInt_t', 0, "tag index")
        self.tree.addBranch('tag_q',              'Int_t', 0, "tag charge")
        
        self.tree.addBranch('tag_probe_dr',     'Float_t', 0, "dR between tag and probe using offline momentum")
        self.tree.addBranch('tag_probe_l1_dr',  'Float_t', 0, "dR between tag and probe using L1 momentum")
        
        self.tree.addBranch('probe_hlt_pt',     'Float_t', 0, "probe HLT trigger object pt")
        self.tree.addBranch('probe_l1_pt',      'Float_t', 0, "probe L1 trigger object pt")
        self.tree.addBranch('probe_l1_eta',     'Float_t', 0, "probe L1 trigger object eta")
        self.tree.addBranch('probe_l1_phi',     'Float_t', 0, "probe L1 trigger object phi")
        self.tree.addBranch('probe_l1_quality',   'Int_t', 0, "probe L1 trigger object quality")

        self.tree.addBranch('tag_hlt_pt',       'Float_t', 0, "tag HLT trigger object pt")
        self.tree.addBranch('tag_l1_pt',        'Float_t', 0, "tag L1 trigger object pt")
        self.tree.addBranch('tag_l1_eta',       'Float_t', 0, "tag L1 trigger object eta")
        self.tree.addBranch('tag_l1_phi',       'Float_t', 0, "tag L1 trigger object phi")
        self.tree.addBranch('tag_l1_quality',     'Int_t', 0, "tag L1 trigger object quality")

        for trigger in self.job_info['tag_triggers']:
            if trigger not in self.triggers:
                self.triggers.append(trigger)
                
        for trigger in self.triggers:
            self.tree.addBranch(trigger, 'UInt_t', 0)
            self.tree.addBranch("%s_ps" % trigger, 'Int_t', -1, "Prescale. 0 - Off, -1 - no information")
            
    def _fill_tree(self, tag, probe):
        self.tree.reset()

        ## event info
        self.tree['run']              = self.event.run
        self.tree['ls']               = self.event.luminosityBlock
        self.tree['evt']              = self.event.event
        self.tree['PV_npvsGood']      = self.event.PV_npvsGood
        
        nl1 = 0
        nl1_SQ = 0
        for index, l1_pt in enumerate(self.event.MuonId_l1_pt):
            if l1_pt >= 0:
                nl1 += 1
                if self.event.MuonId_l1_quality[index] == 12:
                    nl1_SQ += 1
                    
        self.tree['nl1']              = nl1
        self.tree['nl1_SQ']           = nl1_SQ

        self.tree['probe_pt']         = self.event.Muon_pt[probe]
        self.tree['probe_eta']        = self.event.Muon_eta[probe]
        self.tree['probe_phi']        = self.event.Muon_phi[probe]
        self.tree['probe_index']      = probe
        self.tree['probe_q']          = self.event.Muon_charge[probe]

        if tag >=0:
            self.tree['tag_pt']           = self.event.Muon_pt[tag]
            self.tree['tag_eta']          = self.event.Muon_eta[tag]
            self.tree['tag_phi']          = self.event.Muon_phi[tag]
            self.tree['tag_index']        = tag
            self.tree['tag_q']            = self.event.Muon_charge[tag]

            self.tree['tag_probe_dr']     = deltaR(self.event.Muon_eta[tag],
                                                  self.event.Muon_phi[tag],
                                                  self.event.Muon_eta[probe],
                                                  self.event.Muon_phi[probe])

            self.tree['tag_probe_l1_dr']  = deltaR(self.event.MuonId_l1_eta[tag],
                                                  self.event.MuonId_l1_phi[tag],
                                                  self.event.MuonId_l1_eta[probe],
                                                  self.event.MuonId_l1_phi[probe])
        
            self.tree['tag_hlt_pt']       = self.event.MuonId_hlt_pt[tag]
            self.tree['tag_l1_pt']        = self.event.MuonId_l1_pt[tag]
            self.tree['tag_l1_eta']       = self.event.MuonId_l1_eta[tag]
            self.tree['tag_l1_phi']       = self.event.MuonId_l1_phi[tag]
            self.tree['tag_l1_quality']   = self.event.MuonId_l1_quality[tag]
        
        self.tree['probe_hlt_pt']     = self.event.MuonId_hlt_pt[probe]
        self.tree['probe_l1_pt']      = self.event.MuonId_l1_pt[probe]
        self.tree['probe_l1_eta']     = self.event.MuonId_l1_eta[probe]
        self.tree['probe_l1_phi']     = self.event.MuonId_l1_phi[probe]
        self.tree['probe_l1_quality'] = self.event.MuonId_l1_quality[probe]

        for trigger in self.triggers:
            if hasattr(self.event, trigger):
                self.tree[trigger] = getattr(self.event, trigger)
            if hasattr(self.event, "prescale_" + trigger):
                self.tree[trigger + "_ps"] = getattr(self.event, "prescale_" + trigger)
                
        self.tree.fill()

if __name__ == "__main__":

    ### create a test job

    # prefix = "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/SingleMuon+Run2018D-12Nov2019_UL2018-v8+MINIAOD/"
    # prefix = "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/SingleMuon+Run2017F-09Aug2019_UL2017-v1+MINIAOD//"
    prefix = "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/518/EGamma+Run2018D-12Nov2019_UL2018-v4+MINIAOD/"
    job = {
        "input": [
            # SingleMuon+Run2018D
            # prefix + "3E0EF47D-4FF1-A440-81BC-064B0F05D611.root",
            # prefix + "00006384-4AE6-8F4C-9E59-3C78C2461244.root",
            # prefix + "000A7806-F703-9549-9021-424A2553EF38.root",
            # prefix + "001F4E39-30E2-8645-AC79-8345B516DAA3.root",
            # prefix + "0023E28F-4442-B345-855E-95D9E78FFFFD.root",
            # prefix + "003D83E3-BE6E-D84B-9C6E-BCE65CAFE110.root",
            # prefix + "003E9885-071B-384E-B099-39A3BCA96B87.root",
            # prefix + "00530ADF-CCD1-584F-AB96-4B564E5BBAC5.root",
            # prefix + "0053ADEA-081D-8E40-AD76-0BD96DAD2BBB.root",
            # prefix + "0063B8CC-09EA-A543-B920-32245C701DB4.root",

            # SingleMuon+Run2017F
            # prefix + "0014A39D-5620-F143-8EB5-51F760200B4F.root",
            # prefix + "00256B90-A7CE-864B-8A17-0749AAB5F73D.root",
            # prefix + "003CFCAB-0D2F-114A-8099-D1D2C9761E63.root",
            # prefix + "006AC703-3C92-9347-AEEE-1FA6BF3D08EA.root",
            # prefix + "0074A30E-2BD6-8643-B768-41D89B7868EF.root",
            # prefix + "00C0A3BD-11BC-1146-B3D2-E0FA2EEC7A63.root",
            # prefix + "00C641C3-E36F-CC4F-A316-310891E816BB.root",
            # prefix + "00DB73CB-10A5-E748-BD63-EE5D286F70AC.root",
            # prefix + "00E0CBDF-7A40-EC47-862D-42D035B57C34.root",
            # prefix + "00E44F96-2FA5-8644-A698-97F3BD588315.root",

            prefix + "6749DD07-E613-3E44-A0AC-B590CC6C2FB5.root",
            prefix + "0005843C-A8E4-BA44-8791-FB30945D33D4.root",
            prefix + "00141923-3629-944A-88F5-FC5B61374B67.root",
            prefix + "00159164-6FFC-8040-A9CB-704AA15CE165.root",
            prefix + "0018D05F-E404-8347-9BB3-BB5DA36CEAA3.root",
            prefix + "00194ADB-8C7E-E64D-99CB-5A8A3CBA691F.root",
            prefix + "001A4023-03FC-2843-80D9-7CC4A28B7EC6.root",
            prefix + "00254B67-252A-FE4C-B5C3-CBED8C6937A0.root",
            prefix + "002AC038-9B89-6D4B-B034-676C93DAD974.root",
            prefix + "002D84A3-DAD8-6E4B-827C-9CFE2AFC1E13.root",
            prefix + "002D9E80-76EB-8B4E-819A-E29B603503AC.root",
        ],
        "tree_name" : "muons",
        "require_muon_tag" : False, 
        # "tag_triggers" : ["HLT_IsoMu24", "HLT_IsoMu27"],
        "tag_triggers" : ["L1_Mu6_DoubleEG10er2p5", "L1_Mu6_DoubleEG17er2p5",
                          "L1_Mu6_DoubleEG12er2p5", "L1_Mu6_DoubleEG15er2p5",
                          "L1_SingleMu3", "L1_SingleMu5", "L1_SingleMu7",
                          "L1_Mu6_DoubleEG10", "L1_Mu6_DoubleEG17"
                          ],
        "cut" : "Muon_isTracker && Muon_isGlobal && Muon_softMva > 0.45 && Muon_pt > 4.0 && abs(Muon_eta) < 1.5 && Muon_pt<20"
      }  

    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = FlatNtupleForMLFit("/tmp/dmytro/9081fd38604b24f4c6035628a89e18ed.job")
    p = FlatNtupleForTrigEfficiency(file_name)
    print p.__dict__
        
    p.process()
