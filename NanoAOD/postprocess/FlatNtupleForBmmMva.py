from PostProcessingBase import Processor

import os, re, sys, time, subprocess, math
from Bmm5.MVA.mtree import MTree
import multiprocessing
from datetime import datetime
import hashlib
from ROOT import TFile, TTree

class FlatNtupleForBmmMva(Processor):
    """Produce flat ROOT ntuples for Bmm5 MVA training"""
    def _process(self):
        # check for missing information
        for parameter in ['input', 'signal_only', 'tree_name', 'blind']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

        # process files
        results = []
        t0 = time.clock()
        n_events = 0
        for f in self.job_info['input']:
            result, n = self.process_file(f)
            n_events += n
            results.append(result)
        print n//(time.clock()-t0), "Hz"
        
        print "Merging output."
        # merge results
        good_files = []
        for rfile in results:
            if rfile: good_files.append(rfile)
        command = "hadd -f %s %s" % (self.job_output_tmp, " ".join(good_files))
        status = subprocess.call(command, shell=True)

        if status==0:
            print "Merged output."
            for file in good_files:
                os.remove(file)
        else:
            raise Exception("Merge failed")

    def process_file(self, input_file):
        print "Processing file: %s" % input_file
        match = re.search("([^\/]+)\.root$", input_file)
        if match:
            output_filename = "%s/%s_processed.root" % (self.tmp_dir, match.group(1))
        else:
            raise Exception("Unexpected input ROOT file name:\n%s" % input_file)
        
        statisitics = dict()

        fout = TFile(output_filename, 'recreate')
        self.tree = MTree(self.job_info['tree_name'], '')
        self._configure_output_tree()

        fin = TFile(input_file)
        input_tree = fin.Get("Events")
        nevents = 0

        for event in input_tree:
            # if event.eventAuxiliary().event()!=69338226: continue
            # print "Event: %d \tLumi: %d \tRun: %d"% (event.event,
            #                                          event.luminosityBlock,
            #                                          event.run)
            nevents += 1
            self.event = event
            candidates = self._select_candidates()
            # _analyze_selection(info,statisitics)
            for cand in self._good_candidates(candidates):
                if self.job_info['signal_only'] and self.event.mm_gen_pdgId[cand]==0: continue
                self._fill_tree(cand)

        nout = self.tree.tree.GetEntries()
        print 'Selected %d / %d entries from %s (%.2f%%)' % (nout, nevents, input_file, 100.*nout/nevents if nevents else 0)
        fout.Write()
        fout.Close()
        # if collect_statistics:
        #    print statisitics
        return output_filename, nevents

    def _select_candidates(self):
        candidates = []
        if not self.event.nmm>0: return candidates
        for i in range(self.event.nmm):
            selection = []
            selection.append( ('good_muon1', self._good_muon(self.event.mm_mu1_index[i])) )
            selection.append( ('good_muon2', self._good_muon(self.event.mm_mu2_index[i])) )
            selection.append( ('mass_cut',  abs(self.event.mm_kin_mass[i]-5.4)<0.5) )
            selection.append( ('valid_fit', self.event.mm_kin_valid[i]) )
            selection.append( ('decay_length_significance', self.event.mm_kin_sl3d[i]>4) )
            selection.append( ('vertex_chi2', self.event.mm_kin_vtx_chi2dof[i]<5) )
            keeper = True
            if not hasattr(self.event, 'mm_gen_pdgId') and blind_signal_region:
                keeper = False
                if keep_right_sideband:
                    if (self.event.mm_kin_mass[i]-5.3)>0.2: keeper = True
                if keep_left_sideband:
                    if (self.event.mm_kin_mass[i]-5.3)<-0.2: keeper = True
            selection.append( ('blinding', keeper) )
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
        # self.tree.addBranch( 'evt_pu',    'UInt_t', 0)
        self.tree.addBranch( 'evt_nvtx',          'UInt_t',  0, "Number of good reconstructed primary vertices")
        self.tree.addBranch( 'evt_met',           'Float_t', 0, "Missing transverse energy")

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
        self.tree.addBranch('mm_kin_spvip',       'Float_t', 0, "Significance of impact parameter wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_kin_spvlip',      'Float_t', 0, "Significance of longitudinal impact parameter wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_iso',             'Float_t', 0, "B isolation the way it's done in Bmm4")
        self.tree.addBranch('mm_kin_vtx_chi2dof', 'Float_t', 0, "Normalized chi2 of the dimuon vertex kinematic fit")
        self.tree.addBranch('mm_doca',            'Float_t', 0, "dimuon DOCA")
        self.tree.addBranch('mm_docatrk',         'Float_t', 0, "Distance of closest approach of a track to the vertex")
        self.tree.addBranch('mm_closetrk',        'UInt_t',  0, "Number of tracks compatible with the vertex by DOCA (0.3mm)")
        self.tree.addBranch('mm_m1iso',           'Float_t', 0, "Muon isolation the way it's done in Bmm4")
        self.tree.addBranch('mm_m2iso',           'Float_t', 0, "Muon isolation the way it's done in Bmm4")
        self.tree.addBranch('mm_os',              'UInt_t',  0, "Opposite charge or not")

        self.tree.addBranch('mm_kin_alphaXY',     'Float_t', 0, "Cosine of pointing angle in XY wrt BS")
        self.tree.addBranch('mm_nBMTrks',         'UInt_t',  0, "Number of tracks more compatible with the mm vertex than with PV by doca significance")
        self.tree.addBranch('mm_nDisTrks',        'UInt_t',  0, "Number of displaced tracks compatible with the vertex by vertex probability")
        self.tree.addBranch('mm_closetrks1',      'UInt_t',  0, "Number of tracks compatible with the vertex by DOCA(0.3mm) and signifance less than 1")
        self.tree.addBranch('mm_closetrks2',      'UInt_t',  0, "Number of tracks compatible with the vertex by DOCA(0.3mm) and signifance less than 2")
        self.tree.addBranch('mm_closetrks3',      'UInt_t',  0, "Number of tracks compatible with the vertex by DOCA(0.3mm) and signifance less than 3")
        self.tree.addBranch('mm_otherVtxMaxProb', 'Float_t', 0, "Max vertexing probability of one of the muons with a random track with minPt=0.5GeV")
        self.tree.addBranch('mm_otherVtxMaxProb1','Float_t', 0, "Max vertexing probability of one of the muons with a random track with minPt=1.0GeV")
        self.tree.addBranch('mm_otherVtxMaxProb2','Float_t', 0, "Max vertexing probability of one of the muons with a random track with minPt=2.0GeV")

        self.tree.addBranch('HLT_DoubleMu4_3_Bs', 'UInt_t',  0, "Main analysis trigger")
        self.tree.addBranch('mm_kin_pvip',        'Float_t', 0, "Impact parameter wrt Primary Vertex in 3D")
        self.tree.addBranch('mm_kin_pvlip',       'Float_t', 0, "Longitudinal impact parameter wrt Primary Vertex")
        self.tree.addBranch('mm_met',             'Float_t', 0, "MET projected on dimuon direction")

        ## gen info
        self.tree.addBranch('mm_gen_pdgId',       'Int_t',   0, "Gen match: dimuon pdg Id")

        ## BDT info
        self.tree.addBranch('mm_bdt',             'Float_t', 0, "Bmm4 BDT")

    def _fill_tree(self, cand):
        self.tree.reset()

        ## event info
        self.tree['evt_run']   = self.event.run
        self.tree['evt_lumi']  = self.event.luminosityBlock
        self.tree['evt_event'] = self.event.event
        self.tree['evt_nvtx']  = self.event.PV_npvsGood
        self.tree['evt_met']   = self.event.MET_pt

        ## mm info
        self.tree['mm_kin_mass']    = self.event.mm_kin_mass[cand]
        self.tree['mm_kin_eta']     = self.event.mm_kin_eta[cand]
        self.tree['mm_kin_pt']      = self.event.mm_kin_pt[cand]
        self.tree['mm_mu1_pt']      = self.event.Muon_pt[self.event.mm_mu1_index[cand]]
        self.tree['mm_mu2_pt']      = self.event.Muon_pt[self.event.mm_mu2_index[cand]]
        self.tree['mm_mu1_eta']     = self.event.Muon_eta[self.event.mm_mu1_index[cand]]
        self.tree['mm_mu2_eta']     = self.event.Muon_eta[self.event.mm_mu2_index[cand]]
        self.tree['mm_mu1_phi']     = self.event.Muon_phi[self.event.mm_mu1_index[cand]]
        self.tree['mm_mu2_phi']     = self.event.Muon_phi[self.event.mm_mu2_index[cand]]
        self.tree['mm_kin_l3d']     = self.event.mm_kin_l3d[cand]
        self.tree['mm_kin_sl3d']    = self.event.mm_kin_sl3d[cand]
        self.tree['mm_kin_slxy']    = self.event.mm_kin_slxy[cand]
        self.tree['mm_kin_alpha']   = self.event.mm_kin_alpha[cand]
        self.tree['mm_kin_spvip']   = self.event.mm_kin_pvip[cand]/self.event.mm_kin_pvipErr[cand]
        self.tree['mm_kin_spvlip']  = self.event.mm_kin_pvlip[cand]/self.event.mm_kin_pvlipErr[cand]
        self.tree['mm_kin_pvip']    = self.event.mm_kin_pvip[cand]
        self.tree['mm_kin_pvlip']   = self.event.mm_kin_pvlip[cand]
        self.tree['mm_iso']         = self.event.mm_iso[cand]
        self.tree['mm_kin_vtx_chi2dof'] = self.event.mm_kin_vtx_chi2dof[cand]
        self.tree['mm_docatrk']     = self.event.mm_docatrk[cand]
        self.tree['mm_closetrk']    = self.event.mm_closetrk[cand]
        self.tree['mm_closetrks1']  = self.event.mm_closetrks1[cand]
        self.tree['mm_closetrks2']  = self.event.mm_closetrks2[cand]
        self.tree['mm_closetrks3']  = self.event.mm_closetrks3[cand]
        self.tree['mm_otherVtxMaxProb']  = self.event.mm_otherVtxMaxProb[cand]
        self.tree['mm_otherVtxMaxProb1'] = self.event.mm_otherVtxMaxProb1[cand]
        self.tree['mm_otherVtxMaxProb2'] = self.event.mm_otherVtxMaxProb2[cand]
        self.tree['mm_m1iso']       = self.event.mm_m1iso[cand]
        self.tree['mm_m2iso']       = self.event.mm_m2iso[cand]
        self.tree['mm_os'] = 1 if (self.event.Muon_charge[self.event.mm_mu1_index[cand]] * self.event.Muon_charge[self.event.mm_mu2_index[cand]] == -1) else 0
        self.tree['mm_kin_alphaXY'] = self.event.mm_kin_cosAlphaXY[cand]
        self.tree['mm_nBMTrks']     = self.event.mm_nBMTrks[cand]
        self.tree['mm_nDisTrks']    = self.event.mm_nDisTrks[cand]
        self.tree['mm_doca']        = self.event.mm_doca[cand]
        self.tree['mm_met']         = self.event.MET_pt*math.cos(self.event.mm_kin_phi[cand]-self.event.MET_phi)

        ## gen info
        if hasattr(self.event, 'mm_gen_pdgId'):
            self.tree['mm_gen_pdgId'] = self.event.mm_gen_pdgId[cand]
        else:
            self.tree['mm_gen_pdgId'] = 0

        self.tree['HLT_DoubleMu4_3_Bs'] = self.event.HLT_DoubleMu4_3_Bs

        ## old BDT
        self.tree['mm_bdt']       = self.event.mm_bdt[cand]

        self.tree.fill()

    def _good_muon(self, index):
        if index < 0: return False
        if not (abs(self.event.Muon_eta[index]) < 1.4): return False
        if not (self.event.Muon_pt[index] > 4): return False
        if not self.event.Muon_softId[index]: return False
        return True


    ##
## User Input
##
blind_signal_region = True
keep_right_sideband = True
keep_left_sideband = True
tree_name = "mva"
force_recreation = False
n_events_limit = None 
# n_events_limit = 5000
collect_statistics = True

def _get_time():
    return datetime.now().strftime("%H:%M:%S")

def _formated_print(data):
    print "[%d-%s] %s" % (multiprocessing.current_process().pid,
                          _get_time(),data)
    sys.stdout.flush()

def _is_good_file(f):
    tf = TFile.Open(f)
    if not tf: return False
    if tf.IsZombie() or tf.TestBit(TFile.kRecovered): return False
    if not tf.Get(tree_name): return False
    return True

def output_is_already_available(filename):
    if os.path.exists(filename):
        good = _is_good_file(filename)
        message = "Output file %s already exists." % filename
        if force_recreation:
            os.remove(filename)
        elif good:
            print message, " It appears to be good. Skip it. If you want to recreate it, please remove the existing file" 
            return True
        else:
            print message," File is corrupted. Will recreate."
            os.remove(filename)
    return False


# def _select_candidates(event):
#     good_candidates = []
#     if not event.nmm>0: return good_candidates
#     for i in range(event.nmm):
#         if not _good_muon(event,event.mm_mu1_index[i]): continue
#         if not _good_muon(event,event.mm_mu2_index[i]): continue
#         if not event.mm_kin_valid[i]: continue
#         if not event.mm_kin_sl3d[i]>4: continue
#         if abs(event.mm_kin_mass[i]-5.4)>0.5: continue
#         if event.mm_kin_vtx_chi2dof[i]>5: continue 
#         if blind_signal_region and abs(event.mm_kin_mass[i]-5.3)<0.2: continue
#         good_candidates.append(i)
#     return good_candidates


# def _analyze_selection(candidates,statistics):
#     if not collect_statistics: return
#     max_selections = 0
#     best_cand = None
#     for i in range(len(candidates)):
#         n = 0
#         for (name,passed) in candidates[i]:
#             if passed:
#                 n+=1
#             else:
#                 break
#         if n > max_selections:
#             best_cand = i
#             max_selections = n
#     if best_cand!=None:
#         for (selection_name,passed) in candidates[best_cand]:
#             if passed:
#                 if not selection_name in statistics:
#                     statistics[selection_name] = 1
#                 else:
#                     statistics[selection_name] += 1
#             else:
#                 break


def unit_test():
    p = FlatNtupleForBmmMva("/eos/cms/store/group/phys_muon/dmytro/tmp/mva_ntuple_test/ffcb48b1f605ca5a5f2614b94830e806.job")
    # print p.__dict__
    p.process()

    # check single file processing
    # p.job_output_tmp = '/tmp/'
    # p.job_info['tree_name'] = 'mva'
    # p.job_info['signal_only'] = 'True'
    # p.process_file('root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/510/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/02F7319D-D4D6-8340-B94F-2D882775B406.root')

if __name__ == "__main__":
    unit_test()
