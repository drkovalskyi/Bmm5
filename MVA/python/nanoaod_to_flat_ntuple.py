#!/bin/env python
#
# Produce flat ROOT ntuples for MVA training
#

import os, re, sys, time, subprocess
from Bmm5.MVA.mtree import MTree
import multiprocessing
from datetime import datetime
import hashlib
from ROOT import TFile, TTree

##
## User Input
##
blind_signal_region = True
blinding_right_sideband = True
blinding_left_sideband = False
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

def _configure_output_tree(tree):
    ## event info
    tree.addBranch( 'evt_run',   'UInt_t', 0)
    tree.addBranch( 'evt_lumi',  'UInt_t', 0)
    tree.addBranch( 'evt_event', 'ULong64_t', 0)
    # tree.addBranch( 'evt_pu',    'UInt_t', 0)
    # tree.addBranch( 'evt_nvtx',  'UInt_t', 0)

    ## mm info
    tree.addBranch('mm_kin_mass',   'Float_t', 0, "Vertex constrained dimuon mass")
    tree.addBranch('mm_kin_eta',    'Float_t', 0, "Vertex constrained dimuon mass")
    tree.addBranch('mm_mu1_pt',     'Float_t', 0)
    tree.addBranch('mm_mu2_pt',     'Float_t', 0)
    tree.addBranch('mm_mu1_eta',    'Float_t', 0)
    tree.addBranch('mm_mu2_eta',    'Float_t', 0)
    tree.addBranch('mm_kin_l3d',    'Float_t', 0, "Decay length wrt Primary Vertex in 3D")
    tree.addBranch('mm_kin_sl3d',   'Float_t', 0, "Decay length significance wrt Primary Vertex in 3D")
    tree.addBranch('mm_kin_alpha',  'Float_t', 0, "Cosine of ointing angle in 3D wrt PV")
    tree.addBranch('mm_kin_spvip',  'Float_t', 0, "Significance of impact parameter wrt Primary Vertex in 3D")
    tree.addBranch('mm_iso',        'Float_t', 0, "B isolation the way it's done in Bmm4")
    tree.addBranch('mm_kin_vtx_chi2dof', 'Float_t', 0, "")
    tree.addBranch('mm_docatrk',    'Float_t', 0, "Distance of closest approach of a track to the vertex")
    tree.addBranch('mm_closetrk',   'UInt_t',  0, "Number of tracks compatible with the vertex by DOCA (0.3mm)")
    tree.addBranch('mm_m1iso',      'Float_t', 0, "Muon isolation the way it's done in Bmm4")
    tree.addBranch('mm_m2iso',      'Float_t', 0, "Muon isolation the way it's done in Bmm4")
    tree.addBranch('mm_os',         'UInt_t',  0, "Opposite charge or not")

    tree.addBranch('mm_kin_alphaXY','Float_t', 0, "Cosine of pointing angle in XY wrt BS")
    tree.addBranch('mm_nBMTrks',     'UInt_t', 0, "Number of tracks more compatible with the mm vertex than with PV by doca significance")
    tree.addBranch('mm_nDisTrks',    'UInt_t', 0, "Number of displaced tracks compatible with the vertex by vertex probability")
    tree.addBranch('mm_closetrks1',  'UInt_t', 0, "Number of tracks compatible with the vertex by DOCA(0.3mm) and signifance less than 1")
    tree.addBranch('mm_closetrks2',  'UInt_t', 0, "Number of tracks compatible with the vertex by DOCA(0.3mm) and signifance less than 2")
    tree.addBranch('mm_closetrks3',  'UInt_t', 0, "Number of tracks compatible with the vertex by DOCA(0.3mm) and signifance less than 3")
    
    
    
    ## gen info
    tree.addBranch('mm_gen_pdgId', 'Int_t',  0, "Gen match: dimuon pdg Id")
    
    ## BDT info
    tree.addBranch('mm_bdt',      'Float_t', 0, "Bmm4 BDT")

def _good_muon(event,index):
    if not (abs(event.Muon_eta[index])<1.4): return False
    if not (event.Muon_pt[index]>4): return False
    if not (event.Muon_softId[index]): return False
    return True

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

def _apply_selection(event):
    candidates = []
    if not event.nmm>0: return candidates
    for i in range(event.nmm):
        selection = []
        selection.append( ('good_muon1',_good_muon(event,event.mm_mu1_index[i])) )
        selection.append( ('good_muon2',_good_muon(event,event.mm_mu2_index[i])) )
        selection.append( ('mass_cut',  abs(event.mm_kin_mass[i]-5.4)<0.5) )
        selection.append( ('valid_fit', event.mm_kin_valid[i]) )
        selection.append( ('decay_length_significance', event.mm_kin_sl3d[i]>4) )
        selection.append( ('vertex_chi2', event.mm_kin_vtx_chi2dof[i]<5) )
        keeper = True
        if blind_signal_region:
            keeper = False
            if blinding_right_sideband:
                if (event.mm_kin_mass[i]-5.3)>0.2: keeper = True
            if blinding_left_sideband:
                if (event.mm_kin_mass[i]-5.3)<-0.2: keeper = True
        selection.append( ('blinding', keeper) )
        candidates.append(selection)
    return candidates

def _good_candidates(candidates):
    good_candidates = []
    for i in range(len(candidates)):
        good_candidate = True
        for selection in candidates[i]:
            if not selection[1]:
                good_candidate = False
                break
        if good_candidate: good_candidates.append(i)
    return good_candidates

def _fill_tree(tree,event,cand):
    tree.reset()

    ## event info
    tree['evt_run']   = event.run
    tree['evt_lumi']  = event.luminosityBlock
    tree['evt_event'] = event.event

    ## mm info
    tree['mm_kin_mass']  = event.mm_kin_mass[cand]
    tree['mm_kin_eta']   = event.mm_kin_eta[cand]
    tree['mm_mu1_pt']    = event.Muon_pt[event.mm_mu1_index[cand]]
    tree['mm_mu2_pt']    = event.Muon_pt[event.mm_mu2_index[cand]]
    tree['mm_mu1_eta']   = event.Muon_eta[event.mm_mu1_index[cand]]
    tree['mm_mu2_eta']   = event.Muon_eta[event.mm_mu1_index[cand]]
    tree['mm_kin_l3d']   = event.mm_kin_l3d[cand]
    tree['mm_kin_sl3d']  = event.mm_kin_sl3d[cand]
    tree['mm_kin_alpha'] = event.mm_kin_alpha[cand]
    tree['mm_kin_spvip'] = event.mm_kin_pvip[cand]/event.mm_kin_pvipErr[cand]
    tree['mm_iso']       = event.mm_iso[cand]
    tree['mm_kin_vtx_chi2dof'] = event.mm_kin_vtx_chi2dof[cand]
    tree['mm_docatrk']   = event.mm_docatrk[cand]
    tree['mm_closetrk']  = event.mm_closetrk[cand]
    tree['mm_closetrks1']  = event.mm_closetrks1[cand]
    tree['mm_closetrks2']  = event.mm_closetrks2[cand]
    tree['mm_closetrks3']  = event.mm_closetrks3[cand]
    tree['mm_m1iso']     = event.mm_m1iso[cand]
    tree['mm_m2iso']     = event.mm_m2iso[cand]
    tree['mm_os']        = 1 if (event.Muon_charge[event.mm_mu1_index[cand]] * event.Muon_charge[event.mm_mu2_index[cand]] == -1) else 0
    tree['mm_kin_alphaXY'] = event.mm_kin_cosAlphaXY[cand]
    tree['mm_nBMTrks']     = event.mm_nBMTrks[cand]
    tree['mm_nDisTrks']    = event.mm_nDisTrks[cand]
    
    ## gen info
    if hasattr(event, 'mm_gen_pdgId'):
        tree['mm_gen_pdgId'] = event.mm_gen_pdgId[cand]
    else:
        tree['mm_gen_pdgId'] = 0
        

    ## old BDT
    tree['mm_bdt']       = event.mm_bdt[cand]

    tree.fill()

def _analyze_selection(candidates,statistics):
    if not collect_statistics: return
    max_selections = 0
    best_cand = None
    for i in range(len(candidates)):
        n = 0
        for (name,passed) in candidates[i]:
            if passed:
                n+=1
            else:
                break
        if n > max_selections:
            best_cand = i
            max_selections = n
    if best_cand!=None:
        for (selection_name,passed) in candidates[best_cand]:
            if passed:
                if not selection_name in statistics:
                    statistics[selection_name] = 1
                else:
                    statistics[selection_name] += 1
            else:
                break

def process_file(input_file,output_filename=None,signal_only=False):
    _formated_print("Processing file: %s" % input_file)
    statisitics = dict()
    if not output_filename:
        output_filename = "/tmp/%s.root" % hashlib.md5(input_file).hexdigest()
    if output_is_already_available(output_filename):
        return output_filename
    try:
        fout = TFile( output_filename, 'recreate' )
        output_tree = MTree( tree_name, '' )
        _configure_output_tree(output_tree)

        fin = TFile( input_file )
        input_tree = fin.Get("Events")
        nevents = 0

        for event in input_tree:
            # if event.eventAuxiliary().event()!=69338226: continue
            if n_events_limit and nevents>=n_events_limit: break
            # print "Event: %d \tLumi: %d \tRun: %d"% (event.event,
            #                                          event.luminosityBlock,
            #                                          event.run)
            if n_events_limit==None and (nevents)%10000==0: print "[%d] Event: %d" % (multiprocessing.current_process().pid,
                                                                                      nevents)
            nevents += 1
            info = _apply_selection(event)
            _analyze_selection(info,statisitics)
            for cand in _good_candidates(info):
                if signal_only and event.mm_gen_pdgId[cand]==0: continue
                _fill_tree(output_tree,event,cand)
         
        fout.Write()
        fout.Close()
        if collect_statistics:
            print statisitics
        return output_filename
    except:
        print "[%d-%s] process failed" % (multiprocessing.current_process().pid,_get_time())
        raise 
    return None

def process_file_mc_mathed(input_file):
    return process_file(input_file,None,True)

if __name__ == "__main__":
    # print process_file('/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/505/BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0D690850-542C-874D-8553-2454517F6E54.root')
    # print process_file('/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/505/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM_Skim.root')
    print process_file('/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/505/Charmonium+Run2018A-17Sep2018-v1+MINIAOD/0024D10F-B9D6-3E46-BEE1-413765D77D91.root')

