#!/bin/env python
#
# Tool to produce training samples using local resources
#

from nanoaod_to_flat_ntuple import *

##
## User Input
##
# force_recreation = False
remote_read = False

# nProcesses = 32
nProcesses = 16
# require_trigger_match = True

# reduce memory consumption
files_per_pool = 100

# Input can be provided for each samples in the following formats:
# - 'file_with_lfns' - file name containing input logical file names
# - 'file_with_pfns' - file name containing input physical file names
# - 'lfns'           - list of input logical file names
# - 'pfns'           - list of input physical file names
samples = {
   'Charmonium+Run2016B_1':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016B-17Jul2018_ver1-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016B_2':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016C':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016C-17Jul2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016D':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016D-17Jul2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016E':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016E-17Jul2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016F':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016F-17Jul2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016G':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016G-17Jul2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2016H':{
       'file_with_pfns':'skim-lists/Charmonium+Run2016H-17Jul2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2017B':{
       'file_with_pfns':'skim-lists/Charmonium+Run2017B-31Mar2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2017C':{
       'file_with_pfns':'skim-lists/Charmonium+Run2017C-31Mar2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2017D':{
      'file_with_pfns':'skim-lists/Charmonium+Run2017D-31Mar2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2017E':{
       'file_with_pfns':'skim-lists/Charmonium+Run2017E-31Mar2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2017F':{
       'file_with_pfns':'skim-lists/Charmonium+Run2017F-31Mar2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2018A':{
       'file_with_pfns':'skim-lists/Charmonium+Run2018A-17Sep2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2018B':{
        'file_with_pfns':'skim-lists/Charmonium+Run2018B-17Sep2018-v1+MINIAOD.txt'
   },
    'Charmonium+Run2018C':{
       'file_with_pfns':'skim-lists/Charmonium+Run2018C-17Sep2018-v1+MINIAOD.txt'
   },
   'Charmonium+Run2018D':{
       'file_with_pfns':'skim-lists/Charmonium+Run2018D-PromptReco-v2+MINIAOD.txt'
   },
    # 'BsToMuMu':{
    #     'file_with_pfns':'skim-lists/BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM.txt',
    #     'match':True
    # },
    'BsToMuMu_BMuonFilter_RunIIAutumn18MiniAOD':{
        'file_with_pfns':'skim-lists/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt',
        'match':True
    },
    'BsToMuMu_BMuonFilter_RunIIFall17MiniAODv2':{
        'file_with_pfns':'skim-lists/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM.txt',
        'match':True
    },
    # 'QCD_HT100to200_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'skim-lists/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_HT50to100_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/MVA/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_HT100to200_BGenFilter_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/MVA/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_HT200to300_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/MVA/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_Pt-80to120_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'skim-lists/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-30to50_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'skim-lists/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-50to80_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'skim-lists/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-30to50_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'skim-lists/QCD_Pt_30to50_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    'QCD_Pt-30toInf_BmmGenFilter':{
        'file_with_pfns':'skim-lists/QCD_Pt-30toInf_BmmGenFilter-NanoAOD.txt'
    }
}

# import os, re, sys, time, ROOT, subprocess
# from ROOT import TFile, TTree
# from DataFormats.FWLite import Events, Handle
# import multiprocessing
# from datetime import datetime
# import tempfile
# from array import array
# from AnalysisTools.Studies.mtree import MTree
# import hashlib


def getPFNs(lfns):
    files = []
    for file in lfns:
        if not input_is_LFNs:
            files.append(file)
        else:
            if remote_read:
                files.append("root://cms-xrd-global.cern.ch/%s"%file)
            else:
                fullpath = "/eos/cms/" + file
                if os.path.exists(fullpath):
                    files.append(fullpath)
                else:
                    raise Exception("File not found: %s" % fullpath)
    return files

def chunks(input_list, n):
    """Yield successive n-sized chunks from input_list"""
    for i in range(0, len(input_list), n):
        yield input_list[i:i + n]

def parellel_process(input_files,output_file_name,match=False):
    file_chunks = list(chunks(input_files,files_per_pool))
    results = []
    for flist in file_chunks:
        pool = multiprocessing.Pool(nProcesses) 
        if match:
            results.extend(pool.map(process_file_mc_mathed, flist))
        else:
            results.extend(pool.map(process_file, flist))
        pool.close()
    print "Multiprocessing is done. Merging output."

    good_files = []
    for rfile in results:
        if rfile: good_files.append(rfile)
    status = subprocess.call("hadd -f %s %s" % (output_file_name," ".join(good_files)),shell=True)
    if status==0:
        print "Merged output."
        for file in good_files:
            os.remove(file)
    else:
        print "Merge failed"

for sample,info in samples.items():
    print "Processing %s" % sample
    output_filename = "%s.root"%sample
    if output_is_already_available(output_filename):
        continue
    input_files = []
    if 'file_with_pfns' in info:
        with open(info['file_with_pfns']) as f:
            for pfn in f:
                input_files.append(pfn.strip())
    elif 'file_with_lfns' in info:
        with open(info['file_with_pfns']) as f:
            for lfn in f:
                input_files.append(lfn.strip())
        input_files = getPFNs(input_files)
    elif 'pfns' in info:
        input_files = info['pfns']
    elif 'lfns' in info:
        input_files = getPFNs(info['pfns'])
    else:
        print "Input is not specified. Skip the sample"
        continue
    if len(input_files)==0:
        print "Nothing to process. Skip the sample"
        continue

    match = False
    if 'match' in info:
        match = info['match']
    parellel_process(input_files,output_filename,match)

