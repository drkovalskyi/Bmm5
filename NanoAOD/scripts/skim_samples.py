#!/usr/bin/env python
# 
# Skim and merge NanoAOD samples
#
# This script requires https://github.com/cms-nanoAOD/nanoAOD-tools
# 
import os, sys, pprint, multiprocessing, subprocess, re
from ROOT import TFile,TTree
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import hashlib

cut = "Muon_softId[mm_mu1_index]&&abs(Muon_eta[mm_mu1_index])<1.4&&Muon_pt[mm_mu1_index]>4&&Muon_softId[mm_mu2_index]&&abs(Muon_eta[mm_mu2_index])<1.4&&Muon_pt[mm_mu2_index]>4&&abs(mm_kin_mass-5.4)<0.5&&mm_kin_sl3d>4&&mm_kin_vtx_chi2dof<5"

output_dir = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/505"

nProcesses = 16

force_recreation = False

# access data by LFN via global redirector
use_xrootd = False

# Input can be provided for each samples in the following formats:
# - 'file_with_lfns' - file name containing input logical file names
# - 'file_with_pfns' - file name containing input physical file names
# - 'lfns'           - list of input logical file names
# - 'pfns'           - list of input physical file names
samples = {
   # 'Charmonium+Run2016B_1':{
   #     'file_with_pfns':'lists/Charmonium+Run2016B-17Jul2018_ver1-v1+MINIAOD.txt',
   #     'files_per_job':100
   # },
   # 'Charmonium+Run2016B_2':{
   #     'file_with_pfns':'lists/Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD.txt',
   #     'files_per_job':100
   # },
   # 'Charmonium+Run2016C':{
   #     'file_with_pfns':'lists/Charmonium+Run2016C-17Jul2018-v1+MINIAOD.txt',
   #     'files_per_job':100
   # },
   # 'Charmonium+Run2016D':{
   #     'file_with_pfns':'lists/Charmonium+Run2016D-17Jul2018-v1+MINIAOD.txt',
   #     'files_per_job':100
   # },
   # 'Charmonium+Run2016E':{
   #     'file_with_pfns':'lists/Charmonium+Run2016E-17Jul2018-v1+MINIAOD.txt',
   #     'files_per_job':100
   # },
   # 'Charmonium+Run2016F':{
   #     'file_with_pfns':'lists/Charmonium+Run2016F-17Jul2018-v1+MINIAOD.txt',
   #     'files_per_job':100
   # },
   # 'Charmonium+Run2016G':{
   #     'file_with_pfns':'lists/Charmonium+Run2016G-17Jul2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2016H':{
   #     'file_with_pfns':'lists/Charmonium+Run2016H-17Jul2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2017B':{
   #     'file_with_pfns':'lists/Charmonium+Run2017B-31Mar2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2017C':{
   #     'file_with_pfns':'lists/Charmonium+Run2017C-31Mar2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2017D':{
   #    'file_with_pfns':'lists/Charmonium+Run2017D-31Mar2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2017E':{
   #     'file_with_pfns':'lists/Charmonium+Run2017E-31Mar2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2017F':{
   #     'file_with_pfns':'lists/Charmonium+Run2017F-31Mar2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2018A':{
   #     'file_with_pfns':'lists/Charmonium+Run2018A-17Sep2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   # },
   # 'Charmonium+Run2018B':{
   #      'file_with_pfns':'lists/Charmonium+Run2018B-17Sep2018-v1+MINIAOD.txt',
   #      'files_per_job':100,
   # },
   #   'Charmonium+Run2018C':{
   #      'file_with_pfns':'lists/Charmonium+Run2018C-17Sep2018-v1+MINIAOD.txt',
   #     'files_per_job':100,
   #  },
   # 'Charmonium+Run2018D':{
   #     'file_with_pfns':'lists/Charmonium+Run2018D-PromptReco-v2+MINIAOD.txt',
   #     'files_per_job':100,
   # },
    # 'BsToMuMu':{
    #     'file_with_pfns':'lists/BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM.txt',
    #     'match':True
    # },
    # 'BsToMuMu_BMuonFilter_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'lists/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt',
    #     'files_per_job':10,
    # },
    'BsToMuMu_BMuonFilter_RunIIFall17MiniAODv2':{
        'file_with_pfns':'lists/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM.txt',
        'files_per_job':10,
    },
    # 'QCD_HT100to200_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'lists/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
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
    #     'file_with_pfns':'lists/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-30to50_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'lists/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-50to80_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'lists/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-30to50_RunIIAutumn18MiniAOD':{
    #     'file_with_pfns':'lists/QCD_Pt_30to50_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    'QCD_Pt-30toInf_BmmGenFilter':{
        'file_with_pfns':'lists/CD_Pt-30toInf_BmmGenFilter-NanoAOD.txt',
        'files_per_job':10,
    }
}

def getPFNs(lfns):
    files = []
    for file in lfns:
        if not input_is_LFNs:
            files.append(file)
        else:
            if use_xrootd:
                files.append("root://cms-xrd-global.cern.ch/%s"%file)
            else:
                fullpath = "/eos/cms/" + file
                if os.path.exists(fullpath):
                    files.append(fullpath)
                else:
                    raise Exception("File not found: %s" % fullpath)
    return files

def chunks(input_list, n, m):
    result = []
    n_elements = len(input_list)
    for i in range(0, n_elements, n*m):
        sub_range = []
        for j in range(n):
            if i+m*j>=n_elements:
                break
            sub_range.append(input_list[i + m*j : i + m*j + m])
        result.append(sub_range)
    return result

def _is_good_file(f):
    tf = TFile.Open(f)
    if not tf: return False
    if tf.IsZombie() or tf.TestBit(TFile.kRecovered): return False
    if not tf.Get("Events"): return False
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

def process_job(input_files):
    if len(input_files)==0: return None
    subdir = ""
    match = re.search("\/([^\/]+)\/[^\/]+$",input_files[0])
    if match:
        subdir = match.group(1)
    job_output_dir = "%s/%s" % (output_dir,subdir)
    output_filename = "%s/%s.root" % (job_output_dir,hashlib.md5(",".join(input_files)).hexdigest())
    if output_is_already_available(output_filename):
        return output_filename
    # skim first
    p = PostProcessor(job_output_dir,input_files,
            cut = cut,
            compression = "LZMA:4",   # need to tune it
            postfix = "_skim")
    p.run()
    # merge
    skimmed_files = []
    for f in input_files:
        skimmed_files.append(os.path.join(job_output_dir, os.path.basename(f).replace(".root","_skim.root")))
    if len(skimmed_files)>1:
        subprocess.call("haddnano.py %s %s" % (output_filename, " ".join(skimmed_files)),shell=True)
        for f in skimmed_files:
            os.remove(f)
    elif len(skimmed_files)==1:
        os.rename(skimmed_files[0],output_filename)
    else:
        raise Exception("Fatal Error: no skimmped input")
    if _is_good_file(output_filename):
        return output_filename
    else:
        return None

def parellel_process(input_files,files_per_job):
    queue_depth = 10 # keep it small to avoid memory leaks
    r = len(input_files)/nProcesses
    if r<files_per_job:
        if r>0: 
            files_per_job = r
        else:
            files_per_job = 1
        
    files = chunks(input_files,queue_depth*nProcesses,files_per_job)
    # pprint.pprint(files)

    results = []
    for job_list in files:
        # pprint.pprint(job_list)
        pool = multiprocessing.Pool(nProcesses) 
        results.extend(pool.map(process_job, job_list))
        pool.close()
    print "Multiprocessing is done."

    # good_files = []
    # for rfile in results:
    #     if rfile: good_files.append(rfile)
    # status = subprocess.call("hadd -f %s %s" % (output_file_name," ".join(good_files)),shell=True)
    # if status==0:
    #     print "Merged output."
    #     for file in good_files:
    #         os.remove(file)
    # else:
    #     print "Merge failed"

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
    parellel_process(input_files,info['files_per_job'])

