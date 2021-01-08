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
from functools import partial
import shutil


version = 510

output_eos_dir = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD-skims/"
# output_eos_dir = "/user/dmytro/NanoAOD-skims/"

output_tmp_dir = "/tmp/dmytro/NanoAOD-skims/"
# output_tmp_dir = "/user/dmytro/tmp/"

nProcesses = 16

force_recreation = False

single_thread = False # debuging mode

# Access mode for LFNs
#
# fuse - eos fuse mount
# xrootd - direct access to eoscms
# AAA - remote access via global redirector
access_mode = "xrootd"

skims = {
    'mm':{
        'cuts':[
            # Trigger like cuts
            "mm_mu1_index>=0", "mm_mu2_index>=0",
            "abs(Muon_eta[mm_mu1_index])<1.4", "Muon_pt[mm_mu1_index]>4",
            "abs(Muon_eta[mm_mu2_index])<1.4", "Muon_pt[mm_mu2_index]>4",
            "mm_kin_mass>4.9", "mm_kin_mass<6.0", "mm_kin_sl3d>4", "mm_kin_vtx_chi2dof<5",
            "HLT_DoubleMu4_3_Bs"
        ]
    },
    'bkmm':{
        'cuts':[
            "mm_mu1_index[bkmm_mm_index]>=0", "mm_mu2_index[bkmm_mm_index]>=0",
            "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4",
            "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4",
            "mm_kin_sl3d[bkmm_mm_index]>4", "mm_kin_vtx_chi2dof[bkmm_mm_index]<5",
            "abs(bkmm_jpsimc_mass-5.4)<0.5", "bkmm_jpsimc_vtx_chi2dof<5","bkmm_jpsimc_alpha<0.2"
        ]
    },
    'ks':{
        'cuts':[
            "ks_kin_cosAlphaXY>0.98"
        ]
    },
    'phi':{
        'cuts':[
            "phi_kin_mass>0"
        ]
    },
    'd0':{
        'cuts':[
            "d0_kin_mass>0"
        ]
    },
    'lambda':{
        'cuts':[
            "lambda_kin_mass>0"
        ]
    },

}

# Input can be provided for each samples in the following formats:
# - 'file_with_lfns' - file name containing input logical file names
# - 'file_with_pfns' - file name containing input physical file names
# - 'lfns'           - list of input logical file names
# - 'pfns'           - list of input physical file names
samples = {
   # 'Charmonium+Run2016B_1':{
   #     'file_with_lfns':'lists/Charmonium+Run2016B-17Jul2018_ver1-v1+MINIAOD.txt',
   #     'files_per_job':100,
   #     'skims':['mm']
   # },
   'Charmonium+Run2016B_2':{
       'file_with_lfns':'lists/Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2016C':{
       'file_with_lfns':'lists/Charmonium+Run2016C-17Jul2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2016D':{
       'file_with_lfns':'lists/Charmonium+Run2016D-17Jul2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2016E':{
       'file_with_lfns':'lists/Charmonium+Run2016E-17Jul2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2016F':{
       'file_with_lfns':'lists/Charmonium+Run2016F-17Jul2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2016G':{
       'file_with_lfns':'lists/Charmonium+Run2016G-17Jul2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2016H':{
       'file_with_lfns':'lists/Charmonium+Run2016H-17Jul2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2017B':{
       'file_with_lfns':'lists/Charmonium+Run2017B-31Mar2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2017C':{
       'file_with_lfns':'lists/Charmonium+Run2017C-31Mar2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2017D':{
      'file_with_lfns':'lists/Charmonium+Run2017D-31Mar2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2017E':{
       'file_with_lfns':'lists/Charmonium+Run2017E-31Mar2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2017F':{
       'file_with_lfns':'lists/Charmonium+Run2017F-31Mar2018-v1+MINIAOD.txt',
       'files_per_job':100,
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2018A':{
       'file_with_lfns':'lists/Charmonium+Run2018A-17Sep2018-v1+MINIAOD.txt',
       'files_per_job':100,
       # 'skims':['mm', 'bkmm', 'ks']
       'skims':['ks','d0','phi','lambda']
   },
   'Charmonium+Run2018B':{
        'file_with_lfns':'lists/Charmonium+Run2018B-17Sep2018-v1+MINIAOD.txt',
        'files_per_job':100,
       # 'skims':['mm', 'bkmm', 'ks']
       'skims':['ks','d0','phi','lambda']
   },
    'Charmonium+Run2018C':{
        'file_with_lfns':'lists/Charmonium+Run2018C-17Sep2018-v1+MINIAOD.txt',
        'files_per_job':100,
        # 'skims':['mm', 'bkmm', 'ks']
        'skims':['ks','d0','phi','lambda']
    },
    'Charmonium+Run2018D':{
        'file_with_lfns':'lists/Charmonium+Run2018D-PromptReco-v2+MINIAOD.txt',
        'files_per_job':100,
        # 'skims':['mm', 'bkmm', 'ks']
         'skims':['ks','d0','phi','lambda']
    },
    # 'BsToMuMu':{
    #     'file_with_lfns':'lists/BsToMuMu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM.txt',
    #     'match':True
    # },
    # 'BsToMuMu_BMuonFilter_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt',
    #     'files_per_job':10,
    #     'skims':['ks','d0','phi','lambda']
    # },
    # 'BsToMuMu_BMuonFilter_RunIIFall17MiniAODv2':{
    #     'file_with_lfns':'lists/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM.txt',
    #     'files_per_job':10,
    #     'skims':['ks','d0','phi','lambda']
    # },
    # 'QCD_HT100to200_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_HT50to100_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/MVA/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_HT100to200_BGenFilter_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/MVA/QCD_HT100to200_BGenFilter_TuneCP5_13TeV-madgraph-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_HT200to300_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'/afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/MVA/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_Pt-80to120_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-30to50_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-50to80_MuEnrichedPt5_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3+MINIAODSIM.txt',
    # },
    # 'QCD_Pt-30to50_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/QCD_Pt_30to50_TuneCP5_13TeV_pythia8+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM.txt'
    # },
    # 'QCD_Pt-30toInf_BmmGenFilter':{
    #     'file_with_lfns':'lists/CD_Pt-30toInf_BmmGenFilter-NanoAOD.txt',
    #     'files_per_job':10,
    # }
    # 'BuToJpsiK_BMuonFilter_RunIIAutumn18MiniAOD':{
    #     'file_with_lfns':'lists/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM.txt',
    #     'files_per_job':10,
    #     'skims':['bkmm','ks']
    # },
    # 'BuToJpsiK_BMuonFilter_RunIIFall17MiniAODv2':{
    #     'file_with_lfns':'lists/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3+MINIAODSIM.txt',
    #     'files_per_job':10,
    #     'skims':['bkmm','ks']
    # },
    # 'BuToJpsiK_BMuonFilter_RunIISummer16MiniAODv2':{
    #     'file_with_lfns':'lists/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen+RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1+MINIAODSIM.txt',
    #     'files_per_job':10,
    #     'skims':['bkmm','ks']
    # },
}

def getPFNs(lfns):
    """Get physical file name from the logical ones"""
    files = []
    for file in lfns:
        if access_mode == "AAA":
            files.append("root://cms-xrd-global.cern.ch/%s"%file)
        elif access_mode == "fuse":
            fullpath = "/eos/cms/" + file
            if os.path.exists(fullpath):
                files.append(fullpath)
            else:
                raise Exception("File not found: %s" % fullpath)
        elif access_mode == "xrootd":
            files.append("root://eoscms.cern.ch://eos/cms%s"%file)
        else:
            raise Exception("Unsupport access mode: %s" % access_mode)
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

def filename(input_files):
    """Generate unique file name based on the hash of input file names"""
    return hashlib.md5(",".join(input_files)).hexdigest()


def save_file_list(input_files, dir, name):
    """Save input files used in the skim to prevent duplicates"""

    file_temp = "%s/%s.temp" % (dir, name)
    file_list = "%s/%s.list" % (dir, name)
    
    if os.path.exists(file_list):
        return

    with open(file_temp, "w") as f:
        for file in input_files:
            f.write("%s\n" % file)

    # make sure that only complete files are stored
    os.rename(file_temp, file_list)


def get_output_dir(path, file, skim):
    """Get director name based on input file name, 
    skim and processing version"""

    subdir = ""
    match = re.search("\/([^\/]+)\/[^\/]+$",file)
    if match:
        subdir = match.group(1)
    return "%s/%s/%s/%s" % (path, version, skim, subdir)


def process_job(input_files, skim):
    """Skim files and merge output"""

    if len(input_files)==0: return None

    job_output_eos_dir = get_output_dir(output_eos_dir, input_files[0], skim)
    job_output_tmp_dir = get_output_dir(output_tmp_dir, input_files[0], skim)

    fname = filename(input_files)

    save_file_list(input_files, job_output_eos_dir, fname)
    
    job_output_eos_filename = "%s/%s.root" % (job_output_eos_dir, fname)
    job_output_tmp_filename = "%s/%s.root" % (job_output_tmp_dir, fname)

    if output_is_already_available(job_output_eos_filename):
        return job_output_eos_filename

    # Skim- data
    cut = "&&".join(skims[skim]['cuts'])

    processor = PostProcessor(job_output_tmp_dir,
                              input_files,
                              cut=cut,
                              compression="LZMA:4",   # need to tune it
                              postfix="_%s" % skim)
    processor.run()

    # merge skimed data
    skimmed_files = []
    for f in input_files:
        skimmed_files.append(os.path.join(job_output_tmp_dir, 
                                          os.path.basename(f).replace(".root","_%s.root" % skim)))
    if len(skimmed_files) > 1:
        subprocess.call("haddnano.py %s %s" % (job_output_tmp_filename, " ".join(skimmed_files)),shell=True)
        for f in skimmed_files:
            os.remove(f)
        shutil.move(job_output_tmp_filename, job_output_eos_filename)
    elif len(skimmed_files) == 1:
        shutil.move(skimmed_files[0], job_output_eos_filename)
    else:
        raise Exception("Fatal Error: no skimmped input")
    if _is_good_file(job_output_eos_filename):
        return job_output_eos_filename
    else:
        return None

def parallel_process(input_files, files_per_job, skims):
    """Split input in jobs and run them in parallel.

    First check if any input files are already assigned to specific
    jobs by checking *.list file content. For existing assignments
    resubmit jobs with the same input lists as they were used during
    first processing attempt. If the jobs are not processable, please
    remove the corresponding list files. All unassigned input files
    will be used to form new jobs.

    Input variable skim is a list of skim names to run for the input
    files. Running time is linearly increasing with a number of skims
    to run since each skim is processed separately.

    Warning: output file name may differ between skims if they were
    not processed with one set of inputs. It shouldn't be an issue,
    but if you want matching name reskim from scratch.
    """
    queue_depth = 10 # keep it small to avoid memory leaks

    for skim in skims:
        print "Skim name: %s" % skim
        new_files = []
        skim_output_eos_dir = get_output_dir(output_eos_dir, input_files[0], skim)

        if not os.path.exists(skim_output_eos_dir):
            os.makedirs(skim_output_eos_dir)

        job_sets = []

        # Find already assigned files
        assigned_files = []
        if os.path.exists(skim_output_eos_dir):
            list_files = subprocess.check_output(
                'find %s -type f -name "*.list"' % skim_output_eos_dir,
                shell = True
            ).split("\n")

            # Put already existing jobs in the first job set. It can lead
            # to a long queue, but we expect these jobs to be complete and
            # nothing, but a simple check will be performed. If that's not
            # the case we can run out of memory using very long job queue.
            for list_file in list_files:
                if list_file != "":
                    # check if we have valid root file
                    root_file = re.sub("\.list$", ".root", list_file)
                    if not _is_good_file(root_file):
                        print "WARNING: no valid ROOT file is found for %s. Reprocessing" % list_file
                        subprocess.call("rm -v %s" % list_file, shell=True)
                    else:
                        with open(list_file) as f:
                            files = []
                            for entry in f:
                                file = entry.strip("\n")
                                if file != "":
                                    assigned_files.append(file)
                            if len(job_sets) == 0:
                                job_sets.append([])
                            job_sets[0].append(files)

            print "Number of assigned files: %u" % len(assigned_files)
            if len(job_sets) > 0:
                print "Number of existing jobs: %u" % len(job_sets[0])

        for f in input_files:
            if not f in assigned_files:
                new_files.append(f)
            
        print "Number of new files: %u" % len(new_files)
                
        r = len(input_files)/nProcesses
        if r<files_per_job:
            if r>0: 
                files_per_job = r
            else:
                files_per_job = 1

        job_sets.extend(chunks(new_files, queue_depth*nProcesses, files_per_job))
        # pprint.pprint(files)

        results = []
        for job_list in job_sets:
            # pprint.pprint(job_list)
            if not single_thread:
                pool = multiprocessing.Pool(nProcesses)
                results.extend(pool.map(partial(process_job,skim=skim), job_list))
                pool.close()
            else:
                for list_of_files in job_list:
                    process_job(list_of_files,skim=skim)

    print "Processing is done."

for sample,info in samples.items():
    print "Processing %s" % sample
    output_filename = "%s.root" % sample
    if output_is_already_available(output_filename):
        continue
    input_files = []
    if 'file_with_pfns' in info:
        with open(info['file_with_pfns']) as f:
            for pfn in f:
                input_files.append(pfn.strip())
    elif 'file_with_lfns' in info:
        with open(info['file_with_lfns']) as f:
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
    parallel_process(input_files,info['files_per_job'],info['skims'])

