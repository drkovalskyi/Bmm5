import re
import tempfile
import os
import subprocess, subprocess
import json
import time
import platform
import shutil
import postprocessing_cfg as cfg
from pprint import pprint
import hashlib
import fcntl
import sys

from Bmm5.MVA.mtree import MTree
from ROOT import TFile, TTree
import numpy as np

class Processor(object):
    """Base class for processors"""
    def __init__(self, job_filename, take_ownership=False):
        """Set up job"""
        self.job_filename = job_filename
        self.take_ownership = take_ownership
        self.limit = -1
        
        # Load job information
        self.job_info = json.load(open(job_filename))

        # Get directory and job names
        match = re.search("^(.*?)\/([^\/]+)\.job$", job_filename)
        if match:
            self.job_output_dir = match.group(1)
            self.job_name = match.group(2)
            fname = "%s/%s"% (self.job_output_dir, self.job_name)
            self.job_ouput = fname + ".root"
            self.job_lock  = fname + ".lock"
            self.job_log   = fname + ".log"
        else:
            raise Exception("Incorrect input name:\n%s" % job_filename)

    def _prepare(self):
        print("processing %s at %s " % (self.job_filename, platform.node()))
        
        # Create a lock
        self._update_lock(self.take_ownership)

        # Create a temporary directory
        self.tmp_dir = tempfile.mkdtemp(prefix=cfg.tmp_prefix)
        self.job_output_tmp = "%s/%s.root" % (self.tmp_dir, self.job_name)

    def _update_lock(self, take_ownership=False):
        """Create and update job lock"""
        # check if lock exists
        if os.path.exists(self.job_lock):
            info = json.load(open(self.job_lock))

            # check if we own the lock
            if info['pid'] != os.getpid():
                if not take_ownership:
                    raise Exception("The job is locked. Ownership information:\n" + str(info))
        # update
        info = {'pid':os.getpid(), 'node':platform.node(), 'lastupdate':time.time()}
        json.dump(info, open(self.job_lock,'w'))

    def _release_lock(self):
        """Release lock"""
        info = json.load(open(self.job_lock))

        # check if we own the lock
        if info['pid'] == os.getpid():
            if re.search('^\/eos\/', self.job_lock):
                subprocess.call("eos rm %s " % self.job_lock, shell=True)
            else:
                subprocess.call("rm -v %s " % self.job_lock, shell=True)
        else:
            raise Exception("The job is locked. Ownership information:\n" + str(info))
            

    def _finalize(self):
        """Finish processing and clean up"""
        info = json.load(open(self.job_lock))

        # check if we own the lock
        if info['pid'] != os.getpid():
            raise Exception("The job is locked. Ownership information:\n" + str(info))

        # FIXME: should use xrootd for the transfer to EOS
        sys.stdout.flush()
        # shutil.move(self.job_output_tmp, self.job_ouput)
        subprocess.call("mv -v %s %s" % (self.job_output_tmp, self.job_ouput), shell=True)
        # os.rmdir(self.tmp_dir)
        subprocess.call("rm -v -d %s " % self.tmp_dir, shell=True)

    def _process(self):
        """Abstract interface to implement specific processing actions in derived classes"""
        pass

    def process(self):
        """Process job and clean up"""
        self._prepare()
        self._process()
        self._finalize()
        self._release_lock()

class FlatNtupleBase(Processor):
    """Flat ROOT ntuple producer for Bmm5 analysis"""

    def __init__(self, job_filename, take_ownership=False):
        self.n_gen_all = None
        self.n_gen_passed = None
        super(FlatNtupleBase, self).__init__(job_filename, take_ownership)
    

    def _validate_inputs(self):
        """Task specific input validation"""
        raise Exception("Not implemented")


    def _is_certified(self, event, type):
        # load certification information
        if type not in self.goodruns:
            self.goodruns[type] = dict()
            for f in os.listdir('certification/%s' % type):
                self.goodruns[type].update(json.load(open('certification/%s/%s' % (type, f))))
            print("Number of runs in the %s certification: %u" % (type, len(self.goodruns[type])))

        # run number is a string for some reason
        run = str(event.run)

        if run not in self.goodruns[type]:
            return False
        
        for min_lumi, max_lumi in self.goodruns[type][run]:
            if event.luminosityBlock >= min_lumi and event.luminosityBlock <= max_lumi:
                return True
        return False


    def _process(self):
        """Process input files, merge output and report performance"""
        self._validate_inputs()

        # process files
        results = []
        t0 = time.perf_counter()
        n_events = 0
        for f in self.job_info['input']:
            result, n = self.process_file(f)
            n_events += n
            results.append(result)
        print(n_events//(time.perf_counter() - t0), "Hz")

        print("Merging output.")
        # merge results
        good_files = []
        for rfile in results:
            if rfile: good_files.append(rfile)
        command = "hadd -f %s %s" % (self.job_output_tmp, " ".join(good_files))
        status = subprocess.call(command, shell=True)

        if status==0:
            print("Merged output.")
            for file in good_files:
                os.remove(file)

            # Store meta data
            f = TFile(self.job_output_tmp, "UPDATE")
            t = TTree("info","Selection information")
            
            n_processed = np.empty((1), dtype="i")
            t.Branch("n_processed", n_processed, "n_processed/I")
            n_processed[0] = n_events

            if self.n_gen_all != None:
                n_gen_all = np.empty((1), dtype="i")
                n_gen_passed = np.empty((1), dtype="i")
                t.Branch("n_gen_all", n_gen_all, "n_gen_all/I")
                t.Branch("n_gen_passed", n_gen_passed, "n_gen_passed/I")
                n_gen_all[0] = self.n_gen_all
                n_gen_passed[0] = self.n_gen_passed
            
            t.Fill()
            f.Write()
            f.Close()
        else:
            raise Exception("Merge failed")

    def _process_events(self):
        raise Exception("Not implemented")

    def process_file(self, input_file):
        """Initialize input and output trees and initiate the event loop"""
        print("Processing file: %s" % input_file)
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

        # GenFilterInfo
        lumis = fin.Get("LuminosityBlocks")
        info = fin.Get("info")
        if lumis:
            for lumi in lumis:
                if hasattr(lumi, 'GenFilter_numEventsPassed'):
                    if self.n_gen_all == None:
                        self.n_gen_all = 0
                        self.n_gen_passed = 0
                    self.n_gen_passed += lumi.GenFilter_numEventsPassed
                    self.n_gen_all    += lumi.GenFilter_numEventsTotal
        elif info:
            for entry in info:
                if hasattr(entry, 'n_gen_all'):
                    if self.n_gen_all == None:
                        self.n_gen_all = 0
                        self.n_gen_passed = 0
                    self.n_gen_passed += entry.n_gen_passed
                    self.n_gen_all    += entry.n_gen_all

        self.input_tree = fin.Get("Events")
        nevents = self.input_tree.GetEntries()

        self._process_events()

        nout = self.tree.tree.GetEntries()
        print('Selected %d / %d entries from %s (%.2f%%)' % (nout, nevents, input_file, 100.*nout/nevents if nevents else 0))
        fout.Write()
        fout.Close()
        return output_filename, nevents

    def get_cut(self):
        """Parse cut for keywords and add placeholders for tree and index"""

        # load branch info
        branches = dict()

        # collect branch names, corresponding leafs and their event
        # counter for arrays
        for br in self.input_tree.GetListOfBranches():
            name = br.GetName()
            leaf = br.GetLeaf(name)
            if leaf.GetLeafCount():
                branches[name] = leaf.GetLeafCount().GetName()
            else:
                branches[name] = ""

        # process cut

        cut = self.job_info['cut']

        # replace ROOT style AND with the one that is acceptable for python
        cut = re.sub('\&\&', ' and ', cut)

        # tokenize the cut string into elements so that we can add the tree and element index
        cut_list = re.split('([^\w\_]+)', cut)

        parsed_cut = ""
        for i in range(len(cut_list)):
            if not re.search('^[\w\_]+$', cut_list[i]) or cut_list[i] not in branches:
                # element is not a branch name, store and move on
                parsed_cut += cut_list[i]
            else:
                # element is a branch name

                # make an index formater for arrays
                index = branches[cut_list[i]]
                if index != "":
                    index = "[{" + index + "}]"

                if i < len(cut_list)-1 and re.search('^\[', cut_list[i+1]):
                    parsed_cut += "{tree}." + cut_list[i]
                else:
                    parsed_cut += "{tree}." + cut_list[i] + index
        return parsed_cut


class ResourceHandler(object):
    """Base class for resource handlers"""
    def __init__(self):
        self.active_jobs = set() # keep track of submitted jobs by the handler
        
    def _processor_name(self, job):
        job_info = json.load(open(job))
        return job_info['processor']
        
    def submit_job(self, job):
        """Submit and keep track of a job"""
        self.active_jobs.add(job)
        self._submit_job(job)

    def _submit_job(self, job):
        raise Exception("Not implemented")

    def number_of_free_slots(self):
        raise Exception("Not implemented")

    def name(self):
        pass

    def get_running_jobs(self):
        jobs = self._get_running_jobs()
        self.active_jobs = self.active_jobs.intersection(set(jobs))
        return jobs
    
    def _get_running_jobs(self):
        raise Exception("Not implemented")

    def number_of_running_jobs(self, owned=True):
        """Get number of running jobs. Can be restricted to only owned jobs"""
        jobs = self.get_running_jobs()
        if not owned:
            return len(jobs)
        else:
            return len(self.active_jobs)

    def number_of_free_slots(self):
        n = self.max_njobs - self.number_of_running_jobs()
        if n<0: n=0
        return n

    def kill_all_jobs(self):
        raise Exception("Not implemented")

    def clean_up(self):
        """Remove temporary directories for failed jobs"""
        raise Exception("Not implemented")
