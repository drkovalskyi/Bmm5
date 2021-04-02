import re
import tempfile
import os
import subprocess, commands
import json
import time
import platform
import shutil
import postprocessing_cfg as cfg
from pprint import pprint
import hashlib
import fcntl
import sys

class Processor(object):
    """Base class for processors"""
    def __init__(self, job_filename, take_ownership=False):
        """Set up job"""
        self.job_filename = job_filename
        self.take_ownership = take_ownership
        
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
        print "processing %s at %s " % (self.job_filename, platform.node())
        
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
            subprocess.call("rm -v %s " % self.job_lock, shell=True) 

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
