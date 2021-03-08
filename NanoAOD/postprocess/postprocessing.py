import re
import tempfile
import os
import subprocess, commands
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import json
import time
import platform
import shutil
import postprocessing_cfg as cfg
from pprint import pprint
import hashlib

class Processor:
    """Base class for processors"""
    def __init__(self, job_filename, take_ownership=False):
        """Set up job"""
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

        # Create a lock
        self._update_lock(take_ownership)

        # Create a temporary directory
        self.tmp_dir = tempfile.mkdtemp()
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
            os.remove(self.job_lock)

    def _finalize(self):
        """Finish processing and clean up"""
        info = json.load(open(self.job_lock))

        # check if we own the lock
        if info['pid'] != os.getpid():
            raise Exception("The job is locked. Ownership information:\n" + str(info))

        # FIXME: should use xrootd for the transfer to EOS
        shutil.move(self.job_output_tmp, self.job_ouput)
        os.rmdir(self.tmp_dir)

    def _process(self):
        """Abstract interface to implement specific processing actions in derived classes"""
        pass

    def process(self):
        """Process job and clean up"""
        self._process()
        self._finalize()
        self._release_lock()


class Skimmer(Processor):
    """Processor to Skim files and Merge output"""
    def _process(self):
        # check for missing information
        for parameter in ['cut', 'input']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

        input_files = self.job_info['input']

        # skim data
        postfix = "_Skim"
        processor = PostProcessor(self.tmp_dir,
                                  input_files,
                                  cut=self.job_info['cut'],
                                  compression="LZMA:9",                
                                  postfix=postfix)
        processor.run()

        # merge skimed data
        skimmed_files = []
        for f in input_files:
            skimmed_files.append(os.path.join(self.tmp_dir,
                                              os.path.basename(f).replace(".root","%s.root" % postfix)))
        if len(skimmed_files) > 1:
            subprocess.call("haddnano.py %s %s" % (self.job_output_tmp, " ".join(skimmed_files)),shell=True)
        # clean up
        for f in skimmed_files:
            os.remove(f)

class ResourceHandler:
    """Base class for resource handlers"""
    def _processor_name(self, job):
        job_info = json.load(open(job))
        return job_info['processor']
        
    def submit_job(self, job_filename):
        raise Exception("Not implemented")

    def number_of_free_slots(self):
        raise Exception("Not implemented")

    def name(self):
        pass

class SSHResourceHandler(ResourceHandler):
    """Resource handler for ssh-based job execution"""
    def __init__(self, site, max_number_of_jobs_running):
        self.site = site
        self.max_njobs = max_number_of_jobs_running

    def submit_job(self, job_filename):
        path = cfg.workdir
        log = re.sub('\.job$', '\.log', job_filename)
        command = "ssh %s '(cd %s; cmsenv; aklog; python remote_job_starter.py %s %s >& %s) >& /dev/null &'"
        print "submitting %s" % job_filename
        subprocess.call(command % (self.site, path, self._processor_name(job_filename),
                                   job_filename, log),
                        shell=True)
            
        # print "Connecting to %s" % self.site
        # process = subprocess.Popen("ssh -T %s" % self.site, shell=True,
        #                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        # if process:
        #     out, err = process.communicate("""
        #     cd /afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/postprocess
        #     cmsenv
        #     aklog
        #     tokens >& /tmp/dmytro/tokens.log
        #     python submission_test.py %s %s >& /eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.log &
        #     """ % (processor_name, job_filename))
        #     print out

    def number_of_running_jobs(self):
        command = "ssh %s \"ps -Af | grep '[r]emote_job_starter.py'|wc -l\"" % self.site
        code, response = commands.getstatusoutput(command)
        if code not in [0, 256]:
            # 256 is 1 from grep when it doesn't match
            raise Exception("Exit code: %s for command:\n%s" % (code, command))
        return int(response)
        
    def number_of_free_slots(self):
        n = self.max_njobs - self.number_of_running_jobs()
        if n<0: n=0
        return n

    def name(self):
        return self.site

class LocalResourceHandler(ResourceHandler):
    """Resource handler for local job execution"""
    def __init__(self, max_number_of_jobs_running):
        self.max_njobs = max_number_of_jobs_running

    def submit_job(self, job_filename):
        path = cfg.workdir
        log = re.sub('\.job$', '\.log', job_filename)
        command = "python remote_job_starter.py %s %s >& %s &"
        print "submitting %s" % job_filename
        subprocess.call(command % (self._processor_name(job_filename), job_filename, log),
                        shell=True)
        
    def number_of_running_jobs(self):
        return int(subprocess.check_output("ps -Af | grep '[r]emote_job_starter.py'|wc -l", shell=True))
        
    def number_of_free_slots(self):
        n = self.max_njobs - self.number_of_running_jobs()
        if n<0: n=0
        return n

    def name(self):
        return "localhost"

def chunks(input_list, n):
    result = []
    n_elements = len(input_list)
    for i in range(0, n_elements, n):
        result.append(input_list[i : i + n])
    return result

def filename(input_files):
    """Generate unique file name based on the hash of input file names"""
    return hashlib.md5(",".join(input_files)).hexdigest()

class JobCreator:
    """Create jobs according to the specifications in the config file"""
    def __init__(self):
        self.all_inputs_by_datasets = dict()
        self.files_in_use_by_task_and_dataset = dict()
        
    def load_existing_jobs(self):
        """Find existings jobs and store their input"""
        command = 'find -L %s -type f -name "*job" -path "*/%u/*"' % (cfg.output_location, cfg.version)
        all_jobs = subprocess.check_output(command, shell=True).splitlines()
        print "Found %u jobs" % len(all_jobs)
        njobs = 0
        for job in all_jobs:
            match = re.search("([^\/]+)\/(\d+)\/([^\/]+)\/([^\/]+)", job)
            if match:
                if int(match.group(2)) != cfg.version:
                    continue
                task_id = "%s-%s" % (match.group(1), match.group(3))
                if task_id not in self.files_in_use_by_task_and_dataset:
                    self.files_in_use_by_task_and_dataset[task_id] = dict()
                dataset = match.group(4)
                if dataset not in self.files_in_use_by_task_and_dataset[task_id]:
                    self.files_in_use_by_task_and_dataset[task_id][dataset] = dict()
                job_info = json.load(open(job))
                njobs += 1
                for file in job_info['input']:
                    # get eos file name without xrootd prefix
                    file = re.sub("^.*?\/eos\/cms\/store\/", "/eos/cms/store/", file)
                    if file in self.files_in_use_by_task_and_dataset[task_id][dataset]:
                        raise Exception("Found the same input file in different jobs\n" + 
                                        "Task: %s\nDataset: %s\nFile: %s" % (task_id, dataset, file))
                    else:
                        self.files_in_use_by_task_and_dataset[task_id][dataset][file] = 1
        print "Number of processed valid jobs: %u" % njobs

    def find_all_inputs(self):
        """Find all files and splits them in datasets"""
        # look for files that were not modified at least for 30 mins to avoid interferences with transfers
        command = 'find -L %s/%s -mmin +30 -type f -name "*root"' % (cfg.input_location, cfg.version)
        all_inputs = subprocess.check_output(command, shell=True).splitlines()
        all_inputs.sort()
        print "Total number of input file: %u" % len(all_inputs)
        for entry in all_inputs:
            match = re.search("([^\/]+)\/[^\/]+\.root", entry)
            if match:
                dataset = match.group(1)
                if dataset not in self.all_inputs_by_datasets:
                    self.all_inputs_by_datasets[dataset] = []
                self.all_inputs_by_datasets[dataset].append(entry)
            else:
                raise Exception("Failed to get dataset name for file %s" % entry)
        print "Number of datasets: %u" % len(self.all_inputs_by_datasets)
        
    def create_new_jobs(self):
        """Find new files and create jobs"""
        for task in cfg.tasks:
            task_id = "%s-%s" % (task['type'], task['name'])
            print "Processing task %s" % task_id

            for dataset, ds_inputs in self.all_inputs_by_datasets.items():
                if not re.search(task['input_pattern'], dataset): continue
                print "  Processing dataset %s" % dataset
                # find new inputs
                new_inputs = []
                for input in ds_inputs:
                    if task_id in self.files_in_use_by_task_and_dataset:
                        if dataset in self.files_in_use_by_task_and_dataset[task_id]:
                            if input in self.files_in_use_by_task_and_dataset[task_id][dataset]:
                                continue
                    new_inputs.append(input)
                print "    Number of new input files %u" % len(new_inputs)

                # create jobs
                n_elements = len(new_inputs)
                n = task['files_per_job']
                njobs = 0
                for i in range(0, n_elements, n):
                    # get input
                    inputs = new_inputs[i : i + n]
                    inputs.sort()

                    # create job unique id
                    job_id = filename(inputs)

                    # prepare job information
                    job_info = dict()
                    job_info['cut'] = task['cut']
                    job_info['input'] = []
                    job_info['processor'] = task['processor']
                    for entry in inputs:
                        job_info['input'].append(cfg.xrootd_prefix + entry)
                    
                    # prepare output
                    job_dir = "%s/%s/%s/%s/%s" % (cfg.output_location, task['type'], cfg.version, 
                                                  task['name'], dataset)
                    if not os.path.exists(job_dir):
                        os.makedirs(job_dir)
                    job_filename = "%s/%s.job" % (job_dir, job_id)
                    
                    # save job
                    json.dump(job_info, open(job_filename, "w"))
                    njobs += 1
                print "    Number of new jobs created %u" % njobs

class JobDispatcher:
    """Job scheduling"""
    def __init__(self, lifetime=3600):
        self.lock = None
        self.start = time.time()
        self.resources = []
        self.sleep = 60
        self.all_jobs = None
        self.jobs_by_status = None

        # for name,info in cfg.resources.items():
        for resource in cfg.resources:
            # self.resources.append(eval(info['type'])(*info['args']))
            self.resources.append(eval(resource))
    
    def show_resource_availability(self):
        for resource in self.resources:
            print "%s - free slots: %u" % (resource.name(), resource.number_of_free_slots())

    def get_job_status(self, job):
        match = re.search("^(.*?)\.job$", job)
        if match:
            fname = match.group(1)
            output = fname + ".root"
            lock = fname + ".lock"
            log = fname + ".log"
            if os.path.exists(lock):
                # FIXME: this is may also be a failure
                return "Processing"
            elif os.path.exists(output):
                return "Done"
            elif os.path.exists(log):
                return "Failed"
            else: 
                return "New"
        else:
            raise Exception("Incorrect job name:\n%s" % job)
            
    def load_existing_jobs(self):
        """Find existings jobs and store their input"""
        command = 'find -L %s -type f -name "*job" -path "*/%u/*"' % (cfg.output_location, cfg.version)
        self.all_jobs = subprocess.check_output(command, shell=True).splitlines()
        print "Found %u jobs" % len(self.all_jobs)

    def update_status_of_jobs(self):
        self.jobs_by_status = {}
        for job in self.all_jobs:
            status = self.get_job_status(job)
            if status not in self.jobs_by_status:
                self.jobs_by_status[status] = []
            self.jobs_by_status[status].append(job)

        for status in self.jobs_by_status:
            print "\t%s: %u" % (status, len(self.jobs_by_status[status]))

    def process_jobs(self):
        self.update_status_of_jobs()
        if 'New' not in self.jobs_by_status or self.jobs_by_status['New'] == 0:
            print "No new jobs to submit"
            return
        for resource in self.resources:
            n_slots = resource.number_of_free_slots()
            print "%s has %u free slots" % (resource.name(), n_slots)
            for i in range(n_slots):
                resource.submit_job(self.jobs_by_status['New'].pop())
            
            
        
if __name__ == "__main__":
    # p = Skimmer("/eos/cms/store/group/phys_bphys/bmm/bmm5/tmp/NanoAOD-skims/510/ks/EGamma+Run2018B-17Sep2018-v1+MINIAOD/02d352594241e74ddcdbfeb18e2be0d0.job")
    # print p.__dict__
    # p.process()

    # vocms001 = SSHResourceHandler('vocms001.cern.ch', 16)
    # vocms001.submit_job("Skimmer","/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job")
    # time.sleep(10)
    # print vocms001.number_of_free_slots()

    # lh = LocalResourceHandler(16)
    # lh.submit_job("/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job")
    # print lh.number_of_free_slots()

    # jc = JobCreator()
    # jc.find_all_inputs()
    # jc.load_existing_jobs()
    # jc.create_new_jobs()

    jd = JobDispatcher()
    jd.show_resource_availability()
    jd.load_existing_jobs()
    jd.update_status_of_jobs()
    # jd.process_jobs()
