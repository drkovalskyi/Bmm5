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
import fcntl
import sys

class Processor(object):
    """Base class for processors"""
    def __init__(self, job_filename, take_ownership=False):
        """Set up job"""
        # Load job information
        self.job_info = json.load(open(job_filename))

        print "processing %s at %s " % (job_filename, platform.node())
        
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
        sys.stdout.flush()
        processor.run()

        # merge skimed data
        skimmed_files = []
        for f in input_files:
            skimmed_files.append(os.path.join(self.tmp_dir,
                                              os.path.basename(f).replace(".root","%s.root" % postfix)))
        sys.stdout.flush()
        # FIXME: check the logic
        if len(skimmed_files) > 1:
            subprocess.call("haddnano.py %s %s" % (self.job_output_tmp, " ".join(skimmed_files)), shell=True)
        # clean up
        for f in skimmed_files:
            os.remove(f)

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

    
class SSHResourceHandler(ResourceHandler):
    """Resource handler for ssh-based job execution"""
    def __init__(self, site, max_number_of_jobs_running):
        super(SSHResourceHandler, self).__init__()
        self.site = site
        self.max_njobs = max_number_of_jobs_running
        self.proc = subprocess.Popen("ssh -T -x %s 'bash -l'" % site,
                                     shell=True,
                                     stdin  = subprocess.PIPE,
                                     stdout = subprocess.PIPE
                                     )
        # prepare working area
        self._send_command_and_get_response("""
        cd /afs/cern.ch/work/d/dmytro/projects/RunII-NanoAODv6/src/Bmm5/NanoAOD/postprocess
        eval `scramv1 runtime -sh`
        """)
        
    def _non_block_read(self):
        fd = self.proc.stdout.fileno()
        fl = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
        try:
            return self.proc.stdout.read()
        except:
            return ""
        
    def _send_command_and_get_response(self, command):
        end_of_transmission = "end_of_transmission"

        # send command with the end of transimission echo
        self.proc.stdin.write(command + "\n")
        self.proc.stdin.write("echo '%s'\n" % end_of_transmission)
        self.proc.stdin.flush()

        response = ""
        while True:
            line = self._non_block_read()
            if line:
                response += line
                if re.search(end_of_transmission, line):
                    break
        # remove end_of_transmission_pattern
        response = re.sub("%s.*$" % end_of_transmission, "", response) 
        return response.rstrip()
        
    def check_server_status(self):
        # need to check if the server is responding properly within given amount of time
        pass
            
    def _submit_job(self, job):
        log = re.sub('\.job$', '\.log', job)
        command = "nice python job_starter.py %s %s &> %s &\n" % (self._processor_name(job),
                                                                         job, log)
        print "submitting %s" % job
            
        print self._send_command_and_get_response(command),

    def _get_running_jobs(self):
        jobs = []
        response = self._send_command_and_get_response("ps -Af | grep '[j]ob_starter.py'")
        for line in response.splitlines():
            match = re.search('(\S+\.job)', line)
            if match:
                jobs.append(match.group(1))
        return jobs

    def wait_for_jobs_to_finish(self):
        time.sleep(10)
        while self.number_of_running_jobs(owned=True) > 0:
            time.sleep(5)
    
    def name(self):
        return self.site

    def kill_all_jobs(self):
        print "Killing jobs at %s" % self.name()
        response = self._send_command_and_get_response("ps -Af")
        for line in response.splitlines():
            match = re.search('^\S+\s+(\S+).*?job_starter.py', line)
            if match:
                self._send_command_and_get_response("kill %s" % match.group(1))

class LocalResourceHandler(ResourceHandler):
    """Resource handler for local job execution"""
    def __init__(self, max_number_of_jobs_running):
        super(LocalResourceHandler, self).__init__()
        self.max_njobs = max_number_of_jobs_running

    def _submit_job(self, job):
        log = re.sub('\.job$', '\.log', job)
        command = "nice python job_starter.py %s %s >& %s &"
        print "submitting %s" % job
        subprocess.call(command % (self._processor_name(job), job, log),
                        shell=True)
        
    def _get_running_jobs(self):
        jobs = []
        response = subprocess.check_output("ps -Af", shell=True, env={})
        for line in response.splitlines():
            if not re.search('job_starter.py', line): continue
            match = re.search('job_starter.py.*?(\S+\.job)', line)
            if match:
                jobs.append(match.group(1))
        return jobs

    def name(self):
        return "localhost"

    def kill_all_jobs(self):
        print "Killing jobs at %s" % self.name()
        response = subprocess.check_output("ps -Af", shell=True)
        for line in response.splitlines():
            match = re.search('^\S+\s+(\S+).*?job_starter.py', line)
            if match:
                subprocess.call("kill %s" % match.group(1), shell=True)
                
def chunks(input_list, n):
    result = []
    n_elements = len(input_list)
    for i in range(0, n_elements, n):
        result.append(input_list[i : i + n])
    return result

def filename(input_files):
    """Generate unique file name based on the hash of input file names"""
    return hashlib.md5(",".join(input_files)).hexdigest()

class JobCreator(object):
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

class JobDispatcher(object):
    """Job scheduling"""
    def __init__(self, lifetime=36000):
        """Initialization"""
        self.lock = None
        self.end_time = time.time() + lifetime
        self._resources = None
        self.sleep = 60
        self.all_jobs = None
        self.jobs_by_status = None
        self.running_jobs = []

        self._load_existing_jobs()

    def _init_resources(self):
        """On demand initialization of resource handlers"""
        self._resources = []
        print "Initializing the resources"
        for resource in cfg.resources:
            print "\t", resource
            self._resources.append(eval(resource))

    def resources(self):
        """Resource accessor. Will trigger initialization if necessary"""
        if self._resources == None:
            self._init_resources()
        return self._resources

    def show_resource_availability(self):
        """Current statust of resources"""
        for resource in self.resources():
            print "%s - free slots: %u" % (resource.name(), resource.number_of_free_slots())

    def update_running_jobs(self):
        """Collection information about running jobs from resource handlers"""
        self.running_jobs = []            
        if self._resources == None:
            return
        for resource in self.resources():
            self.running_jobs.extend(resource.get_running_jobs())

    def get_job_status(self, job):
        """Determine job status based on existance of associated files and running information"""
        match = re.search("^(.*?)\.job$", job)
        if match:
            fname = match.group(1)
            output = fname + ".root"
            lock = fname + ".lock"
            log = fname + ".log"
            if os.path.exists(lock):
                if job not in self.running_jobs:
                    return "Failed"
                else:
                    return "Running"
            elif os.path.exists(output):
                return "Done"
            elif os.path.exists(log):
                if job not in self.running_jobs:
                    return "Failed"
                else:
                    return "Running"
            else: 
                return "New"
        else:
            raise Exception("Incorrect job name:\n%s" % job)
            
    def _load_existing_jobs(self):
        """Find existings jobs and store their input"""
        command = 'find -L %s -type f -name "*job" -path "*/%u/*"' % (cfg.output_location, cfg.version)
        self.all_jobs = subprocess.check_output(command, shell=True).splitlines()
        print "Found %u jobs" % len(self.all_jobs)

    def update_status_of_jobs(self):
        """Classify jobs by status and store in corresponding lists"""
        self.jobs_by_status = {}
        self.update_running_jobs()
        for job in self.all_jobs:
            status = self.get_job_status(job)
            if status not in self.jobs_by_status:
                self.jobs_by_status[status] = []
            self.jobs_by_status[status].append(job)

        for status in self.jobs_by_status:
            print "\t%s: %u" % (status, len(self.jobs_by_status[status]))

    def number_of_running_jobs(self):
        """Get total number of running jobs on all resources"""
        n_running = 0
        for resource in self.resources():
            n_running += resource.number_of_running_jobs()
        return n_running
            
    def process_jobs(self):
        """Process available jobs"""

        # submit jobs
        print "Submitting new jobs"
        while time.time() < self.end_time:
            self.update_status_of_jobs()
            if 'New' not in self.jobs_by_status or self.jobs_by_status['New'] == 0:
                print "No new jobs to submit. Reloading to check if new jobs where injected"
                self._load_existing_jobs()                
                self.update_status_of_jobs()
                if 'New' not in self.jobs_by_status or self.jobs_by_status['New'] == 0:
                    print "No new jobs is found."
                    break
            for resource in self.resources():
                n_slots = resource.number_of_free_slots()
                print "%s has %u free slots" % (resource.name(), n_slots)
                for i in range(n_slots):
                    if len(self.jobs_by_status['New']) > 0:
                        resource.submit_job(self.jobs_by_status['New'].pop())
                    else:
                        break
            time.sleep(60)

        # finalize running jobs
        print "Finalizing running jobs"
        while True:
            n_running = self.number_of_running_jobs(owned=True)
            print "Number of running jobs: %u" % n_running
            if n_running == 0:
                break
            time.sleep(60)
            
    def show_failures(self, detailed=False):
        self.update_status_of_jobs()
        failures = {}

        if 'Failed' in self.jobs_by_status:
            for job in self.jobs_by_status['Failed']:
                match = re.search("^(.*?)\.job$", job)
                if match:
                    fname = match.group(1)
                    output = fname + ".root"
                    lock = fname + ".lock"
                    log = fname + ".log"

                    failure_type = None
                    if os.path.exists(output):
                        failure_type = 'Output is available'
                    elif os.path.exists(lock):
                        failure_type = 'Locked without output'
                    elif os.path.exists(log):
                        failure_type = 'Only log'
                        if detailed:
                            subprocess.call("tail %s" % log, shell=True)

                    if failure_type:
                        if failure_type not in failures:
                            failures[failure_type] = []
                        failures[failure_type].append(job)

        for failure_type, jobs in failures.items():
            print "Failure type: %s" % failure_type
            print "\tNumber of jobs: %u" % len(jobs)

    def reset_failures(self):
        self.update_status_of_jobs()

        # find all jobs that need to be reset
        jobs_to_reset = []
        if 'Failed' in self.jobs_by_status:
            jobs_to_reset.extend(self.jobs_by_status['Failed'])

        # back up existing output
        for job in jobs_to_reset:
            match = re.search("^(.*?)\.job$", job)
            if match:
                fname = match.group(1)
                output = fname + ".root"
                lock = fname + ".lock"
                log = fname + ".log"
    
                for f in [output, lock, log]:
                    # remove previous backups
                    if os.path.exists(f + ".failed"):
                        os.remove(f + ".failed")
                    # backup latest failure
                    if os.path.exists(f):
                        shutil.move(f, f + ".failed")

    def kill_all_jobs(self):
        for resource in self.resources():
            resource.kill_all_jobs()
                
if __name__ == "__main__":
    # p = Skimmer("/eos/cms/store/group/phys_bphys/bmm/bmm5/tmp/NanoAOD-skims/510/ks/EGamma+Run2018B-17Sep2018-v1+MINIAOD/02d352594241e74ddcdbfeb18e2be0d0.job")
    # print p.__dict__
    # p.process()

    test_job = "/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job"
    vocms001 = SSHResourceHandler('vocms0500.cern.ch', 16)
    vocms001.submit_job(test_job)
    vocms001.wait_for_jobs_to_finish()

    # lh = LocalResourceHandler(16)
    # lh.submit_job("/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job")
    # print lh.get_running_jobs()
    # print lh.number_of_free_slots()

    # jc = JobCreator()
    # jc.find_all_inputs()
    # jc.load_existing_jobs()
    # jc.create_new_jobs()
    
    # jd = JobDispatcher()
    # jd.show_resource_availability()
    # jd.show_failures(True)
    # jd.reset_failures()
    # jd.update_status_of_jobs()
    # jd.process_jobs()
    # jd.kill_all_jobs()
