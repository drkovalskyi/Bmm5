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

    def get_running_jobs(self):
        raise Exception("Not implemented")

    def number_of_running_jobs(self):
        return len(self.get_running_jobs())

class SSHResourceHandler(ResourceHandler):
    """Resource handler for ssh-based job execution"""
    def __init__(self, site, max_number_of_jobs_running):
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
        self.proc.stdin.write(command)
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
        return response.rsplit("\n",2)[0]
        
        
    def check_server_status(self):
        # need to check if the server is responding properly within given amount of time
        pass
    
        
    def submit_job(self, job_filename):
        path = cfg.workdir
        log = re.sub('\.job$', '\.log', job_filename)
        command = "nohup python remote_job_starter.py %s %s &> %s &\n" % (self._processor_name(job_filename),
                                                                    job_filename, log)
        print "submitting %s" % job_filename
            
        print self._send_command_and_get_response(command)

    def get_running_jobs(self):
        jobs = []
        response = self._send_command_and_get_response("ps -Af | grep '[r]emote_job_starter.py'\n")
        for line in response.splitlines():
            match = re.search('(\S+\.job)', line)
            if match:
                jobs.append(match.group(1))
        return jobs

    def wait_for_all_jobs_completion(self):
        time.sleep(10)
        while self.number_of_running_jobs()>0:
            time.sleep(5)
    
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
        
    def get_running_jobs(self):
        jobs = []
        response = subprocess.check_output("ps -Af | grep '[r]emote_job_starter.py'|wc -l", shell=True)
        for line in response.splitlines():
            match = re.search('(\S+\.job)', line)
            if match:
                jobs.append(match.group(1))
        return jobs

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
    def __init__(self, lifetime=36000):
        self.lock = None
        self.end_time = time.time() + lifetime
        self.resources = []
        self.sleep = 60
        self.all_jobs = None
        self.jobs_by_status = None
        self.running_jobs = dict()

        print "Initializing the resources"
        # for name,info in cfg.resources.items():
        for resource in cfg.resources:
            print "\t", resource
            # self.resources.append(eval(info['type'])(*info['args']))
            self.resources.append(eval(resource))
    
    def show_resource_availability(self):
        for resource in self.resources:
            print "%s - free slots: %u" % (resource.name(), resource.number_of_free_slots())

    def update_running_jobs(self):
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

    def number_of_running_jobs(self):
        n_running = 0
        for resource in self.resources:
            n_running += resource.number_of_running_jobs()
        return n_running
            
    def process_jobs(self):
        """Process available jobs"""

        # submit jobs
        print "Submitting new jobs"
        while time.time() < self.end_time:
            self.update_status_of_jobs()
            if 'New' not in self.jobs_by_status or self.jobs_by_status['New'] == 0:
                print "No new jobs to submit"
                break
            for resource in self.resources:
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
            n_running = self.number_of_running_jobs()
            print "Number of running jobs: %u" % n_running
            if n_running == 0:
                break
            time.sleep(60)
            
    def list_failures(self):
        self.update_status_of_jobs()
        if 'Failed' in self.jobs_by_status:
            print "\nFailed to start"
            pprint(self.jobs_by_status['Failed'])
        if self.number_of_running_jobs() == 0:
            if 'Processing' in self.jobs_by_status:
                print "\nFailed in processing"
                pprint(self.jobs_by_status['Processing'])

    def reset_failures(self):
        self.update_status_of_jobs()

        # find all jobs that need to be reset
        jobs_to_reset = []
        if 'Failed' in self.jobs_by_status:
            jobs_to_reset.extend(self.jobs_by_status['Failed'])
        if self.number_of_running_jobs() == 0:
            if 'Processing' in self.jobs_by_status:
                jobs_to_reset.extend(self.jobs_by_status['Processing'])

        # back up existing output
        for job in jobs_to_reset:
            match = re.search("^(.*?)\.job$", job)
            if match:
                fname = match.group(1)
                output = fname + ".root"
                lock = fname + ".lock"
                log = fname + ".log"
                if os.path.exists(output):
                    shutil.move(output, output + ".failed")
                if os.path.exists(lock):
                    shutil.move(lock, lock + ".failed")
                if os.path.exists(log):
                    shutil.move(log, log + ".failed")

                        
if __name__ == "__main__":
    # p = Skimmer("/eos/cms/store/group/phys_bphys/bmm/bmm5/tmp/NanoAOD-skims/510/ks/EGamma+Run2018B-17Sep2018-v1+MINIAOD/02d352594241e74ddcdbfeb18e2be0d0.job")
    # print p.__dict__
    # p.process()

    # vocms001 = SSHResourceHandler('vocms001.cern.ch', 16)
    # print vocms001.number_of_free_slots()
    # vocms001.submit_job("/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job")
    # vocms001.wait_for_all_jobs_completion()
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
    # jd.load_existing_jobs()
    # jd.list_failures()
    # jd.reset_failures()
    # jd.update_status_of_jobs()
    # jd.process_jobs()
