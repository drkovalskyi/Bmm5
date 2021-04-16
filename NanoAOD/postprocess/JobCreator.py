import postprocessing_cfg as cfg
import subprocess, re, json, hashlib, os

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
        
    def filename(self, input_files):
        """Generate unique file name based on the hash of input file names"""
        return hashlib.md5(",".join(input_files)).hexdigest()

    def create_new_jobs(self, make_small_jobs=False):
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
                    if not make_small_jobs and i + n >= n_elements:
                        break
                    # get input
                    inputs = new_inputs[i : i + n]
                    inputs.sort()

                    # create job unique id
                    job_id = self.filename(inputs)

                    # prepare job information
                    job_info = dict()
                    for key, value in task.items():
                        if key in ['files_per_job', 'type', 'name', 'input_pattern']:
                            continue
                        job_info[key] = value
                    job_info['input'] = []
                    for entry in inputs:
                        job_info['input'].append(cfg.xrootd_prefix + entry)
                    
                    # prepare output
                    job_dir = "%s/%s/%s/%s/%s" % (cfg.output_location, task['type'], cfg.version, 
                                                  task['name'], dataset)
                    if not os.path.exists(job_dir):
                        # os.makedirs(job_dir)
                        subprocess.call("mkdir -p %s" % job_dir, shell=True)
                    job_filename = "%s/%s.job" % (job_dir, job_id)
                    
                    # save job
                    json.dump(job_info, open(job_filename, "w"))
                    njobs += 1
                print "    Number of new jobs created %u" % njobs
                            
if __name__ == "__main__":

    jc = JobCreator()
    jc.find_all_inputs()
    jc.load_existing_jobs()
    jc.create_new_jobs(True)
