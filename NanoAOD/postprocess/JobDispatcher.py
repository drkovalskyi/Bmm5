import postprocessing_cfg as cfg
import time, subprocess, re, os, shutil, json
from LocalResourceHandler import LocalResourceHandler
from SSHResourceHandler import SSHResourceHandler
from Skimmer import Skimmer
from FlatNtupleForBmmMva import FlatNtupleForBmmMva
from FlatNtupleForMLFit import FlatNtupleForMLFit
from pprint import pprint

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

    def _job_info(self, job):
        match = re.search("^(.*?)\.job$", job)
        info = dict()
        if match:
            fname = match.group(1)
            info['output'] = fname + ".root"
            info['lock'] = fname + ".lock"
            info['log'] = fname + ".log"
            info['summary'] = fname + ".summary"
        else:
            raise Exception("Incorrect job name:\n%s" % job)
        return info
            
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
        job_info = self._job_info(job)
        if os.path.exists(job_info['lock']):
            if job not in self.running_jobs:
                return "Failed"
            else:
                return "Running"
        elif os.path.exists(job_info['output']):
            if cfg.require_log_for_success:
                if not os.path.exists(job_info['log']):
                    return 'Failed'
            return "Done"
        elif os.path.exists(job_info['log']):
            if job not in self.running_jobs:
                return "Failed"
            else:
                return "Running"
        else: 
            return "New"
            
    def _load_existing_jobs(self):
        """Find existings jobs and store their input"""
        print "Loading existing jobs..."
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

    def number_of_running_jobs(self, owned=False):
        """Get total number of running jobs on all resources"""
        n_running = 0
        for resource in self.resources():
            n_running += resource.number_of_running_jobs(owned)
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
                if detailed:
                    print "\n", job
                     
                job_info = self._job_info(job)

                failure_type = None
                if os.path.exists(job_info['output']):
                    failure_type = 'Output is available'
                elif os.path.exists(job_info['lock']):
                    failure_type = 'Locked without output'
                elif os.path.exists(job_info['log']):
                    failure_type = 'Only log'
                else:
                    failure_type = 'Missing log'

                if detailed and os.path.exists(job_info['log']):
                    subprocess.call("tail %s" % job_info['log'], shell=True)
                    
                if failure_type:
                    if failure_type not in failures:
                        failures[failure_type] = []
                    failures[failure_type].append(job)

        for failure_type, jobs in failures.items():
            print "Failure type: %s" % failure_type
            print "\tNumber of jobs: %u" % len(jobs)
            pprint(jobs)

    def reset_failures(self):
        self.update_status_of_jobs()

        # find all jobs that need to be reset
        jobs_to_reset = []
        if 'Failed' in self.jobs_by_status:
            jobs_to_reset.extend(self.jobs_by_status['Failed'])

        # back up existing output
        for job in jobs_to_reset:
            job_info = self._job_info(job)
    
            for f in [job_info['output'], job_info['lock'], job_info['log']]:
                # remove previous backups
                if os.path.exists(f + ".failed"):
                    os.remove(f + ".failed")
                # backup latest failure
                if os.path.exists(f):
                    shutil.move(f, f + ".failed")

    def kill_all_jobs(self):
        for resource in self.resources():
            resource.kill_all_jobs()

    def clean_up(self):
        for resource in self.resources():
            resource.clean_up()

    def get_job_summary(self, job, reanalyze=False):
        """Analyze job log and produce summary"""
        job_info = self._job_info(job)
        if (reanalyze or not os.path.exists(job_info['summary'])):
            n_selected = 0
            n_processed = 0
            rate = None
            # Skimmer specific analysis for now
            result = subprocess.check_output("grep -E 'Selected|Hz' %s" % job_info['log'],
                                             shell=True)
            for line in result.splitlines():
                match = re.search("^Selected\s+(\d+)[\s\/]+(\d+)\s+entries", line)
                if match:
                    n_selected += int(match.group(1))
                    n_processed += int(match.group(2))
                match = re.search("^([\d\.]+)\s+Hz", line)
                if match:
                    rate = float(match.group(1))
            report = {
                'n_selected':n_selected, 'n_processed':n_processed, 'rate':rate
            }
            json.dump(report, open(job_info['summary'], 'w'))
            return report
        else:
            report = json.load(open(job_info['summary']))
            return report

    def format_data_for_datatable(self, data, column_names, column_types, formatting):
        result = dict()
        result['cols'] = []
        result['rows'] = []
        for i in range(len(column_names)):
            result['cols'].append({"label":column_names[i], "type":column_types[i]})
        for row in data:
            cell = []
            for i in range(len(row)):
                if formatting[i]:
                    cell.append({'v':row[i], 'f':formatting[i] % row[i]})
                else:
                    cell.append({'v':row[i]})
            result['rows'].append({'c':cell})
        return result
        
    def job_report(self):
        """Extract job progress and efficiency information and publish it"""
        self.update_status_of_jobs()

        report = dict()
        if 'Done' in self.jobs_by_status:
            for job in self.jobs_by_status['Done']:
                match = re.search("([^\/]+)\/(\d+)\/([^\/]+)\/([^\/]+)", job)
                if match:
                    task_type = match.group(1)
                    if task_type not in report:
                        report[task_type] = dict()
                    version = match.group(2)
                    if version not in report[task_type]:
                        report[task_type][version] = dict()
                    task_name = match.group(3)
                    if task_name not in report[task_type][version]:
                        report[task_type][version][task_name] = dict()
                    dataset = match.group(4)
                    if dataset not in report[task_type][version][task_name]:
                        report[task_type][version][task_name][dataset] = []
                    try:
                        report[task_type][version][task_name][dataset].append(self.get_job_summary(job))
                    except:
                        pass

        report_template ="""
<html>
  <head>
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
    <script type='text/javascript' src='https://www.google.com/jsapi'></script>
    <style>
     .google-visualization-table-table *  { font-size:11px; } 
    </style>
    <script src="https://www.gstatic.com/charts/loader.js"></script>
    <script>
     google.charts.load('current', {
       callback: draw,
       packages:['table']
     });
     function draw() {
       $.getJSON("%s", 
		 function drawChart(js_obj) {
                    var data = new google.visualization.DataTable(js_obj);
                    var table = new google.visualization.Table(document.getElementById('table_div'));
                    table.draw(data);
                 });
       }
    </script>
  </head>
  <body>
    <div id="table_div"></div>
  </body>
</html>
"""

        # pprint(report)
        for task_type in report:
            for version in report[task_type]:
                for task_name in report[task_type][version]:
                    path = "%s/%s/%s/" % (cfg.web_report_path, task_type, version)
                    if not os.path.exists(path):
                        subprocess.call("mkdir -p %s" % path, shell=True)
                    # with open("%s/%s.txt" % (path, task_name), "w") as f:
                    data = []
                    for dataset, info in report[task_type][version][task_name].items():
                        n_processed = 0
                        n_selected = 0
                        n_rate = 0
                        sum_rate = 0
                        for entry in info:
                            if 'n_processed' in entry and 'n_selected' in entry:
                                n_processed += entry['n_processed']
                                n_selected+= entry['n_selected']
                            if 'rate' in entry:
                                sum_rate += entry['rate']
                                n_rate += 1 
                        if n_processed > 0:
                            average_rate = 0
                            if n_rate > 0:
                                average_rate = float(sum_rate)/n_rate
                            # f.write("%s %u %0.1f%% %0.1f Hz\n" % (dataset, n_selected, 100.*n_selected/n_processed, average_rate))
                            data.append([dataset, n_selected, 100.*n_selected/n_processed, average_rate])
                    data_table = self.format_data_for_datatable(data,
                                                                ["Dataset","Selected", "Efficiency", "Rate, Hz"],
                                                                ["string", "number", "number", "number"],
                                                                [None, None, "%0.1f%%", "%0.1f"])
                    json.dump(data_table, open("%s/%s.json" % (path, task_name), "w"))
                    with open("%s/%s.html" % (path, task_name), "w") as f:
                        f.write(report_template % (task_name + ".json"))
                            
if __name__ == "__main__":
    jd = JobDispatcher()
    # jd.kill_all_jobs()
    # jd.show_resource_availability()
    # jd.show_failures(True)
    # jd.clean_up()
    # jd.update_status_of_jobs()
    # jd.job_report()
    
    jd.show_failures()

    # jd.reset_failures()
    # jd.process_jobs()
