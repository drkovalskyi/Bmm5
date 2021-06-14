from PostProcessingBase import ResourceHandler
import sys, os, subprocess, re, fcntl, time, tempfile
import postprocessing_cfg as cfg

class LocalResourceHandler(ResourceHandler):
    """Resource handler for local job execution"""
    def __init__(self, max_number_of_jobs_running):
        super(LocalResourceHandler, self).__init__()
        self.max_njobs = max_number_of_jobs_running

    def _submit_job(self, job):
        log = re.sub('\.job$', '\.log', job)
        command = "nice bash job_starter.sh %s %s >& /dev/null &" % (self._processor_name(job), job)
        print "submitting %s" % job
        subprocess.call(command, shell=True)
        
    def _get_running_jobs(self):
        jobs = []
        response = subprocess.check_output("ps -Af", shell=True, env={})
        for line in response.splitlines():
            if not re.search('job_starter.sh', line): continue
            match = re.search('job_starter.sh.*?(\S+\.job)', line)
            if match:
                jobs.append(match.group(1))
        return jobs

    def name(self):
        return "localhost"

    def kill_all_jobs(self):
        print "Killing jobs at %s" % self.name()
        response = subprocess.check_output("ps -Af", shell=True)
        for line in response.splitlines():
            match = re.search('^\S+\s+(\S+).*?job_starter', line)
            if match:
                subprocess.call("kill %s" % match.group(1), shell=True)
                
    def clean_up(self):
        # command = "find %s -path '*%s*' -exec rm -rfv {} \;" % (tempfile.gettempdir(), cfg.tmp_prefix)
        command = "find /tmp/ -type d -name '*%s*' -exec rm -rfv {} \;" % (cfg.tmp_prefix)
        subprocess.call(command, shell=True)
                            
def unit_test():
    lh = LocalResourceHandler(16)
    lh.submit_job("/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job")
    # lh.submit_job("/eos/cms/store/group/phys_muon/dmytro/tmp/mva_ntuple_test/ffcb48b1f605ca5a5f2614b94830e806.job")
    print lh.get_running_jobs()
    print lh.number_of_free_slots()
            
if __name__ == "__main__":
    unit_test()
