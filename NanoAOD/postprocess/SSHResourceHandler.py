from PostProcessingBase import ResourceHandler
import sys, os, subprocess, re, fcntl, time, tempfile
import postprocessing_cfg as cfg

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
        # print "%s %s " % (self.name(), command)
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
        command = "nice bash job_starter.sh %s %s &> /dev/null &" % (self._processor_name(job), job)
        print "submitting %s" % job
            
        print self._send_command_and_get_response(command),

    def _get_running_jobs(self):
        jobs = []
        response = self._send_command_and_get_response("ps -Af | grep '[j]ob_starter.sh'")
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
            match = re.search('^\S+\s+(\S+).*?job_starter.sh', line)
            if match:
                self._send_command_and_get_response("kill %s" % match.group(1))

    def clean_up(self):
        command = "find %s -path '*%s*' -exec rm -rfv {} \;" % (tempfile.gettempdir(), cfg.tmp_prefix)
        self._send_command_and_get_response(command)

def unit_test():
    test_job = "/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job"
    res = SSHResourceHandler('vocms0314.cern.ch', 16)
    res.submit_job(test_job)
    res.wait_for_jobs_to_finish()

if __name__ == "__main__":
    unit_test()
