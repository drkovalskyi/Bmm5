from PostProcessingBase import Processor
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import sys, os, subprocess, json

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

        sys.stdout.flush()

        # merge skimed data
        skimmed_files = []
        for f in input_files:
            skimmed_files.append(os.path.join(self.tmp_dir,
                                              os.path.basename(f).replace(".root","%s.root" % postfix)))
        if len(skimmed_files) > 1:
            subprocess.call("haddnano.py %s %s" % (self.job_output_tmp, " ".join(skimmed_files)), shell=True)
            # clean up
            for f in skimmed_files:
                os.remove(f)
        else:
            if len(skimmed_files) == 0:
                raise Exception("There should be more than one skimmed file")
            subprocess.call("mv -v %s %s" % (skimmed_files[0], self.job_output_tmp), shell=True)

        sys.stdout.flush()


def unit_test():
    ### create a test job
    job = {
        "input": [
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD//515/EGamma+Run2018B-17Sep2018-v1+MINIAOD/EC8FE87A-ADBA-8F40-A4F4-C6C28431BB91.root"
        ],
        "cut": "ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
        "processor": "Skimmer"
    }
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    p = Skimmer(file_name)
    print p.__dict__
    p.process()

if __name__ == "__main__":
    unit_test()

