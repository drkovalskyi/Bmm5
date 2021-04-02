from PostProcessingBase import Processor
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import sys, os, subprocess

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

def unit_test():
    p = Skimmer("/eos/cms/store/group/phys_muon/dmytro/tmp/skim-test/1960fd1c81fb0d8371a3899fcf5cd36a.job")
    print p.__dict__
    p.process()

if __name__ == "__main__":
    unit_test()

