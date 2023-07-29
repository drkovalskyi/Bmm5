from PostProcessingBase import Processor
from Bmm5.NanoAOD.postprocessing.postprocessor import PostProcessor
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
                                  compression="ZLIB:1",                
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
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022D-PromptReco-v1+MINIAOD/308d7ea2-c25d-47c2-a567-f7c4abd117af.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/60ef0541-b8f5-479b-becd-4fdbd0e0599b.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/610c8283-348f-485d-8ece-efe360e4a342.root"
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/610c8283-348f-485d-8ece-efe360e4a342.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/6158d285-e9a5-4a80-bb48-46e99684bc50.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/6181b52f-b6eb-4997-8ff9-db7ab0cfaac0.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61a97910-d8aa-4c1d-b42b-61f8f870778f.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61abcdb2-f5b3-4a9e-b41b-bb9435f99be0.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61b12c0c-7555-416b-874e-bd933fe24d7d.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61b17623-a7f8-4e6d-b28c-c24fd633e8ed.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61c53730-fdff-4ab7-b1f0-fa8cdf82bbc5.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61ca2301-d356-4df7-94be-27f648856583.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61d1e05b-b3c1-4795-aa40-1ce162be32df.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61d23b37-9970-4973-bd9d-fdfd4b739ae0.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61d37cc7-cf04-4df3-af64-70c0213fb8cf.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61e9f9f5-f1a7-436e-a04b-a636dc609f83.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61f17f66-f6e3-490e-ad96-8c2a01ef7ecb.root",
        ],
        "cut": "ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
        "processor": "Skimmer"
    }
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    p = Skimmer(file_name)
    print(p.__dict__)
    p.process()

if __name__ == "__main__":
    unit_test()

