from PostProcessingBase import Processor
import time
import json
from ROOT import TChain, RDataFrame, TFile, std
import sys
import re


class SimpleSkimmer(Processor):
    """Processor to Skim and Slim files."""

    def _get_common_branches(self):
        # Find common branches preserving their order
        common_branches = None

        input_files = self.job_info['input']
        for file_name in input_files:
            root_file = TFile.Open(file_name)
            tree = root_file.Get("Events")
            branches = tree.GetListOfBranches()

            current_branches = []
            for j in range(branches.GetEntries()):
                branch_name = branches.At(j).GetName()
                current_branches.append(branch_name)

            if common_branches == None:
                common_branches = current_branches
            else:
                common_branches = [branch for branch in common_branches if branch in current_branches]

            root_file.Close()

        return common_branches


    def _process(self):
        t0 = time.time()

        # check for missing information
        for parameter in ['cut', 'input', 'keep']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

        ## get a list of common branches
        common_branches = self._get_common_branches()
        filtered_list = std.vector('string')()
        for column in common_branches:
            if re.search(self.job_info['keep'], str(column)):
                filtered_list.push_back(column)
        print(filtered_list)

        # setup input TChain
        input_files = self.job_info['input']
        chain = TChain("Events")
        for file in input_files:
            chain.Add(file)

        sys.stdout.flush()

        df = RDataFrame(chain)
        n_events = df.Count().GetValue()
        print("Number of events to process: %d" % n_events)

        df2 = df.Define("goodCandidates", self.job_info['cut'])
        dfFinal = df2.Filter("Sum(goodCandidates) > 0", "Event has good candidates")
        report = dfFinal.Report()

        dfFinal.Snapshot("Events", self.job_output_tmp, filtered_list)
        report.Print()

        print("Total time %.1f sec. to process %i events. Rate = %.1f Hz." % ((time.time() - t0), n_events, n_events / (time.time() - t0)))

        if 'verbose' in self.job_info and self.job_info['verbose']:
            file_size_input = 0
            file_size_output = 0

            for file in input_files:
                f = TFile.Open(file)
                file_size_input += f.GetSize()
                f.Close()

            f = TFile.Open(self.job_output_tmp)
            file_size_output += f.GetSize()
            f.Close()

            print("Input data size: %0.3f MB" % (file_size_input / 1e6))
            print("Output data size: %0.3f MB" % (file_size_output / 1e6))
            if file_size_output > 0:
                print("Reduction is size: %0.1f" % (float(file_size_input) / file_size_output))

        sys.stdout.flush()


def unit_test():
    standard_branches = 'PV_npvs|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock'

    ### create a test job
    job = {
        "input": [
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v4+MINIAODSIM/eeff8699-4ec6-4a6c-93a9-6df3db3992f8.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/522/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/48c83780-3bd9-44e2-848f-c05797f3d474.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022D-PromptReco-v1+MINIAOD/308d7ea2-c25d-47c2-a567-f7c4abd117af.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/60ef0541-b8f5-479b-becd-4fdbd0e0599b.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/610c8283-348f-485d-8ece-efe360e4a342.root"
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/610c8283-348f-485d-8ece-efe360e4a342.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/6158d285-e9a5-4a80-bb48-46e99684bc50.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/6181b52f-b6eb-4997-8ff9-db7ab0cfaac0.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61a97910-d8aa-4c1d-b42b-61f8f870778f.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61abcdb2-f5b3-4a9e-b41b-bb9435f99be0.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61b12c0c-7555-416b-874e-bd933fe24d7d.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61b17623-a7f8-4e6d-b28c-c24fd633e8ed.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61c53730-fdff-4ab7-b1f0-fa8cdf82bbc5.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61ca2301-d356-4df7-94be-27f648856583.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61d1e05b-b3c1-4795-aa40-1ce162be32df.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61d23b37-9970-4973-bd9d-fdfd4b739ae0.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61d37cc7-cf04-4df3-af64-70c0213fb8cf.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61e9f9f5-f1a7-436e-a04b-a636dc609f83.root",
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/521/EGamma+Run2022C-PromptReco-v1+MINIAOD/61f17f66-f6e3-490e-ad96-8c2a01ef7ecb.root",
        ],
        # "cut": "ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
        # "cut": "dstar_dm_pv > 0 ",
        "cut": "Muon_pt>0",
        "processor": "SimpleSkimmer",
        # "keep": "^(MuonId_.*|nMuonId|Muon_.*|nMuon)$",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|" + standard_branches + ")$",
        "verbose": True

    }

    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = SimpleSkimmer(file_name)
    p = SimpleSkimmer("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/524/em/ParkingBPH1+Run2018B-UL2018_MiniAODv2-v1+MINIAOD/039c7188cd7498a349a68a58f8ee392f.job")
    print(p.__dict__)
    p.process()


if __name__ == "__main__":
    unit_test()
