from PostProcessingBase import Processor
import time
import json
from ROOT import TChain, RDataFrame, TFile, std, TTree
import sys
import re
import numpy as np


class SimpleSkimmer(Processor):
    """Processor to Skim and Slim files."""

    def __init__(self, job_filename, take_ownership=False):
        self.n_gen_all = None
        self.n_gen_passed = None
        self.common_branches = None
        super(SimpleSkimmer, self).__init__(job_filename, take_ownership)

    def _preprocess(self):
        input_files = self.job_info['input']
        for file_name in input_files:
            root_file = TFile.Open(file_name)

            # GenFilterInfo
            lumis = root_file.Get("LuminosityBlocks")
            if lumis:
                for lumi in lumis:
                    if hasattr(lumi, 'GenFilter_numEventsPassed'):
                        if self.n_gen_all == None:
                            self.n_gen_all = 0
                            self.n_gen_passed = 0
                        self.n_gen_passed += lumi.GenFilter_numEventsPassed
                        self.n_gen_all    += lumi.GenFilter_numEventsTotal

            # Find common branches preserving their order
            tree = root_file.Get("Events")
            branches = tree.GetListOfBranches()

            current_branches = []
            for j in range(branches.GetEntries()):
                branch_name = branches.At(j).GetName()
                current_branches.append(branch_name)

            if self.common_branches == None:
                self.common_branches = current_branches
            else:
                self.common_branches = [branch for branch in self.common_branches if branch in current_branches]

            root_file.Close()


    def _process(self):
        t0 = time.time()

        # set default value to keep old jobs functional
        if 'candidate_loop' not in self.job_info:
            self.job_info['candidate_loop'] = True

        # check for missing information
        for parameter in ['cut', 'input', 'keep', 'candidate_loop']:
            if parameter not in self.job_info:
                raise Exception("Missing input '%s'" % parameter)

        ## get a list of common branches
        self._preprocess()
        filtered_list = std.vector('string')()
        for column in self.common_branches:
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

        if self.job_info['candidate_loop'] == True:
            df2 = df.Define("goodCandidates", self.job_info['cut'])
            dfFinal = df2.Filter("Sum(goodCandidates) > 0", "Event has good candidates")
        else:
            dfFinal = df.Filter(self.job_info['cut'], "Passed selection")
            
        report = dfFinal.Report()

        if 'keep_only_common_branches' in self.job_info and self.job_info['keep_only_common_branches']:
            print("WARNING: keeping only common branches in the output. May lead to data loss")
            dfFinal.Snapshot("Events", self.job_output_tmp, filtered_list)
        else:
            dfFinal.Snapshot("Events", self.job_output_tmp, self.job_info['keep'])
            
        report.Print()

        # reports

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

        # Update accounting information

        # Store meta data
        f = TFile(self.job_output_tmp, "UPDATE")
        t = TTree("info","Selection information")

        n_processed = np.empty((1), dtype="i")
        t.Branch("n_processed", n_processed, "n_processed/I")
        n_processed[0] = n_events

        if self.n_gen_all != None:
            n_gen_all = np.empty((1), dtype="u8")
            n_gen_passed = np.empty((1), dtype="u8")
            t.Branch("n_gen_all", n_gen_all, "n_gen_all/l")
            t.Branch("n_gen_passed", n_gen_passed, "n_gen_passed/l")
            n_gen_all[0] = self.n_gen_all
            n_gen_passed[0] = self.n_gen_passed

        t.Fill()
        f.Write()
        f.Close()


def unit_test():
    standard_branches = 'PV_npvsGood|PV_npvs|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock'

    ### create a test job
    job = {
        "input": [
            '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/526/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/5114a593-fa24-4fc4-a475-77e6312c0362.root'
        ],
        # "cut": "ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
        # "cut": "dstar_dm_pv > 0 ",
        "cut": "HLT_Mu4_L1DoubleMu",
        "processor": "SimpleSkimmer",
        # "keep": "^(MuonId_.*|nMuonId|Muon_.*|nMuon)$",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|HLT_.*|" + standard_branches + ")$",
        "candidate_loop": False,
        "verbose": True

    }

    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = SimpleSkimmer(file_name)
    p = SimpleSkimmer("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/526/trig/Muon+Run2022E-PromptReco-v1+MINIAOD/04011e9bd966e9590d0df1f91f3c5b40.job")
    print(p.__dict__)
    p.process()


if __name__ == "__main__":
    unit_test()
