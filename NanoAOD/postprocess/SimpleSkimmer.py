from PostProcessingBase import Processor
import time
import json
import ROOT
from ROOT import TChain, RDataFrame, TFile, std, TTree
import sys
import re
import numpy as np
import os

class SimpleSkimmer(Processor):
    """Processor to Skim and Slim files."""

    def __init__(self, job_filename, take_ownership=False):
        self.n_gen_all = None
        self.n_gen_passed = None
        self.common_branches = None
        self.valid_files = []
        super(SimpleSkimmer, self).__init__(job_filename, take_ownership)

    @staticmethod
    def __get_branches(root_file):
        tree = root_file.Get("Events")
        branches = tree.GetListOfBranches()

        branch_names = []
        for j in range(branches.GetEntries()):
            branch_names.append(branches.At(j).GetName())

        return branch_names

    def _preprocess(self):
        input_files = self.job_info['input']

        lumi_mask_type = None
        if 'lumi_mask' in self.job_info:
            lumi_mask_type = self.job_info['lumi_mask']

        for file_name in input_files:
            root_file = TFile.Open(file_name)
            skip_file = False
            if lumi_mask_type:
                skip_file = True

            # GenFilterInfo
            lumis = root_file.Get("LuminosityBlocks")
            if lumis:
                for lumi in lumis:
                    if lumi_mask_type:
                        if not self._is_certified_run_lumi(lumi.run, lumi.luminosityBlock, lumi_mask_type):
                            continue
                        skip_file = False
                    if hasattr(lumi, 'GenFilter_numEventsPassed'):
                        if self.n_gen_all == None:
                            self.n_gen_all = 0
                            self.n_gen_passed = 0
                        self.n_gen_passed += lumi.GenFilter_numEventsPassed
                        self.n_gen_all    += lumi.GenFilter_numEventsTotal

            if not skip_file:
                # Find common branches preserving their order
                current_branches = self.__get_branches(root_file)

                if self.common_branches == None:
                    self.common_branches = current_branches
                else:
                    self.common_branches = [
                        branch
                        for branch in self.common_branches
                        if branch in current_branches
                    ]
                self.valid_files.append(file_name)
            else:
                print(f"Ignore {file_name} - not certified")

            if 'save_branches' in self.job_info and self.job_info['save_branches'] == True:
                branch_info_file_name = re.sub(r'\.root$', '.branches', file_name)
                branch_info_file_name = re.sub(r'^.*?\/eos\/cms', '/eos/cms', branch_info_file_name)
                branch_names = self.__get_branches(root_file)
                branch_names.sort()
                with open(branch_info_file_name, 'w') as file:
                    json.dump(branch_names, file, indent=4)
                    
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

        # preprocess
        self._preprocess()
        if len(self.valid_files) == 0:
            print("No valid input selected. Write empty ROOT file.")
            f = TFile.Open(self.job_output_tmp, "recreate")
            f.Close()
            return
        
        ## get a list of common branches
        filtered_list = std.vector('string')()
        for column in self.common_branches:
            if re.search(self.job_info['keep'], str(column)):
                filtered_list.push_back(column)
        print(filtered_list)

        # setup input TChain
        chain = TChain("Events")
        for file in self.valid_files:
            chain.Add(file)

        sys.stdout.flush()

        df = RDataFrame(chain)
        if 'lumi_mask' in self.job_info:
            self._declare_lumi_mask_code()
            lumi_mask = self._get_lumi_mask(self.job_info['lumi_mask'])
            ROOT.gInterpreter.ProcessLine(f'lumi_mask_string = "{lumi_mask}";')
            df = df.Define("certified", "passed_lumi_mask(run, luminosityBlock)")
            df = df.Filter("certified == 1", "Passed data certification")
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
    # job = {
    #     "input": [
    #         '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/526/ParkingDoubleMuonLowMass0+Run2022C-PromptReco-v1+MINIAOD/5114a593-fa24-4fc4-a475-77e6312c0362.root'
    #     ],
    #     # "cut": "ks_kin_sipPV<3 && ks_kin_slxy>3 && ks_trk1_sip>5 && ks_trk2_sip>5 && ks_kin_cosAlphaXY>0.999",
    #     # "cut": "dstar_dm_pv > 0 ",
    #     "cut": "HLT_Mu4_L1DoubleMu",
    #     "processor": "SimpleSkimmer",
    #     # "keep": "^(MuonId_.*|nMuonId|Muon_.*|nMuon)$",
    #     "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|HLT_.*|" + standard_branches + ")$",
    #     "candidate_loop": False,
    #     "verbose": True

    # }
    job = {
        # "save_branches": true,
        "lumi_mask": "muon",
        "processor": "SimpleSkimmer",
        "cut": "mm_mass > 0",
        "keep": "^(mm_.*|nmm|Muon_.*|nMuon|MuonId_.*|nMuonId|npvs|pvs_.*|HLT_DoubleMu4_3_LowMass|HLT_Mu0_L1DoubleMu|L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4|L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6|L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6|L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5|L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4|L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4|L1_DoubleMu4p5_SQ_OS_dR_Max1p2|L1_DoubleMu4_SQ_OS_dR_Max1p2|PV_npvs|PV_npvsGood|Pileup_nTrueInt|Pileup_nPU|run|event|luminosityBlock)$",
        "input": [
            # "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/HLTPhysics+Run2022F-PromptReco-v1+MINIAOD/b89c0db7-0d8f-420a-9888-87391958106e.root",
            "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/HLTPhysics+Run2022F-PromptReco-v1+MINIAOD/b8aafd1d-54ca-4c4f-bcc0-d901d0b632a0.root"
        ]
    }
    
    file_name = "/tmp/dmytro/test.job"
    json.dump(job, open(file_name, "w"))

    # p = SimpleSkimmer(file_name)
    p = SimpleSkimmer("/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/529/mm/HLTPhysics+Run2022F-PromptReco-v1+MINIAOD/f257b748b8665d04d75551aa22fa7691.job")
    print(p.__dict__)
    p.process()


if __name__ == "__main__":
    unit_test()
