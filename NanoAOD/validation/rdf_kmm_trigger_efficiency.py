"""Trigger efficiency study

The script measures the relative and absolute trigger efficiency of
the dimuon trigger for Ks->mumu analysis in MC and Data.

The input data must be in Bmm5 NanoAOD format.

The script uses RDataFrame for optimal performance, but it needs to
read very large volumes of data to find events from pre-scaled
triggers. Therefore it is advisable to create a trigger skim for data
and use only small number of files to test on raw NanoAOD samples.

The results are stored in two two ways:
- json file with efficiency measurements
- histograms of PV distributions to track conditions of triggers with 
  dynamic prescales

The json file is used in an update mode, i.e. new results override old
results for specific studies keeping results from other studies
unaffected.

"""

import sys
import ROOT
if not hasattr(ROOT.RDataFrame, "DefaultValueFor"):
    sys.exit(f"Error: ROOT.RDataFrame.DefaultValueFor is not available in ROOT {ROOT.gROOT.GetVersion()}. Please use ROOT version 6.34 or later.")
from collections import defaultdict
from pprint import pprint
import os
import re
import glob
import json
import math

ROOT.ROOT.EnableImplicitMT()
# max_files = 10
max_files = 999999
sample_names = ["K0sToMuMu", "ParkingDoubleMuonLowMass", "Muon", "ZeroBias",
                "EGamma", "ParkingDoubleElectronLowMass", "HLTPhysics"]
# sample_names = ["K0sToMuMu"]
recompute_results = True

path   = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/"
path1  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/532/trig/"
path2  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing-NEW/Skims/532/mm/"
path3  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/529/mm"
output = "results/kmm_summary.json"

def load_data(file_patterns):
    chain = ROOT.TChain("Events")
    files = []
    for file_pattern in file_patterns:
        files.extend(glob.glob(file_pattern))
    print("Number of files:", len(files))
    if len(files) == 0:
        return None
    for file in files[:max_files]:
        chain.Add(file)
    n_files = chain.GetListOfFiles().GetEntries()
    print("Number of files in the chain:", n_files)
    if n_files == 0:
        return None
    return chain


def process_sample(sample, era):
    proxies = dict()
    
    file_pattern = []
    for sample_era, pattern in sample["files"].items():
        if sample_era != era:
            continue
        file_pattern.extend(pattern)
    if len(file_pattern) == 0:
        return dict()

    run = None
    if sample["Data"]:
        for i_run, i_eras in eras.items():
            if era in i_eras:
                run = i_run
                break
    else:
        for i_run, i_campaigns in campaigns.items():
            if era in i_campaigns:
                run = i_run
                break
    
    pprint(file_pattern)
    
    chain = load_data(file_pattern)
    if chain == None:
        return dict()
            
    rdf = ROOT.RDataFrame(chain)

    # Make sure all triggers have a default value
    triggers = [
        # HLT
        "HLT_DoubleMu4_3_LowMass", "HLT_Mu4_L1DoubleMu", "HLT_Mu0_L1DoubleMu",
        "HLT_Mu3_PFJet40", "HLT_Mu0_L1DoubleMu", "HLT_Mu8",
        # Run2022
        "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", "L1_DoubleMu4_SQ_OS_dR_Max1p2",
        # Run2023
        "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6",
        "L1_DoubleMu4_SQ_OS_dR_Max1p2",
        # Run2024
        "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
        "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"
    ]
    for trigger in triggers:
        rdf = rdf.DefaultValueFor(trigger, False)

    # Extract information from other branches
    rdf = rdf.Define("mm_mu1_mediumId", "Take(Muon_mediumId,             mm_mu1_index)")
    rdf = rdf.Define("mm_mu1_npixels",  "Take(MuonId_nPixels,            mm_mu1_index)")
    rdf = rdf.Define("mm_mu2_mediumId", "Take(Muon_mediumId,             mm_mu2_index)")
    rdf = rdf.Define("mm_mu2_npixels",  "Take(MuonId_nPixels,            mm_mu2_index)")
        
    # Offline selection
    rdf = rdf.Define("cands", "mm_mu1_index>=0 && mm_mu2_index>=0 && "
                     "mm_mu1_pt > 4 && mm_mu2_pt > 3 && "
                     "(mm_mu1_npixels > 0 && mm_mu2_npixels > 0) &&" # trigger requirement
                     # " mm_kin_lxy > 1 &&"
                     # "mm_mu1_pt > 4 && mm_mu2_pt > 4 &&"
                     # "abs(mm_mu1_eta)<1.4 && abs(mm_mu2_eta)<1.4 &&"
                     # "(abs(mm_mu1_eta)>1.4 || abs(mm_mu2_eta)>1.4) &&"
                     "mm_mu1_mediumId && mm_mu2_mediumId && "
                     "mm_mu1_pdgId * mm_mu2_pdgId == -169 && "
                     "abs(mm_kin_mass-0.50)<0.15 && mm_kin_vtx_prob>0.01")
    
    rdf = rdf.Filter("Sum(cands)>0")
    
    proxies["All"] = rdf.Count()

    if run == "Run2022":
        l1_selection_unprescaled = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
        l1_selection_all = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
    elif run == "Run2023":
        l1_selection_unprescaled = "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4||L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6" + \
            "||L1_DoubleMu4_SQ_OS_dR_Max1p2"
        # l1_selection_all = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6" + \
        #     "||L1_DoubleMu5_SQ_OS_dR_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
        l1_selection_all = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
    elif run == "Run2024":
        l1_selection_unprescaled = "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4||L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6" + \
            "||L1_DoubleMu4p5_SQ_OS_dR_Max1p2"
        # l1_selection_all = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6" + \
        #     "||L1_DoubleMu5_SQ_OS_dR_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
        l1_selection_all = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
    
    ref_trigger = "HLT_Mu4_L1DoubleMu"
    if sample["triggers"] == None or ref_trigger in sample["triggers"]:
        rdf1 = rdf.Filter(ref_trigger)
        proxies[ref_trigger] = rdf1.Count()
        proxies[f"{ref_trigger}&&HLT_DoubleMu4_3_LowMass"] = rdf1.Filter("HLT_DoubleMu4_3_LowMass").Count()
        
    ref_trigger = "HLT_Mu3_PFJet40"
    if sample["triggers"] == None or ref_trigger in sample["triggers"]:
        rdf2 = rdf.Filter(ref_trigger)
        proxies[ref_trigger] = rdf2.Count()
        proxies[f"{ref_trigger}&&HLT_Mu0_L1DoubleMu"] = rdf2.Filter("HLT_Mu0_L1DoubleMu").Count()
        proxies[f"{ref_trigger}&&HLT_Mu0_L1DoubleMu&&L1_DoubleMu_Mix_Unprescaled"] = \
            rdf2.Filter(f"HLT_Mu0_L1DoubleMu&&({l1_selection_unprescaled})").Count()

    ref_trigger = "HLT_Mu8"
    if sample["triggers"] == None or ref_trigger in sample["triggers"]:
        rdf3 = rdf.Filter(ref_trigger)
        proxies[ref_trigger] = rdf3.Count()
        proxies[f"{ref_trigger}&&HLT_Mu0_L1DoubleMu"] = rdf3.Filter("HLT_Mu0_L1DoubleMu").Count()
        proxies[f"{ref_trigger}&&HLT_Mu0_L1DoubleMu&&L1_DoubleMu_Mix_Unprescaled"] = \
            rdf3.Filter(f"HLT_Mu0_L1DoubleMu&&({l1_selection_unprescaled})").Count()

    # no reference
    if sample["triggers"] == None:
        proxies[f"HLT_Mu0_L1DoubleMu"] = rdf.Filter("HLT_Mu0_L1DoubleMu").Count()
        proxies[f"HLT_DoubleMu4_3_LowMass"] = rdf.Filter("HLT_DoubleMu4_3_LowMass").Count()

        # L1
        l1_unprescaled_rdf = rdf.Filter(l1_selection_unprescaled)
        proxies["L1_DoubleMu_Mix_Unprescaled"] = l1_unprescaled_rdf.Count()
        proxies[f"L1_DoubleMu_Mix_Unprescaled&&HLT_DoubleMu4_3_LowMass"] = \
            l1_unprescaled_rdf.Filter("HLT_DoubleMu4_3_LowMass").Count()
        
        l1_all_rdf = rdf.Filter(l1_selection_all)
        proxies["L1_DoubleMu_Mix_All"] = l1_all_rdf.Count()
        proxies[f"L1_DoubleMu_Mix_All&&HLT_DoubleMu4_3_LowMass"] = \
            l1_all_rdf.Filter("HLT_DoubleMu4_3_LowMass").Count()
    
    report = dict()
    for name, proxy in proxies.items():
        report[name] = proxy.GetValue()
    return report


def aggregate_results(run):
    report = dict()
    for sample_name, info in samples.items():
        if sample_name not in results:
            continue
        report[sample_name] = defaultdict(int)
        if info["Data"]:
            for era in results[sample_name]:
                if era in eras[run]:
                    for selection, n in results[sample_name][era].items():
                        report[sample_name][selection] += n
        else:
            for campaign in results[sample_name]:
                if campaign in campaigns[run]:
                    for selection, n in results[sample_name][campaign].items():
                        report[sample_name][selection] += n
    return report

def get_efficiency(n, n0):
    if n0 == 0:
        return ""
    eff = n / float(n0)
    eff_err = math.sqrt(eff * (1 - eff) / n0)

    return  f"{eff * 100:0.1f} \pm {eff_err * 100:0.1f}"
    
        
def print_results(run):
    report = aggregate_results(run)

    eff_report = dict()
    eff_report = {
        "Absolute efficiency": {
            "Data": {},
            "MC": {
                "K0sToMuMu": get_efficiency(report["K0sToMuMu"]["HLT_DoubleMu4_3_LowMass"],
                                            report["K0sToMuMu"]["All"]),
            }
        },
        "Absolute HLT efficiency": {
            "Data": {},
            "MC": {
                "K0sToMuMu": get_efficiency(report["K0sToMuMu"]["L1_DoubleMu_Mix_Unprescaled&&HLT_DoubleMu4_3_LowMass"],
                                            report["K0sToMuMu"]["L1_DoubleMu_Mix_Unprescaled"]),
            }
        },
        "Relative HLT efficiency": {
            "Data": {},
            "MC": {
                "K0sToMuMu": get_efficiency(report["K0sToMuMu"]["HLT_Mu4_L1DoubleMu&&HLT_DoubleMu4_3_LowMass"],
                                            report["K0sToMuMu"]["HLT_Mu4_L1DoubleMu"]),
            }
        },
        "Absolute L1 efficiency": {
            "Data": {},
            "MC": {
                "K0sToMuMu": get_efficiency(report["K0sToMuMu"]["L1_DoubleMu_Mix_Unprescaled"],
                                            report["K0sToMuMu"]["All"]),
            }
        },
        "Relative L1 efficiency": {
            "Data": {},
            "MC": {
                "K0sToMuMu": get_efficiency(report["K0sToMuMu"]["HLT_Mu8&&HLT_Mu0_L1DoubleMu"],
                                            report["K0sToMuMu"]["HLT_Mu8"]),
            }
        },
        
    }
    for pd in ["HLTPhysics", "EGamma"]:
        eff_report["Absolute HLT efficiency"]["Data"][pd] = \
            get_efficiency(report[pd]["L1_DoubleMu_Mix_Unprescaled&&HLT_DoubleMu4_3_LowMass"],
                           report[pd]["L1_DoubleMu_Mix_Unprescaled"])
    for pd in ["ParkingDoubleMuonLowMass"]:
        eff_report["Relative HLT efficiency"]["Data"][pd] = \
            get_efficiency(report[pd]["HLT_Mu4_L1DoubleMu&&HLT_DoubleMu4_3_LowMass"],
                           report[pd]["HLT_Mu4_L1DoubleMu"])
    for pd in ["EGamma", "ZeroBias", "ParkingDoubleElectronLowMass"]:
        eff_report["Absolute L1 efficiency"]["Data"][pd] = \
            get_efficiency(report[pd]["L1_DoubleMu_Mix_Unprescaled"],
                           report[pd]["All"])
    for pd in ["Muon"]:
        eff_report["Relative L1 efficiency"]["Data"][pd] = \
            get_efficiency(report[pd]["HLT_Mu8&&HLT_Mu0_L1DoubleMu"],
                           report[pd]["HLT_Mu8"])
    for pd in ["EGamma"]:
        eff_report["Absolute efficiency"]["Data"][pd] = \
            get_efficiency(report[pd]["HLT_DoubleMu4_3_LowMass"],
                           report[pd]["All"])
            
    
    
    pprint(report)
    pprint(eff_report)

        

####################################################################

results = dict()
if not os.path.exists('results'):
    os.mkdir('results')
if os.path.exists(output):
    results = json.load(open(output))
	
samples = {
    "ParkingDoubleMuonLowMass": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": ["HLT_Mu4_L1DoubleMu"]
    },
    "Muon": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": ["HLT_Mu3_PFJet40", "HLT_Mu8"]
    },
    "SingleMuon": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": ["HLT_Mu3_PFJet40"]
    },
    "DoubleMuon": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": ["HLT_Mu8"]
    },
    "ZeroBias": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": None
    },
    "EGamma": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": None
    },
    "ParkingDoubleElectronLowMass": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": None
    },
    "HLTPhysics": {
        "Data": True,
        "files": defaultdict(list),
        "triggers": None
    },
}

defaultdict(list)
eras = {
    "Run2022": {
        "Run2022C": 1,
        "Run2022D": 2,
        "Run2022E": 1,
        "Run2022F": 1,
    },
    "Run2023": {
        "Run2023C": 4,
        "Run2023D": 2,
    },
    "Run2024": {
        "Run2024C": 1,
        "Run2024D": 1,
        "Run2024E": 2,
        "Run2024F": 1,
        "Run2024G": 1,
        "Run2024H": 1,
        "Run2024I": 2,
    }    
}

campaigns = {
    "Run2022": {
        "Run3Summer22": ["Run2022C", "Run2022D"],
        "Run3Summer22EE": ["Run2022E", "Run2022F"]
    },
    "Run2023": {
        "Run3Summer23": ["Run2023C"],
        "Run3Summer23BPix": ["Run2023D"]
    },
    "Run2024": {
        "RunIII2024Summer24": ["Run2024C", "Run2024D", "Run2024E",
                               "Run2024F", "Run2024G", "Run2024H", "Run2024I"]
    }
}

# Data
for run in eras:
    for era, n_versions in eras[run].items():
        # Data
        for version in range(1, n_versions + 1):
            # ParkingDoubleMuonLowMass
            for i in range(8):
                samples["ParkingDoubleMuonLowMass"]["files"][era].append(
                    path1 + f"/ParkingDoubleMuonLowMass{i}+{era}-PromptReco-v{version}+MINIAOD/*root"
                )

            # Muon
            muon_pds = ["Muon0", "Muon1"]
            if re.search("Run2022", era):
                muon_pds = ["Muon"]
                if era == "Run2022C":
                    muon_pds = []
            for pd in muon_pds:
                samples["Muon"]["files"][era].append(
                    path1 + f"/{pd}+{era}-PromptReco-v{version}+MINIAOD/*root",
                )

            # EGamma
            egamma_pds = ["EGamma0", "EGamma1"]
            if re.search("Run2022", era):
                egamma_pds = ["EGamma"]
            for pd in egamma_pds:
                samples["EGamma"]["files"][era].append(
                    path2 + f"/{pd}+{era}-PromptReco-v{version}+MINIAOD/*root",
                )

            # ZeroBias
            samples["ZeroBias"]["files"][era].append(
                path2 + f"/ZeroBias+{era}-PromptReco-v{version}+MINIAOD/*root",
            )
            
            # HLTPhysics
            samples["HLTPhysics"]["files"][era].append(
                path2 + f"/HLTPhysics+{era}-PromptReco-v{version}+MINIAOD/*root",
            )

            # ParkingDoubleElectronLowMass
            if re.search("Run2022", era):
                for i in range(6):
                    samples["ParkingDoubleElectronLowMass"]["files"][era].append(
                        path3 + f"/ParkingDoubleElectronLowMass{i}+{era}-PromptReco-v{version}+MINIAOD/*root"
                    )
            elif re.search("Run2023", era):
                samples["ParkingDoubleElectronLowMass"]["files"][era].append(
                    path3 + f"/ParkingDoubleElectronLowMass+{era}-PromptReco-v{version}+MINIAOD/*root"
                )

                
samples["SingleMuon"]["files"]["Run2022C"].append(
    path + "/SingleMuon+Run2022C-PromptReco-v1+MINIAOD/*root",
)
samples["DoubleMuon"]["files"]["Run2022C"].append(
    path + "/DoubleMuon+Run2022C-PromptReco-v1+MINIAOD/*root",
)


# MC
samples["K0sToMuMu"] = {
        "Data": False,
        "files": defaultdict(list),
        "triggers": None
}

samples["K0sToMuMu"]["files"]["Run3Summer22"].append(
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*.root"
)
samples["K0sToMuMu"]["files"]["Run3Summer22EE"].append(
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v1+MINIAODSIM/*.root"
)
samples["K0sToMuMu"]["files"]["Run3Summer23"].append(
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4+MINIAODSIM/*.root"
)
samples["K0sToMuMu"]["files"]["Run3Summer23BPix"].append(
    path + "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4+MINIAODSIM/*.root"
)
samples["K0sToMuMu"]["files"]["RunIII2024Summer24"].append(
    path + "/K0sToMuMu_Fil-K0s_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v4+MINIAODSIM/*.root"
)


# pprint(samples)

# process samples
for sample_name in sample_names:
    if not recompute_results and sample_name in results:
        continue
    results[sample_name] = dict()
    for era in samples[sample_name]["files"]:
        results[sample_name][era] = process_sample(samples[sample_name], era)

# pprint(results)

print_results("Run2022")
print_results("Run2023")
print_results("Run2024")

## Save results
json.dump(results, open(output, 'w'), indent=4)
