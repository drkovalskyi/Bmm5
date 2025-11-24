"""Trigger efficiency study

The script measures the absolute trigger efficiency of the dimuon
trigger for Jpsi->mumu events as a function of the Jpsi pt and eta in
MC and Data . These results are used to calculate trigger efficiency
corrections for B to Jpsi X decays and their ratios.

The input data must be in Bmm5 NanoAOD format.

The script uses RDataFrame for optimal performance, but it needs to
read very large volumes of data to collect all information. To speed
up processing it first extracts histograms and stores them in a root
file, which will be used as a source unless instructed to recreate.

For more information about the methods check AN-24-186.

"""

#
# force_recreate options
# - [ ] - processing only missing samples
# - [ 'All' ] - force reprocessing for all samples
# - [ 'Sample1', 'Sample2' ] - force reprocessing of selected samples
#
# force_recreate = ['All']
# force_recreate = ["EGammaLowPt", "BuToJpsiK", "BsToJpsiPhi"]
# force_recreate = ["EGamma", "EGammaLowPt"]
force_recreate = []

debug = True
perform_method_validation = False
compute_run_average_efficiencies = True
make_plots = True
compute_jpsi_efficiency = False
compute_corrections = True

process_only = "RunIII2024Summer24|Run(2024|2025)"
# process_only = None

jpsi_pt_bins = [x / 2 for x in range(10, 21)]
jpsi_pt_bins.extend(range(11, 31))
jpsi_pt_bins.extend(range(32, 51,2))
b_pt_bins = [x / 2 for x in range(16, 41)]
b_pt_bins.extend([25, 30, 40, 50])

jpsi_low_pt_bins = [x / 2 for x in range(6, 21)]
jpsi_low_pt_bins.extend(range(11, 21))
b_low_pt_bins = [x / 2 for x in range(14, 25)]
b_low_pt_bins.extend([13, 14, 15, 16, 18, 20, 23])

mu_pt_bins = [3, 3.5, 4, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 20, 50]
eta_bins = [-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]
dphi_bins = [x / 5 for x in range(0,5)]
dphi_bins.extend([1.5, 3.2])

if debug:
    print(jpsi_pt_bins)
    print(b_pt_bins)
    print(jpsi_low_pt_bins)
    print(b_low_pt_bins)
    print(mu_pt_bins)
    print(eta_bins)
                    
# max_files = 10
max_files = 999999
# sample_names = ["BuToJpsiK", "BsToJpsiPhi",
#                 "ParkingDoubleMuonLowMass", "ParkingDoubleElectronLowMass",
#                 "Muon", "ZeroBias",
#                 "EGamma", "HLTPhysics", 'ZeroBiasExclusive', 'EGammaExclusive']
sample_names = [
    ## Monte Carlo
    "BsToJpsiPhi", 
    "BuToJpsiK",
    "BdToJpsiKstar",
    ## Data
    "EGammaLowPt",
    "EGamma",
    ## Cross checks for EGamma
    # "ZeroBias", 
    # "ParkingDoubleElectronLowMass",
    # "ZeroBiasExclusive", # HLT_ZeroBias events
    # "EGammaExclusive"    # HLT_Ele30_WPTight_Gsf events
]

path   = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/535/"
path1  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/532/trig/"
path2  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/535/mm/"
path3  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/529/mm"
path4   = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/534/"
# output_path = "/eos/home-d/dmytro/www/plots/2025/fsfu-trigger_efficiency/"
output_path = "/eos/home-d/dmytro/www/plots/tmp/2025/fsfu-trigger_efficiency/"

import sys
import ROOT
if not hasattr(ROOT.RDataFrame, "DefaultValueFor"):
    sys.exit(f"Error: ROOT.RDataFrame.DefaultValueFor is not available in ROOT "
             f"{ROOT.gROOT.GetVersion()}. Please use ROOT version 6.34 or later.")
from collections import defaultdict
from pprint import pprint
import os
import re
import glob
import json
import math
from array import array

ROOT.ROOT.EnableImplicitMT()

ROOT.gInterpreter.Declare('''
    using namespace ROOT::VecOps;
    RVec<double> dphi(const RVec<double> &phi1, const RVec<double> &phi2) {
       RVec<double> res;
       for (std::size_t i=0; i<phi1.size(); ++i) {
         res.emplace_back(acos(cos(phi2[i]-phi1[i])));
       }
       return res;
    }
    double dimuon_pt(double m1pt, double m1phi, double m2pt, double m2phi) {
       double pt2 = m1pt*m1pt + m2pt*m2pt + 2.0*m1pt*m2pt*cos(m2phi-m1phi);
       return pt2 > 0 ? sqrt(pt2) : 0;
    }

    double dimuon_eta(double pt1, double eta1, double phi1,
                      double pt2, double eta2, double phi2,
                      double m_mu = 0.105658)  // GeV
    {
       TLorentzVector v1, v2;
       v1.SetPtEtaPhiM(pt1, eta1, phi1, m_mu);
       v2.SetPtEtaPhiM(pt2, eta2, phi2, m_mu);
       return (v1 + v2).Eta();
    }
    RVec<double> b_jpsikstar_mass(const RVec<double> &b_kpi_mass, const RVec<double> &b_pik_mass,
                                  const RVec<double> &kpi_mass,   const RVec<double> &pik_mass) {
       const double kstar_pdg_mass = 0.89167;
       RVec<double> res;
       for (std::size_t i=0; i < kpi_mass.size(); ++i) {
           double b_mass = 0;
           if (fabs(kpi_mass[i] - kstar_pdg_mass) < fabs(pik_mass[i] - kstar_pdg_mass)) {
              if (kpi_mass[i] > 0.846 && kpi_mass[i] < 0.946) {
                 b_mass = b_kpi_mass[i];
              }
           } else {
              if (pik_mass[i] > 0.846 && pik_mass[i] < 0.946) {
                 b_mass = b_kpi_mass[i];
              }
           }
           res.emplace_back(b_mass);
       }
       return res;
    }
''')

class DataProcessor:
    def __init__(self):
        self.report = dict()

        # TH2 histograms
        self.histos = dict()

        self.hist_file = "rdf_jpsimm_trigger_efficiency.root"

        self.load_histograms()

        self.l1seeds = [
            # Run2022
            "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", "L1_DoubleMu4_SQ_OS_dR_Max1p2",
            # Run2023
            "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6",
            "L1_DoubleMu4_SQ_OS_dR_Max1p2",
            # Run2024
            "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
            "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"
        ]

    @staticmethod
    def histo_key(sample_name, era, histo_name):
        return f"{sample_name}___{era}___{histo_name}";
    
    @staticmethod
    def unpack_histo_key(key):
        parts = key.split("___")
        if len(parts) != 3:
            return None
        return tuple(parts)

    
    def load_histograms(self):
        """Load all histograms found in the file"""
        
        if not os.path.exists(self.hist_file):
            return 

        f = ROOT.TFile(self.hist_file)
        if f.IsZombie():
            raise RuntimeError(f"Could not open {self.hist_file} for reading")

        histos = {}
        try:
            keys = f.GetListOfKeys()
            if not keys:
                return histos

            for k in keys:
                hist = k.ReadObj()
                hist_types = ["TH1", "TH2", "TH3", "THnBase"]
                supported = False
                for hist_type in hist_types:
                    if hist.InheritsFrom(hist_type):
                        supported = True
                        break
                if supported:
                    name = hist.GetName()
                    f.Remove(hist)
                    # hist.SetDirectory(0)
                    info = self.unpack_histo_key(name)
                    if info:
                        (sample_name, era, hist_name) = info
                        if sample_name not in self.histos:
                            self.histos[sample_name] = dict()
                        if era not in self.histos[sample_name]:
                            self.histos[sample_name][era] = dict()
                    self.histos[sample_name][era][hist_name] = hist
        finally:
            f.Close()
            
    
    def define_samples(self):
    
        self.samples = {
            # "ParkingDoubleMuonLowMass": {
            #     "Data": True,
            #     "files": defaultdict(list),
            #     "triggers": ["HLT_Mu4_L1DoubleMu"]
            # },
            # "Muon": {
            #     "Data": True,
            #     "files": defaultdict(list),
            #     "triggers": ["HLT_Mu3_PFJet40", "HLT_Mu8"]
            # },
            # "SingleMuon": {
            #     "Data": True,
            #     "files": defaultdict(list),
            #     "triggers": ["HLT_Mu3_PFJet40"]
            # },
            # "DoubleMuon": {
            #     "Data": True,
            #     "files": defaultdict(list),
            #     "triggers": ["HLT_Mu8"]
            # },
            "ZeroBias": {
                "Data": True,
                "files": defaultdict(list),
                "triggers": None
            },
            "ZeroBiasExclusive": {
                "Data": True,
                "files": defaultdict(list),
                "triggers": None,
                "preselection": "HLT_ZeroBias"
            },
            "EGamma": {
                "Data": True,
                "files": defaultdict(list),
                "triggers": None,
            },
            "EGammaExclusive": {
                "Data": True,
                "files": defaultdict(list),
                "triggers": None,
                "preselection": "HLT_Ele30_WPTight_Gsf"
            },
            "EGammaLowPt": {
                "Data": True,
                "files": defaultdict(list),
                "triggers": None,
            },
            "ParkingDoubleElectronLowMass": {
                "Data": True,
                "files": defaultdict(list),
                "triggers": None
            },
            # "HLTPhysics": {
            #     "Data": True,
            #     "files": defaultdict(list),
            #     "triggers": None
            # },
        }

        # defaultdict(list)
        self.eras = {
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
            },
            "Run2025": {
                "Run2025C": 2,
                "Run2025D": 1,
                "Run2025E": 1,
                "Run2025F": 2,
                "Run2025G": 1,
            }    
        }

        self.campaigns = {
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
            },
            "Run2025": {
                "RunIII2024Summer24": ["Run2025C", "Run2025D", "Run2025E",
                                       "Run2025F", "Run2025G"]
            }
        }

        # Data
        for run in self.eras:
            for era, n_versions in self.eras[run].items():
                # Data
                for version in range(1, n_versions + 1):
                    # # ParkingDoubleMuonLowMass
                    # for i in range(8):
                    #     self.samples["ParkingDoubleMuonLowMass"]["files"][era].append(
                    #         path1 + f"/ParkingDoubleMuonLowMass{i}+{era}-PromptReco-v{version}+MINIAOD/*root"
                    #     )

                    # # Muon
                    # muon_pds = ["Muon0", "Muon1"]
                    # if re.search("Run2022", era):
                    #     muon_pds = ["Muon"]
                    #     if era == "Run2022C":
                    #         muon_pds = []
                    # for pd in muon_pds:
                    #     self.samples["Muon"]["files"][era].append(
                    #         path1 + f"/{pd}+{era}-PromptReco-v{version}+MINIAOD/*root",
                    #     )

                    # EGamma
                    egamma_pds = ["EGamma0", "EGamma1"]
                    if re.search("Run2022", era):
                        egamma_pds = ["EGamma"]
                    elif re.search("Run2022", era):
                        egamma_pds = ["EGamma0", "EGamma1", "EGamma2"]
                        
                    for pd in egamma_pds:
                        self.samples["EGamma"]["files"][era].append(
                            path2 + f"/{pd}+{era}-PromptReco-v{version}+MINIAOD/*root",
                        )
                        self.samples["EGammaExclusive"]["files"][era].append(
                            path + f"/{pd}+{era}-PromptReco-v{version}+MINIAOD/*root",
                        )
                        if re.search("Run20(24|25)", era):
                            self.samples["EGammaLowPt"]["files"][era].append(
                                path2 + f"/{pd}+{era}-PromptReco-v{version}+MINIAOD/*root",
                            )


                    # ZeroBias
                    self.samples["ZeroBias"]["files"][era].append(
                        path2 + f"/ZeroBias+{era}-PromptReco-v{version}+MINIAOD/*root",
                    )
                    self.samples["ZeroBiasExclusive"]["files"][era].append(
                        path + f"/ZeroBias+{era}-PromptReco-v{version}+MINIAOD/*root",
                    )

                    # # HLTPhysics
                    # self.samples["HLTPhysics"]["files"][era].append(
                    #     path2 + f"/HLTPhysics+{era}-PromptReco-v{version}+MINIAOD/*root",
                    # )

                    # ParkingDoubleElectronLowMass
                    if re.search("Run2022", era):
                        for i in range(6):
                            self.samples["ParkingDoubleElectronLowMass"]["files"][era].append(
                                path3 + f"/ParkingDoubleElectronLowMass{i}+{era}-PromptReco-v{version}+MINIAOD/*root"
                            )
                    elif re.search("Run2023", era):
                        self.samples["ParkingDoubleElectronLowMass"]["files"][era].append(
                            path3 + f"/ParkingDoubleElectronLowMass+{era}-PromptReco-v{version}+MINIAOD/*root"
                        )


        # self.samples["SingleMuon"]["files"]["Run2022C"].append(
        #     path + "/SingleMuon+Run2022C-PromptReco-v1+MINIAOD/*root",
        # )
        # self.samples["DoubleMuon"]["files"]["Run2022C"].append(
        #     path + "/DoubleMuon+Run2022C-PromptReco-v1+MINIAOD/*root",
        # )


        ###### MC

        ### BuToJpsiK
        
        self.samples["BuToJpsiK"] = {
                "Data": False,
                "files": defaultdict(list),
                "triggers": None
        }
        self.samples["BuToJpsiK"]["files"]["Run3Summer22"].append(
            # path + "/ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root"
            path + "/BuToJpsiK_JpsiToMuMu_MuFilter_Pt-2_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root"
        )
        self.samples["BuToJpsiK"]["files"]["Run3Summer22EE"].append(
            # path + "/ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root"
            # path + "/BuToJpsiK_JpsiToMuMu_MuFilter_Pt-2_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root"
            path + "/ButoJpsiK_Jpsito2Mu_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root"
        )
        self.samples["BuToJpsiK"]["files"]["Run3Summer23"].append(
            # path + "/ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/*.root"
            path + "/BuToJpsiK_JpsiToMuMu_MuFilter_Pt-2_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2+MINIAODSIM/*.root"
        )
        self.samples["BuToJpsiK"]["files"]["Run3Summer23BPix"].append(
            # path + "/ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/*.root"
            path + "/ButoJpsiK_Jpsito2Mu_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/*.root"
        )
        self.samples["BuToJpsiK"]["files"]["RunIII2024Summer24"].append(
            path + "/BuToJpsiK-JpsiToMuMu_Fil-MuPt2_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/*.root"
        )
        
        ### BsToJpsiPhi
        
        self.samples["BsToJpsiPhi"] = {
                "Data": False,
                "files": defaultdict(list),
                "triggers": None
        }
        self.samples["BsToJpsiPhi"]["files"]["Run3Summer22"] = [
            path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root",
            path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2+MINIAODSIM/*.root"
        ]
        self.samples["BsToJpsiPhi"]["files"]["Run3Summer22EE"] = [
            path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root",
            path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2+MINIAODSIM/*.root"
        ]
        self.samples["BsToJpsiPhi"]["files"]["Run3Summer23"].append(
            path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/*.root"
        )
        self.samples["BsToJpsiPhi"]["files"]["Run3Summer23BPix"].append(
            path + "/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/*.root"
        )
        self.samples["BsToJpsiPhi"]["files"]["RunIII2024Summer24"].append(
            path + "/BsToJPsiPhi-JpsiToMuMu-PhiToKK_Fil-MuPt2_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/*.root"
        )

        ### BdToJpsiKstar
        
        self.samples["BdToJpsiKstar"] = {
                "Data": False,
                "files": defaultdict(list),
                "triggers": None
        }
        self.samples["BdToJpsiKstar"]["files"]["RunIII2024Summer24"].append(
            path + "/BdToJpsiKstar-JpsiToMuMu_Fil-MuPt2_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/*.root"
        )
        
        
    @staticmethod
    def load_data(file_patterns, tree_name="Events"):
        assert isinstance(file_patterns, list), "expected a list"
        chain = ROOT.TChain(tree_name)
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

    def book_histo_ND(self, rdf_histos, rdf, sample_name, era, histo_name, bins, vars):
        full_histo_name = self.histo_key(sample_name, era, histo_name)
        if era not in rdf_histos:
            rdf_histos[era] = dict()
        nbins = [len(e) - 1 for e in bins]
        model = ROOT.RDF.THnDModel(full_histo_name, "", len(bins), nbins, bins)
        rdf_histos[era][histo_name] = rdf.HistoND(model, vars)
        
    def book_histo(self, rdf_histos, rdf, sample_name, era, histo_name, bins, vars):
        full_histo_name = self.histo_key(sample_name, era, histo_name)
        if era not in rdf_histos:
            rdf_histos[era] = dict()

        if len(bins) == 2:
            rdf_histos[era][histo_name] = \
                rdf.Histo2D((full_histo_name, "",
                             len(bins[0]) - 1, array('d', bins[0]), 
                             len(bins[1]) - 1, array('d', bins[1])), 
                            vars[0], vars[1])
        if len(bins) == 3:
            rdf_histos[era][histo_name] = \
                rdf.Histo3D((full_histo_name, "",
                             len(bins[0]) - 1, array('d', bins[0]), 
                             len(bins[1]) - 1, array('d', bins[1]), 
                             len(bins[2]) - 1, array('d', bins[2])), 
                            vars[0], vars[1], vars[1])
        elif len(bins) > 3: 
            nbins = [len(e) - 1 for e in bins]
            model = ROOT.RDF.THnDModel(full_histo_name, "", len(bins), nbins, bins)
            rdf_histos[era][histo_name] = rdf.HistoND(model, vars)
        
    def book_histo_2D(self, rdf_histos, rdf, sample_name, era, histo_name):
        self.book_histo(rdf_histos, rdf, sample_name, era, histo_name,
                        [jpsi_pt_bins, eta_bins],
                        ["jpsi_pt", "jpsi_eta"])

    def book_histo_3D(self, rdf_histos, rdf, sample_name, era, histo_name):
        self.book_histo(rdf_histos, rdf, sample_name, era, histo_name,
                        [b_pt_bins, jpsi_pt_bins, eta_bins],
                        ["pt", "jpsi_pt", "jpsi_eta"])

    def process_sample(self, sample_name, era):
        if process_only != None and not re.search(process_only, era):
            return

        sample = self.samples[sample_name]
        rdf_histos = dict()
        
        ### Prepare RDataFrame

        file_pattern = []
        for sample_era, pattern in sample["files"].items():
            if sample_era != era:
                continue
            file_pattern.extend(pattern)
        if len(file_pattern) == 0:
            return dict()

        run = None
        if sample["Data"]:
            for i_run, i_eras in self.eras.items():
                if era in i_eras:
                    run = i_run
                    break
        else:
            for i_run, i_campaigns in self.campaigns.items():
                if era in i_campaigns:
                    run = i_run
                    break

        pprint(file_pattern)

        chain = self.load_data(file_pattern)
        if chain == None:
            return dict()

        rdf = ROOT.RDataFrame(chain)
        # Make sure all triggers have a default value
        triggers = [
            # HLT
            "HLT_DoubleMu4_3_LowMass", "HLT_Mu4_L1DoubleMu", "HLT_Mu0_L1DoubleMu",
            "HLT_Mu3_PFJet40", "HLT_Mu8", "HLT_DoubleMu2_Jpsi_LowPt", "HLT_ZeroBias"
        ]
        triggers.extend(self.l1seeds)
        for trigger in triggers:
            rdf = rdf.DefaultValueFor(trigger, False)

        if "preselection" in sample:
            rdf = rdf.Filter(sample["preselection"])

        ### Extract information from other branches
        
        # mm
        rdf = rdf.Define("mm_mu1_mediumId", "Take(Muon_mediumId,  mm_mu1_index)")
        rdf = rdf.Define("mm_mu1_isGlobal", "Take(Muon_isGlobal,  mm_mu1_index)")
        rdf = rdf.Define("mm_mu1_npixels",  "Take(MuonId_nPixels, mm_mu1_index)")
        rdf = rdf.Define("mm_mu2_mediumId", "Take(Muon_mediumId,  mm_mu2_index)")
        rdf = rdf.Define("mm_mu2_isGlobal", "Take(Muon_isGlobal,  mm_mu2_index)")
        rdf = rdf.Define("mm_mu2_npixels",  "Take(MuonId_nPixels, mm_mu2_index)")
        rdf = rdf.Define("mm_dphi",         "dphi(mm_mu1_phi,     mm_mu2_phi)")

        # Offline selection
        rdf = rdf.Define("nom_cands_loose", "mm_mu1_index>=0 && mm_mu2_index>=0 && "
                         "mm_mu1_pt > 4 && mm_mu2_pt > 3 && "
                         "mm_mu1_mediumId && mm_mu2_mediumId && "
                         "mm_mu1_isGlobal && mm_mu2_isGlobal && "
                         "mm_mu1_pdgId * mm_mu2_pdgId == -169 && "
                         "abs(mm_kin_mass-3.1)<0.2 && mm_kin_vtx_prob>0.01")
        rdf = rdf.Define("nom_cands", "nom_cands_loose and abs(mm_kin_mass-3.09)<0.10")

        rdf = rdf.Define("lowpt_cands", "mm_mu1_index>=0 && mm_mu2_index>=0 && "
                         "(mm_mu1_pt < 4 || mm_mu2_pt < 3) && "
                         "mm_mu1_mediumId && mm_mu2_mediumId && "
                         "mm_mu1_isGlobal && mm_mu2_isGlobal && "
                         "mm_mu1_pdgId * mm_mu2_pdgId == -169 && "
                         "abs(mm_kin_mass-3.09)<0.10 && mm_kin_vtx_prob>0.01")

        rdf_nom_loose = rdf.Filter("Sum(nom_cands_loose)>0")
        rdf_nom_loose = rdf_nom_loose.Define("jpsi_mass",   "mm_kin_mass[nom_cands_loose]")
        rdf_nom = rdf.Filter("Sum(nom_cands)>0")
        rdf_nom = rdf_nom.Define("jpsi_pt",     "mm_kin_pt[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_eta",    "mm_kin_eta[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_m1_pt",  "mm_mu1_pt[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_m1_eta", "mm_mu1_eta[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_m2_pt",  "mm_mu2_pt[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_m2_eta", "mm_mu2_eta[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_dphi",   "mm_dphi[nom_cands]")
        rdf_nom = rdf_nom.Define("jpsi_mass",   "mm_kin_mass[nom_cands]")
        rdf_lowpt = rdf.Filter("Sum(lowpt_cands)>0")
        rdf_lowpt = rdf_lowpt.Define("jpsi_pt",  "mm_kin_pt[lowpt_cands]")
        rdf_lowpt = rdf_lowpt.Define("jpsi_eta", "mm_kin_eta[lowpt_cands]")

        low_pt_bins = [jpsi_low_pt_bins, eta_bins]
        low_pt_vars = ["jpsi_pt", "jpsi_eta"]
        b_bins = [b_pt_bins, eta_bins, jpsi_pt_bins, eta_bins]
        b_vars = ["b_pt", "b_eta", "jpsi_pt", "jpsi_eta"]
        b_lowpt_bins = [b_low_pt_bins, eta_bins, jpsi_low_pt_bins, eta_bins]
        b_lowpt_vars = ["b_pt", "b_eta", "jpsi_pt", "jpsi_eta"]
        bins = [jpsi_pt_bins, eta_bins, mu_pt_bins]
        vars = ["jpsi_pt", "jpsi_eta", "jpsi_m2_pt"]
        bins4D = [jpsi_pt_bins, eta_bins, mu_pt_bins, eta_bins]
        vars4D = ["jpsi_pt", "jpsi_eta", "jpsi_m2_pt", "jpsi_m2_eta"]
        bins4D2 = [mu_pt_bins, eta_bins, mu_pt_bins, eta_bins]
        vars4D2 = ["jpsi_m1_pt", "jpsi_m1_eta", "jpsi_m2_pt", "jpsi_m2_eta"]
        bins4D3 = [jpsi_pt_bins, eta_bins, mu_pt_bins, mu_pt_bins]
        vars4D3 = ["jpsi_pt", "jpsi_eta", "jpsi_m1_pt", "jpsi_m2_pt"]
        bins4D4 = [mu_pt_bins, mu_pt_bins, eta_bins, dphi_bins]
        vars4D4 = ["jpsi_m1_pt", "jpsi_m2_pt", "jpsi_eta", "jpsi_dphi"]
        bins5D = [mu_pt_bins, eta_bins, mu_pt_bins, eta_bins, dphi_bins]
        vars5D = ["jpsi_m1_pt", "jpsi_m1_eta", "jpsi_m2_pt", "jpsi_m2_eta", "jpsi_dphi"]

        self.book_histo_2D(rdf_histos, rdf_nom, sample_name, era, "All")
        self.book_histo_ND(rdf_histos, rdf_nom, sample_name, era, "All3D", bins, vars)
        self.book_histo_ND(rdf_histos, rdf_nom.Filter("HLT_DoubleMu4_3_LowMass"),
                           sample_name, era, "HLT_DoubleMu4_3_LowMass_3D", bins, vars)
        self.book_histo_ND(rdf_histos, rdf_nom, sample_name, era, "All4D", bins4D, vars4D)
        self.book_histo_ND(rdf_histos, rdf_nom.Filter("HLT_DoubleMu4_3_LowMass"),
                           sample_name, era, "HLT_DoubleMu4_3_LowMass_4D", bins4D, vars4D)
        self.book_histo_ND(rdf_histos, rdf_nom, sample_name, era, "All4D2", bins4D2, vars4D2)
        self.book_histo_ND(rdf_histos, rdf_nom.Filter("HLT_DoubleMu4_3_LowMass"),
                           sample_name, era, "HLT_DoubleMu4_3_LowMass_4D2", bins4D2, vars4D2)
        self.book_histo_ND(rdf_histos, rdf_nom, sample_name, era, "All4D3", bins4D3, vars4D3)
        self.book_histo_ND(rdf_histos, rdf_nom.Filter("HLT_DoubleMu4_3_LowMass"),
                           sample_name, era, "HLT_DoubleMu4_3_LowMass_4D3", bins4D3, vars4D3)
        self.book_histo_ND(rdf_histos, rdf_nom, sample_name, era, "All4D4", bins4D4, vars4D4)
        self.book_histo_ND(rdf_histos, rdf_nom.Filter("HLT_DoubleMu4_3_LowMass"),
                           sample_name, era, "HLT_DoubleMu4_3_LowMass_4D4", bins4D4, vars4D4)
        self.book_histo_ND(rdf_histos, rdf_nom, sample_name, era, "All5D", bins5D, vars5D)
        self.book_histo_ND(rdf_histos, rdf_nom.Filter("HLT_DoubleMu4_3_LowMass"),
                           sample_name, era, "HLT_DoubleMu4_3_LowMass_5D", bins5D, vars5D)
        full_histo_name = self.histo_key(sample_name, era, "jpsi_mass")
        rdf_histos[era]["jpsi_mass"] = rdf_nom_loose.Histo1D((full_histo_name,";mass, GeV", 100, 2.9, 3.3), "jpsi_mass")
        
        if sample_name in ["EGammaLowPt", "BsToJpsiPhi", "BuToJpsiK"]:
            self.book_histo(rdf_histos, rdf_lowpt, sample_name, era, "AllLowPt", low_pt_bins, low_pt_vars)
            self.book_histo(rdf_histos, rdf_lowpt.Filter("HLT_DoubleMu2_Jpsi_LowPt"),
                            sample_name, era, "HLT_DoubleMu2_Jpsi_LowPt_AllLowPt", low_pt_bins, low_pt_vars)

        
        if sample_name in ["BsToJpsiPhi", "BdToJpsiKstar"]:
            # bkkmm
            rdf = rdf.Define("bkkmm_mu1_mediumId",    "Take(mm_mu1_mediumId, bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu1_isGlobal",    "Take(mm_mu1_isGlobal, bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu1_pt",          "Take(mm_mu1_pt,       bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu1_eta",         "Take(mm_mu1_eta,      bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu2_mediumId",    "Take(mm_mu2_mediumId, bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu2_isGlobal",    "Take(mm_mu2_isGlobal, bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu2_pt",          "Take(mm_mu2_pt,       bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mu2_eta",         "Take(mm_mu2_eta,      bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mm_kin_pt",       "Take(mm_kin_pt,       bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mm_kin_eta",      "Take(mm_kin_eta,      bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mm_kin_vtx_prob", "Take(mm_kin_vtx_prob, bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_mm_dphi",         "Take(mm_dphi,         bkkmm_mm_index)")
            rdf = rdf.Define("bkkmm_jpsikstar_mass",
                             "b_jpsikstar_mass(bkkmm_jpsikpi_mass, bkkmm_jpsipik_mass, " + \
                             "bkkmm_jpsikpi_hh_mass, bkkmm_jpsipik_hh_mass)")

            selection = "bkkmm_mu1_pt > 4 && bkkmm_mu2_pt > 3" + \
                "&& bkkmm_mu1_mediumId && bkkmm_mu2_mediumId" + \
                "&& bkkmm_mu1_isGlobal && bkkmm_mu2_isGlobal && bkkmm_mm_kin_vtx_prob>0.01" + \
                "&& bkkmm_jpsikk_vtx_prob>0.025" + \
                "&& bkkmm_jpsikk_sl3d>3 && abs(bkkmm_jpsikk_alpha) < 0.1"

            if sample_name == "BsToJpsiPhi":
                selection += "&& abs(bkkmm_jpsikk_mass-5.37)<0.07 && abs(bkkmm_kk_mass-1.02)<0.01"
            elif sample_name == "BdToJpsiKstar":
                selection += "&& abs(bkkmm_jpsikstar_mass-5.3)<0.1"
            
            rdf = rdf.Define("nom_b_cands", selection)
            
            rdf_nom_bsd = rdf.Filter("Sum(nom_b_cands)>0")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_pt",     "bkkmm_mm_kin_pt[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_eta",    "bkkmm_mm_kin_eta[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_m1_pt",  "bkkmm_mu1_pt[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_m1_eta", "bkkmm_mu1_eta[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_m2_pt",  "bkkmm_mu2_pt[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_m2_eta", "bkkmm_mu2_eta[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("jpsi_dphi",   "bkkmm_mm_dphi[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("b_pt",        "bkkmm_jpsikk_pt[nom_b_cands]")
            rdf_nom_bsd = rdf_nom_bsd.Define("b_eta",       "bkkmm_jpsikk_eta[nom_b_cands]")

            # TODO: continue
            self.book_histo_2D(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All")
            self.book_histo_2D(rdf_histos, rdf_nom_bsd.Filter("HLT_DoubleMu4_3_LowMass"),
                               sample_name, era, "B_HLT_DoubleMu4_3_LowMass")

            # self.book_histo_3D(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All3D")
            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All3D",  bins,    vars)
            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All4D",  bins4D,  vars4D)
            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All4D2", bins4D2, vars4D2)
            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All4D3", bins4D3, vars4D3)
            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All4D4", bins4D4, vars4D4)
            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B_All5D",  bins5D,  vars5D)

            self.book_histo_ND(rdf_histos, rdf_nom_bsd, sample_name, era, "B",  b_bins,    b_vars)

            if era == "RunIII2024Summer24" and sample_name == "BsToJpsiPhi":
                rdf = rdf.Define("lowpt_bs_cands",
                                 "(bkkmm_mu1_pt < 4 || bkkmm_mu2_pt < 3) && "
                                 "bkkmm_mu1_mediumId && bkkmm_mu2_mediumId && "
                                 "bkkmm_mu1_isGlobal && bkkmm_mu2_isGlobal && "
                                 "bkkmm_mm_kin_vtx_prob>0.01 && bkkmm_jpsikk_vtx_prob>0.025 && "
                                 "bkkmm_jpsikk_sl3d>3 && abs(bkkmm_jpsikk_alpha) < 0.1 && "
                                 "abs(bkkmm_jpsikk_mass-5.37)<0.07 && abs(bkkmm_kk_mass-1.02)<0.01")

                rdf_lowpt_bs = rdf.Filter("Sum(lowpt_bs_cands)>0")
                rdf_lowpt_bs = rdf_lowpt_bs.Define("jpsi_pt",     "bkkmm_mm_kin_pt[lowpt_bs_cands]")
                rdf_lowpt_bs = rdf_lowpt_bs.Define("jpsi_eta",    "bkkmm_mm_kin_eta[lowpt_bs_cands]")
                rdf_lowpt_bs = rdf_lowpt_bs.Define("b_pt",        "bkkmm_jpsikk_pt[lowpt_bs_cands]")
                rdf_lowpt_bs = rdf_lowpt_bs.Define("b_eta",       "bkkmm_jpsikk_eta[lowpt_bs_cands]")

                self.book_histo(rdf_histos, rdf_lowpt_bs, sample_name, era, "B_AllLowPt",
                                low_pt_bins, low_pt_vars)
                self.book_histo(rdf_histos, rdf_lowpt_bs.Filter("HLT_DoubleMu2_Jpsi_LowPt"),
                                sample_name, era, "B_HLT_DoubleMu2_Jpsi_LowPt",
                                low_pt_bins, low_pt_vars)
                self.book_histo_ND(rdf_histos, rdf_lowpt_bs, sample_name, era, "BsLowPt",
                                   b_lowpt_bins, b_lowpt_vars)
            
        if sample_name == "BuToJpsiK":
            # bkmm
            rdf = rdf.Define("bkmm_mu1_mediumId",    "Take(mm_mu1_mediumId, bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu1_isGlobal",    "Take(mm_mu1_isGlobal, bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu1_pt",          "Take(mm_mu1_pt,       bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu1_eta",         "Take(mm_mu1_eta,      bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu2_mediumId",    "Take(mm_mu2_mediumId, bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu2_isGlobal",    "Take(mm_mu2_isGlobal, bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu2_pt",          "Take(mm_mu2_pt,       bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mu2_eta",         "Take(mm_mu2_eta,      bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mm_kin_pt",       "Take(mm_kin_pt,       bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mm_kin_eta",      "Take(mm_kin_eta,      bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mm_kin_vtx_prob", "Take(mm_kin_vtx_prob, bkmm_mm_index)")
            rdf = rdf.Define("bkmm_mm_dphi",         "Take(mm_dphi,         bkmm_mm_index)")

            rdf = rdf.Define("nom_bu_cands",
                             "bkmm_mu1_pt > 4 && bkmm_mu2_pt > 3 && "
                             "bkmm_mu1_mediumId && bkmm_mu2_mediumId && "
                             "bkmm_mu1_isGlobal && bkmm_mu2_isGlobal && "
                             "bkmm_mm_kin_vtx_prob>0.01 && bkmm_jpsimc_vtx_prob>0.025 && "
                             "bkmm_jpsimc_sl3d>3 && abs(bkmm_jpsimc_alpha) < 0.1 && "
                             "bkmm_kaon_pt>2.3 && abs(bkmm_jpsimc_mass-5.3)<0.1")
            
            rdf_nom_bu = rdf.Filter("Sum(nom_bu_cands)>0")
            rdf_nom_bu = rdf_nom_bu.Define("jpsi_pt",     "bkmm_mm_kin_pt[nom_bu_cands]")
            rdf_nom_bu = rdf_nom_bu.Define("jpsi_eta",    "bkmm_mm_kin_eta[nom_bu_cands]")
            rdf_nom_bu = rdf_nom_bu.Define("b_pt",        "bkmm_jpsimc_pt[nom_bu_cands]")
            rdf_nom_bu = rdf_nom_bu.Define("b_eta",       "bkmm_jpsimc_eta[nom_bu_cands]")
        
            self.book_histo_ND(rdf_histos, rdf_nom_bu, sample_name, era, "Bu",  b_bins,    b_vars)
            
            if era == "RunIII2024Summer24":
                rdf = rdf.Define("lowpt_bu_cands",
                                 "(bkmm_mu1_pt < 4 || bkmm_mu2_pt < 3) && "
                                 "bkmm_mu1_mediumId && bkmm_mu2_mediumId && "
                                 "bkmm_mu1_isGlobal && bkmm_mu2_isGlobal && "
                                 "bkmm_mm_kin_vtx_prob>0.01 && bkmm_jpsimc_vtx_prob>0.025 && "
                                 "bkmm_jpsimc_sl3d>3 && abs(bkmm_jpsimc_alpha) < 0.1 && "
                                 "bkmm_kaon_pt>2.3 && abs(bkmm_jpsimc_mass-5.3)<0.1")
            
                rdf_lowpt_bu = rdf.Filter("Sum(lowpt_bu_cands)>0")
                rdf_lowpt_bu = rdf_lowpt_bu.Define("jpsi_pt",     "bkmm_mm_kin_pt[lowpt_bu_cands]")
                rdf_lowpt_bu = rdf_lowpt_bu.Define("jpsi_eta",    "bkmm_mm_kin_eta[lowpt_bu_cands]")
                rdf_lowpt_bu = rdf_lowpt_bu.Define("b_pt",        "bkmm_jpsimc_pt[lowpt_bu_cands]")
                rdf_lowpt_bu = rdf_lowpt_bu.Define("b_eta",       "bkmm_jpsimc_eta[lowpt_bu_cands]")
                
                self.book_histo_ND(rdf_histos, rdf_lowpt_bu, sample_name, era, "BuLowPt",
                                   b_lowpt_bins, b_lowpt_vars)

            
        if run == "Run2022":
            l1_selection_unprescaled = "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6||L1_DoubleMu4_SQ_OS_dR_Max1p2"
            l1_selection_real_mix = l1_selection_unprescaled
        elif run == "Run2023":
            l1_selection_unprescaled = "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4||L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6"+ \
                "||L1_DoubleMu4_SQ_OS_dR_Max1p2"
            l1_selection_real_mix = l1_selection_unprescaled
        elif run in ["Run2024", "Run2025"]:
            l1_selection_unprescaled = "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4||L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6"+ \
                "||L1_DoubleMu4p5_SQ_OS_dR_Max1p2"
            l1_selection_real_mix = "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"

        # ref_trigger = "HLT_Mu4_L1DoubleMu"
        # if sample["triggers"] == None or ref_trigger in sample["triggers"]:
        #     name = f"{ref_trigger}"
        #     rdf_nom1 = rdf_nom.Filter(ref_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom1, sample_name, era, name)

        #     test_trigger = "HLT_DoubleMu4_3_LowMass"
        #     name += f"_{test_trigger}"
        #     rdf_nom1 = rdf_nom1.Filter(test_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom1, sample_name, era, name)

        # ref_trigger = "HLT_Mu3_PFJet40"
        # if sample["triggers"] == None or ref_trigger in sample["triggers"]:
        #     name = f"{ref_trigger}"
        #     rdf_nom2 = rdf_nom.Filter(ref_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom2, sample_name, era, name)

        #     test_trigger = "HLT_Mu0_L1DoubleMu"
        #     name += f"_{test_trigger}"
        #     rdf_nom2 = rdf_nom2.Filter(test_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom2, sample_name, era, name)

        #     test_trigger = l1_selection_unprescaled
        #     name += f"_L1_DoubleMu_Mix_Unprescaled"
        #     rdf_nom2 = rdf_nom2.Filter(test_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom2, sample_name, era, name)

        # ref_trigger = "HLT_Mu8"
        # if sample["triggers"] == None or ref_trigger in sample["triggers"]:
        #     name = f"{ref_trigger}"
        #     rdf_nom3 = rdf_nom.Filter(ref_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom3, sample_name, era, name)

        #     test_trigger = "HLT_Mu0_L1DoubleMu"
        #     name += f"_{test_trigger}"
        #     rdf_nom3 = rdf_nom3.Filter(test_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom3, sample_name, era, name)

        #     test_trigger = l1_selection_unprescaled
        #     name += f"_L1_DoubleMu_Mix_Unprescaled"
        #     rdf_nom3 = rdf_nom3.Filter(test_trigger)
        #     self.book_histo_2D(rdf_histos, rdf_nom3, sample_name, era, name)

        # no reference
        if sample["triggers"] == None:
            test_triggers = ["HLT_Mu0_L1DoubleMu", "HLT_DoubleMu4_3_LowMass"]
            for test_trigger in test_triggers:
                self.book_histo_2D(rdf_histos, rdf_nom.Filter(test_trigger), sample_name, era, test_trigger)

            test_trigger = "HLT_DoubleMu2_Jpsi_LowPt"
            self.book_histo_2D(rdf_histos, rdf_lowpt.Filter(test_trigger), sample_name, era, f"lowpt_{test_trigger}")

            # L1
            rdf_nom_l1_unprescaled = rdf_nom.Filter(l1_selection_unprescaled)
            name = "HLT_DoubleMu4_3_LowMass_L1_DoubleMu_Mix_Unprescaled"
            self.book_histo_2D(rdf_histos, rdf_nom_l1_unprescaled.Filter("HLT_DoubleMu4_3_LowMass"),
                               sample_name, era, name)

            rdf_nom_l1_real_mix = rdf_nom.Filter(l1_selection_real_mix)
            name = "HLT_DoubleMu4_3_LowMass_L1_DoubleMu_Real_Mix"
            self.book_histo_2D(rdf_histos, rdf_nom_l1_real_mix.Filter("HLT_DoubleMu4_3_LowMass"),
                               sample_name, era, name)

            for l1seed in self.l1seeds:
                rdf_nom_l1 = rdf_nom.Filter(l1seed)
                name = f"HLT_DoubleMu4_3_LowMass_{l1seed}"
                self.book_histo_2D(rdf_histos, rdf_nom_l1.Filter("HLT_DoubleMu4_3_LowMass"),
                                   sample_name, era, name)

        if sample_name not in self.histos:
            self.histos[sample_name] = dict()
        for era, hist_dict in rdf_histos.items():
            if era not in self.histos[sample_name]:
                self.histos[sample_name][era] = dict()
            for hist_name, rdf_hist in hist_dict.items():
                # Materialize histogram
                hist =  rdf_hist.GetValue()
                # hist.SetDirectory(0)
                full_hist_name = self.histo_key(sample_name, era, hist_name)                    
                self.histos[sample_name][era][hist_name] = hist

                # save the histogram
                mode = "UPDATE" if os.path.exists(self.hist_file) else "RECREATE"
                of = ROOT.TFile(self.hist_file, mode)
                hist.Write(hist.GetName(), ROOT.TObject.kOverwrite)
                of.Close()


    def make_report(self):
        for sample_name, era_hist_dict in self.histos.items():
            for era, hist_dict in era_hist_dict.items():
                for hist_name, hist in hist_dict.items():
                    full_hist_name = self.histo_key(sample_name, era, hist_name)
                    if hist.InheritsFrom("THnBase"):
                        self.report[full_hist_name] = hist.Integral(True)
                    else:
                        self.report[full_hist_name] = hist.Integral()


    def process_samples(self, sample_names, force_recreate=[]):
        
        self.define_samples()
        recreate_all = False
        if "All" in force_recreate:
            recreate_all = True
        
        if recreate_all and os.path.exists(self.hist_file):
            os.remove(self.hist_file)
        
        # process samples
        for sample_name in sample_names:
            if not recreate_all and sample_name not in force_recreate and sample_name in self.histos:
                continue
            for era in self.samples[sample_name]["files"]:
                # if not re.search('^(Run3Summer23|Run2023)', era):
                #     continue
                print(f"Processing {sample_name} - {era}")
                self.process_sample(sample_name, era)

        self.make_report()


def weighted_efficiency_thn_allbins(h_num, h_den, h_targ):
    # Basic type and size checks
    for h in (h_num, h_den, h_targ):
        assert h.InheritsFrom("THnBase")
    nlin = h_den.GetNbins()
    assert nlin == h_num.GetNbins() == h_targ.GetNbins()

    # First pass: sum target over bins with valid denominator
    Tsum = 0.0
    # skip 0 bin since python calls a wrong method (0 is interpreted a s pointer I guess)
    for ib in range(1,nlin):
        d = h_den.GetBinContent(ib)
        if d > 0:
            Tsum += h_targ.GetBinContent(ib)
    if Tsum <= 0:
        return None, None

    # Second pass: E = Σ p e, Var(E) = Σ p² e(1−e)/d
    E = 0.0
    varE = 0.0
    for ib in range(1,nlin):
        d = h_den.GetBinContent(ib)
        if d <= 0:
            continue
        n = h_num.GetBinContent(ib)
        e = n / d
        p = h_targ.GetBinContent(ib) / Tsum
        E += p * e
        varE += (p * p) * e * (1.0 - e) / d

    return E, math.sqrt(varE)


def weighted_efficiency(h_num, h_den, h_target):

    # Histograms are expected to contain non-weighted event counts

    if h_num.InheritsFrom("THnBase"):
        return weighted_efficiency_thn_allbins(h_num, h_den, h_target)
    
    # Basic checks
    assert isinstance(h_num, ROOT.TH2) and isinstance(h_den, ROOT.TH2) and isinstance(h_target, ROOT.TH2)
    assert h_num.GetNbinsX()==h_den.GetNbinsX()==h_target.GetNbinsX()
    assert h_num.GetNbinsY()==h_den.GetNbinsY()==h_target.GetNbinsY()
    nx, ny = h_den.GetNbinsX(), h_den.GetNbinsY()

    # Use only bins with den>0 and target>0, then renormalize target over those bins
    valid = []
    Tsum = 0.0
    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            d = h_den.GetBinContent(ix, iy)
            t = h_target.GetBinContent(ix, iy)
            if d > 0 and t > 0:
                valid.append((ix, iy))
                Tsum += t
    if not valid or Tsum <= 0:
        return None, None

    # Compute E = sum p_i e_i with Var(E) = sum p_i^2 * e_i(1-e_i)/d_i
    E = 0.0
    var_ctrl = 0.0
    p = {}
    for ix, iy in valid:
        d = h_den.GetBinContent(ix, iy)
        n = h_num.GetBinContent(ix, iy)
        e = n / d
        pi = h_target.GetBinContent(ix, iy) / Tsum
        p[(ix, iy)] = pi
        E += pi * e
        var_ctrl += (pi * pi) * e * (1.0 - e) / d

    return E, math.sqrt(var_ctrl)

def add_histos(hists, out_name):
    if not hists:
        return None

    hsum = hists[0].Clone(out_name)

    if hsum.InheritsFrom("TH1") or hsum.InheritsFrom("TH2") or hsum.InheritsFrom("TH3"):
        hsum.SetDirectory(0)

    for h in hists[1:]:
        hsum.Add(h)

    return hsum

def weighted_efficiency_1D(h_num, h_den, h_target, axis="x", out_name=None):
    """
    For each bin along `axis`, average over the other axis with
    p = T(ix,iy) / sum_y T(ix,iy)  [or swapped if axis='y'].
    Errors use Var(E) = sum p^2 * e*(1-e)/d with e=n/d.
    """

    assert isinstance(h_num, ROOT.TH2) and isinstance(h_den, ROOT.TH2) and isinstance(h_target, ROOT.TH2)
    assert h_num.GetNbinsX()==h_den.GetNbinsX()==h_target.GetNbinsX()
    assert h_num.GetNbinsY()==h_den.GetNbinsY()==h_target.GetNbinsY()

    nx, ny = h_den.GetNbinsX(), h_den.GetNbinsY()
    ax = axis.lower()
    if ax not in ("x", "y"):
        raise ValueError("axis must be 'x' or 'y'")

    if ax == "x":
        h1 = h_den.ProjectionX(out_name or "eff_vs_x", 1, ny, "e")
        h1.SetDirectory(0)
        h1.Reset("ICES")
        h1.SetMinimum(0)
        h1.SetMaximum(1.)

        for ix in range(1, nx+1):
            # find valid y in this x-slice
            valid_y = []
            Tsum = 0.0
            for iy in range(1, ny+1):
                d = h_den.GetBinContent(ix, iy)
                t = h_target.GetBinContent(ix, iy)
                if d > 0 and t > 0:
                    valid_y.append(iy)
                    Tsum += t
            if not valid_y or Tsum <= 0:
                continue

            E = 0.0
            var_ctrl = 0.0
            for iy in valid_y:
                d = h_den.GetBinContent(ix, iy)
                n = h_num.GetBinContent(ix, iy)
                e = n / d
                p = h_target.GetBinContent(ix, iy) / Tsum
                E += p * e
                var_ctrl += (p * p) * e * (1.0 - e) / d

            h1.SetBinContent(ix, E)
            h1.SetBinError(ix, math.sqrt(var_ctrl))
        return h1

    # axis == "y"
    h1 = h_den.ProjectionY(out_name or "eff_vs_y", 1, nx, "e")
    h1.SetDirectory(0)
    h1.Reset("ICES")
    h1.SetMinimum(0)
    h1.SetMaximum(1.)
        
    for iy in range(1, ny+1):
        valid_x = []
        Tsum = 0.0
        for ix in range(1, nx+1):
            d = h_den.GetBinContent(ix, iy)
            t = h_target.GetBinContent(ix, iy)
            if d > 0 and t > 0:
                valid_x.append(ix)
                Tsum += t
        if not valid_x or Tsum <= 0:
            continue

        E = 0.0
        var_ctrl = 0.0
        for ix in valid_x:
            d = h_den.GetBinContent(ix, iy)
            n = h_num.GetBinContent(ix, iy)
            e = n / d
            p = h_target.GetBinContent(ix, iy) / Tsum
            E += p * e
            var_ctrl += (p * p) * e * (1.0 - e) / d

        h1.SetBinContent(iy, E)
        h1.SetBinError(iy, math.sqrt(var_ctrl))
    return h1

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print(f"{path}/{output_name_without_extention}.png")
    canvas.Print(f"{path}/{output_name_without_extention}.pdf")
    canvas.Print(f"{path}/{output_name_without_extention}.root")

def get_integral(hist):
    # Fixme: would be nice to be consisten
    if hist.InheritsFrom("THnBase"):
        # integrall with underflow and overflow entries
        return hist.Integral(True)
    else:
        # integrall without underflow and overflow entries
        return hist.Integral()
    
def compute_efficiency(histos, pd, era, h_name_num, h_name_denom, prefix):
    h_num = processor.histos[pd][era][h_name_num]
    h_den = processor.histos[pd][era][h_name_denom]
    n_num = get_integral(h_num)
    n_den = get_integral(h_den)
    eff   = n_num / n_den
    eff_err = math.sqrt(eff * (1 - eff)/n_den)
    print(f"{prefix} efficiency using {pd} {era}: {eff:0.4f} +/- {eff_err:0.4f}")

def aggregate_histograms(histos, pd, era_pattern, h_name, h_agg_name="h_agg"):
    if pd in histos:
        hs = []
        for h_era in histos[pd]:
            if not re.search(f'^{era_pattern}', h_era):
                continue
            for name, hist in histos[pd][h_era].items():
                if name == h_name:
                    hs.append(hist)

        h = add_histos(hs, h_name)
        return h
    else:
        if debug:
            print(f"{pd} not found in histograms")
        return None

def aggregate_two_sets_of_histograms(histos, pd, era_pattern, h_name_num, h_name_denom):
    h_num = aggregate_histograms(histos, pd, era_pattern, h_name_num, "h_num")
    h_den = aggregate_histograms(histos, pd, era_pattern, h_name_denom, "h_den")
    return (h_num, h_den)
        
def compute_reweighted_efficiency(histos, pd, era_pattern, h_name_num, h_name_denom, h_target, plot=False):

    h_num, h_den = aggregate_two_sets_of_histograms(histos, pd, era_pattern, h_name_num, h_name_denom)
    if h_num:
        eff, eff_err = weighted_efficiency(h_num, h_den, h_target)
        print(f"{pd} {era_pattern} efficiency: {eff:.4f} +/- {eff_err:.4f}")
        if plot:
            h_eff_pt = weighted_efficiency_1D(h_num, h_den, h_target, axis="x")
            h_eff_pt.Draw()
            print_canvas(f"{h_name_num}_wrt_{h_name_denom}-{pd}_{era_pattern}_eff_vs_pt", output_path)
            h_eff_eta = weighted_efficiency_1D(h_num, h_den, h_target, axis="y")
            h_eff_eta.Draw()
            print_canvas(f"{h_name_num}_wrt_{h_name_denom}-{pd}_{era_pattern}_eff_vs_eta", output_path)
        return eff, eff_err
    else:
        return None, None

def _collect_bins_and_norms(h_den, h_tA, h_tB):
    nx, ny = h_den.GetNbinsX(), h_den.GetNbinsY()
    bins = []
    TsumA = 0.0
    TsumB = 0.0
    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            d = h_den.GetBinContent(ix, iy)
            if d <= 0:
                continue
            tA = h_tA.GetBinContent(ix, iy)
            tB = h_tB.GetBinContent(ix, iy)
            bins.append((ix, iy, d, tA, tB))
            if tA > 0:
                TsumA += tA
            if tB > 0:
                TsumB += tB
    return bins, TsumA, TsumB

def weighted_efficiency_ratio(h_num, h_den, h_tA, h_tB):
    # checks
    for h in (h_num, h_den, h_tA, h_tB):
        assert isinstance(h, ROOT.TH2)
    assert h_num.GetNbinsX()==h_den.GetNbinsX()==h_tA.GetNbinsX()==h_tB.GetNbinsX()
    assert h_num.GetNbinsY()==h_den.GetNbinsY()==h_tA.GetNbinsY()==h_tB.GetNbinsY()

    bins, TsumA, TsumB = _collect_bins_and_norms(h_den, h_tA, h_tB)
    if not bins or TsumA <= 0 or TsumB <= 0:
        return None, None

    EA = EB = 0.0
    VA = VB = 0.0
    CAB = 0.0

    for ix, iy, d, tA, tB in bins:
        n = h_num.GetBinContent(ix, iy)
        e = n / d
        ve = e * (1.0 - e) / d  # binomial variance

        pA = (tA / TsumA) if tA > 0 else 0.0
        pB = (tB / TsumB) if tB > 0 else 0.0

        EA += pA * e
        EB += pB * e
        VA += (pA*pA) * ve
        VB += (pB*pB) * ve
        CAB += (pA*pB) * ve

    if EB <= 0:
        return None, None

    # Delta-method variance of the ratio
    varR = VA/(EB*EB) + (EA*EA)*VB/(EB**4) - 2.0*EA*CAB/(EB**3)
    dR = math.sqrt(max(0.0, varR))

    R = EA / EB
    return R, dR

    
##########################################################

# def aggregate_results(run):
#     report = dict()
#     for sample_name, info in samples.items():
#         if sample_name not in results:
#             continue
#         report[sample_name] = defaultdict(int)
#         if info["Data"]:
#             for era in results[sample_name]:
#                 if era in eras[run]:
#                     for selection, n in results[sample_name][era].items():
#                         report[sample_name][selection] += n
#         else:
#             for campaign in results[sample_name]:
#                 if campaign in campaigns[run]:
#                     for selection, n in results[sample_name][campaign].items():
#                         report[sample_name][selection] += n
#     return report

def get_efficiency(n, n0):
    if n0 == 0:
        return ""
    eff = n / float(n0)
    eff_err = math.sqrt(eff * (1 - eff) / n0)

    return  f"{eff * 100:0.1f} \pm {eff_err * 100:0.1f}"
    
        
# def print_results(run):
#     report = aggregate_results(run)

#     eff_report = dict()
#     eff_report = {
#         "Absolute efficiency": {
#             "Data": {},
#             "MC": {
#                 "BuToJpsiK": get_efficiency(report["BuToJpsiK"]["HLT_DoubleMu4_3_LowMass"],
#                                             report["BuToJpsiK"]["All"]),
#             }
#         },
#         "Absolute HLT efficiency": {
#             "Data": {},
#             "MC": {
#                 "BuToJpsiK": get_efficiency(report["BuToJpsiK"]["L1_DoubleMu_Mix_Unprescaled&&HLT_DoubleMu4_3_LowMass"],
#                                             report["BuToJpsiK"]["L1_DoubleMu_Mix_Unprescaled"]),
#             }
#         },
#         "Relative HLT efficiency": {
#             "Data": {},
#             "MC": {
#                 "BuToJpsiK": get_efficiency(report["BuToJpsiK"]["HLT_Mu4_L1DoubleMu&&HLT_DoubleMu4_3_LowMass"],
#                                             report["BuToJpsiK"]["HLT_Mu4_L1DoubleMu"]),
#             }
#         },
#         "Absolute L1 efficiency": {
#             "Data": {},
#             "MC": {
#                 "BuToJpsiK": get_efficiency(report["BuToJpsiK"]["L1_DoubleMu_Mix_Unprescaled"],
#                                             report["BuToJpsiK"]["All"]),
#             }
#         },
#         "Relative L1 efficiency": {
#             "Data": {},
#             "MC": {
#                 "BuToJpsiK": get_efficiency(report["BuToJpsiK"]["HLT_Mu8&&HLT_Mu0_L1DoubleMu"],
#                                             report["BuToJpsiK"]["HLT_Mu8"]),
#             }
#         },
        
#     }
#     for pd in ["HLTPhysics", "EGamma"]:
#         eff_report["Absolute HLT efficiency"]["Data"][pd] = \
#             get_efficiency(report[pd]["L1_DoubleMu_Mix_Unprescaled&&HLT_DoubleMu4_3_LowMass"],
#                            report[pd]["L1_DoubleMu_Mix_Unprescaled"])
#     for pd in ["ParkingDoubleMuonLowMass"]:
#         eff_report["Relative HLT efficiency"]["Data"][pd] = \
#             get_efficiency(report[pd]["HLT_Mu4_L1DoubleMu&&HLT_DoubleMu4_3_LowMass"],
#                            report[pd]["HLT_Mu4_L1DoubleMu"])
#     for pd in ["EGamma", "ZeroBias", "ParkingDoubleElectronLowMass"]:
#         eff_report["Absolute L1 efficiency"]["Data"][pd] = \
#             get_efficiency(report[pd]["L1_DoubleMu_Mix_Unprescaled"],
#                            report[pd]["All"])
#     for pd in ["Muon"]:
#         eff_report["Relative L1 efficiency"]["Data"][pd] = \
#             get_efficiency(report[pd]["HLT_Mu8&&HLT_Mu0_L1DoubleMu"],
#                            report[pd]["HLT_Mu8"])
#     for pd in ["EGamma"]:
#         eff_report["Absolute efficiency"]["Data"][pd] = \
#             get_efficiency(report[pd]["HLT_DoubleMu4_3_LowMass"],
#                            report[pd]["All"])
            
    
    
#     pprint(report)
#     pprint(eff_report)

def thnd_to_th2(hn, x_axis, y_axis, name="h2",
                include_overflow=False, bins=None, ranges=None):
    """
    bins : dict[int, tuple[int,int]] or dict[int, int], optional
        Per-axis bin selection. Values are either a single 1-based bin (fix) or
        a (first,last) inclusive bin range. You may use 0 and nbins+1 to include UF/OF.
    ranges : dict[int, tuple[float,float]], optional
        Per-axis coordinate range (xmin,xmax), applied with SetRangeUser.
        Ignored for axes present in `bins`.
    """
    ndim = hn.GetNdimensions()
    bins  = bins  or {}
    ranges = ranges or {}

    # Save and then set per-axis ranges
    axes = [hn.GetAxis(i) for i in range(ndim)]
    saved = [(ax.GetFirst(), ax.GetLast()) for ax in axes]

    for i, ax in enumerate(axes):
        nb = ax.GetNbins()

        if i in bins:
            sel = bins[i]
            if isinstance(sel, int):
                first = last = int(sel)
            else:
                first, last = map(int, sel)
            ax.SetRange(first, last)
            continue

        if i in ranges:
            xmin, xmax = map(float, ranges[i])
            ax.SetRangeUser(xmin, xmax)
            continue

        # No explicit constraint: full range, optionally with UF/OF
        if include_overflow:
            ax.SetRange(0, nb + 1)   # includes UF and OF
        else:
            ax.SetRange(1, nb)       # in-range only

    # Do the projection. PyROOT accepts a Python list of axis indices.
    h2 = hn.Projection(y_axis, x_axis)
    h2.SetName(name)
    h2.SetDirectory(0)  # detach; Projection attaches to gDirectory

    # Restore original axis ranges
    for ax, (f, l) in zip(axes, saved):
        ax.SetRange(f, l)

    return h2




####################################################################

# if os.path.exists(output):
#     results = json.load(open(output))
	




# # pprint(results)

# # print_results("Run2022")
# print_results("Run2023")
# # print_results("Run2024")

if __name__ == "__main__":
    
    #################################
    ## SETUP
    #################################
    
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetLabelSize(0.045, "X")
    ROOT.gStyle.SetLabelSize(0.045, "Y")
    ROOT.gStyle.SetTitleSize(0.045, "X")
    ROOT.gStyle.SetTitleSize(0.045, "Y")
    ROOT.gStyle.SetTitleOffset(1.2, "X")
    ROOT.gStyle.SetTitleOffset(1.2, "Y")
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetOptFit(1)

    c = ROOT.TCanvas("c", "", 600, 600)

    # histos = defaultdict(dict)
    
    #######################################################
    ## Extract information for Jpsi efficiency measurement
    ##
    ## Note: it's slow process - don't force recreation of
    ## histograms unless you change data or selection
    #######################################################
    
    processor = DataProcessor()
    # pprint(processor.samples)
    processor.process_samples(sample_names, force_recreate)
    # pprint(processor.report)

    #######################################################
    ## Jpsi efficiency summary
    #######################################################
    print("\nAbsolute HLT_DoubleMu4_3_LowMass efficiency")

    ##############################################################
    ## Method validation
    ##############################################################

    if perform_method_validation:
        print("\nMethod validation")
        # campaign = 'RunIII2024Summer24'
        campaign = 'Run3Summer22'
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All']
        h_target3D = processor.histos['BsToJpsiPhi'][campaign]['B_All3D']

        compute_efficiency(processor.histos, "BsToJpsiPhi", campaign,
                           "B_HLT_DoubleMu4_3_LowMass", "B_All", "\tBsToJpsiPhi")
        compute_efficiency(processor.histos, "BuToJpsiK", campaign,
                           "HLT_DoubleMu4_3_LowMass", "All", "\tJpsi")
        compute_efficiency(processor.histos, "BsToJpsiPhi", campaign,
                           "HLT_DoubleMu4_3_LowMass", "All", "\tJpsi")
        
        print("Reweighted efficiency using BsToJpsiPhi Run3Summer22 as target:")
        print("2D")
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", "Run3Summer22", "HLT_DoubleMu4_3_LowMass", "All", h_target)
        compute_reweighted_efficiency(processor.histos, "BsToJpsiPhi", "Run3Summer22", "HLT_DoubleMu4_3_LowMass", "All", h_target)
                
        print("3D")
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All3D']
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", campaign, "HLT_DoubleMu4_3_LowMass_3D", "All3D", h_target)
        print("4D")
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All4D']
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", campaign, "HLT_DoubleMu4_3_LowMass_4D", "All4D", h_target)
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All4D2']
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", campaign, "HLT_DoubleMu4_3_LowMass_4D2", "All4D2", h_target)
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All4D3']
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", campaign, "HLT_DoubleMu4_3_LowMass_4D3", "All4D3", h_target)
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All4D4']
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", campaign, "HLT_DoubleMu4_3_LowMass_4D4", "All4D4", h_target)
        h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All5D']
        compute_reweighted_efficiency(processor.histos, "BuToJpsiK", campaign, "HLT_DoubleMu4_3_LowMass_5D", "All5D", h_target)

        print("Data validation")
        h_target = processor.histos['BsToJpsiPhi']['Run3Summer22']['B_All']
        # pds = ['ZeroBias', 'ZeroBiasExclusive', 'EGamma', 'EGammaExclusive', 'ParkingDoubleElectronLowMass']
        pds = ['ZeroBias', 'EGamma', 'ParkingDoubleElectronLowMass']
        for pd in pds:
            for run in ["Run2022", "Run2023", "Run2024", "Run2025"]:
                if run != "Run2022" and pd == "ParkingDoubleElectronLowMass":
                    continue
                # for era in processor.eras[run]:
                #     compute_efficiency(processor.histos, pd, era,
                #                        "HLT_DoubleMu4_3_LowMass", "All", "\tJpsi")
                #     print("\t2D ", end="")
                #     compute_reweighted_efficiency(processor.histos, pd, era,
                #                                   "HLT_DoubleMu4_3_LowMass", "All", h_target)
                #     print("\t3D ", end="")
                #     compute_reweighted_efficiency(processor.histos, pd, era,
                #                                   "HLT_DoubleMu4_3_LowMass_3D", "All3D", h_target3D)
                print("\t", end="")
                compute_reweighted_efficiency(processor.histos, pd, run,
                                              "HLT_DoubleMu4_3_LowMass", "All", h_target)
                # h_num_data, h_den_data = aggregate_two_sets_of_histograms(processor.histos, pd, run,
                #                                               "HLT_DoubleMu4_3_LowMass", "All")
                # h_den_data.GetXaxis().SetTitle("p_{T}, GeV")
                # h_den_data.GetYaxis().SetTitle("\eta")
                # h_den_data.Draw("colz")
                # print_canvas(f"jpsi_all_-{pd}_{run}", output_path)
                # h_jpsi_mass = aggregate_histograms(processor.histos, pd, run, "jpsi_mass", "jpsi_mass")
                # if h_jpsi_mass:
                #     h_jpsi_mass.SetLineWidth(2)
                #     h_jpsi_mass.Draw()
                #     print_canvas(f"jpsi_mass_-{pd}_{run}", output_path)


    if compute_run_average_efficiencies:
        print("\nAverage efficiency of HLT_DoubleMu4_3_LowMass for nominal BsToJpsiPhi selection")
        for campaign in ['Run3Summer22', 'Run3Summer23', 'RunIII2024Summer24']:
            if campaign not in processor.histos['BsToJpsiPhi']:
                continue
            compute_efficiency(processor.histos, "BsToJpsiPhi", campaign,
                               "B_HLT_DoubleMu4_3_LowMass", "B_All", "BsToJpsiPhi")
        
        pd = 'EGamma'
        for run,campaign in [("Run2022", 'Run3Summer22'),
                             ("Run2023", 'Run3Summer23'),
                             ("Run2024", 'RunIII2024Summer24'),
                             ("Run2025", 'RunIII2024Summer24')]:
            if campaign not in processor.histos['BsToJpsiPhi']:
                continue
            h_target = processor.histos['BsToJpsiPhi'][campaign]['B_All']
            compute_reweighted_efficiency(processor.histos, pd, run,
                                          "HLT_DoubleMu4_3_LowMass", "All", h_target)

        # Low Pt
        print("\nAverage efficiency of HLT_DoubleMu2_Jpsi_LowPt for low pt BsToJpsiPhi selection")
        pd = 'EGammaLowPt'
        h_target = processor.histos['BsToJpsiPhi']['RunIII2024Summer24']['B_AllLowPt']
        compute_efficiency(processor.histos, "BsToJpsiPhi", "RunIII2024Summer24",
                           "B_HLT_DoubleMu2_Jpsi_LowPt", "B_AllLowPt", "BsToJpsiPhi")
        compute_reweighted_efficiency(processor.histos, pd, "Run2024",
                                      "HLT_DoubleMu2_Jpsi_LowPt_AllLowPt", "AllLowPt", h_target)
        

    ##############################################################
    ## Plot distributions
    ##############################################################

    if make_plots:
        for sample_name in processor.histos:
            print(f"Sample name: {sample_name}")
            for run in ["Run2022", "Run2023", "Run2024", "Run2025", "Run3Summer22",
                        "Run3Summer23", "RunIII2024Summer24"]:
                h2 = aggregate_histograms(processor.histos, sample_name, run, "All")
                if h2:
                    h2.GetXaxis().SetTitle("p_{T}, GeV")
                    h2.GetYaxis().SetTitle("\eta")
                    h2.Draw("colz")
                    print_canvas(f"jpsi_all-{sample_name}_{run}", output_path)
                h_jpsi_mass = aggregate_histograms(processor.histos, sample_name, run, "jpsi_mass")
                if h_jpsi_mass:
                    h_jpsi_mass.SetLineWidth(2)
                    h_jpsi_mass.Draw()
                    print_canvas(f"jpsi_mass-{sample_name}_{run}", output_path)
    
            
    # Compute jpsi efficiency
    if compute_jpsi_efficiency:
        print("\nAbsolute HLT_DoubleMu4_3_LowMass efficiency for Jpsi selection")
        
        ##############################################################
        ## General jpsi efficiency
        ##############################################################
        h_target = processor.histos['BuToJpsiK']['RunIII2024Summer24']['All']

        pds = ['ZeroBias', 'ZeroBiasExclusive', 'EGamma', 'ParkingDoubleElectronLowMass']
        for pd in pds:
            for run in ["Run2022", "Run2023", "Run2024", "Run2025"]:
                if run != "Run2022" and pd == "ParkingDoubleElectronLowMass":
                    continue
                compute_reweighted_efficiency(processor.histos, pd, run,
                                   "HLT_DoubleMu4_3_LowMass", "All", h_target, pd=="EGamma")

        # mc22_n_num = processor.histos['BuToJpsiK']['Run3Summer22']['HLT_DoubleMu4_3_LowMass'].Integral()
        # mc22_n_den = processor.histos['BuToJpsiK']['Run3Summer22']['All'].Integral()
        # mc22_eff = mc22_n_num / mc22_n_den
        # mc22_eff_err = math.sqrt(mc22_eff*(1-mc22_eff)/mc22_n_den)
        # print(f"MC efficiency BuToJpsiK Run3Summer22 (HLT_DoubleMu4_3_LowMass): {mc22_eff:0.4f} +/- {mc22_eff_err:0.4f}")
        # mc22ee_n_num = processor.histos['BuToJpsiK']['Run3Summer22EE']['HLT_DoubleMu4_3_LowMass'].Integral()
        # mc22ee_n_den = processor.histos['BuToJpsiK']['Run3Summer22EE']['All'].Integral()
        # mc22ee_eff = mc22ee_n_num / mc22ee_n_den
        # mc22ee_eff_err = math.sqrt(mc22ee_eff*(1-mc22ee_eff)/mc22ee_n_den)
        # print(f"MC efficiency BuToJpsiK Run3Summer22EE (HLT_DoubleMu4_3_LowMass): {mc22ee_eff:0.4f} +/- {mc22ee_eff_err:0.4f}")
        
        # mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_DoubleMu4_3_LowMass'].Integral()
        # mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_DoubleMu4_3_LowMass_L1_DoubleMu_Mix_Unprescaled'].Integral()
        # mc_eff /= h_target.Integral()
        # print(f"MC efficiency BuToJpsiK RunIII2024Summer24 (HLT_DoubleMu4_3_LowMass_L1_DoubleMu_Mix_Unprescaled): {mc_eff:0.4f}")

        h_target = processor.histos['BuToJpsiK']['RunIII2024Summer24']['All']

        print("\nAbsolute HLT_Mu0_L1DoubleMu efficiency")
        pds = ['ZeroBias', 'EGamma']
        for pd in pds:
            for run in ["Run2022", "Run2023", "Run2024", "Run2025"]: 
                compute_reweighted_efficiency(processor.histos, pd, run,
                                   "HLT_Mu0_L1DoubleMu", "All", h_target, pd=="EGamma")

        mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_Mu0_L1DoubleMu'].Integral()
        mc_eff /= h_target.Integral()
        print(f"MC efficiency: {mc_eff:0.4f}")


        print("\nHLT_DoubleMu4_3_LowMass wrt HLT_Mu4_L1DoubleMu efficiency")
        pds = ['ParkingDoubleMuonLowMass']
        for pd in pds:
            for run in ["Run2022", "Run2023", "Run2024", "Run2025"]: 
                compute_reweighted_efficiency(processor.histos, pd, run,
                                   "HLT_Mu4_L1DoubleMu_HLT_DoubleMu4_3_LowMass", "HLT_Mu4_L1DoubleMu", h_target, True)

        mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_Mu4_L1DoubleMu_HLT_DoubleMu4_3_LowMass'].Integral()
        mc_eff /= processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_Mu4_L1DoubleMu'].Integral()
        print(f"MC efficiency: {mc_eff:0.4f}")


        print("\nHLT_DoubleMu4_3_LowMass wrt HLT_Mu0_L1DoubleMu efficiency")
        pds = ['HLTPhysics']
        for pd in pds:
            for run in ["Run2022", "Run2023", "Run2024", "Run2025"]: 
                compute_reweighted_efficiency(processor.histos, pd, run,
                                   "HLT_DoubleMu4_3_LowMass", "HLT_Mu0_L1DoubleMu", h_target, False)

        mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_DoubleMu4_3_LowMass'].Integral()
        mc_eff /= processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_Mu0_L1DoubleMu'].Integral()
        print(f"MC efficiency: {mc_eff:0.4f}")

        print("\nHLT_Mu0_L1DoubleMu wrt HLT_Mu8 efficiency")
        pds = ['Muon']
        for pd in pds:
            for run in ["Run2023", "Run2024", "Run2025"]: 
                compute_reweighted_efficiency(processor.histos, pd, run,
                                   "HLT_Mu8_HLT_Mu0_L1DoubleMu", "HLT_Mu8", h_target, False)

        mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_Mu8_HLT_Mu0_L1DoubleMu'].Integral()
        mc_eff /= processor.histos['BuToJpsiK']['RunIII2024Summer24']['HLT_Mu8'].Integral()
        print(f"MC efficiency: {mc_eff:0.4f}")


        print("\nAbsolute HLT_DoubleMu2_Jpsi_LowPt efficiency")
        pds = ['EGamma']
        for pd in pds:
            for run in ["Run2024", "Run2025"]:
                compute_reweighted_efficiency(processor.histos, pd, run,
                                   "lowpt_HLT_DoubleMu2_Jpsi_LowPt", "lowpt_All", h_target, False)

        mc_eff = processor.histos['BuToJpsiK']['RunIII2024Summer24']['lowpt_HLT_DoubleMu2_Jpsi_LowPt'].Integral()
        mc_eff /= processor.histos['BuToJpsiK']['RunIII2024Summer24']['lowpt_All'].Integral()
        print(f"MC efficiency: {mc_eff:0.4f}")


    #######################################################
    ## F-ratio correction 
    #######################################################

    if compute_corrections:
        for run, campaign in [ ("Run2022", "Run3Summer22"),
                               ("Run2023", "Run3Summer23"),
                               ("Run2024", "RunIII2024Summer24"),
                               ("Run2025", "RunIII2024Summer24"),
                              ]:
            if campaign not in processor.histos['BuToJpsiK']:
                continue
            h4_bu = processor.histos['BuToJpsiK'][campaign]['Bu']
            h4_bs = processor.histos['BsToJpsiPhi'][campaign]['B']
            h4_bd = processor.histos['BdToJpsiKstar'][campaign]['B']

            pd = "EGamma"
            mc_pd = "BuToJpsiK"


            # "HLT_DoubleMu4_3_LowMass_L1_DoubleMu_Real_Mix", "All", h2_bs)
            h_num_data, h_den_data = aggregate_two_sets_of_histograms(processor.histos, pd, run,
                                                          "HLT_DoubleMu4_3_LowMass", "All")
            h_num_mc, h_den_mc = aggregate_two_sets_of_histograms(processor.histos, mc_pd, campaign,
                                                      "HLT_DoubleMu4_3_LowMass", "All")

            print(f"pT scan for HLT_DoubleMu4_3_LowMass using {pd}-{run} and {mc_pd}-{campaign}")
            latex_report = ""
            results = {"Rsu":[], "Rdu":[]}
            results_for_analysis = {"Rsu":[], "Rdu":[]}
            for b_pt_bin in range(len(b_pt_bins) - 1):
                if debug:
                    print(b_pt_bins[b_pt_bin])

                h2_bs = thnd_to_th2(h4_bs, 2, 3, "h2_bs", False, {0: b_pt_bin + 1})
                bs_data_eff = weighted_efficiency(h_num_data, h_den_data, h2_bs)
                bs_mc_eff = weighted_efficiency(h_num_mc, h_den_mc, h2_bs)

                h2_bd = thnd_to_th2(h4_bd, 2, 3, "h2_bd", False, {0: b_pt_bin + 1})
                bd_data_eff = weighted_efficiency(h_num_data, h_den_data, h2_bd)
                bd_mc_eff = weighted_efficiency(h_num_mc, h_den_mc, h2_bd)
                
                h2_bu = thnd_to_th2(h4_bu, 2, 3, "h2_bu", False, {0: b_pt_bin + 1})
                bu_data_eff = weighted_efficiency(h_num_data, h_den_data, h2_bu)
                bu_mc_eff = weighted_efficiency(h_num_mc, h_den_mc, h2_bu)

                if debug:
                    print(f"Data efficiency for Bu: {bu_data_eff[0]:.4f} +/- {bu_data_eff[1]:.4f}")
                    print(f"MC efficiency for Bu:   {bu_mc_eff[0]:.4f} +/- {bu_mc_eff[1]:.4f}")
                    print(f"MC efficiency for Bs:   {bs_mc_eff[0]:.4f} +/- {bs_mc_eff[1]:.4f}")
                    print(f"Data efficiency for Bs: {bs_data_eff[0]:.4f} +/- {bs_data_eff[1]:.4f}")
                    print(f"(bs_data_eff/bs_mc_eff)/(bu_data_eff/bu_mc_eff): {(bs_data_eff[0]/bs_mc_eff[0])/(bu_data_eff[0]/bu_mc_eff[0])}")
                    print(f"MC efficiency for Bd:   {bd_mc_eff[0]:.4f} +/- {bd_mc_eff[1]:.4f}")
                    print(f"Data efficiency for Bd: {bd_data_eff[0]:.4f} +/- {bd_data_eff[1]:.4f}")
                    print(f"(bd_data_eff/bd_mc_eff)/(bu_data_eff/bu_mc_eff): {(bd_data_eff[0]/bd_mc_eff[0])/(bu_data_eff[0]/bu_mc_eff[0])}")

                r_bs_data, r_bs_data_err = weighted_efficiency_ratio(h_num_data, h_den_data, h2_bs, h2_bu)
                r_bd_data, r_bd_data_err = weighted_efficiency_ratio(h_num_data, h_den_data, h2_bd, h2_bu)
                r_bs_mc, r_bs_mc_err = weighted_efficiency_ratio(h_num_mc, h_den_mc, h2_bs, h2_bu)
                r_bd_mc, r_bd_mc_err = weighted_efficiency_ratio(h_num_mc, h_den_mc, h2_bd, h2_bu)
                corr_bs = r_bs_data / r_bs_mc
                corr_bs_err = corr_bs * math.sqrt((r_bs_data_err / r_bs_data)**2 + (r_bs_mc_err / r_bs_mc)**2)
                corr_bd = r_bd_data / r_bd_mc
                corr_bd_err = corr_bd * math.sqrt((r_bd_data_err / r_bd_data)**2 + (r_bd_mc_err / r_bd_mc)**2)
                if debug:
                    print(f"Data efficiency ratio (Bs/Bu): {r_bs_data:.4f} +/- {r_bs_data_err:.4f}")
                    print(f"MC efficiency ratio (Bs/Bu): {r_bs_mc:.4f} +/- {r_bs_mc_err:.4f}")
                    print(f"Rs correction (Bs/Bu): {corr_bs:.4f} +/- {corr_bs_err:.4f}")
                    print(f"Data efficiency ratio (Bd/Bu): {r_bd_data:.4f} +/- {r_bd_data_err:.4f}")
                    print(f"MC efficiency ratio (Bd/Bu): {r_bd_mc:.4f} +/- {r_bd_mc_err:.4f}")
                    print(f"Rs correction (Bd/Bu): {corr_bd:.4f} +/- {corr_bd_err:.4f}")
                latex_report += f"$[{b_pt_bins[b_pt_bin]:4.1f},{b_pt_bins[b_pt_bin+1]:4.1f}]$" + \
                    f"& ${bs_data_eff[0]*100:.2f}\pm{bs_data_eff[1]*100:.2f}$" + \
                    f"& ${bs_mc_eff[0]*100:.2f}\pm{bs_mc_eff[1]*100:.2f}$" + \
                    f"& ${bu_data_eff[0]*100:.2f}\pm{bu_data_eff[1]*100:.2f}$" + \
                    f"& ${bu_mc_eff[0]*100:.2f}\pm{bu_mc_eff[1]*100:.2f}$" + \
                    f"& ${corr_bs:.4f}\pm{corr_bs_err:.4f}$ \\\\\n"
                results["Rsu"].append((
                    (b_pt_bins[b_pt_bin]+b_pt_bins[b_pt_bin+1])/2, corr_bs,
                    (b_pt_bins[b_pt_bin+1]-b_pt_bins[b_pt_bin])/2, corr_bs_err))
                results["Rdu"].append((
                    (b_pt_bins[b_pt_bin]+b_pt_bins[b_pt_bin+1])/2, corr_bd,
                    (b_pt_bins[b_pt_bin+1]-b_pt_bins[b_pt_bin])/2, corr_bd_err))
                results_for_analysis["Rsu"].append(((
                    b_pt_bins[b_pt_bin], b_pt_bins[b_pt_bin+1]), (corr_bs, corr_bs_err)))
                results_for_analysis["Rdu"].append(((
                    b_pt_bins[b_pt_bin], b_pt_bins[b_pt_bin+1]), (corr_bd, corr_bd_err)))

            print(f"Ratio correction using {pd}-{run} and {mc_pd}-{campaign}")
            print(latex_report)
            for btype, result in results.items():
                n = len(result)
                gr = ROOT.TGraphErrors(n)
                for i in range(n):
                    gr.SetPoint(i, result[i][0], result[i][1])
                    gr.SetPointError(i, result[i][2], result[i][3])
                gr.SetMarkerStyle(20)
                gr.SetMarkerSize(1)
                gr.SetLineColor(ROOT.kBlue)
                gr.SetMinimum(0.9)
                gr.SetMaximum(1.1)
                gr.Draw("AP")
                gr.Fit("pol1","EX0","",10,20)
                print_canvas(f"{btype}_correction_nom-{run}", output_path)
                
                with open(f"{output_path}/{btype}_correction_nom-{run}.json", "w") as f:
                    json.dump(results_for_analysis[btype], f, indent=4)

            # Low Pt
            if run in ["Run2024", "Run2025"]:
                print(f"pT scan for HLT_DoubleMu2_Jpsi_LowPt using {pd}-{run} and {mc_pd}-{campaign}")

                h4_bu = processor.histos['BuToJpsiK'][campaign]['BuLowPt']
                h4_bs = processor.histos['BsToJpsiPhi'][campaign]['BsLowPt']

                pd = "EGammaLowPt"
                mc_pd = "BuToJpsiK"

                h_num_data, h_den_data = aggregate_two_sets_of_histograms(processor.histos, pd, run,
                                                              "HLT_DoubleMu2_Jpsi_LowPt_AllLowPt", "AllLowPt")
                h_num_mc, h_den_mc = aggregate_two_sets_of_histograms(processor.histos, mc_pd, campaign,
                                                          "HLT_DoubleMu2_Jpsi_LowPt_AllLowPt", "AllLowPt")
                latex_report = ""
                results_for_graph = []
                results_for_analysis = []
                for b_pt_bin in range(len(b_low_pt_bins) - 1):
                    if debug:
                        print(b_low_pt_bins[b_pt_bin])

                    h2_bs = thnd_to_th2(h4_bs, 2, 3, "h2_bs", False, {0: b_pt_bin + 1})
                    bs_data_eff = weighted_efficiency(h_num_data, h_den_data, h2_bs)
                    bs_mc_eff = weighted_efficiency(h_num_mc, h_den_mc, h2_bs)
                    
                    h2_bu = thnd_to_th2(h4_bu, 2, 3, "h2_bu", False, {0: b_pt_bin + 1})
                    bu_data_eff = weighted_efficiency(h_num_data, h_den_data, h2_bu)
                    bu_mc_eff = weighted_efficiency(h_num_mc, h_den_mc, h2_bu)

                    skip = False
                    if bs_mc_eff[0] == None:
                        print(f"Cannot get MC efficiency for Bs. Continue.")
                        skip = True
                    if bs_data_eff[0] == None:
                        print(f"Cannot get Data efficiency for Bs. Continue.")
                        skip = True
                    if bu_data_eff[0] == None:
                        print(f"Cannot get Data efficiency for Bu. Continue.")
                        skip = True
                    if bu_mc_eff[0] == None:
                        print(f"Cannot get MC efficiency for Bu. Continue.")
                        skip = True

                    if skip:
                        continue

                    if debug:
                        print(f"MC efficiency for Bs:   {bs_mc_eff[0]:.4f} +/- {bs_mc_eff[1]:.4f}")
                        print(f"Data efficiency for Bs: {bs_data_eff[0]:.4f} +/- {bs_data_eff[1]:.4f}")
                        print(f"Data efficiency for Bu: {bu_data_eff[0]:.4f} +/- {bu_data_eff[1]:.4f}")
                        print(f"MC efficiency for Bu:   {bu_mc_eff[0]:.4f} +/- {bu_mc_eff[1]:.4f}")
                        print(f"(bs_data_eff/bs_mc_eff)/(bu_data_eff/bu_mc_eff): {(bs_data_eff[0]/bs_mc_eff[0])/(bu_data_eff[0]/bu_mc_eff[0])}")

                    r_data, r_data_err = weighted_efficiency_ratio(h_num_data, h_den_data, h2_bs, h2_bu)
                    r_mc, r_mc_err = weighted_efficiency_ratio(h_num_mc, h_den_mc, h2_bs, h2_bu)
                    corr = r_data / r_mc
                    corr_err = corr * math.sqrt((r_data_err / r_data)**2 + (r_mc_err / r_mc)**2)
                    if debug:
                        print(f"Data efficiency ratio: {r_data:.4f} +/- {r_data_err:.4f}")
                        print(f"MC efficiency ratio: {r_mc:.4f} +/- {r_mc_err:.4f}")
                        print(f"Rs correction: {corr:.4f} +/- {corr_err:.4f}")
                    latex_report += f"$[{b_low_pt_bins[b_pt_bin]:4.1f},{b_low_pt_bins[b_pt_bin+1]:4.1f}]$" + \
                        f"& ${bs_data_eff[0]*100:.2f}\pm{bs_data_eff[1]*100:.2f}$" + \
                        f"& ${bs_mc_eff[0]*100:.2f}\pm{bs_mc_eff[1]*100:.2f}$" + \
                        f"& ${bu_data_eff[0]*100:.2f}\pm{bu_data_eff[1]*100:.2f}$" + \
                        f"& ${bu_mc_eff[0]*100:.2f}\pm{bu_mc_eff[1]*100:.2f}$" + \
                        f"& ${corr:.4f}\pm{corr_err:.4f}$ \\\\\n"
                    results_for_graph.append((
                        (b_low_pt_bins[b_pt_bin]+b_low_pt_bins[b_pt_bin+1])/2, corr,
                        (b_low_pt_bins[b_pt_bin+1]-b_low_pt_bins[b_pt_bin])/2, corr_err))
                    results_for_analysis.append(((b_low_pt_bins[b_pt_bin], b_low_pt_bins[b_pt_bin+1]), (corr, corr_err)))

                print(f"Ratio correction using {pd}-{run} and {mc_pd}-{campaign}")
                print(latex_report)
                n = len(results_for_graph)
                gr = ROOT.TGraphErrors(n)
                for i in range(n):
                    gr.SetPoint(i, results_for_graph[i][0], results_for_graph[i][1])
                    gr.SetPointError(i, results_for_graph[i][2], results_for_graph[i][3])
                gr.SetMarkerStyle(20)
                gr.SetMarkerSize(1)
                gr.SetLineColor(ROOT.kBlue)
                gr.SetMinimum(0.8)
                gr.SetMaximum(1.2)
                gr.Draw("AP")
                # f = ROOT.TF1("f", "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])*(x-[3])/[4])", 7.0, 24.0);
                gr.Fit("pol3","EX0")
                print_canvas(f"rs_correction_lowpt-{run}", output_path)
                with open(f"{output_path}/rs_correction_lowpt-{run}.json", "w") as f:
                    json.dump(results_for_analysis, f, indent=4)

            
    #     self.samples["BuToJpsiK"]["files"]["RunIII2024Summer24"].append(
    #         path + "/BuToJpsiK_Fil-BMuon_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2+MINIAODSIM/*.root"
    #     )

        # BsToJPsiPhi-JpsiToMuMu-PhiToKK_Fil-MuPt2_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen
    
    
    # Use flat ntuples
    

    # fill_histograms(force=force_recreate)

    # ## Perform efficiency study
    # histos = load_histograms(args.outfile)

    # ## Save results
    # json.dump(results, open(output, 'w'), indent=4)
