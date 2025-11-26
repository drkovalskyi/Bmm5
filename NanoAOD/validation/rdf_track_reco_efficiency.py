"""Track reconstruction efficiency measurement

The script measures the track reconstruction efficiency and its
corrections using the Tag-n-Probe method, with stand-alone muons used
as a probe in Jpsi -> mumu decays.

The input data must be in Bmm5 NanoAOD format.

The script uses RDataFrame for optimal performance, but it needs to
read large volumes of data to collect all information. To speed up
processing, it first extracts histograms and saves them in a root
file. This file is then reused as input for further processing and for
reruns, unless the script is instructed to recreate it.
"""

force_recreate = False
aggregate_eras = False
force_data_fits = True
do_eta_fits = False
do_pt_fits = False
compute_corrections = True
debug = False
# process_only = "Run2024C|RunIII2024Summer24"
process_only = "Run2022|Run3Summer22EE"
# process_only = None
# max_files = 10
max_files = 999999
path   = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/Skims/535/trig/"
path2  = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/535/"
# output_path = "/eos/home-d/dmytro/www/plots/2025/fsfu-track_reco_efficiency/"
output_path = "/eos/home-d/dmytro/www/plots/tmp/2025/track_reco_efficiency/"
file_fit_results = f"{output_path}/fit_results.json"

eta_bins = [-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]
pt_bins  = [3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0, 20.0]
min_mass = 2.0
max_mass = 5.0
min_mass_tag_tag = 2.5
max_mass_tag_tag = 4.0
nbins = 100
mass_bins = [min_mass + i * (max_mass - min_mass) / nbins for i in range(nbins + 1)]
fine_mass_bins = [min_mass + i * (max_mass - min_mass) / nbins / 2 for i in range(2 * nbins + 1)]

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

class DataProcessor:
    def __init__(self):
        self.report = dict()
        self.samples = dict()
        self.histos = dict()

        self.hist_file = "rdf_track_reco_efficiency.root"

        self.load_histograms()

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

        # Data
        data_files = glob.glob(f"{path}/ParkingDoubleMuonLowMass*202[2345]*/*.root")
        pattern = re.compile(r"(Run202\d\w)")
        for file in data_files:
            match = pattern.search(file)
            if match:
                era = match.group(1)
                if "ParkingDoubleMuonLowMass" not in self.samples:
                    self.samples["ParkingDoubleMuonLowMass"] = {
                        "Data": True,
                        "files": defaultdict(list),
                        "triggers": ["HLT_Mu4_L1DoubleMu"]
                    }
                self.samples["ParkingDoubleMuonLowMass"]["files"][era].append(file)

        # MC
        # data_files = glob.glob(f"{path2}/BuToJpsiK*/*.root")
        data_files = glob.glob(f"{path2}/BuToJpsiPi*/*.root")
        pattern = re.compile(r"(Run3Summer22EE|Run3Summer22|Run3Summer23BPix|Run3Summer23|RunIII2024Summer24)")
        for file in data_files:
            match = pattern.search(file)
            if match:
                era = match.group(1)
                if "BuToJpsiK" not in self.samples:
                    self.samples["BuToJpsiK"] = {
                        "Data": False,
                        "files": defaultdict(list),
                        "triggers": []
                    }
                self.samples["BuToJpsiK"]["files"][era].append(file)

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
        self.era2campaign = dict()
        for run, campaigns in self.campaigns.items():
            for campaign, eras in campaigns.items():
                for era in eras:
                    self.era2campaign[era] = campaign
        if debug:
            for sample, info in self.samples.items():
                print(sample)
                for era, files in info["files"].items():
                    print(f"\t{era}: {len(files)}")
        
    @staticmethod
    def load_data(file_list, tree_name="Events"):
        assert isinstance(file_list, list), "expected a list"
        chain = ROOT.TChain(tree_name)
        print("Number of files:", len(file_list))
        if len(file_list) == 0:
            return None
        for file in file_list[:max_files]:
            chain.Add(file)
        n_files = chain.GetListOfFiles().GetEntries()
        print("Number of files in the chain:", n_files)
        if n_files == 0:
            return None
        return chain

    def book_histo_ND(self, rdf_histos, rdf, sample_name, era, histo_name, bins, vars):
        full_histo_name = self.histo_key(sample_name, era, histo_name)
        nbins = [len(e) - 1 for e in bins]
        model = ROOT.RDF.THnDModel(full_histo_name, "", len(bins), nbins, bins)
        rdf_histos[histo_name] = rdf.HistoND(model, vars)
        
    def book_histo(self, rdf_histos, rdf, sample_name, era, histo_name, bins, vars):
        full_histo_name = self.histo_key(sample_name, era, histo_name)

        if len(bins) == 1:
            rdf_histos[histo_name] = \
                rdf.Histo1D((full_histo_name, "",
                             len(bins[0]) - 1, array('d', bins[0])),
                            vars[0])
        elif len(bins) == 2:
            rdf_histos[histo_name] = \
                rdf.Histo2D((full_histo_name, "",
                             len(bins[0]) - 1, array('d', bins[0]), 
                             len(bins[1]) - 1, array('d', bins[1])), 
                            vars[0], vars[1])
        elif len(bins) == 3:
            rdf_histos[histo_name] = \
                rdf.Histo3D((full_histo_name, "",
                             len(bins[0]) - 1, array('d', bins[0]), 
                             len(bins[1]) - 1, array('d', bins[1]), 
                             len(bins[2]) - 1, array('d', bins[2])), 
                            vars[0], vars[1], vars[2])
        elif len(bins) > 3: 
            nbins = [len(e) - 1 for e in bins]
            model = ROOT.RDF.THnDModel(full_histo_name, "", len(bins), nbins, bins)
            rdf_histos[histo_name] = rdf.HistoND(model, vars)
        
    def book_histo_2D(self, rdf_histos, rdf, sample_name, era, histo_name):
        self.book_histo(rdf_histos, rdf, sample_name, era, histo_name,
                        [jpsi_pt_bins, eta_bins],
                        ["jpsi_pt", "jpsi_eta"])

    def book_histo_3D(self, rdf_histos, rdf, sample_name, era, histo_name):
        self.book_histo(rdf_histos, rdf, sample_name, era, histo_name,
                        [mass_bins, pt_bins, eta_bins],
                        ["mass", "probe_pt", "probe_eta"])

    def process_sample(self, sample_name, era):
        sample = self.samples[sample_name]
        rdf_histos = dict()
        
        ### Prepare RDataFrame

        file_list = []
        for sample_era, files in sample["files"].items():
            if sample_era != era:
                continue
            file_list.extend(files)
        if len(file_list) == 0:
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

        chain = self.load_data(file_list)
        if chain == None:
            return dict()

        rdf = ROOT.RDataFrame(chain)
        # Make sure all triggers have a default value
        triggers = [
            "HLT_Mu4_L1DoubleMu"
        ]
        for trigger in triggers:
            rdf = rdf.DefaultValueFor(trigger, False)

        ### Extract information from other branches
        # tnp_probe1_tag2&&MuonId_HLT_Mu4_L1DoubleMu[tnp_mu2_index]
        
        # mm
        rdf = rdf.Define("tnp_probe1_tag2_trig", "Take(MuonId_HLT_Mu4_L1DoubleMu,  tnp_mu2_index)")
        rdf = rdf.Define("tnp_probe2_tag1_trig", "Take(MuonId_HLT_Mu4_L1DoubleMu,  tnp_mu1_index)")

        # Offline selection
        rdf = rdf.Define("probes_1",        "tnp_probe1_tag2 && tnp_probe1_tag2_trig")
        rdf = rdf.Define("probes_1_failed", "tnp_probe1_tag2 && tnp_probe1_tag2_trig && !tnp_mu1_test")
        rdf = rdf.Define("probes_1_sig",
                         "tnp_probe1_tag2 && tnp_probe1_tag2_trig && tnp_mu1_tag && abs(tnp_mass-3.09)<0.05")
        rdf = rdf.Define("probes_2",        "tnp_probe2_tag1 && tnp_probe2_tag1_trig")
        rdf = rdf.Define("probes_2_failed", "tnp_probe2_tag1 && tnp_probe2_tag1_trig && !tnp_mu2_test")
        rdf = rdf.Define("probes_2_sig",
                         "tnp_probe2_tag1 && tnp_probe2_tag1_trig && tnp_mu2_tag && abs(tnp_mass-3.09)<0.05")
        rdf = rdf.Define("tag_tag",         "tnp_mu1_tag && tnp_mu2_tag")

        rdf1 = rdf.Filter("Sum(probes_1)>0")
        rdf1 = rdf1.Define("probe_pt",     "tnp_mu1_sta_pt[probes_1]")
        rdf1 = rdf1.Define("probe_eta",    "tnp_mu1_eta[probes_1]")
        rdf1 = rdf1.Define("probe_phi",    "tnp_mu1_phi[probes_1]")
        rdf1 = rdf1.Define("mass",         "tnp_probe1_tag2_mass[probes_1]")
        
        rdf1_failed = rdf.Filter("Sum(probes_1_failed)>0")
        rdf1_failed = rdf1_failed.Define("probe_pt",     "tnp_mu1_sta_pt[probes_1_failed]")
        rdf1_failed = rdf1_failed.Define("probe_eta",    "tnp_mu1_eta[probes_1_failed]")
        rdf1_failed = rdf1_failed.Define("probe_phi",    "tnp_mu1_phi[probes_1_failed]")
        rdf1_failed = rdf1_failed.Define("mass",         "tnp_probe1_tag2_mass[probes_1_failed]")

        rdf1_sig = rdf.Filter("Sum(probes_1_sig)>0")
        rdf1_sig = rdf1_sig.Define("probe_pt",     "tnp_mu1_sta_pt[probes_1_sig]")
        rdf1_sig = rdf1_sig.Define("probe_eta",    "tnp_mu1_eta[probes_1_sig]")
        rdf1_sig = rdf1_sig.Define("probe_phi",    "tnp_mu1_phi[probes_1_sig]")
        rdf1_sig = rdf1_sig.Define("mass",         "tnp_probe1_tag2_mass[probes_1_sig]")
        
        rdf2 = rdf.Filter("Sum(probes_2)>0")
        rdf2 = rdf2.Define("probe_pt",     "tnp_mu2_sta_pt[probes_2]")
        rdf2 = rdf2.Define("probe_eta",    "tnp_mu2_eta[probes_2]")
        rdf2 = rdf2.Define("probe_phi",    "tnp_mu2_phi[probes_2]")
        rdf2 = rdf2.Define("mass",         "tnp_probe2_tag1_mass[probes_2]")

        rdf2_failed = rdf.Filter("Sum(probes_2_failed)>0")
        rdf2_failed = rdf2_failed.Define("probe_pt",     "tnp_mu2_sta_pt[probes_2_failed]")
        rdf2_failed = rdf2_failed.Define("probe_eta",    "tnp_mu2_eta[probes_2_failed]")
        rdf2_failed = rdf2_failed.Define("probe_phi",    "tnp_mu2_phi[probes_2_failed]")
        rdf2_failed = rdf2_failed.Define("mass",         "tnp_probe2_tag1_mass[probes_2_failed]")
        
        rdf2_sig = rdf.Filter("Sum(probes_2_sig)>0")
        rdf2_sig = rdf2_sig.Define("probe_pt",     "tnp_mu2_sta_pt[probes_2_sig]")
        rdf2_sig = rdf2_sig.Define("probe_eta",    "tnp_mu2_eta[probes_2_sig]")
        rdf2_sig = rdf2_sig.Define("probe_phi",    "tnp_mu2_phi[probes_2_sig]")
        rdf2_sig = rdf2_sig.Define("mass",         "tnp_probe2_tag1_mass[probes_2_sig]")
        
        rdf3 = rdf.Filter("Sum(tag_tag)>0")
        rdf3 = rdf3.Define("mass",         "tnp_mass[tag_tag]")
        
        self.book_histo_3D(rdf_histos, rdf1, sample_name, era, "probe1")
        self.book_histo_3D(rdf_histos, rdf1_failed, sample_name, era, "probe1_failed")
        self.book_histo_3D(rdf_histos, rdf1_sig, sample_name, era, "probe1_sig")
        self.book_histo_3D(rdf_histos, rdf2, sample_name, era, "probe2")
        self.book_histo_3D(rdf_histos, rdf2_failed, sample_name, era, "probe2_failed")
        self.book_histo_3D(rdf_histos, rdf2_sig, sample_name, era, "probe2_sig")
        self.book_histo(rdf_histos, rdf3, sample_name, era, "tag_tag", [fine_mass_bins], ["mass"])
        
        if sample_name not in self.histos:
            self.histos[sample_name] = dict()
        if era not in self.histos[sample_name]:
            self.histos[sample_name][era] = dict()
        for hist_name, rdf_hist in rdf_histos.items():
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


    def process_samples(self):
        
        self.define_samples()
        
        if force_recreate and os.path.exists(self.hist_file):
            os.remove(self.hist_file)
            
        # process samples
        for sample_name in self.samples:
            if not force_recreate and sample_name in self.histos:
                continue
            for era in self.samples[sample_name]["files"]:
                if process_only and not re.search(process_only, era):
                    continue
                print(f"Processing {sample_name} - {era}")
                self.process_sample(sample_name, era)

        self.make_report()

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    """Print canvas in different formats"""

    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

def make_dataset(name, mass_var, hist):
    """Build RooDataSet from a histogram"""
    
    data = ROOT.RooDataHist(name, "", ROOT.RooArgList(mass_var), hist)

    return data

def build_model(mass_var):
    """Build fit model and save it in a workspace"""

    peak = 3.096
    peak2 = 3.686
    delta = ROOT.RooConstVar("delta", "mass difference", peak2 - peak)
    search_width = 0.1
    
    jpsi_mean   = ROOT.RooRealVar("jpsi_mean", "mu", peak, peak - search_width, peak + search_width)
    jpsi_lambda = ROOT.RooRealVar("jpsi_lambda", "lambda", 0.1, 0.01, 0.3) # sigma
    jpsi_gamma  = ROOT.RooRealVar("jpsi_gamma", "gamma", 1, -5, 5) # skewness
    jpsi_delta  = ROOT.RooRealVar("jpsi_delta", "delta", 1, 0.1, 3) # larger value smaller tails
    # jpsi_john = ROOT.RooJohnson("jpsi_john", "signal", mass_var, jpsi_mean, jpsi_lambda, jpsi_gamma, jpsi_delta)
    jpsi = ROOT.RooJohnson("jpsi", "signal", mass_var, jpsi_mean, jpsi_lambda, jpsi_gamma, jpsi_delta)
    # jpsi_G2_sigma = ROOT.RooRealVar("jpsi_G2_sigma", "", 0.03, 0.01, 0.50)
    # jpsi_G2_mean = ROOT.RooRealVar("jpsi_G2_mean", "", peak, peak - search_width, peak + search_width)
    # jpsi_G2       = ROOT.RooGaussian("jpsi_G2", "", mass_var, jpsi_G2_mean, jpsi_G2_sigma)
    # jpsi_G2_frac = ROOT.RooRealVar("jpsi_G2_frac","",0.3,0.0,1.0)
    # jpsi = ROOT.RooAddPdf("jpsi"," ", ROOT.RooArgList(jpsi_G2,jpsi_john), ROOT.RooArgList(jpsi_G2_frac))

    psi2S_mean = ROOT.RooFormulaVar("psi2S_mean", "@0 + @1", ROOT.RooArgList(jpsi_mean, delta))
    psi2S = ROOT.RooJohnson("psi2S", "", mass_var, psi2S_mean, jpsi_lambda, jpsi_gamma, jpsi_delta)
    # psi2S_john = ROOT.RooJohnson("psi2S_john", "", mass_var, psi2S_mean, jpsi_lambda, jpsi_gamma, jpsi_delta)
    # psi2S_G2 = ROOT.RooGaussian("psi2S_G2", "", mass_var, psi2S_mean, jpsi_G2_sigma)
    # psi2S = ROOT.RooAddPdf("psi2S"," ", ROOT.RooArgList(psi2S_G2, psi2S_john), ROOT.RooArgList(jpsi_G2_frac))
    
    psi2S_frac = ROOT.RooRealVar("psi2S_frac","",0.03,0.0,0.1)
    sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(psi2S, jpsi), ROOT.RooArgList(psi2S_frac))
    
    # # multi-gaussian
    # G1_mean  = ROOT.RooRealVar("sig_G1_mean",  "", peak, peak - search_width, peak + search_width)
    # G1_sigma = ROOT.RooRealVar("sig_G1_sigma", "", 0.03, 0.001, 0.20)
    # G2_scale = ROOT.RooRealVar("sig_G2_scale", "", 2.5, 0.2, 7.5)
    # G3_scale = ROOT.RooRealVar("sig_G3_scale", "", 3.0, 0.5, 6.7)
    # G2_sigma = ROOT.RooProduct("sig_G2_sigma", "", ROOT.RooArgList(G1_sigma,G2_scale))
    # G3_sigma = ROOT.RooProduct("sig_G3_sigma", "", ROOT.RooArgList(G1_sigma,G3_scale))
    # G1 = ROOT.RooGaussian("sig_G1", "", mass_var, G1_mean, G1_sigma)
    # G2 = ROOT.RooGaussian("sig_G2", "", mass_var, G1_mean, G2_sigma)
    # G3 = ROOT.RooGaussian("sig_G3", "", mass_var, G1_mean, G3_sigma)
    
    # G2_frac = ROOT.RooRealVar("sig_G2_frac","",0.3,0.0,1.0)
    # G3_frac = ROOT.RooRealVar("sig_G3_frac","",0.2,0.0,1.0)
    # # sig = ROOT.RooAddPdf("sig"," ", ROOT.RooArgList(G2,G1), ROOT.RooArgList(G2_frac))
    # sig  = ROOT.RooAddPdf("sig"," ",ROOT.RooArgList(G3,G2,G1),ROOT.RooArgList(G2_frac,G3_frac))

    ## Combinatorial background
    
    # a0    = ROOT.RooRealVar("a0", "a0", 0.0, -1.0,  1.0)
    # a1    = ROOT.RooRealVar("a1", "a1", 0.0, -0.3, 0.3)
    # bkg   = ROOT.RooChebychev("bkg", "Background", mass_var, ROOT.RooArgList(a0, a1))
    
    c0 = ROOT.RooRealVar("c0", "c0", 1.0, 0.0, 10.0)
    c1 = ROOT.RooRealVar("c1", "c1", 1.0, 0.0, 10.0)
    c2 = ROOT.RooRealVar("c2", "c2", 0.0, 0.0, 1.0)
    bkg = ROOT.RooBernstein("bkg", "Bernstein background", mass_var, ROOT.RooArgList(c0, c1, c2))
    # bkg = ROOT.RooBernstein("bkg", "Bernstein background", mass_var, ROOT.RooArgList(c0, c1))
    
    Nsig  = ROOT.RooRealVar("Nsig", "Nsig", 1000, 0, 1e9)
    # Nsig2  = ROOT.RooRealVar("Nsig2", "Nsig2", 1000, 0, 1e9)
    Nbkg  = ROOT.RooRealVar("Nbkg", "Nbkg", 0, 0, 1e9)
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sig,bkg), ROOT.RooArgList(Nsig,Nbkg))

    ws = ROOT.RooWorkspace("ws","")
    getattr(ws,'import')(model)

    return ws

def fit_data(ws, filename_prefix, subdirectory,
             hist_probe, hist_failed, hist_ref=None, hist_tag_tag=None,
             fix_shape_to_ref_all_probes=True, fix_shape_to_ref_failed_probes=True):
    results = dict()
    mass_var = ws.var("m")
    
    data_probe = make_dataset("data", mass_var, hist_probe)
    data_failed = make_dataset("data_failed", mass_var, hist_failed)
    
    model = ws.pdf("model")
    latex_list = []
    
    psi2S_frac = None
    ## fit tag-tag to get psi2S fraction
    if hist_tag_tag != None:
        data_tag_tag = make_dataset("data_tag_tag", mass_var, hist_tag_tag)
        ws.var("jpsi_lambda").setConstant(False)
        ws.var("jpsi_lambda").setRange(0.05, 0.2)
        ws.var("jpsi_lambda").setVal(0.01)
        ws.var("jpsi_gamma").setConstant(False)
        ws.var("jpsi_delta").setConstant(False) 
        ws.var("psi2S_frac").setVal(0.002)
        ws.var("psi2S_frac").setConstant(False)
        fit_result = model.fitTo(data_tag_tag, ROOT.RooFit.Save(),
                                 ROOT.RooFit.Range(min_mass_tag_tag, max_mass_tag_tag))
        fit_result.Print()
        # plot results
        frame = mass_var.frame()
        data_tag_tag.plotOn(frame)
        model.plotOn(frame)
        frame.Draw()
        print_canvas(f"{filename_prefix}_tag-tag", f"{output_path}/{subdirectory}")
        psi2S_frac = ws.var("psi2S_frac").getVal()
        latex_list.append("f_{\psi(2S)} = %0.3f \pm %0.3f" % (ws.var("psi2S_frac").getVal(), ws.var("psi2S_frac").getError()))
        # reset
        ws.var("jpsi_lambda").setRange(0.01, 0.3)
    
    ## prefit
    if hist_ref != None:
        # fit signal shape
        data_ref = make_dataset("data_ref", mass_var, hist_ref)
        jpsi = ws.pdf("jpsi")
        ws.var("jpsi_lambda").setConstant(False)
        ws.var("jpsi_gamma").setConstant(False)
        ws.var("jpsi_delta").setConstant(False)
        # ws.var("jpsi_G2_mean").setConstant(False)
        # ws.var("jpsi_G2_sigma").setConstant(False)
        # ws.var("jpsi_G2_frac").setConstant(False)
        ws.var("psi2S_frac").setVal(0.0)
        ws.var("psi2S_frac").setConstant(True)
        fit_result = jpsi.fitTo(data_ref, ROOT.RooFit.Save())
        fit_result.Print()
        # plot results
        frame = mass_var.frame()
        data_ref.plotOn(frame)
        jpsi.plotOn(frame)
        frame.Draw()
        print_canvas(f"{filename_prefix}_ref", f"{output_path}/{subdirectory}")
        
        # fix shape
        ws.var("jpsi_lambda").setConstant(True)
        ws.var("jpsi_gamma").setConstant(True)
        ws.var("jpsi_delta").setConstant(True)
        if hist_tag_tag != None:
            ws.var("psi2S_frac").setVal(psi2S_frac)
            ws.var("psi2S_frac").setConstant(True)
        else:
            ws.var("psi2S_frac").setConstant(False)
        
        ws.var("Nsig").setVal(hist_probe.Integral() * 0.7)
        ws.var("Nbkg").setVal(hist_probe.Integral() * 0.3)
        # ws.var("c2").setConstant(True)
        # ws.var("c2").setVal(0.0)
        fit_result = model.fitTo(data_probe, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Save())
        nDOF = get_number_of_free_parameters(model, mass_var)
        chi2ndof = frame.chiSquare(nDOF)

        latex_list.append("N^{fix}_{sig} = %0.0f \pm %0.0f" % (ws.var("Nsig").getVal(), ws.var("Nsig").getError()))
        latex_list.append("N^{fix}_{bkg} = %0.0f \pm %0.0f" % (ws.var("Nbkg").getVal(), ws.var("Nbkg").getError()))
        latex_list.append("#chi^{2}_{fix}/nDOF = %0.1f" % (chi2ndof))
        
        # ws.var("c2").setConstant(False)
    else:
        ws.var("Nsig").setVal(hist_probe.Integral() * 0.7)
        ws.var("Nbkg").setVal(hist_probe.Integral() * 0.3)
        fit_result = model.fitTo(data_probe, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Save())
        
    if not fix_shape_to_ref_all_probes:
        ws.var("jpsi_lambda").setConstant(False)
        ws.var("jpsi_gamma").setConstant(False)
        ws.var("jpsi_delta").setConstant(False)

        fit_result = model.fitTo(data_probe, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Save())
        nDOF = get_number_of_free_parameters(model, mass_var)
        chi2ndof = frame.chiSquare(nDOF)
        print(f"chiSquare: {chi2ndof:0.2f} nDOF: {nDOF}")
        latex_list.append("N^_{sig} = %0.0f \pm %0.0f" % (ws.var("Nsig").getVal(), ws.var("Nsig").getError()))
        latex_list.append("N_{bkg} = %0.0f \pm %0.0f" % (ws.var("Nbkg").getVal(), ws.var("Nbkg").getError()))
        latex_list.append("#chi^{2}/nDOF = %0.1f" % (chi2ndof))
    fit_result.Print()
    
    ## Plot results
    frame = mass_var.frame()
    data_probe.plotOn(frame)
    model.plotOn(frame, ROOT.RooFit.Components("bkg"),  ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame)
    frame.GetYaxis().SetTitleOffset(1.6)
    
    frame.Draw()
    latex = ROOT.TLatex()
    latex.SetTextSize(0.02)
    latex.SetTextFont(42) # Helvetica-like font, Greek works
    for i, entry in enumerate(latex_list):
        latex.DrawLatexNDC(0.65, 0.85 - i * 0.05, entry) 
    
    print_canvas(f"{filename_prefix}_probe", f"{output_path}/{subdirectory}")
    results["all_probe_fit"] = {
        "Nsig": ws.var("Nsig").getVal(),
        "Nsig_err": ws.var("Nsig").getError(),
        "chiSquare": chi2ndof
    }

    ## failed probe fit
    latex_list.clear()
    frame = mass_var.frame()
    ws.var("Nsig").setVal(ws.var("Nsig").getVal() * 0.05)
    
    ws.var("jpsi_lambda").setConstant(True)
    ws.var("jpsi_gamma").setConstant(True)
    ws.var("jpsi_delta").setConstant(True)
    ws.var("psi2S_frac").setConstant(True)
    fit_result = model.fitTo(data_failed, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Save())
    fit_result.Print()
    latex_list.append("N^{fix}_{sig} = %0.0f \pm %0.0f" % (ws.var("Nsig").getVal(), ws.var("Nsig").getError()))
    latex_list.append("N^{fix}_{bkg} = %0.0f \pm %0.0f" % (ws.var("Nbkg").getVal(), ws.var("Nbkg").getError()))
    latex_list.append("#chi_{fix}^2/nDOF = %0.1f" % (chi2ndof))
    
    if not fix_shape_to_ref_failed_probes:
        ws.var("jpsi_lambda").setConstant(False)
        ws.var("jpsi_gamma").setConstant(False)
        ws.var("jpsi_delta").setConstant(False)
        ws.var("c0").setConstant(True)
        ws.var("c1").setConstant(True)
        ws.var("c2").setConstant(True)

        fit_result = model.fitTo(data_failed, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Save())
        fit_result.Print()
        ws.var("c0").setConstant(False)
        ws.var("c1").setConstant(False)
        ws.var("c2").setConstant(False)
        
        fit_result = model.fitTo(data_failed, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Save())
        fit_result.Print()
        nDOF = get_number_of_free_parameters(model, mass_var)
        chi2ndof = frame.chiSquare(nDOF)
        print(f"chiSquare: {chi2ndof:0.2f} nDOF: {nDOF}")
        latex_list.append("N^_{sig} = %0.0f \pm %0.0f" % (ws.var("Nsig").getVal(), ws.var("Nsig").getError()))
        latex_list.append("N_{bkg} = %0.0f \pm %0.0f" % (ws.var("Nbkg").getVal(), ws.var("Nbkg").getError()))
        latex_list.append("#chi^{2}/nDOF = %0.1f" % (chi2ndof))
    
    data_failed.plotOn(frame)
    model.plotOn(frame, ROOT.RooFit.Components("bkg"),  ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame)
    nDOF = get_number_of_free_parameters(model, mass_var)
    chi2ndof = frame.chiSquare(nDOF)
    print(f"chiSquare: {chi2ndof:0.2f} nDOF: {nDOF}")
    frame.Draw()
    for i, entry in enumerate(latex_list):
        latex.DrawLatexNDC(0.65, 0.85 - i * 0.05, entry) 
    
    print_canvas(f"{filename_prefix}_failed", f"{output_path}/{subdirectory}")
    results["failed_probe_fit"] = {
        "Nsig": ws.var("Nsig").getVal(),
        "Nsig_err": ws.var("Nsig").getError(),
        "chiSquare": chi2ndof
    }
    ws.var("jpsi_lambda").setConstant(False)
    ws.var("jpsi_gamma").setConstant(False)
    ws.var("jpsi_delta").setConstant(False)
    ws.var("psi2S_frac").setConstant(False)

    return results

def get_number_of_free_parameters(pdf, mass_var):

    # mass_var is observable
    params = pdf.getParameters(ROOT.RooArgSet(mass_var))

    # Loop and check which are floating
    free_params = []
    for i in range(params.getSize()):
        p = params[i]
        if not p.isConstant():
            free_params.append(p.GetName())

    return len(free_params)

def compute_ratio(data_result, mc_result):
    if data_result["all_probe_fit"]["Nsig"] > 0 and mc_result["all_probe_fit"]["Nsig"] > 0 :
        failed_data = data_result["failed_probe_fit"]["Nsig"]
        failed_data_err = data_result["failed_probe_fit"]["Nsig_err"]
        total_data = data_result["all_probe_fit"]["Nsig"]

        failed_mc = mc_result["failed_probe_fit"]["Nsig"]
        failed_mc_err = mc_result["failed_probe_fit"]["Nsig_err"]
        total_mc = mc_result["all_probe_fit"]["Nsig"]

        # Efficiencies
        eff_data = 1.0 - failed_data / total_data
        eff_mc = 1.0 - failed_mc / total_mc

        # Correction
        correction = eff_data / eff_mc

        # Propagated uncertainties (only failed_* carry uncertainty)
        eff_data_err = failed_data_err / total_data
        eff_mc_err = failed_mc_err / total_mc

        # Uncertainty on correction
        correction_err = correction * ((eff_data_err / eff_data) ** 2 + (eff_mc_err / eff_mc) ** 2) ** 0.5

        return correction, correction_err
    else:
        return None, None

def aggregate_corrections(name, hist_name, hist_title, projection_type, eras):
    if projection_type not in correction_results:
        return
    if projection_type == "pt":
        projection_bins = pt_bins
    elif projection_type == "eta":
        projection_bins = eta_bins
    else:
        raise Exception(f"Uknown projection: {projection_type}")
    
    h_eff_corr = ROOT.TH1F(hist_name, hist_title, len(projection_bins) - 1, array("d", projection_bins))
    h_eff_corr.SetMinimum(0.9)
    h_eff_corr.SetMaximum(1.05)
    h_eff_corr.SetLineWidth(2)
    h_eff_corr.SetMarkerStyle(20)
    h_eff_corr.SetMarkerSize(1.5)
    for i in range(1, len(projection_bins)):
        sum_weight_val = 0
        sum_weights = 0
        for era in eras:
            if era not in correction_results[projection_type]:
                print(f"{era} is not in correction_results for type: {projection_type}. Abort aggregation.")
                return
            correction, correction_err =  correction_results[projection_type][era][i - 1]
            weight = 1/(correction_err**2)
            sum_weight_val += weight * correction
            sum_weights += weight
            
        h_eff_corr.SetBinContent(i, sum_weight_val / sum_weights)
        h_eff_corr.SetBinError(i, 1 / math.sqrt(sum_weights))

    h_eff_corr.Draw("E1")
    c.SetGridy()
    c.Update()
    print_canvas(f"correction_vs_{projection_type}_{name}", output_path)
    c.SetGridy(0)

        
def corrections_with_bands(name, hist_name, x_title, projection_type, eras):
    if projection_type not in correction_results:
        return
    if projection_type == "pt":
        bins = pt_bins
    elif projection_type == "eta":
        bins = eta_bins
    else:
        raise Exception(f"Uknown projection: {projection_type}")

    x = array("d")
    y = array("d")
    exl = array("d")
    exh = array("d")
    eyl = array("d")
    eyh = array("d")
    nbins = len(bins) - 1
    for ibin in range(nbins):
        central_vals = []

        xlow  = bins[ibin]
        xhigh = bins[ibin + 1]
        xcenter = 0.5 * (xlow + xhigh)
        bin_width = xhigh - xlow

        for era in eras:
            val, err =  correction_results[projection_type][era][ibin]
            central_vals.append(val)
        
        vmin = min(central_vals)
        vmax = max(central_vals)
        ycenter = 0.5 * (vmin + vmax)

        x.append(xcenter)
        y.append(ycenter)
        exl.append(0.5 * bin_width)
        exh.append(0.5 * bin_width)
        eyl.append(ycenter - vmin)
        eyh.append(vmax - ycenter)

    band = ROOT.TGraphAsymmErrors(nbins, x, y, exl, exh, eyl, eyh)
    band.SetFillColorAlpha(ROOT.kGray, 1.0)
    band.SetLineColor(ROOT.kGray)
    band.SetLineWidth(1)
    band.SetMarkerStyle(0)
    band.SetFillStyle(1001)

    band.Draw("A2")  # A = axes, 2 = filled band

    c.SetGridy()
    
    # set y limits
    c.Update() # make sure the frame histogram exists
    hframe = band.GetHistogram()
    hframe.GetXaxis().SetTitle(x_title)
    hframe.GetYaxis().SetTitle("#varepsilon_{data}/#varepsilon_{MC}")
    hframe.GetYaxis().SetRangeUser(0.9, 1.07)
    c.Modified()

    # legend
    leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    leg.AddEntry(band, "Envelope (min to max)", "f")
    
    # add results
    hist_cache = []
    for i, era in enumerate(eras):
        x = array("d")
        y = array("d")
        ex = array("d")
        ey = array("d")
        for ibin, (v, e) in enumerate(correction_results[projection_type][era]):
            xlow = bins[ibin]
            xhigh = bins[ibin + 1]
            xc = 0.5 * (xlow + xhigh)
            x.append(xc)
            y.append(v)
            ex.append(0.0)        # no horizontal error for the points
            ey.append(e)          

        gr = ROOT.TGraphErrors(nbins, x, y, ex, ey)
        gr.SetTitle(era)
        gr.SetLineWidth(2)
        gr.SetMarkerStyle(20 + i)
        gr.SetMarkerSize(1.5)
        gr.Draw("PE same")  # points + error bars + line
        hist_cache.append(gr)

        leg.AddEntry(gr, era, "pl")
        
    leg.Draw()
    
    c.Update()
    print_canvas(f"correction_vs_{projection_type}_{name}", output_path)
    c.SetGridy(0)
    
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

    c = ROOT.TCanvas("c", "", 1600, 1200)
    c2 = ROOT.TCanvas("c2", "", 1600, 800)
    c.cd()
    
    # histos = defaultdict(dict)
    
    #######################################################
    ## Extract information for Jpsi efficiency measurement
    ##
    ## Note: it's slow process - don't force recreation of
    ## histograms unless you change data or selection
    #######################################################
    
    processor = DataProcessor()
    processor.process_samples()
    # pprint(processor.samples)

    mass = ROOT.RooRealVar("m", "",  min_mass, max_mass)
    
    if force_data_fits or not os.path.exists(file_fit_results):
 
        print_level = 0

        ws = build_model(mass)
        fit_results = dict()

        # prepare aggregated histograms by era
        if aggregate_eras:
            for sample in processor.histos:
                agg_histos = dict()
                for era in processor.histos[sample]:
                    agg_era = None
                    if era in ["Run2022C", "Run2022D", "Run2022E", "Run2022F"]:
                        agg_era = "Run2022"
                    if era in ["Run2023C", "Run2023D"]:
                        agg_era = "Run2023"
                    if era in ["Run2024C", "Run2024D"]:
                        agg_era = "Run2024CD"
                    if agg_era != None:
                        if agg_era not in agg_histos:
                            agg_histos[agg_era] = dict()
                            agg_histos[agg_era]["probe1"] = \
                                processor.histos[sample][era]["probe1"].Clone(
                                    processor.histo_key(sample, agg_era, "probe1"))
                            agg_histos[agg_era]["probe1"].SetDirectory(0)
                            agg_histos[agg_era]["probe2"] = \
                                processor.histos[sample][era]["probe2"].Clone(
                                    processor.histo_key(sample, agg_era, "probe2"))
                            agg_histos[agg_era]["probe2"].SetDirectory(0)
                            agg_histos[agg_era]["probe1_sig"] = \
                                processor.histos[sample][era]["probe1_sig"].Clone(
                                    processor.histo_key(sample, agg_era, "probe1_sig"))
                            agg_histos[agg_era]["probe1_sig"].SetDirectory(0)
                            agg_histos[agg_era]["probe2_sig"] = \
                                processor.histos[sample][era]["probe2_sig"].Clone(
                                    processor.histo_key(sample, agg_era, "probe2_sig"))
                            agg_histos[agg_era]["probe2_sig"].SetDirectory(0)
                            agg_histos[agg_era]["probe1_failed"] = \
                                processor.histos[sample][era]["probe1_failed"].Clone(
                                    processor.histo_key(sample, agg_era, "probe1_failed"))
                            agg_histos[agg_era]["probe1_failed"].SetDirectory(0)
                            agg_histos[agg_era]["probe2_failed"] = \
                                processor.histos[sample][era]["probe2_failed"].Clone(
                                    processor.histo_key(sample, agg_era, "probe2_failed"))
                            agg_histos[agg_era]["probe2_failed"].SetDirectory(0)
                        else:
                            agg_histos[agg_era]["probe1"].Add(processor.histos[sample][era]["probe1"])
                            agg_histos[agg_era]["probe2"].Add(processor.histos[sample][era]["probe2"])
                            agg_histos[agg_era]["probe1_sig"].Add(processor.histos[sample][era]["probe1_sig"])
                            agg_histos[agg_era]["probe2_sig"].Add(processor.histos[sample][era]["probe2_sig"])
                            agg_histos[agg_era]["probe1_failed"].Add(processor.histos[sample][era]["probe1_failed"])
                            agg_histos[agg_era]["probe2_failed"].Add(processor.histos[sample][era]["probe2_failed"])
                for agg_era, histos in agg_histos.items():
                    if agg_era not in processor.histos[sample]:
                        processor.histos[sample][agg_era] = dict()
                    for name, hist in histos.items():
                        processor.histos[sample][agg_era][name] = hist

        for sample in processor.histos:
            # if sample != "BuToJpsiK": continue
            for era in processor.histos[sample]:
                if process_only and not re.search(process_only, era):
                    continue
                # if era not in ["Run2022D"]: continue
                if era in ["Run2022G"]: continue

                hist_tag_tag = processor.histos[sample][era]["tag_tag"]
                hist_tag_tag.SetLineWidth(2)
                hist_tag_tag.Draw()
                print_canvas(f"tag_tag_mass_{era}", output_path)
                
                # Integrated Fit
                hist3D_probe = processor.histos[sample][era]["probe1"].Clone("hist3D_probe")
                hist3D_probe.Add(processor.histos[sample][era]["probe2"])
                hist_probe = hist3D_probe.ProjectionX("hist_probe")

                hist3D_ref = processor.histos[sample][era]["probe1_sig"].Clone("hist3D_probe")
                hist3D_ref.Add(processor.histos[sample][era]["probe2_sig"])
                hist_ref = hist3D_ref.ProjectionX("hist_ref")

                hist3D_failed = processor.histos[sample][era]["probe1_failed"].Clone("hist3D_failed")
                hist3D_failed.Add(processor.histos[sample][era]["probe2_failed"])
                hist_failed = hist3D_failed.ProjectionX("hist_failed")

                if sample not in fit_results:
                    fit_results[sample] = dict()
                if era not in fit_results[sample]:
                    fit_results[sample][era] = dict()
                fix_shape_to_ref_all_probes = False
                fix_shape_to_ref_failed_probes = True
                # if re.search("Summer", era):
                #     fix_shape_to_ref_failed_probes = False
                fit_results[sample][era]["average"] = fit_data(ws, era, "fits",
                                                               hist_probe, hist_failed, hist_ref, hist_tag_tag,
                                                               fix_shape_to_ref_all_probes,
                                                               fix_shape_to_ref_failed_probes)
                # Eta Fit
                if do_eta_fits:
                    h_ineff_eta = ROOT.TH1F("h_ineff_eta", "", len(eta_bins) - 1, array("d", eta_bins))
                    h_ineff_eta.SetMinimum(0)
                    h_ineff_eta.SetLineWidth(2)
                    h_ineff_eta.SetMarkerStyle(20)
                    h_ineff_eta.SetMarkerSize(1.5)
                    for i in range(1, len(eta_bins)):
                        hist_probe = hist3D_probe.ProjectionX("hist_probe", 0, -1, i, i)
                        hist_ref = hist3D_ref.ProjectionX("hist_ref", 0, -1, i, i)
                        hist_failed = hist3D_failed.ProjectionX("hist_failed", 0, -1, i, i)

                        result = fit_data(ws, f"{era}_eta_{eta_bins[i-1]}_{eta_bins[i]}", era,
                                          hist_probe, hist_failed, hist_ref)
                        if "eta" not in fit_results[sample][era]:
                            fit_results[sample][era]["eta"] = list()
                        result['bin'] = f"[{eta_bins[i-1]},{eta_bins[i]}]"
                        fit_results[sample][era]["eta"].append(result)

                        if result["all_probe_fit"]["Nsig"] > 0:
                            inefficiency = result["failed_probe_fit"]["Nsig"] / result["all_probe_fit"]["Nsig"]
                            inefficiency_err = result["failed_probe_fit"]["Nsig_err"] / result["all_probe_fit"]["Nsig"]
                            h_ineff_eta.SetBinContent(i, inefficiency)
                            h_ineff_eta.SetBinError(i, inefficiency_err)

                    h_ineff_eta.Draw("E1")
                    c.SetGridy()
                    c.Update()
                    print_canvas(f"inefficiency_vs_eta_{era}", output_path)
                    c.SetGridy(0)

                # Pt Fit
                if do_pt_fits:
                    h_ineff_pt = ROOT.TH1F("h_ineff_pt", "", len(pt_bins) - 1, array("d", pt_bins))
                    h_ineff_pt.SetMinimum(0)
                    h_ineff_eta.SetLineWidth(2)
                    h_ineff_pt.SetMarkerStyle(20)
                    h_ineff_pt.SetMarkerSize(1.5)
                    for i in range(1, len(pt_bins)):
                        hist_probe = hist3D_probe.ProjectionX("hist_probe", i, i)
                        hist_ref = hist3D_ref.ProjectionX("hist_ref", i, i)
                        hist_failed = hist3D_failed.ProjectionX("hist_failed", i, i)

                        result = fit_data(ws, f"{era}_pt_{pt_bins[i-1]}_{pt_bins[i]}", era,
                                          hist_probe, hist_failed, hist_ref)
                        if "pt" not in fit_results[sample][era]:
                            fit_results[sample][era]["pt"] = list()
                        result['bin'] = f"[{pt_bins[i-1]},{pt_bins[i]}]"
                        fit_results[sample][era]["pt"].append(result)

                        if result["all_probe_fit"]["Nsig"] > 0:
                            inefficiency = result["failed_probe_fit"]["Nsig"] / result["all_probe_fit"]["Nsig"]
                            inefficiency_err = result["failed_probe_fit"]["Nsig_err"] / result["all_probe_fit"]["Nsig"]
                            h_ineff_pt.SetBinContent(i, inefficiency)
                            h_ineff_pt.SetBinError(i, inefficiency_err)

                    h_ineff_pt.Draw("E1")
                    c.SetGridy()
                    c.Update()
                    print_canvas(f"inefficiency_vs_pt_{era}", output_path)
                    c.SetGridy(0)

        # pprint(fit_results)
        with open(f"{output_path}/fit_results.json", "w") as f:
            json.dump(fit_results, f, indent=4)
    else:
        with open(f"{output_path}/fit_results.json") as f:
            fit_results = json.load(f)

    for sample in fit_results:
        n = len(fit_results[sample])
        h = ROOT.TH1D(f"h","", n, 0, n)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(1.5)
        # h.SetMarkerSize(0.5)
        h.SetMinimum(0)
        for i, (era, results) in enumerate(fit_results[sample].items()):
            # skip aggregated eras
            if era not in processor.samples[sample]:
                continue
            if "average" in results:
                result = results["average"]
                inefficiency = result["failed_probe_fit"]["Nsig"] / result["all_probe_fit"]["Nsig"]
                inefficiency_err = result["failed_probe_fit"]["Nsig_err"] / result["all_probe_fit"]["Nsig"]
                h.SetBinContent(i + 1, inefficiency)
                h.SetBinError(i + 1, inefficiency_err)
                h.GetXaxis().SetBinLabel(i + 1, era)

        h.Draw()
        print_canvas(f"inefficiency_vs_era_{sample}", output_path)

    if compute_corrections:
        # Compute corrections
        c2.cd()
        c2.SetGridy()
        mc_sample = "BuToJpsiK"
        correction_results = dict()
        for sample, info in processor.samples.items():
            print(f"Processing {sample}")
            if not info["Data"]: continue
            if sample not in fit_results: continue

            for era in fit_results[sample]:
                if era not in processor.era2campaign: continue
                campaign = processor.era2campaign[era]

                if mc_sample not in fit_results: continue
                if campaign not in fit_results[mc_sample]: continue

                # eta
                if "eta" in fit_results[sample][era]:
                    h_eff_corr_eta = ROOT.TH1F("h_eff_corr_eta", ";#eta;#varepsilon_{data}/#varepsilon_{MC}",
                                               len(eta_bins) - 1, array("d", eta_bins))
                    h_eff_corr_eta.SetMinimum(0.9)
                    h_eff_corr_eta.SetMaximum(1.05)
                    h_eff_corr_eta.SetLineWidth(2)
                    h_eff_corr_eta.SetMarkerStyle(20)
                    h_eff_corr_eta.SetMarkerSize(1.5)
                    correction_result = []
                    for i in range(1, len(eta_bins)):
                        data_result = fit_results[sample][era]["eta"][i - 1]
                        mc_result = fit_results[mc_sample][campaign]["eta"][i - 1]

                        correction, correction_err = compute_ratio(data_result, mc_result)
                        correction_result.append((correction, correction_err))
                        if correction != None:
                            h_eff_corr_eta.SetBinContent(i, correction)
                            h_eff_corr_eta.SetBinError(i, correction_err)

                    h_eff_corr_eta.Draw("E1")
                    print_canvas(f"efficiency_correction_vs_eta_{era}", output_path)
                    if "eta" not in correction_results:
                        correction_results["eta"] = dict()
                    correction_results["eta"][era] = correction_result

                # pt
                if "pt" in fit_results[sample][era]:
                    h_eff_corr_pt = ROOT.TH1F("h_eff_corr_pt", ";p_{T}, GeV;#varepsilon_{data}/#varepsilon_{MC}",
                                              len(pt_bins) - 1, array("d", pt_bins))
                    h_eff_corr_pt.SetMinimum(0.9)
                    h_eff_corr_pt.SetMaximum(1.05)
                    h_eff_corr_pt.SetLineWidth(2)
                    h_eff_corr_pt.SetMarkerStyle(20)
                    h_eff_corr_pt.SetMarkerSize(1.5)
                    correction_result = []
                    for i in range(1, len(pt_bins)):
                        data_result = fit_results[sample][era]["pt"][i - 1]
                        mc_result = fit_results[mc_sample][campaign]["pt"][i - 1]

                        correction, correction_err = compute_ratio(data_result, mc_result)
                        correction_result.append((correction, correction_err))
                        if correction != None:
                            h_eff_corr_pt.SetBinContent(i, correction)
                            h_eff_corr_pt.SetBinError(i, correction_err)
                    if "pt" not in correction_results:
                        correction_results["pt"] = dict()
                    correction_results["pt"][era] = correction_result

                    h_eff_corr_pt.Draw("E1")
                    print_canvas(f"efficiency_correction_vs_pt_{era}", output_path)

                # average
                data_result = fit_results[sample][era]["average"]
                mc_result = fit_results[mc_sample][campaign]["average"]

                if "average" not in correction_results:
                    correction_results["average"] = dict()
                correction_results["average"][era] = compute_ratio(data_result, mc_result)

        # Aggregated corrections
        pprint(correction_results)
        with open(f"{output_path}/correction_results.json", "w") as f:
            json.dump(correction_results, f, indent=4)

        result = correction_results["average"]
        h_eff_corr = ROOT.TH1F("h_eff_corr", ";;#varepsilon_{data}/#varepsilon_{MC}", len(result), 0, len(result))
        h_eff_corr.SetMinimum(0.9)
        h_eff_corr.SetLineWidth(2)
        h_eff_corr.SetMaximum(1.05)
        h_eff_corr.SetMarkerSize(1.5)
        h_eff_corr.SetMarkerStyle(20)
        for i, (era, (correction, correction_err)) in enumerate(result.items()):
            h_eff_corr.SetBinContent(i + 1, correction)
            h_eff_corr.SetBinError(i + 1, correction_err)
            h_eff_corr.GetXaxis().SetBinLabel(i + 1, era)

        h_eff_corr.Draw("E1")
        print_canvas(f"correction_vs_era", output_path)

        # Rs plot (different study)
        result = {
            "Run2022C": (0.1108,0.0022),
            "Run2022D": (0.1120,0.0031),
            "Run2022E": (0.1085,0.0020),
            "Run2022F": (0.1073,0.0015),
            "Run2023C": (0.1095,0.0013),
            "Run2023D": (0.1103,0.0019),
            "Run2024C": (0.1096,0.0015),
            "Run2024D": (0.1088,0.0015),
            "Run2024E": (0.1074,0.0013),
            "Run2024F": (0.1072,0.0009),
            "Run2024G": (0.1076,0.0008),
            "Run2024H": (0.1064,0.0020),
            "Run2024I": (0.1046,0.0014),
        }
        h_eff_corr.SetMinimum(0.099)
        h_eff_corr.SetMaximum(0.116)
        h_eff_corr.GetYaxis().SetTitle("Rs")
        for i, (era, (correction, correction_err)) in enumerate(result.items()):
            h_eff_corr.SetBinContent(i + 1, correction)
            h_eff_corr.SetBinError(i + 1, correction_err)
            h_eff_corr.GetXaxis().SetBinLabel(i + 1, era)

        h_eff_corr.Draw("E1")
        print_canvas(f"Rs_vs_era", output_path)
        
        # Run2022
        aggregate_corrections("Run2022", "h_eff_corr_pt", ";p_{T}, GeV;#varepsilon_{data}/#varepsilon_{MC}",
                              "pt", ["Run2022C", "Run2022D", "Run2022E", "Run2022F"])
        aggregate_corrections("Run2022", "h_eff_corr_eta", ";#eta;#varepsilon_{data}/#varepsilon_{MC}", "eta",
                              ["Run2022C", "Run2022D", "Run2022E", "Run2022F"])
        
        aggregate_corrections("Run2023", "h_eff_corr_pt", ";p_{T}, GeV;#varepsilon_{data}/#varepsilon_{MC}",
                              "pt", ["Run2023C", "Run2023D"])
        aggregate_corrections("Run2023", "h_eff_corr_eta", ";#eta;#varepsilon_{data}/#varepsilon_{MC}", "eta",
                              ["Run2023C", "Run2023D"])
        
        aggregate_corrections("Run2024CDE", "h_eff_corr_pt", ";p_{T}, GeV;#varepsilon_{data}/#varepsilon_{MC}",
                              "pt", ["Run2024C", "Run2024D", "Run2024E"])
        aggregate_corrections("Run2024CDE", "h_eff_corr_eta", ";#eta;#varepsilon_{data}/#varepsilon_{MC}",
                              "eta", ["Run2024C", "Run2024D", "Run2024E"])
        
        # Correction Bands
        corrections_with_bands("Run2024", "h_eff_corr_eta", "#eta", "eta", processor.eras["Run2024"])
        corrections_with_bands("Run2025", "h_eff_corr_eta", "#eta", "eta", processor.eras["Run2025"])
        c.cd()
        
# WARNINGs

n_fits = 0
n_warnings = 0

def check_fit(message, result, fit_type, report_bin=False):
    bin = ""
    global n_fits
    global n_warnings
    if report_bin:
        bin = f"-{result['bin']}"
    n_fits += 1
    if result[fit_type]["chiSquare"] > 5:
        n_warnings += 1
        print(f"WARNING: {message}{bin} {fit_type} chi2ndof: {result[fit_type]['chiSquare']:0.2f}")
        
        
for sample, sample_info in fit_results.items():
    # if sample != "ParkingDoubleMuonLowMass": continue
    for era, era_info in sample_info.items():
        for result_type, result in era_info.items():
            if isinstance(result, list):
                for entry in result:
                    check_fit(f"{sample}-{era}-{result_type}", entry, "all_probe_fit", True)
                    check_fit(f"{sample}-{era}-{result_type}", entry, "failed_probe_fit", True)
            else:    
                check_fit(f"{sample}-{era}-{result_type}", result, "all_probe_fit")
                check_fit(f"{sample}-{era}-{result_type}", result, "failed_probe_fit")

print(f"n_fits: {n_fits}")
print(f"n_warnings: {n_warnings}")
