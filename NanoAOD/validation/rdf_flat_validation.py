import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob
from collections import defaultdict
import pprint
from math import sqrt
import numpy as np
import json

ROOT.ROOT.EnableImplicitMT()

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/"
# input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/533/fit-bkkmm2"
# input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/fit-bkkmm"
# input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/533/fit-bkmm2"
# input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/fit-bkmm"
# input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/fit-bkstarmm/"
# input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/533/fit-bkstarmm2/"
output_path = "/eos/home-d/dmytro/www/plots/tmp/flat-validation/"
# path_pattern = "ParkingDoubleMuonLowMass0.*Run2024"
path_pattern = "ParkingDoubleMuonLowMass.*Run2024"

min_mass = 4.9
max_mass = 5.9
nbins = 100

histos = defaultdict(dict)

histograms = {
    "mass": [";Mass, GeV", 100, 4.9, 5.9, "m"],
    "pt": [";p_{T}, GeV",   20,   0,  10, "pt"]
}
studies = {
    "bs2jpsiphi": {
        "samples": {
            "2024" : f"{input_path}/532/fit-bkkmm",
            "2024_lowpt" : f"{input_path}/533/fit-bkkmm2"
        },
        "tree": "bspsiphiData",
    },
    "bu2jpsik": {
        "samples": {
            "2024" : f"{input_path}/532/fit-bkmm",
            "2024_lowpt" : f"{input_path}/533/fit-bkmm2"
        },
        "tree": "bupsikData",
    },
    "bd2jpsikstar": {
        "samples": {
            "2024" : f"{input_path}/532/fit-bkstarmm",
            "2024_lowpt" : f"{input_path}/533/fit-bkstarmm2"
        },
        "tree": "bdpsikstarData",
    }
}

def fill_histograms(rdf, prefix):
    rdf = rdf.Filter("pt<10")
    for name, params in histograms.items():
        h_name = f"{prefix}_{name}"
        histos[h_name] = rdf.Histo1D((f"h_{h_name}", *params[:-1]), params[-1])
    h_list = []
    for name, h in histos.items():
        # h.SetLineColor(ROOT.kBlack)
        # h.SetLineWidth(2)
        h.SetMinimum(0)
        h_list.append(h)
    ROOT.RDF.RunGraphs(h_list)
    n_events = rdf.Count().GetValue()
    print(n_events)
    return n_events


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print(f"{path}/{output_name_without_extention}.png")
    canvas.Print(f"{path}/{output_name_without_extention}.pdf")
    canvas.Print(f"{path}/{output_name_without_extention}.root")


def make_plots():
    for name, h in histos.items():
        h.Draw("hist")
        print_canvas(f"{name}", output_path)


def add_files(chain, path, file_pattern = "^[^\/]+.root$"):
    nfiles = 0
    for root_dir, _, files in os.walk(path):
        if path_pattern != None:
            if not re.search(f"{path_pattern}", root_dir):
                continue
        n = 0
        for file in files:
            if re.search(f"{file_pattern}", file):
                full_path = os.path.join(root_dir, file)
                # print(full_path)
                n += 1
                nfiles += 1
                chain.Add(full_path)
        if n > 0:
            print(f"root_dir: {root_dir}, n:{n}")
    print(f"Number of files: {nfiles}")

    
####################################################################################

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
ROOT.gStyle.SetPadBottomMargin(0.15)

c = ROOT.TCanvas("c", "", 600, 600)

for study_name, info in studies.items():
    for sample_name, input_pattern in info["samples"].items():
        chain = ROOT.TChain(info["tree"])
        add_files(chain, input_pattern)
        rdf = ROOT.RDataFrame(chain)
        fill_histograms(rdf, f"{study_name}_{sample_name}")
make_plots()
