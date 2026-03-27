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

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/533/fit-bkkmm2"
output_path = "/eos/home-d/dmytro/www/plots/tmp/mass-BsJpsi/"
path_pattern = "ParkingDoubleMuonLowMass0.*Run2024"

min_mass = 4.9
max_mass = 5.9
nbins = 100

histos = defaultdict(dict)


def book_histos(rdf):
    histos["mass"] = rdf.Histo1D(("h_mass",";Mass, GeV", \
                                  nbins, min_mass, max_mass), "m")
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


def make_plots(study_name = ""):
    for name, h in histos.items():
        file_name = name
        if study_name != "":
            file_name += "_" + study_name
        h.Draw("hist")
        print_canvas(f"{study_name}", output_path)


def add_files(chain, path, file_pattern = "^[^\/]+.root$"):
    nfiles = 0
    for root_dir, _, files in os.walk(path):
        if path_pattern != None:
            if not re.search(f"{path_pattern}", root_dir):
                continue
        print(f"root_dir: {root_dir}, n:{len(files)}")
        for file in files:
            if re.search(f"{file_pattern}", file):
                full_path = os.path.join(root_dir, file)
                # print(full_path)
                nfiles += 1
                chain.Add(full_path)
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

# chain = ROOT.TChain("ksmmData")
chain = ROOT.TChain("bspsiphiData")
add_files(chain, f"{input_path}/", )

rdf = ROOT.RDataFrame(chain)
book_histos(rdf)

make_plots("2024-lowpt")
