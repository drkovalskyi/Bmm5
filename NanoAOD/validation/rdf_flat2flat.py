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

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/ksmm_sb"

def add_files(chain, path, pattern = "^[^\/]+.root$"):
    nfiles = 0
    for root_dir, _, files in os.walk(path):
        for file in files:
            if re.search(f"{pattern}", file):
                full_path = os.path.join(root_dir, file)
                # print(full_path)
                nfiles += 1
                chain.Add(full_path)
    print(f"Number of files: {nfiles}")

    
####################################################################################

ROOT.gROOT.SetBatch(True)

treename = "ksmmData"
chain = ROOT.TChain(treename)
add_files(chain, f"{input_path}/")

rdf = ROOT.RDataFrame(chain)
rdf.Snapshot(treename,'/tmp/dmytro/test.root','^m$')
