import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob

ROOT.ROOT.EnableImplicitMT()

input_path = '/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/'
output_path = "/eos/home-d/dmytro/www/plots/tmp/ratio_bjpsik_vs_bsjpsiphi"
# xrootd_server = "eoscms.cern.ch:1094"
xrootd_server = None
histos = dict()
chains = []
kaon_pt_values = [1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 3.0]
muon_selection = "m1pt > 4 && m2pt > 3"
muon_selection_name = "43"
# muon_selection_name = "22_m1pt.lt.4.or.m2pt.lt.3"
nbins = 22
min_pt = 8
max_pt = 30
# nfiles = 1
nfiles = 999999
trigger = "HLT_DoubleMu4_3_LowMass"
trigger_name = trigger
# trigger = "HLT_Dimuon0_Jpsi"
# trigger = "HLT_DoubleMu2_Jpsi_LowPt"

def xrd_ls(path):
    if xrootd_server:
        command = f"xrdfs {xrootd_server} ls {path}"
    else:
        command = f"find {path} -type f -name \"*.root\""
    output = subprocess.check_output(command, text=True, shell=True)
    files = []
    for line in output.splitlines():
        if re.search('\.root$', line):
            files.append(line)
    return files[:nfiles]


def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))


def get_data(datasets, tree_name="Events"):
    chain = ROOT.TChain(tree_name)
    chains.append(chain)
    n_files = 0
    for dataset in datasets:
        for file in xrd_ls(f"{input_path}/{dataset}"):
            if xrootd_server:
                n_files += chain.Add(f"root://{xrootd_server}/{file}")
            else:
                n_files += chain.Add(f"{file}")
    print(f"Number of files: {n_files}")
    # print(f"Number of events: {chain.GetEntries()}")
    return ROOT.RDataFrame(chain)


def process_BsToJPsiPhi_data():
    datasets = [
        "fit-bkkmm/ParkingDoubleMuonLowMass4+Run2024C-PromptReco-v1+MINIAOD/",
        
    ]
    rdf = get_data(datasets, "bspsiphiData")
    rdf = rdf.Filter(muon_selection)
    h_reco = rdf.Histo1D(("h_reco_BsToJPsiPhi", "", nbins, min_pt, max_pt), "pt")
    histos['h_reco_BsToJPsiPhi'] = h_reco
    
    h_trig = rdf.Filter(trigger).Histo1D(("h_trig_BsToJPsiPhi", "", nbins, min_pt, max_pt), "pt")
    h_trig.Sumw2()
    histos['h_trig_BsToJPsiPhi'] = h_trig
    ROOT.RDF.RunGraphs([h_reco, h_trig])
    n_events = rdf.Count().GetValue()
    return n_events

def process_BuToJpsiK_data():
    datasets = [
        'fit-bkmm/ParkingDoubleMuonLowMass4+Run2024C-PromptReco-v1+MINIAOD/',
    ]
    rdf = get_data(datasets, "bupsikData")

    rdfs = []
    histos_to_process = []
    for i, pt in enumerate(kaon_pt_values):
        rdfs.append(rdf.Filter(f"kaon_pt>{pt} && {muon_selection}"))
        h_reco = rdfs[i].Histo1D((f"h_reco_BuToJpsiK_{pt}", "", nbins, min_pt, max_pt), "pt")
        histos[f"h_reco_BuToJpsiK_{pt}"] = h_reco
        histos_to_process.append(h_reco)
        
        h_trig = rdfs[i].Filter(trigger).Histo1D((f"h_trig_BuToJpsiK_{pt}", "", nbins, min_pt, max_pt), "pt")
        h_trig.Sumw2()
        histos[f"h_trig_BuToJpsiK_{pt}"] = h_trig
        histos_to_process.append(h_trig)
        
    ROOT.RDF.RunGraphs(histos_to_process)
    n_events = rdf.Count().GetValue()
    return n_events


def plot_ratio(name, h1, h2, r_min=0.5, r_max=1.5):
    hist = h1.Clone(name)
    hist.Divide(h1, h2)
    hist.SetMinimum(r_min)
    hist.SetMaximum(r_max)
    hist.SetMarkerColor(ROOT.kBlue)
    hist.SetMarkerStyle(20)
    hist.GetXaxis().SetTitle("p^{B}_{T}, [GeV]")
    # hist.GetXaxis().SetRangeUser(10, 30)
    hist.Draw("e0")
    print_canvas(name, output_path)

n_BsToJPsiPhi = process_BsToJPsiPhi_data()
print(f"Number of BsToJPsiPhi events: {n_BsToJPsiPhi}")
n_BuToJpsiK = process_BuToJpsiK_data()
print(f"Number of BuToJpsiK events: {n_BuToJpsiK}")

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

c1 = TCanvas("c1", "c1", 600, 600)
ROOT.gPad.SetGrid()

histos[f'h_reco_BsToJPsiPhi'].Draw()
# print_canvas(f'h{muon_selection_name}_reco_BsToJPsiPhi', output_path)
histos[f'h_trig_BsToJPsiPhi'].Draw()
print_canvas(f'h{muon_selection_name}_trig_{trigger_name}_BsToJPsiPhi', output_path)

for i, pt in enumerate(kaon_pt_values):
    histos[f'h_reco_BuToJpsiK_{pt}'].Draw()
    # print_canvas(f'h{muon_selection_name}_reco_BuToJpsiK_{pt}', output_path)
    histos[f'h_trig_BuToJpsiK_{pt}'].Draw()
    print_canvas(f'h{muon_selection_name}_trig_{trigger_name}_BuToJpsiK_{pt}', output_path)
    
    plot_ratio(f'h{muon_selection_name}_trig_{trigger_name}_ratio_k{pt}',
               histos['h_trig_BsToJPsiPhi'].GetPtr(), histos[f'h_trig_BuToJpsiK_{pt}'].GetPtr(),
               0, 1)
