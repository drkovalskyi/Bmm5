import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob
from collections import defaultdict

ROOT.ROOT.EnableImplicitMT()

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/"
output_path = "/eos/home-d/dmytro/www/plots/ksmm/reco_input_ksmm_vs_kspipi/"

histos = defaultdict(dict)

def hh_histos(rdf, h_name, line_color=ROOT.kBlack, fill_color=None):
    rdf = rdf.Define("cands", "abs(hh_gen_pdgId)==310&&hh_had1_pt>4&&hh_had2_pt>3&&hh_kin_vtx_prob>0.001")
    rdf = rdf.Filter("Sum(cands)>0")

    histos["m"][h_name]  = rdf.Define("m", "hh_kin_mass[cands]").\
        Histo1D(("h_m",";m, [GeV]", 100, 0.35, 0.65), "m")
    histos["pt"][h_name] = rdf.Define("pt", "hh_kin_pt[cands]").\
        Histo1D(("h_pt",";p_{T}, [GeV]", 100, 5, 25), "pt")
    histos["eta"][h_name] = rdf.Define("eta", "hh_kin_eta[cands]").\
        Histo1D(("h_eta",";#eta", 100, -2.5, 2.5), "eta")
    histos["d1pt"][h_name] = rdf.Define("d1pt", "hh_had1_pt[cands]").\
        Histo1D(("h_m1pt",";p_{T}, [GeV]", 100, 0, 20), "d1pt")
    histos["d2pt"][h_name] = rdf.Define("d2pt", "hh_had2_pt[cands]").\
        Histo1D(("h_m2pt",";p_{T}, [GeV]", 100, 0, 20), "d2pt")

    h_list = []
    for h_type, h_dict in histos.items():
        h = h_dict[h_name]
        h.SetLineColor(line_color)
        h.SetLineWidth(2)
        if fill_color != None:
            h.SetFillColor(fill_color)
        h_list.append(h)
    ROOT.RDF.RunGraphs(h_list)
    
    n_events = rdf.Count().GetValue()
    return n_events


def mm_histos(rdf, h_name, line_color=ROOT.kBlack, fill_color=None):
    rdf = rdf.Define("cands", "abs(mm_gen_pdgId)==310&&mm_mu1_pt>4&&mm_mu2_pt>3&&mm_kin_vtx_prob>0.001")
    rdf = rdf.Filter("Sum(cands)>0")

    histos["m"][h_name]  = rdf.Define("m", "mm_kin_mass[cands]").\
        Histo1D(("h_m",";m, [GeV]", 100, 0.35, 0.65), "m")
    histos["pt"][h_name] = rdf.Define("pt", "mm_kin_pt[cands]").\
        Histo1D(("h_pt",";p_{T}, [GeV]", 100, 5, 25), "pt")
    histos["eta"][h_name] = rdf.Define("eta", "mm_kin_eta[cands]").\
        Histo1D(("h_eta",";#eta", 100, -2.5, 2.5), "eta")
    histos["d1pt"][h_name] = rdf.Define("d1pt", "mm_mu1_pt[cands]").\
        Histo1D(("h_m1pt",";p_{T}, [GeV]", 100, 0, 20), "d1pt")
    histos["d2pt"][h_name] = rdf.Define("d2pt", "mm_mu2_pt[cands]").\
        Histo1D(("h_m2pt",";p_{T}, [GeV]", 100, 0, 20), "d2pt")

    h_list = []
    for h_type, h_dict in histos.items():
        h = h_dict[h_name]
        h.SetLineColor(line_color)
        h.SetLineWidth(2)
        if fill_color != None:
            h.SetFillColor(fill_color)
        h_list.append(h)
    ROOT.RDF.RunGraphs(h_list)
    
    n_events = rdf.Count().GetValue()
    return n_events

    
def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))


def make_plots(prefix):
    for h_name, h_dict in histos.items():
        legend = ROOT.TLegend(0.60, 0.75, 0.87, 0.87)
        legend.SetShadowColor(ROOT.kWhite)
        legend.SetLineColor(ROOT.kWhite)
        max_val = 0

        # loop over histrograms of the given type
        for name, h in h_dict.items():
            n_bins = h.GetNbinsX()

            # Normalize to unity
            h.Scale(1./h.Integral(0, n_bins + 1))

            # Include overflow
            h.GetXaxis().SetRange(0, n_bins + 1)
        
            # Set maximum
            if max_val < h.GetMaximum():
                max_val = h.GetMaximum()

            legend.AddEntry(h.GetPtr(), name)
            
        for i, (name, h) in enumerate(h_dict.items()):
            h.SetMaximum(max_val * 1.25)
            if i == 0:
                h.Draw("hist")
            else:
                h.Draw("same hist")
                
        legend.Draw()
        print_canvas(f"{prefix}_{h_name}", output_path)
    

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

chain1 = ROOT.TChain("Events")
chain1.Add(f"{input_path}/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*root")
rdf1 = ROOT.RDataFrame(chain1)
n1 = mm_histos(rdf1, "K_{S}#rightarrow#mu#mu", ROOT.kBlue, ROOT.kMagenta)
print(f"n1: {n1}")

chain2 = ROOT.TChain("Events")
chain2.Add(f"{input_path}/K0sToPiPi_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*root")
rdf2 = ROOT.RDataFrame(chain2)
n2 = hh_histos(rdf2, "K_{S}#rightarrow#pi#pi")
print(f"n2: {n2}")

make_plots("sig")
