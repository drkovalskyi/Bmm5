import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob
from collections import defaultdict

ROOT.ROOT.EnableImplicitMT()

input_path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/PostProcessing/FlatNtuples/532/"
output_path = "/eos/home-d/dmytro/www/plots/ksmm/MVA_input_ksmm_vs_kspipimm/"

histos = defaultdict(dict)

def book_histos(rdf, h_name, selection = "", line_color=ROOT.kBlack, fill_color=None):
    rdf = rdf.Filter("displacement_lxy>1")
    if selection != "":
        rdf = rdf.Filter(selection)

    histos["m"][h_name]  = rdf.Histo1D(("h_m",";m, [GeV]", 100, 0.35, 0.65), "m")
    histos["pt"][h_name] = rdf.Histo1D(("h_pt",";p_{T}, [GeV]", 100, 5, 25), "pt")
    histos["eta"][h_name] = rdf.Histo1D(("h_eta",";#eta", 100, -2.5, 2.5), "eta")
    histos["d1pt"][h_name] = rdf.Histo1D(("h_m1pt",";p_{T}, [GeV]", 100, 0, 20), "d1pt")
    histos["d2pt"][h_name] = rdf.Histo1D(("h_m2pt",";p_{T}, [GeV]", 100, 0, 20), "d2pt")
    histos["d1nLostHits"][h_name] = rdf.Histo1D(("h_m1nLostHits",";nLostHitsInner", 5, 0, 5), "vtx_qual_d1_nLostHitsInner")
    histos["d2nLostHits"][h_name] = rdf.Histo1D(("h_m2nLostHits",";nLostHitsInner", 5, 0, 5), "vtx_qual_d2_nLostHitsInner")
    histos["d1nVeto"][h_name] = rdf.Histo1D(("h_m1Veto",";Veto", 5, 0, 5), "vtx_qual_d1_veto")
    histos["d2nVeto"][h_name] = rdf.Histo1D(("h_m2Veto",";Veto", 5, 0, 5), "vtx_qual_d2_veto")
    histos["vtx_prob"][h_name] = rdf.Histo1D(("h_vtx_prob",";Probability", 50, 0, 1), "vtx_qual_prob")
    histos["vtx_qual_doca"][h_name] = rdf.Histo1D(("h_vtx_qual_doca",";Distance of closest approach, [cm]", 100, 0, 0.02), "vtx_qual_doca")
    histos["displacement_l3d"][h_name] = rdf.Histo1D(("h_displacement_l3d",";Vertex displacement (3D), [cm]", 100, 0, 50), "displacement_l3d")
    histos["displacement_sl3d"][h_name] = rdf.Histo1D(("h_displacement_sl3d",";Vertex displacement significance (3D)", 100, 0, 1000), "displacement_sl3d")
    histos["displacement_lxy"][h_name] = rdf.Histo1D(("h_displacement_lxy",";Vertex displacement (transverse), [cm]", 100, 0, 50), "displacement_lxy")
    histos["displacement_slxy"][h_name] = rdf.Histo1D(("h_displacement_slxy",";Vertex displacement significance (transverse)", 100, 0, 1000), "displacement_slxy")
    
    histos["vtx_qual_flightdist_err"][h_name] = rdf.Histo1D(("h_vtx_qual_flightdist_err",";Flight length uncertainty", 100, 0, 1), "vtx_qual_flightdist_err")
    histos["vtx_qual_flightdist2D_err"][h_name] = rdf.Histo1D(("vtx_qual_flightdist2D_err",";Flight length uncertainty (2D)", 100, 0, 0.1), "vtx_qual_flightdist2D_err")

    histos["displacement_pseudo_decaytime"][h_name] = rdf.Histo1D(("h_displacement_pseudo_decaytime",";Decay time with mixed mass, [ps]", 100, 0, 100), "displacement_pseudo_decaytime")
    histos["iso_vtx_tight_ntracks"][h_name] = rdf.Histo1D(("h_iso_vtx_tight_ntracks",";Vertex isolation - ntracks (tight)", 100, 0, 10), "iso_vtx_tight_ntracks")
    histos["iso_vtx_loose_ntracks"][h_name] = rdf.Histo1D(("iso_vtx_loose_ntracks",";Vertex isolation - ntracks (loose)", 100, 0, 10), "iso_vtx_loose_ntracks")

# *Br   26 :iso_vtx_tight_sum_pt : iso_vtx_tight_sum_pt/F                      *
# *Br   28 :iso_vtx_loose_sum_pt : iso_vtx_loose_sum_pt/F                      *
# *............................................................................*
    histos["iso_tight_ntracks"][h_name] = rdf.Histo1D(("iso_tight_ntracks",";Isolation - ntracks (tight)", 100, 0, 10), "iso_tight_ntracks")
    histos["iso_loose_ntracks"][h_name] = rdf.Histo1D(("iso_loose_ntracks",";Isolation - ntracks (loose)", 100, 0, 10), "iso_loose_ntracks")

# *Br   30 :iso_mu_tight_sum_pt : iso_mu_tight_sum_pt/F                        *
# *Br   32 :iso_mu_loose_sum_pt : iso_mu_loose_sum_pt/F                        *

    histos["prod_alpha"][h_name] = rdf.Histo1D(("prod_alpha",";Pointing angle (3D)", 100, 0, 0.01), "prod_alpha")
    histos["prod_alpha_significance"][h_name] = rdf.Histo1D(("prod_alpha_significance",";Pointing angle significance (3D)", 100, 0, 10), "prod_alpha_significance")
    histos["prod_alphaBS"][h_name] = rdf.Histo1D(("prod_alphaBS",";Pointing angle (2D)", 100, 0, 0.01), "prod_alphaBS")
    histos["prod_alphaBS_significance"][h_name] = rdf.Histo1D(("prod_alphaBS_significance",";Pointing angle significance (2D)", 100, 0, 10), "prod_alphaBS_significance")
    histos["prod_pvip"][h_name] = rdf.Histo1D(("prod_pvip",";Impact parameter (3D)", 100, 0, 0.1), "prod_pvip")
    histos["prod_spvip"][h_name] = rdf.Histo1D(("prod_spvip",";Impact parameter significance (3D)", 100, 0, 10), "prod_spvip")

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
    # print(n_events)
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

chain1 = ROOT.TChain("ksmmMc")
chain1.Add(f"{input_path}/ksmm/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1+MINIAODSIM/*.root")
chain1.Add(f"{input_path}/ksmm/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v1+MINIAODSIM/*.root")
chain1.Add(f"{input_path}/ksmm/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4+MINIAODSIM/*.root")
chain1.Add(f"{input_path}/ksmm/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15_ext1-v2+MINIAODSIM/*.root")
chain1.Add(f"{input_path}/ksmm/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4+MINIAODSIM/*.root")
rdf1 = ROOT.RDataFrame(chain1)
n1 = book_histos(rdf1, "K_{S}#rightarrow#mu#mu", "", ROOT.kBlue, ROOT.kMagenta)
print(f"n1: {n1}")

chain2 = ROOT.TChain("kspipimmMc")
chain2.Add(f"{input_path}/ksmm/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2+MINIAODSIM/*.root")
chain2.Add(f"{input_path}/ksmm/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/*.root")
chain2.Add(f"{input_path}/ksmm/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2+MINIAODSIM/*.root")
chain2.Add(f"{input_path}/ksmm/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15_ext1-v2+MINIAODSIM/*.root")
chain2.Add(f"{input_path}/ksmm/K0sTo2PiTo2Mu_K0sFilter_MuFilter_PiLifetime0p02_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2+MINIAODSIM/*.root")
rdf2 = ROOT.RDataFrame(chain2)
n2 = book_histos(rdf2, "K_{S}#rightarrow#pi#pi#rightarrow2#mu2#nu", "m>0.45 and m<0.465")
print(f"n2: {n2}")

make_plots("peak")
