#!/bin/env python
import os, re, ROOT, sys, time
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/validation-bmm5-roc/"

# preselection = "mm_gen_mass>0&&abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4&&Muon_pt[mm_mu1_index]>4 && Muon_pt[mm_mu2_index]>4&&mm_kin_valid&&Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45&&mm_kin_sl3d>3"
preselection = "abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4&&Muon_pt[mm_mu1_index]>4 && Muon_pt[mm_mu2_index]>4&&mm_kin_valid&&Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45&&mm_kin_sl3d>4&&mm_kin_mass<30&&abs(mm_kin_pvip/mm_kin_pvipErr)<5&&mm_kin_vtx_chi2dof<5"

samples = [
    {
        'name':'BsToMuMu_BMuonFilter',
        'type':'sig',
        'path': '/afs/cern.ch/work/d/dmytro/projects/NanoAOD/src/BsToMuMu_BMuonFilter.root',
        'extra_cuts':''
    },
    {
        'name':'QCD_Pt-50to80_MuEnrichedPt5',
        'type':'bkg',
        'path':'/eos/cms/store/user/dmytro/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_NanoAOD_QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8_Bmm5-01/191014_051624/0000/*.root',
        'extra_cuts':''
        
    },
    {
        'name':'QCD_Pt-50to80_EMEnriched',
        'type':'bkg',
        'path':'/afs/cern.ch/work/d/dmytro/projects/NanoAOD/src/QCD_Pt-50to80_EMEnriched.root',
        'extra_cuts':''
    },
    # {
    #     'name':'LambdaBToPMuNu',
    #     'type':'bkg',
    #     'path':'/afs/cern.ch/work/d/dmytro/projects/NanoAOD/src/LambdaBToPMuNu.root',
    #     'extra_cuts':''
    # }
    {
        'name':'QCD_Pt-15to20_MuEnrichedPt5',
        'type':'bkg',
        'path':'/eos/cms/store/user/dmytro/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_NanoAOD_QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8_Bmm5-01/191014_050926/0000/*.root',
        'extra_cuts':''
        
    },
    {
        'name':'QCD_Pt-120to170_MuEnrichedPt5',
        'type':'bkg',
        'path':'/eos/cms/store/user/dmytro/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_NanoAOD_QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8_Bmm5-01/191014_050337/0000/*.root',
        'extra_cuts':''
        
    },
]

# colors = ['#17DFDF','#8800ff','#0000aa','#00ff00','#ffff00','#ff0000']
# ROOT.TColor.GetColor
colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kOrange+5, ROOT.kGreen+3]

def get_nevents_processed(path_to_root_files):
    chain = ROOT.TChain("Runs")
    chain.Add(path_to_root_files)
    nEvents = 0
    for run in chain:
        nEvents += run.genEventCount
    return nEvents

def process_sample(sample):
    # Run info
    chain_Runs = ROOT.TChain("Runs")
    chain_Runs.Add(sample['path'])
    nEvents = 0
    for run in chain_Runs:
        nEvents += run.genEventCount
    sample['nEvents'] = nEvents
    print nEvents

    # Event info
    chain_Events = ROOT.TChain("Events")
    chain_Events.Add(sample['path'])
    h_bdt = ROOT.TH1D("h_bdt",sample['name'],30,-1.5,1.5)
    chain_Events.Draw("mm_bdt>>h_bdt")
    sample['h_bdt'] = h_bdt
    h_bdt.SetDirectory(0)
    print h_bdt.GetEntries()
    print h_bdt.Integral()

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    # canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    # canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def plot_results(hist_title,file_name,sig_sample,bkg_samples,relative_efficiency,hist_name='h_bdt'):
    sig_eff = array( "f" )
    bkg_eff = []
    n_bkg_samples = len(bkg_samples)
    for i in range(n_bkg_samples):
        bkg_eff.append(array( "f" ))
    nbins = sig_sample[hist_name].GetNbinsX()

    for bin in range(nbins):
        if relative_efficiency:
            sig_eff.append(sig_sample[hist_name].Integral(bin,nbins)/sig_sample[hist_name].Integral(0,nbins))
        else:
            sig_eff.append(sig_sample[hist_name].Integral(bin,nbins)/sig_sample['nEvents'])
        for i in range(n_bkg_samples):
            if relative_efficiency:
                bkg_eff[i].append(bkg_samples[i][hist_name].Integral(bin,nbins)/bkg_samples[i][hist_name].Integral(0,nbins))
            else:
                bkg_eff[i].append(bkg_samples[i][hist_name].Integral(bin,nbins)/bkg_samples[i]['nEvents'])
                
    max_bkg_eff = 0.0
    for i in range(n_bkg_samples):
        if max(bkg_eff[i])>max_bkg_eff:
            max_bkg_eff = max(bkg_eff[i])
            
    print max_bkg_eff

    first_graph = True

    legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
    legend.SetFillStyle(0)
    legend.SetLineWidth(0)
    for i in range(n_bkg_samples):
        graph = ROOT.TGraph(len(sig_eff), sig_eff, bkg_eff[i])
        bkg_samples[i]['graph'] = graph
        icolor = i
        if icolor >= len(colors): icolor = 0
        graph.SetLineColor(colors[icolor])
        graph.SetLineWidth(3)
        if first_graph:
            graph.SetTitle(hist_title)
            graph.SetMinimum(0.00)
            graph.SetMaximum(max_bkg_eff)
            if relative_efficiency:
                graph.GetXaxis().SetTitle("Relative Signal Efficiency")
                graph.GetYaxis().SetTitle("Relative Background Efficiency")
            else:
                graph.GetXaxis().SetTitle("Absolute Signal Efficiency")
                graph.GetYaxis().SetTitle("Absolute Background Efficiency")
            graph.Draw("AC")
            first_graph = False
        else:
            graph.Draw("C same")
            
        legend.AddEntry(graph,bkg_samples[i]['name'])
    legend.Draw()

    print_canvas(file_name, output_path)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
# ROOT.TH1.AddDirectory(False)
c1 = TCanvas("c1","c1",800,800)

# for study in studies:

#     fBmm4_sig = TFile.Open(study['bmm4_sig_file_name'])
#     assert(fBmm4_sig)
#     bmm4_sig_events = fBmm4_sig.Get("candAnaMuMu/events")
#     assert(bmm4_sig_events)

#     fBmm4_bkg = TFile.Open(study['bmm4_bkg_file_name'])
#     assert(fBmm4_bkg)
#     bmm4_bkg_events = fBmm4_bkg.Get("candAnaMuMu/events")
#     assert(bmm4_bkg_events)

#     fBmm5_sig = TFile.Open(study['bmm5_sig_file_name'])
#     assert(fBmm5_sig)
#     bmm5_sig_events = fBmm5_sig.Get("Events")
#     assert(bmm5_sig_events)

#     fBmm5_bkg = TFile.Open(study['bmm5_bkg_file_name'])
#     assert(fBmm5_bkg)
#     bmm5_bkg_events = fBmm5_bkg.Get("Events")
#     assert(bmm5_bkg_events)

#     plot_data(study['name'], "%s_absolute"%study['hist_name'],"bdt","mm_bdt",30,-1.5,1.5,False)
#     plot_data(study['name'], "%s_relative"%study['hist_name'],"bdt","mm_bdt",30,-1.5,1.5,True)

# print get_nevents_processed("/afs/cern.ch/work/d/dmytro/projects/NanoAOD/src/BsToMuMu_BMuonFilter.root")
# print get_nevents_processed("/eos/cms/store/user/dmytro/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_NanoAOD_QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8_Bmm5-01/191014_051624/0000/*.root")

for sample in samples:
    process_sample(sample)
    print sample

sig_samples = []
bkg_samples = []

for sample in samples:
    if sample['type']=='sig': sig_samples.append(sample)
    if sample['type']=='bkg': bkg_samples.append(sample)

for sig_sample in sig_samples:
    name = sig_sample['name']
    plot_results(name,"%s_relative"%name,sig_sample,bkg_samples,True)
    plot_results(name,"%s_absolute"%name,sig_sample,bkg_samples,False)
        
