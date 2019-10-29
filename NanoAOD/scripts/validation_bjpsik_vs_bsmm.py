#!/bin/env python
import os, re, ROOT, sys, time
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_dev/bjpsik_vs_bm/"

mm_selection = "mm_kin_sl3d>4&&mm_kin_vtx_prob>0.1"
mm_selection += "&&mm_kin_mass<6.0"
mm_selection += "&&abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4&&Muon_pt[mm_mu1_index]>4 && Muon_pt[mm_mu2_index]>4&&mm_kin_valid>0&&Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45"

bkmm_selection = "bkmm_kaon_sdxy_bs>5&&bkmm_kaon_pt>1&&abs(bkmm_kaon_eta)<1.4" 
bkmm_selection += "&&bkmm_jpsimc_vtx_prob>0.1&&bkmm_jpsimc_sl3d>4"
bkmm_selection += "&&bkmm_jpsimc_mass<6.0"
bkmm_selection += "&&abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4&&Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4&&Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45"

samples = {
    'BdToMuMu_RunIIAutumn18MiniAOD_102X':{
        'files':[
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/04F83F9C-6B10-034F-992D-DC47B06E39F3.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/07973B48-B70B-A846-A7B2-67A2A9980511.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/08CDAA0D-57AD-8C42-8F70-56D2BC06E1C5.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0C6C2C09-D676-204D-AEF9-F84A1C3FD48D.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0F7FC5E6-350F-DA42-9D6D-EF5EFEB06282.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/17ADAFA2-4BBB-9442-8A6E-54A1677B9A1A.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/1C8BBD81-5198-CF4A-80A7-1F44CE6EB7FE.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/1E6EE082-2A94-5B46-822F-8FD0B927AA3B.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/2574F45A-8254-6C4B-A07E-977A99C38776.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502//BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/2E3D44E4-EF97-154C-86E4-3622D06DF682.root',
        ],
        'color':ROOT.kBlue,
        'type':'bmm',
    },
    'BuToJpsiK_RunIIAutumn18MiniAOD_102X':{
        'files':[
            '/afs/cern.ch/work/d/dmytro/projects/NanoAOD-CentOS7/src/BuToJpsiK_BMuonFilter_RunIIAutumn18MiniAOD.root'
        ],
        'color':ROOT.kBlack,
        'type':'bjpsik'
    }
}

# read list of files from a file instead

# with open("Charmonium+Run2018D-PromptReco-v2+MINIAOD.list") as f:
#     samples['Charmonium_Run2018D_PromptReco']['files'] = []
#     for file in f:
#         samples['Charmonium_Run2018D_PromptReco']['files'].append(file.strip())
        
# Load data

mask = '\w\d\_\.'
for name,sample in samples.items():
    if re.search("[^%s]"%mask,name):
        raise Exception("Illigal symbol used for sample name. Allowed %s" % mask) 
    chain = ROOT.TChain("Events")
    for entry in sample['files']:
        chain.Add(entry)
    sample['events'] = chain
    
print "Number of events:"
for name,sample in samples.items():
    sample['nAll'] = sample['events'].GetEntries()
    print "\t%s: \t%u" % (name,sample['nAll'])

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    # canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    # canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def plot_generic_1D(selections,hist_title,file_name,vars,nbins=100,xmin=0,xmax=100):
    c1 = TCanvas("c1","c1",800,800)
    max_value = 0
    for name,sample in samples.items():
        selection = selections[sample['type']]
        var = vars[sample['type']]
        hist = ROOT.TH1D("hist",hist_title,nbins,xmin,xmax)
        hist.SetLineColor(sample['color'])
        hist.SetLineWidth(2)
        sample['nSelected'] = sample['events'].Draw("%s>>hist"%var,selection)
        print_canvas("%s_%s"%(file_name,name), output_path)
        if hist.GetEntries()>0:
            hist.Scale(1/hist.GetEntries()) # normalize
        if max_value < hist.GetMaximum(): max_value = hist.GetMaximum()
        hist.SetDirectory(0)
        sample['hist'] = hist
    
    legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
    legend.SetFillStyle(0)
    legend.SetLineWidth(0)

    first_plot = True
    for name,sample in samples.items():
        sample['hist'].SetMinimum(0)
        sample['hist'].SetMaximum(max_value*1.2)
        legend.AddEntry(sample['hist'],name)
        if first_plot:
            sample['hist'].Draw("hist")
            first_plot = False
        else:
            sample['hist'].Draw("hist same")
    legend.Draw()
    print_canvas(file_name, output_path)
    print "Number of selected events:"
    for name,sample in samples.items():
        print "\t%s: \t%u out of %u" % (name,sample['nSelected'],sample['nAll'])

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

selections = {
    'bmm':mm_selection,
    'bjpsik':bkmm_selection
}

plot_generic_1D(selections,"BDT Matched", "09_bdt_matched",{'bmm':'mm_bdt','bjpsik':'bkmm_bmm_bdt'},100,-1.5,1.5)

plot_generic_1D(selections,"BDT Raw", "09_bdt_raw",{'bmm':'mm_bdt','bjpsik':'mm_bdt[bkmm_mm_index]'},100,-1.5,1.5)
