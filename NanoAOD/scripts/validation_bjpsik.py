#!/bin/env python
import os, re, ROOT, sys, time
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv5-1-502/bjpsik/"

mm_selection = "abs(Muon_eta[mm_mu1_index])<1.4 && abs(Muon_eta[mm_mu2_index])<1.4&&Muon_pt[mm_mu1_index]>4 && Muon_pt[mm_mu2_index]>4&&mm_kin_valid>0&&Muon_softMva[mm_mu1_index]>0.45 && Muon_softMva[mm_mu2_index]>0.45&&mm_kin_sl3d>4&&mm_kin_mass<6.0&&mm_kin_vtx_chi2dof<5"

bkmm_selection = "bkmm_kaon_sdxy_bs>5&&bkmm_kaon_pt>1&&abs(bkmm_kaon_eta)<1.4&&bkmm_jpsimc_vtx_prob>0.1&&bkmm_jpsimc_sl3d>4&&bkmm_jpsimc_mass<6.0&&abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4&&Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4&&Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45"

samples = {
    'BuToJpsiK_RunIIAutumn18MiniAOD_102X':{
        'files':[
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0057C0BD-B380-D74A-8634-BA648F4FF669.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0151F40E-240B-374C-AA69-44EBCE093268.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/017F2874-3AD9-074C-BB1E-E061D429E3F1.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0440EA45-411F-E04F-92A6-B3DDA0AA9DDD.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/05DDBC2E-90BE-044A-88A9-0049686243DF.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0939BA48-B40B-444B-9665-1743D6017488.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0C8090D3-3B84-A142-B22B-F81541642824.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0DF86428-7817-0B40-AF6D-421E65B64647.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0F3A06FA-7BDF-0B49-9EC8-959D2257A614.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/10EAD7D6-C405-D144-9308-0FA1416FBD72.root',
        ],
        'color':ROOT.kBlue,
    },
    'Charmonium_Run2018D_PromptReco':{
        'files':[
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/0067F915-7DB8-E811-A3B5-FA163E3C94E9.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/00A68F5A-8AAF-E811-800C-FA163E32FB27.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/00EE8907-CA9E-E811-A3C9-FA163E1858E5.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/02093F10-FA9D-E811-99F8-FA163EA55BB8.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/021B4E0E-F8E7-E84B-8D73-708675E13E77.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/022E11A6-56A0-E811-AD2B-FA163E205A29.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/023EC343-EDAE-E811-A17C-FA163EDAA78B.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/0268C30A-CEAF-E811-9604-02163E015F3C.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/02995A20-B0B3-E811-9F1A-FA163EB2B58A.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/502/Charmonium+Run2018D-PromptReco-v2+MINIAOD/02AEB5E0-CFAE-E811-BF92-FA163EB9E8EA.root',
            
        ],
        'color':ROOT.kBlack,
    }
}

# # read list of files from a file instead

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

def plot_generic_1D(selection,hist_title,file_name,var,nbins=100,xmin=0,xmax=100):
    c1 = TCanvas("c1","c1",800,800)
    max_value = 0
    for name,sample in samples.items():
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

plot_generic_1D(mm_selection,  "Dimuon vertex constrained;Mass",           "02_mass","mm_kin_mass",100,2.9,3.3)
plot_generic_1D(bkmm_selection,"BtoJ/#psiK Kaon_{pt}>1GeV;Mass",          "02_mass_jpsik","bkmm_jpsimc_mass",100,4.5,6.0)
plot_generic_1D(bkmm_selection+"&&bkmm_kaon_pt>2","BtoJ/#psiK Kaon_{pt}>2GeV;Mass",          
                "02_mass_jpsik_kaon2","bkmm_jpsimc_mass",100,4.5,6.0)
plot_generic_1D(bkmm_selection,"J/#psi from BtoJ/#psiK Kaon_{pt}>1GeV;Mass",              
                "02_mass_jpsi","mm_kin_mass[bkmm_mm_index]",100,2.9,3.3)
plot_generic_1D(bkmm_selection+"&&bkmm_kaon_pt>2","J/#psi from BtoJ/#psiK Kaon_{pt}>2GeV;Mass",              
                "02_mass_jpsi_kaon2","mm_kin_mass[bkmm_mm_index]",100,2.9,3.3)
plot_generic_1D(bkmm_selection+"&&abs(bkmm_jpsimc_mass-5.29)<0.05", "J/#psi from BtoJ/#psiK (B_{M}#in[5.24,5.34]);Mass", 
                "02_mass_jpsi_bmasscut","mm_kin_mass[bkmm_mm_index]",100,2.9,3.3)

plot_generic_1D(bkmm_selection+"&&abs(bkmm_jpsimc_mass-5.29)<0.05","BDT", "09_bdt","mm_bdt[bkmm_mm_index]",100,-1.5,1.5)

# plot_generic_1D("Decay length 3D significance",    "01_decay_length_3D_significance","mm_kin_sl3d",60,0,120)
# plot_generic_1D("Decay length 3D",                 "01_decay_length_3D","mm_kin_l3d",100,0,1)
# plot_generic_1D("Decay length XY significance",    "01_decay_length_XY_significance","mm_kin_slxy",60,0,120)
# plot_generic_1D("Pointing angle;#alpha",           "01_alpha","mm_kin_alpha",100,0.0,0.3)

# plot_generic_1D("Bs vertex constrained;Mass",        "02_mass","mm_kin_mass",70,5.0,5.7)
# plot_generic_1D("Leading muon;P_T,[GeV]",            "02_m1pt","mm_kin_mu1pt",100,0.0,50)
# plot_generic_1D("Trailing muon;P_T,[GeV]",           "02_m2pt","mm_kin_mu2pt",100,0.0,50)
# plot_generic_1D("Leading muon;#eta",                 "02_m1eta","mm_kin_mu1eta",50,-2.5,2.5)
# plot_generic_1D("Trailing muon;#eta",                "02_m2eta","mm_kin_mu2eta",50,-2.5,2.5)
# plot_generic_1D("Leading muon;#phi",                 "02_m1phi","mm_kin_mu1phi",32,-3.2,3.2)
# plot_generic_1D("Trailing muon;#phi",                "02_m2phi","mm_kin_mu2phi",32,-3.2,3.2)
# plot_generic_1D("Bs vertex constrained;P_T,[GeV]",   "02_pt","mm_kin_pt",100,0,50)
# plot_generic_1D("Bs vertex constrained;#eta",        "02_eta","mm_kin_eta",50,-2.5,2.5)

# plot_generic_1D("Muon distance of closest approach", "03_mmdoca","mm_doca",100,0.0,0.03)

# plot_generic_1D("Normalized #chi^{2} vertex fit;#frac{#chi^{2}}{N_{dof}}", "04_chi2dof","mm_kin_vtx_chi2dof",100,0.0,2)
# plot_generic_1D("Vertex fit probability;Prob",                             "04_vtx_prob","mm_kin_vtx_prob",100,0.0,1)

# plot_generic_1D("Impact parameter 3D wrt PV",                            "05_pvip","mm_kin_pvip",100,0.0,0.1)
# plot_generic_1D("Significance of Impact parameter 3D wrt PV",            "05_pvips","mm_kin_pvip/mm_kin_pvipErr",100,0.0,10)
# plot_generic_1D("Longitudinal Impact parameter wrt PV",                  "05_pvlip","mm_kin_pvlip",100,-0.1,0.1)
# plot_generic_1D("Significance of Longitudinal Impact parameter wrt PV",  "05_pvlips","mm_kin_pvlip/mm_kin_pvlipErr",100,-10.0,10)
# plot_generic_1D("Longitudinal Impact parameter wrt PV2",                 "05_pv2lip","mm_kin_pv2lip",100,-1,1)
# plot_generic_1D("Significance of Longitudinal Impact parameter wrt PV2", "05_pv2lips","mm_kin_pv2lip/mm_kin_pv2lipErr",100,-100.0,100.0)

# plot_generic_1D("closetrk",   "06_closetrk","mm_closetrk",100,0,100)
# plot_generic_1D("closetrks1", "06_closetrks1","mm_closetrks1",100,0,100)
# plot_generic_1D("closetrks2", "06_closetrks2","mm_closetrks2",100,0,100)
# plot_generic_1D("closetrks3", "06_closetrks3","mm_closetrks3",100,0,100)
# plot_generic_1D("docatrk",    "06_docatrk","mm_docatrk",100,0,0.03)

# plot_generic_1D("m1iso",   "07_m1iso","mm_m1iso",100,0,1.5)
# plot_generic_1D("m2iso",   "07_m2iso","mm_m2iso",100,0,1.5)
# plot_generic_1D("iso",     "07_iso","mm_iso",100,0,1.5)

# plot_generic_1D("othervtx pt>0.5", "08_othervtx0.5","mm_otherVtxMaxProb-mm_kin_vtx_prob",100,-1.5,1.5)
# plot_generic_1D("othervtx pt>1.0", "08_othervtx1.0","mm_otherVtxMaxProb1-mm_kin_vtx_prob",100,-1.5,1.5)
# plot_generic_1D("othervtx pt>2.0", "08_othervtx2.0","mm_otherVtxMaxProb2-mm_kin_vtx_prob",100,-1.5,1.5)
# plot_generic_1D("othervtx prob pt>2.0", "08_othervtx2.0prob","mm_otherVtxMaxProb2",100,-1.5,1.5)



