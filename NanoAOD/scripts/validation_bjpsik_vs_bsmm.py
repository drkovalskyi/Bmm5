#!/bin/env python
import os, re, ROOT, sys, time
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-507/mc_bjpsik_vs_bsmm"
quick_check = False # quick check uses only the first file for each sample

mm_cuts = [
    # skimming requirements
    "abs(Muon_eta[mm_mu1_index])<1.4", "Muon_pt[mm_mu1_index]>4",
    "abs(Muon_eta[mm_mu2_index])<1.4", "Muon_pt[mm_mu2_index]>4",
    "abs(mm_kin_mass-5.4)<0.5", "mm_kin_sl3d>4", "mm_kin_vtx_chi2dof<5",
    # muon id
    "Muon_softMva[mm_mu1_index]>0.45 ", " Muon_softMva[mm_mu2_index]>0.45",
    # gen matching
    "mm_gen_pdgId!=0",
    # tuning
    # "mm_kin_pt>10.*5.3/3.1", "mm_kin_pt<20*5.3/3.1",
    "mm_kin_alpha<0.2",
    
]
mm_selection = "&&".join(mm_cuts)

bkmm_cuts = [
    # skimming requirements
    "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4",
    "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4",
    "mm_kin_sl3d[bkmm_mm_index]>4", "mm_kin_vtx_chi2dof[bkmm_mm_index]<5",
    "abs(bkmm_jpsimc_mass-5.4)<0.5", 
    "bkmm_jpsimc_vtx_chi2dof<5", # need to get rid of it
    "bkmm_jpsimc_alpha<0.2",
    # muon id
    "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 ", " Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45",
    # gen matching
    "bkmm_gen_pdgId!=0",
    # kaon
    "bkmm_kaon_pt<1.5",
    # tuning
    # "mm_kin_pt[bkmm_mm_index]>10", "mm_kin_pt[bkmm_mm_index]<20",
]
bkmm_selection = "&&".join(bkmm_cuts)


# "bkmm_kaon_sdxy_bs>5&&bkmm_kaon_pt>1&&abs(bkmm_kaon_eta)<1.4" 

samples = {
    'BsToMuMu_RunIIAutumn18MiniAOD_102X':{
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/0ea357ae822b928036f590a1d8fd4281.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/1d6e11b1cfa2bb1e43dc0ba074e3c50b.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/35947643ae86e0eaad2e807ba9a0c6ae.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/3fa621813e82b51afb8326c1e605ac1c.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/455d900fb5d0cc3923aafac00f8e4955.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/4d5679a310048519e04e48331d412be5.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/5b8ceb58a5ee10032c92fcadeb0ca4f2.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/677935795b0818c9f24b237db4a3a29b.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/a23f15b66e4a231d2e22df8fa770c9dc.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/ac3944ff78c688f1b7e8b02b9201f677.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/ae79cec75ea35246052674e2ec90685e.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/b66781c76391873502d36b494dcdb8fe.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/ba7abeb9c53b701de336a5a1e6276141.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/bee26883b8c8fb2da395b151f298fabc.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/bef22b11b7b342a6d0bdae36245f8dd5.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/507/BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/d2d4f8c590beb69dac8f5fbdb5460697.root",
        ],
        'color':ROOT.kBlue,
        'type':'bmm',
    },
    'BuToJpsiK_RunIIAutumn18MiniAOD_102X':{
        'files':[
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0bb34079aa40de8c3fd5615d4816875f.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0e73043348cb10591fd0695444457eb0.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/11f35da73cfdc279a1c462ff4545171d.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/34879d46edcc4751275237a9a5f5e005.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/4a098cd804e128ff3884d364a83ca001.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/5ee0b86054c0cd06cdd8321aa7bd6ad7.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/69256771e201020cf2cb897dcf8a8c50.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/83cca0d22ac9b759c372b35e26383e75.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/858c7ecc88235b88b48d2781029c589b.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/859639c395e8b78e5961a529c8c40e57.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/95f810c82fcd5d0f28c797d434567fd0.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/96ef9d3d833bc619b9302acdaaad4bae.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ae3371239732054142154e65ac6704da.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/bf3334e398d5fe68e5b3b90902b98211.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/c1fccea24bea988bf44999dcbe392944.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/c4829a75fb3d015222a50537e11cf5c5.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/d0f387b0f422c1a5fc9d27473d2aac84.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ec72869cbf788f4b6f7ff70009f1ec71.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ee13334ead0829b10ce5e7451ee11245.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/f9156ede12c49a38bb79033cb8ee80f1.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/feb8d6652fdd5e9037dda99f1a01cbe4.root",
            "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ff104c1d47f6100eee007c5f7d7e348b.root",
        ],
        'color':ROOT.kRed,
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
        if quick_check: break
    sample['events'] = chain
    
print "Number of events:"
for name,sample in samples.items():
    sample['nAll'] = sample['events'].GetEntries()
    print "\t%s: \t%u" % (name,sample['nAll'])

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))

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
        # print_canvas("%s_%s"%(file_name,name), output_path)
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
    'bjpsik':bkmm_selection,
}

# Variables used in MVA
# * bkmm_jpsimc_alpha
# * bkmm_jpsimc_cosAlphaXY
# * bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr
# * bkmm_jpsimc_pvip
# * bkmm_bmm_m1iso
# * bkmm_bmm_m2iso
# * bkmm_bmm_iso
# * bkmm_bmm_nBMTrks
# * bkmm_bmm_otherVtxMaxProb1
# * bkmm_bmm_otherVtxMaxProb2
# * mm_kin_vtx_chi2dof
# * mm_kin_sl3d
#
# Should also consider
# * bkmm_jpsimc_vtx_chi2dof
# * bkmm_jpsimc_sl3d

plot_generic_1D(selections, "#mu#mu;P_{T}, [GeV]", "01_mm_pt",
                {'bmm':'mm_kin_pt', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'}, 100, 0, 100)

plot_generic_1D(selections, "#mu#mu vertex displacement significance;#sigma", "02_mm_sl3d",
                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]'}, 100, 0, 100)
plot_generic_1D(selections, "#mu#muK vertex displacement significance;#sigma", "02_kmm_sl3d",
                {'bmm':'mm_kin_sl3d', 'bjpsik':'bkmm_jpsimc_sl3d'}, 100, 0, 100)
# plot_generic_1D(selections, "#mu#mu vertex displacement significance (#mu#muK is scaled by 1.5);#sigma", "02_mm_sl3d_scaled",
#                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]*1.5'}, 100, 0, 100)
scale = 1.60
plot_generic_1D({'bmm':mm_selection + "&&mm_kin_sl3d>4*%s" % scale, 'bjpsik':bkmm_selection},
                "#mu#mu vertex displacement significance (#mu#muK is scaled by %s with matched seleciton);#sigma" % scale, 
                "02_mm_sl3d_scaled_and_matched_selection",
                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]*%s' % scale}, 100, 0, 100)

plot_generic_1D(selections, "Pointing angle 3D;#alpha", "03_alpha",
                {'bmm':'mm_kin_alpha', 'bjpsik':'bkmm_jpsimc_alpha'}, 100, 0, 0.3)
plot_generic_1D(selections, "Pointing angle XY;cos(#alpha)", "03_cosAlphaXY",
                {'bmm':'mm_kin_cosAlphaXY', 'bjpsik':'bkmm_jpsimc_cosAlphaXY'}, 110, 0.990, 1.001)

plot_generic_1D(selections, "Impact parameter significance", "04_spvip",
                {'bmm':'mm_kin_pvip/mm_kin_pvipErr', 'bjpsik':'bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr'},
                100, 0, 5)
plot_generic_1D(selections, "Impact parameter 3D", "04_pvip",
                {'bmm':'mm_kin_pvip', 'bjpsik':'bkmm_jpsimc_pvip'}, 100, 0, 0.02)

plot_generic_1D(selections, "#mu#mu isolation", "05_iso",
                {'bmm':'mm_iso', 'bjpsik':'bkmm_bmm_iso'}, 120, 0, 1.2)
plot_generic_1D(selections, "#mu1 isolation", "05_m1iso",
                {'bmm':'mm_m1iso', 'bjpsik':'bkmm_bmm_m1iso'}, 120, 0, 1.2)
plot_generic_1D(selections, "#mu2 isolation", "05_m2iso",
                {'bmm':'mm_m2iso', 'bjpsik':'bkmm_bmm_m2iso'}, 120, 0, 1.2)

plot_generic_1D(selections, "#chi/nDof for #mu#mu vertex", "06_mm_chi2dof",
                {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'mm_kin_vtx_chi2dof[bkmm_mm_index]'}, 100, 0, 5)
plot_generic_1D(selections, "#chi/nDof for #mu#muK vertex", "06_mmK_chi2dof",
                {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'bkmm_jpsimc_vtx_chi2dof'}, 100, 0, 5)

plot_generic_1D(selections, "nBMTrks", "07_nBMTrks",
                {'bmm':'mm_nBMTrks','bjpsik':'min(bkmm_bmm_nBMTrks,9)'}, 10, 0, 10)
plot_generic_1D(selections, "otherVtxMaxProb1", "07_otherVtxMaxProb1",
                {'bmm':'mm_otherVtxMaxProb1', 'bjpsik':'bkmm_bmm_otherVtxMaxProb1'}, 120, 0, 1.2)
plot_generic_1D(selections, "otherVtxMaxProb2", "07_otherVtxMaxProb2",
                {'bmm':'mm_otherVtxMaxProb2', 'bjpsik':'bkmm_bmm_otherVtxMaxProb2'}, 120, 0, 1.2)

plot_generic_1D(selections, "BDT Matched", "09_bdt_matched",
                {'bmm':'mm_bdt', 'bjpsik':'bkmm_bmm_bdt'}, 100, -1.5, 1.5)
plot_generic_1D(selections, "BDT Raw", "09_bdt_raw",
                {'bmm':'mm_bdt', 'bjpsik':'mm_bdt[bkmm_mm_index]'}, 100, -1.5, 1.5)

plot_generic_1D(selections,"MVA Matched", "09_mva_matched",
                {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)


# plot_generic_1D(selections,"BDT Raw", "02_mm_pt",{'bmm':'mm_pt','bjpsik':'mm_pt[bkmm_mm_index]'},100,0,100)
