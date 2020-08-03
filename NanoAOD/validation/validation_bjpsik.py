#!/bin/env python
import os, re, ROOT, sys, time
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-507/bjpsik/"

bkmm_cuts = [
    # skimming requirements
    "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu1_index[bkmm_mm_index]]>4",
    "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4", "Muon_pt[mm_mu2_index[bkmm_mm_index]]>4",
    "mm_kin_sl3d[bkmm_mm_index]>4", "mm_kin_vtx_chi2dof[bkmm_mm_index]<5",
    "abs(bkmm_jpsimc_mass-5.4)<0.5", "bkmm_jpsimc_vtx_chi2dof<5","bkmm_jpsimc_alpha<0.2",
    # trigger requirements
    
    # muon id
    "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 ", " Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45",
]
bkmm_selection = "&&".join(bkmm_cuts)

save_root = True
save_pdf  = False

samples = {
    'BuToJpsiK_RunIIAutumn18MiniAOD_102X':{
        'files':[
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0bb34079aa40de8c3fd5615d4816875f.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0e73043348cb10591fd0695444457eb0.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/11f35da73cfdc279a1c462ff4545171d.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/34879d46edcc4751275237a9a5f5e005.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/4a098cd804e128ff3884d364a83ca001.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/5ee0b86054c0cd06cdd8321aa7bd6ad7.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/69256771e201020cf2cb897dcf8a8c50.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/83cca0d22ac9b759c372b35e26383e75.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/858c7ecc88235b88b48d2781029c589b.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/859639c395e8b78e5961a529c8c40e57.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/95f810c82fcd5d0f28c797d434567fd0.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/96ef9d3d833bc619b9302acdaaad4bae.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ae3371239732054142154e65ac6704da.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/bf3334e398d5fe68e5b3b90902b98211.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/c1fccea24bea988bf44999dcbe392944.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/c4829a75fb3d015222a50537e11cf5c5.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/d0f387b0f422c1a5fc9d27473d2aac84.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ec72869cbf788f4b6f7ff70009f1ec71.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ee13334ead0829b10ce5e7451ee11245.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/f9156ede12c49a38bb79033cb8ee80f1.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/feb8d6652fdd5e9037dda99f1a01cbe4.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/ff104c1d47f6100eee007c5f7d7e348b.root',
        ],
        'color':ROOT.kBlue,
    },
    'Charmonium_Run2018D_PromptReco':{
        'files':[
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/0598c756298cd841ea1ba52aa6205eaa.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/093134082bc6cd87be174850e72892c2.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/142c84fb8e8d03eea34e0cf3277a0656.root',
            '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/171c485619d24e7b2de271e713eca657.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/3a2e7df0fa80cecb95f8f3adcf209df1.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/3f4c01f6a59a852150e3a7a30319a8b3.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/41c385f0cbe95de10b190b3d64ad7bdc.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/4743df6d0cc32407bf84a5e83f2404a8.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/6cba9feb4a653616c839744db0f2a04b.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/71cf6c7a252fa917ab6a61b3c60d1b05.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/72080df52e6a36f3e85c6b960daaf510.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/729046729c91eeaef1ebfa1e320ab2d2.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/8e58a2b22e745af953fb0ce4763d8453.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/92c411e0b3a00b2b0606a6b12f009aa8.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/974e71c94fc6b826fc5a57d0a1e11351.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/b9a1a2f74e56da3db7871bd16cff191b.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/d038d09060a5b5487cdf4de572362e2f.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/d606c6838c2bfe67273e6f1bb7f8ecff.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/ddf50979ae49179f598e49ab87a8ceeb.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/e6fc8ef8028c22b1b9489374220f1d4e.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/e88d8f1444f66a53cf4098070036100f.root',
            # '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/f016487d3c89be36327f2575e958c9bf.root',
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
    canvas.Print("%s/%s.png"%(path, output_name_without_extention))
    if save_pdf:
        canvas.Print("%s/%s.pdf"%(path, output_name_without_extention))
    if save_root:
        canvas.Print("%s/%s.root"%(path, output_name_without_extention))

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

def plot_resolution(selection, hist_title, file_name, 
                    profile_var, profile_bins,
                    observable, nbins=100, xmin=-5.0, xmax=5.0):
    """Make a profile plot of the resolution with respect to the profile
    variable."""

    max_value = 0
    n = len(profile_bins)
    if n < 2:
        raise Exception("Too few profile bins: %d" % (n-1))
    elif n < 4:
        c1 = TCanvas("c1", "c1", 1000, 500)
        c1.Divide(2, 1)
    elif n < 6:
        c1 = TCanvas("c1", "c1", 1000, 1000)
        c1.Divide(2, 2)
    elif n < 8:
        c1 = TCanvas("c1", "c1", 1500, 1000)
        c1.Divide(3, 2)
    elif n < 11:
        c1 = TCanvas("c1", "c1", 1500, 1500)
        c1.Divide(3, 3)
    elif n < 14:
        c1 = TCanvas("c1", "c1", 2000, 1500)
        c1.Divide(4, 3)
    elif n < 18:
        c1 = TCanvas("c1", "c1", 2000, 2000)
        c1.Divide(4, 4)
    else:
        raise Exception("Too many profile bins: %d" % (n-1))

    for name,sample in samples.items():
        results = []
        sample['hists'] = []
        x = array("f")
        x_err = array("f")
        y = array("f")
        y_err = array("f")
        for i in range(n-1):
            c1.cd(i+1)
            x.append((profile_bins[i+1] + profile_bins[i]) / 2)
            x_err.append((profile_bins[i+1] - profile_bins[i]) / 2)
            hist = ROOT.TH1D("hist",hist_title,nbins,xmin,xmax)
            hist.SetLineColor(sample['color'])
            hist.SetLineWidth(2)
            cuts = selection + "&&%s>=%s&&%s<%s" % (profile_var, profile_bins[i], profile_var, profile_bins[i+1])
            sample['events'].Draw("%s>>hist" % observable, cuts)
            result = hist.Fit("gaus","S")
            results.append([result.Parameter(2),result.ParError(2)])
            y.append(result.Parameter(2))
            y_err.append(result.ParError(2))
            hist.SetDirectory(0)
            sample['hists'].append(hist)
            
        print results
        print_canvas("%s_%s"%(file_name,name), output_path, c1)

        gr = ROOT.TGraphErrors(len(x), x, y, x_err, y_err)
        # gr.SetLineColor(sample['color'])
        gr.SetLineWidth(2)
        gr.SetTitle(hist_title)
        c2 = TCanvas("c2", "c2", 800, 800)
        gr.Draw("AP")
        print_canvas("%s_graph_%s"%(file_name,name), output_path, c2)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# plot_generic_1D(bkmm_selection, "BtoJ/#psiK (kaon P_T>1 GeV);Mass",
#                 "02_mass_jpsik", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_kaon_pt>2", "BtoJ/#psiK (kaon P_T>2 GeV);Mass",          
#                 "02_mass_jpsik_kaon2", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_kaon_pt>3", "BtoJ/#psiK (kaon P_T>3 GeV);Mass",          
#                 "02_mass_jpsik_kaon3", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_kaon_pt>5", "BtoJ/#psiK (kaon P_T>5 GeV);Mass",          
#                 "02_mass_jpsik_kaon5", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_bmm_mva>0.5", "BtoJ/#psiK (modified MVA>0.5);Mass",
#                 "02_mass_jpsik_mva_gt0p5", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_bmm_mva>0.9", "BtoJ/#psiK (modified MVA>0.9);Mass",
#                 "02_mass_jpsik_mva_gt0p9", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_kaon_pt>5&&bkmm_bmm_mva>0.9", "BtoJ/#psiK (modified MVA>0.9, kaon P_T>5 GeV);Mass",
#                 "02_mass_jpsik_mva_gt0p9_kaon5", "bkmm_jpsimc_mass", 100, 4.5, 6.0)
# plot_generic_1D(bkmm_selection + "&&bkmm_kaon_pt>5&&bkmm_bmm_mva>0.9", "BtoJ/#psiK (modified MVA>0.9, kaon P_T>5 GeV);Mass",
#                 "02_mass_jpsik_mva_gt0p9_kaon5_zoom", "bkmm_jpsimc_mass", 100, 5.15, 5.55)
plot_resolution(bkmm_selection + 
                "&&bkmm_kaon_pt>5&&bkmm_bmm_mva>0.9&&!TMath::IsNaN(bkmm_jpsimc_massErr)&&abs(bkmm_jpsimc_mass-5.28)<0.13", 
                "BtoJ/#psiK (modified MVA>0.9, kaon P_T>5 GeV, [5.15,5.41]);#sigma_M;pull width", 
                "10_resolution_vs_masserr_jpsik_mva_gt0p9_kaon5",
                "bkmm_jpsimc_massErr", [0, 0.013, 0.015, 0.017, 0.020, 0.025, 0.05],
                "(bkmm_jpsimc_mass-5.2792601)/bkmm_jpsimc_massErr")
plot_resolution(bkmm_selection + 
                "&&bkmm_kaon_pt>5&&bkmm_bmm_mva>0.9&&!TMath::IsNaN(bkmm_jpsimc_massErr)&&abs(bkmm_jpsimc_mass-5.28)<0.13", 
                "BtoJ/#psiK (modified MVA>0.9, kaon P_T>5 GeV, [5.15,5.41]);max(|#eta|);pull width", 
                "10_resolution_vs_eta_jpsik_mva_gt0p9_kaon5",
                "max(abs(Muon_eta[mm_mu1_index[bkmm_mm_index]]),abs(Muon_eta[mm_mu2_index[bkmm_mm_index]]))", 
                [0, 0.25, 0.50, 0.75, 1.00, 1.20, 1.40],
                "(bkmm_jpsimc_mass-5.2792601)/bkmm_jpsimc_massErr")
                
# plot_generic_1D(bkmm_selection, "J/#psi from BtoJ/#psiK Kaon_{pt}>1GeV;Mass",              
#                 "02_mass_jpsi","mm_kin_mass[bkmm_mm_index]", 100, 2.9, 3.3)
# plot_generic_1D(bkmm_selection + "&&bkmm_kaon_pt>2", "J/#psi from BtoJ/#psiK (kaon P_T>2 GeV);Mass",              
#                 "02_mass_jpsi_kaon2", "mm_kin_mass[bkmm_mm_index]", 100, 2.9, 3.3)
# plot_generic_1D(bkmm_selection + "&&abs(bkmm_jpsimc_mass-5.29)<0.05", "J/#psi from BtoJ/#psiK (B_{M}#in[5.24,5.34]);Mass", 
#                 "02_mass_jpsi_bmasscut", "mm_kin_mass[bkmm_mm_index]", 100, 2.9, 3.3)

# plot_generic_1D(bkmm_selection + "&&abs(bkmm_jpsimc_mass-5.29)<0.05", "BtoJ/#psiK (B_{M}#in[5.24,5.34], kaon P_T>1 GeV);BDT",
#                 "09_bdt", "mm_bdt[bkmm_mm_index]", 100, -1.5, 1.5)

# plot_generic_1D(bkmm_selection + "&&abs(bkmm_jpsimc_mass-5.29)<0.05&&bkmm_kaon_pt>2", "BtoJ/#psiK (B_{M}#in[5.24,5.34], kaon P_T>2 GeV);BDT", 
#                 "09_bdt_kaon2", "mm_bdt[bkmm_mm_index]", 100, -1.5, 1.5)
# plot_generic_1D(bkmm_selection + "&&abs(bkmm_jpsimc_mass-5.29)<0.05&&bkmm_kaon_pt>5", "BtoJ/#psiK (B_{M}#in[5.24,5.34], kaon P_T>5 GeV);BDT", 
#                 "09_bdt_kaon5", "mm_bdt[bkmm_mm_index]", 100, -1.5, 1.5)
# plot_generic_1D(bkmm_selection + "&&abs(bkmm_jpsimc_mass-5.29)<0.05&&bkmm_kaon_pt>5", "BtoJ/#psiK (B_{M}#in[5.24,5.34], kaon P_T>5 GeV);MVA", 
#                 "09_mva_kaon5", "mm_mva[bkmm_mm_index]", 120, -0.1, 1.1)
# plot_generic_1D(bkmm_selection + "&&abs(bkmm_jpsimc_mass-5.29)<0.05&&bkmm_kaon_pt>5", "BtoJ/#psiK (B_{M}#in[5.24,5.34], kaon P_T>5 GeV);MVA (modified)", 
#                 "09_mva_mod_kaon5", "bkmm_bmm_mva", 120, -0.1, 1.1)


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



