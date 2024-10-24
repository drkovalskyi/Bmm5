import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob

ROOT.ROOT.EnableImplicitMT()

output_path = "/eos/home-d/dmytro/www/plots/mc_bjpsik_vs_bsjpsiphi"
cache_path = "cache/"
rebuild_cache = False
# use_cache = False
use_cache = True
n_max = 999
# n_max = 1

# preselection must use one branch only (no cross links)
bkkmm_preselection = \
    "bkkmm_jpsikk_vtx_prob>0.025 && bkkmm_jpsikk_sl3d>3" \
    "&& bkkmm_jpsikk_alpha<0.1 && abs(bkkmm_jpsikk_mass-5.4)<0.5" \
    "&& abs(bkkmm_kk_mass-1.02)<0.03" \
    "&& bkkmm_gen_pdgId!=0"

bkkmm_selection = bkkmm_preselection + \
    "&& mm_mu1_index[bkkmm_mm_index]>=0 && mm_mu2_index[bkkmm_mm_index]>=0" \
    "&& abs(mm_mu1_eta[bkkmm_mm_index])<1.4 && mm_mu1_pt[bkkmm_mm_index]>4" \
    "&& abs(mm_mu2_eta[bkkmm_mm_index])<1.4 && mm_mu2_pt[bkkmm_mm_index]>3" \
    "&& mm_kin_vtx_prob[bkkmm_mm_index]>0.01"

# preselection must use one branch only (no cross links)
bkmm_preselection = \
    "bkmm_jpsimc_vtx_prob>0.025 && bkmm_jpsimc_sl3d>3" \
    "&& bkmm_jpsimc_alpha<0.1 && abs(bkmm_jpsimc_mass-5.4)<0.5" \
    "&& bkmm_gen_pdgId!=0"

bkmm_selection = bkmm_preselection + \
    "&& mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0" \
    "&& abs(mm_mu1_eta[bkmm_mm_index])<1.4 && mm_mu1_pt[bkmm_mm_index]>4" \
    "&& abs(mm_mu2_eta[bkmm_mm_index])<1.4 && mm_mu2_pt[bkmm_mm_index]>3" \
    "&& mm_kin_vtx_prob[bkmm_mm_index]>0.01"

samples = {
    'BsToJpsiPhi':{
        'datasets':[
            "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2+MINIAODSIM/",
            "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2+MINIAODSIM/",
            "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
            "BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        ],
        'color':ROOT.kBlue,
        'type':'mmkk',
        'legend':'B_{s} #rightarrow J/#psi#phi',
        'skim_selection': bkkmm_preselection,
        'skim_keep_pattern': '^(nbkkmm|bkkmm_.*|nmm|mm_.*)$'
    },
    'BuToJpsiK':{
        'datasets':[
            "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",
            "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
            "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        ],
        'color':ROOT.kRed,
        'type':'mmk',
        'legend':'B #rightarrow J/#psiK',
        'skim_selection': bkmm_preselection,
        'skim_keep_pattern': '^(nbkmm|bkmm_.*|nmm|mm_.*)$'
    }
}

directory = '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/'

# Load data
mask = '\w\d\_\.'
for sample_name, sample_info in samples.items():
    
    # make sure that the sample name is valid
    if re.search(f"[^{mask}]", sample_name):
        raise Exception(f"Illigal symbol used in sample name {sample_name}. Allowed regexp: '{mask}'")
    
    cache_file = f"{cache_path}/{sample_name}.root"
    if not use_cache or rebuild_cache or not os.path.exists(cache_file):
        sample_info['files'] = []
        for ds in sample_info['datasets']:
            sample_info['files'].extend(glob.glob(f"{directory}/{ds}/*.root"))
        print(f"Number of files available for {sample_name}: {len(sample_info['files'])}")
    
        chain = ROOT.TChain("Events")
        for entry in sample_info['files'][:n_max]:
            chain.Add(entry)
        sample_info['events'] = chain
        
        if use_cache and 'skim_selection' in sample_info:
            df = ROOT.RDataFrame(chain)
            df = df.Define("goodCandidates", sample_info['skim_selection'])
            df = df.Filter("Sum(goodCandidates) > 0", "Event has good candidates")
            df.Snapshot("Events", cache_file, sample_info['skim_keep_pattern'])

    if use_cache and os.path.exists(cache_file):
        chain = ROOT.TChain("Events")
        chain.Add(cache_file)
        sample_info['events'] = chain
        print(f"Using cache for {sample_name}")
    
print(f"Number of events using not more than {n_max} files:")
for name, sample in samples.items():
    sample['nAll'] = sample['events'].GetEntries()
    print(f"\t{name}: \t{sample['nAll']:,}")

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def plot_generic_1D(selections, hist_title, file_name,
                    vars, nbins=100, xmin=0, xmax=100, legend_left=True):
    c1 = TCanvas("c1", "c1", 600, 600)
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

    if legend_left:
        legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
    else:
        legend = ROOT.TLegend(0.5,0.75,0.85,0.87)
    legend.SetFillStyle(0)
    legend.SetLineWidth(0)

    first_plot = True
    for name,sample in samples.items():
        sample['hist'].SetMinimum(0)
        sample['hist'].SetMaximum(max_value * 1.2)
        legend.AddEntry(sample['hist'], sample['legend'])
        if first_plot:
            sample['hist'].Draw("hist")
            first_plot = False
        else:
            sample['hist'].Draw("hist same")
    legend.Draw()
    print_canvas(file_name, output_path)
    print("Number of selected events:")
    for name,sample in samples.items():
        print(f"\t{name}: \t{sample['nSelected']} out of {sample['nAll']}")

def plot_ratio_1D(selections, hist_title, file_name, sample_name1, sample_name2,
                  vars, nbins=100, xmin=0, xmax=100):
    c1 = TCanvas("c1", "c1", 600, 600)

    for name in [sample_name1, sample_name2]:
        sample = samples[name]
        selection = selections[sample['type']]
        var = vars[sample['type']]
        hist = ROOT.TH1D("hist", hist_title, nbins, xmin, xmax)
        hist.Sumw2()
        hist.SetMarkerColor(sample['color'])
        hist.SetMarkerStyle(20)
        sample['nSelected'] = sample['events'].Draw(f"{var}>>hist", selection)
        # print_canvas("%s_%s"%(file_name,name), output_path)
        if hist.GetEntries()>0:
            hist.Scale(1 / hist.GetEntries()) # normalize
        hist.SetDirectory(0)
        sample['hist'] = hist
    
    samples[sample_name1]['hist'].Divide(samples[sample_name2]['hist'])
    samples[sample_name1]['hist'].Draw("e0")
    print_canvas(file_name, output_path)
    print("Number of selected events:")
    for name in [sample_name1, sample_name2]:
        sample = samples[name]
        print(f"\t{name}: \t{sample['nSelected']} out of {sample['nAll']}")


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

plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection},
                "Dimuon p_{T};p^{#mu#mu}_{T}, [GeV]", "01_mm_pt",
                {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>1.5"},
                "Dimuon p_{T};p^{#mu#mu}_{T}, [GeV]", "01_mm_pt_kaon1.5",
                {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>2.0"},
                "Dimuon p_{T};p^{#mu#mu}_{T}, [GeV]", "01_mm_pt_kaon2.0",
                {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>2.5"},
                "Dimuon p_{T};p^{#mu#mu}_{T}, [GeV]", "01_mm_pt_kaon2.5",
                {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>3.0"},
                "Dimuon p_{T};p^{#mu#mu}_{T}, [GeV]", "01_mm_pt_kaon3.0",
                {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
                50, 0, 50, False)

plot_ratio_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection},
              "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio",
              "BsToJpsiPhi", "BuToJpsiK",
              {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
              30, 0, 30)
plot_ratio_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>1.5"},
              "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon1.5",
              "BsToJpsiPhi", "BuToJpsiK",
              {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
              30, 0, 30)
plot_ratio_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>2.0"},
              "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon2.0",
              "BsToJpsiPhi", "BuToJpsiK",
              {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
              30, 0, 30)
plot_ratio_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>2.5"},
              "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon2.5",
              "BsToJpsiPhi", "BuToJpsiK",
              {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
              30, 0, 30)
plot_ratio_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>3.0"},
              "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon3.0",
              "BsToJpsiPhi", "BuToJpsiK",
              {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
              30, 0, 30)

plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection},
                "B meson momentum;p^{B}_{T}, [GeV]", "01_B_pt",
                {'mmkk':'bkkmm_jpsikk_pt', 'mmk':'bkmm_jpsimc_pt'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>1.5"},
                "B meson momentum;p^{B}_{T}, [GeV]", "01_B_pt_kaon1.5",
                {'mmkk':'bkkmm_jpsikk_pt', 'mmk':'bkmm_jpsimc_pt'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>2.0"},
                "B meson momentum;p^{B}_{T}, [GeV]", "01_B_pt_kaon2.0",
                {'mmkk':'bkkmm_jpsikk_pt', 'mmk':'bkmm_jpsimc_pt'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>2.5"},
                "B meson momentum;p^{B}_{T}, [GeV]", "01_B_pt_kaon2.5",
                {'mmkk':'bkkmm_jpsikk_pt', 'mmk':'bkmm_jpsimc_pt'},
                50, 0, 50, False)
plot_generic_1D({'mmkk':bkkmm_selection, 'mmk':bkmm_selection + "&& bkmm_kaon_pt>3.0"},
                "B meson momentum;p^{B}_{T}, [GeV]", "01_B_pt_kaon3.0",
                {'mmkk':'bkkmm_jpsikk_pt', 'mmk':'bkmm_jpsimc_pt'},
                50, 0, 50, False)

