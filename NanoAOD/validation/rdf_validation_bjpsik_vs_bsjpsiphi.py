import os, re, ROOT, sys, time, subprocess
from ROOT import TFile,TTree,TH1,TROOT,TDirectory,TPad,TCanvas,TColor
from array import array
import glob

output_path = "/eos/home-d/dmytro/www/plots/tmp/mc_bjpsik_vs_bsjpsiphi"
n_max = 999
n_max = 1

ROOT.ROOT.EnableImplicitMT(8)

bkkmm_cuts = """
mm_mu1_index[bkkmm_mm_index]>=0 && mm_mu2_index[bkkmm_mm_index]>=0 &&
abs(Muon_eta[mm_mu1_index[bkkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkkmm_mm_index]]>4 &&
abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkkmm_mm_index]]>3 &&
mm_kin_vtx_prob[bkkmm_mm_index]>0.01 && bkkmm_jpsikk_vtx_prob>0.025 &&
bkkmm_jpsikk_sl3d>3 && bkkmm_jpsikk_alpha<0.1 &&
abs(bkkmm_jpsikk_mass-5.4)<0.5 && abs(bkkmm_kk_mass-1.02)<0.03 &&
bkkmm_gen_pdgId!=0
"""

bkmm_cuts = """
mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 &&
abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 &&
abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>3 &&
mm_kin_vtx_prob[bkmm_mm_index]>0.01 && bkmm_jpsimc_vtx_prob>0.025 &&
bkmm_jpsimc_sl3d>3 && bkmm_jpsimc_alpha<0.1 &&
abs(bkmm_jpsimc_mass-5.4)<0.5 &&
bkmm_gen_pdgId!=0
"""

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
    },
    'BuToJpsiK':{
        'datasets':[
            "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2+MINIAODSIM/",
            "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3+MINIAODSIM/",
            "ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3+MINIAODSIM/",
        ],
        'color':ROOT.kRed,
        'type':'mmk',
        'legend':'B #rightarrow J/#psiK'
    }
}

directory = '/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/529/'

# Find root files
for sample, info in samples.items():
    info['files'] = []
    for ds in info['datasets']:
        info['files'].extend(glob.glob(f"{directory}/{ds}/*.root"))
    print(f"Number of files available for {sample}: {len(info['files'])}")

# Load data

mask = '\w\d\_\.'
for name, sample in samples.items():
    if re.search(f"[^{mask}]", name):
        raise Exception(f"Illigal symbol used for sample name. Allowed regexp: '{mask}'") 
    chain = ROOT.TChain("Events")
    for entry in sample['files'][:n_max]:
        chain.Add(entry)
    sample['events'] = chain
    sample['rdf'] = ROOT.RDataFrame(chain)
    
print(f"Number of events using not more than {n_max} files:")
for name,sample in samples.items():
    sample['nAll'] = sample['events'].GetEntries()
    print(f"\t{name}: \t{sample['nAll']:,}")

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png"%(path,output_name_without_extention))
    canvas.Print("%s/%s.pdf"%(path,output_name_without_extention))
    canvas.Print("%s/%s.root"%(path,output_name_without_extention))

def book_generic_1D(selections, hist_title, file_name,
                    vars, nbins=100, xmin=0, xmax=100, legend_left=True):
    for name, sample in samples.items():
        selection = selections[sample['type']]
        var = vars[sample['type']]
        hist_name = f"{file_name}_{name}"
        var_name = f"var_{file_name}_{name}"
        if 'histos' not in sample:
            sample['histos'] = dict()
        rd = sample['rdf'].Define("idx",selection)
        sample['histos'][hist_name] = \
            rd.Define(var_name,var).Histo1D((hist_name, hist_title, nbins, xmin, xmax), var_name)

def plot_generic_1D(selections, hist_title, file_name,
                    vars, nbins=100, xmin=0, xmax=100, legend_left=True):
    c1 = TCanvas("c1", "c1", 800, 800)
    max_value = 0
    for name, sample in samples.items():
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
        legend = ROOT.TLegend(0.5,0.75,0.75,0.87)
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
    c1 = TCanvas("c1", "c1", 800, 800)

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

# def integrate(source_hist, destination_hist, left_to_right = True):
#     assert(source_hist.GetNbinsX() == destination_hist.GetNbinsX())
#     assert(source_hist.Integral(0,-1) > 0)
#     for bin in range(source_hist.GetNbinsX() + 1):
#         if left_to_right:
#             destination_hist.SetBinContent(bin, source_hist.Integral(bin,-1)/source_hist.Integral(0,-1))
#         else:
#             destination_hist.SetBinContent(bin, 1 - source_hist.Integral(bin,-1)/source_hist.Integral(0,-1))

# def plot_integral_1D(selections, hist_title, file_name, vars, 
#                      nbins=100, xmin=0, xmax=100, 
#                      left_to_right = True):
#     """Make disitributions of the target variables and make intergral
#        plots. By default the integral is from x till xmax.
#     """
#     c1 = TCanvas("c1", "c1", 800, 800)
#     max_value = 0
#     for name,sample in samples.items():
#         selection = selections[sample['type']]
#         var = vars[sample['type']]
#         hist = ROOT.TH1D("hist", hist_title, nbins, xmin, xmax)
#         hist_integrated = ROOT.TH1D("hist_integrated", hist_title, nbins, xmin, xmax)
#         hist_integrated.SetLineColor(sample['color'])
#         hist_integrated.SetLineWidth(2)
#         sample['nSelected'] = sample['events'].Draw("%s>>hist" % var, selection)
#         integrate(hist, hist_integrated, left_to_right)
#         # print_canvas("%s_%s"%(file_name,name), output_path)
#         if max_value < hist_integrated.GetMaximum(): max_value = hist_integrated.GetMaximum()
#         hist_integrated.SetDirectory(0)
#         hist.SetDirectory(0)
#         sample['hist'] = hist_integrated
    
#     legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
#     legend.SetFillStyle(0)
#     legend.SetLineWidth(0)

#     first_plot = True
#     for name,sample in samples.items():
#         sample['hist'].SetMinimum(0)
#         sample['hist'].SetMaximum(max_value * 1.2)
#         legend.AddEntry(sample['hist'], sample['legend'])
#         if first_plot:
#             sample['hist'].Draw("hist")
#             first_plot = False
#         else:
#             sample['hist'].Draw("hist same")
#     legend.Draw()
#     print_canvas(file_name, output_path)
#     print "Number of selected events:"
#     for name,sample in samples.items():
#         print "\t%s: \t%u out of %u" % (name,sample['nSelected'],sample['nAll'])

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

book_generic_1D({'mmkk': bkkmm_cuts, 'mmk': bkmm_cuts},
                "Dimuon p^{T} spectrum;p^{#mu#mu}_{T}, [GeV]", "01_mm_pt",
                {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
                50, 0, 50, False)

histos = []
for name, sample in samples.items():
    if 'histos' in sample:
        for h_name, h_ref in sample['histos'].items():
            print(h_name)
            histos.append(h_ref)
print(len(histos))
ROOT.RDF.RunGraphs(histos)
print("done with processing")

# plot_generic_1D({'mmkk': bkkmm_cuts, 'mmk': bkmm_cuts},
#                 "Dimuon p^{T} spectrum;p^{#mu#mu}_{T}, [GeV]", "01_mm_pt",
#                 {'mmkk':'mm_kin_pt[bkkmm_mm_index]', 'mmk':'mm_kin_pt[bkmm_mm_index]'},
#                 50, 0, 50, False)


# plot_generic_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>1.5"},
#                 "#mu#mu;P_{T}, [GeV]", "01_mm_pt_kaon1.5",
#                 {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#                 50, 0, 50, False)
# plot_generic_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>2.0"},
#                 "#mu#mu;P_{T}, [GeV]", "01_mm_pt_kaon2.0",
#                 {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#                 50, 0, 50, False)
# plot_generic_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>2.5"},
#                 "#mu#mu;P_{T}, [GeV]", "01_mm_pt_kaon2.5",
#                 {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#                 50, 0, 50, False)
# plot_generic_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>3.0"},
#                 "#mu#mu;P_{T}, [GeV]", "01_mm_pt_kaon3.0",
#                 {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#                 50, 0, 50, False)
# # plot_ratio_1D(selections, "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio",
# #                 "BsToJpsiPhi", "BuToJpsiK",
# #                 {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
# #                 30, 0, 30)

# plot_ratio_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>1.5"},
#               "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon1.5",
#               "BsToJpsiPhi", "BuToJpsiK",
#               {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#               30, 0, 30)
# plot_ratio_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>2.0"},
#               "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon2.0",
#               "BsToJpsiPhi", "BuToJpsiK",
#               {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#               30, 0, 30)
# plot_ratio_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>2.5"},
#               "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon2.5",
#               "BsToJpsiPhi", "BuToJpsiK",
#               {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#               30, 0, 30)
# plot_ratio_1D({'bsjpsiphi':bkkmm_selection, 'bjpsik':bkmm_selection + "&& bkmm_kaon_pt>3.0"},
#               "#mu#mu;P_{T}, [GeV]", "01_mm_pt_ratio_kaon3.0",
#               "BsToJpsiPhi", "BuToJpsiK",
#               {'bsjpsiphi':'mm_kin_pt[bkkmm_mm_index]', 'bjpsik':'mm_kin_pt[bkmm_mm_index]'},
#               30, 0, 30)


# # plot_generic_1D(selections, "#mu#mu vertex displacement significance;#sigma", "02_mm_sl3d",
# #                 {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]'}, 100, 0, 100)
# # plot_generic_1D(selections, "#mu#muK vertex displacement significance;#sigma", "02_kmm_sl3d",
# #                 {'bmm':'mm_kin_sl3d', 'bjpsik':'bkmm_jpsimc_sl3d'}, 100, 0, 100)
# # # plot_generic_1D(selections, "#mu#mu vertex displacement significance (#mu#muK is scaled by 1.5);#sigma", "02_mm_sl3d_scaled",
# # #                {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]*1.5'}, 100, 0, 100)
# # scale = 1.60
# # plot_generic_1D({'bmm':mm_selection + "&&mm_kin_sl3d>4*%s" % scale, 'bjpsik':bkmm_selection},
# #                 "#mu#mu vertex displacement significance (#mu#muK is scaled by %s);#sigma" % scale, 
# #                 "02_mm_sl3d_scaled_and_matched_selection",
# #                 {'bmm':'mm_kin_sl3d', 'bjpsik':'mm_kin_sl3d[bkmm_mm_index]*%s' % scale}, 100, 0, 100)

# # plot_generic_1D(selections, "Pointing angle 3D;#alpha_{3D}", "03_alpha",
# #                 {'bmm':'mm_kin_alpha', 'bjpsik':'bkmm_jpsimc_alpha'}, 100, 0, 0.2)
# # plot_generic_1D(selections, "Pointing angle BS;#alpha_{BS}", "03_alphaBS",
# #                 {'bmm':'mm_kin_alphaBS', 'bjpsik':'bkmm_jpsimc_alphaBS'}, 110, 0, 0.2)

# # plot_generic_1D(selections, "Impact parameter significance", "04_spvip",
# #                 {'bmm':'mm_kin_pvip/mm_kin_pvipErr', 'bjpsik':'bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr'},
# #                 100, 0, 5)
# # plot_generic_1D(selections, "Impact parameter 3D", "04_pvip",
# #                 {'bmm':'mm_kin_pvip', 'bjpsik':'bkmm_jpsimc_pvip'}, 100, 0, 0.02)

# # plot_generic_1D(selections, "#mu#mu isolation", "05_iso",
# #                 {'bmm':'mm_iso', 'bjpsik':'bkmm_bmm_iso'}, 120, 0, 1.2)
# # plot_generic_1D(selections, "#mu1 isolation", "05_m1iso",
# #                 {'bmm':'mm_m1iso', 'bjpsik':'bkmm_bmm_m1iso'}, 120, 0, 1.2)
# # plot_generic_1D(selections, "#mu2 isolation", "05_m2iso",
# #                 {'bmm':'mm_m2iso', 'bjpsik':'bkmm_bmm_m2iso'}, 120, 0, 1.2)

# # plot_generic_1D(selections, "#chi/nDof for #mu#mu vertex", "06_mm_chi2dof",
# #                 {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'mm_kin_vtx_chi2dof[bkmm_mm_index]'}, 100, 0, 5)
# # plot_generic_1D(selections, "#chi/nDof for #mu#muK vertex", "06_mmK_chi2dof",
# #                 {'bmm':'mm_kin_vtx_chi2dof', 'bjpsik':'bkmm_jpsimc_vtx_chi2dof'}, 100, 0, 5)

# # plot_generic_1D(selections, "nBMTrks", "07_nBMTrks",
# #                 {'bmm':'mm_nBMTrks','bjpsik':'min(bkmm_bmm_nBMTrks,9)'}, 10, 0, 10)
# # plot_generic_1D(selections, "otherVtxMaxProb1", "07_otherVtxMaxProb1",
# #                 {'bmm':'mm_otherVtxMaxProb1', 'bjpsik':'bkmm_bmm_otherVtxMaxProb1'}, 120, 0, 1.2)
# # plot_generic_1D(selections, "otherVtxMaxProb2", "07_otherVtxMaxProb2",
# #                 {'bmm':'mm_otherVtxMaxProb2', 'bjpsik':'bkmm_bmm_otherVtxMaxProb2'}, 120, 0, 1.2)

# # # plot_generic_1D(selections, "BDT Matched", "09_bdt_matched",
# # #                 {'bmm':'mm_bdt', 'bjpsik':'bkmm_bmm_bdt'}, 100, -1.5, 1.5)
# # # plot_generic_1D(selections, "BDT Raw", "09_bdt_raw",
# # #                 {'bmm':'mm_bdt', 'bjpsik':'mm_bdt[bkmm_mm_index]'}, 100, -1.5, 1.5)

# # plot_generic_1D(selections,"MVA Matched", "09_mva_matched",
# #                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)
# # plot_generic_1D(selections,"MVA Matched", "09_mva_matched_zoomed",
# #                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 101, 0.9, 1.01)
# # plot_generic_1D(selections,"MVA Raw", "09_mva_raw",
# #                 {'bmm':'mm_mva', 'bjpsik':'mm_mva[bkmm_mm_index]'}, 110, 0, 1.1)

# # plot_integral_1D(selections,"MVA Matched Efficiency", "09_mva_matched_eff",
# #                 {'bmm':'mm_mva', 'bjpsik':'bkmm_bmm_mva'}, 110, 0, 1.1)

