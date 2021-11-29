import os, re
import ROOT
import math 
import sweights
import tdrstyle
from ROOT.RooFit import Binning
from ROOT import RooRealVar

# Set the TDR style
tdrstyle.setTDRStyle()

version = 516

# output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv8-%u/bjpsik-splots_binned/" % version;
output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/AN/bjpsik-splots_binned/";


# bool silent_roofit = true;
# bool store_projections = false;
# bool use_mc_truth_matching = true;

# struct variable{
#   string name, branch;
#   unsigned int nbins;
#   double xmin, xmax;
# };

mass = ROOT.RooRealVar("mm_kin_mass", "", 5.15, 5.45)

trigger = ROOT.RooCategory("trigger", "trigger")
trigger.defineType("failed", 0)
trigger.defineType("passed", 1)

# WARNING: variable ranges are effectively cuts!
# Option: mass.setRange("signal",4.,6.);

variables = [
    ( RooRealVar("trigger",             "", 0),                       None,                   None   ), 
    ( RooRealVar("mm_kin_pt",           "p_{T}(#mu#mu)",              0, 100), Binning(40, 0, 40),     "Right" ),
    ( RooRealVar("mm_mu1_pt",           "p_{T}(#mu_{1})",             0,  50), Binning(25, 0, 25),     "Right" ),
    ( RooRealVar("mm_mu2_pt",           "p_{T}(#mu_{2})",             0,  50), Binning(25, 0, 25),     "Right" ),
    ( RooRealVar("mm_kin_eta",          "#eta(#mu#mu)",            -1.5, 1.5), Binning(60, -1.5, 1.5), "TopRight" ),
    ( RooRealVar("mm_kin_alpha",        "#alpha_{3D}",                0, 0.4), Binning(50, 0, 0.1),    "Right" ),
    ( RooRealVar("mm_kin_alphaSig",     "#alpha_{3D}/#sigma_{#alpha}", 0, 10), Binning(40, 0, 4),      "Right" ),
    ( RooRealVar("mm_kin_alphaBS",      "#alpha_{BS}",                0, 0.4), Binning(50, 0, 0.1),    "Right" ),
    ( RooRealVar("mm_kin_alphaBSSig",   "#alpha_{BS}/#sigma_{#alpha}", 0, 10), Binning(40, 0, 4),      "Right" ),
    ( RooRealVar("mm_kin_spvip",        "#delta_{3D}/#sigma_{#delta}", 0, 10), Binning(40, 0, 4),      "Right" ),
    ( RooRealVar("mm_kin_pvip",         "#delta_{3D}",               0, 0.02), Binning(50, 0, 0.02),   "Right" ),
    ( RooRealVar("mm_iso",              "iso(#mu#mu)",               0,  1.1), Binning(35, 0.4, 1.1),   "Left" ),
    ( RooRealVar("mm_m1iso",            "iso(#mu_{1})",               0, 1.1), Binning(35, 0.4, 1.1),   "Left" ),
    ( RooRealVar("mm_m2iso",            "iso(#mu_{2})",               0, 1.1), Binning(35, 0.4, 1.1),   "Left" ),
    ( RooRealVar("mm_kin_sl3d",         "L_{3D}/#sigma_{L}",          0, 100), Binning(50, 0, 100),    "Right" ),
    ( RooRealVar("mm_mva",              "MVA",                      0., 1.01), Binning(50, 0.9, 1.0),   "Left" ),
    ( RooRealVar("mm_kin_vtx_chi2dof",  "#chi^{2}/dof",                 0, 5), Binning(50, 0, 5),      "Right" ),
    # ( RooRealVar("mm_kin_vtx_prob",     "vtx prob",              0.025, 1.01), Binning(50, 0.025, 1),  "Left" ),
    ( RooRealVar("mm_nBMTrks",          "N_{trk}",                     0, 10), Binning(10, 0, 10),     "Right" ),
    ( RooRealVar("mm_otherVtxMaxProb1", "Other vetex probability 1", -0.01,1.01), Binning(51, -0.01,1.01), "Right" ),
    ( RooRealVar("mm_otherVtxMaxProb2", "Other vetex probability 2", -0.01,1.01), Binning(51, -0.01,1.01), "Right" ),
    ( RooRealVar("evt_nvtx",            "N_{PV}",                     0, 100), Binning(70, 0, 70),      "Right" ),
]

roovars = []
for var in variables:
    roovars.append(var[0])

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/%u/bmm_mva_jpsik/" % version
mc_2018 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
mc_2017 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root"
mc_2016BF = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root"
mc_2016GH = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
    
datasets = {

    "Run2018":{
        "mc":   mc_2018,
        "data": [
            data_path + "Charmonium+Run2018A-12Nov2019_UL2018_rsb-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2018B-12Nov2019_UL2018-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2018C-12Nov2019_UL2018_rsb_v3-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2018D-12Nov2019_UL2018-v1+MINIAOD/*.root"
        ]
    },
    "Run2017":{
        "mc":   mc_2017,
        "data": [
            data_path + "Charmonium+Run2017B-09Aug2019_UL2017-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2017C-09Aug2019_UL2017-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2017D-09Aug2019_UL2017-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2017E-09Aug2019_UL2017-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2017F-09Aug2019_UL2017-v1+MINIAOD/*.root"
        ]
    },
    "Run2016BF":{
        "mc":   mc_2016BF,
        "data": [
            data_path + "Charmonium+Run2016B-21Feb2020_ver2_UL2016_HIPM-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2016C-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2016D-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2016E-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2016F-21Feb2020_UL2016_HIPM-v1+MINIAOD/*.root",
        ]
    },
    "Run2016GH":{
        "mc":   mc_2016GH,
        "data": [
            data_path + "Charmonium+Run2016G-21Feb2020_UL2016-v1+MINIAOD/*.root",
            data_path + "Charmonium+Run2016H-21Feb2020_UL2016-v1+MINIAOD/*.root",
        ]
    },
}

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))


ROOT.gROOT.SetBatch(True)
# c2 = ROOT.TCanvas("c2","c2", 800, 800)
# c2_1 = ROOT.TPad("c2_1", "", 0.0, 0.25, 1.0, 1.0)
# c2_2 = ROOT.TPad("c2_2", "", 0.0, 0.0, 1.0, 0.25)
c1 = ROOT.TCanvas("c1","c1", 800, 800)

for dataset, info in datasets.items():

    name = dataset
    selection = "trigger>0"
    
    chain_mc = ROOT.TChain("mva")
    chain_mc.Add(info['mc'])
    print chain_mc.GetEntries()
    ds_mc = sweights.make_dataset(chain_mc, "mc", mass, roovars, selection)

    # for i in range(10):
    #     ds_mc.get(i).Print("V")
    # break

    chain_data = ROOT.TChain("mva")
    for pattern in info['data']:
        chain_data.Add(pattern)

    # get reference signal mass distribution
    ref_hist = ROOT.TH1F("ref_hist", "", 120, 5.15, 5.45)
    chain_mc.Draw("mm_kin_mass>>ref_hist", selection)

    # ws = sweights.get_workspace_with_weights_for_jpsik(chain_data, "data", mass, roovars, selection)
    ws = sweights.get_workspace_with_weights_for_jpsik(chain_data, "data", mass, roovars,
                                                       cuts=selection, ref_hist=ref_hist)

    ### Fit Validation
    
    ds_data = ws.data("data")
    model = ws.pdf("model")
    # mass = ws.var("mm_kin_mass")

    frame = mass.frame()
    ds_data.plotOn(frame)
    model.plotOn(frame)
    # # model.plotOn(frame, ROOT.RooFit.Components(bkg), ROOT.RooFit.LineStyle(ROOT.kDashed))

    model.paramOn(frame, ROOT.RooFit.Layout(0.6, 0.85, 0.85))
    # # model->paramOn(frame, Layout(0.55));
    frame.getAttText().SetTextSize(0.02)
    frame.Draw()
    print_canvas(name + "_mass_fit", output_path)

    ### Plots results

    dataw_sig = ROOT.RooDataSet(ds_data.GetName(), ds_data.GetTitle(), ds_data, ds_data.get(), "", "Nsig_sw")

    ROOT.gStyle.SetOptFit(1)
    
    for var, binning, legend_position in variables:
        c1.cd()
        if var.GetName()=="trigger":
            continue
        # nbins = 50
        # if v.GetName() == 'mm_mva':
        #     nbins = 11
        h_data = ROOT.RooAbsData.createHistogram(dataw_sig, 'sig', var, binning)
        h_data.SetLineWidth(2)
        # h_data.SetMarkerStyle(21)
        if var.GetName() == 'mm_kin_alphaSig':
            h_data.Fit("gaus")
        h_data.Draw("e")
        # print_canvas(name + "_splot_" + var.GetName(), output_path)
        print_canvas(var.GetName() + "_splot_" + name, output_path)
        if var.GetName() == 'mm_kin_alphaSig':
            h_data.GetListOfFunctions().Clear()
            
        h_mc = ROOT.RooAbsData.createHistogram( ds_mc, 'mc', var, binning)
        h_mc.SetMarkerColor(ROOT.kRed)
        h_mc.SetMarkerStyle(20)
        if var.GetName() == 'mm_kin_alphaSig':
            h_mc.Fit("gaus")
        h_mc.Draw()
        # print_canvas(name + "_mc_" + v.GetName(), output_path)
        print_canvas(var.GetName() + "_mc_" + name, output_path)
        if var.GetName() == 'mm_kin_alphaSig':
            h_mc.GetListOfFunctions().Clear()

        max_value = h_data.GetMaximum()
        h_mc.Scale(h_data.Integral()/h_mc.Integral())
        if h_mc.GetMaximum() > max_value:
            max_value = h_mc.GetMaximum()
        if re.search("Top", legend_position):
            max_value *= 1.4
        h_data.SetMaximum(max_value * 1.1)
        h_data.SetMinimum(0)
        h_mc.SetMaximum(max_value * 1.1)
        h_data.Draw("hist")
        h_mc.Draw("e same")

        # right handside legend
        if re.search("Right", legend_position):
            legend = ROOT.TLegend(0.60,0.75,0.85,0.87)
        else:
            legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
        legend.SetFillStyle(0)
        legend.SetLineWidth(0)
        legend.SetBorderSize(1)

        legend.AddEntry(h_data, "Data")
        legend.AddEntry(h_mc, "MC")
        legend.Draw()

        # print_canvas(name + "_all_" + v.GetName(), output_path)
        print_canvas(var.GetName() + "_all_" + name, output_path)

        ratio_plot = ROOT.TRatioPlot(h_data, h_mc)
        # ratio_plot.SetH1DrawOpt("hist e")
        ratio_plot.Draw()
        ratio_plot.SetSeparationMargin(0.03)
        ratio_plot.GetLowerRefGraph().SetMinimum(0.4)
        ratio_plot.GetLowerRefGraph().SetMaximum(1.6)
        # ratio_plot.GetXaxis().SetTitleSize()
        # SetBottomMargin(2.0)
        # c1.SetBottomMargin(2.0)
        # rp->GetLowerRefYaxis()->SetRange(...)
        # rp->SetH1DrawOpt("E");
        ratio_plot.GetLowerRefYaxis().SetTitle("Data/MC")
        
        ratio_plot.GetLowerRefYaxis().SetTitleSize()
        ratio_plot.GetLowerRefYaxis().SetTitleOffset(1.1)
        ratio_plot.GetLowerRefYaxis().SetLabelSize(0.035)
        ratio_plot.GetLowYaxis().SetNdivisions(503)
        
        ratio_plot.GetLowerRefXaxis().SetTitleSize()
        ratio_plot.GetLowerRefXaxis().SetTitleOffset()
        ratio_plot.GetLowerRefXaxis().SetLabelSize(0.035)
        
        ratio_plot.GetUpperRefYaxis().SetTitle("")
        ratio_plot.GetUpperRefYaxis().SetTitleSize()
        ratio_plot.GetUpperRefYaxis().SetTitleOffset()
        ratio_plot.GetUpperRefYaxis().SetLabelSize(0.035)
        
        ratio_plot.GetUpperRefXaxis().SetTitleSize()
        ratio_plot.GetUpperRefXaxis().SetTitleOffset(0)
        ratio_plot.GetUpperRefXaxis().SetLabelSize(0.035)

        c1.Update()
        c1.cd()
        legend.Draw()
        print_canvas(var.GetName() + "_ratio_" + name, output_path)

        
        # c2_1.cd()
        # h_data.Draw()
        # h_mc.Draw("same")
        # legend.Draw()
        
        # c2_2.cd()
        # h_ratio = h_data.Clone("h_ratio")
        # h_ratio.Divide(h_data, h_mc)
        # h_ratio.SetMinimum(0.5)
        # h_ratio.SetMaximum(1.5)
        # h_ratio.Draw()
        # c2.cd()
        # c2_1.Draw()
        # c2_2.Draw()

        # print_canvas(name + "_ratio_" + v.GetName(), output_path)
        print_canvas(var.GetName() + "_ratio_" + name, output_path)

        
        
        
    ws.Delete() # Cleanup
