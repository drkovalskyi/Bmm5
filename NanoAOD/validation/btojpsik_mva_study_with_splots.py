import os
import ROOT
import math 
import sweights
import tdrstyle

# Set the TDR style
tdrstyle.setTDRStyle()

version = 514

output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-%u/bjpsik-splots/" % version;

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

variables = [
    ROOT.RooRealVar("mm_kin_pt", "", 0, 30),
    ROOT.RooRealVar("mm_kin_eta", "", -1.5, 1.5),
    ROOT.RooRealVar("mm_kin_alpha", "", 0, 0.05),
    ROOT.RooRealVar("mm_kin_alphaSig", "", 0, 10),
    ROOT.RooRealVar("mm_kin_alphaBS", "", 0, 0.05),
    ROOT.RooRealVar("mm_kin_alphaBSSig", "", 0, 10),
    ROOT.RooRealVar("mm_kin_spvip", "", 0, 10),
    ROOT.RooRealVar("mm_kin_pvip", "", 0, 0.02),
    ROOT.RooRealVar("mm_iso", "", 0, 1),
    ROOT.RooRealVar("mm_m1iso", "", 0, 1),
    ROOT.RooRealVar("mm_m2iso", "", 0, 1),
    ROOT.RooRealVar("mm_kin_sl3d", "", 0, 100),
    ROOT.RooRealVar("mm_mva", "", 0.9, 1.01),
    ROOT.RooRealVar("mm_kin_vtx_chi2dof", "", 0, 5),
    ROOT.RooRealVar("mm_nBMTrks", "", 0, 10),
    ROOT.RooRealVar("mm_otherVtxMaxProb1", "", -0.01, 1.01),
    ROOT.RooRealVar("mm_otherVtxMaxProb2", "", -0.01, 1.01),
    ROOT.RooRealVar("trigger","", 0),

]

# struct dataset{
#   string name, mc_files, data_files, trigger;
# };

data_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/%u/bmm_mva_jpsik/" % version
# mc_2018 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root"
# mc_2017 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3+MINIAODSIM/*.root"
# mc_2016 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen+RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2+MINIAODSIM/*.root"
mc_2018 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2+MINIAODSIM/*.root"
mc_2017 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL17MiniAOD-106X_mc2017_realistic_v6-v2+MINIAODSIM/*.root"
mc_2016BF = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1+MINIAODSIM/*.root"
mc_2016GH = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM/*.root"
    
datasets = {
    # "Run2018A":{
    #     "mc":   mc_2018,
    #     "data": data_path + "Charmonium+Run2018A-17Sep2018-v1+MINIAOD/*.root"
    # },
    # "Run2018B":{
    #     "mc":   mc_2018,
    #     "data": data_path + "Charmonium+Run2018B-17Sep2018-v1+MINIAOD/*.root"
    # },
    # "Run2018C":{
    #     "mc":   mc_2018,
    #     "data": data_path + "Charmonium+Run2018C-17Sep2018-v1+MINIAOD/*.root"
    # },
    # "Run2018D":{
    #     "mc":   mc_2018,
    #     "data": data_path + "Charmonium+Run2018D-PromptReco-v2+MINIAOD/*.root"
    # },
    # "Run2018":{
    #     "mc":   mc_2018,
    #     "data": [
    #         data_path + "Charmonium+Run2018A-17Sep2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2018B-17Sep2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2018C-17Sep2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2018D-PromptReco-v2+MINIAOD/*.root",
    #     ]
    # },
    # "Run2017":{
    #     "mc":   mc_2017,
    #     "data": [
    #         data_path + "Charmonium+Run2017B-31Mar2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2017C-31Mar2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2017D-31Mar2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2017E-31Mar2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2017F-31Mar2018-v1+MINIAOD/*.root"
    #     ]
    # },
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
    # "Run2016":{
    #     "mc":   mc_2016,
    #     "data": [
    #         data_path + "Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016C-17Jul2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016D-17Jul2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016E-17Jul2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016F-17Jul2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016G-17Jul2018-v1+MINIAOD/*.root",
    #         data_path + "Charmonium+Run2016H-17Jul2018-v1+MINIAOD/*.root"
    #     ]
    # },
    # "Run2017B":{
    #     "mc":   mc_2017,
    #     "data": data_path + "Charmonium+Run2017B-31Mar2018-v1+MINIAOD/*.root"
    # },
    # "Run2017C":{
    #     "mc":   mc_2017
    #     "data": data_path + "Charmonium+Run2017C-31Mar2018-v1+MINIAOD/*.root"
    # },
    # "Run2017D":{
    #     "mc":   mc_2017,
    #     "data": data_path + "Charmonium+Run2017D-31Mar2018-v1+MINIAOD/*.root"
    # },
    # "Run2017E":{
    #     "mc":   mc_2017,
    #     "data": data_path + "Charmonium+Run2017E-31Mar2018-v1+MINIAOD/*.root"
    # },
    # "Run2017F":{
    #     "mc":   mc_2017,
    #     "data": data_path + "Charmonium+Run2017F-31Mar2018-v1+MINIAOD/*.root"
    # },
    # "Run2016B":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD/*.root"
    # },
    # "Run2016C":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016C-17Jul2018-v1+MINIAOD/*.root"
    # },
    # "Run2016D":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016D-17Jul2018-v1+MINIAOD/*.root"
    # },
    # "Run2016E":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016E-17Jul2018-v1+MINIAOD/*.root"
    # }
    # "Run2016F":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016F-17Jul2018-v1+MINIAOD/*.root"
    # },
    # "Run2016G":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016G-17Jul2018-v1+MINIAOD/*.root"
    # },
    # "Run2016H":{
    #     "mc":   mc_2016,
    #     "data": data_path + "Charmonium+Run2016H-17Jul2018-v1+MINIAOD/*.root"
    # }
}

def print_canvas(output_name_without_extention, path, canvas=ROOT.gPad):
    if not os.path.exists(path):
        os.makedirs(path)
    canvas.Print("%s/%s.png" % (path, output_name_without_extention))
    canvas.Print("%s/%s.pdf" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.root" % (path, output_name_without_extention))
    # canvas.Print("%s/%s.C" % (path, output_name_without_extention))

# void add_data(TFile* f, TChain* mc_bkmm, TChain* data, const char* trigger){

#   TCut mc_match;
#   if (use_mc_truth_matching)
#     // mc_match = "bkmm_gen_pdgId!=0";
#     mc_match = "abs(1-bkmm_kaon_pt/genbmm_kaon1_pt[0])<0.1";

#   TCut base_cut(trigger);
#   // TCut base_cut("HLT_DoubleMu4_Jpsi_NoVertexing");
#   base_cut += "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4";
#   base_cut += "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4";
#   base_cut += "mm_kin_sl3d[bkmm_mm_index]>4 && mm_kin_vtx_chi2dof[bkmm_mm_index]<5";
#   base_cut += "abs(bkmm_jpsimc_mass-5.4)<0.5 && bkmm_jpsimc_vtx_chi2dof<5";
#   base_cut += "bkmm_kaon_pt<1.5";
#   base_cut += "bkmm_jpsimc_alpha<0.2"; 
#   base_cut += "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45";
#   base_cut += "!TMath::IsNaN(bkmm_jpsimc_massErr)";

#   // Monte Carlo
  
#   TH1D* h_mc_bkmm_base = new TH1D("h_mc_bkmm_base", "", 50, 5.17, 5.45);
#   mc_bkmm->Draw("bkmm_jpsimc_mass>>h_mc_bkmm_base", base_cut + mc_match, "goff");
#   h_mc_bkmm_base->Write();

#   for (auto const& var : variables){
#     string h_mass_name = "h_mc_bkmm_" + var.name + "_vs_mass";
#     string h_npv_name  = "h_mc_bkmm_" + var.name + "_vs_npv";
#     string command_mass = var.branch + ":bkmm_jpsimc_mass>>" + h_mass_name;
#     string command_npv  = var.branch + ":PV_npvs>>" + h_npv_name;
#     TH2D* h_mass = new TH2D(h_mass_name.c_str(), "", 50, 5.17, 5.45, var.nbins, var.xmin, var.xmax);
#     TH2D* h_npv  = new TH2D(h_npv_name.c_str(),  "", 50,    0,  100, var.nbins, var.xmin, var.xmax);
#     mc_bkmm->Draw(command_mass.c_str(), base_cut + mc_match, "goff");
#     mc_bkmm->Draw(command_npv.c_str(),  base_cut + mc_match, "goff");
#     h_mass->Write();
#     h_npv->Write();
#   }

#   // Data

#   TH1D* h_data_base = new TH1D("h_data_base", "", 50, 5.17, 5.45);
#   data->Draw("bkmm_jpsimc_mass>>h_data_base", base_cut, "goff");
#   h_data_base->Write();

#   for (auto const& var : variables){
#     string hist_name = "h_data_" + var.name + "_vs_mass";
#     string draw_command = var.branch + ":bkmm_jpsimc_mass>>" + hist_name;
#     TH2D* h_data = new TH2D(hist_name.c_str(), "", 
# 			    50, 5.17, 5.45, 
# 			    var.nbins, var.xmin, var.xmax);
#     data->Draw(draw_command.c_str(), base_cut, "goff");
#     h_data->Write();
#   }

#   TH1D* h_data_npv = new TH1D("h_data_npv", "", 50, 0, 100);
#   data->Draw("PV_npvs>>h_data_npv", base_cut, "goff");
#   h_data_npv->Write();

# }

# struct Result {
#   double sig;
#   double sig_err;
#   double bkg;
#   double bkg_err;
#   Result(const RooRealVar& nsig, const RooRealVar& nbkg):
#     sig(nsig.getVal()), sig_err(nsig.getError()),
#     bkg(nbkg.getVal()), bkg_err(nbkg.getError()){}
#   Result():
#     sig(0), sig_err(0),
#     bkg(0), bkg_err(0){}
# };

# void process_variables(TFile* f, string prefix){
#   TH2* h_mc_bkmm_mva_vs_npv = (TH2*)f->Get("h_mc_bkmm_mva_vs_npv");
#   TH1* h_mc_npv = h_mc_bkmm_mva_vs_npv->ProjectionX();
#   TH1* h_data_npv = ((TH1*)f->Get("h_data_npv"));
#   h_mc_npv->SetLineWidth(2);
#   h_mc_npv->SetLineColor(kBlue);
#   h_mc_npv->Draw();
#   print_canvas(prefix + "-" + "mc_npv", output_path, gPad);

#   h_data_npv->SetLineWidth(2);
#   h_data_npv->SetLineColor(kBlue);
#   h_data_npv->Draw();
#   print_canvas(prefix + "-" + "data_npv", output_path, gPad);

#   for (auto const& var: variables){
#     unsigned int n = var.nbins;
#     vector<Double_t> mc_eff;
#     vector<Double_t> data_eff;
#     vector<Double_t> data_over_mc_eff;
#     vector<Double_t> mc_eff_reweighted;
#     vector<Double_t> data_over_mc_eff_reweighted;
#     double total_data(-1);
#     double total_mc(-1);
#     // vector<Double_t> mc_eff_err;
#     // vector<Double_t> data_over_mc_eff_err;

#     printf("processing %s\n", var.name.c_str());
#     // Scan efficiency in MC and data as a function of the cut on the variable

#     auto h_mc_x_vs_npv    = (TH2*)f->Get(("h_mc_bkmm_"   + var.name + "_vs_npv" ).c_str()); 
#     auto h_mc_x_vs_mass   = (TH2*)f->Get(("h_mc_bkmm_"   + var.name + "_vs_mass").c_str()); 
#     auto h_data_x_vs_mass = (TH2*)f->Get(("h_data_" + var.name + "_vs_mass").c_str());
#     auto h_mc_x_reweighted = reweight_histogram(h_mc_x_vs_npv, h_mc_npv, h_data_npv);
    
#     for (unsigned int i=0; i <= n; ++i){
#       auto mass_mc   = h_mc_x_vs_mass->ProjectionX(  "mass_mc",   i, n+1);
#       auto mass_data = h_data_x_vs_mass->ProjectionX("mass_data", i, n+1);
#       auto result = fitHistogram(mass_mc, mass_data);
      
#       if (store_projections)
# 	print_canvas(prefix + "-" + var.name + "_proj" + to_string(i), output_path, gPad);

#       if (i==0){
# 	total_data = result.sig;
# 	total_mc = mass_mc->Integral();
#       }

#       assert(total_data > 0);
#       data_eff.push_back(result.sig / total_data);
#       assert(total_mc>0);
#       mc_eff.push_back(mass_mc->Integral() / total_mc);
#       data_over_mc_eff.push_back(mc_eff.back()>0?data_eff.back()/mc_eff.back():0);

#       mc_eff_reweighted.push_back(h_mc_x_reweighted->Integral(i,-1) / h_mc_x_reweighted->Integral(0,-1));
#       data_over_mc_eff_reweighted.push_back(mc_eff_reweighted.back()>0?data_eff.back()/mc_eff_reweighted.back():0);
#     }
#     gPad->SetGridx();
#     gPad->SetGridy();
#     auto gr = new TGraph(n, &mc_eff[0], &data_over_mc_eff[0]);
#     gr->SetMinimum(0);
#     gr->SetMaximum(1.5);
#     gr->GetXaxis()->SetLimits(0.,1.);
#     gr->GetXaxis()->SetTitle("MC efficiency");
#     gr->GetYaxis()->SetTitle("Data/MC efficiency ratio");
#     gr->SetTitle(var.name.c_str());
#     gr->Draw("AP*");

#     print_canvas(prefix + "-" + var.name + "_eff", output_path, gPad);

#     auto gr_reweighted = new TGraph(n, &mc_eff_reweighted[0], &data_over_mc_eff_reweighted[0]);
#     gr_reweighted->SetMinimum(0);
#     gr_reweighted->SetMaximum(1.5);
#     gr_reweighted->GetXaxis()->SetLimits(0.,1.);
#     gr_reweighted->GetXaxis()->SetTitle("MC efficiency");
#     gr_reweighted->GetYaxis()->SetTitle("Data/MC efficiency ratio");
#     gr_reweighted->SetTitle(var.name.c_str());
#     gr_reweighted->Draw("AP*");
#     print_canvas(prefix + "-" + var.name + "_eff_reweighted", output_path, gPad);
#     gPad->SetGridx(0);
#     gPad->SetGridy(0);
#   }
# }

ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas("c1","c1", 800, 600)

for dataset, info in datasets.items():

    name = dataset
    
    chain_mc = ROOT.TChain("mva")
    chain_mc.Add(info['mc'])
    print chain_mc.GetEntries()
    ds_mc = sweights.make_dataset(chain_mc, "mc", mass, variables, "trigger>0")
    # for i in range(10):
    #     ds_mc.get(i).Print("V")
    # break

    chain_data = ROOT.TChain("mva")
    for pattern in info['data']:
        chain_data.Add(pattern)
    ws = sweights.get_workspace_with_weights_for_jpsik(chain_data, "data", mass, variables, "trigger>0")

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

    for v in variables:
        if v.GetName()=="trigger":
            continue
        h_data = ROOT.RooAbsData.createHistogram( dataw_sig, 'sig', v, ROOT.RooFit.Binning(50))
        h_data.Draw()
        print_canvas(name + "_splot_" + v.GetName(), output_path)

        h_mc = ROOT.RooAbsData.createHistogram( ds_mc, 'mc', v, ROOT.RooFit.Binning(50))
        h_mc.SetMarkerColor(ROOT.kRed)
        h_mc.Draw()
        print_canvas(name + "_mc_" + v.GetName(), output_path)

        max_value = h_data.GetMaximum()
        h_mc.Scale(h_data.Integral()/h_mc.Integral())
        if h_mc.GetMaximum() > max_value:
            max_value = h_mc.GetMaximum()
        h_data.SetMaximum(max_value * 1.1)
        h_mc.SetMaximum(max_value * 1.1)
        h_data.Draw()
        h_mc.Draw("same")

        legend = ROOT.TLegend(0.15,0.75,0.5,0.87)
        legend.SetFillStyle(0)
        legend.SetLineWidth(0)
        legend.SetBorderSize(1)

        legend.AddEntry(h_data, "Data")
        legend.AddEntry(h_mc, "MC")
        legend.Draw()

        print_canvas(name + "_all_" + v.GetName(), output_path)

        h_data.Divide(h_data, h_mc)
        h_data.SetMinimum(0.5)
        h_data.SetMaximum(1.5)
        h_data.Draw()
        print_canvas(name + "_ratio_" + v.GetName(), output_path)
    
    ws.Delete() # Cleanup
