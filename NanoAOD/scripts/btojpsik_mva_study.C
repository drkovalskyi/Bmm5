#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooNovosibirsk.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooGaussModel.h"
#include "RooHistPdf.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include <tuple>

using namespace RooFit;
using namespace std;

string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-507/bjpsik-fit/";

bool silent_roofit = true;
bool store_projections = false;
bool use_mc_truth_matching = true;

struct variable{
  string name, branch;
  unsigned int nbins;
  double xmin, xmax;
};

// WARNING: signal purity must increase as x increases from xmin to xmax
vector<variable> variables {
  {"mva", "bkmm_bmm_mva", 100, 0, 1},
  {"alpha", "(0.2-bkmm_jpsimc_alpha)", 100, 0, 0.2},
  {"alphaXY", "bkmm_jpsimc_cosAlphaXY", 100, 0.999, 1},
  {"spvip", "(10-bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr)", 100, 0, 10},
  {"pvip", "(0.02-bkmm_jpsimc_pvip)", 100, 0, 0.02},
  {"iso", "bkmm_bmm_iso", 100, 0, 1},
  {"m1iso", "bkmm_bmm_m1iso", 100, 0, 1},
  {"m2iso", "bkmm_bmm_m2iso", 100, 0, 1},
  {"sl3d", "mm_kin_sl3d[bkmm_mm_index]", 100, 0, 100},
  {"bdt", "bkmm_bmm_bdt", 50, -1, 1},
};

struct dataset{
  string name, mc_files, data_files, trigger;
};

string data_path = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/";
string mc_2018 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root";
string mc_2017 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3+MINIAODSIM/*.root";
string mc_2016 = data_path + "BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen+RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1+MINIAODSIM/*.root";
string trigger_2018 = "HLT_DoubleMu4_3_Jpsi";
string trigger_2017 = "HLT_DoubleMu4_Jpsi_Displaced";
// string trigger_2017 = "HLT_DoubleMu4_Jpsi_NoVertexing";
string trigger_2016 = "HLT_DoubleMu4_3_Jpsi_Displaced";
// string trigger_2016 = "HLT_Dimuon6_Jpsi_NoVertexing";
    
vector<dataset> datasets{
  {"Run2018A", mc_2018,	data_path + "Charmonium+Run2018A-17Sep2018-v1+MINIAOD/*.root", trigger_2018},
  {"Run2018B", mc_2018, data_path + "Charmonium+Run2018B-17Sep2018-v1+MINIAOD/*.root", trigger_2018},
  {"Run2018C", mc_2018,	data_path + "Charmonium+Run2018C-17Sep2018-v1+MINIAOD/*.root", trigger_2018},
  {"Run2018D", mc_2018,	data_path + "Charmonium+Run2018D-PromptReco-v2+MINIAOD/*.root",	trigger_2018},

  {"Run2017B", mc_2017,	data_path + "Charmonium+Run2017B-31Mar2018-v1+MINIAOD/*.root", trigger_2017},
  {"Run2017C", mc_2017,	data_path + "Charmonium+Run2017C-31Mar2018-v1+MINIAOD/*.root", trigger_2017},
  {"Run2017D", mc_2017,	data_path + "Charmonium+Run2017D-31Mar2018-v1+MINIAOD/*.root", trigger_2017},
  {"Run2017E", mc_2017, data_path + "Charmonium+Run2017E-31Mar2018-v1+MINIAOD/*.root", trigger_2017},
  {"Run2017F", mc_2017,	data_path + "Charmonium+Run2017F-31Mar2018-v1+MINIAOD/*.root", trigger_2017},

  {"Run2016B", mc_2016, data_path + "Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD/*.root", trigger_2016},
  {"Run2016C", mc_2016,	data_path + "Charmonium+Run2016C-17Jul2018-v1+MINIAOD/*.root", trigger_2016},
  {"Run2016D", mc_2016,	data_path + "Charmonium+Run2016D-17Jul2018-v1+MINIAOD/*.root", trigger_2016},
  {"Run2016E", mc_2016, data_path + "Charmonium+Run2016E-17Jul2018-v1+MINIAOD/*.root", trigger_2016},
  {"Run2016F", mc_2016,	data_path + "Charmonium+Run2016F-17Jul2018-v1+MINIAOD/*.root", trigger_2016},
  {"Run2016G", mc_2016,	data_path + "Charmonium+Run2016G-17Jul2018-v1+MINIAOD/*.root", trigger_2016},
  {"Run2016H", mc_2016,	data_path + "Charmonium+Run2016H-17Jul2018-v1+MINIAOD/*.root", trigger_2016},
};

void print_canvas(string output_name_without_extention, 
		  string path, 
		  TVirtualPad* canvas){
  if (gSystem->AccessPathName(path.c_str())){
    gSystem->Exec(("mkdir -p " + path).c_str());
  }
  string s = path + "/" + output_name_without_extention;
  canvas->Print((s + ".png").c_str());
  canvas->Print((s + ".pdf").c_str());
  canvas->Print((s + ".root").c_str());
}

void add_data(TFile* f, TChain* mc_bkmm, TChain* data, const char* trigger){

  TCut mc_match;
  if (use_mc_truth_matching)
    // mc_match = "bkmm_gen_pdgId!=0";
    mc_match = "abs(1-bkmm_kaon_pt/genbmm_kaon1_pt[0])<0.1";

  TCut base_cut(trigger);
  // TCut base_cut("HLT_DoubleMu4_Jpsi_NoVertexing");
  base_cut += "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4";
  base_cut += "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4";
  base_cut += "mm_kin_sl3d[bkmm_mm_index]>4 && mm_kin_vtx_chi2dof[bkmm_mm_index]<5";
  base_cut += "abs(bkmm_jpsimc_mass-5.4)<0.5 && bkmm_jpsimc_vtx_chi2dof<5";
  base_cut += "bkmm_kaon_pt<1.5";
  base_cut += "bkmm_jpsimc_alpha<0.2"; 
  base_cut += "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45";
  base_cut += "!TMath::IsNaN(bkmm_jpsimc_massErr)";

  // Monte Carlo
  
  TH1D* h_mc_bkmm_base = new TH1D("h_mc_bkmm_base", "", 50, 5.17, 5.45);
  mc_bkmm->Draw("bkmm_jpsimc_mass>>h_mc_bkmm_base", base_cut + mc_match, "goff");
  h_mc_bkmm_base->Write();

  for (auto const& var : variables){
    string h_mass_name = "h_mc_bkmm_" + var.name + "_vs_mass";
    string h_npv_name  = "h_mc_bkmm_" + var.name + "_vs_npv";
    string command_mass = var.branch + ":bkmm_jpsimc_mass>>" + h_mass_name;
    string command_npv  = var.branch + ":PV_npvs>>" + h_npv_name;
    TH2D* h_mass = new TH2D(h_mass_name.c_str(), "", 50, 5.17, 5.45, var.nbins, var.xmin, var.xmax);
    TH2D* h_npv  = new TH2D(h_npv_name.c_str(),  "", 50,    0,  100, var.nbins, var.xmin, var.xmax);
    mc_bkmm->Draw(command_mass.c_str(), base_cut + mc_match, "goff");
    mc_bkmm->Draw(command_npv.c_str(),  base_cut + mc_match, "goff");
    h_mass->Write();
    h_npv->Write();
  }

  // Data

  TH1D* h_data_base = new TH1D("h_data_base", "", 50, 5.17, 5.45);
  data->Draw("bkmm_jpsimc_mass>>h_data_base", base_cut, "goff");
  h_data_base->Write();

  for (auto const& var : variables){
    string hist_name = "h_data_" + var.name + "_vs_mass";
    string draw_command = var.branch + ":bkmm_jpsimc_mass>>" + hist_name;
    TH2D* h_data = new TH2D(hist_name.c_str(), "", 
			    50, 5.17, 5.45, 
			    var.nbins, var.xmin, var.xmax);
    data->Draw(draw_command.c_str(), base_cut, "goff");
    h_data->Write();
  }

  TH1D* h_data_npv = new TH1D("h_data_npv", "", 50, 0, 100);
  data->Draw("PV_npvs>>h_data_npv", base_cut, "goff");
  h_data_npv->Write();

}

struct Result {
  double sig;
  double sig_err;
  double bkg;
  double bkg_err;
  Result(const RooRealVar& nsig, const RooRealVar& nbkg):
    sig(nsig.getVal()), sig_err(nsig.getError()),
    bkg(nbkg.getVal()), bkg_err(nbkg.getError()){}
  Result():
    sig(0), sig_err(0),
    bkg(0), bkg_err(0){}
};

TH1* reweight_histogram(TH2* h_mc_x_vs_npv, TH1* h_mc_npv, TH1* h_data_npv){
  /// Pileup reweighting
  /// - Input is a 2D histrogram of the target observable as a function of npv
  /// - Output is scaled to a unit area including overflow bins
  
  // TODO
  // - do proper error propagation
  
  assert(h_mc_npv->GetNbinsX() == h_data_npv->GetNbinsX());

  TH1* h_mc_x_reweighted = h_mc_x_vs_npv->ProjectionY("h_mc_mva_reweighted");
  h_mc_x_reweighted->SetDirectory(0);
  // h_mc_x_reweighted->Sumw2();

  for (int x_i=0; x_i <= h_mc_x_reweighted->GetNbinsX() + 1; ++x_i){
    double integral(0);
    for (int npv_i=0; npv_i <= h_data_npv->GetNbinsX() + 1; ++npv_i){
      if (h_mc_npv->GetBinContent(npv_i)>0)
    	integral += h_mc_x_vs_npv->GetBinContent(npv_i, x_i) / 
    	  h_mc_npv->GetBinContent(npv_i) * h_data_npv->GetBinContent(npv_i);
    }
    h_mc_x_reweighted->SetBinContent(x_i, integral);
  }
  h_mc_x_reweighted->Scale(1/h_mc_x_reweighted->Integral(0,-1));
  return h_mc_x_reweighted;
}


Result fitHistogram(TH1* h_ref, TH1* h_test){
  if (h_ref->Integral() < 1 or h_test->Integral() < 1) return Result();

  if (silent_roofit)
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooRealVar mass("mass", "BtoJpisK mass", 5.29);
  
  // signal pdf
  RooRealVar bias("bias", "bias", 0, -0.1, 0.1) ;
  RooRealVar sigma("sigma", "sigma", 0.0001, 0., 0.01);
  RooGaussModel gaussM("gaussM", "signal pdf", mass, bias, sigma) ;
  RooDataHist ref_data("ref_data", "", mass, h_ref);
  RooHistPdf ref_pdf("ref_pdf", "theoretical lineshape", mass, ref_data, 2);
  mass.setBins(10000, "fft");
  RooFFTConvPdf sig("sig", "smeared distribution", mass, ref_pdf, gaussM);

  // background pdf
  RooRealVar a0("a0", "a0", 0.0, -1., 1.) ;
  // RooRealVar a1("a1", "a1", 0.0, -0.2, 0.2) ;
  // RooRealVar a2("a2", "a2", 0.0, -1., 1.) ;
  // RooChebychev bkg("bkg", "Background", mass, RooArgSet(a0, a1));
  RooChebychev bkg("bkg", "Background", mass, RooArgSet(a0));
  RooRealVar Nsig("Nsig", "Nsig", h_test->Integral(), 0, h_test->Integral());
  RooRealVar Nbkg("Nbkg", "Nbkg", 0, 0, h_test->Integral());
  RooAddPdf model("model", "", RooArgList(sig,bkg), RooArgList(Nsig,Nbkg));

  // test dataset
  RooDataHist test_data("test_data", "", mass, h_test);

  // freeze signal shape to get background shape right 
  sigma.setConstant(true);
  bias.setConstant(true);
  if (silent_roofit)
    model.fitTo(test_data, PrintEvalErrors(-1), PrintLevel(-1));
  else
    model.fitTo(test_data);
    
  sigma.setConstant(false);
  bias.setConstant(false);
  if (silent_roofit)
    model.fitTo(test_data, PrintEvalErrors(-1), PrintLevel(-1));
  else
    model.fitTo(test_data);

  RooPlot* frame = mass.frame();
  test_data.plotOn(frame);
  model.plotOn(frame);
  model.plotOn(frame, RooFit::Components(bkg), RooFit::LineStyle(kDashed));
  // data.statOn(frame, Layout(0.55, 0.99, 0.8));
  // model->paramOn(frame, Parameters(params), Layout(0.6,0.9,0.9) ) ;

  model.paramOn(frame, Layout(0.6, 0.85, 0.85));
  // model->paramOn(frame, Layout(0.55));
  frame->getAttText()->SetTextSize(0.02);
  frame->Draw();

  return Result(Nsig,Nbkg);
}

void process_variables(TFile* f, string prefix){
  TH2* h_mc_bkmm_mva_vs_npv = (TH2*)f->Get("h_mc_bkmm_mva_vs_npv");
  TH1* h_mc_npv = h_mc_bkmm_mva_vs_npv->ProjectionX();
  TH1* h_data_npv = ((TH1*)f->Get("h_data_npv"));
  h_mc_npv->SetLineWidth(2);
  h_mc_npv->SetLineColor(kBlue);
  h_mc_npv->Draw();
  print_canvas(prefix + "-" + "mc_npv", output_path, gPad);

  h_data_npv->SetLineWidth(2);
  h_data_npv->SetLineColor(kBlue);
  h_data_npv->Draw();
  print_canvas(prefix + "-" + "data_npv", output_path, gPad);

  for (auto const& var: variables){
    unsigned int n = var.nbins;
    vector<Double_t> mc_eff;
    vector<Double_t> data_eff;
    vector<Double_t> data_over_mc_eff;
    vector<Double_t> mc_eff_reweighted;
    vector<Double_t> data_over_mc_eff_reweighted;
    double total_data(-1);
    double total_mc(-1);
    // vector<Double_t> mc_eff_err;
    // vector<Double_t> data_over_mc_eff_err;

    printf("processing %s\n", var.name.c_str());
    // Scan efficiency in MC and data as a function of the cut on the variable

    auto h_mc_x_vs_npv    = (TH2*)f->Get(("h_mc_bkmm_"   + var.name + "_vs_npv" ).c_str()); 
    auto h_mc_x_vs_mass   = (TH2*)f->Get(("h_mc_bkmm_"   + var.name + "_vs_mass").c_str()); 
    auto h_data_x_vs_mass = (TH2*)f->Get(("h_data_" + var.name + "_vs_mass").c_str());
    auto h_mc_x_reweighted = reweight_histogram(h_mc_x_vs_npv, h_mc_npv, h_data_npv);
    
    for (unsigned int i=0; i <= n; ++i){
      auto mass_mc   = h_mc_x_vs_mass->ProjectionX(  "mass_mc",   i, n+1);
      auto mass_data = h_data_x_vs_mass->ProjectionX("mass_data", i, n+1);
      auto result = fitHistogram(mass_mc, mass_data);
      
      if (store_projections)
	print_canvas(prefix + "-" + var.name + "_proj" + to_string(i), output_path, gPad);

      if (i==0){
	total_data = result.sig;
	total_mc = mass_mc->Integral();
      }

      assert(total_data > 0);
      data_eff.push_back(result.sig / total_data);
      assert(total_mc>0);
      mc_eff.push_back(mass_mc->Integral() / total_mc);
      data_over_mc_eff.push_back(mc_eff.back()>0?data_eff.back()/mc_eff.back():0);

      mc_eff_reweighted.push_back(h_mc_x_reweighted->Integral(i,-1) / h_mc_x_reweighted->Integral(0,-1));
      data_over_mc_eff_reweighted.push_back(mc_eff_reweighted.back()>0?data_eff.back()/mc_eff_reweighted.back():0);
    }
    gPad->SetGridx();
    gPad->SetGridy();
    auto gr = new TGraph(n, &mc_eff[0], &data_over_mc_eff[0]);
    gr->SetMinimum(0);
    gr->SetMaximum(1.5);
    gr->GetXaxis()->SetLimits(0.,1.);
    gr->GetXaxis()->SetTitle("MC efficiency");
    gr->GetYaxis()->SetTitle("Data/MC efficiency ratio");
    gr->SetTitle(var.name.c_str());
    gr->Draw("AP*");

    print_canvas(prefix + "-" + var.name + "_eff", output_path, gPad);

    auto gr_reweighted = new TGraph(n, &mc_eff_reweighted[0], &data_over_mc_eff_reweighted[0]);
    gr_reweighted->SetMinimum(0);
    gr_reweighted->SetMaximum(1.5);
    gr_reweighted->GetXaxis()->SetLimits(0.,1.);
    gr_reweighted->GetXaxis()->SetTitle("MC efficiency");
    gr_reweighted->GetYaxis()->SetTitle("Data/MC efficiency ratio");
    gr_reweighted->SetTitle(var.name.c_str());
    gr_reweighted->Draw("AP*");
    print_canvas(prefix + "-" + var.name + "_eff_reweighted", output_path, gPad);
    gPad->SetGridx(0);
    gPad->SetGridy(0);
  }
}

void btojpsik_mva_study(){
  for (auto const& ds: datasets){
    string file_name = "btojpsik_mva_study-" + ds.name + ".root";
    // Check if we already have histograms prepared
    if (gSystem->AccessPathName(file_name.c_str())){
      printf("Recreating histograms for %s\n", ds.name.c_str());
      TFile* f = TFile::Open(file_name.c_str(), "RECREATE");
      TChain* mc = new TChain("Events");
      mc->Add(ds.mc_files.c_str());
      TChain* data = new TChain("Events");
      data->Add(ds.data_files.c_str());
      add_data(f, mc, data, (ds.trigger.c_str()));    
      f->Close();
    }
    // continue;

    TFile* f = TFile::Open(file_name.c_str());
  
    TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
    auto result_all = fitHistogram((TH1*)f->Get("h_mc_bkmm_base"), 
				   (TH1*)f->Get("h_data_base"));
    
    print_canvas(ds.name + "base", output_path, c1);
      
    process_variables(f, ds.name);
  }

  return;

  // unsigned int n = 50;

  // // MVA
  // vector<tuple<Result, unsigned int>> results_mva;

  // for (unsigned int i=1; i <= n; ++i){
  //   auto proj_mc = ((TH2*)f->Get("h_mc_bkmm_mva"))->ProjectionX("proj_mc", i);
  //   auto proj_data = ((TH2*)f->Get("h_data_mva"))->ProjectionX("proj_data", i);
  //   auto result = fitHistogram(proj_mc, proj_data);
  //   results_mva.push_back(make_tuple(result, proj_mc->Integral()));
  //   print_canvas("mva_proj"+to_string(i), output_path, c1);
  // }

  // // BDT
  // vector<tuple<Result, unsigned int>> results_bdt;

  // for (unsigned int i=1; i <= n; ++i){
  //   auto proj_mc = ((TH2*)f->Get("h_mc_bkmm_bdt"))->ProjectionX("proj_mc", i);
  //   auto proj_data = ((TH2*)f->Get("h_data_bdt"))->ProjectionX("proj_data", i);
  //   auto result = fitHistogram(proj_mc, proj_data);
  //   results_bdt.push_back(make_tuple(result, proj_mc->Integral()));
  //   print_canvas("bdt_proj"+to_string(i), output_path, c1);
  // }


  // // Pileup treatment
  // TH1* h_mc_npv = ((TH2*)f->Get("h_mc_bkmm_mva_vs_npv"))->ProjectionX();
  // TH1* h_data_npv = ((TH1*)f->Get("h_data_npv"));
  // h_mc_npv->SetLineWidth(2);
  // h_mc_npv->SetLineColor(kBlue);
  // h_mc_npv->Draw();
  // print_canvas("mc_npv", output_path, c1);

  // h_data_npv->SetLineWidth(2);
  // h_data_npv->SetLineColor(kBlue);
  // h_data_npv->Draw();
  // print_canvas("data_npv", output_path, c1);
  
  // // Results
  // printf("Total: %0.1f+/-%0.1f\n", result_all.sig, result_all.sig_err);
  // printf("MVA:\n");
  // for (unsigned int i=0; i < n; ++i){
  //   double eff_data = get<0>(results_mva[i]).sig/get<0>(results_mva[0]).sig;
  //   double eff_mc = get<1>(results_mva[i])/float(get<1>(results_mva[0]));
  //   printf("\t[%u] \tmva_data: %0.1f+/-%0.1f \tmva_mc:%u \teff_data: %0.1f%% \teff_mc: %0.1f%% \teff_mc/eff_data: %0.2f\n", 
  // 	   i+1,
  // 	   get<0>(results_mva[i]).sig, 
  // 	   get<0>(results_mva[i]).sig_err, 
  // 	   get<1>(results_mva[i]),
  // 	   eff_data*100.,
  // 	   eff_mc*100.,
  // 	   eff_data>0?eff_mc/eff_data:0
  // 	   ); 
  // }
  // printf("BDT:\n");
  // for (unsigned int i=0; i < n; ++i){
  //   double eff_data = get<0>(results_bdt[i]).sig/get<0>(results_bdt[0]).sig;
  //   double eff_mc = get<1>(results_bdt[i])/float(get<1>(results_bdt[0]));
  //   printf("\t[%u] \tbdt_data: %0.1f+/-%0.1f \tbdt_mc:%u \teff_data: %0.1f%% \teff_mc: %0.1f%% \teff_mc/eff_data: %0.2f\n", 
  // 	   i+1,
  // 	   get<0>(results_bdt[i]).sig, 
  // 	   get<0>(results_bdt[i]).sig_err, 
  // 	   get<1>(results_bdt[i]),
  // 	   eff_data*100.,
  // 	   eff_mc*100.,
  // 	   eff_data>0?eff_mc/eff_data:0
  // 	   ); 
  // }
}
