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

std::string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-507/bjpsik-fit/";

// name, branch name, number of bins, xmin, xmax
// WARNING: signal purity must increase as x increases from xmin to xmax
std::vector<std::tuple<std::string, std::string, unsigned int, float, float> > variables;
void add_variables(){
  variables.push_back(std::make_tuple("mva", "bkmm_bmm_mva", 100, 0, 1));
  // variables.push_back(std::make_tuple("alpha", "(0.2-bkmm_jpsimc_alpha)", 100, 0, 0.2));
  // variables.push_back(std::make_tuple("alphaXY", "bkmm_jpsimc_cosAlphaXY", 100, 0.999, 1));
  // variables.push_back(std::make_tuple("spvip", "(10-bkmm_jpsimc_pvip/bkmm_jpsimc_pvipErr)", 100, 0, 10));
  // variables.push_back(std::make_tuple("pvip", "(0.02-bkmm_jpsimc_pvip)", 100, 0, 0.02));
  // variables.push_back(std::make_tuple("iso", "bkmm_bmm_iso", 100, 0, 1));
  // variables.push_back(std::make_tuple("m1iso", "bkmm_bmm_m1iso", 100, 0, 1));
  // variables.push_back(std::make_tuple("m2iso", "bkmm_bmm_m2iso", 100, 0, 1));
  // variables.push_back(std::make_tuple("sl3d", "mm_kin_sl3d[bkmm_mm_index]", 100, 0, 100));
  // variables.push_back(std::make_tuple("bdt", "bkmm_bmm_bdt", 50, -1, 1));
}

void print_canvas(std::string output_name_without_extention, 
		  std::string path, 
		  TVirtualPad* canvas){
  if (gSystem->AccessPathName(path.c_str())){
    gSystem->Exec(("mkdir -p " + path).c_str());
  }
  std::string s = path + "/" + output_name_without_extention;
  canvas->Print((s + ".png").c_str());
  canvas->Print((s + ".pdf").c_str());
  // canvas->Print((s + ".root").c_str());
}

void add_data(TFile* f){

  TCut base_cut("HLT_DoubleMu4_3_Jpsi");
  // TCut base_cut("HLT_DoubleMu4_Jpsi_NoVertexing");
  base_cut += "abs(Muon_eta[mm_mu1_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4";
  base_cut += "abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4";
  base_cut += "mm_kin_sl3d[bkmm_mm_index]>4 && mm_kin_vtx_chi2dof[bkmm_mm_index]<5";
  base_cut += "abs(bkmm_jpsimc_mass-5.4)<0.5 && bkmm_jpsimc_vtx_chi2dof<5";
  base_cut += "bkmm_kaon_pt>1";
  base_cut += "bkmm_jpsimc_alpha<0.2"; 
  base_cut += "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45";
  base_cut += "!TMath::IsNaN(bkmm_jpsimc_massErr)";

  TChain* mc_bkmm = new TChain("Events");
  // mc_bkmm->Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0bb34079aa40de8c3fd5615d4816875f.root");
  mc_bkmm->Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root");

  TH1D* h_mc_bkmm_base = new TH1D("h_mc_bkmm_base", "", 50, 5.17, 5.45);
  mc_bkmm->Draw("bkmm_jpsimc_mass>>h_mc_bkmm_base", base_cut + "bkmm_gen_pdgId!=0", "goff");
  h_mc_bkmm_base->Write();

  TH1D* h_mc_bkmm_kaon5 = new TH1D("h_mc_bkmm_kaon5", "", 50, 5.17, 5.45);
  mc_bkmm->Draw("bkmm_jpsimc_mass>>h_mc_bkmm_kaon5", base_cut + "bkmm_gen_pdgId!=0&&bkmm_kaon_pt>5", "goff");
  h_mc_bkmm_kaon5->Write();

  for (auto var: variables){
    std::string hist_name = "h_mc_bkmm_" + std::get<0>(var);
    std::string draw_command = std::get<1>(var) + ":bkmm_jpsimc_mass>>h_mc_bkmm_" + std::get<0>(var);
    TH2D* h_mc_bkmm = new TH2D(hist_name.c_str(), "", 
			       50, 5.17, 5.45, 
			       std::get<2>(var), std::get<3>(var), std::get<4>(var));
    mc_bkmm->Draw(draw_command.c_str(), base_cut + "bkmm_gen_pdgId!=0", "goff");
    h_mc_bkmm->Write();
  }

  TH2D* h_mc_bkmm_mva_vs_npv = new TH2D("h_mc_bkmm_mva_vs_npv", "", 50, 0, 100, 50, 0.0, 1.0);
  mc_bkmm->Draw("bkmm_bmm_mva:PV_npvs>>h_mc_bkmm_mva_vs_npv", base_cut + "bkmm_gen_pdgId!=0", "goff");
  h_mc_bkmm_mva_vs_npv->Write();

  TH2D* h_mc_bkmm_bdt_vs_npv = new TH2D("h_mc_bkmm_bdt_vs_npv", "", 50, 0, 100, 50, -1.0, 1.0);
  mc_bkmm->Draw("bkmm_bmm_bdt:PV_npvs>>h_mc_bkmm_bdt_vs_npv", base_cut + "bkmm_gen_pdgId!=0", "goff");
  h_mc_bkmm_bdt_vs_npv->Write();

  TChain* data = new TChain("Events");
  data->Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/*.root");
  // data->Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/Charmonium+Run2018D-PromptReco-v2+MINIAOD/142c84fb8e8d03eea34e0cf3277a0656.root");

  TH1D* h_data_base = new TH1D("h_data_base", "", 50, 5.17, 5.45);
  data->Draw("bkmm_jpsimc_mass>>h_data_base", base_cut, "goff");
  h_data_base->Write();

  TH1D* h_data_kaon5 = new TH1D("h_data_kaon5", "", 50, 5.17, 5.45);
  data->Draw("bkmm_jpsimc_mass>>h_data_kaon5", base_cut + "bkmm_kaon_pt>5", "goff");
  h_data_kaon5->Write();

  for (auto var: variables){
    std::string hist_name = "h_data_" + std::get<0>(var);
    std::string draw_command = std::get<1>(var) + ":bkmm_jpsimc_mass>>h_data_" + std::get<0>(var);
    TH2D* h_data = new TH2D(hist_name.c_str(), "", 
			    50, 5.17, 5.45, 
			    std::get<2>(var), std::get<3>(var), std::get<4>(var));
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

Result fitHistogram(TH1* h_ref, TH1* h_test){
  if (h_ref->Integral() < 1 or h_test->Integral() < 1) return Result();

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
  RooRealVar a1("a1", "a1", 0.0, -0.2, 0.2) ;
  // RooRealVar a2("a2", "a2", 0.0, -1., 1.) ;
  RooChebychev bkg("bkg", "Background", mass, RooArgSet(a0, a1));
  RooRealVar Nsig("Nsig", "Nsig", h_test->Integral(), 0, h_test->Integral());
  RooRealVar Nbkg("Nbkg", "Nbkg", 0, 0, h_test->Integral());
  RooAddPdf model("model", "", RooArgList(sig,bkg), RooArgList(Nsig,Nbkg));

  // test dataset
  RooDataHist test_data("test_data", "", mass, h_test);

  // freeze signal shape to get background shape right 
  sigma.setConstant(true);
  bias.setConstant(true);
  model.fitTo(test_data);

  sigma.setConstant(false);
  bias.setConstant(false);
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

void process_variables(TFile* f){
  for (auto var: variables){
    unsigned int n = std::get<2>(var);
    std::vector<Double_t> mc_eff;
    std::vector<Double_t> data_eff;
    std::vector<Double_t> data_over_mc_eff;
    double data_total(-1);
    double mc_total(-1);
    // std::vector<Double_t> mc_eff_err;
    // std::vector<Double_t> data_over_mc_eff_err;

    for (unsigned int i=0; i <= n; ++i){
      printf("processing %s %u\n",  std::get<0>(var).c_str(), i);
      auto proj_mc = ((TH2*)f->Get(("h_mc_bkmm_" + std::get<0>(var)).c_str()))->ProjectionX("proj_mc", i, n+1);
      auto proj_data = ((TH2*)f->Get(("h_data_" + std::get<0>(var)).c_str()))->ProjectionX("proj_data", i, n+1);
      // if (i==0){
      // 	proj_mc->Draw();
      // 	print_canvas(std::get<0>(var) + "_h_mc_proj" + std::to_string(i), output_path, gPad);
      // 	proj_data->Draw();
      // 	print_canvas(std::get<0>(var) + "_h_data_proj" + std::to_string(i), output_path, gPad);
      // }
      auto result = fitHistogram(proj_mc, proj_data);
      print_canvas(std::get<0>(var) + "_proj" + std::to_string(i), output_path, gPad);
      if (i==0){
	data_total = result.sig;
	mc_total = proj_mc->Integral();
      }
      assert(data_total>0);
      data_eff.push_back(result.sig/data_total);
      assert(mc_total>0);
      mc_eff.push_back(proj_mc->Integral()/mc_total);
      data_over_mc_eff.push_back(mc_eff.back()>0?data_eff.back()/mc_eff.back():0);
    }
    auto gr = new TGraph(n, &mc_eff[0], &data_over_mc_eff[0]);
    gr->SetMinimum(0);
    gr->SetMaximum(1.5);
    gr->GetXaxis()->SetLimits(0.,1.);
    gr->GetXaxis()->SetTitle("MC efficiency");
    gr->GetYaxis()->SetTitle("Data/MC efficiency ratio");
    gr->SetTitle(std::get<0>(var).c_str());
    gr->Draw("AP*");
    print_canvas(std::get<0>(var) + "_eff", output_path, gPad);
  }
  
}

void btojpsik_mva_study(){
  add_variables();

  const char* file_name = "btojpsik_mva_study.root";
  // Check if we already have histograms prepared
  if (gSystem->AccessPathName(file_name)){
    printf("Recreating histograms\n");
    TFile* f = TFile::Open(file_name,"RECREATE");
    add_data(f);    
    f->Close();
  }
  TFile* f = TFile::Open(file_name);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  auto result_all = fitHistogram((TH1*)f->Get("h_mc_bkmm_base"), 
				 (TH1*)f->Get("h_data_base"));

  print_canvas("base", output_path, c1);
  
  process_variables(f);

  return;

  unsigned int n = 50;

  // MVA
  std::vector<std::tuple<Result, unsigned int>> results_mva;

  for (unsigned int i=1; i <= n; ++i){
    auto proj_mc = ((TH2*)f->Get("h_mc_bkmm_mva"))->ProjectionX("proj_mc", i);
    auto proj_data = ((TH2*)f->Get("h_data_mva"))->ProjectionX("proj_data", i);
    auto result = fitHistogram(proj_mc, proj_data);
    results_mva.push_back(std::make_tuple(result, proj_mc->Integral()));
    print_canvas("mva_proj"+std::to_string(i), output_path, c1);
  }

  // BDT
  std::vector<std::tuple<Result, unsigned int>> results_bdt;

  for (unsigned int i=1; i <= n; ++i){
    auto proj_mc = ((TH2*)f->Get("h_mc_bkmm_bdt"))->ProjectionX("proj_mc", i);
    auto proj_data = ((TH2*)f->Get("h_data_bdt"))->ProjectionX("proj_data", i);
    auto result = fitHistogram(proj_mc, proj_data);
    results_bdt.push_back(std::make_tuple(result, proj_mc->Integral()));
    print_canvas("bdt_proj"+std::to_string(i), output_path, c1);
  }


  // Pileup treatment
  TH1* h_mc_npv = ((TH2*)f->Get("h_mc_bkmm_mva_vs_npv"))->ProjectionX();
  TH1* h_data_npv = ((TH1*)f->Get("h_data_npv"));
  h_mc_npv->SetLineWidth(2);
  h_mc_npv->SetLineColor(kBlue);
  h_mc_npv->Draw();
  print_canvas("mc_npv", output_path, c1);

  h_data_npv->SetLineWidth(2);
  h_data_npv->SetLineColor(kBlue);
  h_data_npv->Draw();
  print_canvas("data_npv", output_path, c1);
  
  // Results
  printf("Total: %0.1f+/-%0.1f\n", result_all.sig, result_all.sig_err);
  printf("MVA:\n");
  for (unsigned int i=0; i < n; ++i){
    double eff_data = std::get<0>(results_mva[i]).sig/std::get<0>(results_mva[0]).sig;
    double eff_mc = std::get<1>(results_mva[i])/float(std::get<1>(results_mva[0]));
    printf("\t[%u] \tmva_data: %0.1f+/-%0.1f \tmva_mc:%u \teff_data: %0.1f%% \teff_mc: %0.1f%% \teff_mc/eff_data: %0.2f\n", 
	   i+1,
	   std::get<0>(results_mva[i]).sig, 
	   std::get<0>(results_mva[i]).sig_err, 
	   std::get<1>(results_mva[i]),
	   eff_data*100.,
	   eff_mc*100.,
	   eff_data>0?eff_mc/eff_data:0
	   ); 
  }
  printf("BDT:\n");
  for (unsigned int i=0; i < n; ++i){
    double eff_data = std::get<0>(results_bdt[i]).sig/std::get<0>(results_bdt[0]).sig;
    double eff_mc = std::get<1>(results_bdt[i])/float(std::get<1>(results_bdt[0]));
    printf("\t[%u] \tbdt_data: %0.1f+/-%0.1f \tbdt_mc:%u \teff_data: %0.1f%% \teff_mc: %0.1f%% \teff_mc/eff_data: %0.2f\n", 
	   i+1,
	   std::get<0>(results_bdt[i]).sig, 
	   std::get<0>(results_bdt[i]).sig_err, 
	   std::get<1>(results_bdt[i]),
	   eff_data*100.,
	   eff_mc*100.,
	   eff_data>0?eff_mc/eff_data:0
	   ); 
  }
}
