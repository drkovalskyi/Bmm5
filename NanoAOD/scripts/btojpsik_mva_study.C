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

void print_canvas(std::string output_name_without_extention, 
		  std::string path, 
		  TCanvas* canvas){
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
  mc_bkmm->Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/0bb34079aa40de8c3fd5615d4816875f.root");
  // mc_bkmm->Add("/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/bkmm/507/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root");

  TH1D* h_mc_bkmm_base = new TH1D("h_mc_bkmm_base", "", 50, 5.17, 5.45);
  mc_bkmm->Draw("bkmm_jpsimc_mass>>h_mc_bkmm_base", base_cut + "bkmm_gen_pdgId!=0", "goff");
  h_mc_bkmm_base->Write();

  TH1D* h_mc_bkmm_kaon5 = new TH1D("h_mc_bkmm_kaon5", "", 50, 5.17, 5.45);
  mc_bkmm->Draw("bkmm_jpsimc_mass>>h_mc_bkmm_kaon5", base_cut + "bkmm_gen_pdgId!=0&&bkmm_kaon_pt>5", "goff");
  h_mc_bkmm_kaon5->Write();

  TH2D* h_mc_bkmm_mva = new TH2D("h_mc_bkmm_mva", "", 50, 5.17, 5.45, 50, 0.0, 1.0);
  mc_bkmm->Draw("bkmm_bmm_mva:bkmm_jpsimc_mass>>h_mc_bkmm_mva", base_cut + "bkmm_gen_pdgId!=0", "goff");
  h_mc_bkmm_mva->Write();

  TH2D* h_mc_bkmm_bdt = new TH2D("h_mc_bkmm_bdt", "", 50, 5.17, 5.45, 50, -1.0, 1.0);
  mc_bkmm->Draw("bkmm_bmm_bdt:bkmm_jpsimc_mass>>h_mc_bkmm_bdt", base_cut + "bkmm_gen_pdgId!=0", "goff");
  h_mc_bkmm_bdt->Write();

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

  TH2D* h_data_mva = new TH2D("h_data_mva", "", 50, 5.17, 5.45, 50, 0.0, 1.0);
  data->Draw("bkmm_bmm_mva:bkmm_jpsimc_mass>>h_data_mva", base_cut, "goff");
  h_data_mva->Write();

  TH2D* h_data_bdt = new TH2D("h_data_bdt", "", 50, 5.17, 5.45, 50, -1.0, 1.0);
  data->Draw("bkmm_bmm_bdt:bkmm_jpsimc_mass>>h_data_bdt", base_cut, "goff");
  h_data_bdt->Write();

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
};

Result fitHistogram(TH1* h_ref, TH1* h_test){
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
  RooRealVar a1("a1", "a1", 0.0, -1., 1.) ;
  // RooRealVar a2("a2", "a2", 0.0, -1., 1.) ;
  RooChebychev bkg("bkg", "Background", mass, RooArgSet(a0, a1));
  RooRealVar Nsig("Nsig", "Nsig", 1000, 0, 1e9);
  RooRealVar Nbkg("Nbkg", "Nbkg", 1, 0, 1e9);
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

void btojpsik_mva_study(){
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
