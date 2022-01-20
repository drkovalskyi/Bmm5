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
#include "RooDataSet.h"
#include <tuple>
#include "TStopwatch.h"
#include "RooProdPdf.h"
#include "RooAcceptReject.h"
#include "TStyle.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"

using namespace RooFit;
using namespace std;

string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv8-517/dimuon-Bs/";

const bool remake_input_workspaces = false;
const bool pre_processing_only = false; // Just produce samples 
const bool fit_samples = false; // no effect on results, just extra info on plots

const double mm_mass_min = 5.2;
const double mm_mass_max = 5.6;
const double mm_mass_err_min = 0.005;
const double mm_mass_err_max = 0.100;

// Running options
const bool silent_roofit = true;

string path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/PostProcessing/FlatNtuples/517/fit/";


vector<string> files =
  {
    path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1+MINIAODSIM/391cf45093fbe4ff6ab26a77246c11e2.root"
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


void dimuon_fit(){
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

  RooRealVar bdt("bdt", "bdt", 0);
  RooRealVar mass("m", "mass", mm_mass_min, mm_mass_max);
  RooRealVar mass_err("me", "mass_err", mm_mass_err_min, mm_mass_err_max);

  TChain* chain = new TChain("bsmmMc");
  for (auto f: files)
    chain->Add(f.c_str());
  Long64_t n = chain->GetEntries();
  printf("Total number of events: %lld\n", n);
  
  // RooDataSet data("data", "data", RooArgSet(mass, mass_err), Import(*chain), Cut("int(me*10000)-125*int(me*80)==0"));
  RooDataSet data("data", "data", RooArgSet(mass, mass_err, bdt), Import(*chain), Cut("bdt>0.995 && int(me*10000)-100*int(me*100)==0"));
  data.Print("V");
  // return;
  
  // Build pdf to be used as a conditinal pdf (mass|mass_err)
  RooRealVar mean("mean", "mean", (mm_mass_max + mm_mass_min)/2, mm_mass_min, mm_mass_max);
  RooRealVar s("s", "resolution", 1, 0, 20);
  RooFormulaVar sigma("sigma", "", "@0*@1", RooArgList(s, mass_err));
  RooRealVar cb_alpha("cb_alpha", "Crystal-Ball gaussian cutoff", 1, 0.1, 10);
  RooRealVar cb_n("cb_n", "Crystal-Ball exponent of the power-low tail", 5, 0, 10000);
  RooCBShape cond_pdf("cond_pdf", "Conditional PDF", mass, mean, sigma, cb_alpha, cb_n);

  // Fit signal model to extract parameters
  cond_pdf.fitTo(data, ConditionalObservables(mass_err), NumCPU(16), Timer(true));
  cb_alpha.setConstant(true);
  cb_n.setConstant(true);
  s.setConstant(true);
  cond_pdf.fitTo(data, ConditionalObservables(mass_err), NumCPU(16), Timer(true));
    
  // Plot results of the fit
  RooPlot* frame = mass.frame();
  data.plotOn(frame);
  // full_pdf.plotOn(frame);
  cond_pdf.plotOn(frame, ProjWData(mass_err, data));
  // pdf_bsmm.plotOn(frame_bsmm, ProjWData(*mass_err, *data_bsmm));
  // pdf_bsmm.plotOn(frame_bsmm, ProjWData(*mass_err, *data_bsmm), NumCPU(16), Normalization(data_bsmm->sumEntries(), RooAbsReal::NumEvent));
  cond_pdf.paramOn(frame, Layout(0.55, 0.85, 0.85));
  frame->getAttText()->SetTextSize(0.03);
  frame->Draw();
  print_canvas("cond_pdf", output_path, gPad);
}
