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

using namespace RooFit;
using namespace std;

string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/sensitivity/";

const bool do_bsmm_significance = false;
const bool recompute_scale_factors = false;
const bool remake_input_workspaces = false;
const bool update_scale_factors = false;
const bool update_cross_sections = false;
const bool produce_2D_conditional_projections = false;

const double mm_mass_min = 4.9;
const double mm_mass_max = 5.7;
const double mm_mass_err_min = 0.005;
const double mm_mass_err_max = 0.100;

bool silent_roofit = true;
bool store_projections = false;
// bool use_mc_truth_matching = true;

struct Sample{
  string name, files, selection;
  float cross_section, scale_factor;
  bool truth_match, blind;
};

const float luminosity = 140e3; // [1/pb]

string storage_path = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/508/";
    
// Don't use symbols in the name
vector<Sample> samples{
  { "bsmm", 
      // storage_path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/35C867FA-837D-1F4C-82B7-694DBC862D16.root",
      storage_path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/*.root",
      "", 2.98E-02, 1.0, true, false
      },
  { "bmm", 
      // storage_path + "BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/4BC0052B-D9BF-6E46-9E0B-C56BC7BA107A.root",
      storage_path + "BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root",
      "", 3.26E-03, 1.0, true, false
      },
  { "bkpi", 
      // storage_path + "BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/0E1049D2-8F42-E811-9DC2-7845C4FC37A9.root",
      storage_path + "BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1+MINIAODSIM/*.root",
      "", 5.76E+02, 2*2e-3*2e-3, // muon fakes
      true, false
      },
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
  // canvas->Print((s + ".root").c_str());
}

// void process_sample(TChain* chain, string name, string selection=""){
const RooWorkspace* process_sample(Sample& sample){
  string file_name = "sensitivity_study-" + sample.name + ".root";
  
  // Check if we already have histograms prepared
  
  if (not remake_input_workspaces and not gSystem->AccessPathName(file_name.c_str())){
    TFile *f = TFile::Open(file_name.c_str()) ;
    if (not f) 
      throw std::runtime_error( "File is not found" );
    const RooWorkspace* ws = dynamic_cast<const RooWorkspace*>(f->Get("workspace"));
    if (not ws) 
      throw std::runtime_error( "Workspace is not found" );
    return ws;
      
  }
  
  printf("Recreating histograms for %s\n", sample.name.c_str());
    
  RooRealVar mass("mass", "mass", mm_mass_min, mm_mass_max);
  RooRealVar mass_err("mass_err", "mass_err", mm_mass_err_min, mm_mass_err_max);
  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooDataSet data("data", "data", RooArgSet(mass, mass_err));
  printf("RooDataSet storage type: %u\n", RooAbsData::getDefaultStorageType());
  TChain* chain = new TChain("Events");
  chain->Add(sample.files.c_str());
  Long64_t n = chain->GetEntries();
  printf("Total number of events: %lld\n", n);

  // Setup variables and branches to read data efficiently
  // mm
  UInt_t     nmm;
  TBranch* b_nmm;
  chain->SetBranchAddress("nmm", &nmm, &b_nmm);
  Int_t      mm_gen_pdgId[200];
  TBranch* b_mm_gen_pdgId;
  chain->SetBranchAddress("mm_gen_pdgId", mm_gen_pdgId, &b_mm_gen_pdgId);
  Int_t      mm_mu1_index[200];
  TBranch* b_mm_mu1_index;
  chain->SetBranchAddress("mm_mu1_index", mm_mu1_index, &b_mm_mu1_index);
  Int_t      mm_mu2_index[200];
  TBranch* b_mm_mu2_index;
  chain->SetBranchAddress("mm_mu2_index", mm_mu2_index, &b_mm_mu2_index);
  Float_t    mm_mu1_pt[200];
  TBranch* b_mm_mu1_pt;
  chain->SetBranchAddress("mm_mu1_pt", mm_mu1_pt, &b_mm_mu1_pt);
  Float_t    mm_mu2_pt[200];
  TBranch* b_mm_mu2_pt;
  chain->SetBranchAddress("mm_mu2_pt", mm_mu2_pt, &b_mm_mu2_pt);
  Float_t    mm_mu1_eta[200];
  TBranch* b_mm_mu1_eta;
  chain->SetBranchAddress("mm_mu1_eta", mm_mu1_eta, &b_mm_mu1_eta);
  Float_t    mm_mu2_eta[200];
  TBranch* b_mm_mu2_eta;
  chain->SetBranchAddress("mm_mu2_eta", mm_mu2_eta, &b_mm_mu2_eta);
  Float_t    mm_kin_sl3d[200];
  TBranch* b_mm_kin_sl3d;
  chain->SetBranchAddress("mm_kin_sl3d", mm_kin_sl3d, &b_mm_kin_sl3d);
  Float_t    mm_kin_mass[200];
  TBranch* b_mm_kin_mass;
  chain->SetBranchAddress("mm_kin_mass", mm_kin_mass, &b_mm_kin_mass);
  Float_t    mm_kin_massErr[200];
  TBranch* b_mm_kin_massErr;
  chain->SetBranchAddress("mm_kin_massErr", mm_kin_massErr, &b_mm_kin_massErr);
  Float_t    mm_kin_vtx_chi2dof[200];
  TBranch* b_mm_kin_vtx_chi2dof;
  chain->SetBranchAddress("mm_kin_vtx_chi2dof", mm_kin_vtx_chi2dof, &b_mm_kin_vtx_chi2dof);
  // Muon
  Float_t    muon_softMva[200];
  TBranch* b_muon_softMva;
  chain->SetBranchAddress("Muon_softMva", muon_softMva, &b_muon_softMva);

  Long64_t n_good_cands(0);
  Long64_t n_good_events(0);
    
  for (unsigned int i=0; i < n; ++i){
    // Load a proper try and get relative even index
    Long64_t localEntry = chain->LoadTree(i);
    b_nmm->GetEntry(localEntry);
    if (nmm == 0) continue;
    // Get relevant branches
    b_mm_gen_pdgId->GetEntry(localEntry);
    b_mm_mu1_index->GetEntry(localEntry);
    b_mm_mu2_index->GetEntry(localEntry);
    b_mm_mu1_eta->GetEntry(localEntry);
    b_mm_mu2_eta->GetEntry(localEntry);
    b_mm_mu1_pt->GetEntry(localEntry);
    b_mm_mu2_pt->GetEntry(localEntry);
    b_mm_kin_mass->GetEntry(localEntry);
    b_mm_kin_massErr->GetEntry(localEntry);
    b_mm_kin_sl3d->GetEntry(localEntry);
    b_mm_kin_vtx_chi2dof->GetEntry(localEntry);
    b_muon_softMva->GetEntry(localEntry);
    // Loop over candidates
    bool good_event = false;
    for (unsigned int cand=0; cand < nmm; ++cand){
      // std::cout << bkmm_jpsimc_mass[cand] << ", " << bkmm_jpsimc_massErr[cand] << ", " << 
      //   bkmm_mm_index[cand] << ", " << bkmm_jpsimc_vtx_chi2dof[cand] << ", " << bkmm_jpsimc_alpha[cand] << ", " <<
      //   mm_mu1_index[bkmm_mm_index[cand]] << ", " << mm_mu2_index[bkmm_mm_index[cand]] << ", " << 
      //   mm_kin_sl3d[bkmm_mm_index[cand]] << ", " << mm_kin_vtx_chi2dof[bkmm_mm_index[cand]] << ", " << 
      //   muon_eta[mm_mu1_index[bkmm_mm_index[cand]]] << ", " <<  muon_pt[mm_mu1_index[bkmm_mm_index[cand]]] << ", " <<
      //   muon_softMva[mm_mu1_index[bkmm_mm_index[cand]]] << 
      //   std::endl;
      if (sample.truth_match and not mm_gen_pdgId[cand]) continue;
      if (fabs(mm_mu1_eta[cand]) > 1.4) continue;
      if (fabs(mm_mu2_eta[cand]) > 1.4) continue;
      if (mm_mu1_pt[cand] < 4.0) continue;
      if (mm_mu2_pt[cand] < 4.0) continue;
      // if (muon_softMva[mm_mu1_index[bkmm_mm_index[cand]]] < 0.45) continue;
      // if (muon_softMva[mm_mu2_index[bkmm_mm_index[cand]]] < 0.45) continue;
      if (mm_kin_sl3d[cand] < 4.0) continue;
      if (mm_kin_vtx_chi2dof[cand] > 5.0) continue;
      if (mm_kin_mass[cand] < mm_mass_min) continue;
      if (mm_kin_mass[cand] > mm_mass_max) continue;
      if (TMath::IsNaN(mm_kin_massErr[cand])) continue;
      if (mm_kin_massErr[cand] < mm_mass_err_min) continue;
      if (mm_kin_massErr[cand] > mm_mass_err_max) continue;
      good_event = true;
      n_good_cands++;
      mass = mm_kin_mass[cand];
      mass_err = mm_kin_massErr[cand];
      data.add(RooArgSet(mass,mass_err));
    }
    if (good_event) n_good_events++;
    // if (n_good_cands > 10000) break;
  }
  printf("Number of selected events: %lld\n", n_good_events);
  printf("Number of selected candidates: %lld\n", n_good_cands);
  // data.Print("v");

  // Create a new empty workspace
  RooWorkspace *ws = new RooWorkspace("workspace","workspace");
  ws->import(data) ;
  RooRealVar eff("eff", "Event selection efficiency", double(n_good_events)/n);
  ws->import(eff);
  RooRealVar xsec("xsec", "Effective cross-section before event selection, [pb]", sample.cross_section);
  ws->import(xsec);
  RooRealVar scale("scale", "Correction scale factor", sample.scale_factor);
  ws->import(scale);
  ws->writeToFile(file_name.c_str()) ;
  return ws;
  
  

  //   printf("Number of events: %lld\n", chain->GetEntries());
  
  //   TCut mc_match("mm_gen_pdgId!=0");
  //   // mc_match = "abs(1-bkmm_kaon_pt/genbmm_kaon1_pt[0])<0.1";

  //   TCut base_cut;
  //   // TCut base_cut("HLT_DoubleMu4_Jpsi_NoVertexing");
  //   base_cut += "abs(mm_mu1_eta)<1.4 && mm_mu1_pt>4";
  //   base_cut += "abs(mm_mu2_eta)<1.4 && mm_mu2_pt>4";
  //   base_cut += "abs(mm_kin_mass-5.5)<0.5";
  //   base_cut += "!TMath::IsNaN(mm_kin_massErr)";
  //   base_cut += "mm_kin_sl3d>4 && mm_kin_vtx_chi2dof<5";
    
  //   // base_cut += "Muon_softMva[mm_mu1_index[bkmm_mm_index]]>0.45 && Muon_softMva[mm_mu2_index[bkmm_mm_index]]>0.45";

  // string h_mass_name = "h_" + sample.name + "_mass";
  // string h_mass_err_name = "h_" + sample.name + "_mass_err";

  // TH1D* h_mass = new TH1D(h_mass_name.c_str(), "", 100, 5.0, 5.5);
  // TH1D* h_mass_err = new TH1D(h_mass_err_name.c_str(), "", 100, 0.0, 0.1);
  // chain->Draw(("mm_kin_mass>>" + h_mass_name).c_str(), base_cut + mc_match, "goff");
  // chain->Draw(("mm_kin_massErr>>" + h_mass_err_name).c_str(), base_cut + mc_match, "goff");
  // printf("Selected events: %0.0f (efficiency: %0.1f%%)\n", h_mass->Integral(), 
  // 	 100.*h_mass->Integral()/chain->GetEntries());
  // h_mass->SetDirectory(0);
  // h_mass_err->SetDirectory(0);
  // sample.h_mass = h_mass;
  // sample.h_mass_err = h_mass_err;
  // sample.efficiency = h_mass->Integral()/chain->GetEntries();
  // if (recompute_scale_factors){
  //   unsigned int n = chain->GetEntries(base_cut + mc_match);
  //   float mu_eff = n>0? float(chain->GetEntries(base_cut + mc_match + "mm_mu1_index>=0&&mm_mu2_index>=0"))/n : 1;
  //   sample.scale_factor = mu_eff;
  // }

  // // Make a dataset
  // RooRealVar* mass = workspace.var("mm_kin_mass");
  // if (not mass){
  //   workspace.import(RooRealVar("mm_kin_mass", "mass", 5.0, 6.0));
  //   mass = workspace.var("mm_kin_mass");
  //   assert(mass);
  // }
  // RooRealVar* mass_err = workspace.var("mm_kin_massErr");
  // if (not mass_err){
  //   workspace.import(RooRealVar("mm_kin_massErr", "mass uncertainty", 0.0, 0.1));
  //   mass_err = workspace.var("mm_kin_massErr");
  //   assert(mass_err);
  // }

  // TTree* tree = chain->CopyTree(base_cut + mc_match);
  // assert(tree);
  // printf("Number of events selected: %lld\n", tree->GetEntries());
  // RooDataSet data(("data_" + sample.name).c_str(), "data", RooArgSet(*mass, *mass_err), Import(*tree));
  // workspace.import(data);

  // for (auto const& var : variables){
  //   string h_mass_name = "h_mc_bkmm_" + var.name + "_vs_mass";
  //   string h_npv_name  = "h_mc_bkmm_" + var.name + "_vs_npv";
  //   string command_mass = var.branch + ":bkmm_jpsimc_mass>>" + h_mass_name;
  //   string command_npv  = var.branch + ":PV_npvs>>" + h_npv_name;
  //   TH2D* h_mass = new TH2D(h_mass_name.c_str(), "", 50, 5.17, 5.45, var.nbins, var.xmin, var.xmax);
  //   TH2D* h_npv  = new TH2D(h_npv_name.c_str(),  "", 50,    0,  100, var.nbins, var.xmin, var.xmax);
  //   mc_bkmm->Draw(command_mass.c_str(), base_cut + mc_match, "goff");
  //   mc_bkmm->Draw(command_npv.c_str(),  base_cut + mc_match, "goff");
  //   h_mass->Write();
  //   h_npv->Write();
  // }

  // // Data

  // TH1D* h_data_base = new TH1D("h_data_base", "", 50, 5.17, 5.45);
  // data->Draw("bkmm_jpsimc_mass>>h_data_base", base_cut, "goff");
  // h_data_base->Write();

  // for (auto const& var : variables){
  //   string hist_name = "h_data_" + var.name + "_vs_mass";
  //   string draw_command = var.branch + ":bkmm_jpsimc_mass>>" + hist_name;
  //   TH2D* h_data = new TH2D(hist_name.c_str(), "", 
  // 			    50, 5.17, 5.45, 
  // 			    var.nbins, var.xmin, var.xmax);
  //   data->Draw(draw_command.c_str(), base_cut, "goff");
  //   h_data->Write();
  // }

  // TH1D* h_data_npv = new TH1D("h_data_npv", "", 50, 0, 100);
  // data->Draw("PV_npvs>>h_data_npv", base_cut, "goff");
  // h_data_npv->Write();

  // return h_mass;

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

struct ToyStudy{
  TH1* yield_bmm;
  TH1* significance_bmm;
  TH1* yield_bsmm;
  TH1* significance_bsmm;
};

ToyStudy toy_study(const RooWorkspace& ws_ref, string gen_model_name, 
		   string fit_model_name, unsigned int trials){
  if (silent_roofit){
    RooMsgService::instance().saveState();
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  }
  ToyStudy result;
  vector<float> yield_bmm;
  vector<float> yield_bsmm;
  vector<float> significance_bmm;
  vector<float> significance_bsmm;

  TStopwatch stop_watch;
  TStopwatch stop_watch_fit;
  stop_watch_fit.Stop();
  TStopwatch stop_watch_generate;
  stop_watch_generate.Stop();
  TStopwatch stop_watch_init;
  stop_watch_init.Stop();
  for (unsigned int i=0; i<trials; ++i){
    stop_watch_init.Start(false);
    RooWorkspace ws;
    ws.import(*ws_ref.pdf(gen_model_name.c_str()));
    if (gen_model_name != fit_model_name)
      ws.import(*ws_ref.pdf(fit_model_name.c_str()));
    auto mass = ws.var("mass");
    auto n_bmm = ws.var("n_bmm");
    auto n_bsmm = ws.var("n_bsmm");
    auto gen_model = ws.pdf(gen_model_name.c_str());
    auto fit_model = ws.pdf(fit_model_name.c_str());
    auto gen_observables = const_cast<RooWorkspace*>(&ws_ref)->set(gen_model_name.c_str());
    stop_watch_init.Stop();

    stop_watch_generate.Start(false);
    
    RooDataSet *data = gen_model->generate(*gen_observables, Extended(kTRUE));
    stop_watch_generate.Stop();

    stop_watch_fit.Start(false);
    RooFitResult* def_res(nullptr);
    if (silent_roofit){
      def_res = fit_model->fitTo(*data, PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
    } else{
      def_res = fit_model->fitTo(*data, Save(kTRUE));
    }
    stop_watch_fit.Stop();

    yield_bmm.push_back(n_bmm->getVal());
    yield_bsmm.push_back(n_bsmm->getVal());
    n_bmm->setVal(0);
    n_bmm->setConstant(true);

    stop_watch_fit.Start(false);
    auto nobmm_res = fit_model->fitTo(*data, PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
    stop_watch_fit.Stop();
    significance_bmm.push_back(sqrt(max(0.,nobmm_res->minNll() - def_res->minNll())*2.));

    if (do_bsmm_significance){
      n_bmm->setConstant(false);
      n_bsmm->setVal(0);
      n_bsmm->setConstant(true);
      stop_watch_fit.Start(false);
      auto nobsmm_res = fit_model->fitTo(*data, PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
      stop_watch_fit.Stop();
      significance_bsmm.push_back(sqrt(max(0.,nobsmm_res->minNll() - def_res->minNll())*2.));
    }
    delete data;
  }
  stop_watch.Stop();
  stop_watch_fit.Stop();
  stop_watch_generate.Stop();
  printf("Total time per toy study: %0.3f sec\n", stop_watch.RealTime()/trials);
  printf("Generation time per toy study: %0.3f sec\n", stop_watch_generate.RealTime()/trials);
  printf("Fitting time per toy study: %0.3f sec\n", stop_watch_fit.RealTime()/trials);
  printf("Init time per toy study: %0.3f sec\n", stop_watch_init.RealTime()/trials);

  result.yield_bsmm = new TH1F("yield_bsmm", "Toy Study for BsToMuMu Yeild", 100, 0, ws_ref.var("n_bsmm")->getVal()*2);
  result.yield_bsmm->SetDirectory(0);
  for (float n: yield_bsmm){
    result.yield_bsmm->Fill(n);
  }
  result.yield_bmm = new TH1F("yield_bmm", "Toy Study for BdToMuMu Yeild", 100, 0, ws_ref.var("n_bmm")->getVal()*2);
  result.yield_bmm->SetDirectory(0);
  for (float n: yield_bmm){
    result.yield_bmm->Fill(n);
  }
  result.significance_bmm = new TH1F("significance_bmm", "Toy Study for BdToMuMu Significance", 100, 0, 15);
  result.significance_bmm->SetDirectory(0);
  for (float s: significance_bmm){
    result.significance_bmm->Fill(s);
  }
  result.significance_bsmm = new TH1F("significance_bsmm", "Toy Study for BsToMuMu Significance", 100, 0, 100);
  result.significance_bsmm->SetDirectory(0);
  for (float s: significance_bsmm){
    result.significance_bsmm->Fill(s);
  }

  if (silent_roofit)
    RooMsgService::instance().restoreState();
  return result;
}

float get_expected_event_yield(const RooWorkspace& workspace, string name){
  const Sample* sample(nullptr);
  for (auto s: samples)
    if (s.name == name)
      sample = &s;
  if (not sample)
    throw std::runtime_error("Cannot find sample with name: " + name);
    
  float cross_section = update_cross_sections ? sample->cross_section : workspace.var(("xsec_"  + name).c_str())->getVal();
  float scale_factor  = update_scale_factors ? sample->scale_factor : workspace.var(("scale_" + name).c_str())->getVal();
  float efficiency = workspace.var(("eff_"  + name).c_str())->getVal();
  return luminosity * cross_section * scale_factor * efficiency;
}

void build_model_1D(RooWorkspace& workspace){
  const char* model_name = "model_1D";
  RooArgList pdfs;
  RooArgList yields;
  auto mass = workspace.var("mass");
  if (not mass)
    throw std::runtime_error("Cannot get mass");
  
  for (auto& sample: samples){
    // auto ds = workspace.data(("data_" + sample.name).c_str());
    // if (not ds)
    //   throw std::runtime_error("Cannot get data_"  + sample.name);
    auto data = dynamic_cast<const RooDataHist*>(workspace.data(("data_hist_" + sample.name).c_str()));
    if (not data)
      throw std::runtime_error("Cannot get data_hist_" + sample.name);

    // Use 2nd order interpolation to smooth the shape
    RooHistPdf* pdf = new RooHistPdf(("pdf_" + sample.name + "_1D").c_str(), "", *mass, *data, 2);
    pdf->Print();
    pdfs.add(*pdf);
    printf("pdfs.getSize: %u\n", pdfs.getSize());
    float n_expected = get_expected_event_yield(workspace, sample.name);
    RooRealVar* n = new RooRealVar(("n_" + sample.name).c_str(), "blah", n_expected, 0., 1e9);
    n->Print();
    yields.add(*n);
    printf("yields.getSize: %u\n", yields.getSize());
  }
  
  RooAddPdf model(model_name, model_name, pdfs, yields);
  
  workspace.import(model);
  // Specify observables
  workspace.defineSet(model_name, "mass");
  workspace.Print();

  auto frame = mass->frame() ;

  model.plotOn(frame);
  // if (workspace.pdf("pdf_bkpi")){
  //   model.plotOn(frame, Components(*(workspace.pdf("pdf_bkpi"))), FillColor(kMagenta), FillStyle(1001), DrawOption("F"));
  //   model.plotOn(frame, Components(RooArgSet(*(workspace.pdf("pdf_bmm")), *(workspace.pdf("pdf_bkpi")))), LineStyle(kDashed));
  // } else {
  //   model.plotOn(frame, Components(*(workspace.pdf("pdf_bmm"))), LineStyle(kDashed));
  // }
  frame->Draw();

  print_canvas(model_name, output_path, gPad);
}

void build_model_2D(RooWorkspace& workspace){
  auto mass = workspace.var("mass");
  if (not mass)
    throw std::runtime_error("Cannot get mass var");
  auto mass_err = workspace.var("mass_err");
  if (not mass_err)
    throw std::runtime_error("Cannot get mass_err var");

  RooRealVar mean_bsmm("mean_bsmm", "mean_bsmm", 5.35, 5.32, 5.40);
  RooRealVar mean_bmm( "mean_bmm",  "mean_bmm",  5.29, 5.25, 5.31);
  RooRealVar mean_bkpi("mean_bkpi", "mean_bkpi", 5.23, 5.20, 5.25);

  RooRealVar s_mm("s_mm","resolution", 1, 0, 20);
  RooFormulaVar sigma_mm("sigma_mm", "", "@0*@1", RooArgList(s_mm, *mass_err));

  RooRealVar s_kpi("s_kpi","resolution", 1, 0, 20);
  RooFormulaVar sigma_kpi("sigma_kpi", "", "@0*@1", RooArgList(s_kpi, *mass_err));

  // RooGaussian pdf_bsmm("pdf_bsmm", "pdf_bsmm", *mass, mean_bsmm, sigma_mm); 
  // RooGaussian pdf_bsmm("pdf_bsmm", "pdf_bsmm", *mass, mean_bsmm, *mass_err); 
  RooRealVar cb_alpha_mm("cb_alpha_mm", "Crystal-Ball gaussian cutoff for mm signal", 1, 0.1, 10);
  RooRealVar cb_n_mm("cb_n_mm", "Crystal-Ball exponent of the power-low tail for mm signal", 5, 0, 10000);
  RooCBShape pdf_bsmm("pdf_bsmm", "pdf_bsmm", *mass, mean_bsmm, sigma_mm, cb_alpha_mm, cb_n_mm);
  RooCBShape pdf_bmm( "pdf_bmm",  "pdf_bmm",  *mass, mean_bmm,  sigma_mm, cb_alpha_mm, cb_n_mm);

  RooRealVar cb_alpha_kpi("cb_alpha_kpi", "Crystal-Ball gaussian cutoff for kpi events", 1, 0.1, 10);
  RooRealVar cb_n_kpi("cb_n_kpi", "Crystal-Ball exponent of the power-low tail for kpi events", 5, 0, 10000);
  RooCBShape pdf_bkpi("pdf_bkpi", "pdf_bkpi", *mass, mean_bkpi, sigma_kpi, cb_alpha_kpi, cb_n_kpi);

  // RooGaussian pdf_bmm( "pdf_bmm",  "pdf_bmm",  *mass, mean_bmm,  sigma_mm); 
  // RooGaussian pdf_bkpi("pdf_bkpi", "pdf_bkpi", *mass, mean_bkpi, sigma_kpi);
  

  // Process BsToMuMu

  auto data_bsmm = workspace.data("data_bsmm");
  if (not data_bsmm)
    throw std::runtime_error("Cannot get data_bsmm");

  // Get Mass Uncertainty PDF
  // RooDataHist* expHistDterr = expDataDterr->binnedClone() ;
  RooDataHist data_hist_mass_err_bsmm("data_hist_mass_err_bsmm", "", *mass_err, *data_bsmm);
  RooHistPdf pdf_mass_err_bsmm("pdf_mass_err_bsmm", "", *mass_err, data_hist_mass_err_bsmm, 2);
  RooPlot* frame_mass_err_bsmm = mass_err->frame();
  pdf_mass_err_bsmm.plotOn(frame_mass_err_bsmm);
  frame_mass_err_bsmm->Draw();
  print_canvas("mass_err_bsmm", output_path, gPad);
  
  // Fit signal model to extract parameters
  pdf_bsmm.fitTo(*data_bsmm, ConditionalObservables(*mass_err), NumCPU(16), Timer(true));

  // Make a full PDF
  RooProdPdf pdf2_bsmm("pdf2_bsmm","pdf2_bsmm", pdf_mass_err_bsmm, Conditional(pdf_bsmm, *mass));

  // Plot results of the fit
  RooPlot* frame_bsmm = mass->frame();
  data_bsmm->plotOn(frame_bsmm);
  pdf2_bsmm.plotOn(frame_bsmm);
  // pdf_bsmm.plotOn(frame_bsmm, ProjWData(*mass_err, *data_bsmm));
  // pdf_bsmm.plotOn(frame_bsmm, ProjWData(*mass_err, *data_bsmm), NumCPU(16), Normalization(data_bsmm->sumEntries(), RooAbsReal::NumEvent));
  frame_bsmm->Draw();
  print_canvas("gaus2D_bsmm", output_path, gPad);

  cb_alpha_mm.setConstant(true);
  cb_n_mm.setConstant(true);
  s_mm.setConstant(true);

  // Process BdToMuMu

  auto data_bmm = workspace.data("data_bmm");
  if (not data_bmm)
    throw std::runtime_error("Cannot get data_bmm");

  // Fit Bmm signal model to extract parameters
  pdf_bmm.fitTo(*data_bmm, ConditionalObservables(*mass_err), NumCPU(16), Timer(true));

  // Get Mass Uncertainty PDF
  RooDataHist data_hist_mass_err_bmm("data_hist_mass_err_bmm", "", *mass_err, *data_bmm);
  RooHistPdf pdf_mass_err_bmm("pdf_mass_err_bmm", "", *mass_err, data_hist_mass_err_bmm, 2);
  RooPlot* frame_mass_err_bmm = mass_err->frame();
  pdf_mass_err_bmm.plotOn(frame_mass_err_bmm);
  frame_mass_err_bmm->Draw();
  print_canvas("mass_err_bmm", output_path, gPad);

  // Make a full PDF
  RooProdPdf pdf2_bmm("pdf2_bmm","pdf2_bmm", pdf_mass_err_bmm, Conditional(pdf_bmm, *mass));

  RooPlot* frame_bmm = mass->frame();
  data_bmm->plotOn(frame_bmm);
  pdf2_bmm.plotOn(frame_bmm);
  // pdf_bmm.plotOn(frame_bmm, ProjWData(*mass_err, *data_bmm));
  frame_bmm->Draw();
  print_canvas("gaus2D_bmm", output_path, gPad);

  // Process BToKPi

  auto data_bkpi = workspace.data("data_bkpi");
  if (not data_bkpi)
    throw std::runtime_error("Cannot get data_bkpi");
  
  // Fit KPi model to extract parameters
  pdf_bkpi.fitTo(*data_bkpi, ConditionalObservables(*mass_err), NumCPU(16), Timer(true));

  // Get Mass Uncertainty PDF
  RooDataHist data_hist_mass_err_bkpi("data_hist_mass_err_bkpi", "", *mass_err, *data_bkpi);
  RooHistPdf pdf_mass_err_bkpi("pdf_mass_err_bkpi", "", *mass_err, data_hist_mass_err_bkpi, 2);
  RooPlot* frame_mass_err_bkpi = mass_err->frame();
  pdf_mass_err_bkpi.plotOn(frame_mass_err_bkpi);
  frame_mass_err_bkpi->Draw();
  print_canvas("mass_err_bkpi", output_path, gPad);

  // Make a full PDF
  RooProdPdf pdf2_bkpi("pdf2_bkpi","pdf2_bkpi", pdf_mass_err_bkpi, Conditional(pdf_bkpi, *mass));

  RooPlot* frame_bkpi = mass->frame();
  data_bkpi->plotOn(frame_bkpi);
  pdf2_bkpi.plotOn(frame_bkpi);
  // pdf_bkpi.plotOn(frame_bkpi, ProjWData(*mass_err, *data_bkpi));
  frame_bkpi->Draw();
  print_canvas("gaus2D_bkpi", output_path, gPad);


  RooRealVar n_bsmm("n_bsmm", "n_bsmm", get_expected_event_yield(workspace, "bsmm"), 0., 1e9); 
  RooRealVar n_bmm( "n_bmm",  "n_bmm",  get_expected_event_yield(workspace, "bmm"),  0., 1e9); 
  RooRealVar n_bkpi("n_bkpi", "n_bsmm", get_expected_event_yield(workspace, "bkpi"), 0., 1e9); 

  RooAddPdf model("model","model", RooArgList(pdf2_bsmm, pdf2_bmm, pdf2_bkpi), RooArgList(n_bsmm, n_bmm, n_bkpi));
  
  workspace.import(model);
  workspace.Print();
}

// Result fitHistogram(TH1* h_ref, TH1* h_test){
//   if (h_ref->Integral() < 1 or h_test->Integral() < 1) return Result();

//   if (silent_roofit)
//     RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

//   // signal pdf
//   RooRealVar bias("bias", "bias", 0, -0.1, 0.1) ;
//   RooRealVar sigma("sigma", "sigma", 0.0001, 0., 0.01);
//   RooGaussModel gaussM("gaussM", "signal pdf", mass, bias, sigma) ;
//   mass.setBins(10000, "fft");
//   RooFFTConvPdf sig("sig", "smeared distribution", mass, ref_pdf, gaussM);

//   // background pdf
//   RooRealVar a0("a0", "a0", 0.0, -1., 1.) ;
//   // RooRealVar a1("a1", "a1", 0.0, -0.2, 0.2) ;
//   // RooRealVar a2("a2", "a2", 0.0, -1., 1.) ;
//   // RooChebychev bkg("bkg", "Background", mass, RooArgSet(a0, a1));
//   RooChebychev bkg("bkg", "Background", mass, RooArgSet(a0));
//   RooRealVar Nsig("Nsig", "Nsig", h_test->Integral(), 0, h_test->Integral());
//   RooRealVar Nbkg("Nbkg", "Nbkg", 0, 0, h_test->Integral());
//   RooAddPdf model("model", "", RooArgList(sig,bkg), RooArgList(Nsig,Nbkg));

//   // test dataset
//   RooDataHist test_data("test_data", "", mass, h_test);

//   // freeze signal shape to get background shape right 
//   sigma.setConstant(true);
//   bias.setConstant(true);
//   if (silent_roofit)
//     model.fitTo(test_data, PrintEvalErrors(-1), PrintLevel(-1));
//   else
//     model.fitTo(test_data);
    
//   sigma.setConstant(false);
//   bias.setConstant(false);
//   if (silent_roofit)
//     model.fitTo(test_data, PrintEvalErrors(-1), PrintLevel(-1));
//   else
//     model.fitTo(test_data);

//   RooPlot* frame = mass.frame();
//   test_data.plotOn(frame);
//   model.plotOn(frame);
//   model.plotOn(frame, RooFit::Components(bkg), RooFit::LineStyle(kDashed));
//   // data.statOn(frame, Layout(0.55, 0.99, 0.8));
//   // model->paramOn(frame, Parameters(params), Layout(0.6,0.9,0.9) ) ;

//   model.paramOn(frame, Layout(0.6, 0.85, 0.85));
//   // model->paramOn(frame, Layout(0.55));
//   frame->getAttText()->SetTextSize(0.02);
//   frame->Draw();

//   return Result(Nsig,Nbkg);
// }

// void process_variables(TFile* f, string prefix){
//   TH2* h_mc_bkmm_mva_vs_npv = (TH2*)f->Get("h_mc_bkmm_mva_vs_npv");
//   TH1* h_mc_npv = h_mc_bkmm_mva_vs_npv->ProjectionX();
//   TH1* h_data_npv = ((TH1*)f->Get("h_data_npv"));
//   h_mc_npv->SetLineWidth(2);
//   h_mc_npv->SetLineColor(kBlue);
//   h_mc_npv->Draw();
//   print_canvas(prefix + "-" + "mc_npv", output_path, gPad);

//   h_data_npv->SetLineWidth(2);
//   h_data_npv->SetLineColor(kBlue);
//   h_data_npv->Draw();
//   print_canvas(prefix + "-" + "data_npv", output_path, gPad);

//   for (auto const& var: variables){
//     unsigned int n = var.nbins;
//     vector<Double_t> mc_eff;
//     vector<Double_t> data_eff;
//     vector<Double_t> data_over_mc_eff;
//     vector<Double_t> mc_eff_reweighted;
//     vector<Double_t> data_over_mc_eff_reweighted;
//     double total_data(-1);
//     double total_mc(-1);
//     // vector<Double_t> mc_eff_err;
//     // vector<Double_t> data_over_mc_eff_err;

//     printf("processing %s\n", var.name.c_str());
//     // Scan efficiency in MC and data as a function of the cut on the variable

//     auto h_mc_x_vs_npv    = (TH2*)f->Get(("h_mc_bkmm_"   + var.name + "_vs_npv" ).c_str()); 
//     auto h_mc_x_vs_mass   = (TH2*)f->Get(("h_mc_bkmm_"   + var.name + "_vs_mass").c_str()); 
//     auto h_data_x_vs_mass = (TH2*)f->Get(("h_data_" + var.name + "_vs_mass").c_str());
//     auto h_mc_x_reweighted = reweight_histogram(h_mc_x_vs_npv, h_mc_npv, h_data_npv);
    
//     for (unsigned int i=0; i <= n; ++i){
//       auto mass_mc   = h_mc_x_vs_mass->ProjectionX(  "mass_mc",   i, n+1);
//       auto mass_data = h_data_x_vs_mass->ProjectionX("mass_data", i, n+1);
//       auto result = fitHistogram(mass_mc, mass_data);
      
//       if (store_projections)
// 	print_canvas(prefix + "-" + var.name + "_proj" + to_string(i), output_path, gPad);

//       if (i==0){
// 	total_data = result.sig;
// 	total_mc = mass_mc->Integral();
//       }

//       assert(total_data > 0);
//       data_eff.push_back(result.sig / total_data);
//       assert(total_mc>0);
//       mc_eff.push_back(mass_mc->Integral() / total_mc);
//       data_over_mc_eff.push_back(mc_eff.back()>0?data_eff.back()/mc_eff.back():0);

//       mc_eff_reweighted.push_back(h_mc_x_reweighted->Integral(i,-1) / h_mc_x_reweighted->Integral(0,-1));
//       data_over_mc_eff_reweighted.push_back(mc_eff_reweighted.back()>0?data_eff.back()/mc_eff_reweighted.back():0);
//     }
//     gPad->SetGridx();
//     gPad->SetGridy();
//     auto gr = new TGraph(n, &mc_eff[0], &data_over_mc_eff[0]);
//     gr->SetMinimum(0);
//     gr->SetMaximum(1.5);
//     gr->GetXaxis()->SetLimits(0.,1.);
//     gr->GetXaxis()->SetTitle("MC efficiency");
//     gr->GetYaxis()->SetTitle("Data/MC efficiency ratio");
//     gr->SetTitle(var.name.c_str());
//     gr->Draw("AP*");

//     print_canvas(prefix + "-" + var.name + "_eff", output_path, gPad);

//     auto gr_reweighted = new TGraph(n, &mc_eff_reweighted[0], &data_over_mc_eff_reweighted[0]);
//     gr_reweighted->SetMinimum(0);
//     gr_reweighted->SetMaximum(1.5);
//     gr_reweighted->GetXaxis()->SetLimits(0.,1.);
//     gr_reweighted->GetXaxis()->SetTitle("MC efficiency");
//     gr_reweighted->GetYaxis()->SetTitle("Data/MC efficiency ratio");
//     gr_reweighted->SetTitle(var.name.c_str());
//     gr_reweighted->Draw("AP*");
//     print_canvas(prefix + "-" + var.name + "_eff_reweighted", output_path, gPad);
//     gPad->SetGridx(0);
//     gPad->SetGridy(0);
//   }
// }

void import_var(RooWorkspace& target_ws, const RooWorkspace& source_ws, const char* var_name, const char* new_name){
  auto var = source_ws.var(var_name);
    if (not var) 
      throw std::runtime_error(string(var_name) +" is missing");
    RooRealVar new_var(*var);
    new_var.SetName(new_name);
    target_ws.import(new_var);
}

void sensitivity_study(){
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  RooWorkspace workspace("study_workspace","");
  for (auto& sample: samples){
    // Process sample
    const RooWorkspace* sample_workspace = process_sample(sample);

    // Import sample workspace into the main workspace
    auto sample_data = sample_workspace->data("data");
    // sample_data = sample_data->reduce("mass_err<0.035");
    if (not sample_data) 
      throw std::runtime_error("RooDataSet \"data\" is missing");
    // Tree is easier to use for debugging purposes
    sample_data->convertToTreeStore();
    workspace.import(*sample_data, Rename(("data_" + sample.name).c_str()));

    import_var(workspace, *sample_workspace, "eff",   ("eff_"   + sample.name).c_str());
    import_var(workspace, *sample_workspace, "xsec",  ("xsec_"  + sample.name).c_str());
    import_var(workspace, *sample_workspace, "scale", ("scale_" + sample.name).c_str());

    TTree* tree = const_cast<TTree*>(sample_data->tree());
    if (tree){
      TH1D h_mass(("h_mass_" + sample.name).c_str(), "", 100, mm_mass_min, mm_mass_max);
      tree->Draw(("mass>>h_mass_" + sample.name).c_str());
      RooDataHist rdh(("data_hist_" + sample.name).c_str(), "", *workspace.var("mass"), &h_mass);
      workspace.import(rdh);
      print_canvas("mass_" + sample.name, output_path, c1);

      tree->Draw("mass_err");
      print_canvas("mass_err_" + sample.name, output_path, c1);
    } else {
      if (RooRealVar* mass = sample_workspace->var("mass")){
	RooPlot* frame = mass->frame() ;
	sample_workspace->data("data")->plotOn(frame);
	frame->Draw();
	print_canvas("mass_" + sample.name, output_path, c1);
      }
    }
  }
  workspace.Print("V");

  build_model_1D(workspace);
  // build_model_2D(workspace);

  /*  
  auto mass = workspace.var("mass");
  auto model = workspace.pdf("model");
  
  auto frame = mass->frame() ;

  model->plotOn(frame);
  if (workspace.pdf("pdf_bkpi")){
    model->plotOn(frame, Components(*(workspace.pdf("pdf_bkpi"))), FillColor(kMagenta), FillStyle(1001), DrawOption("F"));
    model->plotOn(frame, Components(RooArgSet(*(workspace.pdf("pdf_bmm")), *(workspace.pdf("pdf_bkpi")))), LineStyle(kDashed));
  } else {
    model->plotOn(frame, Components(*(workspace.pdf("pdf_bmm"))), LineStyle(kDashed));
  }
  frame->Draw();

  print_canvas("model", output_path, c1);
  */

  // Toy study
  auto result = toy_study(workspace, "model_1D", "model_1D", 100);
  
  // result.yield_bsmm->Draw();
  // print_canvas("toy_n_bsmm_2D", output_path, c1);
  // result.yield_bmm->Draw();
  // print_canvas("toy_n_bmm_2D", output_path, c1);
  // result.significance_bmm->Draw();
  // print_canvas("toy_significance_bmm_2D", output_path, c1);
  // result.significance_bsmm->Draw();
  // print_canvas("toy_significance_bsmm_2D", output_path, c1);

  result.yield_bsmm->Draw();
  print_canvas("toy_n_bsmm", output_path, c1);
  result.yield_bmm->Draw();
  print_canvas("toy_n_bmm", output_path, c1);
  result.significance_bmm->Draw();
  print_canvas("toy_significance_bmm", output_path, c1);
  result.significance_bsmm->Draw();
  print_canvas("toy_significance_bsmm", output_path, c1);


  //   // continue;

  //   // TFile* f = TFile::Open(file_name.c_str());
  
  //   // TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  //   // auto result_all = fitHistogram((TH1*)f->Get("h_mc_bkmm_base"), 
  //   // 				   (TH1*)f->Get("h_data_base"));
    


      
    // process_variables(f, ds.name);


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
