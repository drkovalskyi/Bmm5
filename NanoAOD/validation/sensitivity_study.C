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

using namespace RooFit;
using namespace std;

string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/sensitivity_mm_only/";

const bool do_bsmm_significance = false;
const bool recompute_scale_factors = false;
const bool remake_input_workspaces = false;
const bool update_scale_factors = false;
const bool update_cross_sections = false;
const bool produce_2D_conditional_projections = false;
const bool no_bhh_in_toys = true;

const double mm_mass_min = 4.9;
const double mm_mass_max = 5.7;
const double mm_mass_err_min = 0.005;
const double mm_mass_err_max = 0.100;

const bool silent_roofit = true;
const bool plot_each_toy = false; // debugging option 

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
  string plot_name = "toy_gen-" + gen_model_name + "_fit-" + fit_model_name;

  bool conditional_fit = fit_model_name == "model_2D_cond";

  TStopwatch stop_watch;
  TStopwatch stop_watch_fit;
  stop_watch_fit.Stop();
  TStopwatch stop_watch_generate;
  stop_watch_generate.Stop();
  TStopwatch stop_watch_init;
  stop_watch_init.Stop();
  printf("Peforming %u toy MC studies\n", trials);
  int i_permille_old = 0;
  for (unsigned int i=0; i < trials; ++i){
    int i_permille = (int)floor(100. * i / trials);
    if (i_permille != i_permille_old) {
      printf("\015\033[32m ---> \033[1m\033[31m%d%%"
             "\033[0m\033[32m <---\033[0m\015", i_permille);
      fflush(stdout);
      i_permille_old = i_permille;
    }
    stop_watch_init.Start(false);
    RooWorkspace ws;
    ws.import(*ws_ref.pdf(gen_model_name.c_str()));
    if (gen_model_name != fit_model_name)
      ws.import(*ws_ref.pdf(fit_model_name.c_str()));
    auto mass = ws.var("mass");
    auto mass_err = ws.var("mass_err");
    auto n_bmm = ws.var("n_bmm");
    auto n_bkpi = ws.var("n_bkpi");
    if (no_bhh_in_toys){
      n_bkpi->setVal(0);
      n_bkpi->setConstant(true);
    }
    auto n_bsmm = ws.var("n_bsmm");
    auto gen_model = ws.pdf(gen_model_name.c_str());
    auto fit_model = ws.pdf(fit_model_name.c_str());
    auto gen_observables = const_cast<RooWorkspace*>(&ws_ref)->set(gen_model_name.c_str());
    stop_watch_init.Stop();

    stop_watch_generate.Start(false);
    
    // RooDataSet *data = gen_model->generate(*gen_observables, Extended(kTRUE));
    // https://sft.its.cern.ch/jira/browse/ROOT-8489
    RooDataSet *data = new RooDataSet("data", "", *gen_observables);
    RooAcceptReject generator(*gen_model, *gen_observables, *gen_model->getGeneratorConfig(), false, 0);
    Double_t resample = 0.;
    int n = n_bsmm->getVal() + n_bmm->getVal() + n_bkpi->getVal();
    for (int i=0; i<n; i++)
      data->add(*generator.generateEvent(n-i,resample));

    stop_watch_generate.Stop();

    stop_watch_fit.Start(false);
    RooFitResult* def_res(nullptr);
    if (silent_roofit){
      if (conditional_fit)
	def_res = fit_model->fitTo(*data, ConditionalObservables(*mass_err), PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
      else
	def_res = fit_model->fitTo(*data, PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
    } else{
      if (conditional_fit)
	def_res = fit_model->fitTo(*data, ConditionalObservables(*mass_err), Save(kTRUE));
      else
	def_res = fit_model->fitTo(*data, Save(kTRUE));
    }
    stop_watch_fit.Stop();

    yield_bmm.push_back(n_bmm->getVal());
    yield_bsmm.push_back(n_bsmm->getVal());
    
    if (plot_each_toy){
      auto frame = mass->frame() ;
      data->plotOn(frame);
      fit_model->plotOn(frame);
      frame->Draw();
      print_canvas(plot_name + "_toy" + to_string(i), output_path, gPad);
      if ( gen_observables->find("mass_err") ){
	auto mass_err = ws.var("mass_err");
	auto frame = mass_err->frame() ;
	data->plotOn(frame);
	// fit_model->plotOn(frame);
	frame->Draw();
	print_canvas(plot_name + "_mass_err_toy" + to_string(i), output_path, gPad);
      }
    }

    n_bmm->setVal(0);
    n_bmm->setConstant(true);

    stop_watch_fit.Start(false);
    RooFitResult* nobmm_res(nullptr);
    if (conditional_fit)
      nobmm_res = fit_model->fitTo(*data, ConditionalObservables(*mass_err), PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
    else
      nobmm_res = fit_model->fitTo(*data, PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
    stop_watch_fit.Stop();
    significance_bmm.push_back(sqrt(max(0.,nobmm_res->minNll() - def_res->minNll())*2.));

    if (do_bsmm_significance){
      n_bmm->setConstant(false);
      n_bsmm->setVal(0);
      n_bsmm->setConstant(true);
      stop_watch_fit.Start(false);
      RooFitResult* nobsmm_res(nullptr);
      if (conditional_fit)
	nobsmm_res = fit_model->fitTo(*data, ConditionalObservables(*mass_err), PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
      else
	nobsmm_res = fit_model->fitTo(*data, PrintEvalErrors(-1), PrintLevel(-1), Save(kTRUE));
      stop_watch_fit.Stop();
      significance_bsmm.push_back(sqrt(max(0.,nobsmm_res->minNll() - def_res->minNll())*2.));
      if (nobsmm_res) delete nobsmm_res;
    }
    delete data;
    if (def_res) delete def_res;
    if (nobmm_res) delete nobmm_res;
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
  ws_ref.Print();

  result.yield_bsmm->Draw();
  print_canvas(plot_name + "_n_bsmm", output_path, gPad);
  result.yield_bmm->Draw();
  print_canvas(plot_name + "_n_bmm", output_path, gPad);
  result.significance_bmm->Draw();
  print_canvas(plot_name + "_significance_bmm", output_path, gPad);
  result.significance_bsmm->Draw();
  print_canvas(plot_name + "_significance_bsmm", output_path, gPad);

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
  string model_name = "model_2D";
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
  print_canvas(model_name + "_mass_err_bsmm", output_path, gPad);
  
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
  print_canvas(model_name + "_bsmm", output_path, gPad);

  cb_alpha_mm.setConstant(true);
  cb_n_mm.setConstant(true);
  s_mm.setConstant(true);
  mean_bsmm.setConstant(true);

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
  print_canvas(model_name + "_mass_err_bmm", output_path, gPad);

  // Make a full PDF
  RooProdPdf pdf2_bmm("pdf2_bmm","pdf2_bmm", pdf_mass_err_bmm, Conditional(pdf_bmm, *mass));

  RooPlot* frame_bmm = mass->frame();
  data_bmm->plotOn(frame_bmm);
  pdf2_bmm.plotOn(frame_bmm);
  // pdf_bmm.plotOn(frame_bmm, ProjWData(*mass_err, *data_bmm));
  frame_bmm->Draw();
  print_canvas(model_name + "_bmm", output_path, gPad);

  mean_bmm.setConstant(true);

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
  print_canvas(model_name + "_mass_err_bkpi", output_path, gPad);

  // Make a full PDF
  RooProdPdf pdf2_bkpi("pdf2_bkpi","pdf2_bkpi", pdf_mass_err_bkpi, Conditional(pdf_bkpi, *mass));

  RooPlot* frame_bkpi = mass->frame();
  data_bkpi->plotOn(frame_bkpi);
  pdf2_bkpi.plotOn(frame_bkpi);
  // pdf_bkpi.plotOn(frame_bkpi, ProjWData(*mass_err, *data_bkpi));
  frame_bkpi->Draw();
  print_canvas(model_name +"_bkpi", output_path, gPad);

  cb_alpha_kpi.setConstant(true);
  cb_n_kpi.setConstant(true);
  s_kpi.setConstant(true);
  mean_bkpi.setConstant(true);

  RooRealVar n_bsmm("n_bsmm", "n_bsmm", get_expected_event_yield(workspace, "bsmm"), 0., 1e9); 
  RooRealVar n_bmm( "n_bmm",  "n_bmm",  get_expected_event_yield(workspace, "bmm"),  0., 1e9); 
  RooRealVar n_bkpi("n_bkpi", "n_bsmm", get_expected_event_yield(workspace, "bkpi"), 0., 1e9); 

  RooAddPdf model(model_name.c_str(), model_name.c_str(), 
		  RooArgList(pdf2_bsmm, pdf2_bmm, pdf2_bkpi), 
		  RooArgList(n_bsmm, n_bmm, n_bkpi));
  
  workspace.import(model);

  // Specify observables
  workspace.defineSet(model_name.c_str(), "mass,mass_err");
  workspace.Print();
  workspace.var("n_bmm")->Print();

  auto frame = mass->frame() ;
  model.plotOn(frame);
  frame->Draw();
  print_canvas(model_name, output_path, gPad);

  // Conditional model
  RooAddPdf model_cond((model_name + "_cond").c_str(), "",
		       RooArgList(pdf_bsmm, pdf_bmm, pdf_bkpi), 
		       RooArgList(n_bsmm, n_bmm, n_bkpi));
  
  workspace.import(model_cond, RenameConflictNodes("_cond"));


}

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
  build_model_2D(workspace);

  // Toy study
  //  auto result = 
  unsigned int n_toys = 10000;
  // toy_study(workspace, "model_1D", "model_1D", n_toys);
  toy_study(workspace, "model_2D", "model_2D_cond", n_toys);
  // toy_study(workspace, "model_2D", "model_2D", n_toys);
  // toy_study(workspace, "model_2D", "model_1D", n_toys);

  return;
}
