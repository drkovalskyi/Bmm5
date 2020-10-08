#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "MRooCBShape.h"
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

// string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/sensitivity_test/";
// string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/sensitivity_kpi_trigger_0.75/";
string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/sensitivity_muon_mva/";

const bool exclude_bhh = true;
const bool exclude_data = true;
const bool do_bsmm_significance = true;
const bool compute_scale_factors = true;    // Compute muon selection efficiency for hadrons
const bool remake_input_workspaces = false;
const bool update_scale_factors = false;     // Ignore scale factors in pre-processed samples and get them from input provided by user in Samples
const bool update_cross_sections = false;
const bool produce_2D_conditional_projections = false;
const bool pre_processing_only = false; // Just produce samples 
// const string final_selection = "mva>0.99&&muonid==3"; // additional selection requirements applied to RooDatasets
const string final_selection = "muonid==3";

const int muon_id = 2; // 1 - loose, 2 - medium, 3 - mva
const double mm_mass_min = 4.95;
const double mm_mass_max = 5.95;
const double mm_mass_blind_min = 5.15;
const double mm_mass_blind_max = 5.50;
const double mm_mass_err_min = 0.005;
const double mm_mass_err_max = 0.100;
// const double mm_min_mva = 0.5;
const double mm_min_mva = -1;
const double trigger_efficiency = 0.6; // applied when sample cannot be used with an explicit trigger requirement

const bool silent_roofit = true;
const bool plot_each_toy = false; // debugging option 

struct Sample{
  string name; 
  vector<string> files;
  float cross_section, scale_factor;
  bool trigger, truth_match, blind, exclusive;
};

const float luminosity = 140e3; // [1/pb]

// string storage_path = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/508/";
string storage_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/508/";
string skim_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD-skims/mm/508/";
    
// Don't use symbols in the name
vector<Sample> samples;
void add_samples(){
  samples.push_back(
			{ "bsmm", 
			  {
			    storage_path + "BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/*.root"
			  },
			  2.98E-02, 1.0, 
			  // trigger, truth_match, blind, exclusive 
			  true,  true, false, true
			});
  samples.push_back(
			{ "bmm", 
			  {
			    storage_path + "BdToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root"
			  },
			  3.26E-03, 1.0,
			  // trigger, truth_match, blind, exclusive 
			  true,  true, false, true
			});
  if (not exclude_data){
    samples.push_back(
		      { "data", 
			{   
			  // Run2016 
			  skim_path + "Charmonium+Run2016B-17Jul2018_ver2-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2016C-17Jul2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2016D-17Jul2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2016E-17Jul2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2016F-17Jul2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2016G-17Jul2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2016H-17Jul2018-v1+MINIAOD/*.root",
			  // Run2017 
			  skim_path + "Charmonium+Run2017B-31Mar2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2017C-31Mar2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2017D-31Mar2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2017E-31Mar2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2017F-31Mar2018-v1+MINIAOD/*.root",
			  // Run2018
			  skim_path + "Charmonium+Run2018A-17Sep2018-v1+MINIAOD/*.root", 
			  skim_path + "Charmonium+Run2018B-17Sep2018-v1+MINIAOD/*.root",
			  skim_path + "Charmonium+Run2018C-17Sep2018-v1+MINIAOD/*.root", 
			  skim_path + "Charmonium+Run2018D-PromptReco-v2+MINIAOD/*.root" 
			},
			1.0, 1.0,
			// trigger, truth_match, blind, exclusive 
			true,  false, true, false
		      });
  }
  if (not exclude_bhh){
    samples.push_back( 
			  { "bkpi", 
			    {
			      storage_path + "BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root"
			    },
			    5.76E+02, 8.64E-06, // muon fakes
			    // trigger, truth_match, blind, exclusive 
			    false,  true, false, true
			  });
  }
}
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


const RooWorkspace* get_workspace(const char* file_name){
  // Check if we already have histograms prepared
  if (not remake_input_workspaces and not gSystem->AccessPathName(file_name)){
    TFile *f = TFile::Open(file_name) ;
    if (not f) 
      throw std::runtime_error( "File is not found" );
    const RooWorkspace* ws = dynamic_cast<const RooWorkspace*>(f->Get("workspace"));
    if (not ws) 
      throw std::runtime_error( "Workspace is not found" );
    return ws;
  }
  return nullptr;
}

void create_workspace(Sample& sample){
  string file_name = "sensitivity_study-" + sample.name + ".root";

  printf("Recreating histograms for %s\n", sample.name.c_str());
    
  RooRealVar mass("mass", "mass", mm_mass_min, mm_mass_max);
  RooRealVar mass_err("mass_err", "mass_err", mm_mass_err_min, mm_mass_err_max);
  RooRealVar mva("mva", "mva", 0.0, 1.0);
  RooCategory muonid("muonid", "muonid");
  // both muons satistify the requirement
  muonid.defineType("mva",    3);
  muonid.defineType("medium", 2);
  muonid.defineType("loose",  1);
  muonid.defineType("failed", 0); 
  RooCategory eta_bin("eta_bin", "eta_bin");
  eta_bin.defineType("barrel", 0);
  eta_bin.defineType("endcap", 1);
  // RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooDataSet data("data", "data", RooArgSet(mass, mass_err, muonid, mva));
  printf("RooDataSet storage type: %u\n", RooAbsData::getDefaultStorageType());
  TChain* chain = new TChain("Events");
  for (auto f: sample.files)
    chain->Add(f.c_str());
  Long64_t n = chain->GetEntries();
  printf("Total number of events: %lld\n", n);

  // Setup variables and branches to read data efficiently
  // mm
  UInt_t     nmm;
  TBranch* b_nmm;
  chain->SetBranchAddress("nmm", &nmm, &b_nmm);
  Int_t      mm_gen_pdgId[200];
  TBranch* b_mm_gen_pdgId;
  if (sample.truth_match)
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
  Float_t    mm_mva[200];
  TBranch* b_mm_mva;
  chain->SetBranchAddress("mm_mva", mm_mva, &b_mm_mva);
  // Muon
  Float_t    muon_softMva[200];
  TBranch* b_muon_softMva;
  chain->SetBranchAddress("Muon_softMva", muon_softMva, &b_muon_softMva);
  Bool_t     muon_mediumId[200];
  TBranch* b_muon_mediumId;
  chain->SetBranchAddress("Muon_mediumId", muon_mediumId, &b_muon_mediumId);
  Bool_t     muon_softMvaId[200];
  TBranch* b_muon_softMvaId;
  chain->SetBranchAddress("Muon_softMvaId", muon_softMvaId, &b_muon_softMvaId);
  // Trigger
  Bool_t     HLT_DoubleMu4_3_Bs;
  TBranch* b_HLT_DoubleMu4_3_Bs;
  chain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
    
  Long64_t n_saved_cands(0);
  Long64_t n_cands_passed_muon_id(0);
  Long64_t n_saved_events(0);
    
  int i_permille_old = 0;
  for (unsigned int i=0; i < n; ++i){
    int i_permille = (int)floor(100. * i / n);
    if (i_permille != i_permille_old) {
      printf("\015\033[32m ---> \033[1m\033[31m%d%%"
             "\033[0m\033[32m <---\033[0m\015", i_permille);
      fflush(stdout);
      i_permille_old = i_permille;
    }
    // Load a proper try and get relative even index
    Long64_t localEntry = chain->LoadTree(i);
    b_HLT_DoubleMu4_3_Bs->GetEntry(localEntry);
    if (sample.trigger && not HLT_DoubleMu4_3_Bs) continue;
    b_nmm->GetEntry(localEntry);
    if (nmm == 0) continue;
    // Get relevant branches
    if (sample.truth_match)
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
    b_mm_mva->GetEntry(localEntry);
    b_muon_softMva->GetEntry(localEntry);
    b_muon_mediumId->GetEntry(localEntry);
    b_muon_softMvaId->GetEntry(localEntry);

    // Loop over candidates
    bool save_event = false;
    for (unsigned int cand=0; cand < nmm; ++cand){
      if (sample.truth_match and not mm_gen_pdgId[cand]) continue;
      if (mm_kin_mass[cand] < mm_mass_min) continue;
      if (mm_kin_mass[cand] > mm_mass_max) continue;
      if (sample.blind){
	if (mm_kin_mass[cand] > mm_mass_blind_min and 
	    mm_kin_mass[cand] < mm_mass_blind_max) continue;
      }
      if (fabs(mm_mu1_eta[cand]) > 1.4) continue;
      if (fabs(mm_mu2_eta[cand]) > 1.4) continue;
      if (mm_mu1_pt[cand] < 4.0) continue;
      if (mm_mu2_pt[cand] < 4.0) continue;
      if (mm_kin_sl3d[cand] < 4.0) continue;
      if (mm_kin_vtx_chi2dof[cand] > 5.0) continue;
      if (TMath::IsNaN(mm_kin_massErr[cand])) continue;
      if (mm_kin_massErr[cand] < mm_mass_err_min) continue;
      if (mm_kin_massErr[cand] > mm_mass_err_max) continue;
      if (mm_mva[cand] < mm_min_mva) continue;
      save_event = true;
      n_saved_cands++;
      mass = mm_kin_mass[cand];
      mass_err = mm_kin_massErr[cand];
      muonid.setIndex(0); 
      if (mm_mu1_index[cand]>=0 and mm_mu2_index[cand]>=0){
	muonid.setIndex(1);
	if (muon_mediumId[mm_mu1_index[cand]] and muon_mediumId[mm_mu2_index[cand]]){
	  muonid.setIndex(2);
	  if (muon_softMvaId[mm_mu1_index[cand]] and muon_softMvaId[mm_mu2_index[cand]]){
	    muonid.setIndex(3);
	  }
	}
      }

      // &&Muon_softMvaId[mm_mu1_index])||(mm_mu2_index>=0
      // and muon_softMva[mm_mu1_index[cand]] > 0.58 and
      // 	  mm_mu2_index[cand]>=0 and muon_softMva[mm_mu2_index[cand]] > 0.58)
      // 	muonid.setIndex(1);
      //       else
      // 	muonid.setIndex(0);
      
      if (muonid.getIndex() >= muon_id)
	n_cands_passed_muon_id++;

      mva = mm_mva[cand];
      data.add(RooArgSet(mass, mass_err, muonid, mva));
    }
    if (save_event) n_saved_events++;
    // if (n_saved_cands > 10000) break;
  }
  printf("Number of saved events: %lld\n", n_saved_events);
  printf("Number of saved candidates: %lld\n", n_saved_cands);
  printf("Number of saved candidates that passed muon id: %lld\n", n_cands_passed_muon_id);
  // data.Print("v");

  // Create a new empty workspace
  RooWorkspace ws("workspace","workspace");
  ws.import(data) ;
  double efficiency = double(n_saved_events)/n;
  if (not sample.trigger)
    efficiency *= trigger_efficiency;
  RooRealVar eff("eff", "Event selection efficiency", efficiency);
  ws.import(eff);
  RooRealVar xsec("xsec", "Effective cross-section before event selection, [pb]", sample.cross_section);
  ws.import(xsec);
  RooRealVar scale("scale", "Correction scale factor", sample.scale_factor);
  if (compute_scale_factors)
    scale.setVal(double(n_cands_passed_muon_id)/n_saved_cands);
  ws.import(scale);
  ws.writeToFile(file_name.c_str()) ;
}

const RooWorkspace* process_sample(Sample& sample){
  string file_name = "sensitivity_study-" + sample.name + ".root";

  auto ws = get_workspace(file_name.c_str());

  if (ws) return ws;

  create_workspace(sample);
  
  return get_workspace(file_name.c_str());
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
      ws.import(*ws_ref.pdf(fit_model_name.c_str()), RenameConflictNodes("_cond"));
    auto mass = ws.var("mass");
    auto mass_err = ws.var("mass_err");
    auto n_bmm = ws.var("n_bmm");
    //  auto n_bkpi = ws.var("n_bkpi");
    // if (exclude_bhh){
    //   n_bkpi->setVal(0);
    //   n_bkpi->setConstant(true);
    // }
    auto n_bsmm = ws.var("n_bsmm");
    auto gen_model = ws.pdf(gen_model_name.c_str());
    auto fit_model = ws.pdf(fit_model_name.c_str());
    auto gen_observables = const_cast<RooWorkspace*>(&ws_ref)->set(gen_model_name.c_str());
    stop_watch_init.Stop();

    stop_watch_generate.Start(false);

    RooDataSet *data = gen_model->generate(*gen_observables, Extended(kTRUE));
    // WARNING:
    // RooCBRooCBShape won't work due to a bug reported at
    // https://sft.its.cern.ch/jira/browse/ROOT-8489
    // use MRooCBRooCBShape

    //// A workaround for the bug - works only for simple PDFs.
    // RooDataSet *data = new RooDataSet("data", "", *gen_observables);
    // RooAcceptReject generator(*gen_model, *gen_observables, *gen_model->getGeneratorConfig(), false, 0);

    // Double_t resample = 0.;
    // int n = 0;
    // for (auto s: samples){
    //   auto n_x = ws.var(("n_" + s.name).c_str());
    //   printf("%s = %f\n", ("n_" + s.name).c_str(), n_x->getVal());
    //   n += n_x->getVal();
    // }
    // printf("n = %d\n", n);
    // for (int i=0; i<n; i++)
    //   data->add(*generator.generateEvent(n-i,resample));

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

  gStyle->SetOptFit();

  result.yield_bsmm->Fit("gaus");
  result.yield_bsmm->Draw();
  print_canvas(plot_name + "_n_bsmm", output_path, gPad);
  result.yield_bmm->Fit("gaus");
  result.yield_bmm->Draw();
  print_canvas(plot_name + "_n_bmm", output_path, gPad);
  result.significance_bmm->Fit("gaus");
  result.significance_bmm->Draw();
  print_canvas(plot_name + "_significance_bmm", output_path, gPad);
  result.significance_bsmm->Fit("gaus");
  result.significance_bsmm->Draw();
  print_canvas(plot_name + "_significance_bsmm", output_path, gPad);

  gStyle->SetOptFit(0);

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
  
  // for (auto& sample: samples){
  for (auto& sample: samples){
    string pdf_name = "pdf_" + sample.name + "_1D";
    if (sample.exclusive){
      auto data = dynamic_cast<const RooDataHist*>(workspace.data(("data_hist_" + sample.name).c_str()));
      if (not data)
	throw std::runtime_error("Cannot get data_hist_" + sample.name);

      // Use 2nd order interpolation to smooth the shape
      RooHistPdf pdf(pdf_name.c_str(), "", *mass, *data, 2);
      workspace.import(pdf);
      auto pdf_ptr = workspace.pdf(pdf_name.c_str());
      pdfs.add(*pdf_ptr);

      float n_expected = get_expected_event_yield(workspace, sample.name);
      const char* n_name = ("n_" + sample.name).c_str();
      RooRealVar* n = new RooRealVar(n_name, "blah", n_expected, 0., 1e9);
      printf("%s = %f\n", n_name, n->getVal());
      yields.add(*n);
    } else {
      auto data = workspace.data(("data_" + sample.name).c_str());
      if (not data)
	throw std::runtime_error("Cannot get data_" + sample.name);
      RooRealVar exp_c_left("exp_c_left","exp_c_left", -1, -1000., 0.);
      RooRealVar exp_c_right("exp_c_right","exp_c_right", -1, -1000., 0.);
      RooRealVar exp_frac("exp_frac","", 0.5, 0, 1);
      RooExponential pdf_left((pdf_name + "_left").c_str(), "Background", *mass, exp_c_left);
      RooExponential pdf_right((pdf_name + "_right").c_str(), "Background", *mass, exp_c_right);
      RooAddPdf pdf(pdf_name.c_str(), "", RooArgList(pdf_left, pdf_right), exp_frac);
      mass->setRange("LeftSideBand", mm_mass_min, mm_mass_blind_min);
      mass->setRange("RightSideBand", mm_mass_blind_max, mm_mass_max);

      pdf_right.fitTo(*data, Range("RightSideBand"));
      exp_c_right.setConstant(true);
      pdf.fitTo(*data, Range("LeftSideBand,RightSideBand"));
      auto fraction = pdf.createIntegral(*mass, *mass, "LeftSideBand,RightSideBand");
      float n_expected = data->numEntries()*fraction->getVal();
      const char* n_name = ("n_" + sample.name).c_str();
      RooRealVar* n = new RooRealVar(n_name, "blah", n_expected, 0., 1e9);
      printf("%s = %f\n", n_name, n->getVal());
      yields.add(*n);
      
      auto frame = mass->frame();
      data->plotOn(frame);
      pdf.plotOn(frame);
      frame->Draw();
      print_canvas(pdf_name, output_path, gPad);
      mass->removeRange("LeftSideBand");
      mass->removeRange("RightSideBand");
      
      workspace.import(pdf);
      auto pdf_ptr = workspace.pdf(pdf_name.c_str());
      pdfs.add(*pdf_ptr);
    }
  }
  
  RooAddPdf model(model_name, model_name, pdfs, yields);
  
  workspace.import(model);
  // Specify observables
  workspace.defineSet(model_name, "mass");
  workspace.Print();

  auto frame = mass->frame() ;

  model.plotOn(frame);
  // } else {
  if (exclude_bhh and workspace.pdf("pdf_bmm_1D")){
    model.plotOn(frame, Components(*(workspace.pdf("pdf_bmm_1D"))), LineStyle(kDashed));
  }
  if (not exclude_bhh and workspace.pdf("pdf_bmm_1D") and workspace.pdf("pdf_bkpi_1D")){
    model.plotOn(frame, Components(*(workspace.pdf("pdf_bkpi_1D"))), FillColor(kMagenta), FillStyle(1001), DrawOption("F"));
    model.plotOn(frame, Components(RooArgSet(*(workspace.pdf("pdf_bmm_1D")), *(workspace.pdf("pdf_bkpi_1D")))), LineStyle(kDashed));
  }

  frame->Draw();

  print_canvas(model_name, output_path, gPad);
}

void build_model_2D(RooWorkspace& workspace){
  string model_name = "model_2D";
  auto mass = workspace.var("mass");
  if (not mass) throw std::runtime_error("Cannot get mass var");
  auto mass_err = workspace.var("mass_err");
  if (not mass_err) throw std::runtime_error("Cannot get mass_err var");

  // Conditional PDFs (mass|mass_err)
  RooArgList cond_pdfs;

  // Full 2D PDFs (mass,mass_err)
  RooArgList full_pdfs;

  RooArgList yields;
  // std::cout << __FILE__ << ":" << __LINE__ << std::endl;

  for (auto sample: samples){
    if (sample.name == "data")
      throw std::runtime_error("2D model is not ready for combinatorial background");
    printf("Building pdfs for %s\n", sample.name.c_str());
    auto data = workspace.data(("data_" + sample.name).c_str());
    if (not data)
      throw std::runtime_error("Cannot get data for " + sample.name);

    // Get Mass Error distribution and make 1D RooHistPdf for it
    RooDataHist data_hist_mass_err(("data_hist_mass_err_" + sample.name).c_str(), "", *mass_err, *data);
    RooHistPdf pdf_mass_err(("pdf_mass_err_" + sample.name).c_str(), "", *mass_err, data_hist_mass_err, 2);
    RooPlot* frame_mass_err = mass_err->frame();
    pdf_mass_err.plotOn(frame_mass_err);
    frame_mass_err->Draw();
    print_canvas(model_name + "_mass_err_" + sample.name, output_path, gPad);
  
    // Build pdf to be used as a conditinal pdf (mass|mass_err)
    // NOTE: it only works for peaking backgrounds around B/Bs mass region
    RooRealVar mean(("mean_" + sample.name).c_str(), "mean", 5.30, 5.20, 5.40);
    RooRealVar s(("s_" + sample.name).c_str(), "resolution", 1, 0, 20);
    RooFormulaVar sigma(("sigma_" + sample.name).c_str(), "", "@0*@1", RooArgList(s, *mass_err));
    RooRealVar cb_alpha(("cb_alpha_" + sample.name).c_str(), "Crystal-Ball gaussian cutoff", 1, 0.1, 10);
    RooRealVar cb_n(("cb_n_" + sample.name).c_str(), "Crystal-Ball exponent of the power-low tail", 5, 0, 10000);
    MRooCBShape cond_pdf(("cond_pdf_" + sample.name).c_str(), "Conditional PDF", *mass, mean, sigma, cb_alpha, cb_n);
    // Fit signal model to extract parameters
    cond_pdf.fitTo(*data, ConditionalObservables(*mass_err), NumCPU(16), Timer(true));
    
    // Fix model parameters
    cb_alpha.setConstant(true);
    cb_n.setConstant(true);
    s.setConstant(true);
    mean.setConstant(true);

    workspace.import(cond_pdf);
    auto cond_pdf_ptr = workspace.pdf(("cond_pdf_" + sample.name).c_str());
    cond_pdfs.add(*cond_pdf_ptr);

    // Make a full 2D PDF
    RooProdPdf full_pdf(("full_pdf_" + sample.name).c_str(), "Full PDF", pdf_mass_err, Conditional(*cond_pdf_ptr, *mass));
    workspace.import(full_pdf);
    full_pdfs.add(*workspace.pdf(("full_pdf_" + sample.name).c_str()));

    // Plot results of the fit
    RooPlot* frame = mass->frame();
    data->plotOn(frame);
    full_pdf.plotOn(frame);
    // pdf_bsmm.plotOn(frame_bsmm, ProjWData(*mass_err, *data_bsmm));
    // pdf_bsmm.plotOn(frame_bsmm, ProjWData(*mass_err, *data_bsmm), NumCPU(16), Normalization(data_bsmm->sumEntries(), RooAbsReal::NumEvent));
    frame->Draw();
    print_canvas(model_name + "_" + sample.name, output_path, gPad);

    // Yield
    RooRealVar* n = new RooRealVar(("n_" + sample.name).c_str(), "n", get_expected_event_yield(workspace, sample.name), 0., 1e9); 

    yields.add(*n);
  }

  // Build models
  RooAddPdf cond_model((model_name + "_cond").c_str(), "Conditional",  cond_pdfs, yields);
  workspace.import(cond_model);

  RooAddPdf full_model((model_name + "_full").c_str(), "Full", full_pdfs, yields);

  // Specify observables
  workspace.defineSet((model_name + "_full").c_str(), "mass,mass_err");

  auto frame = mass->frame() ;
  full_model.plotOn(frame);
  frame->Draw();
  print_canvas(model_name, output_path, gPad);

  workspace.import(full_model, RenameConflictNodes("_full"));
}

void build_bkg_model_1D(RooWorkspace& workspace, const Sample& sample){
  const char* model_name = "bkg_1D";
  RooArgList pdfs;
  RooArgList yields;
  auto mass = workspace.var("mass");
  if (not mass)
    throw std::runtime_error("Cannot get mass");
  
  auto data = workspace.data(("data_" + sample.name).c_str());
  if (not data)
    throw std::runtime_error("Cannot get data_" + sample.name);

  // RooRealVar br1("br1", "", 0.5, 0. , 1);
  // RooRealVar br2("br2", "", 0.5, 0. , 1);
  // RooBernstein bkg("bkg", "Background", *mass, RooArgList(br1, br2));

  RooRealVar exp_c("exp_c","exp_c", -1, -1000., 0.);
  RooExponential bkg(model_name, "Background", *mass, exp_c);
  // RooRealVar n_data("n_data", "", 1, 0, 1e9) ; 
  // RooExtendPdf ebkg("ebkg","ebkg", bkg, n_data) ;  
  mass->setRange("LeftSideBand", mm_mass_min, mm_mass_blind_min);
  mass->setRange("RightSideBand", mm_mass_blind_max, mm_mass_max);

  bkg.fitTo(*data, Range("LeftSideBand,RightSideBand"));
  printf("fraction: %0.3f\n", bkg.createIntegral(*mass, *mass, "LeftSideBand,RightSideBand")->getVal()) ;
  
  auto frame = mass->frame();
  data->plotOn(frame);
  bkg.plotOn(frame);
  frame->Draw();
  print_canvas(model_name, output_path, gPad);
  workspace.import(bkg);

//     // Use 2nd order interpolation to smooth the shape
//     RooHistPdf* pdf = new RooHistPdf(("pdf_" + sample.name + "_1D").c_str(), "", *mass, *data, 2);
//     pdf->Print();
//     pdfs.add(*pdf);
//     printf("pdfs.getSize: %u\n", pdfs.getSize());
//     float n_expected = get_expected_event_yield(workspace, sample.name);
//     RooRealVar* n = new RooRealVar(("n_" + sample.name).c_str(), "blah", n_expected, 0., 1e9);
//     n->Print();
//     yields.add(*n);
//     printf("yields.getSize: %u\n", yields.getSize());
//   }
  
//   RooAddPdf model(model_name, model_name, pdfs, yields);
  
//   workspace.import(model);
//   // Specify observables
//   workspace.defineSet(model_name, "mass");
//   workspace.Print();

//   auto frame = mass->frame() ;

//   model.plotOn(frame);
//   // } else {
//   if (exclude_bhh and workspace.pdf("pdf_bmm_1D")){
//     model.plotOn(frame, Components(*(workspace.pdf("pdf_bmm_1D"))), LineStyle(kDashed));
//   }
//   if (not exclude_bhh and workspace.pdf("pdf_bmm_1D") and workspace.pdf("pdf_bkpi_1D")){
//     model.plotOn(frame, Components(*(workspace.pdf("pdf_bkpi_1D"))), FillColor(kMagenta), FillStyle(1001), DrawOption("F"));
//     model.plotOn(frame, Components(RooArgSet(*(workspace.pdf("pdf_bmm_1D")), *(workspace.pdf("pdf_bkpi_1D")))), LineStyle(kDashed));
//   }

//   frame->Draw();

//   print_canvas(model_name, output_path, gPad);
}

void import_var(RooWorkspace& target_ws, const RooWorkspace& source_ws, const char* var_name, const char* new_name, float scale = 1.0){
  auto var = source_ws.var(var_name);
    if (not var) 
      throw std::runtime_error(string(var_name) +" is missing");
    RooRealVar new_var(*var);
    new_var.SetName(new_name);
    new_var.setVal(new_var.getVal()*scale);
    target_ws.import(new_var);
}

void sensitivity_study(){
  add_samples();

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  RooWorkspace workspace("study_workspace","");
  for (auto& sample: samples){
    // Process sample
    const RooWorkspace* sample_workspace = process_sample(sample);
    if (pre_processing_only) continue;

    // Import sample workspace into the main workspace
    auto sample_data = sample_workspace->data("data");
    if (not sample_data) 
      throw std::runtime_error("RooDataSet \"data\" is missing");
    double n_0 = sample_data->numEntries();
    if (final_selection != "")
      sample_data = sample_data->reduce(final_selection.c_str());

    workspace.import(*sample_data, Rename(("data_" + sample.name).c_str()));

    import_var(workspace, *sample_workspace, "eff",   ("eff_"   + sample.name).c_str(), sample_data->numEntries()/n_0);
    import_var(workspace, *sample_workspace, "xsec",  ("xsec_"  + sample.name).c_str());
    import_var(workspace, *sample_workspace, "scale", ("scale_" + sample.name).c_str());

    RooRealVar* mass = sample_workspace->var("mass");
    if (mass){
      RooPlot* frame = mass->frame() ;
      sample_data->plotOn(frame);
      frame->Draw();
      print_canvas("mass_" + sample.name, output_path, c1);

      RooDataHist rdh(("data_hist_" + sample.name).c_str(), "", *workspace.var("mass"), *sample_data);
      workspace.import(rdh);
    }

    RooRealVar* mass_err = sample_workspace->var("mass_err");
    if (mass_err){
      RooPlot* frame = mass_err->frame() ;
      sample_data->plotOn(frame);
      frame->Draw();
      print_canvas("mass_err_" + sample.name, output_path, c1);
    }


    // TTree* tree = const_cast<TTree*>(sample_data->tree());
    // if (tree){
    //   TH1D h_mass(("h_mass_" + sample.name).c_str(), "", 100, mm_mass_min, mm_mass_max);
    //   tree->Draw(("mass>>h_mass_" + sample.name).c_str());
    //   RooDataHist rdh(("data_hist_" + sample.name).c_str(), "", *workspace.var("mass"), &h_mass);
    //   workspace.import(rdh);
    //   print_canvas("mass_" + sample.name, output_path, c1);

    //   tree->Draw("mass_err");
    //   print_canvas("mass_err_" + sample.name, output_path, c1);
    // } else {
    //   if (RooRealVar* mass = sample_workspace->var("mass")){
    // 	RooPlot* frame = mass->frame() ;
    // 	sample_workspace->data("data")->plotOn(frame);
    // 	frame->Draw();
    // 	print_canvas("mass_" + sample.name, output_path, c1);
    //   }
    // }
  }
  if (pre_processing_only) return;

  build_model_1D(workspace);

  build_model_2D(workspace);

  workspace.Print("V");

  // Toy studies
  unsigned int n_toys = 1000;
  toy_study(workspace, "model_1D", "model_1D", n_toys);
  toy_study(workspace, "model_2D_full", "model_2D_cond", n_toys);
  toy_study(workspace, "model_2D_full", "model_1D", n_toys);
  // // toy_study(workspace, "model_2D", "model_2D", n_toys);


}
				//                   
