#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"

/////////////////////////////////////////////////////
//
//  Objective: estimate hadron to muon fake rate
//
//  Algorithm:
//    - loop over all hadrons in genbmm block and match them by dR to Muons
//    - measure fake rate for different selectors and hadron types
//
////////////////////////////////////////////////////

using namespace std;

string output_path = "/afs/cern.ch/user/d/dmytro/www/public_html/plots/bmm5_NanoAODv6-508/muon_fakerate_test/";


// Study parameters

struct Sample{
  string name; 
  vector<string> files;
};

string storage_path = "/eos/cms/store/group/phys_bphys/bmm/bmm5/NanoAOD/508/";

// Don't use symbols in the name
vector<Sample> samples = {
  {
    "bkpi", 
    { storage_path + "BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root" }
    // storage_path + "BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/D1F627D2-EEAD-C142-AC31-4575EB6E0416.root"
  },
  {
    "bpipi", 
    { storage_path + "BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root" }
  },
  {
    "soup",
    {
      storage_path + "BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root",
      storage_path + "BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root",
      storage_path + "BsToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1+MINIAODSIM/*.root",
      storage_path + "BdToPiMuNu_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root",
      storage_path + "LambdaBToPPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root",
      storage_path + "BdToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root",
      storage_path + "BsToKK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/*.root"
    }
  }

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

int match_muon(const Float_t& had_pt, const Float_t& had_eta, const Float_t& had_phi,
	       const UInt_t& nmuons, const Float_t* muon_pt, const Float_t* muon_eta,
	       const Float_t* muon_phi)
{
  for (unsigned int i=0; i < nmuons; ++i){
    double dr2 = (had_eta - muon_eta[i])*(had_eta - muon_eta[i]) + pow(acos(cos(had_phi - muon_phi[i])),2);
    if (dr2 < 0) continue;
    if (sqrt(dr2) > 0.1) continue;
    if ( fabs(muon_pt[i]/had_pt - 1) > 0.03 ) continue;
    return i;
  }
  return -1;
}

void process_sample(Sample& sample){
  printf("Processing %s\n", sample.name.c_str());

  TChain* chain = new TChain("Events");
  for (auto f: sample.files)
    chain->Add(f.c_str());
  Long64_t n = chain->GetEntries();
  printf("Total number of events: %lld\n", n);

  // Setup variables and branches to read data efficiently
  // genbmm
  UInt_t     ngenbmm;
  TBranch* b_ngenbmm;
  chain->SetBranchAddress("ngenbmm", &ngenbmm, &b_ngenbmm);
  Int_t      genbmm_mu1_pdgId[200];
  TBranch* b_genbmm_mu1_pdgId;
  chain->SetBranchAddress("genbmm_mu1_pdgId", genbmm_mu1_pdgId, &b_genbmm_mu1_pdgId);
  Int_t      genbmm_mu2_pdgId[200];
  TBranch* b_genbmm_mu2_pdgId;
  chain->SetBranchAddress("genbmm_mu2_pdgId", genbmm_mu2_pdgId, &b_genbmm_mu2_pdgId);
  Float_t    genbmm_mu1_pt[200];
  TBranch* b_genbmm_mu1_pt;
  chain->SetBranchAddress("genbmm_mu1_pt", genbmm_mu1_pt, &b_genbmm_mu1_pt);
  Float_t    genbmm_mu2_pt[200];
  TBranch* b_genbmm_mu2_pt;
  chain->SetBranchAddress("genbmm_mu2_pt", genbmm_mu2_pt, &b_genbmm_mu2_pt);
  Float_t    genbmm_mu1_eta[200];
  TBranch* b_genbmm_mu1_eta;
  chain->SetBranchAddress("genbmm_mu1_eta", genbmm_mu1_eta, &b_genbmm_mu1_eta);
  Float_t    genbmm_mu2_eta[200];
  TBranch* b_genbmm_mu2_eta;
  chain->SetBranchAddress("genbmm_mu2_eta", genbmm_mu2_eta, &b_genbmm_mu2_eta);
  Float_t    genbmm_mu1_phi[200];
  TBranch* b_genbmm_mu1_phi;
  chain->SetBranchAddress("genbmm_mu1_phi", genbmm_mu1_phi, &b_genbmm_mu1_phi);
  Float_t    genbmm_mu2_phi[200];
  TBranch* b_genbmm_mu2_phi;
  chain->SetBranchAddress("genbmm_mu2_phi", genbmm_mu2_phi, &b_genbmm_mu2_phi);

  // Muon
  UInt_t     nmuons;
  TBranch* b_nmuons;
  chain->SetBranchAddress("nMuon", &nmuons, &b_nmuons);
  Float_t    muon_pt[200];
  TBranch* b_muon_pt;
  chain->SetBranchAddress("Muon_pt", muon_pt, &b_muon_pt);
  Float_t    muon_eta[200];
  TBranch* b_muon_eta;
  chain->SetBranchAddress("Muon_eta", muon_eta, &b_muon_eta);
  Float_t    muon_phi[200];
  TBranch* b_muon_phi;
  chain->SetBranchAddress("Muon_phi", muon_phi, &b_muon_phi);
  Float_t    muon_softMva[200];
  TBranch* b_muon_softMva;
  chain->SetBranchAddress("Muon_softMva", muon_softMva, &b_muon_softMva);
  Bool_t     muon_mediumId[200];
  TBranch* b_muon_mediumId;
  chain->SetBranchAddress("Muon_mediumId", muon_mediumId, &b_muon_mediumId);
  Bool_t     muon_softMvaId[200];
  TBranch* b_muon_softMvaId;
  chain->SetBranchAddress("Muon_softMvaId", muon_softMvaId, &b_muon_softMvaId);
  Bool_t     muon_isTracker[200];
  TBranch* b_muon_isTracker;
  chain->SetBranchAddress("Muon_isTracker", muon_isTracker, &b_muon_isTracker);
  Bool_t     muon_isGlobal[200];
  TBranch* b_muon_isGlobal;
  chain->SetBranchAddress("Muon_isGlobal", muon_isGlobal, &b_muon_isGlobal);

  Double_t pt_bins[5] = {4.0, 6.0, 10.0, 18.0, 34.0};
  Double_t eta_bins[4] = {0.0, 0.4, 0.9, 1.4};

  TH1D h_pion_pt_all("h_pion_pt_all;Pt", "All pions from relevant B decays", 4, pt_bins);
  h_pion_pt_all.Sumw2();
  h_pion_pt_all.SetDirectory(0);
  TH1D h_pion_eta_all("h_pion_eta_all;Pt", "All pions from relevant B decays", 3, eta_bins);
  h_pion_eta_all.Sumw2();
  h_pion_eta_all.SetDirectory(0);

  TH1D h_pion_pt_trigger("h_pion_pt_trigger;Pt", "Pions that passed Tracker&Global", 4, pt_bins);
  h_pion_pt_trigger.Sumw2();
  h_pion_pt_trigger.SetDirectory(0);
  TH1D h_pion_pt_trigger_fake_rate("h_pion_pt_trigger_fake_rate;Pt;Fake Rate", "Pion fake rate for Tracker&Global", 4, pt_bins);
  h_pion_pt_trigger_fake_rate.Sumw2();
  h_pion_pt_trigger_fake_rate.SetDirectory(0);

  TH1D h_pion_eta_trigger("h_pion_eta_trigger;Pt", "Pions that passed Tracker&Global", 3, eta_bins);
  h_pion_eta_trigger.Sumw2();
  h_pion_eta_trigger.SetDirectory(0);
  TH1D h_pion_eta_trigger_fake_rate("h_pion_trigger_fake_rate;Pt;Fake Rate", "Pion fake rate for Tracker&Global", 3, eta_bins);
  h_pion_eta_trigger_fake_rate.Sumw2();
  h_pion_eta_trigger_fake_rate.SetDirectory(0);
  
  TH1D h_pion_pt_mediumMuId("h_pion_pt_mediumMuId;Pt", "Pions that passed medium muon ID", 4, pt_bins);
  h_pion_pt_mediumMuId.Sumw2();
  h_pion_pt_mediumMuId.SetDirectory(0);
  TH1D h_pion_pt_mediumMuId_fake_rate("h_pion_pt_mediumMuId_fake_rate;Pt;Fake Rate", "Pion fake rate for medium muon ID", 4, pt_bins);
  h_pion_pt_mediumMuId_fake_rate.Sumw2();
  h_pion_pt_mediumMuId_fake_rate.SetDirectory(0);

  TH1D h_pion_eta_mediumMuId("h_pion_eta_mediumMuId;Pt", "Pions that passed medium muon ID", 3, eta_bins);
  h_pion_eta_mediumMuId.Sumw2();
  h_pion_eta_mediumMuId.SetDirectory(0);
  TH1D h_pion_eta_mediumMuId_fake_rate("h_pion_eta_mediumMuId_fake_rate;Pt;Fake Rate", "Pion fake rate for medium muon ID", 3, eta_bins);
  h_pion_eta_mediumMuId_fake_rate.Sumw2();
  h_pion_eta_mediumMuId_fake_rate.SetDirectory(0);
  
  TH1D h_pion_pt_softMvaId("h_pion_pt_softMvaId;Pt", "Pions that passed soft MVA muon ID", 4, pt_bins);
  h_pion_pt_softMvaId.Sumw2();
  h_pion_pt_softMvaId.SetDirectory(0);
  TH1D h_pion_pt_softMvaId_fake_rate("h_pion_pt_softMvaId_fake_rate;Pt;Fake Rate", "Pion fake rate for soft MVA muon ID", 4, pt_bins);
  h_pion_pt_softMvaId_fake_rate.Sumw2();
  h_pion_pt_softMvaId_fake_rate.SetDirectory(0);

  TH1D h_pion_eta_softMvaId("h_pion_eta_softMvaId;Pt", "Pions that passed soft MVA muon ID", 3, eta_bins);
  h_pion_eta_softMvaId.Sumw2();
  h_pion_eta_softMvaId.SetDirectory(0);
  TH1D h_pion_eta_softMvaId_fake_rate("h_pion_eta_softMvaId_fake_rate;Pt;Fake Rate", "Pion fake rate for soft MVA muon ID", 3, eta_bins);
  h_pion_eta_softMvaId_fake_rate.Sumw2();
  h_pion_eta_softMvaId_fake_rate.SetDirectory(0);

  TH1D h_kaon_pt_all("h_kaon_pt_all;Pt", "All kaons from relevant B decays", 4, pt_bins);
  h_kaon_pt_all.Sumw2();
  h_kaon_pt_all.SetDirectory(0);
  TH1D h_kaon_eta_all("h_kaon_eta_all;Pt", "All kaons from relevant B decays", 3, eta_bins);
  h_kaon_eta_all.Sumw2();
  h_kaon_eta_all.SetDirectory(0);

  TH1D h_kaon_pt_trigger("h_kaon_pt_trigger;Pt", "Kaons that passed Tracker&Global", 4, pt_bins);
  h_kaon_pt_trigger.Sumw2();
  h_kaon_pt_trigger.SetDirectory(0);
  TH1D h_kaon_pt_trigger_fake_rate("h_kaon_pt_trigger_fake_rate;Pt;Fake Rate", "Kaon fake rate for Tracker&Global", 4, pt_bins);
  h_kaon_pt_trigger_fake_rate.Sumw2();
  h_kaon_pt_trigger_fake_rate.SetDirectory(0);

  TH1D h_kaon_eta_trigger("h_kaon_eta_trigger;Pt", "Kaons that passed Tracker&Global", 3, eta_bins);
  h_kaon_eta_trigger.Sumw2();
  h_kaon_eta_trigger.SetDirectory(0);
  TH1D h_kaon_eta_trigger_fake_rate("h_kaon_trigger_fake_rate;Pt;Fake Rate", "Kaon fake rate for Tracker&Global", 3, eta_bins);
  h_kaon_eta_trigger_fake_rate.Sumw2();
  h_kaon_eta_trigger_fake_rate.SetDirectory(0);
  
  TH1D h_kaon_pt_mediumMuId("h_kaon_pt_mediumMuId;Pt", "Kaons that passed medium muon ID", 4, pt_bins);
  h_kaon_pt_mediumMuId.Sumw2();
  h_kaon_pt_mediumMuId.SetDirectory(0);
  TH1D h_kaon_pt_mediumMuId_fake_rate("h_kaon_pt_mediumMuId_fake_rate;Pt;Fake Rate", "Kaon fake rate for medium muon ID", 4, pt_bins);
  h_kaon_pt_mediumMuId_fake_rate.Sumw2();
  h_kaon_pt_mediumMuId_fake_rate.SetDirectory(0);

  TH1D h_kaon_eta_mediumMuId("h_kaon_eta_mediumMuId;Pt", "Kaons that passed medium muon ID", 3, eta_bins);
  h_kaon_eta_mediumMuId.Sumw2();
  h_kaon_eta_mediumMuId.SetDirectory(0);
  TH1D h_kaon_eta_mediumMuId_fake_rate("h_kaon_eta_mediumMuId_fake_rate;Pt;Fake Rate", "Kaon fake rate for medium muon ID", 3, eta_bins);
  h_kaon_eta_mediumMuId_fake_rate.Sumw2();
  h_kaon_eta_mediumMuId_fake_rate.SetDirectory(0);
  
  TH1D h_kaon_pt_softMvaId("h_kaon_pt_softMvaId;Pt", "Kaons that passed soft MVA muon ID", 4, pt_bins);
  h_kaon_pt_softMvaId.Sumw2();
  h_kaon_pt_softMvaId.SetDirectory(0);
  TH1D h_kaon_pt_softMvaId_fake_rate("h_kaon_pt_softMvaId_fake_rate;Pt;Fake Rate", "Kaon fake rate for soft MVA muon ID", 4, pt_bins);
  h_kaon_pt_softMvaId_fake_rate.Sumw2();
  h_kaon_pt_softMvaId_fake_rate.SetDirectory(0);

  TH1D h_kaon_eta_softMvaId("h_kaon_eta_softMvaId;Pt", "Kaons that passed soft MVA muon ID", 3, eta_bins);
  h_kaon_eta_softMvaId.Sumw2();
  h_kaon_eta_softMvaId.SetDirectory(0);
  TH1D h_kaon_eta_softMvaId_fake_rate("h_kaon_eta_softMvaId_fake_rate;Pt;Fake Rate", "Kaon fake rate for soft MVA muon ID", 3, eta_bins);
  h_kaon_eta_softMvaId_fake_rate.Sumw2();
  h_kaon_eta_softMvaId_fake_rate.SetDirectory(0);
  
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

    b_ngenbmm->GetEntry(localEntry);
    if (ngenbmm == 0) continue;
    // Get relevant branches
    b_genbmm_mu1_pdgId->GetEntry(localEntry);
    b_genbmm_mu2_pdgId->GetEntry(localEntry);
    b_genbmm_mu1_pt->GetEntry(localEntry);
    b_genbmm_mu2_pt->GetEntry(localEntry);
    b_genbmm_mu1_eta->GetEntry(localEntry);
    b_genbmm_mu2_eta->GetEntry(localEntry);
    b_genbmm_mu1_phi->GetEntry(localEntry);
    b_genbmm_mu2_phi->GetEntry(localEntry);
    
    b_nmuons->GetEntry(localEntry);
    b_muon_pt->GetEntry(localEntry);
    b_muon_eta->GetEntry(localEntry);
    b_muon_phi->GetEntry(localEntry);
    b_muon_softMva->GetEntry(localEntry);
    b_muon_mediumId->GetEntry(localEntry);
    b_muon_softMvaId->GetEntry(localEntry);
    b_muon_isTracker->GetEntry(localEntry);
    b_muon_isGlobal->GetEntry(localEntry);

    // Loop over hadrons that may fake muons
    
    for (unsigned int bhad=0; bhad < ngenbmm; ++bhad){
      if (fabs(genbmm_mu1_eta[bhad]) < 1.4 and genbmm_mu1_pt[bhad] > 4.0)
	{
	  if (abs(genbmm_mu1_pdgId[bhad]) == 211 or abs(genbmm_mu1_pdgId[bhad]) == 321){
	    // look for a matching muon
	    int imatch = match_muon(genbmm_mu1_pt[bhad], genbmm_mu1_eta[bhad], genbmm_mu1_phi[bhad],
				    nmuons, muon_pt, muon_eta, muon_phi);
	    if (abs(genbmm_mu1_pdgId[bhad]) == 211) // Pion
	      {
		h_pion_pt_all.Fill(genbmm_mu1_pt[bhad]);
		h_pion_eta_all.Fill(fabs(genbmm_mu1_eta[bhad]));
		if (imatch >= 0){
		  if (muon_isGlobal[imatch] and muon_isTracker[imatch]){
		    h_pion_pt_trigger.Fill(genbmm_mu1_pt[bhad]);
		    h_pion_eta_trigger.Fill(fabs(genbmm_mu1_eta[bhad]));
		  }
		  if (muon_mediumId[imatch]){
		    h_pion_pt_mediumMuId.Fill(genbmm_mu1_pt[bhad]);
		    h_pion_eta_mediumMuId.Fill(fabs(genbmm_mu1_eta[bhad]));
		  }
		  if (muon_softMvaId[imatch]){
		    h_pion_pt_softMvaId.Fill(genbmm_mu1_pt[bhad]);
		    h_pion_eta_softMvaId.Fill(fabs(genbmm_mu1_eta[bhad]));
		  }
		}
	      }
	    if (abs(genbmm_mu1_pdgId[bhad]) == 321) // Kaon
	      {
		h_kaon_pt_all.Fill(genbmm_mu1_pt[bhad]);
		h_kaon_eta_all.Fill(fabs(genbmm_mu1_eta[bhad]));
		if (imatch >= 0){
		  if (muon_isGlobal[imatch] and muon_isTracker[imatch]){
		    h_kaon_pt_trigger.Fill(genbmm_mu1_pt[bhad]);
		    h_kaon_eta_trigger.Fill(fabs(genbmm_mu1_eta[bhad]));
		  }
		  if (muon_mediumId[imatch]){
		    h_kaon_pt_mediumMuId.Fill(genbmm_mu1_pt[bhad]);
		    h_kaon_eta_mediumMuId.Fill(fabs(genbmm_mu1_eta[bhad]));
		  }
		  if (muon_softMvaId[imatch]){
		    h_kaon_pt_softMvaId.Fill(genbmm_mu1_pt[bhad]);
		    h_kaon_eta_softMvaId.Fill(fabs(genbmm_mu1_eta[bhad]));
		  }
		}
	      }
	  }
	}
      if (fabs(genbmm_mu2_eta[bhad]) < 1.4 and genbmm_mu2_pt[bhad] > 4.0)
	{
	  if (abs(genbmm_mu2_pdgId[bhad]) == 211 or abs(genbmm_mu2_pdgId[bhad]) == 321){
	    // look for a matching muon
	    int imatch = match_muon(genbmm_mu2_pt[bhad], genbmm_mu2_eta[bhad], genbmm_mu2_phi[bhad],
				    nmuons, muon_pt, muon_eta, muon_phi);
	    if (abs(genbmm_mu2_pdgId[bhad]) == 211) // Pion
	      {
		h_pion_pt_all.Fill(genbmm_mu2_pt[bhad]);
		h_pion_eta_all.Fill(fabs(genbmm_mu2_eta[bhad]));
		if (imatch >= 0){
		  if (muon_isGlobal[imatch] and muon_isTracker[imatch]){
		    h_pion_pt_trigger.Fill(genbmm_mu2_pt[bhad]);
		    h_pion_eta_trigger.Fill(fabs(genbmm_mu2_eta[bhad]));
		  }
		  if (muon_mediumId[imatch]){
		    h_pion_pt_mediumMuId.Fill(genbmm_mu2_pt[bhad]);
		    h_pion_eta_mediumMuId.Fill(fabs(genbmm_mu2_eta[bhad]));
		  }
		  if (muon_softMvaId[imatch]){
		    h_pion_pt_softMvaId.Fill(genbmm_mu2_pt[bhad]);
		    h_pion_eta_softMvaId.Fill(fabs(genbmm_mu2_eta[bhad]));
		  }
		}
	      }
	    if (abs(genbmm_mu2_pdgId[bhad]) == 321) // Kaon
	      {
		h_kaon_pt_all.Fill(genbmm_mu2_pt[bhad]);
		h_kaon_eta_all.Fill(fabs(genbmm_mu2_eta[bhad]));
		if (imatch >= 0){
		  if (muon_isGlobal[imatch] and muon_isTracker[imatch]){
		    h_kaon_pt_trigger.Fill(genbmm_mu2_pt[bhad]);
		    h_kaon_eta_trigger.Fill(fabs(genbmm_mu2_eta[bhad]));
		  }
		  if (muon_mediumId[imatch]){
		    h_kaon_pt_mediumMuId.Fill(genbmm_mu2_pt[bhad]);
		    h_kaon_eta_mediumMuId.Fill(fabs(genbmm_mu2_eta[bhad]));
		  }
		  if (muon_softMvaId[imatch]){
		    h_kaon_pt_softMvaId.Fill(genbmm_mu2_pt[bhad]);
		    h_kaon_eta_softMvaId.Fill(fabs(genbmm_mu2_eta[bhad]));
		  }
		}
	      }
	  }
	}
    }
  }
  h_pion_pt_all.Draw();
  print_canvas(sample.name + "_pion_pt_all", output_path, gPad);
  h_pion_pt_trigger.Draw();
  print_canvas(sample.name + "_pion_pt_trigger", output_path, gPad);
  h_pion_pt_trigger.Draw();
  print_canvas(sample.name + "_pion_pt_mediumMuId", output_path, gPad);
  h_pion_pt_softMvaId.Draw();
  print_canvas(sample.name + "_pion_pt_softMvaId", output_path, gPad);
  h_pion_eta_all.Draw();
  print_canvas(sample.name + "_pion_eta_all", output_path, gPad);
  h_pion_eta_trigger.Draw();
  print_canvas(sample.name + "_pion_eta_trigger", output_path, gPad);
  h_pion_eta_trigger.Draw();
  print_canvas(sample.name + "_pion_eta_mediumMuId", output_path, gPad);
  h_pion_eta_softMvaId.Draw();
  print_canvas(sample.name + "_pion_eta_softMvaId", output_path, gPad);

  h_pion_pt_trigger_fake_rate.Divide(&h_pion_pt_trigger, &h_pion_pt_all, 1, 1);
  h_pion_pt_trigger_fake_rate.SetMinimum(0);
  h_pion_pt_trigger_fake_rate.SetMaximum(0.002);
  h_pion_pt_trigger_fake_rate.SetMarkerStyle(20);
  h_pion_pt_trigger_fake_rate.SetMarkerColor(kBlack);
  h_pion_pt_trigger_fake_rate.GetXaxis()->SetTitle("Pt");
  h_pion_pt_trigger_fake_rate.Draw("e0");
  print_canvas(sample.name + "_pion_pt_trigger_fake_rate", output_path, gPad);
  h_pion_eta_trigger_fake_rate.Divide(&h_pion_eta_trigger, &h_pion_eta_all, 1, 1);
  h_pion_eta_trigger_fake_rate.SetMinimum(0);
  h_pion_eta_trigger_fake_rate.SetMarkerStyle(20);
  h_pion_eta_trigger_fake_rate.SetMarkerColor(kBlack);
  h_pion_eta_trigger_fake_rate.GetXaxis()->SetTitle("|#eta|");
  h_pion_eta_trigger_fake_rate.Draw("e0");
  print_canvas(sample.name + "_pion_eta_trigger_fake_rate", output_path, gPad);

  h_pion_pt_mediumMuId_fake_rate.Divide(&h_pion_pt_mediumMuId, &h_pion_pt_all, 1, 1);
  h_pion_pt_mediumMuId_fake_rate.SetMinimum(0);
  h_pion_pt_mediumMuId_fake_rate.SetMaximum(0.002);
  h_pion_pt_mediumMuId_fake_rate.SetMarkerStyle(20);
  h_pion_pt_mediumMuId_fake_rate.SetMarkerColor(kBlue);
  h_pion_pt_mediumMuId_fake_rate.GetXaxis()->SetTitle("Pt");
  h_pion_pt_mediumMuId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_pion_pt_mediumMuId_fake_rate", output_path, gPad);
  h_pion_eta_mediumMuId_fake_rate.Divide(&h_pion_eta_mediumMuId, &h_pion_eta_all, 1, 1);
  h_pion_eta_mediumMuId_fake_rate.SetMinimum(0);
  h_pion_eta_mediumMuId_fake_rate.SetMarkerStyle(20);
  h_pion_eta_mediumMuId_fake_rate.SetMarkerColor(kBlue);
  h_pion_eta_mediumMuId_fake_rate.GetXaxis()->SetTitle("|#eta|");
  h_pion_eta_mediumMuId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_pion_eta_mediumMuId_fake_rate", output_path, gPad);

  h_pion_pt_softMvaId_fake_rate.Divide(&h_pion_pt_softMvaId, &h_pion_pt_all, 1, 1);
  h_pion_pt_softMvaId_fake_rate.SetMinimum(0);
  h_pion_pt_softMvaId_fake_rate.SetMaximum(0.002);
  h_pion_pt_softMvaId_fake_rate.SetMarkerStyle(20);
  h_pion_pt_softMvaId_fake_rate.SetMarkerColor(kRed);
  h_pion_pt_softMvaId_fake_rate.GetXaxis()->SetTitle("Pt");
  h_pion_pt_softMvaId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_pion_pt_softMvaId_fake_rate", output_path, gPad);
  h_pion_eta_softMvaId_fake_rate.Divide(&h_pion_eta_softMvaId, &h_pion_eta_all, 1, 1);
  h_pion_eta_softMvaId_fake_rate.SetMinimum(0);
  h_pion_eta_softMvaId_fake_rate.SetMarkerStyle(20);
  h_pion_eta_softMvaId_fake_rate.SetMarkerColor(kRed);
  h_pion_eta_softMvaId_fake_rate.GetXaxis()->SetTitle("|#eta|");
  h_pion_eta_softMvaId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_pion_eta_softMvaId_fake_rate", output_path, gPad);

  h_kaon_pt_all.Draw();
  print_canvas(sample.name + "_kaon_pt_all", output_path, gPad);
  h_kaon_pt_trigger.Draw();
  print_canvas(sample.name + "_kaon_pt_trigger", output_path, gPad);
  h_kaon_pt_trigger.Draw();
  print_canvas(sample.name + "_kaon_pt_mediumMuId", output_path, gPad);
  h_kaon_pt_softMvaId.Draw();
  print_canvas(sample.name + "_kaon_pt_softMvaId", output_path, gPad);
  h_kaon_eta_all.Draw();
  print_canvas(sample.name + "_kaon_eta_all", output_path, gPad);
  h_kaon_eta_trigger.Draw();
  print_canvas(sample.name + "_kaon_eta_trigger", output_path, gPad);
  h_kaon_eta_trigger.Draw();
  print_canvas(sample.name + "_kaon_eta_mediumMuId", output_path, gPad);
  h_kaon_eta_softMvaId.Draw();
  print_canvas(sample.name + "_kaon_eta_softMvaId", output_path, gPad);

  h_kaon_pt_trigger_fake_rate.Divide(&h_kaon_pt_trigger, &h_kaon_pt_all, 1, 1);
  h_kaon_pt_trigger_fake_rate.SetMinimum(0);
  h_kaon_pt_trigger_fake_rate.SetMaximum(0.0032);
  h_kaon_pt_trigger_fake_rate.SetMarkerStyle(20);
  h_kaon_pt_trigger_fake_rate.SetMarkerColor(kBlack);
  h_kaon_pt_trigger_fake_rate.GetXaxis()->SetTitle("Pt");
  h_kaon_pt_trigger_fake_rate.Draw("e0");
  print_canvas(sample.name + "_kaon_pt_trigger_fake_rate", output_path, gPad);
  h_kaon_eta_trigger_fake_rate.Divide(&h_kaon_eta_trigger, &h_kaon_eta_all, 1, 1);
  h_kaon_eta_trigger_fake_rate.SetMinimum(0);
  h_kaon_eta_trigger_fake_rate.SetMarkerStyle(20);
  h_kaon_eta_trigger_fake_rate.SetMarkerColor(kBlack);
  h_kaon_eta_trigger_fake_rate.GetXaxis()->SetTitle("|#eta|");
  h_kaon_eta_trigger_fake_rate.Draw("e0");
  print_canvas(sample.name + "_kaon_eta_trigger_fake_rate", output_path, gPad);

  h_kaon_pt_mediumMuId_fake_rate.Divide(&h_kaon_pt_mediumMuId, &h_kaon_pt_all, 1, 1);
  h_kaon_pt_mediumMuId_fake_rate.SetMinimum(0);
  h_kaon_pt_mediumMuId_fake_rate.SetMaximum(0.0032);
  h_kaon_pt_mediumMuId_fake_rate.SetMarkerStyle(20);
  h_kaon_pt_mediumMuId_fake_rate.SetMarkerColor(kBlue);
  h_kaon_pt_mediumMuId_fake_rate.GetXaxis()->SetTitle("Pt");
  h_kaon_pt_mediumMuId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_kaon_pt_mediumMuId_fake_rate", output_path, gPad);
  h_kaon_eta_mediumMuId_fake_rate.Divide(&h_kaon_eta_mediumMuId, &h_kaon_eta_all, 1, 1);
  h_kaon_eta_mediumMuId_fake_rate.SetMinimum(0);
  h_kaon_eta_mediumMuId_fake_rate.SetMarkerStyle(20);
  h_kaon_eta_mediumMuId_fake_rate.SetMarkerColor(kBlue);
  h_kaon_eta_mediumMuId_fake_rate.GetXaxis()->SetTitle("|#eta|");
  h_kaon_eta_mediumMuId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_kaon_eta_mediumMuId_fake_rate", output_path, gPad);

  h_kaon_pt_softMvaId_fake_rate.Divide(&h_kaon_pt_softMvaId, &h_kaon_pt_all, 1, 1);
  h_kaon_pt_softMvaId_fake_rate.SetMinimum(0);
  h_kaon_pt_softMvaId_fake_rate.SetMaximum(0.0032);
  h_kaon_pt_softMvaId_fake_rate.SetMarkerStyle(20);
  h_kaon_pt_softMvaId_fake_rate.SetMarkerColor(kRed);
  h_kaon_pt_softMvaId_fake_rate.GetXaxis()->SetTitle("Pt");
  h_kaon_pt_softMvaId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_kaon_pt_softMvaId_fake_rate", output_path, gPad);
  h_kaon_eta_softMvaId_fake_rate.Divide(&h_kaon_eta_softMvaId, &h_kaon_eta_all, 1, 1);
  h_kaon_eta_softMvaId_fake_rate.SetMinimum(0);
  h_kaon_eta_softMvaId_fake_rate.SetMarkerStyle(20);
  h_kaon_eta_softMvaId_fake_rate.SetMarkerColor(kRed);
  h_kaon_eta_softMvaId_fake_rate.GetXaxis()->SetTitle("|#eta|");
  h_kaon_eta_softMvaId_fake_rate.Draw("e0");
  print_canvas(sample.name + "_kaon_eta_softMvaId_fake_rate", output_path, gPad);

}

void fake_muon_study(){
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  for (auto& sample: samples){
    process_sample(sample);
  }
}

