#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>

// 
// GenBmmProducer gen info extraction for Bmm5 analysis

typedef reco::Candidate::LorentzVector LorentzVector;
namespace {
  const float muon_mass_    = 0.10565837;
}

class GenBmmProducer : public edm::EDProducer {
public:
explicit GenBmmProducer(const edm::ParameterSet &iConfig);
~GenBmmProducer() override {};
    
private:
virtual void produce(edm::Event&, const edm::EventSetup&);
// ----------member data ---------------------------
edm::EDGetTokenT<std::vector<reco::GenParticle> >   prunedGenToken_;
edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
};

namespace {
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* candidate){
    if (ancestor == candidate) return true;
    for (unsigned int iMother=0; iMother < candidate->numberOfMothers(); ++iMother){
      if (isAncestor(ancestor,candidate->mother(iMother)))
	return true;
    }
    return false;
  }
  bool order_by_pt(const reco::Candidate* a, const reco::Candidate* b)
  {
    return a->pt() > b->pt();
  }

}

GenBmmProducer::GenBmmProducer(const edm::ParameterSet &iConfig):
  prunedGenToken_( consumes<std::vector<reco::GenParticle>> ( edm::InputTag( "prunedGenParticles" ) ) ),
  packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( edm::InputTag( "packedGenParticles" ) ) )
{
  produces<pat::CompositeCandidateCollection>("genbmm");
  produces<pat::CompositeCandidateCollection>("gensummary");
}

void GenBmmProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Output collection
  auto genbmm = std::make_unique<pat::CompositeCandidateCollection>();
  auto genbinfo = std::make_unique<pat::CompositeCandidateCollection>();
    
  edm::Handle<std::vector<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  edm::Handle<std::vector<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);

  // Pruned GenParticles - subset of original GenParticles keeping
  // only "interesting" objects, i.e. leptons, neutrinos, B-hadrons,
  // vector bosons etc
  //
  // Packed particles are all the status 1 with reduced content.
  //
  // The navigation from packed to pruned is possible (the other
  // direction should be done manually)

  std::vector<const reco::Candidate*> b_hadrons;
  
  unsigned int n_bs(0);
  unsigned int n_anti_bs(0);
  for (auto const & cand: *pruned){
    if (cand.status()>=21 and cand.status()<=29){
      if (cand.pdgId()==5) n_bs++;
      if (cand.pdgId()==-5) n_anti_bs++;
    }

    // keep only interesting b-hadrons
    if ( abs(cand.pdgId()) != 511 and  // B0
	 abs(cand.pdgId()) != 521 and  // B+/-
	 abs(cand.pdgId()) != 531 and  // Bs 
	 abs(cand.pdgId()) != 541 )  // Bc 
      continue;

    // check direct daughter infor for signs of neutral B oscilations
    bool final_b = true;
    for (auto const& daughter: cand.daughterRefVector()){
      if (daughter->pdgId() == -cand.pdgId()){
	final_b = false;
	break;
      }
    }
    if (not final_b) continue;

    // find final state daughters and compute signature
    // note: ignore photos
    long long int signature = 1;
    std::vector<const reco::Candidate*> final_state_particles;
    LorentzVector radiation;

    // Loop over final state particles stored in the packed gen collection
    for (auto const& dau: *packed){
      auto mother = dau.mother(0);
      if (mother and isAncestor(&cand,mother)){
	if (dau.pdgId()!=22){
	  signature *= dau.pdgId();
	  if (abs(dau.pdgId()) == 13  or // mu+/-
	      abs(dau.pdgId()) == 211 or // pi+/-
	      abs(dau.pdgId()) == 321 or // K+/-
	      abs(dau.pdgId()) == 14  or // nu_mu
	      dau.pdgId() == 111)        // pi0
	    final_state_particles.push_back(&dau);
	} else {
	  radiation += dau.p4();
	}
      }
    }
    std::sort(final_state_particles.begin(), final_state_particles.end(), order_by_pt);

    // Interpret the signature
    // 
    // We will store directly acceessible the following information
    // - B pdg id
    // - B p3
    // - "dimuon" pair - most likely combination of daughters that may fake dimuon signature
    // - decay channel id based on the signature
    // 
    // The rest of information can be extracted from other daughters
    const reco::Candidate* mu_cand1(nullptr);
    const reco::Candidate* mu_cand2(nullptr);
    const reco::Candidate* dau3(nullptr);
    const reco::Candidate* dau4(nullptr);
    
    // skip irrelevant signatures
    // if ( (signature != 13*13*321*321)      // kkmm
    // 	 and (abs(signature) != 13*13*321) // kmm
    // 	 and (signature != -13*13) )        // mm
    //   continue;

    switch (signature){
    case -13*13:    // mm
    case -211*211:  // pi+pi-
    case -321*321:  // K+K-
    case -321*211:  // Kpi
    case -321*2212: // Kp
    case -211*2212: // pi p
      mu_cand1 = final_state_particles.at(0);
      mu_cand2 = final_state_particles.at(1);
      break;

    case 211*13*14:
    case -211*13*14: // pi mu nu
    case 321*13*14:
    case -321*13*14: // K mu nu
    case 2212*13*14:
    case -2212*13*14: // p mu nu
      if ( abs(final_state_particles.at(0)->pdgId()) == 14 ){
	dau3     = final_state_particles.at(0);
	mu_cand1 = final_state_particles.at(1);
	mu_cand2 = final_state_particles.at(2);
      }
      if ( abs(final_state_particles.at(1)->pdgId()) == 14 ){
	dau3     = final_state_particles.at(1);
	mu_cand1 = final_state_particles.at(0);
	mu_cand2 = final_state_particles.at(2);
      }
      if ( abs(final_state_particles.at(2)->pdgId()) == 14 ){
	dau3     = final_state_particles.at(2);
	mu_cand1 = final_state_particles.at(0);
	mu_cand2 = final_state_particles.at(1);
      }
      break;

    case -111*13*13:  // pi0 mu mu
      if ( final_state_particles.at(0)->pdgId() == 111 ){
	dau3     = final_state_particles.at(0);
	mu_cand1 = final_state_particles.at(1);
	mu_cand2 = final_state_particles.at(2);
      }
      if ( final_state_particles.at(1)->pdgId() == 111 ){
	dau3     = final_state_particles.at(1);
	mu_cand1 = final_state_particles.at(0);
	mu_cand2 = final_state_particles.at(2);
      }
      if ( final_state_particles.at(2)->pdgId() == 111 ){
	dau3     = final_state_particles.at(2);
	mu_cand1 = final_state_particles.at(0);
	mu_cand2 = final_state_particles.at(1);
      }
      break;

    case 321*321*13*13:  // KKmm
      for (auto dau: final_state_particles){
	if (abs(dau->pdgId())==13){
	  if (not mu_cand1)
	    mu_cand1 = dau;
	  else
	    mu_cand2 = dau;
	} else {
	  if (not dau3)
	    dau3 = dau;
	  else
	    dau4 = dau;
	}
      }
      break;

    default:       // unknown signature
      continue;
    }

    // fill information
    pat::CompositeCandidate bCand;
    bCand.addUserInt("signature",  signature);
    bCand.addUserInt("pdgId",  cand.pdgId());
    bCand.addUserFloat("pt",   cand.pt());
    bCand.addUserFloat("eta",  cand.eta());
    bCand.addUserFloat("phi",  cand.phi());
    bCand.addUserFloat("mass", cand.mass());

    // dimuon info
    bCand.addUserInt(  "mu1_pdgId", mu_cand1->pdgId());
    bCand.addUserFloat("mu1_pt",    mu_cand1->pt());
    bCand.addUserFloat("mu1_eta",   mu_cand1->eta());
    bCand.addUserFloat("mu1_phi",   mu_cand1->phi());
    bCand.addUserInt(  "mu2_pdgId", mu_cand2->pdgId());
    bCand.addUserFloat("mu2_pt",    mu_cand2->pt());
    bCand.addUserFloat("mu2_eta",   mu_cand2->eta());
    bCand.addUserFloat("mu2_phi",   mu_cand2->phi());
    auto dimuon_p2 = 
      pow(mu_cand1->px() + mu_cand2->px(),2) +
      pow(mu_cand1->py() + mu_cand2->py(),2) +
      pow(mu_cand1->pz() + mu_cand2->pz(),2);
    auto dimuon_e = 
      sqrt(pow(mu_cand1->p(),2) + pow(muon_mass_,2)) +
      sqrt(pow(mu_cand2->p(),2) + pow(muon_mass_,2));
    
    bCand.addUserFloat("dimuon_mass", sqrt(pow(dimuon_e,2) - dimuon_p2));
    
    // other daughters
    bCand.addUserInt(  "dau3_pdgId", dau3?dau3->pdgId():0);
    bCand.addUserFloat("dau3_pt",    dau3?dau3->pt():0);
    bCand.addUserFloat("dau3_eta",   dau3?dau3->eta():0);
    bCand.addUserFloat("dau3_phi",   dau3?dau3->phi():0);
    bCand.addUserInt(  "dau4_pdgId", dau4?dau4->pdgId():0);
    bCand.addUserFloat("dau4_pt",    dau4?dau4->pt():0);
    bCand.addUserFloat("dau4_eta",   dau4?dau4->eta():0);
    bCand.addUserFloat("dau4_phi",   dau4?dau4->phi():0);
    // bCand.addUserFloat("kk_mass",   kaons.size()>1?(kaons[0]->p4()+kaons[1]->p4()).mass():0);

    // radiation
    bCand.addUserFloat("rad_p",   radiation.P());    
    bCand.addUserFloat("rad_pt",  radiation.pt());
    bCand.addUserFloat("rad_eta", radiation.eta());
    bCand.addUserFloat("rad_phi", radiation.phi());
      
    genbmm->push_back(bCand);

  }
  // fill information
  pat::CompositeCandidate bgen;
  bgen.addUserInt("n_b",  n_bs);
  bgen.addUserInt("n_anti_b",  n_anti_bs);
  
  unsigned int processType(0);
  if (n_bs==0 and n_anti_bs==0){
    processType = 42; // gluon splitting
  }else{
    if (n_bs==0 or n_anti_bs==0){
      processType = 41; // flavor excitation
    }else{
      processType = 40; // gluon fusion
    }
  }
  bgen.addUserInt("process_type",processType);
  genbinfo->push_back(bgen);

  iEvent.put(std::move(genbmm),"genbmm");
  iEvent.put(std::move(genbinfo),"gensummary");
}
DEFINE_FWK_MODULE(GenBmmProducer);
