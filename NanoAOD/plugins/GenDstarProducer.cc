#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>

// 
// GenDstarProducer gen info extraction for Bmm5 analysis

typedef reco::Candidate::LorentzVector LorentzVector;
namespace {
  const float muon_mass_    = 0.10565837;
}

class GenDstarProducer : public edm::stream::EDProducer<> {
public:
explicit GenDstarProducer(const edm::ParameterSet &iConfig);
~GenDstarProducer() override {};
    
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  int match_muon(const reco::Candidate* gen_muon);
  const pat::PackedCandidate* match_tack(const reco::Candidate* gen_cand);
  bool isGoodMuon(const pat::Muon& muon);
  bool isGoodTrack(const pat::PackedCandidate& trk);
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::GenParticle> >   prunedGenToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;

  edm::Handle<std::vector<pat::Muon>> muonHandle_;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_;
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

  void fill_final_state_daugheters(std::vector<const reco::GenParticle*>& daus,
				   const reco::GenParticle& cand)
  {
    for (auto const& dau: cand.daughterRefVector()){
      if (dau->daughterRefVector().empty())
	daus.push_back(&*dau);
      else
	fill_final_state_daugheters(daus, *dau);
    }
  }

  bool dr_match(const LorentzVector& reco , const LorentzVector& gen){
    if (fabs(reco.pt()-gen.pt())/gen.pt()<0.1 and deltaR(reco,gen)<0.02)
      return true;
    return false;
  }

}

bool GenDstarProducer::isGoodMuon(const pat::Muon& muon){
  if ( not muon.isLooseMuon() ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  return true;
}

bool GenDstarProducer::isGoodTrack(const pat::PackedCandidate& trk){
  if ( not trk.bestTrack()->quality(reco::Track::highPurity)) return false;
  return true;
}

GenDstarProducer::GenDstarProducer(const edm::ParameterSet &iConfig):
  prunedGenToken_( consumes<std::vector<reco::GenParticle>> ( edm::InputTag( "prunedGenParticles" ) ) ),
  packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( edm::InputTag( "packedGenParticles" ) ) ),
  muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  pfCandToken_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) )
{
  produces<pat::CompositeCandidateCollection>("gendstar");
}

int GenDstarProducer::match_muon(const reco::Candidate* gen_muon){
  int index = -1;
  for (unsigned int i=0; i < muonHandle_->size(); ++i){
    if (dr_match(muonHandle_->at(i).p4(), gen_muon->p4())){
      return i;
    }
  }
  return index;
}

const pat::PackedCandidate* GenDstarProducer::match_tack(const reco::Candidate* gen_cand){
  if (not gen_cand) return nullptr;
  for (auto const& pfCand: *pfCandHandle_){
    if (pfCand.charge() == 0 ) continue;
    if (not pfCand.hasTrackDetails()) continue;
    if (dr_match(pfCand.p4(), gen_cand->p4())){
      return &pfCand;
    }
  }
  return nullptr;
}

void GenDstarProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Output collection
  auto gendstar = std::make_unique<pat::CompositeCandidateCollection>();
    
  edm::Handle<std::vector<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  edm::Handle<std::vector<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);

  iEvent.getByToken(muonToken_, muonHandle_);
  iEvent.getByToken(pfCandToken_, pfCandHandle_);
  
  // Pruned GenParticles - subset of original GenParticles keeping
  // only "interesting" objects, i.e. leptons, neutrinos, B-hadrons,
  // vector bosons etc
  //
  // Packed particles are all the status 1 with reduced content.
  //
  // The navigation from packed to pruned is possible (the other
  // direction should be done manually)

  // std::vector<const reco::Candidate*> c_hadrons;
  
  for (auto const & cand: *pruned){
    // keep only interesting c-hadrons
    if ( abs(cand.pdgId()) != 413 )  // D*
      continue;

    // check direct daughter infor for signs of oscilations
    bool final_c = true;
    for (auto const& daughter: cand.daughterRefVector()){
      if (daughter->pdgId() == -cand.pdgId()){
	final_c = false;
	break;
      }
    }
    if (not final_c) continue;

    // find final state daughters and compute signature
    // note: ignore photos
    long long int signature = 1;
    std::vector<const reco::GenParticle*> daughters;
    fill_final_state_daugheters(daughters, cand);
    
    std::vector<const reco::GenParticle*> final_state_particles;
    LorentzVector radiation;

    for (auto const* dau: daughters){
      if (dau->pdgId()!=22){
	signature *= dau->pdgId();
	if (abs(dau->pdgId()) == 13  or // mu+/-
	    abs(dau->pdgId()) == 211 or // pi+/-
	    abs(dau->pdgId()) == 321 or // K+/-
	    abs(dau->pdgId()) == 14  or // nu_mu
	    abs(dau->pdgId()) == 2212 or // protons
	    dau->pdgId() == 111)        // pi0
	  final_state_particles.push_back(dau);
      } else {
	radiation += dau->p4();
      }
    }
    std::sort(final_state_particles.begin(), final_state_particles.end(), order_by_pt);

    // Loop over daughters
    
    // Interpret the signature
    // 
    // We will store directly acceessible the following information
    // - D pdg id
    // - D p3
    // - decay channel id based on the signature
    // 
    // The rest of information can be extracted from other daughters
    const reco::Candidate* dau1(nullptr);
    const reco::Candidate* dau2(nullptr);
    const reco::Candidate* dau3(nullptr);
    
    switch (signature){
    case  -13 *  13 * 211:  // mmpi
    case   13 *  13 * 211:  // mmpi
    case -211 * 211 * 211:  // pipipi
    case  211 * 211 * 211:  // pipipi
    case -321 * 321 * 211:  // KKpi
    case  321 * 321 * 211:  // KKpi
    case -321 * 211 * 211:  // Kpipi
    case  321 * 211 * 211:  // Kpipi
      dau1 = final_state_particles.at(0);
      dau2 = final_state_particles.at(1);
      dau3 = final_state_particles.at(2);
      break;

    default:       // unknown signature
      continue;
    }

    // fill information
    pat::CompositeCandidate cCand;
    cCand.addUserInt("signature",  signature);
    cCand.addUserInt("pdgId",  cand.pdgId());
    cCand.addUserInt("mpdgId", cand.mother()?cand.mother()->pdgId():0);
    cCand.addUserFloat("pt",   cand.pt());
    cCand.addUserFloat("eta",  cand.eta());
    cCand.addUserFloat("phi",  cand.phi());
    cCand.addUserFloat("mass", cand.mass());

    cCand.addUserInt(  "dau1_pdgId", dau1?dau1->pdgId():0);
    cCand.addUserInt(  "dau1_mpdgId",dau1&&dau1->mother()?dau1->mother()->pdgId():0);
    cCand.addUserFloat("dau1_pt",    dau1?dau1->pt():0);
    cCand.addUserFloat("dau1_eta",   dau1?dau1->eta():0);
    cCand.addUserFloat("dau1_phi",   dau1?dau1->phi():0);
    auto dau1_reco = match_tack(dau1);
    cCand.addUserFloat("dau1_reco_pt",    dau1_reco?dau1_reco->pt():0);
    cCand.addUserFloat("dau1_reco_eta",   dau1_reco?dau1_reco->eta():0);
    cCand.addUserFloat("dau1_reco_phi",   dau1_reco?dau1_reco->phi():0);
    cCand.addUserInt(  "dau1_reco_id",    dau1_reco?isGoodTrack(*dau1_reco):0);

    cCand.addUserInt(  "dau2_pdgId", dau2?dau2->pdgId():0);
    cCand.addUserInt(  "dau2_mpdgId",dau2&&dau2->mother()?dau2->mother()->pdgId():0);
    cCand.addUserFloat("dau2_pt",    dau2?dau2->pt():0);
    cCand.addUserFloat("dau2_eta",   dau2?dau2->eta():0);
    cCand.addUserFloat("dau2_phi",   dau2?dau2->phi():0);
    auto dau2_reco = match_tack(dau2);
    cCand.addUserFloat("dau2_reco_pt",    dau2_reco?dau2_reco->pt():0);
    cCand.addUserFloat("dau2_reco_eta",   dau2_reco?dau2_reco->eta():0);
    cCand.addUserFloat("dau2_reco_phi",   dau2_reco?dau2_reco->phi():0);
    cCand.addUserInt(  "dau2_reco_id",    dau2_reco?isGoodTrack(*dau2_reco):0);

    cCand.addUserInt(  "dau3_pdgId", dau3?dau3->pdgId():0);
    cCand.addUserInt(  "dau3_mpdgId",dau3&&dau3->mother()?dau3->mother()->pdgId():0);
    cCand.addUserFloat("dau3_pt",    dau3?dau3->pt():0);
    cCand.addUserFloat("dau3_eta",   dau3?dau3->eta():0);
    cCand.addUserFloat("dau3_phi",   dau3?dau3->phi():0);
    auto dau3_reco = match_tack(dau3);
    cCand.addUserFloat("dau3_reco_pt",    dau3_reco?dau3_reco->pt():0);
    cCand.addUserFloat("dau3_reco_eta",   dau3_reco?dau3_reco->eta():0);
    cCand.addUserFloat("dau3_reco_phi",   dau3_reco?dau3_reco->phi():0);
    cCand.addUserInt(  "dau3_reco_id",    dau3_reco?isGoodTrack(*dau3_reco):0);
      
    // radiation
    cCand.addUserFloat("rad_p",   radiation.P());    
    cCand.addUserFloat("rad_pt",  radiation.pt());
    cCand.addUserFloat("rad_eta", radiation.eta());
    cCand.addUserFloat("rad_phi", radiation.phi());
      
    gendstar->push_back(cCand);

  }
  iEvent.put(std::move(gendstar),"gendstar");
}
DEFINE_FWK_MODULE(GenDstarProducer);
