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
  auto b2kkmm = std::make_unique<pat::CompositeCandidateCollection>();
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
    if( abs(cand.pdgId()) != 521     // B+/-
	and abs(cand.pdgId()) != 511 // B
	and abs(cand.pdgId()) != 531 )// Bs 
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
    std::vector<const reco::Candidate*> muons;
    std::vector<const reco::Candidate*> kaons;
    LorentzVector radiation;

    for (auto const& cand: *packed){
      auto mother = cand.mother(0);
      if (mother and isAncestor(&cand,mother)){
	if (cand.pdgId()!=22){
	  signature *= cand.pdgId();
	  if (abs(cand.pdgId()) == 13)  muons.push_back(&cand);
	  if (abs(cand.pdgId()) == 321) kaons.push_back(&cand);
	} else {
	  radiation += cand.p4();
	}
      }
    }
    
    // skip irrelevant signatures
    if ( (signature != 13*13*321*321)      // kkmm
	 and (abs(signature) != 13*13*321) // kmm
	 and (signature != 13*13) )        // mm
      continue;

    // fill information
    pat::CompositeCandidate bCand;
    bCand.addUserInt("pdgId",  cand.pdgId());
    bCand.addUserFloat("pt",   cand.pt());
    bCand.addUserFloat("eta",  cand.eta());
    bCand.addUserFloat("phi",  cand.phi());
    bCand.addUserFloat("mass", cand.mass());

    // dimuon info
    std::sort(muons.begin(), muons.end(), order_by_pt);
    bCand.addUserInt(  "mu1_pdgId", muons[0]->pdgId());
    bCand.addUserFloat("mu1_pt",    muons[0]->pt());
    bCand.addUserFloat("mu1_eta",   muons[0]->eta());
    bCand.addUserFloat("mu1_phi",   muons[0]->phi());
    bCand.addUserInt(  "mu2_pdgId", muons[1]->pdgId());
    bCand.addUserFloat("mu2_pt",    muons[1]->pt());
    bCand.addUserFloat("mu2_eta",   muons[1]->eta());
    bCand.addUserFloat("mu2_phi",   muons[1]->phi());
    bCand.addUserFloat("dimuon_mass", (muons[0]->p4()+muons[1]->p4()).mass());
    
    // kaon(s) info
    std::sort(kaons.begin(), kaons.end(), order_by_pt);
    bCand.addUserInt(  "kaon1_pdgId", kaons.size()>0?kaons[0]->pdgId():0);
    bCand.addUserFloat("kaon1_pt",    kaons.size()>0?kaons[0]->pt():0);
    bCand.addUserFloat("kaon1_eta",   kaons.size()>0?kaons[0]->eta():0);
    bCand.addUserFloat("kaon1_phi",   kaons.size()>0?kaons[0]->phi():0);
    bCand.addUserInt(  "kaon2_pdgId", kaons.size()>1?kaons[1]->pdgId():0);
    bCand.addUserFloat("kaon2_pt",    kaons.size()>1?kaons[1]->pt():0);
    bCand.addUserFloat("kaon2_eta",   kaons.size()>1?kaons[1]->eta():0);
    bCand.addUserFloat("kaon2_phi",   kaons.size()>1?kaons[1]->phi():0);
    bCand.addUserFloat("kk_mass",   kaons.size()>1?(kaons[0]->p4()+kaons[1]->p4()).mass():0);

    // radiation
    bCand.addUserFloat("rad_p",   radiation.P());    
    bCand.addUserFloat("rad_pt",  radiation.pt());
    bCand.addUserFloat("rad_eta", radiation.eta());
    bCand.addUserFloat("rad_phi", radiation.phi());
      
    b2kkmm->push_back(bCand);

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

  iEvent.put(std::move(b2kkmm),"genbmm");
  iEvent.put(std::move(genbinfo),"gensummary");
}
DEFINE_FWK_MODULE(GenBmmProducer);
