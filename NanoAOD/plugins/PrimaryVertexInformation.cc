#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>
#include <algorithm>

// 
// PrimaryVertexInformation is designed for Bs/d->mumu analysis
//

using namespace std;
typedef reco::Candidate::LorentzVector LorentzVector;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

// offlineSlimmedPrimaryVertices
// floatedmValueMap_offlineSlimmedPrimaryVertices__PAT


class PrimaryVertexInformation : public edm::stream::EDProducer<> {
    
public:
    
  explicit PrimaryVertexInformation(const edm::ParameterSet &iConfig);
    
  ~PrimaryVertexInformation() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> vertexScoreToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  bool isMC_;
};

PrimaryVertexInformation::PrimaryVertexInformation(const edm::ParameterSet &iConfig):
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  vertexScoreToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("vertexScores"))),
  pfCandToken_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
  isMC_( iConfig.getParameter<bool>( "isMC" ) )
{
    produces<pat::CompositeCandidateCollection>("pvs");
}

void PrimaryVertexInformation::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  // Get primary vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);
  unsigned int nvtx = vertices->size();

  // Get vertex scores
  edm::Handle<edm::ValueMap<float>> vertexScores;
  iEvent.getByToken(vertexScoreToken_, vertexScores);

  // Get collection containing tracks
  edm::Handle<edm::View<pat::PackedCandidate>> pfCands;
  iEvent.getByToken(pfCandToken_, pfCands);

  // Output collection
  auto pvs = make_unique<pat::CompositeCandidateCollection>();

  std::vector<float> sumpt(nvtx, 0.0);
  std::vector<float> sumpt2(nvtx, 0.0);
  std::vector<int> ntrks(nvtx, 0);

  // Collect track information
  for (const auto& pfCand: *pfCands.product()){
    // keep only the tracks used in the PV fit
    if (pfCand.pvAssociationQuality() != pat::PackedCandidate::UsedInFitTight) continue;
    unsigned int index = pfCand.vertexRef().key();
    if (index < nvtx){
      sumpt[index] += pfCand.pt();
      sumpt2[index] += pfCand.pt() * pfCand.pt();
      ntrks[index] += 1;
    }
  }

  // Get PV info
  for (unsigned int i = 0; i < nvtx; i++) {
    const auto& vertex = vertices->at(i);
    reco::VertexRef vertexRef(vertices, i);
    float score = (*vertexScores)[vertexRef];

    pat::CompositeCandidate cand;

    cand.addUserFloat("score", score);
    cand.addUserFloat("ndof",  vertex.ndof());
    cand.addUserFloat("chi2",  vertex.chi2());
    cand.addUserFloat("x",     vertex.x());
    cand.addUserFloat("xErr",  vertex.xError());
    cand.addUserFloat("y",     vertex.y());
    cand.addUserFloat("yErr",  vertex.yError());
    cand.addUserFloat("z",     vertex.z());
    cand.addUserFloat("zErr",  vertex.zError());

    cand.addUserFloat("sumpt", sumpt[i]);
    cand.addUserFloat("sumpt2",sumpt2[i]);
    cand.addUserInt(  "ntrks", ntrks[i]);

      // if (isMC_){
      // 	mu_cand.addUserInt("simType", muon.simType());
      // 	mu_cand.addUserInt("simExtType", muon.simExtType());
      // 	mu_cand.addUserInt("simPdgId", muon.simPdgId());
      // 	mu_cand.addUserInt("simMotherPdgId", muon.simMotherPdgId());
      // 	mu_cand.addUserFloat("simProdRho", muon.simProdRho());
      // 	mu_cand.addUserFloat("simProdZ", muon.simProdZ());
      // }
      
      pvs->push_back(cand);
    }
    
    iEvent.put(move(pvs), "pvs");
}

DEFINE_FWK_MODULE(PrimaryVertexInformation);

//  LocalWords:  vertices
