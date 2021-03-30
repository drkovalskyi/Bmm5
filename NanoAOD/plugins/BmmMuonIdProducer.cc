#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/PointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TMVA/Reader.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>
#include <algorithm>

// 
// BmmMuonIdProducer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;

using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class BmmMuonIdProducer : public edm::EDProducer {
    
public:
    
  explicit BmmMuonIdProducer(const edm::ParameterSet &iConfig);
    
  ~BmmMuonIdProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);
  // GenMatchInfo getGenMatchInfo( const pat::PackedCandidate& track1,
  // const pat::PackedCandidate& track2);

  // ----------member data ---------------------------
    
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  bool isMC_;
};

BmmMuonIdProducer::BmmMuonIdProducer(const edm::ParameterSet &iConfig):
muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( iConfig.getParameter<edm::InputTag>( "packedGenParticleCollection" ) ) ),
packedGenParticles_(nullptr),
isMC_(             iConfig.getParameter<bool>( "isMC" ) )
{
    produces<pat::CompositeCandidateCollection>("muons");
}

void BmmMuonIdProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<std::vector<pat::Muon> > muonHandle;
    iEvent.getByToken(muonToken_, muonHandle);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      packedGenParticles_ = nullptr;
    }

    // Output collection
    auto muons = std::make_unique<pat::CompositeCandidateCollection>();

    for ( const auto& muon: *muonHandle.product()){
      pat::CompositeCandidate mu_cand;
      mu_cand.addUserFloat("trkKink", muon.combinedQuality().trkKink);
      mu_cand.addUserFloat("glbTrackProbability", muon.combinedQuality().glbTrackProbability);
      muons->push_back(mu_cand);
    }
    
    iEvent.put(std::move(muons), "muons");
}

DEFINE_FWK_MODULE(BmmMuonIdProducer);

//  LocalWords:  vertices
