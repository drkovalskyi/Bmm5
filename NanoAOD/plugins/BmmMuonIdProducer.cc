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

#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>
#include <algorithm>

// 
// BmmMuonIdProducer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;
typedef std::pair<const reco::MuonChamberMatch*, const reco::MuonSegmentMatch*> MatchPair;
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
  void fillMatchInfo(pat::CompositeCandidate& cand, const pat::Muon& muon);
  
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
      mu_cand.addUserFloat("chi2LocalPosition", muon.combinedQuality().chi2LocalPosition);
      if (muon.isGlobalMuon())
	mu_cand.addUserFloat("glbNormChi2", muon.globalTrack()->normalizedChi2());
      else
	mu_cand.addUserFloat("glbNormChi2", 9999.);

      if (muon.isTrackerMuon() or muon.isGlobalMuon()){
	mu_cand.addUserFloat("trkValidFrac",  muon.innerTrack()->validFraction());
	
	mu_cand.addUserInt("nPixels",         muon.innerTrack()->hitPattern().numberOfValidPixelHits());
	mu_cand.addUserInt("nValidHits",      muon.innerTrack()->hitPattern().numberOfValidTrackerHits());
	mu_cand.addUserInt("nLostHitsInner",  muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS));
	mu_cand.addUserInt("nLostHitsOn",     muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::TRACK_HITS));
	mu_cand.addUserInt("nLostHitsOuter",  muon.innerTrack()->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS));
	
	mu_cand.addUserInt("trkLayers",           muon.innerTrack()->hitPattern().trackerLayersWithMeasurement());
	mu_cand.addUserInt("trkLostLayersInner",  muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS));
	mu_cand.addUserInt("trkLostLayersOn",     muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));
	mu_cand.addUserInt("trkLostLayersOuter",  muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_OUTER_HITS));

	mu_cand.addUserInt("highPurity",   muon.innerTrack()->quality(reco::Track::highPurity));

      } else {
	mu_cand.addUserFloat("trkValidFrac",  0);
	
	mu_cand.addUserInt("nPixels",         0);
	mu_cand.addUserInt("nValidHits",      0);
	mu_cand.addUserInt("nLostHitsInner",  0);
	mu_cand.addUserInt("nLostHitsOn",     0);
	mu_cand.addUserInt("nLostHitsOutter", 0);
	
	mu_cand.addUserInt("trkLayers",           0);
	mu_cand.addUserInt("trkLostLayersInner",  0);
	mu_cand.addUserInt("trkLostLayersOn",     0);
	mu_cand.addUserInt("trkLostLayersOuter",  0);

	mu_cand.addUserInt("highPurity",   0);

      }
	
      fillMatchInfo(mu_cand, muon);
      
      if (isMC_){
	mu_cand.addUserInt("simType", muon.simType());
	mu_cand.addUserInt("simExtType", muon.simExtType());
	mu_cand.addUserInt("simPdgId", muon.simPdgId());
	mu_cand.addUserInt("simMotherPdgId", muon.simMotherPdgId());
	mu_cand.addUserFloat("simProdRho", muon.simProdRho());
	mu_cand.addUserFloat("simProdZ", muon.simProdZ());
      }
      muons->push_back(mu_cand);
    }
    
    iEvent.put(std::move(muons), "muons");
}


const MatchPair&
getBetterMatch(const MatchPair& match1, const MatchPair& match2){

  // Prefer DT over CSC simply because it's closer to IP
  // and will have less multiple scattering (at least for
  // RB1 vs ME1/3 case). RB1 & ME1/2 overlap is tiny
  if (match2.first->detector() == MuonSubdetId::DT and
      match1.first->detector() != MuonSubdetId::DT)
    return match2;

  // For the rest compare local x match. We expect that
  // segments belong to the muon, so the difference in
  // local x is a reflection on how well we can measure it
  if ( abs(match1.first->x - match1.second->x) >
       abs(match2.first->x - match2.second->x) )
    return match2;
    
  return match1;
}

float dX(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return (match.first->x - match.second->x);
  else
    return 9999.;
}

float pullX(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return dX(match) /
      sqrt(std::pow(match.first->xErr, 2) + std::pow(match.second->xErr, 2));
  else
    return 9999.;
}

float pullDxDz(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return (match.first->dXdZ - match.second->dXdZ) /
           sqrt(std::pow(match.first->dXdZErr, 2) + std::pow(match.second->dXdZErr, 2));
  else
    return 9999.;
}

float dY(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return (match.first->y - match.second->y);
  else
    return 9999.;
}

float pullY(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return dY(match) /
      sqrt(std::pow(match.first->yErr, 2) + std::pow(match.second->yErr, 2));
  else
    return 9999.;
}

float pullDyDz(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return (match.first->dYdZ - match.second->dYdZ) /
           sqrt(std::pow(match.first->dYdZErr, 2) + std::pow(match.second->dYdZErr, 2));
  else
    return 9999.;
}

void fillMatchInfoForStation(string prefix,
			     pat::CompositeCandidate& cand,
			     const MatchPair& match){
  cand.addUserFloat(prefix + "_dX",       dX(match));
  cand.addUserFloat(prefix + "_pullX",    pullX(match));
  cand.addUserFloat(prefix + "_pullDxDz", pullDxDz(match));
  cand.addUserFloat(prefix + "_dY",       dY(match));
  cand.addUserFloat(prefix + "_pullY",    pullY(match));
  cand.addUserFloat(prefix + "_pullDyDz", pullDyDz(match));
}

void BmmMuonIdProducer::fillMatchInfo(pat::CompositeCandidate& cand,
				      const pat::Muon& muon){
  // Initiate containter for results
  const int n_stations = 2;
  std::vector<MatchPair> matches;
  for (unsigned int i=0; i < n_stations; ++i)
    matches.push_back(std::pair(nullptr, nullptr));

  // Find best matches
  for (auto& chamberMatch : muon.matches()){
    unsigned int station = chamberMatch.station() - 1;
    if (station >= n_stations) continue;

    // Find best segment match.
    // We could consider all segments, but we will restrict to segments
    // that match to this candidate better than to other muon candidates
    for (auto& segmentMatch : chamberMatch.segmentMatches){
      if ( not segmentMatch.isMask(reco::MuonSegmentMatch::BestInStationByDR) ||
	   not segmentMatch.isMask(reco::MuonSegmentMatch::BelongsToTrackByDR) )
	continue;

      // Multiple segment matches are possible in different
      // chambers that are either overlapping or belong to
      // different detectors. We need to select one.
      auto match_pair = MatchPair(&chamberMatch, &segmentMatch);
      
      if (matches[station].first)
	matches[station] = getBetterMatch(matches[station], match_pair);
      else
	matches[station] = match_pair;
    }
  }

  // Fill matching information
  fillMatchInfoForStation("match1", cand, matches[0]);
  fillMatchInfoForStation("match2", cand, matches[1]);
}


DEFINE_FWK_MODULE(BmmMuonIdProducer);

//  LocalWords:  vertices
