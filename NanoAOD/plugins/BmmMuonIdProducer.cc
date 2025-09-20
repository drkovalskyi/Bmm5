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
#include <algorithm>

#include "Bmm5/NanoAOD/interface/XGBooster.h"
#include "Bmm5/NanoAOD/interface/CommonTools.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

// 
// BmmMuonIdProducer is designed for Bs/d->mumu analysis
//

using namespace std;
typedef reco::Candidate::LorentzVector LorentzVector;
typedef pair<const reco::MuonChamberMatch*, const reco::MuonSegmentMatch*> MatchPair;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class BmmMuonIdProducer : public edm::stream::EDProducer<> {
    
public:
    
  explicit BmmMuonIdProducer(const edm::ParameterSet &iConfig);
    
  ~BmmMuonIdProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);
  // GenMatchInfo getGenMatchInfo( const pat::PackedCandidate& track1,
  // const pat::PackedCandidate& track2);
  void fillMatchInfo(pat::CompositeCandidate& cand, const pat::Muon& muon);
  void fillSoftMva(pat::CompositeCandidate& mu_cand);
  const l1t::Muon* getL1Muon(const reco::Candidate& cand, float max_dr=0.5);
  
  // ----------member data ---------------------------
    
  edm::EDGetTokenT<vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<vector<pat::PackedGenParticle> >   packedGenToken_;
  edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone> > triggerInfoToken_;
  edm::EDGetTokenT<BXVector<l1t::Muon> > l1Token_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  
  const vector<pat::PackedGenParticle>* packedGenParticles_;
  bool isMC_;
  string triggerCollection_;
  vector<string> triggers_;
  vector<string> features_;
  vector<string> xgboost_models_;
  vector<string> xgboost_variable_names_;
  vector<XGBooster> softMuonMva_;
  edm::Handle<BXVector<l1t::Muon> >   l1Handle_;
};

BmmMuonIdProducer::BmmMuonIdProducer(const edm::ParameterSet &iConfig):
muonToken_( consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
packedGenToken_( consumes<vector<pat::PackedGenParticle>> ( iConfig.getParameter<edm::InputTag>( "packedGenParticleCollection" ) ) ),
triggerInfoToken_( consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigger") ) ),
l1Token_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1Src"))),
vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
pfCandToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
packedGenParticles_(nullptr),
isMC_( iConfig.getParameter<bool>( "isMC" ) ),
triggerCollection_( iConfig.getParameter<string>( "triggerCollection" ) ),
triggers_( iConfig.getParameter<vector<string>>( "triggers" ) ),
xgboost_models_( iConfig.getParameter<vector<string>>( "xgboost_models" ) ),
xgboost_variable_names_( iConfig.getParameter<vector<string>>( "xgboost_variable_names" ) )
{
    produces<pat::CompositeCandidateCollection>("muons");
    for (auto model: xgboost_models_)
      softMuonMva_.push_back(XGBooster(edm::FileInPath("Bmm5/NanoAOD/data/muon_mva/" + model + ".model").fullPath(),
				       edm::FileInPath("Bmm5/NanoAOD/data/muon_mva/" + model + ".features").fullPath()));
}

void BmmMuonIdProducer::fillSoftMva(pat::CompositeCandidate& mu_cand){
  // "match2_pullDyDz"
  // "match1_pullDyDz"
  // "match2_pullDxDz"
  // "match1_pullDxDz"
  for (unsigned int i=0; i < softMuonMva_.size(); ++i) {
    softMuonMva_.at(i).set("pt",                  mu_cand.pt());
    softMuonMva_.at(i).set("eta",                 mu_cand.eta());
    softMuonMva_.at(i).set("trkValidFrac",        mu_cand.userFloat("trkValidFrac"));
    softMuonMva_.at(i).set("glbTrackProbability", mu_cand.userFloat("glbTrackProbability"));
    softMuonMva_.at(i).set("nLostHitsInner",      mu_cand.userInt(  "nLostHitsInner"));
    softMuonMva_.at(i).set("nLostHitsOuter",      mu_cand.userInt(  "nLostHitsOuter"));
    softMuonMva_.at(i).set("trkKink",             mu_cand.userFloat("trkKink"));
    softMuonMva_.at(i).set("chi2LocalPosition",   mu_cand.userFloat("chi2LocalPosition"));
    softMuonMva_.at(i).set("match2_dX",           mu_cand.userFloat("match2_dX"));
    softMuonMva_.at(i).set("match2_pullX",        mu_cand.userFloat("match2_pullX"));
    softMuonMva_.at(i).set("match1_dX",           mu_cand.userFloat("match1_dX"));
    softMuonMva_.at(i).set("match1_pullX",        mu_cand.userFloat("match1_pullX"));
    softMuonMva_.at(i).set("nPixels",             mu_cand.userInt(  "nPixels"));
    softMuonMva_.at(i).set("nValidHits",          mu_cand.userInt(  "nValidHits"));
    softMuonMva_.at(i).set("nLostHitsOn",         mu_cand.userInt(  "nLostHitsOn"));
    softMuonMva_.at(i).set("match2_dY",           mu_cand.userFloat("match2_dY"));
    softMuonMva_.at(i).set("match2_pullY",        mu_cand.userFloat("match2_pullY"));
    softMuonMva_.at(i).set("match1_dY",           mu_cand.userFloat("match1_dY"));
    softMuonMva_.at(i).set("match1_pullY",        mu_cand.userFloat("match1_pullY"));
    softMuonMva_.at(i).set("match2_pullDyDz",     mu_cand.userFloat("match2_pullDyDz"));
    softMuonMva_.at(i).set("match1_pullDyDz",     mu_cand.userFloat("match1_pullDyDz"));
    softMuonMva_.at(i).set("match2_pullDxDz",     mu_cand.userFloat("match2_pullDxDz"));
    softMuonMva_.at(i).set("match1_pullDxDz",     mu_cand.userFloat("match1_pullDxDz"));
    softMuonMva_.at(i).set("glbNormChi2",         mu_cand.userFloat("glbNormChi2"));
    softMuonMva_.at(i).set("trkLayers",           mu_cand.userInt("trkLayers"));
    softMuonMva_.at(i).set("highPurity",          mu_cand.userInt("highPurity"));
    
    mu_cand.addUserFloat("xgb_" + xgboost_variable_names_.at(i), softMuonMva_.at(i).predict());
  }
}


const l1t::Muon* BmmMuonIdProducer::getL1Muon( const reco::Candidate& cand, float max_dr){
  const l1t::Muon* match = nullptr;
  float best_dr = 999.;
  // Loop over L1 candidates from BX 0 only
  for (auto it = l1Handle_->begin(0); it != l1Handle_->end(0); it++){
    float dr = deltaR(it->etaAtVtx(), it->phiAtVtx(), cand.eta(), cand.phi());
    if (dr > max_dr) continue;
    if (match == nullptr or dr < best_dr){
      best_dr = dr;
      match = &*it;
    }
  }
  return match;
}


void BmmMuonIdProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<vector<pat::Muon> > muonHandle;
    iEvent.getByToken(muonToken_, muonHandle);
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(pfCandToken_, pfCandHandle);
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByToken(vertexToken_, pvHandle);

    edm::Handle<vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      packedGenParticles_ = nullptr;
    }

    edm::Handle<vector<pat::TriggerObjectStandAlone> > trigger_info;
    iEvent.getByToken(triggerInfoToken_, trigger_info);

    iEvent.getByToken(l1Token_, l1Handle_);

    // Output collection
    auto muons = make_unique<pat::CompositeCandidateCollection>();

    for ( const auto& muon: *muonHandle.product()){
      pat::CompositeCandidate mu_cand(reco::CompositeCandidate(muon.charge(), muon.p4()));
      mu_cand.addUserFloat("trkKink",             muon.combinedQuality().trkKink);
      mu_cand.addUserFloat("glbTrackProbability", muon.combinedQuality().glbTrackProbability);
      mu_cand.addUserFloat("chi2LocalPosition",   muon.combinedQuality().chi2LocalPosition);
      mu_cand.addUserFloat("chi2LocalMomentum",   muon.combinedQuality().chi2LocalMomentum);
      mu_cand.addUserFloat("trkRelChi2",          muon.combinedQuality().trkRelChi2);
      mu_cand.addUserFloat("staRelChi2",          muon.combinedQuality().staRelChi2);
      mu_cand.addUserFloat("mvaId",               muon.mvaIDValue());
      
      if (muon.isGlobalMuon()){
	mu_cand.addUserFloat("glbNormChi2", muon.globalTrack()->normalizedChi2());
	mu_cand.addUserInt("chargeProduct", muon.outerTrack()->charge() * muon.innerTrack()->charge());
      } else {
	mu_cand.addUserFloat("glbNormChi2", 9999.);
	mu_cand.addUserInt("chargeProduct", 0);
      }
      
      if (muon.isStandAloneMuon()) {
	mu_cand.addUserFloat("staNormChi2", muon.outerTrack()->normalizedChi2());
	mu_cand.addUserFloat("staNdof", muon.outerTrack()->ndof());
	mu_cand.addUserInt(  "staValidHits", muon.outerTrack()->hitPattern().muonStationsWithValidHits());
      } else {
	mu_cand.addUserFloat("staNormChi2", 9999.);
	mu_cand.addUserFloat("staNdof", 0);
	mu_cand.addUserInt(  "staValidHits", 0);
      }
      
      if (muon.isTrackerMuon() or muon.isGlobalMuon()){
	bmm::fill_track_info(mu_cand, muon.innerTrack().get());
      } else {
	bmm::fill_track_info(mu_cand, nullptr);
      }
	
      fillMatchInfo(mu_cand, muon);
      fillSoftMva(mu_cand);
      
      if (isMC_){
	mu_cand.addUserInt("simType", muon.simType());
	mu_cand.addUserInt("simExtType", muon.simExtType());
	mu_cand.addUserInt("simPdgId", muon.simPdgId());
	mu_cand.addUserInt("simMotherPdgId", muon.simMotherPdgId());
	mu_cand.addUserFloat("simProdRho", muon.simProdRho());
	mu_cand.addUserFloat("simProdZ", muon.simProdZ());
      }

      // PF info
      const pat::PackedCandidate* pfmuon = nullptr;
      if (muon.bestTrack()) {
	for (const auto &pfcand: *pfCandHandle.product()) {
	  if (abs(pfcand.pdgId()) != 13) continue;
	  if (deltaR(*muon.bestTrack(), pfcand) < 0.01) {
	    pfmuon = &pfcand;
	    break;
	  }
	}
      }

      if (pfmuon and pfmuon->vertexRef().key() < pvHandle->size()) {
	int pv_index = int(pfmuon->vertexRef().key());
	mu_cand.addUserInt("pv_index", pv_index);
	mu_cand.addUserFloat("pv_z", pvHandle->at(pv_index).z());
      } else {
	mu_cand.addUserInt("pv_index", -1);
	mu_cand.addUserFloat("pv_z", 9999.);
      }
      
      /////////////////// Trigger info

      ////// HLT
      const pat::TriggerObjectStandAlone* best_trigger_object(nullptr);
      double best_dr(9999.);
      for (const auto &trigger_object : *trigger_info) {
	if (not trigger_object.hasTriggerObjectType(trigger::TriggerMuon)) continue;
	  // Restict muon collection to L3 muons
	  if (trigger_object.collection().find(triggerCollection_ + ":") == string::npos) continue;

	  double dr = deltaR(trigger_object, muon);
	  if (dr < best_dr){
	    best_dr = dr;
	    best_trigger_object = &trigger_object;
	  }
      }
    
      if (best_dr < 0.02){
	mu_cand.addUserFloat("hlt_pt", best_trigger_object->pt());
	mu_cand.addUserFloat("hlt_dr", best_dr);
	for ( auto const &trigger: triggers_ )
	  mu_cand.addUserInt(trigger, muon.triggered((trigger + "_v*").c_str()));
      } else {
	mu_cand.addUserFloat("hlt_pt", -1);
	mu_cand.addUserFloat("hlt_dr", -1);
	for ( auto const &trigger: triggers_ )
	  mu_cand.addUserInt(trigger, 0);
      }

      ///// L1

      const l1t::Muon* l1_muon = getL1Muon(muon);
      if (l1_muon){
	float dr = deltaR(l1_muon->etaAtVtx(), l1_muon->phiAtVtx(), muon.eta(), muon.phi());
	mu_cand.addUserFloat("l1_pt", l1_muon->pt());
	mu_cand.addUserFloat("l1_eta", l1_muon->eta());
	mu_cand.addUserFloat("l1_phi", l1_muon->phi());
	mu_cand.addUserInt("l1_quality", l1_muon->hwQual());
	mu_cand.addUserFloat("l1_phiAtVtx", l1_muon->phiAtVtx());
	mu_cand.addUserFloat("l1_etaAtVtx", l1_muon->etaAtVtx());
	mu_cand.addUserFloat("l1_dr", dr);
      } else {
	// defaults
	mu_cand.addUserFloat("l1_pt", -1);
	mu_cand.addUserFloat("l1_eta", 0);
	mu_cand.addUserFloat("l1_phi", 0);
	mu_cand.addUserInt("l1_quality", 0);
	mu_cand.addUserFloat("l1_phiAtVtx", 0);
	mu_cand.addUserFloat("l1_etaAtVtx", 0);
	mu_cand.addUserFloat("l1_dr", -1);
      }
      
      muons->push_back(mu_cand);
    }
    
    iEvent.put(move(muons), "muons");
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
      sqrt(pow(match.first->xErr, 2) + pow(match.second->xErr, 2));
  else
    return 9999.;
}

float pullDxDz(const MatchPair& match){
  if (match.first and match.second->hasPhi())
    return (match.first->dXdZ - match.second->dXdZ) /
           sqrt(pow(match.first->dXdZErr, 2) + pow(match.second->dXdZErr, 2));
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
      sqrt(pow(match.first->yErr, 2) + pow(match.second->yErr, 2));
  else
    return 9999.;
}

float pullDyDz(const MatchPair& match){
  if (match.first and match.second->hasZed())
    return (match.first->dYdZ - match.second->dYdZ) /
           sqrt(pow(match.first->dYdZErr, 2) + pow(match.second->dYdZErr, 2));
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
  vector<MatchPair> matches;
  for (unsigned int i=0; i < n_stations; ++i)
    matches.push_back(pair(nullptr, nullptr));

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
