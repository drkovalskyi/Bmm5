// -*- C++ -*-
//
// Package:    Bmm5/BmmProdFilter
// Class:      BmmProdFilter
// 
/**\class BmmProdFilter BmmProdFilter.cc Bmm5/NanoAOD/plugins/BmmProdFilter.cc

 Description: Bmm5 analysis filter for skimming of AOD and MiniAOD data

 Implementation:
    - dimuon reco signature requirements
    - optional trigger requirements
*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Thu, 10 Oct 2019 09:27:23 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"


//
// class declaration
//

class BmmProdFilter : public edm::stream::EDFilter<> {
public:
  explicit BmmProdFilter(const edm::ParameterSet&);
  ~BmmProdFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool 
  isGoodMuon(const reco::Muon& muon);

  float 
  distanceOfClosestApproach( const reco::Track* track1,
			     const reco::Track* track2 );

  virtual bool filter(edm::Event&, const edm::EventSetup&) override;

  //virtual void endStream() override;
  //virtual void beginStream(edm::StreamID) override;
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;
  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  double maxTwoTrackDOCA_;
  double ptMinMu_;
  double etaMaxMu_;
};

bool 
BmmProdFilter::isGoodMuon(const reco::Muon& muon){
  if ( not muon::isLooseMuon(muon) ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < ptMinMu_ || fabs(muon.eta()) > etaMaxMu_ ) return false;
  return true;
}

float 
BmmProdFilter::distanceOfClosestApproach( const reco::Track* track1,
					     const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}


BmmProdFilter::BmmProdFilter(const edm::ParameterSet& iConfig):
  muonToken_( consumes<edm::View<reco::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
  ptMinMu_(         iConfig.getParameter<double>( "MuonMinPt" ) ),
  etaMaxMu_(        iConfig.getParameter<double>( "MuonMaxEta" ) )
{}


BmmProdFilter::~BmmProdFilter(){}

bool
BmmProdFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder_);
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  auto nMuons = muonHandle->size();
  if ( nMuons > 1 ){
    for (unsigned int i = 0; i < nMuons-1; ++i) {
      const reco::Muon & muon1 = muonHandle->at(i);
      if (not isGoodMuon(muon1)) continue;
      for (unsigned int j = i+1; j < nMuons; ++j) {
	const reco::Muon & muon2 = muonHandle->at(j);
	if (not isGoodMuon(muon2)) continue;
	auto mm_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						 muon2.innerTrack().get());
	if (maxTwoTrackDOCA_>0 and mm_doca > maxTwoTrackDOCA_) continue;
	return true;
      }
    }
  }
  return false;
}

/*
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
BmmProdFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
BmmProdFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
void
BmmProdFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
 
// ------------ method called when ending the processing of a run  ------------
void
BmmProdFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 
// ------------ method called when starting to processes a luminosity block  ------------
void
BmmProdFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
 
// ------------ method called when ending the processing of a luminosity block  ------------
void
BmmProdFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BmmProdFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(BmmProdFilter);
