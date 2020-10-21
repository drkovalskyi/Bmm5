// -*- C++ -*-
//
// Package:    Bmm5/MuonFakeFilter
// Class:      MuonFakeFilter
// 
/**\class MuonFakeFilter MuonFakeFilter.cc Bmm5/NanoAOD/plugins/MuonFakeFilter.cc

 Description: MiniAOD filter to select interesting fake muon events.

 Implementation:
    - find clean muons matching pions, kaon and protons by dR and dPt
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"


//
// class declaration
//

class MuonFakeFilter : public edm::stream::EDFilter<> {
public:
  explicit MuonFakeFilter(const edm::ParameterSet&);
  ~MuonFakeFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool 
  isGoodMuon(const pat::Muon& muon);

  virtual bool filter(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenToken_;
  double maxDR_, maxDPt_, ptMinMu_, etaMaxMu_;
};

bool 
MuonFakeFilter::isGoodMuon(const pat::Muon& muon){
  if ( not muon::isLooseMuon(muon) ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.isGlobalMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < ptMinMu_ || fabs(muon.eta()) > etaMaxMu_ ) return false;
  return true;
}

MuonFakeFilter::MuonFakeFilter(const edm::ParameterSet& iConfig):
  muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  prunedGenToken_( consumes<std::vector<reco::GenParticle>> ( iConfig.getParameter<edm::InputTag>( "prunedGenParticleCollection" ) ) ),
  maxDR_(           iConfig.getParameter<double>( "maxDR" ) ),
  maxDPt_(          iConfig.getParameter<double>( "maxDPt" ) ),
  ptMinMu_(         iConfig.getParameter<double>( "minPt" ) ),
  etaMaxMu_(        iConfig.getParameter<double>( "maxEta" ) )
{}

MuonFakeFilter::~MuonFakeFilter(){}

bool
MuonFakeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  bool hasInterstingMuon = false;
  for (const auto& muon: *muonHandle){
    if (not isGoodMuon(muon)) continue;
    hasInterstingMuon = true;
    break;
  }
  if (not hasInterstingMuon) return false;

  edm::Handle<std::vector<reco::GenParticle> > prunedGenParticleHandle;
  iEvent.getByToken(prunedGenToken_,prunedGenParticleHandle);
  auto prunedGenParticles = prunedGenParticleHandle.product();
  
  for (const auto& muon: *muonHandle){
    if (not isGoodMuon(muon)) continue;
    for (auto const & genParticle: *prunedGenParticles){
      if ( abs(genParticle.pdgId()) != 321 and  // Kaon
	   abs(genParticle.pdgId()) != 211 and  // Pion
	   abs(genParticle.pdgId()) != 2212 )  // Proton
	continue;
      if (deltaR(genParticle, muon) > maxDR_) continue;
      if (fabs(muon.pt()/genParticle.pt()-1) > maxDPt_) continue;
      return true;
    }
  }
  return false;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonFakeFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonFakeFilter);
