// -*- C++ -*-
//
//
/*

 Description: Copy pat muon collection updating soft muon MVA id and value

 Implementation:
     This producer is meant to be used to make an updated version of 
     slimmedMuons recomputing Soft Muon MVA.
*/
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "PhysicsTools/PatAlgos/interface/SoftMuonMvaEstimator.h"
#include "PhysicsTools/PatAlgos/interface/MuonMvaIDEstimator.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class MuonWithSoftMvaProducer : public edm::stream::EDProducer<> {
public:
  explicit MuonWithSoftMvaProducer(const edm::ParameterSet&);
  ~MuonWithSoftMvaProducer();

private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  std::unique_ptr<const pat::SoftMuonMvaEstimator> softMuonMvaEstimator_;
  std::unique_ptr<const pat::MuonMvaIDEstimator> muonMvaIDEstimator_;
  bool computeMuonIDMVA_ = false;
};


MuonWithSoftMvaProducer::MuonWithSoftMvaProducer(const edm::ParameterSet& iConfig)
{
  muonToken_ = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("input"));

  edm::FileInPath softMvaTrainingFile = iConfig.getParameter<edm::FileInPath>("softMvaTrainingFile");
  softMuonMvaEstimator_ = std::make_unique<pat::SoftMuonMvaEstimator>(softMvaTrainingFile);

  // Muon POG MVA Id for pt > 10 GeV
  computeMuonIDMVA_ = iConfig.getParameter<bool>("computeMuonIDMVA");
  if ( computeMuonIDMVA_ ) {
    edm::FileInPath mvaIDTrainingFile = iConfig.getParameter<edm::FileInPath>("mvaIDTrainingFile");
    muonMvaIDEstimator_ = std::make_unique<pat::MuonMvaIDEstimator>(mvaIDTrainingFile);
  }
  
  produces<std::vector<pat::Muon>>();
}


MuonWithSoftMvaProducer::~MuonWithSoftMvaProducer()
{
}


void
MuonWithSoftMvaProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);
  
  auto updated_muons  = std::make_unique<std::vector<pat::Muon>>();
  for (auto muon: *muons){
    float mva = softMuonMvaEstimator_->computeMva(muon);
    muon.setSoftMvaValue(mva);
    muon.setSelector(reco::Muon::SoftMvaId,  muon.softMvaValue() >   0.58  ); //WP choose for bmm4

    // Muon POG MVA id for muons with pt > 10
    float mvaID = 0.0;
    if (computeMuonIDMVA_) {
      if (muon.isLooseMuon()) {
        mvaID = muonMvaIDEstimator_->computeMVAID(muon)[1];
      }
    }
    muon.setMvaIDValue(mvaID);
    
    updated_muons->push_back(muon);
  }
  iEvent.put(std::move(updated_muons));
}


void
MuonWithSoftMvaProducer::beginStream(edm::StreamID)
{
}


void
MuonWithSoftMvaProducer::endStream()
{
}


DEFINE_FWK_MODULE(MuonWithSoftMvaProducer);
