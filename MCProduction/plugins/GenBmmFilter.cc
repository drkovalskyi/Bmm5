// -*- C++ -*-
//
// Package:    Bmm5/GenBmmFilter
// Class:      GenBmmFilter
// 
/**\class GenBmmFilter GenBmmFilter.cc Bmm5/NanoAOD/plugins/GenBmmFilter.cc

 Description: Bmm5 analysis gen filter for selecting two near by muons in QCD events

 Implementation:
    - muon kinematic cuts (pt, eta)
    - dimuon doca
    - dimuon mass
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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// #include "DataFormats/MuonReco/interface/Muon.h"
// #include "DataFormats/MuonReco/interface/MuonSelectors.h"
// #include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
// #include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

//
// class declaration
//

class GenBmmFilter : public edm::stream::EDFilter<> {
public:
  explicit GenBmmFilter(const edm::ParameterSet&);
  ~GenBmmFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool 
  isGoodMuon(const HepMC::GenParticle&);

  float 
  distanceOfClosestApproach( const HepMC::GenParticle* track1,
			     const HepMC::GenParticle* track2);


  virtual bool filter(edm::Event&, const edm::EventSetup&) override;

  int charge(int pdg_id);

  //virtual void endStream() override;
  //virtual void beginStream(edm::StreamID) override;
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::HepMCProduct> hep_mc_token_;
  double max_doca_;
  double min_mu_pt_;
  double max_mu_eta_;
  double min_mm_mass_;
  double max_mm_mass_;
  const MagneticField* bField_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
};

bool 
GenBmmFilter::isGoodMuon(const HepMC::GenParticle& cand){
  if ( abs(cand.pdg_id()) != 13 ) return false;
  if ( cand.status() != 1 ) return false; // keep only stable particles
  if ( cand.momentum().perp() < min_mu_pt_) return false;
  if ( fabs(cand.momentum().eta()) > max_mu_eta_) return false;
  return true;
}

int 
GenBmmFilter::charge(int pdg_id){
  if (abs(pdg_id)==11 or abs(pdg_id)==13 or abs(pdg_id)==15){
    if (pdg_id>0) 
      return -1;
    else
      return +1;
  } else {
    throw cms::Exception("Fatal") << "Unsupported pdg_id";
  }
}

float 
GenBmmFilter::distanceOfClosestApproach( const HepMC::GenParticle* track1,
					 const HepMC::GenParticle* track2)
{
  TwoTrackMinimumDistance md;

  const HepMC::ThreeVector& trk1_point(track1->production_vertex()->point3d());
  GlobalPoint trk1_pos(trk1_point.x()/10, trk1_point.y()/10, trk1_point.z()/10); // to cm
  GlobalVector trk1_mom(track1->momentum().px(),track1->momentum().py(),track1->momentum().pz());
  GlobalTrajectoryParameters trk1(trk1_pos,trk1_mom,charge(track1->pdg_id()),bField_);

  const HepMC::ThreeVector& trk2_point(track2->production_vertex()->point3d());
  GlobalPoint trk2_pos(trk2_point.x()/10, trk2_point.y()/10, trk2_point.z()/10); // to cm
  GlobalVector trk2_mom(track2->momentum().px(),track2->momentum().py(),track2->momentum().pz());
  GlobalTrajectoryParameters trk2(trk2_pos,trk2_mom,charge(track2->pdg_id()),bField_);

  if ( not md.calculate( trk1, trk2 ) ) return -1.0;
  return -1.0;
}


GenBmmFilter::GenBmmFilter(const edm::ParameterSet& iConfig):
  hep_mc_token_( consumes<edm::HepMCProduct>(edm::InputTag("generator", "unsmeared")) ),
  max_doca_(   iConfig.getParameter<double>( "max_doca" ) ),
  min_mu_pt_(  iConfig.getParameter<double>( "min_mu_pt" ) ),
  max_mu_eta_( iConfig.getParameter<double>( "max_mu_eta" ) ),
  min_mm_mass_(  iConfig.getParameter<double>( "min_mm_mass" ) ),
  max_mm_mass_(  iConfig.getParameter<double>( "max_mm_mass" ) ),
  bField_(nullptr),
  bFieldToken_(esConsumes())
{}


GenBmmFilter::~GenBmmFilter(){}

namespace{
  double compute_mass(const HepMC::FourVector& v1, const HepMC::FourVector& v2){
    double m2 = pow(v1.e()+v2.e(),2)-pow(v1.px()+v2.px(),2)-pow(v1.py()+v2.py(),2)-pow(v1.pz()+v2.pz(),2);
    return m2>0?sqrt(m2):-1;
  }
}

bool
GenBmmFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bField_ = &iSetup.getData(bFieldToken_);
  edm::Handle<edm::HepMCProduct> hep_mc_handle;
  iEvent.getByToken(hep_mc_token_, hep_mc_handle);
  
  std::vector<const HepMC::GenParticle*> good_muons;

  for (auto cand = hep_mc_handle->GetEvent()->particles_begin(); 
       cand != hep_mc_handle->GetEvent()->particles_end(); ++cand){
    if (isGoodMuon(**cand)) good_muons.push_back(*cand);
  }

  unsigned int nMuons = good_muons.size();
  if ( nMuons > 1 ){
    for (unsigned int i = 0; i < nMuons-1; ++i) {
      for (unsigned int j = i+1; j < nMuons; ++j) {
	double mass = compute_mass(good_muons.at(i)->momentum(),good_muons.at(j)->momentum());
	if ( mass > max_mm_mass_ or mass < min_mm_mass_) continue;
  	float mm_doca = distanceOfClosestApproach(good_muons.at(i),good_muons.at(j));
  	if (max_doca_>0 and mm_doca > max_doca_) continue;
  	return true;
      }
    }
  }
  return false;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenBmmFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenBmmFilter);
