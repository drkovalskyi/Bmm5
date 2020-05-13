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
// V0Producer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;

namespace {
  const float muon_mass_     = 0.10565837;
  const float kaon_mass_     = 0.493677;
  const float mass_err_      = 1.6e-5;
  const float pion_mass_     = 0.139570;
  const float jpsi_mass_     = 3.0969;
};

struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, cosAlpha;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), lxyErr(-1.0), sigLxy(-1.0), cosAlpha(-999.)
  {}

  bool valid() const {
    return treeIsValid and vertexIsValid;
  }

  void postprocess(const reco::BeamSpot& beamSpot)
  {
    if ( not valid() ) return;
    // displacement information
    TVector v(2);
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();

    TMatrix errVtx(2,2);
    errVtx(0,0) = refitVertex->error().cxx();
    errVtx(0,1) = refitVertex->error().matrix()(0,1);
    errVtx(1,0) = errVtx(0,1);
    errVtx(1,1) = refitVertex->error().cyy();

    TMatrix errBS(2,2);
    errBS(0,0) = beamSpot.covariance()(0,0);
    errBS(0,1) = beamSpot.covariance()(0,1);
    errBS(1,0) = beamSpot.covariance()(1,0);
    errBS(1,1) = beamSpot.covariance()(1,1);
    
    lxy = sqrt(v.Norm2Sqr());
    lxyErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / lxy;
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
    
    // compute cosAlpha 2D wrt BeamSpot
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();
    TVector w(2);
    w[0] = refitMother->currentState().globalMomentum().x();
    w[1] = refitMother->currentState().globalMomentum().y();
    cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());

  }
  
  float mass() const
  {
    if ( not valid() ) return -1.0;
    return refitMother->currentState().mass();
  }

  float refit_mass(unsigned int i, unsigned int j) const
  {
    if ( not valid() ) return -1.0;
    if (i >= refitDaughters.size()) return -2.0;
    if (j >= refitDaughters.size()) return -3.0;
    if (refitDaughters.at(i)->currentState().globalMomentum().mag2()<0) return -4.0;
    if (refitDaughters.at(j)->currentState().globalMomentum().mag2()<0) return -5.0;
    auto momentum = refitDaughters.at(i)->currentState().globalMomentum() + 
      refitDaughters.at(j)->currentState().globalMomentum();
    auto energy1 = sqrt(refitDaughters.at(i)->currentState().globalMomentum().mag2() + 
			pow(refitDaughters.at(i)->currentState().mass(),2));
    auto energy2 = sqrt(refitDaughters.at(j)->currentState().globalMomentum().mag2() + 
			pow(refitDaughters.at(j)->currentState().mass(),2));
    return sqrt(pow(energy1+energy2,2)-momentum.mag2());
  }

  GlobalVector p3() const
  {
    if ( not valid() ) return GlobalVector();
    return refitMother->currentState().globalMomentum();
  }

  GlobalVector dau_p3(unsigned int i) const
  {
    if ( not valid() or i>=refitDaughters.size() ) return GlobalVector();
    return refitDaughters.at(i)->currentState().globalMomentum();
  }

  float massErr() const
  {
    if ( not valid() ) return -1.0;
    return sqrt(refitMother->currentState().kinematicParametersError().matrix()(6,6));
  }

  float chi2() const
  {
    if ( not valid() ) return -1.0;
    return refitVertex->chiSquared();
  }

  float ndof() const
  {
    return refitVertex->degreesOfFreedom();
  }

  float vtxProb() const
  {
    if ( not valid() ) return -1.0;
    return TMath::Prob((double)refitVertex->chiSquared(), int(rint(refitVertex->degreesOfFreedom())));
  }
  
};

struct GenKsMatchInfo{
  const pat::PackedGenParticle* mc_trk1;
  const pat::PackedGenParticle* mc_trk2;
  const reco::Candidate*   match;
  GenKsMatchInfo():mc_trk1(0), mc_trk2(0), match(0)
  {}
};

using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class V0Producer : public edm::EDProducer {
    
public:
    
  explicit V0Producer(const edm::ParameterSet &iConfig);
    
  ~V0Producer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isGoodMuon(const pat::Muon& muon);
  bool isGoodTrack(const pat::PackedCandidate& track);
  bool isDisplacedTrack(const pat::PackedCandidate& track);
  bool isGoodMuonProbe(const pat::PackedCandidate& track);
  bool isGoodPion(const pat::PackedCandidate& track);
  bool isGoodPair(const pat::PackedCandidate& track1,
		  const pat::PackedCandidate& track2);
    
  KinematicFitResult 
  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
   				 reco::BeamSpot beamSpot);
  GenKsMatchInfo getGenKsMatchInfo( const pat::PackedCandidate& track1,
				    const pat::PackedCandidate& track2);
  // Two track DOCA
  float 
  distanceOfClosestApproach( const reco::Track* track1,
			     const reco::Track* track2 );
  float 
  distanceOfClosestApproach( const pat::PackedGenParticle* track1,
			     const pat::PackedGenParticle* track2);
  // Track to vertex DOCA
  Measurement1D
  distanceOfClosestApproach( const reco::Track* track,
			     RefCountedKinematicVertex vertex);
  Measurement1D 
  distanceOfClosestApproach( const reco::Track* track,
			     const reco::Vertex& vertex);

  KinematicFitResult 
  fillKsInfo(pat::CompositeCandidate& ksCand,
	     const edm::Event& iEvent,
	     const pat::PackedCandidate & pion1,
	     const pat::PackedCandidate & pion2); 

  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;

  bool isMC_;

  double minMuonPt_;
  double maxMuonEta_;
  double minPionPt_;
  double maxPionEta_;

  double minKsPreselectMass_;
  double maxKsPreselectMass_;
  double minKsMass_;
  double maxKsMass_;
  double maxTwoTrackDOCA_;
  double minDisplaceTrackSignificance_;
  double maxLxy_;
  double minSigLxy_;
  double minCosAlpha_;
  double minVtxProb_;
};

V0Producer::V0Producer(const edm::ParameterSet &iConfig):
beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
beamSpot_(nullptr),
muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
pfCandToken_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( iConfig.getParameter<edm::InputTag>( "packedGenParticleCollection" ) ) ),
packedGenParticles_(nullptr),
impactPointExtrapolator_(0),
isMC_(             iConfig.getParameter<bool>( "isMC" ) ),
minMuonPt_(         iConfig.getParameter<double>( "minMuonPt" ) ),
maxMuonEta_(        iConfig.getParameter<double>( "maxMuonEta" ) ),
minPionPt_(         iConfig.getParameter<double>( "minPionPt" ) ),
maxPionEta_(        iConfig.getParameter<double>( "maxPionEta" ) ),
minKsPreselectMass_(     iConfig.getParameter<double>( "minKsPreselectMass" ) ),
maxKsPreselectMass_(     iConfig.getParameter<double>( "maxKsPreselectMass" ) ),
minKsMass_(     iConfig.getParameter<double>( "minKsMass" ) ),
maxKsMass_(     iConfig.getParameter<double>( "maxKsMass" ) ),
maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
minDisplaceTrackSignificance_( iConfig.getParameter<double>( "minDisplaceTrackSignificance" ) ),
maxLxy_( iConfig.getParameter<double>( "maxLxy" ) ),
minSigLxy_( iConfig.getParameter<double>( "minSigLxy" ) ),
minCosAlpha_( iConfig.getParameter<double>( "minCosAlpha" ) ),
minVtxProb_( iConfig.getParameter<double>( "minVtxProb" ) )
{
    produces<pat::CompositeCandidateCollection>("Ks");
}

bool V0Producer::isGoodTrack(const pat::PackedCandidate& track){
  if (track.charge() == 0 ) return false;
  if ( not track.hasTrackDetails() ) return false;
  if ( not track.bestTrack()->quality(reco::Track::highPurity) ) return false; 
  return true;
}

bool V0Producer::isDisplacedTrack(const pat::PackedCandidate& track){
  double sigDxy = track.bestTrack()->dxyError()>0 ? fabs(track.bestTrack()->dxy(*beamSpot_))/track.bestTrack()->dxyError():0.0;
  return sigDxy > minDisplaceTrackSignificance_;
}

bool V0Producer::isGoodMuonProbe(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxMuonEta_ and track.pt()>minMuonPt_;
}

bool V0Producer::isGoodPion(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxPionEta_ and track.pt()>minPionPt_;
}

bool V0Producer::isGoodPair(const pat::PackedCandidate& track1,
			    const pat::PackedCandidate& track2){
  return (isGoodMuonProbe(track1) and isGoodPion(track2)) or 
    (isGoodMuonProbe(track2) and isGoodPion(track1));
}

namespace {
  void addFitInfo( pat::CompositeCandidate& cand, const KinematicFitResult& fit, std::string name ){
    cand.addUserInt(   name+"_valid",       fit.valid() );
    cand.addUserFloat( name+"_vtx_prob",    fit.vtxProb() );
    cand.addUserFloat( name+"_vtx_chi2dof", fit.chi2()>0?fit.chi2()/fit.ndof():-1);
    cand.addUserFloat( name+"_mass",        fit.mass() );
    cand.addUserFloat( name+"_massErr",     fit.massErr() );
    cand.addUserFloat( name+"_lxy",         fit.lxy );
    cand.addUserFloat( name+"_sigLxy",      fit.sigLxy );
    cand.addUserFloat( name+"_cosAlphaXY",  fit.cosAlpha );
    cand.addUserFloat( name+"_pt",          fit.p3().perp() );
    cand.addUserFloat( name+"_eta",         fit.p3().eta() );
    cand.addUserFloat( name+"_phi",         fit.p3().phi() );
  }    
}

KinematicFitResult 
V0Producer::fillKsInfo(pat::CompositeCandidate& ksCand,
		       const edm::Event& iEvent,
		       const pat::PackedCandidate & pion1,
		       const pat::PackedCandidate & pion2) 
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( pion1.bestTrack() );
  masses.push_back(pion_mass_);
  trks.push_back( pion2.bestTrack() );
  masses.push_back(pion_mass_);
  auto vtxFit = vertexWithKinematicFitter(trks,masses);
  vtxFit.postprocess(*beamSpot_);

  // printf("vtxFit (x,y,z): (%7.3f,%7.3f,%7.3f)\n", 
  // 	 vtxFit.refitVertex->position().x(),
  // 	 vtxFit.refitVertex->position().y(),
  // 	 vtxFit.refitVertex->position().z());
  addFitInfo(ksCand, vtxFit, "kin");
  
  if (isMC_){
    auto gen_tt = getGenKsMatchInfo( pion1, pion2 );
    if (gen_tt.mc_trk1){
      ksCand.addUserInt(  "gen_trk1_pdgId",   gen_tt.mc_trk1->pdgId() );
      ksCand.addUserInt(  "gen_trk1_mpdgId",  gen_tt.mc_trk1->numberOfMothers()>0?gen_tt.mc_trk1->mother(0)->pdgId():0 );
      ksCand.addUserFloat("gen_trk1_pt",      gen_tt.mc_trk1->pt() );
    } else {
      ksCand.addUserInt(  "gen_trk1_pdgId",   0 );
      ksCand.addUserInt(  "gen_trk1_mpdgId",  0 );
      ksCand.addUserFloat("gen_trk1_pt",      0 );
    }
    if (gen_tt.mc_trk2){
      ksCand.addUserInt(  "gen_trk2_pdgId",   gen_tt.mc_trk2->pdgId() );
      ksCand.addUserInt(  "gen_trk2_mpdgId",  gen_tt.mc_trk2->numberOfMothers()>0?gen_tt.mc_trk2->mother(0)->pdgId():0 );
      ksCand.addUserFloat("gen_trk2_pt",      gen_tt.mc_trk2->pt() );
    } else {
      ksCand.addUserInt(  "gen_trk2_pdgId",   0 );
      ksCand.addUserInt(  "gen_trk2_mpdgId",  0 );
      ksCand.addUserFloat("gen_trk2_pt",      0 );
    }
    if (gen_tt.match){
      ksCand.addUserFloat("gen_mass",         gen_tt.match->mass() );
      ksCand.addUserFloat("gen_pt",           gen_tt.match->pt() );
      ksCand.addUserInt(  "gen_pdgId",        gen_tt.match->pdgId() );
    } else {
      ksCand.addUserFloat("gen_mass",         0 );
      ksCand.addUserFloat("gen_pt",           0 );
      ksCand.addUserInt(  "gen_pdgId",        0 );
    }
  }
  return vtxFit;
}

namespace{
  int match_to_muon(const pat::PackedCandidate& pfCand,
		    const std::vector<pat::Muon>& muons){
    if ( abs(pfCand.pdgId())!=13 ) return -1;
    for ( unsigned int i=0; i<muons.size(); ++i ){
      if ( deltaR(muons.at(i), pfCand) < 0.01 ) return i;
    }
    return -1;
  }
}

void V0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder_);

    AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle_.product());
    impactPointExtrapolator_ = &extrapolator;

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("V0Producer") << "No beam spot available from EventSetup" ;
    }
    
    beamSpot_ = beamSpotHandle.product();

    edm::Handle<std::vector<pat::Muon>> muonHandle;
    
    iEvent.getByToken(muonToken_, muonHandle);
    iEvent.getByToken(pfCandToken_, pfCandHandle_);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      packedGenParticles_ = nullptr;
    }

    // auto nMuons = muonHandle->size();
    auto nPFCands = pfCandHandle_->size();
    
    // Output collection
    auto kss  = std::make_unique<pat::CompositeCandidateCollection>();
    AddFourMomenta addP4;

    // Build V0 candidates first

    if ( nPFCands > 1 ){
      for ( unsigned int i = 0; i < nPFCands-1; ++i ) {
	pat::PackedCandidate pfCand1( (*pfCandHandle_)[i] );
	if ( not isGoodTrack(pfCand1) ) continue;
	for ( unsigned int j = i+1; j < nPFCands; ++j ) {
	  pat::PackedCandidate pfCand2( (*pfCandHandle_)[j] );
	  if ( pfCand1.charge()*pfCand2.charge() >= 0 ) continue;
	  if ( not isGoodTrack(pfCand2) ) continue;
	  if ( not isGoodPair(pfCand1,pfCand2) ) continue;

	  // KsToPiPi
	  if ( not isDisplacedTrack(pfCand1) or 
	       not isDisplacedTrack(pfCand2) ) continue;

	  pat::CompositeCandidate ksCand;
	  ksCand.addDaughter( pfCand1 , "pion1" );
	  ksCand.addDaughter( pfCand2 , "pion2" );
	  addP4.set( ksCand );

	  if ( ksCand.mass() < minKsPreselectMass_ or  
	       ksCand.mass() > maxKsPreselectMass_ ) continue;

	  auto tt_doca = distanceOfClosestApproach(pfCand1.bestTrack(),
						   pfCand2.bestTrack());
	  if ( maxTwoTrackDOCA_>0 and tt_doca > maxTwoTrackDOCA_ ) continue;
	  ksCand.addUserFloat( "doca", tt_doca);
	  ksCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
	  ksCand.addUserFloat( "trk1_eta", pfCand1.eta() );
	  ksCand.addUserFloat( "trk1_phi", pfCand1.phi() );
	  ksCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
	  ksCand.addUserFloat( "trk2_eta", pfCand2.eta() );
	  ksCand.addUserFloat( "trk2_phi", pfCand2.phi() );
	  ksCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle));
	  ksCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle));

	  auto ksVtxFit = fillKsInfo(ksCand,iEvent,pfCand1,pfCand2);
	  
	  if ( not ksVtxFit.valid() ) continue;
	  if ( ksVtxFit.vtxProb() < minVtxProb_) continue;
	  if ( ksVtxFit.mass() < minKsMass_ or  
	       ksVtxFit.mass() > maxKsMass_ ) continue;
	  if ( ksVtxFit.lxy > maxLxy_ ) continue;
	  if ( ksVtxFit.sigLxy < minSigLxy_ ) continue;
	  
	  kss->push_back(ksCand);
	}
      }
    }
    
    iEvent.put(std::move(kss), "Ks");
}

KinematicFitResult 
V0Producer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
					    std::vector<float> masses)
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
  if ( trks.size() != masses.size() ) 
    throw cms::Exception("Error") << "number of tracks and number of masses should match";

  std::vector<reco::TransientTrack> transTrks;

  KinematicParticleFactoryFromTransientTrack factory;
  KinematicParticleVertexFitter fitter;
    
  std::vector<RefCountedKinematicParticle> particles;

  double chi = 0.;
  double ndf = 0.;
  float mass_err(mass_err_);
  for (unsigned int i=0; i<trks.size(); ++i){
    transTrks.push_back((*theTTBuilder_).build(trks[i]));
    particles.push_back(factory.particle(transTrks.back(),masses[i],chi,ndf,mass_err));
  }

  RefCountedKinematicTree vertexFitTree = fitter.fit(particles);
  KinematicFitResult result;
    
  if ( !vertexFitTree->isValid()) return result;
  
  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;
  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();
  
  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}


pair<double,double> V0Producer::computeDCA(const pat::PackedCandidate &kaon,
                                                 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle_));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}
namespace{

  bool dr_match(const LorentzVector& reco , const LorentzVector& gen){
    if (fabs(reco.pt()-gen.pt())/gen.pt()<0.1 and deltaR(reco,gen)<0.02)
      return true;
    return false;
  }

  std::vector<unsigned int> 
    get_depth_from_permutation(const std::vector<unsigned int>& elements){
    std::vector<unsigned int> result;
    unsigned int counter(0);
    for (auto element: elements){
      if (element==0){
	counter++;
      } else {
	result.push_back(counter);
	counter = 0;
      }
    }
    result.push_back(counter);
    return result;
  }

  bool is_acceptable(const reco::Candidate* cand){
    if ( not cand) return false; 
    // skip quarks
    if ( abs(cand->pdgId())<10 ) return false;
    // skip protons
    if ( abs(cand->pdgId())==2212 ) return false;
    // skip gluons
    if ( abs(cand->pdgId())==21 ) return false;
    return true;
  }

  // depth 0 - first mother

  const reco::Candidate* get_mother(const reco::Candidate* cand, unsigned int depth){
    if (not cand) return 0;
    const reco::Candidate* mother = cand->mother();
    unsigned int i = 0;
    while ( is_acceptable(mother) and i<depth ){
      i++;
      mother = mother->mother();
    }
    if (is_acceptable(mother))
      return mother;
    else
      return 0;
  }

  const reco::Candidate* 
    find_common_ancestor(const std::vector<const reco::Candidate*>& particles, 
			 unsigned int max_depth=10){
    auto n = particles.size();
    for (unsigned int depth=0; depth<max_depth; ++depth){
      // make a list of elements (0) and separators (1) and
      // find all possible permutations of the elements
      std::vector<unsigned int> elements;
      for (unsigned int i=0; i<depth; ++i)
	elements.push_back(0);
      for (unsigned int i=0; i<n-1; ++i)
	elements.push_back(1);
      do {
	auto depth_vector = get_depth_from_permutation(elements);
	const reco::Candidate* common_mother(0);
	for (unsigned int i=0; i<n; ++i){
	  auto mother = get_mother(particles[i],depth_vector[i]);
	  if (not mother) {
	    common_mother = 0;
	    break;
	  }
	  if (not common_mother) common_mother = mother;
	  if (common_mother != mother) {
	    common_mother = 0;
	    break;
	  }	  
	}
	if (common_mother) return common_mother;
      } while(std::next_permutation(elements.begin(), elements.end()));
    }
    return 0;
  }

}

GenKsMatchInfo V0Producer::getGenKsMatchInfo( const pat::PackedCandidate& track1,
					      const pat::PackedCandidate& track2 )
{
  GenKsMatchInfo result;
  std::vector<const reco::Candidate*> daughters;
  for (auto const & genParticle: *packedGenParticles_){
    if (dr_match(track1.p4(),genParticle.p4())){
      result.mc_trk1 = &genParticle;
      daughters.push_back(result.mc_trk1);
    }
    if (dr_match(track2.p4(),genParticle.p4())){
      result.mc_trk2 = &genParticle;
      daughters.push_back(result.mc_trk2);
    }
  }
  if (daughters.size()==2){
    const auto* mother = find_common_ancestor(daughters);
    if (mother) result.match        = mother;
  }
  return result;
}

float V0Producer::distanceOfClosestApproach( const reco::Track* track1,
					     const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
V0Producer::distanceOfClosestApproach( const reco::Track* track,
					     RefCountedKinematicVertex vertex)
{
  if (not vertex->vertexIsValid()) return Measurement1D(-1.0,-1.0);
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), vertex->position());
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
  return doca;
}

Measurement1D 
V0Producer::distanceOfClosestApproach( const reco::Track* track,
					     const reco::Vertex& vertex)
{
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), GlobalPoint(Basic3DVector<float>(vertex.position())));
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex);
  return doca;
}


DEFINE_FWK_MODULE(V0Producer);

//  LocalWords:  vertices
