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
// BmmV0Producer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;

namespace {
  const float muon_mass_     = 0.10565837;
  const float kaon_mass_     = 0.493677;
  const float mass_err_      = 1.6e-5;
  const float pion_mass_     = 0.139570;
  const float jpsi_mass_     = 3.0969;
  const float proton_mass_   = 0.9382;
};

struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, cosAlpha, ipSigBS, ipSigPV;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), lxyErr(-1.0), sigLxy(-1.0), cosAlpha(-999.),
		       ipSigBS(999.), ipSigPV(999.)
  {}

  bool valid() const {
    return treeIsValid and vertexIsValid;
  }

  void postprocess(const reco::BeamSpot& beamSpot,
		   const reco::VertexCollection& vertices)
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

    // comput impact parameter
    auto transientTrack = refitMother->refittedTransientTrack();
    ipSigBS = transientTrack.stateAtBeamLine().transverseImpactParameter().significance();

    ipSigPV = 999.;
    for (const auto& vertex: vertices){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(transientTrack, vertex);
      if (impactParameter3D.first)
	if (ipSigPV > impactParameter3D.second.significance())
	  ipSigPV = impactParameter3D.second.significance();
    }
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

struct GenMatchInfo{
  const pat::PackedGenParticle* mc_trk1;
  const pat::PackedGenParticle* mc_trk2;
  const reco::Candidate*   match;
  GenMatchInfo():mc_trk1(0), mc_trk2(0), match(0)
  {}
};

using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class BmmV0Producer : public edm::EDProducer {
    
public:
    
  explicit BmmV0Producer(const edm::ParameterSet &iConfig);
    
  ~BmmV0Producer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isGoodMuon(const pat::Muon& muon);
  bool isGoodTrack(const pat::PackedCandidate& track);
  bool displacedTrack(const pat::PackedCandidate& track);
  bool isGoodMuonProbe(const pat::PackedCandidate& track);
  bool isGoodPion(const pat::PackedCandidate& track);
  bool isGoodPair(const pat::PackedCandidate& track1,
		  const pat::PackedCandidate& track2);
    
  float
  trackImpactParameterSignificance( const pat::PackedCandidate& track);
  
  pat::CompositeCandidate
  getKsToPiPi(const edm::Event& iEvent,
	      const pat::PackedCandidate& pfCand1,
	      const pat::PackedCandidate& pfCand2);
  pat::CompositeCandidate
  getD0ToKPi(const edm::Event& iEvent,
	     const pat::PackedCandidate& kaonCand,
	     const pat::PackedCandidate& pion);
  pat::CompositeCandidate
  getPhiToKK(const edm::Event& iEvent,
	     const pat::PackedCandidate& pfCand1,
	     const pat::PackedCandidate& pfCand2);
  pat::CompositeCandidate
  getLambdaToPPi(const edm::Event& iEvent,
		 const pat::PackedCandidate& pfCand1,
		 const pat::PackedCandidate& pfCand2);

  KinematicFitResult 
  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
   				 reco::BeamSpot beamSpot);
  GenMatchInfo getGenMatchInfo( const pat::PackedCandidate& track1,
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
  fillInfo(pat::CompositeCandidate& ksCand,
	   const edm::Event& iEvent,
	   const pat::PackedCandidate & cand1,
	   const pat::PackedCandidate & cand2, 
	   string cand1_name = "trk1",
	   string cand2_name = "trk2");

  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const reco::VertexCollection* primaryVertices_;

  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_;
  edm::Handle<std::vector<pat::Muon>> muonHandle_;

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
  double minPhiPreselectMass_;
  double maxPhiPreselectMass_;
  double minPhiMass_;
  double maxPhiMass_;
  double minD0PreselectMass_;
  double maxD0PreselectMass_;
  double minD0Mass_;
  double maxD0Mass_;
  double minLambdaPreselectMass_;
  double maxLambdaPreselectMass_;
  double minLambdaMass_;
  double maxLambdaMass_;
  double maxTwoTrackDOCA_;
  double minDisplaceTrackSignificance_;
  double maxLxy_;
  double minSigLxy_;
  double minCosAlpha_;
  double minVtxProb_;
};

BmmV0Producer::BmmV0Producer(const edm::ParameterSet &iConfig):
beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
beamSpot_(nullptr),
vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
primaryVertices_(nullptr),
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
minPhiPreselectMass_(     iConfig.getParameter<double>( "minPhiPreselectMass" ) ),
maxPhiPreselectMass_(     iConfig.getParameter<double>( "maxPhiPreselectMass" ) ),
minPhiMass_(     iConfig.getParameter<double>( "minPhiMass" ) ),
maxPhiMass_(     iConfig.getParameter<double>( "maxPhiMass" ) ),
minD0PreselectMass_(     iConfig.getParameter<double>( "minD0PreselectMass" ) ),
maxD0PreselectMass_(     iConfig.getParameter<double>( "maxD0PreselectMass" ) ),
minD0Mass_(     iConfig.getParameter<double>( "minD0Mass" ) ),
maxD0Mass_(     iConfig.getParameter<double>( "maxD0Mass" ) ),
minLambdaPreselectMass_(     iConfig.getParameter<double>( "minLambdaPreselectMass" ) ),
maxLambdaPreselectMass_(     iConfig.getParameter<double>( "maxLambdaPreselectMass" ) ),
minLambdaMass_(     iConfig.getParameter<double>( "minLambdaMass" ) ),
maxLambdaMass_(     iConfig.getParameter<double>( "maxLambdaMass" ) ),
maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
minDisplaceTrackSignificance_( iConfig.getParameter<double>( "minDisplaceTrackSignificance" ) ),
maxLxy_( iConfig.getParameter<double>( "maxLxy" ) ),
minSigLxy_( iConfig.getParameter<double>( "minSigLxy" ) ),
minCosAlpha_( iConfig.getParameter<double>( "minCosAlpha" ) ),
minVtxProb_( iConfig.getParameter<double>( "minVtxProb" ) )
{
    produces<pat::CompositeCandidateCollection>("Ks");
    produces<pat::CompositeCandidateCollection>("D0");
    produces<pat::CompositeCandidateCollection>("Phi");
    produces<pat::CompositeCandidateCollection>("Lambda");
}

bool BmmV0Producer::isGoodTrack(const pat::PackedCandidate& track){
  if ( track.charge() == 0 ) return false;
  if ( not track.hasTrackDetails() ) return false;
  if ( not track.bestTrack()->quality(reco::Track::highPurity) ) return false; 
  return true;
}

float BmmV0Producer::trackImpactParameterSignificance(const pat::PackedCandidate& track){
  return track.bestTrack()->dxyError()>0 ? fabs(track.bestTrack()->dxy(*beamSpot_))/track.bestTrack()->dxyError():0.0;
}

bool BmmV0Producer::displacedTrack(const pat::PackedCandidate& track){
  return trackImpactParameterSignificance(track) > minDisplaceTrackSignificance_;
}

bool BmmV0Producer::isGoodMuonProbe(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxMuonEta_ and track.pt()>minMuonPt_;
}

bool BmmV0Producer::isGoodPion(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxPionEta_ and track.pt()>minPionPt_;
}

bool BmmV0Producer::isGoodPair(const pat::PackedCandidate& track1,
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
    cand.addUserFloat( name+"_sipBS",       fit.ipSigBS );
    cand.addUserFloat( name+"_sipPV",       fit.ipSigPV );
    cand.addUserFloat( name+"_pt",          fit.p3().perp() );
    cand.addUserFloat( name+"_eta",         fit.p3().eta() );
    cand.addUserFloat( name+"_phi",         fit.p3().phi() );
  }    
}

KinematicFitResult 
BmmV0Producer::fillInfo(pat::CompositeCandidate& v0Cand,
			const edm::Event& iEvent,
			const pat::PackedCandidate & cand1,
			const pat::PackedCandidate & cand2,
			string cand1_name,
			string cand2_name)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back(cand1.bestTrack());
  trks.push_back(cand2.bestTrack());
  masses.push_back(cand1.mass());
  masses.push_back(cand2.mass());
  auto vtxFit = vertexWithKinematicFitter(trks, masses);
  vtxFit.postprocess(*beamSpot_, *primaryVertices_);

  // printf("vtxFit (x,y,z): (%7.3f,%7.3f,%7.3f)\n", 
  // 	 vtxFit.refitVertex->position().x(),
  // 	 vtxFit.refitVertex->position().y(),
  // 	 vtxFit.refitVertex->position().z());
  addFitInfo(v0Cand, vtxFit, "kin");
  
  if (isMC_){
    auto gen_tt = getGenMatchInfo( cand1, cand2 );
    if (gen_tt.mc_trk1){
      v0Cand.addUserInt(  "gen_" + cand1_name + "_pdgId",   gen_tt.mc_trk1->pdgId() );
      v0Cand.addUserInt(  "gen_" + cand1_name + "_mpdgId",  gen_tt.mc_trk1->numberOfMothers()>0?gen_tt.mc_trk1->mother(0)->pdgId():0 );
      v0Cand.addUserFloat("gen_" + cand1_name + "_pt",      gen_tt.mc_trk1->pt() );
    } else {
      v0Cand.addUserInt(  "gen_" + cand1_name + "_pdgId",   0 );
      v0Cand.addUserInt(  "gen_" + cand1_name + "_mpdgId",  0 );
      v0Cand.addUserFloat("gen_" + cand1_name + "_pt",      0 );
    }
    if (gen_tt.mc_trk2){
      v0Cand.addUserInt(  "gen_" + cand2_name + "_pdgId",   gen_tt.mc_trk2->pdgId() );
      v0Cand.addUserInt(  "gen_" + cand2_name + "_mpdgId",  gen_tt.mc_trk2->numberOfMothers()>0?gen_tt.mc_trk2->mother(0)->pdgId():0 );
      v0Cand.addUserFloat("gen_" + cand2_name + "_pt",      gen_tt.mc_trk2->pt() );
    } else {
      v0Cand.addUserInt(  "gen_" + cand2_name + "_pdgId",   0 );
      v0Cand.addUserInt(  "gen_" + cand2_name + "_mpdgId",  0 );
      v0Cand.addUserFloat("gen_" + cand2_name + "_pt",      0 );
    }
    if (gen_tt.match){
      v0Cand.addUserFloat("gen_mass",         gen_tt.match->mass() );
      v0Cand.addUserFloat("gen_pt",           gen_tt.match->pt() );
      v0Cand.addUserInt(  "gen_pdgId",        gen_tt.match->pdgId() );
    } else {
      v0Cand.addUserFloat("gen_mass",         0 );
      v0Cand.addUserFloat("gen_pt",           0 );
      v0Cand.addUserInt(  "gen_pdgId",        0 );
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

void BmmV0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder_);

    AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle_.product());
    impactPointExtrapolator_ = &extrapolator;

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BmmV0Producer") << "No beam spot available from EventSetup" ;
    }
    beamSpot_ = beamSpotHandle.product();

    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByToken(vertexToken_, pvHandle);
    primaryVertices_ = pvHandle.product();
    
    iEvent.getByToken(muonToken_, muonHandle_);
    iEvent.getByToken(pfCandToken_, pfCandHandle_);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      packedGenParticles_ = nullptr;
    }

    // auto nMuons = muonHandle_->size();
    auto nPFCands = pfCandHandle_->size();
    
    // Output collection
    auto kss  = std::make_unique<pat::CompositeCandidateCollection>();
    auto d0s  = std::make_unique<pat::CompositeCandidateCollection>();
    auto phis = std::make_unique<pat::CompositeCandidateCollection>();
    auto lambdas = std::make_unique<pat::CompositeCandidateCollection>();

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

	  auto tt_doca = distanceOfClosestApproach(pfCand1.bestTrack(),
						   pfCand2.bestTrack());
	  if ( maxTwoTrackDOCA_>0 and tt_doca > maxTwoTrackDOCA_ )
	    continue;
	  
	  // PhiToKK
	  auto phiCand = getPhiToKK(iEvent, pfCand1, pfCand2);
	  if (phiCand.numberOfDaughters() > 0){
	    phiCand.addUserFloat( "doca", tt_doca);
	    phis->push_back(phiCand);
	  }

	  // Look for V0s built from displaced tracks
	  if ( not displacedTrack(pfCand1) or 
	       not displacedTrack(pfCand2) ) continue;

	  // KsToPiPi
	  auto ksCand = getKsToPiPi(iEvent, pfCand1, pfCand2);
	  if (ksCand.numberOfDaughters() > 0){
	    ksCand.addUserFloat( "doca", tt_doca);
	    kss->push_back(ksCand);
	  }

	  // D0ToKPi
	  auto d0Cand1 = getD0ToKPi(iEvent, pfCand1, pfCand2);
	  if (d0Cand1.numberOfDaughters() > 0){
	    d0Cand1.addUserFloat( "doca", tt_doca);
	    d0s->push_back(d0Cand1);
	  }
	  auto d0Cand2 = getD0ToKPi(iEvent, pfCand2, pfCand1);
	  if (d0Cand2.numberOfDaughters() > 0){
	    d0Cand2.addUserFloat( "doca", tt_doca);
	    d0s->push_back(d0Cand2);
	  }
	  
	  // LambdaToPPi
	  auto lambdaCand1 = getLambdaToPPi(iEvent, pfCand1, pfCand2);
	  if (lambdaCand1.numberOfDaughters() > 0){
	    lambdaCand1.addUserFloat( "doca", tt_doca);
	    lambdas->push_back(lambdaCand1);
	  }
	  auto lambdaCand2 = getLambdaToPPi(iEvent, pfCand2, pfCand1);
	  if (lambdaCand2.numberOfDaughters() > 0){
	    lambdaCand2.addUserFloat( "doca", tt_doca);
	    lambdas->push_back(lambdaCand2);
	  }

	}
      }
    }
    
    iEvent.put(std::move(kss), "Ks");
    iEvent.put(std::move(d0s), "D0");
    iEvent.put(std::move(phis), "Phi");
    iEvent.put(std::move(lambdas), "Lambda");
}

pat::CompositeCandidate
BmmV0Producer::getKsToPiPi(const edm::Event& iEvent,
			   const pat::PackedCandidate& ipfCand1,
			   const pat::PackedCandidate& ipfCand2)
{
  pat::CompositeCandidate ksCand;
  pat::PackedCandidate pfCand1(ipfCand1);
  pfCand1.setMass(pion_mass_);
  pat::PackedCandidate pfCand2(ipfCand2);
  pfCand2.setMass(pion_mass_);
  ksCand.addDaughter( pfCand1 , "trk1" );
  ksCand.addDaughter( pfCand2 , "trk2" );
  AddFourMomenta addP4;
  addP4.set( ksCand );

  if ( ksCand.mass() < minKsPreselectMass_ or ksCand.mass() > maxKsPreselectMass_ )
    return pat::CompositeCandidate();

  ksCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
  ksCand.addUserFloat( "trk1_eta", pfCand1.eta() );
  ksCand.addUserFloat( "trk1_phi", pfCand1.phi() );
  ksCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
  ksCand.addUserFloat( "trk2_eta", pfCand2.eta() );
  ksCand.addUserFloat( "trk2_phi", pfCand2.phi() );
  ksCand.addUserFloat( "trk1_sip", trackImpactParameterSignificance(pfCand1) );
  ksCand.addUserFloat( "trk2_sip", trackImpactParameterSignificance(pfCand2) );
  ksCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  ksCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto ksVtxFit = fillInfo(ksCand, iEvent, pfCand1, pfCand2);
	  
  if ( not ksVtxFit.valid()  or 
       ksVtxFit.vtxProb() < minVtxProb_ or
       ksVtxFit.mass() < minKsMass_ or  
       ksVtxFit.mass() > maxKsMass_ or
       ksVtxFit.lxy > maxLxy_ or 
       ksVtxFit.sigLxy < minSigLxy_ ) return pat::CompositeCandidate();

  return ksCand;
}

pat::CompositeCandidate
BmmV0Producer::getD0ToKPi(const edm::Event& iEvent,
		       const pat::PackedCandidate& ikaon,
		       const pat::PackedCandidate& ipion)
{
  pat::CompositeCandidate d0Cand;
  pat::PackedCandidate kaon(ikaon);
  kaon.setMass(kaon_mass_);
  pat::PackedCandidate pion(ipion);
  pion.setMass(pion_mass_);
  d0Cand.addDaughter( kaon , "kaon" );
  d0Cand.addDaughter( pion , "pion" );
  AddFourMomenta addP4;
  addP4.set( d0Cand );
  
  if ( d0Cand.mass() < minD0PreselectMass_ or d0Cand.mass() > maxD0PreselectMass_ )
    return pat::CompositeCandidate();

  d0Cand.addUserFloat( "kaon_pt",  kaon.pt() );
  d0Cand.addUserFloat( "kaon_eta", kaon.eta() );
  d0Cand.addUserFloat( "kaon_phi", kaon.phi() );
  d0Cand.addUserFloat( "pion_pt",  pion.pt() );
  d0Cand.addUserFloat( "pion_eta", pion.eta() );
  d0Cand.addUserFloat( "pion_phi", pion.phi() );
  d0Cand.addUserFloat( "kaon_sip", trackImpactParameterSignificance(kaon) );
  d0Cand.addUserFloat( "pion_sip", trackImpactParameterSignificance(pion) );
  d0Cand.addUserInt( "kaon_mu_index", match_to_muon(kaon, *muonHandle_));
  d0Cand.addUserInt( "pion_mu_index", match_to_muon(pion, *muonHandle_));

  auto d0VtxFit = fillInfo(d0Cand, iEvent, kaon, pion, "kaon", "pion");
	  
  if ( not d0VtxFit.valid()  or 
       d0VtxFit.vtxProb() < minVtxProb_ or
       d0VtxFit.mass() < minD0Mass_ or  
       d0VtxFit.mass() > maxD0Mass_ or
       d0VtxFit.lxy > maxLxy_ or 
       d0VtxFit.sigLxy < minSigLxy_ ) return pat::CompositeCandidate();

  return d0Cand;
}

pat::CompositeCandidate
BmmV0Producer::getPhiToKK(const edm::Event& iEvent,
			  const pat::PackedCandidate& ipfCand1,
			  const pat::PackedCandidate& ipfCand2)
{
  pat::CompositeCandidate phiCand;
  pat::PackedCandidate pfCand1(ipfCand1);
  pfCand1.setMass(kaon_mass_);
  pat::PackedCandidate pfCand2(ipfCand2);
  pfCand2.setMass(kaon_mass_);
  phiCand.addDaughter( pfCand1 , "trk1" );
  phiCand.addDaughter( pfCand2 , "trk2" );
  AddFourMomenta addP4;
  addP4.set( phiCand );

  if ( phiCand.mass() < minPhiPreselectMass_ or phiCand.mass() > maxPhiPreselectMass_ )
    return pat::CompositeCandidate();

  phiCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
  phiCand.addUserFloat( "trk1_eta", pfCand1.eta() );
  phiCand.addUserFloat( "trk1_phi", pfCand1.phi() );
  phiCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
  phiCand.addUserFloat( "trk2_eta", pfCand2.eta() );
  phiCand.addUserFloat( "trk2_phi", pfCand2.phi() );
  phiCand.addUserFloat( "trk1_sip", trackImpactParameterSignificance(pfCand1) );
  phiCand.addUserFloat( "trk2_sip", trackImpactParameterSignificance(pfCand2) );
  phiCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  phiCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto ksVtxFit = fillInfo(phiCand, iEvent, pfCand1, pfCand2);
	  
  if ( not ksVtxFit.valid()  or 
       ksVtxFit.vtxProb() < minVtxProb_ or
       ksVtxFit.mass() < minPhiMass_ or  
       ksVtxFit.mass() > maxPhiMass_ ) return pat::CompositeCandidate();

  return phiCand;
}

pat::CompositeCandidate
BmmV0Producer::getLambdaToPPi(const edm::Event& iEvent,
			   const pat::PackedCandidate& iproton,
			   const pat::PackedCandidate& ipion)
{
  pat::CompositeCandidate lambdaCand;
  pat::PackedCandidate proton( iproton );
  proton.setMass( proton_mass_ );
  pat::PackedCandidate pion( ipion );
  pion.setMass( pion_mass_ );
  lambdaCand.addDaughter( proton, "proton" );
  lambdaCand.addDaughter( pion, "pion" );
  AddFourMomenta addP4;
  addP4.set( lambdaCand );

  if ( lambdaCand.mass() < minLambdaPreselectMass_ or lambdaCand.mass() > maxLambdaPreselectMass_ )
    return pat::CompositeCandidate();

  lambdaCand.addUserFloat( "proton_pt",  proton.pt() );
  lambdaCand.addUserFloat( "proton_eta", proton.eta() );
  lambdaCand.addUserFloat( "proton_phi", proton.phi() );
  lambdaCand.addUserFloat( "pion_pt",  pion.pt() );
  lambdaCand.addUserFloat( "pion_eta", pion.eta() );
  lambdaCand.addUserFloat( "pion_phi", pion.phi() );
  lambdaCand.addUserFloat( "proton_sip", trackImpactParameterSignificance(proton) );
  lambdaCand.addUserFloat( "pion_sip", trackImpactParameterSignificance(pion) );
  lambdaCand.addUserInt( "proton_mu_index", match_to_muon(proton,*muonHandle_));
  lambdaCand.addUserInt( "pion_mu_index", match_to_muon(pion,*muonHandle_));

  auto ksVtxFit = fillInfo(lambdaCand, iEvent, proton, pion, "proton", "pion");
	  
  if ( not ksVtxFit.valid()  or 
       ksVtxFit.vtxProb() < minVtxProb_ or
       ksVtxFit.mass() < minLambdaMass_ or  
       ksVtxFit.mass() > maxLambdaMass_ or
       ksVtxFit.lxy > maxLxy_ or 
       ksVtxFit.sigLxy < minSigLxy_
       ) return pat::CompositeCandidate();

  return lambdaCand;
}

KinematicFitResult 
BmmV0Producer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
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


pair<double,double> BmmV0Producer::computeDCA(const pat::PackedCandidate &kaon,
                                                 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle_));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}
namespace{

  const pat::PackedGenParticle*
    gen_match(const std::vector<pat::PackedGenParticle>& packedGenParticles,
	      const LorentzVector& reco)
  {
    const pat::PackedGenParticle* best_match(nullptr);
    float min_delta = 1e9;
    for (auto const & gen: packedGenParticles){
      if (deltaR(reco,gen) > 0.15) continue;
      float delta = fabs(reco.pt()-gen.pt())/gen.pt();
      if (delta > min_delta) continue;
      min_delta = delta;
      best_match = &gen;
    }
    return best_match;
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

GenMatchInfo BmmV0Producer::getGenMatchInfo( const pat::PackedCandidate& track1,
					     const pat::PackedCandidate& track2 )
{
  GenMatchInfo result;
  std::vector<const reco::Candidate*> daughters;
  result.mc_trk1 = gen_match(*packedGenParticles_, track1.p4());
  if (result.mc_trk1)
    daughters.push_back(result.mc_trk1);
  result.mc_trk2 = gen_match(*packedGenParticles_, track2.p4());
  if (result.mc_trk2)
    daughters.push_back(result.mc_trk2);
  if (daughters.size()==2){
    const auto* mother = find_common_ancestor(daughters);
    if (mother) result.match        = mother;
  }
  return result;
}

float BmmV0Producer::distanceOfClosestApproach( const reco::Track* track1,
					     const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
BmmV0Producer::distanceOfClosestApproach( const reco::Track* track,
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
BmmV0Producer::distanceOfClosestApproach( const reco::Track* track,
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


DEFINE_FWK_MODULE(BmmV0Producer);

//  LocalWords:  vertices
