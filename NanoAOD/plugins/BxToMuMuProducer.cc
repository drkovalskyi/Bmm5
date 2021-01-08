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

#include "Bmm5/NanoAOD/interface/XGBooster.h"

// 
// BxToMuMuProducer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;

namespace {
  const float MuonMass_    = 0.10565837;
  const float MuonMassErr_ = 3.5*1e-9;
  const float KaonMass_    = 0.493677;
  const float KaonMassErr_ = 1.6e-5;
  const float pionMass_    = 0.139570;
  const float pionMassErr_ = 3.5e-7;
  const float JPsiMass_    = 3.0969;
  const float JPsiMassErr_ = 92.9e-6;
};


struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, cosAlphaXY;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), lxyErr(-1.0), sigLxy(-1.0), cosAlphaXY(-999.)
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
    cosAlphaXY = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());

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

struct KalmanVertexFitResult{
  float vtxProb;
  bool  valid;
  std::vector<LorentzVector> refitVectors;
  GlobalPoint position;
  GlobalError err;
  float lxy, lxyErr, sigLxy;

  KalmanVertexFitResult():vtxProb(-1.0),valid(false),lxy(-1.0),lxyErr(-1.0),sigLxy(-1.0){}

  float mass() const
  {
    if (not valid) return -1.0;
    LorentzVector p4;
    for (auto v: refitVectors)
      p4 += v;
    return p4.mass();
  }
  
  void postprocess(const reco::BeamSpot& bs)
  {
    if (not valid) return;
    // position of the beam spot at a given z value (it takes into account the dxdz and dydz slopes)
    reco::BeamSpot::Point bs_at_z(bs.position(position.z()));
    GlobalPoint xy_displacement(position.x() - bs_at_z.x(),
				position.y() - bs_at_z.y(),
				0);
    lxy = xy_displacement.perp();
    lxyErr = sqrt(err.rerr(xy_displacement));
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
  }
};

struct DisplacementInformationIn3D{
  double decayLength, decayLengthErr, decayLength2, decayLength2Err, 
    distaceOfClosestApproach, distaceOfClosestApproachErr, distaceOfClosestApproachSig,
    distaceOfClosestApproach2, distaceOfClosestApproach2Err, distaceOfClosestApproach2Sig,
    longitudinalImpactParameter, longitudinalImpactParameterErr, longitudinalImpactParameterSig,
    longitudinalImpactParameter2, longitudinalImpactParameter2Err,longitudinalImpactParameter2Sig,
    cosAlpha, cosAlphaXY, decayTime, decayTimeError, decayTimeXY, decayTimeXYError;
  const reco::Vertex *pv,*pv2;
  int pvIndex,pv2Index;
  DisplacementInformationIn3D():decayLength(-1.0), decayLengthErr(0.), decayLength2(-1.0), decayLength2Err(0.),
				distaceOfClosestApproach(-1.0), distaceOfClosestApproachErr(0.0), distaceOfClosestApproachSig(0.0),
				distaceOfClosestApproach2(-1.0), distaceOfClosestApproach2Err(0.0), distaceOfClosestApproach2Sig(0.0),
				longitudinalImpactParameter(0.0), longitudinalImpactParameterErr(0.), longitudinalImpactParameterSig(0.),
				longitudinalImpactParameter2(0.0), longitudinalImpactParameter2Err(0.), longitudinalImpactParameter2Sig(0.),
				cosAlpha(-999.), cosAlphaXY(-999.), decayTime(-999.), decayTimeError(-999.),
				decayTimeXY(-999.), decayTimeXYError(-999.),
				pv(0), pv2(0),
				pvIndex(-1), pv2Index(-1)
  {};
};

LorentzVector makeLorentzVectorFromPxPyPzM(double px, double py, double pz, double m){
  double p2 = px*px+py*py+pz*pz;
  return LorentzVector(px,py,pz,sqrt(p2+m*m));
}

struct GenMatchInfo{
  int mu1_pdgId, mu1_motherPdgId, mu2_pdgId, mu2_motherPdgId, kaon1_pdgId, kaon1_motherPdgId,
    kaon2_pdgId, kaon2_motherPdgId, mm_pdgId, mm_motherPdgId, kmm_pdgId, kkmm_pdgId;
  float mu1_pt, mu2_pt, kaon1_pt, kaon2_pt, mm_mass, mm_pt, kmm_mass, kkmm_mass, kmm_pt, kkmm_pt;
  math::XYZPoint mm_prod_vtx, mm_vtx, kmm_prod_vtx, kkmm_prod_vtx;
  const reco::Candidate* mc_mu1;
  const reco::Candidate* mc_mu2;
  const reco::Candidate* mc_kaon1;
  const reco::Candidate* mc_kaon2;
  const reco::Candidate* match;
  const reco::Candidate* common_mother;
  GenMatchInfo():mu1_pdgId(0), mu1_motherPdgId(0), mu2_pdgId(0), mu2_motherPdgId(0), 
		 kaon1_pdgId(0), kaon1_motherPdgId(0), kaon2_pdgId(0), kaon2_motherPdgId(0),
		 mm_pdgId(0), mm_motherPdgId(0), 
		 kmm_pdgId(0), kkmm_pdgId(0), mu1_pt(0), mu2_pt(0), 
		 kaon1_pt(0), kaon2_pt(0), mm_mass(0), mm_pt(0), 
		 kmm_mass(0), kkmm_mass(0), kmm_pt(0), kkmm_pt(0),
		 mc_mu1(0), mc_mu2(0), mc_kaon1(0), mc_kaon2(0),
		 match(0), common_mother(0)
  {}
  const reco::GenParticle* gen_mu1(){
    return dynamic_cast<const reco::GenParticle*>(mc_mu1);
  }
  const reco::GenParticle* gen_mu2(){
    return dynamic_cast<const reco::GenParticle*>(mc_mu2);
  }
    

};

struct GenEventInfo{};

struct CloseTrack{
  float svDoca, svDocaErr, svProb,
    pvDoca, pvDocaErr,
    impactParameterSignificanceBS;
  const pat::PackedCandidate* pfCand;
  CloseTrack(): svDoca(-1), svDocaErr(-1), svProb(-1),
		pvDoca(-1), pvDocaErr(-1),
		impactParameterSignificanceBS(-1),
		pfCand(0)
  {};
};


struct CloseTrackInfo{
  std::vector<CloseTrack> tracks;
  unsigned int nTracksByVertexProbability(double minProb = 0.1, 
					  double minIpSignificance = -1,
					  int pvIndex = -1,
					  const pat::PackedCandidate* ignoreTrack1 = 0)
  {
    unsigned int n = 0;
    for (auto track: tracks){
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (minIpSignificance>0 and track.impactParameterSignificanceBS<minIpSignificance) continue;
      if (track.svProb<minProb) continue;
      if (pvIndex >= 0 and int(track.pfCand->vertexRef().key())!=pvIndex) continue;
      n++;
    }
    return n;
  }
  unsigned int nTracksByDisplacementSignificance(double max_svDoca = 0.03, 
						 double maxSignificance = -1,
						 int pvIndex = -1,
						 const pat::PackedCandidate* ignoreTrack1 = 0)
  {
    unsigned int n = 0;
    for (auto track: tracks){
      if (track.svDoca>max_svDoca) continue;
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (maxSignificance>0 and (track.svDocaErr<=0 or 
				 track.svDoca/track.svDocaErr > maxSignificance) ) continue;
      if (pvIndex >= 0 and int(track.pfCand->vertexRef().key())!=pvIndex) continue;
      n++;
    }
    return n;
  }
  unsigned int nTracksByBetterMatch(double max_svDoca = 0.03, 
				    double maxSignificance = 2,
				    int pvIndex = -1,
				    const pat::PackedCandidate* ignoreTrack1 = 0)
  {
    unsigned int n = 0;
    for (auto track: tracks){
      if (track.svDoca>max_svDoca) continue;
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (maxSignificance>0 and (track.svDocaErr<=0 or 
				 track.svDoca/track.svDocaErr > maxSignificance) ) continue;
      if (track.svDocaErr<=0 or (track.pvDocaErr>0 and track.svDoca/track.svDocaErr > track.pvDoca/track.pvDocaErr) ) continue;
      if (pvIndex >= 0 and int(track.pfCand->vertexRef().key())!=pvIndex) continue;
      n++;
    }
    return n;
  }
  float minDoca(double max_svDoca = 0.03, 
		int pvIndex = -1,
		const pat::PackedCandidate* ignoreTrack1 = 0)
  {
    float doca = 99.;
    for (auto track: tracks){
      if (track.svDoca>max_svDoca) continue;
      if (ignoreTrack1 and track.pfCand==ignoreTrack1) continue;
      if (pvIndex >= 0 and int(track.pfCand->vertexRef().key())!=pvIndex) continue;
      if (doca>track.svDoca) doca = track.svDoca;
    }
    return doca;
  }

  void fillCandInfo(pat::CompositeCandidate& cand, int pvIndex, std::string name)
  {
    if (name!="") name += "_";
    cand.addUserInt(   name + "nTrks",       nTracksByVertexProbability(0.1,-1.0,pvIndex) );
    cand.addUserInt(   name + "nBMTrks",     nTracksByBetterMatch() );
    cand.addUserInt(   name + "nDisTrks",    nTracksByVertexProbability(0.1, 2.0,pvIndex) );
    cand.addUserInt(   name + "closetrk",    nTracksByDisplacementSignificance(0.03, -1, pvIndex) );
    cand.addUserInt(   name + "closetrks1",  nTracksByDisplacementSignificance(0.03, 1, pvIndex) );
    cand.addUserInt(   name + "closetrks2",  nTracksByDisplacementSignificance(0.03, 2, pvIndex) );
    cand.addUserInt(   name + "closetrks3",  nTracksByDisplacementSignificance(0.03, 3, pvIndex) );
    cand.addUserFloat( name + "docatrk",     minDoca(0.03, pvIndex) );
  }
};

struct BdtReaderData {
  float fls3d, alpha, pvips, iso, chi2dof, docatrk, closetrk, m1iso, m2iso, eta, m;
};

namespace {
  // Muon container to hold muons and hadrons that may decays to muons
  // Index is the position of the muon in the original muon collection
  // For hadrons Index is -1

  class MuonCand: public pat::Muon{
  public:
    MuonCand(const pat::Muon& muon, int index):
      pat::Muon(muon), index_(index)
    {
    }
    MuonCand(const pat::PackedCandidate& hadron):
      pat::Muon(reco::Muon(hadron.charge(), hadron.p4())),
      index_(-1)
    {
      std::vector<reco::Track> tracks;
      assert(hadron.hasTrackDetails());
      tracks.push_back(*hadron.bestTrack());
      setInnerTrack(reco::TrackRef(&tracks,0));
      embedTrack();
      setPdgId(hadron.pdgId());
    }
    int index() const { return index_; }
  private:
    int index_;
  };
}


using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class BxToMuMuProducer : public edm::EDProducer {
    
public:
    
  explicit BxToMuMuProducer(const edm::ParameterSet &iConfig);
    
  ~BxToMuMuProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isGoodMuon(const pat::Muon& muon);
    
  KalmanVertexFitResult 
  vertexWithKalmanFitter( std::vector<const reco::Track*> trks, 
			 std::vector<float> masses);

  KalmanVertexFitResult 
  vertexMuonsWithKalmanFitter(const pat::Muon& muon1,
			      const pat::Muon& muon2);

  KinematicFitResult 
  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  KinematicFitResult 
  vertexMuonsWithKinematicFitter(const pat::Muon& muon1,
				 const pat::Muon& muon2);

  KinematicFitResult 
  vertexWithKinematicFitter(const pat::Muon& muon1,
			    const pat::Muon& muon2,
			    const pat::PackedCandidate& pfCand);

  /// BToKJPsiMuMuFitResult
  KinematicFitResult
  fitBToKJPsiMuMu( RefCountedKinematicParticle jpsi,
		   const pat::PackedCandidate& kaon,
		   bool applyJpsiMassConstraint);

  KinematicFitResult
  fitBToKJPsiMuMuNew( RefCountedKinematicTree jpsi,
		      const pat::PackedCandidate& kaon,
		      bool applyJpsiMassConstraint);
  KinematicFitResult
  fitBToKKMuMu( RefCountedKinematicTree jpsi,
		const pat::PackedCandidate& kaon1,
		const pat::PackedCandidate& kaon2,
		bool applyJpsiMassConstraint);

  KinematicFitResult
  vertexMuonsWithPointingConstraint( const pat::Muon& muon1,
				     const pat::Muon& muon2,
				     const reco::Vertex& primaryVertex);

  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
   				 reco::BeamSpot beamSpot);
  GenMatchInfo getGenMatchInfo( const pat::Muon& muon1,
				const pat::Muon& muon2,
				const pat::PackedCandidate* kaon1 = 0,
				const pat::PackedCandidate* kaon2 = 0 );
  const reco::Candidate* getGenParticle(const reco::Candidate* cand);

  // Two track DOCA
  float 
  distanceOfClosestApproach( const reco::Track* track1,
			     const reco::Track* track2 );
  float 
  distanceOfClosestApproach( const reco::GenParticle* track1,
			     const reco::GenParticle* track2);
  // Track to vertex DOCA
  Measurement1D
  distanceOfClosestApproach( const reco::Track* track,
			     RefCountedKinematicVertex vertex);
  Measurement1D 
  distanceOfClosestApproach( const reco::Track* track,
			     const reco::Vertex& vertex);

  DisplacementInformationIn3D 
  compute3dDisplacement(const KinematicFitResult& fit,
			const reco::VertexCollection& vertices,
			bool closestIn3D = true);

  CloseTrackInfo 
  findTracksCompatibleWithTheVertex(const pat::Muon& muon1,
				    const pat::Muon& muon2,
				    const KinematicFitResult& fit,
				    double maxDoca=0.03,
				    std::vector<const pat::PackedCandidate*> ignoreTracks = 
				    std::vector<const pat::PackedCandidate*>());
  float
  computeTrkMuonIsolation(const pat::Muon& muon, 
			  const pat::Muon& the_other_muon,
			  unsigned int primaryVertexIndex,
			  float minPt=0.5, float dR=0.5,
			  std::vector<const pat::PackedCandidate*> ignoreTracks = 
			  std::vector<const pat::PackedCandidate*>());
  float
  computeTrkMuMuIsolation(const pat::Muon& muon1, 
			  const pat::Muon& muon2,
			  unsigned int primaryVertexIndex,
			  float minPt=0.9, float dR=0.7,
			  std::vector<const pat::PackedCandidate*> ignoreTracks = 
			  std::vector<const pat::PackedCandidate*>());

  float
  otherVertexMaxProb(const pat::Muon& muon1, 
		     const pat::Muon& muon2,
		     float min_pt = 0.5,
		     float max_doca = 0.1,
		     std::vector<const pat::PackedCandidate*> ignoreTracks = 
		     std::vector<const pat::PackedCandidate*>());

  void 
  fillBtoJpsiKInfo(pat::CompositeCandidate& bCand,
		   const edm::Event& iEvent,
		   const KinematicFitResult& kinematicMuMuVertexFit,
		   const pat::Muon& muon1,
		   const pat::Muon& muon2,
		   const pat::PackedCandidate & kaon); 
  void
  fillBstoJpsiKKInfo(pat::CompositeCandidate& bCand,
		     const edm::Event& iEvent,
		     const KinematicFitResult& kinematicMuMuVertexFit,
		     const pat::Muon& muon1,
		     const pat::Muon& muon2,
		     const pat::PackedCandidate & kaon1,
		     const pat::PackedCandidate & kaon2
		    ); 
  
  void 
  fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(pat::CompositeCandidate& btokmmCand,
				    const pat::CompositeCandidate& dimuonCand,
				    const edm::Event& iEvent,
				    const KinematicFitResult& kinematicMuMuVertexFit,
				    const pat::Muon& muon1,
				    const pat::Muon& muon2,
				    const pat::PackedCandidate & kaon);
 
  KinematicFitResult 
  fillMuMuInfo(pat::CompositeCandidate& dimuonCand,
	       const edm::Event& iEvent,
	       const pat::Muon& muon1,
	       const pat::Muon& muon2); 

  void 
  injectHadronsThatMayFakeMuons(std::vector<MuonCand>& good_muon_candidates);

  float  computeAnalysisBDT(unsigned int event_idx);
  
  void setupTmvaReader(TMVA::Reader& reader, std::string file);


  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >   prunedGenToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<reco::GenParticle>* prunedGenParticles_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;

  bool isMC_;

  double ptMinMu_;
  double etaMaxMu_;
  double ptMinKaon_;

  double etaMaxKaon_;
  double DCASigMinKaon_;
  bool   diMuonCharge_;
  double minBKmmMass_;
  double maxBKmmMass_;
  double minBKKmmMass_;
  double maxBKKmmMass_;
  double maxTwoTrackDOCA_;
  BdtReaderData bdtData_;
  TMVA::Reader  bdtReader0_;
  TMVA::Reader  bdtReader1_;
  TMVA::Reader  bdtReader2_;
  std::vector<XGBooster> xgBoosters_;
};

BxToMuMuProducer::BxToMuMuProducer(const edm::ParameterSet &iConfig):
beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
beamSpot_(nullptr),
vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
pfCandToken_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
prunedGenToken_( consumes<std::vector<reco::GenParticle>> ( iConfig.getParameter<edm::InputTag>( "prunedGenParticleCollection" ) ) ),
packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( edm::InputTag( "packedGenParticles" ) ) ),
prunedGenParticles_(nullptr),
packedGenParticles_(nullptr),
impactPointExtrapolator_(0),
isMC_(            iConfig.getParameter<bool>( "isMC" ) ),
ptMinMu_(         iConfig.getParameter<double>( "MuonMinPt" ) ),
etaMaxMu_(        iConfig.getParameter<double>( "MuonMaxEta" ) ),
ptMinKaon_(       iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_(      iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_(   iConfig.getParameter<double>( "KaonMinDCASig" ) ),
diMuonCharge_(    iConfig.getParameter<bool>(   "DiMuonChargeCheck" ) ),
minBKmmMass_(     iConfig.getParameter<double>( "minBKmmMass" ) ),
maxBKmmMass_(     iConfig.getParameter<double>( "maxBKmmMass" ) ),
minBKKmmMass_(    iConfig.getParameter<double>( "minBKKmmMass" ) ),
maxBKKmmMass_(    iConfig.getParameter<double>( "maxBKKmmMass" ) ),
maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
bdtReader0_("!Color:Silent"),
bdtReader1_("!Color:Silent"),
bdtReader2_("!Color:Silent")
{
    produces<pat::CompositeCandidateCollection>("DiMuon");
    produces<pat::CompositeCandidateCollection>("BToKmumu");
    produces<pat::CompositeCandidateCollection>("BToKKmumu");
    setupTmvaReader(bdtReader0_,(iConfig.getParameter<edm::FileInPath>("bdtEvent0")).fullPath());
    setupTmvaReader(bdtReader1_,(iConfig.getParameter<edm::FileInPath>("bdtEvent1")).fullPath());
    setupTmvaReader(bdtReader2_,(iConfig.getParameter<edm::FileInPath>("bdtEvent2")).fullPath());

    xgBoosters_.push_back(XGBooster(iConfig.getParameter<edm::FileInPath>("xgbEvent0").fullPath()));
    xgBoosters_.push_back(XGBooster(iConfig.getParameter<edm::FileInPath>("xgbEvent1").fullPath()));
    xgBoosters_.push_back(XGBooster(iConfig.getParameter<edm::FileInPath>("xgbEvent2").fullPath()));

    // XGBooster
    std::vector<std::string> features = {"mm_kin_alpha", "mm_kin_alphaXY", "mm_kin_spvip", "mm_kin_pvip", 
					 "mm_iso", "mm_m1iso", "mm_m2iso", "mm_kin_sl3d", "mm_kin_vtx_chi2dof", 
					 "mm_nBMTrks", "mm_otherVtxMaxProb1", "mm_otherVtxMaxProb2"};
    for (auto& xgBooster: xgBoosters_)
      for (const auto& feature: features)
	xgBooster.addFeature(feature);
}

bool BxToMuMuProducer::isGoodMuon(const pat::Muon& muon){
  if ( not muon.isLooseMuon() ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < ptMinMu_ || fabs(muon.eta()) > etaMaxMu_ ) return false;
  return true;
}

namespace {
  void addFitInfo(pat::CompositeCandidate& cand, const KinematicFitResult& fit, std::string name, 
		  const DisplacementInformationIn3D& displacement3d = DisplacementInformationIn3D(),
		  int firstMuonDaughterIndex = -1, int secondMuonDaughterIndex = -1,
		  int firstKaonDaughterIndex = -1, int secondKaonDaughterIndex = -1 ){
    cand.addUserInt(   name+"_valid",       fit.valid() );
    cand.addUserFloat( name+"_vtx_prob",    fit.vtxProb() );
    cand.addUserFloat( name+"_vtx_chi2dof", fit.chi2()>0?fit.chi2()/fit.ndof():-1);
    cand.addUserFloat( name+"_mass",        fit.mass() );
    cand.addUserFloat( name+"_massErr",     fit.massErr() );
    cand.addUserFloat( name+"_lxy",         fit.lxy );
    cand.addUserFloat( name+"_sigLxy",      fit.sigLxy );
    cand.addUserFloat( name+"_cosAlphaXY",  fit.cosAlphaXY );
    cand.addUserFloat( name+"_alpha",       fabs(displacement3d.cosAlpha)<=1?acos(displacement3d.cosAlpha):-999. );
    cand.addUserFloat( name+"_alphaXY",     fabs(fit.cosAlphaXY)<=1?acos(fit.cosAlphaXY):-999. );
    cand.addUserFloat( name+"_vtx_x",       fit.valid()?fit.refitVertex->position().x():0 );
    cand.addUserFloat( name+"_vtx_xErr",    fit.valid()?sqrt(fit.refitVertex->error().cxx()):0 );
    cand.addUserFloat( name+"_vtx_y",       fit.valid()?fit.refitVertex->position().y():0 );
    cand.addUserFloat( name+"_vtx_yErr",    fit.valid()?sqrt(fit.refitVertex->error().cyy()):0 );
    cand.addUserFloat( name+"_vtx_z",       fit.valid()?fit.refitVertex->position().z():0 );
    cand.addUserFloat( name+"_vtx_zErr",    fit.valid()?sqrt(fit.refitVertex->error().czz()):0 );
    cand.addUserFloat( name+"_pt",          fit.p3().perp() );
    cand.addUserFloat( name+"_eta",         fit.p3().eta() );
    cand.addUserFloat( name+"_phi",         fit.p3().phi() );
    
    // IP info
    cand.addUserFloat( name+"_l3d",         displacement3d.decayLength);
    cand.addUserFloat( name+"_sl3d",        displacement3d.decayLengthErr>0?displacement3d.decayLength/displacement3d.decayLengthErr:0);
    cand.addUserFloat( name+"_pv_z",        displacement3d.pv?displacement3d.pv->position().z():0);
    cand.addUserFloat( name+"_pv_zErr",     displacement3d.pv?displacement3d.pv->zError():0);
    cand.addUserFloat( name+"_pvip",        displacement3d.distaceOfClosestApproach);
    cand.addUserFloat( name+"_pvipSig",     displacement3d.distaceOfClosestApproachSig);
    cand.addUserFloat( name+"_pvipErr",     displacement3d.distaceOfClosestApproachErr);
    cand.addUserFloat( name+"_pv2ip",       displacement3d.distaceOfClosestApproach2);
    cand.addUserFloat( name+"_pv2ipSig",    displacement3d.distaceOfClosestApproach2Sig);
    cand.addUserFloat( name+"_pv2ipErr",    displacement3d.distaceOfClosestApproach2Err);
    cand.addUserFloat( name+"_pvlip",       displacement3d.longitudinalImpactParameter);
    cand.addUserFloat( name+"_pvlipSig",    displacement3d.longitudinalImpactParameterSig);
    cand.addUserFloat( name+"_pvlipErr",    displacement3d.longitudinalImpactParameterErr);
    cand.addUserFloat( name+"_pv2lip",      displacement3d.longitudinalImpactParameter2);
    cand.addUserFloat( name+"_pv2lipSig",   displacement3d.longitudinalImpactParameter2Sig);
    cand.addUserFloat( name+"_pv2lipErr",   displacement3d.longitudinalImpactParameter2Err);
    cand.addUserInt(   name+"_pvIndex",     displacement3d.pvIndex);

    // DecayTime
    cand.addUserFloat( name+"_tau",         displacement3d.decayTime);
    cand.addUserFloat( name+"_taue",        displacement3d.decayTimeError);
    cand.addUserFloat( name+"_tauxy",       displacement3d.decayTimeXY);
    cand.addUserFloat( name+"_tauxye",      displacement3d.decayTimeXYError);

    // Refitted daughter information
    if (firstMuonDaughterIndex>=0){
      cand.addUserFloat( name+"_mu1pt",       fit.dau_p3(firstMuonDaughterIndex).perp() );
      cand.addUserFloat( name+"_mu1eta",      fit.dau_p3(firstMuonDaughterIndex).eta() );
      cand.addUserFloat( name+"_mu1phi",      fit.dau_p3(firstMuonDaughterIndex).phi() );
    }
    if (secondMuonDaughterIndex>=0){
      cand.addUserFloat( name+"_mu2pt",       fit.dau_p3(secondMuonDaughterIndex).perp() );
      cand.addUserFloat( name+"_mu2eta",      fit.dau_p3(secondMuonDaughterIndex).eta() );
      cand.addUserFloat( name+"_mu2phi",      fit.dau_p3(secondMuonDaughterIndex).phi() );
    }
    if (firstKaonDaughterIndex>=0){
      cand.addUserFloat( name+"_kaon1pt",     fit.dau_p3(firstKaonDaughterIndex).perp() );
      cand.addUserFloat( name+"_kaon1eta",    fit.dau_p3(firstKaonDaughterIndex).eta() );
      cand.addUserFloat( name+"_kaon1phi",    fit.dau_p3(firstKaonDaughterIndex).phi() );
    }
    if (secondKaonDaughterIndex>=0){
      cand.addUserFloat( name+"_kaon2pt",     fit.dau_p3(secondKaonDaughterIndex).perp() );
      cand.addUserFloat( name+"_kaon2eta",    fit.dau_p3(secondKaonDaughterIndex).eta() );
      cand.addUserFloat( name+"_kaon2phi",    fit.dau_p3(secondKaonDaughterIndex).phi() );
    }
  }
}

CloseTrackInfo 
BxToMuMuProducer::findTracksCompatibleWithTheVertex(const pat::Muon& muon1,
						    const pat::Muon& muon2,
						    const KinematicFitResult& fit, 
						    double maxDoca,
						    std::vector<const pat::PackedCandidate*> ignoreTracks)
{
  CloseTrackInfo result;
  if (not fit.valid()) return result;
  for (const auto& pfCand: *pfCandHandle_.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (deltaR(*trk, pfCand) < 0.01){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    
    if (deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    double mu1_kaon_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						     pfCand.bestTrack());
    double mu2_kaon_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
						     pfCand.bestTrack());
    if (mu1_kaon_doca>maxDoca or mu2_kaon_doca>maxDoca) continue;
    
    CloseTrack track;
    track.pfCand = &pfCand;
    auto doca = distanceOfClosestApproach(pfCand.bestTrack(),fit.refitVertex);
    track.svDoca = doca.value();
    track.svDocaErr = doca.error();

    // add PV doca
    if (pfCand.vertexRef().key()<pvHandle_->size()){
      doca = distanceOfClosestApproach(pfCand.bestTrack(),pvHandle_->at(pfCand.vertexRef().key()) );
      track.pvDoca = doca.value();
      track.pvDocaErr = doca.error();
    }
    
    auto fit_result = vertexWithKinematicFitter(muon1, muon2, pfCand);
    if (fit_result.valid()){
      track.svProb = fit_result.vtxProb();
      track.impactParameterSignificanceBS = pfCand.bestTrack()->dxyError()>0 ? fabs(pfCand.bestTrack()->dxy(*beamSpot_))/pfCand.bestTrack()->dxyError():0.0;
    }
    result.tracks.push_back(track);
  }

  return result;
}

float
BxToMuMuProducer::computeTrkMuonIsolation(const pat::Muon& the_muon, const pat::Muon& the_other_muon, 
					  unsigned int primaryVertexIndex,
					  float minPt, float dR,
					  std::vector<const pat::PackedCandidate*> ignoreTracks)
{
  float sumPt(0);
  for (const auto& pfCand: *pfCandHandle_.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(the_muon, pfCand) < 0.01 || deltaR(the_other_muon, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    if (pfCand.pt()<minPt) continue;
    if (pfCand.vertexRef().key()!=primaryVertexIndex) continue;
    if (deltaR(the_muon, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }

  return the_muon.pt()/(the_muon.pt()+sumPt);
}

float
BxToMuMuProducer::computeTrkMuMuIsolation(const pat::Muon& muon1, const pat::Muon& muon2, 
					  unsigned int primaryVertexIndex,
					  float minPt, float dR,
					  std::vector<const pat::PackedCandidate*> ignoreTracks)
{
  float sumPt(0);
  auto b_p4 = muon1.p4()+muon2.p4();
  for (const auto& pfCand: *pfCandHandle_.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    if (pfCand.pt()<minPt) continue;
    if (pfCand.vertexRef().key()!=primaryVertexIndex) continue;
    if (deltaR(b_p4, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }

  return b_p4.pt()/(b_p4.pt()+sumPt);
}


float
BxToMuMuProducer::otherVertexMaxProb(const pat::Muon& muon1, 
				     const pat::Muon& muon2,
				     float minPt,
				     float max_doca,
				     std::vector<const pat::PackedCandidate*> ignoreTracks){
  float bestMu1Vtx = 0;
  float bestMu2Vtx = 0;
  KalmanVertexFitter kvf;
  std::vector<reco::TransientTrack> transTrksForMu1Vertex;
  transTrksForMu1Vertex.push_back((*theTTBuilder_).build(muon1.innerTrack().get()));
  std::vector<reco::TransientTrack> transTrksForMu2Vertex;
  transTrksForMu2Vertex.push_back((*theTTBuilder_).build(muon2.innerTrack().get()));


  for (const auto& pfCand: *pfCandHandle_.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    if (pfCand.pt()<minPt) continue;
    double mu1_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						pfCand.bestTrack());
    double mu2_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
						pfCand.bestTrack());
    if (mu1_doca < max_doca and mu1_doca < mu2_doca){
      // first  muon is closer - check vertex probability
      transTrksForMu1Vertex.push_back((*theTTBuilder_).build(pfCand.bestTrack()));
      TransientVertex tv = kvf.vertex(transTrksForMu1Vertex);
      if ( tv.isValid() ){
	float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
	if (vtxProb > bestMu1Vtx) bestMu1Vtx = vtxProb;
      }
      transTrksForMu1Vertex.pop_back();
    }
    if (mu2_doca < max_doca and mu2_doca < mu1_doca){
      // second  muon is closer - check vertex probability
      transTrksForMu2Vertex.push_back((*theTTBuilder_).build(pfCand.bestTrack()));
      TransientVertex tv = kvf.vertex(transTrksForMu2Vertex);
      if ( tv.isValid() ){
	float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
	if (vtxProb > bestMu2Vtx) bestMu2Vtx = vtxProb;
      }
      transTrksForMu2Vertex.pop_back();
    }
  }
  return max(bestMu1Vtx,bestMu2Vtx);
}

namespace {
  math::XYZPoint getProductionVertex( const reco::Candidate* cand){
    if (not cand) return math::XYZPoint();
    const reco::Candidate* primary = cand;
    // handle oscillation and radiation
    while (primary->mother() and abs(primary->pdgId())==abs(primary->mother()->pdgId()))
      primary = primary->mother();
    return primary->vertex();
  }

  double computeDecayTime( const GenMatchInfo& info ){
    if (not info.match) return -1.0;
    auto prod_vtx = getProductionVertex(info.match);
    if (prod_vtx.r()<1e-12) return -2.0;
    return (prod_vtx-info.mm_vtx).r()/TMath::Ccgs()*info.match->mass()/info.match->p();
  }
}

KinematicFitResult 
BxToMuMuProducer::fillMuMuInfo(pat::CompositeCandidate& dimuonCand,
			       const edm::Event& iEvent,
			       const pat::Muon& muon1,
			       const pat::Muon& muon2
			       ) 
{
  auto kinematicMuMuVertexFit = vertexMuonsWithKinematicFitter(muon1, muon2);
  kinematicMuMuVertexFit.postprocess(*beamSpot_);

  // printf("kinematicMuMuVertexFit (x,y,z): (%7.3f,%7.3f,%7.3f)\n", 
  // 	 kinematicMuMuVertexFit.refitVertex->position().x(),
  // 	 kinematicMuMuVertexFit.refitVertex->position().y(),
  // 	 kinematicMuMuVertexFit.refitVertex->position().z());
  auto displacement3D = compute3dDisplacement(kinematicMuMuVertexFit, *pvHandle_.product(),true);
  addFitInfo(dimuonCand, kinematicMuMuVertexFit, "kin", displacement3D,0,1);
  
  if (isMC_){
    auto gen_mm = getGenMatchInfo(muon1,muon2);
    // int mu1_pdgId, mu1_motherPdgId, mu2_pdgId, mu2_motherPdgId, kaon_pdgId, kaon_motherPdgId,
    // mm_pdgId, mm_motherPdgId, kmm_pdgId;
    // float mu1_pt, mu2_pt, kaon_pt, mm_mass, mm_pt, kmm_mass, kmm_pt;
    dimuonCand.addUserInt(  "gen_mu1_pdgId",   gen_mm.mu1_pdgId);
    dimuonCand.addUserInt(  "gen_mu1_mpdgId",  gen_mm.mu1_motherPdgId);
    dimuonCand.addUserFloat("gen_mu1_pt",      gen_mm.mu1_pt);
    dimuonCand.addUserInt(  "gen_mu2_pdgId",   gen_mm.mu2_pdgId);
    dimuonCand.addUserInt(  "gen_mu2_mpdgId",  gen_mm.mu2_motherPdgId);
    dimuonCand.addUserFloat("gen_mu2_pt",      gen_mm.mu2_pt);
    dimuonCand.addUserFloat("gen_mass",        gen_mm.mm_mass);
    dimuonCand.addUserFloat("gen_pt",          gen_mm.mm_pt);
    dimuonCand.addUserInt(  "gen_pdgId",       gen_mm.mm_pdgId);
    dimuonCand.addUserInt(  "gen_mpdgId",      gen_mm.mm_motherPdgId);
    dimuonCand.addUserInt(  "gen_cpdgId",      gen_mm.common_mother?gen_mm.common_mother->pdgId():0);
    dimuonCand.addUserFloat("gen_prod_x",      gen_mm.mm_prod_vtx.x());
    dimuonCand.addUserFloat("gen_prod_y",      gen_mm.mm_prod_vtx.y());
    dimuonCand.addUserFloat("gen_prod_z",      gen_mm.mm_prod_vtx.z());
    dimuonCand.addUserFloat("gen_vtx_x",       gen_mm.mm_vtx.x());
    dimuonCand.addUserFloat("gen_vtx_y",       gen_mm.mm_vtx.y());
    dimuonCand.addUserFloat("gen_vtx_z",       gen_mm.mm_vtx.z());
    dimuonCand.addUserFloat("gen_l3d",         (gen_mm.mm_prod_vtx-gen_mm.mm_vtx).r());
    dimuonCand.addUserFloat("gen_lxy",         (gen_mm.mm_prod_vtx-gen_mm.mm_vtx).rho());
    dimuonCand.addUserFloat("gen_tau",         computeDecayTime(gen_mm));
    
    double mm_doca = -1;
    if (gen_mm.gen_mu1() and gen_mm.gen_mu2())
      mm_doca = distanceOfClosestApproach(gen_mm.gen_mu1(), gen_mm.gen_mu2());
    dimuonCand.addUserFloat("gen_doca",        mm_doca);
    
    
  }

  int pvIndex = displacement3D.pvIndex;

  // Look for additional tracks compatible with the dimuon vertex
  auto closeTracks = findTracksCompatibleWithTheVertex(muon1,muon2,kinematicMuMuVertexFit);
  closeTracks.fillCandInfo(dimuonCand, pvIndex, "");

  dimuonCand.addUserFloat( "m1iso",     computeTrkMuonIsolation(muon1,muon2,pvIndex,0.5,0.5));
  dimuonCand.addUserFloat( "m2iso",     computeTrkMuonIsolation(muon2,muon1,pvIndex,0.5,0.5));
  dimuonCand.addUserFloat( "iso",       computeTrkMuMuIsolation(muon2,muon1,pvIndex,0.9,0.7));
  dimuonCand.addUserFloat( "otherVtxMaxProb", otherVertexMaxProb(muon1,muon2,0.5));
  dimuonCand.addUserFloat( "otherVtxMaxProb1", otherVertexMaxProb(muon1,muon2,1.0));
  dimuonCand.addUserFloat( "otherVtxMaxProb2", otherVertexMaxProb(muon1,muon2,2.0));

  // BDT
  bdtData_.fls3d    = dimuonCand.userFloat("kin_sl3d");
  bdtData_.alpha    = dimuonCand.userFloat("kin_alpha");
  bdtData_.pvips    = dimuonCand.userFloat("kin_pvipErr")>0?dimuonCand.userFloat("kin_pvip")/dimuonCand.userFloat("kin_pvipErr"):999.;
  bdtData_.iso      = dimuonCand.userFloat("iso");
  bdtData_.chi2dof  = dimuonCand.userFloat("kin_vtx_chi2dof");
  bdtData_.docatrk  = dimuonCand.userFloat("docatrk");
  bdtData_.closetrk = dimuonCand.userInt(  "closetrk");
  bdtData_.m1iso    = dimuonCand.userFloat("m1iso");
  bdtData_.m2iso    = dimuonCand.userFloat("m2iso");
  bdtData_.eta      = dimuonCand.userFloat("kin_eta");	  
  bdtData_.m        = dimuonCand.userFloat("kin_mass");	  

  dimuonCand.addUserFloat("bdt",computeAnalysisBDT(iEvent.eventAuxiliary().event()%3));

  // XGBoost
  unsigned int xg_index = iEvent.eventAuxiliary().event()%3;
  xgBoosters_.at(xg_index).set("mm_kin_alpha",       dimuonCand.userFloat("kin_alpha"));
  xgBoosters_.at(xg_index).set("mm_kin_alphaXY",     dimuonCand.userFloat("kin_cosAlphaXY"));
  xgBoosters_.at(xg_index).set("mm_kin_spvip",       dimuonCand.userFloat("kin_pvipErr")>0?dimuonCand.userFloat("kin_pvip")/dimuonCand.userFloat("kin_pvipErr"):999.);
  xgBoosters_.at(xg_index).set("mm_kin_pvip",        dimuonCand.userFloat("kin_pvip"));
  xgBoosters_.at(xg_index).set("mm_iso",             dimuonCand.userFloat("iso"));
  xgBoosters_.at(xg_index).set("mm_m1iso",           dimuonCand.userFloat("m1iso"));
  xgBoosters_.at(xg_index).set("mm_m2iso",           dimuonCand.userFloat("m2iso"));
  xgBoosters_.at(xg_index).set("mm_kin_sl3d",        dimuonCand.userFloat("kin_sl3d"));
  xgBoosters_.at(xg_index).set("mm_kin_vtx_chi2dof", dimuonCand.userFloat("kin_vtx_chi2dof"));
  xgBoosters_.at(xg_index).set("mm_nBMTrks",         dimuonCand.userInt(  "nBMTrks"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb1", dimuonCand.userFloat(  "otherVtxMaxProb1"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb2", dimuonCand.userFloat(  "otherVtxMaxProb2"));
  
  dimuonCand.addUserFloat("mva", xgBoosters_.at(xg_index).predict());

  // std::cout << "\n\nmm_kin_alpha: " <<        dimuonCand.userFloat("kin_alpha") << std::endl;
  // std::cout << "mm_kin_alphaXY: " <<      dimuonCand.userFloat("kin_cosAlphaXY") << std::endl;
  // std::cout << "mm_kin_spvip: " <<        (dimuonCand.userFloat("kin_pvipErr")>0?dimuonCand.userFloat("kin_pvip")/dimuonCand.userFloat("kin_pvipErr"):999.) << std::endl;
  // std::cout << "kin_pvipErr: " <<         dimuonCand.userFloat("kin_pvipErr") << std::endl;
  // std::cout << "mm_kin_pvip: " <<         dimuonCand.userFloat("kin_pvip") << std::endl;
  // std::cout << "mm_iso: " <<              dimuonCand.userFloat("iso") << std::endl;
  // std::cout << "mm_m1iso: " <<            dimuonCand.userFloat("m1iso") << std::endl;
  // std::cout << "mm_m2iso: " <<            dimuonCand.userFloat("m2iso") << std::endl;
  // std::cout << "mm_kin_sl3d: " <<         dimuonCand.userFloat("kin_sl3d") << std::endl;
  // std::cout << "mm_nBMTrks: " <<          dimuonCand.userInt(  "nBMTrks") << std::endl;
  // std::cout << "mm_kin_vtx_chi2dof: " <<  dimuonCand.userFloat("kin_vtx_chi2dof") << std::endl;
  // std::cout << "mm_closetrks1: " <<       dimuonCand.userInt(  "closetrk") << std::endl;
  // std::cout << "mm_nDisTrks: " <<         dimuonCand.userInt(  "nDisTrks") << std::endl;
  // std::cout << "mva: " <<  xgBoosters_.at(xg_index).predict() << std::endl;

  // Refit with pointing constraint
  auto bToMuMu_PC = vertexMuonsWithPointingConstraint(muon1,muon2,*displacement3D.pv);
  addFitInfo(dimuonCand, bToMuMu_PC, "kinpc");
  
  return kinematicMuMuVertexFit;
}

void BxToMuMuProducer::fillBtoJpsiKInfo(pat::CompositeCandidate& btokmmCand,
					const edm::Event& iEvent,
					const KinematicFitResult& kinematicMuMuVertexFit,
					const pat::Muon& muon1,
					const pat::Muon& muon2,
					const pat::PackedCandidate & kaon
					) 
{
  btokmmCand.addUserFloat("kaon_pt",     kaon.pt());
  btokmmCand.addUserFloat("kaon_eta",    kaon.eta());
  btokmmCand.addUserFloat("kaon_phi",    kaon.phi());
  btokmmCand.addUserFloat("kaon_dxy_bs", kaon.bestTrack()->dxy(*beamSpot_));
  btokmmCand.addUserFloat("kaon_sdxy_bs", 
			  kaon.bestTrack()->dxyError()>0 ? fabs(kaon.bestTrack()->dxy(*beamSpot_))/kaon.bestTrack()->dxyError():0.0);
  btokmmCand.addUserInt("kaon_charge", kaon.charge());
  if (isMC_){
    auto gen_kmm = getGenMatchInfo(muon1,muon2,&kaon);
    btokmmCand.addUserInt(  "gen_kaon_pdgId",  gen_kmm.kaon1_pdgId);
    btokmmCand.addUserInt(  "gen_kaon_mpdgId", gen_kmm.kaon1_motherPdgId);
    btokmmCand.addUserFloat("gen_kaon_pt",     gen_kmm.kaon1_pt);
    btokmmCand.addUserFloat("gen_mass",        gen_kmm.kmm_mass);
    btokmmCand.addUserFloat("gen_pt",          gen_kmm.kmm_pt);
    btokmmCand.addUserInt(  "gen_pdgId",       gen_kmm.kmm_pdgId);
    btokmmCand.addUserFloat("gen_prod_x",      gen_kmm.kmm_prod_vtx.x());
    btokmmCand.addUserFloat("gen_prod_y",      gen_kmm.kmm_prod_vtx.y());
    btokmmCand.addUserFloat("gen_prod_z",      gen_kmm.kmm_prod_vtx.z());
    btokmmCand.addUserFloat("gen_l3d",         (gen_kmm.kmm_prod_vtx-gen_kmm.mm_vtx).r());
    btokmmCand.addUserFloat("gen_lxy",         (gen_kmm.kmm_prod_vtx-gen_kmm.mm_vtx).rho());
    btokmmCand.addUserFloat("gen_tau",         computeDecayTime(gen_kmm));
    btokmmCand.addUserFloat("gen_cpdgId",      gen_kmm.common_mother?gen_kmm.common_mother->pdgId():0);
  }

  // if (kaon.genParticle()){
  // 	btokmmCand.addUserInt("kaon_mc_pdgId", kaon.genParticle().pdgId());
  // } else {
  // 	btokmmCand.addUserInt("kaon_mc_pdgId", 0);
  // }
  
  auto bToKJPsiMuMu_NoMassConstraint = fitBToKJPsiMuMu(kinematicMuMuVertexFit.refitMother, kaon, false);
  bToKJPsiMuMu_NoMassConstraint.postprocess(*beamSpot_);
  auto bToKJPsiMuMu_NoMassConstraint_displacement = compute3dDisplacement(bToKJPsiMuMu_NoMassConstraint, *pvHandle_.product(),true);
  addFitInfo(btokmmCand, bToKJPsiMuMu_NoMassConstraint, "nomc", bToKJPsiMuMu_NoMassConstraint_displacement,-1,-1,1);
  
  // worse performing option
  // auto bToKJPsiMuMuWithMassConstraint = fitBToKJPsiMuMu(kinematicMuMuVertexFit.refitMother, kaon, true);
  // bToKJPsiMuMuWithMassConstraint.postprocess(beamSpot);
  // addFitInfo(btokmmCand, bToKJPsiMuMuWithMassConstraint, "jpsimc");
  
  auto bToKJPsiMuMu_MassConstraint = fitBToKJPsiMuMuNew(kinematicMuMuVertexFit.refitTree, kaon, true);
  bToKJPsiMuMu_MassConstraint.postprocess(*beamSpot_);
  auto bToKJPsiMuMu_MassConstraint_displacement = compute3dDisplacement(bToKJPsiMuMu_MassConstraint, *pvHandle_.product(),true);
  addFitInfo(btokmmCand, bToKJPsiMuMu_MassConstraint, "jpsimc", bToKJPsiMuMu_MassConstraint_displacement,-1,-1,1);
  
  // broken pointing constraint
  // auto bToKJPsiMuMu_MC_PC = refitWithPointingConstraint(bToKJPsiMuMu_MC.refitTree, primaryVertex);
  // bToKJPsiMuMu_MC_PC.postprocess(beamSpot);
  // addFitInfo(btokmmCand, bToKJPsiMuMu_MC_PC, "mcpc");
}

void BxToMuMuProducer::fillBstoJpsiKKInfo(pat::CompositeCandidate& bCand,
					  const edm::Event& iEvent,
					  const KinematicFitResult& kinematicMuMuVertexFit,
					  const pat::Muon& muon1,
					  const pat::Muon& muon2,
					  const pat::PackedCandidate & kaon1,
					  const pat::PackedCandidate & kaon2
					) 
{
  bCand.addUserFloat("kaon1_pt",     kaon1.pt());
  bCand.addUserFloat("kaon1_eta",    kaon1.eta());
  bCand.addUserFloat("kaon1_phi",    kaon1.phi());
  bCand.addUserFloat("kaon1_dxy_bs", kaon1.bestTrack()->dxy(*beamSpot_));
  bCand.addUserFloat("kaon1_sdxy_bs", 
			  kaon1.bestTrack()->dxyError()>0 ? fabs(kaon1.bestTrack()->dxy(*beamSpot_))/kaon1.bestTrack()->dxyError():0.0);
  bCand.addUserInt("kaon1_charge", kaon1.charge());

  bCand.addUserFloat("kaon2_pt",     kaon2.pt());
  bCand.addUserFloat("kaon2_eta",    kaon2.eta());
  bCand.addUserFloat("kaon2_phi",    kaon2.phi());
  bCand.addUserFloat("kaon2_dxy_bs", kaon2.bestTrack()->dxy(*beamSpot_));
  bCand.addUserFloat("kaon2_sdxy_bs", 
			  kaon2.bestTrack()->dxyError()>0 ? fabs(kaon2.bestTrack()->dxy(*beamSpot_))/kaon2.bestTrack()->dxyError():0.0);
  bCand.addUserInt("kaon2_charge", kaon2.charge());
  bCand.addUserFloat("kk_mass",    (kaon1.p4()+kaon2.p4()).mass());

  if (isMC_){
    auto gen_info = getGenMatchInfo(muon1,muon2,&kaon1,&kaon2);
    bCand.addUserInt(  "gen_kaon1_pdgId",  gen_info.kaon1_pdgId);
    bCand.addUserInt(  "gen_kaon1_mpdgId", gen_info.kaon1_motherPdgId);
    bCand.addUserFloat("gen_kaon1_pt",     gen_info.kaon1_pt);
    bCand.addUserInt(  "gen_kaon2_pdgId",  gen_info.kaon2_pdgId);
    bCand.addUserInt(  "gen_kaon2_mpdgId", gen_info.kaon2_motherPdgId);
    bCand.addUserFloat("gen_kaon2_pt",     gen_info.kaon2_pt);
    bCand.addUserFloat("gen_mass",         gen_info.kkmm_mass);
    bCand.addUserFloat("gen_pt",           gen_info.kkmm_pt);
    bCand.addUserInt(  "gen_pdgId",        gen_info.kkmm_pdgId);
    bCand.addUserFloat("gen_prod_x",       gen_info.kkmm_prod_vtx.x());
    bCand.addUserFloat("gen_prod_y",       gen_info.kkmm_prod_vtx.y());
    bCand.addUserFloat("gen_prod_z",       gen_info.kkmm_prod_vtx.z());
    bCand.addUserFloat("gen_l3d",         (gen_info.kkmm_prod_vtx-gen_info.mm_vtx).r());
    bCand.addUserFloat("gen_lxy",         (gen_info.kkmm_prod_vtx-gen_info.mm_vtx).rho());
    bCand.addUserFloat("gen_tau",          computeDecayTime(gen_info));
    bCand.addUserFloat("gen_cpdgId",       gen_info.common_mother?gen_info.common_mother->pdgId():0);
  }

  auto bToKKJPsiMuMu = fitBToKKMuMu(kinematicMuMuVertexFit.refitTree, kaon1, kaon2, true);
  bToKKJPsiMuMu.postprocess(*beamSpot_);
  auto bToKKJPsiMuMu_displacement = compute3dDisplacement(bToKKJPsiMuMu, *pvHandle_.product(),true);
  addFitInfo(bCand, bToKKJPsiMuMu, "jpsikk", bToKKJPsiMuMu_displacement,-1,-1,1,2);
  bCand.addUserFloat("jpsikk_kk_mass",    bToKKJPsiMuMu.refit_mass(1,2));
}


void 
BxToMuMuProducer::fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(pat::CompositeCandidate& mmK,
							       const pat::CompositeCandidate& mm,
							       const edm::Event& iEvent,
							       const KinematicFitResult& kinematicMuMuVertexFit,
							       const pat::Muon& muon1,
							       const pat::Muon& muon2,
							       const pat::PackedCandidate & kaon
					) 
{
  ///////
  //// Treat B->JpsiK as B->mm for signal studies in data
  //
  std::vector<const pat::PackedCandidate*> ignoreTracks;
  ignoreTracks.push_back(&kaon);

  int pvIndex = mmK.userInt("jpsimc_pvIndex");

  // Look for additional tracks compatible with the dimuon vertex
  auto closeTracks = findTracksCompatibleWithTheVertex(muon1,muon2,kinematicMuMuVertexFit,0.03,ignoreTracks);
  closeTracks.fillCandInfo(mmK, pvIndex, "bmm");

  mmK.addUserFloat( "bmm_m1iso",     computeTrkMuonIsolation(muon1,muon2,pvIndex,0.5,0.5,ignoreTracks));
  mmK.addUserFloat( "bmm_m2iso",     computeTrkMuonIsolation(muon2,muon1,pvIndex,0.5,0.5,ignoreTracks));
  mmK.addUserFloat( "bmm_iso",       computeTrkMuMuIsolation(muon2,muon1,pvIndex,0.9,0.7,ignoreTracks));
  mmK.addUserFloat( "bmm_otherVtxMaxProb",  otherVertexMaxProb(muon1,muon2,0.5,0.1,ignoreTracks));
  mmK.addUserFloat( "bmm_otherVtxMaxProb1", otherVertexMaxProb(muon1,muon2,1.0,0.1,ignoreTracks));
  mmK.addUserFloat( "bmm_otherVtxMaxProb2", otherVertexMaxProb(muon1,muon2,2.0,0.1,ignoreTracks));

  // BDT
  bdtData_.fls3d    = mm.userFloat("kin_sl3d");
  bdtData_.alpha    = mmK.userFloat("jpsimc_alpha");
  bdtData_.pvips    = mmK.userFloat("jpsimc_pvipErr")>0?mmK.userFloat("jpsimc_pvip")/mmK.userFloat("jpsimc_pvipErr"):999;
  // One can use bkmm without mass constraint, but it doesn't help
  // bdtData_.alpha    = mmK.userFloat("nomc_alpha");
  // bdtData_.pvips    = mmK.userFloat("nomc_pvip")/mmK.userFloat("nomc_pvipErr");
  bdtData_.iso      = mmK.userFloat("bmm_iso");
  bdtData_.chi2dof  = mm.userFloat("kin_vtx_chi2dof");
  bdtData_.docatrk  = mmK.userFloat("bmm_docatrk");
  bdtData_.closetrk = mmK.userInt(  "bmm_closetrk");
  bdtData_.m1iso    = mmK.userFloat("bmm_m1iso");
  bdtData_.m2iso    = mmK.userFloat("bmm_m2iso");
  bdtData_.eta      = mmK.userFloat("jpsimc_eta");	  
  bdtData_.m        = mmK.userFloat("jpsimc_mass");	  

  mmK.addUserFloat("bmm_bdt",computeAnalysisBDT(iEvent.eventAuxiliary().event()%3));

  // XGBoost
  unsigned int xg_index = iEvent.eventAuxiliary().event()%3;
  // Pointing angle - mmK
  xgBoosters_.at(xg_index).set("mm_kin_alpha",       mmK.userFloat("jpsimc_alpha"));
  xgBoosters_.at(xg_index).set("mm_kin_alphaXY",     mmK.userFloat("jpsimc_cosAlphaXY"));
  // PV matching - mmK
  xgBoosters_.at(xg_index).set("mm_kin_spvip",       mmK.userFloat("jpsimc_pvipErr")>0?
			       mmK.userFloat("jpsimc_pvip")/mmK.userFloat("jpsimc_pvipErr"):
			       999.);
  xgBoosters_.at(xg_index).set("mm_kin_pvip",        mmK.userFloat("jpsimc_pvip"));
  // Isolation and extra track variables need to be recomputed ignoring kaon
  xgBoosters_.at(xg_index).set("mm_iso",             mmK.userFloat("bmm_iso"));
  xgBoosters_.at(xg_index).set("mm_m1iso",           mmK.userFloat("bmm_m1iso"));
  xgBoosters_.at(xg_index).set("mm_m2iso",           mmK.userFloat("bmm_m2iso"));
  xgBoosters_.at(xg_index).set("mm_nBMTrks",         mmK.userInt(  "bmm_nBMTrks"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb1", mmK.userFloat(  "bmm_otherVtxMaxProb1"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb2", mmK.userFloat(  "bmm_otherVtxMaxProb2"));
  // Vertexing - mm
  xgBoosters_.at(xg_index).set("mm_kin_vtx_chi2dof", mm.userFloat("kin_vtx_chi2dof"));
  // Flight length significance - mm
  xgBoosters_.at(xg_index).set("mm_kin_sl3d",        mm.userFloat("kin_sl3d"));
  
  mmK.addUserFloat("bmm_mva", xgBoosters_.at(xg_index).predict());


}
namespace {
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* candidate){
    if (ancestor == candidate) return true;
    for (unsigned int iMother=0; iMother < candidate->numberOfMothers(); ++iMother){
      if (isAncestor(ancestor,candidate->mother(iMother)))
	return true;
    }
    return false;
  }
}

void 
BxToMuMuProducer::injectHadronsThatMayFakeMuons(std::vector<MuonCand>& good_muon_candidates){
  // Loop over gen info and find interesting events
  for (auto const & cand: *prunedGenParticles_){
    // keep only interesting b-hadrons
    if ( abs(cand.pdgId()) != 511 and  // B0
	 abs(cand.pdgId()) != 521 and  // B+/-
	 abs(cand.pdgId()) != 531 and  // Bs 
	 abs(cand.pdgId()) != 541 )  // Bc 
      continue;

    // check direct daughter infor for signs of neutral B oscilations
    bool final_b = true;
    for (auto const& daughter: cand.daughterRefVector()){
      if (daughter->pdgId() == -cand.pdgId()){
	final_b = false;
	break;
      }
    }
    if (not final_b) continue;

    long long int signature = 1;
    std::vector<const reco::Candidate*> final_state_particles;

    // Loop over packed gen particles that represent final state
    // particles at gen level and compute the decay signature ignoring
    // photons
    for (auto const& dau: *packedGenParticles_){
      auto mother = dau.mother(0);
      if (mother and isAncestor(&cand,mother)){
	if (dau.pdgId()!=22){
	  signature *= dau.pdgId();
	  if (abs(dau.pdgId()) == 211 or // pions
	      abs(dau.pdgId()) == 321 or // kaons
	      abs(dau.pdgId()) == 2212)  // protons 
	    final_state_particles.push_back(&dau);
	}
      }
    }

    // select relevent backgound events
    static const std::vector<long long int> relevant_backgrounds = {
      -211*211,  // pi+pi-
      -321*321,  // K+K-
      -321*211,  // Kpi
      -321*2212, // Kp
      -211*2212, // pi p
      211*13*14, // pi mu nu
      211*13*14, // pi mu nu
      321*13*14, // K mu nu
     -321*13*14, // K mu nu
     2212*13*14, // p mu nu
    -2212*13*14, // p mu nu
     -111*13*13 // pi0 mu mu
    };
    bool relevant_event = false;
    for (auto bkg_signature: relevant_backgrounds)
      if (bkg_signature == signature){
	relevant_event = true;
	break;
      }
    if (not relevant_event) continue;
    
    // find reco tracks matching selected gen level hadrons
    for (const auto& pfCand: *pfCandHandle_.product()){
      if (pfCand.charge() == 0 ) continue;
      if (not pfCand.hasTrackDetails()) continue;
      for (auto hadron: final_state_particles){
	if (deltaR(*hadron, pfCand) > 0.01) continue;
	// check if the hadron is matching one of the selected muons
	bool good_candidate = true;
	for (const auto& good_muon_candidate: good_muon_candidates){
	  if (deltaR(*hadron, good_muon_candidate) < 0.01) {
	    good_candidate=false;
	    break;
	  }
	}
	if (good_candidate)
	  good_muon_candidates.push_back(MuonCand(pfCand));
      }
    }
  }
}

void BxToMuMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder_);

    AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle_.product());
    impactPointExtrapolator_ = &extrapolator;

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BxToMuMuProducer") << "No beam spot available from EventSetup" ;
    }
    
    beamSpot_ = beamSpotHandle.product();
    
    iEvent.getByToken(vertexToken_, pvHandle_);
    // const reco::Vertex & primaryVertex = vertexHandle->front();

    edm::Handle<std::vector<pat::Muon>> muonHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> lostTrackHandle;
    
    iEvent.getByToken(muonToken_, muonHandle);
    iEvent.getByToken(pfCandToken_, pfCandHandle_);
    
    edm::Handle<std::vector<reco::GenParticle> > prunedGenParticleHandle;
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(prunedGenToken_,prunedGenParticleHandle);
      prunedGenParticles_ = prunedGenParticleHandle.product();
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      prunedGenParticles_ = nullptr;
      packedGenParticles_ = nullptr;
    }

    auto nMuons = muonHandle->size();
    auto nPFCands = pfCandHandle_->size();
    // unsigned int lostTrackNumber = useLostTracks_ ? lostTrackHandle->size() : 0;
    
    // Output collection
    auto dimuon  = std::make_unique<pat::CompositeCandidateCollection>();
    auto btokmm  = std::make_unique<pat::CompositeCandidateCollection>();
    auto btokkmm = std::make_unique<pat::CompositeCandidateCollection>();
    AddFourMomenta addP4;

    // good muon candidates
    std::vector<MuonCand> good_muon_candidates;
    for (unsigned int i = 0; i < nMuons; ++i) {
      const pat::Muon & muon = muonHandle->at(i);
      if (not isGoodMuon(muon)) continue;
      good_muon_candidates.push_back(MuonCand(muon, i));
    }
    
    // FIXME: add a switch to disable this feature
    if ( isMC_ ) {
      injectHadronsThatMayFakeMuons(good_muon_candidates);
    }

    // Build dimuon candidates
    if ( good_muon_candidates.size() > 1 ){
      for (unsigned int i = 0; i < good_muon_candidates.size(); ++i) {
	const MuonCand & muon1 = good_muon_candidates.at(i);
	for (unsigned int j = 0; j < good_muon_candidates.size(); ++j) {
	  if (i==j) continue;
	  const MuonCand & muon2 = good_muon_candidates.at(j);
	  // Ensure that muon1.pt > muon2.pt
	  if (muon2.pt() > muon1.pt()) continue;
	  
	  auto mm_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
						   muon2.innerTrack().get());
	  if (maxTwoTrackDOCA_>0 and mm_doca > maxTwoTrackDOCA_) continue;
	  if (diMuonCharge_ && muon1.charge()*muon2.charge()>0) continue;

	  pat::CompositeCandidate dimuonCand;
	  dimuonCand.addDaughter( muon1 , "muon1");
	  dimuonCand.addDaughter( muon2 , "muon2");
	  addP4.set( dimuonCand );

	  dimuonCand.addUserInt("mu1_index", muon1.index());
	  dimuonCand.addUserInt("mu1_pdgId", muon1.pdgId());
	  dimuonCand.addUserFloat("mu1_pt",  muon1.pt());
	  dimuonCand.addUserFloat("mu1_eta", muon1.eta());
	  dimuonCand.addUserFloat("mu1_phi", muon1.phi());
	  dimuonCand.addUserInt("mu2_index", muon2.index());
	  dimuonCand.addUserInt("mu2_pdgId", muon2.pdgId());
	  dimuonCand.addUserFloat("mu2_pt",  muon2.pt());
	  dimuonCand.addUserFloat("mu2_eta", muon2.eta());
	  dimuonCand.addUserFloat("mu2_phi", muon2.phi());
	  dimuonCand.addUserFloat( "doca", mm_doca);

	  // Kalman Vertex Fit
	  auto kalmanMuMuVertexFit = vertexMuonsWithKalmanFitter(muon1, muon2);
	  kalmanMuMuVertexFit.postprocess(*beamSpot_);
	  dimuonCand.addUserInt(   "kalman_valid",    kalmanMuMuVertexFit.valid);
	  dimuonCand.addUserFloat( "kalman_vtx_prob", kalmanMuMuVertexFit.vtxProb);
	  dimuonCand.addUserFloat( "kalman_mass",     kalmanMuMuVertexFit.mass() );
	  dimuonCand.addUserFloat( "kalman_lxy",      kalmanMuMuVertexFit.lxy );
	  dimuonCand.addUserFloat( "kalman_sigLxy",   kalmanMuMuVertexFit.sigLxy );
	  
	  // Kinematic Fits
	  auto kinematicMuMuVertexFit = fillMuMuInfo(dimuonCand,iEvent,muon1,muon2);

	  auto imm = dimuon->size();
	  for (unsigned int k = 0; k < nPFCands; ++k) {
	    pat::PackedCandidate kaonCand1((*pfCandHandle_)[k]);
	    kaonCand1.setMass(KaonMass_);
	    if (deltaR(muon1, kaonCand1) < 0.01 || deltaR(muon2, kaonCand1) < 0.01) continue;
	    if (kaonCand1.charge() == 0 ) continue;
	    if (!kaonCand1.hasTrackDetails()) continue;
	    double mu1_kaon_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
							     kaonCand1.bestTrack());
	    double mu2_kaon_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
							     kaonCand1.bestTrack());
	    // BtoJpsiK
	    bool goodBtoJpsiK = true;
	    if (fabs(kinematicMuMuVertexFit.mass()-3.1)>0.2) goodBtoJpsiK = false;
	    if (abs(kaonCand1.pdgId())!=211) goodBtoJpsiK = false; //Charged hadrons
	    if (kaonCand1.pt()<ptMinKaon_ || abs(kaonCand1.eta())>etaMaxKaon_) goodBtoJpsiK = false;
	    if (maxTwoTrackDOCA_>0 and mu1_kaon_doca> maxTwoTrackDOCA_) goodBtoJpsiK = false;	      
	    if (maxTwoTrackDOCA_>0 and mu2_kaon_doca> maxTwoTrackDOCA_) goodBtoJpsiK = false;	      

	    // BToJpsiKK
	    bool goodBtoJpsiKK = goodBtoJpsiK;

	    double kmm_mass = (muon1.p4()+muon2.p4()+kaonCand1.p4()).mass();
	    if ( kmm_mass<minBKmmMass_ || kmm_mass>maxBKmmMass_ ) goodBtoJpsiK = false;

	    
	    if (goodBtoJpsiK){
	      // fill BtoJpsiK candidate info
	    
	      pat::CompositeCandidate btokmmCand;
	      btokmmCand.addUserInt("mm_index", imm);
	      btokmmCand.addUserFloat("kaon_mu1_doca", mu1_kaon_doca);
	      btokmmCand.addUserFloat("kaon_mu2_doca", mu2_kaon_doca);
	      
	      fillBtoJpsiKInfo(btokmmCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1);
	      fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokmmCand,dimuonCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1);

	      btokmm->push_back(btokmmCand);
	    }
	    // else{
	    // New methods of counting tracks
	    // if (maxTwoTrackDOCA_>0 and mu1_kaon_doca<maxTwoTrackDOCA_ and mu2_kaon_doca<maxTwoTrackDOCA_){
	    // 	auto fit_result = vertexWithKinematicFitter(muon1, muon2, kaonCand1);
	    // 	if ( fit_result.vtxProb()>0.1 ) {
	    // 	  nTracksCompatibleWithTheMuMuVertex++;
	    // 	  double sigDxy = kaonCand1.bestTrack()->dxyError()>0 ? fabs(kaonCand1.bestTrack()->dxy(beamSpot))/kaonCand1.bestTrack()->dxyError():0.0;
	    // 	  if (sigDxy>2.0) nDisplacedTracksCompatibleWithTheMuMuVertex++;
	    // 	}
	    // }	
	    //
	    //  continue;
	    // }

	    if (goodBtoJpsiKK){ // good candidate to consider for JpsiKK
	      for (unsigned int k2 = k+1; k2 < nPFCands; ++k2) { // only works if selection requirements for both kaons are identical
		pat::PackedCandidate kaonCand2((*pfCandHandle_)[k2]);
		kaonCand2.setMass(KaonMass_);
		if (deltaR(muon1, kaonCand2) < 0.01 || deltaR(muon2, kaonCand2) < 0.01) continue;
		if (kaonCand2.charge() == 0 ) continue;
		if (!kaonCand2.hasTrackDetails()) continue;
		double mu1_kaon2_doca = distanceOfClosestApproach(muon1.innerTrack().get(),
							     kaonCand2.bestTrack());
		double mu2_kaon2_doca = distanceOfClosestApproach(muon2.innerTrack().get(),
							     kaonCand2.bestTrack());
	      
		if (abs(kaonCand2.pdgId())!=211) goodBtoJpsiKK = false; //Charged hadrons
		if (kaonCand2.pt()<ptMinKaon_ || abs(kaonCand2.eta())>etaMaxKaon_) goodBtoJpsiKK = false;
		if (maxTwoTrackDOCA_>0 and mu1_kaon2_doca> maxTwoTrackDOCA_) goodBtoJpsiKK = false;	      
		if (maxTwoTrackDOCA_>0 and mu2_kaon2_doca> maxTwoTrackDOCA_) goodBtoJpsiKK = false;	      

		double kkmm_mass = (muon1.p4()+muon2.p4()+kaonCand1.p4()+kaonCand2.p4()).mass();
		if ( kkmm_mass<minBKKmmMass_ || kkmm_mass>maxBKKmmMass_ ) goodBtoJpsiKK = false;
		
		if (goodBtoJpsiKK){
		  // fill BtoJpsiKK candidate info
		  
		  pat::CompositeCandidate btokkmmCand;
		  btokkmmCand.addUserInt("mm_index", imm);
		  int ikmm = -1;
		  if (goodBtoJpsiK) ikmm = btokmm->size()-1;
		  btokkmmCand.addUserInt("kmm_index", ikmm);
		  btokkmmCand.addUserFloat("kaon_mu1_doca", mu1_kaon2_doca);
		  btokkmmCand.addUserFloat("kaon_mu2_doca", mu2_kaon2_doca);
	      
		  fillBstoJpsiKKInfo(btokkmmCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1,kaonCand2);
		  // FIXME
		  // fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokkmmCand,dimuonCand,iEvent,kinematicMuMuVertexFit,muon1,muon2,kaonCand1);

		  btokkmm->push_back(btokkmmCand);
		}
	      }
	    }
	  }                  

	  dimuon->push_back(dimuonCand);
	}
      }
    }
    
    iEvent.put(std::move(dimuon), "DiMuon");
    iEvent.put(std::move(btokmm), "BToKmumu");
    iEvent.put(std::move(btokkmm),"BToKKmumu");
    
}

KalmanVertexFitResult 
BxToMuMuProducer::vertexWithKalmanFitter(std::vector<const reco::Track*> trks, 
					 std::vector<float> masses){
  if (trks.size()!=masses.size()) 
    throw cms::Exception("Error") << "number of tracks and number of masses should match";
  KalmanVertexFitResult results;
  std::vector<reco::TransientTrack> transTrks;
  for (auto trk: trks){
    transTrks.push_back((*theTTBuilder_).build(trk));
  }
  KalmanVertexFitter kvf(true);
  TransientVertex tv = kvf.vertex(transTrks);

  if ( tv.isValid() ){
    results.vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
    results.valid = true;
    results.position = tv.position();
    results.err = tv.positionError();
    if (tv.hasRefittedTracks()){
      assert(tv.refittedTracks().size()==transTrks.size());
      for (unsigned int i=0; i<transTrks.size(); ++i){
	// Is it safe to assume that the order hasn't changed?
	GlobalVector gvP = tv.refittedTracks()[i].trajectoryStateClosestToPoint(tv.position()).momentum();
	// GlobalVector gvP = tv.originalTracks()[i].trajectoryStateClosestToPoint(tv.position()).momentum();
	results.refitVectors.push_back(makeLorentzVectorFromPxPyPzM(gvP.x(), gvP.y(), gvP.z(), masses[i]));
      }
    }
  }
  return results;
}

KalmanVertexFitResult 
BxToMuMuProducer::vertexMuonsWithKalmanFitter(const pat::Muon& muon1,
					      const pat::Muon& muon2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  return vertexWithKalmanFitter(trks,masses);
}


KinematicFitResult 
BxToMuMuProducer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
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
  float muonMassErr(MuonMassErr_);
  for (unsigned int i=0; i<trks.size(); ++i){
    transTrks.push_back((*theTTBuilder_).build(trks[i]));
    particles.push_back(factory.particle(transTrks.back(),masses[i],chi,ndf,muonMassErr));
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

KinematicFitResult
BxToMuMuProducer::vertexMuonsWithKinematicFitter(const pat::Muon& muon1,
						 const pat::Muon& muon2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  return vertexWithKinematicFitter(trks,masses);
}

KinematicFitResult
BxToMuMuProducer::vertexWithKinematicFitter(const pat::Muon& muon1,
					    const pat::Muon& muon2,
					    const pat::PackedCandidate& pion)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( pion.bestTrack() );
  masses.push_back(pionMass_);
  return vertexWithKinematicFitter(trks,masses);
}

KinematicFitResult
BxToMuMuProducer::fitBToKJPsiMuMu( RefCountedKinematicParticle refitMuMu,
				   const pat::PackedCandidate &kaon,
				   bool applyJpsiMassConstraint)
{
  const reco::TransientTrack mmTT = refitMuMu->refittedTransientTrack();
  const reco::TransientTrack kaonTT = theTTBuilder_->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  float MuMu_mass = refitMuMu->currentState().mass();
  float MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

  if ( applyJpsiMassConstraint ){
    MuMu_mass = JPsiMass_;
    MuMu_mass_err = JPsiMassErr_;
  }

  BToKMuMuParticles.push_back(partFactory.particle(mmTT,MuMu_mass,chi,ndf,MuMu_mass_err));
  float kaonMassErr(KaonMassErr_);
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,kaonMassErr));

  KinematicFitResult result; 
  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(BToKMuMuParticles);
  } catch (const std::exception& e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout << "Fit failed. Result is invalid." << std::endl;
    return result;
  }

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

KinematicFitResult
BxToMuMuProducer::fitBToKJPsiMuMuNew( RefCountedKinematicTree jpsiTree,
				      const pat::PackedCandidate& kaon,
				      bool applyJpsiMassConstraint)
{
  KinematicFitResult result; 
  if ( !jpsiTree->isValid()) return result;

  KinematicConstraint* jpsi_mc(0);
  if (applyJpsiMassConstraint){
    ParticleMass jpsi = JPsiMass_;
    // jpsi mass constraint fit
    KinematicParticleFitter csFitter;
    float jp_m_sigma = JPsiMassErr_;
    // FIXME: memory leak
    jpsi_mc = new MassKinematicConstraint(jpsi, jp_m_sigma);
    jpsiTree = csFitter.fit(jpsi_mc, jpsiTree);
  }

  const reco::TransientTrack kaonTT = theTTBuilder_->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  jpsiTree->movePointerToTheTop();
  BToKMuMuParticles.push_back(jpsiTree->currentParticle());
  float kaonMassErr(KaonMassErr_);
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,kaonMassErr));

  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(BToKMuMuParticles);
  } catch (const std::exception& e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cout << "Fit failed. Result is invalid." << std::endl;
    return result;
  }

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

KinematicFitResult
BxToMuMuProducer::fitBToKKMuMu( RefCountedKinematicTree jpsiTree,
				const pat::PackedCandidate& kaon1,
				const pat::PackedCandidate& kaon2,
				bool applyJpsiMassConstraint)
{
  KinematicFitResult result; 
  if ( !jpsiTree->isValid()) return result;

  KinematicConstraint* jpsi_mc(0);
  if (applyJpsiMassConstraint){
    ParticleMass jpsi = JPsiMass_;
    // jpsi mass constraint fit
    KinematicParticleFitter csFitter;
    float jp_m_sigma = JPsiMassErr_;
    // FIXME: memory leak
    jpsi_mc = new MassKinematicConstraint(jpsi, jp_m_sigma);
    jpsiTree = csFitter.fit(jpsi_mc, jpsiTree);
  }

  const reco::TransientTrack kaonTT1 = theTTBuilder_->build(kaon1.bestTrack());
  const reco::TransientTrack kaonTT2 = theTTBuilder_->build(kaon2.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  jpsiTree->movePointerToTheTop();
  float kaonMassErr(KaonMassErr_);
  BToKKMuMuParticles.push_back(jpsiTree->currentParticle());
  BToKKMuMuParticles.push_back(partFactory.particle(kaonTT1,KaonMass_,chi,ndf,kaonMassErr));
  BToKKMuMuParticles.push_back(partFactory.particle(kaonTT2,KaonMass_,chi,ndf,kaonMassErr));

  RefCountedKinematicTree vertexFitTree = fitter.fit(BToKKMuMuParticles);

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

KinematicFitResult
BxToMuMuProducer::vertexMuonsWithPointingConstraint( const pat::Muon& muon1,
						     const pat::Muon& muon2,
						     const reco::Vertex& primaryVertex)
{
  auto kinematicMuMuVertexFit = vertexMuonsWithKinematicFitter(muon1, muon2);
  kinematicMuMuVertexFit.postprocess(*beamSpot_);
  if ( !kinematicMuMuVertexFit.valid()) return KinematicFitResult();
  auto tree = kinematicMuMuVertexFit.refitTree;
  if ( !tree->isValid()) return KinematicFitResult();

  GlobalPoint pv(primaryVertex.position().x(), 
		 primaryVertex.position().y(), 
		 primaryVertex.position().z());
  // FIXIT: potential memory leak
  PointingKinematicConstraint* pointing_constraint = new PointingKinematicConstraint(pv);

  KinematicParticleFitter fitter;

  tree->movePointerToTheTop(); // not sure it's needed

  RefCountedKinematicTree refittedTree;
  try {
    refittedTree = fitter.fit(pointing_constraint,tree);
  } catch (const VertexException &e) {
    // fit failed
    return KinematicFitResult();
  }

  if ( !refittedTree->isValid()) return KinematicFitResult();

  KinematicFitResult result; 
  result.treeIsValid = true;

  refittedTree->movePointerToTheTop();
  result.refitVertex = refittedTree->currentDecayVertex();
  result.refitMother = refittedTree->currentParticle();
  result.refitTree   = refittedTree;

  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  refittedTree->movePointerToTheTop();

  if ( refittedTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(refittedTree->currentParticle());
    } while (refittedTree->movePointerToTheNextChild());
  }
  return result;
}

pair<double,double> BxToMuMuProducer::computeDCA(const pat::PackedCandidate &kaon,
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

const reco::Candidate* BxToMuMuProducer::getGenParticle(const reco::Candidate* cand)
{
  if (not cand) return nullptr;
  auto muon = dynamic_cast<const pat::Muon*>(cand);
  if (muon and muon->genParticle()) return muon->genParticle();
  for (auto const & genParticle: *packedGenParticles_){
      if (dr_match(cand->p4(), genParticle.p4()))
	return &genParticle;
  }
  return nullptr;
}

GenMatchInfo BxToMuMuProducer::getGenMatchInfo( const pat::Muon& muon1,
						const pat::Muon& muon2,
						const pat::PackedCandidate* kaon1,
						const pat::PackedCandidate* kaon2 )
{
  auto result = GenMatchInfo();
  const reco::Candidate*   mm_mother(0);
  assert(prunedGenParticles_);
  assert(packedGenParticles_);
  std::vector<const reco::Candidate*> daughters;

  result.mc_mu1 = getGenParticle(&muon1);
  if (result.mc_mu1){
    result.mu1_pdgId = result.mc_mu1->pdgId();
    result.mu1_pt    = result.mc_mu1->pt();
    if (result.mc_mu1->mother()){
      result.mu1_motherPdgId = result.mc_mu1->mother()->pdgId();
    }
    daughters.push_back(result.mc_mu1);
  }

  result.mc_mu2 = getGenParticle(&muon2);
  if (result.mc_mu2){
    result.mu2_pdgId = result.mc_mu2->pdgId();
    result.mu2_pt    = result.mc_mu2->pt();
    if (result.mc_mu2->mother()){
      result.mu2_motherPdgId = result.mc_mu2->mother()->pdgId();
    }
    daughters.push_back(result.mc_mu2);
  }

  if ( result.mc_mu1 and result.mc_mu2 ){
    if ( (result.mc_mu1->vertex()-result.mc_mu2->vertex()).r() < 1e-4)
      result.mm_vtx    = result.mc_mu1->vertex();
    if ( result.mc_mu1->mother() and result.mc_mu1->mother() == result.mc_mu2->mother() ){
      mm_mother = result.mc_mu1->mother();
      result.match = result.mc_mu1->mother();
      result.mm_mass      = mm_mother->mass();
      result.mm_pt        = mm_mother->pt();
      result.mm_pdgId     = mm_mother->pdgId();
      if (mm_mother->mother()) result.mm_motherPdgId = mm_mother->mother()->pdgId();
      result.mm_prod_vtx = getProductionVertex(mm_mother);
    }
  }
  
  if (kaon1){
    for (auto const & genParticle: *packedGenParticles_){
      if (dr_match(kaon1->p4(),genParticle.p4())){
	result.mc_kaon1 = &genParticle;
	daughters.push_back(result.mc_kaon1);
	result.kaon1_pdgId = genParticle.pdgId();
	result.kaon1_pt    = genParticle.pt();
	if (genParticle.mother(0)){
	  result.kaon1_motherPdgId = genParticle.mother(0)->pdgId();
	}
	break;
      }
    }
    if (daughters.size()==3){
      const auto* mother = find_common_ancestor(daughters);
      if (mother){
	result.match        = mother;
	result.kmm_pdgId    = mother->pdgId();
	result.kmm_mass     = mother->mass();
	result.kmm_pt       = mother->pt();
	result.kmm_prod_vtx = getProductionVertex(mother);
      }
    }
  }
  if (kaon2){
    for (auto const & genParticle: *packedGenParticles_){
      if (dr_match(kaon2->p4(),genParticle.p4())){
	result.mc_kaon2 = &genParticle;
	daughters.push_back(result.mc_kaon2);
	result.kaon2_pdgId = genParticle.pdgId();
	result.kaon2_pt    = genParticle.pt();
	if (genParticle.mother(0)){
	  result.kaon2_motherPdgId = genParticle.mother(0)->pdgId();
	}
	break;
      }
    }
    if (daughters.size()==4){
      const auto* mother = find_common_ancestor(daughters);
      if (mother){
	result.match         = mother;
	result.kkmm_pdgId    = mother->pdgId();
	result.kkmm_mass     = mother->mass();
	result.kkmm_pt       = mother->pt();
	result.kkmm_prod_vtx = getProductionVertex(mother);
      }
    }
  }
  if (daughters.size() > 1){
    const auto* mother = find_common_ancestor(daughters); 
    if (mother){ 
      result.common_mother = mother;
    }
  }

  return result;
}

float BxToMuMuProducer::distanceOfClosestApproach( const reco::GenParticle* track1,
						   const reco::GenParticle* track2)
{
  TwoTrackMinimumDistance md;
  GlobalPoint trk1_pos(track1->vertex().x(), 
		       track1->vertex().y(), 
		       track1->vertex().z());
  GlobalVector trk1_mom(track1->px(),track1->py(),track1->pz());

  GlobalTrajectoryParameters trk1(trk1_pos,trk1_mom,track1->charge(),bFieldHandle_.product());
  GlobalPoint trk2_pos(track2->vertex().x(), 
		       track2->vertex().y(), 
		       track2->vertex().z());
  GlobalVector trk2_mom(track2->px(),track2->py(),track2->pz());

  GlobalTrajectoryParameters trk2(trk2_pos,trk2_mom,track2->charge(),bFieldHandle_.product());
  if ( not md.calculate( trk1, trk2 ) ) return -1.0;
  return md.distance();
}

float BxToMuMuProducer::distanceOfClosestApproach( const reco::Track* track1,
						   const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
BxToMuMuProducer::distanceOfClosestApproach( const reco::Track* track,
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
BxToMuMuProducer::distanceOfClosestApproach( const reco::Track* track,
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

namespace {
  typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > cov33_t;
  typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> > cov66_t;
  typedef ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > cov77_t;
  typedef ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > cov99_t;
  typedef ROOT::Math::SVector<double,9> jac9_t;
  
  cov33_t GlobalError2SMatrix_33(GlobalError m_in) 
  {
    cov33_t m_out;
    for (int i=0; i<3; i++) {
      for (int j=i; j<3; j++)  {
	m_out(i,j) = m_in.matrix()(i,j);
      }
    }
    return m_out;
  }
  
  cov99_t makeCovarianceMatrix(const cov33_t cov_vtx1,
			       const cov77_t cov_vtx2) 
  {
    cov99_t cov;
    cov.Place_at(cov_vtx1,0,0);
    cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
    return cov;
  }

  jac9_t makeJacobianVector3d(const AlgebraicVector3 &vtx1, 
			      const AlgebraicVector3 &vtx2, 
			      const AlgebraicVector3 &momentum) 
  {
    jac9_t jac;
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double factor2 = 1. / ROOT::Math::Mag2(momentum);
    const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
    jac.Place_at(-momentum*factor2,0);
    jac.Place_at( momentum*factor2,3);
    jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
    return jac;
  }
  
  jac9_t makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			      ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum) 
  {
    return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
				AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
				AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
  }

  jac9_t makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
			      const AlgebraicVector3 &momentum) {
    jac9_t jac;
    const double momentumMag = ROOT::Math::Mag(momentum);
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double distMag = ROOT::Math::Mag(dist);
    const double factorPositionComponent = 1./(distMag*momentumMag);
    const double factorMomentumComponent = 1./pow(momentumMag,3);
    jac(0)=-dist(0)*factorPositionComponent;
    jac(1)=-dist(1)*factorPositionComponent;
    jac(3)= dist(0)*factorPositionComponent;
    jac(4)= dist(1)*factorPositionComponent;
    jac(6)= momentum(0)*factorMomentumComponent;
    jac(7)= momentum(1)*factorMomentumComponent;
    return jac;
  }
  
  jac9_t makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			      ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
    return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
				AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
				AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
  }
}

DisplacementInformationIn3D BxToMuMuProducer::compute3dDisplacement(const KinematicFitResult& fit,
								    const reco::VertexCollection& vertices,
								    bool closestIn3D)
{
  DisplacementInformationIn3D result;
  if (not fit.valid()) return result;

  // Potential issue: tracks used to build the candidate could 
  // also be used in the primary vertex fit. One can refit the vertices
  // excluding tracks from the cadndidate. It's not done at the moment 
  // due to non-trivial linkig between primary vertex and its tracks 
  // in MiniAOD. Also not all muons associated to a vertex are really 
  // used in the fit, so the potential bias most likely small.
  
  auto candTransientTrack = fit.refitMother->refittedTransientTrack();

  const reco::Vertex* bestVertex(0);
  int bestVertexIndex(-1);
  double minDistance(999.);

  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    if (closestIn3D){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
      if (impactParameter3D.first and impactParameter3D.second.value() < minDistance){
	minDistance = impactParameter3D.second.value();
	bestVertex = &vertex;
	bestVertexIndex = i;
      }
    } else{
      auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), vertex);
      double distance = fabs(impactParameterZ.second.value());
      if (impactParameterZ.first and distance < minDistance){
	minDistance = distance;
	bestVertex = &vertex;
	bestVertexIndex = i;
      }
    }
  }

  // find second best vertex
  const reco::Vertex* bestVertex2(0);
  int bestVertexIndex2(-1);
  double minDistance2(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    if (closestIn3D){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
      if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance){
	minDistance2 = impactParameter3D.second.value();
	bestVertex2 = &vertex;
	bestVertexIndex2 = i;
      }
    } else{
      auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), vertex);
      double distance = fabs(impactParameterZ.second.value());
      if (impactParameterZ.first and distance < minDistance2 and distance > minDistance){
	minDistance2 = distance;
	bestVertex2 = &vertex;
	bestVertexIndex2 = i;
      }
    }
  }

  if (! bestVertex) return result;

  auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex);
  auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex);
  result.pv = bestVertex;
  result.pvIndex = bestVertexIndex;
  if (impactParameterZ.first) {
    result.longitudinalImpactParameter    = impactParameterZ.second.value();
    result.longitudinalImpactParameterSig = impactParameterZ.second.significance();
    result.longitudinalImpactParameterErr = impactParameterZ.second.error();
  }
  if (impactParameter3D.first) {
    result.distaceOfClosestApproach       = impactParameter3D.second.value();
    result.distaceOfClosestApproachSig    = impactParameter3D.second.value();
    result.distaceOfClosestApproachErr    = impactParameter3D.second.error();
  }

  // compute decay length
  VertexDistance3D distance3D;
  auto dist = distance3D.distance(*bestVertex, fit.refitVertex->vertexState() );
  result.decayLength    = dist.value();
  result.decayLengthErr = dist.error();
  
  VertexDistanceXY distanceXY;
  auto distXY = distanceXY.distance(*bestVertex, fit.refitVertex->vertexState() );

  if (bestVertex2){
    auto impactParameter3D2 = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex2);
    auto impactParameterZ2  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex2);
    result.pv2 = bestVertex2;
    result.pv2Index = bestVertexIndex2;
    if (impactParameterZ2.first) {
      result.longitudinalImpactParameter2    = impactParameterZ2.second.value();
      result.longitudinalImpactParameter2Sig = impactParameterZ2.second.significance();
      result.longitudinalImpactParameter2Err = impactParameterZ2.second.error();
    }
    if (impactParameter3D2.first) {
      result.distaceOfClosestApproach2       = impactParameter3D2.second.value();
      result.distaceOfClosestApproach2Sig    = impactParameter3D2.second.value();
      result.distaceOfClosestApproach2Err    = impactParameter3D2.second.error();
    }

    // compute decay length
    VertexDistance3D distance3D;
    auto dist = distance3D.distance(*bestVertex2, fit.refitVertex->vertexState() );
    result.decayLength2    = dist.value();
    result.decayLength2Err = dist.error();

  }

  // cosAlpha

  TVector3 plab(fit.refitMother->currentState().globalMomentum().x(),
		fit.refitMother->currentState().globalMomentum().y(),
		fit.refitMother->currentState().globalMomentum().z());
  TVector3 p1(bestVertex->x(), bestVertex->y(), bestVertex->z());
  TVector3 p2(fit.refitVertex->vertexState().position().x(), 
	      fit.refitVertex->vertexState().position().y(), 
	      fit.refitVertex->vertexState().position().z());
  TVector3 pDiff = p2-p1;
  TVector3 pDiffXY = TVector3(pDiff.X(), pDiff.Y(), 0.);
  TVector3 ptrans  = TVector3(plab.X(), plab.Y(), 0.);
  result.cosAlpha  = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
  result.cosAlphaXY  = ptrans.Dot(pDiffXY) / (ptrans.Mag() * pDiffXY.Mag());

  // compute decayTime information

  const double massOverC = fit.mass()/TMath::Ccgs();

  // get covariance matrix for error propagation in decayTime calculation
  auto vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(bestVertex->error()),
					     fit.refitMother->currentState().kinematicParametersError().matrix());
  auto vtxDistanceJac3d = makeJacobianVector3d(bestVertex->position(), fit.refitVertex->vertexState().position(), plab);
  auto vtxDistanceJac2d = makeJacobianVector2d(bestVertex->position(), fit.refitVertex->vertexState().position(), plab);

  result.decayTime = dist.value() / plab.Mag() * result.cosAlpha * massOverC;
  result.decayTimeError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;

  result.decayTimeXY = distXY.value() / plab.Perp() * result.cosAlphaXY * massOverC;
  result.decayTimeXYError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;
    
  return result;
}


// GenInfoSummary HeavyFlavDileptonNtupleMakerMiniAOD::getGenInfoSummary(const edm::Event& iEvent){
//   GenInfoSummary summary;
  
//   edm::Handle<std::vector<reco::GenParticle> > pruned;
//   iEvent.getByToken(pruned_gen_token,pruned);

//   // Packed particles are all the status 1, so usable to remake jets
//   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
//   edm::Handle<edm::View<pat::PackedGenParticle> > packed;
//   iEvent.getByToken(packed_gen_token,packed);

//   //let's try to find all status1 originating directly from a B meson decay
//   for (auto const & bMeson: *pruned){
//     if(abs(bMeson.pdgId()) != 521 and abs(bMeson.pdgId()) != 511) continue;
//     std::vector<const reco::Candidate*> b_muons;
//     std::vector<const reco::Candidate*> b_hadrons;
//     std::vector<const reco::Candidate*> other_muons;
//     for (auto const& pa: *packed){ 
//       //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
//       const reco::Candidate * mother = pa.mother(0) ;
//       // FIXME: may need to check particle status to avoid intermediate gen particles
//       if(mother != nullptr && isAncestor( &bMeson , mother)){
// 	if (abs(pa.pdgId())==13)
// 	  b_muons.push_back(&pa);
// 	if (abs(pa.pdgId())==211 or abs(pa.pdgId())==321)
// 	  b_hadrons.push_back(&pa);
//       }else{
// 	if (abs(pa.pdgId())==13)
// 	  other_muons.push_back(&pa);
//       }
//     }
//     if (b_muons.size()==2 && b_hadrons.size()>0){
//       summary.b = &bMeson;
//       summary.b_muons = b_muons;
//       summary.b_hadrons = b_hadrons;
//       std::sort(other_muons.begin(),other_muons.end(), [](const reco::Candidate* a,const reco::Candidate* b){
// 	  return a->pt() > b->pt();
// 	});
//       summary.other_muons = other_muons;
//     }
//   }

//   return summary;
// }

void BxToMuMuProducer::setupTmvaReader(TMVA::Reader& reader, std::string file){
  reader.AddVariable("fls3d",    & bdtData_.fls3d);
  reader.AddVariable("alpha",    & bdtData_.alpha);
  reader.AddVariable("pvips",    & bdtData_.pvips);
  reader.AddVariable("iso",      & bdtData_.iso);
  reader.AddVariable("chi2dof",  & bdtData_.chi2dof);
  reader.AddVariable("docatrk",  & bdtData_.docatrk);
  reader.AddVariable("closetrk", & bdtData_.closetrk);
  reader.AddVariable("m1iso",    & bdtData_.m1iso);
  reader.AddVariable("m2iso",    & bdtData_.m2iso);
  reader.AddVariable("eta",      & bdtData_.eta);
  reader.AddSpectator("m",       & bdtData_.m);
  reader.BookMVA("BDTG",file);
}

float
BxToMuMuProducer::computeAnalysisBDT(unsigned int event_idx)
{
  switch (event_idx){
  case 0:
    return bdtReader0_.EvaluateMVA("BDTG");
    break;
  case 1:
    return bdtReader1_.EvaluateMVA("BDTG");
    break;
  case 2:
    return bdtReader2_.EvaluateMVA("BDTG");
    break;
  default:
    throw cms::Exception("FatalError") << "event index must be in [0-3] range\n";
  }
}


DEFINE_FWK_MODULE(BxToMuMuProducer);

//  LocalWords:  vertices
