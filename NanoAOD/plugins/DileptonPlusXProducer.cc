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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

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
#include "RecoVertex/KinematicFit/interface/ColinearityKinematicConstraintT.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitterT.h"
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
#include <math.h>

#include "Bmm5/NanoAOD/interface/XGBooster.h"
#include "Bmm5/NanoAOD/interface/KinFitUtils.h"
#include "Bmm5/NanoAOD/interface/KinematicFitResult.h"
#include "Bmm5/NanoAOD/interface/Displacement.h"
#include "Bmm5/NanoAOD/interface/CommonTools.h"
#include "Bmm5/NanoAOD/interface/Candidate.h"

// 
// DileptonPlusXProducer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;

using namespace bmm;

namespace {
  bool overlap(const reco::Track* t1, const reco::Track* t2){
    assert(t1);
    assert(t2);
    return deltaR(*t1, *t2) < 0.01;
  }

  bool overlap(const bmm::Candidate& lep, const pat::PackedCandidate& can){
    return overlap(lep.track(), can.bestTrack());
  }

  bool overlap(const bmm::Candidate& lep1, const bmm::Candidate& lep2){
    return overlap(lep1.track(), lep2.track());
  }

  float momentum_resolution(const pat::Photon& photon){
    if (fabs(photon.eta()) < 0.4) return photon.p() * 0.025;
    if (fabs(photon.eta()) < 0.8) return photon.p() * 0.025;
    if (fabs(photon.eta()) < 1.2) return photon.p() * 0.031;
    if (fabs(photon.eta()) < 1.6) return photon.p() * 0.040;
    if (fabs(photon.eta()) < 2.0) return photon.p() * 0.045;
    return photon.p() * 0.044;
  }
  
  struct GenMatchInfo{
    int l1_pdgId{0}, l1_motherPdgId{0}, l2_pdgId{0}, l2_motherPdgId{0},
      kaon1_pdgId{0}, kaon1_motherPdgId{0}, kaon2_pdgId{0}, kaon2_motherPdgId{0},
      photon_pdgId{0}, photon_motherPdgId{0},
      ll_pdgId{0}, ll_motherPdgId{0}, kll_pdgId{0}, kll_motherPdgId{0},
      kkll_pdgId{0}, llg_pdgId{0};
    float l1_pt{0}, l2_pt{0}, kaon1_pt{0}, kaon2_pt{0}, photon_pt{0}, ll_mass{0}, ll_pt{0},
      kll_mass{0}, kkll_mass{0}, llg_mass{0}, kll_pt{0}, kkll_pt{0}, llg_pt{0};
    math::XYZPoint ll_prod_vtx, ll_vtx, kll_prod_vtx, kkll_prod_vtx, llg_prod_vtx;
    const reco::Candidate* mc_l1{nullptr}, *mc_l2{nullptr}, *mc_kaon1{nullptr},
      *mc_kaon2{nullptr}, *mc_photon{nullptr}, *match{nullptr}, *common_mother{nullptr};
    
    const reco::GenParticle* gen_l1(){
      return dynamic_cast<const reco::GenParticle*>(mc_l1);
    }
    const reco::GenParticle* gen_l2(){
      return dynamic_cast<const reco::GenParticle*>(mc_l2);
    }
  };

  struct GenEventInfo{};

  struct BdtReaderData {
    float fls3d, alpha, pvips, iso, chi2dof, docatrk, closetrk, m1iso, m2iso, eta, m;
  };

}


using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class DileptonPlusXProducer : public edm::stream::EDProducer<> {
    
public:
    
  explicit DileptonPlusXProducer(const edm::ParameterSet &iConfig);
    
  ~DileptonPlusXProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool preprocess(pat::CompositeCandidate& candidate,
		  const edm::Event& iEvent,
		  const bmm::Candidate& lepton1,
		  const bmm::Candidate& lepton2);

  void buildLLXCandidates(pat::CompositeCandidateCollection& llk,
			  pat::CompositeCandidateCollection& llkk,
			  const edm::Event& iEvent,
			  const KinematicFitResult& kinematicLLVertexFit,
			  const pat::CompositeCandidate& dileptonCand,
			  int ll_index,
			  const bmm::Candidate& lepton1,
			  const bmm::Candidate& lepton2);
  void
  fillDstarInfo(pat::CompositeCandidateCollection& dstar_collection,
		const edm::Event& iEvent,
		const KinematicFitResult& d0VertexFit,
		const pat::CompositeCandidate& d0Cand,
		const pat::PackedCandidate& soft_pion,
		int mm_index,
		int hh_index,
		const bmm::Candidate& daughter1,
		const bmm::Candidate& daughter2);

  void
  buildDstarCandidates(pat::CompositeCandidateCollection& dstar_collection,
		       pat::CompositeCandidateCollection& hh_collection,
		       const edm::Event& iEvent,
		       const pat::PackedCandidate& had1,
		       const pat::PackedCandidate& had2);


  bool isGoodMuon(const pat::Muon& muon);
  bool isGoodTrack(const pat::PackedCandidate& cand);
  bool isGoodMuonCandidateFromTrack(const pat::PackedCandidate& cand);
  bool isGoodHadron(const pat::PackedCandidate& cand);
  bool isGoodElectron(const pat::Electron& el);
    
  KalmanVertexFitResult 
  vertexWithKalmanFitter( std::vector<const reco::Track*> trks, 
			  std::vector<float> masses);

  KalmanVertexFitResult 
  vertexLeptonsWithKalmanFitter(const bmm::Candidate& lepton1,
				const bmm::Candidate& lepton2);

  KinematicFitResult 
  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  KinematicFitResult 
  vertexLeptonsWithKinematicFitter(const bmm::Candidate& lepton1,
				   const bmm::Candidate& lepton2);

  KinematicFitResult 
  vertexKaonsWithKinematicFitter(const pat::PackedCandidate& pfCand1,
				 const pat::PackedCandidate& pfCand2);
  
  KinematicFitResult 
  vertexWithKinematicFitter(const bmm::Candidate& lepton1,
			    const bmm::Candidate& lepton2,
			    const pat::PackedCandidate& pfCand);

  template <typename T> std::vector<const reco::Track*>
  getGoodTracksToRefitPV(int pvIndex, T track_to_ignore);
  
  template <typename T> std::vector<const reco::Track*>
  getGoodTracksToRefitPV(int pvIndex, std::vector<T> ignoreTracks);

  std::pair<KinematicFitResult, KinematicFitResult>
  refitWithVertexConstraint(const reco::Track& track,
			    int pvIndex);

  KinematicFitResult
  fitBToKLL(const bmm::Candidate& lepton1,
	    const bmm::Candidate& lepton2,
	    const pat::PackedCandidate& kaon,
	    float mass_constraint=-1.0);

  KinematicFitResult
  fitBToKKLL( const bmm::Candidate& lepton1,
	      const bmm::Candidate& lepton2,
	      const pat::PackedCandidate& kaon1,
	      const pat::PackedCandidate& kaon2,
	      float ll_mass_constraint=-1.0,
	      float kk_mass_constraint=-1.0);
  
  KinematicFitResult
  fitDstar(const bmm::Candidate& lepton1,
	   const bmm::Candidate& lepton2,
	   const pat::PackedCandidate& pion,
	   float mass_constraint);

  KinematicFitResult
  fitLLGamma( RefCountedKinematicTree tree,
	      const pat::Photon& photon,
	      float mass_constraint=-1.0);

  KinematicFitResult
  fitLLGammaConv( RefCountedKinematicTree tree,
		  const pat::CompositeCandidate& photon,
		  float mass_constraint=-1.0);

  KinematicFitResult
  vertexLeptonsWithPointingConstraint( const bmm::Candidate& lepton1,
				       const bmm::Candidate& lepton2,
				       const reco::Vertex& primaryVertex);

  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
   				 reco::BeamSpot beamSpot);
  GenMatchInfo getGenMatchInfo( const bmm::Candidate& lepton1,
				const bmm::Candidate& lepton2,
				const pat::PackedCandidate* kaon1 = 0,
				const pat::PackedCandidate* kaon2 = 0,
				const reco::Candidate* photon = 0);
  
  const reco::Candidate* getGenParticle(const bmm::Candidate& cand);

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
			     const VertexState& vertex_state);
  Measurement1D 
  distanceOfClosestApproach( const reco::Track* track,
			     const reco::Vertex& vertex);

  bmm::Displacements
  compute3dDisplacement(const KinematicFitResult& fit,
			bool closestIn3D = true);

  CloseTrackInfo 
  findTracksCompatibleWithTheVertex(const bmm::Candidate& lepton1,
				    const bmm::Candidate& lepton2,
				    const KinematicFitResult& fit,
				    double maxDoca=0.03,
				    std::vector<const pat::PackedCandidate*> ignoreTracks = 
				    std::vector<const pat::PackedCandidate*>());
  float
  computeTrkLeptonIsolation(const bmm::Candidate& lepton, 
			    const bmm::Candidate& the_other_lepton,
			    unsigned int primaryVertexIndex,
			    float minPt=0.5, float dR=0.5,
			    std::vector<const pat::PackedCandidate*> ignoreTracks = 
			    std::vector<const pat::PackedCandidate*>());
  float
  computeTrkDileptonIsolation(const bmm::Candidate& lepton1, 
			      const bmm::Candidate& lepton2,
			      unsigned int primaryVertexIndex,
			      float minPt=0.9, float dR=0.7,
			      std::vector<const pat::PackedCandidate*> ignoreTracks = 
			      std::vector<const pat::PackedCandidate*>());

  float
  otherVertexMaxProb(const bmm::Candidate& lepton1, 
		     const bmm::Candidate& lepton2,
		     float min_pt = 0.5,
		     float max_doca = 0.1,
		     std::vector<const pat::PackedCandidate*> ignoreTracks = 
		     std::vector<const pat::PackedCandidate*>());

  void 
  fillBtoKllInfo(pat::CompositeCandidate& bCand,
		 const edm::Event& iEvent,
		 const bmm::Candidate& muon1,
		 const bmm::Candidate& muon2,
		 const pat::PackedCandidate & kaon); 
  void
  fillBtoKKllInfo(pat::CompositeCandidate& bCand,
		  const edm::Event& iEvent,
		  const bmm::Candidate& muon1,
		  const bmm::Candidate& muon2,
		  const pat::PackedCandidate & kaon1,
		  const pat::PackedCandidate & kaon2); 
  void
  fillLLGammaGenInfo(pat::CompositeCandidate& llgCand,
		     const edm::Event& iEvent,
		     const bmm::Candidate& muon1,
		     const bmm::Candidate& muon2,
		     const reco::Candidate & photon);
  void
  fillLLGammaInfo(pat::CompositeCandidate& mmgCand,
		  const edm::Event& iEvent,
		  const bmm::Candidate& muon1,
		  const bmm::Candidate& muon2,
		  const pat::Photon & photon);
  void
  fillLLGammaConvInfo(pat::CompositeCandidate& mmgCand,
		      const edm::Event& iEvent,
		      const bmm::Candidate& muon1,
		      const bmm::Candidate& muon2,
		      const pat::CompositeCandidate& photon);
  void 
  fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(pat::CompositeCandidate& btokllCand,
					       const pat::CompositeCandidate& dimuonCand,
					       const edm::Event& iEvent,
					       const KinematicFitResult& kinematicLLVertexFit,
					       const bmm::Candidate& muon1,
					       const bmm::Candidate& muon2,
					       const pat::PackedCandidate & kaon);
 
  void 
  fillMvaInfoForLLGamma(pat::CompositeCandidate& mmg,
			const pat::CompositeCandidate& ll,
			const edm::Event& iEvent,
			const KinematicFitResult& kinematicLLVertexFit,
			const bmm::Candidate& muon1,
			const bmm::Candidate& muon2,
			const pat::Photon & photon);
  
  KinematicFitResult 
  fillDileptonInfo(pat::CompositeCandidate& dileptonCand,
		   const edm::Event& iEvent,
		   const bmm::Candidate& lepton1,
		   const bmm::Candidate& lepton2); 

  void 
  injectHadronsThatMayFakeMuonsInMC(std::vector<bmm::Candidate>& good_lepton_candidates);

  void 
  injectBhhHadrons(std::vector<bmm::Candidate>& good_lepton_candidates);

  void 
  injectJpsiTracks(std::vector<bmm::Candidate>& good_lepton_candidates);

  float  computeAnalysisBDT(unsigned int event_idx);
  
  void setupTmvaReader(TMVA::Reader& reader, std::string file);


  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
  edm::EDGetTokenT<std::vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> conversionToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >   prunedGenToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<reco::GenParticle>* prunedGenParticles_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  const TransientTrackBuilder* theTTBuilder_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBuilderToken_;
  const MagneticField* bField_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;

  bool isMC_;

  double ptMinMu_;
  double etaMaxMu_;
  double ptMinKaon_;
  double etaMaxKaon_;
  double ptMinEl_;
  double etaMaxEl_;

  double DCASigMinKaon_;
  bool   diLeptonCharge_;
  double minBKllMass_;
  double maxBKllMass_;
  double minLLGammaMass_;
  double maxLLGammaMass_;
  double minGammaPt_;
  double minBKKllMass_;
  double maxBKKllMass_;
  double maxTwoTrackDOCA_;
  bool   injectMatchedBtohh_;
  bool   injectBtohh_;
  bool   injectJpsiTracks_;
  bool   recoElElX_;
  bool   recoMuMuGamma_;
  bool   recoMuMuGammaConv_;
  bool   recoDstar_;
  bool   recoD0pipi_;
  bool   recoD0Kpi_;
  double minBhhTrkPt_;
  double minDhhTrkPt_;
  double minJpsiHadronPt_;
  double minBhhMass_;
  double maxBhhMass_;
  double maxBhhTrkEta_;
  double maxDhhTrkEta_;
  double minBhhSigLxy_;
  double minBhhVtxProb_;
  double minD0Mass_;
  double maxD0Mass_;
  double min_dm_;
  double max_dm_;
  
  BdtReaderData bdtData_;
  TMVA::Reader  bdtReader0_;
  TMVA::Reader  bdtReader1_;
  TMVA::Reader  bdtReader2_;
  std::vector<XGBooster> xgBoosters_;
};

DileptonPlusXProducer::DileptonPlusXProducer(const edm::ParameterSet &iConfig):
  beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
  beamSpot_(nullptr),
  vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
  muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  electronToken_( consumes<std::vector<pat::Electron>> ( iConfig.getParameter<edm::InputTag>( "electronCollection" ) ) ),
  photonToken_( consumes<std::vector<pat::Photon>> ( iConfig.getParameter<edm::InputTag>( "photonCollection" ) ) ),
  conversionToken_( consumes<pat::CompositeCandidateCollection> ( iConfig.getParameter<edm::InputTag>( "conversionCollection" ) ) ),
  pfCandToken_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
  prunedGenToken_( consumes<std::vector<reco::GenParticle>> ( iConfig.getParameter<edm::InputTag>( "prunedGenParticleCollection" ) ) ),
  packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( edm::InputTag( "packedGenParticles" ) ) ),
  prunedGenParticles_(nullptr),
  packedGenParticles_(nullptr),
  theTTBuilder_(nullptr),
  theTTBuilderToken_(esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"})),
  bField_(nullptr),
  bFieldToken_(esConsumes()),
  impactPointExtrapolator_(0),
  isMC_(            iConfig.getParameter<bool>( "isMC" ) ),
  ptMinMu_(         iConfig.getParameter<double>( "MuonMinPt" ) ),
  etaMaxMu_(        iConfig.getParameter<double>( "MuonMaxEta" ) ),
  ptMinKaon_(       iConfig.getParameter<double>( "KaonMinPt" ) ),
  etaMaxKaon_(      iConfig.getParameter<double>( "KaonMaxEta" ) ),
  ptMinEl_(         iConfig.getParameter<double>( "ElectronMinPt" ) ),
  etaMaxEl_(        iConfig.getParameter<double>( "ElectronMaxEta" ) ),
  DCASigMinKaon_(   iConfig.getParameter<double>( "KaonMinDCASig" ) ),
  diLeptonCharge_(    iConfig.getParameter<bool>( "DiLeptonChargeCheck" ) ),
  minBKllMass_(     iConfig.getParameter<double>( "minBKllMass" ) ),
  maxBKllMass_(     iConfig.getParameter<double>( "maxBKllMass" ) ),
  minLLGammaMass_(     iConfig.getParameter<double>( "minLLGammaMass" ) ),
  maxLLGammaMass_(     iConfig.getParameter<double>( "maxLLGammaMass" ) ),
  minGammaPt_(     iConfig.getParameter<double>( "minGammaPt" ) ),
  minBKKllMass_(    iConfig.getParameter<double>( "minBKKllMass" ) ),
  maxBKKllMass_(    iConfig.getParameter<double>( "maxBKKllMass" ) ),
  maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
  injectMatchedBtohh_( iConfig.getParameter<bool>( "injectMatchedBtohh" ) ),
  injectBtohh_(        iConfig.getParameter<bool>( "injectBtohh" ) ),
  injectJpsiTracks_(  iConfig.getParameter<bool>( "injectJpsiTracks" ) ),
  recoElElX_(      iConfig.getParameter<bool>( "recoElElX" ) ),
  recoMuMuGamma_(      iConfig.getParameter<bool>( "recoMuMuGamma" ) ),
  recoMuMuGammaConv_(  iConfig.getParameter<bool>( "recoMuMuGammaConv" ) ),
  recoDstar_(      iConfig.getParameter<bool>( "recoDstar" ) ),
  recoD0pipi_(     iConfig.getParameter<bool>( "recoD0pipi" ) ),
  recoD0Kpi_(      iConfig.getParameter<bool>( "recoD0Kpi" ) ),
  minBhhTrkPt_(        iConfig.getParameter<double>( "minBhhHadronPt" ) ),
  minDhhTrkPt_(        iConfig.getParameter<double>( "minDhhHadronPt" ) ),
  minJpsiHadronPt_( iConfig.getParameter<double>( "minJpsiHadronPt" ) ),
  minBhhMass_(      iConfig.getParameter<double>( "minBhhMass" ) ),
  maxBhhMass_(      iConfig.getParameter<double>( "maxBhhMass" ) ),
  maxBhhTrkEta_(       iConfig.getParameter<double>( "maxBhhHadronEta" ) ),
  maxDhhTrkEta_(       iConfig.getParameter<double>( "maxDhhHadronEta" ) ),
  minBhhSigLxy_(    iConfig.getParameter<double>( "minBhhSigLxy" ) ),
  minBhhVtxProb_(   iConfig.getParameter<double>( "minBhhVtxProb" ) ),
  minD0Mass_(      iConfig.getParameter<double>( "minD0Mass" ) ),
  maxD0Mass_(      iConfig.getParameter<double>( "maxD0Mass" ) ),
  min_dm_(      iConfig.getParameter<double>( "minDm" ) ),
  max_dm_(      iConfig.getParameter<double>( "maxDm" ) ),
  bdtReader0_("!Color:Silent"),
  bdtReader1_("!Color:Silent"),
  bdtReader2_("!Color:Silent")
{
  produces<pat::CompositeCandidateCollection>("MuMu");
  produces<pat::CompositeCandidateCollection>("ElEl");
  produces<pat::CompositeCandidateCollection>("HH");
  produces<pat::CompositeCandidateCollection>("BToKmumu");
  produces<pat::CompositeCandidateCollection>("BToKee");
  produces<pat::CompositeCandidateCollection>("BToKKmumu");
  produces<pat::CompositeCandidateCollection>("BToKKee");
  produces<pat::CompositeCandidateCollection>("BToMuMuGamma");
  produces<pat::CompositeCandidateCollection>("Dstar");
    
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

bool DileptonPlusXProducer::isGoodMuon(const pat::Muon& muon){
  if ( not muon.isLooseMuon() ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < ptMinMu_ || fabs(muon.eta()) > etaMaxMu_ ) return false;
  return true;
}

bool DileptonPlusXProducer::isGoodTrack(const pat::PackedCandidate& cand){
  if (cand.charge() == 0) return false;
  if (not cand.hasTrackDetails()) return false;
  if (not cand.bestTrack()->quality(reco::Track::highPurity)) return false;
  if (isnan(cand.pt())) return false;
  return true;
}

bool DileptonPlusXProducer::isGoodHadron(const pat::PackedCandidate& cand)
{
  if (not isGoodTrack(cand)) return false;
  if (abs(cand.pdgId()) != 211) return false;
  return true;
}  

bool DileptonPlusXProducer::isGoodMuonCandidateFromTrack(const pat::PackedCandidate& cand)
{
  if (not isGoodTrack(cand)) return false; 
  if (cand.pt() < minJpsiHadronPt_) return false;
  if (abs(cand.pdgId()) != 211) return false;
  return true;
}  

bool DileptonPlusXProducer::isGoodElectron(const pat::Electron& el){
  if ( el.pt() < ptMinEl_ || fabs(el.eta()) > etaMaxEl_ ) return false;
  return true;
}

namespace {
  void addFitInfo(pat::CompositeCandidate& cand, const KinematicFitResult& fit, std::string name, 
		  const bmm::Displacements& displacements = bmm::Displacements(),
		  int firstMuonDaughterIndex = -1, int secondMuonDaughterIndex = -1,
		  int firstKaonDaughterIndex = -1, int secondKaonDaughterIndex = -1 ){
    cand.addUserInt(   name+"_valid",       fit.valid() );
    cand.addUserFloat( name+"_vtx_prob",    fit.vtxProb() );
    cand.addUserFloat( name+"_vtx_chi2dof", fit.chi2()>0?fit.chi2()/fit.ndof():-1);
    cand.addUserFloat( name+"_mass",        fit.mass() );
    cand.addUserFloat( name+"_massErr",     fit.massErr() );
    cand.addUserFloat( name+"_lxy",         fit.lxy() );
    cand.addUserFloat( name+"_sigLxy",      fit.sigLxy() );
    cand.addUserFloat( name+"_alphaBS",     fit.alphaBS() );
    cand.addUserFloat( name+"_alphaBSErr",  fit.alphaBSErr() );
    auto vtx_position = fit.vtx_position();
    auto vtx_error = fit.vtx_error();
    cand.addUserFloat( name+"_vtx_x",       vtx_position.x() );
    cand.addUserFloat( name+"_vtx_xErr",    sqrt(vtx_error.cxx()) );
    cand.addUserFloat( name+"_vtx_y",       vtx_position.y() );
    cand.addUserFloat( name+"_vtx_yErr",    sqrt(vtx_error.cyy()) );
    cand.addUserFloat( name+"_vtx_z",       vtx_position.z() );
    cand.addUserFloat( name+"_vtx_zErr",    sqrt(vtx_error.czz()) );
    cand.addUserFloat( name+"_pt",          fit.p3().perp() );
    cand.addUserFloat( name+"_eta",         fit.p3().eta() );
    cand.addUserFloat( name+"_phi",         fit.p3().phi() );
    
    // displacement information
    
    for (const auto& displacement: displacements){
      std::string prefix = name;
      if (displacement.name() != "pv")
	prefix += displacement.name();
      cand.addUserFloat( prefix + "_alpha",       displacement.alpha());
      cand.addUserFloat( prefix + "_alphaErr",    displacement.alphaErr());

      // IP info
      cand.addUserFloat( prefix + "_l3d",         displacement.decayLength());
      cand.addUserFloat( prefix + "_sl3d",        displacement.decayLengthSig());
      cand.addUserFloat( prefix + "_pv_z",        displacement.prodVertex().z());
      cand.addUserFloat( prefix + "_pv_zErr",     displacement.prodVertex().zError());
      cand.addUserFloat( prefix + "_pvip",        displacement.distaceOfClosestApproach());
      cand.addUserFloat( prefix + "_spvip",       displacement.distaceOfClosestApproachSig());
      cand.addUserFloat( prefix + "_pvipErr",     displacement.distaceOfClosestApproachErr());
      // cand.addUserFloat( prefix + "_pv2ip",       displacement.distaceOfClosestApproach2());
      // cand.addUserFloat( prefix + "_spv2ip",      displacement.distaceOfClosestApproach2Sig());
      // cand.addUserFloat( prefix + "_pv2ipErr",    displacement.distaceOfClosestApproach2Err());
      cand.addUserFloat( prefix + "_pvlip",       displacement.longitudinalImpactParameter());
      cand.addUserFloat( prefix + "_pvlipSig",    displacement.longitudinalImpactParameterSig());
      cand.addUserFloat( prefix + "_pvlipErr",    displacement.longitudinalImpactParameterErr());
      // cand.addUserFloat( prefix + "_pv2lip",      displacement.longitudinalImpactParameter2());
      // cand.addUserFloat( prefix + "_pv2lipSig",   displacement.longitudinalImpactParameter2Sig());
      // cand.addUserFloat( prefix + "_pv2lipErr",   displacement.longitudinalImpactParameter2Err());
      cand.addUserInt(   prefix + "_pvIndex",     displacement.pvIndex());

      // DecayTime
      cand.addUserFloat( prefix + "_tau",         displacement.decayTime());
      cand.addUserFloat( prefix + "_taue",        displacement.decayTimeError());
      cand.addUserFloat( prefix + "_tauxy",       displacement.decayTimeXY());
      cand.addUserFloat( prefix + "_tauxye",      displacement.decayTimeXYError());
    }

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
      cand.addUserFloat( name+"_kaon1_pt",     fit.dau_p3(firstKaonDaughterIndex).perp() );
      cand.addUserFloat( name+"_kaon1_eta",    fit.dau_p3(firstKaonDaughterIndex).eta() );
      cand.addUserFloat( name+"_kaon1_phi",    fit.dau_p3(firstKaonDaughterIndex).phi() );
    }
    if (secondKaonDaughterIndex>=0){
      cand.addUserFloat( name+"_kaon2_pt",     fit.dau_p3(secondKaonDaughterIndex).perp() );
      cand.addUserFloat( name+"_kaon2_eta",    fit.dau_p3(secondKaonDaughterIndex).eta() );
      cand.addUserFloat( name+"_kaon2_phi",    fit.dau_p3(secondKaonDaughterIndex).phi() );
    }
  }
}

CloseTrackInfo 
DileptonPlusXProducer::findTracksCompatibleWithTheVertex(const bmm::Candidate& lepton1,
							 const bmm::Candidate& lepton2,
							 const KinematicFitResult& fit, 
							 double maxDoca,
							 std::vector<const pat::PackedCandidate*> ignoreTracks)
{
  CloseTrackInfo result;
  if (not fit.valid()) return result;
  for (const auto& pfCand: *pfCandHandle_.product()){
    if (not isGoodTrack(pfCand)) continue; 
    if (overlap(lepton1, pfCand) || overlap(lepton2, pfCand)) continue;

    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (deltaR(*trk, pfCand) < 0.01){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    
    double mu1_kaon_doca = distanceOfClosestApproach(lepton1.track(),
						     pfCand.bestTrack());
    double mu2_kaon_doca = distanceOfClosestApproach(lepton2.track(),
						     pfCand.bestTrack());
    if (mu1_kaon_doca>maxDoca or mu2_kaon_doca>maxDoca) continue;
    
    CloseTrack track;
    track.pfCand = &pfCand;
    auto doca = distanceOfClosestApproach(pfCand.bestTrack(),fit.vtx_state());
    track.svDoca = doca.value();
    track.svDocaErr = doca.error();

    // add PV doca
    if (pfCand.vertexRef().key()<pvHandle_->size()){
      doca = distanceOfClosestApproach(pfCand.bestTrack(),pvHandle_->at(pfCand.vertexRef().key()) );
      track.pvDoca = doca.value();
      track.pvDocaErr = doca.error();
    }
    
    auto fit_result = vertexWithKinematicFitter(lepton1, lepton2, pfCand);
    if (fit_result.valid()){
      track.svProb = fit_result.vtxProb();
      track.impactParameterSignificanceBS = pfCand.bestTrack()->dxyError()>0 ? fabs(pfCand.bestTrack()->dxy(*beamSpot_))/pfCand.bestTrack()->dxyError():0.0;
    }
    result.tracks.push_back(track);
  }

  return result;
}

float
DileptonPlusXProducer::computeTrkLeptonIsolation(const bmm::Candidate& theLepton,
						 const bmm::Candidate& theOtherLepton, 
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
    if (not isGoodTrack(pfCand)) continue; 
    if (pfCand.pt() < minPt) continue;
    if (pfCand.vertexRef().key() != primaryVertexIndex) continue;
    if (overlap(theLepton, pfCand) || overlap(theOtherLepton, pfCand)) continue;
    if (deltaR(theLepton, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }

  return theLepton.pt() / (theLepton.pt() + sumPt);
}

float
DileptonPlusXProducer::computeTrkDileptonIsolation(const bmm::Candidate& lepton1,
						   const bmm::Candidate& lepton2, 
						   unsigned int primaryVertexIndex,
						   float minPt, float dR,
						   std::vector<const pat::PackedCandidate*> ignoreTracks)
{
  float sumPt(0);
  auto b_p4 = lepton1.p4() + lepton2.p4();
  for (const auto& pfCand: *pfCandHandle_.product()){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (not isGoodTrack(pfCand)) continue;
    if (pfCand.pt() < minPt) continue;
    if (pfCand.vertexRef().key() != primaryVertexIndex) continue;
    if (overlap(lepton1, pfCand) || overlap(lepton2, pfCand)) continue;
    if (deltaR(b_p4, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }

  return b_p4.pt()/(b_p4.pt()+sumPt);
}


float
DileptonPlusXProducer::otherVertexMaxProb(const bmm::Candidate& lepton1, 
					  const bmm::Candidate& lepton2,
					  float minPt,
					  float max_doca,
					  std::vector<const pat::PackedCandidate*> ignoreTracks){
  float bestLep1Vtx = 0;
  float bestLep2Vtx = 0;
  KalmanVertexFitter kvf;
  std::vector<reco::TransientTrack> transTrksForLep1Vertex;
  transTrksForLep1Vertex.push_back((*theTTBuilder_).build(lepton1.track()));
  std::vector<reco::TransientTrack> transTrksForLep2Vertex;
  transTrksForLep2Vertex.push_back((*theTTBuilder_).build(lepton2.track()));


  for (const auto& pfCand: *pfCandHandle_.product()){
    if (not isGoodTrack(pfCand)) continue; 
    if (pfCand.pt() < minPt) continue;
    if (overlap(lepton1, pfCand) || overlap(lepton2, pfCand)) continue;

    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;

    double lep1_doca = distanceOfClosestApproach(lepton1.track(),
						 pfCand.bestTrack());
    double lep2_doca = distanceOfClosestApproach(lepton2.track(),
						 pfCand.bestTrack());
    if (lep1_doca < max_doca and lep1_doca < lep2_doca){
      // first  lepton is closer - check vertex probability
      transTrksForLep1Vertex.push_back((*theTTBuilder_).build(pfCand.bestTrack()));
      try {
	TransientVertex tv = kvf.vertex(transTrksForLep1Vertex);
	if ( tv.isValid() ){
	  float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
	  if (vtxProb > bestLep1Vtx) bestLep1Vtx = vtxProb;
	}
      } catch (const std::exception& e) {}
      transTrksForLep1Vertex.pop_back();
    }
    if (lep2_doca < max_doca and lep2_doca < lep1_doca){
      // second  lepton is closer - check vertex probability
      transTrksForLep2Vertex.push_back((*theTTBuilder_).build(pfCand.bestTrack()));
      try {
	TransientVertex tv = kvf.vertex(transTrksForLep2Vertex);
	if ( tv.isValid() ){
	  float vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
	  if (vtxProb > bestLep2Vtx) bestLep2Vtx = vtxProb;
	}
      } catch (const std::exception& e) {}
      transTrksForLep2Vertex.pop_back();
    }
  }
  return max(bestLep1Vtx,bestLep2Vtx);
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
    if (prod_vtx.r() < 1e-12) return -2.0;
    return (prod_vtx - info.ll_vtx).r()/TMath::Ccgs() * info.match->mass() / info.match->p();
  }
}

KinematicFitResult 
DileptonPlusXProducer::fillDileptonInfo(pat::CompositeCandidate& dileptonCand,
					const edm::Event& iEvent,
					const bmm::Candidate& lepton1,
					const bmm::Candidate& lepton2 ) 
{
  auto kinematicLLVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  kinematicLLVertexFit.postprocess(*beamSpot_);
  
  auto displacements = compute3dDisplacement(kinematicLLVertexFit);
  addFitInfo(dileptonCand, kinematicLLVertexFit, "kin", displacements, 0, 1);
  
  if (isMC_){
    auto gen_ll = getGenMatchInfo(lepton1, lepton2);
    dileptonCand.addUserInt(  "gen_" + lepton1.name() + "1_pdgId",   gen_ll.l1_pdgId);
    dileptonCand.addUserInt(  "gen_" + lepton1.name() + "1_mpdgId",  gen_ll.l1_motherPdgId);
    dileptonCand.addUserFloat("gen_" + lepton1.name() + "1_pt",      gen_ll.l1_pt);
    dileptonCand.addUserInt(  "gen_" + lepton2.name() + "2_pdgId",   gen_ll.l2_pdgId);
    dileptonCand.addUserInt(  "gen_" + lepton2.name() + "2_mpdgId",  gen_ll.l2_motherPdgId);
    dileptonCand.addUserFloat("gen_" + lepton2.name() + "2_pt",      gen_ll.l2_pt);
    dileptonCand.addUserFloat("gen_mass",       gen_ll.ll_mass);
    dileptonCand.addUserFloat("gen_pt",         gen_ll.ll_pt);
    dileptonCand.addUserInt(  "gen_pdgId",      gen_ll.ll_pdgId);
    dileptonCand.addUserInt(  "gen_mpdgId",     gen_ll.ll_motherPdgId);
    dileptonCand.addUserInt(  "gen_cpdgId",     gen_ll.common_mother?gen_ll.common_mother->pdgId():0);
    dileptonCand.addUserFloat("gen_prod_x",     gen_ll.ll_prod_vtx.x());
    dileptonCand.addUserFloat("gen_prod_y",     gen_ll.ll_prod_vtx.y());
    dileptonCand.addUserFloat("gen_prod_z",     gen_ll.ll_prod_vtx.z());
    dileptonCand.addUserFloat("gen_vtx_x",      gen_ll.ll_vtx.x());
    dileptonCand.addUserFloat("gen_vtx_y",      gen_ll.ll_vtx.y());
    dileptonCand.addUserFloat("gen_vtx_z",      gen_ll.ll_vtx.z());
    dileptonCand.addUserFloat("gen_l3d",        (gen_ll.ll_prod_vtx-gen_ll.ll_vtx).r());
    dileptonCand.addUserFloat("gen_lxy",        (gen_ll.ll_prod_vtx-gen_ll.ll_vtx).rho());
    dileptonCand.addUserFloat("gen_tau",        computeDecayTime(gen_ll));
    if (gen_ll.match and kinematicLLVertexFit.valid()){
      dileptonCand.addUserFloat("gen_alpha_p_phi", kinematicLLVertexFit.p3().phi() - gen_ll.match->phi());
      dileptonCand.addUserFloat("gen_alpha_p_theta", kinematicLLVertexFit.p3().theta() - gen_ll.match->theta());
      TVector3 p_gen(gen_ll.match->px(),
		     gen_ll.match->py(),
		     gen_ll.match->pz());
      TVector3 ip_reco(displacements.get("pv").prodVertex().x(),
		       displacements.get("pv").prodVertex().y(),
		       displacements.get("pv").prodVertex().z());
      TVector3 ip_gen(gen_ll.ll_prod_vtx.x(),
		      gen_ll.ll_prod_vtx.y(),
		      gen_ll.ll_prod_vtx.z());
      TVector3 vtx_reco(kinematicLLVertexFit.vtx_position().x(), 
			kinematicLLVertexFit.vtx_position().y(), 
			kinematicLLVertexFit.vtx_position().z());
      TVector3 vtx_gen(gen_ll.ll_vtx.x(),
		       gen_ll.ll_vtx.y(),
		       gen_ll.ll_vtx.z());
      float cosAlpha_ip  = p_gen.Dot(vtx_gen - ip_reco) / (p_gen.Mag() * (vtx_gen - ip_reco).Mag());
      float cosAlpha_vtx = p_gen.Dot(vtx_reco - ip_gen) / (p_gen.Mag() * (vtx_reco - ip_gen).Mag());
      
      dileptonCand.addUserFloat("gen_alpha_ip", acos(cosAlpha_ip));
      dileptonCand.addUserFloat("gen_alpha_vtx", acos(cosAlpha_vtx));
    } else {
      dileptonCand.addUserFloat("gen_alpha_p_phi", 999);
      dileptonCand.addUserFloat("gen_alpha_p_theta", 999);
      dileptonCand.addUserFloat("gen_alpha_ip", 999);
      dileptonCand.addUserFloat("gen_alpha_vtx", 999);
    }
    
    double ll_doca = -1;

    if (gen_ll.gen_l1() and gen_ll.gen_l2())
      ll_doca = distanceOfClosestApproach(gen_ll.gen_l1(), gen_ll.gen_l2());
    dileptonCand.addUserFloat("gen_doca",        ll_doca);
    
    
  }

  int pvIndex = displacements.get("pv").pvIndex();

  // Look for additional tracks compatible with the dilepton vertex
  auto closeTracks = findTracksCompatibleWithTheVertex(lepton1,lepton2,kinematicLLVertexFit);
  closeTracks.fillCandInfo(dileptonCand, pvIndex, "");

  dileptonCand.addUserFloat( "m1iso",     computeTrkLeptonIsolation(lepton1,lepton2,pvIndex,0.5,0.5));
  dileptonCand.addUserFloat( "m2iso",     computeTrkLeptonIsolation(lepton2,lepton1,pvIndex,0.5,0.5));
  dileptonCand.addUserFloat( "iso",       computeTrkDileptonIsolation(lepton2,lepton1,pvIndex,0.9,0.7));
  dileptonCand.addUserFloat( "otherVtxMaxProb", otherVertexMaxProb(lepton1,lepton2,0.5));
  dileptonCand.addUserFloat( "otherVtxMaxProb1", otherVertexMaxProb(lepton1,lepton2,1.0));
  dileptonCand.addUserFloat( "otherVtxMaxProb2", otherVertexMaxProb(lepton1,lepton2,2.0));

  // BDT
  bdtData_.fls3d    = dileptonCand.userFloat("kin_sl3d");
  bdtData_.alpha    = dileptonCand.userFloat("kin_alpha");
  bdtData_.pvips    = dileptonCand.userFloat("kin_pvipErr")>0?dileptonCand.userFloat("kin_pvip")/dileptonCand.userFloat("kin_pvipErr"):999.;
  bdtData_.iso      = dileptonCand.userFloat("iso");
  bdtData_.chi2dof  = dileptonCand.userFloat("kin_vtx_chi2dof");
  bdtData_.docatrk  = dileptonCand.userFloat("docatrk");
  bdtData_.closetrk = dileptonCand.userInt(  "closetrk");
  bdtData_.m1iso    = dileptonCand.userFloat("m1iso");
  bdtData_.m2iso    = dileptonCand.userFloat("m2iso");
  bdtData_.eta      = dileptonCand.userFloat("kin_eta");	  
  bdtData_.m        = dileptonCand.userFloat("kin_mass");	  

  dileptonCand.addUserFloat("bdt",computeAnalysisBDT(iEvent.eventAuxiliary().event()%3));

  // XGBoost
  unsigned int xg_index = iEvent.eventAuxiliary().event()%3;

  xgBoosters_.at(xg_index).set("mm_kin_alpha",       dileptonCand.userFloat("kin_alpha"));
  xgBoosters_.at(xg_index).set("mm_kin_alphaXY",     cos(dileptonCand.userFloat("kin_alphaBS"))); // FIXME - need new training
  xgBoosters_.at(xg_index).set("mm_kin_spvip",       dileptonCand.userFloat("kin_spvip"));
  xgBoosters_.at(xg_index).set("mm_kin_pvip",        dileptonCand.userFloat("kin_pvip"));
  xgBoosters_.at(xg_index).set("mm_iso",             dileptonCand.userFloat("iso"));
  xgBoosters_.at(xg_index).set("mm_m1iso",           dileptonCand.userFloat("m1iso"));
  xgBoosters_.at(xg_index).set("mm_m2iso",           dileptonCand.userFloat("m2iso"));
  xgBoosters_.at(xg_index).set("mm_kin_sl3d",        dileptonCand.userFloat("kin_sl3d"));
  xgBoosters_.at(xg_index).set("mm_kin_vtx_chi2dof", dileptonCand.userFloat("kin_vtx_chi2dof"));
  xgBoosters_.at(xg_index).set("mm_nBMTrks",         dileptonCand.userInt(  "nBMTrks"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb1", dileptonCand.userFloat(  "otherVtxMaxProb1"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb2", dileptonCand.userFloat(  "otherVtxMaxProb2"));

  dileptonCand.addUserFloat("mva", xgBoosters_.at(xg_index).predict());

  // Refit with pointing constraint
  auto bToLL_PC = vertexLeptonsWithPointingConstraint(lepton1,lepton2,displacements.get("pv").prodVertex());
  addFitInfo(dileptonCand, bToLL_PC, "kinpc");
  
  return kinematicLLVertexFit;
}

void DileptonPlusXProducer::fillBtoKllInfo(pat::CompositeCandidate& btokllCand,
					   const edm::Event& iEvent,
					   const bmm::Candidate& lepton1,
					   const bmm::Candidate& lepton2,
					   const pat::PackedCandidate & kaon) 
{
  btokllCand.addUserFloat("kaon_pt",     kaon.pt());
  btokllCand.addUserFloat("kaon_eta",    kaon.eta());
  btokllCand.addUserFloat("kaon_phi",    kaon.phi());
  btokllCand.addUserFloat("kaon_dxy_bs", kaon.bestTrack()->dxy(*beamSpot_));
  btokllCand.addUserFloat("kaon_sdxy_bs", 
			  kaon.bestTrack()->dxyError()>0 ? fabs(kaon.bestTrack()->dxy(*beamSpot_))/kaon.bestTrack()->dxyError() : 0.0);
  btokllCand.addUserInt("kaon_charge", kaon.charge());
  if (isMC_){
    auto gen_kll = getGenMatchInfo(lepton1,lepton2,&kaon);
    btokllCand.addUserInt(  "gen_kaon_pdgId",  gen_kll.kaon1_pdgId);
    btokllCand.addUserInt(  "gen_kaon_mpdgId", gen_kll.kaon1_motherPdgId);
    btokllCand.addUserFloat("gen_kaon_pt",     gen_kll.kaon1_pt);
    btokllCand.addUserFloat("gen_mass",        gen_kll.kll_mass);
    btokllCand.addUserFloat("gen_pt",          gen_kll.kll_pt);
    btokllCand.addUserInt(  "gen_pdgId",       gen_kll.kll_pdgId);
    btokllCand.addUserFloat("gen_prod_x",      gen_kll.kll_prod_vtx.x());
    btokllCand.addUserFloat("gen_prod_y",      gen_kll.kll_prod_vtx.y());
    btokllCand.addUserFloat("gen_prod_z",      gen_kll.kll_prod_vtx.z());
    btokllCand.addUserFloat("gen_l3d",         (gen_kll.kll_prod_vtx - gen_kll.ll_vtx).r());
    btokllCand.addUserFloat("gen_lxy",         (gen_kll.kll_prod_vtx - gen_kll.ll_vtx).rho());
    btokllCand.addUserFloat("gen_tau",         computeDecayTime(gen_kll));
    btokllCand.addUserFloat("gen_cpdgId",      gen_kll.common_mother ? gen_kll.common_mother->pdgId() : 0);
  }

  // if (kaon.genParticle()){
  // 	btokllCand.addUserInt("kaon_mc_pdgId", kaon.genParticle().pdgId());
  // } else {
  // 	btokllCand.addUserInt("kaon_mc_pdgId", 0);
  // }

  // Inclusive
  
  auto bToKJPsiLL_NoMassConstraint = fitBToKLL(lepton1, lepton2, kaon, -1.0);
  bToKJPsiLL_NoMassConstraint.postprocess(*beamSpot_);
  auto bToKJPsiLL_NoMassConstraint_displacement = compute3dDisplacement(bToKJPsiLL_NoMassConstraint);
  addFitInfo(btokllCand, bToKJPsiLL_NoMassConstraint, "nomc", bToKJPsiLL_NoMassConstraint_displacement, -1, -1, 1);

  // JpsiK
  KinematicFitResult bToKJPsiLL_MassConstraint;
  if (fabs((lepton1.p4() + lepton2.p4()).mass()-3.1) < 0.2) {
    bToKJPsiLL_MassConstraint = fitBToKLL(lepton1, lepton2, kaon, JPsiMass_);
    bToKJPsiLL_MassConstraint.postprocess(*beamSpot_);
  }
  auto bToKJPsiLL_MassConstraint_displacement = compute3dDisplacement(bToKJPsiLL_MassConstraint);
  addFitInfo(btokllCand, bToKJPsiLL_MassConstraint, "jpsimc", bToKJPsiLL_MassConstraint_displacement,-1,-1,1);

  // Psi(2S)K
  KinematicFitResult bToKPsi2SLL_MassConstraint;
  if (fabs((lepton1.p4() + lepton2.p4()).mass()-3.7) < 0.2) {
    bToKPsi2SLL_MassConstraint = fitBToKLL(lepton1, lepton2, kaon, Psi2SMass_);
    bToKPsi2SLL_MassConstraint.postprocess(*beamSpot_);
  }
  auto bToKPsi2SLL_MassConstraint_displacement = compute3dDisplacement(bToKPsi2SLL_MassConstraint);
  addFitInfo(btokllCand, bToKPsi2SLL_MassConstraint, "psimc", bToKPsi2SLL_MassConstraint_displacement,-1,-1,1);
  
  // broken pointing constraint
  // auto bToKJPsiLL_MC_PC = refitWithPointingConstraint(bToKJPsiLL_MC.refitTree, primaryVertex);
  // bToKJPsiLL_MC_PC.postprocess(beamSpot);
  // addFitInfo(btokllCand, bToKJPsiLL_MC_PC, "mcpc");
}


void
DileptonPlusXProducer::fillDstarInfo(pat::CompositeCandidateCollection& dstar_collection,
				     const edm::Event& iEvent,
				     const KinematicFitResult& d0VertexFit,
				     const pat::CompositeCandidate& d0Cand,
				     const pat::PackedCandidate& soft_pion,
				     int mm_index,
				     int hh_index,
				     const bmm::Candidate& daughter1,
				     const bmm::Candidate& daughter2)
{
  pat::CompositeCandidate dstarCand;
  dstarCand.addUserInt("mm_index", mm_index);
  dstarCand.addUserInt("hh_index", hh_index);

  // Dstar preselection is performed using original tracks. The best
  // estimate of D0 parameters is achieved by vertexing its decay
  // products. Dstar can be prompt and non-prompt. Therefore we can
  // have the following dm estimates
  // * dm_raw - dm computed without kinematic fits
  // * dm_prompt - D0 vertexed and soft_pion is restricted to PV
  // * dm_prompt2 - same as dm_prompt, but using second best PV
  // * dm_free - D0 vertexed and soft_pion is not refit
  //
  // Dstar mass may have similiar to dm definitions, but since it's no
  // expected to be used directly, we store only one mass varian: raw

  // soft pion raw information
  dstarCand.addUserFloat("pion_pt",     soft_pion.pt());
  dstarCand.addUserFloat("pion_eta",    soft_pion.eta());
  dstarCand.addUserFloat("pion_phi",    soft_pion.phi());
  dstarCand.addUserFloat("pion_dxy_bs", soft_pion.bestTrack()->dxy(*beamSpot_));
  auto pion_sdxy_bs = 0;
  if (soft_pion.bestTrack()->dxyError() > 0)
    pion_sdxy_bs = fabs(soft_pion.bestTrack()->dxy(*beamSpot_))/soft_pion.bestTrack()->dxyError();
  dstarCand.addUserFloat("pion_sdxy_bs", pion_sdxy_bs);
  dstarCand.addUserInt("pion_charge", soft_pion.charge());

  // gen info
  // FIXME
  if (isMC_){
    auto gen_info = getGenMatchInfo(daughter1, daughter2, &soft_pion);
    dstarCand.addUserInt(  "gen_pion_pdgId",  gen_info.kaon1_pdgId);
    dstarCand.addUserInt(  "gen_pion_mpdgId", gen_info.kaon1_motherPdgId);
    dstarCand.addUserFloat("gen_pion_pt",     gen_info.kaon1_pt);
    dstarCand.addUserFloat("gen_mass",        gen_info.kll_mass);
    dstarCand.addUserFloat("gen_pt",          gen_info.kll_pt);
    dstarCand.addUserInt(  "gen_pdgId",       gen_info.kll_pdgId);
    dstarCand.addUserInt(  "gen_mpdgId",      gen_info.kll_motherPdgId);
    dstarCand.addUserFloat("gen_prod_x",      gen_info.kll_prod_vtx.x());
    dstarCand.addUserFloat("gen_prod_y",      gen_info.kll_prod_vtx.y());
    dstarCand.addUserFloat("gen_prod_z",      gen_info.kll_prod_vtx.z());
    dstarCand.addUserFloat("gen_l3d",         (gen_info.kll_prod_vtx - gen_info.ll_vtx).r());
    dstarCand.addUserFloat("gen_lxy",         (gen_info.kll_prod_vtx - gen_info.ll_vtx).rho());
    dstarCand.addUserFloat("gen_tau",         computeDecayTime(gen_info));
    dstarCand.addUserInt(  "gen_cpdgId",      gen_info.common_mother ? gen_info.common_mother->pdgId() : 0);
  }

  // store mass information
  auto d0_p4(makeLorentzVectorFromPxPyPzM(d0VertexFit.p3().x(),
					  d0VertexFit.p3().y(),
					  d0VertexFit.p3().z(),
					  d0VertexFit.mass()));
  float dstar_mass_raw = (soft_pion.p4() + daughter1.p4() + daughter2.p4()).mass();
  dstarCand.addUserFloat("mass", dstar_mass_raw);
  dstarCand.addUserFloat("dm_raw", dstar_mass_raw - (daughter1.p4() + daughter2.p4()).mass());
  dstarCand.addUserFloat("dm_free", (soft_pion.p4() + d0_p4).mass() - d0_p4.mass());
  float dm_prompt = 0;

  // refit soft_pion with PV constraint
  int pvIndex = d0Cand.userInt("kin_pvIndex");
  
  pat::PackedCandidate soft_pion_refit(soft_pion);
  double pv_prob(0), pv_with_pion_prob(0), pv_sum_pt(0), pv_sum_pt2(0);
  int pv_ntrks(0);
  if (pvIndex >= 0){
    auto fit_results = refitWithVertexConstraint(*soft_pion.bestTrack(), pvIndex);
    auto& pv_refit = fit_results.first;
    auto& pv_refit_with_soft_pion = fit_results.second;
    
    if (pv_refit_with_soft_pion.valid()){
      auto pion_index = pv_refit_with_soft_pion.number_of_daughters() - 1;
      auto soft_pion_refit_p3 = pv_refit_with_soft_pion.dau_p3(pion_index);
      auto refit_p4(makeLorentzVectorFromPxPyPzM(soft_pion_refit_p3.x(),
						 soft_pion_refit_p3.y(),
						 soft_pion_refit_p3.z(),
						 PionMass_));
      dm_prompt = (refit_p4 + d0_p4).mass() - d0_p4.mass();

      pv_with_pion_prob = pv_refit_with_soft_pion.vtxProb();
    }

    if (pv_refit.valid()){
      pv_prob = pv_refit.vtxProb();
      pv_sum_pt = pv_refit.sumPt();
      pv_sum_pt2 = pv_refit.sumPt2();
      pv_ntrks = pv_refit.number_of_daughters();
    }
  }
  
  dstarCand.addUserFloat("dm_prompt", dm_prompt);
  dstarCand.addUserFloat("pv_prob", pv_prob);
  dstarCand.addUserFloat("pv_sum_pt", pv_sum_pt);
  dstarCand.addUserFloat("pv_sum_pt2", pv_sum_pt2);
  dstarCand.addUserFloat("pv_with_pion_prob", pv_with_pion_prob);
  dstarCand.addUserInt("pv_ntrks", pv_ntrks);
  
  dstar_collection.push_back(dstarCand);
}

void DileptonPlusXProducer::fillBtoKKllInfo(pat::CompositeCandidate& bCand,
					    const edm::Event& iEvent,
					    const bmm::Candidate& lepton1,
					    const bmm::Candidate& lepton2,
					    const pat::PackedCandidate & kaon1,
					    const pat::PackedCandidate & kaon2) 
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
    auto gen_info = getGenMatchInfo(lepton1,lepton2,&kaon1,&kaon2);
    bCand.addUserInt(  "gen_kaon1_pdgId",  gen_info.kaon1_pdgId);
    bCand.addUserInt(  "gen_kaon1_mpdgId", gen_info.kaon1_motherPdgId);
    bCand.addUserFloat("gen_kaon1_pt",     gen_info.kaon1_pt);
    bCand.addUserInt(  "gen_kaon2_pdgId",  gen_info.kaon2_pdgId);
    bCand.addUserInt(  "gen_kaon2_mpdgId", gen_info.kaon2_motherPdgId);
    bCand.addUserFloat("gen_kaon2_pt",     gen_info.kaon2_pt);
    bCand.addUserFloat("gen_mass",         gen_info.kkll_mass);
    bCand.addUserFloat("gen_pt",           gen_info.kkll_pt);
    bCand.addUserInt(  "gen_pdgId",        gen_info.kkll_pdgId);
    bCand.addUserFloat("gen_prod_x",       gen_info.kkll_prod_vtx.x());
    bCand.addUserFloat("gen_prod_y",       gen_info.kkll_prod_vtx.y());
    bCand.addUserFloat("gen_prod_z",       gen_info.kkll_prod_vtx.z());
    bCand.addUserFloat("gen_l3d",         (gen_info.kkll_prod_vtx - gen_info.ll_vtx).r());
    bCand.addUserFloat("gen_lxy",         (gen_info.kkll_prod_vtx - gen_info.ll_vtx).rho());
    bCand.addUserFloat("gen_tau",          computeDecayTime(gen_info));
    bCand.addUserFloat("gen_cpdgId",       gen_info.common_mother ? gen_info.common_mother->pdgId() : 0);
  }

  // Inclusive
  auto bToKKll = fitBToKKLL(lepton1, lepton2, kaon1, kaon2);
  bToKKll.postprocess(*beamSpot_);
  auto bToKKll_displacement = compute3dDisplacement(bToKKll);
  addFitInfo(bCand, bToKKll, "kin", bToKKll_displacement,-1,-1,1,2);
  bCand.addUserFloat("kin_kk_mass",    bToKKll.refit_mass(1,2));

  // Jpsi KK
  
  KinematicFitResult bToJpsiKK;
  if (fabs((lepton1.p4() + lepton2.p4()).mass() - JPsiMass_) < 0.2) { 
    bToJpsiKK = fitBToKKLL(lepton1, lepton2, kaon1, kaon2, JPsiMass_);
    bToJpsiKK.postprocess(*beamSpot_);
  }
  auto bToJpsiKK_displacement = compute3dDisplacement(bToJpsiKK);
  addFitInfo(bCand, bToJpsiKK, "jpsikk", bToJpsiKK_displacement, -1, -1, 1, 2);
  bCand.addUserFloat("jpsikk_kk_mass",    bToJpsiKK.refit_mass(1,2));

  // Phi ll
  
  KinematicFitResult bToPhill;
  if (fabs((kaon1.p4() + kaon2.p4()).mass() - PhiMass_) < 0.01) { 
    bToPhill = fitBToKKLL(lepton1, lepton2, kaon1, kaon2, -1, PhiMass_);
    bToPhill.postprocess(*beamSpot_);
  }
  auto bToPhill_displacement = compute3dDisplacement(bToPhill);
  addFitInfo(bCand, bToPhill, "phill", bToPhill_displacement,-1,-1,1,2);
  
}

void DileptonPlusXProducer::fillLLGammaGenInfo(pat::CompositeCandidate& llgCand,
					       const edm::Event& iEvent,
					       const bmm::Candidate& lepton1,
					       const bmm::Candidate& lepton2,
					       const reco::Candidate & photon) 
{
  if (isMC_){
    auto gen_llg = getGenMatchInfo(lepton1,lepton2,nullptr,nullptr,&photon);
    llgCand.addUserInt(  "gen_ph_pdgId",    gen_llg.photon_pdgId);
    llgCand.addUserInt(  "gen_ph_mpdgId",   gen_llg.photon_motherPdgId);
    llgCand.addUserFloat("gen_ph_pt",       gen_llg.photon_pt);
    llgCand.addUserFloat("gen_mass",        gen_llg.llg_mass);
    llgCand.addUserFloat("gen_pt",          gen_llg.llg_pt);
    llgCand.addUserInt(  "gen_pdgId",       gen_llg.llg_pdgId);
    llgCand.addUserFloat("gen_prod_x",      gen_llg.llg_prod_vtx.x());
    llgCand.addUserFloat("gen_prod_y",      gen_llg.llg_prod_vtx.y());
    llgCand.addUserFloat("gen_prod_z",      gen_llg.llg_prod_vtx.z());
    llgCand.addUserFloat("gen_l3d",         (gen_llg.llg_prod_vtx-gen_llg.ll_vtx).r());
    llgCand.addUserFloat("gen_lxy",         (gen_llg.llg_prod_vtx-gen_llg.ll_vtx).rho());
    llgCand.addUserFloat("gen_tau",         computeDecayTime(gen_llg));
    llgCand.addUserFloat("gen_cpdgId",      gen_llg.common_mother?gen_llg.common_mother->pdgId():0);
  }
}

void DileptonPlusXProducer::fillLLGammaInfo(pat::CompositeCandidate& llgCand,
					    const edm::Event& iEvent,
					    const bmm::Candidate& lepton1,
					    const bmm::Candidate& lepton2,
					    const pat::Photon & photon) 
{
  // Rebuild ll vertex to ensure that the KinematicTree remains self
  // consistent and no elements get out of scope or get deleted
  // when the tree is used in subsequent fits
  auto llVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  
  llVertexFit.tree()->movePointerToTheTop();
  
  llgCand.addUserFloat("ph_pt",     photon.pt());
  llgCand.addUserFloat("ph_eta",    photon.eta());
  llgCand.addUserFloat("ph_phi",    photon.phi());
  // FIXME: add photon id

  fillLLGammaGenInfo(llgCand, iEvent, lepton1, lepton2, photon);

  KinematicFitResult jpsiGamma_NoMassConstraint;
  jpsiGamma_NoMassConstraint = fitLLGamma(llVertexFit.tree(), photon);
  if (lepton1.track()) 
    jpsiGamma_NoMassConstraint.tracks.push_back(lepton1.track());
  if (lepton2.track()) 
    jpsiGamma_NoMassConstraint.tracks.push_back(lepton2.track());
  jpsiGamma_NoMassConstraint.postprocess(*beamSpot_);
  auto jpsiGamma_NoMassConstraint_displacement = compute3dDisplacement(jpsiGamma_NoMassConstraint);
  addFitInfo(llgCand, jpsiGamma_NoMassConstraint, "nomc", jpsiGamma_NoMassConstraint_displacement, -1, -1, 1);

  KinematicFitResult jpsiGamma_MassConstraint;
  if (fabs(llVertexFit.mass()-3.1) < 0.2) {
    jpsiGamma_MassConstraint = fitLLGamma(llVertexFit.tree(), photon, JPsiMass_);
    if (lepton1.track()) 
      jpsiGamma_MassConstraint.tracks.push_back(lepton1.track());
    if (lepton2.track()) 
      jpsiGamma_MassConstraint.tracks.push_back(lepton2.track());
    jpsiGamma_MassConstraint.postprocess(*beamSpot_);
  }
  auto jpsiGamma_MassConstraint_displacement = compute3dDisplacement(jpsiGamma_MassConstraint);
  addFitInfo(llgCand, jpsiGamma_MassConstraint, "jpsimc", jpsiGamma_MassConstraint_displacement, -1, -1, 1);

  // need bestVertex

  // FIXME add relevant info that normally is added by addFitInfo
  // - not all of it is needed, because everything vertex related should be
  //   extracted from mm vertex
  // - most of the variables are etracted using bmm::Displacements info
  //   computed with compute3dDisplacement
  // - massErr - this is something non-trivial
  // - pt,eta,phi
  
  // auto alpha = getAlpha(fit.refitVertex->vertexState().position(),
  // 			fit.refitVertex->vertexState().error(),
  // 			GlobalPoint(Basic3DVector<float>(bestVertex->position())),
  // 			GlobalError(bestVertex->covariance()),
  // 			fit.refitMother->currentState().globalMomentum());

  // auto alphaXY = getAlpha(fit.refitVertex->vertexState().position(),
  // 			  fit.refitVertex->vertexState().error(),
  // 			  GlobalPoint(Basic3DVector<float>(bestVertex->position())),
  // 			  GlobalError(bestVertex->covariance()),
  // 			  fit.refitMother->currentState().globalMomentum(),
  // 			  true);

  // result.alpha    = alpha.first;
  // result.alphaErr = alpha.second;

  // result.alphaXY    = alphaXY.first;
  // result.alphaXYErr = alphaXY.second;

}

void DileptonPlusXProducer::fillLLGammaConvInfo(pat::CompositeCandidate& llgCand,
						const edm::Event& iEvent,
						const bmm::Candidate& lepton1,
						const bmm::Candidate& lepton2,
						const pat::CompositeCandidate& photon) 
{
  // Rebuild ll vertex to ensure that the KinematicTree remains self
  // consistent and no elements get out of scope or get deleted
  // when the tree is used in subsequent fits
  auto llVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  
  llVertexFit.tree()->movePointerToTheTop();
  
  llgCand.addUserFloat("ph_pt",     photon.pt());
  llgCand.addUserFloat("ph_eta",    photon.eta());
  llgCand.addUserFloat("ph_phi",    photon.phi());

  // FIXME: add photon id

  fillLLGammaGenInfo(llgCand, iEvent, lepton1, lepton2, photon);

  // Kinematic fits
  
  KinematicFitResult jpsiGamma_NoMassConstraint;
  jpsiGamma_NoMassConstraint = fitLLGammaConv(llVertexFit.tree(), photon);
  if (lepton1.track()) 
    jpsiGamma_NoMassConstraint.tracks.push_back(lepton1.track());
  if (lepton2.track()) 
    jpsiGamma_NoMassConstraint.tracks.push_back(lepton2.track());
  jpsiGamma_NoMassConstraint.postprocess(*beamSpot_);
  auto jpsiGamma_NoMassConstraint_displacement = compute3dDisplacement(jpsiGamma_NoMassConstraint);
  addFitInfo(llgCand, jpsiGamma_NoMassConstraint, "nomc", jpsiGamma_NoMassConstraint_displacement, -1, -1, 1);

  KinematicFitResult jpsiGamma_MassConstraint;
  if (fabs(llVertexFit.mass()-3.1) < 0.2) {
    jpsiGamma_MassConstraint = fitLLGammaConv(llVertexFit.tree(), photon, JPsiMass_);
    if (lepton1.track())
      jpsiGamma_MassConstraint.tracks.push_back(lepton1.track());
    if (lepton2.track())
      jpsiGamma_MassConstraint.tracks.push_back(lepton2.track());
    jpsiGamma_MassConstraint.postprocess(*beamSpot_);
  }
  auto jpsiGamma_MassConstraint_displacement = compute3dDisplacement(jpsiGamma_MassConstraint);
  addFitInfo(llgCand, jpsiGamma_MassConstraint, "jpsimc", jpsiGamma_MassConstraint_displacement, -1, -1, 1);

}


void 
DileptonPlusXProducer::fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(pat::CompositeCandidate& llK,
								    const pat::CompositeCandidate& ll,
								    const edm::Event& iEvent,
								    const KinematicFitResult& kinematicLLVertexFit,
								    const bmm::Candidate& lepton1,
								    const bmm::Candidate& lepton2,
								    const pat::PackedCandidate & kaon) 
{
  ///////
  //// Treat B->JpsiK as B->mm for signal studies in data
  //
  std::vector<const pat::PackedCandidate*> ignoreTracks;
  ignoreTracks.push_back(&kaon);

  int pvIndex = llK.userInt("jpsimc_pvIndex");

  // Look for additional tracks compatible with the dimuon vertex
  auto closeTracks = findTracksCompatibleWithTheVertex(lepton1,lepton2,kinematicLLVertexFit,0.03,ignoreTracks);
  closeTracks.fillCandInfo(llK, pvIndex, "bmm");

  llK.addUserFloat( "bmm_m1iso",     computeTrkLeptonIsolation(lepton1,lepton2,pvIndex,0.5,0.5,ignoreTracks));
  llK.addUserFloat( "bmm_m2iso",     computeTrkLeptonIsolation(lepton2,lepton1,pvIndex,0.5,0.5,ignoreTracks));
  llK.addUserFloat( "bmm_iso",       computeTrkDileptonIsolation(lepton2,lepton1,pvIndex,0.9,0.7,ignoreTracks));
  llK.addUserFloat( "bmm_otherVtxMaxProb",  otherVertexMaxProb(lepton1,lepton2,0.5,0.1,ignoreTracks));
  llK.addUserFloat( "bmm_otherVtxMaxProb1", otherVertexMaxProb(lepton1,lepton2,1.0,0.1,ignoreTracks));
  llK.addUserFloat( "bmm_otherVtxMaxProb2", otherVertexMaxProb(lepton1,lepton2,2.0,0.1,ignoreTracks));

  // BDT
  bdtData_.fls3d    = ll.userFloat("kin_sl3d");
  bdtData_.alpha    = llK.userFloat("jpsimc_alpha");
  bdtData_.pvips    = llK.userFloat("jpsimc_pvipErr")>0?llK.userFloat("jpsimc_pvip")/llK.userFloat("jpsimc_pvipErr"):999;
  // One can use bkmm without mass constraint, but it doesn't help
  // bdtData_.alpha    = mmK.userFloat("nomc_alpha");
  // bdtData_.pvips    = mmK.userFloat("nomc_pvip")/mmK.userFloat("nomc_pvipErr");
  bdtData_.iso      = llK.userFloat("bmm_iso");
  bdtData_.chi2dof  = ll.userFloat("kin_vtx_chi2dof");
  bdtData_.docatrk  = llK.userFloat("bmm_docatrk");
  bdtData_.closetrk = llK.userInt(  "bmm_closetrk");
  bdtData_.m1iso    = llK.userFloat("bmm_m1iso");
  bdtData_.m2iso    = llK.userFloat("bmm_m2iso");
  bdtData_.eta      = llK.userFloat("jpsimc_eta");	  
  bdtData_.m        = llK.userFloat("jpsimc_mass");	  

  llK.addUserFloat("bmm_bdt",computeAnalysisBDT(iEvent.eventAuxiliary().event()%3));

  // XGBoost
  unsigned int xg_index = iEvent.eventAuxiliary().event()%3;
  // Pointing angle - mmK
  xgBoosters_.at(xg_index).set("mm_kin_alpha",       llK.userFloat("jpsimc_alpha"));
  xgBoosters_.at(xg_index).set("mm_kin_alphaXY",     cos(llK.userFloat("jpsimc_alphaBS"))); // FIXME - need new training
  // PV matching - mmK
  xgBoosters_.at(xg_index).set("mm_kin_spvip",       llK.userFloat("jpsimc_spvip"));
  xgBoosters_.at(xg_index).set("mm_kin_pvip",        llK.userFloat("jpsimc_pvip"));
  // Isolation and extra track variables need to be recomputed ignoring kaon
  xgBoosters_.at(xg_index).set("mm_iso",             llK.userFloat("bmm_iso"));
  xgBoosters_.at(xg_index).set("mm_m1iso",           llK.userFloat("bmm_m1iso"));
  xgBoosters_.at(xg_index).set("mm_m2iso",           llK.userFloat("bmm_m2iso"));
  xgBoosters_.at(xg_index).set("mm_nBMTrks",         llK.userInt(  "bmm_nBMTrks"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb1", llK.userFloat(  "bmm_otherVtxMaxProb1"));
  xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb2", llK.userFloat(  "bmm_otherVtxMaxProb2"));
  // Vertexing - mm
  xgBoosters_.at(xg_index).set("mm_kin_vtx_chi2dof", ll.userFloat("kin_vtx_chi2dof"));
  // Flight length significance - mm
  xgBoosters_.at(xg_index).set("mm_kin_sl3d",        ll.userFloat("kin_sl3d")*1.6);

  llK.addUserFloat("bmm_mva", xgBoosters_.at(xg_index).predict());
}


// // FIXME: need to add info similar to addFitInfo for mmg
void 
DileptonPlusXProducer::fillMvaInfoForLLGamma(pat::CompositeCandidate& llg,
					     const pat::CompositeCandidate& ll,
					     const edm::Event& iEvent,
					     const KinematicFitResult& kinematicLLVertexFit,
					     const bmm::Candidate& lepton1,
					     const bmm::Candidate& lepton2,
					     const pat::Photon & photon) 
{
  int pvIndex = ll.userInt("kin_pvIndex");

  // Look for additional tracks compatible with the dimuon vertex
  auto closeTracks = findTracksCompatibleWithTheVertex(lepton1, lepton2, kinematicLLVertexFit, 0.03);
  closeTracks.fillCandInfo(llg, pvIndex, "");

  llg.addUserFloat( "_m1iso",     computeTrkLeptonIsolation(lepton1, lepton2, pvIndex, 0.5, 0.5));
  llg.addUserFloat( "_m2iso",     computeTrkLeptonIsolation(lepton2, lepton1, pvIndex, 0.5, 0.5));
  llg.addUserFloat( "_iso",       computeTrkDileptonIsolation(lepton2, lepton1, pvIndex, 0.9, 0.7));
  llg.addUserFloat( "_otherVtxMaxProb",  otherVertexMaxProb(lepton1, lepton2, 0.5, 0.1));
  llg.addUserFloat( "_otherVtxMaxProb1", otherVertexMaxProb(lepton1, lepton2, 1.0, 0.1));
  llg.addUserFloat( "_otherVtxMaxProb2", otherVertexMaxProb(lepton1, lepton2, 2.0, 0.1));

  // // BDT
  // bdtData_.fls3d    = mm.userFloat("kin_sl3d");
  // bdtData_.alpha    = mmK.userFloat("jpsimc_alpha");
  // bdtData_.pvips    = mmK.userFloat("jpsimc_pvipErr")>0?mmK.userFloat("jpsimc_pvip")/mmK.userFloat("jpsimc_pvipErr"):999;
  // // One can use bkmm without mass constraint, but it doesn't help
  // // bdtData_.alpha    = mmK.userFloat("nomc_alpha");
  // // bdtData_.pvips    = mmK.userFloat("nomc_pvip")/mmK.userFloat("nomc_pvipErr");
  // bdtData_.iso      = mmK.userFloat("bmm_iso");
  // bdtData_.chi2dof  = mm.userFloat("kin_vtx_chi2dof");
  // bdtData_.docatrk  = mmK.userFloat("bmm_docatrk");
  // bdtData_.closetrk = mmK.userInt(  "bmm_closetrk");
  // bdtData_.m1iso    = mmK.userFloat("bmm_m1iso");
  // bdtData_.m2iso    = mmK.userFloat("bmm_m2iso");
  // bdtData_.eta      = mmK.userFloat("jpsimc_eta");	  
  // bdtData_.m        = mmK.userFloat("jpsimc_mass");	  

  //   mmK.addUserFloat("bmm_bdt",computeAnalysisBDT(iEvent.eventAuxiliary().event()%3));

  //   // XGBoost
  //   unsigned int xg_index = iEvent.eventAuxiliary().event()%3;
  //   // Pointing angle - mmK
  //   xgBoosters_.at(xg_index).set("mm_kin_alpha",       mmK.userFloat("jpsimc_alpha"));
  //   xgBoosters_.at(xg_index).set("mm_kin_alphaXY",     cos(mmK.userFloat("jpsimc_alphaBS"))); // FIXME - need new training
  //   // PV matching - mmK
  //   xgBoosters_.at(xg_index).set("mm_kin_spvip",       mmK.userFloat("jpsimc_spvip"));
  //   xgBoosters_.at(xg_index).set("mm_kin_pvip",        mmK.userFloat("jpsimc_pvip"));
  //   // Isolation and extra track variables need to be recomputed ignoring kaon
  //   xgBoosters_.at(xg_index).set("mm_iso",             mmK.userFloat("bmm_iso"));
  //   xgBoosters_.at(xg_index).set("mm_m1iso",           mmK.userFloat("bmm_m1iso"));
  //   xgBoosters_.at(xg_index).set("mm_m2iso",           mmK.userFloat("bmm_m2iso"));
  //   xgBoosters_.at(xg_index).set("mm_nBMTrks",         mmK.userInt(  "bmm_nBMTrks"));
  //   xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb1", mmK.userFloat(  "bmm_otherVtxMaxProb1"));
  //   xgBoosters_.at(xg_index).set("mm_otherVtxMaxProb2", mmK.userFloat(  "bmm_otherVtxMaxProb2"));
  //   // Vertexing - mm
  //   xgBoosters_.at(xg_index).set("mm_kin_vtx_chi2dof", mm.userFloat("kin_vtx_chi2dof"));
  //   // Flight length significance - mm
  //   xgBoosters_.at(xg_index).set("mm_kin_sl3d",        mm.userFloat("kin_sl3d")*1.6);

  //   mmK.addUserFloat("bmm_mva", xgBoosters_.at(xg_index).predict());
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
DileptonPlusXProducer::injectHadronsThatMayFakeMuonsInMC(std::vector<bmm::Candidate>& good_lepton_candidates){
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
    std::vector<const pat::PackedGenParticle*> final_state_particles;

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
	for (const auto& good_muon_candidate: good_lepton_candidates){
	  if (deltaR(*hadron, good_muon_candidate) < 0.01) {
	    good_candidate=false;
	    break;
	  }
	}
	if (good_candidate)
	  good_lepton_candidates.push_back(bmm::Candidate(pfCand, hadron));
      }
    }
  }
}

void 
DileptonPlusXProducer::injectBhhHadrons(std::vector<bmm::Candidate>& good_lepton_candidates){

  // find reco tracks matching preselection requirements
  for (const auto& pfCand: *pfCandHandle_.product()){
    if (not isGoodTrack(pfCand)) continue; 
    if (pfCand.pt() < minBhhTrkPt_) continue;
    if (abs(pfCand.eta()) > maxBhhTrkEta_) continue;
    good_lepton_candidates.push_back(bmm::Candidate(pfCand));
  }
}

// void 
// DileptonPlusXProducer::injectDhhHadrons(std::vector<bmm::Candidate>& good_lepton_candidates){
//   std::vector<const pat::PackedCandidate*> selected_candidates;

//   // find reco tracks matching preselection requirements
//   for (const auto& pfCand: *pfCandHandle_.product()){
//     if (pfCand.charge() == 0) continue;
//     if (pfCand.pt() < minDhhTrkPt_) continue;
//     if (abs(pfCand.eta()) > maxDhhTrkEta_) continue;
//     if (not pfCand.hasTrackDetails()) continue;
//     if (not pfCand.bestTrack()->quality(reco::Track::highPurity)) continue;
//     selected_candidates.push_back(&pfCand);
//   }
  
//   if (selected_candidates.size() < 1) return;
//   for (unsigned int i=0; i < selected_candidates.size() - 1; ++i){
//     auto kaon1 bmm::Candidate(*selected_candidates.at(i));
//     kaon1.setType(KaonMass_, "had");
//     auto pion1 bmm::Candidate(*selected_candidates.at(i));
//     pion1.setType(PionMass_, "had");
//     for (unsigned int j=i + 1; j < selected_candidates.size(); ++j){
//       auto kaon2 bmm::Candidate(*selected_candidates.at(j));
//       kaon2.setType(KaonMass_, "had");
//       auto pion2 bmm::Candidate(*selected_candidates.at(j));
//       pion2.setType(PionMass_, "had");
      



//   // good_lepton_candidates.push_back(bmm::Candidate(pfCand));


// }


void 
DileptonPlusXProducer::injectJpsiTracks(std::vector<bmm::Candidate>& good_lepton_candidates){
  // Look for Jpsi to mu mu candidates even if muon is too soft to be
  // constructed as a muon

  // find reco tracks matching preselection requirements
  const auto& pf_candidates = *pfCandHandle_.product();
    
  std::vector<bmm::Candidate> muon_candidates;
  muon_candidates.assign(good_lepton_candidates.begin(), good_lepton_candidates.end());

  for (unsigned int i=0; i < pf_candidates.size(); ++i){
    if (not isGoodMuonCandidateFromTrack(pf_candidates.at(i))) continue;
    bool new_mu = true;
    for (const auto& lep: muon_candidates){
      if (overlap(lep, pf_candidates.at(i))) new_mu = false;
    }
    if (not new_mu) continue;
    muon_candidates.push_back(bmm::Candidate(pf_candidates.at(i)));
    muon_candidates.back().setType(MuonMass_, "mu");
  }
  
  if (muon_candidates.size() - good_lepton_candidates.size() < 2 ) return;

  for (unsigned int i=0; i < muon_candidates.size()-1; ++i){
    for (unsigned int j=i+1; j < muon_candidates.size(); ++j){
      if (muon_candidates.at(i).charge() == muon_candidates.at(j).charge()) continue;
      if (fabs((muon_candidates.at(i).p4() + muon_candidates.at(j).p4()).mass() - 3.1) > 0.2) continue;
      auto ll_doca = distanceOfClosestApproach(muon_candidates.at(i).track(), muon_candidates.at(j).track());
      if (maxTwoTrackDOCA_>0 and ll_doca > maxTwoTrackDOCA_) continue;
      bool new_mu1 = true;
      bool new_mu2 = true;
      for (const auto& lep: good_lepton_candidates){
	if (overlap(lep, muon_candidates.at(i))) new_mu1 = false;
	if (overlap(lep, muon_candidates.at(j))) new_mu2 = false;
      }

      if (new_mu1)
	good_lepton_candidates.push_back(muon_candidates.at(i));
      if (new_mu2)
	good_lepton_candidates.push_back(muon_candidates.at(j));
    }
  }
}


// namespace {
//   void dump_track(const reco::Track* trk){
//     std::cout << "pt(): " << trk->pt() << std::endl;
//     std::cout << "eta(): " << trk->eta() << std::endl;
//     std::cout << "charge(): " << trk->charge() << std::endl;
//     std::cout << "chi2(): " << trk->chi2() << std::endl;
//     std::cout << "ndof(): " << trk->ndof() << std::endl;
//     std::cout << "vertex(): " << trk->vertex().x() << ", " << trk->vertex().y() << ", " << trk->vertex().z() << std::endl;
//     std::cout << "momentum(): " << trk->momentum().x() << ", " << trk->momentum().y() << ", " << trk->momentum().z() << std::endl;
//     std::cout << "covariance matrix: " << std::endl;
//     // for (unsigned int i=0; i<5; ++i)
//     //   for (unsigned int j=0; j<5; ++j)
//     //     std::cout << "\t(" << i << "," << j << "): " << trk->covariance(i,j) << std::endl;
//     for (unsigned int i=0; i<15; ++i)
//       std::cout << "\t" << trk->covariance().Array()[i] << ", ";
//     std::cout << std::endl;
//   }

//   void dump_photon(const reco::Photon* photon){
//     std::cout << "pt(): " << photon->pt() << std::endl;
//     std::cout << "eta(): " << photon->eta() << std::endl;
//     std::cout << "caloPosition: " << photon->caloPosition().x() << ", " << photon->caloPosition().y() << ", " << photon->caloPosition().z() << std::endl;
//     std::cout << "XYZT: " << photon->px() << ", " << photon->py() << ", " << photon->pz() << ", " << photon->energy() << std::endl;
//   }
// }



bool DileptonPlusXProducer::preprocess(pat::CompositeCandidate& candidate,
				       const edm::Event& iEvent,
				       const bmm::Candidate& cand1,
				       const bmm::Candidate& cand2)
{
  auto ll_doca = distanceOfClosestApproach(cand1.track(), cand2.track());
  
  if (maxTwoTrackDOCA_>0 and ll_doca > maxTwoTrackDOCA_) return false;
  if (diLeptonCharge_ && cand1.charge() * cand2.charge() > 0) return false;

  // We have a decent candidate. Let's fill some information
  candidate.addUserInt(   cand1.name() + "1_index", cand1.index());
  candidate.addUserInt(   cand1.name() + "1_pdgId", cand1.pdgId());
  candidate.addUserFloat( cand1.name() + "1_pt",    cand1.pt());
  candidate.addUserFloat( cand1.name() + "1_eta",   cand1.eta());
  candidate.addUserFloat( cand1.name() + "1_phi",   cand1.phi());
  candidate.addUserInt(   cand2.name() + "2_index", cand2.index());
  candidate.addUserInt(   cand2.name() + "2_pdgId", cand2.pdgId());
  candidate.addUserFloat( cand2.name() + "2_pt",    cand2.pt());
  candidate.addUserFloat( cand2.name() + "2_eta",   cand2.eta());
  candidate.addUserFloat( cand2.name() + "2_phi",   cand2.phi());
  candidate.addUserFloat( "doca",                     ll_doca);

  return true;
}

void
DileptonPlusXProducer::buildLLXCandidates(pat::CompositeCandidateCollection& llk,
					  pat::CompositeCandidateCollection& llkk,
					  const edm::Event& iEvent,
					  const KinematicFitResult& kinematicLLVertexFit,
					  const pat::CompositeCandidate& dileptonCand,
					  int ll_index,
					  const bmm::Candidate& lepton1,
					  const bmm::Candidate& lepton2)
{
  auto nPFCands = pfCandHandle_->size();
  for (unsigned int k = 0; k < nPFCands; ++k) {
    pat::PackedCandidate kaonCand1((*pfCandHandle_)[k]);
    kaonCand1.setMass(KaonMass_);
    if (not isGoodTrack(kaonCand1)) continue; 
    if (abs(kaonCand1.pdgId()) != 211) continue;
    if (kaonCand1.pt() < ptMinKaon_ or abs(kaonCand1.eta()) > etaMaxKaon_) continue;
    if (overlap(lepton1, kaonCand1) || overlap(lepton2, kaonCand1)) continue;
    double l1_kaon_doca = distanceOfClosestApproach(lepton1.track(),
						    kaonCand1.bestTrack());
    double l2_kaon_doca = distanceOfClosestApproach(lepton2.track(),
						    kaonCand1.bestTrack());
    if (maxTwoTrackDOCA_ > 0 and l1_kaon_doca > maxTwoTrackDOCA_) continue;
    if (maxTwoTrackDOCA_ > 0 and l2_kaon_doca > maxTwoTrackDOCA_) continue;
	      
    bool goodBtoLLK = true;

    double kll_mass = (lepton1.p4() + lepton2.p4() + kaonCand1.p4()).mass();
    if (kll_mass < minBKllMass_ || kll_mass > maxBKllMass_) goodBtoLLK = false;
	    
    // fill BtoLLK candidate info
    if (goodBtoLLK){
      pat::CompositeCandidate btokllCand;
      btokllCand.addUserInt(dileptonCand.name() + "_index", ll_index);
      btokllCand.addUserFloat("kaon_l1_doca", l1_kaon_doca);
      btokllCand.addUserFloat("kaon_l2_doca", l2_kaon_doca);
      
      fillBtoKllInfo(btokllCand, iEvent, lepton1, lepton2, kaonCand1);
      fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokllCand, dileptonCand, iEvent, 
						   kinematicLLVertexFit, lepton1, lepton2, kaonCand1);

      llk.push_back(btokllCand);
      
    }

    // Build BsToKKll
    for (unsigned int k2 = k+1; k2 < nPFCands; ++k2) { 
      // only works if selection requirements for both kaons are identical
      pat::PackedCandidate kaonCand2((*pfCandHandle_)[k2]);
      kaonCand2.setMass(KaonMass_);
      if (not isGoodTrack(kaonCand2)) continue; 
      if (abs(kaonCand2.pdgId()) != 211) continue;
      if (kaonCand2.pt() < ptMinKaon_ || abs(kaonCand2.eta()) > etaMaxKaon_) continue;
      if (overlap(lepton1, kaonCand2) || overlap(lepton2, kaonCand2)) continue;
      double l1_kaon2_doca = distanceOfClosestApproach(lepton1.track(),
						       kaonCand2.bestTrack());
      double l2_kaon2_doca = distanceOfClosestApproach(lepton2.track(),
						       kaonCand2.bestTrack());
      if (maxTwoTrackDOCA_>0 and l1_kaon2_doca > maxTwoTrackDOCA_) continue;
      if (maxTwoTrackDOCA_>0 and l2_kaon2_doca > maxTwoTrackDOCA_) continue;
      
      bool goodBtoLLKK = true;
      
      double kkll_mass = (lepton1.p4() + lepton2.p4() + kaonCand1.p4() + kaonCand2.p4()).mass();
      if ( kkll_mass < minBKKllMass_ || kkll_mass > maxBKKllMass_ ) goodBtoLLKK = false;
		  
      // fill BtoLLKK candidate info
      if (goodBtoLLKK){
	pat::CompositeCandidate btokkllCand;
	btokkllCand.addUserInt(dileptonCand.name() + "_index", ll_index);
	btokkllCand.addUserFloat("kaon1_l1_doca", l1_kaon_doca);
	btokkllCand.addUserFloat("kaon1_l2_doca", l2_kaon_doca);
	btokkllCand.addUserFloat("kaon2_l1_doca", l1_kaon2_doca);
	btokkllCand.addUserFloat("kaon2_l2_doca", l2_kaon2_doca);
	
	fillBtoKKllInfo(btokkllCand, iEvent, lepton1, lepton2, kaonCand1, kaonCand2);
	// FIXME
	// fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokkllCand,dileptonCand,iEvent,kinematicLLVertexFit,lepton1,lepton2,kaonCand1);

	llkk.push_back(btokkllCand);
	
      }
    }
  }
}

void
DileptonPlusXProducer::buildDstarCandidates(pat::CompositeCandidateCollection& dstar_collection,
					    pat::CompositeCandidateCollection& hh_collection,
					    const edm::Event& iEvent,
					    const pat::PackedCandidate& had1,
					    const pat::PackedCandidate& had2)
{
  if (had1.pt() < minDhhTrkPt_ || fabs(had1.eta()) > maxDhhTrkEta_) return;
  if (had2.pt() < minDhhTrkPt_ || fabs(had2.eta()) > maxDhhTrkEta_) return;
  AddFourMomenta addP4;
  auto nPFCands = pfCandHandle_->size();
  for (unsigned int k = 0; k < nPFCands; ++k) {
    pat::PackedCandidate soft_pion((*pfCandHandle_)[k]);
    if ( not isGoodHadron(soft_pion) ) continue;
    
    if (overlap(had1, soft_pion) || overlap(had2, soft_pion)) continue;

    soft_pion.setMass(PionMass_);

    bmm::Candidate pion1(had1);
    pion1.setType(PionMass_, "had", 211 * had1.charge());
    bmm::Candidate pion2(had2);
    pion2.setType(PionMass_, "had", 211 * had2.charge());

    bmm::Candidate kaon1(had1);
    kaon1.setType(KaonMass_, "had", 321 * had1.charge());
    bmm::Candidate kaon2(had2);
    kaon2.setType(KaonMass_, "had", 321 * had2.charge());

    // D0->pipi
    if (recoD0pipi_){
      double d0_mass = (pion1.p4() + pion2.p4()).mass();
      double dstar_mass = (pion1.p4() + pion2.p4() + soft_pion.p4()).mass();

      if (d0_mass > minD0Mass_ && d0_mass < maxD0Mass_ &&
	  (dstar_mass - d0_mass) > min_dm_ && (dstar_mass - d0_mass) < max_dm_){
	
	pat::CompositeCandidate d0Cand(std::string("hh"));
	d0Cand.addDaughter( pion1 , "pion1");
	d0Cand.addDaughter( pion2 , "pion2");
	addP4.set( d0Cand );
	
	if (preprocess(d0Cand, iEvent, pion1, pion2)){
	  // Kinematic Fits
	  auto d0VertexFit = fillDileptonInfo(d0Cand, iEvent, pion1, pion2);
	  int hh_index = hh_collection.size();
	  hh_collection.push_back(d0Cand);
	  fillDstarInfo(dstar_collection, iEvent, d0VertexFit, d0Cand, soft_pion,
			-1, hh_index, pion1, pion2);
	}
      }
    }

    // D0->Kpi
    if (recoD0Kpi_){
      const bmm::Candidate *daughter1(nullptr), *daughter2(nullptr);
      
      if (pion2.charge() == soft_pion.charge()){
	// Kpi case
	daughter1 = &kaon1;
	daughter2 = &pion2;
      } else {
	// piK case
	daughter1 = &pion1;
	daughter2 = &kaon2;
      }
	
      double d0_mass = (daughter1->p4() + daughter2->p4()).mass();
      double dstar_mass = (daughter1->p4() + daughter2->p4() + soft_pion.p4()).mass();

      if (d0_mass > minD0Mass_ && d0_mass < maxD0Mass_ &&
	  (dstar_mass - d0_mass) > min_dm_ && (dstar_mass - d0_mass) < max_dm_){
	
	pat::CompositeCandidate d0Cand(std::string("hh"));
	d0Cand.addDaughter( *daughter1 , "had1");
	d0Cand.addDaughter( *daughter2 , "had2");
	addP4.set( d0Cand );
	
	if (preprocess(d0Cand, iEvent, *daughter1, *daughter2)){
	  // Kinematic Fits
	  auto d0VertexFit = fillDileptonInfo(d0Cand, iEvent, *daughter1, *daughter2);
	  int hh_index = hh_collection.size();
	  hh_collection.push_back(d0Cand);
	  fillDstarInfo(dstar_collection, iEvent, d0VertexFit, d0Cand, soft_pion,
			-1, hh_index, *daughter1, *daughter2);
	}
      }
    }

  }
}


void DileptonPlusXProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  bField_ = &iSetup.getData(bFieldToken_);
  theTTBuilder_ = &iSetup.getData(theTTBuilderToken_);

  AnalyticalImpactPointExtrapolator extrapolator(bField_);
  impactPointExtrapolator_ = &extrapolator;

  edm::Handle<reco::BeamSpot> beamSpotHandle;
    
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("DileptonPlusXProducer") << "No beam spot available from EventSetup" ;
  }
    
  beamSpot_ = beamSpotHandle.product();
    
  iEvent.getByToken(vertexToken_, pvHandle_);
  // const reco::Vertex & primaryVertex = vertexHandle->front();

  edm::Handle<std::vector<pat::Muon>> muonHandle;
  edm::Handle<std::vector<pat::Electron>> electronHandle;
  edm::Handle<std::vector<pat::Photon>> photonHandle;
  edm::Handle<pat::CompositeCandidateCollection> conversionHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> lostTrackHandle;
    
  iEvent.getByToken(muonToken_, muonHandle);
  iEvent.getByToken(electronToken_, electronHandle);
  iEvent.getByToken(photonToken_, photonHandle);
  iEvent.getByToken(conversionToken_, conversionHandle);
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

  auto nMuons   = muonHandle->size();
  auto nPhotons = photonHandle->size();
  auto nConversions = conversionHandle->size();
  auto nPFCands = pfCandHandle_->size();
  // unsigned int lostTrackNumber = useLostTracks_ ? lostTrackHandle->size() : 0;
    
  // Output collection
  auto mm_collection = std::make_unique<pat::CompositeCandidateCollection>();
  auto ee_collection = std::make_unique<pat::CompositeCandidateCollection>();
  auto hh_collection = std::make_unique<pat::CompositeCandidateCollection>();
  auto btokmm  = std::make_unique<pat::CompositeCandidateCollection>();
  auto btokee  = std::make_unique<pat::CompositeCandidateCollection>();
  auto btokkmm = std::make_unique<pat::CompositeCandidateCollection>();
  auto btokkee = std::make_unique<pat::CompositeCandidateCollection>();
  auto btommg  = std::make_unique<pat::CompositeCandidateCollection>();
  auto dstar_collection = std::make_unique<pat::CompositeCandidateCollection>();
  AddFourMomenta addP4;

  // Build input lists

  // Good muons
  std::vector<bmm::Candidate> good_muon_candidates;
  for (unsigned int i = 0; i < nMuons; ++i) {
    const pat::Muon & muon = muonHandle->at(i);
    if (not isGoodMuon(muon)) continue;
    good_muon_candidates.push_back(bmm::Candidate(muon, i));
  }
    
  // Good electrons
  std::vector<bmm::Candidate> good_electron_candidates;
  if ( recoElElX_ ){
    for (unsigned int i = 0; i < electronHandle->size(); ++i) {
      const pat::Electron & el = electronHandle->at(i);
      if (not isGoodElectron(el)) continue;
      good_electron_candidates.push_back(bmm::Candidate(el, i));
    }
  }
    
  // std::vector<bmm::Candidate> good_hadron_candidates;
  // // Inject B to hh candidates where hadrons are explicitely matched
  // // to gen level decays 
  // if ( injectMatchedBtohh_ and isMC_ ) {
  //   injectHadronsThatMayFakeMuonsInMC(good_hadron_candidates);
  // }
  // // Inject reco B to hh candidates
  // if ( injectBtohh_ ) {
  //   injectBhhHadrons(good_hadron_candidates);
  // }
  // // Inject Jpsi to mumu based on charged tracks
  // if ( injectJpsiTracks_ ) {
  //   injectJpsiTracks(good_hadron_candidates);
  // }

  // Build dimuon candidates
  if ( good_muon_candidates.size() > 1 ){
    for (unsigned int i = 0; i < good_muon_candidates.size(); ++i) {
      const bmm::Candidate & muon1 = good_muon_candidates.at(i);
      for (unsigned int j = 0; j < good_muon_candidates.size(); ++j) {
	const bmm::Candidate & muon2 = good_muon_candidates.at(j);
	if (i==j) continue;

	// Enforce momentum ordering
	if (muon2.pt() > muon1.pt()) continue;

	pat::CompositeCandidate dimuonCand(std::string("mm"));
	dimuonCand.addDaughter( muon1 , "muon1");
	dimuonCand.addDaughter( muon2 , "muon2");
	addP4.set( dimuonCand );

	if (not preprocess(dimuonCand, iEvent, muon1, muon2)) continue;
	  
	// Kinematic Fits
	auto kinematicLLVertexFit = fillDileptonInfo(dimuonCand, iEvent, muon1, muon2);

	// dimuon + X
	int mm_index = mm_collection->size();
	    
	auto dimuon_p4(makeLorentzVectorFromPxPyPzM(kinematicLLVertexFit.p3().x(),
						    kinematicLLVertexFit.p3().y(),
						    kinematicLLVertexFit.p3().z(),
						    kinematicLLVertexFit.mass()));
	// MuMuGamma
	if (recoMuMuGamma_ && kinematicLLVertexFit.valid()){
	  for (unsigned int k=0; k < nPhotons; ++k){
	    auto photon(photonHandle->at(k));
	    if (photon.pt() < minGammaPt_) continue;
	    const auto & vtx_point = kinematicLLVertexFit.vtx_position();
	    photon.setVertex(reco::Photon::Point(vtx_point.x(), vtx_point.y(), vtx_point.z()));
	    double mmg_mass = (dimuon_p4 + photon.p4()).mass();
	    if (mmg_mass >= minLLGammaMass_ and mmg_mass <= maxLLGammaMass_){
	      // fill BtoLLPhoton candidate info
	      pat::CompositeCandidate mmgCand;
	      mmgCand.addUserInt("mm_index", mm_index);
	      mmgCand.addUserInt("ph_index", k);
	      mmgCand.addUserFloat("mass", mmg_mass);
		
	      fillLLGammaInfo(mmgCand,iEvent,muon1,muon2,photon);
	      fillMvaInfoForLLGamma(mmgCand,dimuonCand,iEvent,kinematicLLVertexFit,muon1,muon2,photon);
		
	      btommg->push_back(mmgCand);
	    }
	  }
	}

	// MuMuGamma with photon conversion
	if (recoMuMuGammaConv_ && kinematicLLVertexFit.valid()){
	  for (unsigned int k=0; k < nConversions; ++k){
	    auto conversion(conversionHandle->at(k));
	    if (conversion.pt() < minGammaPt_) continue;
	    double mmg_mass = (dimuon_p4 + conversion.p4()).mass();
	    if (mmg_mass >= minLLGammaMass_ and mmg_mass <= maxLLGammaMass_){
	      // fill BtoLLPhoton candidate info
	      pat::CompositeCandidate mmgCand;
	      mmgCand.addUserInt("mm_index", mm_index);
	      mmgCand.addUserInt("ph_index", -1);
	      mmgCand.addUserFloat("mass", mmg_mass);
		  
	      fillLLGammaConvInfo(mmgCand,iEvent, muon1, muon2, conversion);
	      // fillMvaInfoForLLGamma(mmgCand,dileptonCand,iEvent,kinematicLLVertexFit,lepton1,lepton2,photon);

	      btommg->push_back(mmgCand);
	    }
	  }
	}

	// mmK and mmKK
	buildLLXCandidates(*btokmm, *btokkmm, iEvent, kinematicLLVertexFit, dimuonCand, 
			   mm_index, muon1, muon2);

	// Dstar->D0pi->mmpi
	if (recoDstar_){
	  for (unsigned int k = 0; k < nPFCands; ++k) {
	    pat::PackedCandidate soft_pion((*pfCandHandle_)[k]);
	    if ( not isGoodHadron(soft_pion) ) continue;
    
	    if (overlap(muon1, soft_pion) || overlap(muon2, soft_pion)) continue;

	    soft_pion.setMass(PionMass_);

	    double d0_mass = (muon1.p4() + muon2.p4()).mass();
	    double dstar_mass = (muon1.p4() + muon2.p4() + soft_pion.p4()).mass();

	    if (d0_mass > minD0Mass_ && d0_mass < maxD0Mass_ &&
		(dstar_mass - d0_mass) > min_dm_ && (dstar_mass - d0_mass) < max_dm_){
		
	      fillDstarInfo(*dstar_collection, iEvent, kinematicLLVertexFit, dimuonCand, 
			    soft_pion, mm_index, -1, muon1, muon2);
	    }
	  }
	}

	// save dimuon
	mm_collection->push_back(dimuonCand);
      }
    } 
  }

  // Build dielectron candidates
  if ( good_electron_candidates.size() > 1 ){
    for (unsigned int i = 0; i < good_electron_candidates.size(); ++i) {
      const bmm::Candidate & electron1 = good_electron_candidates.at(i);
      for (unsigned int j = 0; j < good_electron_candidates.size(); ++j) {
	const bmm::Candidate & electron2 = good_electron_candidates.at(j);
	if (i==j) continue;

	// Enforce momentum ordering
	if (electron2.pt() > electron1.pt()) continue;

	pat::CompositeCandidate dielectronCand(std::string("ee"));
	dielectronCand.addDaughter( electron1 , "electron1");
	dielectronCand.addDaughter( electron2 , "electron2");
	addP4.set( dielectronCand );

	if (not preprocess(dielectronCand, iEvent, electron1, electron2)) continue;
	  
	// Kinematic Fits
	auto kinematicLLVertexFit = fillDileptonInfo(dielectronCand, iEvent, electron1, electron2);

	// dielectron + X
	int ee_index = ee_collection->size();
	    
	// ElElK and ElElKK
	buildLLXCandidates(*btokee, *btokkee, iEvent, kinematicLLVertexFit, dielectronCand, ee_index, electron1, electron2);

	// save dielectron
	ee_collection->push_back(dielectronCand);
	  
      }
    } 
  }

  // Build hh candidates
  // - loop over all hh combinations
  // - let individual studies fill hh_collection
  // - no check for duplicate entries
  if (nPFCands > 1){
    for (unsigned int i=0; i < nPFCands - 1; ++i){
      const auto& had1 = pfCandHandle_->at(i);
      if (not isGoodHadron(had1)) continue;
      for (unsigned int j=i + 1; j < nPFCands; ++j){
	const auto& had2 = pfCandHandle_->at(j);
	if (not isGoodHadron(had2)) continue;
	  
	if (had1.charge() * had2.charge() > 0) continue;
	  
	buildDstarCandidates(*dstar_collection, *hh_collection, iEvent, had1, had2);

	// // reco Btohh
	// if (dileptonCand.name() == "hh")
	//   if (not lepton1.from_gen() and not lepton2.from_gen())
	//     {
	// 	if (not kinematicLLVertexFit.valid()) continue;
	// 	if (kinematicLLVertexFit.vtxProb() < minBhhVtxProb_) continue;
	// 	if (kinematicLLVertexFit.sigLxy < minBhhSigLxy_) continue;
	//     }

      }
    }
  }
  // // dilepton + X
  // if (dileptonCand.name() == "mm" || dileptonCand.name() == "ee"){
  //   int ll_index = -1;
  //   if (dileptonCand.name() == "mm") ll_index = dimuon->size();
  //   if (dileptonCand.name() == "ee") ll_index = dielectron->size();
	    
  //   auto dilepton_p4(makeLorentzVectorFromPxPyPzM(kinematicLLVertexFit.p3().x(),
  // 						  kinematicLLVertexFit.p3().y(),
  // 						  kinematicLLVertexFit.p3().z(),
  // 						  kinematicLLVertexFit.mass()));
  //   // MuMuGamma
  //   if (recoMuMuGamma_ && dileptonCand.name() == "mm" && kinematicLLVertexFit.valid()){
  //     for (unsigned int k=0; k < nPhotons; ++k){
  // 	auto photon(photonHandle->at(k));
  // 	if (photon.pt() < minGammaPt_) continue;
  // 	const auto & vtx_point = kinematicLLVertexFit.refitVertex->vertexState().position();
  // 	photon.setVertex(reco::Photon::Point(vtx_point.x(), vtx_point.y(), vtx_point.z()));
  // 	double mmg_mass = (dilepton_p4 + photon.p4()).mass();
  // 	if (mmg_mass >= minLLGammaMass_ and mmg_mass <= maxLLGammaMass_){
  // 	  // fill BtoLLPhoton candidate info
  // 	  pat::CompositeCandidate mmgCand;
  // 	  mmgCand.addUserInt("mm_index", ll_index);
  // 	  mmgCand.addUserInt("ph_index", k);
  // 	  mmgCand.addUserFloat("mass", mmg_mass);
		  
  // 	  fillLLGammaInfo(mmgCand,iEvent,lepton1,lepton2,photon);
  // 	  fillMvaInfoForLLGamma(mmgCand,dileptonCand,iEvent,kinematicLLVertexFit,lepton1,lepton2,photon);

  // 	  btommg->push_back(mmgCand);
  // 	}
  //     }
  //   }

  //   // MuMuGamma with photon conversion
  //   if (recoMuMuGammaConv_ && dileptonCand.name() == "mm" && kinematicLLVertexFit.valid()){
  //     for (unsigned int k=0; k < nConversions; ++k){
  // 	auto conversion(conversionHandle->at(k));
  // 	if (conversion.pt() < minGammaPt_) continue;
  // 	double mmg_mass = (dilepton_p4 + conversion.p4()).mass();
  // 	if (mmg_mass >= minLLGammaMass_ and mmg_mass <= maxLLGammaMass_){
  // 	  // fill BtoLLPhoton candidate info
  // 	  pat::CompositeCandidate mmgCand;
  // 	  mmgCand.addUserInt("mm_index", ll_index);
  // 	  mmgCand.addUserInt("ph_index", -1);
  // 	  mmgCand.addUserFloat("mass", mmg_mass);
		  
  // 	  fillLLGammaConvInfo(mmgCand,iEvent, lepton1, lepton2, conversion);
  // 	  // fillMvaInfoForLLGamma(mmgCand,dileptonCand,iEvent,kinematicLLVertexFit,lepton1,lepton2,photon);

  // 	  btommg->push_back(mmgCand);
  // 	}
  //     }
  //   }

  //   // Dstar to D0 pi, D0 to MuMu
  //   if (recoDstar_ && dileptonCand.name() == "mm" && kinematicLLVertexFit.valid() &&
  // 	dilepton_p4.mass() > minD0Mass_ && dilepton_p4.mass() < maxD0Mass_){
  //     for (unsigned int k = 0; k < nPFCands; ++k) {
  // 	pat::PackedCandidate pionCand((*pfCandHandle_)[k]);
  // 	if ( abs(pionCand.pdgId()) != 211 ) continue;
  // 	if (pionCand.charge() == 0 ) continue;
  // 	if ( not pionCand.hasTrackDetails() ) continue;
  // 	pionCand.setMass(PionMass_);
  // 	if (overlap(lepton1, pionCand) || overlap(lepton2, pionCand)) continue;
  // 	double l1_pion_doca = distanceOfClosestApproach(lepton1.track(),
  // 							pionCand.bestTrack());
  // 	double l2_pion_doca = distanceOfClosestApproach(lepton2.track(),
  // 							pionCand.bestTrack());
  // 	if (maxTwoTrackDOCA_ > 0 and l1_pion_doca > maxTwoTrackDOCA_) continue;
  // 	if (maxTwoTrackDOCA_ > 0 and l2_pion_doca > maxTwoTrackDOCA_) continue;
		
  // 	double mmpi_mass = (dilepton_p4 + pionCand.p4()).mass();
  // 	if (mmpi_mass < minDstarMass_ or mmpi_mass > maxDstarMass_) continue;
		
  // 	pat::CompositeCandidate dstarCand;
  // 	dstarCand.addUserInt(dileptonCand.name() + "_index", ll_index);
  // 	dstarCand.addUserFloat("pion_l1_doca", l1_pion_doca);
  // 	dstarCand.addUserFloat("pion_l2_doca", l2_pion_doca);
		
  // 	fillDstarInfo(dstarCand, iEvent, lepton1, lepton2, pionCand);

  // 	dstar2mmpi->push_back(dstarCand);
  //     }
  //   }
	    
  //   // llK and llKK
  //   for (unsigned int k = 0; k < nPFCands; ++k) {
  //     pat::PackedCandidate kaonCand1((*pfCandHandle_)[k]);
  //     kaonCand1.setMass(KaonMass_);
  //     if (kaonCand1.charge() == 0 ) continue;
  //     if (!kaonCand1.hasTrackDetails()) continue;
  //     if (abs(kaonCand1.pdgId()) != 211) continue;
  //     if (kaonCand1.pt() < ptMinKaon_ or abs(kaonCand1.eta()) > etaMaxKaon_) continue;
  //     if (overlap(lepton1, kaonCand1) || overlap(lepton2, kaonCand1)) continue;
  //     double l1_kaon_doca = distanceOfClosestApproach(lepton1.track(),
  // 						       kaonCand1.bestTrack());
  //     double l2_kaon_doca = distanceOfClosestApproach(lepton2.track(),
  // 						       kaonCand1.bestTrack());
  //     if (maxTwoTrackDOCA_ > 0 and l1_kaon_doca > maxTwoTrackDOCA_) continue;
  //     if (maxTwoTrackDOCA_ > 0 and l2_kaon_doca > maxTwoTrackDOCA_) continue;
	      
	      
  //     bool goodBtoLLK = true;

  //     double kll_mass = (lepton1.p4() + lepton2.p4() + kaonCand1.p4()).mass();
  //     if (kll_mass < minBKllMass_ || kll_mass > maxBKllMass_) goodBtoLLK = false;
	    
  //     // fill BtoLLK candidate info
  //     if (goodBtoLLK){
  // 	pat::CompositeCandidate btokllCand;
  // 	btokllCand.addUserInt(dileptonCand.name() + "_index", ll_index);
  // 	btokllCand.addUserFloat("kaon_l1_doca", l1_kaon_doca);
  // 	btokllCand.addUserFloat("kaon_l2_doca", l2_kaon_doca);
		
  // 	fillBtoKllInfo(btokllCand, iEvent, lepton1, lepton2, kaonCand1);
  // 	fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokllCand, dileptonCand, iEvent, kinematicLLVertexFit, lepton1, lepton2, kaonCand1);

  // 	if (dileptonCand.name() == "mm")
  // 	  btokmm->push_back(btokllCand);
		
  // 	if (dileptonCand.name() == "ee")
  // 	  btokee->push_back(btokllCand);
  //     }

  //     // Build BsToKKll
  //     for (unsigned int k2 = k+1; k2 < nPFCands; ++k2) { // only works if selection requirements for both kaons are identical
  // 	pat::PackedCandidate kaonCand2((*pfCandHandle_)[k2]);
  // 	kaonCand2.setMass(KaonMass_);
  // 	if (kaonCand2.charge() == 0 ) continue;
  // 	if (!kaonCand2.hasTrackDetails()) continue;
  // 	if (abs(kaonCand2.pdgId()) != 211) continue;
  // 	if (kaonCand2.pt() < ptMinKaon_ || abs(kaonCand2.eta()) > etaMaxKaon_) continue;
  // 	if (overlap(lepton1, kaonCand2) || overlap(lepton2, kaonCand2)) continue;
  // 	double l1_kaon2_doca = distanceOfClosestApproach(lepton1.track(),
  // 							 kaonCand2.bestTrack());
  // 	double l2_kaon2_doca = distanceOfClosestApproach(lepton2.track(),
  // 							 kaonCand2.bestTrack());
  // 	if (maxTwoTrackDOCA_>0 and l1_kaon2_doca > maxTwoTrackDOCA_) continue;
  // 	if (maxTwoTrackDOCA_>0 and l2_kaon2_doca > maxTwoTrackDOCA_) continue;
		
  // 	bool goodBtoLLKK = true;
		      
  // 	double kkll_mass = (lepton1.p4() + lepton2.p4() + kaonCand1.p4() + kaonCand2.p4()).mass();
  // 	if ( kkll_mass < minBKKllMass_ || kkll_mass > maxBKKllMass_ ) goodBtoLLKK = false;
		  
  // 	// fill BtoLLKK candidate info
  // 	if (goodBtoLLKK){
  // 	  pat::CompositeCandidate btokkllCand;
  // 	  btokkllCand.addUserInt(dileptonCand.name() + "_index", ll_index);
  // 	  btokkllCand.addUserFloat("kaon1_l1_doca", l1_kaon_doca);
  // 	  btokkllCand.addUserFloat("kaon1_l2_doca", l2_kaon_doca);
  // 	  btokkllCand.addUserFloat("kaon2_l1_doca", l1_kaon2_doca);
  // 	  btokkllCand.addUserFloat("kaon2_l2_doca", l2_kaon2_doca);
		    
  // 	  fillBtoKKllInfo(btokkllCand, iEvent, lepton1, lepton2, kaonCand1, kaonCand2);
  // 	  // FIXME
  // 	  // fillMvaInfoForBtoJpsiKCandidatesEmulatingBmm(btokkllCand,dileptonCand,iEvent,kinematicLLVertexFit,lepton1,lepton2,kaonCand1);

  // 	  if (dileptonCand.name() == "mm")
  // 	    btokkmm->push_back(btokkllCand);
		  
  // 	  if (dileptonCand.name() == "ee")
  // 	    btokkee->push_back(btokkllCand);
  // 	}
  //     }
  //   }                  
  // }
	  
  iEvent.put(std::move(mm_collection),     "MuMu");
  iEvent.put(std::move(ee_collection), "ElEl");
  iEvent.put(std::move(hh_collection),   "HH");
  iEvent.put(std::move(btokmm), "BToKmumu");
  iEvent.put(std::move(btokee), "BToKee");
  iEvent.put(std::move(btokkmm),"BToKKmumu");
  iEvent.put(std::move(btokkee),"BToKKee");
  iEvent.put(std::move(btommg), "BToMuMuGamma");
  iEvent.put(std::move(dstar_collection), "Dstar");
}

KalmanVertexFitResult 
DileptonPlusXProducer::vertexWithKalmanFitter(std::vector<const reco::Track*> trks, 
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
DileptonPlusXProducer::vertexLeptonsWithKalmanFitter(const bmm::Candidate& lepton1,
						     const bmm::Candidate& lepton2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( lepton1.track() );
  masses.push_back( lepton1.mass() );
  trks.push_back( lepton2.track() );
  masses.push_back( lepton2.mass() );
  return vertexWithKalmanFitter(trks, masses);
}


KinematicFitResult 
DileptonPlusXProducer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
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

  RefCountedKinematicTree vertexFitTree;
  KinematicFitResult result;
  result.tracks = trks;
  try {
    vertexFitTree = fitter.fit(particles);
  } catch (const std::exception& e) {
    return result;
  }
  result.set_tree(vertexFitTree);
  return result;
}


KinematicFitResult
DileptonPlusXProducer::vertexLeptonsWithKinematicFitter(const bmm::Candidate& lepton1,
							const bmm::Candidate& lepton2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( lepton1.track() );
  masses.push_back( lepton1.mass() );
  trks.push_back( lepton2.track() );
  masses.push_back( lepton2.mass() );
  return vertexWithKinematicFitter(trks, masses);
}

KinematicFitResult 
DileptonPlusXProducer::vertexKaonsWithKinematicFitter(const pat::PackedCandidate& pfCand1,
						      const pat::PackedCandidate& pfCand2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( pfCand1.bestTrack() );
  masses.push_back( KaonMass_ );
  trks.push_back( pfCand2.bestTrack() );
  masses.push_back( KaonMass_ );
  return vertexWithKinematicFitter(trks, masses);
}


KinematicFitResult
DileptonPlusXProducer::vertexWithKinematicFitter(const bmm::Candidate& lepton1,
						 const bmm::Candidate& lepton2,
						 const pat::PackedCandidate& pion)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( lepton1.track() );
  masses.push_back( lepton1.mass() );
  trks.push_back( lepton2.track() );
  masses.push_back( lepton2.mass() );
  trks.push_back( pion.bestTrack() );
  masses.push_back(PionMass_);
  return vertexWithKinematicFitter(trks,masses);
}

KinematicFitResult
DileptonPlusXProducer::fitDstar(const bmm::Candidate& lepton1,
				const bmm::Candidate& lepton2,
				const pat::PackedCandidate& pion,
				float mass_constraint)
{
  // Rebuild ll vertex to ensure that the KinematicTree remains self
  // consistent and no elements get out of scope or get deleted
  // when the tree is used in subsequent fits
  auto llVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  auto tree = llVertexFit.tree();
  
  KinematicFitResult result; 
  if (lepton1.track()) result.tracks.push_back(lepton1.track());
  if (lepton2.track()) result.tracks.push_back(lepton2.track());
  if (pion.bestTrack()) result.tracks.push_back(pion.bestTrack());

  if ( not llVertexFit.valid()) return result;

  KinematicConstraint* mc(0);
  if (mass_constraint > 0){
    ParticleMass mass = mass_constraint;
    // mass constraint fit
    KinematicParticleFitter csFitter;
    float mass_sigma = 1e-4;
    // FIXME: potential memory leak
    mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      tree = csFitter.fit(mc, tree);
    } catch (const std::exception& e) {
      return result;
    }
  }

  const reco::TransientTrack pionTT = theTTBuilder_->build(pion.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> particles;
  double chi = 0.;
  double ndf = 0.;

  tree->movePointerToTheTop();
  particles.push_back(tree->currentParticle());
  float pionMassErr(PionMassErr_);
  particles.push_back(partFactory.particle(pionTT,PionMass_,chi,ndf,pionMassErr));

  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(particles);
  } catch (const std::exception& e) {
    return result;
  }

  result.set_tree(vertexFitTree);

  return result;
}


KinematicFitResult
DileptonPlusXProducer::fitBToKLL(const bmm::Candidate& lepton1,
				 const bmm::Candidate& lepton2,
				 const pat::PackedCandidate& kaon,
				 float mass_constraint)
{
  // Rebuild ll vertex to ensure that the KinematicTree remains self
  // consistent and no elements get out of scope or get deleted
  // when the tree is used in subsequent fits
  auto llVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  auto tree = llVertexFit.tree();
  
  KinematicFitResult result;
  if (lepton1.track()) result.tracks.push_back(lepton1.track());
  if (lepton2.track()) result.tracks.push_back(lepton2.track());
  if (kaon.bestTrack()) result.tracks.push_back(kaon.bestTrack());

  if ( not llVertexFit.valid()) return result;

  KinematicConstraint* mc(0);
  if (mass_constraint > 0){
    ParticleMass mass = mass_constraint;
    // mass constraint fit
    KinematicParticleFitter csFitter;
    float mass_sigma = JPsiMassErr_;
    // FIXME: memory leak
    mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      tree = csFitter.fit(mc, tree);
    } catch (const std::exception& e) {
      return result;
    }
  }

  const reco::TransientTrack kaonTT = theTTBuilder_->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> particles;
  double chi = 0.;
  double ndf = 0.;

  tree->movePointerToTheTop();
  particles.push_back(tree->currentParticle());
  float kaonMassErr(KaonMassErr_);
  particles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,kaonMassErr));

  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(particles);
  } catch (const std::exception& e) {
    return result;
  }

  result.set_tree(vertexFitTree);
  return result;
}

KinematicFitResult
DileptonPlusXProducer::fitLLGamma( RefCountedKinematicTree tree,
				   const pat::Photon& photon,
				   float mass_constraint)
{
  KinematicFitResult result; 
  if ( !tree->isValid()) return result;
  tree->movePointerToTheTop();
  
  KinematicConstraint* mc(0);
  if (mass_constraint > 0){
    ParticleMass mass = mass_constraint;
    // mass constraint fit
    KinematicParticleFitter csFitter;
    float mass_sigma = JPsiMassErr_;
    // FIXME: memory leak
    mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      tree = csFitter.fit(mc, tree);
    } catch (const std::exception& e) {
      return result;
    }
  }
  
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> particles;
  float chi2 = 0.;
  float ndf = 0.;

  tree->movePointerToTheTop();
  particles.push_back(tree->currentParticle());
  bmm::KinematicParticleRef
    photon_kp(bmm::build_particle(photon, momentum_resolution(photon),
				  bField_, chi2, ndf));

  particles.push_back(photon_kp);

  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(particles);
  } catch (const std::exception& e) {
    return result;
  }

  result.set_tree(vertexFitTree);
  return result;
}

KinematicFitResult
DileptonPlusXProducer::fitLLGammaConv( RefCountedKinematicTree mmVertexTree,
				       const pat::CompositeCandidate& photon,
				       float mass_constraint)
{
  KinematicFitResult result; 
  if ( !mmVertexTree->isValid()) return result;

  // Refit mm if necessary
  
  mmVertexTree->movePointerToTheTop();
  
  KinematicParticleFitter csFitter;
  
  KinematicConstraint* mm_mc(0);
  if (mass_constraint > 0){
    ParticleMass mass = mass_constraint;
    // mass constraint fit
    float mass_sigma = JPsiMassErr_;
    // FIXME: memory leak
    mm_mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      mmVertexTree = csFitter.fit(mm_mc, mmVertexTree);
    } catch (const std::exception& e) {
      return result;
    }
  }

  // Build ee vertex and refit

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter vtxFitter;

  ColinearityKinematicConstraintT<colinearityKinematic::PhiTheta> constr;
  KinematicConstrainedVertexFitterT<2, 2> kcvFitter(bField_);
  
  std::vector<RefCountedKinematicParticle> particles;
  
  auto tk0 = photon.userData<reco::Track>("track0");
  auto tt0 = theTTBuilder_->build(*tk0);
  auto tk1 = photon.userData<reco::Track>("track1");
  auto tt1 = theTTBuilder_->build(*tk1);

  float ElectronMassErr(ElectronMassErr_);
  
  particles.push_back(partFactory.particle(tt0, ElectronMass_, float(0), float(0), ElectronMassErr));
  particles.push_back(partFactory.particle(tt1, ElectronMass_, float(0), float(0), ElectronMassErr));

  RefCountedKinematicTree photonVertexTree;
  try {
    // photonVertexTree = vtxFitter.fit(particles);
    // default parameters
    // - maxDelta: 0.01
    // - maxNbrOfIterations: 1000
    // - maxReducedChiSq: 225.
    // - minChiSqImprovement: 50.
    // Note: the maximum number of convergence iterations is the only parameter that
    // tends to be changed. 40 is used in EgammaPhotonProducers
    photonVertexTree = kcvFitter.fit(particles, &constr);
  } catch (const std::exception& e) {
    return result;
  }
  
  if ( !photonVertexTree->isValid()) return result;

  // // add mass constraint to the photon candidate

  // const ParticleMass photon_mass(0);
  // float photon_mass_err(1e-6);

  // KinematicConstraint* photon_mc = new MassKinematicConstraint(photon_mass, photon_mass_err);
  // photonVertexTree->movePointerToTheTop();
  // try {
  //   photonVertexTree = csFitter.fit(photon_mc, photonVertexTree);
  // } catch (const std::exception& e) {
  //   return result;
  // }

  // if ( !photonVertexTree->isValid()) return result;

  // common vertex fit
  std::vector<RefCountedKinematicParticle> mmg_particles;

  mmVertexTree->movePointerToTheTop();
  mmg_particles.push_back(mmVertexTree->currentParticle());

  photonVertexTree->movePointerToTheTop(); 
  mmg_particles.push_back(photonVertexTree->currentParticle());
  
  RefCountedKinematicTree mmgVertexTree;
  try {
    mmgVertexTree = vtxFitter.fit(mmg_particles);
  } catch (const std::exception& e) {
    return result;
  }
  
  result.set_tree(mmgVertexTree);
  return result;
}



KinematicFitResult
DileptonPlusXProducer::fitBToKKLL( const bmm::Candidate& lepton1,
				   const bmm::Candidate& lepton2,
				   const pat::PackedCandidate& kaon1,
				   const pat::PackedCandidate& kaon2,
				   float ll_mass_constraint,
				   float kk_mass_constraint)
{
  // Rebuild ll vertex to ensure that the KinematicTree remains self
  // consistent and no elements get out of scope or get deleted
  // when the tree is used in subsequent fits
  auto llVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  auto ll_tree = llVertexFit.tree();

  KinematicFitResult result; 
  if (lepton1.track()) result.tracks.push_back(lepton1.track());
  if (lepton2.track()) result.tracks.push_back(lepton2.track());
  if (kaon1.bestTrack()) result.tracks.push_back(kaon1.bestTrack());
  if (kaon2.bestTrack()) result.tracks.push_back(kaon2.bestTrack());

  if ( not llVertexFit.valid()) return result;

  auto kkVertexFit = vertexKaonsWithKinematicFitter(kaon1, kaon2);
  auto kk_tree = kkVertexFit.tree();

  if ( not kkVertexFit.valid()) return result;
  
  KinematicConstraint* ll_mc(0);
  if (ll_mass_constraint > 0){
    ParticleMass mass = ll_mass_constraint;
    // mass constraint fit
    KinematicParticleFitter csFitter;
    float mass_sigma = 1e-4;
    // FIXME: potential memory leak
    ll_mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      ll_tree = csFitter.fit(ll_mc, ll_tree);
    } catch (const std::exception& e) {
      return result;
    }
  }

  KinematicConstraint* kk_mc(0);
  if (kk_mass_constraint > 0){
    ParticleMass mass = kk_mass_constraint;
    // mass constraint fit
    KinematicParticleFitter csFitter;
    float mass_sigma = 1e-4;
    // FIXME: potential memory leak
    kk_mc = new MassKinematicConstraint(mass, mass_sigma);
    try {
      kk_tree = csFitter.fit(kk_mc, kk_tree);
    } catch (const std::exception& e) {
      return result;
    }
  }

  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> particles;

  ll_tree->movePointerToTheTop();
  particles.push_back(ll_tree->currentParticle());
  kk_tree->movePointerToTheTop();
  particles.push_back(kk_tree->currentParticle());

  RefCountedKinematicTree vertexFitTree;
  try {
    vertexFitTree = fitter.fit(particles);
  } catch (const std::exception& e) {
    return result;
  }

  result.set_tree(vertexFitTree);
  return result;
}

KinematicFitResult
DileptonPlusXProducer::vertexLeptonsWithPointingConstraint( const bmm::Candidate& lepton1,
							    const bmm::Candidate& lepton2,
							    const reco::Vertex& primaryVertex)
{
  KinematicFitResult result; 
  if (lepton1.track()) result.tracks.push_back(lepton1.track());
  if (lepton2.track()) result.tracks.push_back(lepton2.track());

  auto kinematicLLVertexFit = vertexLeptonsWithKinematicFitter(lepton1, lepton2);
  kinematicLLVertexFit.postprocess(*beamSpot_);
  if ( !kinematicLLVertexFit.valid()) return result;
  auto tree = kinematicLLVertexFit.tree();
  if ( !tree->isValid()) return result;

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
    return result;
  }
  result.set_tree(refittedTree);
  return result;
}

pair<double,double> DileptonPlusXProducer::computeDCA(const pat::PackedCandidate &kaon,
						      reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), bField_);

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}

namespace{
}

const reco::Candidate* DileptonPlusXProducer::getGenParticle(const bmm::Candidate& cand)
{
  
  if (cand.genParticle()) return cand.genParticle();

  for (auto const & genParticle: *packedGenParticles_){
    if (dr_match(cand.p4(), genParticle.p4()))
      return &genParticle;
  }
  return nullptr;
}

GenMatchInfo DileptonPlusXProducer::getGenMatchInfo( const bmm::Candidate& lepton1,
						     const bmm::Candidate& lepton2,
						     const pat::PackedCandidate* kaon1,
						     const pat::PackedCandidate* kaon2,
						     const reco::Candidate* photon)
{
  auto result = GenMatchInfo();
  const reco::Candidate*   ll_mother(0);
  assert(prunedGenParticles_);
  assert(packedGenParticles_);
  std::vector<const reco::Candidate*> daughters;

  result.mc_l1 = getGenParticle(lepton1);
  if (result.mc_l1){
    result.l1_pdgId = result.mc_l1->pdgId();
    result.l1_pt    = result.mc_l1->pt();
    if (result.mc_l1->mother()){
      result.l1_motherPdgId = result.mc_l1->mother()->pdgId();
    }
    daughters.push_back(result.mc_l1);
  }

  result.mc_l2 = getGenParticle(lepton2);
  if (result.mc_l2){
    result.l2_pdgId = result.mc_l2->pdgId();
    result.l2_pt    = result.mc_l2->pt();
    if (result.mc_l2->mother()){
      result.l2_motherPdgId = result.mc_l2->mother()->pdgId();
    }
    daughters.push_back(result.mc_l2);
  }

  if ( result.mc_l1 and result.mc_l2 ){
    if ( (result.mc_l1->vertex() - result.mc_l2->vertex()).r() < 1e-4)
      result.ll_vtx    = result.mc_l1->vertex();
    if ( result.mc_l1->mother() and result.mc_l1->mother() == result.mc_l2->mother() ){
      ll_mother = result.mc_l1->mother();
      result.match = result.mc_l1->mother();
      result.ll_mass      = ll_mother->mass();
      result.ll_pt        = ll_mother->pt();
      result.ll_pdgId     = ll_mother->pdgId();
      if (ll_mother->mother()) result.ll_motherPdgId = ll_mother->mother()->pdgId();
      result.ll_prod_vtx = getProductionVertex(ll_mother);
    }
  }
  
  if (kaon1){
    for (auto const & genParticle: *packedGenParticles_){
      if (dr_match(kaon1->p4(), genParticle.p4())){
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
	result.kll_pdgId    = mother->pdgId();
	result.kll_motherPdgId = mother->mother() ? mother->mother()->pdgId() : 0;
	result.kll_mass     = mother->mass();
	result.kll_pt       = mother->pt();
	result.kll_prod_vtx = getProductionVertex(mother);
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
	result.kkll_pdgId    = mother->pdgId();
	result.kkll_mass     = mother->mass();
	result.kkll_pt       = mother->pt();
	result.kkll_prod_vtx = getProductionVertex(mother);
      }
    }
  }
  
  if (photon){
    for (auto const & genParticle: *packedGenParticles_){
      if (dr_match(photon->p4(),genParticle.p4())){
	result.mc_photon = &genParticle;
	daughters.push_back(result.mc_photon);
	result.photon_pdgId = genParticle.pdgId();
	result.photon_pt    = genParticle.pt();
	if (genParticle.mother(0)){
	  result.photon_motherPdgId = genParticle.mother(0)->pdgId();
	}
	break;
      }
    }
    if (daughters.size()==3){
      const auto* mother = find_common_ancestor(daughters);
      if (mother){
	result.match        = mother;
	result.llg_pdgId    = mother->pdgId();
	result.llg_mass     = mother->mass();
	result.llg_pt       = mother->pt();
	result.llg_prod_vtx = getProductionVertex(mother);
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

float DileptonPlusXProducer::distanceOfClosestApproach( const reco::GenParticle* track1,
							const reco::GenParticle* track2)
{
  TwoTrackMinimumDistance md;
  GlobalPoint trk1_pos(track1->vertex().x(), 
		       track1->vertex().y(), 
		       track1->vertex().z());
  GlobalVector trk1_mom(track1->px(),track1->py(),track1->pz());

  GlobalTrajectoryParameters trk1(trk1_pos,trk1_mom,track1->charge(),bField_);
  GlobalPoint trk2_pos(track2->vertex().x(), 
		       track2->vertex().y(), 
		       track2->vertex().z());
  GlobalVector trk2_mom(track2->px(),track2->py(),track2->pz());

  GlobalTrajectoryParameters trk2(trk2_pos,trk2_mom,track2->charge(),bField_);
  if ( not md.calculate( trk1, trk2 ) ) return -1.0;
  return md.distance();
}

float DileptonPlusXProducer::distanceOfClosestApproach( const reco::Track* track1,
							const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
DileptonPlusXProducer::distanceOfClosestApproach( const reco::Track* track,
						  const VertexState& vertex_state)
{
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), vertex_state.position());
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex_state);
  return doca;
}

Measurement1D 
DileptonPlusXProducer::distanceOfClosestApproach( const reco::Track* track,
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

template <typename T> std::vector<const reco::Track*>
DileptonPlusXProducer::getGoodTracksToRefitPV(int pvIndex, T track_to_ignore)
{
  return getGoodTracksToRefitPV(pvIndex, std::vector<T>(1,track_to_ignore));
}

template <typename T> std::vector<const reco::Track*>
DileptonPlusXProducer::getGoodTracksToRefitPV(int pvIndex, std::vector<T> ignoreTracks)
{
  std::vector<const reco::Track*> tracks;

  for (const auto& pfCand: *pfCandHandle_.product()){
    if (not isGoodTrack(pfCand)) continue; 
    if (int(pfCand.vertexRef().key()) != pvIndex) continue;
    // keep only the tracks used in the PV fit
    if (pfCand.pvAssociationQuality() != pat::PackedCandidate::UsedInFitTight) continue;
    bool keep_track = true;
    for (const auto& track: ignoreTracks){
      if (overlap(track, pfCand.bestTrack())) {
	keep_track = false;
	break;
      }
    }
    if (keep_track)
      tracks.push_back(pfCand.bestTrack());
  }
  return tracks;
}


std::pair<KinematicFitResult, KinematicFitResult>
DileptonPlusXProducer::refitWithVertexConstraint(const reco::Track& track,
						 int pvIndex)
{
  std::vector<const reco::Track*> tracks(getGoodTracksToRefitPV(pvIndex, &track));
  std::vector<float> masses(tracks.size(), PionMass_);

  auto pv_refit = vertexWithKinematicFitter(tracks, masses);


  // make sure the track is the last element
  tracks.push_back(&track);
  masses.push_back(PionMass_);
  auto pv_refit_with_track = vertexWithKinematicFitter(tracks, masses);
 
  return std::make_pair(pv_refit, pv_refit_with_track);
}

bmm::Displacements
DileptonPlusXProducer::compute3dDisplacement(const KinematicFitResult& fit, 
					     bool closestIn3D)
{
  // WARNING: all variables need to be filled for even if the fit is not valid

  bmm::Displacements result;

  const reco::Vertex* bestVertex(0);
  int bestVertexIndex(-1);
  const reco::Vertex* bestVertex2(0);
  int bestVertexIndex2(-1);

  if (fit.valid()){

    const auto& vertices = *pvHandle_.product();

    auto candTransientTrack = fit.particle()->refittedTransientTrack();
    
    // find best matching primary vertex
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
  }

  if (bestVertex)
    result.push_back(bmm::Displacement("pv", fit, *bestVertex, bestVertexIndex));
  else
    result.push_back(bmm::Displacement("pv"));

  if (bestVertex2)
    result.push_back(bmm::Displacement("_pv2", fit, *bestVertex2, bestVertexIndex2));
  else
    result.push_back(bmm::Displacement("_pv2"));

  // Refit the primary vertex
  if (bestVertex)
    {
      std::vector<const reco::Track*> tracks(getGoodTracksToRefitPV(bestVertexIndex, fit.tracks));
      std::vector<float> masses(tracks.size(), PionMass_);
      auto pv_refit = vertexWithKinematicFitter(tracks, masses);
      result.push_back(bmm::Displacement("_refit", fit, pv_refit.vertex(), bestVertexIndex));
    }
  else
    result.push_back(bmm::Displacement("_refit"));
    
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

void DileptonPlusXProducer::setupTmvaReader(TMVA::Reader& reader, std::string file){
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
DileptonPlusXProducer::computeAnalysisBDT(unsigned int event_idx)
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

DEFINE_FWK_MODULE(DileptonPlusXProducer);

//  LocalWords:  vertices
