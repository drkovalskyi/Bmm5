#include "RecoVertex/KinematicFitPrimitives/interface/TrackKinematicStatePropagator.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateAccessor.h"
#include "DataFormats/GeometrySurface/interface/Surface.h" 
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrajectoryState/interface/BasicSingleTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToCurvilinear.h"
#include "TrackingTools/TrajectoryState/interface/PerigeeConversions.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/PointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "TMath.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include <iostream>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LorentzVectorM;
typedef ROOT::Math::SMatrix<double, 6, 6, ROOT::Math::MatRepSym<double, 6> > Matrix6S;
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > Matrix3S;
typedef ROOT::Math::SMatrix<double, 3, 3> Matrix33;

class ConstMagneticField : public MagneticField {
public:

  virtual GlobalVector inTesla ( const GlobalPoint& ) const {
    return GlobalVector(0,0,4);
  }

};

LorentzVectorM makeLorentzVector(const GlobalVector &momentum, double mass){
  return LorentzVectorM(momentum.x(), momentum.y(), momentum.z(), mass);
}

void dump_result( const RefCountedKinematicTree& fit ){
  std::cout << "Fit results\n=== Mother ===" << std::endl;
  
  fit->movePointerToTheTop();
  auto kp = fit->currentParticle();
  auto state = kp->currentState();
  
  std::cout << "\tmass: " <<  state.mass() << std::endl;
  std::cout << "\tmass uncertainty: " << sqrt(state.kinematicParametersError().matrix()(6,6)) << std::endl;

  double prob = TMath::Prob((double)fit->currentDecayVertex()->chiSquared(), int(rint(fit->currentDecayVertex()->degreesOfFreedom())));
  std::cout << "\tprobability: " << prob << std::endl;

  auto fts = state.freeTrajectoryState();
  std::cout << fts << std::endl;

  std::cout << "=== Daughters ===" << std::endl;
  
  if ( fit->movePointerToTheFirstChild() ){
    do {
      std::cout << fit->currentParticle()->currentState().freeTrajectoryState() << std::endl;
    } while (fit->movePointerToTheNextChild());
  }

  std::cout << "================" << std::endl;
  fit->movePointerToTheTop();

}

Matrix33 jacobianSphToCart(const reco::Candidate::LorentzVector& p4){
  Matrix33 jac;
  jac(0,0) =   p4.Px() / p4.P();
  jac(0,1) =   p4.Py() / p4.P();
  jac(0,2) =   p4.Pz() / p4.P();
  jac(1,0) =   p4.Px() * p4.Pz() / p4.Pt() / p4.P() / p4.P();
  jac(1,1) =   p4.Py() * p4.Pz() / p4.Pt() / p4.P() / p4.P();
  jac(1,2) = - p4.Pt() / p4.P() / p4.P();
  jac(2,0) = - p4.Py() / p4.Pt() / p4.Pt();
  jac(2,1) =   p4.Px() / p4.Pt() / p4.Pt();
  jac(2,2) =   0;
  
  return jac;
}

int main() {

  float muon_mass(0.10565837);
  float muon_mass_err(3.5 * 1e-9);
  
  MagneticField * field = new ConstMagneticField;

  //
  // Reco particles
  //
  
  // Muon 1
  
  // pt(): 7.33023
  // eta(): -0.362706
  double e1[] = {1.49217e-06, 3.87728e-09, 1.72144e-07, 2.8663e-07, 4.97307e-10,
    1.62012e-07, -2.22939e-06, 7.02559e-09, -1.49479e-06, 1.52818e-05, -7.57048e-08,
    -2.36024e-06, -1.41106e-08, -5.46833e-08, 3.6118e-05};

  reco::TrackBase::CovarianceMatrix cov1(e1, e1 + 15);
  reco::Track::Point v1(0.00545092, 0.0481168, -6.28238);
  reco::Track::Vector p1(-5.55419, -4.78364, -2.71739);
  
  reco::Track trk1(6.15377, 14., v1, p1, 1, cov1);
  reco::TransientTrack tt1(trk1, field);

  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> mu1_p4(p1.x(), p1.y(), p1.z(), muon_mass);
  
  // Muon 2

  // pt(): 4.10482
  // eta(): -0.522372
  double e2[] = {4.64278e-06, -4.88327e-09, 1.92113e-07, 6.29305e-07, -1.07874e-09,
    2.93906e-07, -3.11981e-06, -3.72671e-10, -1.96178e-06, 1.42558e-05,
    8.35985e-08, -1.76157e-06,  2.81995e-08, -1.18262e-07, 1.76285e-05};
  reco::TrackBase::CovarianceMatrix cov2(e2, e2 + 15);
  reco::Track::Point v2(0.0178698, 0.0160669, -6.27348);
  reco::Track::Vector p2(-3.96015, -1.08015, -2.2431);
  
  reco::Track trk2(14.5696, 19., v2, p2, -1, cov2);
  reco::TransientTrack tt2(trk2, field);

  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> mu2_p4(p2.x(), p2.y(), p2.z(), muon_mass);
  
  // Photon
  
  // pt(): 10.271
  // eta(): -0.169703
  reco::Photon photon(reco::Candidate::LorentzVector(-6.90274, -7.60566, -1.7514, 10.4193),
		      reco::Photon::Point(-92.7337, -102.121, -29.8284),
		      reco::PhotonCoreRef());
		      
  std::cout << "photon energy: " << photon.energy() << std::endl;
  
  
  //
  // Build KinematicParticle for the Photon
  //
  
  GlobalTrajectoryParameters photon_gtp(GlobalPoint(photon.caloPosition().x(),
						    photon.caloPosition().y(),
						    photon.caloPosition().z()),
					GlobalVector(photon.px(),
						     photon.py(),
						     photon.pz()),
					0, field);

  Matrix3S ph_cov;
  ph_cov(0,0) = pow(0.18,2);
  ph_cov(1,1) = 1;
  ph_cov(2,2) = 1;
  for (unsigned int i = 0; i < 3; ++i){
      for (unsigned int j = 0; j < 3; ++j)
        std::cout << ph_cov(i,j) << "\t";
      std::cout << std::endl;
  }

  auto jac = jacobianSphToCart(photon.p4());

  AlgebraicSymMatrix66 photon_catesian_cov;
  photon_catesian_cov(0,0) = 1.0; // cm
  photon_catesian_cov(1,1) = 1.0; // cm
  photon_catesian_cov(2,2) = 1.0; // cm
  
  Matrix3S ph_cov2 = ROOT::Math::Similarity(jac, ph_cov);
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
        std::cout << ph_cov2(i,j) << "\t";
	photon_catesian_cov(3 + i, 3 + j) = ph_cov2(i,j);
    }
    std::cout << std::endl;
  }

  // CurvilinearTrajectoryError photon_cte;
  // photon_cte.matrix()(0,0) = 1e-6;
  // photon_cte.matrix()(1,1) = 1;
  // photon_cte.matrix()(2,2) = 1;
  // photon_cte.matrix()(3,3) = 1;
  // photon_cte.matrix()(4,4) = 1;
  // // for (unsigned int i = 0; i < 5; ++i){
  // //     for (unsigned int j = 0; j < 5; ++j)
  // //       std::cout << photon_cte.matrix()(i,j) << "\t";
  // //     std::cout << std::endl;
  // // }
  FreeTrajectoryState photon_fts(photon_gtp, CartesianTrajectoryError(photon_catesian_cov));

  std::cout << photon_fts << std::endl;
  
  auto photon_ce = photon_fts.cartesianError();
  for (unsigned int i = 0; i < 6; ++i){
      for (unsigned int j = 0; j < 6; ++j)
        std::cout << photon_ce.matrix()(i,j) << "\t";
      std::cout << std::endl;
  }
  
  // FreeTrajectoryState photon_fts(photon_gtp);

  // Get KinematicState
  // not sure what to put for mass uncertainty
  
  KinematicState photon_kinstate(photon_fts, 0, 1e-6);
  KinematicConstraint* lastConstraint(nullptr);
  ReferenceCountingPointer<KinematicParticle> previousParticle(nullptr);
  float photon_chi2(1);
  float photon_ndf(1);
  ReferenceCountingPointer<KinematicParticle> photon_kp(new VirtualKinematicParticle(photon_kinstate,
                      photon_chi2, photon_ndf, lastConstraint, previousParticle, nullptr));

  std::cout << "Unfitted mass: " << (mu1_p4 + mu2_p4 + photon.p4()).mass() << std::endl;

  
  //
  // Kinimatic Fit
  //

  KinematicParticleFactoryFromTransientTrack factory;
  KinematicParticleVertexFitter fitter;
    
  double chi = 0.;
  double ndf = 0.;
  std::vector<RefCountedKinematicParticle> particles;
  particles.push_back(factory.particle(tt1, muon_mass, chi, ndf, muon_mass_err));
  particles.push_back(factory.particle(tt2, muon_mass, chi, ndf, muon_mass_err));

  std::cout << "Dimuon vertex fit" << std::endl;
  RefCountedKinematicTree vertexFitTree = fitter.fit(particles);

  if ( !vertexFitTree->isValid() ) {
    std::cout << "Fit failed" << std::endl;
    return 1;
  }

  vertexFitTree->movePointerToTheTop();
  dump_result(vertexFitTree);
  auto mm_kp = vertexFitTree->currentParticle();
  // result.refitVertex = vertexFitTree->currentDecayVertex();
  auto mm_state = mm_kp->currentState();
  
  auto mm_p4 = makeLorentzVector(mm_state.globalMomentum(), mm_state.mass());
  
  std::cout << "dimuon fit + photon mass: " << (mm_p4 + photon.p4()).mass() << std::endl;

  
  //
  // 3-body fit
  //

  std::vector<RefCountedKinematicParticle> bfit;
  bfit.push_back(mm_kp);
  bfit.push_back(photon_kp);

  std::cout << "Dimuon vertex fit" << std::endl;
  RefCountedKinematicTree mmg_fit = fitter.fit(bfit);

  if ( !mmg_fit->isValid() ) {
    std::cout << "Fit failed" << std::endl;
    return 1;
  }

  mmg_fit->movePointerToTheTop();
  dump_result(mmg_fit);

  std::cout << "photon_chi2: " << photon_chi2 << std::endl;
  std::cout << "photon_ndf: "  << photon_ndf  << std::endl;


  return 0;

}
