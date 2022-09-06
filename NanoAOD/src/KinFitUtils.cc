#include "Bmm5/NanoAOD/interface/KinFitUtils.h"

using namespace bmm;

Matrix33
bmm::jacobianSphToCart(const reco::Candidate::LorentzVector& p4)
{
  Matrix33 jacobian;
  
  jacobian(0,0) =   p4.Px() / p4.P();
  jacobian(0,1) =   p4.Py() / p4.P();
  jacobian(0,2) =   p4.Pz() / p4.P();
  jacobian(1,0) =   p4.Px() * p4.Pz() / p4.Pt() / p4.P() / p4.P();
  jacobian(1,1) =   p4.Py() * p4.Pz() / p4.Pt() / p4.P() / p4.P();
  jacobian(1,2) = - p4.Pt() / p4.P() / p4.P();
  jacobian(2,0) = - p4.Py() / p4.Pt() / p4.Pt();
  jacobian(2,1) =   p4.Px() / p4.Pt() / p4.Pt();
  jacobian(2,2) =   0;
  
  return jacobian;
}

FreeTrajectoryState
bmm::getFTS(const reco::Photon& photon,
	    double energyErr,
	    const MagneticField * field)
{
  GlobalTrajectoryParameters gtp(GlobalPoint(photon.caloPosition().x(),
					     photon.caloPosition().y(),
					     photon.caloPosition().z()),
				 GlobalVector(photon.px(),
					      photon.py(),
					      photon.pz()),
				 0, field);
  
  Matrix3S momentum_covariance_spherical;

  // Energy resolution
  momentum_covariance_spherical(0,0) = pow(energyErr, 2);

  // Direction uncertainty
  // - without production vertix photon direction is not known
  // - take the default direction and put large uncertainty
  momentum_covariance_spherical(1,1) = 1; 
  momentum_covariance_spherical(2,2) = 1;
  
  // for (unsigned int i = 0; i < 3; ++i){
  //     for (unsigned int j = 0; j < 3; ++j)
  //       std::cout << momentum_covariance_spherical(i,j) << "\t";
  //     std::cout << std::endl;
  // }

  auto jac = jacobianSphToCart(photon.p4());
  
  AlgebraicSymMatrix66 covariance_catesian;
  covariance_catesian(0,0) = 1.0; // cm
  covariance_catesian(1,1) = 1.0; // cm
  covariance_catesian(2,2) = 1.0; // cm
  
  Matrix3S momentum_covariance_spherical2 = ROOT::Math::Similarity(jac, momentum_covariance_spherical);
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      // std::cout << momentum_covariance_spherical2(i,j) << "\t";
      covariance_catesian(3 + i, 3 + j) = momentum_covariance_spherical2(i,j);
    }
    // std::cout << std::endl;
  }

  return FreeTrajectoryState(gtp, CartesianTrajectoryError(covariance_catesian));
}

KinematicParticleRef
bmm::build_particle(const reco::Photon& photon,
	       double energyErr,
	       const MagneticField * field,
	       float& chi2,
	       float& ndof)
{
  auto fts(getFTS(photon, energyErr, field));
  KinematicState kinstate(fts, 0, 1e-6);
  KinematicConstraint* last_constraint(nullptr);
  KinematicParticleRef previous_particle(nullptr);
  return KinematicParticleRef(new VirtualKinematicParticle(kinstate, chi2, ndof,
							   last_constraint, previous_particle,
							   nullptr));
}
