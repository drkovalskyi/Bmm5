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

float KalmanVertexFitResult::mass() const
{
  if (not valid) return -1.0;
  LorentzVector p4;
  for (auto v: refitVectors)
    p4 += v;
  return p4.mass();
}
  
void KalmanVertexFitResult::postprocess(const reco::BeamSpot& bs)
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

LorentzVector
bmm::makeLorentzVectorFromPxPyPzM(double px, double py, double pz, double m){
  double p2 = px*px+py*py+pz*pz;
  return LorentzVector(px,py,pz,sqrt(p2+m*m));
}

LorentzVector
bmm::makeLorentzVectorFromP3M(const math::XYZVector& p, double m){
  return LorentzVector(p.x(), p.y(), p.z(), sqrt(p.Mag2() + m * m));
}

LorentzVector
bmm::makeLorentzVectorFromP3M(const GlobalVector& p, double m){
  return LorentzVector(p.x(), p.y(), p.z(), sqrt(p.mag2() + m * m));
}
