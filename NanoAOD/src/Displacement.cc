#include "Bmm5/NanoAOD/interface/Displacement.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TVector3.h"

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

void bmm::Displacement::compute_displacement()
{
  if (not fit_.valid()) return;
  auto candTransientTrack = fit_.particle()->refittedTransientTrack();
  auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, prodVertex_);
  auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), prodVertex_);
  
  if (impactParameterZ.first) {
    longitudinalImpactParameter_    = impactParameterZ.second.value();
    longitudinalImpactParameterSig_ = impactParameterZ.second.significance();
    longitudinalImpactParameterErr_ = impactParameterZ.second.error();
  }

  if (impactParameter3D.first and not isnan(impactParameter3D.second.error())) {
    distaceOfClosestApproach_       = impactParameter3D.second.value();
    distaceOfClosestApproachSig_    = impactParameter3D.second.significance();
    distaceOfClosestApproachErr_    = impactParameter3D.second.error();
  }

  // compute decay length
  VertexDistance3D distance3D;
  auto dist = distance3D.distance(prodVertex_, fit_.vtx_state() );
  decayLength_    = dist.value();
  decayLengthErr_ = dist.error();
  
  VertexDistanceXY distanceXY;
  auto distXY = distanceXY.distance(prodVertex_, fit_.vtx_state() );

  //
  // Pointing angle
  //
  auto alpha = getAlpha(fit_.vtx_position(),
			fit_.vtx_error(),
			GlobalPoint(Basic3DVector<float>(prodVertex_.position())),
			GlobalError(prodVertex_.covariance()),
			fit_.p3());

  auto alphaXY = getAlpha(fit_.vtx_position(),
			  fit_.vtx_error(),
			  GlobalPoint(Basic3DVector<float>(prodVertex_.position())),
			  GlobalError(prodVertex_.covariance()),
			  fit_.p3(),
			  true);

  alpha_    = alpha.first;
  alphaErr_ = alpha.second;

  alphaXY_    = alphaXY.first;
  alphaXYErr_ = alphaXY.second;

  
  //
  // Decay time information
  //
  TVector3 plab(fit_.p3().x(), fit_.p3().y(), fit_.p3().z());
  const double massOverC = fit_.mass()/TMath::Ccgs();

  // get covariance matrix for error propagation in decayTime calculation
  auto vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(prodVertex_.error()),
					     fit_.particle()->currentState().kinematicParametersError().matrix());
  auto vtxDistanceJac3d = makeJacobianVector3d(prodVertex_.position(), fit_.vtx_position(), plab);
  auto vtxDistanceJac2d = makeJacobianVector2d(prodVertex_.position(), fit_.vtx_position(), plab);

  decayTime_ = dist.value() / plab.Mag() * cos(alpha_) * massOverC;
  decayTimeError_ = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;

  decayTimeXY_ = distXY.value() / plab.Perp() * cos(alphaXY_) * massOverC;
  decayTimeXYError_ = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;
}
