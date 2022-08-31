#ifndef Bmm5_NanoAOD_KinFitUtils_h
#define Bmm5_NanoAOD_KinFitUtils_h

// #include "RecoVertex/KinematicFitPrimitives/interface/TrackKinematicStatePropagator.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
// #include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h"

// #include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
// #include "TrackingTools/TrajectoryState/interface/TrajectoryStateAccessor.h"
// #include "DataFormats/GeometrySurface/interface/Surface.h" 
// #include "DataFormats/GeometrySurface/interface/BoundPlane.h"
// #include "MagneticField/Engine/interface/MagneticField.h"

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

namespace bmm
{
  
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LorentzVectorM;
  typedef ROOT::Math::SMatrix<double, 6, 6, ROOT::Math::MatRepSym<double, 6> > Matrix6S;
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > Matrix3S;
  typedef ROOT::Math::SMatrix<double, 3, 3> Matrix33;
  typedef ReferenceCountingPointer<KinematicParticle> KinematicParticleRef;
  
  // Jacobian for momentum transformation from the spherical coordinate
  // system to the cartesian one
  Matrix33 jacobianSphToCart(const reco::Candidate::LorentzVector& p4);

  // create trajectory state for a photon
  FreeTrajectoryState getFTS(const reco::Photon& photon,
			     double energyErr,
			     MagneticField* field);

  // build KinematicParticle for a photon
  KinematicParticleRef build_particle(const reco::Photon& photon,
				      double energyErr,
				      MagneticField * field,
				      float& chi2,
				      float& ndof);
}

#endif
