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
  
  const float ElectronMass_    = 0.511e-3;
  const float ElectronMassErr_ = 3.5*1e-9;
  const float MuonMass_        = 0.10565837;
  const float MuonMassErr_     = 3.5*1e-9;
  const float KaonMass_        = 0.493677;
  const float KaonMassErr_     = 1.6e-5;
  const float PionMass_        = 0.139570;
  const float PionMassErr_     = 3.5e-7;
  const float JPsiMass_        = 3.0969;
  const float JPsiMassErr_     = 92.9e-6;
  const float Psi2SMass_       = 3.6861;
  const float D0Mass_          = 1.86484;
  const float D0MassErr_       = 0.05e-3;
  const float PhiMass_         = 1.01946;

  // ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
  typedef reco::Candidate::LorentzVector LorentzVector;
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
			     const MagneticField* field);

  // build KinematicParticle for a photon
  KinematicParticleRef build_particle(const reco::Photon& photon,
				      double energyErr,
				      const MagneticField * field,
				      float& chi2,
				      float& ndof);
  struct KalmanVertexFitResult{
    float vtxProb {-1.0};
    bool  valid {false};
    std::vector<LorentzVector> refitVectors;
    GlobalPoint position;
    GlobalError err;
    float lxy{-1.0}, lxyErr{-1.0}, sigLxy{-1.0};
    
    float mass() const;
    void postprocess(const reco::BeamSpot& bs);
  };

  struct DisplacementInformationIn3D{
    double decayLength{-1.0}, decayLengthErr{0.0}, decayLength2{-1.0}, decayLength2Err{0.0}, 
      distaceOfClosestApproach{-1.0}, distaceOfClosestApproachErr{0.0}, distaceOfClosestApproachSig{0.0},
      distaceOfClosestApproach2{-1.0}, distaceOfClosestApproach2Err{0.0}, distaceOfClosestApproach2Sig{0.0},
      longitudinalImpactParameter{0.0}, longitudinalImpactParameterErr{0.0}, longitudinalImpactParameterSig{0.0},
      longitudinalImpactParameter2{0.0}, longitudinalImpactParameter2Err{0.0}, longitudinalImpactParameter2Sig{0.0},
      decayTime{-999.}, decayTimeError{-999.}, decayTimeXY{-999.}, decayTimeXYError{-999.},
      alpha{-999.}, alphaErr{-999.}, alphaXY{-999.}, alphaXYErr{-999.};
    const reco::Vertex *pv{nullptr}, *pv2{nullptr};
    int pvIndex{-1}, pv2Index{-1};
  };

  LorentzVector makeLorentzVectorFromPxPyPzM(double px, double py, double pz, double m);

  
}

#endif
