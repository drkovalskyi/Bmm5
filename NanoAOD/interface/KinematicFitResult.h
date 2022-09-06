#ifndef Bmm5_NanoAOD_KinematicFitResult_h
#define Bmm5_NanoAOD_KinematicFitResult_h

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

std::pair<float, float>
getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
	 const GlobalPoint& ip_position,  const GlobalError& ip_error,
			     const GlobalVector &momentum,
	 bool transverse = true);

struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, alphaBS, alphaBSErr;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), lxyErr(-1.0), sigLxy(-1.0),
		       alphaBS(-999.), alphaBSErr(-999.)
  {}

  bool valid() const;
  void postprocess(const reco::BeamSpot& beamSpot);
  float mass() const;
  float refit_mass(unsigned int i, unsigned int j) const;
  GlobalVector p3() const;
  GlobalVector dau_p3(unsigned int i) const;
  float massErr() const;
  float chi2() const;
  float ndof() const;
  float vtxProb() const;
};

#endif
