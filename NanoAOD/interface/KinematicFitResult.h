#ifndef Bmm5_NanoAOD_KinematicFitResult_h
#define Bmm5_NanoAOD_KinematicFitResult_h

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

std::pair<float, float>
getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
	 const GlobalPoint& ip_position,  const GlobalError& ip_error,
			     const GlobalVector &momentum,
	 bool transverse = true);

typedef reco::Candidate::LorentzVector LorentzVector;

class KinematicFitResult{
 public:
  KinematicFitResult() : lxy_(-1.0), lxyErr_(-1.0), sigLxy_(-1.0),
    alphaBS_(-999.), alphaBSErr_(-999.),
    treeIsValid(false)
  {}

  void set_tree(RefCountedKinematicTree tree);
  
  bool valid() const;

  // compute displacement with respect to the beam spot
  void postprocess(const reco::BeamSpot& beamSpot);
  
  float mass() const;
  float refit_mass(unsigned int i, unsigned int j) const;
  GlobalVector p3() const;
  LorentzVector p4() const;
  unsigned int number_of_daughters() const
  {
    return refitDaughters.size();
  }
  GlobalVector dau_p3(unsigned int i) const;
  float dau_mass(unsigned int i) const;
  float massErr() const;
  float chi2() const;
  float ndof() const;
  float vtxProb() const;
  GlobalPoint vtx_position() const;
  VertexState vtx_state() const;
  GlobalError vtx_error() const;
  reco::Vertex vertex() const;
  float sumPt() const;
  float sumPt2() const;
  float lxy() const {return lxy_;}
  float lxyErr() const {return lxyErr_;}
  float sigLxy() const {return sigLxy_;}
  float alphaBS() const {return alphaBS_;}
  float alphaBSErr() const {return alphaBSErr_;}
  RefCountedKinematicTree tree(){return refitTree;}
  const RefCountedKinematicParticle particle() const {return refitMother;}

  std::vector<const reco::Track*> tracks;

 private:
  float lxy_, lxyErr_, sigLxy_, alphaBS_, alphaBSErr_;
  bool treeIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
};

#endif
