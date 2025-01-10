#include "Bmm5/NanoAOD/interface/KinematicFitResult.h"
#include "Bmm5/NanoAOD/interface/KinFitUtils.h"
#include <TVector.h>
#include <TMatrix.h>
#include <TMath.h>


std::pair<float, float>
getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
	 const GlobalPoint& ip_position,  const GlobalError& ip_error,
			     const GlobalVector &momentum,
			     bool transverse){
  AlgebraicSymMatrix33 error_matrix(vtx_error.matrix() + ip_error.matrix());
  GlobalVector dir(vtx_position - ip_position);
  if (dir.mag() == 0)
    return std::pair<float, float>(999., 999.);

  GlobalVector p(momentum);
  if (transverse){
    dir = GlobalVector(dir.x(), dir.y(), 0);
    p = GlobalVector(p.x(), p.y(), 0);
  }
  
  double dot_product = dir.dot(p);
  double cosAlpha = dot_product / p.mag() / dir.mag();
  if (cosAlpha > 1) cosAlpha = 1;
  if (cosAlpha < -1) cosAlpha = -1;

  // Error propagation

  double c1 = 1 / dir.mag() / p.mag();
  double c2 = dot_product / pow(dir.mag(), 3) / p.mag();
  
  double dfdx = p.x() * c1 - dir.x() * c2;
  double dfdy = p.y() * c1 - dir.y() * c2;
  double dfdz = p.z() * c1 - dir.z() * c2;

  double err2_cosAlpha =
    pow(dfdx, 2) * error_matrix(0, 0) +
    pow(dfdy, 2) * error_matrix(1, 1) +
    pow(dfdz, 2) * error_matrix(2, 2) +
    2 * dfdx * dfdy * error_matrix(0, 1) +
    2 * dfdx * dfdz * error_matrix(0, 2) +
    2 * dfdy * dfdz * error_matrix(1, 2);

  float err_alpha = fabs(cosAlpha) <= 1 and err2_cosAlpha >=0 ? sqrt(err2_cosAlpha) / sqrt(1-pow(cosAlpha, 2)) : 999;
  float alpha = acos(cosAlpha);
  if (isnan(alpha) or isnan(err_alpha))
    return std::pair<float, float>(999., 999.);
  else
    return std::pair<float, float>(alpha, err_alpha);
}

bool KinematicFitResult::valid() const {
  return treeIsValid and refitVertex->vertexIsValid();
}

void KinematicFitResult::postprocess(const reco::BeamSpot& beamSpot)
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
  
  lxy_ = sqrt(v.Norm2Sqr());
  lxyErr_ = sqrt( v*(errVtx*v) + v*(errBS*v) ) / lxy_;
  if (lxyErr_ > 0) sigLxy_ = lxy_/lxyErr_;
    
  // compute pointing angle wrt BeamSpot (2D)

  // rotatedCovariance3D - is a proper covariance matrix for Beam Spot,
  // which includes the beam spot width, not just uncertainty on the
  // absolute beamspot position
  auto alphaXY = getAlpha(refitVertex->vertexState().position(),
			  refitVertex->vertexState().error(),
			  GlobalPoint(Basic3DVector<float>(beamSpot.position())),
			  GlobalError(beamSpot.rotatedCovariance3D()),
			  refitMother->currentState().globalMomentum(),
			  true);
  alphaBS_    = alphaXY.first;
  alphaBSErr_ = alphaXY.second;
}
  
float KinematicFitResult::mass() const
{
  if ( not valid() ) return -1.0;
  return refitMother->currentState().mass();
}

float KinematicFitResult::refit_mass(unsigned int i, unsigned int j) const
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

GlobalVector KinematicFitResult::p3() const
{
  if ( not valid() ) return GlobalVector();
  return refitMother->currentState().globalMomentum();
}

LorentzVector KinematicFitResult::p4() const
{
  if ( not valid() ) return LorentzVector();
  return bmm::makeLorentzVectorFromP3M(p3(), mass());
}


GlobalVector KinematicFitResult::dau_p3(unsigned int i) const
{
  if ( not valid() or i>=refitDaughters.size() ) return GlobalVector();
  return refitDaughters.at(i)->currentState().globalMomentum();
}

float KinematicFitResult::massErr() const
{
  if ( not valid() ) return -1.0;
  return sqrt(refitMother->currentState().kinematicParametersError().matrix()(6,6));
}

float KinematicFitResult::chi2() const
{
  if ( not valid() ) return -1.0;
  return refitVertex->chiSquared();
}

float KinematicFitResult::ndof() const
{
  return refitVertex->degreesOfFreedom();
}

float KinematicFitResult::vtxProb() const
{
  if ( not valid() ) return -1.0;
  return TMath::Prob((double)refitVertex->chiSquared(), int(rint(refitVertex->degreesOfFreedom())));
}

float KinematicFitResult::sumPt() const
{
  float result(0);
  for (const auto& dau: refitDaughters){
    result += dau->currentState().globalMomentum().mag();
  }
  return result;
}

float KinematicFitResult::sumPt2() const
{
  float result(0);
  for (const auto& dau: refitDaughters){
    result += dau->currentState().globalMomentum().mag2();
  }
  return result;
}

reco::Vertex KinematicFitResult::vertex() const
{
//   if (not valid()) return reco::Vertex();
//   return *refitVertex;
// }
  if (not valid()) return reco::Vertex();

  // Extract position, error, and chi-squared from the KinematicVertex
  const auto& position_global = refitVertex->position();
  math::XYZPoint position(position_global.x(), position_global.y(), position_global.z()); // Convert GlobalPoint to math::XYZPoint
  const auto& error = refitVertex->error().matrix(); // Convert GlobalError to AlgebraicSymMatrix33
  double chi2 = refitVertex->chiSquared();
  double ndof = refitVertex->degreesOfFreedom();
  
  return reco::Vertex(position, error, chi2, ndof, 0);
}

void KinematicFitResult::set_tree(RefCountedKinematicTree tree)
{
  if (not tree->isValid()) return;

  treeIsValid = true;

  // extract the re-fitted tracks  
  if ( tree->movePointerToTheFirstChild() ){
    do {
      refitDaughters.push_back(tree->currentParticle());
    } while (tree->movePointerToTheNextChild());
  }

  tree->movePointerToTheTop();
  refitVertex = tree->currentDecayVertex();
  refitMother = tree->currentParticle();
  refitTree   = tree;
}

GlobalPoint KinematicFitResult::vtx_position() const
{
  if (not valid()) return GlobalPoint();
  return refitVertex->vertexState().position();
}

GlobalError KinematicFitResult::vtx_error() const
{
  if (not valid()) return GlobalError();
  return refitVertex->vertexState().error();
}

VertexState KinematicFitResult::vtx_state() const
{
  if (not valid()) return VertexState();
  return refitVertex->vertexState();
}
