#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"

 typedef ObjectCountFilter<pat::CompositeCandidateCollection>::type BxFilter;

DEFINE_FWK_MODULE( BxFilter );
