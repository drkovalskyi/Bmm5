#ifndef Bmm5_NanoAOD_CommonTools_h
#define Bmm5_NanoAOD_CommonTools_h
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

namespace bmm
{
  typedef reco::Candidate::LorentzVector LorentzVector;
  struct CloseTrack{
    float svDoca{-1.0}, svDocaErr{-1.0}, svProb{-1.0},
      pvDoca{-1.0}, pvDocaErr{-1.0},
      impactParameterSignificanceBS{-1.0};
    const pat::PackedCandidate* pfCand{nullptr};
  };

  struct CloseTrackInfo{
    std::vector<CloseTrack> tracks;
    unsigned int nTracksByVertexProbability(double minProb = 0.1, 
					    double minIpSignificance = -1,
					    int pvIndex = -1,
					    const pat::PackedCandidate* ignoreTrack1 = 0);
    unsigned int nTracksByDisplacementSignificance(double max_svDoca = 0.03, 
						   double maxSignificance = -1,
						   int pvIndex = -1,
						   const pat::PackedCandidate* ignoreTrack1 = 0);
    unsigned int nTracksByBetterMatch(double max_svDoca = 0.03, 
				      double maxSignificance = 2,
				      int pvIndex = -1,
				      const pat::PackedCandidate* ignoreTrack1 = 0);
    float minDoca(double max_svDoca = 0.03, 
		  int pvIndex = -1,
		  const pat::PackedCandidate* ignoreTrack1 = 0);
    void fillCandInfo(pat::CompositeCandidate& cand, int pvIndex, std::string name);
  };


  bool dr_match(const LorentzVector& reco , const LorentzVector& gen);
  
  std::vector<unsigned int> 
  get_depth_from_permutation(const std::vector<unsigned int>& elements);

  bool is_acceptable(const reco::Candidate* cand);

  // depth 0 - first mother
  const reco::Candidate* get_mother(const reco::Candidate* cand, unsigned int depth);

  const reco::Candidate* 
    find_common_ancestor(const std::vector<const reco::Candidate*>& particles, 
			 unsigned int max_depth=10);
  
}

#endif
