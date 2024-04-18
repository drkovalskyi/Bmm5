#ifndef Bmm5_NanoAOD_Candidate_h
#define Bmm5_NanoAOD_Candidate_h
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"

// Canidate is a container to hold a muon, electron or a hadron. The
// class is used to simplify code reuse
//
// index is the position of the corresponding original muon or
// electron in their collection. For hadrons Index is -1
//
// name should be "mu", "el" or "had", but it's not enforced.

namespace bmm
{
  class Candidate: public reco::LeafCandidate{
  public:
    Candidate(const pat::Muon& muon, int index);
    Candidate(const Run3ScoutingMuon& muon, int index);
    // Candidate(const Run3ScoutingTrack& track, int index);
    Candidate(const pat::Electron& elec, int index);
    Candidate(const pat::PackedCandidate& hadron, const pat::PackedGenParticle* = nullptr);
    Candidate(const reco::Track& track, int index);

    int index() const { return index_; }
    bool from_gen() const { return packed_gen_particle_ != nullptr; }
    const std::string& name() const { return name_; }
    const reco::Track* track() const{
      if (track_is_embedded_)
	return &embedded_track_;
      else
	return track_;
    }
    const reco::Track* bestTrack() const { return track(); }
    const reco::GenParticle* genParticle() const { return gen_particle_; }
    const pat::PackedGenParticle* packedGenParticle() const { return packed_gen_particle_; }
    void setType(double mass, std::string name, int pdgId = 0);
      
  private:
    int index_{-1};
    const reco::GenParticle* gen_particle_{nullptr};
    const pat::PackedGenParticle* packed_gen_particle_{nullptr};
    const reco::Track* track_{nullptr};
    reco::Track embedded_track_;
    bool track_is_embedded_ = false;
    std::string name_;
  };
  
}
#endif
