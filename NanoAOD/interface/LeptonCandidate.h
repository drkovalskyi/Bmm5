#ifndef Bmm5_NanoAOD_LeptonCandidate_h
#define Bmm5_NanoAOD_LeptonCandidate_h
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


// LeptonCanidate is a container to hold a muon, electron or a hadron
// that may represent a lepton. Most of the functionallity is available
// via various Candidate types. The class is used mostly to treat all 
// different types of jobs in the same way.

// index is the position of the corresponding original muon or
// electron in their collection. For hadrons Index is -1

// name should be "mu", "el" or "had", but it's not enforced.

namespace bmm
{
  class LeptonCandidate: public reco::LeafCandidate{
  public:
    LeptonCandidate(const pat::Muon& muon, int index);
    LeptonCandidate(const pat::Electron& elec, int index);
    LeptonCandidate(const pat::PackedCandidate& hadron, const pat::PackedGenParticle* = nullptr);

    int index() const { return index_; }
    bool from_gen() const { return packed_gen_particle_ != nullptr; }
    const std::string& name() const { return name_; }
    const reco::Track* track() const{ return track_; }
    const reco::GenParticle* genParticle() const { return gen_particle_; }
    const pat::PackedGenParticle* packedGenParticle() const { return packed_gen_particle_; }
    void setType(double mass, std::string name);
      
  private:
    int index_{-1};
    const reco::GenParticle* gen_particle_{nullptr};
    const pat::PackedGenParticle* packed_gen_particle_{nullptr};
    const reco::Track* track_{nullptr};
    std::string name_;
  };
  
}
#endif
