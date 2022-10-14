#ifndef Bmm5_NanoAOD_LeptonCandidate_h
#define Bmm5_NanoAOD_LeptonCandidate_h
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


// LeptonCanidate is a container to hold a muon, electron or a hadron
// that may decay to a muon. Most of the functionallity is available
// via various Candidate types. The class is used mostly to facilitate
// some non-trivial use cases

// index is the position of the corresponding original muon or
// electron in their collection. For hadrons Index is -1

namespace bmm
{
  class LeptonCandidate: public reco::LeafCandidate{
  public:
    LeptonCandidate(const pat::Muon& muon, int index);
    LeptonCandidate(const pat::Electron& elec, int index);
    LeptonCandidate(const pat::PackedCandidate& hadron, bool from_gen);

    const reco::Track* track() const;

    const pat::Muon*            muon() const     {return muon_;}
    const pat::Electron*        electron() const {return electron_;}
    const pat::PackedCandidate* hadron() const   {return hadron_;}
    
    int index() const { return index_; }
    bool from_gen() const { return gen_; }
    const std::string& name() const { return name_; }
      
  private:
    int index_{-1};
    bool gen_{false};
    const pat::Muon *muon_{nullptr};
    const pat::PackedCandidate *hadron_{nullptr};
    const pat::Electron *electron_{nullptr};
    std::string name_;
  };
  
}
#endif
