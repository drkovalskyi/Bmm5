#include "Bmm5/NanoAOD/interface/Candidate.h"


using namespace bmm;

Candidate::Candidate(const pat::Muon& muon, int index):
  reco::LeafCandidate(muon.charge(), muon.polarP4(), muon.vertex(), muon.pdgId(), muon.status()),
  index_(index), name_("mu")
{
  track_ = muon.innerTrack().get();
  gen_particle_ = muon.genParticle();
}

Candidate::Candidate(const pat::Electron& elec, int index):
  reco::LeafCandidate(elec.charge(), elec.polarP4(), elec.vertex(), elec.pdgId(), elec.status()),
  index_(index), name_("el")
{
  track_ = elec.bestTrack();
  gen_particle_ = elec.genParticle();
}
    
Candidate::Candidate(const pat::PackedCandidate& hadron, const pat::PackedGenParticle* gen):
  reco::LeafCandidate(hadron.charge(), hadron.polarP4(), hadron.vertex(), hadron.pdgId(), hadron.status()),
  packed_gen_particle_(gen), name_("had")
{
  track_ = hadron.bestTrack();
}


void Candidate::setType(double mass, std::string name, int pdgId){
  this->setMass(mass);
  name_ = name;
  if (pdgId != 0)
    this->setPdgId(pdgId);
}
