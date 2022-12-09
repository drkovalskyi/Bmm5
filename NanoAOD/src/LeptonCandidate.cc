#include "Bmm5/NanoAOD/interface/LeptonCandidate.h"


using namespace bmm;

LeptonCandidate::LeptonCandidate(const pat::Muon& muon, int index):
  reco::LeafCandidate(muon.charge(), muon.polarP4(), muon.vertex(), muon.pdgId(), muon.status()),
  index_(index), name_("mu")
{
  track_ = muon.innerTrack().get();
  gen_particle_ = muon.genParticle();
}

LeptonCandidate::LeptonCandidate(const pat::Electron& elec, int index):
  reco::LeafCandidate(elec.charge(), elec.polarP4(), elec.vertex(), elec.pdgId(), elec.status()),
  index_(index), name_("el")
{
  track_ = elec.bestTrack();
  gen_particle_ = elec.genParticle();
}
    
LeptonCandidate::LeptonCandidate(const pat::PackedCandidate& hadron, const pat::PackedGenParticle* gen):
  reco::LeafCandidate(hadron.charge(), hadron.polarP4(), hadron.vertex(), hadron.pdgId(), hadron.status()),
  packed_gen_particle_(gen), name_("had")
{
  track_ = hadron.bestTrack();
}


void LeptonCandidate::setType(double mass, std::string name){
  this->setMass(mass);
  name_ = name;
}
