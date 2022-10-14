#include "Bmm5/NanoAOD/interface/LeptonCandidate.h"


using namespace bmm;

const reco::Track*
LeptonCandidate::track() const
{
  if (muon_) return muon_->innerTrack().get();
  if (hadron_) return hadron_->bestTrack();
  if (electron_) return electron_->bestTrack();
  return nullptr;
}

LeptonCandidate::LeptonCandidate(const pat::Muon& muon, int index):
  reco::LeafCandidate(muon.charge(), muon.polarP4(), muon.vertex(), muon.pdgId(), muon.status()),
  index_(index), muon_(&muon), name_("mu")
{}

LeptonCandidate::LeptonCandidate(const pat::Electron& elec, int index):
  reco::LeafCandidate(elec.charge(), elec.polarP4(), elec.vertex(), elec.pdgId(), elec.status()),
  index_(index), electron_(&elec), name_("el")
{}
    
LeptonCandidate::LeptonCandidate(const pat::PackedCandidate& hadron, bool from_gen):
  reco::LeafCandidate(hadron.charge(), hadron.polarP4(), hadron.vertex(), hadron.pdgId(), hadron.status()),
  gen_(from_gen), hadron_(&hadron), name_("had")
{}
