#include "Bmm5/NanoAOD/interface/Candidate.h"
#include "Bmm5/NanoAOD/interface/ScoutingDataHandling.h"

using namespace bmm;

Candidate::Candidate(const pat::Muon& muon, int index):
  reco::LeafCandidate(muon.charge(), muon.polarP4(), muon.vertex(), muon.pdgId(), muon.status()),
  index_(index), name_("mu")
{
  track_ = muon.innerTrack().get();
  gen_particle_ = muon.genParticle();
}

Candidate::Candidate(const Run3ScoutingMuon& smuon, int index):
  index_(index), name_("mu")
{
  embedded_track_ = makeRecoTrack(smuon);
  track_is_embedded_ = true;

  reco::Candidate::PolarLorentzVector p4(smuon.pt(), smuon.eta(), smuon.phi(),  0.10566);
  setP4(p4);
  setCharge(smuon.charge());
  setPdgId(- 13 * charge());
  setVertex(embedded_track_.vertex());

  // gen_particle_ = muon.genParticle();
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

// I think it's not a good idea to take ownership of a reco::Track

// Candidate::Candidate(const Run3ScoutingTrack& strack, int index):
//   index_(index), name_("had")
// {
//   embedded_track_ = makeRecoTrack(strack);
//   track_ = &embedded_track_;

//   reco::Candidate::PolarLorentzVector p4(strack.tk_pt(), strack.tk_eta(), strack.tk_phi(),  0.139570);
//   this->setP4(p4);
//   this->setCharge(strack.tk_charge());
//   this->setPdgId(211 * this->charge());
//   this->setVertex(embedded_track_.vertex());
// }

Candidate::Candidate(const reco::Track& track, int index):
  index_(index), name_("had")
{
  reco::Candidate::PolarLorentzVector p4(track.pt(), track.eta(), track.phi(),  0.139570);
  setP4(p4);
  setCharge(track.charge());
  setPdgId(211 * charge());
  setVertex(track.vertex());
  track_ = &track;
}
