#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelQuality.h"
#include "CondFormats/DataRecord/interface/SiPixelQualityFromDbRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <tuple>

//
// PixelQualityTableProducer stores per-event pixel ROC quality counts
// for each pixel region, computed from the stuckTBM (auto-masked) and
// base (dead) SiPixelQuality conditions payloads. Counts are constant
// within a conditions IOV and are recomputed only when it changes.
//
// Region indexing (fixed 10-region superset, works for every Run-2/Run-3
// geometry): 0-3 BPix L1-L4; 4-6 FPix minus D1-D3; 7-9 FPix plus D1-D3.
// The Phase-0 detector (2016) has only 3 barrel layers and 2 disks/side,
// so its regions 3, 6, 9 are absent and get nRocsTotal == 0; the Phase-1
// detector (2017+, Run 3) populates all 10. ROCs per module are taken from
// the geometry (Phase-1: all 16; Phase-0: 2/5/6/8/10/16), not assumed to be
// 16; every ROC has 80x52 channels in both.
//
// Sentinels (-1) mark "not available", kept distinct from a real 0:
//   * readPixelQuality=False (MC): every column -1.
//   * readStuckTBM=False (data before the stuckTBM PCL tag exists, i.e.
//     all 2016/2017 and early 2018): masked columns -1, dead/total/active
//     still filled from the base payload.
//   * a region absent in the current geometry: nRocsTotal 0, fRocsMasked -1.
//

using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

namespace {
  constexpr int kNPixRegions = 10;

  int pixelRegionIndex(const TrackerTopology& topo, const DetId& id) {
    if (id.det() != DetId::Tracker) return -1;
    if (id.subdetId() == PixelSubdetector::PixelBarrel)
      return topo.pxbLayer(id) - 1;
    if (id.subdetId() == PixelSubdetector::PixelEndcap)
      return 4 + 3 * (topo.pxfSide(id) - 1) + (topo.pxfDisk(id) - 1);
    return -1;
  }
}

class PixelQualityTableProducer : public edm::stream::EDProducer<> {

public:

  explicit PixelQualityTableProducer(const edm::ParameterSet &iConfig);

  ~PixelQualityTableProducer() override {};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void produce(edm::Event&, const edm::EventSetup&);
  void updatePixelQualityCounts(const edm::EventSetup& iSetup);

  // ----------member data ---------------------------
  struct PixQualCounts {
    std::array<int, kNPixRegions> nRocsTotal{}, nRocsMasked{}, nRocsDead{},
      nRocsBadTotal{}, nModulesMasked{};
  };

  const bool readPixelQuality_;
  const bool readStuckTBM_;
  const std::string name_;
  const bool verbose_;

  edm::ESGetToken<SiPixelQuality, SiPixelQualityFromDbRcd> pixQualStuckTBMToken_;
  edm::ESGetToken<SiPixelQuality, SiPixelQualityFromDbRcd> pixQualBaseToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopoToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
  edm::ESWatcher<SiPixelQualityFromDbRcd> pixQualWatcher_;
  PixQualCounts pixQual_;

};

PixelQualityTableProducer::PixelQualityTableProducer(const edm::ParameterSet &iConfig):
  readPixelQuality_( iConfig.getParameter<bool>("readPixelQuality") ),
  readStuckTBM_(     iConfig.getParameter<bool>("readStuckTBM") ),
  name_(            iConfig.getParameter<std::string>("name") ),
  verbose_(         iConfig.getUntrackedParameter<bool>("verbosePixelQuality") )
{
  if (readPixelQuality_) {
    pixQualBaseToken_ = esConsumes<SiPixelQuality, SiPixelQualityFromDbRcd>(
        edm::ESInputTag("", iConfig.getParameter<std::string>("pixelQualityBaseLabel")));
    trackerTopoToken_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
    trackerGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
    if (readStuckTBM_)
      pixQualStuckTBMToken_ = esConsumes<SiPixelQuality, SiPixelQualityFromDbRcd>(
          edm::ESInputTag("", iConfig.getParameter<std::string>("pixelQualityStuckTBMLabel")));
  }
  produces<nanoaod::FlatTable>();
}

void PixelQualityTableProducer::updatePixelQualityCounts(const edm::EventSetup& iSetup) {
  const TrackerTopology& tTopo = iSetup.getData(trackerTopoToken_);
  const TrackerGeometry& geom = iSetup.getData(trackerGeomToken_);
  pixQual_ = PixQualCounts{};

  // Number of ROCs per module is read from the geometry, not assumed to be 16:
  // Phase-1 (2017+, Run 3) modules are all 2x8=16 ROCs, but Phase-0 (2016) has
  // BPix half-modules and FPix plaquettes with 2/5/6/8/10 ROCs. Every ROC has
  // the same 80x52 pixel channel count in both, so ROC counting is uniform.
  std::map<uint32_t, int> rocsPerDet;
  auto addDet = [&](const GeomDet* det) {
    const auto* pdu = dynamic_cast<const PixelGeomDetUnit*>(det);
    if (pdu == nullptr) return;
    const PixelTopology& t = pdu->specificTopology();
    const int nroc = t.rocsX() * t.rocsY();
    const int region = pixelRegionIndex(tTopo, det->geographicalId());
    if (region < 0) return;
    pixQual_.nRocsTotal[region] += nroc;
    rocsPerDet[det->geographicalId().rawId()] = nroc;
  };
  for (const auto* det : geom.detsPXB()) addDet(det);
  for (const auto* det : geom.detsPXF()) addDet(det);

  std::set<std::tuple<int, uint32_t, int>> unionBad;
  auto countOne = [&](const SiPixelQuality& quality, std::array<int, kNPixRegions>& nRocs,
		      std::array<int, kNPixRegions>* nModules) {
    for (const auto& module : quality.getBadComponentList()) {
      int region = pixelRegionIndex(tTopo, DetId(module.DetID));
      if (region < 0) continue;
      auto it = rocsPerDet.find(module.DetID);
      const int nroc = (it != rocsPerDet.end()) ? it->second : 16;
      // errorType==0 means the whole module is bad -> all its ROCs; otherwise
      // BadRocs is a per-ROC mask. Iterate only the module's real ROCs so a
      // whole-bad Phase-0 sub-16 module is not overcounted as 16.
      const uint16_t mask = (module.errorType == 0) ? 0xFFFF : module.BadRocs;
      int nbad = 0;
      for (int roc = 0; roc < nroc; ++roc)
	if ((mask >> roc) & 0x1) {
	  ++nRocs[region];
	  ++nbad;
	  unionBad.emplace(region, module.DetID, roc);
	}
      if (nModules != nullptr && nbad == nroc) ++(*nModules)[region];
    }
  };
  countOne(iSetup.getData(pixQualBaseToken_), pixQual_.nRocsDead, nullptr);
  if (readStuckTBM_)
    countOne(iSetup.getData(pixQualStuckTBMToken_), pixQual_.nRocsMasked, &pixQual_.nModulesMasked);
  for (const auto& badRoc : unionBad)
    ++pixQual_.nRocsBadTotal[std::get<0>(badRoc)];
}

void PixelQualityTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (readPixelQuality_ && pixQualWatcher_.check(iSetup)) {
    updatePixelQualityCounts(iSetup);
    if (verbose_) {
      int nMaskedTotal = 0;
      edm::LogInfo log("PixelQualityTableProducer");
      log << "SiPixelQuality IOV change at " << iEvent.id().run() << ":" << iEvent.luminosityBlock()
	  << (readStuckTBM_ ? "" : " (stuckTBM off)");
      for (int region = 0; region < kNPixRegions; ++region) {
	log << "\n  region " << region << ": total " << pixQual_.nRocsTotal[region]
	    << " masked " << (readStuckTBM_ ? pixQual_.nRocsMasked[region] : -1)
	    << " dead " << pixQual_.nRocsDead[region]
	    << " union " << pixQual_.nRocsBadTotal[region];
	nMaskedTotal += pixQual_.nRocsMasked[region];
      }
      log << "\n  whole-detector masked ROCs: " << (readStuckTBM_ ? nMaskedTotal : -1);
    }
  }

  vector<int32_t> region(kNPixRegions), nRocsTotal(kNPixRegions, -1),
    nRocsMasked(kNPixRegions, -1), nRocsDead(kNPixRegions, -1),
    nRocsBadTotal(kNPixRegions, -1), nRocsActive(kNPixRegions, -1), nModulesMasked(kNPixRegions, -1);
  vector<float> fRocsMasked(kNPixRegions, -1.f);
  for (int i = 0; i < kNPixRegions; ++i) {
    region[i] = i;
    if (readPixelQuality_) {
      nRocsTotal[i] = pixQual_.nRocsTotal[i];
      nRocsDead[i] = pixQual_.nRocsDead[i];
      nRocsBadTotal[i] = pixQual_.nRocsBadTotal[i];
      nRocsActive[i] = pixQual_.nRocsTotal[i] - pixQual_.nRocsBadTotal[i];
      if (readStuckTBM_) {
	nRocsMasked[i] = pixQual_.nRocsMasked[i];
	nModulesMasked[i] = pixQual_.nModulesMasked[i];
	fRocsMasked[i] = (pixQual_.nRocsTotal[i] > 0)
	  ? float(pixQual_.nRocsMasked[i]) / pixQual_.nRocsTotal[i] : -1.f;
      }
    }
  }

  auto out = std::make_unique<nanoaod::FlatTable>(kNPixRegions, name_, false);
  out->setDoc("Pixel ROC quality counts per region from SiPixelQuality conditions (constant within an IOV; -1 = not available: MC, stuckTBM off, or region absent in this geometry)");
  out->addColumn<int32_t>("region", region,
			  "region index: 0-3 BPix L1-L4, 4-6 FPix- D1-D3, 7-9 FPix+ D1-D3 (Phase-0 lacks 3/6/9)");
  out->addColumn<int32_t>("nRocsTotal", nRocsTotal,
			  "total ROCs in region from geometry (0 if region absent in this geometry, -1 if not read)");
  out->addColumn<int32_t>("nRocsMasked", nRocsMasked,
			  "ROCs auto-masked (stuckTBM payload) at this run:ls, permanent part subtracted (-1 if stuckTBM not read)");
  out->addColumn<int32_t>("nRocsDead", nRocsDead,
			  "ROCs bad in the base quality payload (permanent + persistently low-occupancy)");
  out->addColumn<int32_t>("nRocsBadTotal", nRocsBadTotal,
			  "ROC-level union of masked and dead (equals dead when stuckTBM not read)");
  out->addColumn<int32_t>("nRocsActive", nRocsActive,
			  "nRocsTotal minus nRocsBadTotal");
  out->addColumn<int32_t>("nModulesMasked", nModulesMasked,
			  "modules with all 16 ROCs auto-masked (-1 if stuckTBM not read)");
  out->addColumn<float>("fRocsMasked", fRocsMasked,
			"nRocsMasked / nRocsTotal (geometric proxy, not an efficiency; -1 if undefined)");
  iEvent.put(std::move(out));
}

void PixelQualityTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("readPixelQuality", true);
  desc.add<bool>("readStuckTBM", true);
  desc.add<std::string>("pixelQualityStuckTBMLabel", "stuckTBM");
  desc.add<std::string>("pixelQualityBaseLabel", "");
  desc.add<std::string>("name", "pix");
  desc.addUntracked<bool>("verbosePixelQuality", false);
  descriptions.add("pixelQualityTable", desc);
}

DEFINE_FWK_MODULE(PixelQualityTableProducer);
