import ROOT
import numpy as np
import glob

ROOT.ROOT.EnableImplicitMT()

# Select Events from flat ntuples
flat_files = [
    "/eos/user/a/alen/SWAN_projects/MIT_Ksmumu/Step3/Kspipi_same_sign_IcnlusiveDilepton_full_DNN_max.root"
]
rdf_flat = ROOT.RDataFrame("kspipiMc", flat_files)
rdf_flat = rdf_flat.Filter("my_dnn_val>0.999&&d1mc_pdgId!=0&&d2mc_pdgId!=0")
np_flat = rdf_flat.AsNumpy(["run", "ls", "evt"])  # materializes once
keys = set(zip(np_flat["run"].astype(np.int64),
               np_flat["ls"].astype(np.int64),
               np_flat["evt"].astype(np.int64)))
print(f"Number of events selected: {len(keys)}")

# Find all NanoAOD files
path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/534/"
nano_files = glob.glob(f"{path}/InclusiveDileptonMinBias*/11*.root")
print(f"Number of NanoAOD files to process: {len(nano_files)}")

# build one C++ function with a static unordered_set of constants
triples_src = ",\n    ".join(
    f"Triple{{{r}ULL, {l}ULL, {e}ULL}}" for r, l, e in keys
)

cpp = f"""
#include <unordered_set>
#include <cstdint>
#include <vector>

namespace KeyDB {{

struct Triple {{
  ULong64_t run;
  ULong64_t ls;
  ULong64_t evt;
  bool operator==(const Triple& o) const noexcept {{
    return run==o.run && ls==o.ls && evt==o.evt;
  }}
}};

struct TripleHash {{
  size_t operator()(const Triple& t) const noexcept {{
    auto h1 = std::hash<ULong64_t>{{}}(t.run);
    auto h2 = std::hash<ULong64_t>{{}}(t.ls);
    auto h3 = std::hash<ULong64_t>{{}}(t.evt);
    return h1 ^ (h2<<1) ^ (h3<<2);
  }}
}};

// build once at library load, no first-call lock
static const std::unordered_set<Triple, TripleHash> KS = {{
  {triples_src}
}};

inline bool InKeys(ULong64_t run, ULong64_t ls, ULong64_t evt) {{
  return KS.find(Triple{{run, ls, evt}}) != KS.end();
}}

}} // namespace KeyDB
"""
ROOT.gInterpreter.Declare(cpp)


# Pick events from NanoAOD

rdf_nano = ROOT.RDataFrame("Events", nano_files)
rdf_nano = rdf_nano.Filter("KeyDB::InKeys(run, luminosityBlock, event)")

file = "/tmp/dmytro/selected_events.root"
rdf_nano.Snapshot("Events", file)
print(f"Number of events found in NanoAOD: {rdf_nano.Count().GetValue()}")
