#
# Basic wrapper script to run a processor for a given job
#
# Usage example:
#  python job_starter.py [Processor] [Job file]
#
from FlatNtupleForBmmMva  import FlatNtupleForBmmMva
from FlatNtupleForBmmMvaJpsiK import FlatNtupleForBmmMvaJpsiK
from FlatNtupleForMLFit   import FlatNtupleForMLFit
from FlatNtupleForMuonMVA import FlatNtupleForMuonMVA
from FlatNtupleForTrigEfficiency import FlatNtupleForTrigEfficiency
from FlatNtupleForTrigInfo import FlatNtupleForTrigInfo
from Skimmer import Skimmer
from SimpleSkimmer import SimpleSkimmer
import sys

if len(sys.argv)==3:
    processorType = globals()[sys.argv[1]]
    p = processorType(sys.argv[2])
    p.process()
