from CRABClient.UserUtilities import config
config = config()

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

##
##  General
##

# A name the user gives to it's request/task. Used by CRAB to create a
# project directory (named crab_<requestName>) where files
# corresponding to this particular task will be stored. Defaults to
# <time-stamp>, where the time stamp is of the form
# <YYYYMMDD>_<hhmmss> and corresponds to the submission time. The
# maximum allowed length is 100 characters matching 'a-zA-Z0-9\-_:'
# pattern.
config.General.requestName = 'QCD_Pt-30toInf_BmmGenFilter_v10'

# The area (full or relative path) where to create the CRAB project
# directory. If the area doesn't exist, CRAB will try to create it
# using the mkdir command. Defaults to the current working directory.
config.General.workArea = '../../../crab_privateMCProduction'

# Whether or not to transfer the output files to the storage site. If
# set to False, the output files are discarded and the user can not
# recover them. Defaults to True.
config.General.transferOutputs = True

# Whether or not to copy the jobs log files to the storage site. If
# set to False, the log files are discarded and the user can not
# recover them. Notice however that a short version of the log files
# containing the first 1000 lines and the last 3000 lines are still
# available through the monitoring web pages. Defaults to False.
config.General.transferLogs = False

##
## JobType
##

# Specifies if this task is running an analysis ('Analysis') on an
# existing dataset or is running MC event generation ('PrivateMC').
config.JobType.pluginName = 'PrivateMC'

# Whether to disable or not the automatic recognition of output files
# produced by PoolOutputModule or TFileService in the CMSSW
# parameter-set configuration. If set to True, it becomes the user's
# responsibility to specify in the JobType.outputFiles parameter all
# the output files that need to be collected. Defaults to False.
config.JobType.disableAutomaticOutputCollection = False

# List of output files that need to be collected. If
# disableAutomaticOutputCollection = False (the default), output files
# produced by PoolOutputModule or TFileService in the CMSSW
# parameter-set configuration are automatically recognized by CRAB and
# don't need to be included in this parameter.
config.JobType.outputFiles = []

# Maximum amount of memory (in MB) a job is allowed to use. Defaults
# to 2000.  
config.JobType.maxMemoryMB = 2000

# The name of the CMSSW parameter-set configuration file that should
# be run via cmsRun. Defaults to 'pset.py'.
# config.JobType.psetName = 'fake.py' 
config.JobType.psetName = 'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8_fake.py' 

# List of private input files (and/or directories) needed by the
# jobs. They will be added to the input sandbox. The input sandbox can
# not exceed 120 MB. The input sandbox is shipped with each job. The
# input files will be placed in the working directory where the users'
# application (e.g. cmsRun) is launched regardless of a possible path
# indicated in this parameter (i.e. only the file name at right of
# last / is relevant). Directories are tarred and their subtree
# structure is preserved. 
config.JobType.inputFiles = [
    'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIFall18GS-RunIIAutumn18DRPremix-RunIIAutumn18MiniAOD_scriptExe.sh',
    'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIFall18GS_cfg.py',
    'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-DIGI_cfg.py',
    'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18DRPremix-RECO_cfg.py',
    'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIAutumn18MiniAOD_cfg.py'
]

# A user script that should be run on the worker node instead of the
# default cmsRun.
config.JobType.scriptExe = 'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8-RunIIFall18GS-RunIIAutumn18DRPremix-RunIIAutumn18MiniAOD_scriptExe.sh'

# Number of requested cores per job. Defaults to 1. If you increase
# this value to run multi-threaded cmsRun, you may need to increase
# maxMemoryMB as well
config.JobType.numCores = 1

config.JobType.allowUndistributedCMSSW = True

##
## Data 
##

# When JobType.pluginName = 'PrivateMC', the splitting mode can only
# be 'EventBased'.
config.Data.splitting = 'EventBased'

# Mandatory when Data.splitting is not 'Automatic', suggests (but not
# impose) how many units (i.e. files, luminosity sections or events -
# depending on the splitting mode - see the note about Data.splitting
# below) to include in each job.
config.Data.unitsPerJob = 10
# Mandatory when JobType.pluginName = 'PrivateMC', in which case the
# parameter tells how many events to generate in total.
config.Data.totalUnits = 100000

# The first part of the LFN of the output files. Accepted values are
# /store/user/<username>[/<subdir>*] (the trailing / after <username>
# can not be omitted if a subdir is not given) and
# /store/group/<groupname>[/<subgroupname>*] (and
# /store/local/<dir>[/<subdir>*] if Data.publication =
# False). Defaults to /store/user/<username>/. CRAB creates the
# outLFNDirBase path on the storage site if needed, do not create it
# yourself otherwise the file stage-out may fail due to permissions
# inconsistency.

config.Data.outLFNDirBase = '/store/group/phys_muon/dmytro/'

# Whether to publish or not the EDM output files (i.e. output files
# produced by PoolOutputModule) in DBS. Notice that for publication to
# be possible, the corresponding output files have to be transferred
# to the permanent storage element. Defaults to True.
config.Data.publication = True

# When running an analysis over private input files or running MC
# generation, this parameter specifies the primary dataset name that
# should be used in the LFN of the output/log files and in the
# publication dataset name
config.Data.outputPrimaryDataset = 'QCD_Pt-30toInf_BmmGenFilter_TuneCP5_13TeV_pythia8'

# A custom string used in both, the LFN of the output files (even if
# Data.publication = False) and the publication dataset name (if
# Data.publication = True)
config.Data.outputDatasetTag ='RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15'

##
## Site
##

# Site where the output files should be permanently copied to. See the
# note about storageSite below. The user MUST have write permission in
# the storage site.
config.Site.storageSite = 'T2_CH_CERN'

#config.Site.storageSite = 'T3_US_MIT'

