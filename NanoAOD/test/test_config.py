"""Test process definitions for Bmm5/NanoAOD regression testing.

Reference layout on EOS: REFERENCE_BASE / <tag> / <sample> /
Local output:            /tmp/$USER/bmm_test / <sample> /

Test levels:
  smoke       - 5% of nominal (quick sanity check)
  nominal     - base event count (commit comparison)
  pre-release - 20x nominal, single sample (performance gate)
  release     - 20x nominal (tagged release)
"""

import subprocess
import os
from pathlib import Path

# --- Test levels ---

LEVELS = {
    "smoke":       0.05,
    "nominal":     1.0,
    "pre-release": 20.0,
    "release":     20.0,
}

def get_nevents(sample_nevents, level):
    """Compute event count for a given level."""
    return max(1, int(sample_nevents * LEVELS[level]))

# --- Paths ---

REFERENCE_BASE = "/eos/cms/store/user/dmytro/tmp/bmm_test_references"
LOCAL_OUTPUT = f"/tmp/{os.environ.get('USER', 'nobody')}/bmm_test"

# Sample used for smoke and commit-level tests (single sample)
DEFAULT_SAMPLE = "BuToJpsiK"

# Sample used for pre-release performance validation
PRE_RELEASE_SAMPLE = "BsToMuMu"

# --- Common customizations ---

BMM_CUSTOMISE = [
    "Bmm5/NanoAOD/nano_cff.nanoAOD_customizeDileptonPlusX",
    "Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake",
    "Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBmmMuonId",
]

BMM_CUSTOMISE_COMMANDS = [
    "process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))",
]

PERF_CUSTOMISE = [
    "Configuration/DataProcessing/Utils.addMonitoring",
    "Validation/Performance/TimeMemoryInfo.py",
]

PERF_CUSTOMISE_COMMANDS = [
    "process.Timing.summaryOnly = cms.untracked.bool(True)",
]

ALL_CUSTOMISE = BMM_CUSTOMISE + PERF_CUSTOMISE
ALL_CUSTOMISE_COMMANDS = BMM_CUSTOMISE_COMMANDS + PERF_CUSTOMISE_COMMANDS

# --- Git version detection ---
# All git commands run in the Bmm5 directory (separate repo from CMSSW)

BMM5_DIR = str(Path(__file__).resolve().parent.parent.parent)

def _git(*args):
    """Run a git command in the Bmm5 repo directory."""
    return subprocess.check_output(
        ["git"] + list(args),
        cwd=BMM5_DIR,
        stderr=subprocess.DEVNULL
    ).decode().strip()

def get_git_version():
    """Return git tag (if HEAD is tagged) or short commit hash."""
    try:
        return _git("describe", "--tags", "--exact-match", "HEAD")
    except subprocess.CalledProcessError:
        pass
    try:
        return _git("rev-parse", "--short", "HEAD")
    except subprocess.CalledProcessError:
        return "unknown"

def get_git_commit():
    """Return full commit hash."""
    try:
        return _git("rev-parse", "HEAD")
    except subprocess.CalledProcessError:
        return "unknown"

def get_latest_tag():
    """Return the latest tag reachable from HEAD."""
    try:
        return _git("describe", "--tags", "--abbrev=0")
    except subprocess.CalledProcessError:
        return None

# --- Sample definitions ---

processes = {

    "BsToMuMu": {
        "type": "mc",
        "era": "Run3,run3_nanoAOD_124",
        "conditions": "auto:phase1_2022_realistic_postEE",
        "input": [
            "/store/user/dmytro/tmp/store+mc+Run3Summer22EEMiniAODv3+BsToMuMu_BMuonFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+124X_mcRun3_2022_realistic_postEE_v1-v2+2820000+0096d5dd-88d3-46a0-a8cc-255a3090c71e.root",
        ],
        "customise": ALL_CUSTOMISE + [
            "Bmm5/NanoAOD/nano_cff.run3_nanoAOD_124",
        ],
        "customise_commands": ALL_CUSTOMISE_COMMANDS,
        "nevents": 500,
        "exclude_branches": "",
        "selections":{
            "reco (mc matched)":"mm_kin_vtx_prob>0.01&&mm_mu1_pdgId*mm_mu2_pdgId==-169&&mm_gen_pdgId!=0",
            "reco (mc not matched)":"mm_kin_vtx_prob>0.01&&mm_mu1_pdgId*mm_mu2_pdgId==-169&&mm_gen_pdgId==0",
        },        
        "plots": {
            "B-mass (reco, mc matched)": [
                "mm_kin_mass>>h(100,4.9,5.9)",
                "mm_kin_vtx_prob>0.1&&mm_mu1_pdgId*mm_mu2_pdgId==-169&&abs(mm_kin_mass-5.4)<0.5&&mm_gen_pdgId!=0"],
        }
    },

    "BsToJPsiPhi": {
        "type": "mc",
        "era": "Run3_2024",
        "conditions": "auto:phase1_2024_realistic",
        "input": [
            "file:/eos/cms/store/user/dmytro/tmp/store+mc+RunIII2024Summer24MiniAOD+BsToJPsiPhi-JpsiToMuMu-PhiToKK_Fil-MuPt2_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+140X_mcRun3_2024_realistic_v26_ext1-v3+140000+7055677a-69c3-48af-9a14-c6870cef04d6.root",
        ],
        "customise": ALL_CUSTOMISE,
        "customise_commands": ALL_CUSTOMISE_COMMANDS,
        "nevents": 500,
        "exclude_branches": "",
        "selections":{
            "reco (mc matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge",
            "reco (mc not matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId==0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge",
            "nominal (mc matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge" \
                                   "&&abs(bkkmm_jpsikk_alpha)<0.01&&bkkmm_jpsikk_sl3d>5&&abs(bkkmm_kk_mass-1.02)<0.01",
            "tight (mc matched)":"bkkmm_jpsikk_vtx_prob>0.1&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge" \
                                 "&&abs(bkkmm_jpsikk_alpha)<0.005&&bkkmm_jpsikk_sl3d>5&&abs(bkkmm_kk_mass-1.02)<0.01",
            "tight (mc not matched)":"bkkmm_jpsikk_vtx_prob>0.1&&bkkmm_gen_pdgId==0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge" \
                                     "&&abs(bkkmm_jpsikk_alpha)<0.005&&bkkmm_jpsikk_sl3d>5&&abs(bkkmm_kk_mass-1.02)<0.01"
        },        
        "plots": {
            "B-mass (reco, mc matched)": [
                "bkkmm_jpsikk_mass>>h(100,4.9,5.9)",
                "bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge"],
            "B-mass (nominal, mc matched)": [
                "bkkmm_jpsikk_mass>>h(100,4.9,5.9)",
                "bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge&&" \
                "abs(bkkmm_jpsikk_alpha)<0.01&&bkkmm_jpsikk_sl3d>5&&abs(bkkmm_kk_mass-1.02)<0.01"],
            "B-mass (tight, mc matched)": [
                "bkkmm_jpsikk_mass>>h(100,4.9,5.9)",
                "bkkmm_jpsikk_vtx_prob>0.1&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge&&" \
                "abs(bkkmm_jpsikk_alpha)<0.005&&bkkmm_jpsikk_sl3d>5&&abs(bkkmm_kk_mass-1.02)<0.01"],
        }
    },

    "BuToJpsiK": {
        "type": "mc",
        "era": "Run3_2024",
        "conditions": "auto:phase1_2024_realistic",
        "input": [
            # "file:/tmp/dmytro/dmytro/store+mc+RunIII2024Summer24MiniAOD+BuToJpsiK-JpsiToMuMu_Fil-MuPt2_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+140X_mcRun3_2024_realistic_v26-v2+110000+04ce715e-06f7-4858-b27c-5b7feff16c09.root",
            "/store/user/dmytro/tmp/store+mc+RunIII2024Summer24MiniAOD+BuToJpsiK-JpsiToMuMu_Fil-MuPt2_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+140X_mcRun3_2024_realistic_v26-v2+110000+04ce715e-06f7-4858-b27c-5b7feff16c09.root",
        ],
        "customise": ALL_CUSTOMISE,
        "customise_commands": ALL_CUSTOMISE_COMMANDS,
        "nevents": 500,
        "exclude_branches": "",
        "selections":{
            "reco (mc matched)":"bkmm_gen_pdgId!=0",
            "loose (mc matched)":"bkmm_jpsimc_sl3d>3 && bkmm_jpsimc_vtx_prob>0.01 && bkmm_gen_pdgId!=0",
            "loose (mc not matched)":"bkmm_jpsimc_sl3d>3 && bkmm_jpsimc_vtx_prob>0.01 && bkmm_gen_pdgId==0",
            "tight (mc matched)":"bkmm_jpsimc_sl3d>5 && bkmm_jpsimc_vtx_prob>0.1 && bkmm_gen_pdgId!=0 && abs(bkmm_jpsimc_alpha)<0.005",
            "tight (mc not matched)":"bkmm_jpsimc_sl3d>5 && bkmm_jpsimc_vtx_prob>0.1&&bkmm_gen_pdgId==0&&abs(bkmm_jpsimc_alpha)<0.005",
        },
        "plots": {
            "B-mass (loose, mc matched)": [
                "bkmm_jpsimc_mass>>h(100,4.9,5.9)",
                "bkmm_jpsimc_sl3d>3 && bkmm_jpsimc_vtx_prob>0.01 && bkmm_gen_pdgId!=0"],
            "B vertex probability (loose, mc matched)": [
                "bkmm_jpsimc_vtx_prob>>h(100,0,1)",
                "bkmm_jpsimc_sl3d>3 && bkmm_jpsimc_vtx_prob>0.01 && bkmm_gen_pdgId!=0"],
        }
    },

    "BdToJpsiKstar": {
        "type": "mc",
        "era": "Run3_2024",
        "conditions": "auto:phase1_2024_realistic",
        "input": [
            "/store/user/dmytro/tmp/store+mc+RunIII2024Summer24MiniAOD+BdToJpsiKstar_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+140X_mcRun3_2024_realistic_v26-v2+2550000+65cb1337-c183-496e-8233-0bd7fa4d19a5.root",
        ],
        "customise": ALL_CUSTOMISE,
        "customise_commands": ALL_CUSTOMISE_COMMANDS,
        "nevents": 500,
        "exclude_branches": "",
        "selections":{
            "reco (mc matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge",
            "reco (mc not matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId==0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge",
            "loose (mc matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge" \
                                 "&&(abs(bkkmm_jpsikpi_hh_mass-0.89167)<0.05||abs(bkkmm_jpsipik_hh_mass-0.89167)<0.05)" \
                                 "&&bkkmm_jpsikk_sl3d>3&&abs(bkkmm_jpsikk_alpha)<0.1",
            "loose (mc not matched)":"bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId==0&&bkkmm_kaon1_charge!=bkkmm_kaon2_charge" \
                                     "&&(abs(bkkmm_jpsikpi_hh_mass-0.89167)<0.05||abs(bkkmm_jpsipik_hh_mass-0.89167)<0.05)" \
                                     "&&bkkmm_jpsikk_sl3d>3&&abs(bkkmm_jpsikk_alpha)<0.1",
        },        
        "plots": {
            "B-mass (reco, mc matched, true mass hypothesis)": [
                "abs(bkkmm_gen_kaon1_pdgId)==321?bkkmm_jpsikpi_mass:bkkmm_jpsipik_mass>>h(100,4.9,5.9)",
                "bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge != bkkmm_kaon2_charge"],
            "B-mass (reco, mc matched, best mass hypothesis)": [
                "abs(bkkmm_jpsikpi_hh_mass-0.89167)<abs(bkkmm_jpsipik_hh_mass-0.89167)?bkkmm_jpsikpi_mass:bkkmm_jpsipik_mass>>h(100,4.9,5.9)",
                "bkkmm_jpsikk_vtx_prob>0.01&&bkkmm_gen_pdgId!=0&&bkkmm_kaon1_charge != bkkmm_kaon2_charge&&" \
                "(abs(bkkmm_jpsikpi_hh_mass-0.89167)<0.05||abs(bkkmm_jpsipik_hh_mass-0.89167)<0.05)"],
        }
    },
    "BdToJpsiKShort": {
        "type": "mc",
        "era": "Run3_2024",
        "conditions": "auto:phase1_2024_realistic",
        "input": [
            "file:/eos/cms/store/user/dmytro/tmp/store+mc+RunIII2024Summer24MiniAOD+BdToJpsiKShort-JpsiTo2Mu-KShortTo2Pi_Fil-Jpsi-KShort-Mu_TuneCP5_13p6TeV_pythia8-evtgen+MINIAODSIM+140X_mcRun3_2024_realistic_v26-v2+120000+018e5aeb-b5dc-404b-a279-145d0e65329b.root",
        ],
        "customise": ALL_CUSTOMISE,
        "customise_commands": ALL_CUSTOMISE_COMMANDS,
        "nevents": 500,
        "exclude_branches": "",
        "selections":{
            "reco (mc matched)":"bjpsiks_gen_pdgId!=0",
            "loose (mc matched)":"bjpsiks_kin_sl3d>3 && bjpsiks_kin_vtx_prob>0.001 && bjpsiks_gen_pdgId!=0",
        },
    },

    "Data": {
        "type": "data",
        "era": "Run3",
        "conditions": "140X_dataRun3_Prompt_v4",
        "input": [
            "/store/user/dmytro/tmp/store+data+Run2024G+ParkingDoubleMuonLowMass0+MINIAOD+PromptReco-v1+000+385+620+00000+b8d018fd-5c14-4e02-9d0b-d455cff1e2f9.root",
        ],
        "customise": ALL_CUSTOMISE,
        "customise_commands": ALL_CUSTOMISE_COMMANDS,
        "nevents": 500,
        "exclude_branches": "",
        "selections":{
            "jpsi reco+hlt":"mm_kin_vtx_prob>0.1&&mm_mu1_pdgId*mm_mu2_pdgId==-169&&abs(mm_kin_mass-3.1)<0.1&&HLT_DoubleMu4_3_LowMass",
        },        
        "plots": {
            "jpsi-mass": [
                "mm_kin_mass>>h(100,2.9,3.3)",
                "mm_kin_vtx_prob>0.1&&mm_mu1_pdgId*mm_mu2_pdgId==-169&&abs(mm_kin_mass-3.1)<0.2"],
        }
    },
}
