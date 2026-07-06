from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def pixel_stuckTBM_enabled(globaltag):
    # The SiPixelQuality stuckTBM (auto-masking) PCL tag is in no GlobalTag and
    # its IOVs only start mid-2018 (run 315690). Run-3 is exactly the run range
    # (>= 355k) for which it is guaranteed to be covered, so key the decision on
    # the data GlobalTag: 'dataRun3' -> enable, everything else (Run-2 data, MC)
    # -> disable. This is derived from the GT the job already sets; no manual
    # switch. Run-2 keeps dead/total counts, the geometry-driven denominators,
    # and the per-muon HitPattern flags -- only the event-level masked columns
    # fall back to the -1 sentinel (auto-masking was not recorded before 2018).
    return 'dataRun3' in globaltag

def disable_legacy_tau_reprocessing(process):
    # Run-2 NanoAOD (run2_nanoAOD_106Xv2) re-embeds the legacy tau anti-electron
    # MVA6 discriminant, whose GBRForest payloads no longer exist in the conditions
    # DB for >= 14_0_X releases -- so the job aborts with
    # "No data of type GBRForest ... RecoTauTag_antiElectronMVA_...". Bmm uses no
    # taus, but PATObjectCrossLinker ('linkedObjects', consumed by BmmMuonId) needs
    # valid pat::Tau collections. Fix: drop the anti-electron-MVA6 entries from
    # every tau-ID embedder's tauIDSources, so its producers are never consumed and
    # therefore never scheduled / never touch the missing conditions; and drop the
    # tau tables (their Vars would read the now-absent ID). The deepTau / isolation
    # IDs (conditions present) stay, so finalTaus / finalBoostedTaus and the tau
    # collections are otherwise unchanged. Keyed on module contents, not on
    # release-specific collection names, so it works in 14_0_X and 15_0_X alike;
    # no-op for Run 3 (no anti-electron-rejection producer is scheduled there).
    # Must be invoked AFTER nanoAOD_customizeCommon (which builds the tau
    # reprocessing) -- i.e. from --customise_commands, not a --customise function.
    if not any('ElectronRejection' in n for n in process.producers_()):
        return process

    # Walk a selector's src back to the raw MiniAOD tau collection: follow each
    # reprocessing producer's own 'src' until we reach a label with no producer in
    # this process, or an InputTag that explicitly names the input process (e.g.
    # "@skipCurrentProcess", used in 15_0_X where the embedder shares the MiniAOD
    # collection's name). Release-agnostic -- no hardcoded reprocessed-collection
    # names.
    def raw_src(tag, _seen=None):
        _seen = _seen or set()
        proc = tag.getProcessName()
        if proc and proc != '@currentProcess':
            return tag
        label = tag.getModuleLabel()
        if label in _seen:
            return tag
        _seen.add(label)
        prod = getattr(process, label, None)
        if prod is not None and hasattr(prod, 'src') and isinstance(prod.src, cms.InputTag):
            return raw_src(prod.src, _seen)
        return tag

    for sel, ptcut in (('finalTaus', 'pt > 18'), ('finalBoostedTaus', 'pt > 40')):
        m = getattr(process, sel, None)
        if m is not None and hasattr(m, 'src'):
            m.src = raw_src(m.src)
            m.cut = cms.string(ptcut)
    if hasattr(process, 'nanoTableTaskCommon'):
        for _t in ('tauTablesTask', 'boostedTauTablesTask'):
            if hasattr(process, _t):
                try:
                    process.nanoTableTaskCommon.remove(getattr(process, _t))
                except Exception:
                    pass
    return process

def run3_nanoAOD_124(process):
    process.finalTaus.cut = cms.string("pt > 18 && ((tauID('decayModeFindingNewDMs') > 0.5 && (tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') || (tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum'))<2.5) || tauID('byVVVLooseDeepTau2017v2p1VSjet') || (tauID('byDeepTau2018v2p5VSjetraw') > {}))) || (?isTauIDAvailable('byUTagCHSVSjetraw')?tauID('byUTagCHSVSjetraw'):-1) > {} || (?isTauIDAvailable('byUTagPUPPIVSjetraw')?tauID('byUTagPUPPIVSjetraw'):-1) > {})".format(WORKING_POINTS_v2p5["jet"]["VVVLoose"], 0.05, 0.05))
    
    return process

def nanoAOD_keepLowPtMuons(process):
    process.muonTable.doc = cms.string("slimmedMuons after basic selection (pt > 2 || (pt > 2 && (passed(\'CutBasedIdLoose\') || passed(\'SoftCutBasedId\') || passed(\'SoftMvaId\') || passed(\'CutBasedIdGlobalHighPt\') || passed(\'CutBasedIdTrkHighPt\'))))")

    process.finalMuons.cut = cms.string("pt > 2 || (pt > 2 && (isStandAloneMuon() || passed(\'CutBasedIdLoose\') || passed(\'SoftCutBasedId\') || passed(\'SoftMvaId\') || passed(\'CutBasedIdGlobalHighPt\') || passed(\'CutBasedIdTrkHighPt\')))")

    return process

def nanoAOD_customizeDileptonPlusX(process):

    nanoAOD_keepLowPtMuons(process)
    
    process.load('Bmm5.NanoAOD.DileptonPlusX_cff')
    process.load('Bmm5.NanoAOD.UpdateSlimmedMuons_cff')
    process.load('Bmm5.NanoAOD.PixelQuality_cff')
    process.load('PhysicsTools.NanoAOD.muons_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.DileptonPlusXSequence + process.DileptonPlusXTables + process.PixelQualityTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.DileptonPlusXMcSequence + process.DileptonPlusXMcTables + process.PixelQualityMcTables)

    # SiPixelQuality stuckTBM (auto-masking) tag: not in any GlobalTag, and its
    # IOVs only start mid-2018, so it is loaded + read only for Run-3 data. The
    # decision is taken automatically from the job's GlobalTag (see
    # pixel_stuckTBM_enabled) -- no manual switch, nothing beyond the --conditions
    # and --era cmsDriver already needs. Data only (MC uses pixelQualityMcTable,
    # which does not read conditions). An appended-but-uncovered tag would
    # invalidate the whole SiPixelQualityFromDbRcd record for the job
    # (CondDBESSource intersects IOV validity across a record's labels), which is
    # why Run-2 must not append it; there dead/total counts and the per-muon
    # flags still fill, masked columns get the -1 sentinel.
    if hasattr(process, 'NANOAODoutput'):
        enable_stuckTBM = pixel_stuckTBM_enabled(process.GlobalTag.globaltag.value())
        process.pixelQualityTable.readStuckTBM = cms.bool(enable_stuckTBM)
        if enable_stuckTBM:
            process.GlobalTag.toGet.append(cms.PSet(
                record = cms.string("SiPixelQualityFromDbRcd"),
                tag    = cms.string("SiPixelQuality_byPCL_stuckTBM_v1"),
                label  = cms.untracked.string("stuckTBM"),
            ))
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)

    # NB: the legacy-tau bypass (disable_legacy_tau_reprocessing) is NOT called
    # here -- the tau reprocessing is built by nanoAOD_customizeCommon, which
    # cmsDriver always runs AFTER the user --customise functions, so anything set
    # here would be clobbered. Invoke it via --customise_commands instead (runs
    # last); see the helper's docstring. It is a no-op on Run 3.

    # enforce process name
    # process.load('PhysicsTools.NanoAOD.globals_cff')
    process.genFilterTable.src = cms.InputTag("genFilterEfficiencyProducer")

    # keep all genparticles
    process.finalGenParticles.select = cms.vstring("keep *")

    # save gen particle vertex
    process.genParticleTable.variables = cms.PSet(
        process.genParticleTable.variables,
        vx = Var("vx", "float", doc="x coordinate of vertex position"),
        vy = Var("vy", "float", doc="y coordinate of vertex position"),
        vz = Var("vz", "float", doc="z coordinate of vertex position")
    )

    process.load('FWCore.MessageService.MessageLogger_cfi')
    # Kill all messages from that category on 'cerr'
    process.MessageLogger.cerr.TwoTrackMinimumDistance = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    )
    
    return process

def nanoAOD_customizeV0ForMuonFake(process):
    process.load('Bmm5.NanoAOD.BmmV0ForMuonFake_cff')
    process.load('Bmm5.NanoAOD.UpdateSlimmedMuons_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.V0ForMuonFakeSequence + process.V0ForMuonFakeTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
    return process

def nanoAOD_customizeBmmMuonId(process):
    process.load('Bmm5.NanoAOD.BmmMuonId_cff')
    # Data 
    process.nanoSequence   = cms.Sequence(process.nanoSequence + process.BmmMuonIdSequence + process.BmmMuonIdTables)
    # MC
    process.nanoSequenceMC = cms.Sequence(process.nanoSequenceMC + process.BmmMuonIdMcSequence + process.BmmMuonIdMcTables)

    process.load('PhysicsTools.NanoAOD.muons_cff')
    # process.muonTable.variables.mvaMuID = Var("mvaIDValue()",float,doc="MVA-based ID score ",precision=6)
    return process

