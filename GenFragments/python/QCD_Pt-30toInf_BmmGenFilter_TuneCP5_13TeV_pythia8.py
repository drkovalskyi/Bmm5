import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
        maxEventsToPrint = cms.untracked.int32(1),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        filterEfficiency = cms.untracked.double(0.03844),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        comEnergy = cms.double(13000.0),

        crossSection = cms.untracked.double(2.75842e+06),

        PythiaParameters = cms.PSet(
            pythia8CommonSettingsBlock,
            pythia8CP5SettingsBlock,
            processParameters = cms.vstring(
                        'ParticleDecays:limitTau0 = off',
                        'ParticleDecays:limitCylinder = on',
                        'ParticleDecays:xyMax = 2000',
                        'ParticleDecays:zMax = 4000',
                        'HardQCD:all = on',
                        'PhaseSpace:pTHatMin = 30',
                        'PhaseSpace:pTHatMax = 13000',
                        '130:mayDecay = on',
                        '211:mayDecay = on',
                        '321:mayDecay = on'
            ),
            parameterSets = cms.vstring('pythia8CommonSettings',
                                        'pythia8CP5Settings',
                                        'processParameters',
                                        )
        )
)

bmm_filter = cms.EDFilter("GenBmmFilter",
                          max_doca = cms.double(0.25),
                          min_mu_pt = cms.double(3.5),
                          max_mu_eta = cms.double(1.5),
                          min_mm_mass = cms.double(4.9),
                          max_mm_mass = cms.double(5.9),

)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('\$Revision$'),
    name = cms.untracked.string('\$Source$'),
    annotation = cms.untracked.string('QCD dijet production, pThat 30toInf GeV, with BmmGenFilter, 13 TeV, TuneCP5')
)

ProductionFilterSequence = cms.Sequence(generator*bmm_filter)
