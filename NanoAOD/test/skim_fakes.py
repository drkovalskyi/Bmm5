import FWCore.ParameterSet.Config as cms
process = cms.Process('SKIM')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Bmm5.NanoAOD.MuonFakeFilter_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/5C3622AA-CC64-F349-86B5-98A9308B485E.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B88FA0A2-D90E-B14E-9C3F-C6642B87F56B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/08518D08-B188-8C4F-8E9B-E26E58B17B1E.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/96B75080-F718-C145-8CB4-4950825AB09C.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/F1DAE84D-7842-AB48-B545-F59B847376D4.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/74077EE3-C078-4C49-8BC2-221072F6687B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/25766CE9-CD77-034B-BBD8-608062CB156D.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/32E1ADA1-FF2F-ED4A-B8D0-393CEEFFA399.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/7118966B-AF6C-BB49-9BAC-7EEB3B473053.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2D080D33-DBC2-BB45-A064-1AB1660C4A42.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/E9364FD4-6C61-A243-9F43-5B75052AC651.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/FF46F4FE-92EA-D840-8A0A-1BCD644D9066.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2FB687F1-7369-D547-B19D-87711660B32B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B8FD36E5-749C-3E49-B79C-C93548FE747B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/AE54591B-2E25-B14F-A116-5BCDAD5FECC2.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/C8BD1426-85D0-F842-957A-F2311DA93AAE.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/ACF7CCED-7AD7-AE46-A0B4-07F976ECB660.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B967392F-676A-734F-B89A-F9D52693E4E3.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/6B35D7AF-A112-AA4A-8F6D-6BA23837C2D4.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2CCD1BF1-9649-AC46-9394-59557894F89D.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/6C1476F5-CDAB-4E44-A4DC-8F453B41675E.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/BE3FE6CC-B57B-6541-8D14-794E212693A4.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B02436D1-992E-9D47-B9F6-2EFF277A85D6.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/E114E832-2311-9140-9B10-894BB715B954.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/37C523A2-F919-2D49-9B4A-24063D8D3157.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B1431570-EC92-C041-AD42-1FF0E0B8D1D6.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/33E40D89-63FD-DE41-9D85-FD372342C851.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/43EDAC42-19E4-074C-BF63-E7245B953BD0.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2DC54373-4CE6-CA40-AACF-B827CAAF96CB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/A28B0225-69AF-A34B-9AC4-43D7909619E4.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/0C753584-3BF2-D642-BE10-A6F577CCB4B0.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/7C28278A-FE74-8449-A572-EFD0D8651A41.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/44B8D07F-93F1-FB4A-BC4A-E376687308BB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B284E9D4-9C89-F94A-9FEF-E98069C0BB74.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B9B911A1-7D70-884F-A850-71C0C3958BD5.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/A17CA3CA-D162-D545-9208-5F358C408250.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/D9B62F4B-FF35-724C-B991-E799BADCD348.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/F9AB084E-8F0B-4942-B0E6-5EE4A9F697DA.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/519689EA-0A99-CD49-BF9F-F4917EEAA197.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B395D325-A3F1-464F-AE25-445DA7C8F731.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/4DFFC2E9-6278-3E4A-9D09-FE84655D2C6A.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2EEB961A-317D-604B-8F84-4FE0F523E25D.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/46FE6653-C94E-9444-BF5C-AA257FF2B432.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/AE53B33C-3037-0E41-AA3B-98C0419C3862.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/DFEF7DEF-6816-B74C-A305-F20131D00FAB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/EB720F14-CD63-0347-A1A1-25886EA3E2F5.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/5E39D516-B2AD-0B49-BC45-B0F4443135F4.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/B221648D-5095-E24A-9232-4ECD3BDB29C5.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/D0D4FF16-A2C6-D84D-92C6-D70C9E771A20.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2D091B30-22EA-CA40-ABD1-E4C93F8A304A.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/3FA9DB11-D2FD-4C4D-96B8-8B49054B0AB4.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/73B9E6A6-6D26-2041-A619-F637D230DB30.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/73E8ADBF-9697-C147-8DE9-492CEE0DC25D.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/CDCEE538-EE9B-5D40-BA6A-0293CCDA2CA8.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/5E5EF9C9-2000-1942-9AB9-65B207AB3BFC.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/A9BF8AB9-F6DE-8D4E-B7C2-9B6888A70421.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/E1096FB7-A980-E44A-9F8D-1B5D827EDB0F.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/C2D8E774-1506-7140-B22E-74E67AA2F103.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/E2A76731-22B2-FF4F-BF0F-D18CBB9612AE.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/4EE5A8EA-929F-974A-8E43-369244BC4CD3.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/CCC95E07-CF21-7E4C-8EE8-83B25A3BAE2E.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/9CF11AEB-BE30-3C43-B332-1E2DFB6ABD24.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/AE938C92-56C8-B24E-A4AC-208582050CD6.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/F7C10BA4-910C-4847-BE60-CEA2F3D40820.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/57AF39E3-CCE3-F843-A594-CC5910A5EDEC.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/D567BC0F-A9CF-1F4A-B876-D7AE9517C6B2.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/D1F627D2-EEAD-C142-AC31-4575EB6E0416.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/3E12C3CF-B636-E340-989D-813DE335163C.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/A8C60797-9A09-0A47-A134-0B925E50D5DD.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/2D85FA91-D6C6-3948-ABAB-A19AD0EDE886.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/61C0F105-B70D-2040-BC1D-2E59F976F22B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/9EEE02A9-C881-FB43-90FF-156224301F5A.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/01E77550-F5DA-FC45-8E97-EDF489BBF3B6.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/C79A872F-7BF4-8846-B8B1-09759A99D562.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/3C247C0D-7804-4D49-9067-08CAB165FFEB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/99D86C9B-6399-7941-8A6B-76F1A327EBDB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/3CEC5C26-14EC-204F-9FCB-C7FB4E79B3D5.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/7AE7BE29-234E-A848-9C18-E24124FB3169.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/3447C35E-4F21-2C4C-8AAE-AA7F919C85D6.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/D2E15AB5-86C5-724D-B7C4-B50D40422B8E.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/6182E34F-4AC1-B943-AE9A-C28BB63CCF82.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/719CC2DD-AFE5-CE47-A219-13C39794C5B6.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/05D93C65-717B-2642-A0B0-D66EC70B15A5.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/6B3F23CF-B7EC-D746-AE2B-116E9E546B63.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/93B2314F-843E-D440-8C29-93048DF00BBD.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/A594A3E3-961E-4446-AE73-34359A9F8D1F.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/10000/BA88F258-CE72-4646-8631-9460D99DA355.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/2F958CB5-268F-5D4B-99A8-99C6433B681B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/5DFC459B-1BE3-0D48-8827-54D9350054DC.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/860C9118-A34C-8244-8963-F0E2B6BFE4B5.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/F6635EC9-0685-B249-A129-F18ADF44F9A9.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/7B8C3EC1-7133-D041-A92C-6893E73419DB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/6A367C5A-B386-CB41-9DF6-83C43E4E212B.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/F407B405-EFA9-DF42-A47D-797510208882.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/D94DD567-A4B9-534E-A962-8774D7DFA5BA.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/F5B2827D-5C9F-A94A-A4D3-DC0B7D7856BD.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/F76CE91B-EE80-7041-80C1-0C4C4FB4DA10.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/B4AA5325-3CD7-324E-BF6E-53CD2A0531B7.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/5935654B-97A7-AE4D-8287-CA59C58FF23C.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/1C97EA9E-9483-324E-8184-1A9138741D97.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/BC209D08-9205-D74B-A589-C7E4FF9F5ABB.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/0DEB659F-B9E8-FD4E-9F6B-3230945CACAE.root',
        '/store/mc/RunIIAutumn18MiniAOD/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/280000/6D54F65F-15F4-B94A-9A8D-639BD2C0A884.root',)
)

process.outputPath = cms.EndPath()

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('file:BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_RunIIAutumn18MiniAOD_muon_fake_skim_2.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filterPath'))                              
)

process.filterPath = cms.Path( process.muonFakeFilter )
process.output_step = cms.EndPath( process.output )
process.schedule = cms.Schedule( process.filterPath, process.output_step )

# Spit out filter efficiency at the end.                                                                                                                                         
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
