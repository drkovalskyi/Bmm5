
import ROOT

path = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/523"
samples = [
    {
        'final_state':'mm',
        'name':'\dzmm',
        'files':[
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/089fcee0-0260-470b-8f08-a458129f2c4a.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/164a5846-abdf-44ae-ba68-00ac4669fdff.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1fbe986c-bc25-4a90-a0ea-03d2d4f76001.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/27fb3f48-a6ea-4674-9fdd-af87426e0c0d.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2f2337ac-420c-477c-a781-1065f0c218a8.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2ff197b6-51fc-404f-8109-1db55acddc01.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/3236c958-eb94-4338-95e4-8978dbb3d590.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4264fdbb-4822-43bb-966b-2a3da02bc49b.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5c7f65e6-9a42-4e87-9463-32b437dfd863.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/7ff86924-1e30-4633-ba52-0915badd3f66.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/95982060-e6de-4d62-b0e2-15dd9adaac06.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/a121d5c0-d816-42e6-a415-31706a143da1.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ae31daa2-7105-4bdc-b456-c1c10d9a81ba.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/af37d534-abfe-4ad8-935a-22c6971b6ef8.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/b9c0f04e-75d8-489a-9670-b1e720ce532e.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/bccd752e-2f56-494b-83f0-86fd0d8059c4.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c7ab600c-b132-42d1-a07c-4dab362332bc.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d0400b0b-f754-42ee-ba1e-cd9838c2a4aa.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d3b0ec5b-97c6-4b66-8370-c61de101327a.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/dbf39517-9639-4f68-9fb9-eecf8c21e431.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/de11d9f6-f03f-4b76-806c-65120b70cf36.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e86f4986-bd7d-4739-acd9-c4da17b471ec.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ea54cd31-9792-4d65-8c21-1236921055ff.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/eb9dfbdc-addd-4c63-903d-8edc27c33428.root",
            path + "/DstarToD0Pi_D0To2Mu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f5217816-f36b-441a-8d28-6112aeab01a7.root",
        ],
        'chain':None
    },
    {
        'final_state':'fakes',
        'name':r'\dzpipimm',
        'files':[
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/0022a93d-6eb6-44f8-8ec1-2e39db36ccd0.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/02902f40-c433-46c8-b3be-ab6073863d83.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/02ab694d-62d0-4c1f-b87e-e6665d1d5886.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/07978518-9861-4ae3-8d81-cc022d96ca83.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/0a11c749-72d1-4736-a13a-8578ab5c185f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/0a659d88-a023-4a05-a388-038a9b2f6672.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/0b2914d7-41c4-44a4-8cf9-cc5283deeb4c.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/0d6ea42b-95e5-4cb3-94e7-7c1695c98147.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/11a99e8a-3f55-43d6-8ca7-10d5eb596310.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1325c0dd-56b3-442f-adc7-24c427f64aa6.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/13c5a451-56d8-4c78-9b06-b81b241423e0.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/13dc1195-e020-478c-bfee-447fd9831f94.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/14a3a0c3-a892-4e2e-8670-ed5154f38c15.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/151fe901-9024-40a5-a292-7142e8c6bc5f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/17ebcd48-8e07-4870-9c79-99b78d63839f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/185b40ed-07c5-4eac-8ea4-e49145537703.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/191f4612-7e23-4222-b1ad-9915bab91748.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/19cfd8d2-6424-4ce7-af3e-b3c796c56bf2.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1a33e6d5-353d-4b4f-92eb-c5761f0cb61c.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1c9ae6f7-8cce-4a10-a9d6-467e43f25a57.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1d414cae-5fa0-49cf-af10-a0b2f8b4a40b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/1ff299a8-adbe-4925-a9a2-b76beccbb9ab.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/23ba1069-8ff5-44c0-a7a8-3d1ab294fc19.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/25ba7d5a-5d53-4ffb-b7c1-ebd38920bbbe.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/26218d47-559e-4266-ad14-a524471e7003.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/28c5f86d-0fbc-4522-8fc7-96dede3d40e0.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2aea1c5c-ceab-4f69-8495-abca3a6ec6e5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2bddcea6-c502-4152-923b-e7cc7c436de1.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2e12961d-a786-4d7e-bc3a-43b7470e927a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2ed92edf-af85-4cb5-8aa1-50b58e045910.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2f249f10-b603-41c7-b6b0-2f57ebb40caa.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/2fd21647-3b9c-4e60-b311-af2003528244.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/360799d7-b64d-4259-84f1-06189472f16f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/3c3d2a66-86b6-4f93-8600-84e0f63510e8.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/3ffe70b0-23af-489a-9dff-b2d7f01ed7a4.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4101f721-984d-40c8-99f6-cf0e79300acb.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/426dba8c-b673-4448-99bd-8b277c35b55b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4546c916-3eb3-46dc-b15e-ac750eeddcba.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4d37a46a-86b0-48b7-9fcb-969290af8160.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4e883874-71b0-4903-a3ad-c251f5b6ed48.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4e9ec1f8-35df-498b-b48e-715fc08c5727.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4f564d6b-6529-4416-8296-b8f32c48c2c5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/4f8f4798-647a-44e6-a70c-45a8a560f62a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/50a73b8b-a6de-4aa2-8eb9-ed3f6d8f5376.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5145dd37-9ea1-43c6-9da3-42e6fc912193.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/54c8896b-be43-4144-969a-37ab084d7f85.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/575b1ccb-3883-48d8-8116-5a94b9e1ecab.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/586015fa-3ab1-46fd-ab17-8275a594da81.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5a283a40-f08d-454d-88a3-480571b12cb2.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5aff288f-a126-447d-ad17-14369083358b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5bd1d1ff-3707-44c5-9d47-0dd18b775b70.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5c8ddcd0-257b-4ce7-88a2-76dafd8a7147.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5cc72128-43e5-4605-9d58-c1f9bbd8782e.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5da35fdc-8cc8-41a9-a0e7-75b3c14e27c4.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5e7069e9-073b-4a30-a1bd-28d89f4b2bd7.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/5ef84280-d405-4d3e-982d-e2916c0ef59a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/643f6569-a490-4a7c-9fbb-789632d967d5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/657713fd-cf3e-40e9-9ac0-05b37b62b827.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/66179e35-7674-4113-bc8b-eff6d91df01d.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/682fd7d6-1849-4053-a0b6-8bb508323526.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/68ba506d-7223-40e5-b854-7b620c7bbdb7.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/6915c9ba-ab8a-4a1a-b49f-762246b1498d.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/6b7caf70-e4d7-43ad-b23b-8e553c6b6c1b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/6cf8493b-745b-4de4-9180-cc9a99ba4244.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/6d16b9f5-66ba-49d9-a0a7-e88c8a3358e9.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/72fb4b79-4c27-49ac-8f1c-f3edb377971a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/72fd43ef-90fa-4ad2-89e3-201610eae40c.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/742d5256-6b20-496b-b155-db41aeeb852d.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/7556fc41-7c1d-464b-95a9-338c59b544f8.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/79112b6f-1fc8-4b0e-ad11-cfe27e6603ef.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/7978be0b-0665-42cd-b1e9-5c50c4072b3b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/7ba7fc54-644f-406f-a77b-2e20a7b80012.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/7e198e91-1a99-48c5-8225-dd00af2126bc.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/814110ce-51f4-43a1-9349-fd15b208b340.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/846c0123-6b8b-4664-9d79-4381f7999993.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/84dbe7da-64c2-4de5-87fb-1dcbadc14f0a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/875489ad-5c51-4b85-989f-3e1d1e9fed0d.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/890ac6be-5554-44e8-8931-8373e0557c28.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/89b6886a-0717-462d-b50d-f5c132583d3c.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8a91ab00-23df-4479-bb19-5972a39666f8.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8bd8e27c-0c72-4c4a-a18c-549d851412f1.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8d4a6a02-c65a-4767-b921-a897daf713f4.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8d91c1b2-6b52-4147-aa9b-78581f577a2a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8db14ded-d293-4fad-95e8-09e31ad8295f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8e41e53d-a2eb-4833-af00-357107152f0a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/8f6893c3-b8f0-45f1-a5d1-83427a701aef.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/91fa0802-7d12-4156-8212-1cfded4700d7.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/92c8e716-2a47-4d9b-9eab-2da6ef4b81e8.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/93038767-0aa3-4a5a-8039-4e5a4eb245ba.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/9336607e-127a-4f44-91b6-fd8836bb7504.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/969fec18-a102-419e-b375-964b7e8f3797.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/9a1fd67a-3595-4d58-a6b1-7cac0b420563.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/9b364b11-0060-4874-9221-75cadae450e2.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/9bb8a62b-092c-4514-aacd-6e26742ecfe6.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/9d1c9e14-592c-4ef2-86cc-7d9dc5ecd44f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/9eb2e730-d215-4ce9-a878-447521201575.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/a04d5f3f-51aa-4da0-b891-a8516054d98d.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/a515d95b-8c60-42c8-b210-d8665c99762d.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/a5b10435-689c-481e-863d-2ab7538e68af.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ab6c29b1-8fc1-424e-9a8f-fd4097b9cd2e.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ac1c831f-4810-41fb-b383-ed48a6984516.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/b328bf12-9320-4ded-9ebd-ad1ca504f0bd.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/b4521e3b-8a63-449c-93d8-f4c461004bf5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/b5cbdda7-842e-4409-ab20-cad3f8d12c94.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/b9c18a57-865c-4ee4-9ab9-79dcbdc936f5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ba2316d1-d905-4a07-901e-52653f8bf9e5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ba35ec23-7d36-4869-9753-fcafc6081b25.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/bb26000f-0d1e-4ea3-b947-755246aac747.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c0a77d86-0950-407d-992b-fc22ecca7eac.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c5008768-05b4-4348-b7e2-ad03f505a13f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c510839b-3560-4967-8a24-5418e34e5eff.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c6a4af02-8857-46a3-a3b2-ecbe971db6ff.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c7688dd4-483c-4212-8a46-fde1fda6a4f3.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/c8e5f582-2c91-4205-acd1-78ebd7c23027.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/caad7dd9-12d4-4b7e-b04c-dde576b9f33b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/cc97e058-31bc-4f67-adf3-59ac12e03639.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/cd4bf18c-a4d3-4fb2-9478-e5549967d390.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/cdf190f7-3ec2-4dd4-8c34-e9a8d72c8fc4.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d392edcc-fbaf-4f3d-930d-71752a1e5d65.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d5c88cb3-3865-4b25-8a29-55802664f78f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d6f9af22-db31-4e2a-b7bd-12cb84577b65.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d7ac8f45-3348-4848-9087-c6884c4ffff1.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/d9a0d856-fa0a-4c1e-852c-9e6888ed7b5e.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/dd5b7334-9198-41c8-b8bf-b5d429574de6.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/de6ddb69-34a9-47f1-8fa0-7550981732d2.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e070fc51-eb70-434a-9c80-ee91bfa63770.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e16074e2-b44c-4339-9b7e-8385b7375f9a.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e2ffc78a-ff12-4384-ba7f-b6db8360acef.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e37e2638-77b5-42c5-8524-7fbb395b9685.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e5c8949c-077d-43d9-b118-f972ebab789e.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e5c9c74d-e928-4d1e-b4f9-2d4e8cc9d94c.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e5e25907-81b5-4def-a84f-99c100efd929.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/e7e7b309-b250-4ba6-9ce9-7c0093f0d290.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ea011839-cacd-41c5-9352-61f278584041.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/eb51f912-69a6-4508-9561-2c5f99695ff4.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ec1b7662-19dd-4f28-8744-33d6a9922588.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ed26a478-1db0-461a-9075-0837b34425b7.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ed5ec56f-58ac-48b9-81e7-ef9acfb69f09.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/edde5912-a726-498d-a15b-b326b8a85fcb.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f0ab34aa-ff81-447a-bb37-2f7ef0551d98.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f10243bd-a0d9-45f2-8402-3811c58893fb.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f285854d-c113-440a-9ad2-b3fc6fccf962.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f2e68f68-7777-4fbd-a90f-36d45563ee9f.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f4dd621e-b649-4617-8105-19d75b475c4b.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f530bf95-51e6-44cb-a8c1-42d2de9a03a1.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f5a9d4a3-de52-4c38-bcac-36f1ec1711f8.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f6e3722b-a5de-4fc0-8b5c-3244a2800789.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f806064e-9ea8-4764-92bc-930eddfc2dcd.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/f8e5c2f5-d60d-458d-93cb-3fc70134a776.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/faf3ec01-dbe6-46e9-8420-37c1022c2580.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/fb6a0785-a0c6-4fa1-8785-b5e8eabe47b1.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/fdaf0bed-07c9-4af4-9060-1ade5b0d5784.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/fe6204fc-a918-4d53-9017-951e540af604.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ff35b958-db53-4826-977d-d6e23e8f6fb5.root",
            path + "/DstarToD0Pi_D0To2Pi_PiToMu_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen+Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2+MINIAODSIM/ffa88ee2-bb40-4c51-aabb-0c708dcb701c.root",
        ],
        'chain':None
    },
    # {
    #     'final_state':'mmk',
    #     'name':'\bjpsik',
    #     'files':[
    #         "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
    #     ],
    #     'chain':None
    # },
    # {
    #     'final_state':'mmkk',
    #     'name':'blah',
    #     'files':[
    #         "/eos/cms/store/group/phys_muon/dmytro/tmp/RunIISummer20UL18MiniAOD_BsToJPsiPhi_JPsiToMuMu_PhiToKK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.root"
    #     ],
    #     'chain':None
    # }
]

cuts = [
    {
        'cut':{
            'mm':'dstar_gen_pdgId!=0 && dstar_mm_index>=0',
            'fakes':'dstar_gen_pdgId!=0 && dstar_mm_index>=0 && abs(mm_gen_cpdgId[dstar_mm_index])==421',
            # 'mmk':'mm_mu1_index[bkmm_mm_index]>=0 && mm_mu2_index[bkmm_mm_index]>=0 && bkmm_gen_pdgId!=0 && bkmm_jpsimc_mass>0'
            # + ' && Muon_pt[mm_mu1_index[bkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkmm_mm_index]]>4'
            # + ' && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkmm_mm_index]])<1.4',
            # 'mmkk':'mm_mu1_index[bkkmm_mm_index]>=0 && mm_mu2_index[bkkmm_mm_index]>=0 && bkkmm_gen_pdgId!=0 && bkkmm_jpsikk_mass>0'
            # + ' && Muon_pt[mm_mu1_index[bkkmm_mm_index]]>4 && Muon_pt[mm_mu2_index[bkkmm_mm_index]]>4'
            # + ' && abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4 && abs(Muon_eta[mm_mu2_index[bkkmm_mm_index]])<1.4',
        },
        'name':'MC matched baseline',
    },
    {
        'cut':{
            'mm':'mm_mu1_pt[dstar_mm_index] > 4 && mm_mu2_pt[dstar_mm_index] > 3',
            'fakes':'mm_mu1_pt[dstar_mm_index] > 4 && mm_mu2_pt[dstar_mm_index] > 3',
        },
        'name':r'$\ptmuone>4$, $\ptmutwo>3$',
    },
    {
        'cut':{
            'mm':'HLT_DoubleMu4_3_LowMass',
            'fakes':'HLT_DoubleMu4_3_LowMass',
        },
        'name':r'\verb|HLT_DoubleMu4_3_LowMass|',
    },
    {
        'cut':{
            'mm':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
            'fakes':'dstar_dm_pv>0.140 && dstar_dm_pv<0.155',
        },
        'name':'$\dm$ is in [0.140, 0.155] GeV',
    },
    {
        'cut':{
            'mm':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
            'fakes':'mm_kin_mass[dstar_mm_index]>1.81 && mm_kin_mass[dstar_mm_index]<1.94',
        },
        'name':'$\PDz$ mass is in [1.81, 1.94] GeV',
    },
    {
        'cut':{
            'mm':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
            'fakes':'mm_kin_vtx_prob[dstar_mm_index]>0.01',
        },
        'name':'$\PDz$ vertex probability $> 0.01$',
    },
    {
        'cut':{
            'mm':'dstar_pv_with_pion_prob>0.1',
            'fakes':'dstar_pv_with_pion_prob>0.1',
        },
        'name':'$\PDstpm$ vertex probability $> 0.1$',
    },
    {
        'cut':{
            'mm':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
            'fakes':'Muon_softMva[mm_mu1_index[dstar_mm_index]] > 0.45 && Muon_softMva[mm_mu2_index[dstar_mm_index]] > 0.45',
        },
        'name':'Muon identification',
    },
    {
        'cut':{
            'mm':'mm_kin_sl3d[dstar_mm_index]>3',
            'fakes':'mm_kin_sl3d[dstar_mm_index]>3',
        },
        'name':r'$\fls > 3$',
    },
    {
        'cut':{
            'mm':'mm_kin_alpha[dstar_mm_index]<0.1',
            'fakes':'mm_kin_alpha[dstar_mm_index]<0.1',
        },
        'name':'$\pa < 0.1$',
    },
    
]


total_counts = []
current_counts = []
final_counts = []

current_cuts = [""] * len(samples)

# function to get N and N-1 cuts
def get_final_cut(final_state, cut_to_exclude=None):
    cut = ""
    for entry in cuts:
        if final_state in entry['cut']:
            if cut_to_exclude and cut_to_exclude == entry['cut'][final_state]:
                continue
            if cut != "":  cut += "&&"
            cut += entry['cut'][final_state]
    #return "mm_mu1_index>=0 && mm_mu2_index>=0 && " + cut
    return cut

for sample in samples:
    chain = ROOT.TChain("Events")
    for f in sample['files']:
        chain.Add(f)
    sample['chain'] = chain

    n = chain.GetEntries()
    total_counts.append(n)
    current_counts.append(n)
    
    cut = get_final_cut(sample['final_state'])
    final_counts.append(chain.GetEntries(cut))


for entry in cuts:
    print("%35s " % entry['name'], end=' ')
    for i,sample in enumerate(samples):
        if sample['final_state'] in entry['cut']:
            baseline_cut = True
            if current_cuts[i] != "":
                current_cuts[i] += "&&"
                baseline_cut = False
            current_cuts[i] += entry['cut'][sample['final_state']]
            #print "cut ",current_cuts[i]
            n1 = sample['chain'].GetEntries(current_cuts[i])
            print("& %6.2f & %6.2f " % (100.0 * n1 / total_counts[i], 100.0 * n1 / current_counts[i]), end=' ')
            #print "& %6.2f & %6.2f & %6.2f" % (n1 , total_counts[i], current_counts[i]),
            current_counts[i] = n1
            if not baseline_cut:
                n2 = sample['chain'].GetEntries(get_final_cut(sample['final_state'], entry['cut'][sample['final_state']]))
                #print "cut ",get_final_cut(sample['final_state'], entry['cut'][sample['final_state']])
                print("& %6.2f " % (100.0 * final_counts[i] / n2), end=' ')
                #print "& %6.2f & %6.2f " % ( final_counts[i] , n2),
                
            else:
                print("&        ", end=' ')
        else:
            print("&        &         &        ", end=' ')
        
    print("\\\\")
