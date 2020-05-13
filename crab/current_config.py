from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8'
config.section_('JobType')
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/TreeWriter/python/runTreeWriter.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['dataset=/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', 'user=dmeuser']
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
config.Data.outputDatasetTag = 'v09'
config.Data.publication = False
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_top/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
