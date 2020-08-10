from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8'
config.section_('JobType')
config.JobType.maxJobRuntimeMin = 2400
config.JobType.pyCfgParams = ['dataset=/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM', 'user=dmeuser']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTree/CMSSW_10_2_18/src/TreeWriter/TreeWriter/python/runTreeWriter.py'
config.JobType.inputFiles = ['/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTree/CMSSW_10_2_18/src/TreeWriter/data']
config.section_('Data')
config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'v17'
config.Data.publication = False
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_top/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
