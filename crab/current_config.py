from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8_2018'
config.section_('JobType')
config.JobType.maxJobRuntimeMin = 2400
config.JobType.pyCfgParams = ['dataset=/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM', 'user=dmeuser']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/TreeWriter/python/runTreeWriter_2018.py'
config.JobType.inputFiles = ['/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/data']
config.section_('Data')
config.Data.inputDataset = '/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
config.Data.outputDatasetTag = 'v07'
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_topUL/2018/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
