from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_2018'
config.section_('JobType')
config.JobType.maxJobRuntimeMin = 2400
config.JobType.pyCfgParams = ['dataset=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM', 'user=dmeuser']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_20/src/TreeWriter/TreeWriter/python/runTreeWriter_2018.py'
config.JobType.inputFiles = ['/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_20/src/TreeWriter/data']
config.section_('Data')
config.Data.inputDataset = '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM'
config.Data.outputDatasetTag = 'v06'
config.Data.publication = False
config.Data.unitsPerJob = 100000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_topUL/2018/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
