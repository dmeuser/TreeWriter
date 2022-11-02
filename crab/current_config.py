from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'ST_tW_top_5f_DS_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_2017'
config.section_('JobType')
config.JobType.maxJobRuntimeMin = 2400
config.JobType.pyCfgParams = ['dataset=/ST_tW_top_5f_DS_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM', 'user=dmeuser']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/TreeWriter/python/runTreeWriter_2017.py'
config.JobType.inputFiles = ['/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/data']
config.section_('Data')
config.Data.inputDataset = '/ST_tW_top_5f_DS_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'v04'
config.Data.publication = False
config.Data.unitsPerJob = 5
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_topUL/2017/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
