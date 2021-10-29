from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'EGamma_Run2018D-12Nov2019_UL2018-v8'
config.section_('JobType')
config.JobType.maxJobRuntimeMin = 1800
config.JobType.pyCfgParams = ['dataset=/EGamma/Run2018D-12Nov2019_UL2018-v8/MINIAOD', 'user=dmeuser']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_20/src/TreeWriter/TreeWriter/python/runTreeWriter_2018.py'
config.JobType.inputFiles = ['/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_20/src/TreeWriter/data']
config.section_('Data')
config.Data.inputDataset = '/EGamma/Run2018D-12Nov2019_UL2018-v8/MINIAOD'
config.Data.outputDatasetTag = 'v04'
config.Data.publication = False
config.Data.unitsPerJob = 50
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/home/home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_topUL/2018/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
