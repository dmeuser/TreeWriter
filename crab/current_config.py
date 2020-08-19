from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'SingleElectron_Run2016H-17Jul2018-v1'
config.section_('JobType')
config.JobType.maxJobRuntimeMin = 1200
config.JobType.pyCfgParams = ['dataset=/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD', 'user=dmeuser']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTree/CMSSW_10_2_18/src/TreeWriter/TreeWriter/python/runTreeWriter.py'
config.JobType.inputFiles = ['/.automount/home/home__home4/institut_1b/dmeuser/CMSSW_ttbarTree/CMSSW_10_2_18/src/TreeWriter/data']
config.section_('Data')
config.Data.inputDataset = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
config.Data.outputDatasetTag = 'v18'
config.Data.publication = False
config.Data.unitsPerJob = 50
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.outLFNDirBase = '/store/user/dmeuser/run2_top/'
config.section_('Site')
config.Site.storageSite = 'T2_DE_RWTH'
