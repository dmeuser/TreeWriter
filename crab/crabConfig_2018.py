#!/usr/bin/env python2
import getpass
import os
import subprocess
import CRABClient
from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration,saveConfigurationFile


def searchUserDatasets( name ):
    return subprocess.check_output( 'das_client.py --limit=0 --query="dataset={} instance=prod/phys03"'.format(name), shell=True ).split("\n")[0:-1]

def getLumiMask(cmssw_src):
    """Extracts the Lumi-Mask path from the makefile used to calculate the pileup reweighting"""
    out = ""
    with open(cmssw_src+"TreeWriter/PUreweighting/Makefile") as f:
        lines = f.readlines()
    for l in lines:
        if l.startswith("ANA_JSON2018="):
            out = l[13:-1]
    if not out:
        print "ERROR: could not find lumi mask"
    return out

def getRequestName(dataset, isSim):
    out = ""
    d1, d2, d3, d4 = dataset.split("/")
    if isSim:
        out = d2
        if "ext1" in dataset: out += "_ext"
        elif "ext2" in dataset: out += "_ext2"
        elif "ext3" in dataset: out += "_ext3"
        elif "ext4" in dataset: out += "_ext4"
        elif "ext5" in dataset: out += "_ext5"
        elif "backup" in dataset: out += "_backup"
        out += "_2018"
    else:
        out = "{}_{}".format(d2,d3)
    # CRABClient.ClientExceptions.ConfigurationException: Invalid CRAB configuration: Parameter General.requestName should not be longer than 100 characters.
    out = out[:100]
    return out

cmssw_src = os.environ['CMSSW_BASE']+'/src/'


datasets={}

datasets["DoubleMuon"] = [
    "/DoubleMuon/Run2018A-12Nov2019_UL2018-v2/MINIAOD",
    "/DoubleMuon/Run2018B-12Nov2019_UL2018-v2/MINIAOD",
    "/DoubleMuon/Run2018C-12Nov2019_UL2018-v2/MINIAOD",
    "/DoubleMuon/Run2018D-12Nov2019_UL2018-v3/MINIAOD",
    
    #  ~"/DoubleMuon/Run2018A-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/DoubleMuon/Run2018B-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/DoubleMuon/Run2018C-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/DoubleMuon/Run2018D-UL2018_MiniAODv2-v1/MINIAOD",
]

datasets["EGamma"] = [
    "/EGamma/Run2018A-12Nov2019_UL2018-v2/MINIAOD",
    "/EGamma/Run2018B-12Nov2019_UL2018-v2/MINIAOD",
    "/EGamma/Run2018C-12Nov2019_UL2018-v2/MINIAOD",
    "/EGamma/Run2018D-12Nov2019_UL2018-v4/MINIAOD",
    
    #  ~"/EGamma/Run2018A-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/EGamma/Run2018B-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/EGamma/Run2018C-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/EGamma/Run2018D-UL2018_MiniAODv2-v1/MINIAOD",
]

datasets["MuonEG"] = [
    "/MuonEG/Run2018A-12Nov2019_UL2018_rsb-v1/MINIAOD",
    "/MuonEG/Run2018B-12Nov2019_UL2018-v1/MINIAOD",
    "/MuonEG/Run2018C-12Nov2019_UL2018-v1/MINIAOD",
    "/MuonEG/Run2018D-12Nov2019_UL2018_rsb-v1/MINIAOD",
    
    #  ~"/MuonEG/Run2018A-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/MuonEG/Run2018B-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/MuonEG/Run2018C-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/MuonEG/Run2018D-UL2018_MiniAODv2-v1/MINIAOD",
]

datasets["SingleMuon"] = [
    "/SingleMuon/Run2018A-12Nov2019_UL2018_rsb-v1/MINIAOD",
    "/SingleMuon/Run2018B-12Nov2019_UL2018-v2/MINIAOD",
    "/SingleMuon/Run2018C-12Nov2019_UL2018-v2/MINIAOD",
    "/SingleMuon/Run2018D-12Nov2019_UL2018-v8/MINIAOD",
    
    #  ~"/SingleMuon/Run2018A-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/SingleMuon/Run2018B-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/SingleMuon/Run2018C-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/SingleMuon/Run2018D-UL2018_MiniAODv2-v1/MINIAOD",
]

datasets["MET"] = [
    "/MET/Run2018A-12Nov2019_UL2018-v3/MINIAOD",
    "/MET/Run2018B-12Nov2019_UL2018_rsb-v1/MINIAOD",
    "/MET/Run2018C-12Nov2019_UL2018_rsb-v1/MINIAOD",
    "/MET/Run2018D-12Nov2019_UL2018_rsb-v2/MINIAOD",
    
    #  ~"/MET/Run2018A-UL2018_MiniAODv2-v2/MINIAOD",
    #  ~"/MET/Run2018B-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/MET/Run2018C-UL2018_MiniAODv2-v1/MINIAOD",
    #  ~"/MET/Run2018D-UL2018_MiniAODv2-v1/MINIAOD",
]

datasets["Standard_ttbar"] = [
    "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM"
    
    #  ~"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",      #not finished
    #  ~"/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM"       #not finished
]

datasets["SingleTop"] = [
    "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",
    "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",
    "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    #missing:
    #"ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"
    #"ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"
    
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",        #not finished
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",        #not finished
    #  ~"/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",        #not finished
    #missing:
    #"ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"
    #"ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"
]

datasets["V+Jets"] = [
    "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAOD-Pilot_106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",     #correct with pilot?
    "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",
    "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",
    
    #  ~"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-Pilot_106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",     #correct with pilot?
    #  ~"/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
    #  ~"/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1_ext1-v1/MINIAODSIM",
    #  ~"/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
]

datasets["Diboson"] = [
    "/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    
    #  ~"/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",       #not finished
]

datasets["tt+X"] = [
    "/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM",
    "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",
    "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM",
    
    #  ~"//TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",        #not finished
    #  ~"/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",
    #  ~"/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM",     #not finished
]

datasets["dmeuser"] = []
datasets["dmeuser"] += datasets["DoubleMuon"]
datasets["dmeuser"] += datasets["EGamma"]
datasets["dmeuser"] += datasets["MuonEG"]
datasets["dmeuser"] += datasets["SingleMuon"]
datasets["dmeuser"] += datasets["MET"]
datasets["dmeuser"] += datasets["Standard_ttbar"]
datasets["dmeuser"] += datasets["SingleTop"]
datasets["dmeuser"] += datasets["V+Jets"]
datasets["dmeuser"] += datasets["Diboson"]
datasets["dmeuser"] += datasets["tt+X"]

# call with 'python crabConfig.py'
if __name__ == '__main__':
    user = getpass.getuser()
    
    print "!!!!!!!!!!!!check runTreeWriter for maxEvents!!!!!!!!!!!!!"

    iSub = 0
    failed = []
    for dataset in datasets[user]:

        isSim = 'SIM' in dataset
        isFastSim = "Fast" in dataset

        config = Configuration()

        config.section_("General")
        config.General.requestName = getRequestName(dataset, isSim)
        config.General.transferOutputs = True
        config.General.transferLogs = True

        config.section_("JobType")
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = cmssw_src + 'TreeWriter/TreeWriter/python/runTreeWriter_2018.py'
        config.JobType.pyCfgParams = ["dataset="+dataset,"user="+user]
        config.JobType.inputFiles  = [cmssw_src + "TreeWriter/" + x for x in ["data"]]
        config.JobType.allowUndistributedCMSSW = True
        #  ~config.JobType.maxJobRuntimeMin = 2400
        config.JobType.maxJobRuntimeMin = 1800


        config.section_("Data")
        config.Data.inputDataset = dataset
        config.Data.splitting = 'FileBased' if isSim else 'LumiBased'
        #  ~config.Data.unitsPerJob = 1 if isSim else 50
        config.Data.unitsPerJob = 5 if isSim else 50
        config.Data.publication = False
        config.Data.outputDatasetTag = 'outputDatasetTag'
        config.Data.outLFNDirBase = "outLFNDirBase"

        config.section_("Site")
        config.Site.storageSite = 'T2_DE_RWTH'

        if not isSim:
            config.Data.lumiMask = getLumiMask(cmssw_src)

        if user=="dmeuser":
            config.Data.outputDatasetTag = 'v03'
            config.Data.outLFNDirBase = "/store/user/dmeuser/run2_topUL/2018/"
        else:
            print "you shall not pass!"
            print "(unkown user '%s')"%user
            exit()
         
        saveConfigurationFile(config,"current_config.py")
        try:
            print "submitting",dataset
            # ~crabCommand('submit', config = config)
            subprocess.call("crab submit current_config.py", shell=True)
        except CRABClient.ClientExceptions.ConfigException,e:
            print "  ERROR:",e
            failed.append(dataset)
        except Exception,e:
            print "  ERROR:",e
            failed.append(dataset)
        else:
            iSub+=1
    print "submitted",iSub,"tasks"
    if failed:
        print "failed to submit:"
        print "   "+"\n   ".join(failed)
