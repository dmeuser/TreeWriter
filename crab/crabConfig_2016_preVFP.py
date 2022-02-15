#!/usr/bin/env python2
import getpass
import os
import subprocess
import CRABClient
from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration


def searchUserDatasets( name ):
    return subprocess.check_output( 'das_client.py --limit=0 --query="dataset={} instance=prod/phys03"'.format(name), shell=True ).split("\n")[0:-1]

def getLumiMask(cmssw_src):
    """Extracts the Lumi-Mask path from the makefile used to calculate the pileup reweighting"""
    out = ""
    with open(cmssw_src+"TreeWriter/PUreweighting/Makefile") as f:
        lines = f.readlines()
    for l in lines:
        if l.startswith("ANA_JSON2016="):
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
        out += "_2016_preVFP"
    else:
        out = "{}_{}".format(d2,d3)
    # CRABClient.ClientExceptions.ConfigurationException: Invalid CRAB configuration: Parameter General.requestName should not be longer than 100 characters.
    out = out[:100]
    return out

cmssw_src = os.environ['CMSSW_BASE']+'/src/'


datasets={}

datasets["DoubleMuon"] = [
    # "/DoubleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD",
    "/DoubleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD",
    "/DoubleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/DoubleMuon/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/DoubleMuon/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/DoubleMuon/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",

    # "/DoubleMuon/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleMuon/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleMuon/Run2016C-HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleMuon/Run2016D-HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleMuon/Run2016E-HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleMuon/Run2016F-HIPM_UL2016_MiniAODv2-v1/MINIAOD",

]

datasets["DoubleEG"] = [
    # "/DoubleEG/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD",
    "/DoubleEG/Run2016B-21Feb2020_ver2_UL2016_HIPM-v2/MINIAOD",
    "/DoubleEG/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/DoubleEG/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/DoubleEG/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/DoubleEG/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    
    # "/DoubleEG/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleEG/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v3/MINIAOD",
    # "/DoubleEG/Run2016C-HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleEG/Run2016D-HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleEG/Run2016E-HIPM_UL2016_MiniAODv2-v1/MINIAOD",
    # "/DoubleEG/Run2016F-HIPM_UL2016_MiniAODv2-v1/MINIAOD",    

]
    
datasets["SingleElectron"] = [
    # "/SingleElectron/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD",
    "/SingleElectron/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD",
    "/SingleElectron/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/SingleElectron/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/SingleElectron/Run2016E-21Feb2020_UL2016_HIPM-v2/MINIAOD",
    "/SingleElectron/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",

    # "/SingleElectron/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleElectron/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleElectron/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleElectron/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleElectron/Run2016E-HIPM_UL2016_MiniAODv2-v5/MINIAOD",
    # "/SingleElectron/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD",

]

datasets["MuonEG"] = [
    # "/MuonEG/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD",
    "/MuonEG/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD",
    "/MuonEG/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/MuonEG/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/MuonEG/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/MuonEG/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",

    # "/MuonEG/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MuonEG/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MuonEG/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MuonEG/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MuonEG/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MuonEG/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD",

]

datasets["SingleMuon"] = [
    # "/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD",
    "/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD",
    "/SingleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/SingleMuon/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/SingleMuon/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/SingleMuon/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",

    # "/SingleMuon/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleMuon/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleMuon/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleMuon/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleMuon/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/SingleMuon/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    
]

datasets["MET"] = [
    # "/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD",
    "/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD",
    "/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/MET/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD",
    "/MET/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD",

    # "/MET/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MET/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MET/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MET/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MET/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD",
    # "/MET/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD",

]

datasets["Standard_ttbar"] = [
    "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM"
    
    # "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM" # production
    
]

datasets["SingleTop"] = [
    "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/MINIAODSIM",
    "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/MINIAODSIM",
    
    # "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3/MINIAODSIM",
    # "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3/MINIAODSIM",
    
]

datasets["V+Jets"] = [
    "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",     
    "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
        
    #"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    #"/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    #"/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-FlatPU0to75_106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    #"/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    #"/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM",


]

datasets["Diboson"] = [
    "/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    
    #"/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    #"/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    #"/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
]

datasets["tt+X"] = [
    "/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",

    # "/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM",
    # "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM",
    # "/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM",
    
]

datasets["hdamp_ttbar"] = [
    "/TTTo2L2Nu_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",

    # "/TTTo2L2Nu_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTToHadronic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToHadronic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    
]

datasets["UE_top"] = [
    "/TTTo2L2Nu_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTTo2L2Nu_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    # missing:
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5down_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5up_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_TuneCP5down_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_TuneCP5up_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",

    # "/TTTo2L2Nu_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTTo2L2Nu_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTToHadronic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTToHadronic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",    
    
]

datasets["ColorRec_top"] = [
    "/TTTo2L2Nu_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTTo2L2Nu_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTTo2L2Nu_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    # missing:
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_erdON_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5CR1_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5CR2_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_erdON_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_TuneCP5CR1_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_TuneCP5CR2_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    
    
    # "/TTTo2L2Nu_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTTo2L2Nu_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTTo2L2Nu_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToHadronic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToHadronic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    
]

datasets["TopMass_top"] = [
    "/TTTo2L2Nu_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTTo2L2Nu_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToHadronic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    # missing:
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_mtop1695_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_antitop_5f_NoFullyHadronicDecays_mtop1755_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_mtop1695_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_tW_top_5f_NoFullyHadronicDecays_mtop1755_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_mtop1695_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_antitop_4f_InclusiveDecays_mtop1755_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_mtop1695_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    #  ~"/ST_t-channel_top_4f_InclusiveDecays_mtop1755_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    
    # "/TTTo2L2Nu_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production
    # "/TTTo2L2Nu_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToHadronic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToHadronic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM",
    # "/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", # production

]

datasets["all"] = []
datasets["all"] += datasets["DoubleMuon"]
datasets["all"] += datasets["DoubleEG"]
datasets["all"] += datasets["SingleElectron"]
datasets["all"] += datasets["MuonEG"]
datasets["all"] += datasets["SingleMuon"]
datasets["all"] += datasets["MET"]
datasets["all"] += datasets["Standard_ttbar"]
datasets["all"] += datasets["SingleTop"]
datasets["all"] += datasets["V+Jets"]
datasets["all"] += datasets["Diboson"]
datasets["all"] += datasets["tt+X"]
datasets["all"] += datasets["hdamp_ttbar"]
datasets["all"] += datasets["UE_top"]
datasets["all"] += datasets["ColorRec_top"]
datasets["all"] += datasets["TopMass_top"]

# call with 'python crabConfig.py'
if __name__ == '__main__':
    user = getpass.getuser()
    
    print "!!!!!!!!!!!!check runTreeWriter for maxEvents!!!!!!!!!!!!!"

    iSub = 0
    failed = []

    for dataset in datasets["all"]:
        isSim = 'SIM' in dataset
        isFastSim = "Fast" in dataset

        config = Configuration()

        config.section_("General")
        config.General.requestName = getRequestName(dataset, isSim)
        config.General.transferOutputs = True
        config.General.transferLogs = True

        config.section_("JobType")
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = cmssw_src + 'TreeWriter/TreeWriter/python/runTreeWriter_2016_preVFP.py'
        config.JobType.pyCfgParams = ["dataset="+dataset,"user="+user]
        config.JobType.inputFiles  = [cmssw_src + "TreeWriter/" + x for x in ["data"]]
        config.JobType.allowUndistributedCMSSW = True
        config.JobType.maxJobRuntimeMin = 2400
        #  ~config.JobType.maxJobRuntimeMin = 1800


        config.section_("Data")
        config.Data.inputDataset = dataset
        config.Data.splitting = 'FileBased' if isSim else 'LumiBased'
        config.Data.unitsPerJob = 5 if isSim else 50
        config.Data.publication = False
        config.Data.outputDatasetTag = 'outputDatasetTag'
        config.Data.outLFNDirBase = "outLFNDirBase"

        config.section_("Site")
        config.Site.storageSite = 'T2_DE_RWTH'

        if not isSim:
            config.Data.lumiMask = getLumiMask(cmssw_src)

        from TreeWriter.paths.getPath import getPath
        config.Data.outputDatasetTag = getPath("datasetTag_2016_preVFP", user)
        config.Data.outLFNDirBase = os.path.join(getPath("outLFNDirBase", user),getPath("dir2016_preVFP", user))
        
         
        with open("current_config.py", 'w') as handle:
            handle.write(str(config))
            
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