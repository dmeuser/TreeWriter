import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os,re
import getpass

def guessDatasetFromFileName(filename):
    # This reproduces the dataset roughly if the file is on /store
    # Not reproduced are e.g. the pileup scenario
    # For local files, specify your own rules or run it with the 'dataset' option
    nParts = filename.split("/")
    if "store" in nParts and len(nParts)>6:
        nParts = nParts[nParts.index("store"):]
        return "/{}/{}-{}/{}".format(nParts[3], nParts[2], nParts[5], nParts[4])
    if "user" in nParts:
        return nParts[-1].replace(".root", "")
    return filename

options = VarParsing ('analysis')
options.register ('dataset',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name of the dataset, used to do further settings")
options.register ('user',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name the user. If not set by crab, this script will determine it.")

# defaults
options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/001BCE65-30C3-E811-A357-A4BF0112DD3C.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016D/MuonEG/MINIAOD/17Jul2018-v1/80000/D48E0B3A-A28E-E811-8A12-7CD30AD08EB4.root'
options.outputFile = 'ttbarTree.root'
#~ options.outputFile = 'overlap_lepton_2.root'
options.maxEvents = -1
#~ options.maxEvents = 1000
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

# the actual TreeWriter module
process = cms.Process("TreeWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


# determine global tag here only 2016
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if isRealData:
        process.GlobalTag.globaltag = "94X_dataRun2_v10"
else:
        process.GlobalTag.globaltag = "94X_mcRun2_asymptotic_v3"
        

#############
# ELECTRONS #
#############
# Geometry neccessary to run setupEgammaPostRecoSeq
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                       era='2016-Legacy')
                       

##########################
# MET                    #
##########################
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(
    process,
    isData=isRealData,
)


################################
# MET Filter                   #
################################
#~ process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#~ process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#~ process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")


################################
# Define input and output      #
################################
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
#~ process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), lumisToProcess = cms.untracked.VLuminosityBlockRange("276315:134"))#eventsToProcess = cms.untracked.VEventRange("276315:134:185994346"))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))#, eventsToProcess = cms.untracked.VEventRange("276315:134:186839926"))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))


################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(0),
                                    jet_pT_cut=cms.untracked.double(14), # for all jets
                                    isolatedPhotons=cms.untracked.bool(True), # for all photons in the collection
                                    minNumberElectrons_cut=cms.untracked.uint32(0),
                                    NumberLeptons_cut=cms.untracked.uint32(2),
                                    # physics objects
                                    jets = cms.InputTag("slimmedJets"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    metCorr = cms.InputTag(""),
                                    metCorrCal = cms.InputTag(""),
                                    caloMets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    pileUpSummary = cms.InputTag('slimmedAddPileupInfo'),
                                    lheEventProduct = cms.InputTag('externalLHEProducer'),
                                    packedCandidates=cms.InputTag("packedPFCandidates"),
                                    # electron IDs
                                    electronVetoIdMap   = cms.untracked.string("cutBasedElectronID-Fall17-94X-V1-veto"),
                                    electronLooseIdMap  = cms.untracked.string("cutBasedElectronID-Fall17-94X-V1-loose"),
                                    electronMediumIdMap = cms.untracked.string("cutBasedElectronID-Fall17-94X-V1-medium"),
                                    electronTightIdMap  = cms.untracked.string("cutBasedElectronID-Fall17-94X-V1-tight"),
                                    # met filters to apply
                                    metFilterNames=cms.untracked.vstring(
                                        "Flag_HBHENoiseFilter",
                                        "Flag_HBHENoiseIsoFilter",
                                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                                        "Flag_goodVertices",
                                        "Flag_globalSuperTightHalo2016Filter",
                                        "Flag_BadPFMuonFilter",
                                        "Flag_eeBadScFilter",

                                    ),
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    # choose pileup data
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_Moriond17MC_PoissonOOTPU"),
                                    hardPUveto=cms.untracked.bool(False),
                                    reMiniAOD=cms.untracked.bool(False),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    triggerObjectNames=cms.vstring(),
                                    pfJetIDSelector=cms.PSet(version=cms.string('FIRSTDATA'), quality=cms.string('LOOSE')),
                                    triggerPrescales=cms.vstring(), # also useful to check whether a trigger was run
                                    metCorrected = cms.InputTag("slimmedMETs"),
                                    metCalibrated = cms.InputTag("slimmedMETs")
)

################################
# Modify the TreeWriter module #
################################

process.TreeWriter.hardPUveto=dataset.startswith("/QCD_HT100to200")

if "03Feb2017" in dataset:
    process.TreeWriter.reMiniAOD = True
    process.TreeWriter.mets = cms.InputTag("slimmedMETsMuEGClean", "", "PAT")
    process.TreeWriter.metCorrected = cms.InputTag("slimmedMETsMuEGClean", "", "TreeWriter")
    process.TreeWriter.metCalibrated = cms.InputTag("slimmedMETsMuEGCleanCalibrated", "", "TreeWriter")
    process.TreeWriter.metFilterNames.extend(["Flag_chargedHadronTrackResolutionFilter", "Flag_muonBadTrackFilter"])

    # Now you are creating the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
    from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
    corMETFromMuonAndEG(
        process,
        pfCandCollection = "",
        electronCollection = "slimmedElectronsBeforeGSFix",
        photonCollection = "slimmedPhotonsBeforeGSFix",
        corElectronCollection = "slimmedElectrons",
        corPhotonCollection = "slimmedPhotons",
        allMETEGCorrected = True,
        muCorrection = False,
        eGCorrection = True,
        runOnMiniAOD = True,
        postfix = "MuEGClean"
    )
    process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
    process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
    process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
    process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    del process.slimmedMETsMuEGClean.caloMET

    corMETFromMuonAndEG(
        process,
        pfCandCollection = "",
        electronCollection = "slimmedElectronsBeforeGSFix",
        photonCollection = "slimmedPhotonsBeforeGSFix",
        corElectronCollection = "calibratedPatElectrons",
        corPhotonCollection = "calibratedPatPhotons",
        allMETEGCorrected = True,
        muCorrection = False,
        eGCorrection = True,
        runOnMiniAOD = True,
        postfix = "MuEGCleanCalibrated"
    )
    process.slimmedMETsMuEGCleanCalibrated = process.slimmedMETs.clone()
    process.slimmedMETsMuEGCleanCalibrated.src = cms.InputTag("patPFMetT1MuEGClean")
    process.slimmedMETsMuEGCleanCalibrated.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
    process.slimmedMETsMuEGCleanCalibrated.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    del process.slimmedMETsMuEGCleanCalibrated.caloMET

if not isRealData:
    process.TreeWriter.metFilterNames.remove("Flag_eeBadScFilter")
#~ if "Fast" in dataset:
    #~ process.TreeWriter.metFilterNames.remove("Flag_globalTightHalo2016Filter")
    #~ process.TreeWriter.lheEventProduct = "source"
    #~ if "T5Wg" in dataset or "T6Wg" in dataset:
        #~ process.TreeWriter.minNumberBinos_cut = 1

# determine user if not set by crab
user=options.user or getpass.getuser()
# user settings
if user=="jschulz" or user=="dmeuser":
    process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter", "hltEG165HE10Filter"]
    process.TreeWriter.triggerNames=[
        # ee Channel
        #2016
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Ele27_WPTight_Gsf_v",
        "HLT_IsoMu24_v",
        "HLT_IsoTkMu24_v",
        # emu Channel
        #2016
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        # mumu Chanel
        #2016
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        # HT 2016
        "HLT_PFHT125_v",
        "HLT_PFHT200_v",
        "HLT_PFHT250_v",
        "HLT_PFHT300_v",
        "HLT_PFHT350_v",
        "HLT_PFHT400_v",
        "HLT_PFHT475_v",
        "HLT_PFHT600_v",
        "HLT_PFHT650_v",
        "HLT_PFHT800_v",
        # MET 2016
        "HLT_PFMET110_PFMHT110_IDTight_v",
        "HLT_PFMET120_PFMHT120_IDTight_v",
        "HLT_PFMET170_NoiseCleaned_v",
        "HLT_PFMET170_HBHECleaned_v",
        "HLT_PFMET170_JetIdCleaned_v",
        "HLT_PFMET170_NotCleaned_v",
        "HLT_PFMET300_v",
        "HLT_PFMET400_v",
        "HLT_PFMET500_v",
        "HLT_PFMET600_v",
    ]
    process.TreeWriter.triggerPrescales=[
        # HT 2016
        "HLT_PFHT125_v",
        "HLT_PFHT200_v",
        "HLT_PFHT250_v",
        "HLT_PFHT300_v",
        "HLT_PFHT350_v",
        "HLT_PFHT400_v",
        "HLT_PFHT475_v",
        "HLT_PFHT600_v",
        "HLT_PFHT650_v",
        "HLT_PFHT800_v",
        # MET 2016
        "HLT_PFMET110_PFMHT110_IDTight_v",
        "HLT_PFMET120_PFMHT120_IDTight_v",
        "HLT_PFMET170_NoiseCleaned_v",
        "HLT_PFMET170_HBHECleaned_v",
        "HLT_PFMET170_JetIdCleaned_v",
        "HLT_PFMET170_NotCleaned_v",
        "HLT_PFMET300_v",
        "HLT_PFMET400_v",
        "HLT_PFMET500_v",
        "HLT_PFMET600_v",
    ]
else:
    print "you shall not pass!"
    print "(unkown user '%s')"%user
    exit()

for trig in process.TreeWriter.triggerPrescales:
    assert(trig in process.TreeWriter.triggerNames),"Trigger '"+trig+"' is not used, so prescale cannot be stored!"

####################
#     RUN          #
####################

process.p = cms.Path(
    #~ process.BadPFMuonFilter*
    process.egammaPostRecoSeq
    *process.TreeWriter
)
