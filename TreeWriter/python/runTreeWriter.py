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
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/60000/00E7F059-0BD5-E611-9267-001E67397CB5.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/02DD5E46-7ABE-E611-8F20-0025905B8582.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/BC527183-C0B7-E611-BC15-001E67348055.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/SMS-T5Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/80000/F227DD10-813E-E611-A722-6C3BE5B5C460.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/SMS-TChiWG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/000C22CC-092D-E611-8162-B083FED04D68.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/SMS-TChiNG_BF50N50G_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/0009D56C-8613-E711-BE86-0CC47A4D76A2.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/GMSB_GravitinoLSP_N1decays_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/130000/00E65FA0-5B0C-E711-BF18-0025905B8566.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016E/SinglePhoton/MINIAOD/03Feb2017-v1/80000/74BA03B0-62EB-E611-9F30-0025905A6082.root'
#options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016partBtoH-03Feb2017_MINIAOD.root'
#~ options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016C-03Feb2017-v1_MINIAOD.root'
#~ options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016D-03Feb2017-v1_MINIAOD.root'
#~ options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016E-03Feb2017-v1_MINIAOD.root'
#~ options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016F-03Feb2017-v1_MINIAOD.root'
#~ options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016G-03Feb2017-v1_MINIAOD.root'
#options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016H-03Feb2017_ver2-v1_MINIAOD.root'
#options.inputFiles = 'root:///user/kiesel/root-files/johannes/SinglePhoton_Run2016H-03Feb2017_ver3-v1_MINIAOD.root'
#~ options.inputFiles = 'root:///user/dmeuser/master/root-files/MINIAOD_80X_mcRun2_asymptotic_2016_miniAODv2_Fast_PAT_500_SIM.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/06B8CB51-31EB-E711-B445-0025905A6136.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/000C44EA-F8EC-E711-8145-0242AC130002.root'
options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/0E0EDE89-9AC4-E811-9484-AC1F6B1AF16E.root'
#~ options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/SMS-T5qqqqHg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00078287-CB49-E711-B388-FA163E1C5C83.root'
#~ options.inputFiles = 'root:///net/data_cms1b/user/dmeuser/master/root-files/tching.root'
options.outputFile = 'photonTree.root'
#~ options.outputFile = 'overlap_lepton_2.root'
#~ options.maxEvents = -1
options.maxEvents = 1000
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

# the actual TreeWriter module
process = cms.Process("TreeWriter")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


# determine global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if isRealData:
    if "Run2016H" in dataset:
        process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v16"
    else:
        process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"
else:
    if "80X_mcRun2_asymptotic_2016_TrancheIV_v6" in dataset:
        process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_TrancheIV_v8"
    elif "80X_mcRun2_asymptotic_2016_miniAODv2" in dataset:
        process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_miniAODv2_v1"
    else:
        print "Could not guess correct global tag for", dataset

######################
# PHOTONS, ELECTRONS #
######################
# Regression: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)

process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(81),
        engineName = cms.untracked.string('TRandom3'),
    ),
    calibratedPatPhotons = cms.PSet(
        initialSeed = cms.untracked.uint32(81),
        engineName = cms.untracked.string('TRandom3'),
    )
)

process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt > 5 && abs(eta)<2.5")
)
process.selectedPhotons = cms.EDFilter("PATPhotonSelector",
    src = cms.InputTag("slimmedPhotons"),
    cut = cms.string("pt > 5 && abs(eta)<2.5")
)
process.calibratedPatPhotons.isMC = not isRealData
process.calibratedPatElectrons.isMC = not isRealData
process.calibratedPatPhotons.photons = "selectedPhotons"
process.calibratedPatElectrons.electrons = "selectedElectrons"

process.EGMRegression = cms.Path(process.regressionApplication)
process.EGMSmearerElectrons = cms.Path(process.calibratedPatElectrons)
process.EGMSmearerPhotons = cms.Path(process.calibratedPatPhotons)


# Identification
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import *

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD

# turn on VID producer, indicate data format to be DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer  (process, dataFormat)

# define which IDs we want to produce
el_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
ph_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                 'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
#                 ,'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']

#add them to the VID producer
for idmod in el_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.photonIDValueMapProducer.srcMiniAOD = "calibratedPatPhotons"
#process.photonMVAValueMapProducer.srcMiniAOD = "calibratedPatPhotons"
process.egmPhotonIDs.physicsObjectSrc = "calibratedPatPhotons"
process.egmGsfElectronIDs.physicsObjectSrc = "calibratedPatElectrons"


##########################
# Jet Energy Corrections #
##########################
# where .db files are placed (e.g. for JEC, JER)
# Crab will always be in the $CMSSW_BASE directory, so to run the code locally,
# a symbolic link is added
#~ if not os.path.exists("src"): os.symlink(os.environ["CMSSW_BASE"]+"/src/", "src")
if "Fast" in dataset:
    dbPath = 'sqlite:'

    from CondCore.CondDB.CondDB_cfi import CondDB
    CondDB.__delattr__('connect')

    process.jec = cms.ESSource('PoolDBESSource',
        CondDB,
        connect = cms.string(dbPath+'Spring16_25nsFastSimV1_MC.db'),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsFastSimMC_V1_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
            )
        )
    )
    # Add an ESPrefer to override JEC that might be available from the global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')


jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
if isRealData: jecLevels.append('L2L3Residual')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None')
)


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
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")


################################
# Define input and output      #
################################
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))#, eventsToProcess = cms.untracked.VEventRange("1:29986:5916215"))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))


################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(0),
                                    photon_pT_cut=cms.untracked.double(20), # for leading photon
                                    jet_pT_cut=cms.untracked.double(30), # for all jets
                                    isolatedPhotons=cms.untracked.bool(True), # for all photons in the collection
                                    minNumberPhotons_cut=cms.untracked.uint32(1),
                                    minNumberElectrons_cut=cms.untracked.uint32(0),
                                    minNumberBinos_cut=cms.untracked.uint32(0),
                                    # physics objects
                                    photons = cms.InputTag("calibratedPatPhotons"),
                                    jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("calibratedPatElectrons"),
                                    mets = cms.InputTag("slimmedMETs", "", "TreeWriter"),
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
                                    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
                                    # photon IDs
                                    photonLooseId15Map   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
                                    photonMediumId15Map  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
                                    photonTightId15Map   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
                                    photonLooseIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
                                    photonMediumIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
                                    photonTightIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),
#                                    photonMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
                                    # met filters to apply
                                    metFilterNames=cms.untracked.vstring(
                                        "Flag_HBHENoiseFilter",
                                        "Flag_HBHENoiseIsoFilter",
                                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                                        "Flag_goodVertices",
                                        "Flag_eeBadScFilter",
                                        "Flag_globalTightHalo2016Filter",

                                    ),
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU"),
                                    hardPUveto=cms.untracked.bool(False),
                                    reMiniAOD=cms.untracked.bool(False),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    triggerObjectNames=cms.vstring(),
                                    pfJetIDSelector=cms.PSet(version=cms.string('FIRSTDATA'), quality=cms.string('LOOSE')),
                                    triggerPrescales=cms.vstring(), # also useful to check whether a trigger was run
                                    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
                                    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
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
if "Fast" in dataset:
    process.TreeWriter.metFilterNames.remove("Flag_globalTightHalo2016Filter")
    process.TreeWriter.lheEventProduct = "source"
    if "T5Wg" in dataset or "T6Wg" in dataset:
        process.TreeWriter.minNumberBinos_cut = 1

if "PUMoriond17" in dataset or "GGM_GravitinoLSP_M1" in dataset:
    process.TreeWriter.pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_Moriond17MC_PoissonOOTPU")

# determine user if not set by crab
user=options.user or getpass.getuser()
# user settings
if user=="kiesel":
    process.TreeWriter.HT_cut=500.
    process.TreeWriter.photon_pT_cut=90.
    process.TreeWriter.minNumberPhotons_cut=0
    process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter"]
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT600_v",
        "HLT_Photon90_v",
        "HLT_PFHT600_v",
        "HLT_PFHT800_v",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v",
        "HLT_PFJet450_v",
        "HLT_AK8PFJet450_v",
    ]
    process.TreeWriter.triggerPrescales=process.TreeWriter.triggerNames
    if "SingleElectron" in dataset or "DY" in dataset:
        process.TreeWriter.triggerObjectNames = ["hltEle27erWPLooseGsfTrackIsoFilter", "hltEle27erWPTightGsfTrackIsoFilter", "hltEG90CaloIdLHEFilter"]
        process.TreeWriter.triggerNames = ["HLT_Ele27_eta2p1_WPLoose_Gsf_v", "HLT_Ele27_eta2p1_WPTight_Gsf_v"]
        process.TreeWriter.HT_cut = 0.
        process.TreeWriter.photon_pT_cut = 25.
        process.TreeWriter.minNumberPhotons_cut = 1
    if "Fast" in dataset: # signal scan
        process.TreeWriter.HT_cut = 0.

elif user=="jschulz" or user=="dmeuser":
    process.TreeWriter.photon_pT_cut=100
    process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter", "hltEG165HE10Filter"]
    if "Fast" in dataset:
        process.TreeWriter.minNumberPhotons_cut=0
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT600_v",
        "HLT_Photon90_CaloIdL_PFHT500_v",
        "HLT_PFHT600_v",
        'HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
        'HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        'HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        'HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v",
        "HLT_Photon22_v",
        "HLT_Photon30_v",
        "HLT_Photon36_v",
        "HLT_Photon50_v",
        "HLT_Photon75_v",
        "HLT_Photon90_v",
        "HLT_Photon120_v",
        "HLT_Photon36_R9Id90_HE10_IsoM_v",
        "HLT_Photon165_R9Id90_HE10_IsoM_v",
        "HLT_Photon165_HE10_v",
        "HLT_Photon135_PFMET100_JetIdCleaned_v", # used in early data taking
        "HLT_Photon135_PFMET100_v",              # used in later data taking
        "HLT_Photon135_PFMET100_NoiseCleaned_v", # used for MC
        "HLT_Photon175_v",
        "HLT_Photon500_v",
        "HLT_PFMET170_NoiseCleaned_v",
        "HLT_PFMET170_HBHECleaned_v",
        "HLT_PFMET170_JetIdCleaned_v",
        "HLT_PFMET170_NotCleaned_v",
        "HLT_IsoMu18_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu22_v",
        "HLT_Mu20_v",
        "HLT_Mu45_eta2p1_v",
        "HLT_Mu50_v",
        "HLT_Mu30_TkMu11_v",
        "HLT_DoubleIsoMu17_eta2p1_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_Photon22_CaloIdL_L1ISO_v",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
         # HT
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
        #Lepton
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v",
        "HLT_Mu17_Photon30_CaloIdL_L1ISO_v",
        "HLT_Mu38NoFilterNoVtx_Photon38_CaloIdL",
    ]
    process.TreeWriter.triggerPrescales=[
        "HLT_Photon135_PFMET100_JetIdCleaned_v", # used in early data taking
        "HLT_Photon135_PFMET100_v",              # used in later data taking
        "HLT_Photon135_PFMET100_NoiseCleaned_v", # used for MC
        "HLT_Photon22_v",
        "HLT_Photon30_v",
        "HLT_Photon36_v",
        "HLT_Photon50_v",
        "HLT_Photon75_v",
        "HLT_Photon90_v",
        "HLT_Photon120_v",
        "HLT_Photon36_R9Id90_HE10_IsoM_v",
        "HLT_IsoMu18_v",
        "HLT_Mu20_v",
        "HLT_Mu17_Photon22_CaloIdL_L1ISO_v",
        # HT
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
    process.BadPFMuonFilter
    *process.BadChargedCandidateFilter
    *process.TreeWriter
)
