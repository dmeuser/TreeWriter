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
#  ~options.inputFiles =    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/00000/531C1968-9806-4346-834C-2A1EE1A86AEB.root',
#  ~options.inputFiles =    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/04A0B676-D63A-6D41-B47F-F4CF8CBE7DB8.root',      # miniAODv2
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2018B/MuonEG/MINIAOD/12Nov2019_UL2018-v1/100000/00BE9C7C-F659-EB4C-A6C4-EAC5054243B2.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/270002/B4155943-3E29-6849-BF8D-227F49F72A15.root',
options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2018A/DoubleMuon/MINIAOD/12Nov2019_UL2018-v2/100000/A988574D-C798-6246-AFE0-435CD7A9C836.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAOD/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/00000/0095F60C-3B5B-C44F-A15B-E923A692B6F5.root',
options.outputFile = 'ttbarTree.root'
#~ options.outputFile = 'overlap_lepton_2.root'
options.maxEvents = -1
#  ~options.maxEvents = 1000
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

# the actual TreeWriter module
process = cms.Process("TreeWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# determine global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if isRealData:
        process.GlobalTag.globaltag = "106X_dataRun2_v32"
        #  ~process.GlobalTag.globaltag = "106X_dataRun2_v33"       #miniAODv2 (correct? or only for production?)
else:
        process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v15_L1v1"
        #  ~process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v16_L1v1"     #miniAODv2 (correct? or only for production?)
        
#timing information
process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(True)
)

#############
# ELECTRONS #
#############
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
# Geometry neccessary to runVID
#  ~process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#  ~process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#  ~process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
#  ~process.load("Configuration.Geometry.GeometryECALHCAL_cff")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True) needed for Puppi correction
                       era='2018-UL') 
                       

#############
# MUONS     #
#############
from TreeWriter.TreeWriter.muonPATUserDataRochesterCorrectionAdder_cfi import muonPATUserDataRochesterCorrectionAdder
process.MuonsAddedRochesterCorr = muonPATUserDataRochesterCorrectionAdder.clone(
   src = "slimmedMuons",
   applyEnergyCorrections = False,
   debug = False,
)
process.MuonsAddedRochesterCorr.path = cms.string('TreeWriter/data/2018/RoccoR2018UL.txt')

##########################
# MET                    #
##########################
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(   #update pfMET 
    process,
    isData=isRealData,
    fixEE2017 = False,
    #  ~jetCollUnskimmed="updatedPatJetsUpdatedJECID",
    #  ~reapplyJEC=False,
)
# If using patSmearedJETS is removing too many jets due to lepton matching, use the following to keep all jets: https://indico.cern.ch/event/882528/contributions/3718330/attachments/1974802/3286397/Sync_summary.pdf
# Puppi MET is correct when applying new Puppi Tune

#########################################
# MET Filter (not needed in MiniAODv2   #
#########################################
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Recipe_for_BadPFMuonDz_filter_in
from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
)


################################
# Lepton Scale Factors         #
################################
LeptonFullSimScaleFactorMapPars = cms.PSet(
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2018/EfficienciesStudies_UL2018_DEN_TrackerMuons_rootfiles_Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2018/EfficienciesStudies_UL2018_DEN_TrackerMuons_rootfiles_Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root'),
    
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2018/egammaEffi.txt_Ele_Tight_EGM2D.root'),
    
    dataMCScaleFactorHisto_mu_ID = cms.string('NUM_TightID_DEN_TrackerMuons_abseta_pt'),
    dataMCScaleFactorHisto_mu_Iso = cms.string('NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt'),
    
    dataMCScaleFactorHisto_ele_ID = cms.string('EGamma_SF2D'),
)

################################
# Jets                         #
################################
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(        #JEC
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.PFJetIDTightLepVeto = cms.EDProducer('PatJetIDValueMapProducer',
    src = cms.InputTag("updatedPatJetsUpdatedJEC"),
    filterParams = cms.PSet(
        version = cms.string('RUNIIULCHS'),
        quality = cms.string('TIGHTLEPVETO')
    )
)
process.updatedPatJetsUpdatedJECID = cms.EDProducer('PATJetUserDataEmbedder',
    src = cms.InputTag("updatedPatJetsUpdatedJEC"),
    userInts = cms.PSet(
        PFJetIDTightLepVeto = cms.InputTag('PFJetIDTightLepVeto')
    ),
)
process.jetIDSequence = cms.Sequence(process.PFJetIDTightLepVeto * process.updatedPatJetsUpdatedJECID)

#from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/test/runJERsmearingOnMiniAOD.py
process.updatedPatJetsUpdatedJECIDsmeared = cms.EDProducer('SmearedPATJetProducer',        #JER
       src = cms.InputTag('updatedPatJetsUpdatedJECID'),
       enabled = cms.bool(True),
       rho = cms.InputTag("fixedGridRhoFastjetAll"),
       algo = cms.string('AK4PFchs'),
       algopt = cms.string('AK4PFchs_pt'),
       #resolutionFile = cms.FileInPath('Autumn18_V7_MC_PtResolution_AK4PFchs.txt'),
       #scaleFactorFile = cms.FileInPath('combined_SFs_uncertSources.txt'),

       genJets = cms.InputTag('slimmedGenJets'),
       dRMax = cms.double(0.2),
       dPtMaxFactor = cms.double(3),

       debug = cms.untracked.bool(True),
   # Systematic variation
   # 0: Nominal
   # -1: -1 sigma (down variation)
   # 1: +1 sigma (up variation)
   variation = cms.int32(0),  # If not specified, default to 0
   uncertaintySource = cms.string(""), # If not specified, default to Total
       )

################################
# Puppi Jets                   #
################################
from CommonTools.PileupAlgos.customizePuppiTune_cff import UpdatePuppiTuneV15
UpdatePuppiTuneV15(process, runOnMC=(isRealData==False))

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsPuppi'),
   labelName = 'UpdatedJECPuppi',
   jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
)
process.jecSequencePuppi = cms.Sequence(process.patJetCorrFactorsUpdatedJECPuppi * process.updatedPatJetsUpdatedJECPuppi)

process.PFJetIDTightLepVetoPuppi = process.PFJetIDTightLepVeto.clone()
process.PFJetIDTightLepVetoPuppi.src = cms.InputTag("updatedPatJetsUpdatedJECPuppi")

process.updatedPatJetsUpdatedJECIDPuppi = process.updatedPatJetsUpdatedJECID.clone()
process.updatedPatJetsUpdatedJECIDPuppi.src = cms.InputTag("updatedPatJetsUpdatedJECPuppi")
process.updatedPatJetsUpdatedJECIDPuppi.userInts.PFJetIDTightLepVeto = cms.InputTag("PFJetIDTightLepVetoPuppi")

process.jetIDSequencePuppi = cms.Sequence(process.PFJetIDTightLepVetoPuppi * process.updatedPatJetsUpdatedJECIDPuppi)


################################################
# BTag Event Weights DeepJet loose WP          #
################################################

BTag_DeepJet_looseWP_cut=cms.double(0.0490)

BTagCalibrationReaderPars = cms.PSet(
    measurementType_bJets = cms.string('mujets'),
    measurementType_cJets = cms.string('mujets'),
    measurementType_lightJets = cms.string('incl'),
)

BTagCalibrationPars = cms.PSet(
    FullSimTagger = cms.string('DeepJet'),
    FullSimFileName = cms.string('data/2018/DeepJet_106XUL18SF_WPonly.csv'),       
)

bTagEffMapPars = cms.PSet(
    fullSimFile = cms.string('data/2018/bTagEff_2018.root'),
    bEffFullSimName_ee = cms.string('ee/B_DeepJet_loose'),
    cEffFullSimName_ee = cms.string('ee/C_DeepJet_loose'),
    lightEffFullSimName_ee = cms.string('ee/Light_DeepJet_loose'),
    bEffFullSimName_mumu = cms.string('mumu/B_DeepJet_loose'),
    cEffFullSimName_mumu = cms.string('mumu/C_DeepJet_loose'),
    lightEffFullSimName_mumu = cms.string('mumu/Light_DeepJet_loose'),
    bEffFullSimName_emu = cms.string('emu/B_DeepJet_loose'),
    cEffFullSimName_emu = cms.string('emu/C_DeepJet_loose'),
    lightEffFullSimName_emu = cms.string('emu/Light_DeepJet_loose'),
)

################################################
# BTag Event Weights DeepCSV loose WP          #
################################################

BTag_DeepCSV_looseWP_cut=cms.double(0.1208)

BTagCalibrationPars_DeepCSV = cms.PSet(
    FullSimTagger = cms.string('DeepCSV'),
    FullSimFileName = cms.string('data/2018/DeepCSV_106XUL18SF_WPonly.csv'),       
)

bTagEffMapPars_DeepCSV = cms.PSet(
    fullSimFile = cms.string('data/2018/bTagEff_2018.root'),
    bEffFullSimName_ee = cms.string('ee/B_DeepCSV_loose'),
    cEffFullSimName_ee = cms.string('ee/C_DeepCSV_loose'),
    lightEffFullSimName_ee = cms.string('ee/Light_DeepCSV_loose'),
    bEffFullSimName_mumu = cms.string('mumu/B_DeepCSV_loose'),
    cEffFullSimName_mumu = cms.string('mumu/C_DeepCSV_loose'),
    lightEffFullSimName_mumu = cms.string('mumu/Light_DeepCSV_loose'),
    bEffFullSimName_emu = cms.string('emu/B_DeepCSV_loose'),
    cEffFullSimName_emu = cms.string('emu/C_DeepCSV_loose'),
    lightEffFullSimName_emu = cms.string('emu/Light_DeepCSV_loose'),
)


################################
# Ttbar Gen Info               #
################################
genParticleCollection = "prunedGenParticles"
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
process.initSubset.src = genParticleCollection
process.decaySubset.src = genParticleCollection
process.decaySubset.fillMode = "kME" # Status3, use kStable for Status2
#  ~process.decaySubset.fillMode = "kStable" # Status3, use kStable for Status2
runMode = "Run2"
process.decaySubset.runMode = runMode

################################
# Ttbar Pseudo Info            #
################################
process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
process.pseudoTop = cms.EDProducer("PseudoTopProducer",
   genParticles = cms.InputTag("prunedGenParticles"),
   finalStates = cms.InputTag("packedGenParticles"),
   minLeptonPt = cms.double(20),
   maxLeptonEta = cms.double(2.4),
   minJetPt = cms.double(30),
   maxJetEta = cms.double(2.4),
   leptonConeSize = cms.double(0.1),
   jetConeSize = cms.double(0.4),
   wMass = cms.double(80.4),
   tMass = cms.double(172.5),

   minLeptonPtDilepton = cms.double(20),
   maxLeptonEtaDilepton = cms.double(2.4),
   minDileptonMassDilepton = cms.double(20),
   minLeptonPtSemilepton = cms.double(20),
   maxLeptonEtaSemilepton = cms.double(2.4),
   minVetoLeptonPtSemilepton = cms.double(15),
   maxVetoLeptonEtaSemilepton = cms.double(2.5),
   minMETSemiLepton = cms.double(30),
   minMtWSemiLepton = cms.double(30),
   clusterNeutrinosFromHadronsIntoJetsDilepton = cms.bool(True)
)
if not isRealData: process.TTbarGen = cms.Sequence(process.makeGenEvt * process.pseudoTop)
else: process.TTbarGen = cms.Sequence()

################################
# Define input and output      #
################################
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), lumisToProcess = cms.untracked.VLuminosityBlockRange("321393:46"))#eventsToProcess = cms.untracked.VEventRange("276315:134:185994346"))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), eventsToProcess = cms.untracked.VEventRange("315257:72:49554411"))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), eventsToProcess = cms.untracked.VEventRange("315259:66:43335866"))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), eventsToProcess = cms.untracked.VEventRange("315259:66:43335866"))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))

################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    jet_pT_cut=cms.untracked.double(15), # for all jets
                                    NumberLeptons_cut=cms.untracked.uint32(2),
                                    # physics objects
                                    #  ~jets = cms.InputTag("slimmedJets"),
                                    #  ~jets = cms.InputTag("patSmearedJets"),
                                    #  ~jets = cms.InputTag("updatedPatJetsUpdatedJECIDsmeared"),
                                    jets = cms.InputTag("updatedPatJetsUpdatedJECID"),    #Use unsmeared jets and apply JER locally
                                    jets_puppi = cms.InputTag("updatedPatJetsUpdatedJECIDPuppi"),
                                    #  ~jets_puppi = cms.InputTag("slimmedJetsPuppi"),
                                    muons = cms.InputTag("MuonsAddedRochesterCorr"),
                                    #  ~muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    photons = cms.InputTag("slimmedPhotons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    metsPuppi = cms.InputTag("slimmedMETsPuppi"),
                                    metsNoHF = cms.InputTag("slimmedMETsNoHF"),
                                    caloMets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    pileUpSummary = cms.InputTag('slimmedAddPileupInfo'),
                                    lheEventProduct = cms.InputTag('externalLHEProducer'),
                                    packedCandidates=cms.InputTag("packedPFCandidates"),
                                    # electron IDs
                                    electronVetoIdMap   = cms.untracked.string("cutBasedElectronID-Fall17-94X-V2-veto"),
                                    electronLooseIdMap  = cms.untracked.string("cutBasedElectronID-Fall17-94X-V2-loose"),
                                    electronMediumIdMap = cms.untracked.string("cutBasedElectronID-Fall17-94X-V2-medium"),
                                    electronTightIdMap  = cms.untracked.string("cutBasedElectronID-Fall17-94X-V2-tight"),
                                    # photon IDs
                                    #  ~photonLooseIdMap   = cms.untracked.string("cutBasedPhotonID-Fall17-94X-V2-loose"),
                                    #  ~photonMediumIdMap  = cms.untracked.string("cutBasedPhotonID-Fall17-94X-V2-medium"),
                                    #  ~photonTightIdMap   = cms.untracked.string("cutBasedPhotonID-Fall17-94X-V2-tight"),
                                    # met filters to apply
                                    metFilterNames=cms.untracked.vstring(
                                        "Flag_goodVertices",
                                        "Flag_globalSuperTightHalo2016Filter",
                                        "Flag_HBHENoiseFilter",
                                        "Flag_HBHENoiseIsoFilter",
                                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                                        "Flag_BadPFMuonFilter",
                                        #  ~"Flag_BadPFMuonDzFilter",       #only available in miniAODv2 (for v1 hack applied)
                                        "Flag_eeBadScFilter",
                                        "Flag_ecalBadCalibFilter",

                                    ),
                                    # choose pileup data
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2018_25ns_UltraLegacy_PoissonOOTPU"),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    triggerObjectNames=cms.vstring(),
                                    pfJetIDSelector=cms.PSet(version=cms.string('RUNIIULCHS'), quality=cms.string('TIGHT')),
                                    triggerPrescales=cms.vstring(), # also useful to check whether a trigger was run
                                    LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars,
                                    BTag_DeepJet_looseWP_cut=BTag_DeepJet_looseWP_cut,
                                    BTag_DeepCSV_looseWP_cut=BTag_DeepCSV_looseWP_cut,
                                    BTagCalibration = BTagCalibrationPars,
                                    BTagCalibration_DeepCSV = BTagCalibrationPars_DeepCSV,
                                    BTagCalibrationReader = BTagCalibrationReaderPars,
                                    bTagEfficiencies = bTagEffMapPars,
                                    bTagEfficiencies_DeepCSV = bTagEffMapPars_DeepCSV,
                                    ttbarGenInfo = cms.bool(False),
                                    ttbarPseudoInfo = cms.bool(False),
                                    DYptInfo = cms.bool(False),
)

################################
# Modify the TreeWriter module #
################################

process.TreeWriter.ttbarGenInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/SMS-T"))
process.TreeWriter.ttbarPseudoInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/SMS-T"))
process.TreeWriter.DYptInfo=(dataset.startswith("/DY"))

if not isRealData:
    if "Fast" in dataset:
        process.TreeWriter.metFilterNames.remove("Flag_globalSuperTightHalo2016Filter")
        process.TreeWriter.lheEventProduct = "source"

# determine user if not set by crab
user=options.user or getpass.getuser()
# user settings
if user=="dmeuser":
    process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter", "hltEG165HE10Filter"]
    process.TreeWriter.triggerNames=[
        # ee Channel
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_DoubleEle25_CaloIdL_MW_v",
        "HLT_Ele32_WPTight_Gsf_v",
        # emu Channel
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        # mumu Chanel
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
        "HLT_IsoMu24_v",
        # MET
        "HLT_PFMET200_HBHECleaned_v",
        "HLT_PFMET200_HBHE_BeamHaloCleaned_v",
        "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",
        "HLT_PFMET120_PFMHT120_IDTight_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v",
        "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v",
        "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v",
    ]
    process.TreeWriter.triggerPrescales=[
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
    process.jecSequence
    *process.jetIDSequence
    #  ~*process.updatedPatJetsUpdatedJECIDsmeared
    *process.egammaPostRecoSeq
    *process.MuonsAddedRochesterCorr
    *process.fullPatMetSequence
    *process.puppiSequence
    *process.jecSequencePuppi
    *process.jetIDSequencePuppi
    *process.TTbarGen
    *process.BadPFMuonFilterUpdateDz
    *process.TreeWriter
)
