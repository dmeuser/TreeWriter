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

# input files for local testing
options.inputFiles =    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/00000/054AE840-3924-BD47-880F-12FA2F909901.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2017B/DoubleMuon/MINIAOD/09Aug2019_UL2017-v1/260000/032A7B27-0F31-D348-9B82-9E1B62ED9587.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17MiniAOD/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/05F2529D-1C66-3745-925D-C6B3C1E850ED.root',

# defaults
options.outputFile = 'ttbarTree.root'
#  ~options.maxEvents = -1
options.maxEvents = 100
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

# the actual TreeWriter module
process = cms.Process("TreeWriter")

# message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# determine global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if isRealData:
        process.GlobalTag.globaltag = "106X_dataRun2_v32"
        #  ~process.GlobalTag.globaltag = "106X_dataRun2_v33"       #miniAODv2 (correct? or only for production?)

else:
        process.GlobalTag.globaltag = "106X_mc2017_realistic_v8"
        #  ~process.GlobalTag.globaltag = "106X_mc2017_realistic_v9"     #miniAODv2 (correct? or only for production?)

        
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

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True, #needed for Puppi correction
                       era='2017-UL') 
                       

#############
# MUONS     #
#############
from TreeWriter.TreeWriter.muonPATUserDataRochesterCorrectionAdder_cfi import muonPATUserDataRochesterCorrectionAdder
process.MuonsAddedRochesterCorr = muonPATUserDataRochesterCorrectionAdder.clone(
   src = "slimmedMuons",
   applyEnergyCorrections = False,
   debug = False,
)
process.MuonsAddedRochesterCorr.path = cms.string('TreeWriter/data/2017/RoccoR2017UL.txt')

##########################
# MET                    #
##########################
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(   #update pfMET 
    process,
    isData=isRealData,
    fixEE2017 = False,
)
# If using patSmearedJETS is removing too many jets due to lepton matching, use the following to keep all jets: https://indico.cern.ch/event/882528/contributions/3718330/attachments/1974802/3286397/Sync_summary.pdf
# Puppi MET is correct when applying new Puppi Tune

#########################################
# MET Filter (not needed in MiniAODv2)  #
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
# Jets                         #
################################
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#JEC
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
)
process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

#TightLeptonVeto ID
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.PFJetIDTightLepVeto = cms.EDProducer('PatJetIDValueMapProducer',
    src = cms.InputTag("updatedPatJetsUpdatedJEC"),
    filterParams = cms.PSet(
        version = cms.string('RUN2ULCHS'),
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

#JER (currently not used, but applied in local framework)
#from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/test/runJERsmearingOnMiniAOD.py
process.updatedPatJetsUpdatedJECIDsmeared = cms.EDProducer('SmearedPATJetProducer',
       src = cms.InputTag('updatedPatJetsUpdatedJECID'),
       enabled = cms.bool(True),
       rho = cms.InputTag("fixedGridRhoFastjetAll"),
       algo = cms.string('AK4PFchs'),
       algopt = cms.string('AK4PFchs_pt'),

       genJets = cms.InputTag('slimmedGenJets'),
       dRMax = cms.double(0.2),
       dPtMaxFactor = cms.double(3),

       debug = cms.untracked.bool(False),
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
# Update to new Puppi tune
UpdatePuppiTuneV15(process, runOnMC=(isRealData==False))

#JEC
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsPuppi'),
   labelName = 'UpdatedJECPuppi',
   jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
)
process.jecSequencePuppi = cms.Sequence(process.patJetCorrFactorsUpdatedJECPuppi * process.updatedPatJetsUpdatedJECPuppi)

#TightLeptonVeto ID
process.PFJetIDTightLepVetoPuppi = process.PFJetIDTightLepVeto.clone()
process.PFJetIDTightLepVetoPuppi.src = cms.InputTag("updatedPatJetsUpdatedJECPuppi")
process.PFJetIDTightLepVetoPuppi.filterParams.version = cms.string("RUN2ULPUPPI")

process.updatedPatJetsUpdatedJECIDPuppi = process.updatedPatJetsUpdatedJECID.clone()
process.updatedPatJetsUpdatedJECIDPuppi.src = cms.InputTag("updatedPatJetsUpdatedJECPuppi")
process.updatedPatJetsUpdatedJECIDPuppi.userInts.PFJetIDTightLepVeto = cms.InputTag("PFJetIDTightLepVetoPuppi")

process.jetIDSequencePuppi = cms.Sequence(process.PFJetIDTightLepVetoPuppi * process.updatedPatJetsUpdatedJECIDPuppi)

################################
# Ttbar Gen Info               #
################################
genParticleCollection = "prunedGenParticles"
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
process.initSubset.src = genParticleCollection
process.decaySubset.src = genParticleCollection
process.decaySubset.fillMode = "kME" # Status3, use kStable for Status2
runMode = "Run2"
process.decaySubset.runMode = runMode

####################################
# Ttbar Pseudo Info(particle level)#
####################################
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
# BFragmenation Weights        #
################################
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                    inputPruned = cms.InputTag("prunedGenParticles"),
                    inputPacked = cms.InputTag("packedGenParticles"),
)
from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
process.genParticles2HepMC = genParticles2HepMC.clone(genParticles = cms.InputTag("mergedGenParticles"))
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
process.particleLevel.excludeNeutrinosFromJetClustering = False
process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')
if not isRealData: process.bFragWgtSequence = cms.Sequence(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer)
else: process.bFragWgtSequence = cms.Sequence()

################################
# Define input and output      #
################################
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), lumisToProcess = cms.untracked.VLuminosityBlockRange("273302:322"))#eventsToProcess = cms.untracked.VEventRange("276315:134:185994346"))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), eventsToProcess = cms.untracked.VEventRange("279694:2133:3845835844"))
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
                                    jets = cms.InputTag("updatedPatJetsUpdatedJECID"),    #Use unsmeared jets and apply JER locally
                                    jets_puppi = cms.InputTag("updatedPatJetsUpdatedJECIDPuppi"),
                                    muons = cms.InputTag("MuonsAddedRochesterCorr"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
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
                                        "Flag_ecalBadCalibFilter"

                                    ),
                                    # choose pileup data
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2017_25ns_UltraLegacy_PoissonOOTPU"),
                                    pileupHistogramNameUp=cms.untracked.string("pileupWeightUp_mix_2017_25ns_UltraLegacy_PoissonOOTPU"),
                                    pileupHistogramNameDown=cms.untracked.string("pileupWeightDown_mix_2017_25ns_UltraLegacy_PoissonOOTPU"),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    triggerObjectNames=cms.vstring(),
                                    pfJetIDSelector=cms.PSet(version=cms.string('RUN2ULCHS'), quality=cms.string('TIGHT')),
                                    triggerPrescales=cms.vstring(), # also useful to check whether a trigger was run
                                    ttbarGenInfo = cms.bool(False),
                                    ttbarPseudoInfo = cms.bool(False),
                                    DYptInfo = cms.bool(False),
                                    bFragInfo = cms.bool(False),
)

################################
# Modify the TreeWriter module #
################################

# set booleans depending on sample name
process.TreeWriter.ttbarGenInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/SMS-T"))
process.TreeWriter.ttbarPseudoInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/SMS-T"))
process.TreeWriter.bFragInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/ST"))
process.TreeWriter.DYptInfo=(dataset.startswith("/DY"))

# set triggers
process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter", "hltEG165HE10Filter"]
process.TreeWriter.triggerNames=[
    # ee Channel
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_DoubleEle25_CaloIdL_MW_v",
    "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",
    # emu Channel
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    # mumu Chanel
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
    "HLT_IsoMu27_v",
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

for trig in process.TreeWriter.triggerPrescales:
    assert(trig in process.TreeWriter.triggerNames),"Trigger '"+trig+"' is not used, so prescale cannot be stored!"

# determine user if not set by crab
user=options.user or getpass.getuser()
####################
#     RUN          #
####################

process.p = cms.Path(
    process.jecSequence
    *process.jetIDSequence
    *process.egammaPostRecoSeq
    *process.MuonsAddedRochesterCorr
    *process.fullPatMetSequence
    *process.puppiSequence
    *process.jecSequencePuppi
    *process.jetIDSequencePuppi
    *process.TTbarGen
    *process.BadPFMuonFilterUpdateDz
    *process.bFragWgtSequence
    *process.TreeWriter
)
