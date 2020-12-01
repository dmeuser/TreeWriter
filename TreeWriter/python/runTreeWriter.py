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
#  ~options.inputFiles =    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/80EBB916-A6C4-E811-9A55-A4BF011254E0.root'
#  ~options.inputFiles =    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/100000/AA72240A-9EC2-E811-9D0A-008CFA1C6414.root'
#  ~options.inputFiles =    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/FEF559B1-95DF-E811-9528-D4AE52900EF9.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/A6781AB0-A774-E911-8648-0425C5903030.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/90000/D2E079F1-29F5-E811-B34D-0CC47AF9B2CA.root'
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleMuon/MINIAOD/17Jul2018_ver2-v1/40000/A4DF1FC2-E28B-E811-8D40-509A4C838D01.root'
options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016H/DoubleEG/MINIAOD/17Jul2018-v1/80000/A896F673-5C8D-E811-AE20-0242AC1C0502.root'
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016H/MET/MINIAOD/17Jul2018-v2/270000/FC265DF1-71B9-E811-AF11-34E6D7BEAF0E.root'
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/ZZTo2L2Nu_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/40000/00D06C25-A823-E911-AE44-001E67792800.root',
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/60000/F40DF1A0-E937-E911-AB5B-AC1F6BAC7D12.root'
#  ~options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/10000/009FB1D8-0117-E911-B425-A4BF01125E66.root'
options.outputFile = 'ttbarTree.root'
#~ options.outputFile = 'overlap_lepton_2.root'
options.maxEvents = -1
#  ~options.maxEvents = 100
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

# the actual TreeWriter module
process = cms.Process("TreeWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# determine global tag here only 2016
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if isRealData:
        process.GlobalTag.globaltag = "94X_dataRun2_v10"
else:
        process.GlobalTag.globaltag = "94X_mcRun2_asymptotic_v3"
#  ~if isRealData:
        #  ~process.GlobalTag.globaltag = "102X_dataRun2_v13"
#  ~else:
        #  ~process.GlobalTag.globaltag = "102X_mcRun2_asymptotic_v8"
        
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
# https://hypernews.cern.ch/HyperNews/CMS/get/egamma/2204/1/1.html because of PUPPI MET: added phoIDModules=[]
# Geometry neccessary to run setupEgammaPostRecoSeq
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                       era='2016-Legacy', phoIDModules=[])
                       

#############
# MUONS     #
#############
from TreeWriter.TreeWriter.muonPATUserDataRochesterCorrectionAdder_cfi import muonPATUserDataRochesterCorrectionAdder
process.MuonsAddedRochesterCorr = muonPATUserDataRochesterCorrectionAdder.clone(
   src = "slimmedMuons",
   #  ~applyEnergyCorrections = True,
   applyEnergyCorrections = False,
   debug = False,
)
#  ~process.MuonsAddedRochesterCorr.path = '${CMSSW_BASE}/src/TreeWriter/data/2016/RoccoR2016.txt'
process.MuonsAddedRochesterCorr.path = cms.string('TreeWriter/data/2016/RoccoR2016.txt')

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

runMetCorAndUncFromMiniAOD(   #new MET collection with BJet regression corrected jets
    process,
    isData=isRealData,
    fixEE2017 = False,
    jetCollUnskimmed = "updatedJetsBJetRegression",
    postfix = "BJetRegressionCorr",
    reapplyJEC = False,    #updatedJetsBJetRegression have already latest JEC included
)

runMetCorAndUncFromMiniAOD(   #new MET collection with BJet regression corrected jets (loose selection for jets)
    process,
    isData=isRealData,
    fixEE2017 = False,
    jetCollUnskimmed = "updatedJetsBJetRegressionLoose",
    postfix = "BJetRegressionCorrLoose",
    reapplyJEC = False,    #updatedJetsBJetRegression have already latest JEC included
)

from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True);
runMetCorAndUncFromMiniAOD(process,    #update PuppiMET 
                           isData=isRealData,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           )
process.puppiNoLep.useExistingWeights = True
process.puppi.useExistingWeights = True
process.puppiMETSequence = cms.Sequence(process.egmPhotonIDSequence * process.puppiMETSequence * process.fullPatMetSequencePuppi)

#Add DeepMET
from RecoMET.METPUSubtraction.deepMETProducer_cfi import deepMETProducer
process.deepMETProducer = deepMETProducer.clone()
process.deepMETProducer.graph_path="RecoMET/METPUSubtraction/data/tf_models/deepmet_2016.pb"


################################
# MET Filter                   #
################################
#~ process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#~ process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#~ process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

################################
# Lepton Scale Factors         #
################################
LeptonFullSimScaleFactorMapPars2016 = cms.PSet(
    #~ dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2016/RunBCDEF_SF_ID.root'),
    #~ dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2016/RunBCDEF_SF_ISO.root'),
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2016/MuonSF_ID_merged.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2016/MuonSF_ISO_merged.root'),
    
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2016/2016LegacyReReco_ElectronTight_Fall17V2.root'),
    
    dataMCScaleFactorHisto_mu_ID = cms.string('NUM_TightID_DEN_genTracks_eta_pt'),
    dataMCScaleFactorHisto_mu_Iso = cms.string('NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt'),
    
    dataMCScaleFactorHisto_ele_ID = cms.string('EGamma_SF2D'),
)

################################
# Jets                         #
################################
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
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
        version = cms.string('WINTER16'),
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

################################
# Puppi Jets                   #
################################
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


################################
# BTag Event Weights           #
################################
BTagCalibrationReaderPars2016 = cms.PSet(
    measurementType_bJets = cms.string('mujets'),
    measurementType_cJets = cms.string('mujets'),
    measurementType_lightJets = cms.string('incl'),
)

BTagCalibrationPars2016 = cms.PSet(
    CSVFullSimTagger = cms.string('DeepCSV'),
    CSVFullSimFileName = cms.string('data/2016/DeepCSV_2016LegacySF_WP_V1.csv'),       
)

bTagEffMapPars2016 = cms.PSet(
    fullSimFile = cms.string('${CMSSW_BASE}/src/TreeWriter/data/2016/btageff__ttbar_powheg_pythia8_25ns_Moriond17_deepCSV.root'),
    bEffFullSimName = cms.string('h2_BTaggingEff_csv_loose_Eff_b'),
    cEffFullSimName = cms.string('h2_BTaggingEff_csv_loose_Eff_c'),
    lightEffFullSimName = cms.string('h2_BTaggingEff_csv_loose_Eff_udsg'),
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
# Bjet Regression              #
################################
# https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/jets_cff.py
# https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/BJetRegression

collection="updatedPatJetsUpdatedJEC"

process.bJetVars = cms.EDProducer("JetRegressionVarProducer",
   pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
   src = cms.InputTag(collection),
   svsrc = cms.InputTag("slimmedSecondaryVertices"),
   gpsrc = cms.InputTag("prunedGenParticles"),
)


#  ~if opts['sample']['era'] == '2016':
wgtFile = cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2016.pb")
#  ~wgtFile = cms.FileInPath('TreeWriter/data/2016/breg_training_2016_JECv11_Oct_2019.pb')
outputFormulas = cms.vstring(["at(0)*0.31976690888404846+1.047176718711853","0.5*(at(2)-at(1))*0.31976690888404846"])
wgtFile_c = cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2016.pb")
outputFormulas_c = cms.vstring(["at(0)*0.28862622380256653+0.9908722639083862","0.5*(at(2)-at(1))*0.28862622380256653"])
#  ~elif opts['sample']['era'] == '2017':
   #  ~wgtFile = cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2017.pb")
   #  ~outputFormulas = cms.vstring(["at(0)*0.28225210309028625+1.055067777633667","0.5*(at(2)-at(1))*0.28225210309028625"])
   #  ~wgtFile_c = cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2017.pb")
   #  ~outputFormulas_c = cms.vstring(["at(0)*0.24718524515628815+0.9927206635475159","0.5*(at(2)-at(1))*0.24718524515628815"])
#  ~elif opts['sample']['era'] == '2018':
   #  ~wgtFile = cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2018.pb")
   #  ~outputFormulas = cms.vstring(["at(0)*0.27912887930870056+1.0545977354049683","0.5*(at(2)-at(1))*0.27912887930870056"])
   #  ~wgtFile_c = cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2018.pb")
   #  ~outputFormulas_c = cms.vstring(["at(0)*0.24325256049633026+0.993854820728302","0.5*(at(2)-at(1))*0.24325256049633026"])

process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag(collection),
    userFloats = cms.PSet(
        leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
        leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
        leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
        leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
        leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
        leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
        leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
        leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
        leptonPt = cms.InputTag("bJetVars:leptonPt"),
        vtxPt = cms.InputTag("bJetVars:vtxPt"),
        vtxMass = cms.InputTag("bJetVars:vtxMass"),
        vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
        vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
        ptD = cms.InputTag("bJetVars:ptD"),
        genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
        ),
    userInts = cms.PSet(
       vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
       leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
    ),
)
#  ~if opts['sample']['era'] == '2016':
    #  ~process.updatedJetsWithUserData.userInts = cms.PSet(
       #  ~vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
       #  ~leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
    #  ~)
collection = 'updatedJetsWithUserData'

process.bjetNN = cms.EDProducer("BJetEnergyRegressionMVA",
   backend = cms.string("TF"),
   src = cms.InputTag(collection),
   pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
   svsrc = cms.InputTag("slimmedSecondaryVertices"),
   rhosrc = cms.InputTag("fixedGridRhoFastjetAll"),

   weightFile =  wgtFile,
   name = cms.string("JetRegNN"),
   isClassifier = cms.bool(False),
   variablesOrder = cms.vstring(["Jet_pt","Jet_eta","rho","Jet_mt","Jet_leadTrackPt","Jet_leptonPtRel","Jet_leptonDeltaR","Jet_neHEF",
                                 "Jet_neEmEF","Jet_vtxPt","Jet_vtxMass","Jet_vtx3dL","Jet_vtxNtrk","Jet_vtx3deL",
                                 "Jet_numDaughters_pt03","Jet_energyRing_dR0_em_Jet_rawEnergy","Jet_energyRing_dR1_em_Jet_rawEnergy",
                                 "Jet_energyRing_dR2_em_Jet_rawEnergy","Jet_energyRing_dR3_em_Jet_rawEnergy","Jet_energyRing_dR4_em_Jet_rawEnergy",
                                 "Jet_energyRing_dR0_neut_Jet_rawEnergy","Jet_energyRing_dR1_neut_Jet_rawEnergy","Jet_energyRing_dR2_neut_Jet_rawEnergy",
                                 "Jet_energyRing_dR3_neut_Jet_rawEnergy","Jet_energyRing_dR4_neut_Jet_rawEnergy","Jet_energyRing_dR0_ch_Jet_rawEnergy",
                                 "Jet_energyRing_dR1_ch_Jet_rawEnergy","Jet_energyRing_dR2_ch_Jet_rawEnergy","Jet_energyRing_dR3_ch_Jet_rawEnergy",
                                 "Jet_energyRing_dR4_ch_Jet_rawEnergy","Jet_energyRing_dR0_mu_Jet_rawEnergy","Jet_energyRing_dR1_mu_Jet_rawEnergy",
                                 "Jet_energyRing_dR2_mu_Jet_rawEnergy","Jet_energyRing_dR3_mu_Jet_rawEnergy","Jet_energyRing_dR4_mu_Jet_rawEnergy",
                                 "Jet_chHEF","Jet_chEmEF","Jet_leptonPtRelInv","isEle","isMu","isOther","Jet_mass","Jet_ptd"]),
   variables = cms.PSet(
   Jet_pt = cms.string("pt*jecFactor('Uncorrected')"),
   Jet_mt = cms.string("mt*jecFactor('Uncorrected')"),
   Jet_eta = cms.string("eta"),
   Jet_mass = cms.string("mass*jecFactor('Uncorrected')"),
   Jet_ptd = cms.string("userFloat('ptD')"),
   Jet_leadTrackPt = cms.string("userFloat('leadTrackPt')"),
   Jet_vtxNtrk = cms.string("userInt('vtxNtrk')"),
   Jet_vtxMass = cms.string("userFloat('vtxMass')"),
   Jet_vtx3dL = cms.string("userFloat('vtx3dL')"),
   Jet_vtx3deL = cms.string("userFloat('vtx3deL')"),
   Jet_vtxPt = cms.string("userFloat('vtxPt')"),
   Jet_leptonPtRel = cms.string("userFloat('leptonPtRelv0')"),
   Jet_leptonPtRelInv = cms.string("userFloat('leptonPtRelInvv0')*jecFactor('Uncorrected')"),
   Jet_leptonDeltaR = cms.string("userFloat('leptonDeltaR')"),
   Jet_neHEF = cms.string("neutralHadronEnergyFraction()"),
   Jet_neEmEF = cms.string("neutralEmEnergyFraction()"),
   Jet_chHEF = cms.string("chargedHadronEnergyFraction()"),
   Jet_chEmEF = cms.string("chargedEmEnergyFraction()"),
   isMu = cms.string("?abs(userInt('leptonPdgId'))==13?1:0"),
   isEle = cms.string("?abs(userInt('leptonPdgId'))==11?1:0"),
   isOther = cms.string("?userInt('leptonPdgId')==0?1:0"),
   ),
    inputTensorName = cms.string("ffwd_inp"),
    outputTensorName = cms.string("ffwd_out/BiasAdd"),
    outputNames = cms.vstring(["corr","res"]),
    outputFormulas = outputFormulas,
    nThreads = cms.uint32(1),
    singleThreadPool = cms.string("no_threads"),
)

process.UpdatedPatJetsAK04CHSBCRegressionUserData = cms.EDProducer('PATJetUserDataEmbedder',
   src = cms.InputTag(collection),
   userFloats = cms.PSet(
      BJetEnergyCorrFactor       = cms.InputTag('bjetNN:corr'),
      BJetEnergyCorrResolution       = cms.InputTag('bjetNN:res'),
   ),
)

from TreeWriter.TreeWriter.JetPATBJetRegressionCorrectionAdder_cfi import JetPATBJetRegressionCorrectionAdder
process.updatedJetsBJetRegression = JetPATBJetRegressionCorrectionAdder.clone(   #New Jet collection with b jet regression corrected Jets
   src = "UpdatedPatJetsAK04CHSBCRegressionUserData",
   debug = False,
   applyLooseSelection = False,
)

process.updatedJetsBJetRegressionLoose = JetPATBJetRegressionCorrectionAdder.clone(   #New Jet collection with b jet regression corrected Jets (if pt>15 and eta<2.5)
   src = "UpdatedPatJetsAK04CHSBCRegressionUserData",
   debug = False,
   applyLooseSelection = True,
)


################################
# Define input and output      #
################################
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), lumisToProcess = cms.untracked.VLuminosityBlockRange("273302:322"))#eventsToProcess = cms.untracked.VEventRange("276315:134:185994346"))
#  ~process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles), eventsToProcess = cms.untracked.VEventRange("1:260919:52183687"))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))


################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(0),
                                    jet_pT_cut=cms.untracked.double(30), # for all jets
                                    isolatedPhotons=cms.untracked.bool(True), # for all photons in the collection
                                    minNumberElectrons_cut=cms.untracked.uint32(0),
                                    NumberLeptons_cut=cms.untracked.uint32(2),
                                    # physics objects
                                    #  ~jets = cms.InputTag("slimmedJets"),
                                    #  ~jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
                                    jets = cms.InputTag("updatedPatJetsUpdatedJECID"),
                                    #  ~jets = cms.InputTag("UpdatedPatJetsAK04CHSBCRegressionUserData"),
                                    jets_puppi = cms.InputTag("updatedPatJetsUpdatedJECIDPuppi"),
                                    muons = cms.InputTag("MuonsAddedRochesterCorr"),
                                    #  ~muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    photons = cms.InputTag("slimmedPhotons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    metsPuppi = cms.InputTag("slimmedMETsPuppi"),
                                    metsNoHF = cms.InputTag("slimmedMETsNoHF"),
                                    metsDeepMET = cms.InputTag("deepMETProducer"),
                                    metCorr = cms.InputTag(""),
                                    metCorrCal = cms.InputTag(""),
                                    caloMets = cms.InputTag("slimmedMETs"),
                                    mets_BJetRegression = cms.InputTag("slimmedMETsBJetRegressionCorr"),
                                    mets_BJetRegressionLoose = cms.InputTag("slimmedMETsBJetRegressionCorrLoose"),
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
                                    metCalibrated = cms.InputTag("slimmedMETs"),
                                    LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars2016,
                                    BTagCalibration = BTagCalibrationPars2016,
                                    BTagCalibrationReader = BTagCalibrationReaderPars2016,
                                    bTagEfficiencies = bTagEffMapPars2016,
                                    ttbarGenInfo = cms.bool(False),
                                    ttbarPseudoInfo = cms.bool(False),
                                    DYptInfo = cms.bool(False),
)

################################
# Modify the TreeWriter module #
################################

process.TreeWriter.hardPUveto=dataset.startswith("/QCD_HT100to200")
process.TreeWriter.ttbarGenInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/SMS-T"))
process.TreeWriter.ttbarPseudoInfo=(dataset.startswith("/TT") or dataset.startswith("/tt") or dataset.startswith("/SMS-T"))
process.TreeWriter.DYptInfo=(dataset.startswith("/DY"))

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
        process.TreeWriter.metFilterNames.remove("Flag_globalSuperTightHalo2016Filter")
        process.TreeWriter.lheEventProduct = "source"

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
        # MET 2016 (top recomm.)
        "HLT_PFHT300_PFMET110_v",
        "HLT_PFMET120_PFMHT120_IDTight_v",
        "HLT_PFMET170_HBHECleaned_v",
        "HLT_PFMET300_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
        "HLT_MET200_v",
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
    process.jecSequence
    *process.jetIDSequence
    *process.jecSequencePuppi
    *process.jetIDSequencePuppi
    #  ~*process.bJetVars
    #  ~*process.updatedJetsWithUserData
    #  ~*process.bjetNN
    #  ~*process.UpdatedPatJetsAK04CHSBCRegressionUserData
    #  ~*process.updatedJetsBJetRegression
    #  ~*process.updatedJetsBJetRegressionLoose
    *process.egammaPostRecoSeq
    *process.MuonsAddedRochesterCorr
    *process.fullPatMetSequence
    #  ~*process.fullPatMetSequenceBJetRegressionCorr
    #  ~*process.fullPatMetSequenceBJetRegressionCorrLoose
    *process.puppiMETSequence
    *process.TTbarGen
    # ~*process.generatorTopFilter
    #  ~*process.deepMETProducer
    *process.TreeWriter
)
