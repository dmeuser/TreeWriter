#ifndef TREEWRITER_HPP__
#define TREEWRITER_HPP__

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <regex>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/Common/interface/EDCollection.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "AnalysisDataFormats/TopObjects/interface/TopGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TMath.h"

#include "TreeParticles.hpp"

#include "MT2Functor.h"
#include <TreeWriter/TreeWriter/plugins/LeptonFullSimScaleFactorMapFunctor.h>


typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> EcalRecHitCollection;


//
// class declaration
//

class TreeWriter : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks> {
public:
   explicit TreeWriter(const edm::ParameterSet&);
   ~TreeWriter() {};

   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
   virtual void beginJob() override {};
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override {};

   virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
   virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override {};

   TH1F* createCutFlowHist(std::string modelName = "");

   // ----------member data ---------------------------
   double dHT_cut_;
   double dJet_pT_cut_;
   unsigned minNumberElectrons_cut_;
   unsigned NumberLeptons_cut_;

   bool newLumiBlock_;

   edm::EDGetTokenT<reco::VertexCollection>    vtxToken_;
   edm::EDGetTokenT<pat::JetCollection>        jetCollectionToken_;
   edm::EDGetTokenT<reco::GenJetCollection>    genJetCollectionToken_;
   edm::EDGetTokenT<pat::MuonCollection>       muonCollectionToken_;
   edm::EDGetTokenT<edm::View<pat::Electron>>  electronCollectionToken_;
   edm::EDGetTokenT<edm::View<pat::Photon>>    photonCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metPuppiCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metNoHFCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCorrectedCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCalibratedCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metDeepCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metBJetRegressionCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metBJetRegressionLooseCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        caloMetCollectionToken_;
   edm::EDGetTokenT<double>                    rhoToken_;
   edm::EDGetTokenT<EcalRecHitCollection>      ebRecHitsToken_;
   edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedGenToken_;
   edm::EDGetTokenT<PileupSummaryInfoCollection>  pileUpSummaryToken_;
   edm::EDGetTokenT<LHEEventProduct>           LHEEventToken_;

   edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidateToken_;

   // electron id
   std::string electronVetoIdMapToken_;
   std::string electronLooseIdMapToken_;
   std::string electronMediumIdMapToken_;
   std::string electronTightIdMapToken_;
   
   // photon id
   std::string photonLooseIdMapToken_;
   std::string photonMediumIdMapToken_;
   std::string photonTightIdMapToken_;

   // met filters to apply
   const std::vector<std::string> metFilterNames_;

   // from photon ID value map producer
   edm::EDGetTokenT<edm::ValueMap<float>> phoWorstChargedIsolationToken_;

   const std::string pileupHistogramName_;
   const bool hardPUveto_;
   const bool reMiniAOD_; // 03Feb2017 campaign

   PFJetIDSelectionFunctor jetIdSelector;

   // === TREE DATA ===
   TTree *eventTree_;

   Int_t   nPV_;   // number of reconsrtucted primary vertices
   Int_t   true_nPV_;   // true number of reconsrtucted primary vertices
   Int_t   nGoodVertices_;
   Int_t   nTracksPV_;
   Float_t rho_;
   Float_t caloMetPt_;

   Float_t pu_weight_;
   Char_t  mc_weight_; // +1 or -1 event weights (take care when reading with python, this is a character!)
   std::vector<float> vPdf_weights_;
   
   Float_t mll_;
   Float_t Ht_;
   Float_t genHt_;
   Float_t puPtHat_;
   Float_t EWKinoPairPt_;
   Float_t MT2_;
   Float_t genMT2_;
   Float_t genMT2neutrino_;
   Float_t dyPt_;


   ULong64_t evtNo_;
   UInt_t    runNo_;
   UInt_t    lumNo_;

   Bool_t particleFlowEGammaGSFixed_dupECALClusters_;
   Bool_t ecalMultiAndGSGlobalRecHitEB_hitsNotReplaced_;

   UShort_t signal_m1_; // usually mass of first particle in decay chain
   UShort_t signal_m2_; // usually neutarlino mass
   UShort_t signal_nBinos_; // 2 for T5gg, 1 for T5Wg, 0 for T5WW
   UShort_t signal_nNeutralinoDecays_;

   // Trigger decisions
   std::vector<std::string>      triggerNames_;
   std::vector<std::string>      triggerPrescales_;
   std::map<std::string, Bool_t> triggerDecision_;
   std::map<std::string, int>    triggerPrescale_;
   std::map<std::string, int>    triggerIndex_;
   std::vector<std::string>      triggerObjectNames_;

   // physics Objects
   std::vector<tree::Jet>      vJets_;
   std::vector<tree::Particle> vGenJets_;
   std::vector<tree::Electron> vElectrons_;
   std::vector<tree::Muon>     vMuons_;
   std::vector<tree::Electron> vElectrons_add_;
   std::vector<tree::Muon>     vMuons_add_;
   std::vector<tree::Photon>   vPhotons_;
   tree::MET                   met_;
   tree::MET                   metCalo_;
   tree::MET                   metPuppi_;
   tree::MET                   metNoHF_;
   tree::MET                   metCorrected_;
   tree::MET                   metCalibrated_;
   tree::MET                   met_raw_;
   tree::MET                   met_gen_;
   tree::MET                   met_JESu_;
   tree::MET                   met_JESd_;
   tree::MET                   met_JERu_;
   tree::MET                   met_JERd_;
   tree::MET                   metDeep_;
   tree::MET                   metBJetRegression_;
   tree::MET                   metBJetRegressionLoose_;
   std::map<std::string,std::vector<tree::Particle>> triggerObjectMap_;

   std::vector<tree::GenParticle> vGenParticles_;
   std::vector<tree::IntermediateGenParticle> vIntermediateGenParticles_;
   
   // TrackIsolation
   Bool_t electronTrackIsoVeto;
   Bool_t muonTrackIsoVeto;
   Bool_t pionTrackIsoVeto;
   
   // Scale factors
   LeptonFullSimScaleFactorMapFunctor fctLeptonFullSimScaleFactors_;
   Float_t lepton1SF_;
   Float_t lepton2SF_;
   Float_t lepton1SF_unc_;
   Float_t lepton2SF_unc_;

   // File Service to store things to a root file
   edm::Service<TFileService> fs_;

   // histogram to store #evts after each "cut"
   TH1F* hCutFlow_;
   std::map<std::string,TH1F*> hCutFlowMap_;

   // Pileup histogram(s)
   TH1F hPU_;
   
   // Dilepton Event Type (ee,mumu,emu)
   bool mumu_;
   bool ee_;
   bool emu_;
   
   //MT2 functor
   MT2Functor fctMT2_;
   double pa[3];
   double pb[3];
   double pmiss[3];
   
   //Ttbar gen Event
   const bool ttbarGenInfo_;
   int ttbarProductionMode_;
   int ttbarDecayMode_;
   TLorentzVector genTop_;
   TLorentzVector genAntiTop_;
   TLorentzVector genLepton_;
   TLorentzVector genAntiLepton_;
   TLorentzVector genTau_;
   TLorentzVector genAntiTau_;
   int genLeptonPdgId_;
   int genAntiLeptonPdgId_;
   TLorentzVector genB_;
   TLorentzVector genAntiB_;
   TLorentzVector genNeutrino_;
   TLorentzVector genAntiNeutrino_;
   TLorentzVector genWMinus_;
   TLorentzVector genWPlus_;
   TLorentzVector genNeutrinoSum_;
   
   //Ttbar pseudo gen Event
   const bool pseudoTopInfo_;
   int ttbarPseudoDecayMode_;
   TLorentzVector pseudoTop_;
   TLorentzVector pseudoAntiTop_;
   TLorentzVector pseudoLepton_;
   TLorentzVector pseudoAntiLepton_;
   TLorentzVector pseudoTau_;
   TLorentzVector pseudoAntiTau_;
   int pseudoLeptonPdgId_;
   int pseudoAntiLeptonPdgId_;
   TLorentzVector pseudoBJet_;
   TLorentzVector pseudoAntiBJet_;
   TLorentzVector pseudoNeutrino_;
   TLorentzVector pseudoAntiNeutrino_;
   TLorentzVector pseudoWMinus_;
   TLorentzVector pseudoWPlus_;
   TLorentzVector pseudoNeutrinoSum_;
   //~std::vector<tree::Particle> v_allPseudoJet_;
   //~std::vector<tree::Particle> v_allPseudoLepton_;
   
   //DY pt Info
   const bool dyPtInfo_;
};

#endif /* TREEWRITER_HPP__ */
