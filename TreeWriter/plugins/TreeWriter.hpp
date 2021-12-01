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
#include "FWCore/Framework/interface/Run.h"

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
#include <TreeWriter/TreeWriter/plugins/XYMETCorrection.h>


typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> EcalRecHitCollection;


//
// class declaration
//

class TreeWriter : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks,edm::one::WatchRuns> {
public:
   explicit TreeWriter(const edm::ParameterSet&);
   ~TreeWriter() {};

   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
   virtual void beginJob() override {};
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override {};
   
   virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
   virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

   virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
   virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override {};

   TH1D* createCutFlowHist(std::string modelName = "");
   TH1D* createSystMCWeightHist(const std::string &histName, const int &nBins);

   // ----------member data ---------------------------
   double dJet_pT_cut_;
   unsigned NumberLeptons_cut_;

   bool newLumiBlock_;

   edm::EDGetTokenT<reco::VertexCollection>    vtxToken_;
   edm::EDGetTokenT<pat::JetCollection>        jetCollectionToken_;
   edm::EDGetTokenT<pat::JetCollection>        jetPuppiCollectionToken_;
   edm::EDGetTokenT<reco::GenJetCollection>    genJetCollectionToken_;
   edm::EDGetTokenT<pat::MuonCollection>       muonCollectionToken_;
   edm::EDGetTokenT<edm::View<pat::Electron>>  electronCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metPuppiCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metNoHFCollectionToken_;
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

   // met filters to apply
   const std::vector<std::string> metFilterNames_;

   const std::string pileupHistogramName_;
   const std::string pileupHistogramNameUp_;
   const std::string pileupHistogramNameDown_;

   PFJetIDSelectionFunctor jetIdSelector;

   // === TREE DATA ===
   TTree *eventTree_;

   Int_t   nPV_;   // number of reconsrtucted primary vertices
   Int_t   true_nPV_;   // true number of reconsrtucted primary vertices
   Int_t   nGoodVertices_;
   Int_t   nTracksPV_;
   Float_t rho_;
   Float_t caloMetPt_;

   // gen weights
   Float_t pu_weight_;
   Float_t pu_weight_up_;
   Float_t pu_weight_down_;
   Float_t mc_weight_;
   Float_t topPTweight_;
   std::vector<float> vPdf_weights_;
   std::vector<float> vPS_weights_;
   
   // event information
   ULong64_t evtNo_;
   UInt_t    runNo_;
   UInt_t    lumNo_;

   // Trigger decisions
   std::vector<std::string>      triggerNames_;
   std::vector<std::string>      triggerPrescales_;
   std::map<std::string, Bool_t> triggerDecision_;
   std::map<std::string, int>    triggerPrescale_;
   std::map<std::string, int>    triggerIndex_;
   std::vector<std::string>      triggerObjectNames_;

   // physics Objects
   std::vector<tree::Jet>      vJets_;
   std::vector<tree::Jet>      vJetsPuppi_;
   std::vector<tree::Particle> vGenJets_;
   std::vector<tree::Electron> vElectrons_;
   std::vector<tree::Muon>     vMuons_;
   std::vector<tree::Electron> vElectrons_add_;
   std::vector<tree::Muon>     vMuons_add_;
   tree::MET                   met_;
   tree::MET                   met_UnclEu_;
   tree::MET                   met_UnclEd_;
   tree::MET                   metCalo_;
   tree::MET                   metPuppi_;
   tree::MET                   metPuppi_UnclEu_;
   tree::MET                   metPuppi_UnclEd_;
   tree::MET                   met_raw_;
   tree::MET                   met_gen_;
   tree::MET                   metXYcorr_;
   std::map<std::string,std::vector<tree::Particle>> triggerObjectMap_;
   
   // other variables
   Float_t Ht_;
   Float_t genHt_;
   Float_t puPtHat_;
   Float_t EWKinoPairPt_;
   Float_t genMT2_;
   Float_t genMT2neutrino_;
   
   // gen particles
   std::vector<tree::GenParticle> vGenParticles_;
   std::vector<tree::IntermediateGenParticle> vIntermediateGenParticles_;

   // File Service to store things to a root file
   edm::Service<TFileService> fs_;

   // histogram to store #evts after each "cut"
   TH1D* hCutFlow_;
   
   // histograms to store sum of MC weights
   TH1D* hSystMCweight_PS_norm_;
   TH1D* hSystMCweight_PS_;
   TH1D* hSystMCweight_PDF_norm_;
   TH1D* hSystMCweight_PDF_;
   TH1D* hSystMCweight_topPt_norm_;
   TH1D* hSystMCweight_topPt_;
   TH1D* hSystMCweight_bFrag_norm_;
   TH1D* hSystMCweight_bFrag_;
   TH1D* hSystMCweight_PU_norm_;
   TH1D* hSystMCweight_PU_;

   // Pileup histogram(s)
   TH1F hPU_;
   TH1F hPU_up_;
   TH1F hPU_down_;
   
   // Add lepton veto
   bool addLepton_;
   
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
   
   //DY gen Info
   const bool dyPtInfo_;
   Float_t dyPt_;
   std::vector<tree::GenParticle> v_dyLep_;
   
   //BFragmenation
   const bool bFragInfo_;
   edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsBfragToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > bfragWeight_BLCentralToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > bfragWeight_BLUpToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > bfragWeight_BLDownToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > bfragWeight_PetersonToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > bfragWeight_BSemiLepUpToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > bfragWeight_BSemiLepDownToken_;
   
   float fragUpWeight_;
   float fragCentralWeight_;
   float fragDownWeight_;
   float fragPetersonWeight_;
   float semilepbrUpWeight_;
   float semilepbrDownWeight_;
   
};

#endif /* TREEWRITER_HPP__ */
