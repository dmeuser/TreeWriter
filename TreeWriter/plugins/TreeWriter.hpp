#ifndef TREEWRITER_HPP__
#define TREEWRITER_HPP__

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <regex>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

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


typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> EcalRecHitCollection;


//
// class declaration
//

class TreeWriter : public edm::EDAnalyzer {
public:
   explicit TreeWriter(const edm::ParameterSet&);
   ~TreeWriter() {};

   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
   virtual void beginJob() override {};
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override {};

   virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
   virtual void endRun(edm::Run const&, edm::EventSetup const&) override {};
   virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
   //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

   int matchToTruth(const pat::Photon &pho,
                    const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
   int matchToTruthAlternative(const pat::Photon &pho,
                               const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

   void findFirstNonPhotonMother(const reco::Candidate *particle,
                                 int &ancestorPID, int &ancestorStatus);

   TH1F* createCutFlowHist(std::string modelName = "");

   // ----------member data ---------------------------
   double dHT_cut_;
   double dPhoton_pT_cut_;
   double dJet_pT_cut_;
   bool isolatedPhotons_;
   unsigned minNumberPhotons_cut_;
   unsigned minNumberElectrons_cut_;
   unsigned minNumberBinos_cut_;

   bool newLumiBlock_;

   edm::EDGetTokenT<reco::VertexCollection>    vtxToken_;
   edm::EDGetTokenT<edm::View<pat::Photon>>    photonCollectionToken_;
   edm::EDGetTokenT<pat::JetCollection>        jetCollectionToken_;
   edm::EDGetTokenT<reco::GenJetCollection>    genJetCollectionToken_;
   edm::EDGetTokenT<pat::MuonCollection>       muonCollectionToken_;
   edm::EDGetTokenT<edm::View<pat::Electron>>  electronCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCorrectedCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        metCalibratedCollectionToken_;
   edm::EDGetTokenT<pat::METCollection>        caloMetCollectionToken_;
   edm::EDGetTokenT<double>                    rhoToken_;
   edm::EDGetTokenT<EcalRecHitCollection>      ebRecHitsToken_;
   edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedGenToken_;
   edm::EDGetTokenT<PileupSummaryInfoCollection>  pileUpSummaryToken_;
   edm::EDGetTokenT<LHEEventProduct>           LHEEventToken_;

   edm::EDGetTokenT<std::vector<pat::PackedCandidate>> packedCandidateToken_;

   // electron id
   edm::EDGetTokenT<edm::ValueMap<bool>> electronVetoIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> electronLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> electronMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> electronTightIdMapToken_;
   // photon id
   edm::EDGetTokenT<edm::ValueMap<bool>> photonLooseId15MapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> photonMediumId15MapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> photonTightId15MapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> photonLooseIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> photonMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool>> photonTightIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<float>> photonMvaValuesMapToken_;
   edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> phoLooseIdFullInfoMapToken_;

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

   Float_t genHt_;
   Float_t puPtHat_;
   Float_t EWKinoPairPt_;

   Int_t nISR_;

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

   // met filter tokens
   edm::EDGetTokenT<bool> BadChCandFilterToken_;
   edm::EDGetTokenT<bool> BadPFMuonFilterToken_;

   // physics Objects
   std::vector<tree::Photon>   vPhotons_;
   std::vector<tree::Jet>      vJets_;
   std::vector<tree::Particle> vGenJets_;
   std::vector<tree::Electron> vElectrons_;
   std::vector<tree::Muon>     vMuons_;
   tree::MET                   met_;
   tree::MET                   metCorrected_;
   tree::MET                   metCalibrated_;
   tree::MET                   met_raw_;
   tree::MET                   met_gen_;
   tree::MET                   met_JESu_;
   tree::MET                   met_JESd_;
   tree::MET                   met_JERu_;
   tree::MET                   met_JERd_;
   std::vector<tree::Particle> vGenPhotons_;
   std::map<std::string,std::vector<tree::Particle>> triggerObjectMap_;

   std::vector<tree::GenParticle> vGenParticles_;
   std::vector<tree::IntermediateGenParticle> vIntermediateGenParticles_;
   
   // TrackIsolation
   Bool_t electronTrackIsoVeto;
   Bool_t muonTrackIsoVeto;
   Bool_t pionTrackIsoVeto;

   // File Service to store things to a root file
   edm::Service<TFileService> fs_;

   // histogram to store #evts after each "cut"
   TH1F* hCutFlow_;
   std::map<std::string,TH1F*> hCutFlowMap_;

   // Pileup histogram(s)
   TH1F hPU_;
};

#endif /* TREEWRITER_HPP__ */
