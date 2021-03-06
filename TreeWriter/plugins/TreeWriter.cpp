// -*- C++ -*-
//
// Package:    TreeWriter/TreeWriter
// Class:      TreeWriter
//

#include "TreeWriter.hpp"

using namespace std;

// compute HT using RECO objects to "reproduce" the trigger requirements
static double computeHT(const std::vector<tree::Jet>& jets) {
   double HT = 0;
   double pt = 0;
   for (const tree::Jet& jet: jets) {
      pt = jet.p.Pt();
      if (fabs(jet.p.Eta())<3 && pt>30) {
         HT += pt;
      }
   }
   return HT;
}

// checks if a particle has a special mother. Treats anti-particles as particles
bool hasAncestor(int index, const lhef::HEPEUP& info, int searchId) {
   if (index < 2 or index > info.NUP) {
      return false;
   } else if (abs(info.IDUP[index]) == searchId) {
      return true;
   } else {
      auto mothers = info.MOTHUP[index];
      return
         (index != mothers.first-1 and hasAncestor(mothers.first-1, info, searchId))
         or (index != mothers.second-1 and hasAncestor(mothers.second-1, info, searchId));
   }
}

// taken from https://twiki.cern.ch/twiki/bin/view/CMS/SUSYRecommendationsForZinv#Photon_jets_control_region
PromptStatusType getPromptStatus(const reco::GenParticle& p, const edm::Handle<edm::View<reco::GenParticle>>& particles) {
   if (p.status()==1 && p.numberOfMothers() && (fabs(p.mother(0)->pdgId())<=22 || p.mother(0)->pdgId() == 2212)) {
      for (auto& genP : *particles) {
         auto absId = fabs(genP.pdgId());
         if (genP.status()==23 && ROOT::Math::VectorUtil::DeltaR(p.p4(), genP.p4())<0.4) {
             if (absId==21 || absId<6) return FRAGMENTPROMPT;
             if (absId==11 || absId==13 || absId==15) return LEPTONPROMPT;
         }
      }
      return DIRECTPROMPT;
   }
   return NOPROMPT;
};

float seedCrystalEnergyEB(const reco::SuperCluster& sc, const edm::Handle<EcalRecHitCollection>& ebRecHits) {
   float energy = -1;
   DetId detid = sc.seed()->seed();
   const EcalRecHit* rh = NULL;
   if (detid.subdetId() == EcalBarrel) {
      auto rh_i = ebRecHits->find(detid);
      if (rh_i != ebRecHits->end()) {
        rh = &(*rh_i);
      }
   }
   if (rh != NULL) {
     energy = rh->energy();
   }
   return energy;
}

//PFIsolation
double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> const &pfcands,
                        reco::Candidate const &ptcl,
                        double r_iso_min=0.05, double r_iso_max=0.2, double kt_scale=10,
                        bool charged_only=false) {
    // copied from https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY

    if (ptcl.pt()<5.) return 99999.;

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl.isElectron()) {
      if (fabs(ptcl.eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl.isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl.isElectron()) ptThresh = 0;
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl.pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, ptcl);
      if (dr > r_iso) continue;

      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
        /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    double iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl.pt();

    return iso;
}

// Calculate TrackIsolation
void GetTrkIso(edm::Handle<pat::PackedCandidateCollection> pfcands, const unsigned tkInd, float& trkiso, float& activity) {
  if (tkInd>pfcands->size()) {
	  trkiso = -999.;
	  activity = -999.;
	  return;
  }
  trkiso = 0.;
  activity = 0.;
  double r_iso = 0.3;
  for (unsigned int iPF(0); iPF<pfcands->size(); iPF++) {
    const pat::PackedCandidate &pfc = pfcands->at(iPF);
    if (pfc.charge()==0) continue;
    if (iPF==tkInd) continue; // don't count track in its own sum
    float dz_other = pfc.dz();
    if( fabs(dz_other) > 0.1 ) continue;
    double dr = deltaR(pfc, pfcands->at(tkInd));
    // activity annulus
    if (dr >= r_iso && dr <= 0.4) activity += pfc.pt();
    // mini iso cone
    if (dr <= r_iso) trkiso += pfc.pt();
  }
  trkiso = trkiso/pfcands->at(tkInd).pt();
  activity = activity/pfcands->at(tkInd).pt();
}

// Check TrackIsolation
bool TrackIsolation(edm::Handle<pat::PackedCandidateCollection> const &pfcands,edm::Handle<pat::METCollection> const &metColl,
                    edm::Handle<reco::VertexCollection> const &vertices, int const &pdgId_){
    
    //Get MET
    reco::MET::LorentzVector metLorentz(0,0,0,0);
    if (metColl.isValid()){
        metLorentz = metColl->at(0).p4();
    }
    
    //Good Vertices
    bool hasGoodVtx = false;
    if(vertices->size() > 0) hasGoodVtx = true;
    
    //Define MinPt cut and iso cut
    float minPt_ = 5.0;
    float isoCut_ = 0.2;
    if (pdgId_ == 211){
        minPt_ = 10.0;
        isoCut_ = 0.1;
    }
    
    //Define other cuts
    float maxEta_ = 2.5; //?????? correct for all three particle types
    float mTCut_ = 100.;
    float dzcut_ = 0.1;
    
    //Define TrackIso veto
    bool TrackIso = false;
    
    //loop over PFCandidates and calculate the trackIsolation
    for(size_t i=0; i<pfcands->size();i++)
    {
		const pat::PackedCandidate pfCand = (*pfcands)[i];
		
		//calculated mT value
		double dphiMET = fabs(pfCand.phi()-metLorentz.phi());
		double mT = sqrt(2 *metLorentz.pt() * pfCand.pt() * (1 - cos(dphiMET)));
		
		//-------------------------------------------------------------------------------------
		// skip events with no good vertices
		//-------------------------------------------------------------------------------------
		if(!hasGoodVtx) {
            continue;
		}
		//-------------------------------------------------------------------------------------
		// only consider charged tracks
		//-------------------------------------------------------------------------------------
		if (pfCand.charge() == 0) {
            continue;
		}
		//-------------------------------------------------------------------------------------
		// only store PFCandidate values if PFCandidate.pdgId() == pdgId_
		//-------------------------------------------------------------------------------------
		if( pdgId_ != 0 && abs( pfCand.pdgId() ) != pdgId_ ) {
            continue;
		}
		//-------------------------------------------------------------------------------------
		// only store PFCandidate values if pt > minPt
		//-------------------------------------------------------------------------------------
		if(pfCand.pt() <minPt_) {
            continue;
		}
		if(fabs(pfCand.eta()) >maxEta_) {
            continue;
		}
		//-------------------------------------------------------------------------------------
		// cut on mT of track and MET
		//-------------------------------------------------------------------------------------
		if(mTCut_>0.01 && mT>mTCut_) {
            continue;
		}
		//----------------------------------------------------------------------------
		// now make cuts on isolation and dz
		//----------------------------------------------------------------------------
		float trkiso = 0.;
		float activity = 0.;
		GetTrkIso(pfcands, i, trkiso, activity);
		float dz_it = pfCand.dz();

		if( isoCut_>0 && trkiso > isoCut_ ) {
            continue;
		}
		if( std::abs(dz_it) > dzcut_ ) {
            continue;
		}
        
        //Change TrackIso Veto
        TrackIso = true;
        }
    
    return TrackIso;
}

// Caracterize ttbar decay
const TLorentzVector nullP4_(0., 0., 0., 0.);
void GenLorentzVector(const reco::GenParticle* gen, TLorentzVector& TLV)
{
   TLV.SetPtEtaPhiE(gen->pt(), gen->eta(), gen->phi(), gen->energy());
}

void GenLorentzVector_Jet(const reco::GenJet* gen, TLorentzVector& TLV)
{
   TLV.SetPtEtaPhiE(gen->pt(), gen->eta(), gen->phi(), gen->energy());
}

const reco::GenParticle* tauDaughter(const reco::GenParticle* tau)
{
   for(size_t iDaughter = 0; iDaughter < tau->numberOfDaughters(); ++iDaughter){
      const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(tau->daughter(iDaughter));
      if(std::abs(daughter->pdgId())==11 || std::abs(daughter->pdgId())==13) return daughter;
      else if(abs(daughter->pdgId()) == 15) return tauDaughter(daughter);
   }
   return tau;
}

void assignLeptonAndTau(const reco::GenParticle* lepton, TLorentzVector& genLepton, int& pdgId, TLorentzVector& genTau)
{
   const reco::GenParticle* finalLepton;
   if(abs(lepton->pdgId()) == 15){
      GenLorentzVector(lepton,genTau);
      finalLepton = tauDaughter(lepton);
   }
   else{
      genTau = nullP4_;
      finalLepton = lepton;
   }

   if(abs(finalLepton->pdgId()) != 15){
      GenLorentzVector(finalLepton,genLepton);
      pdgId = finalLepton->pdgId();
   }
   else{
      genLepton = nullP4_;
      pdgId = 0;
   }
}

bool isWDaughter(const reco::GenParticle& particle)
{
  return (std::abs(particle.mother()->pdgId()) == 24);
}

bool isTopDaughter(const reco::GenParticle& particle)
{
  return (std::abs(particle.mother()->pdgId()) == 6);
}


template <typename T> int sign(T val) {
   return (T(0) < val) - (val < T(0));
}

TreeWriter::TreeWriter(const edm::ParameterSet& iConfig)
   : dHT_cut_(iConfig.getUntrackedParameter<double>("HT_cut"))
   , dJet_pT_cut_(iConfig.getUntrackedParameter<double>("jet_pT_cut"))
   , minNumberElectrons_cut_(iConfig.getUntrackedParameter<unsigned>("minNumberElectrons_cut"))
   , NumberLeptons_cut_(iConfig.getUntrackedParameter<unsigned>("NumberLeptons_cut"))
   , newLumiBlock_(true)
   , vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
   , jetCollectionToken_     (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
   , jetPuppiCollectionToken_     (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets_puppi")))
   , genJetCollectionToken_  (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")))
   , muonCollectionToken_    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
   , electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons")))
   , photonCollectionToken_  (consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons")))
   , metCollectionToken_     (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
   , metPuppiCollectionToken_     (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi")))
   , metNoHFCollectionToken_     (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF")))
   , metCorrectedCollectionToken_  (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metCorrected")))
   , metCalibratedCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metCalibrated")))
   , metDeepCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsDeepMET")))
   , metBJetRegressionCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_BJetRegression")))
   , metBJetRegressionLooseCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_BJetRegressionLoose")))
   , caloMetCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("caloMets")))
   , rhoToken_               (consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
   , ebRecHitsToken_         (consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRecHits")))
   , prunedGenToken_         (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles")))
   , pileUpSummaryToken_     (consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileUpSummary")))
   , LHEEventToken_          (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProduct")))
   , packedCandidateToken_   (consumes<std::vector<pat::PackedCandidate>> (iConfig.getParameter<edm::InputTag>("packedCandidates")))
   // electron id
   , electronVetoIdMapToken_  (iConfig.getUntrackedParameter<std::string>("electronVetoIdMap"))
   , electronLooseIdMapToken_ (iConfig.getUntrackedParameter<std::string>("electronLooseIdMap"))
   , electronMediumIdMapToken_(iConfig.getUntrackedParameter<std::string>("electronMediumIdMap"))
   , electronTightIdMapToken_ (iConfig.getUntrackedParameter<std::string>("electronTightIdMap"))
   // photon id
   // ~, photonLooseIdMapToken_  (iConfig.getUntrackedParameter<std::string>("photonLooseIdMap"  ))
   // ~, photonMediumIdMapToken_ (iConfig.getUntrackedParameter<std::string>("photonMediumIdMap" ))
   // ~, photonTightIdMapToken_  (iConfig.getUntrackedParameter<std::string>("photonTightIdMap"  ))
   // met filters to apply
   , metFilterNames_(iConfig.getUntrackedParameter<std::vector<std::string>>("metFilterNames"))
   , phoWorstChargedIsolationToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoWorstChargedIsolation")))
   , pileupHistogramName_(iConfig.getUntrackedParameter<std::string>("pileupHistogramName"))
   , hardPUveto_(iConfig.getUntrackedParameter<bool>("hardPUveto"))
   , reMiniAOD_(iConfig.getUntrackedParameter<bool>("reMiniAOD"))
   , jetIdSelector(iConfig.getParameter<edm::ParameterSet>("pfJetIDSelector"))
   , triggerNames_(iConfig.getParameter<std::vector<std::string>>("triggerNames"))
   , triggerPrescales_(iConfig.getParameter<std::vector<std::string>>("triggerPrescales"))
   , triggerObjectNames_(iConfig.getParameter<std::vector<std::string>>("triggerObjectNames"))
   // scale factor map
   , fctLeptonFullSimScaleFactors_(iConfig.getParameter<edm::ParameterSet>("LeptonFullSimScaleFactors"))
   , fctBTagEff_    (iConfig.getParameter<edm::ParameterSet>("bTagEfficiencies") )
   , fctBTagCalibFullSim_    (iConfig.getParameter<edm::ParameterSet>("BTagCalibration").getParameter<std::string>("CSVFullSimTagger"),iConfig.getParameter<edm::ParameterSet>("BTagCalibration").getParameter<std::string>("CSVFullSimFileName") )
   , fctBTagCalibReaderFullSimBJets_    (BTagEntry::OP_LOOSE,"central" )
   , fctBTagCalibReaderFullSimCJets_    (BTagEntry::OP_LOOSE,"central" )
   , fctBTagCalibReaderFullSimLightJets_    (BTagEntry::OP_LOOSE,"central" )
   , fctBTagCalibReaderFullSimBJetsUp_    (BTagEntry::OP_LOOSE,"up" )
   , fctBTagCalibReaderFullSimCJetsUp_    (BTagEntry::OP_LOOSE,"up" )
   , fctBTagCalibReaderFullSimLightJetsUp_    (BTagEntry::OP_LOOSE,"up" )
   , fctBTagCalibReaderFullSimBJetsDown_    (BTagEntry::OP_LOOSE,"down" )
   , fctBTagCalibReaderFullSimCJetsDown_    (BTagEntry::OP_LOOSE,"down" )
   , fctBTagCalibReaderFullSimLightJetsDown_    (BTagEntry::OP_LOOSE,"down" )
   // Ttbar gen Event Info
   , ttbarGenInfo_(iConfig.getParameter<bool>("ttbarGenInfo"))
   , pseudoTopInfo_(iConfig.getParameter<bool>("ttbarPseudoInfo"))
   , dyPtInfo_(iConfig.getParameter<bool>("DYptInfo"))
{
   // declare consumptions that are used "byLabel" in analyze()
   mayConsume<GenLumiInfoHeader,edm::InLumi> (edm::InputTag("generator"));
   consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "",edm::InputTag::kSkipCurrentProcess));
   consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
   consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("slimmedPatTrigger"));
   consumes<bool>(edm::InputTag("particleFlowEGammaGSFixed", "dupECALClusters"));
   consumes<edm::EDCollection<DetId>>(edm::InputTag("ecalMultiAndGSGlobalRecHitEB", "hitsNotReplaced"));
   consumes<TtGenEvent>(edm::InputTag("genEvt"));
   consumes<int>(edm::InputTag("generatorTopFilter","decayMode"));
   consumes<std::vector<reco::GenParticle>>(edm::InputTag("pseudoTop"));
   //~consumes<reco::GenJetCollection>(edm::InputTag("pseudoTop","leptons"));
   //~consumes<reco::GenJetCollection>(edm::InputTag("pseudoTop","jets"));
   
   // load BTagCalibrations
   fctBTagCalibReaderFullSimBJets_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_B,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_bJets"));
   fctBTagCalibReaderFullSimCJets_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_C,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_cJets"));
   fctBTagCalibReaderFullSimLightJets_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_UDSG,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_lightJets"));
   fctBTagCalibReaderFullSimBJetsUp_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_B,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_bJets"));
   fctBTagCalibReaderFullSimCJetsUp_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_C,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_cJets"));
   fctBTagCalibReaderFullSimLightJetsUp_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_UDSG,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_lightJets"));
   fctBTagCalibReaderFullSimBJetsDown_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_B,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_bJets"));
   fctBTagCalibReaderFullSimCJetsDown_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_C,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_cJets"));
   fctBTagCalibReaderFullSimLightJetsDown_.load(fctBTagCalibFullSim_,BTagEntry::FLAV_UDSG,iConfig.getParameter<edm::ParameterSet>("BTagCalibrationReader").getParameter<std::string>("measurementType_lightJets"));

   // setup tree and define branches
   eventTree_ = fs_->make<TTree> ("eventTree", "event data");
   eventTree_->Branch("jets", &vJets_);
   eventTree_->Branch("jetsPuppi", &vJetsPuppi_);
   eventTree_->Branch("genJets", &vGenJets_);
   eventTree_->Branch("electrons", &vElectrons_);
   eventTree_->Branch("muons", &vMuons_);
   // ~eventTree_->Branch("electrons_add", &vElectrons_add_);
   // ~eventTree_->Branch("muons_add", &vMuons_add_);
   // ~eventTree_->Branch("photons", &vPhotons_);
   eventTree_->Branch("met", &met_);
   eventTree_->Branch("metCalo", &metCalo_);
   eventTree_->Branch("metPuppi", &metPuppi_);
   eventTree_->Branch("metNoHF", &metNoHF_);
   // ~eventTree_->Branch("metDeep", &metDeep_);
   eventTree_->Branch("metXYcorr", &metXYcorr_);
   // ~eventTree_->Branch("metBJetRegression", &metBJetRegression_);
   // ~eventTree_->Branch("metBJetRegressionLoose", &metBJetRegressionLoose_);
   eventTree_->Branch("met_raw", &met_raw_);
   eventTree_->Branch("met_gen", &met_gen_);
   eventTree_->Branch("met_JESu", &met_JESu_);
   eventTree_->Branch("met_JESd", &met_JESd_);
   eventTree_->Branch("met_JERu", &met_JERu_);
   eventTree_->Branch("met_JERd", &met_JERd_);
   eventTree_->Branch("genParticles", &vGenParticles_);
   for (const auto& n : triggerObjectNames_) {  //Trigger branches
     triggerObjectMap_[n] = std::vector<tree::Particle>();
     eventTree_->Branch(n.c_str(), &triggerObjectMap_[n]);
   }
   eventTree_->Branch("intermediateGenParticles", &vIntermediateGenParticles_);
   eventTree_->Branch("ee", &ee_);
   eventTree_->Branch("mumu", &mumu_);
   eventTree_->Branch("emu", &emu_);
   eventTree_->Branch("addLepton", &addLepton_);
   eventTree_->Branch("true_nPV", &true_nPV_, "true_nPV/I");
   eventTree_->Branch("nPV", &nPV_, "nPV/I");
   eventTree_->Branch("nGoodVertices" , &nGoodVertices_ , "nGoodVertices/I");
   eventTree_->Branch("nTracksPV", &nTracksPV_, "nTracksPV/I");
   eventTree_->Branch("rho", &rho_, "rho/F");
   eventTree_->Branch("caloMetPt", &caloMetPt_, "caloMetPt/F");

   eventTree_->Branch("pu_weight", &pu_weight_, "pu_weight/F");
   eventTree_->Branch("mc_weight", &mc_weight_, "mc_weight/B");
   eventTree_->Branch("pdf_weights", &vPdf_weights_);

   eventTree_->Branch("mll", &mll_, "mll/F");
   eventTree_->Branch("Ht", &Ht_, "Ht/F");
   eventTree_->Branch("genHt", &genHt_, "genHt/F");
   eventTree_->Branch("EWKinoPairPt", &EWKinoPairPt_, "EWKinoPairPt/F");
   eventTree_->Branch("MT2", &MT2_, "MT2/F");
   eventTree_->Branch("genMT2", &genMT2_, "genMT2/F");
   eventTree_->Branch("genMT2neutrino", &genMT2neutrino_, "genMT2neutrino/F");
   
   eventTree_->Branch("lepton1SF", &lepton1SF_, "lepton1SF/F");
   eventTree_->Branch("lepton2SF", &lepton2SF_, "lepton2SF/F");
   eventTree_->Branch("lepton1SF_unc", &lepton1SF_unc_, "lepton1SF_unc/F");
   eventTree_->Branch("lepton2SF_unc", &lepton2SF_unc_, "lepton2SF_unc/F");
   
   eventTree_->Branch("bTagWeight", &bTagWeight_, "bTagWeight/F");
   eventTree_->Branch("bTagWeightErrHeavy", &bTagWeightErrHeavy_, "bTagWeightErrHeavy/F");
   eventTree_->Branch("bTagWeightErrLight", &bTagWeightErrLight_, "bTagWeightErrLight/F");
   
   eventTree_->Branch("topPTweight", &topPTweight_, "topPTweight/F");

   eventTree_->Branch("evtNo", &evtNo_, "evtNo/l");
   eventTree_->Branch("runNo", &runNo_, "runNo/i");
   eventTree_->Branch("lumNo", &lumNo_, "lumNo/i");
   eventTree_->Branch("gs_duplicate", &particleFlowEGammaGSFixed_dupECALClusters_, "gs_duplicate/O");
   eventTree_->Branch("gs_notReplaced", &ecalMultiAndGSGlobalRecHitEB_hitsNotReplaced_, "gs_notReplaced/O");

   eventTree_->Branch("signal_m1", &signal_m1_, "signal_m1/s");
   eventTree_->Branch("signal_m2", &signal_m2_, "signal_m2/s");
   eventTree_->Branch("signal_nBinos", &signal_nBinos_, "signal_nBinos/s");
   eventTree_->Branch("signal_nNeutralinoDecays", &signal_nNeutralinoDecays_, "signal_nNeutralinoDecays/s");
   
   eventTree_->Branch("electronTrackIsoVeto", &electronTrackIsoVeto);
   eventTree_->Branch("muonTrackIsoVeto", &muonTrackIsoVeto);
   eventTree_->Branch("pionTrackIsoVeto", &pionTrackIsoVeto);
   
   eventTree_->Branch("ttbarProductionMode", &ttbarProductionMode_);
   eventTree_->Branch("ttbarDecayMode", &ttbarDecayMode_);
   eventTree_->Branch("genTop", &genTop_);
   eventTree_->Branch("genAntiTop", &genAntiTop_);
   eventTree_->Branch("genLepton", &genLepton_);
   eventTree_->Branch("genAntiLepton", &genAntiLepton_);
   eventTree_->Branch("genTau", &genTau_);
   eventTree_->Branch("genAntiTau", &genAntiTau_);
   eventTree_->Branch("genLeptonPdgId", &genLeptonPdgId_);
   eventTree_->Branch("genAntiLeptonPdgId", &genAntiLeptonPdgId_);
   eventTree_->Branch("genB", &genB_);
   eventTree_->Branch("genAntiB", &genAntiB_);
   eventTree_->Branch("genNeutrino", &genNeutrino_);
   eventTree_->Branch("genAntiNeutrino", &genAntiNeutrino_);
   eventTree_->Branch("genWMinus", &genWMinus_);
   eventTree_->Branch("genWPlus", &genWPlus_);
   
   eventTree_->Branch("ttbarPseudoDecayMode", &ttbarPseudoDecayMode_);
   eventTree_->Branch("pseudoTop", &pseudoTop_);
   eventTree_->Branch("pseudoAntiTop", &pseudoAntiTop_);
   eventTree_->Branch("pseudoLepton", &pseudoLepton_);
   eventTree_->Branch("pseudoAntiLepton", &pseudoAntiLepton_);
   eventTree_->Branch("pseudoTau", &pseudoTau_);
   eventTree_->Branch("pseudoAntiTau", &pseudoAntiTau_);
   eventTree_->Branch("pseudoLeptonPdgId", &pseudoLeptonPdgId_);
   eventTree_->Branch("pseudoAntiLeptonPdgId", &pseudoAntiLeptonPdgId_);
   eventTree_->Branch("pseudoBJet", &pseudoBJet_);
   eventTree_->Branch("pseudoAntiBJet", &pseudoAntiBJet_);
   eventTree_->Branch("pseudoNeutrino", &pseudoNeutrino_);
   eventTree_->Branch("pseudoAntiNeutrino", &pseudoAntiNeutrino_);
   eventTree_->Branch("pseudoWMinus", &pseudoWMinus_);
   eventTree_->Branch("pseudoWPlus", &pseudoWPlus_);
   //~eventTree_->Branch("allPseudoJets", &v_allPseudoJet_);
   //~eventTree_->Branch("allPseudoLeptons", &v_allPseudoLepton_);
   
   eventTree_->Branch("DYgenPt", &dyPt_, "DYgenPt/F");
   eventTree_->Branch("DYgenLep", &v_dyLep_);

   // Fill trigger maps
   for (const auto& n : triggerNames_) {
      triggerIndex_[n] = -10; // not set and not found
      triggerDecision_[n] = false;
      eventTree_->Branch(n.c_str(), &triggerDecision_[n], (n+"/O").c_str());
   }
   // create branches for prescales
   for (std::string const& n : triggerPrescales_) {
      std::string const name = n + "_pre";
      eventTree_->Branch(name.c_str(), &triggerPrescale_[n], (name+"/I").c_str());
   }

   // get pileup histogram(s)
   std::string cmssw_base_src = getenv("CMSSW_BASE");
   cmssw_base_src += "/src/";
   TFile puFile(TString(cmssw_base_src+"/TreeWriter/PUreweighting/data/puWeights.root"));
   if (puFile.IsZombie()) {
      edm::LogError("File not found") << "create puWeights.root! (see README)";
      std::exit(84);
   } else {
      TH1F* hPU_ptr = (TH1F*)puFile.Get(pileupHistogramName_.c_str());
      if (hPU_ptr) {
         hPU_ = *hPU_ptr;
      } else {
         edm::LogError("Pileup histogram not found") << "recreate puWeights.root! (see README)";
         std::exit(84);
      }
   }
   puFile.Close();
   
}

TH1F* TreeWriter::createCutFlowHist(std::string modelName)
{
   std::string const name("hCutFlow" + modelName);
   std::vector<TString> vCutBinNames{{
      "initial_unweighted",
      "initial_mc_weighted",
      "initial",
      "trigger",
      "METfilters",
      "nGoodVertices",
      "Dilepton",
      "HT",
      "final"}};
   TH1F* h = fs_->make<TH1F>(name.c_str(), name.c_str(), vCutBinNames.size(), 0, vCutBinNames.size());
   for (uint i=0; i<vCutBinNames.size(); i++) { h->GetXaxis()->SetBinLabel(i+1, vCutBinNames.at(i)); }
   return h;
}

// ------------ method called for each event  ------------
void TreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Bool_t isRealData;
   isRealData = iEvent.isRealData();

   hCutFlow_->Fill("initial_unweighted", 1);
   
   ///////////////////
   // PileUp weights//
   ///////////////////
   puPtHat_ = 0;
   if (!isRealData) {
      edm::Handle<PileupSummaryInfoCollection> PupInfo;
      iEvent.getByToken(pileUpSummaryToken_, PupInfo);
      float Tnpv = -1;
      for (auto const& PVI: *PupInfo) {
         int BX = PVI.getBunchCrossing();
         if (BX == 0) {
            Tnpv = PVI.getTrueNumInteractions();
            auto ptHats = PVI.getPU_pT_hats();
            puPtHat_ = ptHats.size() ? *max_element(ptHats.begin(),ptHats.end()) : 0;
            continue;
         }
      }
      true_nPV_ = Tnpv;
      pu_weight_ = hPU_.GetBinContent(hPU_.FindBin(Tnpv));
   } else { // real data
      true_nPV_ = -1;
      pu_weight_ = 1.;
   }
   
   //////////////////////
   // generator weights//
   //////////////////////
   mc_weight_ = 1; // 1 for data
   if (!isRealData) {
      edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
      iEvent.getByLabel("generator", GenEventInfoHandle);
      mc_weight_ = sign(GenEventInfoHandle->weight());
      auto weightsize = GenEventInfoHandle->weights().size();
      if (weightsize < 2) {   // for most SM samples
         edm::Handle<LHEEventProduct> LHEEventProductHandle;
         iEvent.getByToken(LHEEventToken_, LHEEventProductHandle);
         if (LHEEventProductHandle.isValid()) {
            unsigned iMax = 110; // these are 9 scale variations and 100 variation of the first pdf set
            if (iMax>LHEEventProductHandle->weights().size()) iMax = LHEEventProductHandle->weights().size();
            vPdf_weights_ = std::vector<float>(iMax, 1.0);
            for (unsigned i=0; i<iMax; i++) {
               // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface#Retrieving_information_on_LHE_ev
               vPdf_weights_[i] = LHEEventProductHandle->weights()[i].wgt/LHEEventProductHandle->originalXWGTUP();
            }
         }
      } else { // for SMS scans
         unsigned iMax = 110;
         if (iMax>GenEventInfoHandle->weights().size()-1) iMax=GenEventInfoHandle->weights().size()-1;
         vPdf_weights_ = std::vector<float>(iMax, 1.0);
         for (unsigned i=0; i<iMax; i++) {
            // 0 and 1 are the same for 80X scans
            // https://hypernews.cern.ch/HyperNews/CMS/get/susy-interpretations/242/1/1.html
            vPdf_weights_[i] = GenEventInfoHandle->weights()[i+1]/GenEventInfoHandle->weights()[1];
         }
      }
   }

   hCutFlow_->Fill("initial_mc_weighted", mc_weight_);
   hCutFlow_->Fill("initial", mc_weight_*pu_weight_);
   
   ////////////
   // Trigger//
   ////////////
   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   edm::InputTag triggerTag("TriggerResults", "", "HLT");
   edm::InputTag triggerPrescaleTag("patTrigger");
   iEvent.getByLabel(triggerTag, triggerBits);
   iEvent.getByLabel(triggerPrescaleTag, triggerPrescales);

   // for each lumiBlock, re-read the trigger indices (rather changes for new run)
   const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
   if (triggerIndex_.size() && newLumiBlock_) {
      newLumiBlock_ = false;
      // set all trigger indeces to -1 as "not available"-flag
      for (auto& it: triggerIndex_) { it.second = -1; }
      // store the indices of the trigger names that we really find
      for (unsigned i=0; i<triggerNames.size(); i++) {
         for (auto& it : triggerIndex_) {
            if (triggerNames.triggerName(i).find(it.first) == 0) {
               it.second = i;
            }
         }
      } // end trigger names
   } // found indices

   // set trigger decision
   bool anyTriggerFired = false;
   for (auto& it : triggerIndex_) {
      if (it.second != -1) {
         anyTriggerFired |= triggerBits->accept( it.second );
         triggerDecision_[it.first] = triggerBits->accept( it.second );
      }
   }
   if (isRealData && !anyTriggerFired) return;
   hCutFlow_->Fill("trigger", mc_weight_*pu_weight_);

   // store prescales
   for (std::string const& n : triggerPrescales_) {
      int const index = triggerIndex_[n];
      // if the index was not found, store '0': trigger was not run!
      triggerPrescale_[n] = index == -1 ? 0 : triggerPrescales->getPrescaleForIndex(triggerIndex_[n]);
   }

   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   // ~edm::InputTag triggerObjects_("selectedPatTrigger"); //80x
   edm::InputTag triggerObjects_("slimmedPatTrigger");
   iEvent.getByLabel(triggerObjects_, triggerObjects);

   for (const auto& n : triggerObjectNames_) triggerObjectMap_.at(n).clear();
   tree::Particle trObj;
   for (pat::TriggerObjectStandAlone obj: *triggerObjects) {
      obj.unpackFilterLabels(iEvent,*triggerBits);
      for (const auto& n : triggerObjectNames_) {
         if (std::count(obj.filterLabels().begin(), obj.filterLabels().end(), n)) {
            //~ trObj.p.SetPtEtaPhi(obj.pt(), obj.eta(), obj.phi());
            trObj.p.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
            triggerObjectMap_.at(n).push_back(trObj);
         }
      }
   }
   
   ////////////////
   // MET Filters//
   ////////////////
   edm::Handle<edm::TriggerResults> metFilterBits;
   edm::InputTag metFilterTag("TriggerResults", "",edm::InputTag::kSkipCurrentProcess);
   iEvent.getByLabel(metFilterTag, metFilterBits);
   // go through the filters and check if they were passed
   const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
   //for (unsigned i=0; i<allFilterNames.size(); i++) std::cout << allFilterNames.triggerName(i) << std::endl;
   for (std::string const &name: metFilterNames_) {
      const unsigned index = allFilterNames.triggerIndex(name);
      if (index >= allFilterNames.size()) std::cerr << "MET filter '" << name << "' not found!" << std::endl;
      if (!metFilterBits->accept(index)) return; // not passed
   }
   hCutFlow_->Fill("METfilters", mc_weight_*pu_weight_);
   
   ////////////////////
   // Primary Vertex//
   ///////////////////
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   nPV_ = vertices->size();

   reco::Vertex firstGoodVertex;
   nGoodVertices_ = 0;
   for (const auto& vtx: *vertices) {
      // from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_14/doc/html/db/d49/GoodVertexFilter_8cc_source.html
      if (vtx.ndof()>4
         && vtx.position().Rho()<=2.0
         && fabs(vtx.position().Z())<=24.0)
      {
         nGoodVertices_++;
         // first one?
         if (nGoodVertices_ == 1) firstGoodVertex = vtx;
      }
   }
   if (!nGoodVertices_) return;
   hCutFlow_->Fill("nGoodVertices", mc_weight_*pu_weight_);

   // Get rho
   edm::Handle< double > rhoH;
   iEvent.getByToken(rhoToken_, rhoH);
   rho_ = *rhoH;
   
   /////////////////////
   // number of tracks//
   /////////////////////
   edm::Handle<std::vector<pat::PackedCandidate>> packedCandidates;
   iEvent.getByToken(packedCandidateToken_, packedCandidates);
   nTracksPV_ = std::count_if(packedCandidates->begin(),packedCandidates->end(), [] (const pat::PackedCandidate& cand) {
      return cand.pt()>.9 && cand.charge() && cand.pvAssociationQuality() == pat::PackedCandidate::UsedInFitTight && cand.fromPV() == pat::PackedCandidate::PVUsedInFit;});

   // get gen particles before photons for the truth match
   edm::Handle<edm::View<reco::GenParticle>> prunedGenParticles;
   if (!isRealData) { iEvent.getByToken(prunedGenToken_,prunedGenParticles); }

   edm::Handle<edm::ValueMap<float>> phoWorstChargedIsolationMap;
   iEvent.getByToken(phoWorstChargedIsolationToken_, phoWorstChargedIsolationMap);

   edm::Handle<EcalRecHitCollection> ebRecHits;
   iEvent.getByToken(ebRecHitsToken_, ebRecHits);
   
   //////////
   // Muons//
   //////////
   edm::Handle<pat::MuonCollection> muonColl;
   iEvent.getByToken(muonCollectionToken_, muonColl);

   vMuons_.clear();
   vMuons_add_.clear();
   tree::Muon trMuon;
   for (const pat::Muon &mu : *muonColl) {
      if (! (mu.isPFMuon() || mu.isGlobalMuon() || mu.isTrackerMuon())) continue; //(can probably be removed, not sure about the veto with looser lepton)
      // ~if (mu.pt()<20 || abs(mu.eta())>2.4) continue;
      if (abs(mu.eta())>2.4) continue;
      //~ trMuon.p.SetPtEtaPhi(mu.pt(), mu.eta(), mu.phi());
      trMuon.p.SetPtEtaPhiE(mu.pt(), mu.eta(), mu.phi(), mu.energy());
      trMuon.isTight = mu.isTightMuon(firstGoodVertex);
      trMuon.isMedium = mu.isMediumMuon();
      trMuon.isLoose = mu.isLooseMuon();
      auto const& pfIso = mu.pfIsolationR04();
      trMuon.rIso = (pfIso.sumChargedHadronPt + std::max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))/mu.pt();
      trMuon.charge = mu.charge();
      trMuon.PFminiIso = getPFIsolation(packedCandidates, mu);
      math::XYZPoint vtx_point = firstGoodVertex.position();
      trMuon.d0 = mu.bestTrack()->dxy( vtx_point );
      trMuon.dZ = mu.bestTrack()->dz( vtx_point );
      trMuon.rochesterCorrection = mu.hasUserFloat("MuonEnergyCorr") ? mu.userFloat("MuonEnergyCorr") : 1.;

      if (mu.pt()*trMuon.rochesterCorrection>20 && mu.isTightMuon(firstGoodVertex) && trMuon.rIso<0.15) vMuons_.push_back(trMuon); // take only 'tight' muons
      else if (mu.pt()*trMuon.rochesterCorrection>10 && trMuon.isLoose) vMuons_add_.push_back(trMuon); //Save all additional muons, which are at least loose
   } // muon loop
   sort(vMuons_.begin(), vMuons_.end(), tree::PtGreater);
   sort(vMuons_add_.begin(), vMuons_add_.end(), tree::PtGreater);
   
   /////////////
   //Electrons//
   /////////////
   edm::Handle<edm::View<pat::Electron>> electronColl;
   iEvent.getByToken(electronCollectionToken_, electronColl);

   vElectrons_.clear();
   vElectrons_add_.clear();
   tree::Electron trEl;
   for (edm::View<pat::Electron>::const_iterator el = electronColl->begin();el != electronColl->end(); el++) {
      // ~if (el->pt()<20 || abs(el->superCluster()->eta())>2.4 || ((1.4442 < abs(el->superCluster()->eta())) && abs(el->superCluster()->eta()) < 1.5660)) continue;
      if (abs(el->superCluster()->eta())>2.4 || ((1.4442 < abs(el->superCluster()->eta())) && abs(el->superCluster()->eta()) < 1.5660)) continue;
      const edm::Ptr<pat::Electron> elPtr(electronColl, el - electronColl->begin());
      trEl.isLoose = el->electronID(electronLooseIdMapToken_);
      trEl.isMedium = el->electronID(electronMediumIdMapToken_);
      trEl.isTight = el->electronID(electronTightIdMapToken_);
      //~ trEl.p.SetPtEtaPhi(el->pt(), el->superCluster()->eta(), el->superCluster()->phi());
      //~ trEl.p.SetPtEtaPhiE(el->pt(), el->superCluster()->eta(), el->superCluster()->phi(), el->energy());
      trEl.p.SetPtEtaPhiE(el->pt(), el->eta(), el->phi(), el->energy());
      trEl.seedCrystalE = seedCrystalEnergyEB(*el->superCluster(), ebRecHits);
      trEl.charge = el->charge();
      auto const & pfIso = el->pfIsolationVariables();
      trEl.rIso = (pfIso.sumChargedHadronPt + std::max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))/el->pt();
      trEl.r9 = el->r9();
      trEl.SigmaIEtaIEtaFull5x5 = el->full5x5_sigmaIetaIeta();
      trEl.dPhiAtVtx = el->deltaPhiSuperClusterTrackAtVtx();
      trEl.dEtaAtVtx = el->deltaEtaSuperClusterTrackAtVtx();
      trEl.HoverE = el->hcalOverEcal();
      trEl.MissHits = el->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
      trEl.ConvVeto = el->passConversionVeto();
      trEl.PFminiIso = getPFIsolation(packedCandidates, *el);
      math::XYZPoint vtx_point = firstGoodVertex.position();
      trEl.d0 = el->bestTrack()->dxy( vtx_point );
      trEl.dZ = el->bestTrack()->dz( vtx_point );
      trEl.phiSC = el->superCluster()->phi();
      trEl.etaSC = el->superCluster()->eta();
      trEl.corr = el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
            
      // VID calculation of (1/E - 1/p)
      if (el->ecalEnergy() == 0)   trEl.EoverPInv = 1e30;
      else if (!std::isfinite(el->ecalEnergy()))  trEl.EoverPInv = 1e30;
      else trEl.EoverPInv = (1.0 - el->eSuperClusterOverP())/el->ecalEnergy();
      
      //uncorrected electrons
      trEl.pUncorrected.SetXYZ(0,0,0);
      //~ auto itPos = std::distance(electronColl->begin(), el); 80x
      //~ if (itPos<electronCollUncorrected->size()) {
        //~ auto ele = (*electronCollUncorrected).at(itPos);
        //~ trEl.pUncorrected.SetPtEtaPhi(ele.pt(), ele.superCluster()->eta(), ele.superCluster()->phi());
      //~ }
		
      //Apply recommended IP cuts
      if((abs(trEl.etaSC) < 1.4442) && ((abs(trEl.d0) > 0.05) || abs(trEl.dZ) > 0.10)) continue;
      else if((abs(trEl.etaSC) > 1.5660) && ((abs(trEl.d0) > 0.10) || abs(trEl.dZ) > 0.20)) continue;
            
      if (el->pt()*trEl.corr>20 && el->electronID(electronTightIdMapToken_)) vElectrons_.push_back(trEl); // take only 'tight' electrons
      else if (el->pt()*trEl.corr>10 && el->electronID(electronVetoIdMapToken_)) vElectrons_add_.push_back(trEl); //Save all additional electrons, which are at least 'veto' electrons
   }
   sort(vElectrons_.begin(), vElectrons_.end(), tree::PtGreater);
   sort(vElectrons_add_.begin(), vElectrons_add_.end(), tree::PtGreater);
   
   /////////////////////////////
   // Pseudo TTbar Information//
   /////////////////////////////
   if(pseudoTopInfo_){
      edm::Handle<std::vector<reco::GenParticle> > pseudoTopQuarks;
      iEvent.getByLabel("pseudoTop", pseudoTopQuarks);
      if(!pseudoTopQuarks.failedToGet()){
         pseudoTop_ = nullP4_;
         pseudoAntiTop_ = nullP4_;
         pseudoBJet_ = nullP4_;
         pseudoAntiBJet_ = nullP4_;
         pseudoWPlus_ = nullP4_;
         pseudoWMinus_ = nullP4_;
         pseudoLepton_ = nullP4_;
         pseudoAntiLepton_ = nullP4_;
         pseudoLeptonPdgId_ = 0;
         pseudoAntiLeptonPdgId_ = 0;
         pseudoNeutrino_ = nullP4_;
         pseudoAntiNeutrino_ = nullP4_;
         ttbarPseudoDecayMode_ = 0;

         for(std::vector<reco::GenParticle>::const_iterator i_particle = pseudoTopQuarks->begin(); i_particle != pseudoTopQuarks->end(); ++i_particle)
         {
            if(i_particle->pdgId() == 6) GenLorentzVector(&(*i_particle),pseudoTop_);
            if(i_particle->pdgId() == -6) GenLorentzVector(&(*i_particle),pseudoAntiTop_);
            if(i_particle->pdgId() == 5 && isTopDaughter(*i_particle)) GenLorentzVector(&(*i_particle),pseudoBJet_);
            if(i_particle->pdgId() == -5 && isTopDaughter(*i_particle)) GenLorentzVector(&(*i_particle),pseudoAntiBJet_);
            if(i_particle->pdgId() == 24 && isTopDaughter(*i_particle)) GenLorentzVector(&(*i_particle),pseudoWPlus_);
            if(i_particle->pdgId() == -24 && isTopDaughter(*i_particle)) GenLorentzVector(&(*i_particle),pseudoWMinus_);
            if((i_particle->pdgId() == 11 || i_particle->pdgId() == 13) && isWDaughter(*i_particle)) {GenLorentzVector(&(*i_particle),pseudoLepton_); pseudoLeptonPdgId_ = i_particle->pdgId();}
            if((i_particle->pdgId() == -11 || i_particle->pdgId() == -13) && isWDaughter(*i_particle)) {GenLorentzVector(&(*i_particle),pseudoAntiLepton_); pseudoAntiLeptonPdgId_ = i_particle->pdgId();}
            if((i_particle->pdgId() == 12 || i_particle->pdgId() == 14 || i_particle->pdgId() == 16) && isWDaughter(*i_particle)) GenLorentzVector(&(*i_particle),pseudoNeutrino_);
            if((i_particle->pdgId() == -12 || i_particle->pdgId() == -14 || i_particle->pdgId() == -16) && isWDaughter(*i_particle)) GenLorentzVector(&(*i_particle),pseudoAntiNeutrino_);
         }

         if(std::abs(pseudoLeptonPdgId_) == 11 && std::abs(pseudoAntiLeptonPdgId_) == 11) ttbarPseudoDecayMode_ = 1; // ee
         else if(std::abs(pseudoLeptonPdgId_) == 13 && std::abs(pseudoAntiLeptonPdgId_) == 13) ttbarPseudoDecayMode_ = 2; // mumu
         else if(std::abs(pseudoLeptonPdgId_*pseudoAntiLeptonPdgId_) == 143) ttbarPseudoDecayMode_ = 3; // emu

      }
      else {
         std::cerr<<"\nError: no pseudo top?!\n\n";
         ttbarPseudoDecayMode_ = -1;
         pseudoTop_ = nullP4_;
         pseudoAntiTop_ = nullP4_;
         pseudoBJet_ = nullP4_;
         pseudoAntiBJet_ = nullP4_;
         pseudoWPlus_ = nullP4_;
         pseudoWMinus_ = nullP4_;
         pseudoLepton_ = nullP4_;
         pseudoAntiLepton_ = nullP4_;
         pseudoLeptonPdgId_ = 0;
         pseudoAntiLeptonPdgId_ = 0;
         pseudoNeutrino_ = nullP4_;
         pseudoAntiNeutrino_ = nullP4_;
      }
      
      //~// All particle jets
      //~v_allPseudoJet_.clear();
      //~tree::Particle tempParticle_;
      //~edm::Handle<reco::GenJetCollection> pseudoJets;
      //~iEvent.getByLabel("pseudoTop","jets", pseudoJets);
      //~for(std::vector<reco::GenJet>::const_iterator i_jet = pseudoJets->begin(); i_jet != pseudoJets->end(); ++i_jet){
         //~GenLorentzVector_Jet(&(*i_jet),tempParticle_.p);
         //~v_allPseudoJet_.push_back(tempParticle_);
      //~}

      //~// All particle leptons
      //~v_allPseudoLepton_.clear();
      //~edm::Handle<reco::GenJetCollection> pseudoLeptons;
      //~iEvent.getByLabel("pseudoTop","leptons",pseudoLeptons);
      //~for(std::vector<reco::GenJet>::const_iterator i_lepton = pseudoLeptons->begin(); i_lepton != pseudoLeptons->end(); ++i_lepton){
         //~GenLorentzVector_Jet(&(*i_lepton),tempParticle_.p);
         //~v_allPseudoLepton_.push_back(tempParticle_);
      //~}


   }
   
   ///////////////////////
   // Dilepton Selection//
   ///////////////////////
   // Keep only events which satisfy the pseudoTop selection (recommended by the LHCTopWG plus additonal dileptoMass cut) or satisfy a loose reco dilepton selection
   bool pseudoDileptonSelection=true;
   if (ttbarPseudoDecayMode_==0 || pseudoTopInfo_==0) pseudoDileptonSelection=false;
   bool recoDileptonSelection=true;
      
   ee_=false;
   mumu_=false;
   emu_=false;
   mll_=0;
   lepton1SF_=1.;
   lepton2SF_=1.;
   if ((vElectrons_.size()+vMuons_.size())==NumberLeptons_cut_){
   
      if (vElectrons_.size()==2){
         if (vElectrons_[0].charge*vElectrons_[1].charge!=-1) recoDileptonSelection=false;
         else{
            mll_=(vElectrons_[0].p+vElectrons_[1].p).M();
            ee_=true;
            lepton1SF_=(fctLeptonFullSimScaleFactors_(vElectrons_[0],vElectrons_[0].p.Pt(),vElectrons_[0].etaSC))[0];
            lepton1SF_unc_=(fctLeptonFullSimScaleFactors_(vElectrons_[0],vElectrons_[0].p.Pt(),vElectrons_[0].etaSC))[1];
            lepton2SF_=(fctLeptonFullSimScaleFactors_(vElectrons_[1],vElectrons_[1].p.Pt(),vElectrons_[1].etaSC))[0];
            lepton2SF_unc_=(fctLeptonFullSimScaleFactors_(vElectrons_[1],vElectrons_[1].p.Pt(),vElectrons_[1].etaSC))[1];
         }
      }
      else if (vMuons_.size()==2){
         if (vMuons_[0].charge*vMuons_[1].charge!=-1) recoDileptonSelection=false;
         else{
            mll_=(vMuons_[0].p+vMuons_[1].p).M();
            mumu_=true;
            lepton1SF_=(fctLeptonFullSimScaleFactors_(vMuons_[0],vMuons_[0].p.Pt(),vMuons_[0].p.Eta()))[0];
            lepton1SF_unc_=(fctLeptonFullSimScaleFactors_(vMuons_[0],vMuons_[0].p.Pt(),vMuons_[0].p.Eta()))[1];
            lepton2SF_=(fctLeptonFullSimScaleFactors_(vMuons_[1],vMuons_[1].p.Pt(),vMuons_[1].p.Eta()))[0];
            lepton2SF_unc_=(fctLeptonFullSimScaleFactors_(vMuons_[1],vMuons_[1].p.Pt(),vMuons_[1].p.Eta()))[1];
         }
      }
      else {
         if (vMuons_[0].charge*vElectrons_[0].charge!=-1) recoDileptonSelection=false;
         else{
            mll_=(vMuons_[0].p+vElectrons_[0].p).M();
            emu_=true;
            lepton1SF_=(fctLeptonFullSimScaleFactors_(vElectrons_[0],vElectrons_[0].p.Pt(),vElectrons_[0].etaSC))[0];
            lepton1SF_unc_=(fctLeptonFullSimScaleFactors_(vElectrons_[0],vElectrons_[0].p.Pt(),vElectrons_[0].etaSC))[1];
            lepton2SF_=(fctLeptonFullSimScaleFactors_(vMuons_[0],vMuons_[0].p.Pt(),vMuons_[0].p.Eta()))[0];
            lepton2SF_unc_=(fctLeptonFullSimScaleFactors_(vMuons_[0],vMuons_[0].p.Pt(),vMuons_[0].p.Eta()))[1];
         }
      }
   }
   else recoDileptonSelection=false;
   
   if (!recoDileptonSelection && !pseudoDileptonSelection) return;
   
   hCutFlow_->Fill("Dilepton", mc_weight_*pu_weight_);
   
   addLepton_=false;
   if((vElectrons_add_.size()+vMuons_add_.size())>0) addLepton_=true;
   
   /*
   ///////////
   //Photons//
   ///////////
   edm::Handle<edm::View<pat::Photon>> photonColl;
   iEvent.getByToken(photonCollectionToken_, photonColl);

   vPhotons_.clear();
   tree::Photon trPho;
   for (edm::View<pat::Photon>::const_iterator pho = photonColl->begin(); pho != photonColl->end(); pho++) {
      // Kinematics
      if (pho->pt() < 20) continue;

      trPho.p.SetPtEtaPhiE(pho->pt(), pho->superCluster()->eta(), pho->superCluster()->phi(), pho->energy());
      const edm::Ptr<pat::Photon> phoPtr( photonColl, pho - photonColl->begin() );
      trPho.sigmaIetaIeta = pho->full5x5_sigmaIetaIeta();
      trPho.hOverE = pho->hadTowOverEm();
      trPho.hasPixelSeed = pho->hasPixelSeed();
      trPho.passElectronVeto = pho->passElectronVeto();

      trPho.cIso = pho->chargedHadronIso();
      trPho.nIso = pho->neutralHadronIso();
      trPho.pIso = pho->photonIso();

      // check photon working points
      trPho.isLoose = pho->photonID(photonLooseIdMapToken_);
      trPho.isMedium = pho->photonID(photonMediumIdMapToken_);
      trPho.isTight = pho->photonID(photonTightIdMapToken_);
      // write the photon to collection
      vPhotons_.push_back(trPho);
   } // photon loop

   sort(vPhotons_.begin(), vPhotons_.end(), tree::PtGreater);
   */

   
   /////////
   // Jets//
   /////////
   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
   iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", JetCorParColl);
   JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
   JetCorrectionUncertainty jecUnc(JetCorPar);

   JME::JetResolution resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
   JME::JetResolution resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");
   JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

   edm::Handle<pat::JetCollection> jetColl;
   iEvent.getByToken(jetCollectionToken_, jetColl);

   vJets_.clear();
   tree::Jet trJet;
   Ht_=0;
   for (const pat::Jet& jet : *jetColl) {
      if (fabs(jet.eta())>2.4) continue;
      if (jet.pt()<dJet_pT_cut_) continue;
      //~ trJet.p.SetPtEtaPhi(jet.pt(), jet.eta(), jet.phi());
      trJet.p.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      trJet.bTagCSVv2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      trJet.bTagMVAv2 = jet.bDiscriminator("pfCombinedMVAV2BJetTags");
      trJet.bTagDeepCSV = jet.bDiscriminator("pfDeepCSVJetTags:probb")+jet.bDiscriminator("pfDeepCSVJetTags:probbb");
      trJet.bTagDeepJet = jet.bDiscriminator("pfDeepFlavourJetTags:probb")+jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
      trJet.bTagSoftMuon = (jet.bDiscriminator("softPFMuonBJetTags")<0) ? -1. : jet.bDiscriminator("softPFMuonBJetTags");
      trJet.bTagSoftElectron = (jet.bDiscriminator("softPFElectronBJetTags")<0) ? -1. : jet.bDiscriminator("softPFElectronBJetTags");
      trJet.isLoose = jetIdSelector(jet);
      trJet.TightIDlepVeto = jet.hasUserInt("PFJetIDTightLepVeto") ? jet.userInt("PFJetIDTightLepVeto") : 0;
      jecUnc.setJetEta(jet.eta());
      jecUnc.setJetPt(jet.pt());
      trJet.uncert = jecUnc.getUncertainty(true);
      JME::JetParameters parameters = {{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, rho_}};
      trJet.ptRes = resolution_pt.getResolution(parameters);
      trJet.phiRes = resolution_phi.getResolution(parameters);
      trJet.sfRes = resolution_sf.getScaleFactor(parameters);
      trJet.sfResUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
      trJet.sfResDn = resolution_sf.getScaleFactor(parameters, Variation::DOWN);
      trJet.uncorJecFactor = jet.jecFactor(0);
      trJet.chf = jet.chargedHadronEnergyFraction();
      trJet.nhf = jet.neutralHadronEnergyFraction();
      trJet.cef = jet.chargedEmEnergyFraction();
      trJet.nef = jet.neutralEmEnergyFraction();
      trJet.muonf = jet.muonEnergyFraction();
      trJet.electronf = jet.electronEnergyFraction();
      trJet.nch = jet.chargedMultiplicity();
      trJet.nconstituents = jet.numberOfDaughters();
      trJet.bJetRegressionCorr = jet.hasUserFloat("BJetEnergyCorrFactor") ? jet.userFloat("BJetEnergyCorrFactor") : 1.;
      trJet.bJetRegressionRes = jet.hasUserFloat("BJetEnergyCorrResolution") ? jet.userFloat("BJetEnergyCorrResolution") : 1.;
      trJet.hadronFlavour = jet.hadronFlavour();
      
      // object matching
      trJet.hasElectronMatch = false;
      for (tree::Electron const &el: vElectrons_) {
         if (el.isTight && el.p.Pt()>=20 && trJet.p.DeltaR(el.p)<0.4) {
            trJet.hasElectronMatch = true;
            break;
         }
      }
      trJet.hasMuonMatch = false;
      for (tree::Muon const &mu: vMuons_){
         if (mu.isTight && mu.p.Pt()>=20 && trJet.p.DeltaR(mu.p)<0.4){
            trJet.hasMuonMatch = true;
            break;
         }
      }
      Ht_+=trJet.p.Pt();
            
      vJets_.push_back(trJet);
   } // jet loop
   sort(vJets_.begin(), vJets_.end(), tree::PtGreater);
   
   ///////////////
   // Puppi Jets//
   ///////////////
   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl_Puppi;
   iSetup.get<JetCorrectionsRecord>().get("AK4PFPuppi", JetCorParColl_Puppi);
   JetCorrectorParameters const & JetCorPar_Puppi = (*JetCorParColl_Puppi)["Uncertainty"];
   JetCorrectionUncertainty jecUnc_Puppi(JetCorPar_Puppi);

   JME::JetResolution resolution_pt_Puppi = JME::JetResolution::get(iSetup, "AK4PFPuppi_pt");
   JME::JetResolution resolution_phi_Puppi = JME::JetResolution::get(iSetup, "AK4PFPuppi_phi");
   JME::JetResolutionScaleFactor resolution_sf_Puppi = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFPuppi");

   edm::Handle<pat::JetCollection> jetColl_Puppi;
   iEvent.getByToken(jetPuppiCollectionToken_, jetColl_Puppi);

   vJetsPuppi_.clear();
   tree::Jet trJet_Puppi;
   for (const pat::Jet& jet : *jetColl_Puppi) {
      if (fabs(jet.eta())>2.4) continue;
      if (jet.pt()<dJet_pT_cut_) continue;
      //~ trJet_Puppi.p.SetPtEtaPhi(jet.pt(), jet.eta(), jet.phi());
      trJet_Puppi.p.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      trJet_Puppi.bTagCSVv2 = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      trJet_Puppi.bTagMVAv2 = jet.bDiscriminator("pfCombinedMVAV2BJetTags");
      trJet_Puppi.bTagDeepCSV = jet.bDiscriminator("pfDeepCSVJetTags:probb")+jet.bDiscriminator("pfDeepCSVJetTags:probbb");
      trJet_Puppi.bTagSoftMuon = (jet.bDiscriminator("softPFMuonBJetTags")<0) ? -1. : jet.bDiscriminator("softPFMuonBJetTags");
      trJet_Puppi.bTagSoftElectron = (jet.bDiscriminator("softPFElectronBJetTags")<0) ? -1. : jet.bDiscriminator("softPFElectronBJetTags");
      trJet_Puppi.isLoose = jetIdSelector(jet);
      trJet_Puppi.TightIDlepVeto = jet.hasUserInt("PFJetIDTightLepVeto") ? jet.userInt("PFJetIDTightLepVeto") : 0;
      jecUnc.setJetEta(jet.eta());
      jecUnc.setJetPt(jet.pt());
      trJet_Puppi.uncert = jecUnc.getUncertainty(true);
      JME::JetParameters parameters = {{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, rho_}};
      trJet_Puppi.ptRes = resolution_pt_Puppi.getResolution(parameters);
      trJet_Puppi.phiRes = resolution_phi_Puppi.getResolution(parameters);
      trJet_Puppi.sfRes = resolution_sf_Puppi.getScaleFactor(parameters);
      trJet_Puppi.sfResUp = resolution_sf_Puppi.getScaleFactor(parameters, Variation::UP);
      trJet_Puppi.sfResDn = resolution_sf_Puppi.getScaleFactor(parameters, Variation::DOWN);
      trJet_Puppi.uncorJecFactor = jet.jecFactor(0);
      trJet_Puppi.chf = jet.chargedHadronEnergyFraction();
      trJet_Puppi.nhf = jet.neutralHadronEnergyFraction();
      trJet_Puppi.cef = jet.chargedEmEnergyFraction();
      trJet_Puppi.nef = jet.neutralEmEnergyFraction();
      trJet_Puppi.muonf = jet.muonEnergyFraction();
      trJet_Puppi.electronf = jet.electronEnergyFraction();
      trJet_Puppi.nch = jet.chargedMultiplicity();
      trJet_Puppi.nconstituents = jet.numberOfDaughters();
      trJet_Puppi.bJetRegressionCorr = jet.hasUserFloat("BJetEnergyCorrFactor") ? jet.userFloat("BJetEnergyCorrFactor") : 1.;
      trJet_Puppi.bJetRegressionRes = jet.hasUserFloat("BJetEnergyCorrResolution") ? jet.userFloat("BJetEnergyCorrResolution") : 1.;
      trJet_Puppi.hadronFlavour = jet.hadronFlavour();
      
      // object matching
      trJet_Puppi.hasElectronMatch = false;
      for (tree::Electron const &el: vElectrons_) {
         if (el.isTight && el.p.Pt()>=20 && trJet_Puppi.p.DeltaR(el.p)<0.4) {
            trJet_Puppi.hasElectronMatch = true;
            break;
         }
      }
      trJet_Puppi.hasMuonMatch = false;
      for (tree::Muon const &mu: vMuons_){
         if (mu.isTight && mu.p.Pt()>=20 && trJet_Puppi.p.DeltaR(mu.p)<0.4){
            trJet_Puppi.hasMuonMatch = true;
            break;
         }
      }
      vJetsPuppi_.push_back(trJet_Puppi);
   } // jet loop
   sort(vJetsPuppi_.begin(), vJetsPuppi_.end(), tree::PtGreater);
   
   
   //////////////
   //BTagWeight//
   //////////////
   if (!isRealData) {
      float P_MC = 1.;
      float P_Data = 1.;
      float err1 = 0.;
      float err2 = 0.;
      float err3 = 0.;
      float err4 = 0.;
      
      for(std::vector<tree::Jet>::const_iterator it = vJets_.begin(); it != vJets_.end() ; ++it){
         if (it->p.Pt() >=30.0 && fabs(it->p.Eta())<2.4 && !it->hasMuonMatch && !it->hasElectronMatch && it->TightIDlepVeto){
           
            int jetFlavor;
            float eff;
            float SF, SF_up, SF_down;
            float SF_err;
            float temp_pt = std::min(599.,it->p.Pt());
            float temp_eta = it->p.Eta();
            jetFlavor = it->hadronFlavour;
            if (jetFlavor == 5)
            {
               SF = fctBTagCalibReaderFullSimBJets_.eval(BTagEntry::FLAV_B, temp_eta, temp_pt);
               SF_up = fctBTagCalibReaderFullSimBJetsUp_.eval(BTagEntry::FLAV_B, temp_eta, temp_pt);
               SF_down = fctBTagCalibReaderFullSimBJetsDown_.eval(BTagEntry::FLAV_B, temp_eta, temp_pt);   
            }
            else if (jetFlavor == 4)
            {
               SF = fctBTagCalibReaderFullSimCJets_.eval(BTagEntry::FLAV_C, temp_eta, temp_pt);
               SF_up = fctBTagCalibReaderFullSimCJetsUp_.eval(BTagEntry::FLAV_C, temp_eta, temp_pt);
               SF_down = fctBTagCalibReaderFullSimCJetsDown_.eval(BTagEntry::FLAV_C, temp_eta, temp_pt);
            }
            else
            {
               SF = fctBTagCalibReaderFullSimLightJets_.eval(BTagEntry::FLAV_UDSG, temp_eta, temp_pt);
               SF_up = fctBTagCalibReaderFullSimLightJetsUp_.eval(BTagEntry::FLAV_UDSG, temp_eta, temp_pt);
               SF_down = fctBTagCalibReaderFullSimLightJetsDown_.eval(BTagEntry::FLAV_UDSG, temp_eta, temp_pt);
            }

            eff = fctBTagEff_(jetFlavor, temp_pt, fabs(temp_eta));

            SF_err = std::max(fabs(SF_up-SF),fabs(SF_down-SF));
             
            // check if jet is btagged
            bool tagged = it->bTagDeepCSV>0.2217;
           
            if (tagged){
               P_MC = P_MC * eff;
               P_Data = P_Data * SF * eff;
               if (jetFlavor == 5 || jetFlavor == 4) err1 += SF_err/SF;
               else err3 += SF_err/SF;
            }
            else{
               P_MC = P_MC * (1 - eff);
               P_Data = P_Data  * (1 - SF * eff);
               if (jetFlavor == 5 || jetFlavor == 4) err2 += (-eff*SF_err) / (1-eff*SF);
               else err4 += (-eff*SF_err) / (1-eff*SF);
            }
         }
      }
      if (P_MC > 0. && P_Data > 0.){
         bTagWeight_ = P_Data/P_MC;
         bTagWeightErrHeavy_ = (err1+err2)*P_Data/P_MC;
         bTagWeightErrLight_ = (err3+err4)*P_Data/P_MC;
      }
      else{
         bTagWeight_ = 0.;
         bTagWeightErrHeavy_ = 0.;
         bTagWeightErrLight_ = 0.;
      }
   }
   
   ///////////
   // GenJet//
   ///////////
   edm::Handle<reco::GenJetCollection> genJetColl;
   if (!isRealData) {
     iEvent.getByToken(genJetCollectionToken_, genJetColl);
     vGenJets_.clear();
     tree::Particle trGJet;
     for (const reco::GenJet& jet: *genJetColl) {
        if (jet.pt()<dJet_pT_cut_-5) continue;
        //~ trGJet.p.SetPtEtaPhi(jet.pt(), jet.eta(), jet.phi());
        trGJet.p.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.energy());
        vGenJets_.push_back(trGJet);
     }
     sort(vGenJets_.begin(), vGenJets_.end(), tree::PtGreater);
   } // gen-jet loop

   if (hardPUveto_) {
      for (tree::Jet const &j: vJets_) {
         if (j.isLoose) {
            if (j.p.Pt()>300) return;
            break; // only check first loose jet
         }
      }
   }

   double const HT = computeHT(vJets_);
   if (HT<dHT_cut_) return;
   hCutFlow_->Fill("HT", mc_weight_*pu_weight_);

   ////////
   // MET//
   ////////
   edm::Handle<pat::METCollection> metCollCalo;
   iEvent.getByToken(caloMetCollectionToken_, metCollCalo);
   caloMetPt_ = metCollCalo->front().caloMETPt();

   edm::Handle<pat::METCollection> metColl;
   iEvent.getByToken(metCollectionToken_, metColl);

   const pat::MET &met = metColl->front();
   double metPt = met.pt();
   //~ met_.p.SetPtEtaPhi(metPt, met.eta(), met.phi());
   met_.p.SetPtEtaPhiE(metPt, met.eta(), met.phi(), met.energy());

   if( !isRealData ) {
      const reco::GenMET *genMet = met.genMET();
      //~ met_gen_.p.SetPtEtaPhi(genMet->pt(), genMet->eta(), genMet->phi());
      met_gen_.p.SetPtEtaPhiE(genMet->pt(), genMet->eta(), genMet->phi(), genMet->energy());
   }

   // jet resolution shift is set to 0 for 74X
   met_.uncertainty = 0;
   // loop over all up-shifts save for last one (=NoShift)
   for (uint iShift=0; iShift<(pat::MET::METUncertaintySize-1); iShift+=2) {
      // up and down shifts
      const double u = fabs(met.shiftedPt(pat::MET::METUncertainty(iShift))  -metPt);
      const double d = fabs(met.shiftedPt(pat::MET::METUncertainty(iShift+1))-metPt);
      // average
      const double a = .5*(u+d);
      if (isinf(a)) continue;    //few events with problems for shift in Muon Energy
      // add deviations in quadrature
      met_.uncertainty += a*a;
   }
   met_.uncertainty=TMath::Sqrt(met_.uncertainty);

   pat::MET::LorentzVector metShifted;
   metShifted = met.shiftedP4(pat::MET::NoShift, pat::MET::Raw);
   //~ met_raw_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   met_raw_.p.SetPtEtaPhiE(metShifted.pt(), metShifted.eta(), metShifted.phi(), metShifted.energy());
   
   metShifted = met.shiftedP4(pat::MET::JetEnUp);
   //~ met_JESu_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   met_JESu_.p.SetPtEtaPhiE(metShifted.pt(), metShifted.eta(), metShifted.phi(), metShifted.energy());
   metShifted = met.shiftedP4(pat::MET::JetEnDown);
   //~ met_JESd_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   met_JESd_.p.SetPtEtaPhiE(metShifted.pt(), metShifted.eta(), metShifted.phi(), metShifted.energy());

   metShifted = met.shiftedP4(pat::MET::JetResUp);
   //~ met_JERu_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   met_JERu_.p.SetPtEtaPhiE(metShifted.pt(), metShifted.eta(), metShifted.phi(), metShifted.energy());
   metShifted = met.shiftedP4(pat::MET::JetResDown);
   //~ met_JERd_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   met_JERd_.p.SetPtEtaPhiE(metShifted.pt(), metShifted.eta(), metShifted.phi(), metShifted.energy());
   
   metShifted = met.shiftedP4(pat::MET::NoShift, pat::MET::RawCalo);
   //~ met_JERd_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   metCalo_.p.SetPtEtaPhiE(metShifted.pt(), metShifted.eta(), metShifted.phi(), metShifted.energy());
   
   

   met_.sig = met.metSignificance();
   met_raw_.sig = met_.sig;
   met_JESu_.sig = met_.sig;
   met_JESd_.sig = met_.sig;
   met_JERu_.sig = met_.sig;
   met_JERd_.sig = met_.sig;
   metCalo_.sig = met_.sig;
   
   //PuppiMET
   edm::Handle<pat::METCollection> metCollPuppi;
   iEvent.getByToken(metPuppiCollectionToken_, metCollPuppi);

   const pat::MET &metPuppi = metCollPuppi->front();
   double metPt_Puppi = metPuppi.pt();
   metPuppi_.p.SetPtEtaPhiE(metPt_Puppi, metPuppi.eta(), metPuppi.phi(), metPuppi.energy());
   metPuppi_.sig = metPuppi.metSignificance();
   // loop over all up-shifts save for last one (=NoShift)
   metPuppi_.uncertainty = 0;
   for (uint iShift=0; iShift<(pat::MET::METUncertaintySize-1); iShift+=2) {
      // up and down shifts
      const double u = fabs(metPuppi.shiftedPt(pat::MET::METUncertainty(iShift))  -metPt_Puppi);
      const double d = fabs(metPuppi.shiftedPt(pat::MET::METUncertainty(iShift+1))-metPt_Puppi);
      // average
      const double a = .5*(u+d);
      if (isinf(a)) continue;    //few events with problems for shift in Muon Energy
      if (isnan(a)) continue;    //few events with problems for shift in JetRes
      // add deviations in quadrature
      metPuppi_.uncertainty += a*a;
   }
   metPuppi_.uncertainty=TMath::Sqrt(metPuppi_.uncertainty);
      
   //MET no HF
   edm::Handle<pat::METCollection> metCollNoHF;
   iEvent.getByToken(metNoHFCollectionToken_, metCollNoHF);

   const pat::MET &metNoHF = metCollNoHF->front();
   double metPt_NoHF = metNoHF.pt();
   metNoHF_.p.SetPtEtaPhiE(metPt_NoHF, metNoHF.eta(), metNoHF.phi(), metNoHF.energy());
   metNoHF_.sig = metNoHF.metSignificance();
   
   //DeepMET
   // ~edm::Handle<pat::METCollection> metCollDeep;
   // ~iEvent.getByToken(metDeepCollectionToken_, metCollDeep);
   // ~const pat::MET &metDeep = metCollDeep->front();
   // ~double metPt_Deep = metDeep.pt();
   // ~metDeep_.p.SetPtEtaPhiE(metPt_Deep, metDeep.eta(), metDeep.phi(), metDeep.energy());
   // ~metDeep_.sig = metDeep.metSignificance();
   
   //XY Corrected MET
   std::pair<double,double> MET_XYpair = METXYCorr_Met_MetPhi(metPt, met.phi(), iEvent.run(), 2016, !isRealData, nPV_);
   metXYcorr_.p.SetPtEtaPhiE(MET_XYpair.first, met.eta(), MET_XYpair.second, met.energy());
   metXYcorr_.sig = met.metSignificance();
   metXYcorr_.uncertainty = met_.uncertainty;   //should probably be changed
   
   // ~//MET with BJet Regression Jets
   // ~edm::Handle<pat::METCollection> metCollBJetRegr;
   // ~iEvent.getByToken(metBJetRegressionCollectionToken_, metCollBJetRegr);
   // ~const pat::MET &metBJetRegression = metCollBJetRegr->front();
   // ~double metPt_BJetRegression = metBJetRegression.pt();
   // ~metBJetRegression_.p.SetPtEtaPhiE(metPt_BJetRegression, metBJetRegression.eta(), metBJetRegression.phi(), metBJetRegression.energy());
   // ~metBJetRegression_.sig = metBJetRegression.metSignificance();
   
   // ~//MET with BJet Regression Jets (only applied to jets with a loose selection)
   // ~edm::Handle<pat::METCollection> metCollBJetRegrLoose;
   // ~iEvent.getByToken(metBJetRegressionLooseCollectionToken_, metCollBJetRegrLoose);
   // ~const pat::MET &metBJetRegressionLoose = metCollBJetRegrLoose->front();
   // ~double metPt_BJetRegressionLoose = metBJetRegressionLoose.pt();
   // ~metBJetRegressionLoose_.p.SetPtEtaPhiE(metPt_BJetRegressionLoose, metBJetRegressionLoose.eta(), metBJetRegressionLoose.phi(), metBJetRegressionLoose.energy());
   // ~metBJetRegressionLoose_.sig = metBJetRegressionLoose.metSignificance();
   // ~std::cout<<metBJetRegressionLoose_.p.Pt()<<std::endl;
      
   ///////////
   //MT2//////
   ///////////
   
   if (recoDileptonSelection){
      if (emu_){
         pa[0]=vMuons_[0].p.M(); pa[1]=vMuons_[0].p.Px(); pa[2]=vMuons_[0].p.Py();
         pb[0]=vElectrons_[0].p.M(); pb[1]=vElectrons_[0].p.Px(); pb[2]=vElectrons_[0].p.Py();
      }
      else if (mumu_){
         pa[0]=vMuons_[0].p.M(); pa[1]=vMuons_[0].p.Px(); pa[2]=vMuons_[0].p.Py();
         pb[0]=vMuons_[1].p.M(); pb[1]=vMuons_[1].p.Px(); pb[2]=vMuons_[1].p.Py();
      }
      else {
         pa[0]=vElectrons_[0].p.M(); pa[1]=vElectrons_[0].p.Px(); pa[2]=vElectrons_[0].p.Py();
         pb[0]=vElectrons_[1].p.M(); pb[1]=vElectrons_[1].p.Px(); pb[2]=vElectrons_[1].p.Py();
      }
      pmiss[0]=0; pmiss[1]=met_.p.Px(); pmiss[2]=met_.p.Py();
      
      fctMT2_.set_mn(0.);
      fctMT2_.set_momenta(pa,pb,pmiss);
      
      MT2_=static_cast<float>(fctMT2_.get_mt2());
   }
   else MT2_=0;
   
   
   ///////////////////////////
   // generated HT and DY pt//
   ///////////////////////////
   // copied from https://github.com/Aachen-3A/PxlSkimmer/blob/master/Skimming/src/PxlSkimmer_miniAOD.cc#L590
   genHt_ = -1;
   dyPt_ = -1;
   v_dyLep_.clear();
   tree::GenParticle trDYlep;
   if (!isRealData) {
      edm::Handle<LHEEventProduct> lheInfoHandle;
      iEvent.getByToken(LHEEventToken_, lheInfoHandle);
      if (lheInfoHandle.isValid()) {
         lhef::HEPEUP lheParticleInfo = lheInfoHandle->hepeup();
         // get the five vector
         // (Px, Py, Pz, E and M in GeV)
         std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
         std::vector<int> statusCodes = lheParticleInfo.ISTUP;
         TLorentzVector DYptVec(0., 0., 0., 0.);
         TLorentzVector DYptVec_temp(0., 0., 0., 0.);
         genHt_ = 0;
         for (unsigned int i = 0; i < statusCodes.size(); i++) {
            auto absId = abs(lheParticleInfo.IDUP[i]);
            if (statusCodes[i] == 1 && ( absId < 11 || absId > 16 ) && absId != 22 && !hasAncestor(i, lheParticleInfo, 6)) {
               genHt_ += sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2));
            }
            if(dyPtInfo_ && (absId==11 || absId==13 || absId==15)) { //save gen DY pt
               DYptVec_temp.SetPxPyPzE(allParticles[i][0],allParticles[i][1],allParticles[i][2],allParticles[i][3]);
               DYptVec += DYptVec_temp;
               trDYlep.p = DYptVec_temp;
               trDYlep.pdgId =absId;
               v_dyLep_.push_back(trDYlep);
            }
         } // end particle loop
         if(dyPtInfo_) dyPt_ = DYptVec.Pt();
      } else { // if no lheEventProduct is found
        genHt_ = 0;
        for (const auto& genP : *prunedGenParticles) {
          auto absId = abs(genP.pdgId());
          if (genP.status() == 23 and (absId<11 || absId > 16 ) && genP.pdgId() != 22 ) {
            genHt_ += genP.pt();
          }
        } // genParticle loop
      }
   }

   ////////////////////////
   // Generated Particles//
   ////////////////////////
   // status flags: https://indico.cern.ch/event/459797/contribution/2/attachments/1181555/1710844/mcaod-Nov4-2015.pdf
   vGenParticles_.clear();
   vIntermediateGenParticles_.clear();
   tree::GenParticle trP;
   tree::IntermediateGenParticle trIntermP;
   signal_nBinos_ = 0;
   signal_nNeutralinoDecays_ = 0;
   TVector3 p_EWK_temp;
   TVector3 p_EWK_tot;
   p_EWK_tot.SetPtEtaPhi(0.,0.,0.);
   if (!isRealData) {
      // Get generator level info
      // Pruned particles are the one containing "important" stuff
      for (const reco::GenParticle &genP: *prunedGenParticles){
         auto absId = abs(genP.pdgId());

         if (absId==6||absId==24 || absId==22 || absId==23) { // store intermediate tops and w bosons and DY
            int iNdaugh = genP.numberOfDaughters();
            if (iNdaugh>1) { // skip "decays" V->V
               trIntermP.pdgId = genP.pdgId();
               trIntermP.status = genP.status();
               trIntermP.isPrompt = genP.statusFlags().isPrompt();
               //~ trIntermP.p.SetPtEtaPhi(genP.pt(), genP.eta(), genP.phi());
               trIntermP.p.SetPtEtaPhiE(genP.pt(), genP.eta(), genP.phi(), genP.energy());
               trIntermP.daughters.clear();
               for (int i=0; i<iNdaugh; i++) { // store the decay products
                  reco::Candidate const& daugh = *genP.daughter(i);
                  trP.pdgId = daugh.pdgId();
                  trP.status = daugh.status();
                  trP.isPrompt = false;
                  //~ trP.p.SetPtEtaPhi(daugh.pt(), daugh.eta(), daugh.phi());
                  trP.p.SetPtEtaPhiE(daugh.pt(), daugh.eta(), daugh.phi(), daugh.energy());
                  trIntermP.daughters.push_back(trP);
               }
               vIntermediateGenParticles_.push_back(trIntermP);
            }
         }
         
         // save particles
         if (genP.status()==22 || genP.status()==23 || // some generator particles
               (genP.status() == 1 && absId>100000 && genP.numberOfDaughters()!=1)||    // BSM particles
               (genP.status() == 1 && genP.pt()>10 && (absId==22 || (11 <= absId && absId <= 16)))) { // status 1 photons and leptons (including neutrinos)
            if(genP.numberOfDaughters()==1) continue; //skip particles, which decay like V->V
            trP.pdgId = genP.pdgId();
            trP.status = genP.status();
            trP.isPrompt = genP.statusFlags().isPrompt();
            trP.fromHardProcess = genP.statusFlags().fromHardProcess();
            //~ trP.p.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
            trP.p.SetPtEtaPhiE(genP.pt(), genP.eta(), genP.phi(), genP.energy());
            trP.promptStatus = getPromptStatus(genP, prunedGenParticles);
            vGenParticles_.push_back(trP);
         }
         
         //Calculate EWKinoPairPt
         if (genP.status()==22 && absId>=1000001 && (abs(genP.mother(0)->pdgId())>1000039 || abs(genP.mother(0)->pdgId())<1000001)) {
            p_EWK_temp.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
            p_EWK_tot=p_EWK_tot+p_EWK_temp;
         }
            
      }
      sort(vGenParticles_.begin(), vGenParticles_.end(), tree::PtGreater);
      EWKinoPairPt_ = p_EWK_tot.Pt();
   }
   hCutFlow_->Fill("nBinos", mc_weight_*pu_weight_);

   hCutFlow_->Fill("final", mc_weight_*pu_weight_);
   
   //////////////////////////
   // Gen Information Ttbar//
   //////////////////////////
   topPTweight_ = 1.;
   if(ttbarGenInfo_){
      // Gen-level particles of ttbar system
      edm::Handle<TtGenEvent> ttbarGenEvent;
      iEvent.getByLabel("genEvt", ttbarGenEvent);
      if(!ttbarGenEvent.failedToGet()){
         // Which process generates the ttbar event: gluon-gluon-fusion (0), quark-quark-annihilation (1), all other processes (2)
         if(ttbarGenEvent->fromGluonFusion()) ttbarProductionMode_ = 0;
         else if(ttbarGenEvent->fromQuarkAnnihilation()) ttbarProductionMode_ = 1;
         else ttbarProductionMode_ = 2;
         
         // Top quarks, and for dileptonic decays also decay products
         if(ttbarGenEvent->top()) GenLorentzVector(ttbarGenEvent->top(),genTop_);
         else genTop_ = nullP4_;
         if(ttbarGenEvent->topBar()) GenLorentzVector(ttbarGenEvent->topBar(),genAntiTop_);
         else genAntiTop_ = nullP4_;
         if(ttbarGenEvent->lepton()){
            assignLeptonAndTau(ttbarGenEvent->lepton(), genLepton_, genLeptonPdgId_, genTau_);
         }
         else{
            genLepton_ = nullP4_;
            genLeptonPdgId_ = 0;
            genTau_ = nullP4_;
         }
         if(ttbarGenEvent->leptonBar()){
            assignLeptonAndTau(ttbarGenEvent->leptonBar(), genAntiLepton_, genAntiLeptonPdgId_, genAntiTau_);
         }
         else{
            genAntiLepton_ = nullP4_;
            genAntiLeptonPdgId_ = 0;
            genAntiTau_ = nullP4_;
         }
         if(ttbarGenEvent->b()) GenLorentzVector(ttbarGenEvent->b(),genB_);
         else genB_ = nullP4_;
         if(ttbarGenEvent->bBar()) GenLorentzVector(ttbarGenEvent->bBar(),genAntiB_);
         else genAntiB_ = nullP4_;
         if(ttbarGenEvent->neutrino()) GenLorentzVector(ttbarGenEvent->neutrino(),genNeutrino_);
         else genNeutrino_ = nullP4_;
         if(ttbarGenEvent->neutrinoBar()) GenLorentzVector(ttbarGenEvent->neutrinoBar(),genAntiNeutrino_);
         else genAntiNeutrino_ = nullP4_;
         if(ttbarGenEvent->wPlus()) GenLorentzVector(ttbarGenEvent->wPlus(),genWPlus_);
         else genWPlus_ = nullP4_;
         if(ttbarGenEvent->wMinus()) GenLorentzVector(ttbarGenEvent->wMinus(),genWMinus_);
         else genWMinus_ = nullP4_;
         
         // Top decay mode
         ttbarDecayMode_=0;
         if(ttbarGenEvent->isFullLeptonic(WDecay::kElec,WDecay::kElec)) ttbarDecayMode_=1;
         else if(ttbarGenEvent->isFullLeptonic(WDecay::kMuon,WDecay::kMuon)) ttbarDecayMode_=2;
         else if(ttbarGenEvent->isFullLeptonic(WDecay::kElec,WDecay::kMuon)) ttbarDecayMode_=3;
         else if(ttbarGenEvent->isFullLeptonic(WDecay::kTau,WDecay::kElec)) ttbarDecayMode_=5;
         else if(ttbarGenEvent->isFullLeptonic(WDecay::kTau,WDecay::kMuon)) ttbarDecayMode_=6;
         else if(ttbarGenEvent->isFullLeptonic(WDecay::kTau,WDecay::kTau)) ttbarDecayMode_=7;
         
         // PT reweighting
         topPTweight_ = sqrt(exp(0.0615-0.0005*genTop_.Pt())*exp(0.0615-0.0005*genAntiTop_.Pt()));
         
      }
      else{
         std::cerr<<"\nError: no ttbar gen event?!\n\n";
         genTop_ = nullP4_;
         genAntiTop_ = nullP4_;
         genLepton_ = nullP4_;
         genAntiLepton_ = nullP4_;
         genTau_ = nullP4_;
         genAntiTau_ = nullP4_;
         genLeptonPdgId_ = 0;
         genAntiLeptonPdgId_ = 0;
         genB_ = nullP4_;
         genAntiB_ = nullP4_;
         genNeutrino_ = nullP4_;
         genAntiNeutrino_ = nullP4_;
         genWMinus_ = nullP4_;
         genWPlus_ = nullP4_;
         ttbarDecayMode_=0;
      }
      
   }
   
   //////////////////
   // Gen MT2     //
   /////////////////
   pa[0]=genLepton_.M(); pa[1]=genLepton_.Px(); pa[2]=genLepton_.Py();
   pb[0]=genAntiLepton_.M(); pb[1]=genAntiLepton_.Px(); pb[2]=genAntiLepton_.Py();
   
   pmiss[0]=0; pmiss[1]=met_gen_.p.Px(); pmiss[2]=met_gen_.p.Py();
   
   fctMT2_.set_mn(0.);
   fctMT2_.set_momenta(pa,pb,pmiss);
   
   genMT2_=static_cast<float>(fctMT2_.get_mt2());
      
   if(ttbarGenInfo_){
      genNeutrinoSum_=genNeutrino_+genAntiNeutrino_;
      pmiss[0]=0; pmiss[1]=genNeutrinoSum_.Px(); pmiss[2]=genNeutrinoSum_.Py();
      
      fctMT2_.set_mn(0.);
      fctMT2_.set_momenta(pa,pb,pmiss);
      
      genMT2neutrino_=static_cast<float>(fctMT2_.get_mt2());
   }
   
   ///////////////////
   // Event identity//
   ///////////////////
   evtNo_ = iEvent.id().event();
   runNo_ = iEvent.run();
   lumNo_ = iEvent.luminosityBlock();
   edm::Handle<edm::EDCollection<DetId>> unreplacedGSFixedHandle;
   iEvent.getByLabel("ecalMultiAndGSGlobalRecHitEB", "hitsNotReplaced", unreplacedGSFixedHandle);
   ecalMultiAndGSGlobalRecHitEB_hitsNotReplaced_ = reMiniAOD_ && (!unreplacedGSFixedHandle.isValid() || !unreplacedGSFixedHandle->empty());

   edm::Handle<bool> duplicateGSFixedHandle;
   iEvent.getByLabel("particleFlowEGammaGSFixed", "dupECALClusters", duplicateGSFixedHandle);
   particleFlowEGammaGSFixed_dupECALClusters_ = reMiniAOD_ && *duplicateGSFixedHandle;
   
   //////////////////
   //TrackIsolation//
   //////////////////
   electronTrackIsoVeto = TrackIsolation(packedCandidates, metColl, vertices, 11);
   muonTrackIsoVeto = TrackIsolation(packedCandidates, metColl, vertices, 13);
   pionTrackIsoVeto = TrackIsolation(packedCandidates, metColl, vertices, 211);
   
   // write the event
   eventTree_->Fill();
}

void TreeWriter::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&)
{
   newLumiBlock_ = true;

   edm::Handle<GenLumiInfoHeader> gen_header;
   iLumi.getByLabel("generator", gen_header);
   std::string modelName_ = "";
   signal_m1_ = 0;
   signal_m2_ = 0;
   if (gen_header.isValid()) {
      modelName_ = gen_header->configDescription();
      std::smatch sm;
      if (regex_match(modelName_, sm, std::regex(".*_(\\d+)_(\\d+)"))) {
         signal_m1_ = std::stoi(sm[1]);
         signal_m2_ = std::stoi(sm[2]);
      } else if (regex_match(modelName_, sm, std::regex(".*_(\\d+)"))) {
         signal_m1_ = std::stoi(sm[1]);
      } else if (regex_match(modelName_, sm, std::regex(".*_M1(\\d+)_M3(\\d+)"))) {
         signal_m1_ = std::stoi(sm[1]);
         signal_m2_ = std::stoi(sm[2]);
      } else if (regex_match(modelName_, sm, std::regex(".*_M1(\\d+)_M2(\\d+)"))) {
         signal_m1_ = std::stoi(sm[1]);
         signal_m2_ = std::stoi(sm[2]);
      }
   }

   // create the cutflow histogram for the model if not there yet
   // (modelName="" for non-signal samples)
   if (!hCutFlowMap_.count(modelName_)) {
      hCutFlowMap_[modelName_] = createCutFlowHist(modelName_);
   }
   // point to the right cut flow histogram
   hCutFlow_ = hCutFlowMap_.at(modelName_);
}

void TreeWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(TreeWriter);
