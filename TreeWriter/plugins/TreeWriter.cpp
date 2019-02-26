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

// taken from https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
// "recipe" https://indico.cern.ch/event/557678/contributions/2247944/attachments/1311994/1963568/16-07-19_ana_manuelf_isr.pdf
int n_isr_jets(edm::Handle<edm::View<reco::GenParticle>> const &genParticles,
                            std::vector<tree::Jet> const &jets) {
   int nisr(0);
   bool matched;
   int momid;
   TVector3 pGen;
   for (tree::Jet const&jet: jets) {
      if (jet.hasMuonMatch || jet.hasElectronMatch || jet.hasPhotonMatch) continue;
      matched = false;
      for (size_t imc(0); imc < genParticles->size(); imc++) {
         if (matched) break;
         const reco::GenParticle &mc = (*genParticles)[imc];
         if (mc.status()!=23 || abs(mc.pdgId())>5) continue;
         momid = abs(mc.mother()->pdgId());
         if (!(momid==6 || momid==23 || momid==24 || momid==25 || momid>1e6)) continue;
         //check against daughter in case of hard initial splitting
         for (size_t idau(0); idau < mc.numberOfDaughters(); idau++) {
            pGen.SetXYZ(mc.daughter(idau)->px(),mc.daughter(idau)->py(),mc.daughter(idau)->pz());
            if (jet.p.DeltaR(pGen)<0.3) {
               matched = true;
               break;
            }
         }
      } // Loop over MC particles
      if (!matched) {
         nisr++;
      }
   } // Loop over jets
   return nisr;
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

template <typename T> int sign(T val) {
   return (T(0) < val) - (val < T(0));
}

TreeWriter::TreeWriter(const edm::ParameterSet& iConfig)
   : dHT_cut_(iConfig.getUntrackedParameter<double>("HT_cut"))
   , dPhoton_pT_cut_(iConfig.getUntrackedParameter<double>("photon_pT_cut"))
   , dJet_pT_cut_(iConfig.getUntrackedParameter<double>("jet_pT_cut"))
   , isolatedPhotons_(iConfig.getUntrackedParameter<bool>("isolatedPhotons"))
   , minNumberPhotons_cut_(iConfig.getUntrackedParameter<unsigned>("minNumberPhotons_cut"))
   , minNumberElectrons_cut_(iConfig.getUntrackedParameter<unsigned>("minNumberElectrons_cut"))
   , minNumberBinos_cut_(iConfig.getUntrackedParameter<unsigned>("minNumberBinos_cut"))
   , newLumiBlock_(true)
   , vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
   , photonCollectionToken_  (consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons")))
   , jetCollectionToken_     (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
   , genJetCollectionToken_  (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")))
   , muonCollectionToken_    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
   , electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons")))
   , metCollectionToken_     (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
   , metCorrectedCollectionToken_  (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metCorrected")))
   , metCalibratedCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metCalibrated")))
   , caloMetCollectionToken_ (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("caloMets")))
   , rhoToken_               (consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
   , ebRecHitsToken_         (consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRecHits")))
   , prunedGenToken_         (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles")))
   , pileUpSummaryToken_     (consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileUpSummary")))
   , LHEEventToken_          (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProduct")))
   , packedCandidateToken_   (consumes<std::vector<pat::PackedCandidate>> (iConfig.getParameter<edm::InputTag>("packedCandidates")))
   // electron id
   , electronVetoIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"   )))
   , electronLooseIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"  )))
   , electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap" )))
   , electronTightIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"  )))
   // photon id
   , photonLooseId15MapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonLooseId15Map"  )))
   , photonMediumId15MapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMediumId15Map" )))
   , photonTightId15MapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonTightId15Map"  )))
   , photonLooseIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonLooseIdMap"  )))
   , photonMediumIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMediumIdMap" )))
   , photonTightIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonTightIdMap"  )))
//   , photonMvaValuesMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("photonMvaValuesMap")))
   , phoLooseIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("photonLooseIdMap" )))
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
   , BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter")))
   , BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter")))
{
   // declare consumptions that are used "byLabel" in analyze()
   mayConsume<GenLumiInfoHeader,edm::InLumi> (edm::InputTag("generator"));
   consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));
   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", ""));
   consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
   consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
   consumes<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotonsBeforeGSFix"));
   consumes<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotons", "", "PAT"));
   consumes<edm::View<pat::Photon>>(edm::InputTag("slimmedPhotons", "", "RECO"));
   consumes<bool>(edm::InputTag("particleFlowEGammaGSFixed", "dupECALClusters"));
   consumes<edm::EDCollection<DetId>>(edm::InputTag("ecalMultiAndGSGlobalRecHitEB", "hitsNotReplaced"));
   consumes<edm::View<pat::Electron>>(edm::InputTag("slimmedElectrons", "", "PAT"));
   consumes<edm::View<pat::Electron>>(edm::InputTag("slimmedElectrons", "", "RECO"));

   eventTree_ = fs_->make<TTree> ("eventTree", "event data");

   eventTree_->Branch("photons", &vPhotons_);
   eventTree_->Branch("jets", &vJets_);
   eventTree_->Branch("genJets", &vGenJets_);
   eventTree_->Branch("electrons", &vElectrons_);
   eventTree_->Branch("muons", &vMuons_);
   eventTree_->Branch("met", &met_);
   eventTree_->Branch("metCorrected", &metCorrected_);
//   eventTree_->Branch("metCalibrated", &metCalibrated_);
   eventTree_->Branch("met_raw", &met_raw_);
   eventTree_->Branch("met_gen", &met_gen_);
   eventTree_->Branch("met_JESu", &met_JESu_);
   eventTree_->Branch("met_JESd", &met_JESd_);
   eventTree_->Branch("met_JERu", &met_JERu_);
   eventTree_->Branch("met_JERd", &met_JERd_);
   eventTree_->Branch("genParticles", &vGenParticles_);
   for (const auto& n : triggerObjectNames_) {
     triggerObjectMap_[n] = std::vector<tree::Particle>();
     eventTree_->Branch(n.c_str(), &triggerObjectMap_[n]);
   }
   eventTree_->Branch("intermediateGenParticles", &vIntermediateGenParticles_);

   //eventTree_->Branch("nPV", &nPV_, "nPV/I");
   eventTree_->Branch("true_nPV", &true_nPV_, "true_nPV/I");
   eventTree_->Branch("nGoodVertices" , &nGoodVertices_ , "nGoodVertices/I");
   eventTree_->Branch("nTracksPV", &nTracksPV_, "nTracksPV/I");
   eventTree_->Branch("rho", &rho_, "rho/F");
   eventTree_->Branch("caloMetPt", &caloMetPt_, "caloMetPt/F");

   eventTree_->Branch("pu_weight", &pu_weight_, "pu_weight/F");
   eventTree_->Branch("mc_weight", &mc_weight_, "mc_weight/B");
   eventTree_->Branch("pdf_weights", &vPdf_weights_);

   eventTree_->Branch("genHt", &genHt_, "genHt/F");
   eventTree_->Branch("nISR", &nISR_, "nISR/I");
   //eventTree_->Branch("puPtHat", &puPtHat_ , "puPtHat/F");
   eventTree_->Branch("EWKinoPairPt", &EWKinoPairPt_, "EWKinoPairPt/F");

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
      "photons",
      "HT",
      "nBinos",
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

   // PileUp weights
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

   // generator weights
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

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   edm::InputTag triggerTag("TriggerResults", "", "HLT");
   edm::InputTag triggerPrescaleTag("patTrigger");
   iEvent.getByLabel(triggerTag, triggerBits);
   iEvent.getByLabel(triggerPrescaleTag, triggerPrescales);

   // for each lumiBlock, re-read the trigger indices (rather changes for new run)
   if (triggerIndex_.size() && newLumiBlock_) {
      newLumiBlock_ = false;
      // set all trigger indeces to -1 as "not available"-flag
      for (auto& it: triggerIndex_) { it.second = -1; }
      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
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
   edm::InputTag triggerObjects_("selectedPatTrigger");
   iEvent.getByLabel(triggerObjects_, triggerObjects);

   for (const auto& n : triggerObjectNames_) triggerObjectMap_.at(n).clear();
   tree::Particle trObj;
   for (const pat::TriggerObjectStandAlone obj: *triggerObjects) {
      for (const auto& n : triggerObjectNames_) {
         if (std::count(obj.filterLabels().begin(), obj.filterLabels().end(), n)) {
            trObj.p.SetPtEtaPhi(obj.pt(), obj.eta(), obj.phi());
            triggerObjectMap_.at(n).push_back(trObj);
         }
      }
   }

   // MET Filters
   edm::Handle<edm::TriggerResults> metFilterBits;
   edm::InputTag metFilterTag("TriggerResults", "");
   iEvent.getByLabel(metFilterTag, metFilterBits);
   // go through the filters and check if they were passed
   const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
   //for (unsigned i=0; i<allFilterNames.size(); i++) std::cout << allFilterNames.triggerName(i) << std::endl;
   for (std::string const &name: metFilterNames_) {
      const unsigned index = allFilterNames.triggerIndex(name);
      if (index >= allFilterNames.size()) std::cerr << "MET filter '" << name << "' not found!" << std::endl;
      if (!metFilterBits->accept(index)) return; // not passed
   }
   if (!reMiniAOD_) {
       edm::Handle<bool> ifilterbadChCand;
       iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
       if (!*ifilterbadChCand) return;

       edm::Handle<bool> ifilterbadPFMuon;
       iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
       if (!*ifilterbadPFMuon) return;
   }
   hCutFlow_->Fill("METfilters", mc_weight_*pu_weight_);

   // Get PV
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
   
   // number of tracks
   edm::Handle<std::vector<pat::PackedCandidate>> packedCandidates;
   iEvent.getByToken(packedCandidateToken_, packedCandidates);
   nTracksPV_ = std::count_if(packedCandidates->begin(),packedCandidates->end(), [] (const pat::PackedCandidate& cand) {
      return cand.pt()>.9 && cand.charge() && cand.pvAssociationQuality() == pat::PackedCandidate::UsedInFitTight && cand.fromPV() == pat::PackedCandidate::PVUsedInFit;});

   // get gen particles before photons for the truth match
   edm::Handle<edm::View<reco::GenParticle>> prunedGenParticles;
   if (!isRealData) { iEvent.getByToken(prunedGenToken_,prunedGenParticles); }

   edm::Handle<edm::ValueMap<float>> phoWorstChargedIsolationMap;
   iEvent.getByToken(phoWorstChargedIsolationToken_, phoWorstChargedIsolationMap);

   edm::Handle<edm::ValueMap<bool>> loose_id15_dec;
   edm::Handle<edm::ValueMap<bool>> medium_id15_dec;
   edm::Handle<edm::ValueMap<bool>> tight_id15_dec;
   edm::Handle<edm::ValueMap<bool>> loose_id_dec;
   edm::Handle<edm::ValueMap<bool>> medium_id_dec;
   edm::Handle<edm::ValueMap<bool>> tight_id_dec;
//   edm::Handle<edm::ValueMap<float>> mva_value;
   edm::Handle<edm::ValueMap<vid::CutFlowResult>> loose_id_cutflow;
   iEvent.getByToken(photonLooseId15MapToken_, loose_id15_dec);
   iEvent.getByToken(photonMediumId15MapToken_, medium_id15_dec);
   iEvent.getByToken(photonTightId15MapToken_, tight_id15_dec);
   iEvent.getByToken(photonLooseIdMapToken_, loose_id_dec);
   iEvent.getByToken(photonMediumIdMapToken_, medium_id_dec);
   iEvent.getByToken(photonTightIdMapToken_, tight_id_dec);

//   iEvent.getByToken(photonMvaValuesMapToken_,mva_value);
   iEvent.getByToken(phoLooseIdFullInfoMapToken_, loose_id_cutflow);

   edm::Handle<EcalRecHitCollection> ebRecHits;
   iEvent.getByToken(ebRecHitsToken_, ebRecHits);

   // old photon collection
   edm::Handle<edm::View<pat::Photon>> photonCollOld;
   iEvent.getByLabel("slimmedPhotonsBeforeGSFix", photonCollOld);

   // uncorrected photon collection
   edm::Handle<edm::View<pat::Photon>> photonCollUncorrected;
   if (reMiniAOD_ || !isRealData) iEvent.getByLabel(edm::InputTag("slimmedPhotons", "", "PAT"), photonCollUncorrected);
   else iEvent.getByLabel(edm::InputTag("slimmedPhotons", "", "RECO"), photonCollUncorrected);

   // photon collection
   edm::Handle<edm::View<pat::Photon>> photonColl;
   iEvent.getByToken(photonCollectionToken_, photonColl);

   vPhotons_.clear();
   tree::Photon trPho;
   for (edm::View<pat::Photon>::const_iterator pho = photonColl->begin(); pho != photonColl->end(); pho++) {
      // Kinematics
      if (pho->pt() < 15) continue;

      trPho.p.SetPtEtaPhi(pho->pt(), pho->superCluster()->eta(), pho->superCluster()->phi());

      trPho.seedCrystalE = seedCrystalEnergyEB(*pho->superCluster(), ebRecHits);
      const edm::Ptr<pat::Photon> phoPtr( photonColl, pho - photonColl->begin() );
      trPho.sigmaPt = pho->getCorrectedEnergyError(pho->getCandidateP4type())*sin(trPho.p.Theta());
      trPho.sigmaIetaIeta = pho->full5x5_sigmaIetaIeta();
      trPho.sigmaIphiIphi = pho->full5x5_showerShapeVariables().sigmaIphiIphi;
      trPho.hOverE = pho->hadTowOverEm();
      trPho.hasPixelSeed = pho->hasPixelSeed();
      trPho.passElectronVeto = pho->passElectronVeto();
      trPho.r9 = pho->r9();
      trPho.hasGainSwitch = !pho->hasUserInt("hasGainSwitchFlag") || pho->userInt("hasGainSwitchFlag");

      auto itPos = std::distance(photonColl->begin(), pho);
      trPho.pMultifit.SetXYZ(0,0,0);
      trPho.pUncorrected.SetXYZ(0,0,0);
      if (reMiniAOD_ && itPos<photonCollOld->size()) {
        auto p = (*photonCollOld).at(itPos);
        trPho.pMultifit.SetPtEtaPhi(p.pt(), p.superCluster()->eta(), p.superCluster()->phi());
      }
      if (itPos<photonCollUncorrected->size()) {
        auto p = (*photonCollUncorrected).at(itPos);
        trPho.pUncorrected.SetPtEtaPhi(p.pt(), p.superCluster()->eta(), p.superCluster()->phi());
      }

      vid::CutFlowResult cutFlow = (*loose_id_cutflow)[phoPtr];
      trPho.cIso = cutFlow.getValueCutUpon(4);
      trPho.nIso = cutFlow.getValueCutUpon(5);
      trPho.pIso = cutFlow.getValueCutUpon(6);
      trPho.cIsoWorst = (*phoWorstChargedIsolationMap)[phoPtr];

//      trPho.mvaValue=(*mva_value)[phoPtr];

      // MC match
      if (!isRealData) {
         trPho.isTrue = matchToTruth(*pho, prunedGenParticles);
         trPho.isTrueAlternative = matchToTruthAlternative(*pho, prunedGenParticles);
      } else {
         trPho.isTrue = UNMATCHED;
         trPho.isTrueAlternative = UNMATCHED;
      }

      // check photon working points
      trPho.isLoose15 = (*loose_id15_dec) [phoPtr];
      trPho.isMedium15 = (*medium_id15_dec)[phoPtr];
      trPho.isTight15 = (*tight_id15_dec) [phoPtr];
      trPho.isLoose = (*loose_id_dec) [phoPtr];
      trPho.isMedium = (*medium_id_dec)[phoPtr];
      trPho.isTight = (*tight_id_dec) [phoPtr];
      // write the photon to collection
      if (isolatedPhotons_ && !trPho.isLoose15 && !trPho.isLoose) continue;
      vPhotons_.push_back(trPho);
   } // photon loop

   sort(vPhotons_.begin(), vPhotons_.end(), tree::PtGreater);
   if (minNumberPhotons_cut_ && (vPhotons_.size()<minNumberPhotons_cut_|| vPhotons_.at(0).p.Pt()<dPhoton_pT_cut_)) return;
   hCutFlow_->Fill("photons", mc_weight_*pu_weight_);

   // Muons
   edm::Handle<pat::MuonCollection> muonColl;
   iEvent.getByToken(muonCollectionToken_, muonColl);

   vMuons_.clear();
   tree::Muon trMuon;
   for (const pat::Muon &mu : *muonColl) {
      //~ if (!mu.isLooseMuon()) continue;
      if (! (mu.isPFMuon() || mu.isGlobalMuon() || mu.isTrackerMuon())) continue;
      if (mu.pt()<3) continue;
      trMuon.p.SetPtEtaPhi(mu.pt(), mu.eta(), mu.phi());
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
      vMuons_.push_back(trMuon);
   } // muon loop
   sort(vMuons_.begin(), vMuons_.end(), tree::PtGreater);

   // Electrons
   // Get the electron ID data from the event stream
   edm::Handle<edm::ValueMap<bool>> veto_id_decisions;
   edm::Handle<edm::ValueMap<bool>> loose_id_decisions;
   edm::Handle<edm::ValueMap<bool>> medium_id_decisions;
   edm::Handle<edm::ValueMap<bool>> tight_id_decisions;
   iEvent.getByToken(electronVetoIdMapToken_, veto_id_decisions);
   iEvent.getByToken(electronLooseIdMapToken_, loose_id_decisions);
   iEvent.getByToken(electronMediumIdMapToken_, medium_id_decisions);
   iEvent.getByToken(electronTightIdMapToken_, tight_id_decisions);
   
   //Uncorrected electron collection
   edm::Handle<edm::View<pat::Electron>> electronCollUncorrected;
   if (reMiniAOD_ || !isRealData) iEvent.getByLabel(edm::InputTag("slimmedElectrons", "", "PAT"), electronCollUncorrected);
   else iEvent.getByLabel(edm::InputTag("slimmedElectrons", "", "RECO"), electronCollUncorrected);
   
   //Electron collection
   edm::Handle<edm::View<pat::Electron>> electronColl;
   iEvent.getByToken(electronCollectionToken_, electronColl);

   vElectrons_.clear();
   tree::Electron trEl;
   for (edm::View<pat::Electron>::const_iterator el = electronColl->begin();el != electronColl->end(); el++) {
      //~ if (el->pt()<5) continue;
      const edm::Ptr<pat::Electron> elPtr(electronColl, el - electronColl->begin());
      //~ if (!(*veto_id_decisions)[elPtr]) continue; // take only 'veto' electrons
      trEl.isLoose =(*loose_id_decisions) [elPtr];
      trEl.isMedium=(*medium_id_decisions) [elPtr];
      trEl.isTight =(*tight_id_decisions) [elPtr];
      trEl.p.SetPtEtaPhi(el->pt(), el->superCluster()->eta(), el->superCluster()->phi());
      trEl.seedCrystalE = seedCrystalEnergyEB(*el->superCluster(), ebRecHits);
      trEl.charge = el->charge();
      auto const & pfIso = el->pfIsolationVariables();
      trEl.rIso = (pfIso.sumChargedHadronPt + std::max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))/el->pt();
      trEl.r9 = el->r9();
      trEl.SigmaIEtaIEtaFull5x5 = el->full5x5_sigmaIetaIeta();
      trEl.dPhiAtVtx = el->deltaPhiSuperClusterTrackAtVtx();
      trEl.dEtaAtVtx = el->deltaEtaSuperClusterTrackAtVtx();
      trEl.HoverE = el->hcalOverEcal();
      trEl.MissHits = el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
      trEl.ConvVeto = el->passConversionVeto();
      trEl.PFminiIso = getPFIsolation(packedCandidates, *el);
      math::XYZPoint vtx_point = firstGoodVertex.position();
      trEl.d0 = el->bestTrack()->dxy( vtx_point );
      trEl.dZ = el->bestTrack()->dz( vtx_point );
      trEl.phiObj = el->phi();
      trEl.etaObj = el->eta();
      
      // VID calculation of (1/E - 1/p)
      if (el->ecalEnergy() == 0)   trEl.EoverPInv = 1e30;
      else if (!std::isfinite(el->ecalEnergy()))  trEl.EoverPInv = 1e30;
      else trEl.EoverPInv = (1.0 - el->eSuperClusterOverP())/el->ecalEnergy();
      
      //uncorrected electrons
      auto itPos = std::distance(electronColl->begin(), el);
      trEl.pUncorrected.SetXYZ(0,0,0);
      if (itPos<electronCollUncorrected->size()) {
        auto ele = (*electronCollUncorrected).at(itPos);
        trEl.pUncorrected.SetPtEtaPhi(ele.pt(), ele.superCluster()->eta(), ele.superCluster()->phi());
      }
		
      vElectrons_.push_back(trEl);
   }
   sort(vElectrons_.begin(), vElectrons_.end(), tree::PtGreater);
   if (vElectrons_.size()<minNumberElectrons_cut_) return;

   // Jets
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
   for (const pat::Jet& jet : *jetColl) {
      if (fabs(jet.eta())>3) continue;
      if (jet.pt()<dJet_pT_cut_) continue;
      trJet.p.SetPtEtaPhi(jet.pt(), jet.eta(), jet.phi());
      trJet.bDiscriminator = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      trJet.isLoose = jetIdSelector(jet);
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
      trJet.nch = jet.chargedMultiplicity();
      trJet.nconstituents = jet.numberOfDaughters();
      // object matching
      trJet.hasPhotonMatch = false;
      for (tree::Photon const &ph: vPhotons_) {
         if (ph.isLoose && trJet.p.DeltaR(ph.p)<0.4) {
            trJet.hasPhotonMatch = true;
            break;
         }
      }
      trJet.hasElectronMatch = false;
      for (tree::Electron const &el: vElectrons_) {
         if (el.isLoose && el.p.Pt()>=5 && trJet.p.DeltaR(el.p)<0.4) {
            trJet.hasElectronMatch = true;
            break;
         }
      }
      trJet.hasMuonMatch = false;
      for (tree::Muon const &mu: vMuons_){
         if (mu.isLoose && mu.p.Pt()>=5 && trJet.p.DeltaR(mu.p)<0.4){
            trJet.hasMuonMatch = true;
            break;
         }
      }
      vJets_.push_back(trJet);
   } // jet loop
   sort(vJets_.begin(), vJets_.end(), tree::PtGreater);

   // number of ISR jets
   nISR_ = isRealData? 0 : n_isr_jets(prunedGenParticles, vJets_);

   edm::Handle<reco::GenJetCollection> genJetColl;
   if (!isRealData) {
     iEvent.getByToken(genJetCollectionToken_, genJetColl);
     vGenJets_.clear();
     tree::Particle trGJet;
     for (const reco::GenJet& jet: *genJetColl) {
        if (jet.pt()<dJet_pT_cut_-5) continue;
        trGJet.p.SetPtEtaPhi(jet.pt(), jet.eta(), jet.phi());
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

   // MET
   edm::Handle<pat::METCollection> metCollCalo;
   iEvent.getByToken(caloMetCollectionToken_, metCollCalo);
   caloMetPt_ = metCollCalo->front().caloMETPt();

   edm::Handle<pat::METCollection> metCorrectedColl;
   iEvent.getByToken(metCorrectedCollectionToken_, metCorrectedColl);
   metCorrected_.p.SetPtEtaPhi(metCorrectedColl->front().pt(), metCorrectedColl->front().eta(), metCorrectedColl->front().phi());

   edm::Handle<pat::METCollection> metCalibratedColl;
   iEvent.getByToken(metCalibratedCollectionToken_, metCalibratedColl);
   metCalibrated_.p.SetPtEtaPhi(metCalibratedColl->front().pt(), metCalibratedColl->front().eta(), metCalibratedColl->front().phi());

   edm::Handle<pat::METCollection> metColl;
   iEvent.getByToken(metCollectionToken_, metColl);

   const pat::MET &met = metColl->front();
   double metPt = met.pt();
   met_.p.SetPtEtaPhi(metPt, met.eta(), met.phi());

   if( !isRealData ) {
      const reco::GenMET *genMet = met.genMET();
      met_gen_.p.SetPtEtaPhi(genMet->pt(), genMet->eta(), genMet->phi());
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
      // add deviations in quadrature
      met_.uncertainty += a*a;
   }
   met_.uncertainty=TMath::Sqrt(met_.uncertainty);

   pat::MET::LorentzVector metShifted;
   metShifted = met.shiftedP4(pat::MET::NoShift, pat::MET::Raw);
   met_raw_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());

   metShifted = met.shiftedP4(pat::MET::JetEnUp);
   met_JESu_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   metShifted = met.shiftedP4(pat::MET::JetEnDown);
   met_JESd_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());

   metShifted = met.shiftedP4(pat::MET::JetResUp);
   met_JERu_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());
   metShifted = met.shiftedP4(pat::MET::JetResDown);
   met_JERd_.p.SetPtEtaPhi(metShifted.pt(), metShifted.eta(), metShifted.phi());

   met_.sig = met.metSignificance();
   met_raw_.sig = met_.sig;
   met_JESu_.sig = met_.sig;
   met_JESd_.sig = met_.sig;
   met_JERu_.sig = met_.sig;
   met_JERd_.sig = met_.sig;

   // generated HT
   // copied from https://github.com/Aachen-3A/PxlSkimmer/blob/master/Skimming/src/PxlSkimmer_miniAOD.cc#L590
   genHt_ = -1;
   if (!isRealData) {
      edm::Handle<LHEEventProduct> lheInfoHandle;
      iEvent.getByToken(LHEEventToken_, lheInfoHandle);
      if (lheInfoHandle.isValid()) {
         lhef::HEPEUP lheParticleInfo = lheInfoHandle->hepeup();
         // get the five vector
         // (Px, Py, Pz, E and M in GeV)
         std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
         std::vector<int> statusCodes = lheParticleInfo.ISTUP;
         genHt_ = 0;
         for (unsigned int i = 0; i < statusCodes.size(); i++) {
            auto absId = abs(lheParticleInfo.IDUP[i]);
            if (statusCodes[i] == 1 && ( absId < 11 || absId > 16 ) && absId != 22 && !hasAncestor(i, lheParticleInfo, 6)) {
               genHt_ += sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2));
            }
         } // end paricle loop
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


   // Generated Particles
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
      //~ std::cout<<"------------------------------"<<std::endl;
      // Get generator level info
      // Pruned particles are the one containing "important" stuff
      for (const reco::GenParticle &genP: *prunedGenParticles){
         auto absId = abs(genP.pdgId());
         // estimate number of binos
         if (absId == 1000023 && abs(genP.mother(0)->pdgId()) != 1000023) signal_nBinos_++;
         
         // number of neutralino_1 decays
         if (absId == 1000022 && genP.status()==22) signal_nNeutralinoDecays_++;

         if (absId==1000022||absId==1000023||absId==1000025||absId==1000024) { // store intermediate neutralinos and charginos
            int iNdaugh = genP.numberOfDaughters();
            if (iNdaugh>1) { // skip "decays" V->V
               trIntermP.pdgId = genP.pdgId();
               trIntermP.isPrompt = genP.statusFlags().isPrompt();
               trIntermP.p.SetPtEtaPhi(genP.pt(), genP.eta(), genP.phi());
               trIntermP.daughters.clear();
               //~ std::cout<<absId<<std::endl;
               for (int i=0; i<iNdaugh; i++) { // store the decay products
                  reco::Candidate const& daugh = *genP.daughter(i);
                  trP.pdgId = daugh.pdgId();
                  trP.isPrompt = false;
                  trP.p.SetPtEtaPhi(daugh.pt(), daugh.eta(), daugh.phi());
                  trIntermP.daughters.push_back(trP);
                  //~ std::cout<<daugh.pdgId()<<std::endl;
               }
               vIntermediateGenParticles_.push_back(trIntermP);
            }
         }
         
         // save particles
         if (genP.status()==22 || genP.status()==23 || // some generator particles
               (genP.status() == 1 && genP.pt()>20 && (absId==22 || (11 <= absId && absId <= 16)))) { // status 1 photons and leptons (including neutrinos)
            trP.pdgId = genP.pdgId();
            trP.isPrompt = genP.statusFlags().isPrompt();
            trP.fromHardProcess = genP.statusFlags().fromHardProcess();
            trP.p.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
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
   if (signal_nBinos_ < minNumberBinos_cut_) return;
   hCutFlow_->Fill("nBinos", mc_weight_*pu_weight_);

   hCutFlow_->Fill("final", mc_weight_*pu_weight_);
   // store event identity
   evtNo_ = iEvent.id().event();
   runNo_ = iEvent.run();
   lumNo_ = iEvent.luminosityBlock();

   edm::Handle<edm::EDCollection<DetId>> unreplacedGSFixedHandle;
   iEvent.getByLabel("ecalMultiAndGSGlobalRecHitEB", "hitsNotReplaced", unreplacedGSFixedHandle);
   ecalMultiAndGSGlobalRecHitEB_hitsNotReplaced_ = reMiniAOD_ && (!unreplacedGSFixedHandle.isValid() || !unreplacedGSFixedHandle->empty());

   edm::Handle<bool> duplicateGSFixedHandle;
   iEvent.getByLabel("particleFlowEGammaGSFixed", "dupECALClusters", duplicateGSFixedHandle);
   particleFlowEGammaGSFixed_dupECALClusters_ = reMiniAOD_ && *duplicateGSFixedHandle;
   
   //TrackIsolation
   electronTrackIsoVeto = TrackIsolation(packedCandidates, metColl, vertices, 11);
   muonTrackIsoVeto = TrackIsolation(packedCandidates, metColl, vertices, 13);
   pionTrackIsoVeto = TrackIsolation(packedCandidates, metColl, vertices, 211);
   
   // write the event
   eventTree_->Fill();
}

void TreeWriter::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{
   // use this to print the weight indices that are used for muR, muF and PDF variations
   // see https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
   // Please also add "consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"));"
   // in the constructor
   /*
   try {
      edm::Handle<LHERunInfoProduct> run;
      typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

      iRun.getByLabel("externalLHEProducer", run);
      LHERunInfoProduct myLHERunInfoProduct = *(run.product());

      for (headers_const_iterator iter = myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++) {
         if (iter->tag().find("initrwgt")==std::string::npos) continue;
         std::cout << iter->tag() << std::endl;
         std::vector<std::string> lines = iter->lines();
         for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
            std::cout << lines.at(iLine);
         }
      }
   } catch (std::exception &e) {
      std::cout<<"cannot read scale/pdf information from lhe:"<<std::endl;
      std::cout<<e.what()<<std::endl;
   }
   */
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

int TreeWriter::matchToTruth(const pat::Photon &pho, const edm::Handle<edm::View<reco::GenParticle>> &genParticles)
{
   //
   // Explicit loop and geometric matching method
   //

   // Find the closest status 1 gen photon to the reco photon
   double dR = 999;
   const reco::Candidate *closestPhoton = 0;
   for (auto const& particle: *genParticles) {
      // Drop everything that is not photon or not status 1
      if (abs(particle.pdgId()) != 22 || particle.status() != 1) continue;

      double dRtmp = ROOT::Math::VectorUtil::DeltaR(pho.p4(), particle.p4());
      if (dRtmp < dR) {
         dR = dRtmp;
         closestPhoton = &particle;
      }
   }
   // See if the closest photon (if it exists) is close enough.
   // If not, no match found.
   if (!(closestPhoton != 0 && dR < 0.1)) {
      return UNMATCHED;
   }

   // Find ID of the parent of the found generator level photon match
   int ancestorPID = -999;
   int ancestorStatus = -999;
   findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

   // Allowed parens: quarks pdgId 1-5, or a gluon 21
   std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
   if (!(std::find(allowedParents.begin(),
                   allowedParents.end(), ancestorPID)
         != allowedParents.end())) {
      // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not.
      if (abs(ancestorPID) == 111)
         return MATCHED_FROM_PI0;
      else
         return MATCHED_FROM_OTHER_SOURCES;
   }
   return MATCHED_FROM_GUDSCB;
}

void TreeWriter::findFirstNonPhotonMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus)
{
   if (particle == 0) {
      printf("TreeWriter: ERROR! null candidate pointer, this should never happen\n");
      return;
   }

   // Is this the first non-photon parent? If yes, return, otherwise
   // go deeper into recursion
   if (abs(particle->pdgId()) == 22) {
      findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
   } else {
      ancestorPID = particle->pdgId();
      ancestorStatus = particle->status();
   }

   return;
}

int TreeWriter::matchToTruthAlternative(const pat::Photon &pho, const edm::Handle<edm::View<reco::GenParticle>> &genParticles)
{
   // Explicit loop and geometric matching method
   int isMatched = UNMATCHED;
   for (auto const& particle: *genParticles) {
      int pid = particle.pdgId();
      int ancestorPID = -999;
      int ancestorStatus = -999;
      findFirstNonPhotonMother(&particle, ancestorPID, ancestorStatus);
      if (pid ==22 && TMath::Abs(ancestorPID) <= 22) {
         double dr = ROOT::Math::VectorUtil::DeltaR(pho.p4(), particle.p4());
         float dpt = fabs( (pho.pt() - particle.pt() )/particle.pt());
         if (dr < 0.2 && dpt < 0.2) {
            isMatched = MATCHED_FROM_GUDSCB;
            if(ancestorPID == 22) {
               printf("Ancestor of a photon is a photon!\n");
            }
         }
      }
   }
   return isMatched;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeWriter);
