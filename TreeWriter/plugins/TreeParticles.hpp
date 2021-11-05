#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>
#include <array>

enum PhotonMatchType {UNMATCHED = 0,
                      MATCHED_FROM_GUDSCB,
                      MATCHED_FROM_PI0,
                      MATCHED_FROM_OTHER_SOURCES};

enum PromptStatusType {
  DIRECTPROMPT, FRAGMENTPROMPT, LEPTONPROMPT, NOPROMPT
};

namespace tree
{
   struct Particle
   {
      //~ TVector3 p;
      TLorentzVector p;
      bool isEB() { return fabs(p.Eta())<1.4442; }
      bool isEE() {
         auto aEta = fabs(p.Eta());
         return 1.566<aEta && aEta < 2.5;
      }
   };

   struct GenParticle: public Particle
   {
      Int_t pdgId=0;
      Int_t status=0;
      bool isPrompt;
      bool fromHardProcess;
      UChar_t promptStatus;
   };

   struct IntermediateGenParticle: public GenParticle
   {
      std::vector<GenParticle> daughters;
   };


   struct Jet : public Particle
   {
      bool isTight;
      bool TightIDlepVeto;
      bool PileupIDloose;
      // ~bool hasElectronMatch;
      // ~bool hasMuonMatch;
      bool hasElectronMatch_loose;
      bool hasMuonMatch_loose;
      float bTagCSVv2;
      float bTagMVAv2;
      float bTagDeepCSV;
      float bTagDeepJet;
      float bTagSoftMuon;
      float bTagSoftElectron;
      float uncert;
      float chf;
      float nhf;
      float cef;
      float nef;
      float muonf;
      float electronf;
      int nch;
      int nconstituents;
      float uncorJecFactor; // uncorrected jet momentum over corrected jet momentum
      float uncorJecFactor_L1; // L1-corrected jet momentum over corrected jet momentum
      float bJetRegressionCorr;
      float bJetRegressionRes;
      int hadronFlavour; //so far used to derive BTag weight
      TLorentzVector matchedGenJet;
      uint32_t seed;
   };

   struct Muon: public Particle
   {
      Char_t charge; // +/- 1
      bool isTight;
      bool isMedium;
      bool isLoose;
      // PF-based combined relative isolation with Δβ correction:
      // (∑pT(ch.had from PV) + max(0, ∑ET(neut.had) + ∑ET(phot) − 0.5*∑pT(ch.had from PU)))/pT(μ)
      float rIso;
      float d0;
      float dZ;
      float PFminiIso;
      
      float rochesterCorrection;
      std::array<float, 6> corrections;
   };

   struct Electron: public Particle
   {
      Char_t charge; // +/- 1
      bool isLoose;
      bool isMedium;
      bool isTight;
      float rIso;
      float d0;
      float dZ;
      Float_t seedCrystalE;
      float r9;
      float SigmaIEtaIEtaFull5x5;
      float dPhiAtVtx;
      float dEtaAtVtx;
      float HoverE;
      float EoverPInv;
      int MissHits;
      bool ConvVeto;
      float PFminiIso;
      TVector3 pUncorrected;
      float phiSC;     //Supercluster Phi
      float etaSC;     //Supercluster Eta
      
      float corr;
      std::array<float, 7> corrections;
   };
   
   struct Photon : public Particle
   {
      float sigmaIetaIeta; // full 5x5
      float hOverE;
      bool hasPixelSeed;
      bool passElectronVeto;

      float cIso;
      float nIso;
      float pIso;
      
      float corr;

      // IDs
      bool  isLoose;
      bool  isMedium;
      bool  isTight;
   };

   struct MET : public Particle
   {
      Float_t  uncertainty;
      Float_t  sig; // MET significance
   };

   inline bool PtGreater(const tree::Particle p1, const tree::Particle p2) {
      return p1.p.Pt() > p2.p.Pt();
   }
   
   inline bool PtGreaterLorentz(const TLorentzVector p1, const TLorentzVector p2) {
      return p1.Pt() > p2.Pt();
   }

} // end namespace definition
#endif /* TREEPARTICLES_H */
