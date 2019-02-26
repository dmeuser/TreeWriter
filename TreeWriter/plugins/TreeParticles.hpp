#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>

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
      TVector3 p;
      bool isEB() { return fabs(p.Eta())<1.4442; }
      bool isEE() {
         auto aEta = fabs(p.Eta());
         return 1.566<aEta && aEta < 2.5;
      }
   };

   struct GenParticle: public Particle
   {
      Int_t pdgId=0;
      bool isPrompt;
      bool fromHardProcess;
      UChar_t promptStatus;
   };

   struct IntermediateGenParticle: public GenParticle
   {
      std::vector<GenParticle> daughters;
   };

   struct Photon : public Particle
   {
      Float_t sigmaIetaIeta; // full 5x5
      Float_t sigmaIphiIphi;
      Float_t hOverE;
      Bool_t hasPixelSeed;
      Bool_t passElectronVeto;
      Float_t r9;
      Float_t sigmaPt;
      Float_t seedCrystalE;
      Bool_t hasGainSwitch;

      Float_t cIso;
      Float_t nIso;
      Float_t pIso;
      Float_t cIsoWorst;

      Int_t isTrue;
      Int_t isTrueAlternative;
      TVector3 pMultifit;
      TVector3 pUncorrected;

      // IDs
      Bool_t  isLoose15;
      Bool_t  isMedium15;
      Bool_t  isTight15;
      Bool_t  isLoose;
      Bool_t  isMedium;
      Bool_t  isTight;
      Float_t mvaValue;
   };

   struct Jet : public Particle
   {
      bool isLoose;
      bool hasPhotonMatch;
      bool hasElectronMatch;
      bool hasMuonMatch;
      float bDiscriminator;
      float uncert;
      float chf;
      float nhf;
      float cef;
      float nef;
      int nch;
      int nconstituents;
      float ptRes;
      float phiRes;
      float sfRes;
      float sfResUp;
      float sfResDn;
      float uncorJecFactor; // uncorrected jet momentum over corrected jet momentum
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
      float phiObj;     //Since in nominal p the supercluster phi is used
      float etaObj;     //Since in nominal p the supercluster eta is used
   };

   struct MET : public Particle
   {
      Float_t  uncertainty;
      Float_t  sig; // MET significance
   };

   inline bool PtGreater(const tree::Particle p1, const tree::Particle p2) {
      return p1.p.Pt() > p2.p.Pt();
   }

} // end namespace definition
#endif /* TREEPARTICLES_H */
