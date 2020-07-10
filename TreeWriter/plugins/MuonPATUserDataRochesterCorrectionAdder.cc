#ifndef TopAnalysis_ZTopUtils_MuonPATUserDataRochesterCorrectionAdder
#define TopAnalysis_ZTopUtils_MuonPATUserDataRochesterCorrectionAdder

#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <string>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "consumeTemplate.h"

#include "../roccor/RoccoR.cc"

#include "TRandom3.h"

class MuonPATUserDataRochesterCorrectionAdder : public ProducerTemplate {

 public:

  explicit MuonPATUserDataRochesterCorrectionAdder(const edm::ParameterSet&);
  ~MuonPATUserDataRochesterCorrectionAdder() {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 protected:

  virtual void produce(edm::Event&, const edm::EventSetup&);

  const edm::InputTag muons_;

  const bool applyEnergyCorrections_;

  const std::string path_;

  const bool debug_;

  RoccoR muonRochCorr_;

  TRandom3 rand_;
};

MuonPATUserDataRochesterCorrectionAdder::MuonPATUserDataRochesterCorrectionAdder(const edm::ParameterSet& cfg)
 : muons_(cfg.getParameter<edm::InputTag>("src"))
 , applyEnergyCorrections_(cfg.getParameter<bool>("applyEnergyCorrections"))
 , path_(cfg.getParameter<std::string>("path"))
 , debug_(cfg.getParameter<bool>("debug"))
{
  muonRochCorr_.init(edm::FileInPath(path_).fullPath());

  rand_ = TRandom3();

  consumeTemplate<edm::View<pat::Muon> >(muons_);

  produces<pat::MuonCollection>();
}

void MuonPATUserDataRochesterCorrectionAdder::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<pat::Muon> > muons;
  evt.getByLabel(muons_, muons);

  std::unique_ptr<pat::MuonCollection> newMuons(new pat::MuonCollection);
  newMuons->reserve(muons->size());

  for(unsigned int i_muo=0; i_muo<muons->size(); ++i_muo)
  {
    newMuons->emplace_back(muons->at(i_muo));
    pat::Muon& muon = newMuons->back();

    double scaleFactor(1.);
    double scaleFactor_stat_RMS(1.);
    double scaleFactor_Zpt(1.);
    double scaleFactor_Ewk(1.);
    double scaleFactor_deltaM(1.);
    double scaleFactor_Ewk2(1.);

    // Data
    if(evt.isRealData())
    {
      scaleFactor = muonRochCorr_.kScaleDT(muon.charge(), muon.pt(), muon.eta(), muon.phi());
    }
    // MC
    else if(muon.genLepton())
    {
      scaleFactor = muonRochCorr_.kSpreadMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon.genLepton()->pt());

      double rms = 0.;
      for(int iVar=0; iVar<100; ++iVar)
      {
        const double temp_sf = muonRochCorr_.kSpreadMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon.genLepton()->pt(), 1, iVar);

        rms += (temp_sf-scaleFactor)*(temp_sf-scaleFactor);
      }

      scaleFactor_stat_RMS = sqrt(rms/99.);

      scaleFactor_Zpt    = muonRochCorr_.kSpreadMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon.genLepton()->pt(), 2);
      scaleFactor_Ewk    = muonRochCorr_.kSpreadMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon.genLepton()->pt(), 3);
      scaleFactor_deltaM = muonRochCorr_.kSpreadMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon.genLepton()->pt(), 4);
      scaleFactor_Ewk2   = muonRochCorr_.kSpreadMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon.genLepton()->pt(), 5);
    }
    else if(muon.innerTrack().isNonnull())
    {
      rand_.SetSeed(muon.hasUserInt("deterministicSeed") ? ((uint32_t) muon.userInt("deterministicSeed")) : 1);

      const auto rdnm_val = rand_.Rndm();

      const auto muon_TrkLayers = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement();

      scaleFactor = muonRochCorr_.kSmearMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon_TrkLayers, rdnm_val);

      double rms = 0.;
      for(int iVar = 0; iVar < 100; ++iVar)
      {
        const double temp_sf = muonRochCorr_.kSmearMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon_TrkLayers, rdnm_val, 1, iVar);

        rms += (temp_sf-scaleFactor)*(temp_sf-scaleFactor);
      }

      scaleFactor_stat_RMS = sqrt(rms/99.);

      scaleFactor_Zpt    = muonRochCorr_.kSmearMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon_TrkLayers, rdnm_val, 2);
      scaleFactor_Ewk    = muonRochCorr_.kSmearMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon_TrkLayers, rdnm_val, 3);
      scaleFactor_deltaM = muonRochCorr_.kSmearMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon_TrkLayers, rdnm_val, 4);
      scaleFactor_Ewk2   = muonRochCorr_.kSmearMC(muon.charge(), muon.pt(), muon.eta(), muon.phi(), muon_TrkLayers, rdnm_val, 5);
    }

    if(applyEnergyCorrections_)
    {
      if(debug_){ std::cout << muon.pt() << " "; }

      muon.setP4(math::XYZTLorentzVector(muon.px()*scaleFactor, muon.py()*scaleFactor, muon.pz()*scaleFactor, muon.energy()*scaleFactor));

      if(debug_){ std::cout << muon.pt() << std::endl; }
    }

    muon.addUserFloat("MuonEnergyCorr"          , scaleFactor          );
    muon.addUserFloat("MuonEnergyCorr_stat_RMS" , scaleFactor_stat_RMS );
    muon.addUserFloat("MuonEnergyCorr_Zpt"      , scaleFactor_Zpt      );
    muon.addUserFloat("MuonEnergyCorr_Ewk"      , scaleFactor_Ewk      );
    muon.addUserFloat("MuonEnergyCorr_deltaM"   , scaleFactor_deltaM   );
    muon.addUserFloat("MuonEnergyCorr_Ewk2"     , scaleFactor_Ewk2     );
  }

  if(applyEnergyCorrections_)
  {
    std::sort(newMuons->begin(), newMuons->end(), GreaterByPt<pat::Muon>());
  }

  evt.put(std::move(newMuons));

  return;
}

void MuonPATUserDataRochesterCorrectionAdder::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("slimmedMuons"))->setComment("edm::View<pat::Muon>");
  desc.add<bool>("applyEnergyCorrections", false)->setComment("muon energy corrections applied to output collection");
  desc.add<std::string>("path", std::string("TopAnalysis/ZTopUtils/data/MuonRochCorr.txt"))->setComment("path to input file with Rochester corrections (relative to ${CMSSW_BASE}/src)");
  desc.add<bool>("debug", false)->setComment("printout debug information");
  descriptions.add("muonPATUserDataRochesterCorrectionAdder", desc);
}

DEFINE_FWK_MODULE(MuonPATUserDataRochesterCorrectionAdder);

#endif
