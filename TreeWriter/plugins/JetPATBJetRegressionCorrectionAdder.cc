#ifndef TopAnalysis_ZTopUtils_JetPATBJetRegressionCorrectionAdder
#define TopAnalysis_ZTopUtils_JetPATBJetRegressionCorrectionAdder

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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "consumeTemplate.h"


class JetPATBJetRegressionCorrectionAdder : public ProducerTemplate {

 public:

  explicit JetPATBJetRegressionCorrectionAdder(const edm::ParameterSet&);
  ~JetPATBJetRegressionCorrectionAdder() {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 protected:

  virtual void produce(edm::Event&, const edm::EventSetup&);

  const edm::InputTag jets_;

  const bool applyLooseSelection_;
  
  const bool debug_;

};

JetPATBJetRegressionCorrectionAdder::JetPATBJetRegressionCorrectionAdder(const edm::ParameterSet& cfg)
 : jets_(cfg.getParameter<edm::InputTag>("src"))
 , applyLooseSelection_(cfg.getParameter<bool>("applyLooseSelection"))
 , debug_(cfg.getParameter<bool>("debug"))
{

  consumeTemplate<edm::View<pat::Jet> >(jets_);

  produces<pat::JetCollection>();
}

void JetPATBJetRegressionCorrectionAdder::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<edm::View<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  std::unique_ptr<pat::JetCollection> newJets(new pat::JetCollection);
  newJets->reserve(jets->size());

  for(unsigned int i_jet=0; i_jet<jets->size(); ++i_jet)
  {
    newJets->emplace_back(jets->at(i_jet));
    pat::Jet& jet = newJets->back();
    
    
    
    float scaleFactor = jet.hasUserFloat("BJetEnergyCorrFactor") ? jet.userFloat("BJetEnergyCorrFactor") : 1.;
    
    if(debug_){ std::cout << jet.pt() << " "; }
    
    if(applyLooseSelection_ && (jet.pt()<15 || jet.eta()>2.5)) continue;

    jet.setP4(math::XYZTLorentzVector(jet.px()*scaleFactor, jet.py()*scaleFactor, jet.pz()*scaleFactor, jet.energy()*scaleFactor));

    if(debug_){ std::cout << jet.pt() << std::endl; }

  }

  std::sort(newJets->begin(), newJets->end(), GreaterByPt<pat::Jet>());

  evt.put(std::move(newJets));

  return;
}

void JetPATBJetRegressionCorrectionAdder::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("slimmedJets"))->setComment("edm::View<pat::Jet>");
  desc.add<bool>("applyLooseSelection", false)->setComment("only correct Jets with ");
  desc.add<bool>("debug", false)->setComment("printout debug information");
  descriptions.add("JetPATBJetRegressionCorrectionAdder", desc);
}

DEFINE_FWK_MODULE(JetPATBJetRegressionCorrectionAdder);

#endif
