#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/Run.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Utilities/interface/InputTag.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"


using namespace cms::Ort;

class MLSystematicWeightsProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {

  public:
    explicit MLSystematicWeightsProducer(const edm::ParameterSet&, const ONNXRuntime*);
    ~MLSystematicWeightsProducer() override;

    static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(const ONNXRuntime*);

  protected:
    void produce(edm::Event&, const edm::EventSetup&) override;
    const edm::InputTag tag_LHEEventProduct_;
    const bool verbose_;
    const bool doBFragWeights_;
    const bool doHdampWeights_;
    const float defaultVal_;
    const edm::EDPutTokenT<float> weightToken_;
    const edm::EDPutTokenT<float> xbtopToken_;
    const edm::EDPutTokenT<float> xbantitopToken_;
    const edm::EDPutTokenT<float> topBhadPtToken_;
    const edm::EDPutTokenT<float> antitopBhadPtToken_;
    edm::EDGetToken lheEvtToken_;
    std::vector<std::string> input_names_;
    std::vector<std::string> output_names_;
    std::vector<std::vector<int64_t>> input_shapes_;
    FloatArrays data_; // each stream hosts its own data

    const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    const edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken_;
    const edm::EDGetTokenT<std::vector<int> > genBHadFromTopWeakDecayToken_;
    const edm::EDGetTokenT<std::vector<reco::GenParticle> > genBHadPlusMothersToken_;
    const edm::EDGetTokenT<std::vector<std::vector<int> > > genBHadPlusMothersIndicesToken_;
    const edm::EDGetTokenT<std::vector<int> > genBHadIndexToken_;
    std::multimap<int, const reco::GenParticle*> firstBHadFromTop(const std::vector<int>&,
                                                           const std::vector<int>&,
                                                           const std::vector<int>&,
                                                           const std::vector<reco::GenParticle>&);
};

/////////////////////////////////////////////////////////////////////////////////

MLSystematicWeightsProducer::MLSystematicWeightsProducer(const edm::ParameterSet& cfg, const ONNXRuntime* cache):
  tag_LHEEventProduct_(cfg.getParameter<edm::InputTag>("LHEEventProduct")),
  verbose_(cfg.getParameter<bool>("verbose")),
  doBFragWeights_(cfg.getParameter<bool>("doBFragWeights")),
  doHdampWeights_(cfg.getParameter<bool>("doHdampWeights")),
  defaultVal_((float)cfg.getParameter<double>("defaultVal")),
  weightToken_(produces<float>("weight")),
  xbtopToken_(produces<float>("xbtop")),
  xbantitopToken_(produces<float>("xbantitop")),
  topBhadPtToken_(produces<float>("topbhadpt")),
  antitopBhadPtToken_(produces<float>("antitopbhadpt")),
  input_names_(cfg.getParameter<std::vector<std::string>>("input_names")),
  output_names_(cfg.getParameter<std::vector<std::string>>("output_names")),
  input_shapes_(),
  genParticlesToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))),
  genBHadFlavourToken_(consumes<std::vector<int> >(cfg.getParameter<edm::InputTag>("genBHadFlavour"))),
  genBHadFromTopWeakDecayToken_(consumes<std::vector<int> >(cfg.getParameter<edm::InputTag>("genBHadFromTopWeakDecay"))),
  genBHadPlusMothersToken_(consumes<std::vector<reco::GenParticle> >(cfg.getParameter<edm::InputTag>("genBHadPlusMothers"))),
  genBHadPlusMothersIndicesToken_(consumes<std::vector<std::vector<int> > >(cfg.getParameter<edm::InputTag>("genBHadPlusMothersIndices"))),
  genBHadIndexToken_(consumes<std::vector<int> >(cfg.getParameter<edm::InputTag>("genBHadIndex")))
{
  lheEvtToken_ = consumes<LHEEventProduct>(tag_LHEEventProduct_);
  if (doHdampWeights_ && doBFragWeights_) {
      throw cms::Exception("MLSystematicWeightsProducer::produce") << "Both doHdampWeights and doBFragWeights are set to true. Only one option is valid. Setting default weights.";
  }
  if (doHdampWeights_) {
      input_shapes_.push_back({1,2,6});
      data_.emplace_back(12, 0);
  }
  if (doBFragWeights_) {
      input_shapes_.push_back({1,2,2});
      data_.emplace_back(4, 0);
  }
  // produces<float>();
}

/////////////////////////////////////////////////////////////////////////////////

MLSystematicWeightsProducer::~MLSystematicWeightsProducer()
{
}

std::unique_ptr<ONNXRuntime> MLSystematicWeightsProducer::initializeGlobalCache(const edm::ParameterSet& cfg) {
  return std::make_unique<ONNXRuntime>(cfg.getParameter<edm::FileInPath>("model_path").fullPath());
}

void MLSystematicWeightsProducer::globalEndJob(const ONNXRuntime* cache) {}

/////////////////////////////////////////////////////////////////////////////////

void MLSystematicWeightsProducer::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  float weight(defaultVal_);
  if (evt.isRealData()) { return; }

  if (doHdampWeights_) {

    edm::Handle<LHEEventProduct> lhe_info;
    evt.getByToken(lheEvtToken_, lhe_info);

    if (!lhe_info.isValid()) {
      return;
    }

    float topPt(0.);
    float topRapidity(0.);
    float topPhi(0.);
    float topMass(0.);
    float topPDGID(0.);
    float antitopPt(0.);
    float antitopRapidity(0.);
    float antitopPhi(0.);
    float antitopMass(0.);
    float antitopPDGID(0.);
    float hdamp(1.379); //This value is the default value of hdamp divided by 172.5
    float maxM(243.9517); //This value is needed to normalise the mass of the particles in each event and comes from the maximum mass value we had in the training+validation sample

    // Get the necessary information
    // https://twiki.cern.ch/twiki/bin/view/CMS/MLReweighting
    // https://twiki.cern.ch/twiki/pub/CMS/MLReweighting/readNanoAODShower.py.txt
    // LHE level top/antitop: lheParticles idx 0-4 px/py/pz/e/m + pdgID

    ROOT::Math::PxPyPzEVector top4;
    ROOT::Math::PxPyPzEVector antitop4;
    const lhef::HEPEUP& lheEvent = lhe_info->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
    for (size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle) {
      int pdgId = lheEvent.IDUP.at(idxParticle);
      if (pdgId == 6) { // top quark
          top4 = ROOT::Math::PxPyPzEVector(lheParticles.at(idxParticle)[0], lheParticles.at(idxParticle)[1], lheParticles.at(idxParticle)[2], lheParticles.at(idxParticle)[3]);
          topPt = top4.Pt();
          topRapidity = top4.Rapidity();
          topPhi = top4.Phi();
          topMass = top4.M();
          topPDGID = pdgId;
      }
      if (pdgId == -6) { // top antiquark
          antitop4 = ROOT::Math::PxPyPzEVector(lheParticles.at(idxParticle)[0], lheParticles.at(idxParticle)[1], lheParticles.at(idxParticle)[2], lheParticles.at(idxParticle)[3]);
          antitopPt = antitop4.Pt();
          antitopRapidity = antitop4.Rapidity();
          antitopPhi = antitop4.Phi();
          antitopMass = antitop4.M();
          antitopPDGID = pdgId;
      }
    }

    // get ttbar 4 vector
    auto ttbar = top4 + antitop4;
    float ttbarPt(ttbar.Pt());

    // perform necessary transformations
    topPt = std::log10(topPt);
    antitopPt = std::log10(antitopPt);
    topMass = topMass / maxM;
    antitopMass = antitopMass / maxM;
    topPDGID = topPDGID>0 ? 0.1 : 0.2;
    antitopPDGID = antitopPDGID>0 ? 0.1 : 0.2;

    // set input data tensor
    data_[0][0] = topPt;
    data_[0][1] = topRapidity;
    data_[0][2] = topPhi;
    data_[0][3] = topMass;
    data_[0][4] = topPDGID;
    data_[0][5] = hdamp;

    data_[0][6]  = antitopPt;
    data_[0][7]  = antitopRapidity;
    data_[0][8]  = antitopPhi;
    data_[0][9]  = antitopMass;
    data_[0][10] = antitopPDGID;
    data_[0][11] = hdamp;

    // run inference
    std::vector<float> outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
    if (ttbarPt < 1000.) {
      weight = outputs[0] / outputs[1];
    } else {
      weight = 1.0;
    }

    if (verbose_) {
      std::cout << "\n" << "MLSystematicWeightsProducer | Run=" << evt.id().run() << ", LuminosityBlock=" << evt.id().luminosityBlock() << ", Event=" << evt.id().event() << "\n\n";
      std::cout << "input data -> ";
      for (auto &i: data_[0]) {
        std::cout << i << " ";
      }
      std::cout << std::endl << "output data -> ";
      for (auto &i: outputs) {
        std::cout << i << " ";
      }
      std::cout << std::endl;
      std::cout << "final weight = " << weight << std::endl;
    }

  }

  // Bfragmentation uses hadron information
  else if (doBFragWeights_){

    edm::Handle<LHEEventProduct> lhe_info;
    evt.getByToken(lheEvtToken_, lhe_info);

    if (!lhe_info.isValid()) {
      return;
    }

    edm::Handle<std::vector<int> > genBHadIndex;
    evt.getByToken(genBHadIndexToken_, genBHadIndex);

    edm::Handle<std::vector<int> > genBHadFlavour;
    evt.getByToken(genBHadFlavourToken_, genBHadFlavour);

    edm::Handle<std::vector<int> > genBHadFromTopWeakDecay;
    evt.getByToken(genBHadFromTopWeakDecayToken_, genBHadFromTopWeakDecay);

    edm::Handle<std::vector<reco::GenParticle> > genBHadPlusMothers;
    evt.getByToken(genBHadPlusMothersToken_, genBHadPlusMothers);

    edm::Handle<std::vector<std::vector<int> > > genBHadPlusMothersIndices;
    evt.getByToken(genBHadPlusMothersIndicesToken_, genBHadPlusMothersIndices);

    ROOT::Math::PxPyPzEVector bhadtop4;
    ROOT::Math::PxPyPzEVector bhadantitop4;

    //std::vector<const reco::GenParticle*> bHads;
    std::multimap<int, const reco::GenParticle*> bHadsAndPdgIds;

    bool hasTwoBHads = false;

    if (!genBHadIndex.failedToGet()) {
      bHadsAndPdgIds = firstBHadFromTop(*genBHadIndex, *genBHadFlavour, *genBHadFromTopWeakDecay, *genBHadPlusMothers);

      for (auto had = bHadsAndPdgIds.begin(); had != bHadsAndPdgIds.end(); had++) {
        //auto tmpBhad = (&**had)->p4();
        auto tmpBhad = (*had).second->p4();
        if (verbose_) std::cout << "pT of b hadron: " << tmpBhad.pt() << " with pdgId " << (*had).first << std::endl;
        if (had->first == 6) bhadtop4 = ROOT::Math::PxPyPzEVector(tmpBhad.px(), tmpBhad.py(), tmpBhad.pz(), tmpBhad.e());
        else if (had->first == - 6) bhadantitop4 = ROOT::Math::PxPyPzEVector(tmpBhad.px(), tmpBhad.py(), tmpBhad.pz(), tmpBhad.e());
      }
      if (verbose_) std::cout <<"end of event" <<std::endl;
    }

    if (bHadsAndPdgIds.size() == 2 and bhadtop4.Pt() > 0 and bhadantitop4.Pt() > 0) hasTwoBHads = true;

    float top_bhad_pt(-1.0);
    float antitop_bhad_pt(-1.0);

    top_bhad_pt = bhadtop4.Pt();
    antitop_bhad_pt = bhadantitop4.Pt();

    float wplusMass(0.);
    float wminusMass(0.);
    float topMass(0.);
    float antitopMass(0.);
    float xbtop(-1.0);
    float xbantitop(-1.0);

    // Get the necessary information
    // https://twiki.cern.ch/twiki/bin/view/CMS/MLReweighting
    // https://twiki.cern.ch/twiki/pub/CMS/MLReweighting/readNanoAODShower.py.txt
    // LHE level top/antitop: lheParticles idx 0-4 px/py/pz/e/m + pdgID

    ROOT::Math::PxPyPzEVector wplus4;
    ROOT::Math::PxPyPzEVector wminus4;
    ROOT::Math::PxPyPzEVector top4;
    ROOT::Math::PxPyPzEVector antitop4;
    const lhef::HEPEUP& lheEvent = lhe_info->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
    for(size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle){
      int pdgId = lheEvent.IDUP.at(idxParticle);
      if (pdgId == 24) { // W plus
          wplus4 = ROOT::Math::PxPyPzEVector(lheParticles.at(idxParticle)[0], lheParticles.at(idxParticle)[1], lheParticles.at(idxParticle)[2], lheParticles.at(idxParticle)[3]);
          wplusMass = wplus4.M();
      }
      if (pdgId == -24) { // W minus
          wminus4 = ROOT::Math::PxPyPzEVector(lheParticles.at(idxParticle)[0], lheParticles.at(idxParticle)[1], lheParticles.at(idxParticle)[2], lheParticles.at(idxParticle)[3]);
          wminusMass = wminus4.M();
      }
      if (pdgId == 6) { // top quark
          top4 = ROOT::Math::PxPyPzEVector(lheParticles.at(idxParticle)[0], lheParticles.at(idxParticle)[1], lheParticles.at(idxParticle)[2], lheParticles.at(idxParticle)[3]);
          topMass = top4.M();
      }
      if (pdgId == -6) { // top antiquark
          antitop4 = ROOT::Math::PxPyPzEVector(lheParticles.at(idxParticle)[0], lheParticles.at(idxParticle)[1], lheParticles.at(idxParticle)[2], lheParticles.at(idxParticle)[3]);
          antitopMass = antitop4.M();
      }
    }

    float wtop = (wplusMass * wplusMass) / (topMass * topMass);
    float wantitop = (wminusMass * wminusMass) / (antitopMass * antitopMass);
    float xEtop = (2 * bhadtop4).Dot(top4) / (topMass * topMass);
    float xEantitop = (2 * bhadantitop4).Dot(antitop4) / (antitopMass * antitopMass);

    xbtop = xEtop / (1 - wtop);
    xbantitop = xEantitop / (1 - wantitop);

    if (verbose_) {
      std::cout << "W+ " << wplusMass << " t: " << topMass << " W- " << wminusMass<< " tbar: " << antitopMass << std::endl;
      std::cout << "wtop: " << wtop << " wantitop: " << wantitop << std::endl;
      std::cout << "xEtop: " << xEtop << " xEantitop: " << xEantitop << std::endl;
      std::cout << "xb top: " << xbtop << " and xb antitop: " << xbantitop << std::endl;
    }

    // rb: parameter of the Lund-Bowler function of Pythia. 1.056 for CMS nominal CP5, 1.254 for up var.
    // rb should be fixed to 0.855 for input array.

    // set input data tensor
    data_[0][0] = std::min(xbtop, static_cast<float>(1.2));
    //data_[0][0] = xbtop;
    data_[0][1] = 0.855;

    data_[0][2] = std::min(xbantitop, static_cast<float>(1.2));
    //data_[0][2] = xbantitop;
    data_[0][3] = 0.855;

    // run inference
    std::vector<float> outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
    weight = outputs[0] / outputs[1];

    if (verbose_) {
      std::cout << "\n" << "MLSystematicWeightsProducer | Run=" << evt.id().run() << ", LuminosityBlock=" << evt.id().luminosityBlock() << ", Event=" << evt.id().event() << "\n\n";
      std::cout << "input data -> ";
      for (auto &i: data_[0]) {
        std::cout << i << " ";
      }
      std::cout << std::endl << "output data -> ";
      for (auto &i: outputs) {
        std::cout << i << " ";
      }
      std::cout << std::endl;
      std::cout << "final weight = " << weight << std::endl;
    }

    if (!hasTwoBHads) {
      weight = 1.0;
      if (bhadtop4.Pt() <= 0) {
        xbtop = -1.0;
        top_bhad_pt = -1.0;
      }
      if (bhadantitop4.Pt() <= 0) {
        xbantitop = -1.0;
        antitop_bhad_pt = -1.0;
      }
    }
    evt.emplace(xbtopToken_, xbtop);
    evt.emplace(xbantitopToken_, xbantitop);
    evt.emplace(topBhadPtToken_, top_bhad_pt);
    evt.emplace(antitopBhadPtToken_, antitop_bhad_pt);
  }
  // set output of producer
  evt.emplace(weightToken_, weight);
}

/////////////////////////////////

// With help of PhysicsTools/JetMCAlgos/plugins/ttHFGenFilter.cc
std::multimap<int, const reco::GenParticle*>
    MLSystematicWeightsProducer::firstBHadFromTop(const std::vector<int>& genBHadIndex,
                                                  const std::vector<int>& genBHadFlavour,
                                                  const std::vector<int>& genBHadFromTopWeakDecay,
                                                  const std::vector<reco::GenParticle>& genBHadPlusMothers) {

    //std::vector<const reco::GenParticle*> out;
    std::multimap<int, const reco::GenParticle*> out;

    for (unsigned int i = 0; i < genBHadIndex.size(); i++) {
      const reco::GenParticle* bhadron = genBHadIndex[i] >= 0 && genBHadIndex[i] < int(genBHadPlusMothers.size())
                                             ? &(genBHadPlusMothers[genBHadIndex[i]])
                                             : nullptr;
      int motherflav = genBHadFlavour[i];
      int ifAfterTopWeakDecay = genBHadFromTopWeakDecay[i];
      if (abs(motherflav) == 6 and ifAfterTopWeakDecay == 1)
        //out.emplace_back(bhadron);
        out.insert(std::pair<int, const reco::GenParticle*> (motherflav, bhadron));
    }

    if (out.size() > 2)
      std::cout << "There is more than two b hadrons found!" << std::endl;
    else if (out.size() < 2)
      std::cout << "There is less than two b hadrons found!" << std::endl;

    return out;
}

/////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(MLSystematicWeightsProducer);
