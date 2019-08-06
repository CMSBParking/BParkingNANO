// Merges the PF and LowPT collections, sets the isPF and isLowPt 
// UserInt's accordingly

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <limits>
#include <algorithm>

namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }

class PATLowPtElectronSeedingEmbedder : public edm::global::EDProducer<> {
public:
  explicit PATLowPtElectronSeedingEmbedder(const edm::ParameterSet &cfg):
    lowpt_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ptBiased_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("ptbiasedSeeding") )},
    unBiased_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("unbiasedSeeding") )},
    mvaId_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("mvaId") )},
    minBdtUnbiased_{cfg.getParameter<double>("minBdtUnbiased")} {
      produces<pat::ElectronCollection>();
    }

  ~PATLowPtElectronSeedingEmbedder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::ElectronCollection> lowpt_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ptBiased_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> unBiased_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mvaId_src_;
  const double minBdtUnbiased_;
};

void PATLowPtElectronSeedingEmbedder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  edm::Handle<pat::ElectronCollection> lowpt;
  evt.getByToken(lowpt_src_, lowpt);

  edm::Handle<edm::ValueMap<float> > ptBiased;  
  evt.getByToken(ptBiased_src_, ptBiased);
  edm::Handle<edm::ValueMap<float> > unBiased;  
  evt.getByToken(unBiased_src_, unBiased);
  edm::Handle<edm::ValueMap<float> > mvaId;  
  evt.getByToken(mvaId_src_, mvaId);


  std::unique_ptr<pat::ElectronCollection> out(new pat::ElectronCollection);

  for ( size_t iele = 0; iele < lowpt->size(); ++iele ) {
    edm::Ref<pat::ElectronCollection> ref(lowpt,iele);
    const reco::GsfTrackRef gsfTrk = ref->gsfTrack();
    float unbiased_seedBDT = float((*unBiased)[gsfTrk]);
    float ptbiased_seedBDT = float((*ptBiased)[gsfTrk]);
    float mva_id = float((*mvaId)[ref]);

    if(unbiased_seedBDT < minBdtUnbiased_) continue;

    pat::Electron ele = *ref;
    ele.addUserFloat("ptBiased", ptbiased_seedBDT);
    ele.addUserFloat("unBiased", unbiased_seedBDT);
    ele.addUserFloat("mvaId", mva_id);
    out->push_back(ele);
  }

  evt.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATLowPtElectronSeedingEmbedder);
