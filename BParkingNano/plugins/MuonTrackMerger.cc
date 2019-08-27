// Merges the PFPackedCandidates and Lost tracks
// beam spot readout in case dcasig to be calculated wrt beam spot
// currently computed wrt triggeringMuon vertex
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// #include "DataFormats/BeamSpot/interface/BeamSpot.h"
// #include "MagneticField/Engine/interface/MagneticField.h"
// #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
// #include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "helper.h"

class MuonTrackMerger : public edm::global::EDProducer<> {


public:

  //would it be useful to give this a bit more standard structure?
  explicit MuonTrackMerger(const edm::ParameterSet &cfg):
    muonsToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
    muons_ttracks_(consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("muonTransientTracks"))),
    tracksToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    tracks_ttracks_(consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("trackTransientTracks"))),
    muonSelection_(cfg.getParameter<std::string>("muonSelection")),
    trackSelection_(cfg.getParameter<std::string>("trackSelection")),
    sortOutputCollections_(cfg.getParameter<bool>("sortOutputCollections"))
{
    produces<pat::CompositeCandidateCollection>("SelectedMuonsTracks");  
    produces<TransientTrackCollection>("SelectedTransientMuonsTracks");  
}

  ~MuonTrackMerger() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:
  const edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  const edm::EDGetTokenT<TransientTrackCollection> muons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<TransientTrackCollection> tracks_ttracks_;

  const StringCutObjectSelector<pat::CompositeCandidate> muonSelection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> trackSelection_;

  const bool sortOutputCollections_;
};






void MuonTrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {

  //input
  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muonsToken_, muons);

  edm::Handle<TransientTrackCollection> muonsTT;
  evt.getByToken(muons_ttracks_, muonsTT);

  edm::Handle<pat::CompositeCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);

  edm::Handle<TransientTrackCollection> tracksTT;
  evt.getByToken(tracks_ttracks_, tracksTT);



  // output                                                                                                                                                                   
  std::unique_ptr<pat::CompositeCandidateCollection> muonTrack_out  (new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection> trans_muonTrack_out     (new TransientTrackCollection);
  std::vector<float> pts;

  auto muonSize = muons->size();
  for(size_t mu_idx = 0; mu_idx < muonSize; ++mu_idx) {

    const pat::Muon& mu = (*muons)[mu_idx];

    pat::CompositeCandidate pcand;
    pcand.setP4(mu.p4());
    pcand.setCharge(mu.charge());
    pcand.setVertex(mu.vertex());
    pcand.setPdgId(mu.pdgId());
    pcand.addUserInt("isPF", 1);
    pcand.addUserInt("isTrack", 0);
    pcand.addUserFloat("dxy", mu.dB(pat::Muon::PV2D));
    pcand.addUserFloat("dxyS", mu.edB(pat::Muon::PV2D));
    pcand.addUserFloat("dz", mu.dB(pat::Muon::PVDZ));
    pcand.addUserFloat("dzS", mu.edB(pat::Muon::PVDZ));

    if(!muonSelection_(pcand)) continue;

    muonTrack_out->emplace_back(pcand);
    trans_muonTrack_out->emplace_back(muonsTT->at(mu_idx));
    pts.push_back(pcand.pt());
  }

  auto trkSize = tracks->size();
  for(size_t trk_idx = 0; trk_idx < trkSize; ++trk_idx) {

    const pat::CompositeCandidate& trk = (*tracks)[trk_idx];
    if(!trackSelection_(trk)) continue;

    pat::CompositeCandidate pcand;
    pcand.setP4(trk.p4());
    pcand.setCharge(trk.charge());
    pcand.setVertex(trk.vertex());
    pcand.setPdgId(trk.pdgId());
    pcand.addUserInt("isPF", 0);
    pcand.addUserInt("isTrack", 1);
    pcand.addUserFloat("dxy", trk.userFloat("dxy"));
    pcand.addUserFloat("dxyS", trk.userFloat("dxyS"));
    pcand.addUserFloat("dz", trk.userFloat("dz"));
    pcand.addUserFloat("dzS", trk.userFloat("dzS"));

    muonTrack_out->emplace_back(pcand);
    trans_muonTrack_out->emplace_back(muonsTT->at(trk_idx));
    pts.push_back(pcand.pt());
  }

  for(auto ij : pts) std::cout << " pre sort pt = " << ij << std::endl;

  if(sortOutputCollections_){
    std::vector<size_t> sortedIndex = decrease_sorted_indices(pts);

    std::unique_ptr<pat::CompositeCandidateCollection>  out_tmp      (new pat::CompositeCandidateCollection );
    std::unique_ptr<TransientTrackCollection> trans_out_tmp(new TransientTrackCollection);
    out_tmp->reserve(muonTrack_out->size());
    trans_out_tmp->reserve(trans_muonTrack_out->size());

    for (std::size_t i=0; i<sortedIndex.size(); ++i) {
      out_tmp->emplace_back((*muonTrack_out)[sortedIndex[i]]);
      trans_out_tmp->emplace_back((*trans_muonTrack_out)[sortedIndex[i]]);
    }
    muonTrack_out.swap(out_tmp);
    trans_muonTrack_out.swap(trans_out_tmp);
  }

  for(auto ij : *muonTrack_out) std::cout << " post sort pt = " << ij.pt() << std::endl;

  evt.put(std::move(muonTrack_out),       "SelectedMuonsTracks");
  evt.put(std::move(trans_muonTrack_out), "SelectedTransientMuonsTracks");
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonTrackMerger);
