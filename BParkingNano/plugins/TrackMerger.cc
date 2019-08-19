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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "helper.h"

class TrackMerger : public edm::global::EDProducer<> {


public:

  //would it be useful to give this a bit more standard structure?
  explicit TrackMerger(const edm::ParameterSet &cfg):
    beamSpotSrc_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgMuonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    trkPtCut_(cfg.getParameter<double>("trkPtCut")),
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
    dzTrg_cleaning_(cfg.getParameter<double>("dzTrg_cleaning")),
    drTrg_Cleaning_(cfg.getParameter<double>("drTrg_Cleaning")),
    dcaSig_(cfg.getParameter<double>("dcaSig")),
    trkNormChiMin_(cfg.getParameter<int>("trkNormChiMin")),
    trkNormChiMax_(cfg.getParameter<int>("trkNormChiMax")) 
{
    produces<pat::CompositeCandidateCollection>("SelectedTracks");  
    produces<TransientTrackCollection>("SelectedTransientTracks");  
}

  ~TrackMerger() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;

  //selections                                                                 
  const double trkPtCut_;
  const double trkEtaCut_;
  const double dzTrg_cleaning_;
  const double drTrg_Cleaning_;
  const double dcaSig_;
  const int trkNormChiMin_;
  const int trkNormChiMax_;
};






void TrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {

  //input
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from EventSetup" ;
  }  
  reco::BeamSpot beamSpot = *beamSpotHandle;

  edm::ESHandle<MagneticField> bFieldHandle;
  stp.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  edm::Handle<pat::PackedCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  evt.getByToken(lostTracksToken_, lostTracks);
  edm::Handle<pat::MuonCollection> trgMuons;
  evt.getByToken(trgMuonToken_, trgMuons);

  //for lost tracks / pf discrimination
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();

  //ok this was CompositeCandidateCollection 
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out      (new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection>          trans_tracks_out(new TransientTrackCollection);


  //try topreserve same logic avoiding the copy of the full collection
  /*
  //correct logic but a bit convoluted -> changing to smthn simpler
   std::vector<pat::PackedCandidate> totalTracks(*tracks);
   totalTracks.insert(totalTracks.end(),lostTracks->begin(),lostTracks->end());
  */
 
  // for loop is better to be range based - especially for large ensembles  
  for(unsigned int iTrk=0; iTrk<totalTracks; ++iTrk){
    const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];

    //arranging cuts for speed
    if(!trk.hasTrackDetails()) continue;
    if(abs(trk.pdgId()) != 211) continue; //do we want also to keep muons?
    if(trk.pt() < trkPtCut_ ) continue;
    if(fabs(trk.eta()) > trkEtaCut_) continue;

    if( (trk.pseudoTrack().normalizedChi2() < trkNormChiMin_ &&
         trkNormChiMin_>=0 ) ||
        (trk.pseudoTrack().normalizedChi2() > trkNormChiMax_ &&
         trkNormChiMax_>0)  )    continue; 

    bool skipTrack=true;
    for (const pat::Muon & mu: *trgMuons){
      //remove tracks inside trg muons jet
      if(reco::deltaR(trk, mu) < drTrg_Cleaning_ && drTrg_Cleaning_ >0) 
        continue;
      //if dz is negative it is deactivated
      if((fabs(trk.vz() - mu.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0))
        continue;
      skipTrack=false;
      break; // at least for one trg muon to pass this cuts
    }
    // if track is closer to at least a triggering muon keep it
    if (skipTrack) continue;

    // high purity requirment applied only in packedCands
    if( iTrk < nTracks && !trk.trackHighPurity()) continue;
   
    // build transient track
    const reco::TransientTrack trackTT((*(trk.bestTrack())), &(*bFieldHandle));
    if (!trackTT.isValid()) continue;
   
    //distance closest approach in x,y wrt beam spot
    std::pair<double,double> DCA = computeDCA(trackTT, beamSpot);
    float DCABS = DCA.first;
    float DCABSErr = DCA.second;
    float DCASig = (DCABSErr != 0 && float(DCABSErr) == DCABSErr) ? fabs(DCABS/DCABSErr) : -1;
    if (DCASig >  dcaSig_  && dcaSig_ >0) continue;

    pat::CompositeCandidate pcand;
    pcand.setP4(trk.p4());
    pcand.setCharge(trk.charge());
    pcand.setVertex(trk.vertex());
    pcand.setPdgId(trk.pdgId());
    pcand.addUserInt("isPacked", (iTrk < nTracks) ? 1 : 0);
    pcand.addUserInt("isLostTrk", (iTrk < nTracks) ? 0 : 1);      
    pcand.addUserFloat("dxy", trk.dxy());
    pcand.addUserFloat("dxyS", trk.dxy()/trk.dxyError());
    pcand.addUserFloat("dz", trk.dz()); 
    pcand.addUserFloat("dzS", trk.dz()/trk.dzError());
    pcand.addUserFloat("DCASig", DCASig);
    //adding the candidate in the composite stuff for fit (need to test)
    if ( iTrk < nTracks )
      pcand.addUserCand( "cand", edm::Ptr<pat::PackedCandidate> ( tracks, iTrk ));
    else 
      pcand.addUserCand( "cand", edm::Ptr<pat::PackedCandidate> ( lostTracks, iTrk-nTracks ));
    
    tracks_out       -> emplace_back(pcand);
    trans_tracks_out -> emplace_back(trackTT);
   
  }
 
  evt.put(std::move(tracks_out),       "SelectedTracks");
  evt.put(std::move(trans_tracks_out), "SelectedTransientTracks");
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackMerger);
