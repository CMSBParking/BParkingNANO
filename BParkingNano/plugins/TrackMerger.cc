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


class TrackMerger : public edm::global::EDProducer<> {


public:

  //would it be useful to give this a bit more standard structure?
  explicit TrackMerger(const edm::ParameterSet &cfg):
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgMuonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    trkPtCut_(cfg.getParameter<double>("trkPtCut")),
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
    dzTrg_cleaning_(cfg.getParameter<double>("dzTrg_cleaning")),
  //  drTrg_ProbeCleaning_(cfg.getParameter<double>("drTrg_ProbeCleaning")),
  //  drTrg_TagCleaning_(cfg.getParameter<double>("drTrg_TagCleaning")),
  //  dcaSig_probe_(cfg.getParameter<double>("dcaSig_probe")),
  //  dcaSig_tag_(cfg.getParameter<double>("dcaSig_tag")),
    drTrg_Cleaning_(cfg.getParameter<double>("drTrg_Cleaning")),
    dcaSig_(cfg.getParameter<double>("dcaSig")),
    trkNormChiMin_(cfg.getParameter<int>("trkNormChiMin")),
    trkNormChiMax_(cfg.getParameter<int>("trkNormChiMax")) 
{
  // removed until we see the plot
//    produces<pat::CompositeCandidateCollection>("TagSide");
    produces<pat::CompositeCandidateCollection>("SelectedTracks");  
}

  ~TrackMerger() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  std::pair<double,double> computeDCA(const pat::PackedCandidate &pfCand,
				      edm::ESHandle<MagneticField> bFieldHandle,
				      const GlobalPoint& refP) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:

  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;

  //selections                                                                 
  const double trkPtCut_;                const double trkEtaCut_;
  const double dzTrg_cleaning_;          const double drTrg_Cleaning_;
  // untill the study for this finalized, perhaps we keep it as comment. Then we uncomment it
 // const double drTrg_ProbeCleaning_;   const double drTrg_TagCleaning_;
 // const double dcaSig_probe_;         const double dcaSig_tag_;
  const double dcaSig_;                  const int trkNormChiMin_;
  const int trkNormChiMax_;
};






void TrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {

  //input
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

 // do we need all those ints?
//  int nLostTracks = lostTracks->size();      
//  int totalTracks = nTracks + nLostTracks;
 //same as before comment out untill we are sure 
//  std::unique_ptr<pat::CompositeCandidateCollection> outTag(new pat::CompositeCandidateCollection());

//ok this was CompositeCandidateCollection 
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out(new pat::CompositeCandidateCollection);


//correct logic but a bit convoluted -> changing to smthn simpler
 std::vector<pat::PackedCandidate> totalTracks(*tracks);

 totalTracks.insert(totalTracks.end(),lostTracks->begin(),lostTracks->end());

 
 // for loop is better to be range based - especially for large ensembles  
 for( const pat::PackedCandidate & trk: totalTracks){

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
   GlobalPoint trgvtx;
   for (const pat::Muon & mu: *trgMuons){
    //remove tracks inside trg muons jet
    if(reco::deltaR(trk, mu) < drTrg_Cleaning_ && drTrg_Cleaning_ >0) 
      continue;
    //if dz is negative it is deactivated
    if((fabs(trk.vz() - mu.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0))
       continue;
    skipTrack=false;
    trgvtx=GlobalPoint(mu.vx(),mu.vy(),mu.vz());
    break; // at least for one trg muon to pass this cuts
   }
   // if track is closer to at least a triggering muon keep it
   if (skipTrack) continue;

   // high purity requirment applied only in packedCands
   unsigned int itrk=&trk-&totalTracks[0];
   if( itrk < nTracks && !trk.trackHighPurity()) continue;
   
   //distance closest approach in x,y wrt triggeringMuon
   std::pair<double,double> DCA = computeDCA(trk, bFieldHandle, trgvtx);
   float DCABS = DCA.first;
   float DCABSErr = DCA.second;
   float DCASig = DCABS/DCABSErr;
   if (DCASig >  dcaSig_  && dcaSig_ >0) continue;
   pat::CompositeCandidate pcand;
   
   pcand.setP4(trk.p4());
   pcand.setCharge(trk.charge());
   pcand.setVertex(trk.vertex());
   pcand.addUserInt("isPacked", (itrk < nTracks) ? 1 : 0);
   pcand.addUserInt("isLostTrk", (itrk < nTracks) ? 0 : 1);      
   pcand.addUserFloat("dxy", trk.dxy());
   pcand.addUserFloat("dxyS", trk.dxy()/trk.dxyError());
   pcand.addUserFloat("dz", trk.dz()); 
   pcand.addUserFloat("dzS", trk.dz()/trk.dzError());
   pcand.addUserFloat("DCASig", DCASig);
   //adding the candidate in the composite stuff for fit (need to test)
//   pcand.addUserCand("cand",trk.sourceCandidatePtr(0));
    
   tracks_out->emplace_back(pcand);
 }
 
//evt.put(std::move(outTag), "TagSide");
  evt.put(std::move(tracks_out), "SelectedTracks");
}


std::pair<double,double> TrackMerger::computeDCA(const pat::PackedCandidate &pfCand,
						 edm::ESHandle<MagneticField> bFieldHandle,
						 const GlobalPoint& refP) const 
{
  
  const reco::TransientTrack trackTT((*(pfCand.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint(refP);

  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();

  std::pair<double,double> DCA = std::make_pair(DCABS,DCABSErr);
  return DCA;
}



//define this as a plug-in
DEFINE_FWK_MODULE(TrackMerger);
