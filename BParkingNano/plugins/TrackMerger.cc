// Merges the PFPackedCandidates and Lost tracks

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
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "helper.h"


class TrackMerger : public edm::global::EDProducer<> {
public:
  explicit TrackMerger(const edm::ParameterSet &cfg):
    beamSpotSrc_( consumes<reco::BeamSpot> (cfg.getParameter<edm::InputTag>("beamSpot"))),
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgMuonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    trkPtCut_(cfg.getParameter<double>("trkPtCut")),
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
    dzTrg_cleaning_(cfg.getParameter<double>("dzTrg_cleaning")),
    drTrg_cleaning_(cfg.getParameter<double>("drTrg_cleaning")),
    dcaSig_probe_(cfg.getParameter<double>("dcaSig_probe")),
    dcaSig_tag_(cfg.getParameter<double>("dcaSig_tag")),
    trkNormChiMin_(cfg.getParameter<int>("trkNormChiMin")),
    trkNormChiMax_(cfg.getParameter<int>("trkNormChiMax")) 
{
    produces<pat::CompositeCandidateCollection>("TagSide");
    produces<pat::CompositeCandidateCollection>("ProbeSide");
}

  ~TrackMerger() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  std::pair<double,double> computeDCA(const pat::PackedCandidate &pfCand,
				      edm::ESHandle<MagneticField> bFieldHandle,
				      const GlobalPoint& refP) const;

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
  const double drTrg_cleaning_;
  const double dcaSig_probe_;
  const double dcaSig_tag_;
  const int trkNormChiMin_;
  const int trkNormChiMax_;
};

void TrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {
  //get data  
  edm::ESHandle<MagneticField> bFieldHandle;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  stp.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("PFCandProducer") << "No beam spot available from EventSetup" ;
  }
  reco::BeamSpot beamSpot = *beamSpotHandle;

  edm::Handle<pat::PackedCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  evt.getByToken(lostTracksToken_, lostTracks);
  
  edm::Handle<pat::MuonCollection> trgMuon;
  evt.getByToken(trgMuonToken_, trgMuon);

  int nTracks = tracks->size();
  int nLostTracks = lostTracks->size();      
  int totalTracks = nTracks + nLostTracks;
 
  std::unique_ptr<pat::CompositeCandidateCollection> outTag(new pat::CompositeCandidateCollection());
  std::unique_ptr<pat::CompositeCandidateCollection> outProbe(new pat::CompositeCandidateCollection());


  std::vector<int> alreadySaved;
  alreadySaved.resize(totalTracks, 0);

  for(auto muonTrg : *trgMuon) {

    for(int iTrk=0; iTrk<totalTracks; ++iTrk){

      const pat::PackedCandidate* trk;

      if(iTrk < nTracks){
	if(alreadySaved[iTrk]) continue;
	trk = &((*tracks)[iTrk]);
	if(!trk->trackHighPurity()) continue;
	if(abs(trk->pdgId()) != 211 && abs(trk->pdgId()) != 13) continue;
      }
      else{
	if(alreadySaved[iTrk-nTracks]) continue;
	trk = &((*lostTracks)[iTrk-nTracks]);
	if(abs(trk->pdgId()) != 211) continue;
      }
      if(!trk->hasTrackDetails()) continue;
      if(trk->pt() < trkPtCut_ || std::fabs(trk->eta()) > trkEtaCut_) continue;
      if(trk->pseudoTrack().normalizedChi2() < trkNormChiMin_ ||
	 trk->pseudoTrack().normalizedChi2() > trkNormChiMax_ ) continue;
      
      if((std::fabs(trk->vz() - muonTrg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ != -1)) continue;

      bool saved = false;

      double DCABS = -1.;
      double DCABSErr = -1.;

      //probe side DCASig wrt beamspot
      if(reco::deltaR(*trk, muonTrg) > drTrg_cleaning_ || drTrg_cleaning_ == -1){

	std::pair<double,double> DCA = computeDCA(*trk,
						  bFieldHandle,
						  GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
        DCABS = DCA.first;
        DCABSErr = DCA.second;
      }
      else{ //tag side DCASig wrt trigger muon                                                                      
	std::pair<double,double> DCA = computeDCA(*trk,
						  bFieldHandle,
						  GlobalPoint(muonTrg.vx(), muonTrg.vy(), muonTrg.vz()));
	DCABS = DCA.first;
	DCABSErr = DCA.second;
      }

      float DCASig = DCABS/DCABSErr;
      if(DCASig > dcaSig_probe_ || dcaSig_probe_ == -1){
	pat::CompositeCandidate pcand;
	pcand.addDaughter(*trk);
	pcand.addUserInt("isPacked", (iTrk < nTracks) ? 1 : 0);
	pcand.addUserInt("isLostTrk", (iTrk < nTracks) ? 0 : 1);
	pcand.addUserFloat("dxy", trk->dxy());
	pcand.addUserFloat("dxyS", trk->dxy()/trk->dxyError());
	pcand.addUserFloat("dz", trk->dz());
	pcand.addUserFloat("dzS", trk->dz()/trk->dzError());
	pcand.addUserFloat("DCASig", DCASig);
	outProbe->push_back(pcand);
	saved = true;
      }    
      if(DCASig < dcaSig_tag_ || dcaSig_tag_ == -1){
	pat::CompositeCandidate pcand;
	pcand.addDaughter(*trk);
	pcand.addUserInt("isPacked", (iTrk < nTracks) ? 1 : 0);
	pcand.addUserInt("isLostTrk", (iTrk < nTracks) ? 0 : 1);
	pcand.addUserFloat("dxy", trk->dxy());
	pcand.addUserFloat("dxyS", trk->dxy()/trk->dxyError());
	pcand.addUserFloat("dz", trk->dz());
	pcand.addUserFloat("dzS", trk->dz()/trk->dzError());
	pcand.addUserFloat("DCASig", DCASig);
	outTag->push_back(pcand);
	saved = true;
      }
    
      if(saved) alreadySaved[((iTrk < nTracks) ? iTrk : (iTrk - nTracks))] = 1;
    }
  }

  evt.put(std::move(outTag), "TagSide");
  evt.put(std::move(outProbe), "ProbeSide");
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
