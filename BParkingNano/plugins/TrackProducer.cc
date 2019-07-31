#include "TrackProducer.h"

TrackProducer::TrackProducer(const edm::ParameterSet& iConfig):
  tracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  lostTracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
  trgMuonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("trgMuon")))

{
  // get params
     runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
     SkipNoTrkEvt=runParameters.getParameter<bool>("SkipNoTrkEvt");
     TrkPtCut=runParameters.getParameter<double>("TrkPtCut");    
     TrkDzCut=runParameters.getParameter<double>("TrkDzCut");
     TrkEtaCut=runParameters.getParameter<double>("TrkEtaCut");
     TrkNormChiMin=runParameters.getParameter<int>("TrkNormChiMin");
     TrkNormChiMax=runParameters.getParameter<int>("TrkNormChiMax");

     //produced data
     produces<std::vector<pat::PackedCandidate>>("myTracks");
    
  }

TrackProducer::~TrackProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

bool
TrackProducer::filter( edm::Event& iEvent, edm::EventSetup const& iSetup)
{
  //get data  
   edm::Handle<pat::PackedCandidateCollection> tracks;
   iEvent.getByToken(tracksToken_, tracks);
   edm::Handle<pat::PackedCandidateCollection> lostTracks;
   iEvent.getByToken(lostTracksToken_, lostTracks);
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(trgMuonToken_, muons);
      
  // cmssw has a bug in lost high purity flag of lost tracks. keep the number of pf cands and apply high quality only there	
  unsigned int nPFcands = tracks->size();
  //add lost + cands
  std::vector<pat::PackedCandidate> totalTracks(*tracks);
  totalTracks.insert(totalTracks.end(),lostTracks->begin(),lostTracks->end());

  std::unique_ptr<pat::PackedCandidateCollection> tracks_out(new pat::PackedCandidateCollection);
 
  for ( auto trk : totalTracks ){
    unsigned int itrk(&trk-&totalTracks[0]);
    if (!trk.trackHighPurity() && itrk<nPFcands) continue;
    if (!trk.hasTrackDetails()) continue;
    if (fabs(trk.pdgId())!=211 ) continue;
    if ( trk.pt() < TrkPtCut ) continue;
    if ( fabs(trk.eta()) > TrkEtaCut ) continue;   
    if (trk.pseudoTrack().normalizedChi2()<TrkNormChiMin ||
        trk.pseudoTrack().normalizedChi2()>TrkNormChiMax ) continue;
    bool dzmuSkip=true;
    for (const pat::Muon & mu : *muons){
      if (fabs(mu.vz() - trk.vz())>TrkDzCut) continue;
      dzmuSkip=false; break;
    }
    if (dzmuSkip) continue;
    trk.setMass(K_MASS);
    tracks_out->push_back(trk);
  } 
  if (SkipNoTrkEvt && tracks_out->size()==0) return false;   
  iEvent.put(std::move(tracks_out),"myTracks");
  
  return true;
}


void TrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackProducer);
