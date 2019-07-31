#include "MuonProducer.h"

MuonProducer::MuonProducer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  trgMuonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("trgMuon")))

{
  // get params
     runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
     SkipNoMuEvt=runParameters.getParameter<bool>("SkipNoMuEvt");
     MuPtCut=runParameters.getParameter<double>("MuPtCut");    
     MuDzCut=runParameters.getParameter<double>("MuDzCut");
     MuEtaCut=runParameters.getParameter<double>("MuEtaCut");
     MuSoftQ=runParameters.getParameter<bool>("MuSoftQ");

     //produced data
     produces<std::vector<pat::Muon>>("myMuons");
    
  }

MuonProducer::~MuonProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

bool
MuonProducer::filter( edm::Event& iEvent, edm::EventSetup const& iSetup)
{
  //get data  
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   edm::Handle<pat::MuonCollection> trgMuons;
   iEvent.getByToken(trgMuonToken_, trgMuons);
   reco::Vertex goodvtx;
   for (const reco::Vertex &vtx : *vertices) {
     if (  vtx.isFake() || !vtx.isValid() ) continue;
     goodvtx=vtx;
   }
  

  std::unique_ptr<pat::MuonCollection> muons_out(new pat::MuonCollection);
 
  for ( auto mu : *muons ){
    if ( mu.pt() < MuPtCut ) continue;
    if ( fabs(mu.eta()) > MuEtaCut ) continue;   
    if ( mu.isSoftMuon(goodvtx) && MuSoftQ ) continue;
    bool dztrgSkip=true;
    for (const pat::Muon & trg : *trgMuons){
      if (fabs(trg.vz() - mu.vz())>MuDzCut) continue;
      dztrgSkip=false; break;
    }
    if (dztrgSkip) continue;
    muons_out->push_back(mu);
  } 

  if (SkipNoMuEvt && muons_out->size()==0) return false;   
  iEvent.put(std::move(muons_out),"myMuons");
  
  return true;
}


void MuonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonProducer);
