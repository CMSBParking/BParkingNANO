#include "BKllProducer.h"

BKllProducer::BKllProducer(const edm::ParameterSet& iConfig):
   muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>  ("muons"))),
   electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>  ("electrons"))),
   tracksToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
   beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot")))
{
   //produced data
   runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
   SkipNoRecoBEvt=runParameters.getParameter<bool>("SkipNoRecoBEvt");
   LeptonPdgId=runParameters.getParameter<int>("pdgId");
   name=runParameters.getParameter<std::string>("name");
   if (name=="default" && LeptonPdgId==13 ) name="BKmumu";
   else if (name=="default" && LeptonPdgId==11 ) name="BKee";
   produces<std::vector<pat::CompositeCandidate>>(name);
}


BKllProducer::~BKllProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


bool
BKllProducer::filter(edm::Event& iEvent, edm::EventSetup const& iSetup)
{ 
   edm::Handle<std::vector<pat::Electron>> electrons;
   iEvent.getByToken(electronsToken_, electrons); 
   edm::Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken(muonsToken_,muons);
   edm::Handle<pat::CompositeCandidateCollection> tracks;
   iEvent.getByToken(tracksToken_, tracks);
   edm::Handle<reco::BeamSpot> Bspot;
   iEvent.getByToken(beamSpotToken_,Bspot);
   edm::ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

   //make pairs and transient tracks to be faster 
   //create all transient tracks in advance of fit
   if (LeptonPdgId==13){
       BtoXLLTemplateFitter<pat::Muon> Kmumu(*muons,*tracks,*Bspot, 
                                              runParameters,bFieldHandle);

       std::unique_ptr<std::vector<pat::CompositeCandidate>> recoB(
                                                    Kmumu.ReconstructB(K_MASS));

       if (recoB->size()==0){ 
         std::unique_ptr<std::vector<pat::CompositeCandidate>> recoB(new std::vector<pat::CompositeCandidate>());
         iEvent.put(std::move(recoB),name);
         return true;
       } else{
         iEvent.put(std::move(recoB),name);
         return true;
       }
   } else if (LeptonPdgId==11){
       BtoXLLTemplateFitter<pat::Electron> Kee(*electrons,*tracks,*Bspot, 
                                                runParameters,bFieldHandle);

       std::unique_ptr<std::vector<pat::CompositeCandidate>> recoB(
						      Kee.ReconstructB(K_MASS));

       if (recoB->size()==0) 
          return false;
       else{
	  iEvent.put(std::move(recoB),name);
          return true;
       }  
   } else{
       std::cout<<"Lepton can be ONLY 13 or 11. Please fix "<<std::endl;
       return true;
   }
}
void BKllProducer::beginJob(){}

void BKllProducer::endJob() {}

void BKllProducer::beginRun( edm::Run const& run,  edm::EventSetup const& iSetup) {}

void BKllProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BKllProducer);
