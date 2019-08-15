#include "BKllProducer.h"



BKllProducer::BKllProducer( const edm::ParameterSet& iConfig ):

  muonsToken_( consumes<std::vector<pat::Muon>>( iConfig.getParameter
                                            <edm::InputTag>("muons") ) ),
  electronsToken_( consumes<std::vector<pat::Electron>>( iConfig.getParameter
                                            <edm::InputTag>("electrons") ) ),
  tracksToken_( consumes<pat::CompositeCandidateCollection>(iConfig.getParameter
                                            <edm::InputTag>("tracks") ) ),
  beamSpotToken_( consumes<reco::BeamSpot>(iConfig.getParameter 
                                            <edm::InputTag>("beamSpot") ) )
{
   // B parameters  - explanations in comments
   runParameters = iConfig.getParameter<edm::ParameterSet>( "RunParameters" );

   // If an event has no reconstructed B then skip it
   SkipNoRecoBEvt = runParameters.getParameter<bool>( "SkipNoRecoBEvt" );

   // PdgID of the lepton in final state (13 or 11)
   LeptonPdgId = runParameters.getParameter<int>( "FinalLeptonId" );

   // User defined output label
   name = runParameters.getParameter<std::string>( "label" );

   // Default laels
   if ( name=="default" && LeptonPdgId==13 ) name = "BKmumu";
   else if ( name=="default" && LeptonPdgId==11 ) name = "BKee";
   produces<std::vector<pat::CompositeCandidate>>( name );
}


BKllProducer::~BKllProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


void
BKllProducer::produce( edm::Event& iEvent, edm::EventSetup const& iSetup )
{ 
   // get needed objects
   edm::Handle<std::vector<pat::Electron>> electrons;
   iEvent.getByToken( electronsToken_, electrons ); 
   edm::Handle<std::vector<pat::Muon>> muons;
   iEvent.getByToken( muonsToken_, muons );
   edm::Handle<pat::CompositeCandidateCollection> tracks;
   iEvent.getByToken( tracksToken_, tracks );

   edm::Handle<reco::BeamSpot> Bspot;
   iEvent.getByToken( beamSpotToken_, Bspot );
   edm::ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get( bFieldHandle );

    
  // depending on the lepton flavour pass the corresponding data for fitting 
  if ( LeptonPdgId == 13 ) {
     // Call to the templated class function for B -> XLL state . This way same code for is used for mu and e
     BtoXLLTemplateFitter<pat::Muon> Kmumu( *muons, *tracks, *Bspot, 
                                           runParameters, bFieldHandle );
    
     // Reconstruct B; needs mass hypothesis for hadron track as input
     std::unique_ptr< std::vector<pat::CompositeCandidate> > recoB(
                                           Kmumu.ReconstructB( K_MASS ) );
 
     iEvent.put( std::move( recoB ), name );


  } else if ( LeptonPdgId == 11 ){

    BtoXLLTemplateFitter<pat::Electron> Kee( *electrons, *tracks, *Bspot, 
                                                runParameters, bFieldHandle );

    std::unique_ptr< std::vector<pat::CompositeCandidate> > recoB(
						Kee.ReconstructB( K_MASS ) );

    iEvent.put( std::move( recoB ), name );
   
  } else {
       // If pdg id is not 11 or 13 print warning
       std::cout<<"Lepton can be ONLY 13(mu) or 11(e). B will not be reconstructed. Please fix "<<std::endl;
       std::unique_ptr< std::vector<pat::CompositeCandidate> > recoB ( new std::vector<pat::CompositeCandidate> () );

       iEvent.put( std::move( recoB ), name );
   }
}

void BKllProducer::beginJob(){}

void BKllProducer::endJob() {}

void BKllProducer::beginRun( edm::Run const& run, edm::EventSetup const& iSetup) {}

void BKllProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );

}

//define this as a plug-in
DEFINE_FWK_MODULE( BKllProducer );
