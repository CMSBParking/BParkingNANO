//// table to produce hlt bits that we are going to use in the analysis


// system include files
#include <memory>


#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "TString.h"
#include <string>



class TrgBitTableProducer : public edm::stream::EDProducer<> {


public:

 
 explicit TrgBitTableProducer(const edm::ParameterSet &cfg):
    trgresultsToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag> ("triggerresults"))),
    hltpaths_( cfg.getParameter< std::vector<std::string> >( "paths" ) )
{
    produces<nanoaod::FlatTable>();

}

  ~TrgBitTableProducer() override {}

  void produce(edm::Event&, edm::EventSetup const&) override;


  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:

  const edm::EDGetTokenT< edm::TriggerResults > trgresultsToken_;
  // l1 seeds not implemented yet but can be added with litle effort
  const std::vector< std::string > l1seeds_;              
  const std::vector< std::string > hltpaths_;
 

};



void 
TrgBitTableProducer::produce( edm::Event &evt, edm::EventSetup const &stp) 
{

  //input
  edm::Handle< edm::TriggerResults > trigResults;
  evt.getByToken( trgresultsToken_, trigResults);
  // returns uint8 instead of bool, because addCollumnValue<bool> is unsuported. (cmsRun error). the next "economical class is uint8
  std::vector<uint8_t> bits;
  unsigned int Npaths = hltpaths_.size();
  bits.reserve( Npaths );
  edm::TriggerNames trigName = evt.triggerNames( *trigResults );   
 

  if ( trigResults.failedToGet() ){
    for ( unsigned int ibit = 0; ibit < Npaths; ++ibit){
       bits.push_back( 0 );
     }
  } else {
    int Ntrg = trigResults->size();
    for ( auto& hltpath: hltpaths_ ){
      bool fire = false; 
      for( int itrg = 0; itrg < Ntrg; ++itrg ){
        if ( !trigResults->accept( itrg ) ) continue;
        TString TrigPath = trigName.triggerName( itrg );
        if ( TrigPath.Contains( hltpath ) ) fire=true;      
      } 
      
      if( fire ) bits.push_back( 1 );
      else bits.push_back( 0 );
    }
  }
 
 
  auto tab  = std::make_unique<nanoaod::FlatTable>(1,"", true);
  for (unsigned int ipath = 0; ipath <Npaths; ++ipath ){
    tab->addColumnValue<uint8_t> (hltpaths_[ipath], bits[ipath], "hlt path", nanoaod::FlatTable::UInt8Column);
  }
  

    evt.put(std::move(tab));

}


//define this as a plug-in
DEFINE_FWK_MODULE(TrgBitTableProducer);
