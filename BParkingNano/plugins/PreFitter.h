// short classto create pairs and make tracks transient. It applys also cuts on leading object

#ifndef PREFITTER_H
#define PREFITTER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"


//simple class to prepare objects for fitting
class PreFitter{

////class members ////////////////////
//CreateTransientPair: creates a pair of OS transient tracks from any object
//CreateTransientObjects: create a collection of transient tracks given a collection of objects
private:

 std::vector< std::pair<unsigned int, unsigned int> > vindex;
 edm::ESHandle<MagneticField> bField;

/////////////////////////////////CODE//////////////////////////////////////
//written in such way so we can use the template

public:
 PreFitter( edm::ESHandle< MagneticField >& bFieldHandle ):
  bField( bFieldHandle ) {}
 
 ~PreFitter() {};

 template< typename T >
 std::pair< std::vector<reco::TransientTrack>, std::vector<reco::TransientTrack> > CreateTransientPairs( const std::vector< T >& objs, 
                                double ptLep1Cut_ = 0, double bdtEl1Cut_ = 0,
                                double mllMin_    = 0, double mllMax_ = 100,
                                double dzLepLepCut_ = 0
                                )
{

  vindex.clear();       
  // object to return
  std::vector< reco::TransientTrack > goodttrk1, goodttrk2;

  for( typename std::vector< T >::const_iterator obj1 = objs.begin(); 
       obj1 != objs.end(); ++obj1 ){
    for( typename std::vector< T >::const_iterator obj2 = obj1 + 1; 
       obj2 != objs.end(); ++obj2 ){
       // os requirment
      if ( obj1->charge() == obj2->charge() ) continue;
      // leading pT cut
      if ( obj1->pt() < ptLep1Cut_ && obj2->pt() < ptLep1Cut_ ) continue;
      // dz between two objects
      if ( fabs( obj1->vz() - obj2->vz() ) > dzLepLepCut_ ) continue;
      // check if we run on electrons and apply bdt cut on the leading
      if ( obj1->hasUserFloat( "unBiased" ) && obj2->hasUserFloat( "unBiased" ) ) 
         {
          if ( ( obj1->userFloat( "unBiased" ) < bdtEl1Cut_ || 
                 obj1->pt() < ptLep1Cut_ ) && 
               ( obj2->userFloat( "unBiased" ) < bdtEl1Cut_ || 
                 obj2->pt() < ptLep1Cut_ )
               )
                continue;  
         }      
      float mass = ( obj1->p4() + obj2->p4() ).M();
      // check mass
      if ( mass < mllMin_ || mass > mllMax_ ) continue;

      // save good pairs
      vindex.emplace_back( std::make_pair( obj1 - objs.begin(), 
                                           obj2 - objs.begin()
                                           ) );      

      goodttrk1.emplace_back( reco::TransientTrack( *( obj1->bestTrack() ), 
                                                    &( *bField )
                                                    ) );

      goodttrk2.emplace_back( reco::TransientTrack( *( obj2->bestTrack() ), 
                                                    &( *bField )
                                                    ) );
      
    }
  }  
  return std::move( std::make_pair( goodttrk1, goodttrk2 ) );
}


  // make an object collection transient
template< typename T >
std::vector< reco::TransientTrack > CreateTransientObject( 
                                                 const std::vector< T >& objs)
{
  std::vector< reco::TransientTrack > ttks;
  for( const T& obj: objs)
      ttks.emplace_back(reco::TransientTrack( 
                       *( obj.userCand( "cand" )->bestTrack() ), &( *bField )
                        ) );
    

  return std::move(ttks);
}

// return indices
std::vector< std::pair<unsigned int, unsigned int> > GetIndexPairs()
{
    return vindex;
}

};

#endif
