#ifndef TRACKPRODUCER_H
#define TRACKPRODUCER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "helper.h"



class TrackProducer : public edm::one::EDFilter<edm::one::SharedResources>  {

public:
  //constractor
   explicit TrackProducer(const edm::ParameterSet&);
   //destructor
   ~TrackProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
 
  bool filter(edm::Event&, const edm::EventSetup&);

  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_; 
  //selections
  double TrkPtCut = 0;        double TrkEtaCut = 0;      double TrkDzCut  = 0; 
  int TrkNormChiMin = 0;      bool SkipNoTrkEvt = false; int TrkNormChiMax = 0; 
  edm::ParameterSet runParameters;

};
#endif
