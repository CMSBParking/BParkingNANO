#ifndef MUONPRODUCER_H
#define MUONPRODUCER_H

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



class MuonProducer : public edm::one::EDFilter<edm::one::SharedResources>  {

public:
  //constractor
   explicit MuonProducer(const edm::ParameterSet&);
   //destructor
   ~MuonProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
 
  bool filter(edm::Event&, const edm::EventSetup&);
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_; 

  //selections
  double MuPtCut = 0;        double MuEtaCut = 0;      double MuDzCut  = 0; 
  bool SkipNoMuEvt = false;  bool MuSoftQ = false;
  edm::ParameterSet runParameters;

};
#endif
