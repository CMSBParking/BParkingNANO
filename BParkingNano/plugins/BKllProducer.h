//Class to produce the Kee and Kmumu resonances. It reads selected objects and reconstructs B, cuts can be applied from configuration in before and after fit


#ifndef BKLLPRODUCER_H
#define BKLLPRODUCER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include <vector>
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "BtoXLLTemplateFitter.h"


class BKllProducer : public edm::one::EDProducer<edm::one::SharedResources,edm::one::WatchRuns> {

public:

  //constractor
  explicit BKllProducer( const edm::ParameterSet& );

  //destructor
  ~BKllProducer();
  
  static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );



private:

  virtual void beginJob() override;
  void         beginRun( edm::Run const& iEvent, edm::EventSetup const& );
  virtual void produce( edm::Event&, edm::EventSetup const& ) override;
  void         endRun( edm::Run const& iEvent, edm::EventSetup const& ) {};
  virtual void endJob() override;

  //data
  const edm::EDGetTokenT<pat::MuonCollection>               muonsToken_;
  const edm::EDGetTokenT<pat::ElectronCollection>           electronsToken_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<reco::BeamSpot>                    beamSpotToken_;
  
  //selections - user inputs
  int LeptonPdgId = 0;                    bool SkipNoRecoBEvt = true;
  std::string name = "default";          
  edm::ParameterSet runParameters;

};
#endif
