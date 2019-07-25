#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include <TLorentzVector.h>

using namespace std;


float MuonMass_ = 0.10565837;

class TriggeringMuonProducer : public edm::EDProducer {
    
public:
    
    explicit TriggeringMuonProducer(const edm::ParameterSet &iConfig);
    
    ~TriggeringMuonProducer() override {};
    
    
private:

    virtual void produce(edm::Event&, const edm::EventSetup&);


    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    bool debug;
};


TriggeringMuonProducer::TriggeringMuonProducer(const edm::ParameterSet &iConfig):
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) )
{
    produces<pat::CompositeCandidateCollection>();
    debug = false;
}



void TriggeringMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vertexSrc_, vertexHandle);
    const reco::Vertex & PV = vertexHandle->front();

    if(debug) std::cout << " BToKstllProducer::identifyTriggeringMuons " << std::endl;

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    std::vector<std::vector<float>> triggeringMuons;
    //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);
    if(debug) std::cout << "\n TRIGGER OBJECTS " << std::endl;

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackFilterLabels(iEvent, *triggerBits);
      obj.unpackPathNames(names);

      bool isTriggerMuon = false;
      for (unsigned h = 0; h < obj.filterIds().size(); ++h)
	if(obj.filterIds()[h] == 83){ 
	  isTriggerMuon = true; 
	  if(debug) std::cout << "\t   Type IDs:   " << 83;  //83 = muon
	  break;
	} 

      if(!isTriggerMuon) continue; 
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
	std::string filterName = obj.filterLabels()[h];
	if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
	  isTriggerMuon = true;
	  if(debug) std::cout << "\t   Filters:   " << filterName; 
	  break;
	}
	else{ isTriggerMuon = false; }
      }

      if(!isTriggerMuon) continue;
      std::vector<float> localMuon;
      localMuon.push_back(obj.pt());
      localMuon.push_back(obj.eta());
      localMuon.push_back(obj.phi());
      triggeringMuons.push_back(localMuon);
      if(debug){ std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
	// Print trigger object collection and type
	std::cout << "\t   Collection: " << obj.collection() << std::endl;
      }
    }//trigger objects

    if(debug){
      std::cout << "\n total n of triggering muons = " << triggeringMuons.size() << std::endl;
      for(auto ij : triggeringMuons){
	for(auto il : ij) std::cout << " >>> components (pt, eta, phi) = " << il << std::endl;
      }
    }


    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );


    //now check for reco muons matched to triggering muons
    edm::Handle<std::vector<pat::Muon>> muonForTrgHandle;
    iEvent.getByToken(muonSrc_, muonForTrgHandle);

    for(unsigned int iTrg=0; iTrg<muonForTrgHandle->size(); ++iTrg){
      const pat::Muon & muon1 = (*muonForTrgHandle)[iTrg];

      if(!(muon1.isLooseMuon() && muon1.isSoftMuon(PV))) continue;
      TLorentzVector recoMuon(muon1.pt(), muon1.eta(), muon1.phi(), MuonMass_);

      float dRMuonMatching = -1.;
      int muonMatching_index = -1;
      for(unsigned int ij=0; ij<triggeringMuons.size(); ++ij){

	TLorentzVector trgMuon((triggeringMuons.at(ij))[0], (triggeringMuons.at(ij))[1], (triggeringMuons.at(ij))[2], MuonMass_);
	float dR = recoMuon.DeltaR(trgMuon);
	if((dR < dRMuonMatching || dRMuonMatching == -1) && dR < 0.01){
	  dRMuonMatching = dR;
	  muonMatching_index = iTrg;
	  if(debug) std::cout << " dR = " << dR 
			      << " reco = " << recoMuon.Pt() << " " << recoMuon.Eta() << " " << recoMuon.Phi() << " " 
			      << " HLT = " << trgMuon.Pt() << " " << trgMuon.Eta() << " " << trgMuon.Phi()
			      << std::endl;
	}
      }

      //save reco muon 
      // can add p4 of triggering muon in case
      if(muonMatching_index != -1){
	pat::CompositeCandidate recoTriggerMuonCand;
	recoTriggerMuonCand.addDaughter( (*muonForTrgHandle)[muonMatching_index] );
	recoTriggerMuonCand.addUserInt("recoMuonIndex", muonMatching_index);
	result->push_back(recoTriggerMuonCand);
      }
    }
    
    iEvent.put(std::move(result));
}



DEFINE_FWK_MODULE(TriggeringMuonProducer);
