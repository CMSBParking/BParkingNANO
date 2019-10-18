///////////////////////// Code to produce phi candidates ////////////////////////


#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"




class PhiBuilder : public edm::global::EDProducer<> {

  
public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit PhiBuilder(const edm::ParameterSet &cfg):
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pfcands") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )} {

      //output
       produces<pat::CompositeCandidateCollection>();

    }

  ~PhiBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> trk1_selection_; // cuts on leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk2_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_; //input PF cands this is sorted in pT in previous step
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_; //input TTracks of PF cands
};


void PhiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //inputs  
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);
 

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> phi_out(new pat::CompositeCandidateCollection());

  

  // main loop
  for(size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx ){

    edm::Ptr<pat::CompositeCandidate> trk1_ptr( pfcands, trk1_idx );
    if(!trk1_selection_(*trk1_ptr)) continue; 
    
    for(size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx) {

     edm::Ptr<pat::CompositeCandidate> trk2_ptr( pfcands, trk2_idx );
     if (trk1_ptr->charge() == trk2_ptr->charge()) continue; 
     if(!trk2_selection_(*trk2_ptr)) continue;
          
     // create a K* candidate; add first quantities that can be used for pre fit selection
     pat::CompositeCandidate phi_cand;
     auto trk1_p4=trk1_ptr->polarP4();
     auto trk2_p4=trk2_ptr->polarP4();
     trk1_p4.SetM(K_MASS);
     trk2_p4.SetM(K_MASS);

     //adding stuff for pre fit selection
     phi_cand.setP4(trk1_p4 + trk2_p4);
     phi_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));

     // save indices
     phi_cand.addUserInt("trk1_idx", trk1_idx );
     phi_cand.addUserInt("trk2_idx", trk2_idx );

     // save cands      
     phi_cand.addUserCand("trk1", trk1_ptr );
     phi_cand.addUserCand("trk2", trk2_ptr );

     // selection before fit
     if( !pre_vtx_selection_(phi_cand) ) continue;
           
     KinVtxFitter fitter(
       {ttracks->at(trk1_idx), ttracks->at(trk2_idx)},
       { K_MASS, K_MASS },
       {K_SIGMA, K_SIGMA}
        );
      if ( !fitter.success() ) continue;           

      // save quantities after fit
      phi_cand.addUserFloat("sv_chi2", fitter.chi2());
      phi_cand.addUserFloat("sv_ndof", fitter.dof()); 
      phi_cand.addUserFloat("sv_prob", fitter.prob());    
      phi_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
      phi_cand.addUserFloat("fitted_pt", fitter.fitted_candidate().globalMomentum().perp() );
      phi_cand.addUserFloat("fitted_eta", fitter.fitted_candidate().globalMomentum().eta() );
      phi_cand.addUserFloat("fitted_phi", fitter.fitted_candidate().globalMomentum().phi() );

      // after fit selection
      if( !post_vtx_selection_(phi_cand) ) continue;
      phi_out->emplace_back(phi_cand);
      }
  }
  
  evt.put(std::move(phi_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhiBuilder);

