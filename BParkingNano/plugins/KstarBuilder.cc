///////////////////////// Code to produce K* candidates ////////////////////////


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




class KstarBuilder : public edm::global::EDProducer<> {

  
public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  
  explicit KstarBuilder(const edm::ParameterSet &cfg):
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("pfcands") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )} {

      //output
       produces<pat::CompositeCandidateCollection>();

    }

  ~KstarBuilder() override {}
  
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


void KstarBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //inputs  
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);
 

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> kstar_out(new pat::CompositeCandidateCollection());

  

  // main loop
  for(size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx ){

    edm::Ptr<pat::CompositeCandidate> trk1_ptr( pfcands, trk1_idx );
    if(!trk1_selection_(*trk1_ptr)) continue; 
    
    for(size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx) {

     edm::Ptr<pat::CompositeCandidate> trk2_ptr( pfcands, trk2_idx );
     if (trk1_ptr->charge() == trk2_ptr->charge()) continue; 
     if(!trk2_selection_(*trk2_ptr)) continue;
          
     // create a K* candidate; add first quantities that can be used for pre fit selection
     pat::CompositeCandidate kstar_cand;
     auto trk1_p4=trk1_ptr->polarP4();
     auto trk2_p4=trk2_ptr->polarP4();
     trk1_p4.SetM(K_MASS);
     trk2_p4.SetM(PI_MASS);

     //adding stuff for pre fit selection
     kstar_cand.setP4(trk1_p4 + trk2_p4);
     kstar_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));

     // save indices
     kstar_cand.addUserInt("trk1_idx", trk1_idx );
     kstar_cand.addUserInt("trk2_idx", trk2_idx );

     // save cands      
     kstar_cand.addUserCand("trk1", trk1_ptr );
     kstar_cand.addUserCand("trk2", trk2_ptr );

     //second mass hypothesis
     trk1_p4.SetM(PI_MASS);
     trk2_p4.SetM(K_MASS);
     kstar_cand.addUserFloat("barMass", (trk1_p4 + trk2_p4).M() );
     
     // selection before fit
     if( !pre_vtx_selection_(kstar_cand) ) continue;
           
     KinVtxFitter fitter(
       {ttracks->at(trk1_idx), ttracks->at(trk2_idx)},
       { K_MASS, PI_MASS },
       {K_SIGMA, K_SIGMA} //K and PI sigma equal...
        );
      if ( !fitter.success() ) continue;           

      // save quantities after fit
      kstar_cand.addUserFloat("sv_chi2", fitter.chi2());
      kstar_cand.addUserFloat("sv_ndof", fitter.dof()); 
      kstar_cand.addUserFloat("sv_prob", fitter.prob());    
      kstar_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
      kstar_cand.addUserFloat("fitted_pt", fitter.fitted_candidate().globalMomentum().perp() );
      kstar_cand.addUserFloat("fitted_eta", fitter.fitted_candidate().globalMomentum().eta() );
      kstar_cand.addUserFloat("fitted_phi", fitter.fitted_candidate().globalMomentum().phi() );

      // second mass hypothesis
      auto fitted_trk1= fitter.daughter_p4(0);
      auto fitted_trk2= fitter.daughter_p4(1);
      fitted_trk1.SetM(PI_MASS);
      fitted_trk2.SetM(K_MASS);
      kstar_cand.addUserFloat("fitted_barMass", (fitted_trk1+fitted_trk2).M() );
                    
      // after fit selection
      if( !post_vtx_selection_(kstar_cand) ) continue;
      kstar_out->emplace_back(kstar_cand);
      }
  }
  
  evt.put(std::move(kstar_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KstarBuilder);

