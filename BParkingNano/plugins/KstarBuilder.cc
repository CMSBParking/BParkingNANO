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

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<pat::PackedCandidate> CandCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector< std::pair<reco::TransientTrack, reco::TransientTrack> > TransientKstarTrackCollection;


  explicit KstarBuilder(const edm::ParameterSet &cfg):
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pfcands_{consumes<CandCollection>( cfg.getParameter<edm::InputTag>("pfcands") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )} {
       produces<pat::CompositeCandidateCollection>();
       produces< TransientKstarTrackCollection >();
    }

  ~KstarBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::PackedCandidate> trk1_selection_; // cuts on leading cand
  const StringCutObjectSelector<pat::PackedCandidate> trk2_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<CandCollection> pfcands_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_;
};

void KstarBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<CandCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> kstar_out(new pat::CompositeCandidateCollection());

  std::unique_ptr< TransientKstarTrackCollection > kstar_tt_out(new TransientKstarTrackCollection());

  //two mass hypothesis for tracks; we cannot discriminate which trak is K and which is pi
  std::vector<std::pair<float, float>> mhypothesis{ std::make_pair(K_MASS,PI_MASS), std::make_pair(PI_MASS,K_MASS)  };

  for(size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx) {
    edm::Ptr<pat::PackedCandidate> trk1_ptr(pfcands, trk1_idx);
    if(!trk1_selection_(*trk1_ptr)) continue; 
    
    for(size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx) {
      edm::Ptr<pat::PackedCandidate> trk2_ptr(pfcands, trk2_idx);
      if(!trk2_selection_(*trk2_ptr)) continue;
      if (trk1_ptr->charge()== trk2_ptr->charge()) continue;
      
      std::vector <float> trk1_mass;
      std::vector <float> trk2_mass;
      
      for ( auto & masses: mhypothesis){
        // M(K*) selection must be checked for both mass hypothesis
        pat::CompositeCandidate kstar_temp;
        auto mtrk1=masses.first;
        auto mtrk2=masses.second;
        auto trk1_p4=trk1_ptr->polarP4();
        auto trk2_p4=trk2_ptr->polarP4();
        trk1_p4.SetM(mtrk1);
        trk2_p4.SetM(mtrk2);
        kstar_temp.setP4(trk1_p4 + trk2_p4);
        kstar_temp.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));
        if( pre_vtx_selection_(kstar_temp) ){
          trk1_mass.emplace_back(mtrk1);
          trk2_mass.emplace_back(mtrk2);
        }
      }
      if (trk1_mass.size()==0) continue;
      // if at least one hypothesis pass the cuts
      // we now that at least there is one entry
       KinVtxFitter fitter(
        {ttracks->at(trk1_idx), ttracks->at(trk2_idx)},
        { trk1_mass[0],  trk2_mass[0]},
        {K_SIGMA, K_SIGMA} //K and PI sigma are equal... no need to change like masses
        );
      if ( !fitter.success() ) continue;
      
      pat::CompositeCandidate kstar_cand;
      if ( trk1_mass[0]==K_MASS) {
         kstar_cand.addUserInt("trk1isK_trk2isPi",1);
      } else {
         kstar_cand.addUserInt("trk1isK_trk2isPi",0);
         kstar_cand.addUserInt("trk1isPi_trk2isK",1);
      }
      auto trk1_p4=trk1_ptr->polarP4();
      auto trk2_p4=trk2_ptr->polarP4();
      trk1_p4.SetM(trk1_mass[0]);
      trk2_p4.SetM(trk2_mass[0]);
      kstar_cand.setP4(trk1_p4 + trk2_p4);
      kstar_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));

      // save indices
      kstar_cand.addUserInt("trk1_idx", trk1_idx );
      kstar_cand.addUserInt("trk2_idx", trk2_idx );

      // save cands      
      kstar_cand.addUserCand("trk1", trk1_ptr );
      kstar_cand.addUserCand("trk2", trk2_ptr );

      kstar_cand.addUserFloat("sv_chi2", fitter.chi2());
      kstar_cand.addUserFloat("sv_ndof", fitter.dof()); 
      kstar_cand.addUserFloat("sv_prob", fitter.prob());    
      kstar_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );

      //deal with multiple good mass hypothesis
      if ( trk1_mass.size()>1) {
        //now we know that second combination passes the cuts
        trk1_p4.SetM(trk1_mass[1]);
        trk2_p4.SetM(trk2_mass[1]);
        kstar_cand.addUserFloat("bar_mass",(trk1_p4+trk2_p4).M());
        kstar_cand.addUserInt("valid_mhypothesis",trk1_mass.size());  
        kstar_cand.addUserInt("trk1isPi_trk2isK",1);
        auto fitted_trk1= fitter.daughter_p4(0);
        auto fitted_trk2= fitter.daughter_p4(1);
        fitted_trk1.SetM(trk1_mass[1]);
        fitted_trk2.SetM(trk2_mass[1]);
        kstar_cand.addUserFloat("bar_fitted_mass", (fitted_trk1+fitted_trk2).M() );
      }else {
        kstar_cand.addUserFloat("bar_mass",-99);
        kstar_cand.addUserInt("valid_mhypothesis",1); 
        kstar_cand.addUserFloat("bar_fitted_mass", -99 );
        kstar_cand.addUserInt("trk1isPi_trk2isK",0);
      }  
        
      
      // if needed, add here more stuff

      // cut on the SV info
      if( !post_vtx_selection_(kstar_cand) ) continue;
      kstar_out->emplace_back(kstar_cand);
      //produce the transient tracks so not to drag indices all over the code
      kstar_tt_out->emplace_back( std::make_pair(ttracks->at(trk1_idx), ttracks->at(trk2_idx)) );
    }
  }
  
  evt.put(std::move(kstar_out));
  evt.put(std::move(kstar_tt_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KstarBuilder);

