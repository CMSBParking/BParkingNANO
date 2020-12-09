#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class ConstrVtxBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit ConstrVtxBuilder(const edm::ParameterSet &cfg):
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    bmesons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("bmesons") )},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    jpsi_low_{cfg.getParameter<double>("jpsiLow")},
    jpsi_up_{cfg.getParameter<double>("jpsiUp")},
    psi2s_up_{cfg.getParameter<double>("psi2sUp")}
    {
      produces<pat::CompositeCandidateCollection>();
    }

  ~ConstrVtxBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> bmesons_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
  const double jpsi_low_;
  const double jpsi_up_;
  const double psi2s_up_;
};

void ConstrVtxBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<pat::CompositeCandidateCollection> bmesons;
  evt.getByToken(bmesons_, bmesons);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  KinVtxFitter *constr_fitter = 0;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  for(size_t b_idx = 0; b_idx < bmesons->size(); ++b_idx) {
    edm::Ptr<pat::CompositeCandidate> b_ptr(bmesons, b_idx);
    edm::Ptr<reco::Candidate> l1_ptr = b_ptr->userCand("l1");
    edm::Ptr<reco::Candidate> l2_ptr = b_ptr->userCand("l2");
    edm::Ptr<reco::Candidate> k_ptr = b_ptr->userCand("K");

    int l1_idx = b_ptr->userInt("l1_idx");
    int l2_idx = b_ptr->userInt("l2_idx");
    int k_idx = b_ptr->userInt("k_idx");

    pat::CompositeCandidate cand;

    cand.setP4(b_ptr->p4());
    cand.setCharge(b_ptr->charge());
    cand.addUserFloat("original_fitted_mass", b_ptr->userFloat("fitted_mass"));
    cand.addUserFloat("original_fitted_mll", b_ptr->userFloat("fitted_mll"));
    cand.addUserInt("b_idx", b_idx);

    if( !pre_vtx_selection_(cand) ) continue;

    bool constr_sv_OK = false;
    double constr_mass = 0;
    if ((cand.userFloat("original_fitted_mll") > jpsi_low_) && (cand.userFloat("original_fitted_mll") < jpsi_up_)) { 
      constr_mass = JPSI_MASS;
    } else if ((cand.userFloat("original_fitted_mll") > jpsi_up_) && (cand.userFloat("original_fitted_mll") < psi2s_up_)) {
      constr_mass = PSI2S_MASS;
    }
    
    if (constr_mass != 0) {
      constr_fitter = new KinVtxFitter(
          {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
          {l1_ptr->mass(), l2_ptr->mass(), K_MASS},
          {LEP_SIGMA, LEP_SIGMA, K_SIGMA}, //some small sigma for the lepton mass
          constr_mass
          );
      constr_sv_OK = constr_fitter->success();
    
      if (constr_sv_OK) {
        cand.addUserInt("sv_OK" , constr_fitter->success());
        cand.addUserFloat("sv_chi2", constr_fitter->chi2());
        cand.addUserFloat("sv_ndof", constr_fitter->dof()); // float??
        cand.addUserFloat("sv_prob", constr_fitter->prob());
        cand.addUserFloat("fitted_mll" , (constr_fitter->daughter_p4(0) + constr_fitter->daughter_p4(1)).mass());

        auto fit_p4 = constr_fitter->fitted_p4();
        cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
        cand.addUserFloat("fitted_eta" , fit_p4.eta());
        cand.addUserFloat("fitted_phi" , fit_p4.phi());
        cand.addUserFloat("fitted_mass", constr_fitter->fitted_candidate().mass());      
        cand.addUserFloat("fitted_massErr", sqrt(constr_fitter->fitted_candidate().kinematicParametersError().matrix()(6,6)));      
        cand.addUserFloat(
          "fitted_cos_theta_2D", 
          cos_theta_2D(*constr_fitter, *beamspot, fit_p4)
          );
        auto lxy = l_xy(*constr_fitter, *beamspot);
        cand.addUserFloat("l_xy", lxy.value());
        cand.addUserFloat("l_xy_unc", lxy.error());
        cand.addUserFloat("vtx_x", constr_fitter->fitted_vtx().x());
        cand.addUserFloat("vtx_y", constr_fitter->fitted_vtx().y());
        cand.addUserFloat("vtx_z", constr_fitter->fitted_vtx().z());
        cand.addUserFloat("vtx_ex", sqrt(constr_fitter->fitted_vtx_uncertainty().cxx()));
        cand.addUserFloat("vtx_ey", sqrt(constr_fitter->fitted_vtx_uncertainty().cyy()));
        cand.addUserFloat("vtx_ez", sqrt(constr_fitter->fitted_vtx_uncertainty().czz()));

        
        //cand.addUserFloat("fitted_l1_pt" , constr_fitter->daughter_p4(0).pt()); 
        //cand.addUserFloat("fitted_l1_eta", constr_fitter->daughter_p4(0).eta());
        //cand.addUserFloat("fitted_l1_phi", constr_fitter->daughter_p4(0).phi());
        //cand.addUserFloat("fitted_l2_pt" , constr_fitter->daughter_p4(1).pt()); 
        //cand.addUserFloat("fitted_l2_eta", constr_fitter->daughter_p4(1).eta());
        //cand.addUserFloat("fitted_l2_phi", constr_fitter->daughter_p4(1).phi()); 
        cand.addUserFloat("fitted_k_pt"  , constr_fitter->daughter_p4(2).pt()); 
        cand.addUserFloat("fitted_k_eta" , constr_fitter->daughter_p4(2).eta());
        cand.addUserFloat("fitted_k_phi" , constr_fitter->daughter_p4(2).phi());
        

      }
      delete constr_fitter;
    }
    
    if (! constr_sv_OK) continue;
    if( !post_vtx_selection_(cand) ) continue;   

    ret_val->push_back(cand);

  } // for(size_t b_idx = 0; b_idx < bmesons->size(); ++b_idx)

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ConstrVtxBuilder);
