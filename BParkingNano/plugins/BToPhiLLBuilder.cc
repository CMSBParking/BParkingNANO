////////////////////// Code to produce PhiLL candidates /////////////////////////

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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"




class BToPhiLLBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToPhiLLBuilder(const edm::ParameterSet &cfg):
    // selections
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    //inputs
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    phis_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("phis") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    phis_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("phisTransientTracks") )},
    isotracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))},
    isolostTracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))},
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} 
    {
       //output
      produces<pat::CompositeCandidateCollection>();
    }

  ~BToPhiLLBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  // selections
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit


  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> phis_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<TransientTrackCollection> phis_ttracks_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BToPhiLLBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {


  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> phis;
  evt.getByToken(phis_, phis);  
  edm::Handle<TransientTrackCollection> phis_ttracks;
  evt.getByToken(phis_ttracks_, phis_ttracks);   

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  

  for(size_t phi_idx = 0; phi_idx < phis->size(); ++phi_idx) {
    // both phi and lep pair already passed cuts; no need for more preselection
    edm::Ptr<pat::CompositeCandidate> phi_ptr(phis, phi_idx);
    edm::Ptr<reco::Candidate> trk1_ptr= phi_ptr->userCand("trk1");
    edm::Ptr<reco::Candidate> trk2_ptr= phi_ptr->userCand("trk2");
    int trk1_idx = phi_ptr->userInt("trk1_idx");
    int trk2_idx = phi_ptr->userInt("trk2_idx");

    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      edm::Ptr<pat::CompositeCandidate> ll_ptr(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_ptr->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_ptr->userCand("l2");
      int l1_idx = ll_ptr->userInt("l1_idx");
      int l2_idx = ll_ptr->userInt("l2_idx");

      // B0 candidate
      pat::CompositeCandidate cand;
      cand.setP4(ll_ptr->p4() + phi_ptr->p4());
      cand.setCharge( 0 ); //B0 has 0 charge

      // save daughters - unfitted
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("trk1", trk1_ptr);
      cand.addUserCand("trk2", trk2_ptr);
      cand.addUserCand("phi", phi_ptr);
      cand.addUserCand("dilepton", ll_ptr);

      // save indices
      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("trk1_idx", trk1_idx);
      cand.addUserInt("trk2_idx", trk2_idx);
      cand.addUserInt("phi_idx" ,phi_idx);


      auto dr_info = min_max_dr({l1_ptr, l2_ptr, trk1_ptr, trk2_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);


      // check if pass pre vertex cut
      if( !pre_vtx_selection_(cand) ) continue;
        
      KinVtxFitter fitter(
        {phis_ttracks->at(trk1_idx), phis_ttracks->at(trk2_idx), 
         leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx)},
        {K_MASS, K_MASS, l1_ptr->mass(), l2_ptr->mass()},
        { K_SIGMA, K_SIGMA, LEP_SIGMA, LEP_SIGMA}  //K_SIGMA==PI_SIGMA
        );

      if(!fitter.success()) continue; 

      // B0 position
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );

      // vertex vars
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof());
      cand.addUserFloat("sv_prob", fitter.prob());

      // refitted kinematic vars
      cand.addUserFloat("fitted_phi_mass",(fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass() );
      cand.addUserFloat("fitted_phi_pt"  ,(fitter.daughter_p4(0) + fitter.daughter_p4(1)).pt());
      cand.addUserFloat("fitted_phi_eta" ,(fitter.daughter_p4(0) + fitter.daughter_p4(1)).eta());
      cand.addUserFloat("fitted_phi_phi" ,(fitter.daughter_p4(0) + fitter.daughter_p4(1)).phi());
      cand.addUserFloat("fitted_mll"       ,(fitter.daughter_p4(2) + fitter.daughter_p4(3)).mass());

      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fit_p4.mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6))); 

      // refitted daughters (leptons/tracks)     
      std::vector<std::string> dnames{ "trk1", "trk2", "l1", "l2" };
      
      for (size_t idaughter=0; idaughter<dnames.size(); idaughter++){
	cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt" ,fitter.daughter_p4(idaughter).pt() );
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta",fitter.daughter_p4(idaughter).eta() );
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi",fitter.daughter_p4(idaughter).phi() );
      }
      
      // other vars
      cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, cand.p4())
        );
      cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );

      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());
      cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());
      cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));

      // post fit selection
      if( !post_vtx_selection_(cand) ) continue;        
      
      //compute isolation
      float l1_iso03  = 0;
      float l1_iso04  = 0;
      float l2_iso03  = 0;
      float l2_iso04  = 0;
      float trk1_iso03 = 0;
      float trk1_iso04 = 0;
      float trk2_iso03 = 0;
      float trk2_iso04 = 0;
      float b_iso03   = 0;
      float b_iso04   = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon or the pion
        if (trk1_ptr ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        if (trk2_ptr ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1  = deltaR(cand.userFloat("fitted_l1_eta"),  cand.userFloat("fitted_l1_phi"),  trk.eta(), trk.phi());
        float dr_to_l2  = deltaR(cand.userFloat("fitted_l2_eta"),  cand.userFloat("fitted_l2_phi"),  trk.eta(), trk.phi());
        float dr_to_trk1 = deltaR(cand.userFloat("fitted_trk1_eta"),cand.userFloat("fitted_trk1_phi"),trk.eta(), trk.phi());
        float dr_to_trk2 = deltaR(cand.userFloat("fitted_trk2_eta"),cand.userFloat("fitted_trk2_phi"),trk.eta(), trk.phi());
        float dr_to_b   = deltaR(cand.userFloat("fitted_eta"),     cand.userFloat("fitted_phi"),     trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_trk1 < 0.4){
          trk1_iso04 += trk.pt();
          if (dr_to_trk1 < 0.3) trk1_iso03 += trk.pt();
        }
        if (dr_to_trk2 < 0.4){
          trk2_iso04 += trk.pt();
          if (dr_to_trk2 < 0.3) trk2_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      }
      cand.addUserFloat("l1_iso03" , l1_iso03);
      cand.addUserFloat("l1_iso04" , l1_iso04);
      cand.addUserFloat("l2_iso03" , l2_iso03);
      cand.addUserFloat("l2_iso04" , l2_iso04);
      cand.addUserFloat("trk1_iso03", trk1_iso03);
      cand.addUserFloat("trk1_iso04", trk1_iso04);
      cand.addUserFloat("trk2_iso03", trk2_iso03);
      cand.addUserFloat("trk2_iso04", trk2_iso04);
      cand.addUserFloat("b_iso03"  , b_iso03 );
      cand.addUserFloat("b_iso04"  , b_iso04 );
            
      ret_val->push_back(cand);

    } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
   
  } // for(size_t p_idx = 0; p_idx < phis->size(); ++p_idx)
  
  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToPhiLLBuilder);
