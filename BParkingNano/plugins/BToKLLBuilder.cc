#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

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

class BToKLLBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToKLLBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    //isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    isotrkDCACut_(cfg.getParameter<double>("isotrkDCACut")),
    drIso_cleaning_(cfg.getParameter<double>("drIso_cleaning")),
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BToKLLBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  //const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::CompositeCandidate> isotrk_selection_; 
  const double isotrkDCACut_;
  const double drIso_cleaning_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BToKLLBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &iSetup) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  const MagneticField *fMagneticField = fieldHandle.product();
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);

  edm::ESHandle<TransientTrackBuilder> theB ;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  VertexDistance3D a3d;  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  //edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  //evt.getByToken(isolostTracksToken_, iso_lostTracks);

  std::vector<int> used_lep1_id, used_lep2_id, used_trk_id;


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
    edm::Ptr<pat::CompositeCandidate> k_ptr(kaons, k_idx);
    if( !k_selection_(*k_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector k_p4(
      k_ptr->pt(), 
      k_ptr->eta(),
      k_ptr->phi(),
      K_MASS
      );

    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
      int l1_idx = ll_prt->userInt("l1_idx");
      int l2_idx = ll_prt->userInt("l2_idx");
    
      pat::CompositeCandidate cand;
      cand.setP4(ll_prt->p4() + k_p4);
      cand.setCharge(ll_prt->charge() + k_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the lepton passing the corresponding selection
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("K", k_ptr);
      cand.addUserCand("dilepton", ll_prt);

      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("k_idx", k_idx);
    
      auto dr_info = min_max_dr({l1_ptr, l2_ptr, k_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      // TODO add meaningful variables
      
      if( !pre_vtx_selection_(cand) ) continue;
    
      KinVtxFitter fitter(
        {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
        {l1_ptr->mass(), l2_ptr->mass(), K_MASS},
        {LEP_SIGMA, LEP_SIGMA, K_SIGMA} //some small sigma for the lepton mass
        );
      if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );
      used_lep1_id.emplace_back(l1_idx);
      used_lep2_id.emplace_back(l2_idx);
      used_trk_id.emplace_back(k_idx);
      cand.addUserInt("sv_OK" , fitter.success());
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      cand.addUserFloat("sv_prob", fitter.prob());
      cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
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

      cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());
    
      if( !post_vtx_selection_(cand) ) continue;        

      //compute isolation
      float l1_iso03 = 0;
      float l1_iso04 = 0;
      float l2_iso03 = 0;
      float l2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;
      int n_isotrk = 0;

      for(size_t trk_idx = 0; trk_idx < kaons->size(); ++trk_idx) {
        // corss clean kaon
        if (trk_idx == k_idx) continue;
        edm::Ptr<pat::CompositeCandidate> trk_ptr(kaons, trk_idx);

        if( !isotrk_selection_(*trk_ptr) ) continue;

        // cross clean PF (electron and muon)
        unsigned int iTrk = trk_ptr->userInt("keyPacked");
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) {
          continue;
        }

        // cross clean leptons
        // hard to trace the source particles of low-pT electron in B builder
        // use simple dR cut instead
        float dr_to_l1_prefit = deltaR(l1_ptr->eta(), l1_ptr->phi(), trk_ptr->eta(), trk_ptr->phi());
        float dr_to_l2_prefit = deltaR(l2_ptr->eta(), l2_ptr->phi(), trk_ptr->eta(), trk_ptr->phi());
        if ((dr_to_l1_prefit < drIso_cleaning_) || (dr_to_l2_prefit < drIso_cleaning_)) continue;

        TrajectoryStateOnSurface tsos_iso = extrapolator.extrapolate(kaons_ttracks->at(trk_idx).impactPointState(), fitter.fitted_vtx());
        std::pair<bool,Measurement1D> cur3DIP_iso = absoluteImpactParameter(tsos_iso, fitter.fitted_refvtx(), a3d);
        float svip_iso = cur3DIP_iso.second.value();
        if (cur3DIP_iso.first && svip_iso < isotrkDCACut_) {
          // add to final particle iso if dR < cone
          float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk_ptr->eta(), trk_ptr->phi());
          float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk_ptr->eta(), trk_ptr->phi());
          float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk_ptr->eta(), trk_ptr->phi());
          float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk_ptr->eta(), trk_ptr->phi());

          bool isotrk_matched = false;
          if (dr_to_l1 < 0.4){
            isotrk_matched = true;
            l1_iso04 += trk_ptr->pt();
            if ( dr_to_l1 < 0.3) l1_iso03 += trk_ptr->pt();
          }
          if (dr_to_l2 < 0.4){
            isotrk_matched = true;
            l2_iso04 += trk_ptr->pt();
            if (dr_to_l2 < 0.3)  l2_iso03 += trk_ptr->pt();
          }
          if (dr_to_k < 0.4){
            isotrk_matched = true;
            k_iso04 += trk_ptr->pt();
            if (dr_to_k < 0.3) k_iso03 += trk_ptr->pt();
          }
          if (dr_to_b < 0.4){
            isotrk_matched = true;
            b_iso04 += trk_ptr->pt();
            if (dr_to_b < 0.3) b_iso03 += trk_ptr->pt();
          }
          if (isotrk_matched) n_isotrk++;
        }
      }

      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );
      cand.addUserInt("n_isotrk" , n_isotrk );

      ret_val->push_back(cand);
    } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
  } // for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)

  for (auto & cand: *ret_val){
    cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
    cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
    cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
  }

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToKLLBuilder);
