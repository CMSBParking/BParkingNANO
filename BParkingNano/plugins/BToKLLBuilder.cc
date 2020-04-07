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

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


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
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
    vertex_src_{consumes<reco::VertexCollection>( cfg.getParameter<edm::InputTag>("offlinePrimaryVertexSrc") )}
  {
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
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
  const edm::EDGetTokenT<reco::VertexCollection> vertex_src_;
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

  edm::Handle<reco::VertexCollection> pvtxs;
  evt.getByToken(vertex_src_, pvtxs);

  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  const MagneticField *fMagneticField = fieldHandle.product();

  edm::ESHandle<TransientTrackBuilder> theB ;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

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

      float ll_vtx_3d_x = ll_prt->userFloat("vtx_3d_x");
      float ll_vtx_3d_y = ll_prt->userFloat("vtx_3d_y");
      float ll_vtx_3d_z = ll_prt->userFloat("vtx_3d_z");

    
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
      cand.addUserFloat("vtx_3d_x",    ll_vtx_3d_x);
      cand.addUserFloat("vtx_3d_y",    ll_vtx_3d_y);
      cand.addUserFloat("vtx_3d_z",    ll_vtx_3d_z);
    
      if( !post_vtx_selection_(cand) ) continue;        



      Float_t iso_sv = 0;
      Int_t ntracks_sv = 0;	  

      //compute isolation
      float l1_iso03 = 0;
      float l1_iso04 = 0;
      float l2_iso03 = 0;
      float l2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon
        if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
        float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4){
          k_iso04 += trk.pt();
          if (dr_to_k < 0.3) k_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }


	if(trk.hasTrackDetails()){
	  
	  if(trk.trackHighPurity() &&
	     trk.pseudoTrack().hitPattern().numberOfValidPixelHits() >= 0 && 
	     trk.pseudoTrack().hitPattern().numberOfValidHits() >= 3 && 
	     trk.pseudoTrack().normalizedChi2() <= 100 && 
	     TMath::Abs(trk.pdgId())==211 &&
	     TMath::Abs(trk.eta()) < 2.5 &&
	     TMath::Abs(trk.charge())==1 && 
	     trk.vertexRef()->z() == ll_vtx_3d_z){

	    //	    Float_t precut_dz = 
	    //	    if(TMath::Abs(precut_dz) > ) continue;
	    
	    reco::TransientTrack tt_track = (*theB).build(trk.pseudoTrack());

	    TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(tt_track.impactPointState(), fitter.fitted_vtx());
	    
	    
	    VertexDistance3D a3d_pf;  
	    
	    std::pair<bool,Measurement1D> cur3DIP_pf = absoluteImpactParameter(tsos_pf, fitter.fitted_refvtx(), a3d_pf);
	    
	    Float_t pvip_pf = cur3DIP_pf.second.value();
	    
	    if(pvip_pf < 0.03){
	      ntracks_sv +=1;
	      iso_sv += trk.pt();
	    }
	  }
	}



      }

      //      std::cout << ntracks_sv << " " << iso_sv << std::endl;

      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );
      cand.addUserFloat("iso_sv" , iso_sv );
      cand.addUserFloat("iso_ntrack" , ntracks_sv );


      reco::Vertex closestVertex;
      bool isVtxFound = false;

      for( reco::VertexCollection::const_iterator vtx = pvtxs->begin(); vtx != pvtxs->end(); ++vtx){	
	if(vtx->position().z()!=ll_vtx_3d_z) continue;
	closestVertex = *vtx;
	isVtxFound = true;
      }

      if(!isVtxFound){
	std::cout << "Vertex was not found !!!!!!!!!!!!!!!" << std::endl;
	continue;
      }

      particle_cand Bcand = calculateIPvariables(extrapolator, fitter.fitted_particle(), fitter.fitted_refvtx(), closestVertex);

      cand.addUserFloat("vtx_3d_lip",   Bcand.lip);
      cand.addUserFloat("vtx_3d_lips",   Bcand.lips);
      cand.addUserFloat("vtx_3d_pvip",   Bcand.pvip);
      cand.addUserFloat("vtx_3d_pvips",   Bcand.pvips);
      cand.addUserFloat("vtx_3d_fl3d",   Bcand.fl3d);
      cand.addUserFloat("vtx_3d_fls3d",   Bcand.fls3d);
      cand.addUserFloat("vtx_3d_alpha",   Bcand.alpha);



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
