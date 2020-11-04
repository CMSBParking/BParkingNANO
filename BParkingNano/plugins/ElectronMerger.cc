// Merges the PF and LowPT collections, sets the isPF and isLowPt 
// UserInt's accordingly

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <limits>
#include <algorithm>
#include "helper.h"

class ConversionInfo {
public:

  ConversionInfo() {;}
  ~ConversionInfo() {;}
  void reset() { ConversionInfo dummy; *this = dummy; }

  // quality
  bool valid = false;
  float chi2prob = -1.;
  bool quality_high_purity = false;
  bool quality_high_efficiency = false;
  // tracks
  uint ntracks = 0;
  float min_trk_pt = -1.;
  int ilead = -1;
  int itrail = -1;
  // displacement
  float l_xy = -1.;
  float vtx_radius = -1.;
  // invariant mass
  float mass_from_conv = -1.;
  float mass_from_Pin = -1.;
  float mass_before_fit = -1.;
  float mass_after_fit = -1.;
  // hits before vertex
  uint lead_nhits_before_vtx = 0;
  uint trail_nhits_before_vtx = 0;
  uint max_nhits_before_vtx = 0;
  uint sum_nhits_before_vtx = 0;
  int delta_expected_nhits_inner = 0;
  // opening angle
  float delta_cot_from_Pin = -1.;
  // matched tracks
  edm::RefToBase<reco::Track> matched_lead;
  edm::RefToBase<reco::Track> matched_trail;
  
};

class ElectronMerger : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)


public:
  bool debug=false; 

  explicit ElectronMerger(const edm::ParameterSet &cfg):
    triggerMuons_{ consumes<pat::MuonCollection>( cfg.getParameter<edm::InputTag>("trgMuon") )},
    lowpt_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("lowptSrc") )},
    pf_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("pfSrc") )},
    ptBiased_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("ptbiasedSeeding") )},
    unBiased_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("unbiasedSeeding") )},
    mvaId_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("mvaId") )},
    pf_mvaId_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("pfmvaId") )},
    vertexSrc_{ consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>("vertexCollection") )},
    conversions_{ consumes<edm::View<reco::Conversion> > ( cfg.getParameter<edm::InputTag>("conversions") )},
    beamSpot_{ consumes<reco::BeamSpot> ( cfg.getParameter<edm::InputTag>("beamSpot") )},
    drTrg_cleaning_{cfg.getParameter<double>("drForCleaning_wrtTrgMuon")},
    dzTrg_cleaning_{cfg.getParameter<double>("dzForCleaning_wrtTrgMuon")},
    dr_cleaning_{cfg.getParameter<double>("drForCleaning")},
    dz_cleaning_{cfg.getParameter<double>("dzForCleaning")},
    flagAndclean_{cfg.getParameter<bool>("flagAndclean")},
    pf_ptMin_{cfg.getParameter<double>("pf_ptMin")},
    ptMin_{cfg.getParameter<double>("ptMin")},
    etaMax_{cfg.getParameter<double>("etaMax")},
    bdtMin_{cfg.getParameter<double>("bdtMin")},
    use_gsf_mode_for_p4_{cfg.getParameter<bool>("useGsfModeForP4")},
    use_regression_for_p4_{cfg.getParameter<bool>("useRegressionModeForP4")},
    sortOutputCollections_{cfg.getParameter<bool>("sortOutputCollections")},
    saveLowPtE_{cfg.getParameter<bool>("saveLowPtE")}
    {
       produces<pat::ElectronCollection>("SelectedElectrons");
       produces<TransientTrackCollection>("SelectedTransientElectrons");  
    }

  ~ElectronMerger() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

  bool matchToConversionTracks(const edm::Handle<reco::BeamSpot>& beamSpot,
			       const edm::Handle<edm::View<reco::Conversion> >& conversions,
			       const pat::Electron& ele) const;
  
  bool matchToConversionTracks(const edm::Handle<reco::BeamSpot>& beamSpot,
			       const edm::Handle<edm::View<reco::Conversion> >& conversions,
			       const pat::Electron& ele,
			       ConversionInfo& info) const;
  
  float mee(float ipx1, float ipy1, float ipz1, 
	    float ipx2, float ipy2, float ipz2) const;
  
private:
  const edm::EDGetTokenT<pat::MuonCollection> triggerMuons_;
  const edm::EDGetTokenT<pat::ElectronCollection> lowpt_src_;
  const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ptBiased_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> unBiased_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mvaId_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> pf_mvaId_src_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  const edm::EDGetTokenT<edm::View<reco::Conversion> > conversions_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  const double drTrg_cleaning_;
  const double dzTrg_cleaning_;
  const double dr_cleaning_;
  const double dz_cleaning_;
  const bool flagAndclean_;
  const double pf_ptMin_;
  const double ptMin_; //pt min cut
  const double etaMax_; //eta max cut
  const double bdtMin_; //bdt min cut
  const bool use_gsf_mode_for_p4_;
  const bool use_regression_for_p4_;
  const bool sortOutputCollections_;
  const bool saveLowPtE_;

};

void ElectronMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & iSetup) const {

  //input
  edm::Handle<pat::MuonCollection> trgMuon;
  evt.getByToken(triggerMuons_, trgMuon);
  edm::Handle<pat::ElectronCollection> lowpt;
  evt.getByToken(lowpt_src_, lowpt);
  edm::Handle<pat::ElectronCollection> pf;
  evt.getByToken(pf_src_, pf);
  edm::Handle<edm::ValueMap<float> > ptBiased;
  evt.getByToken(ptBiased_src_, ptBiased);
  edm::Handle<edm::ValueMap<float> > unBiased;
  evt.getByToken(unBiased_src_, unBiased);
  edm::Handle<edm::ValueMap<float> > mvaId;  
  evt.getByToken(mvaId_src_, mvaId);
  edm::Handle<edm::ValueMap<float> > pfmvaId;  
  evt.getByToken(pf_mvaId_src_, pfmvaId);
  // 
  edm::ESHandle<TransientTrackBuilder> theB ;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  //
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();
  //
  edm::Handle<edm::View<reco::Conversion> > conversions;
  evt.getByToken(conversions_, conversions);
  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByToken(beamSpot_, beamSpot);

  // output
  std::unique_ptr<pat::ElectronCollection>  ele_out      (new pat::ElectronCollection );
  std::unique_ptr<TransientTrackCollection> trans_ele_out(new TransientTrackCollection);
  std::vector<std::pair<float, float>> pfEtaPhi;
  std::vector<float> pfVz;
  
  // -> changing order of loops ert Arabella's fix this without need for more vectors  
  size_t ipfele=-1;
  for(auto ele : *pf) {
   ipfele++;

   if (debug) std::cout << "ElectronMerger, Event " << (evt.id()).event() 
			<< " => PF: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
			<< ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
			<< ", ele gsf track chi2 = " << ele.gsfTrack()->normalizedChi2()
			<< ", ele.p = " << ele.p() << std::endl;

   //cuts
   if (ele.pt()<ptMin_ || ele.pt() < pf_ptMin_) continue;
   if (fabs(ele.eta())>etaMax_) continue;
   // apply conversion veto unless we want conversions
   if (!ele.passConversionVeto()) continue;

   // Fix the mass to the proper one
   reco::Candidate::PolarLorentzVector p4( 
					  ele.pt(),
					  ele.eta(),
					  ele.phi(),
					  ELECTRON_MASS
					   );
   ele.setP4(p4);     

   // skip electrons inside tag's jet or from different PV
   bool skipEle=true;
   for(const auto & trg : *trgMuon) {
     if(reco::deltaR(ele, trg) < drTrg_cleaning_ && drTrg_cleaning_ > 0)
        continue;
     if(fabs(ele.vz() - trg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0)
        continue;
     skipEle=false;
     break; // one trg muon to pass is enough :)
   }
   // we skip evts without trg muon
   if (skipEle) continue;

   // for PF e we set BDT outputs to much higher number than the max
   edm::Ref<pat::ElectronCollection> ref(pf,ipfele);
   float pf_mva_id = float((*pfmvaId)[ref]);
   ele.addUserInt("isPF", 1);
   ele.addUserInt("isLowPt", 0);
   ele.addUserFloat("ptBiased", 20.);
   ele.addUserFloat("unBiased", 20.);
   ele.addUserFloat("mvaId", 20.);
   ele.addUserFloat("pfmvaId", pf_mva_id);
   ele.addUserFloat("chargeMode", ele.charge());
   ele.addUserInt("isPFoverlap", 0);

   pfEtaPhi.push_back(std::pair<float, float>(ele.eta(), ele.phi()));
   pfVz.push_back(ele.vz());
   ele_out       -> emplace_back(ele);
  }

  unsigned int pfSelectedSize = pfEtaPhi.size();

  if ( saveLowPtE_ ) {
  size_t iele=-1;
  /// add and clean low pT e
  for(auto ele : *lowpt) {
    iele++;

    if (debug) std::cout << "ElectronMerger, Event " << (evt.id()).event() 
			 << " => LPT: ele.superCluster()->rawEnergy() = " << ele.superCluster()->rawEnergy()
			 << ", ele.correctedEcalEnergy() = " << ele.correctedEcalEnergy()
			 << ", ele gsf track chi2 = " << ele.gsfTrack()->normalizedChi2()
			 << ", ele.p = " << ele.p() << std::endl;
   
    //take modes
   if (use_regression_for_p4_) {
     // pt from regression, eta and phi from gsf track mode
     reco::Candidate::PolarLorentzVector p4( ele.pt(),
                                             ele.gsfTrack()->etaMode(),
                                             ele.gsfTrack()->phiMode(),
                                             ELECTRON_MASS    );
     ele.setP4(p4);
   }else if(use_gsf_mode_for_p4_) {
     reco::Candidate::PolarLorentzVector p4( ele.gsfTrack()->ptMode(),
                                             ele.gsfTrack()->etaMode(),
                                             ele.gsfTrack()->phiMode(),
                                             ELECTRON_MASS    );
     ele.setP4(p4);
   } else {
     // Fix the mass to the proper one
     reco::Candidate::PolarLorentzVector p4( 
       ele.pt(),
       ele.eta(),
       ele.phi(),
       ELECTRON_MASS
       );
     ele.setP4(p4);     
   }
  
   //same cuts as in PF
   if (ele.pt()<ptMin_) continue;
   if (fabs(ele.eta())>etaMax_) continue;
   // apply conversion veto?
   if (!ele.passConversionVeto()) continue;

   // Match to tracks in "gsfTracksOpenConversions" collection
   ConversionInfo info;
   bool matched = matchToConversionTracks(beamSpot,conversions,ele,info);
   ele.addUserInt("isConversion", matched?1:0);
   ele.addUserInt("convLeadTrack", info.matched_lead.isNonnull()?1:0);
   ele.addUserInt("convTrailTrack", info.matched_trail.isNonnull()?1:0);

   if (/*debug && */matched) {
     std::cout << "[ElectronMerger::produce]"
	       << " event " << (evt.id()).event()
	       << ", electron: " << iele
	       << ", isConversion: " << (matched?1:0)
	       << ", convLeadTrack: " << (info.matched_lead.isNonnull()?1:0)
	       << ", convTrailTrack: " << (info.matched_trail.isNonnull()?1:0)
	       << std::endl;
   }

   //assigning BDT values
   edm::Ref<pat::ElectronCollection> ref(lowpt,iele);
   float mva_id = float((*mvaId)[ref]);
 //  if ( unbiased_seedBDT <bdtMin_) continue; //extra cut for low pT e on BDT
   if ( mva_id <bdtMin_) continue; //extra cut for low pT e on BDT


   bool skipEle=true;
   for(const auto & trg : *trgMuon) {
     if(reco::deltaR(ele, trg) < drTrg_cleaning_ && drTrg_cleaning_ > 0)
        continue;
     if(fabs(ele.vz() - trg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0)
        continue;
     skipEle=false;
     break;  // one trg muon is enough 
   }
   // same here Do we need evts without trg muon? now we skip them
   if (skipEle) continue;   
   

   //pf cleaning    
   bool clean_out = false;
   for(unsigned int iEle=0; iEle<pfSelectedSize; ++iEle) {

      clean_out |= (
	           fabs(pfVz[iEle] - ele.vz()) < dz_cleaning_ &&
                   reco::deltaR(ele.eta(), ele.phi(), pfEtaPhi[iEle].first, pfEtaPhi[iEle].second) < dr_cleaning_   );

   }
   if(clean_out && flagAndclean_) continue;
   else if(clean_out) ele.addUserInt("isPFoverlap", 1);
   else ele.addUserInt("isPFoverlap", 0);

   const reco::GsfTrackRef gsfTrk = ele.gsfTrack();
   float unbiased_seedBDT = float((*unBiased)[gsfTrk]);
   float ptbiased_seedBDT = float((*ptBiased)[gsfTrk]);
   ele.addUserInt("isPF", 0);
   ele.addUserInt("isLowPt", 1);
   ele.addUserFloat("chargeMode", ele.gsfTrack()->chargeMode());
   ele.addUserFloat("ptBiased", ptbiased_seedBDT);
   ele.addUserFloat("unBiased", unbiased_seedBDT);
   ele.addUserFloat("mvaId", mva_id);
   ele.addUserFloat("pfmvaId", 20.);

   ele_out       -> emplace_back(ele);
  }
}
  if(sortOutputCollections_){

    //sorting increases sligtly the time but improves the code efficiency in the Bcandidate builder
    //easier identification of leading and subleading with smarter loop
    std::sort( ele_out->begin(), ele_out->end(), [] (pat::Electron e1, pat::Electron e2) -> bool {return e1.pt() > e2.pt();}
             );
  }

  // build transient track collection
  for(auto &ele : *ele_out){
    const reco::TransientTrack eleTT =(*theB).buildfromGSF( ele.gsfTrack() );
    trans_ele_out -> emplace_back(eleTT);

    if(ele.userInt("isPF")) continue;
    //compute IP for electrons: need transient track
    //from PhysicsTools/PatAlgos/plugins/LeptonUpdater.cc
    const reco::GsfTrackRef gsfTrk = ele.gsfTrack();
    // PVDZ
    ele.setDB(gsfTrk->dz(PV.position()), std::hypot(gsfTrk->dzError(), PV.zError()), pat::Electron::PVDZ);

    //PV2D
    std::pair<bool, Measurement1D> result = IPTools::signedTransverseImpactParameter(eleTT, GlobalVector(gsfTrk->px(), gsfTrk->py(), gsfTrk->pz()), PV);
    double d0_corr = result.second.value();
    double d0_err = PV.isValid() ? result.second.error() : -1.0;
    ele.setDB(d0_corr, d0_err, pat::Electron::PV2D);

    // PV3D
    result = IPTools::signedImpactParameter3D(eleTT, GlobalVector(gsfTrk->px(), gsfTrk->py(), gsfTrk->pz()), PV);
    d0_corr = result.second.value();
    d0_err = PV.isValid() ? result.second.error() : -1.0;
    ele.setDB(d0_corr, d0_err, pat::Electron::PV3D);
  }
   
  //adding label to be consistent with the muon and track naming
  evt.put(std::move(ele_out),      "SelectedElectrons");
  evt.put(std::move(trans_ele_out),"SelectedTransientElectrons");
}

bool ElectronMerger::matchToConversionTracks(const edm::Handle<reco::BeamSpot>& beamSpot,
					     const edm::Handle<edm::View<reco::Conversion> >& conversions,
					     const pat::Electron& ele) const {
  ConversionInfo info;
  return matchToConversionTracks(beamSpot,conversions,ele,info);
}

bool ElectronMerger::matchToConversionTracks(const edm::Handle<reco::BeamSpot>& beamSpot,
					     const edm::Handle<edm::View<reco::Conversion> >& conversions,
					     const pat::Electron& ele,
					     ConversionInfo& info) const {
  
  bool matched = false;

  // Valid handles?
  if ( !(beamSpot.isValid()) ) {
    edm::LogError("ElectronMerger::matchToConversionTracks")
      << " !(beamSpot.isValid())" << std::endl;
    return false;
  }
  if ( !(conversions.isValid()) ) {
    edm::LogError("ElectronMerger::matchToConversionTracks")
      << " !(conversions.isValid())" << std::endl;
    return false;
  }
  
  // Iterate through conversions and calculate quantities (requirement from Nancy)
  for ( const auto& conv : *conversions ) {

    // Filter
    if ( conv.tracks().size() != 2 ) { continue; }

    // Quality
    info.valid = conv.conversionVertex().isValid(); // (=true)
    info.chi2prob = ChiSquaredProbability(conv.conversionVertex().chi2(),conv.conversionVertex().ndof()); // (<0.005)
    info.quality_high_purity = conv.quality(reco::Conversion::highPurity); // (=true)
    info.quality_high_efficiency = conv.quality(reco::Conversion::highEfficiency); // (none)

    // Tracks
    info.ntracks = conv.tracks().size(); // (=2)
    info.min_trk_pt = -1.; // (>0.5)
    for ( const auto& trk : conv.tracks() ) {
      if ( info.min_trk_pt < 0. || trk->pt() < info.min_trk_pt ) { info.min_trk_pt = trk->pt(); }
    }
    info.ilead = -1; info.itrail = -1;
    if ( conv.tracks().size() == 2 ) {
      edm::RefToBase<reco::Track> trk1 = conv.tracks().front();
      edm::RefToBase<reco::Track> trk2 = conv.tracks().back();
      if      ( trk1->pt() > trk2->pt() ) { info.ilead = 0; info.itrail = 1; }
      else if ( trk1->pt() < trk2->pt() ) { info.ilead = 1; info.itrail = 0; }
    }

    // Transverse displacement (with respect to beamspot) and vertex radius
    math::XYZVectorF p_refitted =  conv.refittedPairMomentum();
    float dx = conv.conversionVertex().x() - beamSpot->x0();
    float dy = conv.conversionVertex().y() - beamSpot->y0();
    info.l_xy = (p_refitted.x()*dx + p_refitted.y()*dy) / p_refitted.rho();
    info.vtx_radius = sqrt(conv.conversionVertex().position().perp2()); // (1.5<r<4.)

    // invariant mass from track pair from conversion
    info.mass_from_conv = conv.pairInvariantMass();
    
    // Invariant mass from Pin before fit to common vertex 
    if ( conv.tracksPin().size() >= 2 ) {
      math::XYZVectorF lead_Pin = conv.tracksPin().at(info.ilead);
      math::XYZVectorF trail_Pin = conv.tracksPin().at(info.itrail);
      info.mass_from_Pin = mee( lead_Pin.x(), lead_Pin.y(), lead_Pin.z(),
				trail_Pin.x(), trail_Pin.y(), trail_Pin.z() );
      // Opening angle
      info.delta_cot_from_Pin = 1. / tan(trail_Pin.theta()) - 1. / tan(lead_Pin.theta());
    }

    // Invariant mass before fit to common vertex
    if ( conv.tracks().size() >= 2 ) {
      edm::RefToBase<reco::Track> lead_before_vtx_fit = conv.tracks().at(info.ilead);
      edm::RefToBase<reco::Track> trail_before_vtx_fit = conv.tracks().at(info.itrail);
      info.mass_before_fit = mee( lead_before_vtx_fit->px(), lead_before_vtx_fit->py(), lead_before_vtx_fit->pz(),
				  trail_before_vtx_fit->px(), trail_before_vtx_fit->py(), trail_before_vtx_fit->pz() );
    }

    // Invariant mass after the fit to common vertex
    if ( conv.conversionVertex().refittedTracks().size() >=2 ) {
      const reco::Track lead_after_vtx_fit = conv.conversionVertex().refittedTracks().at(info.ilead);
      const reco::Track trail_after_vtx_fit = conv.conversionVertex().refittedTracks().at(info.itrail);
      info.mass_after_fit = mee( lead_after_vtx_fit.px(), lead_after_vtx_fit.py(), lead_after_vtx_fit.pz(),
				 trail_after_vtx_fit.px(), trail_after_vtx_fit.py(), trail_after_vtx_fit.pz());
      // Difference in expeted hits
      info.delta_expected_nhits_inner =
	lead_after_vtx_fit.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)
	- trail_after_vtx_fit.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    }
    
    // Hits prior to vertex
    info.lead_nhits_before_vtx  = conv.nHitsBeforeVtx().size() > 1 ? conv.nHitsBeforeVtx().at(info.ilead) : 0;
    info.trail_nhits_before_vtx = conv.nHitsBeforeVtx().size() > 1 ? conv.nHitsBeforeVtx().at(info.itrail) : 0;
    info.max_nhits_before_vtx = conv.nHitsBeforeVtx().size() > 1 ?
      ( conv.nHitsBeforeVtx().at(0) > conv.nHitsBeforeVtx().at(1) ?
	conv.nHitsBeforeVtx().at(0) :
	conv.nHitsBeforeVtx().at(1) ) : 0;
    info.sum_nhits_before_vtx = conv.nHitsBeforeVtx().size() > 1 ?
      conv.nHitsBeforeVtx().at(0) +
      conv.nHitsBeforeVtx().at(1) : 0;
    
    // Attempt to match conversion track to electron
    for ( uint itrk = 0; itrk < conv.tracks().size(); ++itrk ) {

      edm::RefToBase<reco::Track> trk = conv.tracks()[itrk];
      if ( trk.isNull() ) { continue; }
      
      reco::GsfTrackRef gsf = ele.gsfTrack();
      if ( gsf.isNull() ) { continue; }

      if ( gsf.id() == trk.id() && gsf.key() == trk.key() )  { 
	if ( (int)itrk == info.ilead ) {
	  info.matched_lead = trk;
	  matched = true;
	}
	if ( (int)itrk == info.itrail ) {
	  info.matched_trail = trk;
	  matched = true;
	}
      }
      
    } // track loop
    
  } // conversions loop

  return matched;

}

float ElectronMerger::mee(float ipx1, float ipy1, float ipz1, 
			  float ipx2, float ipy2, float ipz2) const {
  const float  m = 0.000511;
  const float px = ipx1+ipx2;
  const float py = ipy1+ipy2;
  const float pz = ipz1+ipz2;
  const float p1 = ipx1*ipx1 + ipy1*ipy1 + ipz1*ipz1;
  const float p2 = ipx2*ipx2 + ipy2*ipy2 + ipz2*ipz2;
  const float  e = sqrt( p1 + m*m ) + sqrt( p2 + m*m );
  const float mass = ( e*e - px*px - py*py - pz*pz );
  return mass > 0. ? sqrt(mass) : -1.;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronMerger);
