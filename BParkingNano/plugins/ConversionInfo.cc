#include "ConversionInfo.h"

////////////////////////////////////////////////////////////////////////////////
// Nancy's baseline selections for conversions
bool ConversionInfo::wpLoose() {
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Nancy's selection for analysis of conversions
bool ConversionInfo::wpTight() {
  return true;
}
    
////////////////////////////////////////////////////////////////////////////////
//
void ConversionInfo::addExtraUserVars() {
}

////////////////////////////////////////////////////////////////////////////////
//
bool ConversionInfo::match(const edm::Handle<reco::BeamSpot>& beamSpot,
			   const edm::Handle<edm::View<reco::Conversion> >& conversions,
			   const pat::Electron& ele) {
  ConversionInfo info;
  return ConversionInfo::match(beamSpot,conversions,ele,info);
}

////////////////////////////////////////////////////////////////////////////////
//
bool ConversionInfo::match(const edm::Handle<reco::BeamSpot>& beamSpot,
			   const edm::Handle<edm::View<reco::Conversion> >& conversions,
			   const pat::Electron& ele,
			   ConversionInfo& info) {
  
  bool matched = false;
  
  // Valid handles?
  if ( !(beamSpot.isValid()) ) {
    edm::LogError("ConversionInfo::match")
      << " !(beamSpot.isValid())" << std::endl;
    return false;
  }
  if ( !(conversions.isValid()) ) {
    edm::LogError("ConversionInfo::match")
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

////////////////////////////////////////////////////////////////////////////////
//
float ConversionInfo::mee(float ipx1, float ipy1, float ipz1, 
			  float ipx2, float ipy2, float ipz2) {
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
