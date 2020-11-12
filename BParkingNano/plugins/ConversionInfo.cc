#include "ConversionInfo.h"

////////////////////////////////////////////////////////////////////////////////
// Matched to any conversion (without selections)
// 
bool ConversionInfo::wpOpen() {
  if ( matched ) { return true; }
  else { return false; }
}

////////////////////////////////////////////////////////////////////////////////
// Nancy's baseline selections for conversions
// Based on: https://github.com/CMSBParking/BParkingNANO/blob/b2664ed/BParkingNano/plugins/ConversionSelector.cc#L253-L300
bool ConversionInfo::wpLoose() {
  if ( wpOpen() &&
       ntracks == 2 &&
       valid == true &&
       quality_high_purity == true &&
       min_trk_pt > 0.5 &&
       chi2prob > 0.0005 ) { return true; }
  else { return false; }
}

////////////////////////////////////////////////////////////////////////////////
// Nancy's selection for analysis of conversions
// Based on: slide 20 of https://indico.cern.ch/event/814402/contributions/3401312/
bool ConversionInfo::wpTight() {
  if ( wpLoose() &&
       sum_nhits_before_vtx <= 1 &&
       l_xy > 0. &&
       // min_trk_pt > 1. &&
       // vtx_radius < 4. &&
       mass_from_conv > 0. && // sanity check
       mass_from_conv < 0.05 ) { return true; }
  else { return false; }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// adds minimal set of flags to electron userData
void ConversionInfo::addUserVars(pat::Electron& ele) {
  ele.addUserInt("convOpen", this->matched?1:0);
  ele.addUserInt("convLoose", this->wpLoose()?1:0);
  ele.addUserInt("convTight", this->wpTight()?1:0);
  ele.addUserInt("convLead", this->matched_lead.isNonnull()?1:0);
  ele.addUserInt("convTrail", this->matched_trail.isNonnull()?1:0);
  //ele.addUserInt("convExtra", 0);
  if ( ele.hasUserInt("convExtra") == false ) { ele.addUserInt("convExtra", 0); }
}

////////////////////////////////////////////////////////////////////////////////
// adds all variables to electron userData
void ConversionInfo::addUserVarsExtra(pat::Electron& ele) {
  
  // Flag that indicates if extra variables are added to electron userData
  ele.addUserInt("convExtra", 1, true); // overwrite
  
  // quality
  ele.addUserInt("convValid", this->valid?1:0);
  ele.addUserFloat("convChi2Prob", this->chi2prob);
  ele.addUserInt("convQualityHighPurity", this->quality_high_purity?1:0);
  ele.addUserInt("convQualityHighEff", this->quality_high_efficiency?1:0);
  
  // tracks
  ele.addUserInt("convTracksN", this->ntracks);
  ele.addUserFloat("convMinTrkPt", this->min_trk_pt);
  ele.addUserInt("convLeadIdx", this->ilead);
  ele.addUserInt("convTrailIdx", this->itrail);

  // displacement
  ele.addUserFloat("convLxy", this->l_xy);
  ele.addUserFloat("convVtxRadius", this->vtx_radius);

  // invariant mass
  ele.addUserFloat("convMass", this->mass_from_conv);
  ele.addUserFloat("convMassFromPin", this->mass_from_Pin);
  ele.addUserFloat("convMassBeforeFit", this->mass_before_fit);
  ele.addUserFloat("convMassAfterFit", this->mass_after_fit);

  // hits before vertex
  ele.addUserInt("convLeadNHitsBeforeVtx", this->lead_nhits_before_vtx);
  ele.addUserInt("convTrailNHitsBeforeVtx", this->trail_nhits_before_vtx);
  ele.addUserInt("convMaxNHitsBeforeVtx", this->max_nhits_before_vtx);
  ele.addUserInt("convSumNHitsBeforeVtx", this->sum_nhits_before_vtx);
  ele.addUserInt("convDeltaExpectedNHitsInner", this->delta_expected_nhits_inner);

  // opening angle
  ele.addUserFloat("convDeltaCotFromPin", this->delta_cot_from_Pin);

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
	  info.matched = true;
	  info.matched_lead = trk;
	}
	if ( (int)itrk == info.itrail ) {
	  info.matched = true;
	  info.matched_trail = trk;
	}
      }
      
    } // track loop
    
  } // conversions loop

  return info.matched;

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
