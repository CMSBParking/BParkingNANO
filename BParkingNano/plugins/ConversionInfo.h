#ifndef ConversionInfo_h
#define ConversionInfo_h

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"

class ConversionInfo {
  
 public:
  
  ConversionInfo() {;}
  ~ConversionInfo() {;}
  
  void reset() { ConversionInfo dummy; *this = dummy; }
  
  bool wpOpen();  // Matched to any conversion (without selections)
  bool wpLoose(); // Nancy's baseline selections for conversions
  bool wpTight(); // Nancy's selection for analysis of conversions
  
  void addUserVars(pat::Electron& ele); // adds minimal set of flags to electron userData
  void addUserVarsExtra(pat::Electron& ele); // adds all variables to electron userData
  
  static bool match(const edm::Handle<reco::BeamSpot>& beamSpot,
		    const edm::Handle<edm::View<reco::Conversion> >& conversions,
		    const pat::Electron& ele);
  
  static bool match(const edm::Handle<reco::BeamSpot>& beamSpot,
		    const edm::Handle<edm::View<reco::Conversion> >& conversions,
		    const pat::Electron& ele,
		    ConversionInfo& info);
  
  static float mee(float ipx1, float ipy1, float ipz1, 
		   float ipx2, float ipy2, float ipz2);

public:

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

  // match?
  bool matched = false;
  edm::RefToBase<reco::Track> matched_lead;
  edm::RefToBase<reco::Track> matched_trail;
  
};

#endif // ConversionInfo_h
