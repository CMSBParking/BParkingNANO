#ifndef PhysicsTools_BParkingNano_helpers
#define PhysicsTools_BParkingNano_helpers

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <memory>

typedef std::vector<reco::TransientTrack> TransientTrackCollection;

constexpr float K_MASS = 0.493677;
constexpr float PI_MASS = 0.139571;
constexpr float LEP_SIGMA = 0.0000001;
constexpr float K_SIGMA = 0.000016;
constexpr float PI_SIGMA = 0.000016;

inline std::pair<float, float> min_max_dr(const std::vector< edm::Ptr<reco::Candidate> > & cands) {
  float min_dr = std::numeric_limits<float>::max();
  float max_dr = 0.;
  for(size_t i = 0; i < cands.size(); ++i) {
    for(size_t j = i+1; j < cands.size(); ++j) {
      float dr = reco::deltaR(*cands.at(i), *cands.at(j));
      min_dr = std::min(min_dr, dr);
      max_dr = std::max(max_dr, dr);
    }
  }
  return std::make_pair(min_dr, max_dr);
}

template<typename FITTER>
inline double cos_theta_2D(const FITTER& fitter, const reco::BeamSpot &bs, const reco::Candidate::LorentzVector& p4) {
  if(!fitter.success()) return -2;
  GlobalPoint point = fitter.fitted_vtx();
  auto bs_pos = bs.position(point.z());
  math::XYZVector delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);
  math::XYZVector pt(p4.px(), p4.py(), 0.);
  double den = (delta.R() * pt.R());
  return (den != 0.) ? delta.Dot(pt)/den : -2;
}

template<typename FITTER>
inline Measurement1D l_xy(const FITTER& fitter, const reco::BeamSpot &bs) {
  if(!fitter.success()) return {-2, -2};
  GlobalPoint point = fitter.fitted_vtx();
  GlobalError err = fitter.fitted_vtx_uncertainty();
  auto bs_pos = bs.position(point.z());
  GlobalPoint delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);  
  return {delta.perp(), err.rerr(delta)};
}

#endif
>>>>>>> expanded common helpers
