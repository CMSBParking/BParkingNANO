#ifndef HELPER_H
#define HELPER_H

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"
#include "Math/LorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


constexpr float K_MASS = 0.493677;
constexpr float PI_MASS = 0.139571;
constexpr float LEP_SIGMA = 0.0000001;
constexpr float K_SIGMA = 0.000016;
constexpr float PI_SIGMA = 0.000016;

typedef std::vector<reco::TransientTrack> TransientTrackCollection;
constexpr float MOUN_MASS = 0.10565837;




inline GlobalPoint FlightDistVector (const reco::BeamSpot & bm, GlobalPoint Bvtx)
{
   GlobalPoint Dispbeamspot(-1*( (bm.x0()-Bvtx.x()) + (Bvtx.z()-bm.z0()) * bm.dxdz()),
			   -1*( (bm.y0()-Bvtx.y()) + (Bvtx.z()-bm.z0()) * bm.dydz()), 
                            0);                    
   return std::move(Dispbeamspot);
}


inline float CosA(GlobalPoint & dist, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> & Bp4)
{
    math::XYZVector vperp(dist.x(),dist.y(),0);
    math::XYZVector pperp(Bp4.Px(),Bp4.Py(),0); 
    return std::move(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
}



 
#endif
