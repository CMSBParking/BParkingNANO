// class that creates common vertex using kinematic fit

#ifndef KINFITTER_H
#define KINFITTER_H

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>
#include <cmath>




class KinFitter{

 public: 

  // several constractors for every need
  KinFitter();
  // if refit option set to true will refit daughters before returning them
  // same for mother. If false is the normal p4 sum
  KinFitter( std::vector< RefCountedKinematicParticle >& vallParticles, bool vrefit );
  // fit with mass constraint
  KinFitter( std::vector< RefCountedKinematicParticle >& vallParticles, 
             bool vrefit,
             ParticleMass vKstar_m 
             );

  // destructor    
  virtual ~KinFitter();

  // returns false if fit is unseccesfull
  bool success() { return m_success; }

  // set tracks to fit 
  void SetKinFitTracks( std::vector< RefCountedKinematicParticle >& vallParticles,
                        bool vrefit
                        );

  // momentum of mother
  GlobalVector Mother_Momentum() {
    if ( refit )
       return bs_state.globalMomentum();
    else 
       return UnfittedMotherMomentum();    
  }  
  // mass of mother
  ParticleMass Mother_Mass() {
    if ( refit )
       return bs_state.mass();
    else
       return UnfittedMotherMass();
}
  // charge
  float Mother_Charge() {
       return bs_state.particleCharge();
  }
  
  // energy
  float Mother_Energy() {
    if ( refit )
      return TMath::Sqrt( bs_state.globalMomentum().mag2()                      
                          + pow( bs_state.mass(), 2 )
                          );
    else
      return TMath::Sqrt( UnfittedMotherMomentum().mag2() 
                          + pow( UnfittedMotherMass(), 2 )
                          );
   }

  // chi 
  float chi() { return b_dec_vertex->chiSquared(); }
  // dof
  float dof() { return b_dec_vertex->degreesOfFreedom(); }
  // prob
  float prob() { return ChiSquaredProbability( b_dec_vertex->chiSquared(),
                                               b_dec_vertex->degreesOfFreedom()
                                               ); }
  // momentum of daughters
  GlobalVector Daughter_Momentum( unsigned int vdaughter );
  // mass
  ParticleMass Daughter_Mass( unsigned int vdaughter );
  // charge
  float Daughter_Charge( unsigned int vdaughter );
  // B decay point
  GlobalPoint Mother_XYZ(){ return b_dec_vertex->position(); }
  // error
  GlobalError Mother_XYZError(){ return b_dec_vertex->error(); }
  // errors in kinematic vars
  double Mother_PtError(){ return bs_track.track().ptError(); }
  double Mother_EtaError(){ return bs_track.track().etaError(); }
  double Mother_PhiError(){ return bs_track.track().phiError(); }


private:
  
  bool m_success = false;        bool refit = false;
  RefCountedKinematicTree bsTree;
  RefCountedKinematicVertex b_dec_vertex; 
  KinematicState bs_state;
  RefCountedKinematicParticle b_s;
  std::vector< RefCountedKinematicParticle > bs_children;
  reco::TransientTrack bs_track; 
  unsigned int m_npart;
  GlobalVector UnfittedMotherMomentum();
  ParticleMass UnfittedMotherMass();

};
#endif
