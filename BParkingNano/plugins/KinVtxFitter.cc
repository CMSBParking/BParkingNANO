#include "KinVtxFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" // MIGHT be useful for Phi->KK?
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           std::vector<float> sigmas):
  n_particles_{masses.size()} {
  
  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> particles;
  for(size_t i = 0; i < tracks.size(); ++i) {
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  KinematicParticleVertexFitter kcv_fitter;    
  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(particles);

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != n_particles_) { 
    success_=false; 
    return;
  }
  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
}

// Constrained vertex fit - sequential fit
KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           std::vector<float> sigmas,
                           double mass_constr, double mass_constr_sigma
                           ):
  n_particles_{masses.size()} {
  
  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> leptons;
  std::vector<RefCountedKinematicParticle> particles;

  // Assume the first two paricles are leptons
  for(size_t i = 0; i < 2; ++i) {
    leptons.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  // Assume the remaining particles are kaons/pions
  for(size_t i = 2; i < tracks.size(); ++i) {
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  KinematicParticleVertexFitter kcv_fitter;    

  // Normal fit on two leptons
  RefCountedKinematicTree jp_tree = kcv_fitter.fit(leptons);
  if (jp_tree->isEmpty() || !jp_tree->isValid() || !jp_tree->isConsistent()) {
    success_ = false; 
    return;
  }
  
  jp_tree->movePointerToTheTop(); 
  fitted_particle_ = jp_tree->currentParticle();
  fitted_vtx_ = jp_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    return;
  }

  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = jp_tree->finalStateParticles();
  if(fitted_children_.size() != 2) { 
    success_ = false; 
    return;
  }

  // Constrained fit on two leptons
  KinematicParticleFitter cs_fitter;
  KinematicConstraint * jpsi_c = new MassKinematicConstraint(mass_constr, mass_constr_sigma);
  jp_tree = cs_fitter.fit(jpsi_c, jp_tree);
 
  if (jp_tree->isEmpty() || !jp_tree->isValid() || !jp_tree->isConsistent()) {
    success_ = false; 
    delete jpsi_c;
    return;
  }
  
  jp_tree->movePointerToTheTop(); 
  fitted_particle_ = jp_tree->currentParticle();
  fitted_vtx_ = jp_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    delete jpsi_c;
    return;
  }

  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = jp_tree->finalStateParticles();
  if(fitted_children_.size() != 2) { 
    success_ = false; 
    delete jpsi_c;
    return;
  }


  // Constrained fit on B meson
  jp_tree->movePointerToTheTop();
  RefCountedKinematicParticle jpsi_part = jp_tree->currentParticle();
  particles.push_back(jpsi_part);
  RefCountedKinematicTree vtx_tree;
  try {
    vtx_tree = kcv_fitter.fit(particles);
  } catch (...) {
    std::cout<<"PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
    success_ = false;
    delete jpsi_c;
    return;
  }

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    delete jpsi_c;
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    delete jpsi_c;
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != (n_particles_-1)) { 
    success_ = false; 
    delete jpsi_c;
    return;
  }
  
  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
  delete jpsi_c;
}


// Constrained vertex fit - global fit
KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           std::vector<float> sigmas,
                           double mass_constr
                           ):
  n_particles_{masses.size()} {

  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> particles;
  for(size_t i = 0; i < tracks.size(); ++i) {
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  KinematicConstrainedVertexFitter kcv_fitter;
  MultiTrackKinematicConstraint* jpsi_constr = new TwoTrackMassKinematicConstraint(mass_constr);
  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(particles, jpsi_constr); // Constraint applies to first two particles

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    delete jpsi_constr;
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    delete jpsi_constr;
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != n_particles_) { 
    success_ = false; 
    delete jpsi_constr;
    return;
  }

  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
  delete jpsi_constr;
  return;
}

