#include "KinFitter.h"


KinFitter::KinFitter() {};

KinFitter::KinFitter( std::vector< RefCountedKinematicParticle >& allParticles,
                      bool vrefit )
{ 

  KinematicConstrainedVertexFitter kcvFitter;    
  bsTree = kcvFitter.fit( allParticles );
  m_npart = allParticles.size();

  // check if fit is succesfull
  if( bsTree->isEmpty() )       { m_success=false; return; }
  if( !bsTree->isValid() )      { m_success=false; return; }
  if( !bsTree->isConsistent() ) { m_success=false; return; }
  // get B
  bsTree->movePointerToTheTop(); 
  b_s = bsTree->currentParticle();
  b_dec_vertex = bsTree->currentDecayVertex();
  // more quality checks
  if( !b_s->currentState().isValid() ) { m_success=false; return; }
  bs_state = b_s->currentState();
  if( !b_dec_vertex->vertexIsValid() ) { m_success=false; return; }
  bs_children = bsTree->finalStateParticles();
  if( bs_children.size() != m_npart )  { m_success=false; return;}
  bs_track = b_s->refittedTransientTrack();
  // if pass set success to true
  m_success = true;
  refit = vrefit;

}

// fit with mass constraint
KinFitter::KinFitter( std::vector< RefCountedKinematicParticle >& allParticles,
                      bool vrefit,
                      ParticleMass Kstar_m
                     )
{
    KinematicConstrainedVertexFitter kcvFitter;
    MultiTrackKinematicConstraint* Kstar_c = new TwoTrackMassKinematicConstraint( Kstar_m );
    bsTree = kcvFitter.fit( allParticles, Kstar_c );
    m_npart = allParticles.size();
    // check if fit is succesfull
    if( bsTree->isEmpty() )       { m_success=false; return; }
    if( !bsTree->isValid() )      { m_success=false; return; }
    if( !bsTree->isConsistent() ) { m_success=false; return; }
    // get B
    bsTree->movePointerToTheTop(); 
    b_s = bsTree->currentParticle();
    b_dec_vertex = bsTree->currentDecayVertex();
    // more quality checks
    if( !b_s->currentState().isValid() ) { m_success=false; return; }
    bs_state = b_s->currentState();
    if( !b_dec_vertex->vertexIsValid() ) { m_success=false; return; }
    bs_children = bsTree->finalStateParticles();
    if( bs_children.size() != m_npart )  { m_success=false; return;}
    bs_track = b_s->refittedTransientTrack();
    // if pass set success to true
    m_success = true;
    refit = vrefit;    
}


void KinFitter::SetKinFitTracks( std::vector< RefCountedKinematicParticle >& allParticles, bool vrefit ){
  KinematicConstrainedVertexFitter kcvFitter;    
  bsTree = kcvFitter.fit( allParticles );
  m_npart = allParticles.size();
  // same checks as above
  if ( bsTree->isEmpty() )       { m_success = false; return; }
  if ( !bsTree->isValid() )      { m_success = false; return; }
  if ( !bsTree->isConsistent() ) { m_success = false; return; }
  bsTree->movePointerToTheTop(); 
  b_s = bsTree->currentParticle();
  b_dec_vertex = bsTree->currentDecayVertex();
  if ( !b_s->currentState().isValid() ){ m_success = false; return; }
  bs_state = b_s->currentState();
  if ( !b_dec_vertex->vertexIsValid() ){ m_success = false; return; }
  bs_children = bsTree->finalStateParticles();
  if ( bs_children.size() != m_npart ){ m_success = false; return; }
  bs_track = b_s->refittedTransientTrack();
  m_success = true;
  refit = vrefit;
}


GlobalVector KinFitter::Daughter_Momentum( unsigned int vdaughter )
{
  KinematicState child_state;
  if( refit ) child_state = bs_children.at( vdaughter )->currentState();
  else child_state = bs_children.at( vdaughter )->initialState();
  return child_state.globalMomentum();
}

ParticleMass KinFitter::Daughter_Mass( unsigned int vdaughter )
{
  KinematicState child_state;
  if( refit ) child_state = bs_children.at( vdaughter )->currentState();
  else child_state = bs_children.at( vdaughter )->initialState();
  return child_state.mass();
}

float KinFitter::Daughter_Charge( unsigned int vdaughter )
{
  KinematicState child_state;
  if( refit ) child_state = bs_children.at( vdaughter )->currentState();
  else child_state = bs_children.at( vdaughter )->initialState();
  return child_state.particleCharge();
}

GlobalVector KinFitter::UnfittedMotherMomentum()
{
  TLorentzVector vMother, v1;
  // loop over initial parts and compute p4 in a loop
  for ( unsigned int i=0; i < m_npart; ++i ){
     KinematicState child_state = bs_children.at( i )->initialState();
     v1.SetPtEtaPhiM( child_state.globalMomentum().perp(), 
                      child_state.globalMomentum().eta(), 
                      child_state.globalMomentum().phi(), 
                      child_state.mass()
                    );
    vMother += v1;
  }
  return GlobalVector( vMother.Px(), vMother.Py(), vMother.Pz() );
}

ParticleMass KinFitter::UnfittedMotherMass()
{
  TLorentzVector vMother, v1;
  for ( unsigned int i=0; i < m_npart; i++ ){
    KinematicState child_state = bs_children.at( i )->initialState();
    v1.SetPtEtaPhiM( child_state.globalMomentum().perp(), 
                     child_state.globalMomentum().eta(), 
                     child_state.globalMomentum().phi(), 
                     child_state.mass()
                    );
    vMother += v1;
  }
  return vMother.M();
}

KinFitter::~KinFitter() {}
