#include "KinFitter.h"


KinFitter::KinFitter() {};

KinFitter::KinFitter(std::vector<RefCountedKinematicParticle>& allParticles)
{ 
  KinematicConstrainedVertexFitter kcvFitter;    
  bsTree = kcvFitter.fit(allParticles);
  m_npart=allParticles.size();
  if (bsTree->isEmpty()) {m_success=false; return;}
  if(!bsTree->isValid()) {m_success=false; return;}
  if (!bsTree->isConsistent()) {m_success=false; return;}
  bsTree->movePointerToTheTop(); b_s=bsTree->currentParticle();
  b_dec_vertex=bsTree->currentDecayVertex();
  if (!b_s->currentState().isValid()){ m_success=false; return;}
  bs_state=b_s->currentState();
  if(!b_dec_vertex->vertexIsValid()){m_success=false; return;}
  bs_children = bsTree->finalStateParticles();
  if(bs_children.size()!=m_npart){ m_success=false; return;}
  bs_track=b_s->refittedTransientTrack();
  m_success=true;
}

KinFitter::KinFitter(std::vector<RefCountedKinematicParticle>& allParticles,ParticleMass Kstar_m)
{
    KinematicConstrainedVertexFitter kcvFitter;
    MultiTrackKinematicConstraint * Kstar_c = new TwoTrackMassKinematicConstraint(Kstar_m);
    bsTree = kcvFitter.fit(allParticles,Kstar_c);
    m_npart=allParticles.size();
}

void KinFitter::SetKinFitTracks(std::vector<RefCountedKinematicParticle>& allParticles){
  KinFit=true;
  KinematicConstrainedVertexFitter kcvFitter;    
  bsTree = kcvFitter.fit(allParticles);
  m_npart=allParticles.size();
  if (bsTree->isEmpty()) {m_success=false; return;}
  if (!bsTree->isValid()) {m_success=false; return;}
  if (!bsTree->isConsistent()) {m_success=false; return;}
  bsTree->movePointerToTheTop(); b_s=bsTree->currentParticle();
  b_dec_vertex=bsTree->currentDecayVertex();
  if (!b_s->currentState().isValid()){ m_success=false; return;}
  bs_state=b_s->currentState();
  if (!b_dec_vertex->vertexIsValid()){m_success=false; return;}
  bs_children = bsTree->finalStateParticles();
  if (bs_children.size()!=m_npart){ m_success=false; return;}
  bs_track=b_s->refittedTransientTrack();
  m_success=true;
}


GlobalVector KinFitter::Daughter_Momentum(unsigned int idaughter,bool refit)
{
  KinematicState child_state;
  if(refit) child_state=bs_children.at(idaughter)->currentState();
  else child_state=bs_children.at(idaughter)->initialState();
  return child_state.globalMomentum();
}

ParticleMass KinFitter::Daughter_Mass(unsigned int idaughter,bool refit)
{
  KinematicState child_state;
  if(refit) child_state=bs_children.at(idaughter)->currentState();
  else child_state=bs_children.at(idaughter)->initialState();
  return child_state.mass();
}

float KinFitter::Daughter_Charge(unsigned int idaughter,bool refit)
{
  KinematicState child_state;
  if(refit) child_state=bs_children.at(idaughter)->currentState();
  else child_state=bs_children.at(idaughter)->initialState();
  return child_state.particleCharge();
}

GlobalVector KinFitter::UnfittedMotherMomentum()
{
  TLorentzVector vMother; TLorentzVector v1;
  for (unsigned int i=0; i<m_npart; i++){
     KinematicState child_state=bs_children.at(i)->initialState();
     v1.SetPtEtaPhiM( child_state.globalMomentum().perp(), child_state.globalMomentum().eta(), child_state.globalMomentum().phi(), child_state.mass());
    vMother+=v1;
  }
  return GlobalVector(vMother.Px(),vMother.Py(),vMother.Pz());
}

ParticleMass KinFitter::UnfittedMotherMass()
{
  TLorentzVector vMother; TLorentzVector v1;
  for (unsigned int i=0; i<m_npart; i++){
    KinematicState child_state=bs_children.at(i)->initialState();
    v1.SetPtEtaPhiM( child_state.globalMomentum().perp(), child_state.globalMomentum().eta(), child_state.globalMomentum().phi(), child_state.mass());
    vMother+=v1;
  }
  return vMother.M();
}

KinFitter::~KinFitter() {}
