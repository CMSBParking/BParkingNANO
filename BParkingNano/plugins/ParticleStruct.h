#ifndef ParticleStruct_H
#define ParticleStruct_H

struct particle_cand 
{
  
  // Longitudinal impact parameter and its significance
  Float_t lip;
  Float_t lips;
  
  // Impact parameter for the PV and its significance
  Float_t pvip;
  Float_t pvips;
  
  // Flight length and its significance
  Float_t fl3d;
  Float_t fls3d;
  
  // opening angle in 3D
  Float_t alpha;
  
};


#endif
