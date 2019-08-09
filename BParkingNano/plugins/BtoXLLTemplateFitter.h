#ifndef BTOXLLTEMPLATEFITTER_H
#define BTOXLLTEMPLATEFITTER_H
#include "Math/LorentzVector.h"
#include "Math/UnaryOperators.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "KinFitter.h"
#include "PreFitter.h"
#include "Math/PxPyPzM4D.h"
#include "helper.h"

template<typename T>
class BtoXLLTemplateFitter{
 private:
  const std::vector<T> & leptons;
  const std::vector<pat::CompositeCandidate> & tracks;
  const reco::BeamSpot & Bspot;
  edm::ESHandle<MagneticField> bField;
  std::pair< std::vector<T>, std::vector<T> > pairs;
  
  double ptLep1Cut_=0;  double mllMin_=0;  double mllMax_=0; 
  double  MBMin_ =0;    double  MBMax_ =0; double BPtCut_=0;
  double drTrkLepCut_ =0; float chi=0, ndf=0, part_sigma=0, kaon_sigma=0;
 public:
  BtoXLLTemplateFitter( const std::vector<T> & vleptons,
                        const std::vector<pat::CompositeCandidate> & vtracks, 
                        const reco::BeamSpot & vBspot,
                        const edm::ParameterSet & runParameters,
                        edm::ESHandle<MagneticField> & bFieldHandle );
  //  t1(t),{}
  ~BtoXLLTemplateFitter(){}
  std::unique_ptr<std::vector<pat::CompositeCandidate>> ReconstructB(float trackMass);
  


};
#endif

template<typename T> 
BtoXLLTemplateFitter<T>::BtoXLLTemplateFitter(
                  const std::vector<T> & vleptons,
                  const std::vector<pat::CompositeCandidate> & vtracks, 
                  const reco::BeamSpot & vBspot,
                  const edm::ParameterSet & runParameters,
                  edm::ESHandle<MagneticField> & vbField ):
leptons(vleptons), tracks(vtracks), Bspot(vBspot), bField(vbField)
{
 ptLep1Cut_ = runParameters.getParameter<double>("ptLep1Cut");
 mllMin_ = runParameters.getParameter<double>("mllMin");
 mllMax_ = runParameters.getParameter<double>("mllMax");
 MBMin_ = runParameters.getParameter<double>("MBMin");
 MBMax_ = runParameters.getParameter<double>("MBMax");
 drTrkLepCut_ = runParameters.getParameter<double>("drTrkLepCut");
}

template <typename T>
std::unique_ptr<std::vector<pat::CompositeCandidate>>  BtoXLLTemplateFitter<T>::ReconstructB(float trackMass)
{
  // this makes transient stuff
  PreFitter prefit(bField);
  auto tpairs = prefit.CreateTransientPairs(leptons,ptLep1Cut_,mllMin_,mllMax_);
  auto ttracks = prefit.CreateTransientObject(tracks);

  //provides also indices of good lep pirs
  auto indices = prefit.GetIndexPairs();

  auto cands = std::make_unique< std::vector<pat::CompositeCandidate> >();
  for (unsigned int ilep = 0; ilep < indices.size(); ++ilep){
    unsigned int ilep1(indices[ilep].first);
    unsigned int ilep2(indices[ilep].second);
    ROOT::Math::LorentzVector lep1 ((leptons[ilep1]).p4());
    ROOT::Math::LorentzVector lep2 ((leptons[ilep2]).p4());


    for ( const pat::CompositeCandidate & trk : tracks){
      ROOT::Math::LorentzVector Kaon (trk.p4());
      //change mass E^2new= E^2old - m^2old + m^2new
      Kaon.SetE(TMath::Sqrt(pow(Kaon.E(),2)+pow(trackMass,2)-pow(Kaon.M(),2)) );

      // changing mass 
      if ( (Kaon + lep1 + lep2 ).M() < MBMin_ ||
           (Kaon + lep1 + lep2 ).M() > MBMax_ ||
           (Kaon + lep1 + lep2 ).Pt()< BPtCut_ )
             continue;

      // dr cut makes sense only on same charge objects
      if ( deltaR(Kaon.Eta(),Kaon.Phi(),lep1.Eta(),lep1.Phi()) < drTrkLepCut_ &&
           trk.charge()==leptons[ilep].charge() )
             continue;
      if ( deltaR(Kaon.Eta(),Kaon.Phi(),lep2.Eta(),lep2.Phi()) < drTrkLepCut_ &&
           trk.charge()==leptons[ilep].charge() )
             continue;

      unsigned int itrk=&trk-&tracks[0];
      KinematicParticleFactoryFromTransientTrack pFactory;
      // cands
      std::vector<RefCountedKinematicParticle> allParticles{
           pFactory.particle(tpairs.first[ilep],lep1.M(),chi,ndf,part_sigma),
           pFactory.particle(tpairs.second[ilep],lep2.M(),chi,ndf,part_sigma),
           pFactory.particle(ttracks[itrk],Kaon.M(),chi,ndf,part_sigma)
           };
      // actual fit    
      KinFitter fitter(allParticles);
      if (!fitter.success()) continue;

      //reconstruct B
      pat::CompositeCandidate cand;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> Bp4(
				   double(fitter.Mother_Momentum(true).x()),
			           double(fitter.Mother_Momentum(true).y()),
				   double(fitter.Mother_Momentum(true).z()),
				   double(fitter.Mother_Energy(true))   );


      cand.setP4(Bp4);
      cand.setCharge(fitter.Mother_Charge());
      cand.setVertex( reco::Candidate::Point(fitter.Mother_XYZ().x(),
                                             fitter.Mother_XYZ().y(),
                                             fitter.Mother_XYZ().z())  );
      cand.addUserFloat("vtx_prob", fitter.prob());
      cand.addUserFloat("vtx_chi2", fitter.chi());
      cand.addUserFloat("vtx_ndof", fitter.dof());
            
      GlobalPoint dist=FlightDistVector(Bspot,fitter.Mother_XYZ());
      cand.addUserFloat("Lxy", dist.perp());
      cand.addUserFloat("eLxy", TMath::Sqrt(fitter.Mother_XYZError().rerr(dist)));
      cand.addUserFloat("cosA",CosA(dist,Bp4));      
      cand.addUserInt("l1_id", ilep1 );
      cand.addUserInt("l2_id", ilep2 );
      cand.addUserInt("k_idx", itrk  );
     // unsigned int idaughter =0;
     /* for ( std::string daughter : {"l1refit_", "l2refit_", "Hadron_refit_" } ){
        cand.addUserFloat(daughter+"pt",fitter.Daughter_Momentum(idaughter,true).perp());
        cand.addUserFloat(daughter+"eta",fitter.Daughter_Momentum(idaughter,true).eta());
        cand.addUserFloat(daughter+"phi",fitter.Daughter_Momentum(idaughter,true).phi());
        idaughter++;
      }*/

      cands->emplace_back(cand);
    }
  }
 
  return std::move(cands);
}


