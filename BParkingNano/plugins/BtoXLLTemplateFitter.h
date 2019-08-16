// Template class that can produce general XYY were Y are opposite sign leptons or tracks and X is a hadron

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

template< typename T >
class BtoXLLTemplateFitter {

 private:
  // containers
  const std::vector<T>& leptons;
  const std::vector<pat::CompositeCandidate>& tracks;
  const reco::BeamSpot& Bspot;
  edm::ESHandle< MagneticField > bField;
  std::pair< std::vector< T >, std::vector< T > > pairs;
  
  // selection
  double ptLep1Cut_   = 0;  double mllMin_ = 0;  double mllMax_ = 0; 
  double MBMin_       = 0;  double MBMax_  = 0;  double BPtCut_ = 0;
  double drTrkLepCut_ = 0;  double bdtEl1Cut_ = 0;
  double BVtxProbCut_ = -1; double dzLepLepCut_ =0;
  double BCosACut_    = -1.01;
  
  // needed for fit 
  float chi = 0, ndf = 0, part_sigma = 0, kaon_sigma = 0;


 public:
  //constructor
  BtoXLLTemplateFitter( const std::vector< T >& vleptons,
                        const std::vector< pat::CompositeCandidate >& vtracks, 
                        const reco::BeamSpot& vBspot,
                        const edm::ParameterSet& runParameters,
                        edm::ESHandle< MagneticField >& bFieldHandle );
  // destructor
  ~BtoXLLTemplateFitter(){}

  //only one me,ber function
  std::unique_ptr< std::vector<pat::CompositeCandidate> > ReconstructB(
                                                          float trackMass );
  
};
#endif

//////////////////////////////////////////////////////////////////////////////

template<typename T> 
BtoXLLTemplateFitter< T >::BtoXLLTemplateFitter(
                  const std::vector< T >& vleptons,
                  const std::vector< pat::CompositeCandidate >& vtracks, 
                  const reco::BeamSpot& vBspot,
                  const edm::ParameterSet& runParameters,
                  edm::ESHandle< MagneticField >& vbField ):
   leptons( vleptons ), tracks( vtracks ), Bspot( vBspot ), bField( vbField )
{

 // get parameters
 ptLep1Cut_   = runParameters.getParameter<double>( "ptLep1Cut" );
 bdtEl1Cut_   = runParameters.getParameter<double>( "bdtEl1Cut" );
 dzLepLepCut_ = runParameters.getParameter<double>( "dzLepLepCut" );
 mllMin_      = runParameters.getParameter<double>( "mllMin" );
 mllMax_      = runParameters.getParameter<double>( "mllMax" );
 MBMin_       = runParameters.getParameter<double>( "MBMin" );
 MBMax_       = runParameters.getParameter<double>( "MBMax" );
 drTrkLepCut_ = runParameters.getParameter<double>( "drTrkLepCut" );
 BPtCut_      = runParameters.getParameter<double>( "BPtCut" );
 BVtxProbCut_ = runParameters.getParameter<double>( "BVtxProbCut" );
 BCosACut_    = runParameters.getParameter<double>( "BCosACut" );

}

template < typename T >
std::unique_ptr< std::vector<pat::CompositeCandidate> >
BtoXLLTemplateFitter< T >::ReconstructB( float trackMass )
{
  // first we make everything transient to avoid invoking the builder many times

  // templated class; creates transient objects and applies cuts
  PreFitter prefit(bField);

  // pair creation returns transient tracks cuts ( ptL1 bdt mllmin mllmax dzll)
  auto tpairs = prefit.CreateTransientPairs( leptons, ptLep1Cut_, bdtEl1Cut_,
                                             mllMin_, mllMax_, dzLepLepCut_ );

  // transient tracks from composite part
  auto ttracks = prefit.CreateTransientObject( tracks );

  // provides also indices of good lep pairs
  auto indices = prefit.GetIndexPairs();

  // output
  auto cands = std::make_unique< std::vector<pat::CompositeCandidate> >();

  // loop on the pairs of indices of leptons
  for ( unsigned int ilep = 0; ilep < indices.size(); ++ilep ){
    unsigned int ilep1( indices[ilep].first );
    unsigned int ilep2( indices[ilep].second );

    // get 4 vector
    ROOT::Math::LorentzVector lep1( (leptons[ilep1]).p4() );
    ROOT::Math::LorentzVector lep2( (leptons[ilep2]).p4() );

    // now loop over tracks and reconstruct  B
    for ( const pat::CompositeCandidate & trk : tracks ){
      
      // get tracks 4 vector
      ROOT::Math::LorentzVector Kaon( trk.p4() );

      //change mass according to E^2new= E^2old - m^2old + m^2new
      Kaon.SetE( TMath::Sqrt( pow(Kaon.E(),2) + pow(trackMass,2) 
                            - pow(Kaon.M(),2) ) );

    
      // cut on candidate mass and pT
      if ( ( Kaon + lep1 + lep2 ).M() < MBMin_ ||
           ( Kaon + lep1 + lep2 ).M() > MBMax_ ||
           ( Kaon + lep1 + lep2 ).Pt() < BPtCut_ )
           continue;
    
      // dr cut (trk,lep) only on same charge objects
      if ( deltaR( Kaon.Eta(), Kaon.Phi(), lep1.Eta(), lep1.Phi() ) < drTrkLepCut_ &&
           trk.charge() == leptons[ilep].charge() )
           continue;
      if ( deltaR( Kaon.Eta(), Kaon.Phi(), lep2.Eta(), lep2.Phi() ) < drTrkLepCut_  &&
           trk.charge() == leptons[ilep].charge() )
           continue;
     
      unsigned int itrk( &trk - &tracks[0] );

      KinematicParticleFactoryFromTransientTrack pFactory;

      // B candidate creation
      std::vector< RefCountedKinematicParticle > allParticles{
           pFactory.particle( tpairs.first[ilep], lep1.M(), chi,
                              ndf, part_sigma ),
           pFactory.particle( tpairs.second[ilep], lep2.M(), chi, 
                              ndf, part_sigma ),
           pFactory.particle( ttracks[itrk], Kaon.M(), chi, ndf, part_sigma )
           };

      // actual fit and track refitting done in this class. bool if set to true refits daughters in the common vtx and B is calculated from refitted daughters
      KinFitter fitter( allParticles, true );

      // skip unseccesfull fits
      if ( !fitter.success() ) continue;

      // apply cuts on probability and cos
      if ( fitter.prob() < BVtxProbCut_ ) continue;
  
      // cos A from Beamspot (Bspot)
      GlobalPoint dist = FlightDistVector( Bspot, fitter.Mother_XYZ() );
      // take 4 vector
      ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<double> > Bp4(
				   double( fitter.Mother_Momentum().x() ),
			           double( fitter.Mother_Momentum().y() ),
				   double( fitter.Mother_Momentum().z() ),
				   double( fitter.Mother_Energy() )   
                                   );

      float cos = CosA( dist, Bp4 ); // cos function in helper.h
      if ( cos < BCosACut_ ) continue;

      //create and save B candidate
      pat::CompositeCandidate cand;

      cand.setP4( Bp4 );
      cand.setCharge( fitter.Mother_Charge() );
      cand.setVertex( reco::Candidate::Point( fitter.Mother_XYZ().x(),
                                              fitter.Mother_XYZ().y(),
                                              fitter.Mother_XYZ().z()
                                              )  );
      // B specific vars
      cand.addUserFloat( "vtx_prob", fitter.prob() );
      cand.addUserFloat( "vtx_chi2", fitter.chi() );
      cand.addUserFloat( "vtx_ndof", fitter.dof() );
      cand.addUserFloat( "Lxy", dist.perp() );
      cand.addUserFloat( "eLxy", TMath::Sqrt( fitter.Mother_XYZError().rerr(
                                              dist
                                              ) ) );
      cand.addUserFloat( "cosA", cos ); 

      // dautghter stuff    
      cand.addUserInt( "l1_idx", ilep1 );
      cand.addUserInt( "l2_idx", ilep2 );
      cand.addUserInt( "k_idx", itrk  );

      unsigned int idaughter = 0;
      for ( std::string daughter: { "l1refit_", "l2refit_", "Hrefit_" } ) {

        cand.addUserFloat( daughter + "pt", 
                           fitter.Daughter_Momentum( idaughter ).perp() );
        cand.addUserFloat( daughter + "eta", 
                           fitter.Daughter_Momentum( idaughter ).eta() );
        cand.addUserFloat( daughter + "phi",
                           fitter.Daughter_Momentum( idaughter ).phi() );
        idaughter++;
      }
      // compute new vectors 
      lep1.SetPx ( fitter.Daughter_Momentum( 0 ).x() );
      lep1.SetPy( fitter.Daughter_Momentum( 0 ).y() );
      lep1.SetPz( fitter.Daughter_Momentum( 0 ).z() );
      lep2.SetPx ( fitter.Daughter_Momentum( 1 ).x() );
      lep2.SetPy( fitter.Daughter_Momentum( 1 ).y() );
      lep2.SetPz( fitter.Daughter_Momentum( 1 ).z() );
      cand.addUserFloat( "mll", ( lep1 + lep2 ).M() );
      cands->emplace_back( cand );
    }
  }
  // candidates 
  return std::move( cands );
}


