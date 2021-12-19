#ifndef ConversionSelector_H
#define ConversionSelector_H
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//
#include <map>
#include <vector>

// forward declarations
class TFile;
class TTree;
/** \class ConversionSelector
 **
 **
 **  $Id: ConversionSelector
 **  \author N.Marinelli - Univ. of Notre Dame
 **
 ***/


class ConversionSelector : public edm::EDProducer {

 public:

  //
  explicit ConversionSelector( const edm::ParameterSet& ) ;
  virtual ~ConversionSelector();

  void produce(edm::Event&, const edm::EventSetup&) ;

 private:
  //

  float mee (float ipx1, float ipy1, float ipz1, float ipx2, float ipy2, float ipz2);

  TTree *tree_;
  TFile *histfile_;
  std::string outputFileName_;

  std::string fName_;
  edm::ESHandle<MagneticField> theMF_;

  int verbosity_;
  int nEvt_;
  int nEntry_;
  int nSimConv_[2];
  int nMatched_;
  int nRecConv_;
  int nRecConvAss_;
  int nRecConvAssWithEcal_;
  int nInvalidPCA_;

  int nAcceptedConvTot;


  edm::ParameterSet parameters_;
  edm::ESHandle<CaloGeometry> theCaloGeom_;
  edm::ESHandle<CaloTopology> theCaloTopo_;

  std::string conversionCollectionProducer_;
  std::string conversionCollection_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionCollectionPr_Token_;

  std::string filteredConvertedPhotonCollection_;

  std::string conversionTrackProducer_;

  std::string photonCollectionProducer_;
  std::string photonCollection_;
  edm::EDGetTokenT<reco::PhotonCollection> photonCollectionPr_Token_;

  std::string lowPtElectronProducer_;
  std::string lowPtElectronCollection_;
  edm::EDGetTokenT<std::vector<pat::Electron> > lowPtElectronCollectionPr_Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> unBiased_;
  edm::EDGetTokenT<edm::ValueMap<float>> mvaId_src_;


  edm::EDGetTokenT<reco::VertexCollection> offline_pvToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  edm::EDGetTokenT<double> rho_token_;
  edm::EDGetTokenT<edm::SimTrackContainer> g4_simTk_Token_;
  edm::EDGetTokenT<edm::SimVertexContainer> g4_simVtx_Token_;
  edm::EDGetTokenT<TrackingParticleRefVector> tpSelForEff_Token_;
  edm::EDGetTokenT<TrackingParticleRefVector> tpSelForFake_Token_;
  edm::EDGetTokenT<edm::HepMCProduct> hepMC_Token_;
  edm::EDGetTokenT<reco::GenJetCollection> genjets_Token_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> trackAssociator_Token_;


  edm::InputTag label_tp_;

  PhotonMCTruthFinder*  thePhotonMCTruthFinder_;

  bool isRunCentrally_;
  double minPtOnTracks_;
  double minPhoEtCut_;
  double trkIsolExtRadius_;
  double trkIsolInnRadius_;
  double trkPtLow_;
  double lip_;
  double ecalIsolRadius_;
  double bcEtLow_;
  double hcalIsolExtRadius_;
  double hcalIsolInnRadius_;
  double hcalHitEtLow_;
  int  numOfTracksInCone_;
  double trkPtSumCut_;
  double ecalEtSumCut_;
  double hcalEtSumCut_;
  bool dCotCutOn_;
  double dCotCutValue_;
  double dCotHardCutValue_;
  bool generalTracksOnly_;
  bool gsfTracksOpenOnly_;
  bool arbitratedMerged_;
  bool arbitratedEcalSeeded_;
  bool ecalalgotracks_;
  bool highPurity_;
  double minProb_;
  uint maxHitsBeforeVtx_;
  double minLxy_;


  /// global variable for the MC photon
  double mcPhi_;
  double mcEta_;
  double mcConvPt_;
  double mcConvR_;
  double mcConvZ_;
  double mcConvY_;
  double mcConvX_;
  double mcConvPhi_;
  double mcConvEta_;
  double mcJetEta_;
  double mcJetPhi_;

  edm::RefVector<TrackingParticleCollection> theConvTP_;
  //std::vector<TrackingParticleRef>    theConvTP_;

  double minPhoPtForEffic;
  double maxPhoEtaForEffic;
  double maxPhoZForEffic;
  double maxPhoRForEffic;
  double minPhoPtForPurity;
  double maxPhoEtaForPurity;
  double maxPhoZForPurity;
  double maxPhoRForPurity;

  double simMinPt_;
  double simMaxPt_;

  /// Global variables for reco Photon
  double recMinPt_;
  double recMaxPt_;

  // variables in the tree per conversion
  double rho_;
  int convAlgo_;
  int convMatchSC_;
  int isTwoSuperClusterConversion_;
  float vtxChiSquaredProb_;
  double lxy_;
  double vtxR_;
  double vtxX_;
  double vtxY_;
  double vtxZ_;


  float etaBeforeFit_;
  float phiBeforeFit_;
  float ptBeforeFit_;
  float etaAfterFit_;
  float phiAfterFit_;
  float ptAfterFit_;


  float invMass_fromDF_;
  float invMass_fromPin_;
  float invMass_beforeFit_;
  float invMass_afterFit_;
  float leadTrkPt_;
  float subleadTrkPt_;
  float leadTrkEta_;
  float subleadTrkEta_;
  float leadTrkTheta_;
  float subleadTrkTheta_;
  float dCotTracksFromPin_;
  float dCotTracksBeforeVtxFit_;
  float dCotTracksAfterVtxFit_ ;
  float maxNHitsBeforeVtx_;
  float leadNHitsBeforeVtx_;
  float trailNHitsBeforeVtx_;
  float distMinAppTracks_;
  float dPhiTracksAtVtx_;
  float dEtaTracksAtECAL_whenTwoSC_;
  float dPhiTracksAtECAL_whenTwoSC_;
  float dEtaSC_whenTwoSC_;
  float dPhiSC_whenTwoSC_;
  float dEtaTracksAtECAL_whenOneSC_;
  float dPhiTracksAtECAL_whenOneSC_;


  // variables in the tree for all tracks
  int nTracks_;
  std::vector<float> trk_pt_beforeVtx_;
  std::vector<float> trk_eta_beforeVtx_;
  std::vector<float> trk_phi_beforeVtx_;


  std::vector<float> trk_ptMode_;
  std::vector<float> trk_etaMode_;
  std::vector<float> trk_phiMode_;



  std::vector<float> trk_pt_afterVtx_;
  std::vector<float> trk_eta_afterVtx_;
  std::vector<float> trk_phi_afterVtx_;

  std::vector<float> trk_found_nhits_;
  std::vector<float> trk_valid_nhits_;
  std::vector<float> trk_dxy_;
  std::vector<float> trk_dxy_err_;
  std::vector<float> trk_pin_;
  std::vector<float> trk_pout_;
  std::vector<float> trk_chi2_;
  std::vector<float> trk_d0_;
  std::vector<float> trk_charge_;
  std::vector<float> trk_chargeMode_;



  std::vector<float> dEtaTrkSCAtVtx_;
  std::vector<float> dPhiTrkSCAtVtx_;
  std::vector<float> dEtaTrkSCAtEcal_;
  std::vector<float> dPhiTrkSCAtEcal_;
  std::vector<float> isMatchedToLowPtSC_;
  
  std::vector<float> EoverP_;
  std::vector<float> BremFrac_;
 
  std::vector<float> elePt_;
  std::vector<float> eleSC_energy_;
  std::vector<float> eleSC_eta_;
  std::vector<float> eleSC_clusterSize_;
  std::vector<float> eleSC_etaWidth_;
  std::vector<float> eleSC_phiWidth_;
  std::vector<float> full5x5_sigmaIetaIeta_;
  std::vector<float> full5x5_sigmaIphiIphi_;
  std::vector<float> full5x5_circularity_;
  std::vector<float> full5x5_e1x5_;
  std::vector<float> full5x5_e2x5Max_;
  std::vector<float> full5x5_e5x5_;
  std::vector<float> full5x5_r9_;
  std::vector<float> full5x5_HoverE_;



};




#endif
