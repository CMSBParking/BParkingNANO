#include <iostream>
//
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
#include "ConversionSelector.h"

//
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
//#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
//#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

//
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionHitChecker.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//
//
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"
#include "TProfile.h"
//
/** \class ConversionSelector
 **
 **
 **  $Id: ConversionSelector
 **  \author N.Marinelli - Univ. of Notre Dame
 **
 ***/

using namespace std;


ConversionSelector::ConversionSelector( const edm::ParameterSet& pset )
  {
    outputFileName_ = pset.getParameter<std::string>("outFileName");
    verbosity_ = pset.getUntrackedParameter<int>("Verbosity");
    parameters_ = pset;

    lowPtElectronCollection_ = pset.getParameter<std::string>("lowPtElectronCollection");
    lowPtElectronProducer_ = pset.getParameter<std::string>("lowPtElectronProducer");
    lowPtElectronCollectionPr_Token_ = consumes<std::vector<pat::Electron> >(edm::InputTag(lowPtElectronProducer_,lowPtElectronCollection_));

    conversionCollectionProducer_ = pset.getParameter<std::string>("convProducer");
    conversionCollection_ = pset.getParameter<std::string>("conversionCollection");
    conversionCollectionPr_Token_ = consumes<reco::ConversionCollection> (
        edm::InputTag(conversionCollectionProducer_,
                      conversionCollection_));

    unBiased_ = consumes<edm::ValueMap<float>>( pset.getParameter<edm::InputTag>("unbiasedSeeding") );
    mvaId_src_ = consumes<edm::ValueMap<float>>( pset.getParameter<edm::InputTag>("mvaId") );

    // use configuration file to setup output collection names
    filteredConvertedPhotonCollection_ = pset.getParameter<std::string>("filteredConversions");
  
   // Register the product
    produces<pat::CompositeCandidateCollection>(filteredConvertedPhotonCollection_);

    minPtOnTracks_ = pset.getParameter<double>("minPtOnTracks");
    minPhoEtCut_ = pset.getParameter<double>("minPhoEtCut");
    generalTracksOnly_ = pset.getParameter<bool>("generalTracksOnly");
    gsfTracksOpenOnly_ = pset.getParameter<bool>("gsfTracksOpenOnly");
    arbitratedMerged_ = pset.getParameter<bool>("arbitratedMerged");
    arbitratedEcalSeeded_ = pset.getParameter<bool>("arbitratedEcalSeeded");
    ecalalgotracks_ = pset.getParameter<bool>("ecalalgotracks");
    highPurity_ = pset.getParameter<bool>("highPurity");
    minProb_ = pset.getParameter<double>("minProb");
    maxHitsBeforeVtx_ = pset.getParameter<uint>("maxHitsBeforeVtx");
    minLxy_           = pset.getParameter<double>("minLxy");
 
    
    offline_pvToken_ = consumes<reco::VertexCollection>(pset.getUntrackedParameter<edm::InputTag> ("offlinePV",
											  edm::InputTag("offlineSlimmedPrimaryVertices")));

    beamspotToken_ = consumes<reco::BeamSpot>(pset.getUntrackedParameter<edm::InputTag> ("beamspot",edm::InputTag("offlineBeamSpot")));
    rho_token_ = consumes<double>(pset.getParameter<string> ("rhoHandle"));
  
 
    nEvt_=0;
    isMatchedToLowPtSC_.reserve(2);

    nAcceptedConvTot=0;

  }



ConversionSelector::~ConversionSelector() {
    

  edm::LogInfo("ConversionSelector") << "Analyzed " << nEvt_  << "\n";
  // std::cout  << "::endJob Analyzed " << nEvt_ << " events " << " with total " << nPho_ << " Photons " << "\n";
  std::cout  << "ConversionSelector::DTOR Analyzed " << nEvt_ << " events " << "\n";



}




void ConversionSelector::produce(edm::Event& e, const edm::EventSetup& esup ) {

 
  using namespace edm;
  //  const float etaPhiDistance=0.01;
  // Fiducial region
  // const float TRK_BARL =0.9;
  const float BARL = 1.4442; // DAQ TDR p.290
  //  const float END_LO = 1.566; // unused
  const float END_HI = 2.5;
  // Electron mass
  //  const Float_t mElec= 0.000511; // unused

  isMatchedToLowPtSC_.clear();

  pat::CompositeCandidateCollection outFilteredConversionCollection;
  auto outFilteredConversionCollection_p = std::make_unique<pat::CompositeCandidateCollection>();


  nEvt_++;
  LogInfo("ConversionSelector") << "ConversionSelector Analyzing event number: " << e.id() << " Global Counter " << nEvt_ <<"\n";
  


  ///// Get the recontructed  conversions
  Handle<reco::ConversionCollection> convHandle;
  e.getByToken(conversionCollectionPr_Token_, convHandle);
  const reco::ConversionCollection convCollection = *(convHandle.product());
  if (!convHandle.isValid()) {
    edm::LogError("ConversionSelector") << "Error! Can't get the  collection "<< std::endl;
    return;
  }


  ///// Get the recontructed low pt gsfElectrons
  Handle<std::vector<pat::Electron> > lowPtEleHandle;
  e.getByToken(lowPtElectronCollectionPr_Token_, lowPtEleHandle);
  std::vector<pat::Electron> lowPtElectronCollection = *(lowPtEleHandle.product());
  if (!lowPtEleHandle.isValid()) {
    edm::LogError("ConversionSelector") << "Error! Can't get the lowPt electron collection "<< std::endl;
    return;
  }

  edm::Handle<edm::ValueMap<float> > unBiased;
  e.getByToken(unBiased_, unBiased);

  edm::Handle<edm::ValueMap<float> > mvaId;  
  e.getByToken(mvaId_src_, mvaId);



  // offline  Primary vertex
  /*  
bool valid_pvtx = false; 
  edm::Handle<reco::VertexCollection> vertexHandle;
  reco::VertexCollection vertexCollection;
  e.getByToken(offline_pvToken_, vertexHandle);
  if (!vertexHandle.isValid()) {
    edm::LogError("ConvAnalysisMiniAod") << "Error! Can't get the product primary Vertex Collection "<< "\n";
  } else {
    vertexCollection = *(vertexHandle.product());
  }

  reco::Vertex pvtx;
  if (!vertexCollection.empty()){
    pvtx = *(vertexCollection.begin());
    //asking for one good vertex
    if (pvtx.isValid() && fabs(pvtx.position().z())<=15 && pvtx.position().Rho()<=2){
      valid_pvtx = true;
    }
  }
  */
 


  edm::Handle<reco::BeamSpot> bsHandle;
  e.getByToken(beamspotToken_, bsHandle);
  if (!bsHandle.isValid()) {
      edm::LogError("TrackerOnlyConversionProducer")
          << "Error! Can't get the product primary Vertex Collection "<< "\n";
      return;
  }
  const reco::BeamSpot &thebs = *bsHandle.product();


  edm::Handle<double> rhoHandle;
  e.getByToken(rho_token_, rhoHandle);
  rho_ = *rhoHandle;

  edm::ESHandle<TransientTrackBuilder> theB ;
  esup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  

  //std::cout << " ConversionSelector Conversion collection size " << convHandle->size() << std::endl;
  int nAcceptedConv=0;
  for (reco::ConversionCollection::const_iterator conv = convHandle->begin();conv!=convHandle->end();++conv) {
    const reco::Conversion aConv = (*conv);

    if ( arbitratedMerged_ && !aConv.quality(reco::Conversion::arbitratedMerged)  ) continue;
    if ( generalTracksOnly_ && !aConv.quality(reco::Conversion::generalTracksOnly) ) continue;
    if ( arbitratedEcalSeeded_ && !aConv.quality(reco::Conversion::arbitratedEcalSeeded)  ) continue;
    if ( highPurity_ && !aConv.quality(reco::Conversion::highPurity) ) continue;

     
    std::vector<edm::RefToBase<reco::Track> > tracks = aConv.tracks();
    RefToBase<reco::Track> tk1 = aConv.tracks().front();
    RefToBase<reco::Track> tk2 = aConv.tracks().back();
    const reco::Track refTk1 = aConv.conversionVertex().refittedTracks().front();
    const reco::Track refTk2 = aConv.conversionVertex().refittedTracks().back();
    const reco::Vertex& vtx = aConv.conversionVertex();
   
 
    //requires two tracks and a valid vertex
    if (tracks.size() !=2 || !(vtx.isValid())) continue;
    if ( tk1->pt() < minPtOnTracks_ || tk2->pt() < minPtOnTracks_) continue;   


    //compute transverse decay length with respect to beamspot
    math::XYZVectorF refittedMom =  aConv.refittedPairMomentum();
    double dbsx = aConv.conversionVertex().x() - thebs.x0();
    double dbsy = aConv.conversionVertex().y() - thebs.y0();
    double lxy = (refittedMom.x()*dbsx + refittedMom.y()*dbsy)/refittedMom.rho();
    lxy_=lxy;
    //    if (lxy<minLxy_) continue;
     

    if ( ecalalgotracks_ && ( !(tk1->algo()==reco::TrackBase::outInEcalSeededConv || tk1->algo()==reco::TrackBase::inOutEcalSeededConv) || !(tk2->algo()==reco::TrackBase::outInEcalSeededConv || tk2->algo()==reco::TrackBase::inOutEcalSeededConv)  )  ) continue;

    nAcceptedConv++;
    
   
   
    //TODO replace it with phi at vertex
    //float  dPhiTracksAtVtx =  aConv.dPhiTracksAtVtx();
    
    ///////////  Quantities per conversion
    
    float invM=aConv.pairInvariantMass();
    invMass_fromDF_=invM;

    float chi2Prob = ChiSquaredProbability( aConv.conversionVertex().chi2(),  aConv.conversionVertex().ndof() );
    uint maxNHitsBeforeVtx = aConv.nHitsBeforeVtx().size()>1 ? max(aConv.nHitsBeforeVtx().at(0),aConv.nHitsBeforeVtx().at(1)) : 0;
    uint sumNHitsBeforeVtx = aConv.nHitsBeforeVtx().size()>1 ? aConv.nHitsBeforeVtx().at(0) + aConv.nHitsBeforeVtx().at(1) : 0;

 
    if (chi2Prob <= 0.0005) continue;
    float vtxR = sqrt(aConv.conversionVertex().position().perp2());
    if (vtxR < 1.5 || vtxR > 4) continue;
    
   

    unsigned int ilead = 0, itrail = 1;
    // use pt of tracks before refitting
    if (tk2->pt() > tk1->pt()) {
      ilead = 1;
      itrail = 0;
    }

    //std::cout << " ilead " << ilead << " itrail " << itrail << std::endl;

    math::XYZVectorF leadPin = aConv.tracksPin().at(ilead);
    math::XYZVectorF trailPin = aConv.tracksPin().at(itrail);

    RefToBase<reco::Track> tkleadBeforeVtxFit = aConv.tracks().at(ilead);
    RefToBase<reco::Track> tktrailBeforeVtxFit = aConv.tracks().at(itrail);

    const reco::Track tklead  = aConv.conversionVertex().refittedTracks().at(ilead);
    const reco::Track tktrail = aConv.conversionVertex().refittedTracks().at(itrail);

    // calculate invariant mass before the fit to common vertex from Pin 
    float invMass_fromPin=mee( leadPin.x(), 
			       leadPin.y(),
			       leadPin.z(),
			       trailPin.x(),
			       trailPin.y(),
			       trailPin.z());
    invMass_fromPin_=invMass_fromPin;

    // calculate invariant mass before the fit to common vertex

    float invMass_beforeFit= mee ( tkleadBeforeVtxFit->px(), 
				   tkleadBeforeVtxFit->py(), 
				   tkleadBeforeVtxFit->pz(), 
				   tktrailBeforeVtxFit->px(), 
				   tktrailBeforeVtxFit->py(), 
				   tktrailBeforeVtxFit->pz() );
    invMass_beforeFit_=invMass_beforeFit;

    // calculate invariant mass after the fit to common vertex
    float invMass_afterFit=mee( tklead.px(), 
				tklead.py(), 
				tklead.pz(), 
				tktrail.px(), 
				tktrail.py(), 
				tktrail.pz());
    invMass_afterFit_ = invMass_afterFit;


    int deltaExpectedHitsInner = tklead.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)
      - tktrail.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    //int leadExpectedHitsInner = tklead.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    uint leadNHitsBeforeVtx = aConv.nHitsBeforeVtx().size()>1 ? aConv.nHitsBeforeVtx().at(ilead) : 0;
    uint trailNHitsBeforeVtx = aConv.nHitsBeforeVtx().size()>1 ? aConv.nHitsBeforeVtx().at(itrail) : 0;
    float dCotFromPin = 1./tan(trailPin.theta())-1./tan(leadPin.theta());
  

    pat::CompositeCandidate cand;
    math::PtEtaPhiMLorentzVector convp4(
					aConv.pairMomentum().perp2(), 
					aConv.pairMomentum().eta(),
					aConv.pairMomentum().phi(),
					invMass_fromPin
					);
    
    cand.setP4(convp4);
    cand.addUserFloat("rho",rho_);
    cand.addUserFloat("vtxChi2Prob",chi2Prob);
    cand.addUserFloat("vtxR",sqrt(aConv.conversionVertex().position().perp2()));
    cand.addUserFloat("vtxX",aConv.conversionVertex().position().x());
    cand.addUserFloat("vtxY",aConv.conversionVertex().position().y());
    cand.addUserFloat("vtxZ",aConv.conversionVertex().position().z());
    cand.addUserFloat("dCotFromPin",dCotFromPin);
    cand.addUserFloat("dMinAppTracks",aConv.distOfMinimumApproach());
    cand.addUserFloat("dPhiTracksAtVtx", aConv.dPhiTracksAtVtx());
    cand.addUserFloat("maxNHitsBeforeVtx",maxNHitsBeforeVtx);
    cand.addUserFloat("leadNHitsBeforeVtx",leadNHitsBeforeVtx);
    cand.addUserFloat("trailNHitsBeforeVtx",trailNHitsBeforeVtx);
    cand.addUserFloat("sumNHitsBeforeVtx",sumNHitsBeforeVtx);
    cand.addUserFloat("deltaExpectedHitsInner",deltaExpectedHitsInner);

    
    
    
    std::vector<float> unbiased_seedBDT;
    std::vector<float>  mvaIdV;

    
    std::vector<reco::GsfElectron> eleV;
    
    std::vector<float> elePt;    
    std::vector<float> eleCh;    
    std::vector<float> eleEta;    
    std::vector<float> elePhi;    
    std::vector<float> ele_seed_dEta;
    std::vector<float> ele_eclu_EoverP;
    std::vector<float> ele_sc_dEta;
    std::vector<float> ele_sc_dPhi;
    std::vector<float> ele_fBrem;
    std::vector<float> ele_shFrInHits;

    std::vector<int> cluSize;    
    std::vector<float> gsfPmode;    
    std::vector<float> gsfChi2;
    std::vector<float> gsfNhits;
    std::vector<float> gsfDR;

    std::vector<float> scEnergy;
    std::vector<float> scEta;
    std::vector<float> scEtaWidth;
    std::vector<float> scPhiWidth;
    std::vector<float> full5x5_r9;    
    std::vector<float> full5x5_HoverE;    
    std::vector<float> kf_Chi2;
    std::vector<float> kf_p;
    std::vector<float> kf_nHits;
    std::vector<float> kf_DR;



    
    eleV.resize(2);    
   
    elePt.resize(2);
    eleCh.resize(2);
    eleEta.resize(2);
    elePhi.resize(2);

    ele_seed_dEta.resize(2);
    ele_eclu_EoverP.resize(2);
    ele_sc_dEta.resize(2);
    ele_sc_dPhi.resize(2);
    ele_fBrem.resize(2);
    ele_shFrInHits.resize(2);
    cluSize.resize(2);
    scEnergy.resize(2);
    scEta.resize(2);
    scEtaWidth.resize(2);
    scPhiWidth.resize(2);
    full5x5_HoverE.resize(2);
    full5x5_r9.resize(2);
    kf_Chi2.resize(2);
    kf_p.resize(2);
    kf_nHits.resize(2);
    kf_DR.resize(2);
    
    gsfPmode.resize(2);
    gsfChi2.resize(2);
    gsfNhits.resize(2);
    gsfDR.resize(2);

    unbiased_seedBDT.resize(2);
    mvaIdV.resize(2);

    // check here the -hard- matching between the tracks and lowPt electron objects. Should be 100% by construction
    int nMatchedTracks=0;
    int iOrd=0;
    for ( unsigned int i=0; i<tracks.size(); i++) {

      if ( i == ilead ) iOrd=ilead;
      else iOrd=itrail;
      //      std::cout << " track pt " << tracks[iOrd]->pt();
      //std::cout << " Track " << iOrd << " tracks[i].id() " <<  tracks[iOrd].id() << " key " <<  tracks[iOrd].key() << std::endl;
      
      isMatchedToLowPtSC_[iOrd]=0;
      // std::cout << " gsfEle size " <<  lowPtElectronCollection.size() << std::endl;
    

      int idx=0;
      for(  std::vector<pat::Electron>::const_iterator iEle = lowPtElectronCollection.begin(); iEle != lowPtElectronCollection.end(); iEle++) {
	edm::Ref<pat::ElectronCollection> ref(lowPtEleHandle,idx);
        idx++;
	
	pat::Electron  myPatEle = pat::Electron(*iEle);
	reco::GsfElectron  myEle = reco::GsfElectron(*iEle);

	reco::GsfTrackRef gsfTrk= myEle.gsfTrack();

        if (gsfTrk.isNull()) continue;
	
	float mva_id = 999.;
	mva_id = float((*mvaId)[ref]);

	
	if ( gsfTrk.id() == tracks[iOrd].id() &&  gsfTrk.key() == tracks[iOrd].key() )  {
	  nMatchedTracks++;

          float clusize = myEle.superCluster()->clustersSize();
	  eleV[iOrd]=myEle;
        
          

	  //
	  elePt[iOrd]=myEle.pt();
	  eleCh[iOrd]=myEle.charge();
	  eleEta[iOrd]=myEle.eta();
	  elePhi[iOrd]=myEle.phi();
	  //
	  ele_seed_dEta[iOrd]= myEle.deltaEtaSeedClusterTrackAtCalo();
	  ele_eclu_EoverP[iOrd]=(1./myEle.ecalEnergy()) - (1./myEle.p());
	  ele_sc_dEta[iOrd]=myEle.deltaEtaSuperClusterTrackAtVtx();
	  ele_sc_dPhi[iOrd]=myEle.deltaPhiSuperClusterTrackAtVtx();
	  ele_fBrem[iOrd]= myEle.fbrem();
	  ele_shFrInHits[iOrd]=myEle.shFracInnerHits();

	  if ( myEle.core().isNonnull() ) {
	    cluSize[iOrd]=myEle.superCluster()->clustersSize();
	    scEta[iOrd] = myEle.superCluster()->eta();
	    scEnergy[iOrd] = myEle.superCluster()->energy();
	    scEtaWidth[iOrd] = myEle.superCluster()->etaWidth();
	    scPhiWidth[iOrd] = myEle.superCluster()->phiWidth();
	  }
	  full5x5_HoverE[iOrd] = myEle.full5x5_hcalOverEcal(); 
	  full5x5_r9[iOrd] = myEle.full5x5_r9();
          gsfPmode[iOrd]=gsfTrk->pMode();
          gsfChi2[iOrd]=gsfTrk->normalizedChi2();
          gsfNhits[iOrd]=(float)gsfTrk->found();

	  TVector3 gsfTV3(0,0,0);
	  gsfTV3.SetPtEtaPhi(gsfTrk->ptMode(), gsfTrk->etaMode(), gsfTrk->phiMode());  
	  TVector3 eleTV3(0,0,0);
	  eleTV3.SetPtEtaPhi(myEle.pt(), myEle.eta(), myEle.phi());
	  gsfDR[iOrd] = eleTV3.DeltaR(gsfTV3);  




	  unbiased_seedBDT[iOrd] =  float((*unBiased)[gsfTrk]);
          mvaIdV[iOrd]=mva_id;           


	  //	  std::cout << " Matching ele to track " << i <<std::endl;	  
          
	  //	  	  std::cout << " ********  gsfTrk id " << gsfTrk.id() << " TrackRef id " << tracks[iOrd].id() << " brem " << " gsfTrk key " <<  gsfTrk.key()  << " TrackRef key " << tracks[iOrd].key() << " brem " <<   myEle.trackFbrem() << std::endl;	
	  
	  
	  isMatchedToLowPtSC_[iOrd]=1;

 
	}      

      }



      
      
    } //loop over tracks


    cand.addUserFloat("tkLead_pt",aConv.tracks()[ilead]->pt());
    cand.addUserFloat("tkTrail_pt",aConv.tracks()[itrail]->pt());
    cand.addUserFloat("tkLead_eta",aConv.tracks()[ilead]->eta());
    cand.addUserFloat("tkTrail_eta",aConv.tracks()[itrail]->eta());
    cand.addUserFloat("tkLead_phi",aConv.tracks()[ilead]->phi());
    cand.addUserFloat("tkTrail_phi",aConv.tracks()[itrail]->phi());
    cand.addUserFloat("tkLead_nHits",aConv.tracks()[ilead]->found());
    cand.addUserFloat("tkTrail_nHits",aConv.tracks()[itrail]->found());
    cand.addUserInt("tkLead_isMatchedToLowPtEle",isMatchedToLowPtSC_[ilead]);
    cand.addUserInt("tkTrail_isMatchedToLowPtEle",isMatchedToLowPtSC_[itrail]);
    /// variables taken from the electron object after checking that is the same ele as in the conversion
    cand.addUserFloat("eleLead_pt", elePt[ilead]);
    cand.addUserFloat("eleTrail_pt", elePt[itrail]);
    cand.addUserFloat("eleLead_ch", eleCh[ilead]);
    cand.addUserFloat("eleTrail_ch", eleCh[itrail]);
    cand.addUserFloat("eleLead_eta", eleEta[ilead]);
    cand.addUserFloat("eleTrail_eta", eleEta[itrail]);
    cand.addUserFloat("eleLead_phi", elePhi[ilead]);
    cand.addUserFloat("eleTrail_phi", elePhi[itrail]);
    //
    cand.addUserFloat("eleLead_scEnergy", scEnergy[ilead]);
    cand.addUserFloat("eleTrail_scEnergy", scEnergy[itrail]);
    cand.addUserFloat("eleLead_scEta", scEta[ilead]);
    cand.addUserFloat("eleTrail_scEta", scEta[itrail]);
    cand.addUserFloat("eleLead_scEtaWidth", scEtaWidth[ilead]);
    cand.addUserFloat("eleTrail_scEtaWidth", scEtaWidth[itrail]);
    cand.addUserFloat("eleLead_scPhiWidth", scPhiWidth[ilead]);
    cand.addUserFloat("eleTrail_scPhiWidth", scPhiWidth[itrail]);
    cand.addUserFloat("eleLead_cluSize",cluSize[ilead]);
    cand.addUserFloat("eleTrail_cluSize",cluSize[itrail]);
    cand.addUserFloat("eleLead_full5x5_HoverE",full5x5_HoverE[ilead]);
    cand.addUserFloat("eleTrail_full5x5_HoverE",full5x5_HoverE[itrail]);
    cand.addUserFloat("eleLead_full5x5_r9",full5x5_r9[ilead]);
    cand.addUserFloat("eleTrail_full5x5_r9",full5x5_r9[itrail]);
    //

    cand.addUserFloat("eleLead_seed_dEta",ele_seed_dEta[ilead]);
    cand.addUserFloat("eleTrail_seed_dEta",ele_seed_dEta[itrail]);
    cand.addUserFloat("eleLead_eclu_EOverP",ele_eclu_EoverP[ilead]);
    cand.addUserFloat("eleTrail_eclu_EOverP",ele_eclu_EoverP[itrail]);
    cand.addUserFloat("eleLead_scEOverP",eleV[ilead].eSuperClusterOverP());
    cand.addUserFloat("eleTrail_scEOverP",eleV[itrail].eSuperClusterOverP());

    cand.addUserFloat("eleLead_sc_dEta",ele_sc_dEta[ilead]);
    cand.addUserFloat("eleTrail_sc_dEta",ele_sc_dEta[itrail]);
    cand.addUserFloat("eleLead_sc_dPhi",ele_sc_dPhi[ilead]);
    cand.addUserFloat("eleTrail_sc_dPhi",ele_sc_dPhi[itrail]);
    //
    cand.addUserFloat("eleLead_fBrem",ele_fBrem[ilead]);
    cand.addUserFloat("eleTrail_fBrem",ele_fBrem[itrail]);
    cand.addUserFloat("eleLead_shFracHits",ele_shFrInHits[ilead]);
    cand.addUserFloat("eleTrail_shFracHits",ele_shFrInHits[itrail]);

    cand.addUserFloat("eleLead_pMode",gsfPmode[ilead]);
    cand.addUserFloat("eleTrail_pMode",gsfPmode[itrail]);
    cand.addUserFloat("eleLead_chi2",gsfChi2[ilead]);
    cand.addUserFloat("eleTrail_chi2",gsfChi2[itrail]);
    cand.addUserFloat("eleLead_nHits",gsfNhits[ilead]);
    cand.addUserFloat("eleTrail_nHits",gsfNhits[itrail]);
    cand.addUserFloat("eleLead_dR",gsfDR[ilead]);
    cand.addUserFloat("eleTrail_dR",gsfDR[itrail]);

    cand.addUserFloat("eleLead_unbiasedSeedBDT",unbiased_seedBDT[ilead]);
    cand.addUserFloat("eleTrail_unbiasedSeedBDT",unbiased_seedBDT[itrail]);
    cand.addUserFloat("eleLead_mvaId",mvaIdV[ilead]);
    cand.addUserFloat("eleTrail_mvaId",mvaIdV[itrail]);
 


    if (  nMatchedTracks ==2  ) {
      
      /*   
      isTwoSuperClusterConversion_=false;
      dEtaTracksAtECAL_whenTwoSC_=-99.;
      dPhiTracksAtECAL_whenTwoSC_=-99.;
      dEtaSC_whenTwoSC_ = -99.;
      dPhiSC_whenTwoSC_ = -99.;
      
      
      
      if ( eleV[ilead].superCluster().id() == eleV[itrail].superCluster().id() && eleV[ilead].superCluster().key() == eleV[itrail].superCluster().key() ) {
	std::cout << " the two electrons have the same supercluster " << std::endl;
      } else {
	std::cout << " the two electrons DO NOT have the same supercluster " << std::endl;
	isTwoSuperClusterConversion_=true;
	dEtaTracksAtECAL_whenTwoSC_ =   eleV[ilead].trackPositionAtCalo().eta() - eleV[itrail].trackPositionAtCalo().eta();       
	dPhiTracksAtECAL_whenTwoSC_ =  reco::deltaPhi( eleV[ilead].trackPositionAtCalo().phi(),eleV[itrail].trackPositionAtCalo().phi());  
	dEtaSC_whenTwoSC_ = eleV[ilead].superClusterPosition().eta() -    eleV[itrail].superClusterPosition().eta();   
	dPhiSC_whenTwoSC_ =  reco::deltaPhi( eleV[ilead].superClusterPosition().phi(),eleV[itrail].superClusterPosition().phi());  
	
      }
      
      */
    }  

    outFilteredConversionCollection.push_back(cand);    
    
   } // loop over reco conversions
  
  
  nAcceptedConvTot=nAcceptedConvTot +  outFilteredConversionCollection.size();
  //std::cout << " Accepted conv per event " << nAcceptedConv << std::endl;
  //std::cout << " All Accepted conv  " << nAcceptedConvTot << std::endl;
  
  


  outFilteredConversionCollection_p ->assign(outFilteredConversionCollection.begin(), outFilteredConversionCollection.end());
  e.put(std::move(outFilteredConversionCollection_p), filteredConvertedPhotonCollection_);
}







float ConversionSelector::mee (float ipx1, float ipy1, float ipz1, float ipx2, float ipy2, float ipz2){
                const float mElec= 0.000511;
                float invMass=-99.;
                const float px = ipx1+ipx2;
                const float py = ipy1+ipy2;
                const float pz = ipz1+ipz2;
                const float mom1 = ipx1*ipx1+ipy1*ipy1+ipz1*ipz1;
                const float mom2 = ipx2*ipx2+ipy2*ipy2+ipz2*ipz2;
                const float e = sqrt( mom1+ mElec*mElec ) + sqrt( mom2 + mElec*mElec );
                invMass= ( e*e - px*px -py*py - pz*pz);
                if ( invMass>0.) invMass = sqrt(invMass);
                else 
                                invMass = -1;
                return invMass;
}

DEFINE_FWK_MODULE(ConversionSelector);
