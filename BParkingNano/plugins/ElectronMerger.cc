// Merges the PF and LowPT collections, sets the isPF and isLowPt 
// UserInt's accordingly

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <limits>
#include <algorithm>
#include "helper.h"

class ElectronMerger : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)


public:

  explicit ElectronMerger(const edm::ParameterSet &cfg):
    triggerMuons_{ consumes<pat::MuonCollection>( cfg.getParameter<edm::InputTag>("trgMuon") )},
    lowpt_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("lowptSrc") )},
    pf_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("pfSrc") )},
    ptBiased_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("ptbiasedSeeding") )},
    unBiased_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("unbiasedSeeding") )},
    mvaId_src_{ consumes<edm::ValueMap<float>>( cfg.getParameter<edm::InputTag>("mvaId") )},
    vertexSrc_{ consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>("vertexCollection") )},
    drTrg_cleaning_{cfg.getParameter<double>("drForCleaning_wrtTrgMuon")},
    dzTrg_cleaning_{cfg.getParameter<double>("dzForCleaning_wrtTrgMuon")},
    dr_cleaning_{cfg.getParameter<double>("drForCleaning")},
    dz_cleaning_{cfg.getParameter<double>("dzForCleaning")},
    flagAndclean_{cfg.getParameter<bool>("flagAndclean")},
    pf_ptMin_{cfg.getParameter<double>("pf_ptMin")},
    ptMin_{cfg.getParameter<double>("ptMin")},
    etaMax_{cfg.getParameter<double>("etaMax")},
    bdtMin_{cfg.getParameter<double>("bdtMin")},
    use_gsf_mode_for_p4_{cfg.getParameter<bool>("useGsfModeForP4")},
    sortOutputCollections_{cfg.getParameter<bool>("sortOutputCollections")} 
    {
       produces<pat::ElectronCollection>("SelectedElectrons");
       produces<TransientTrackCollection>("SelectedTransientElectrons");  
    }

  ~ElectronMerger() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::MuonCollection> triggerMuons_;
  const edm::EDGetTokenT<pat::ElectronCollection> lowpt_src_;
  const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ptBiased_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> unBiased_src_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mvaId_src_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  const double drTrg_cleaning_;
  const double dzTrg_cleaning_;
  const double dr_cleaning_;
  const double dz_cleaning_;
  const bool flagAndclean_;
  const double pf_ptMin_;
  const double ptMin_; //pt min cut
  const double etaMax_; //eta max cut
  const double bdtMin_; //bdt min cut
  const bool use_gsf_mode_for_p4_;
  const bool sortOutputCollections_;
};

void ElectronMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & iSetup) const {

  //input
  edm::Handle<pat::MuonCollection> trgMuon;
  evt.getByToken(triggerMuons_, trgMuon);
  edm::Handle<pat::ElectronCollection> lowpt;
  evt.getByToken(lowpt_src_, lowpt);
  edm::Handle<pat::ElectronCollection> pf;
  evt.getByToken(pf_src_, pf);
  edm::Handle<edm::ValueMap<float> > ptBiased;
  evt.getByToken(ptBiased_src_, ptBiased);
  edm::Handle<edm::ValueMap<float> > unBiased;
  evt.getByToken(unBiased_src_, unBiased);
  edm::Handle<edm::ValueMap<float> > mvaId;  
  evt.getByToken(mvaId_src_, mvaId);
  // 
  edm::ESHandle<TransientTrackBuilder> theB ;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  //
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  // output
  std::unique_ptr<pat::ElectronCollection>  ele_out      (new pat::ElectronCollection );
  std::unique_ptr<TransientTrackCollection> trans_ele_out(new TransientTrackCollection);
  std::vector<std::pair<float, float>> pfEtaPhi;
  std::vector<float> pfVz;

  // -> changing order of loops ert Arabella's fix this without need for more vectors  
  for(auto ele : *pf) {
   //cuts
   if (ele.pt()<ptMin_ || ele.pt() < pf_ptMin_) continue;
   if (fabs(ele.eta())>etaMax_) continue;
   // apply conversion veto unless we want conversions
   if (!ele.passConversionVeto()) continue;

   // Fix the mass to the proper one
   reco::Candidate::PolarLorentzVector p4( 
     ele.pt(),
     ele.eta(),
     ele.phi(),
     ELECTRON_MASS
     );
   ele.setP4(p4);     

   // skip electrons inside tag's jet or from different PV
   bool skipEle=true;
   for(const auto & trg : *trgMuon) {
     if(reco::deltaR(ele, trg) < drTrg_cleaning_ && drTrg_cleaning_ > 0)
        continue;
     if(fabs(ele.vz() - trg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0)
        continue;
     skipEle=false;
     break; // one trg muon to pass is enough :)
   }
   // we skip evts without trg muon
   if (skipEle) continue;

   // for PF e we set BDT outputs to much higher number than the max
   ele.addUserInt("isPF", 1);
   ele.addUserInt("isLowPt", 0);
   ele.addUserFloat("ptBiased", 20.);
   ele.addUserFloat("unBiased", 20.);
   ele.addUserFloat("mvaId", 20);
   ele.addUserFloat("chargeMode", ele.charge());
   ele.addUserInt("isPFoverlap", 0);

   pfEtaPhi.push_back(std::pair<float, float>(ele.eta(), ele.phi()));
   pfVz.push_back(ele.vz());
   ele_out       -> emplace_back(ele);
  }

  unsigned int pfSelectedSize = pfEtaPhi.size();

  size_t iele=-1;
  /// add and clean low pT e
  for(auto ele : *lowpt) {
    iele++;
    //take modes
   if(use_gsf_mode_for_p4_) {
     reco::Candidate::PolarLorentzVector p4( ele.gsfTrack()->ptMode(),
                                             ele.gsfTrack()->etaMode(),
                                             ele.gsfTrack()->phiMode(),
                                             ELECTRON_MASS    );
     ele.setP4(p4);
   } 
   else {
     // Fix the mass to the proper one
     reco::Candidate::PolarLorentzVector p4( 
       ele.pt(),
       ele.eta(),
       ele.phi(),
       ELECTRON_MASS
       );
     ele.setP4(p4);     
   }

   //same cuts as in PF
   if (ele.pt()<ptMin_) continue;
   if (fabs(ele.eta())>etaMax_) continue;
   // apply conversion veto?
   if (!ele.passConversionVeto()) continue;

   //assigning BDT values
   const reco::GsfTrackRef gsfTrk = ele.gsfTrack();
   float unbiased_seedBDT = float((*unBiased)[gsfTrk]);
   float ptbiased_seedBDT = float((*ptBiased)[gsfTrk]);
   if ( unbiased_seedBDT <bdtMin_) continue; //extra cut for low pT e on BDT

   bool skipEle=true;
   for(const auto & trg : *trgMuon) {
     if(reco::deltaR(ele, trg) < drTrg_cleaning_ && drTrg_cleaning_ > 0)
        continue;
     if(fabs(ele.vz() - trg.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ > 0)
        continue;
     skipEle=false;
     break;  // one trg muon is enough 
   }
   // same here Do we need evts without trg muon? now we skip them
   if (skipEle) continue;   
   

   //pf cleaning    
   bool clean_out = false;
   for(unsigned int iEle=0; iEle<pfSelectedSize; ++iEle) {

      clean_out |= (
	           fabs(pfVz[iEle] - ele.vz()) < dz_cleaning_ &&
                   reco::deltaR(ele.eta(), ele.phi(), pfEtaPhi[iEle].first, pfEtaPhi[iEle].second) < dr_cleaning_   );

   }
   if(clean_out && flagAndclean_) continue;
   else if(clean_out) ele.addUserInt("isPFoverlap", 1);
   else ele.addUserInt("isPFoverlap", 0);

   edm::Ref<pat::ElectronCollection> ref(lowpt,iele);
   float mva_id = float((*mvaId)[ref]);
   ele.addUserInt("isPF", 0);
   ele.addUserInt("isLowPt", 1);
   ele.addUserFloat("chargeMode", ele.gsfTrack()->chargeMode());
   ele.addUserFloat("ptBiased", ptbiased_seedBDT);
   ele.addUserFloat("unBiased", unbiased_seedBDT);
   ele.addUserFloat("mvaId", mva_id);

   ele_out       -> emplace_back(ele);
  }

  if(sortOutputCollections_){

    //sorting increases sligtly the time but improves the code efficiency in the Bcandidate builder
    //easier identification of leading and subleading with smarter loop
    std::sort( ele_out->begin(), ele_out->end(), [] (pat::Electron e1, pat::Electron e2) -> bool {return e1.pt() > e2.pt();}
             );
  }

  // build transient track collection
  for(auto &ele : *ele_out){
    const reco::TransientTrack eleTT =(*theB).buildfromGSF( ele.gsfTrack() );
    trans_ele_out -> emplace_back(eleTT);

    if(ele.userInt("isPF")) continue;
    //compute IP for electrons: need transient track
    //from PhysicsTools/PatAlgos/plugins/LeptonUpdater.cc
    const reco::GsfTrackRef gsfTrk = ele.gsfTrack();
    // PVDZ
    ele.setDB(gsfTrk->dz(PV.position()), std::hypot(gsfTrk->dzError(), PV.zError()), pat::Electron::PVDZ);

    //PV2D
    std::pair<bool, Measurement1D> result = IPTools::signedTransverseImpactParameter(eleTT, GlobalVector(gsfTrk->px(), gsfTrk->py(), gsfTrk->pz()), PV);
    double d0_corr = result.second.value();
    double d0_err = PV.isValid() ? result.second.error() : -1.0;
    ele.setDB(d0_corr, d0_err, pat::Electron::PV2D);

    // PV3D
    result = IPTools::signedImpactParameter3D(eleTT, GlobalVector(gsfTrk->px(), gsfTrk->py(), gsfTrk->pz()), PV);
    d0_corr = result.second.value();
    d0_err = PV.isValid() ? result.second.error() : -1.0;
    ele.setDB(d0_corr, d0_err, pat::Electron::PV3D);
  }
   
  //adding label to be consistent with the muon and track naming
  evt.put(std::move(ele_out),      "SelectedElectrons");
  evt.put(std::move(trans_ele_out),"SelectedTransientElectrons");
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronMerger);
