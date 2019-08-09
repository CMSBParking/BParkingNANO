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
#include <limits>
#include <algorithm>

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
    drTrg_cleaning_{cfg.getParameter<double>("drForCleaning_wrtTrgMuon")},
    dzTrg_cleaning_{cfg.getParameter<double>("dzForCleaning_wrtTrgMuon")},
    dr_cleaning_{cfg.getParameter<double>("drForCleaning")},
    dz_cleaning_{cfg.getParameter<double>("dzForCleaning")},
    ptMin_{cfg.getParameter<double>("ptMin")},
    etaMax_{cfg.getParameter<double>("etaMax")},
    bdtMin_{cfg.getParameter<double>("bdtMin")},
    use_gsf_mode_for_p4_{cfg.getParameter<bool>("useGsfModeForP4")} 
    {
       produces<pat::ElectronCollection>("SelectedElectrons");
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
  const double drTrg_cleaning_;
  const double dzTrg_cleaning_;
  const double dr_cleaning_;
  const double dz_cleaning_;
  const double ptMin_; //pt min cut
  const double etaMax_; //eta max cut
  const double bdtMin_; //bdt min cut
  const bool use_gsf_mode_for_p4_;
};

void ElectronMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

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

  // output
  std::unique_ptr<pat::ElectronCollection> ele_out(new pat::ElectronCollection);

  
  // -> changing order of loops ert Arabella's fix this without need for more vectors  
 for(auto ele : *pf) {
   //cuts
   if (ele.pt()<ptMin_) continue;
   if (fabs(ele.eta())>etaMax_) continue;
   // apply conversion veto unless we want conversions
   if (!ele.passConversionVeto()) continue;

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
   ele_out->emplace_back(ele);
 }


 size_t iele=-1;
 /// add and clean low pT e
  for(auto ele : *lowpt) {
    iele++;
    //take modes
   if(use_gsf_mode_for_p4_) {
     reco::Candidate::PolarLorentzVector p4( ele.gsfTrack()->ptMode(),
                                             ele.gsfTrack()->etaMode(),
                                             ele.gsfTrack()->phiMode(),
                                             ele.mass()    );
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
   for(const auto& pfele : *pf) {
      clean_out |= (
	           fabs(pfele.vz() - ele.vz()) < dz_cleaning_ &&
                   reco::deltaR(ele, pfele) < dr_cleaning_   );
   }
   if(clean_out) continue;
   edm::Ref<pat::ElectronCollection> ref(lowpt,iele);
   float mva_id = float((*mvaId)[ref]);
   ele.addUserInt("isPF", 0);
   ele.addUserInt("isLowPt", 1);
   ele.addUserFloat("chargeMode", ele.gsfTrack()->chargeMode());
   ele.addUserFloat("ptBiased", ptbiased_seedBDT);
   ele.addUserFloat("unBiased", unbiased_seedBDT);
   ele.addUserFloat("mvaId", mva_id);
   ele_out->emplace_back(ele);
 }

//is this nescaisery ? because it is additional loop
   /* std::sort( out->begin(), out->end(), [] (pat::Electron e1, pat::Electron e2) -> bool {return e1.pt() > e2.pt();}
	      );
  }*/

  //adding label to be consistent with the muon and track naming
  evt.put(std::move(ele_out),"SelectedElectrons");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronMerger);
