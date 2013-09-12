// -*- C++ -*-
//
// Package:    HltGsfEleAnalyser
// Class:      HltGsfEleAnalyser
// 
/**\class HltGsfEleAnalyser HltGsfEleAnalyser.cc UserCode/HltGsfTrkAnalyser/src/HltGsfEleAnalyser.cc

 Description: Analysing GSF electrons wrt. HLT paths.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Reis
//         Created:  Wed Jul 24 14:40:33 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TH2F.h"
//
// class declaration
//

class HltGsfEleAnalyser : public edm::EDAnalyzer {
   public:
      explicit HltGsfEleAnalyser(const edm::ParameterSet&);
      ~HltGsfEleAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      bool PassHeep(reco::GsfElectron const&, double &);
      // ----------member data ---------------------------
      edm::InputTag electronCollTag_;
      edm::InputTag trgResultsTag_;
      std::vector<edm::ParameterSet> trgVPSet_;

      float deltaEtCut = 15.;
      float etCutStepSize = 0.5;

      HLTConfigProvider hltConfig;
      std::vector<unsigned int> trgIndices;
      std::vector<std::string> trgNames;
      std::vector<bool> trgInvs;
      std::vector<unsigned int> trgMinEles;
      std::vector<std::vector<double> > trgMinEtss;

      // histograms
      TH1F* h_nHeep;
      TH1F* h_total;
      TH1F* h_1PassHeep_total;
      TH1F* h_2PassHeep_total;
      TH1F* h_nPassHeep_total;
      TH1F* h_nPassHeep_trgd;
      TH1F* h_nPassHeep_notTrgd;
      TH1F* h_nPassHeep_nPassEt_total;
      TH1F* h_nPassHeep_nPassEt_trgd;
      TH1F* h_nPassHeep_nPassEt_notTrgd;
      TH1F* h_trgd;
      TH1F* h_notTrgd;
      std::vector<TH1F*> v_h_nPassHeep_nMinus1PassEt_vsEt_total;
      std::vector<TH1F*> v_h_nPassHeep_nMinus1PassEt_vsEt_trgd;

      // ratio histograms
      TH1F* hRatio_nPassHeep_trgdVsTotal;
      TH1F* hRatio_nPassHeep_notTrgdVsTotal;
      TH1F* hRatio_nPassHeep_nPassEt_trgdVsTotal;
      TH1F* hRatio_nPassHeep_nPassEt_notTrgdVsTotal;
      std::vector<TH1F*> v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal;

      TH2F* h2Ratio_nPassHeep_trgd;
      TH2F* h2Ratio_nPassHeep_notTrgd;
      TH2F* h2Ratio_nPassHeep_nPassEt_trgd;
      TH2F* h2Ratio_nPassHeep_nPassEt_notTrgd;
      TH2F* h2Ratio_trgd;
      TH2F* h2Ratio_notTrgd;
      //TODO: make 2D histos wrt. nPassHEEP_total 
};


// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HltGsfEleAnalyser::HltGsfEleAnalyser(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   electronCollTag_ = iConfig.getParameter<edm::InputTag>("electronCollTag");
   trgResultsTag_ = iConfig.getParameter<edm::InputTag>("trgResultsTag");
   trgVPSet_ = iConfig.getUntrackedParameter<std::vector<edm::ParameterSet> >("triggers");

   for (std::vector<edm::ParameterSet>::const_iterator trgIt = trgVPSet_.begin(); trgIt < trgVPSet_.end(); ++trgIt) {
      trgNames.push_back(trgIt->getUntrackedParameter<std::string>("triggerName", "HLTriggerFinalPath")); 
      trgInvs.push_back(trgIt->getUntrackedParameter<bool>("invertBit", 0));
      trgMinEles.push_back(trgIt->getUntrackedParameter<unsigned int>("minEle", 1));
      trgMinEtss.push_back(trgIt->getUntrackedParameter<std::vector<double> >("minEts", std::vector<double> (1, 0.)));
   }

   TH1::SetDefaultSumw2(kTRUE);

   unsigned int trgsSize = trgVPSet_.size();
   edm::Service<TFileService> fs;
   h_total = fs->make<TH1F>("h_total", "Total # events", trgsSize, 0., trgsSize);
   h_1PassHeep_total = fs->make<TH1F>("h_1PassHeep_total", "Total events with one HEEP electron", trgsSize, 0., trgsSize);
   h_2PassHeep_total = fs->make<TH1F>("h_2PassHeep_total", "Total events with two HEEP electrons", trgsSize, 0., trgsSize);
   h_nPassHeep_total = fs->make<TH1F>("h_nPassHeep_total", "Total events with n HEEP electrons (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_trgd = fs->make<TH1F>("h_nPassHeep_trgd", "Triggered events with n HEEP electrons (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_notTrgd = fs->make<TH1F>("h_nPassHeep_notTrgd", "Not triggered events with n HEEP electrons (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_total = fs->make<TH1F>("h_nPassHeep_nPassEt_total", "Total events with n HEEP electrons passing Et cut (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_trgd = fs->make<TH1F>("h_nPassHeep_nPassEt_trgd", "Triggered events with n HEEP electrons passing Et cut (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_notTrgd = fs->make<TH1F>("h_nPassHeep_nPassEt_notTrgd", "Not triggered events with n HEEP electrons passing Et cut (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nHeep = fs->make<TH1F>("h_nHeep", "# of HEEP electrons", 10, 0., 10.);
   h_trgd = fs->make<TH1F>("h_trgd", "# of triggered events", trgsSize, 0., trgsSize);
   h_notTrgd = fs->make<TH1F>("h_notTrgd", "# of not triggered events", trgsSize, 0., trgsSize);

   // ratio histos
   hRatio_nPassHeep_trgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_trgdVsTotal", "# triggered / # total n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_notTrgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_notTrgdVsTotal", "# not triggered / # total n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_trgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_trgdVsTotal", "# triggered / # total n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_notTrgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_notTrgdVsTotal", "# not triggered / # total n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);

   h2Ratio_nPassHeep_trgd = fs->make<TH2F>("h2Ratio_nPassHeep_trgd", "Inter trigger ratio of triggered n-HEEP events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_notTrgd = fs->make<TH2F>("h2Ratio_nPassHeep_notTrgd", "Inter trigger ratio of not triggered n-HEEP events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_nPassEt_trgd = fs->make<TH2F>("h2Ratio_nPassHeep_nPassEt_trgd", "Inter trigger ratio of triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_nPassEt_notTrgd = fs->make<TH2F>("h2Ratio_nPassHeep_nPassEt_notTrgd", "Inter trigger ratio of not triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_trgd = fs->make<TH2F>("h2Ratio_trgd", "Inter trigger ratio of triggered events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_notTrgd = fs->make<TH2F>("h2Ratio_notTrgd", "Inter trigger ratio of not triggered events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);

   for (unsigned int i = 0; i < trgNames.size(); ++i) {
      float lastTrgMinEt = trgMinEtss.at(i).back(); 
      std::string nameString = "h_nPassHeep_nMinus1PassEt_vsEt_total_" + trgNames[i];
      v_h_nPassHeep_nMinus1PassEt_vsEt_total.push_back(fs->make<TH1F>(nameString.data(), "# of events vs. Et cut", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "h_nPassHeep_nMinus1PassEt_vsEt_trgd_" + trgNames[i];
      v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.push_back(fs->make<TH1F>(nameString.data(), "# of triggered events vs. Et cut", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal_" + trgNames[i];
      v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal.push_back(fs->make<TH1F>(nameString.data(), "# triggered / # total n-HEEP events vs. Et cuts", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));

      // set bin labels
      const char* trgName = trgNames[i].data();
      h_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_1PassHeep_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_2PassHeep_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);

      // ratio histos
      hRatio_nPassHeep_trgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_notTrgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_trgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_notTrgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);

      h2Ratio_nPassHeep_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_trgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_notTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_trgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_notTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_trgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_notTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
   }
}


HltGsfEleAnalyser::~HltGsfEleAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HltGsfEleAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // get GSF electron collection
   edm::InputTag inTagGsfEle(electronCollTag_);
   edm::Handle<reco::GsfElectronCollection> gsfEleHandle;
   if (!iEvent.getByLabel(inTagGsfEle, gsfEleHandle)) return;
   const reco::GsfElectronCollection* gsfElectrons = gsfEleHandle.product();

   // get the rho variable
   double rho = 0;
   edm::Handle<double> rhoHandle;
   if (iEvent.getByLabel(edm::InputTag("kt6PFJets:rho"), rhoHandle)) rho = *rhoHandle.product();

   unsigned int nPassHeep = 0;
   std::vector<unsigned int> heepEleIndices;
   for (reco::GsfElectronCollection::const_iterator eleIt = gsfElectrons->begin(); eleIt < gsfElectrons->end(); ++eleIt) {
      if (PassHeep(*eleIt, rho)) {
         ++nPassHeep;
         heepEleIndices.push_back(eleIt - gsfElectrons->begin());
      }
   }
   h_nHeep->Fill(nPassHeep);

   // trigger requirements
   edm::TriggerResultsByName trgRes = iEvent.triggerResultsByName(trgResultsTag_.process());
   for (unsigned int trgIndex = 0; trgIndex < trgNames.size(); ++trgIndex) {
      bool triggered = trgRes.accept(trgIndices[trgIndex]) ^ trgInvs[trgIndex];

      h_total->Fill(trgIndex);
      if (nPassHeep > 0) h_1PassHeep_total->Fill(trgIndex);
      if (nPassHeep > 1) h_2PassHeep_total->Fill(trgIndex);

      if (nPassHeep >= trgMinEles[trgIndex]) {
         h_nPassHeep_total->Fill(trgIndex);
         if (triggered) h_nPassHeep_trgd->Fill(trgIndex);
         else h_nPassHeep_notTrgd->Fill(trgIndex);

         // find heep electrons passing also Et cuts from cfg file
         std::vector<unsigned int> heepEleSearchVector(heepEleIndices);
         std::vector<unsigned int> heepEleSearchVectorForked(heepEleIndices);
         for (unsigned int etCutIndex = 0; etCutIndex < trgMinEtss.at(trgIndex).size(); ++etCutIndex) {
            // find smallest passing electron
            float smallestPassEt = 9999999.;
            unsigned int smallestPassIndex = 9999999;
            for (unsigned int k = 0; k < heepEleSearchVector.size(); ++k) {
               float et = gsfElectrons->at(k).caloEnergy() * sin(gsfElectrons->at(k).p4().theta());
               if (et > trgMinEtss.at(trgIndex).at(etCutIndex) && et < smallestPassEt) {
                  smallestPassEt = et;
                  smallestPassIndex = k;
               }
            }
            // found HEEP electron passing Et cut -> erase electron index from search vector
            if (smallestPassIndex < heepEleSearchVector.size()) heepEleSearchVector.erase(heepEleSearchVector.begin() + smallestPassIndex);

            // fork search vector for et scan
            if (etCutIndex == trgMinEtss.at(trgIndex).size() - 2) heepEleSearchVectorForked = heepEleSearchVector;
         }

         // fill histos in case we found enough HEEP electrons passing Et cuts
         if (heepEleIndices.size() - heepEleSearchVector.size() == trgMinEtss.at(trgIndex).size()) {
            h_nPassHeep_nPassEt_total->Fill(trgIndex);
            if (triggered) h_nPassHeep_nPassEt_trgd->Fill(trgIndex);
            else h_nPassHeep_nPassEt_notTrgd->Fill(trgIndex);
         }

         // fill histograms vs. et of electron
         float lastTrgMinEt = trgMinEtss.at(trgIndex).back(); 
         for (float etCut = lastTrgMinEt-deltaEtCut; etCut < lastTrgMinEt+deltaEtCut; etCut += etCutStepSize) {
            // check the remaining HEEP electrons if one is in the Et range
            int subtractValue = 0;
            for (unsigned int k = 0; k < heepEleSearchVectorForked.size(); ++k) {
               float et = gsfElectrons->at(k).caloEnergy() * sin(gsfElectrons->at(k).p4().theta());
               if (et > etCut && et < etCut + etCutStepSize) {
                  subtractValue = 1;
                  break;
               }
            }

            // fill histos in case we found enough HEEP electrons passing Et cuts
            if (heepEleIndices.size() - (heepEleSearchVectorForked.size() - subtractValue) == trgMinEtss.at(trgIndex).size()) {
               v_h_nPassHeep_nMinus1PassEt_vsEt_total.at(trgIndex)->Fill(etCut);
               if (triggered) v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.at(trgIndex)->Fill(etCut);
            }
         }
      }

      if (triggered) h_trgd->Fill(trgIndex);
      else h_notTrgd->Fill(trgIndex);

   }
}

bool
HltGsfEleAnalyser::PassHeep(reco::GsfElectron const& ele, double &rho) {
   bool pass = 0;

   float et = ele.caloEnergy() * sin(ele.p4().theta());
   float absScEta = fabs(ele.superCluster()->eta());
   float ecalHcal1Iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
   const float ecalHcal1EffArea_bar = 0.28;
   const float ecalHcal1EffArea_end = 0.28;

   bool scBar = absScEta < 1.442;
   bool scEnd = absScEta > 1.56 && absScEta < 2.5;

   bool ecalSeeded = 1;
   bool missingHits_bar = 1;
   bool dxyFirstPv_bar = 1;
   bool dEta_bar = 1;
   bool dPhi_bar = 1;
   bool hOverE_bar = 1;
   bool eXx5overe5x5_bar = 1;
   bool ecalHcal1Iso_bar = 1;
   bool trackIso_bar = 1;

   bool missingHits_end = 1;
   bool dxyFirstPv_end = 1;
   bool dEta_end = 1;
   bool dPhi_end = 1;
   bool hOverE_end = 1;
   bool sigmaIetaIeta_end = 1;
   bool ecalHcal1Iso_end = 1;
   bool trackIso_end = 1;

   if (scBar) {
      ecalSeeded = ele.ecalDrivenSeed();
      missingHits_bar = ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < 2;
      dxyFirstPv_bar = fabs(ele.gsfTrack()->dxy(ele.vertex())) < 0.02;
      dEta_bar = fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.005;
      dPhi_bar = fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.06;
      hOverE_bar = ele.hadronicOverEm() < 0.05;
      eXx5overe5x5_bar = (ele.scE2x5Max() / ele.scE5x5() > 0.94) || (ele.scE1x5() / ele.scE5x5() > 0.83);
      ecalHcal1Iso_bar = ecalHcal1Iso < (2. + 0.03*et + rho*ecalHcal1EffArea_bar);
      trackIso_bar = ele.dr03TkSumPt() < 5.;

      pass = ecalSeeded && missingHits_bar && dxyFirstPv_bar && dEta_bar && dPhi_bar && hOverE_bar && eXx5overe5x5_bar && ecalHcal1Iso_bar && trackIso_bar;
   }
   else if (scEnd) {
      ecalSeeded = ele.ecalDrivenSeed();
      missingHits_end = ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() < 2;
      dxyFirstPv_end = fabs(ele.gsfTrack()->dxy(ele.vertex())) < 0.05;
      dEta_end = fabs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.007;
      dPhi_end = fabs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.06;
      hOverE_end = ele.hadronicOverEm() < 0.05;
      sigmaIetaIeta_end = ele.sigmaIetaIeta() < 0.03;
      if (et < 50.) ecalHcal1Iso_end = ecalHcal1Iso < 2.5 + rho*ecalHcal1EffArea_end;
      else ecalHcal1Iso_end = ecalHcal1Iso < (2.5 + 0.03*(et-50.) + rho*ecalHcal1EffArea_end);
      trackIso_end = ele.dr03TkSumPt() < 5.;

      pass = ecalSeeded && missingHits_end && dxyFirstPv_end && dEta_end && dPhi_end && hOverE_end && sigmaIetaIeta_end && ecalHcal1Iso_end && trackIso_end;
   }

   return pass;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HltGsfEleAnalyser::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HltGsfEleAnalyser::endJob() 
{
   // make ratios
   hRatio_nPassHeep_trgdVsTotal->Divide(h_nPassHeep_trgd, h_nPassHeep_total);
   hRatio_nPassHeep_notTrgdVsTotal->Divide(h_nPassHeep_notTrgd, h_nPassHeep_total);
   hRatio_nPassHeep_nPassEt_trgdVsTotal->Divide(h_nPassHeep_nPassEt_trgd, h_nPassHeep_nPassEt_total);
   hRatio_nPassHeep_nPassEt_notTrgdVsTotal->Divide(h_nPassHeep_nPassEt_notTrgd, h_nPassHeep_nPassEt_total);
   for (unsigned int i = 0; i < trgNames.size(); ++i) {
      v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal.at(i)->Divide(v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.at(i), v_h_nPassHeep_nMinus1PassEt_vsEt_total.at(i));
   }

   TH2F* h2Denom_nPassHeep_trgd = (TH2F*)h2Ratio_nPassHeep_trgd->Clone("h2Denom_nPassHeep_trgd");
   TH2F* h2Denom_nPassHeep_notTrgd = (TH2F*)h2Ratio_nPassHeep_notTrgd->Clone("h2Denom_nPassHeep_notTrgd");
   TH2F* h2Denom_nPassHeep_nPassEt_trgd = (TH2F*)h2Ratio_nPassHeep_nPassEt_trgd->Clone("h2Denom_nPassHeep_nPassEt_trgd");
   TH2F* h2Denom_nPassHeep_nPassEt_notTrgd = (TH2F*)h2Ratio_nPassHeep_nPassEt_notTrgd->Clone("h2Denom_nPassHeep_nPassEt_notTrgd");
   TH2F* h2Denom_trgd = (TH2F*)h2Ratio_trgd->Clone("h2Denom_trgd");
   TH2F* h2Denom_notTrgd = (TH2F*)h2Ratio_notTrgd->Clone("h2Denom_notTrgd");
   for (unsigned int ix = 1; ix <= trgNames.size(); ++ix) {
      for (unsigned int iy = 1; iy <= trgNames.size(); ++iy) {
         h2Ratio_nPassHeep_trgd->SetBinContent(ix, iy, h_nPassHeep_trgd->GetBinContent(ix));
         h2Ratio_nPassHeep_trgd->SetBinError(ix, iy, h_nPassHeep_trgd->GetBinError(ix));
         h2Denom_nPassHeep_trgd->SetBinContent(ix, iy, h_nPassHeep_trgd->GetBinContent(iy));
         h2Denom_nPassHeep_trgd->SetBinError(ix, iy, h_nPassHeep_trgd->GetBinError(iy));
         h2Ratio_nPassHeep_notTrgd->SetBinContent(ix, iy, h_nPassHeep_notTrgd->GetBinContent(ix));
         h2Ratio_nPassHeep_notTrgd->SetBinError(ix, iy, h_nPassHeep_notTrgd->GetBinError(ix));
         h2Denom_nPassHeep_notTrgd->SetBinContent(ix, iy, h_nPassHeep_notTrgd->GetBinContent(iy));
         h2Denom_nPassHeep_notTrgd->SetBinError(ix, iy, h_nPassHeep_notTrgd->GetBinError(iy));

         h2Ratio_nPassHeep_nPassEt_trgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinContent(ix));
         h2Ratio_nPassHeep_nPassEt_trgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinError(ix));
         h2Denom_nPassHeep_nPassEt_trgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinContent(iy));
         h2Denom_nPassHeep_nPassEt_trgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinError(iy));
         h2Ratio_nPassHeep_nPassEt_notTrgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinContent(ix));
         h2Ratio_nPassHeep_nPassEt_notTrgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinError(ix));
         h2Denom_nPassHeep_nPassEt_notTrgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinContent(iy));
         h2Denom_nPassHeep_nPassEt_notTrgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinError(iy));

         h2Ratio_trgd->SetBinContent(ix, iy, h_trgd->GetBinContent(ix));
         h2Ratio_trgd->SetBinError(ix, iy, h_trgd->GetBinError(ix));
         h2Denom_trgd->SetBinContent(ix, iy, h_trgd->GetBinContent(iy));
         h2Denom_trgd->SetBinError(ix, iy, h_trgd->GetBinError(iy));
         h2Ratio_notTrgd->SetBinContent(ix, iy, h_notTrgd->GetBinContent(ix));
         h2Ratio_notTrgd->SetBinError(ix, iy, h_notTrgd->GetBinError(ix));
         h2Denom_notTrgd->SetBinContent(ix, iy, h_notTrgd->GetBinContent(iy));
         h2Denom_notTrgd->SetBinError(ix, iy, h_notTrgd->GetBinError(iy));
      }
   }
   h2Ratio_nPassHeep_trgd->Divide(h2Denom_nPassHeep_trgd);
   h2Ratio_nPassHeep_notTrgd->Divide(h2Denom_nPassHeep_notTrgd);
   h2Ratio_nPassHeep_nPassEt_trgd->Divide(h2Denom_nPassHeep_nPassEt_trgd);
   h2Ratio_nPassHeep_nPassEt_notTrgd->Divide(h2Denom_nPassHeep_nPassEt_notTrgd);
   h2Ratio_trgd->Divide(h2Denom_trgd);
   h2Ratio_notTrgd->Divide(h2Denom_notTrgd);
}

// ------------ method called when starting to processes a run  ------------
void 
HltGsfEleAnalyser::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
   bool changed = 1;
   if (hltConfig.init(iRun, iSetup, trgResultsTag_.process(), changed)) {
      if (changed) {
         trgIndices.clear();
         for (unsigned int i = 0; i < trgNames.size(); ++i) {
            trgIndices.push_back(hltConfig.triggerIndex(trgNames[i]));
         }
      }
   } else {
      edm::LogError("beginRunHlt") << "HLTConfigProvider could not be initialised: '" << trgResultsTag_.label() << "', '" << trgResultsTag_.instance() << "', '" << trgResultsTag_.process() << "'";
   }
}

// ------------ method called when ending the processing of a run  ------------
void 
HltGsfEleAnalyser::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HltGsfEleAnalyser::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HltGsfEleAnalyser::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HltGsfEleAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HltGsfEleAnalyser);
