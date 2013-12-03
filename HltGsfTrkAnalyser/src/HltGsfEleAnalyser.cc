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
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
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
      bool MatchTrgFltr(const edm::Event&, reco::GsfElectron const&, const char*);
      // ----------member data ---------------------------
      edm::InputTag electronCollTag_;
      edm::InputTag trgResultsTag_;
      edm::InputTag trgEventTag_;
      std::vector<edm::ParameterSet> trgVPSet_;

      float deltaEtCut = 15.;
      float etCutStepSize = 0.5;

      HLTConfigProvider hltConfig;
      std::vector<unsigned int> trgIndices;
      std::vector<unsigned int> refTrgIndices;
      std::vector<unsigned int> tpTrgIndices;
      std::vector<std::string> trgNames;
      std::vector<std::string> refTrgNames;
      std::vector<bool> trgInvs;
      std::vector<unsigned int> trgMinEles;
      std::vector<std::vector<double> > trgMinEtss;
      std::vector<std::string> tpTrgNames;
      std::vector<std::string> tagFilterNames;
      std::vector<std::string> probeFilterNames;

      // histograms
      TH1F* h_gsfEle_et;
      TH1F* h_nHeep;
      TH1F* h_total;
      TH1F* h_refTrgd;
      TH1F* h_trgd;
      TH1F* h_notTrgd;
      TH1F* h_notTrgdRefTrgd;
      TH1F* h_1PassHeep_total;
      TH1F* h_1PassHeep_refTrgd;
      TH1F* h_2PassHeep_total;
      TH1F* h_2PassHeep_refTrgd;
      TH1F* h_nPassHeep_total;
      TH1F* h_nPassHeep_refTrgd;
      TH1F* h_nPassHeep_trgd;
      TH1F* h_nPassHeep_notTrgd;
      TH1F* h_nPassHeep_notTrgdRefTrgd;
      TH1F* h_nPassHeep_nPassEt_total;
      TH1F* h_nPassHeep_nPassEt_refTrgd;
      TH1F* h_nPassHeep_nPassEt_trgd;
      TH1F* h_nPassHeep_nPassEt_notTrgd;
      TH1F* h_nPassHeep_nPassEt_notTrgdRefTrgd;
      std::vector<TH1F*> v_h_nPassHeep_nMinus1PassEt_vsEt_total;
      std::vector<TH1F*> v_h_nPassHeep_nMinus1PassEt_vsEt_refTrgd;
      std::vector<TH1F*> v_h_nPassHeep_nMinus1PassEt_vsEt_trgd;

      TH1F* hTp_tags;
      TH1F* hTp_probes;
      TH1F* hTp_allProbesGsf;
      TH1F* hTp_passProbesGsf;
      TH1F* hTp_allProbesGsf_et25;
      TH1F* hTp_passProbesGsf_et25;
      TH1F* hTp_allProbesGsf_et33;
      TH1F* hTp_passProbesGsf_et33;
      TH1F* hTp_allProbesGsf_et35;
      TH1F* hTp_passProbesGsf_et35;
      TH1F* hTp_allProbesGsf_et80;
      TH1F* hTp_passProbesGsf_et80;
      TH1F* hTp_allProbesGsf_et100;
      TH1F* hTp_passProbesGsf_et100;
      TH1F* hTp_allProbesGsf_withEtCut;
      TH1F* hTp_passProbesGsf_withEtCut;
      std::vector<TH1F*> v_hTp_allProbesGsf_et;
      std::vector<TH1F*> v_hTp_passProbesGsf_et;
      TH1F* hTp_allProbesHeep;
      TH1F* hTp_passProbesHeep;
      TH1F* hTp_allProbesHeep_et25;
      TH1F* hTp_passProbesHeep_et25;
      TH1F* hTp_allProbesHeep_et33;
      TH1F* hTp_passProbesHeep_et33;
      TH1F* hTp_allProbesHeep_et35;
      TH1F* hTp_passProbesHeep_et35;
      TH1F* hTp_allProbesHeep_et80;
      TH1F* hTp_passProbesHeep_et80;
      TH1F* hTp_allProbesHeep_et100;
      TH1F* hTp_passProbesHeep_et100;
      TH1F* hTp_allProbesHeep_withEtCut;
      TH1F* hTp_passProbesHeep_withEtCut;
      std::vector<TH1F*> v_hTp_allProbesHeep_et;
      std::vector<TH1F*> v_hTp_passProbesHeep_et;

      // ratio histograms
      TH1F* hRatio_nPassHeep_trgdVsTotal;
      TH1F* hRatio_nPassHeep_notTrgdVsTotal;
      TH1F* hRatio_nPassHeep_nPassEt_trgdVsTotal;
      TH1F* hRatio_nPassHeep_nPassEt_notTrgdVsTotal;
      TH1F* hRatio_trgdVsRefTrgd;
      TH1F* hRatio_notTrgdVsRefTrgd;
      TH1F* hRatio_notTrgdRefTrgdVsRefTrgd;
      TH1F* hRatio_nPassHeep_trgdVsRefTrgd;
      TH1F* hRatio_nPassHeep_notTrgdVsRefTrgd;
      TH1F* hRatio_nPassHeep_notTrgdRefTrgdVsRefTrgd;
      TH1F* hRatio_nPassHeep_nPassEt_trgdVsRefTrgd;
      TH1F* hRatio_nPassHeep_nPassEt_notTrgdVsRefTrgd;
      TH1F* hRatio_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd;
      std::vector<TH1F*> v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal;
      std::vector<TH1F*> v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsRefTrgd;
      TH1F* hTpRatio_passProbesGsf_vs_allProbesGsf;
      TH1F* hTpRatio_passProbesGsf_et25_vs_allProbesGsf_et25;
      TH1F* hTpRatio_passProbesGsf_et33_vs_allProbesGsf_et33;
      TH1F* hTpRatio_passProbesGsf_et35_vs_allProbesGsf_et35;
      TH1F* hTpRatio_passProbesGsf_et80_vs_allProbesGsf_et80;
      TH1F* hTpRatio_passProbesGsf_et100_vs_allProbesGsf_et100;
      TH1F* hTpRatio_passProbesGsf_withEtCut_vs_allProbesGsf_withEtCut;
      std::vector<TH1F*> v_hTpRatio_passProbesGsf_et_vs_allProbesGsf_et;
      TH1F* hTpRatio_passProbesHeep_vs_allProbesHeep;
      TH1F* hTpRatio_passProbesHeep_et25_vs_allProbesHeep_et25;
      TH1F* hTpRatio_passProbesHeep_et33_vs_allProbesHeep_et33;
      TH1F* hTpRatio_passProbesHeep_et35_vs_allProbesHeep_et35;
      TH1F* hTpRatio_passProbesHeep_et80_vs_allProbesHeep_et80;
      TH1F* hTpRatio_passProbesHeep_et100_vs_allProbesHeep_et100;
      TH1F* hTpRatio_passProbesHeep_withEtCut_vs_allProbesHeep_withEtCut;
      std::vector<TH1F*> v_hTpRatio_passProbesHeep_et_vs_allProbesHeep_et;

      TGraphAsymmErrors* gEff_trgdVsRefTrgd;
      TGraphAsymmErrors* gEff_notTrgdRefTrgdVsRefTrgd;
      TGraphAsymmErrors* gEff_nPassHeep_trgdVsRefTrgd;
      TGraphAsymmErrors* gEff_nPassHeep_notTrgdRefTrgdVsRefTrgd;
      TGraphAsymmErrors* gEff_nPassHeep_nPassEt_trgdVsRefTrgd;
      TGraphAsymmErrors* gEff_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd;

      TH2F* h2Ratio_trgd;
      TH2F* h2Ratio_notTrgd;
      TH2F* h2Ratio_notTrgdRefTrgd;
      TH2F* h2Ratio_nPassHeep_trgd;
      TH2F* h2Ratio_nPassHeep_notTrgd;
      TH2F* h2Ratio_nPassHeep_notTrgdRefTrgd;
      TH2F* h2Ratio_nPassHeep_nPassEt_trgd;
      TH2F* h2Ratio_nPassHeep_nPassEt_notTrgd;
      TH2F* h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd;
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
   trgEventTag_ = iConfig.getParameter<edm::InputTag>("trgEventTag");
   trgVPSet_ = iConfig.getUntrackedParameter<std::vector<edm::ParameterSet> >("triggers");

   for (std::vector<edm::ParameterSet>::const_iterator trgIt = trgVPSet_.begin(); trgIt < trgVPSet_.end(); ++trgIt) {
      trgNames.push_back(trgIt->getUntrackedParameter<std::string>("triggerName", "HLTriggerFinalPath")); 
      refTrgNames.push_back(trgIt->getUntrackedParameter<std::string>("refTriggerName", "HLTriggerFinalPath")); 
      trgInvs.push_back(trgIt->getUntrackedParameter<bool>("invertBit", 0));
      trgMinEles.push_back(trgIt->getUntrackedParameter<unsigned int>("minEle", 1));
      trgMinEtss.push_back(trgIt->getUntrackedParameter<std::vector<double> >("minEts", std::vector<double> (1, 0.)));
      tpTrgNames.push_back(trgIt->getUntrackedParameter<std::string>("tpTriggerName", "")); 
      tagFilterNames.push_back(trgIt->getUntrackedParameter<std::string>("tagFilterName", "")); 
      probeFilterNames.push_back(trgIt->getUntrackedParameter<std::string>("probeFilterName", "")); 
   }

   TH1::SetDefaultSumw2(kTRUE);

   unsigned int trgsSize = trgVPSet_.size();
   edm::Service<TFileService> fs;
   h_gsfEle_et = fs->make<TH1F>("h_gsfEle_et", "gsfElectron E_{T}", 100, 0., 500.);
   h_nHeep = fs->make<TH1F>("h_nHeep", "# of HEEP electrons", 10, 0., 10.);
   h_total = fs->make<TH1F>("h_total", "Total # events", trgsSize, 0., trgsSize);
   h_refTrgd = fs->make<TH1F>("h_refTrgd", "Events triggered by reference trigger", trgsSize, 0., trgsSize);
   h_trgd = fs->make<TH1F>("h_trgd", "# of triggered events", trgsSize, 0., trgsSize);
   h_notTrgd = fs->make<TH1F>("h_notTrgd", "# of not triggered events", trgsSize, 0., trgsSize);
   h_notTrgdRefTrgd = fs->make<TH1F>("h_notTrgdRefTrgd", "# of not triggered events triggered by reference trigger", trgsSize, 0., trgsSize);
   h_1PassHeep_total = fs->make<TH1F>("h_1PassHeep_total", "Total events with one HEEP electron", trgsSize, 0., trgsSize);
   h_1PassHeep_refTrgd = fs->make<TH1F>("h_1PassHeep_refTrgd", "Events with one HEEP electron triggered by reference trigger", trgsSize, 0., trgsSize);
   h_2PassHeep_total = fs->make<TH1F>("h_2PassHeep_total", "Total events with two HEEP electrons", trgsSize, 0., trgsSize);
   h_2PassHeep_refTrgd = fs->make<TH1F>("h_2PassHeep_refTrgd", "Events with two HEEP electrons triggered by reference trigger", trgsSize, 0., trgsSize);
   h_nPassHeep_total = fs->make<TH1F>("h_nPassHeep_total", "Total events with n HEEP electrons (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_refTrgd = fs->make<TH1F>("h_nPassHeep_refTrgd", "Events with n HEEP electrons triggered by reference trigger (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_trgd = fs->make<TH1F>("h_nPassHeep_trgd", "Triggered events with n HEEP electrons (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_notTrgd = fs->make<TH1F>("h_nPassHeep_notTrgd", "Not triggered events with n HEEP electrons (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_notTrgdRefTrgd = fs->make<TH1F>("h_nPassHeep_notTrgdRefTrgd", "Not triggered events with n HEEP electrons triggered by reference trigger (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_total = fs->make<TH1F>("h_nPassHeep_nPassEt_total", "Total events with n HEEP electrons passing Et cut (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_refTrgd = fs->make<TH1F>("h_nPassHeep_nPassEt_refTrgd", "Events with n HEEP electrons passing Et cut triggered by reference trigger (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_trgd = fs->make<TH1F>("h_nPassHeep_nPassEt_trgd", "Triggered events with n HEEP electrons passing Et cut (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_notTrgd = fs->make<TH1F>("h_nPassHeep_nPassEt_notTrgd", "Not triggered events with n HEEP electrons passing Et cut (n defined by trigger)", trgsSize, 0., trgsSize);
   h_nPassHeep_nPassEt_notTrgdRefTrgd = fs->make<TH1F>("h_nPassHeep_nPassEt_notTrgdRefTrgd", "Not triggered events with n HEEP electrons passing Et cut triggered by reference trigger (n defined by trigger)", trgsSize, 0., trgsSize);

   hTp_tags = fs->make<TH1F>("hTp_tags", "# tags;# tags", 10, 0., 10.);
   hTp_probes = fs->make<TH1F>("hTp_probes", "# probes;# probes", 10, 0., 10.);
   hTp_allProbesGsf = fs->make<TH1F>("hTp_allProbesGsf", "# GSF probes", trgsSize, 0., trgsSize);
   hTp_passProbesGsf = fs->make<TH1F>("hTp_passProbesGsf", "# passing GSF probes", trgsSize, 0., trgsSize);
   hTp_allProbesGsf_et25 = fs->make<TH1F>("hTp_allProbesGsf_et25", "# GSF probes with E_{T} > 25 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesGsf_et25 = fs->make<TH1F>("hTp_passProbesGsf_et25", "# passing GSF probes with E_{T} > 25 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesGsf_et33 = fs->make<TH1F>("hTp_allProbesGsf_et33", "# GSF probes with E_{T} > 33 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesGsf_et33 = fs->make<TH1F>("hTp_passProbesGsf_et33", "# passing GSF probes with E_{T} > 33 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesGsf_et35 = fs->make<TH1F>("hTp_allProbesGsf_et35", "# GSF probes with E_{T} > 35 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesGsf_et35 = fs->make<TH1F>("hTp_passProbesGsf_et35", "# passing GSF probes with E_{T} > 35 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesGsf_et80 = fs->make<TH1F>("hTp_allProbesGsf_et80", "# GSF probes with E_{T} > 80 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesGsf_et80 = fs->make<TH1F>("hTp_passProbesGsf_et80", "# passing GSF probes with E_{T} > 80 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesGsf_et100 = fs->make<TH1F>("hTp_allProbesGsf_et100", "# GSF probes with E_{T} > 100 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesGsf_et100 = fs->make<TH1F>("hTp_passProbesGsf_et100", "# passing GSF probes with E_{T} > 100 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesGsf_withEtCut = fs->make<TH1F>("hTp_allProbesGsf_withEtCut", "# GSF probes passing E_{T} cut", trgsSize, 0., trgsSize);
   hTp_passProbesGsf_withEtCut = fs->make<TH1F>("hTp_passProbesGsf_withEtCut", "# passing GSF probes passing E_{T} cut", trgsSize, 0., trgsSize);
   hTp_allProbesHeep = fs->make<TH1F>("hTp_allProbesHeep", "# HEEP probes", trgsSize, 0., trgsSize);
   hTp_passProbesHeep = fs->make<TH1F>("hTp_passProbesHeep", "# passing HEEP probes", trgsSize, 0., trgsSize);
   hTp_allProbesHeep_et25 = fs->make<TH1F>("hTp_allProbesHeep_et25", "# HEEP probes with E_{T} > 25 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesHeep_et25 = fs->make<TH1F>("hTp_passProbesHeep_et25", "# passing HEEP probes with E_{T} > 25 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesHeep_et33 = fs->make<TH1F>("hTp_allProbesHeep_et33", "# HEEP probes with E_{T} > 33 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesHeep_et33 = fs->make<TH1F>("hTp_passProbesHeep_et33", "# passing HEEP probes with E_{T} > 33 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesHeep_et35 = fs->make<TH1F>("hTp_allProbesHeep_et35", "# HEEP probes with E_{T} > 35 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesHeep_et35 = fs->make<TH1F>("hTp_passProbesHeep_et35", "# passing HEEP probes with E_{T} > 35 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesHeep_et80 = fs->make<TH1F>("hTp_allProbesHeep_et80", "# HEEP probes with E_{T} > 80 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesHeep_et80 = fs->make<TH1F>("hTp_passProbesHeep_et80", "# passing HEEP probes with E_{T} > 80 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesHeep_et100 = fs->make<TH1F>("hTp_allProbesHeep_et100", "# HEEP probes with E_{T} > 100 GeV", trgsSize, 0., trgsSize);
   hTp_passProbesHeep_et100 = fs->make<TH1F>("hTp_passProbesHeep_et100", "# passing HEEP probes with E_{T} > 100 GeV", trgsSize, 0., trgsSize);
   hTp_allProbesHeep_withEtCut = fs->make<TH1F>("hTp_allProbesHeep_withEtCut", "# HEEP probes passing E_{T} cut", trgsSize, 0., trgsSize);
   hTp_passProbesHeep_withEtCut = fs->make<TH1F>("hTp_passProbesHeep_withEtCut", "# passing HEEP probes passing E_{T} cut", trgsSize, 0., trgsSize);

   // ratio histos
   hRatio_nPassHeep_trgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_trgdVsTotal", "# triggered / # total n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_notTrgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_notTrgdVsTotal", "# not triggered / # total n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_trgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_trgdVsTotal", "# triggered / # total n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_notTrgdVsTotal = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_notTrgdVsTotal", "# not triggered / # total n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);

   hRatio_trgdVsRefTrgd = fs->make<TH1F>("hRatio_trgdVsRefTrgd", "# triggered / # reference triggered events", trgsSize, 0., trgsSize);
   hRatio_notTrgdVsRefTrgd = fs->make<TH1F>("hRatio_notTrgdVsRefTrgd", "# not triggered / # reference triggered events", trgsSize, 0., trgsSize);
   hRatio_notTrgdRefTrgdVsRefTrgd = fs->make<TH1F>("hRatio_notTrgdRefTrgdVsRefTrgd", "# not triggered but reference triggered/ # reference triggered events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_trgdVsRefTrgd = fs->make<TH1F>("hRatio_nPassHeep_trgdVsRefTrgd", "# triggered / # reference triggered n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_notTrgdVsRefTrgd = fs->make<TH1F>("hRatio_nPassHeep_notTrgdVsRefTrgd", "# not triggered / # reference triggered n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_notTrgdRefTrgdVsRefTrgd = fs->make<TH1F>("hRatio_nPassHeep_notTrgdRefTrgdVsRefTrgd", "# not triggered but reference triggered / # reference triggered n-HEEP events", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_trgdVsRefTrgd = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_trgdVsRefTrgd", "# triggered / # reference triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_notTrgdVsRefTrgd = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_notTrgdVsRefTrgd", "# not triggered / # reference triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);
   hRatio_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd = fs->make<TH1F>("hRatio_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd", "# not triggered but reference triggered / # reference triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize);

   hTpRatio_passProbesGsf_vs_allProbesGsf = fs->make<TH1F>("hTpRatio_passProbesGsf_vs_allProbesGsf", "# passing GSF probes / # all GSF probes", trgsSize, 0., trgsSize);
   hTpRatio_passProbesGsf_et25_vs_allProbesGsf_et25 = fs->make<TH1F>("hTpRatio_passProbesGsf_et25_vs_allProbesGsf_et25", "# passing GSF probes / # all GSF probes with E_{T} > 25 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesGsf_et33_vs_allProbesGsf_et33 = fs->make<TH1F>("hTpRatio_passProbesGsf_et33_vs_allProbesGsf_et33", "# passing GSF probes / # all GSF probes with E_{T} > 33 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesGsf_et35_vs_allProbesGsf_et35 = fs->make<TH1F>("hTpRatio_passProbesGsf_et35_vs_allProbesGsf_et35", "# passing GSF probes / # all GSF probes with E_{T} > 35 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesGsf_et80_vs_allProbesGsf_et80 = fs->make<TH1F>("hTpRatio_passProbesGsf_et80_vs_allProbesGsf_et80", "# passing GSF probes / # all GSF probes with E_{T} > 80 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesGsf_et100_vs_allProbesGsf_et100 = fs->make<TH1F>("hTpRatio_passProbesGsf_et100_vs_allProbesGsf_et100", "# passing GSF probes / # all GSF probes with E_{T} > 100 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesGsf_withEtCut_vs_allProbesGsf_withEtCut = fs->make<TH1F>("hTpRatio_passProbesGsf_withEtCut_vs_allProbesGsf_withEtCut", "# passing GSF probes / # all GSF probes passing E_{T} cut", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_vs_allProbesHeep = fs->make<TH1F>("hTpRatio_passProbesHeep_vs_allProbesHeep", "# passing HEEP probes / # all HEEP probes", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_et25_vs_allProbesHeep_et25 = fs->make<TH1F>("hTpRatio_passProbesHeep_et25_vs_allProbesHeep_et25", "# passing HEEP probes / # all HEEP probes with E_{T} > 25 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_et33_vs_allProbesHeep_et33 = fs->make<TH1F>("hTpRatio_passProbesHeep_et33_vs_allProbesHeep_et33", "# passing HEEP probes / # all HEEP probes with E_{T} > 33 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_et35_vs_allProbesHeep_et35 = fs->make<TH1F>("hTpRatio_passProbesHeep_et35_vs_allProbesHeep_et35", "# passing HEEP probes / # all HEEP probes with E_{T} > 35 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_et80_vs_allProbesHeep_et80 = fs->make<TH1F>("hTpRatio_passProbesHeep_et80_vs_allProbesHeep_et80", "# passing HEEP probes / # all HEEP probes with E_{T} > 80 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_et100_vs_allProbesHeep_et100 = fs->make<TH1F>("hTpRatio_passProbesHeep_et100_vs_allProbesHeep_et100", "# passing HEEP probes / # all HEEP probes with E_{T} > 100 GeV", trgsSize, 0., trgsSize);
   hTpRatio_passProbesHeep_withEtCut_vs_allProbesHeep_withEtCut = fs->make<TH1F>("hTpRatio_passProbesHeep_withEtCut_vs_allProbesHeep_withEtCut", "# passing HEEP probes / # all HEEP probes passing E_{T} cut", trgsSize, 0., trgsSize);

   h2Ratio_trgd = fs->make<TH2F>("h2Ratio_trgd", "Inter trigger ratio of triggered events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_notTrgd = fs->make<TH2F>("h2Ratio_notTrgd", "Inter trigger ratio of not triggered events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_notTrgdRefTrgd = fs->make<TH2F>("h2Ratio_notTrgdRefTrgd", "Inter trigger ratio of not triggered events triggered by the reference trigger", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_trgd = fs->make<TH2F>("h2Ratio_nPassHeep_trgd", "Inter trigger ratio of triggered n-HEEP events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_notTrgd = fs->make<TH2F>("h2Ratio_nPassHeep_notTrgd", "Inter trigger ratio of not triggered n-HEEP events", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_notTrgdRefTrgd = fs->make<TH2F>("h2Ratio_nPassHeep_notTrgdRefTrgd", "Inter trigger ratio of not triggered n-HEEP events triggered by the reference trigger", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_nPassEt_trgd = fs->make<TH2F>("h2Ratio_nPassHeep_nPassEt_trgd", "Inter trigger ratio of triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_nPassEt_notTrgd = fs->make<TH2F>("h2Ratio_nPassHeep_nPassEt_notTrgd", "Inter trigger ratio of not triggered n-HEEP events passing Et cuts", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);
   h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd = fs->make<TH2F>("h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd", "Inter trigger ratio of not triggered n-HEEP events passing Et cuts triggered by the reference trigger", trgsSize, 0., trgsSize, trgsSize, 0., trgsSize);

   for (unsigned int i = 0; i < trgNames.size(); ++i) {
      float lastTrgMinEt = trgMinEtss.at(i).back(); 
      std::string nameString = "h_nPassHeep_nMinus1PassEt_vsEt_total_" + trgNames[i];
      v_h_nPassHeep_nMinus1PassEt_vsEt_total.push_back(fs->make<TH1F>(nameString.data(), "# of events vs. Et cut;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "h_nPassHeep_nMinus1PassEt_vsEt_refTrgd_" + trgNames[i];
      v_h_nPassHeep_nMinus1PassEt_vsEt_refTrgd.push_back(fs->make<TH1F>(nameString.data(), "# of reference trigger triggered events vs. Et cut;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "h_nPassHeep_nMinus1PassEt_vsEt_trgd_" + trgNames[i];
      v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.push_back(fs->make<TH1F>(nameString.data(), "# of triggered events vs. Et cut;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal_" + trgNames[i];
      v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal.push_back(fs->make<TH1F>(nameString.data(), "# triggered / # total n-HEEP events vs. Et cuts;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsRefTrgd_" + trgNames[i];
      v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsRefTrgd.push_back(fs->make<TH1F>(nameString.data(), "# triggered / # reference triggered n-HEEP events vs. Et cuts;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));

      nameString = "hTp_allProbesGsf_et_" + trgNames[i];
      v_hTp_allProbesGsf_et.push_back(fs->make<TH1F>(nameString.data(), "E_{T} of GSF probes;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hTp_passProbesGsf_et_" + trgNames[i];
      v_hTp_passProbesGsf_et.push_back(fs->make<TH1F>(nameString.data(), "E_{T} of passing GSF probes;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hTpRatio_passProbesGsf_et_vs_allProbesGsf_et_" + trgNames[i];
      v_hTpRatio_passProbesGsf_et_vs_allProbesGsf_et.push_back(fs->make<TH1F>(nameString.data(), "# pass GSF probes / # all GSF probes;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hTp_allProbesHeep_et_" + trgNames[i];
      v_hTp_allProbesHeep_et.push_back(fs->make<TH1F>(nameString.data(), "E_{T} of HEEP probes;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hTp_passProbesHeep_et_" + trgNames[i];
      v_hTp_passProbesHeep_et.push_back(fs->make<TH1F>(nameString.data(), "E_{T} of passing HEEP probes;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));
      nameString = "hTpRatio_passProbesHeep_et_vs_allProbesHeep_et_" + trgNames[i];
      v_hTpRatio_passProbesHeep_et_vs_allProbesHeep_et.push_back(fs->make<TH1F>(nameString.data(), "# pass HEEP probes / # all HEEP probes;E_{T} (GeV)", 2*deltaEtCut/etCutStepSize, lastTrgMinEt-deltaEtCut, lastTrgMinEt+deltaEtCut));

      // set bin labels
      std::string triggerName = trgNames[i];
      if (triggerName.compare(refTrgNames[i]) == 0) triggerName += " [ref]";
      const char* trgName = triggerName.data();
      h_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_refTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_notTrgdRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_1PassHeep_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_1PassHeep_refTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_2PassHeep_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_2PassHeep_refTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_refTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_notTrgdRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_total->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_refTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h_nPassHeep_nPassEt_notTrgdRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);

      hTp_allProbesGsf->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesGsf_et25->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf_et25->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesGsf_et33->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf_et33->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesGsf_et35->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf_et35->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesGsf_et80->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf_et80->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesGsf_et100->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf_et100->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesGsf_withEtCut->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesGsf_withEtCut->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep_et25->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep_et25->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep_et33->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep_et33->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep_et35->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep_et35->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep_et80->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep_et80->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep_et100->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep_et100->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_allProbesHeep_withEtCut->GetXaxis()->SetBinLabel(i+1, trgName);
      hTp_passProbesHeep_withEtCut->GetXaxis()->SetBinLabel(i+1, trgName);

      // ratio histos
      hRatio_nPassHeep_trgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_notTrgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_trgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_notTrgdVsTotal->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_trgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_notTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_notTrgdRefTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_trgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_notTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_notTrgdRefTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_trgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_notTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      hRatio_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);

      hTpRatio_passProbesGsf_vs_allProbesGsf->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesGsf_et25_vs_allProbesGsf_et25->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesGsf_et33_vs_allProbesGsf_et33->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesGsf_et35_vs_allProbesGsf_et35->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesGsf_et80_vs_allProbesGsf_et80->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesGsf_et100_vs_allProbesGsf_et100->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesGsf_withEtCut_vs_allProbesGsf_withEtCut->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_vs_allProbesHeep->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_et25_vs_allProbesHeep_et25->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_et33_vs_allProbesHeep_et33->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_et35_vs_allProbesHeep_et35->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_et80_vs_allProbesHeep_et80->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_et100_vs_allProbesHeep_et100->GetXaxis()->SetBinLabel(i+1, trgName);
      hTpRatio_passProbesHeep_withEtCut_vs_allProbesHeep_withEtCut->GetXaxis()->SetBinLabel(i+1, trgName);

      h2Ratio_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_trgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_notTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_notTrgdRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_notTrgdRefTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_trgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_notTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_notTrgdRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_notTrgdRefTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_trgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_trgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_notTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_notTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd->GetYaxis()->SetBinLabel(i+1, trgName);
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
      h_gsfEle_et->Fill(eleIt->caloEnergy() * sin(eleIt->p4().theta()));
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
      bool refTriggered = trgRes.accept(refTrgIndices[trgIndex]);
      bool tpTriggered = trgRes.accept(tpTrgIndices[trgIndex]);

      // tag and probe
      if (tpTriggered) {
         std::vector<unsigned int> tagEleIndices;
         std::vector<unsigned int> probeEleIndicesHeep;
         for (unsigned int heepInd = 0; heepInd < heepEleIndices.size(); ++heepInd) {
            float eleEt = gsfElectrons->at(heepInd).caloEnergy() * sin(gsfElectrons->at(heepInd).p4().theta());
            if (eleEt > 35. && tagFilterNames[trgIndex].size() == 0) tagEleIndices.push_back(heepInd);
            else if (eleEt > 35. && MatchTrgFltr(iEvent, gsfElectrons->at(heepInd), tagFilterNames[trgIndex].data())) tagEleIndices.push_back(heepInd);
            if (probeFilterNames[trgIndex].size() == 0) probeEleIndicesHeep.push_back(heepInd);
            else if (MatchTrgFltr(iEvent, gsfElectrons->at(heepInd), probeFilterNames[trgIndex].data())) probeEleIndicesHeep.push_back(heepInd);
         }
         hTp_tags->Fill(tagEleIndices.size());
         hTp_probes->Fill(probeEleIndicesHeep.size());
         for (unsigned int tagInd = 0; tagInd < tagEleIndices.size(); ++tagInd) {
            // loop over all GSF electrons
            for (unsigned int gsfInd = 0; gsfInd < gsfElectrons->size(); ++gsfInd) {
               // continue if it is the tag
               if (gsfInd == tagInd) continue;
               // match probe candidate to trigger filter if defined
               if (!(probeFilterNames[trgIndex].size() == 0 || MatchTrgFltr(iEvent, gsfElectrons->at(gsfInd), probeFilterNames[trgIndex].data()))) continue;
               float eleEt = gsfElectrons->at(gsfInd).caloEnergy() * sin(gsfElectrons->at(gsfInd).p4().theta());
               if (eleEt > 100.) {
                  hTp_allProbesGsf_et100->Fill(trgIndex);
                  hTp_allProbesGsf_et80->Fill(trgIndex);
                  hTp_allProbesGsf_et35->Fill(trgIndex);
                  hTp_allProbesGsf_et33->Fill(trgIndex);
                  hTp_allProbesGsf_et25->Fill(trgIndex);
                  hTp_allProbesGsf->Fill(trgIndex);
               } else if (eleEt > 80.) {
                  hTp_allProbesGsf_et80->Fill(trgIndex);
                  hTp_allProbesGsf_et35->Fill(trgIndex);
                  hTp_allProbesGsf_et33->Fill(trgIndex);
                  hTp_allProbesGsf_et25->Fill(trgIndex);
                  hTp_allProbesGsf->Fill(trgIndex);
               } else if (eleEt > 35.) {
                  hTp_allProbesGsf_et35->Fill(trgIndex);
                  hTp_allProbesGsf_et33->Fill(trgIndex);
                  hTp_allProbesGsf_et25->Fill(trgIndex);
                  hTp_allProbesGsf->Fill(trgIndex);
               } else if (eleEt > 33.) {
                  hTp_allProbesGsf_et33->Fill(trgIndex);
                  hTp_allProbesGsf_et25->Fill(trgIndex);
                  hTp_allProbesGsf->Fill(trgIndex);
               } else if (eleEt > 25.) {
                  hTp_allProbesGsf_et25->Fill(trgIndex);
                  hTp_allProbesGsf->Fill(trgIndex);
               } else {
                  hTp_allProbesGsf->Fill(trgIndex);
               }
               if (eleEt > trgMinEtss.at(trgIndex).front()) hTp_allProbesGsf_withEtCut->Fill(trgIndex);
               v_hTp_allProbesGsf_et.at(trgIndex)->Fill(eleEt);

               // assume that the last filter hast saveTags == true and use it to see if the candidate fired the trigger
               std::vector<std::string> saveTagsModules = hltConfig.saveTagsModules(trgIndices[trgIndex]);
               bool matchLast = MatchTrgFltr(iEvent, gsfElectrons->at(gsfInd), saveTagsModules.back().data());

               if (matchLast) {
                  if (eleEt > 100.) {
                     hTp_passProbesGsf_et100->Fill(trgIndex);
                     hTp_passProbesGsf_et80->Fill(trgIndex);
                     hTp_passProbesGsf_et35->Fill(trgIndex);
                     hTp_passProbesGsf_et33->Fill(trgIndex);
                     hTp_passProbesGsf_et25->Fill(trgIndex);
                     hTp_passProbesGsf->Fill(trgIndex);
                  } else if (eleEt > 80.) {
                     hTp_passProbesGsf_et80->Fill(trgIndex);
                     hTp_passProbesGsf_et35->Fill(trgIndex);
                     hTp_passProbesGsf_et33->Fill(trgIndex);
                     hTp_passProbesGsf_et25->Fill(trgIndex);
                     hTp_passProbesGsf->Fill(trgIndex);
                  } else if (eleEt > 35.) {
                     hTp_passProbesGsf_et35->Fill(trgIndex);
                     hTp_passProbesGsf_et33->Fill(trgIndex);
                     hTp_passProbesGsf_et25->Fill(trgIndex);
                     hTp_passProbesGsf->Fill(trgIndex);
                  } else if (eleEt > 33.) {
                     hTp_passProbesGsf_et33->Fill(trgIndex);
                     hTp_passProbesGsf_et25->Fill(trgIndex);
                     hTp_passProbesGsf->Fill(trgIndex);
                  } else if (eleEt > 25.) {
                     hTp_passProbesGsf_et25->Fill(trgIndex);
                     hTp_passProbesGsf->Fill(trgIndex);
                  } else {
                     hTp_passProbesGsf->Fill(trgIndex);
                  }
                  if (eleEt > trgMinEtss.at(trgIndex).front()) hTp_passProbesGsf_withEtCut->Fill(trgIndex);
                  v_hTp_passProbesGsf_et.at(trgIndex)->Fill(eleEt);
               }

               // Only HEEP electrons may pass here
               if (!PassHeep(gsfElectrons->at(gsfInd), rho)) continue;

               if (eleEt > 100.) {
                  hTp_allProbesHeep_et100->Fill(trgIndex);
                  hTp_allProbesHeep_et80->Fill(trgIndex);
                  hTp_allProbesHeep_et35->Fill(trgIndex);
                  hTp_allProbesHeep_et33->Fill(trgIndex);
                  hTp_allProbesHeep_et25->Fill(trgIndex);
                  hTp_allProbesHeep->Fill(trgIndex);
               } else if (eleEt > 80.) {
                  hTp_allProbesHeep_et80->Fill(trgIndex);
                  hTp_allProbesHeep_et35->Fill(trgIndex);
                  hTp_allProbesHeep_et33->Fill(trgIndex);
                  hTp_allProbesHeep_et25->Fill(trgIndex);
                  hTp_allProbesHeep->Fill(trgIndex);
               } else if (eleEt > 35.) {
                  hTp_allProbesHeep_et35->Fill(trgIndex);
                  hTp_allProbesHeep_et33->Fill(trgIndex);
                  hTp_allProbesHeep_et25->Fill(trgIndex);
                  hTp_allProbesHeep->Fill(trgIndex);
               } else if (eleEt > 33.) {
                  hTp_allProbesHeep_et33->Fill(trgIndex);
                  hTp_allProbesHeep_et25->Fill(trgIndex);
                  hTp_allProbesHeep->Fill(trgIndex);
               } else if (eleEt > 25.) {
                  hTp_allProbesHeep_et25->Fill(trgIndex);
                  hTp_allProbesHeep->Fill(trgIndex);
               } else {
                  hTp_allProbesHeep->Fill(trgIndex);
               }
               if (eleEt > trgMinEtss.at(trgIndex).front()) hTp_allProbesHeep_withEtCut->Fill(trgIndex);
               v_hTp_allProbesHeep_et.at(trgIndex)->Fill(eleEt);

               if (matchLast) {
                  if (eleEt > 100.) {
                     hTp_passProbesHeep_et100->Fill(trgIndex);
                     hTp_passProbesHeep_et80->Fill(trgIndex);
                     hTp_passProbesHeep_et35->Fill(trgIndex);
                     hTp_passProbesHeep_et33->Fill(trgIndex);
                     hTp_passProbesHeep_et25->Fill(trgIndex);
                     hTp_passProbesHeep->Fill(trgIndex);
                  } else if (eleEt > 80.) {
                     hTp_passProbesHeep_et80->Fill(trgIndex);
                     hTp_passProbesHeep_et35->Fill(trgIndex);
                     hTp_passProbesHeep_et33->Fill(trgIndex);
                     hTp_passProbesHeep_et25->Fill(trgIndex);
                     hTp_passProbesHeep->Fill(trgIndex);
                  } else if (eleEt > 35.) {
                     hTp_passProbesHeep_et35->Fill(trgIndex);
                     hTp_passProbesHeep_et33->Fill(trgIndex);
                     hTp_passProbesHeep_et25->Fill(trgIndex);
                     hTp_passProbesHeep->Fill(trgIndex);
                  } else if (eleEt > 33.) {
                     hTp_passProbesHeep_et33->Fill(trgIndex);
                     hTp_passProbesHeep_et25->Fill(trgIndex);
                     hTp_passProbesHeep->Fill(trgIndex);
                  } else if (eleEt > 25.) {
                     hTp_passProbesHeep_et25->Fill(trgIndex);
                     hTp_passProbesHeep->Fill(trgIndex);
                  } else {
                     hTp_passProbesHeep->Fill(trgIndex);
                  }
                  if (eleEt > trgMinEtss.at(trgIndex).front()) hTp_passProbesHeep_withEtCut->Fill(trgIndex);
                  v_hTp_passProbesHeep_et.at(trgIndex)->Fill(eleEt);
               }
            }
         }
      }

      h_total->Fill(trgIndex);
      if (refTriggered) {
         h_refTrgd->Fill(trgIndex);
         if (nPassHeep > 0) h_1PassHeep_refTrgd->Fill(trgIndex);
         if (nPassHeep > 1) h_2PassHeep_refTrgd->Fill(trgIndex);
      }
      if (nPassHeep > 0) h_1PassHeep_total->Fill(trgIndex);
      if (nPassHeep > 1) h_2PassHeep_total->Fill(trgIndex);

      if (nPassHeep >= trgMinEles[trgIndex]) {
         h_nPassHeep_total->Fill(trgIndex);
         if (refTriggered) h_nPassHeep_refTrgd->Fill(trgIndex);
         if (triggered) h_nPassHeep_trgd->Fill(trgIndex);
         else {
            h_nPassHeep_notTrgd->Fill(trgIndex);
            if (refTriggered) h_nPassHeep_notTrgdRefTrgd->Fill(trgIndex);
         }

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
            if (refTriggered) h_nPassHeep_nPassEt_refTrgd->Fill(trgIndex);
            if (triggered) h_nPassHeep_nPassEt_trgd->Fill(trgIndex);
            else {
               h_nPassHeep_nPassEt_notTrgd->Fill(trgIndex);
               if (refTriggered) h_nPassHeep_nPassEt_notTrgdRefTrgd->Fill(trgIndex);
            }
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
               if (refTriggered) v_h_nPassHeep_nMinus1PassEt_vsEt_refTrgd.at(trgIndex)->Fill(etCut);
               if (triggered) v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.at(trgIndex)->Fill(etCut);
            }
         }
      }

      if (triggered) h_trgd->Fill(trgIndex);
      else {
         h_notTrgd->Fill(trgIndex);
         if (refTriggered) h_notTrgdRefTrgd->Fill(trgIndex);
      }
   }
}

bool
HltGsfEleAnalyser::PassHeep(reco::GsfElectron const& ele, double &rho) 
{
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

bool
HltGsfEleAnalyser::MatchTrgFltr(const edm::Event& iEvent, reco::GsfElectron const& ele, const char* filterName) 
{
   bool pass = 0;

   edm::Handle<trigger::TriggerEvent> trgEvent; 
   if (!iEvent.getByLabel(trgEventTag_, trgEvent)) return pass;

   trigger::size_type filterIndex = trgEvent->filterIndex(edm::InputTag(filterName, "", trgEventTag_.process())); 
   if (filterIndex < trgEvent->sizeFilters()){ 
      const trigger::Keys& trgKeys = trgEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection& trgObjColl(trgEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt = trgKeys.begin(); keyIt != trgKeys.end(); ++keyIt) { 
         const trigger::TriggerObject& obj = trgObjColl[*keyIt];
         if (deltaR(ele.eta(), ele.phi(), obj.eta(), obj.phi()) < 0.5) {
            pass = 1;
            break;
         }
      }   
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

   hRatio_trgdVsRefTrgd->Divide(h_trgd, h_refTrgd);
   hRatio_notTrgdVsRefTrgd->Divide(h_notTrgd, h_refTrgd);
   hRatio_notTrgdRefTrgdVsRefTrgd->Divide(h_notTrgdRefTrgd, h_refTrgd);
   hRatio_nPassHeep_trgdVsRefTrgd->Divide(h_nPassHeep_trgd, h_nPassHeep_refTrgd);
   hRatio_nPassHeep_notTrgdVsRefTrgd->Divide(h_nPassHeep_notTrgd, h_nPassHeep_refTrgd);
   hRatio_nPassHeep_notTrgdRefTrgdVsRefTrgd->Divide(h_nPassHeep_notTrgdRefTrgd, h_nPassHeep_refTrgd);
   hRatio_nPassHeep_nPassEt_trgdVsRefTrgd->Divide(h_nPassHeep_nPassEt_trgd, h_nPassHeep_nPassEt_refTrgd);
   hRatio_nPassHeep_nPassEt_notTrgdVsRefTrgd->Divide(h_nPassHeep_nPassEt_notTrgd, h_nPassHeep_nPassEt_refTrgd);
   hRatio_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd->Divide(h_nPassHeep_nPassEt_notTrgdRefTrgd, h_nPassHeep_nPassEt_refTrgd);
   hTpRatio_passProbesGsf_vs_allProbesGsf->Divide(hTp_passProbesGsf, hTp_allProbesGsf);
   hTpRatio_passProbesGsf_et25_vs_allProbesGsf_et25->Divide(hTp_passProbesGsf_et25, hTp_allProbesGsf_et25);
   hTpRatio_passProbesGsf_et33_vs_allProbesGsf_et33->Divide(hTp_passProbesGsf_et33, hTp_allProbesGsf_et33);
   hTpRatio_passProbesGsf_et35_vs_allProbesGsf_et35->Divide(hTp_passProbesGsf_et35, hTp_allProbesGsf_et35);
   hTpRatio_passProbesGsf_et80_vs_allProbesGsf_et80->Divide(hTp_passProbesGsf_et80, hTp_allProbesGsf_et80);
   hTpRatio_passProbesGsf_et100_vs_allProbesGsf_et100->Divide(hTp_passProbesGsf_et100, hTp_allProbesGsf_et100);
   hTpRatio_passProbesGsf_withEtCut_vs_allProbesGsf_withEtCut->Divide(hTp_passProbesGsf_withEtCut, hTp_allProbesGsf_withEtCut);
   hTpRatio_passProbesHeep_vs_allProbesHeep->Divide(hTp_passProbesHeep, hTp_allProbesHeep);
   hTpRatio_passProbesHeep_et25_vs_allProbesHeep_et25->Divide(hTp_passProbesHeep_et25, hTp_allProbesHeep_et25);
   hTpRatio_passProbesHeep_et33_vs_allProbesHeep_et33->Divide(hTp_passProbesHeep_et33, hTp_allProbesHeep_et33);
   hTpRatio_passProbesHeep_et35_vs_allProbesHeep_et35->Divide(hTp_passProbesHeep_et35, hTp_allProbesHeep_et35);
   hTpRatio_passProbesHeep_et80_vs_allProbesHeep_et80->Divide(hTp_passProbesHeep_et80, hTp_allProbesHeep_et80);
   hTpRatio_passProbesHeep_et100_vs_allProbesHeep_et100->Divide(hTp_passProbesHeep_et100, hTp_allProbesHeep_et100);
   hTpRatio_passProbesHeep_withEtCut_vs_allProbesHeep_withEtCut->Divide(hTp_passProbesHeep_withEtCut, hTp_allProbesHeep_withEtCut);
   for (unsigned int i = 0; i < trgNames.size(); ++i) {
      v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsTotal.at(i)->Divide(v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.at(i), v_h_nPassHeep_nMinus1PassEt_vsEt_total.at(i));
      v_hRatio_nPassHeep_nMinus1PassEt_vsEt_trgdVsRefTrgd.at(i)->Divide(v_h_nPassHeep_nMinus1PassEt_vsEt_trgd.at(i), v_h_nPassHeep_nMinus1PassEt_vsEt_refTrgd.at(i));
      v_hTpRatio_passProbesGsf_et_vs_allProbesGsf_et.at(i)->Divide(v_hTp_passProbesGsf_et.at(i), v_hTp_allProbesGsf_et.at(i));
      v_hTpRatio_passProbesHeep_et_vs_allProbesHeep_et.at(i)->Divide(v_hTp_passProbesHeep_et.at(i), v_hTp_allProbesHeep_et.at(i));
   }

   edm::Service<TFileService> fs;
   gEff_trgdVsRefTrgd = fs->make<TGraphAsymmErrors>(h_trgd, h_refTrgd, "cp");
   gEff_notTrgdRefTrgdVsRefTrgd = fs->make<TGraphAsymmErrors>(h_notTrgdRefTrgd, h_refTrgd, "cp");
   gEff_nPassHeep_trgdVsRefTrgd = fs->make<TGraphAsymmErrors>(h_nPassHeep_trgd, h_nPassHeep_refTrgd, "cp");
   gEff_nPassHeep_notTrgdRefTrgdVsRefTrgd = fs->make<TGraphAsymmErrors>(h_nPassHeep_notTrgdRefTrgd, h_nPassHeep_refTrgd, "cp");
   gEff_nPassHeep_nPassEt_trgdVsRefTrgd = fs->make<TGraphAsymmErrors>(h_nPassHeep_nPassEt_trgd, h_nPassHeep_nPassEt_refTrgd, "cp");
   gEff_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd = fs->make<TGraphAsymmErrors>(h_nPassHeep_nPassEt_notTrgdRefTrgd, h_nPassHeep_nPassEt_refTrgd, "cp");

   gEff_trgdVsRefTrgd->SetNameTitle("gEff_trgdVsRefTrgd", "# triggered / # reference triggered events");
   gEff_notTrgdRefTrgdVsRefTrgd->SetNameTitle("gEff_notTrgdRefTrgdVsRefTrgd", "# not triggered but reference triggered / # reference triggered events");
   gEff_nPassHeep_trgdVsRefTrgd->SetNameTitle("gEff_nPassHeep_trgdVsRefTrgd", "# triggered / # reference triggered n-HEEP events");
   gEff_nPassHeep_notTrgdRefTrgdVsRefTrgd->SetNameTitle("gEff_nPassHeep_notTrgdRefTrgdVsRefTrgd", "# not triggered but reference triggered / # reference triggered n-HEEP events");
   gEff_nPassHeep_nPassEt_trgdVsRefTrgd->SetNameTitle("gEff_nPassHeep_nPassEt_trgdVsRefTrgd", "# triggered / # reference triggered n-HEEP events passing Et cuts");
   gEff_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd->SetNameTitle("gEff_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd", "# not triggered but reference triggered / # reference triggered n-HEEP events passing Et cuts");
   for (unsigned int i = 0; i < trgNames.size(); ++i) {
      std::string triggerName = trgNames[i];
      if (triggerName.compare(refTrgNames[i]) == 0) triggerName += " [ref]";
      const char* trgName = triggerName.data();
      gEff_trgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      gEff_notTrgdRefTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      gEff_nPassHeep_trgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      gEff_nPassHeep_notTrgdRefTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      gEff_nPassHeep_nPassEt_trgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
      gEff_nPassHeep_nPassEt_notTrgdRefTrgdVsRefTrgd->GetXaxis()->SetBinLabel(i+1, trgName);
   }

   TH2F* h2Denom_trgd = (TH2F*)h2Ratio_trgd->Clone("h2Denom_trgd");
   TH2F* h2Denom_notTrgd = (TH2F*)h2Ratio_notTrgd->Clone("h2Denom_notTrgd");
   TH2F* h2Denom_notTrgdRefTrgd = (TH2F*)h2Ratio_notTrgdRefTrgd->Clone("h2Denom_notTrgdRefTrgd");
   TH2F* h2Denom_nPassHeep_trgd = (TH2F*)h2Ratio_nPassHeep_trgd->Clone("h2Denom_nPassHeep_trgd");
   TH2F* h2Denom_nPassHeep_notTrgd = (TH2F*)h2Ratio_nPassHeep_notTrgd->Clone("h2Denom_nPassHeep_notTrgd");
   TH2F* h2Denom_nPassHeep_notTrgdRefTrgd = (TH2F*)h2Ratio_nPassHeep_notTrgdRefTrgd->Clone("h2Denom_nPassHeep_notTrgdRefTrgd");
   TH2F* h2Denom_nPassHeep_nPassEt_trgd = (TH2F*)h2Ratio_nPassHeep_nPassEt_trgd->Clone("h2Denom_nPassHeep_nPassEt_trgd");
   TH2F* h2Denom_nPassHeep_nPassEt_notTrgd = (TH2F*)h2Ratio_nPassHeep_nPassEt_notTrgd->Clone("h2Denom_nPassHeep_nPassEt_notTrgd");
   TH2F* h2Denom_nPassHeep_nPassEt_notTrgdRefTrgd = (TH2F*)h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd->Clone("h2Denom_nPassHeep_nPassEt_notTrgdRefTrgd");
   for (unsigned int ix = 1; ix <= trgNames.size(); ++ix) {
      for (unsigned int iy = 1; iy <= trgNames.size(); ++iy) {
         h2Ratio_trgd->SetBinContent(ix, iy, h_trgd->GetBinContent(ix));
         h2Ratio_trgd->SetBinError(ix, iy, h_trgd->GetBinError(ix));
         h2Denom_trgd->SetBinContent(ix, iy, h_trgd->GetBinContent(iy));
         h2Denom_trgd->SetBinError(ix, iy, h_trgd->GetBinError(iy));
         h2Ratio_notTrgd->SetBinContent(ix, iy, h_notTrgd->GetBinContent(ix));
         h2Ratio_notTrgd->SetBinError(ix, iy, h_notTrgd->GetBinError(ix));
         h2Denom_notTrgd->SetBinContent(ix, iy, h_notTrgd->GetBinContent(iy));
         h2Denom_notTrgd->SetBinError(ix, iy, h_notTrgd->GetBinError(iy));
         h2Ratio_notTrgdRefTrgd->SetBinContent(ix, iy, h_notTrgdRefTrgd->GetBinContent(ix));
         h2Ratio_notTrgdRefTrgd->SetBinError(ix, iy, h_notTrgdRefTrgd->GetBinError(ix));
         h2Denom_notTrgdRefTrgd->SetBinContent(ix, iy, h_notTrgdRefTrgd->GetBinContent(iy));
         h2Denom_notTrgdRefTrgd->SetBinError(ix, iy, h_notTrgdRefTrgd->GetBinError(iy));

         h2Ratio_nPassHeep_trgd->SetBinContent(ix, iy, h_nPassHeep_trgd->GetBinContent(ix));
         h2Ratio_nPassHeep_trgd->SetBinError(ix, iy, h_nPassHeep_trgd->GetBinError(ix));
         h2Denom_nPassHeep_trgd->SetBinContent(ix, iy, h_nPassHeep_trgd->GetBinContent(iy));
         h2Denom_nPassHeep_trgd->SetBinError(ix, iy, h_nPassHeep_trgd->GetBinError(iy));
         h2Ratio_nPassHeep_notTrgd->SetBinContent(ix, iy, h_nPassHeep_notTrgd->GetBinContent(ix));
         h2Ratio_nPassHeep_notTrgd->SetBinError(ix, iy, h_nPassHeep_notTrgd->GetBinError(ix));
         h2Denom_nPassHeep_notTrgd->SetBinContent(ix, iy, h_nPassHeep_notTrgd->GetBinContent(iy));
         h2Denom_nPassHeep_notTrgd->SetBinError(ix, iy, h_nPassHeep_notTrgd->GetBinError(iy));
         h2Ratio_nPassHeep_notTrgdRefTrgd->SetBinContent(ix, iy, h_nPassHeep_notTrgdRefTrgd->GetBinContent(ix));
         h2Ratio_nPassHeep_notTrgdRefTrgd->SetBinError(ix, iy, h_nPassHeep_notTrgdRefTrgd->GetBinError(ix));
         h2Denom_nPassHeep_notTrgdRefTrgd->SetBinContent(ix, iy, h_nPassHeep_notTrgdRefTrgd->GetBinContent(iy));
         h2Denom_nPassHeep_notTrgdRefTrgd->SetBinError(ix, iy, h_nPassHeep_notTrgdRefTrgd->GetBinError(iy));

         h2Ratio_nPassHeep_nPassEt_trgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinContent(ix));
         h2Ratio_nPassHeep_nPassEt_trgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinError(ix));
         h2Denom_nPassHeep_nPassEt_trgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinContent(iy));
         h2Denom_nPassHeep_nPassEt_trgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_trgd->GetBinError(iy));
         h2Ratio_nPassHeep_nPassEt_notTrgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinContent(ix));
         h2Ratio_nPassHeep_nPassEt_notTrgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinError(ix));
         h2Denom_nPassHeep_nPassEt_notTrgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinContent(iy));
         h2Denom_nPassHeep_nPassEt_notTrgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_notTrgd->GetBinError(iy));
         h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_notTrgdRefTrgd->GetBinContent(ix));
         h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_notTrgdRefTrgd->GetBinError(ix));
         h2Denom_nPassHeep_nPassEt_notTrgdRefTrgd->SetBinContent(ix, iy, h_nPassHeep_nPassEt_notTrgdRefTrgd->GetBinContent(iy));
         h2Denom_nPassHeep_nPassEt_notTrgdRefTrgd->SetBinError(ix, iy, h_nPassHeep_nPassEt_notTrgdRefTrgd->GetBinError(iy));
      }
   }
   h2Ratio_trgd->Divide(h2Denom_trgd);
   h2Ratio_notTrgd->Divide(h2Denom_notTrgd);
   h2Ratio_notTrgdRefTrgd->Divide(h2Denom_notTrgdRefTrgd);
   h2Ratio_nPassHeep_trgd->Divide(h2Denom_nPassHeep_trgd);
   h2Ratio_nPassHeep_notTrgd->Divide(h2Denom_nPassHeep_notTrgd);
   h2Ratio_nPassHeep_notTrgdRefTrgd->Divide(h2Denom_nPassHeep_notTrgdRefTrgd);
   h2Ratio_nPassHeep_nPassEt_trgd->Divide(h2Denom_nPassHeep_nPassEt_trgd);
   h2Ratio_nPassHeep_nPassEt_notTrgd->Divide(h2Denom_nPassHeep_nPassEt_notTrgd);
   h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd->Divide(h2Denom_nPassHeep_nPassEt_notTrgdRefTrgd);
}

// ------------ method called when starting to processes a run  ------------
void 
HltGsfEleAnalyser::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
   bool changed = 1;
   if (hltConfig.init(iRun, iSetup, trgResultsTag_.process(), changed)) {
      if (changed) {
         trgIndices.clear();
         refTrgIndices.clear();
         for (unsigned int i = 0; i < trgNames.size(); ++i) {
            trgIndices.push_back(hltConfig.triggerIndex(trgNames[i]));
            refTrgIndices.push_back(hltConfig.triggerIndex(refTrgNames[i]));
            tpTrgIndices.push_back(hltConfig.triggerIndex(tpTrgNames[i]));
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
