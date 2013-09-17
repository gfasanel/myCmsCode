// -*- C++ -*-
//
// Package:    HltGsfTrkAnalyser
// Class:      HltGsfTrkAnalyser
// 
/**\class HltGsfTrkAnalyser HltGsfTrkAnalyser.cc UserCode/HltGsfTrkAnalyser/src/HltGsfTrkAnalyser.cc

 Description: Comparing of gsf tracking algorithms at HLT level.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Reis
//         Created:  Wed Apr 10 16:56:33 CEST 2013
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

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
//
// class declaration
//

class HltGsfTrkAnalyser : public edm::EDAnalyzer {
   public:
      explicit HltGsfTrkAnalyser(const edm::ParameterSet&);
      ~HltGsfTrkAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      double TrackMatcher(reco::GsfTrack const&, reco::GsfTrack const &);
      unsigned int GetNSharedHits(const reco::GsfTrack &t1, const reco::GsfTrack &t2);
      // ----------member data ---------------------------
      edm::InputTag referenceTrackCollTag_;
      edm::InputTag testTrackCollTag_;
      edm::InputTag referenceTrackVarsTag_;
      edm::InputTag testTrackVarsTag_;
      edm::InputTag ecalCandTagForDEtaR_;
      edm::InputTag ecalCandTagForDEtaT_;
      edm::InputTag ecalCandTagForDPhiR_;
      edm::InputTag ecalCandTagForDPhiT_;
      double dEtaCutR_;
      double dPhiCutR_;
      double dEtaCutT_;
      double dPhiCutT_;
      unsigned int matchMethod_;
      double discrCut_;
      bool largerThan_;
      double minPtR_;
      double minPtT_;
      edm::ParameterSet refTrgPSet_;
      edm::ParameterSet testTrgPSet_;
      std::string refTrgName_;
      bool refTrgInv_;
      std::string testTrgName_;
      bool testTrgInv_;

      TH1D* hR_size;
      TH1D* hR_p;
      TH1D* hR_pt;
      TH1D* hR_eta;
      TH1D* hR_phi;
      TH1D* hR_ptError;
      TH1D* hR_dPtOverPt;
      TH1D* hR_etaError;
      TH1D* hR_phiError;
      TH1D* hR_outerP;
      TH1D* hR_outerPt;
      TH1D* hR_outerEta;
      TH1D* hR_outerPhi;
      TH1D* hR_outerPOverP;
      TH1D* hR_outerPtOverPt;
      TH1D* hR_chi2;
      TH1D* hR_ndof;
      TH1D* hR_normalizedChi2;
      TH1D* hR_charge;
      TH1D* hR_dxy;
      TH1D* hR_dz;
      TH1D* hR_dxyError;
      TH1D* hR_dzError;
      TH1D* hR_numberOfValidHits;
      TH1D* hR_numberOfLostHits;
      TH1D* hR_validFraction;
      TH2D* hR2d_pt_eta;
      TH2D* hR2d_pt_dpt;
      TH2D* hR2d_p_outerP;
      TH2D* hR2d_pt_outerPt;

      TH1D* hR_ptCut_size;
      TH1D* hR_ptCut_p;
      TH1D* hR_ptCut_pt;
      TH1D* hR_ptCut_eta;
      TH1D* hR_ptCut_phi;
      TH1D* hR_ptCut_ptError;
      TH1D* hR_ptCut_dPtOverPt;
      TH1D* hR_ptCut_etaError;
      TH1D* hR_ptCut_phiError;
      TH1D* hR_ptCut_outerP;
      TH1D* hR_ptCut_outerPt;
      TH1D* hR_ptCut_outerEta;
      TH1D* hR_ptCut_outerPhi;
      TH1D* hR_ptCut_outerPOverP;
      TH1D* hR_ptCut_outerPtOverPt;
      TH1D* hR_ptCut_chi2;
      TH1D* hR_ptCut_ndof;
      TH1D* hR_ptCut_normalizedChi2;
      TH1D* hR_ptCut_charge;
      TH1D* hR_ptCut_dxy;
      TH1D* hR_ptCut_dz;
      TH1D* hR_ptCut_dxyError;
      TH1D* hR_ptCut_dzError;
      TH1D* hR_ptCut_numberOfValidHits;
      TH1D* hR_ptCut_numberOfLostHits;
      TH1D* hR_ptCut_validFraction;

      TH1D* hT_size;
      TH1D* hT_p;
      TH1D* hT_pt;
      TH1D* hT_eta;
      TH1D* hT_phi;
      TH1D* hT_ptError;
      TH1D* hT_dPtOverPt;
      TH1D* hT_etaError;
      TH1D* hT_phiError;
      TH1D* hT_outerP;
      TH1D* hT_outerPt;
      TH1D* hT_outerEta;
      TH1D* hT_outerPhi;
      TH1D* hT_outerPOverP;
      TH1D* hT_outerPtOverPt;
      TH1D* hT_chi2;
      TH1D* hT_ndof;
      TH1D* hT_normalizedChi2;
      TH1D* hT_charge;
      TH1D* hT_dxy;
      TH1D* hT_dz;
      TH1D* hT_dxyError;
      TH1D* hT_dzError;
      TH1D* hT_numberOfValidHits;
      TH1D* hT_numberOfLostHits;
      TH1D* hT_validFraction;

      TH1D* hT_ptCut_size;
      TH1D* hT_ptCut_p;
      TH1D* hT_ptCut_pt;
      TH1D* hT_ptCut_eta;
      TH1D* hT_ptCut_phi;
      TH1D* hT_ptCut_ptError;
      TH1D* hT_ptCut_dPtOverPt;
      TH1D* hT_ptCut_etaError;
      TH1D* hT_ptCut_phiError;
      TH1D* hT_ptCut_outerP;
      TH1D* hT_ptCut_outerPt;
      TH1D* hT_ptCut_outerEta;
      TH1D* hT_ptCut_outerPhi;
      TH1D* hT_ptCut_outerPOverP;
      TH1D* hT_ptCut_outerPtOverPt;
      TH1D* hT_ptCut_chi2;
      TH1D* hT_ptCut_ndof;
      TH1D* hT_ptCut_normalizedChi2;
      TH1D* hT_ptCut_charge;
      TH1D* hT_ptCut_dxy;
      TH1D* hT_ptCut_dz;
      TH1D* hT_ptCut_dxyError;
      TH1D* hT_ptCut_dzError;
      TH1D* hT_ptCut_numberOfValidHits;
      TH1D* hT_ptCut_numberOfLostHits;
      TH1D* hT_ptCut_validFraction;

      TH1D* hMatch_size;
      TH1D* hMatch_dr;
      TH1D* hMatch_nSharedHits;
      TH1D* hMatch_nNonSharedHitsRef;
      TH1D* hMatch_nNonSharedHitsTest;
      TH1D* hMatch_diff_p;
      TH1D* hMatch_diff_pt;
      TH1D* hMatch_diff_eta;
      TH1D* hMatch_diff_phi;
      TH1D* hMatch_diff_ptError;
      TH1D* hMatch_diff_dPtOverPt;
      TH1D* hMatch_diff_etaError;
      TH1D* hMatch_diff_phiError;
      TH1D* hMatch_diff_outerP;
      TH1D* hMatch_diff_outerPt;
      TH1D* hMatch_diff_outerEta;
      TH1D* hMatch_diff_outerPhi;
      TH1D* hMatch_diff_outerPOverP;
      TH1D* hMatch_diff_outerPtOverPt;
      TH1D* hMatch_diff_chi2;
      TH1D* hMatch_diff_ndof;
      TH1D* hMatch_diff_normalizedChi2;
      TH1D* hMatch_diff_charge;
      TH1D* hMatch_diff_dxy;
      TH1D* hMatch_diff_dz;
      TH1D* hMatch_diff_dxyError;
      TH1D* hMatch_diff_dzError;
      TH1D* hMatch_diff_numberOfValidHits;
      TH1D* hMatch_diff_numberOfLostHits;
      TH1D* hMatch_diff_validFraction;

      TH1D* hMatch_pull_p;
      TH1D* hMatch_pull_pt;
      TH1D* hMatch_pull_eta;
      TH1D* hMatch_pull_ptError;
      TH1D* hMatch_pull_dPtOverPt;
      TH1D* hMatch_pull_etaError;
      TH1D* hMatch_pull_phiError;
      TH1D* hMatch_pull_outerP;
      TH1D* hMatch_pull_outerPt;
      TH1D* hMatch_pull_outerEta;
      TH1D* hMatch_pull_outerPOverP;
      TH1D* hMatch_pull_outerPtOverPt;
      TH1D* hMatch_pull_chi2;
      TH1D* hMatch_pull_ndof;
      TH1D* hMatch_pull_normalizedChi2;
      TH1D* hMatch_pull_dxy;
      TH1D* hMatch_pull_dz;
      TH1D* hMatch_pull_dxyError;
      TH1D* hMatch_pull_dzError;
      TH1D* hMatch_pull_numberOfValidHits;
      TH1D* hMatch_pull_numberOfLostHits;
      TH1D* hMatch_pull_validFraction;

      TH2D* hMatch2d_dp_dr;
      TH2D* hMatch2d_dpt_dr;
      TH2D* hMatch2d_nSharedHits_dr;
      TH2D* hMatch2d_dp_dz;
      TH2D* hMatch2d_dxy_dz;

      TH1D* hR_matched_p;
      TH1D* hR_matched_pt;
      TH1D* hR_matched_eta;
      TH1D* hR_matched_phi;
      TH1D* hR_matched_ptError;
      TH1D* hR_matched_dPtOverPt;
      TH1D* hR_matched_etaError;
      TH1D* hR_matched_phiError;
      TH1D* hR_matched_outerP;
      TH1D* hR_matched_outerPt;
      TH1D* hR_matched_outerEta;
      TH1D* hR_matched_outerPhi;
      TH1D* hR_matched_outerPOverP;
      TH1D* hR_matched_outerPtOverPt;
      TH1D* hR_matched_chi2;
      TH1D* hR_matched_ndof;
      TH1D* hR_matched_normalizedChi2;
      TH1D* hR_matched_charge;
      TH1D* hR_matched_dxy;
      TH1D* hR_matched_dz;
      TH1D* hR_matched_dxyError;
      TH1D* hR_matched_dzError;
      TH1D* hR_matched_numberOfValidHits;
      TH1D* hR_matched_numberOfLostHits;
      TH1D* hR_matched_validFraction;

      TH1D* hT_matched_p;
      TH1D* hT_matched_pt;
      TH1D* hT_matched_eta;
      TH1D* hT_matched_phi;
      TH1D* hT_matched_ptError;
      TH1D* hT_matched_dPtOverPt;
      TH1D* hT_matched_etaError;
      TH1D* hT_matched_phiError;
      TH1D* hT_matched_outerP;
      TH1D* hT_matched_outerPt;
      TH1D* hT_matched_outerEta;
      TH1D* hT_matched_outerPhi;
      TH1D* hT_matched_outerPOverP;
      TH1D* hT_matched_outerPtOverPt;
      TH1D* hT_matched_chi2;
      TH1D* hT_matched_ndof;
      TH1D* hT_matched_normalizedChi2;
      TH1D* hT_matched_charge;
      TH1D* hT_matched_dxy;
      TH1D* hT_matched_dz;
      TH1D* hT_matched_dxyError;
      TH1D* hT_matched_dzError;
      TH1D* hT_matched_numberOfValidHits;
      TH1D* hT_matched_numberOfLostHits;
      TH1D* hT_matched_validFraction;

      TH1D* hMatch_ptCut_size;
      TH1D* hMatch_ptCut_dr;
      TH1D* hMatch_ptCut_nSharedHits;
      TH1D* hMatch_ptCut_nNonSharedHitsRef;
      TH1D* hMatch_ptCut_nNonSharedHitsTest;
      TH1D* hMatch_ptCut_diff_p;
      TH1D* hMatch_ptCut_diff_pt;
      TH1D* hMatch_ptCut_diff_eta;
      TH1D* hMatch_ptCut_diff_phi;
      TH1D* hMatch_ptCut_diff_ptError;
      TH1D* hMatch_ptCut_diff_dPtOverPt;
      TH1D* hMatch_ptCut_diff_etaError;
      TH1D* hMatch_ptCut_diff_phiError;
      TH1D* hMatch_ptCut_diff_outerP;
      TH1D* hMatch_ptCut_diff_outerPt;
      TH1D* hMatch_ptCut_diff_outerEta;
      TH1D* hMatch_ptCut_diff_outerPhi;
      TH1D* hMatch_ptCut_diff_outerPOverP;
      TH1D* hMatch_ptCut_diff_outerPtOverPt;
      TH1D* hMatch_ptCut_diff_chi2;
      TH1D* hMatch_ptCut_diff_ndof;
      TH1D* hMatch_ptCut_diff_normalizedChi2;
      TH1D* hMatch_ptCut_diff_charge;
      TH1D* hMatch_ptCut_diff_dxy;
      TH1D* hMatch_ptCut_diff_dz;
      TH1D* hMatch_ptCut_diff_dxyError;
      TH1D* hMatch_ptCut_diff_dzError;
      TH1D* hMatch_ptCut_diff_numberOfValidHits;
      TH1D* hMatch_ptCut_diff_numberOfLostHits;
      TH1D* hMatch_ptCut_diff_validFraction;

      TH1D* hMatch_ptCut_pull_p;
      TH1D* hMatch_ptCut_pull_pt;
      TH1D* hMatch_ptCut_pull_eta;
      TH1D* hMatch_ptCut_pull_phi;
      TH1D* hMatch_ptCut_pull_ptError;
      TH1D* hMatch_ptCut_pull_dPtOverPt;
      TH1D* hMatch_ptCut_pull_etaError;
      TH1D* hMatch_ptCut_pull_phiError;
      TH1D* hMatch_ptCut_pull_outerP;
      TH1D* hMatch_ptCut_pull_outerPt;
      TH1D* hMatch_ptCut_pull_outerEta;
      TH1D* hMatch_ptCut_pull_outerPhi;
      TH1D* hMatch_ptCut_pull_outerPOverP;
      TH1D* hMatch_ptCut_pull_outerPtOverPt;
      TH1D* hMatch_ptCut_pull_chi2;
      TH1D* hMatch_ptCut_pull_ndof;
      TH1D* hMatch_ptCut_pull_normalizedChi2;
      TH1D* hMatch_ptCut_pull_dxy;
      TH1D* hMatch_ptCut_pull_dz;
      TH1D* hMatch_ptCut_pull_dxyError;
      TH1D* hMatch_ptCut_pull_dzError;
      TH1D* hMatch_ptCut_pull_numberOfValidHits;
      TH1D* hMatch_ptCut_pull_numberOfLostHits;
      TH1D* hMatch_ptCut_pull_validFraction;

      // Deta and Dphi
      TH1D* hR_dEta_mapSize;
      TH1D* hR_dEta_mapSizeMatch;
      TH1D* hR_dEta_dEta;
      TH1D* hR_dEta_nPassCut;
      TH1D* hR_dEta_recoEcalCandsSize;
      TH1D* hR_dEta_recoEcalCandsSizeMatch;
      TH1D* hR_dEta_dEta_trgdEcalCands;
      TH1D* hR_dEta_nPassCut_trgdEcalCands;

      TH1D* hR_dPhi_mapSize;
      TH1D* hR_dPhi_mapSizeMatch;
      TH1D* hR_dPhi_dPhi;
      TH1D* hR_dPhi_nPassCut;
      TH1D* hR_dPhi_recoEcalCandsSize;
      TH1D* hR_dPhi_recoEcalCandsSizeMatch;
      TH1D* hR_dPhi_dPhi_trgdEcalCands;
      TH1D* hR_dPhi_nPassCut_trgdEcalCands;

      TH1D* hT_dEta_mapSize;
      TH1D* hT_dEta_mapSizeMatch;
      TH1D* hT_dEta_dEta;
      TH1D* hT_dEta_nPassCut;
      TH1D* hT_dEta_recoEcalCandsSize;
      TH1D* hT_dEta_recoEcalCandsSizeMatch;
      TH1D* hT_dEta_dEta_trgdEcalCands;
      TH1D* hT_dEta_nPassCut_trgdEcalCands;

      TH1D* hT_dPhi_mapSize;
      TH1D* hT_dPhi_mapSizeMatch;
      TH1D* hT_dPhi_dPhi;
      TH1D* hT_dPhi_nPassCut;
      TH1D* hT_dPhi_recoEcalCandsSize;
      TH1D* hT_dPhi_recoEcalCandsSizeMatch;
      TH1D* hT_dPhi_dPhi_trgdEcalCands;
      TH1D* hT_dPhi_nPassCut_trgdEcalCands;

      //TTree* tree;
      //unsigned int trkR_size;
      //float pt[50];
};


// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HltGsfTrkAnalyser::HltGsfTrkAnalyser(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   referenceTrackCollTag_ = iConfig.getParameter<edm::InputTag>("referenceTrackCollTag");
   testTrackCollTag_ = iConfig.getParameter<edm::InputTag>("testTrackCollTag");
   referenceTrackVarsTag_ = iConfig.getParameter<edm::InputTag>("referenceTrackVarsTag");
   testTrackVarsTag_ = iConfig.getParameter<edm::InputTag>("testTrackVarsTag");
   ecalCandTagForDEtaR_ = iConfig.getParameter<edm::InputTag>("referenceEcalCandTagForDEta");
   ecalCandTagForDEtaT_ = iConfig.getParameter<edm::InputTag>("testEcalCandTagForDEta");
   ecalCandTagForDPhiR_ = iConfig.getParameter<edm::InputTag>("referenceEcalCandTagForDPhi");
   ecalCandTagForDPhiT_ = iConfig.getParameter<edm::InputTag>("testEcalCandTagForDPhi");
   dEtaCutR_ = iConfig.getUntrackedParameter<double>("referenceDEtaCut", 100.);
   dPhiCutR_ = iConfig.getUntrackedParameter<double>("referenceDPhiCut", 100.);
   dEtaCutT_ = iConfig.getUntrackedParameter<double>("testDEtaCut", 100.);
   dPhiCutT_ = iConfig.getUntrackedParameter<double>("testDPhiCut", 100.);
   matchMethod_ = iConfig.getUntrackedParameter<unsigned int>("matchMethod", 0);
   discrCut_ = iConfig.getUntrackedParameter<double>("discrCutForMatch", 10.);
   largerThan_ = iConfig.getUntrackedParameter<bool>("largerThan", false);
   minPtR_ = iConfig.getUntrackedParameter<double>("referenceMinPt", 0.);
   minPtT_ = iConfig.getUntrackedParameter<double>("testMinPt", 0.);

   edm::ParameterSet defTrgPSet;
   defTrgPSet.addUntrackedParameter<std::string>("triggerName", "HLTriggerFinalPath");
   defTrgPSet.addUntrackedParameter<bool>("invertBit", 0);
   refTrgPSet_ = iConfig.getUntrackedParameter<edm::ParameterSet>("referenceTrigger", defTrgPSet);
   testTrgPSet_ = iConfig.getUntrackedParameter<edm::ParameterSet>("testTrigger", defTrgPSet);
   refTrgName_ = refTrgPSet_.getUntrackedParameter<std::string>("triggerName", "HLTriggerFinalPath");
   refTrgInv_ = refTrgPSet_.getUntrackedParameter<bool>("invertBit", 0);
   testTrgName_ = testTrgPSet_.getUntrackedParameter<std::string>("triggerName", "HLTriggerFinalPath");
   testTrgInv_ = testTrgPSet_.getUntrackedParameter<bool>("invertBit", 0);

   edm::Service<TFileService> fs;
   hR_size = fs->make<TH1D>("hR_size", "GsfTrack size", 50, 0., 50.);
   hR_p = fs->make<TH1D>("hR_p", "GsfTrack p", 100, 0., 250.);
   hR_pt = fs->make<TH1D>("hR_pt", "GsfTrack p_{T}", 100, 0., 250.);
   hR_eta = fs->make<TH1D>("hR_eta", "GsfTrack #eta", 100, -3., 3.);
   hR_phi = fs->make<TH1D>("hR_phi", "GsfTrack #phi", 100, -3.2, 3.2);
   hR_ptError = fs->make<TH1D>("hR_ptError", "GsfTrack p_{T}-error", 100, 0., 15.);
   hR_dPtOverPt = fs->make<TH1D>("hR_dPtOverPt", "GsfTrack dp_{T}/p_{T}", 100, 0., 1.);
   hR_etaError = fs->make<TH1D>("hR_etaError", "GsfTrack #eta-error", 100, 0., 0.1);
   hR_phiError = fs->make<TH1D>("hR_phiError", "GsfTrack #phi-error", 100, 0., 0.1);
   hR_outerP = fs->make<TH1D>("hR_outerP", "GsfTrack p_{out}", 100, 0., 250.);
   hR_outerPt = fs->make<TH1D>("hR_outerPt", "GsfTrack p_{T}^{out}", 100, 0., 250.);
   hR_outerEta = fs->make<TH1D>("hR_outerEta", "GsfTrack #eta_{out}", 100, -3., 3.);
   hR_outerPhi = fs->make<TH1D>("hR_outerPhi", "GsfTrack #phi_{out}", 100, -3.2, 3.2);
   hR_outerPOverP = fs->make<TH1D>("hR_outerPOverP", "GsfTrack p_{out}/p", 100, 0., 2.);
   hR_outerPtOverPt = fs->make<TH1D>("hR_outerPtOverPt", "GsfTrack p_{T}^{out}/p_{T}", 100, 0., 2.);
   hR_chi2 = fs->make<TH1D>("hR_chi2", "GsfTrack #chi^{2}", 100, 0., 100.);
   hR_ndof = fs->make<TH1D>("hR_ndof", "GsfTrack ndof", 45, 0., 45.);
   hR_normalizedChi2 = fs->make<TH1D>("hR_normalizedChi2", "GsfTrack #chi^{2}_{normalized}", 100, 0., 50.);
   hR_charge = fs->make<TH1D>("hR_charge", "GsfTrack charge", 3, -1., 2.);
   hR_dxy = fs->make<TH1D>("hR_dxy", "GsfTrack d_{xy}", 100, -25., 25.);
   hR_dz = fs->make<TH1D>("hR_dz", "GsfTrack d_{z}", 100, -100., 100.);
   hR_dxyError = fs->make<TH1D>("hR_dxyError", "GsfTrack d_{xy}-error", 100, 0., 5.);
   hR_dzError = fs->make<TH1D>("hR_dzError", "GsfTrack d_{z}-error", 100, 0., 10.);
   hR_numberOfValidHits = fs->make<TH1D>("hR_numberOfValidHits", "GsfTrack #hits_{valid}", 25, 0., 25.);
   hR_numberOfLostHits = fs->make<TH1D>("hR_numberOfLostHits", "GsfTrack #hits_{lost}", 5, 0., 5.);
   hR_validFraction = fs->make<TH1D>("hR_validFraction", "GsfTrack valid fraction", 100, 0., 1.);
   hR2d_pt_eta = fs->make<TH2D>("hR2d_pt_eta", "GsfTrack p_{T}:#eta", 500, 0., 500., 121, -3., 3.);
   hR2d_pt_dpt = fs->make<TH2D>("hR2d_pt_dpt", "GsfTrack p_{T}:dp_{T}", 500, 0., 500., 100, 0., 100.);
   hR2d_p_outerP = fs->make<TH2D>("hR2d_p_outerP", "GsfTrack p:p_{out}", 500, 0., 500., 500, 0., 500.);
   hR2d_pt_outerPt = fs->make<TH2D>("hR2d_pt_outerPt", "GsfTrack p_{T}:p_{T}^{out}", 500, 0., 500., 500, 0., 500.);

   hR_ptCut_size = fs->make<TH1D>("hR_ptCut_size", "GsfTrack size", 50, 0., 50.);
   hR_ptCut_p = fs->make<TH1D>("hR_ptCut_p", "GsfTrack p", 100, 0., 250.);
   hR_ptCut_pt = fs->make<TH1D>("hR_ptCut_pt", "GsfTrack p_{T}", 100, 0., 250.);
   hR_ptCut_eta = fs->make<TH1D>("hR_ptCut_eta", "GsfTrack #eta", 100, -3., 3.);
   hR_ptCut_phi = fs->make<TH1D>("hR_ptCut_phi", "GsfTrack #phi", 100, -3.2, 3.2);
   hR_ptCut_ptError = fs->make<TH1D>("hR_ptCut_ptError", "GsfTrack p_{T}-error", 100, 0., 15.);
   hR_ptCut_dPtOverPt = fs->make<TH1D>("hR_ptCut_dPtOverPt", "GsfTrack dp_{T}/p_{T}", 100, 0., 1.);
   hR_ptCut_etaError = fs->make<TH1D>("hR_ptCut_etaError", "GsfTrack #eta-error", 100, 0., 0.1);
   hR_ptCut_phiError = fs->make<TH1D>("hR_ptCut_phiError", "GsfTrack #phi-error", 100, 0., 0.1);
   hR_ptCut_outerP = fs->make<TH1D>("hR_ptCut_outerP", "GsfTrack p_{out}", 100, 0., 250.);
   hR_ptCut_outerPt = fs->make<TH1D>("hR_ptCut_outerPt", "GsfTrack p_{T}^{out}", 100, 0., 250.);
   hR_ptCut_outerEta = fs->make<TH1D>("hR_ptCut_outerEta", "GsfTrack #eta_{out}", 100, -3., 3.);
   hR_ptCut_outerPhi = fs->make<TH1D>("hR_ptCut_outerPhi", "GsfTrack #phi_{out}", 100, -3.2, 3.2);
   hR_ptCut_outerPOverP = fs->make<TH1D>("hR_ptCut_outerPOverP", "GsfTrack p_{out}/p", 100, 0., 2.);
   hR_ptCut_outerPtOverPt = fs->make<TH1D>("hR_ptCut_outerPtOverPt", "GsfTrack p_{T}^{out}/p_{T}", 100, 0., 2.);
   hR_ptCut_chi2 = fs->make<TH1D>("hR_ptCut_chi2", "GsfTrack #chi^{2}", 100, 0., 100.);
   hR_ptCut_ndof = fs->make<TH1D>("hR_ptCut_ndof", "GsfTrack ndof", 45, 0., 45.);
   hR_ptCut_normalizedChi2 = fs->make<TH1D>("hR_ptCut_normalizedChi2", "GsfTrack #chi^{2}_{normalized}", 100, 0., 50.);
   hR_ptCut_charge = fs->make<TH1D>("hR_ptCut_charge", "GsfTrack charge", 3, -1., 2.);
   hR_ptCut_dxy = fs->make<TH1D>("hR_ptCut_dxy", "GsfTrack d_{xy}", 100, -25., 25.);
   hR_ptCut_dz = fs->make<TH1D>("hR_ptCut_dz", "GsfTrack d_{z}", 100, -100., 100.);
   hR_ptCut_dxyError = fs->make<TH1D>("hR_ptCut_dxyError", "GsfTrack d_{xy}-error", 100, 0., 5.);
   hR_ptCut_dzError = fs->make<TH1D>("hR_ptCut_dzError", "GsfTrack d_{z}-error", 100, 0., 10.);
   hR_ptCut_numberOfValidHits = fs->make<TH1D>("hR_ptCut_numberOfValidHits", "GsfTrack #hits_{valid}", 25, 0., 25.);
   hR_ptCut_numberOfLostHits = fs->make<TH1D>("hR_ptCut_numberOfLostHits", "GsfTrack #hits_{lost}", 5, 0., 5.);
   hR_ptCut_validFraction = fs->make<TH1D>("hR_ptCut_validFraction", "GsfTrack valid fraction", 100, 0., 1.);

   hT_size = fs->make<TH1D>("hT_size", "GsfTrack size", 50, 0., 50.);
   hT_p = fs->make<TH1D>("hT_p", "GsfTrack p", 100, 0., 250.);
   hT_pt = fs->make<TH1D>("hT_pt", "GsfTrack p_{T}", 100, 0., 250.);
   hT_eta = fs->make<TH1D>("hT_eta", "GsfTrack #eta", 100, -3., 3.);
   hT_phi = fs->make<TH1D>("hT_phi", "GsfTrack #phi", 100, -3.2, 3.2);
   hT_ptError = fs->make<TH1D>("hT_ptError", "GsfTrack p_{T}-error", 100, 0., 15.);
   hT_dPtOverPt = fs->make<TH1D>("hT_dPtOverPt", "GsfTrack dp_{T}/p_{T}", 100, 0., 1.);
   hT_etaError = fs->make<TH1D>("hT_etaError", "GsfTrack #eta-error", 100, 0., 0.1);
   hT_phiError = fs->make<TH1D>("hT_phiError", "GsfTrack #phi-error", 100, 0., 0.1);
   hT_outerP = fs->make<TH1D>("hT_outerP", "GsfTrack p_{out}", 100, 0., 250.);
   hT_outerPt = fs->make<TH1D>("hT_outerPt", "GsfTrack p_{T}^{out}", 100, 0., 250.);
   hT_outerEta = fs->make<TH1D>("hT_outerEta", "GsfTrack #eta_{out}", 100, -3., 3.);
   hT_outerPhi = fs->make<TH1D>("hT_outerPhi", "GsfTrack #phi_{out}", 100, -3.2, 3.2);
   hT_outerPOverP = fs->make<TH1D>("hT_outerPOverP", "GsfTrack p_{out}/p", 100, 0., 2.);
   hT_outerPtOverPt = fs->make<TH1D>("hT_outerPtOverPt", "GsfTrack p_{T}^{out}/p_{T}", 100, 0., 2.);
   hT_chi2 = fs->make<TH1D>("hT_chi2", "GsfTrack #chi^{2}", 100, 0., 100.);
   hT_ndof = fs->make<TH1D>("hT_ndof", "GsfTrack ndof", 45, 0., 45.);
   hT_normalizedChi2 = fs->make<TH1D>("hT_normalizedChi2", "GsfTrack #chi^{2}_{normalized}", 100, 0., 50.);
   hT_charge = fs->make<TH1D>("hT_charge", "GsfTrack charge", 3, -1., 2.);
   hT_dxy = fs->make<TH1D>("hT_dxy", "GsfTrack d_{xy}", 100, -25., 25.);
   hT_dz = fs->make<TH1D>("hT_dz", "GsfTrack d_{z}", 100, -100., 100.);
   hT_dxyError = fs->make<TH1D>("hT_dxyError", "GsfTrack d_{xy}-error", 100, 0., 5.);
   hT_dzError = fs->make<TH1D>("hT_dzError", "GsfTrack d_{z}-error", 100, 0., 10.);
   hT_numberOfValidHits = fs->make<TH1D>("hT_numberOfValidHits", "GsfTrack #hits_{valid}", 25, 0., 25.);
   hT_numberOfLostHits = fs->make<TH1D>("hT_numberOfLostHits", "GsfTrack #hits_{lost}", 5, 0., 5.);
   hT_validFraction = fs->make<TH1D>("hT_validFraction", "GsfTrack valid fraction", 100, 0., 1.);

   hT_ptCut_size = fs->make<TH1D>("hT_ptCut_size", "GsfTrack size", 50, 0., 50.);
   hT_ptCut_p = fs->make<TH1D>("hT_ptCut_p", "GsfTrack p", 100, 0., 250.);
   hT_ptCut_pt = fs->make<TH1D>("hT_ptCut_pt", "GsfTrack p_{T}", 100, 0., 250.);
   hT_ptCut_eta = fs->make<TH1D>("hT_ptCut_eta", "GsfTrack #eta", 100, -3., 3.);
   hT_ptCut_phi = fs->make<TH1D>("hT_ptCut_phi", "GsfTrack #phi", 100, -3.2, 3.2);
   hT_ptCut_ptError = fs->make<TH1D>("hT_ptCut_ptError", "GsfTrack p_{T}-error", 100, 0., 15.);
   hT_ptCut_dPtOverPt = fs->make<TH1D>("hT_ptCut_dPtOverPt", "GsfTrack dp_{T}/p_{T}", 100, 0., 1.);
   hT_ptCut_etaError = fs->make<TH1D>("hT_ptCut_etaError", "GsfTrack #eta-error", 100, 0., 0.1);
   hT_ptCut_phiError = fs->make<TH1D>("hT_ptCut_phiError", "GsfTrack #phi-error", 100, 0., 0.1);
   hT_ptCut_outerP = fs->make<TH1D>("hT_ptCut_outerP", "GsfTrack p_{out}", 100, 0., 250.);
   hT_ptCut_outerPt = fs->make<TH1D>("hT_ptCut_outerPt", "GsfTrack p_{T}^{out}", 100, 0., 250.);
   hT_ptCut_outerEta = fs->make<TH1D>("hT_ptCut_outerEta", "GsfTrack #eta_{out}", 100, -3., 3.);
   hT_ptCut_outerPhi = fs->make<TH1D>("hT_ptCut_outerPhi", "GsfTrack #phi_{out}", 100, -3.2, 3.2);
   hT_ptCut_outerPOverP = fs->make<TH1D>("hT_ptCut_outerPOverP", "GsfTrack p_{out}/p", 100, 0., 2.);
   hT_ptCut_outerPtOverPt = fs->make<TH1D>("hT_ptCut_outerPtOverPt", "GsfTrack p_{T}^{out}/p_{T}", 100, 0., 2.);
   hT_ptCut_chi2 = fs->make<TH1D>("hT_ptCut_chi2", "GsfTrack #chi^{2}", 100, 0., 100.);
   hT_ptCut_ndof = fs->make<TH1D>("hT_ptCut_ndof", "GsfTrack ndof", 45, 0., 45.);
   hT_ptCut_normalizedChi2 = fs->make<TH1D>("hT_ptCut_normalizedChi2", "GsfTrack #chi^{2}_{normalized}", 100, 0., 50.);
   hT_ptCut_charge = fs->make<TH1D>("hT_ptCut_charge", "GsfTrack charge", 3, -1., 2.);
   hT_ptCut_dxy = fs->make<TH1D>("hT_ptCut_dxy", "GsfTrack d_{xy}", 100, -25., 25.);
   hT_ptCut_dz = fs->make<TH1D>("hT_ptCut_dz", "GsfTrack d_{z}", 100, -100., 100.);
   hT_ptCut_dxyError = fs->make<TH1D>("hT_ptCut_dxyError", "GsfTrack d_{xy}-error", 100, 0., 5.);
   hT_ptCut_dzError = fs->make<TH1D>("hT_ptCut_dzError", "GsfTrack d_{z}-error", 100, 0., 10.);
   hT_ptCut_numberOfValidHits = fs->make<TH1D>("hT_ptCut_numberOfValidHits", "GsfTrack #hits_{valid}", 25, 0., 25.);
   hT_ptCut_numberOfLostHits = fs->make<TH1D>("hT_ptCut_numberOfLostHits", "GsfTrack #hits_{lost}", 5, 0., 5.);
   hT_ptCut_validFraction = fs->make<TH1D>("hT_ptCut_validFraction", "GsfTrack valid fraction", 100, 0., 1.);

   hMatch_size = fs->make<TH1D>("hMatch_size", "GsfTrack_{matched}: size", 35, 0., 35.);
   hMatch_dr = fs->make<TH1D>("hMatch_dr", "GsfTrack_{matched}: #DeltaR", 100, 0., 0.5);
   hMatch_nSharedHits = fs->make<TH1D>("hMatch_nSharedHits", "GsfTrack_{matched}: # hits_{shared}", 30, 0., 30);
   hMatch_nNonSharedHitsRef = fs->make<TH1D>("hMatch_nNonSharedHitsRef", "GsfTrack_{matched}^{ref}: # hits_{not shared}", 10, 0., 10.);
   hMatch_nNonSharedHitsTest = fs->make<TH1D>("hMatch_nNonSharedHitsTest", "GsfTrack_{matched}^{test}: # hits_{not shared}", 10, 0., 10.);
   hMatch_diff_p = fs->make<TH1D>("hMatch_diff_p", "GsfTrack_{matched}: p_{test} - p_{ref}", 100, -250., 250.);
   hMatch_diff_pt = fs->make<TH1D>("hMatch_diff_pt", "GsfTrack_{matched}: p_{T}^{test} - p_{T}^{ref}", 100, -250., 250.);
   hMatch_diff_eta = fs->make<TH1D>("hMatch_diff_eta", "GsfTrack_{matched}: #eta_{test} -  #eta_{ref}", 100, -1., 1.);
   hMatch_diff_phi = fs->make<TH1D>("hMatch_diff_phi", "GsfTrack_{matched}: #phi_{test} - #phi_{ref}", 100, -1., 1.);
   hMatch_diff_ptError = fs->make<TH1D>("hMatch_diff_ptError", "GsfTrack_{matched}: p_{T}-error_{test} - p_{T}-error_{ref}", 100, -20., 20.);
   hMatch_diff_dPtOverPt = fs->make<TH1D>("hMatch_diff_dPtOverPt", "GsfTrack_{matched}: dp_{T}/p_{T}^{test} - dp_{T}/p_{T}^{ref}", 100, -1., 1.);
   hMatch_diff_etaError = fs->make<TH1D>("hMatch_diff_etaError", "GsfTrack_{matched}: #eta-error_{test} - #eta-error_{ref}", 100, -0.1, 0.1);
   hMatch_diff_phiError = fs->make<TH1D>("hMatch_diff_phiError", "GsfTrack_{matched}: #phi-error_{test} - #phi-error_{ref}", 100, -0.1, 0.1);
   hMatch_diff_outerP = fs->make<TH1D>("hMatch_diff_outerP", "GsfTrack_{matched}: p_{out}^{test} - p_{out}^{ref}", 100, -250., 250.);
   hMatch_diff_outerPt = fs->make<TH1D>("hMatch_diff_outerPt", "GsfTrack_{matched}: (p_{T}^{out})_{test} - (p_{T}^{out})_{ref}", 100, -250., 250.);
   hMatch_diff_outerEta = fs->make<TH1D>("hMatch_diff_outerEta", "GsfTrack_{matched}: #eta_{out}^{test} - #eta_{out}^{ref}", 100, -1., 1.);
   hMatch_diff_outerPhi = fs->make<TH1D>("hMatch_diff_outerPhi", "GsfTrack_{matched}: #phi_{out}^{test} - #phi_{out}^{ref}", 100, -1., 1.);
   hMatch_diff_outerPOverP = fs->make<TH1D>("hMatch_diff_outerPOverP", "GsfTrack_{matched}: (p_{out}/p)_{test} - (p_{out}/p)_{ref}", 100, -4., 4.);
   hMatch_diff_outerPtOverPt = fs->make<TH1D>("hMatch_diff_outerPtOverPt", "GsfTrack_{matched}: (p_{T}^{out}/p_{T})_{test} - (p_{T}^{out}/p_{T})_{ref}", 100, -4., 4.);
   hMatch_diff_chi2 = fs->make<TH1D>("hMatch_diff_chi2", "GsfTrack_{matched}: #chi^{2}_{test} - #chi^{2}_{ref}", 100, -150., 150.);
   hMatch_diff_ndof = fs->make<TH1D>("hMatch_diff_ndof", "GsfTrack_{matched}: ndof_{test} - ndof_{ref}", 40, -20., 20.);
   hMatch_diff_normalizedChi2 = fs->make<TH1D>("hMatch_diff_normalizedChi2", "GsfTrack_{matched}: (#chi^{2}_{normalized})_{test} - (#chi^{2}_{normalized})_{ref}", 100, -50., 50.);
   hMatch_diff_charge = fs->make<TH1D>("hMatch_diff_charge", "GsfTrack_{matched}: charge_{test} - charge_{ref}", 5, -2., 3.);
   hMatch_diff_dxy = fs->make<TH1D>("hMatch_diff_dxy", "GsfTrack_{matched}: d_{xy}^{test} - d_{xy}^{ref}", 100, -10., 10.);
   hMatch_diff_dz = fs->make<TH1D>("hMatch_diff_dz", "GsfTrack_{matched}: d_{z}^{test} - d_{z}^{ref}", 100, -50., 50.);
   hMatch_diff_dxyError = fs->make<TH1D>("hMatch_diff_dxyError", "GsfTrack_{matched}: d_{xy}-error_{test} - d_{xy}-error_{ref}", 100, -5., 5.);
   hMatch_diff_dzError = fs->make<TH1D>("hMatch_diff_dzError", "GsfTrack_{matched}: d_{z}-error_{test} - d_{z}-error_{ref}", 100, -10., 10.);
   hMatch_diff_numberOfValidHits = fs->make<TH1D>("hMatch_diff_numberOfValidHits", "GsfTrack_{matched}: #hits_{valid}^{test} - #hits_{valid}^{ref}", 20, -10., 10.);
   hMatch_diff_numberOfLostHits = fs->make<TH1D>("hMatch_diff_numberOfLostHits", "GsfTrack_{matched}: #hits_{lost}^{test} - #hits_{lost}^{ref}", 10, -5., 5.);
   hMatch_diff_validFraction = fs->make<TH1D>("hMatch_diff_validFraction", "GsfTrack_{matched}: valid fraction_{test} - valid fraction_{ref}", 100, -1., 1.);

   hMatch_pull_p = fs->make<TH1D>("hMatch_pull_p", "GsfTrack_{matched}: |p_{test} - p_{ref}|/p_{ref}", 100, 0., 30.);
   hMatch_pull_pt = fs->make<TH1D>("hMatch_pull_pt", "GsfTrack_{matched}: |p_{T}^{test} - p_{T}^{ref}|/p_{T}^{ref}", 100, 0., 20.);
   hMatch_pull_eta = fs->make<TH1D>("hMatch_pull_eta", "GsfTrack_{matched}: |#eta_{test} -  #eta_{ref}|/#eta_{ref}", 100, 0., 1.);
   hMatch_pull_ptError = fs->make<TH1D>("hMatch_pull_ptError", "GsfTrack_{matched}: |p_{T}-error_{test} - p_{T}-error_{ref}|/p_{T}-error_{ref}", 100, 0., 20.);
   hMatch_pull_dPtOverPt = fs->make<TH1D>("hMatch_pull_dPtOverPt", "GsfTrack_{matched}: |dp_{T}/p_{T}^{test} - dp_{T}/p_{T}^{ref}|/(dp_{T}/p_{T}^{ref})", 100, 0., 3.);
   hMatch_pull_etaError = fs->make<TH1D>("hMatch_pull_etaError", "GsfTrack_{matched}: |#eta-error_{test} - #eta-error_{ref}|/#eta-error_{ref}", 100, 0., 3.);
   hMatch_pull_phiError = fs->make<TH1D>("hMatch_pull_phiError", "GsfTrack_{matched}: |#phi-error_{test} - #phi-error_{ref}|/#phi-error_{ref}", 100, 0., 3.);
   hMatch_pull_outerP = fs->make<TH1D>("hMatch_pull_outerP", "GsfTrack_{matched}: |p_{out}^{test} - p_{out}^{ref}|/p_{out}^{ref}", 100, 0., 25.);
   hMatch_pull_outerPt = fs->make<TH1D>("hMatch_pull_outerPt", "GsfTrack_{matched}: |(p_{T}^{out})_{test} - (p_{T}^{out})_{ref}|/(p_{T}^{out})_{ref}", 100, 0., 25.);
   hMatch_pull_outerEta = fs->make<TH1D>("hMatch_pull_outerEta", "GsfTrack_{matched}: |#eta_{out}^{test} - #eta_{out}^{ref}|/#eta_{out}^{ref}", 100, 0., 1.);
   hMatch_pull_outerPOverP = fs->make<TH1D>("hMatch_pull_outerPOverP", "GsfTrack_{matched}: |(p_{out}/p)_{test} - (p_{out}/p)_{ref}|/(p_{out}/p)_{ref}", 100, 0., 8.);
   hMatch_pull_outerPtOverPt = fs->make<TH1D>("hMatch_pull_outerPtOverPt", "GsfTrack_{matched}: |(p_{T}^{out}/p_{T})_{test} - (p_{T}^{out}/p_{T})_{ref}|/(p_{T}^{out}/p_{T})_{ref}", 100, 0., 8.);
   hMatch_pull_chi2 = fs->make<TH1D>("hMatch_pull_chi2", "GsfTrack_{matched}: |#chi^{2}_{test} - #chi^{2}_{ref}|/#chi^{2}_{ref}", 100, 0., 20.);
   hMatch_pull_ndof = fs->make<TH1D>("hMatch_pull_ndof", "GsfTrack_{matched}: |ndof_{test} - ndof_{ref}|/ndof_{ref}", 100, 0., 3.);
   hMatch_pull_normalizedChi2 = fs->make<TH1D>("hMatch_pull_normalizedChi2", "GsfTrack_{matched}: |(#chi^{2}_{normalized})_{test} - (#chi^{2}_{normalized})_{ref}|/(#chi^{2}_{normalized})_{ref}", 100, 0., 10.);
   hMatch_pull_dxy = fs->make<TH1D>("hMatch_pull_dxy", "GsfTrack_{matched}: |d_{xy}^{test} - d_{xy}^{ref}|/|d_{xy}^{ref}|", 100, 0., 10.);
   hMatch_pull_dz = fs->make<TH1D>("hMatch_pull_dz", "GsfTrack_{matched}: |d_{z}^{test} - d_{z}^{ref}|/|d_{z}^{ref}|", 100, 0., 10.);
   hMatch_pull_dxyError = fs->make<TH1D>("hMatch_pull_dxyError", "GsfTrack_{matched}: |d_{xy}-error_{test} - d_{xy}-error_{ref}|/d_{xy}-error_{ref}", 100, 0., 4.);
   hMatch_pull_dzError = fs->make<TH1D>("hMatch_pull_dzError", "GsfTrack_{matched}: |d_{z}-error_{test} - d_{z}-error_{ref}|/d_{z}-error_{ref}", 100, 0., 4.);
   hMatch_pull_numberOfValidHits = fs->make<TH1D>("hMatch_pull_numberOfValidHits", "GsfTrack_{matched}: (#hits_{valid}^{test} - #hits_{valid}^{ref}|/#hits_{valid}^{ref}", 100, 0., 5.);
   hMatch_pull_numberOfLostHits = fs->make<TH1D>("hMatch_pull_numberOfLostHits", "GsfTrack_{matched}: |#hits_{lost}^{test} - #hits_{lost}^{ref}|/#hits_{lost}^{ref}", 100, 0., 5.);
   hMatch_pull_validFraction = fs->make<TH1D>("hMatch_pull_validFraction", "GsfTrack_{matched}: |valid fraction_{test} - valid fraction_{ref}|/valid fraction_{ref}", 100, 0., 1.);

   hMatch2d_dp_dr = fs->make<TH2D>("hMatch2d_dp_dr", "GsfTrack_{matched}: #Deltap:#DeltaR", 500, -250., 250., 100, 0., 0.5);
   hMatch2d_dpt_dr = fs->make<TH2D>("hMatch2d_dpt_dr", "GsfTrack_{matched}: #Deltap_{T}:#DeltaR", 500, -250., 250., 100, 0., 0.5);
   hMatch2d_nSharedHits_dr = fs->make<TH2D>("hMatch2d_nSharedHits_dr", "GsfTrack_{matched}: # hits_{shared}:#DeltaR", 30, 0., 30., 100, 0., 0.5);
   hMatch2d_dp_dz = fs->make<TH2D>("hMatch2d_dp_dz", "GsfTrack_{matched}: #Deltap:d_{z}", 500, -250., 250., 200, -100., 100.);
   hMatch2d_dxy_dz = fs->make<TH2D>("hMatch2d_dxy_dz", "GsfTrack_{matched}: d_{xy}:d_{z}", 101, -25., 25., 200, -100., 100.);

   hR_matched_p = fs->make<TH1D>("hR_matched_p", "GsfTrack p", 100, 0., 250.);
   hR_matched_pt = fs->make<TH1D>("hR_matched_pt", "GsfTrack p_{T}", 100, 0., 250.);
   hR_matched_eta = fs->make<TH1D>("hR_matched_eta", "GsfTrack #eta", 100, -3., 3.);
   hR_matched_phi = fs->make<TH1D>("hR_matched_phi", "GsfTrack #phi", 100, -3.2, 3.2);
   hR_matched_ptError = fs->make<TH1D>("hR_matched_ptError", "GsfTrack p_{T}-error", 100, 0., 15.);
   hR_matched_dPtOverPt = fs->make<TH1D>("hR_matched_dPtOverPt", "GsfTrack dp_{T}/p_{T}", 100, 0., 1.);
   hR_matched_etaError = fs->make<TH1D>("hR_matched_etaError", "GsfTrack #eta-error", 100, 0., 0.1);
   hR_matched_phiError = fs->make<TH1D>("hR_matched_phiError", "GsfTrack #phi-error", 100, 0., 0.1);
   hR_matched_outerP = fs->make<TH1D>("hR_matched_outerP", "GsfTrack p_{out}", 100, 0., 250.);
   hR_matched_outerPt = fs->make<TH1D>("hR_matched_outerPt", "GsfTrack p_{T}^{out}", 100, 0., 250.);
   hR_matched_outerEta = fs->make<TH1D>("hR_matched_outerEta", "GsfTrack #eta_{out}", 100, -3., 3.);
   hR_matched_outerPhi = fs->make<TH1D>("hR_matched_outerPhi", "GsfTrack #phi_{out}", 100, -3.2, 3.2);
   hR_matched_outerPOverP = fs->make<TH1D>("hR_matched_outerPOverP", "GsfTrack p_{out}/p", 100, 0., 2.);
   hR_matched_outerPtOverPt = fs->make<TH1D>("hR_matched_outerPtOverPt", "GsfTrack p_{T}^{out}/p_{T}", 100, 0., 2.);
   hR_matched_chi2 = fs->make<TH1D>("hR_matched_chi2", "GsfTrack #chi^{2}", 100, 0., 100.);
   hR_matched_ndof = fs->make<TH1D>("hR_matched_ndof", "GsfTrack ndof", 45, 0., 45.);
   hR_matched_normalizedChi2 = fs->make<TH1D>("hR_matched_normalizedChi2", "GsfTrack #chi^{2}_{normalized}", 100, 0., 50.);
   hR_matched_charge = fs->make<TH1D>("hR_matched_charge", "GsfTrack charge", 3, -1., 2.);
   hR_matched_dxy = fs->make<TH1D>("hR_matched_dxy", "GsfTrack d_{xy}", 100, -25., 25.);
   hR_matched_dz = fs->make<TH1D>("hR_matched_dz", "GsfTrack d_{z}", 100, -100., 100.);
   hR_matched_dxyError = fs->make<TH1D>("hR_matched_dxyError", "GsfTrack d_{xy}-error", 100, 0., 5.);
   hR_matched_dzError = fs->make<TH1D>("hR_matched_dzError", "GsfTrack d_{z}-error", 100, 0., 10.);
   hR_matched_numberOfValidHits = fs->make<TH1D>("hR_matched_numberOfValidHits", "GsfTrack #hits_{valid}", 25, 0., 25.);
   hR_matched_numberOfLostHits = fs->make<TH1D>("hR_matched_numberOfLostHits", "GsfTrack #hits_{lost}", 5, 0., 5.);
   hR_matched_validFraction = fs->make<TH1D>("hR_matched_validFraction", "GsfTrack valid fraction", 100, 0., 1.);

   hT_matched_p = fs->make<TH1D>("hT_matched_p", "GsfTrack p", 100, 0., 250.);
   hT_matched_pt = fs->make<TH1D>("hT_matched_pt", "GsfTrack p_{T}", 100, 0., 250.);
   hT_matched_eta = fs->make<TH1D>("hT_matched_eta", "GsfTrack #eta", 100, -3., 3.);
   hT_matched_phi = fs->make<TH1D>("hT_matched_phi", "GsfTrack #phi", 100, -3.2, 3.2);
   hT_matched_ptError = fs->make<TH1D>("hT_matched_ptError", "GsfTrack p_{T}-error", 100, 0., 15.);
   hT_matched_dPtOverPt = fs->make<TH1D>("hT_matched_dPtOverPt", "GsfTrack dp_{T}/p_{T}", 100, 0., 1.);
   hT_matched_etaError = fs->make<TH1D>("hT_matched_etaError", "GsfTrack #eta-error", 100, 0., 0.1);
   hT_matched_phiError = fs->make<TH1D>("hT_matched_phiError", "GsfTrack #phi-error", 100, 0., 0.1);
   hT_matched_outerP = fs->make<TH1D>("hT_matched_outerP", "GsfTrack p_{out}", 100, 0., 250.);
   hT_matched_outerPt = fs->make<TH1D>("hT_matched_outerPt", "GsfTrack p_{T}^{out}", 100, 0., 250.);
   hT_matched_outerEta = fs->make<TH1D>("hT_matched_outerEta", "GsfTrack #eta_{out}", 100, -3., 3.);
   hT_matched_outerPhi = fs->make<TH1D>("hT_matched_outerPhi", "GsfTrack #phi_{out}", 100, -3.2, 3.2);
   hT_matched_outerPOverP = fs->make<TH1D>("hT_matched_outerPOverP", "GsfTrack p_{out}/p", 100, 0., 2.);
   hT_matched_outerPtOverPt = fs->make<TH1D>("hT_matched_outerPtOverPt", "GsfTrack p_{T}^{out}/p_{T}", 100, 0., 2.);
   hT_matched_chi2 = fs->make<TH1D>("hT_matched_chi2", "GsfTrack #chi^{2}", 100, 0., 100.);
   hT_matched_ndof = fs->make<TH1D>("hT_matched_ndof", "GsfTrack ndof", 45, 0., 45.);
   hT_matched_normalizedChi2 = fs->make<TH1D>("hT_matched_normalizedChi2", "GsfTrack #chi^{2}_{normalized}", 100, 0., 50.);
   hT_matched_charge = fs->make<TH1D>("hT_matched_charge", "GsfTrack charge", 3, -1., 2.);
   hT_matched_dxy = fs->make<TH1D>("hT_matched_dxy", "GsfTrack d_{xy}", 100, -25., 25.);
   hT_matched_dz = fs->make<TH1D>("hT_matched_dz", "GsfTrack d_{z}", 100, -100., 100.);
   hT_matched_dxyError = fs->make<TH1D>("hT_matched_dxyError", "GsfTrack d_{xy}-error", 100, 0., 5.);
   hT_matched_dzError = fs->make<TH1D>("hT_matched_dzError", "GsfTrack d_{z}-error", 100, 0., 10.);
   hT_matched_numberOfValidHits = fs->make<TH1D>("hT_matched_numberOfValidHits", "GsfTrack #hits_{valid}", 25, 0., 25.);
   hT_matched_numberOfLostHits = fs->make<TH1D>("hT_matched_numberOfLostHits", "GsfTrack #hits_{lost}", 5, 0., 5.);
   hT_matched_validFraction = fs->make<TH1D>("hT_matched_validFraction", "GsfTrack valid fraction", 100, 0., 1.);

   hMatch_ptCut_size = fs->make<TH1D>("hMatch_ptCut_size", "GsfTrack_{matched}: size", 35, 0., 35.);
   hMatch_ptCut_dr = fs->make<TH1D>("hMatch_ptCut_dr", "GsfTrack_{matched}: #DeltaR", 100, 0., 0.5);
   hMatch_ptCut_nSharedHits = fs->make<TH1D>("hMatch_ptCut_nSharedHits", "GsfTrack_{matched}: # hits_{shared}", 30, 0., 30);
   hMatch_ptCut_nNonSharedHitsRef = fs->make<TH1D>("hMatch_ptCut_nNonSharedHitsRef", "GsfTrack_{matched}^{ref}: # hits_{not shared}", 10, 0., 10.);
   hMatch_ptCut_nNonSharedHitsTest = fs->make<TH1D>("hMatch_ptCut_nNonSharedHitsTest", "GsfTrack_{matched}^{test}: # hits_{not shared}", 10, 0., 10.);
   hMatch_ptCut_diff_p = fs->make<TH1D>("hMatch_ptCut_diff_p", "GsfTrack_{matched}: p_{test} - p_{ref}", 100, -250., 250.);
   hMatch_ptCut_diff_pt = fs->make<TH1D>("hMatch_ptCut_diff_pt", "GsfTrack_{matched}: p_{T}^{test} - p_{T}^{ref}", 100, -250., 250.);
   hMatch_ptCut_diff_eta = fs->make<TH1D>("hMatch_ptCut_diff_eta", "GsfTrack_{matched}: #eta_{test} -  #eta_{ref}", 100, -1., 1.);
   hMatch_ptCut_diff_phi = fs->make<TH1D>("hMatch_ptCut_diff_phi", "GsfTrack_{matched}: #phi_{test} - #phi_{ref}", 100, -1., 1.);
   hMatch_ptCut_diff_ptError = fs->make<TH1D>("hMatch_ptCut_diff_ptError", "GsfTrack_{matched}: p_{T}-error_{test} - p_{T}-error_{ref}", 100, -20., 20.);
   hMatch_ptCut_diff_dPtOverPt = fs->make<TH1D>("hMatch_ptCut_diff_dPtOverPt", "GsfTrack_{matched}: dp_{T}/p_{T}^{test} - dp_{T}/p_{T}^{ref}", 100, -1., 1.);
   hMatch_ptCut_diff_etaError = fs->make<TH1D>("hMatch_ptCut_diff_etaError", "GsfTrack_{matched}: #eta-error_{test} - #eta-error_{ref}", 100, -0.1, 0.1);
   hMatch_ptCut_diff_phiError = fs->make<TH1D>("hMatch_ptCut_diff_phiError", "GsfTrack_{matched}: #phi-error_{test} - #phi-error_{ref}", 100, -0.1, 0.1);
   hMatch_ptCut_diff_outerP = fs->make<TH1D>("hMatch_ptCut_diff_outerP", "GsfTrack_{matched}: p_{out}^{test} - p_{out}^{ref}", 100, -250., 250.);
   hMatch_ptCut_diff_outerPt = fs->make<TH1D>("hMatch_ptCut_diff_outerPt", "GsfTrack_{matched}: (p_{T}^{out})_{test} - (p_{T}^{out})_{ref}", 100, -250., 250.);
   hMatch_ptCut_diff_outerEta = fs->make<TH1D>("hMatch_ptCut_diff_outerEta", "GsfTrack_{matched}: #eta_{out}^{test} - #eta_{out}^{ref}", 100, -1., 1.);
   hMatch_ptCut_diff_outerPhi = fs->make<TH1D>("hMatch_ptCut_diff_outerPhi", "GsfTrack_{matched}: #phi_{out}^{test} - #phi_{out}^{ref}", 100, -1., 1.);
   hMatch_ptCut_diff_outerPOverP = fs->make<TH1D>("hMatch_ptCut_diff_outerPOverP", "GsfTrack_{matched}: (p_{out}/p)_{test} - (p_{out}/p)_{ref}", 100, -4., 4.);
   hMatch_ptCut_diff_outerPtOverPt = fs->make<TH1D>("hMatch_ptCut_diff_outerPtOverPt", "GsfTrack_{matched}: (p_{T}^{out}/p_{T})_{test} - (p_{T}^{out}/p_{T})_{ref}", 100, -4., 4.);
   hMatch_ptCut_diff_chi2 = fs->make<TH1D>("hMatch_ptCut_diff_chi2", "GsfTrack_{matched}: #chi^{2}_{test} - #chi^{2}_{ref}", 100, -150., 150.);
   hMatch_ptCut_diff_ndof = fs->make<TH1D>("hMatch_ptCut_diff_ndof", "GsfTrack_{matched}: ndof_{test} - ndof_{ref}", 40, -20., 20.);
   hMatch_ptCut_diff_normalizedChi2 = fs->make<TH1D>("hMatch_ptCut_diff_normalizedChi2", "GsfTrack_{matched}: (#chi^{2}_{normalized})_{test} - (#chi^{2}_{normalized})_{ref}", 100, -50., 50.);
   hMatch_ptCut_diff_charge = fs->make<TH1D>("hMatch_ptCut_diff_charge", "GsfTrack_{matched}: charge_{test} - charge_{ref}", 5, -2., 3.);
   hMatch_ptCut_diff_dxy = fs->make<TH1D>("hMatch_ptCut_diff_dxy", "GsfTrack_{matched}: d_{xy}^{test} - d_{xy}^{ref}", 100, -10., 10.);
   hMatch_ptCut_diff_dz = fs->make<TH1D>("hMatch_ptCut_diff_dz", "GsfTrack_{matched}: d_{z}^{test} - d_{z}^{ref}", 100, -50., 50.);
   hMatch_ptCut_diff_dxyError = fs->make<TH1D>("hMatch_ptCut_diff_dxyError", "GsfTrack_{matched}: d_{xy}-error_{test} - d_{xy}-error_{ref}", 100, -5., 5.);
   hMatch_ptCut_diff_dzError = fs->make<TH1D>("hMatch_ptCut_diff_dzError", "GsfTrack_{matched}: d_{z}-error_{test} - d_{z}-error_{ref}", 100, -10., 10.);
   hMatch_ptCut_diff_numberOfValidHits = fs->make<TH1D>("hMatch_ptCut_diff_numberOfValidHits", "GsfTrack_{matched}: #hits_{valid}^{test} - #hits_{valid}^{ref}", 20, -10., 10.);
   hMatch_ptCut_diff_numberOfLostHits = fs->make<TH1D>("hMatch_ptCut_diff_numberOfLostHits", "GsfTrack_{matched}: #hits_{lost}^{test} - #hits_{lost}^{ref}", 10, -5., 5.);
   hMatch_ptCut_diff_validFraction = fs->make<TH1D>("hMatch_ptCut_diff_validFraction", "GsfTrack_{matched}: valid fraction_{test} - valid fraction_{ref}", 100, -1., 1.);

   hMatch_ptCut_pull_p = fs->make<TH1D>("hMatch_ptCut_pull_p", "GsfTrack_{matched}: |p_{test} - p_{ref}|/p_{ref}", 100, 0., 30.);
   hMatch_ptCut_pull_pt = fs->make<TH1D>("hMatch_ptCut_pull_pt", "GsfTrack_{matched}: |p_{T}^{test} - p_{T}^{ref}|/p_{T}^{ref}", 100, 0., 20.);
   hMatch_ptCut_pull_eta = fs->make<TH1D>("hMatch_ptCut_pull_eta", "GsfTrack_{matched}: |#eta_{test} -  #eta_{ref}|/#eta_{ref}", 100, 0., 1.);
   hMatch_ptCut_pull_ptError = fs->make<TH1D>("hMatch_ptCut_pull_ptError", "GsfTrack_{matched}: |p_{T}-error_{test} - p_{T}-error_{ref}|/p_{T}-error_{ref}", 100, 0., 20.);
   hMatch_ptCut_pull_dPtOverPt = fs->make<TH1D>("hMatch_ptCut_pull_dPtOverPt", "GsfTrack_{matched}: |dp_{T}/p_{T}^{test} - dp_{T}/p_{T}^{ref}|/(dp_{T}/p_{T}^{ref})", 100, 0., 3.);
   hMatch_ptCut_pull_etaError = fs->make<TH1D>("hMatch_ptCut_pull_etaError", "GsfTrack_{matched}: |#eta-error_{test} - #eta-error_{ref}|/#eta-error_{ref}", 100, 0., 3.);
   hMatch_ptCut_pull_phiError = fs->make<TH1D>("hMatch_ptCut_pull_phiError", "GsfTrack_{matched}: |#phi-error_{test} - #phi-error_{ref}|/#phi-error_{ref}", 100, 0., 3.);
   hMatch_ptCut_pull_outerP = fs->make<TH1D>("hMatch_ptCut_pull_outerP", "GsfTrack_{matched}: |p_{out}^{test} - p_{out}^{ref}|/p_{out}^{ref}", 100, 0., 25.);
   hMatch_ptCut_pull_outerPt = fs->make<TH1D>("hMatch_ptCut_pull_outerPt", "GsfTrack_{matched}: |(p_{T}^{out})_{test} - (p_{T}^{out})_{ref}|/(p_{T}^{out})_{ref}", 100, 0., 25.);
   hMatch_ptCut_pull_outerEta = fs->make<TH1D>("hMatch_ptCut_pull_outerEta", "GsfTrack_{matched}: |#eta_{out}^{test} - #eta_{out}^{ref}|/#eta_{out}^{ref}", 100, 0., 1.);
   hMatch_ptCut_pull_outerPOverP = fs->make<TH1D>("hMatch_ptCut_pull_outerPOverP", "GsfTrack_{matched}: |(p_{out}/p)_{test} - (p_{out}/p)_{ref}|/(p_{out}/p)_{ref}", 100, 0., 8.);
   hMatch_ptCut_pull_outerPtOverPt = fs->make<TH1D>("hMatch_ptCut_pull_outerPtOverPt", "GsfTrack_{matched}: |(p_{T}^{out}/p_{T})_{test} - (p_{T}^{out}/p_{T})_{ref}|/(p_{T}^{out}/p_{T})_{ref}", 100, 0., 8.);
   hMatch_ptCut_pull_chi2 = fs->make<TH1D>("hMatch_ptCut_pull_chi2", "GsfTrack_{matched}: |#chi^{2}_{test} - #chi^{2}_{ref}|/#chi^{2}_{ref}", 100, 0., 20.);
   hMatch_ptCut_pull_ndof = fs->make<TH1D>("hMatch_ptCut_pull_ndof", "GsfTrack_{matched}: |ndof_{test} - ndof_{ref}|/ndof_{ref}", 100, 0., 3.);
   hMatch_ptCut_pull_normalizedChi2 = fs->make<TH1D>("hMatch_ptCut_pull_normalizedChi2", "GsfTrack_{matched}: |(#chi^{2}_{normalized})_{test} - (#chi^{2}_{normalized})_{ref}|/(#chi^{2}_{normalized})_{ref}", 100, 0., 10.);
   hMatch_ptCut_pull_dxy = fs->make<TH1D>("hMatch_ptCut_pull_dxy", "GsfTrack_{matched}: |d_{xy}^{test} - d_{xy}^{ref}|/|d_{xy}^{ref}|", 100, 0., 10.);
   hMatch_ptCut_pull_dz = fs->make<TH1D>("hMatch_ptCut_pull_dz", "GsfTrack_{matched}: |d_{z}^{test} - d_{z}^{ref}|/|d_{z}^{ref}|", 100, 0., 10.);
   hMatch_ptCut_pull_dxyError = fs->make<TH1D>("hMatch_ptCut_pull_dxyError", "GsfTrack_{matched}: |d_{xy}-error_{test} - d_{xy}-error_{ref}|/d_{xy}-error_{ref}", 100, 0., 4.);
   hMatch_ptCut_pull_dzError = fs->make<TH1D>("hMatch_ptCut_pull_dzError", "GsfTrack_{matched}: |d_{z}-error_{test} - d_{z}-error_{ref}|/d_{z}-error_{ref}", 100, 0., 4.);
   hMatch_ptCut_pull_numberOfValidHits = fs->make<TH1D>("hMatch_ptCut_pull_numberOfValidHits", "GsfTrack_{matched}: (#hits_{valid}^{test} - #hits_{valid}^{ref}|/#hits_{valid}^{ref}", 100, 0., 5.);
   hMatch_ptCut_pull_numberOfLostHits = fs->make<TH1D>("hMatch_ptCut_pull_numberOfLostHits", "GsfTrack_{matched}: |#hits_{lost}^{test} - #hits_{lost}^{ref}|/#hits_{lost}^{ref}", 100, 0., 5.);
   hMatch_ptCut_pull_validFraction = fs->make<TH1D>("hMatch_ptCut_pull_validFraction", "GsfTrack_{matched}: |valid fraction_{test} - valid fraction_{ref}|/valid fraction_{ref}", 100, 0., 1.);

   hR_dEta_mapSize = fs->make<TH1D>("hR_dEta_mapSize", "GsfTrack-SC: mapSize_{#Delta#eta}", 20, 0., 20.);
   hR_dEta_mapSizeMatch = fs->make<TH1D>("hR_dEta_mapSizeMatch", "GsfTrack-SC: mapSize^{matched}_{#Delta#eta}", 20, 0., 20.);
   hR_dEta_dEta = fs->make<TH1D>("hR_dEta_dEta", "GsfTrack-SC: #Delta#eta", 100, 0., 1.);
   hR_dEta_nPassCut = fs->make<TH1D>("hR_dEta_nPassCut", "GsfTrack-SC: #passCut_{#Delta#eta}", 20, 0., 20.);
   hR_dEta_recoEcalCandsSize = fs->make<TH1D>("hR_dEta_recoEcalCandsSize", "GsfTrack-SC: recoEcalCandsSize_{#Delta#eta}", 5, 0., 5.);
   hR_dEta_recoEcalCandsSizeMatch = fs->make<TH1D>("hR_dEta_recoEcalCandsSizeMatch", "GsfTrack-SC: recoEcalCandsSize^{matched}_{#Delta#eta}", 5, 0., 5.);
   hR_dEta_dEta_trgdEcalCands = fs->make<TH1D>("hR_dEta_dEta_trgdEcalCands", "GsfTrack-SC: #Delta#eta_{passing recoEcalCands}", 100, 0., 0.1);
   hR_dEta_nPassCut_trgdEcalCands = fs->make<TH1D>("hR_dEta_nPassCut_trgdEcalCands", "GsfTrack-SC: #passCut^{passing recoEcalCands}_{#Delta#eta}", 5, 0., 5.);

   hR_dPhi_mapSize = fs->make<TH1D>("hR_dPhi_mapSize", "GsfTrack-SC: mapSize_{#Delta#phi}", 20, 0., 20.);
   hR_dPhi_mapSizeMatch = fs->make<TH1D>("hR_dPhi_mapSizeMatch", "GsfTrack-SC: mapSize^{matched}_{#Delta#phi}", 20, 0., 20.);
   hR_dPhi_dPhi = fs->make<TH1D>("hR_dPhi_dPhi", "GsfTrack-SC: #Delta#phi", 100, 0., 1.);
   hR_dPhi_nPassCut = fs->make<TH1D>("hR_dPhi_nPassCut", "GsfTrack-SC: #passCut_{#Delta#phi}", 20, 0., 20.);
   hR_dPhi_recoEcalCandsSize = fs->make<TH1D>("hR_dPhi_recoEcalCandsSize", "GsfTrack-SC: recoEcalCandsSize_{#Delta#phi}", 5, 0., 5.);
   hR_dPhi_recoEcalCandsSizeMatch = fs->make<TH1D>("hR_dPhi_recoEcalCandsSizeMatch", "GsfTrack-SC: recoEcalCandsSize^{matched}_{#Delta#phi}", 5, 0., 5.);
   hR_dPhi_dPhi_trgdEcalCands = fs->make<TH1D>("hR_dPhi_dPhi_trgdEcalCands", "GsfTrack-SC: #Delta#phi_{passing recoEcalCands}", 100, 0., 0.6);
   hR_dPhi_nPassCut_trgdEcalCands = fs->make<TH1D>("hR_dPhi_nPassCut_trgdEcalCands", "GsfTrack-SC: #passCut^{passing recoEcalCands}_{#Delta#phi}", 5, 0., 5.);

   hT_dEta_mapSize = fs->make<TH1D>("hT_dEta_mapSize", "GsfTrack-SC: mapSize_{#Delta#eta}", 20, 0., 20.);
   hT_dEta_mapSizeMatch = fs->make<TH1D>("hT_dEta_mapSizeMatch", "GsfTrack-SC: mapSize^{matched}_{#Delta#eta}", 20, 0., 20.);
   hT_dEta_dEta = fs->make<TH1D>("hT_dEta_dEta", "GsfTrack-SC: #Delta#eta", 100, 0., 1.);
   hT_dEta_nPassCut = fs->make<TH1D>("hT_dEta_nPassCut", "GsfTrack-SC: #passCut_{#Delta#eta}", 20, 0., 20.);
   hT_dEta_recoEcalCandsSize = fs->make<TH1D>("hT_dEta_recoEcalCandsSize", "GsfTrack-SC: recoEcalCandsSize_{#Delta#eta}", 5, 0., 5.);
   hT_dEta_recoEcalCandsSizeMatch = fs->make<TH1D>("hT_dEta_recoEcalCandsSizeMatch", "GsfTrack-SC: recoEcalCandsSize^{matched}_{#Delta#eta}", 5, 0., 5.);
   hT_dEta_dEta_trgdEcalCands = fs->make<TH1D>("hT_dEta_dEta_trgdEcalCands", "GsfTrack-SC: #Delta#eta_{passing recoEcalCands}", 100, 0., 0.1);
   hT_dEta_nPassCut_trgdEcalCands = fs->make<TH1D>("hT_dEta_nPassCut_trgdEcalCands", "GsfTrack-SC: #passCut^{passing recoEcalCands}_{#Delta#eta}", 5, 0., 5.);

   hT_dPhi_mapSize = fs->make<TH1D>("hT_dPhi_mapSize", "GsfTrack-SC: mapSize_{#Delta#phi}", 20, 0., 20.);
   hT_dPhi_mapSizeMatch = fs->make<TH1D>("hT_dPhi_mapSizeMatch", "GsfTrack-SC: mapSize^{matched}_{#Delta#phi}", 20, 0., 20.);
   hT_dPhi_dPhi = fs->make<TH1D>("hT_dPhi_dPhi", "GsfTrack-SC: #Delta#phi", 100, 0., 1.);
   hT_dPhi_nPassCut = fs->make<TH1D>("hT_dPhi_nPassCut", "GsfTrack-SC: #passCut_{#Delta#phi}", 20, 0., 20.);
   hT_dPhi_recoEcalCandsSize = fs->make<TH1D>("hT_dPhi_recoEcalCandsSize", "GsfTrack-SC: recoEcalCandsSize_{#Delta#phi}", 5, 0., 5.);
   hT_dPhi_recoEcalCandsSizeMatch = fs->make<TH1D>("hT_dPhi_recoEcalCandsSizeMatch", "GsfTrack-SC: recoEcalCandsSize^{matched}_{#Delta#phi}", 5, 0., 5.);
   hT_dPhi_dPhi_trgdEcalCands = fs->make<TH1D>("hT_dPhi_dPhi_trgdEcalCands", "GsfTrack-SC: #Delta#phi_{passing recoEcalCands}", 100, 0., 0.6);
   hT_dPhi_nPassCut_trgdEcalCands = fs->make<TH1D>("hT_dPhi_nPassCut_trgdEcalCands", "GsfTrack-SC: #passCut^{passing recoEcalCands}_{#Delta#phi}", 5, 0., 5.);

   //tree = fs->make<TTree>("tree","tree");
   //tree->Branch("trkSize", &trkSize, "trkSize/i");
   //tree->Branch("pt", pt, "pt/F");
}


HltGsfTrkAnalyser::~HltGsfTrkAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HltGsfTrkAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // trigger requirements
   edm::TriggerResultsByName trgRes = iEvent.triggerResultsByName(testTrackCollTag_.process());
   //for (unsigned int i = 0; i < trgRes.size(); ++i) std::cout << trgRes.index(i) << " " << trgRes.triggerName(i) << " " << trgRes.accept(i) << std::endl;
   bool refTrgAcc = trgRes.accept(refTrgName_);
   bool testTrgAcc = trgRes.accept(testTrgName_);
   if (!(refTrgAcc ^ refTrgInv_)) return;
   if (!(testTrgAcc ^ testTrgInv_)) return;

   // get the gsftrack collections
   edm::InputTag inTagR(referenceTrackCollTag_);
   edm::InputTag inTagT(testTrackCollTag_);
   edm::Handle<reco::GsfTrackCollection > trkHandleR;
   if (!iEvent.getByLabel(inTagR, trkHandleR)) return;
   const reco::GsfTrackCollection* trksR = trkHandleR.product();
   edm::Handle<reco::GsfTrackCollection > trkHandleT;
   if (!iEvent.getByLabel(inTagT, trkHandleT)) return;
   const reco::GsfTrackCollection* trksT = trkHandleT.product();

   //trkRSize = 0;
   unsigned int trksPassingPtCut = 0;
   hR_size->Fill(trksR->size());
   for (reco::GsfTrackCollection::const_iterator it = trksR->begin(); it < trksR->end(); ++it) {
      hR_p->Fill(it->p());
      hR_pt->Fill(it->pt());
      hR_eta->Fill(it->eta());
      hR_phi->Fill(it->phi());
      hR_ptError->Fill(it->ptError());
      hR_dPtOverPt->Fill(it->ptError() / it->pt());
      hR_etaError->Fill(it->etaError());
      hR_phiError->Fill(it->phiError());
      hR_outerP->Fill(it->outerP());
      hR_outerPt->Fill(it->outerPt());
      hR_outerEta->Fill(it->outerEta());
      hR_outerPhi->Fill(it->outerPhi());
      hR_outerPOverP->Fill(it->outerP() / it->p());
      hR_outerPtOverPt->Fill(it->outerPt() / it->pt());
      hR_chi2->Fill(it->chi2());
      hR_ndof->Fill(it->ndof());
      hR_normalizedChi2->Fill(it->normalizedChi2());
      hR_charge->Fill(it->charge());
      hR_dxy->Fill(it->dxy());
      hR_dz->Fill(it->dz());
      hR_dxyError->Fill(it->dxyError());
      hR_dzError->Fill(it->dzError());
      hR_numberOfValidHits->Fill(it->numberOfValidHits());
      hR_numberOfLostHits->Fill(it->numberOfLostHits());
      hR_validFraction->Fill(it->validFraction());
      hR2d_pt_eta->Fill(it->pt(), it->eta());
      hR2d_pt_dpt->Fill(it->pt(), it->ptError());
      hR2d_p_outerP->Fill(it->p(), it->outerP());
      hR2d_pt_outerPt->Fill(it->pt(), it->outerPt());

      if (it->pt() > minPtR_) {
         ++trksPassingPtCut;
         hR_ptCut_p->Fill(it->p());
         hR_ptCut_pt->Fill(it->pt());
         hR_ptCut_eta->Fill(it->eta());
         hR_ptCut_phi->Fill(it->phi());
         hR_ptCut_ptError->Fill(it->ptError());
         hR_ptCut_dPtOverPt->Fill(it->ptError() / it->pt());
         hR_ptCut_etaError->Fill(it->etaError());
         hR_ptCut_phiError->Fill(it->phiError());
         hR_ptCut_outerP->Fill(it->outerP());
         hR_ptCut_outerPt->Fill(it->outerPt());
         hR_ptCut_outerEta->Fill(it->outerEta());
         hR_ptCut_outerPhi->Fill(it->outerPhi());
         hR_ptCut_outerPOverP->Fill(it->outerP() / it->p());
         hR_ptCut_outerPtOverPt->Fill(it->outerPt() / it->pt());
         hR_ptCut_chi2->Fill(it->chi2());
         hR_ptCut_ndof->Fill(it->ndof());
         hR_ptCut_normalizedChi2->Fill(it->normalizedChi2());
         hR_ptCut_charge->Fill(it->charge());
         hR_ptCut_dxy->Fill(it->dxy());
         hR_ptCut_dz->Fill(it->dz());
         hR_ptCut_dxyError->Fill(it->dxyError());
         hR_ptCut_dzError->Fill(it->dzError());
         hR_ptCut_numberOfValidHits->Fill(it->numberOfValidHits());
         hR_ptCut_numberOfLostHits->Fill(it->numberOfLostHits());
         hR_ptCut_validFraction->Fill(it->validFraction());
      }

      //pt[trkRSize++] = it->pt();
   }
   hR_ptCut_size->Fill(trksPassingPtCut);
   trksPassingPtCut = 0;

   hT_size->Fill(trksT->size());
   for (reco::GsfTrackCollection::const_iterator it = trksT->begin(); it < trksT->end(); ++it) {
      hT_p->Fill(it->p());
      hT_pt->Fill(it->pt());
      hT_eta->Fill(it->eta());
      hT_phi->Fill(it->phi());
      hT_ptError->Fill(it->ptError());
      hT_dPtOverPt->Fill(it->ptError() / it->pt());
      hT_etaError->Fill(it->etaError());
      hT_phiError->Fill(it->phiError());
      hT_outerP->Fill(it->outerP());
      hT_outerPt->Fill(it->outerPt());
      hT_outerEta->Fill(it->outerEta());
      hT_outerPhi->Fill(it->outerPhi());
      hT_outerPOverP->Fill(it->outerP() / it->p());
      hT_outerPtOverPt->Fill(it->outerPt() / it->pt());
      hT_chi2->Fill(it->chi2());
      hT_ndof->Fill(it->ndof());
      hT_normalizedChi2->Fill(it->normalizedChi2());
      hT_charge->Fill(it->charge());
      hT_dxy->Fill(it->dxy());
      hT_dz->Fill(it->dz());
      hT_dxyError->Fill(it->dxyError());
      hT_dzError->Fill(it->dzError());
      hT_numberOfValidHits->Fill(it->numberOfValidHits());
      hT_numberOfLostHits->Fill(it->numberOfLostHits());
      hT_validFraction->Fill(it->validFraction());

      if (it->pt() > minPtT_) {
         ++trksPassingPtCut;
         hT_ptCut_p->Fill(it->p());
         hT_ptCut_pt->Fill(it->pt());
         hT_ptCut_eta->Fill(it->eta());
         hT_ptCut_phi->Fill(it->phi());
         hT_ptCut_ptError->Fill(it->ptError());
         hT_ptCut_dPtOverPt->Fill(it->ptError() / it->pt());
         hT_ptCut_etaError->Fill(it->etaError());
         hT_ptCut_phiError->Fill(it->phiError());
         hT_ptCut_outerP->Fill(it->outerP());
         hT_ptCut_outerPt->Fill(it->outerPt());
         hT_ptCut_outerEta->Fill(it->outerEta());
         hT_ptCut_outerPhi->Fill(it->outerPhi());
         hT_ptCut_outerPOverP->Fill(it->outerP() / it->p());
         hT_ptCut_outerPtOverPt->Fill(it->outerPt() / it->pt());
         hT_ptCut_chi2->Fill(it->chi2());
         hT_ptCut_ndof->Fill(it->ndof());
         hT_ptCut_normalizedChi2->Fill(it->normalizedChi2());
         hT_ptCut_charge->Fill(it->charge());
         hT_ptCut_dxy->Fill(it->dxy());
         hT_ptCut_dz->Fill(it->dz());
         hT_ptCut_dxyError->Fill(it->dxyError());
         hT_ptCut_dzError->Fill(it->dzError());
         hT_ptCut_numberOfValidHits->Fill(it->numberOfValidHits());
         hT_ptCut_numberOfLostHits->Fill(it->numberOfLostHits());
         hT_ptCut_validFraction->Fill(it->validFraction());
      }
   }
   hT_ptCut_size->Fill(trksPassingPtCut);

   // track matching
   typedef edm::AssociationMap<edm::OneToOne<reco::GsfTrackCollection, reco::GsfTrackCollection> > TrackToTrackMap;
   TrackToTrackMap trkTrkMatchMap;
   // make RefVectors for both track collections
   reco::GsfTrackRefVector trkRefsR;
   reco::GsfTrackRefVector trkRefsT;
   for (std::size_t i = 0; i < trksR->size(); ++i) {
      reco::GsfTrackRef trkRefR(trksR, i);
      trkRefsR.push_back(trkRefR);
   }
   for (std::size_t i = 0; i < trksT->size(); ++i) {
      reco::GsfTrackRef trkRefT(trksT, i);
      trkRefsT.push_back(trkRefT);
   }

   // loop over reference tracks
   reco::GsfTrackRefVector::const_iterator it = trkRefsR.begin();
   unsigned int matchTrkCtr = 0;
   unsigned int ambTrkCtr = 0;
   unsigned int nJumps = 0;
   while (it < trkRefsR.end()) {
      //std::cout << trkRefsR.size() << " " << trkRefsT.size() << std::endl;
      if (trkRefsR.size() * trkRefsT.size() == 0) break;
      if (nJumps > trkRefsR.size() + trkRefsT.size()) {
         std::cout << "Catched infinit loop. Break." << std::endl;
         nJumps = 0;
         break; // have to implement something other than break. Could take one track out of the loop or all?
      }
      double bestDiscrFwd = 100. * discrCut_;
      double bestDiscrBwd = 100. * discrCut_;
      if (largerThan_) {
         bestDiscrFwd = 0.;
         bestDiscrBwd = 0.;
      }
      reco::GsfTrackRefVector::iterator bestItFwd = trkRefsT.begin();
      reco::GsfTrackRefVector::iterator bestItBwd = trkRefsR.begin();
      // search forward for best match
      for (reco::GsfTrackRefVector::const_iterator itT = trkRefsT.begin(); itT < trkRefsT.end(); ++itT) { 
         double trkTrkDiscr = TrackMatcher(*(it->get()), *(itT->get()));
         if (!largerThan_ && (trkTrkDiscr < bestDiscrFwd)) {
            bestDiscrFwd = trkTrkDiscr;
            bestItFwd = itT;
         }
         else if (largerThan_ && (trkTrkDiscr > bestDiscrFwd)) {
            bestDiscrFwd = trkTrkDiscr;
            bestItFwd = itT;
         }
         //std::cout << "Fwd: " << trkTrkDiscr << "  " << bestDiscrFwd << "  " << it->get()->pt() << "  " << itT->get()->pt() << "  " << bestItFwd->get()->pt() << std::endl;
      }
      // search backward for best match
      for (reco::GsfTrackRefVector::const_iterator itR = trkRefsR.begin(); itR < trkRefsR.end(); ++itR) { 
         double trkTrkDiscr = TrackMatcher(*(bestItFwd->get()), *(itR->get()));
         if (!largerThan_ && (trkTrkDiscr < bestDiscrBwd)) {
            bestDiscrBwd = trkTrkDiscr;
            bestItBwd = itR;
         }
         else if (largerThan_ && (trkTrkDiscr > bestDiscrBwd)) {
            bestDiscrBwd = trkTrkDiscr;
            bestItBwd = itR;
         }
         //std::cout << "Bwd: " << trkTrkDiscr << "  " << bestDiscrBwd << "  " << itR->get()->pt() << "  " << bestItFwd->get()->pt() << "  " << bestItBwd->get()->pt() << std::endl;
      }
      // do we have a match?
      if (bestDiscrFwd == bestDiscrBwd) {
         //std::cout << "Match: " << bestDiscrFwd << "  " << bestItFwd->get()->pt() << "  " << bestItBwd->get()->pt() << std::endl;
         // fill AssociationMap and then erase the matched TrackRefs from RefVectors
         if (!largerThan_ && (bestDiscrFwd < discrCut_)) trkTrkMatchMap.insert(reco::GsfTrackRef(bestItBwd->product(), bestItBwd->index()), reco::GsfTrackRef(bestItFwd->product(), bestItFwd->index()));
         else if (largerThan_ && (bestDiscrFwd > discrCut_)) trkTrkMatchMap.insert(reco::GsfTrackRef(bestItBwd->product(), bestItBwd->index()), reco::GsfTrackRef(bestItFwd->product(), bestItFwd->index()));
         trkRefsR.erase(bestItBwd);
         trkRefsT.erase(bestItFwd);
         it = trkRefsR.begin();
         nJumps = 0;
         ++matchTrkCtr;
      } else {
         //std::cout << "No match: " << bestDiscrFwd << "  " << bestDiscrBwd << std::endl;
         it = bestItBwd;
         nJumps += 2;
         ++ambTrkCtr;
      }
   }
   //std::cout << "Tracks matched: " << matchTrkCtr << ", Ambiguous tracks: " << ambTrkCtr << ", TrackToTrackMap size: " << trkTrkMatchMap.size() << std::endl;

   // test matches
   reco::GsfTrackCollection::const_iterator tIt = trksR->begin();
   for (std::size_t i = 0; i < trksR->size(); ++i, ++tIt) {
      reco::GsfTrackRef trkRefKey(trksR, i);
      TrackToTrackMap::const_iterator mIt = trkTrkMatchMap.find(trkRefKey);
      if (mIt == trkTrkMatchMap.end()) {
         //std::cout << "No Match found for this track" << std::endl;
         continue;
      }
      reco::GsfTrackRef trkRefVal = mIt->val;
      //std::cout << trkRefKey.get()->pt() << " " << trkRefVal->pt() << std::endl;

      double dR = reco::deltaR(*trkRefVal, *trkRefKey.get());
      double nSharedHits = GetNSharedHits(*trkRefVal, *trkRefKey.get());

      //if (trkRefVal->ptError() / trkRefVal->pt() > 1.) continue;
      //if (trkRefKey.get()->ptError() / trkRefKey.get()->pt() > 1.) continue;

      hMatch_size->Fill(trkTrkMatchMap.size());
      hMatch_dr->Fill(dR);
      hMatch_nSharedHits->Fill(nSharedHits);
      hMatch_nNonSharedHitsRef->Fill(trkRefKey.get()->numberOfValidHits() - nSharedHits);
      hMatch_nNonSharedHitsTest->Fill(trkRefVal->numberOfValidHits() - nSharedHits);
      hMatch_diff_p->Fill(trkRefVal->p() - trkRefKey.get()->p());
      hMatch_diff_pt->Fill(trkRefVal->pt() - trkRefKey.get()->pt());
      hMatch_diff_eta->Fill(trkRefVal->eta() - trkRefKey.get()->eta());
      double phiDiff = trkRefVal->phi() - trkRefKey.get()->phi();
      if (phiDiff > ROOT::Math::Pi()) phiDiff -= 2. * ROOT::Math::Pi();
      else if (phiDiff < -1*ROOT::Math::Pi()) phiDiff += 2. * ROOT::Math::Pi();
      hMatch_diff_phi->Fill(phiDiff);
      hMatch_diff_ptError->Fill(trkRefVal->ptError() - trkRefKey.get()->ptError());
      hMatch_diff_dPtOverPt->Fill(trkRefVal->ptError() / trkRefVal->pt() - trkRefKey.get()->ptError() / trkRefKey.get()->pt());
      hMatch_diff_etaError->Fill(trkRefVal->etaError() - trkRefKey.get()->etaError());
      hMatch_diff_phiError->Fill(trkRefVal->phiError() - trkRefKey.get()->phiError());
      hMatch_diff_outerP->Fill(trkRefVal->outerP() - trkRefKey.get()->outerP());
      hMatch_diff_outerPt->Fill(trkRefVal->outerPt() - trkRefKey.get()->outerPt());
      hMatch_diff_outerEta->Fill(trkRefVal->outerEta() - trkRefKey.get()->outerEta());
      double outerPhiDiff = trkRefVal->phi() - trkRefKey.get()->phi();
      if (outerPhiDiff > ROOT::Math::Pi()) outerPhiDiff -= 2. * ROOT::Math::Pi();
      else if (outerPhiDiff < -1*ROOT::Math::Pi()) outerPhiDiff += 2. * ROOT::Math::Pi();
      hMatch_diff_outerPhi->Fill(outerPhiDiff);
      hMatch_diff_outerPOverP->Fill(trkRefVal->outerP() / trkRefVal->p() - trkRefKey.get()->outerP() / trkRefKey.get()->p());
      hMatch_diff_outerPtOverPt->Fill(trkRefVal->outerPt() / trkRefVal->pt() - trkRefKey.get()->outerPt() / trkRefKey.get()->pt());
      hMatch_diff_chi2->Fill(trkRefVal->chi2() - trkRefKey.get()->chi2());
      hMatch_diff_ndof->Fill(trkRefVal->ndof() - trkRefKey.get()->ndof());
      hMatch_diff_normalizedChi2->Fill(trkRefVal->normalizedChi2() - trkRefKey.get()->normalizedChi2());
      hMatch_diff_charge->Fill(trkRefVal->charge() - trkRefKey.get()->charge());
      hMatch_diff_dxy->Fill(trkRefVal->dxy() - trkRefKey.get()->dxy());
      hMatch_diff_dz->Fill(trkRefVal->dz() - trkRefKey.get()->dz());
      hMatch_diff_dxyError->Fill(trkRefVal->dxyError() - trkRefKey.get()->dxyError());
      hMatch_diff_dzError->Fill(trkRefVal->dzError() - trkRefKey.get()->dzError());
      hMatch_diff_numberOfValidHits->Fill(trkRefVal->numberOfValidHits() - trkRefKey.get()->numberOfValidHits());
      hMatch_diff_numberOfLostHits->Fill(trkRefVal->numberOfLostHits() - trkRefKey.get()->numberOfLostHits());
      hMatch_diff_validFraction->Fill(trkRefVal->validFraction() - trkRefKey.get()->validFraction());

      hMatch_pull_p->Fill(fabs(trkRefVal->p() - trkRefKey.get()->p()) / trkRefKey.get()->p());
      hMatch_pull_pt->Fill(fabs(trkRefVal->pt() - trkRefKey.get()->pt()) / trkRefKey.get()->pt());
      hMatch_pull_eta->Fill(fabs(trkRefVal->eta() - trkRefKey.get()->eta()) / trkRefKey.get()->eta());
      hMatch_pull_ptError->Fill(fabs(trkRefVal->ptError() - trkRefKey.get()->ptError()) / trkRefKey.get()->ptError());
      hMatch_pull_dPtOverPt->Fill(fabs(trkRefVal->ptError() / trkRefVal->pt() - trkRefKey.get()->ptError() / trkRefKey.get()->pt()) / (trkRefKey.get()->ptError() / trkRefKey.get()->pt()));
      hMatch_pull_etaError->Fill(fabs(trkRefVal->etaError() - trkRefKey.get()->etaError()) / trkRefKey.get()->etaError());
      hMatch_pull_phiError->Fill(fabs(trkRefVal->phiError() - trkRefKey.get()->phiError()) / trkRefKey.get()->phiError());
      hMatch_pull_outerP->Fill(fabs(trkRefVal->outerP() - trkRefKey.get()->outerP()) / trkRefKey.get()->outerP());
      hMatch_pull_outerPt->Fill(fabs(trkRefVal->outerPt() - trkRefKey.get()->outerPt()) / trkRefKey.get()->outerPt());
      hMatch_pull_outerEta->Fill(fabs(trkRefVal->outerEta() - trkRefKey.get()->outerEta()) / trkRefKey.get()->outerEta());
      hMatch_pull_outerPOverP->Fill(fabs(trkRefVal->outerP() / trkRefVal->p() - trkRefKey.get()->outerP() / trkRefKey.get()->p()) / (trkRefKey.get()->outerP() / trkRefKey.get()->p()));
      hMatch_pull_outerPtOverPt->Fill(fabs(trkRefVal->outerPt() / trkRefVal->pt() - trkRefKey.get()->outerPt() / trkRefKey.get()->pt()) / (trkRefKey.get()->outerPt() / trkRefKey.get()->pt()));
      hMatch_pull_chi2->Fill(fabs(trkRefVal->chi2() - trkRefKey.get()->chi2()) / trkRefKey.get()->chi2());
      hMatch_pull_ndof->Fill(fabs(trkRefVal->ndof() - trkRefKey.get()->ndof()) / (float)trkRefKey.get()->ndof());
      hMatch_pull_normalizedChi2->Fill(fabs(trkRefVal->normalizedChi2() - trkRefKey.get()->normalizedChi2()) / trkRefKey.get()->normalizedChi2());
      hMatch_pull_dxy->Fill(fabs(trkRefVal->dxy() - trkRefKey.get()->dxy()) / fabs(trkRefKey.get()->dxy()));
      hMatch_pull_dz->Fill(fabs(trkRefVal->dz() - trkRefKey.get()->dz()) / fabs(trkRefKey.get()->dz()));
      hMatch_pull_dxyError->Fill(fabs(trkRefVal->dxyError() - trkRefKey.get()->dxyError()) / trkRefKey.get()->dxyError());
      hMatch_pull_dzError->Fill(fabs(trkRefVal->dzError() - trkRefKey.get()->dzError()) / trkRefKey.get()->dzError());
      if (trkRefKey.get()->numberOfValidHits() != 0) hMatch_pull_numberOfValidHits->Fill(fabs(trkRefVal->numberOfValidHits() - trkRefKey.get()->numberOfValidHits()) / (float)trkRefKey.get()->numberOfValidHits());
      if (trkRefKey.get()->numberOfLostHits() != 0) hMatch_pull_numberOfLostHits->Fill(fabs(trkRefVal->numberOfLostHits() - trkRefKey.get()->numberOfLostHits()) / (float)trkRefKey.get()->numberOfLostHits());
      hMatch_pull_validFraction->Fill(fabs(trkRefVal->validFraction() - trkRefKey.get()->validFraction()) / trkRefKey.get()->validFraction());

      hMatch2d_dp_dr->Fill(trkRefVal->p() - trkRefKey.get()->p(), dR);
      hMatch2d_dpt_dr->Fill(trkRefVal->pt() - trkRefKey.get()->pt(), dR);
      hMatch2d_nSharedHits_dr->Fill(nSharedHits, dR);
      hMatch2d_dp_dz->Fill(trkRefVal->p() - trkRefKey.get()->p(), trkRefVal->dz() - trkRefKey.get()->dz());
      hMatch2d_dxy_dz->Fill(trkRefVal->dxy() - trkRefKey.get()->dxy(), trkRefVal->dz() - trkRefKey.get()->dz());

      hR_matched_p->Fill(trkRefKey.get()->p());
      hR_matched_pt->Fill(trkRefKey.get()->pt());
      hR_matched_eta->Fill(trkRefKey.get()->eta());
      hR_matched_phi->Fill(trkRefKey.get()->phi());
      hR_matched_ptError->Fill(trkRefKey.get()->ptError());
      hR_matched_dPtOverPt->Fill(trkRefKey.get()->ptError() / trkRefKey.get()->pt());
      hR_matched_etaError->Fill(trkRefKey.get()->etaError());
      hR_matched_phiError->Fill(trkRefKey.get()->phiError());
      hR_matched_outerP->Fill(trkRefKey.get()->outerP());
      hR_matched_outerPt->Fill(trkRefKey.get()->outerPt());
      hR_matched_outerEta->Fill(trkRefKey.get()->outerEta());
      hR_matched_outerPhi->Fill(trkRefKey.get()->outerPhi());
      hR_matched_outerPOverP->Fill(trkRefKey.get()->outerP() / trkRefKey.get()->p());
      hR_matched_outerPtOverPt->Fill(trkRefKey.get()->outerPt() / trkRefKey.get()->pt());
      hR_matched_chi2->Fill(trkRefKey.get()->chi2());
      hR_matched_ndof->Fill(trkRefKey.get()->ndof());
      hR_matched_normalizedChi2->Fill(trkRefKey.get()->normalizedChi2());
      hR_matched_charge->Fill(trkRefKey.get()->charge());
      hR_matched_dxy->Fill(trkRefKey.get()->dxy());
      hR_matched_dz->Fill(trkRefKey.get()->dz());
      hR_matched_dxyError->Fill(trkRefKey.get()->dxyError());
      hR_matched_dzError->Fill(trkRefKey.get()->dzError());
      hR_matched_numberOfValidHits->Fill(trkRefKey.get()->numberOfValidHits());
      hR_matched_numberOfLostHits->Fill(trkRefKey.get()->numberOfLostHits());
      hR_matched_validFraction->Fill(trkRefKey.get()->validFraction());

      hT_matched_p->Fill(trkRefVal->p());
      hT_matched_pt->Fill(trkRefVal->pt());
      hT_matched_eta->Fill(trkRefVal->eta());
      hT_matched_phi->Fill(trkRefVal->phi());
      hT_matched_ptError->Fill(trkRefVal->ptError());
      hT_matched_dPtOverPt->Fill(trkRefVal->ptError() / trkRefVal->pt());
      hT_matched_etaError->Fill(trkRefVal->etaError());
      hT_matched_phiError->Fill(trkRefVal->phiError());
      hT_matched_outerP->Fill(trkRefVal->outerP());
      hT_matched_outerPt->Fill(trkRefVal->outerPt());
      hT_matched_outerEta->Fill(trkRefVal->outerEta());
      hT_matched_outerPhi->Fill(trkRefVal->outerPhi());
      hT_matched_outerPOverP->Fill(trkRefVal->outerP() / trkRefVal->p());
      hT_matched_outerPtOverPt->Fill(trkRefVal->outerPt() / trkRefVal->pt());
      hT_matched_chi2->Fill(trkRefVal->chi2());
      hT_matched_ndof->Fill(trkRefVal->ndof());
      hT_matched_normalizedChi2->Fill(trkRefVal->normalizedChi2());
      hT_matched_charge->Fill(trkRefVal->charge());
      hT_matched_dxy->Fill(trkRefVal->dxy());
      hT_matched_dz->Fill(trkRefVal->dz());
      hT_matched_dxyError->Fill(trkRefVal->dxyError());
      hT_matched_dzError->Fill(trkRefVal->dzError());
      hT_matched_numberOfValidHits->Fill(trkRefVal->numberOfValidHits());
      hT_matched_numberOfLostHits->Fill(trkRefVal->numberOfLostHits());
      hT_matched_validFraction->Fill(trkRefVal->validFraction());

      if (trkRefKey.get()->pt() > minPtR_ && trkRefVal->pt() > minPtT_) {
         hMatch_ptCut_size->Fill(trkTrkMatchMap.size());
         hMatch_ptCut_dr->Fill(dR);
         hMatch_ptCut_nSharedHits->Fill(nSharedHits);
         hMatch_ptCut_nNonSharedHitsRef->Fill(trkRefKey.get()->numberOfValidHits() - nSharedHits);
         hMatch_ptCut_nNonSharedHitsTest->Fill(trkRefVal->numberOfValidHits() - nSharedHits);
         hMatch_ptCut_diff_p->Fill(trkRefVal->p() - trkRefKey.get()->p());
         hMatch_ptCut_diff_pt->Fill(trkRefVal->pt() - trkRefKey.get()->pt());
         hMatch_ptCut_diff_eta->Fill(trkRefVal->eta() - trkRefKey.get()->eta());
         hMatch_ptCut_diff_phi->Fill(phiDiff);
         hMatch_ptCut_diff_ptError->Fill(trkRefVal->ptError() - trkRefKey.get()->ptError());
         hMatch_ptCut_diff_dPtOverPt->Fill(trkRefVal->ptError() / trkRefVal->pt() - trkRefKey.get()->ptError() / trkRefKey.get()->pt());
         hMatch_ptCut_diff_etaError->Fill(trkRefVal->etaError() - trkRefKey.get()->etaError());
         hMatch_ptCut_diff_phiError->Fill(trkRefVal->phiError() - trkRefKey.get()->phiError());
         hMatch_ptCut_diff_outerP->Fill(trkRefVal->outerP() - trkRefKey.get()->outerP());
         hMatch_ptCut_diff_outerPt->Fill(trkRefVal->outerPt() - trkRefKey.get()->outerPt());
         hMatch_ptCut_diff_outerEta->Fill(trkRefVal->outerEta() - trkRefKey.get()->outerEta());
         hMatch_ptCut_diff_outerPhi->Fill(outerPhiDiff);
         hMatch_ptCut_diff_outerPOverP->Fill(trkRefVal->outerP() / trkRefVal->p() - trkRefKey.get()->outerP() / trkRefKey.get()->p());
         hMatch_ptCut_diff_outerPtOverPt->Fill(trkRefVal->outerPt() / trkRefVal->pt() - trkRefKey.get()->outerPt() / trkRefKey.get()->pt());
         hMatch_ptCut_diff_chi2->Fill(trkRefVal->chi2() - trkRefKey.get()->chi2());
         hMatch_ptCut_diff_ndof->Fill(trkRefVal->ndof() - trkRefKey.get()->ndof());
         hMatch_ptCut_diff_normalizedChi2->Fill(trkRefVal->normalizedChi2() - trkRefKey.get()->normalizedChi2());
         hMatch_ptCut_diff_charge->Fill(trkRefVal->charge() - trkRefKey.get()->charge());
         hMatch_ptCut_diff_dxy->Fill(trkRefVal->dxy() - trkRefKey.get()->dxy());
         hMatch_ptCut_diff_dz->Fill(trkRefVal->dz() - trkRefKey.get()->dz());
         hMatch_ptCut_diff_dxyError->Fill(trkRefVal->dxyError() - trkRefKey.get()->dxyError());
         hMatch_ptCut_diff_dzError->Fill(trkRefVal->dzError() - trkRefKey.get()->dzError());
         hMatch_ptCut_diff_numberOfValidHits->Fill(trkRefVal->numberOfValidHits() - trkRefKey.get()->numberOfValidHits());
         hMatch_ptCut_diff_numberOfLostHits->Fill(trkRefVal->numberOfLostHits() - trkRefKey.get()->numberOfLostHits());
         hMatch_ptCut_diff_validFraction->Fill(trkRefVal->validFraction() - trkRefKey.get()->validFraction());

         hMatch_ptCut_pull_p->Fill(fabs(trkRefVal->p() - trkRefKey.get()->p()) / trkRefKey.get()->p());
         hMatch_ptCut_pull_pt->Fill(fabs(trkRefVal->pt() - trkRefKey.get()->pt()) / trkRefKey.get()->pt());
         hMatch_ptCut_pull_eta->Fill(fabs(trkRefVal->eta() - trkRefKey.get()->eta()) / trkRefKey.get()->eta());
         hMatch_ptCut_pull_ptError->Fill(fabs(trkRefVal->ptError() - trkRefKey.get()->ptError()) / trkRefKey.get()->ptError());
         hMatch_ptCut_pull_dPtOverPt->Fill(fabs(trkRefVal->ptError() / trkRefVal->pt() - trkRefKey.get()->ptError() / trkRefKey.get()->pt()) / (trkRefKey.get()->ptError() / trkRefKey.get()->pt()));
         hMatch_ptCut_pull_etaError->Fill(fabs(trkRefVal->etaError() - trkRefKey.get()->etaError()) / trkRefKey.get()->etaError());
         hMatch_ptCut_pull_phiError->Fill(fabs(trkRefVal->phiError() - trkRefKey.get()->phiError()) / trkRefKey.get()->phiError());
         hMatch_ptCut_pull_outerP->Fill(fabs(trkRefVal->outerP() - trkRefKey.get()->outerP()) / trkRefKey.get()->outerP());
         hMatch_ptCut_pull_outerPt->Fill(fabs(trkRefVal->outerPt() - trkRefKey.get()->outerPt()) / trkRefKey.get()->outerPt());
         hMatch_ptCut_pull_outerEta->Fill(fabs(trkRefVal->outerEta() - trkRefKey.get()->outerEta()) / trkRefKey.get()->outerEta());
         hMatch_ptCut_pull_outerPOverP->Fill(fabs(trkRefVal->outerP() / trkRefVal->p() - trkRefKey.get()->outerP() / trkRefKey.get()->p()) / (trkRefKey.get()->outerP() / trkRefKey.get()->p()));
         hMatch_ptCut_pull_outerPtOverPt->Fill(fabs(trkRefVal->outerPt() / trkRefVal->pt() - trkRefKey.get()->outerPt() / trkRefKey.get()->pt()) / (trkRefKey.get()->outerPt() / trkRefKey.get()->pt()));
         hMatch_ptCut_pull_chi2->Fill(fabs(trkRefVal->chi2() - trkRefKey.get()->chi2()) / trkRefKey.get()->chi2());
         hMatch_ptCut_pull_ndof->Fill(fabs(trkRefVal->ndof() - trkRefKey.get()->ndof()) / (float)trkRefKey.get()->ndof());
         hMatch_ptCut_pull_normalizedChi2->Fill(fabs(trkRefVal->normalizedChi2() - trkRefKey.get()->normalizedChi2()) / trkRefKey.get()->normalizedChi2());
         hMatch_ptCut_pull_dxy->Fill(fabs(trkRefVal->dxy() - trkRefKey.get()->dxy()) / fabs(trkRefKey.get()->dxy()));
         hMatch_ptCut_pull_dz->Fill(fabs(trkRefVal->dz() - trkRefKey.get()->dz()) / fabs(trkRefKey.get()->dz()));
         hMatch_ptCut_pull_dxyError->Fill(fabs(trkRefVal->dxyError() - trkRefKey.get()->dxyError()) / trkRefKey.get()->dxyError());
         hMatch_ptCut_pull_dzError->Fill(fabs(trkRefVal->dzError() - trkRefKey.get()->dzError()) / trkRefKey.get()->dzError());
         if (trkRefKey.get()->numberOfValidHits() != 0) hMatch_ptCut_pull_numberOfValidHits->Fill(fabs(trkRefVal->numberOfValidHits() - trkRefKey.get()->numberOfValidHits()) / (float)trkRefKey.get()->numberOfValidHits());
         if (trkRefKey.get()->numberOfLostHits() != 0) hMatch_ptCut_pull_numberOfLostHits->Fill(fabs(trkRefVal->numberOfLostHits() - trkRefKey.get()->numberOfLostHits()) / (float)trkRefKey.get()->numberOfLostHits());
         hMatch_ptCut_pull_validFraction->Fill(fabs(trkRefVal->validFraction() - trkRefKey.get()->validFraction()) / trkRefKey.get()->validFraction());
      }
   }

   // Retrieving Deta and Dphi
   edm::InputTag dEtaMapTagR(referenceTrackVarsTag_.label(), "Deta", referenceTrackVarsTag_.process());
   edm::InputTag dPhiMapTagR(referenceTrackVarsTag_.label(), "Dphi", referenceTrackVarsTag_.process());
   edm::InputTag dEtaMapTagT(testTrackVarsTag_.label(), "Deta", testTrackVarsTag_.process());
   edm::InputTag dPhiMapTagT(testTrackVarsTag_.label(), "Dphi", testTrackVarsTag_.process());
 
   edm::Handle<reco::RecoEcalCandidateIsolationMap> ecalMatchHandleDEtaR;
   if (!iEvent.getByLabel(dEtaMapTagR, ecalMatchHandleDEtaR)) return;
   const reco::RecoEcalCandidateIsolationMap* ecalMatchMapDEtaR = ecalMatchHandleDEtaR.product();

   edm::Handle<reco::RecoEcalCandidateIsolationMap> ecalMatchHandleDPhiR;
   if (!iEvent.getByLabel(dPhiMapTagR, ecalMatchHandleDPhiR)) return;
   const reco::RecoEcalCandidateIsolationMap* ecalMatchMapDPhiR = ecalMatchHandleDPhiR.product();

   edm::Handle<reco::RecoEcalCandidateIsolationMap> ecalMatchHandleDEtaT;
   if (!iEvent.getByLabel(dEtaMapTagT, ecalMatchHandleDEtaT)) return;
   const reco::RecoEcalCandidateIsolationMap* ecalMatchMapDEtaT = ecalMatchHandleDEtaT.product();

   edm::Handle<reco::RecoEcalCandidateIsolationMap> ecalMatchHandleDPhiT;
   if (!iEvent.getByLabel(dPhiMapTagT, ecalMatchHandleDPhiT)) return;
   const reco::RecoEcalCandidateIsolationMap* ecalMatchMapDPhiT = ecalMatchHandleDPhiT.product();

   unsigned int goodMatch = 0;
   unsigned int nMatch = 0;
   reco::RecoEcalCandidateIsolationMap::const_iterator mapIt;
   //std::cout << "Map size: " << ecalMatchMapDEtaR->size() << std::endl;
   hR_dEta_mapSize->Fill(ecalMatchMapDEtaR->size());
   for (mapIt = ecalMatchMapDEtaR->begin(); mapIt != ecalMatchMapDEtaR->end(); ++mapIt) {
      float value = mapIt->val;
      //std::cout << "Deta: " <<  value << std::endl;
      hR_dEta_dEta->Fill(value);
      if (value < 999999.) ++goodMatch;
      if (value < dEtaCutR_) ++nMatch;
   }
   hR_dEta_mapSizeMatch->Fill(goodMatch);
   hR_dEta_nPassCut->Fill(nMatch);

   goodMatch = 0;
   nMatch = 0;
   hR_dPhi_mapSize->Fill(ecalMatchMapDPhiR->size());
   for (mapIt = ecalMatchMapDPhiR->begin(); mapIt != ecalMatchMapDPhiR->end(); ++mapIt) {
      float value = mapIt->val;
      //std::cout << "Dphi: " << value << std::endl;
      hR_dPhi_dPhi->Fill(value);
      if (value < 999999.) ++goodMatch;
      if (value < dPhiCutR_) ++nMatch;
   }
   hR_dPhi_mapSizeMatch->Fill(goodMatch);
   hR_dPhi_nPassCut->Fill(nMatch);

   goodMatch = 0;
   nMatch = 0;
   hT_dEta_mapSize->Fill(ecalMatchMapDEtaT->size());
   for (mapIt = ecalMatchMapDEtaT->begin(); mapIt != ecalMatchMapDEtaT->end(); ++mapIt) {
      float value = mapIt->val;
      //std::cout << "Deta: " <<  value << std::endl;
      hT_dEta_dEta->Fill(value);
      if (value < 999999.) ++goodMatch;
      if (value < dEtaCutT_) ++nMatch;
   }
   hT_dEta_mapSizeMatch->Fill(goodMatch);
   hT_dEta_nPassCut->Fill(nMatch);

   goodMatch = 0;
   nMatch = 0;
   hT_dPhi_mapSize->Fill(ecalMatchMapDPhiT->size());
   for (mapIt = ecalMatchMapDPhiT->begin(); mapIt != ecalMatchMapDPhiT->end(); ++mapIt) {
      float value = mapIt->val;
      //std::cout << "Dphi: " << value << std::endl;
      hT_dPhi_dPhi->Fill(value);
      if (value < 999999.) ++goodMatch;
      if (value < dPhiCutT_) ++nMatch;
   }
   hT_dPhi_mapSizeMatch->Fill(goodMatch);
   hT_dPhi_nPassCut->Fill(nMatch);

   // now look only at the ones that pass the previous filters
   // mimic behaviour of hltEgammaGenericFilter
   edm::Handle<trigger::TriggerFilterObjectWithRefs> prevFilterOutput;
   std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > recoEcalCands;
   edm::Ref<reco::RecoEcalCandidateCollection> ecalCandRef;
   if (iEvent.getByLabel (ecalCandTagForDEtaR_, prevFilterOutput)) {
      prevFilterOutput->getObjects(trigger::TriggerCluster, recoEcalCands);
      if(recoEcalCands.empty()) prevFilterOutput->getObjects(trigger::TriggerPhoton, recoEcalCands);
      goodMatch = 0;
      nMatch = 0;
      hR_dEta_recoEcalCandsSize->Fill(recoEcalCands.size());
      for (unsigned int i = 0; i < recoEcalCands.size(); ++i) {
         ecalCandRef = recoEcalCands[i];
         mapIt = ecalMatchMapDEtaR->find(ecalCandRef);
         if (mapIt == ecalMatchMapDEtaR->end()) continue;
         float value = mapIt->val;
         hR_dEta_dEta_trgdEcalCands->Fill(value);
         if (value < 999999.) ++goodMatch;
         if (value < dEtaCutR_) ++nMatch;
      }
      hR_dEta_recoEcalCandsSizeMatch->Fill(goodMatch);
      hR_dEta_nPassCut_trgdEcalCands->Fill(nMatch);
   }

   if (iEvent.getByLabel (ecalCandTagForDPhiR_, prevFilterOutput)) {
      prevFilterOutput->getObjects(trigger::TriggerCluster, recoEcalCands);
      if(recoEcalCands.empty()) prevFilterOutput->getObjects(trigger::TriggerPhoton, recoEcalCands);
      goodMatch = 0;
      nMatch = 0;
      hR_dPhi_recoEcalCandsSize->Fill(recoEcalCands.size());
      for (unsigned int i = 0; i < recoEcalCands.size(); ++i) {
         ecalCandRef = recoEcalCands[i];
         mapIt = ecalMatchMapDPhiR->find(ecalCandRef);
         if (mapIt == ecalMatchMapDPhiR->end()) continue;
         float value = mapIt->val;
         hR_dPhi_dPhi_trgdEcalCands->Fill(value);
         if (value < 999999.) ++goodMatch;
         if (value < dPhiCutR_) ++nMatch;
      }
      hR_dPhi_recoEcalCandsSizeMatch->Fill(goodMatch);
      hR_dPhi_nPassCut_trgdEcalCands->Fill(nMatch);
   }

   if (iEvent.getByLabel (ecalCandTagForDEtaT_, prevFilterOutput)) {
      prevFilterOutput->getObjects(trigger::TriggerCluster, recoEcalCands);
      if(recoEcalCands.empty()) prevFilterOutput->getObjects(trigger::TriggerPhoton, recoEcalCands);
      goodMatch = 0;
      nMatch = 0;
      hT_dEta_recoEcalCandsSize->Fill(recoEcalCands.size());
      for (unsigned int i = 0; i < recoEcalCands.size(); ++i) {
         ecalCandRef = recoEcalCands[i];
         mapIt = ecalMatchMapDEtaT->find(ecalCandRef);
         if (mapIt == ecalMatchMapDEtaT->end()) continue;
         float value = mapIt->val;
         hT_dEta_dEta_trgdEcalCands->Fill(value);
         if (value < 999999.) ++goodMatch;
         if (value < dEtaCutT_) ++nMatch;
      }
      hT_dEta_recoEcalCandsSizeMatch->Fill(goodMatch);
      hT_dEta_nPassCut_trgdEcalCands->Fill(nMatch);
   }

   if (iEvent.getByLabel (ecalCandTagForDPhiT_, prevFilterOutput)) {
      prevFilterOutput->getObjects(trigger::TriggerCluster, recoEcalCands);
      if(recoEcalCands.empty()) prevFilterOutput->getObjects(trigger::TriggerPhoton, recoEcalCands);
      goodMatch = 0;
      nMatch = 0;
      hT_dPhi_recoEcalCandsSize->Fill(recoEcalCands.size());
      for (unsigned int i = 0; i < recoEcalCands.size(); ++i) {
         ecalCandRef = recoEcalCands[i];
         mapIt = ecalMatchMapDPhiT->find(ecalCandRef);
         if (mapIt == ecalMatchMapDPhiT->end()) continue;
         float value = mapIt->val;
         hT_dPhi_dPhi_trgdEcalCands->Fill(value);
         if (value < 999999.) ++goodMatch;
         if (value < dPhiCutT_) ++nMatch;
      }
      hT_dPhi_recoEcalCandsSizeMatch->Fill(goodMatch);
      hT_dPhi_nPassCut_trgdEcalCands->Fill(nMatch);
   }

   //std::cout << std::endl;

   //tree->Fill();
}

double 
HltGsfTrkAnalyser::TrackMatcher(const reco::GsfTrack &t1, const reco::GsfTrack &t2) 
{
   double discr = -1.;

   if (matchMethod_ == 0) discr = reco::deltaR(t1, t2);
   else if (matchMethod_ == 1) discr = fabs(t1.pt()-t2.pt());
   else if (matchMethod_ == 2) discr = fabs(t1.dz()-t2.dz());
   else if (matchMethod_ == 3) discr = (double)GetNSharedHits(t1, t2);

   return discr;
}

unsigned int 
HltGsfTrkAnalyser::GetNSharedHits(const reco::GsfTrack &t1, const reco::GsfTrack &t2) 
{
   unsigned int nSharedHits = 0;
   trackingRecHit_iterator nhit=t1.recHitsBegin();
   trackingRecHit_iterator nhit_end=t1.recHitsEnd();

   for (; nhit != nhit_end; ++nhit) {
      if ((*nhit)->isValid()) {

         trackingRecHit_iterator ihit=t2.recHitsBegin();
         trackingRecHit_iterator ihit_end=t2.recHitsEnd();

         for (; ihit != ihit_end; ++ihit) {
            if ((*ihit)->isValid() && (*nhit)->sharesInput(&*(*ihit), TrackingRecHit::all)) {
               ++nSharedHits;
               break;
            }
         }
      }
   }
   //std::cout << "found " << nSharedHits << " shared hits" << std::endl;

   return nSharedHits;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HltGsfTrkAnalyser::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HltGsfTrkAnalyser::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
HltGsfTrkAnalyser::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HltGsfTrkAnalyser::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HltGsfTrkAnalyser::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HltGsfTrkAnalyser::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HltGsfTrkAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HltGsfTrkAnalyser);
