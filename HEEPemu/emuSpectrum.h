//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 16 17:19:37 2013 by ROOT version 5.32/00
// from TTree tree/tree
// found on file: /user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root
//////////////////////////////////////////////////////////

#ifndef emuSpectrum_h
#define emuSpectrum_h

#include "TFile.h"
#include "TDCacheFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TPaveLabel.h"
#include "TParameter.h"
#include "THashList.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TStopwatch.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxnWasRun = 1;
const Int_t kMaxnAccept = 1;
const Int_t kMaxnErrors = 1;

class EmuSpectrum {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          runnumber;
   UInt_t          eventnumber;
   UInt_t          luminosityBlock;
   Int_t           HLT_Mu22_Photon22_CaloIdL;
   Int_t           HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Mu22_Photon22_CaloIdL;
   Int_t           prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Float_t         rho;
   Float_t         pfmet;
   Int_t           pvsize;
   Int_t           JetColl_size;
   Float_t         Jet_pt[100];   //[JetColl_size]
   Int_t           muon_size;
   Float_t         muon_pt[100];   //[muon_size]
   Float_t         muon_ptError[100];   //[muon_size]
   Float_t         muon_eta[100];   //[muon_size]
   Float_t         muon_phi[100];   //[muon_size]
   Int_t           muon_charge[100];   //[muon_size]
   Int_t           muon_nhitspixel[100];   //[muon_size]
   Int_t           muon_nhitstrack[100];   //[muon_size]
   Int_t           muon_nhitsmuons[100];   //[muon_size]
   Int_t           muon_nlayerswithhits[100];   //[muon_size]
   Int_t           muon_nSegmentMatch[100];   //[muon_size]
   Bool_t          muon_isTrackerMuon[100];   //[muon_size]
   Float_t         muon_normChi2[100];   //[muon_size]
   Float_t         muon_dz_beamSpot[100];   //[muon_size]
   Float_t         muon_dz_firstPVtx[100];   //[muon_size]
   Float_t         muon_dxy_cmsCenter[100];   //[muon_size]
   Float_t         muon_dxy_beamSpot[100];   //[muon_size]
   Float_t         muon_dxy_firstPVtx[100];   //[muon_size]
   Float_t         muon_trackIso03[100];   //[muon_size]
   Float_t         muon_emIso03[100];   //[muon_size]
   Float_t         muon_hadIso03[100];   //[muon_size]
   Int_t           gsf_size;
   Float_t         gsf_eta[100];   //[gsf_size]
   Float_t         gsf_phi[100];   //[gsf_size]
   Float_t         gsf_theta[100];   //[gsf_size]
   Int_t           gsf_charge[100];   //[gsf_size]
   Float_t         gsf_sigmaetaeta[100];   //[gsf_size]
   Float_t         gsf_sigmaIetaIeta[100];   //[gsf_size]
   Float_t         gsf_dxy_firstPVtx[100];   //[gsf_size]
   Float_t         gsf_dz_beamSpot[100];   //[gsf_size]
   Float_t         gsf_dz_firstPVtx[100];   //[gsf_size]
   Int_t           gsf_nLostInnerHits[100];   //[gsf_size]
   Float_t         gsf_deltaeta[100];   //[gsf_size]
   Float_t         gsf_deltaphi[100];   //[gsf_size]
   Float_t         gsf_hovere[100];   //[gsf_size]
   Float_t         gsf_trackiso[100];   //[gsf_size]
   Float_t         gsf_ecaliso[100];   //[gsf_size]
   Float_t         gsf_hcaliso1[100];   //[gsf_size]
   Float_t         gsf_hcaliso2[100];   //[gsf_size]
   Bool_t          gsf_isecaldriven[100];   //[gsf_size]
   Float_t         gsfsc_e[100];   //[gsf_size]
   Float_t         gsfsc_eta[100];   //[gsf_size]
   Float_t         gsfsc_phi[100];   //[gsf_size]
   Float_t         gsf_e2x5overe5x5[100];   //[gsf_size]
   Float_t         gsf_e1x5overe5x5[100];   //[gsf_size]
   Float_t         gsf_gsfet[100];   //[gsf_size]
   Float_t         genPair_mass;
   Int_t           trueNVtx;

   // List of branches
   TBranch        *b_runnumber;   //!
   TBranch        *b_eventnumber;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_HLT_Mu22_Photon22_CaloIdL;   //!
   TBranch        *b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Mu22_Photon22_CaloIdL;   //!
   TBranch        *b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_pvsize;   //!
   TBranch        *b_JetColl_size;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_muon_size;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_ptError;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_nhitspixel;   //!
   TBranch        *b_muon_nhitstrack;   //!
   TBranch        *b_muon_nhitsmuons;   //!
   TBranch        *b_muon_nlayerswithhits;   //!
   TBranch        *b_muon_nSegmentMatch;   //!
   TBranch        *b_muon_isTrackerMuon;   //!
   TBranch        *b_muon_normChi2;   //!
   TBranch        *b_muon_dz_beamSpot;   //!
   TBranch        *b_muon_dz_firstPVtx;   //!
   TBranch        *b_muon_dxy_cmsCenter;   //!
   TBranch        *b_muon_dxy_beamSpot;   //!
   TBranch        *b_muon_dxy_firstPVtx;   //!
   TBranch        *b_muon_trackIso03;   //!
   TBranch        *b_muon_emIso03;   //!
   TBranch        *b_muon_hadIso03;   //!
   TBranch        *b_gsf_size;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_charge;   //!
   TBranch        *b_gsf_sigmaetaeta;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
   TBranch        *b_gsf_dz_beamSpot;   //!
   TBranch        *b_gsf_dz_firstPVtx;   //!
   TBranch        *b_gsf_nLostInnerHits;   //!
   TBranch        *b_gsf_deltaeta;   //!
   TBranch        *b_gsf_deltaphi;   //!
   TBranch        *b_gsf_hovere;   //!
   TBranch        *b_gsf_trackiso;   //!
   TBranch        *b_gsf_ecaliso;   //!
   TBranch        *b_gsf_hcaliso1;   //!
   TBranch        *b_gsf_hcaliso2;   //!
   TBranch        *b_gsf_isecaldriven;   //!
   TBranch        *b_gsfsc_e;   //!
   TBranch        *b_gsfsc_eta;   //!
   TBranch        *b_gsfsc_phi;   //!
   TBranch        *b_gsf_e2x5overe5x5;   //!
   TBranch        *b_gsf_e1x5overe5x5;   //!
   TBranch        *b_gsf_gsfet;   //!
   TBranch        *b_genPair_mass;   //!
   TBranch        *b_trueNVtx;   //!

   EmuSpectrum(TTree *tree=0);
   virtual ~EmuSpectrum();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

protected :
   pair<unsigned int, unsigned int> runs_HLT_Mu22_Photon22_CaloIdL;
   pair<unsigned int, unsigned int> runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   pair<unsigned int, unsigned int> runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;

   // HEEP v4.1
   // barrel
   float bar_et;
   float bar_hoE;
   float bar_DEta;
   float bar_DPhi;
   float bar_e2x5e5x5;
   float bar_e1x5e5x5;
   float bar_isoEcalHcal1_1;
   float bar_isoEcalHcal1_2;
   float bar_isoEcalHcalRho;
   float bar_isoTrack;
   float bar_dxy;
   int bar_missInnerHits;
   // endcap
   float end_et;
   float end_hoE;
   float end_DEta;
   float end_DPhi;
   float end_sigmaietaieta;
   float end_isoEcalHcal1_1_1;
   float end_isoEcalHcal1_1_2;
   float end_isoEcalHcal1_2;
   float end_isoEcalHcalRho;
   float end_isoTrack;
   float end_dxy;
   int end_missInnerHits;

   //MUON selection
   float muon_pt_min;
   float muon_dptOverPt;
   float muon_etaMax;
   int muon_nHitsMinGlobal;
   int muon_nHitsMinPixel;
   int muon_nHitsMinMuon;
   int muon_nLayersMin; 
   float muon_impactParamMaxXY;   // in cm
   //float muon_impactParamMaxZ;   // in cm , not used
   int muon_nSegMatchMin;
   float muon_relIsoCutMax;

   // strings for histogram names
   vector<TString> suffix;

   // Functions
   bool PassHEEP(const int &n);
   bool PassHighPtMu(const int &n);
   int Trigger(int &prescale, unsigned int *trig, const int &selector = 0);
};

#endif

#ifdef emuSpectrum_cxx
EmuSpectrum::EmuSpectrum(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
   //   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
   //   if (!f || !f->IsOpen()) {
   //      f = new TFile("/user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
   //   }
   //   TDirectory * dir = (TDirectory*)f->Get("/user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root:/gsfcheckerjob");
   //   dir->GetObject("tree",tree);

   //}
   } else {
      Init(tree);
   }
   runs_HLT_Mu22_Photon22_CaloIdL = make_pair(999999999, 0);
   runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = make_pair(999999999, 0);
   runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = make_pair(999999999, 0);
   
   // HEEP v4.1
   // barrel
   bar_et = 35.;
   bar_hoE = 0.05;
   bar_DEta = 0.005;
   bar_DPhi = 0.06;
   bar_e2x5e5x5 = 0.94;
   bar_e1x5e5x5 = 0.83;
   bar_isoEcalHcal1_1 = 2.;
   bar_isoEcalHcal1_2 = 0.03;
   bar_isoEcalHcalRho = 0.28;
   bar_isoTrack = 5.;
   bar_dxy = 0.02;
   bar_missInnerHits = 1;
   // endcap
   end_et = 35.;
   end_hoE = 0.05;
   end_DEta = 0.007 ;
   end_DPhi = 0.06;
   end_sigmaietaieta = 0.03;
   end_isoEcalHcal1_1_1 = 2.5;
   end_isoEcalHcal1_1_2 = 1.;
   end_isoEcalHcal1_2 = 0.03;
   end_isoEcalHcalRho = 0.28;
   end_isoTrack = 5.;
   end_dxy = 0.05;
   end_missInnerHits = 1;

   //MUON selection
   muon_pt_min = 35.;
   muon_dptOverPt = 0.3;
   muon_etaMax = 2.4;
   muon_nHitsMinGlobal = 0;
   muon_nHitsMinPixel = 1;
   muon_nHitsMinMuon = 1;
   muon_nLayersMin = 6; 
   muon_impactParamMaxXY = 0.2;   // in cm
   //muon_impactParamMaxZ = 0.5;   // in cm , not used
   muon_nSegMatchMin = 2;
   muon_relIsoCutMax = 0.1;

   suffix.push_back("data");
   suffix.push_back("ttbar");
   suffix.push_back("ttbar700to1000");
   suffix.push_back("ttbar1000up");
   suffix.push_back("ztautau");
   suffix.push_back("ww");
   suffix.push_back("wz");
   suffix.push_back("tw");
   suffix.push_back("wjets");
   suffix.push_back("zmumu");
   suffix.push_back("zee");
   suffix.push_back("zz");
   suffix.push_back("ttbarto2l");
   suffix.push_back("ttbarto1l1jet");
   suffix.push_back("ttbarw");
   suffix.push_back("ttbarww");
}

EmuSpectrum::~EmuSpectrum()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EmuSpectrum::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EmuSpectrum::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EmuSpectrum::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("eventnumber", &eventnumber, &b_eventnumber);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &HLT_Mu22_Photon22_CaloIdL, &b_HLT_Mu22_Photon22_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Mu22_Photon22_CaloIdL", &prescale_HLT_Mu22_Photon22_CaloIdL, &b_prescale_HLT_Mu22_Photon22_CaloIdL);
   fChain->SetBranchAddress("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("pvsize", &pvsize, &b_pvsize);
   fChain->SetBranchAddress("JetColl_size", &JetColl_size, &b_JetColl_size);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("muon_size", &muon_size, &b_muon_size);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_ptError", muon_ptError, &b_muon_ptError);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_nhitspixel", muon_nhitspixel, &b_muon_nhitspixel);
   fChain->SetBranchAddress("muon_nhitstrack", muon_nhitstrack, &b_muon_nhitstrack);
   fChain->SetBranchAddress("muon_nhitsmuons", muon_nhitsmuons, &b_muon_nhitsmuons);
   fChain->SetBranchAddress("muon_nlayerswithhits", muon_nlayerswithhits, &b_muon_nlayerswithhits);
   fChain->SetBranchAddress("muon_nSegmentMatch", muon_nSegmentMatch, &b_muon_nSegmentMatch);
   fChain->SetBranchAddress("muon_isTrackerMuon", muon_isTrackerMuon, &b_muon_isTrackerMuon);
   fChain->SetBranchAddress("muon_normChi2", muon_normChi2, &b_muon_normChi2);
   fChain->SetBranchAddress("muon_dz_beamSpot", muon_dz_beamSpot, &b_muon_dz_beamSpot);
   fChain->SetBranchAddress("muon_dz_firstPVtx", muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
   fChain->SetBranchAddress("muon_dxy_cmsCenter", muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
   fChain->SetBranchAddress("muon_dxy_beamSpot", muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
   fChain->SetBranchAddress("muon_dxy_firstPVtx", muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
   fChain->SetBranchAddress("muon_trackIso03", muon_trackIso03, &b_muon_trackIso03);
   fChain->SetBranchAddress("muon_emIso03", muon_emIso03, &b_muon_emIso03);
   fChain->SetBranchAddress("muon_hadIso03", muon_hadIso03, &b_muon_hadIso03);
   fChain->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
   fChain->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_charge", gsf_charge, &b_gsf_charge);
   fChain->SetBranchAddress("gsf_sigmaetaeta", gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   fChain->SetBranchAddress("gsf_dz_beamSpot", gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
   fChain->SetBranchAddress("gsf_dz_firstPVtx", gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
   fChain->SetBranchAddress("gsf_nLostInnerHits", gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   fChain->SetBranchAddress("gsf_deltaeta", gsf_deltaeta, &b_gsf_deltaeta);
   fChain->SetBranchAddress("gsf_deltaphi", gsf_deltaphi, &b_gsf_deltaphi);
   fChain->SetBranchAddress("gsf_hovere", gsf_hovere, &b_gsf_hovere);
   fChain->SetBranchAddress("gsf_trackiso", gsf_trackiso, &b_gsf_trackiso);
   fChain->SetBranchAddress("gsf_ecaliso", gsf_ecaliso, &b_gsf_ecaliso);
   fChain->SetBranchAddress("gsf_hcaliso1", gsf_hcaliso1, &b_gsf_hcaliso1);
   fChain->SetBranchAddress("gsf_hcaliso2", gsf_hcaliso2, &b_gsf_hcaliso2);
   fChain->SetBranchAddress("gsf_isecaldriven", gsf_isecaldriven, &b_gsf_isecaldriven);
   fChain->SetBranchAddress("gsfsc_e", gsfsc_e, &b_gsfsc_e);
   fChain->SetBranchAddress("gsfsc_eta", gsfsc_eta, &b_gsfsc_eta);
   fChain->SetBranchAddress("gsfsc_phi", gsfsc_phi, &b_gsfsc_phi);
   fChain->SetBranchAddress("gsf_e2x5overe5x5", gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
   fChain->SetBranchAddress("gsf_e1x5overe5x5", gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
   fChain->SetBranchAddress("gsf_gsfet", gsf_gsfet, &b_gsf_gsfet);
   fChain->SetBranchAddress("genPair_mass", &genPair_mass, &b_genPair_mass);
   fChain->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);

   // enable only used branches
   fChain->SetBranchStatus("*", 0);
   fChain->SetBranchStatus("runnumber", 1);
   fChain->SetBranchStatus("eventnumber", 1);
   fChain->SetBranchStatus("luminosityBlock", 1);
   fChain->SetBranchStatus("HLT_Mu22_Photon22_CaloIdL", 1);
   fChain->SetBranchStatus("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
   fChain->SetBranchStatus("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
   fChain->SetBranchStatus("prescale_HLT_Mu22_Photon22_CaloIdL", 1);
   fChain->SetBranchStatus("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
   fChain->SetBranchStatus("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
   fChain->SetBranchStatus("rho", 1);
   fChain->SetBranchStatus("pfmet", 1);
   fChain->SetBranchStatus("pvsize", 1);
   fChain->SetBranchStatus("JetColl_size", 1);
   fChain->SetBranchStatus("Jet_pt", 1);
   fChain->SetBranchStatus("muon_size", 1);
   fChain->SetBranchStatus("muon_pt", 1);
   fChain->SetBranchStatus("muon_ptError", 1);
   fChain->SetBranchStatus("muon_eta", 1);
   fChain->SetBranchStatus("muon_phi", 1);
   fChain->SetBranchStatus("muon_charge", 1);
   fChain->SetBranchStatus("muon_nhitspixel", 1);
   fChain->SetBranchStatus("muon_nhitstrack", 1);
   fChain->SetBranchStatus("muon_nhitsmuons", 1);
   fChain->SetBranchStatus("muon_nlayerswithhits", 1);
   fChain->SetBranchStatus("muon_nSegmentMatch", 1);
   fChain->SetBranchStatus("muon_isTrackerMuon", 1);
   fChain->SetBranchStatus("muon_normChi2", 1);
   fChain->SetBranchStatus("muon_dz_beamSpot", 1);
   fChain->SetBranchStatus("muon_dz_firstPVtx", 1);
   fChain->SetBranchStatus("muon_dxy_cmsCenter", 1);
   fChain->SetBranchStatus("muon_dxy_beamSpot", 1);
   fChain->SetBranchStatus("muon_dxy_firstPVtx", 1);
   fChain->SetBranchStatus("muon_trackIso03", 1);
   fChain->SetBranchStatus("muon_emIso03", 1);
   fChain->SetBranchStatus("muon_hadIso03", 1);
   fChain->SetBranchStatus("gsf_size", 1);
   fChain->SetBranchStatus("gsf_eta", 1);
   fChain->SetBranchStatus("gsf_phi", 1);
   fChain->SetBranchStatus("gsf_theta", 1);
   fChain->SetBranchStatus("gsf_charge", 1);
   fChain->SetBranchStatus("gsf_sigmaetaeta", 1);
   fChain->SetBranchStatus("gsf_sigmaIetaIeta", 1);
   fChain->SetBranchStatus("gsf_dxy_firstPVtx", 1);
   fChain->SetBranchStatus("gsf_dz_beamSpot", 1);
   fChain->SetBranchStatus("gsf_dz_firstPVtx", 1);
   fChain->SetBranchStatus("gsf_nLostInnerHits", 1);
   fChain->SetBranchStatus("gsf_deltaeta", 1);
   fChain->SetBranchStatus("gsf_deltaphi", 1);
   fChain->SetBranchStatus("gsf_hovere", 1);
   fChain->SetBranchStatus("gsf_trackiso", 1);
   fChain->SetBranchStatus("gsf_ecaliso", 1);
   fChain->SetBranchStatus("gsf_hcaliso1", 1);
   fChain->SetBranchStatus("gsf_hcaliso2", 1);
   fChain->SetBranchStatus("gsf_isecaldriven", 1);
   fChain->SetBranchStatus("gsfsc_e", 1);
   fChain->SetBranchStatus("gsfsc_eta", 1);
   fChain->SetBranchStatus("gsfsc_phi", 1);
   fChain->SetBranchStatus("gsf_e2x5overe5x5", 1);
   fChain->SetBranchStatus("gsf_e1x5overe5x5", 1);
   fChain->SetBranchStatus("gsf_gsfet", 1);
   fChain->SetBranchStatus("genPair_mass", 1);
   fChain->SetBranchStatus("trueNVtx", 1);

   Notify();
}

Bool_t EmuSpectrum::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EmuSpectrum::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EmuSpectrum::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef emuSpectrum_cxx
