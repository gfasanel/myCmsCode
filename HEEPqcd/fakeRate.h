#ifndef fakeRate_h
#define fakeRate_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <vector>
#include <stdio.h>

#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"

class FakeRate {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TFile* _outF;
   virtual bool IsVBTF(float pt,float eta, float dEta,float dPhi, float tckiso, float ecaliso, float hcaliso1,float hcaliso2,float sigmaietaieta, float hovere); 
   virtual bool IsHEEP(float pt,float eta, float dEta,float dPhi, float tckiso, float ecaliso, float hcaliso1,float hcaliso2,float sigmaietaieta, float hovere, float  e2x5overe5x5,float e1x5overe5x5 );
   virtual void SingleEff(vector<TH1F *> histoB, vector<TH1F *> histobisB, vector<TH1F *> histoE, vector<TH1F *> histobisE,TString canvasname, TString Title, TString XTitle,TString YTitle );
   // virtual void Nmin1Eff(vector<TH1F *> histoB ,vector<TH1F *> histobisB,vector<TH1F *> histoE,vector<TH1F *> histobisE,TString canvasname, TString Title, TString XTitle,TString YTitle );

   // Declaration of leaf types
   Int_t           L1trigger_size;
   Int_t           L1trigger_bool[64];   //[L1trigger_size]
   Int_t           PhysDecl_bool;
   Int_t           HLTriggers[200];
   Int_t           nJetsAKT_pt15;
   Int_t           nJetsIC5_pt15;
   Float_t         calomet;
   Float_t         met;
   Int_t           jetIC5_size;
   Float_t         jetIC5_pt[100];   //[jetIC5_size]
   Float_t         jetIC5_eta[100];   //[jetIC5_size]
   Float_t         jetIC5_phi[100];   //[jetIC5_size]
   Float_t         jetIC5_em[100];   //[jetIC5_size]
   Int_t           muon_size;
   Float_t         muon_pt[100];   //[muon_size]
   Float_t         muon_ptError[100];   //[muon_size]
   Float_t         muon_eta[100];   //[muon_size]
   Float_t         muon_etaError[100];   //[muon_size]
   Float_t         muon_theta[100];   //[muon_size]
   Float_t         muon_thetaError[100];   //[muon_size]
   Float_t         muon_phi[100];   //[muon_size]
   Float_t         muon_phiError[100];   //[muon_size]
   Float_t         muon_outerPt[100];   //[muon_size]
   Float_t         muon_outerEta[100];   //[muon_size]
   Float_t         muon_outerPhi[100];   //[muon_size]
   Float_t         muon_outerTheta[100];   //[muon_size]
   Float_t         muon_px[100];   //[muon_size]
   Float_t         muon_py[100];   //[muon_size]
   Float_t         muon_pz[100];   //[muon_size]
   Int_t           muon_charge[100];   //[muon_size]
   Int_t           muon_nhitstotal[100];   //[muon_size]
   Int_t           muon_nhitstrack[100];   //[muon_size]
   Int_t           muon_nhitsmuons[100];   //[muon_size]
   Int_t           muon_nlosthits[100];   //[muon_size]
   Float_t         muon_chi2[100];   //[muon_size]
   Int_t           muon_ndof[100];   //[muon_size]
   Float_t         muon_normChi2[100];   //[muon_size]
   Float_t         muon_d0[100];   //[muon_size]
   Float_t         muon_d0Error[100];   //[muon_size]
   Float_t         muon_dz_cmsCenter[100];   //[muon_size]
   Float_t         muon_dz_beamSpot[100];   //[muon_size]
   Float_t         muon_dz_firstPVtx[100];   //[muon_size]
   Float_t         muon_dzError[100];   //[muon_size]
   Float_t         muon_dxy_cmsCenter[100];   //[muon_size]
   Float_t         muon_dxy_beamSpot[100];   //[muon_size]
   Float_t         muon_dxy_firstPVtx[100];   //[muon_size]
   Float_t         muon_dxyError[100];   //[muon_size]
   Float_t         muon_innerPosx[100];   //[muon_size]
   Float_t         muon_innerPosy[100];   //[muon_size]
   Float_t         muon_innerPosz[100];   //[muon_size]
   Float_t         muon_trackIso03[100];   //[muon_size]
   Float_t         muon_trackIso05[100];   //[muon_size]
   Float_t         muon_trackIso03_ptInVeto[100];   //[muon_size]
   Float_t         muon_trackIso05_ptInVeto[100];   //[muon_size]
   Int_t           runnumber;
   Int_t           eventnumber;
   Int_t           eventcounter;
   Int_t           processid;
   Float_t         pthat;
   Float_t         alphaqcd;
   Float_t         alphaqed;
   Float_t         qscale;
   Float_t         weight;
   Float_t         sigmaZ;
   Float_t         sigmaZ0Error;
   Float_t         sq;
   Float_t         bsposx;
   Float_t         bsposy;
   Float_t         bsposz;
   Int_t           pvsize;
   Float_t         pvx[10];   //[pvsize]
   Float_t         pvy[10];   //[pvsize]
   Float_t         pvz[10];   //[pvsize]
   Int_t           scsize;
   Float_t         scenergy[100];   //[scsize]
   Float_t         sceta[100];   //[scsize]
   Float_t         scetacorr[100];   //[scsize]
   Float_t         sctheta[100];   //[scsize]
   Float_t         scthetacorr[100];   //[scsize]
   Float_t         scet[100];   //[scsize]
   Float_t         scphi[100];   //[scsize]
   Float_t         scpx[100];   //[scsize]
   Float_t         scpy[100];   //[scsize]
   Float_t         scpz[100];   //[scsize]
   Float_t         scx[100];   //[scsize]
   Float_t         scy[100];   //[scsize]
   Float_t         scz[100];   //[scsize]
   Float_t         scgsfmatched[100];   //[scsize]
   Float_t         genelec_e_branch;
   Float_t         genelec_eta_branch;
   Float_t         genelec_phi_branch;
   Float_t         genelec_et_branch;
   Float_t         genposi_e_branch;
   Float_t         genposi_eta_branch;
   Float_t         genposi_phi_branch;
   Float_t         genposi_et_branch;
   Int_t           genelec_hassc_branch;
   Int_t           genposi_hassc_branch;
   Float_t         unstablegenelec_e_branch;
   Float_t         unstablegenelec_eta_branch;
   Float_t         unstablegenelec_phi_branch;
   Float_t         unstablegenelec_et_branch;
   Float_t         unstablegenposi_e_branch;
   Float_t         unstablegenposi_eta_branch;
   Float_t         unstablegenposi_phi_branch;
   Float_t         unstablegenposi_et_branch;
   Float_t         genboson_m_branch;
   Float_t         genboson_eta_branch;
   Float_t         genboson_phi_branch;
   Float_t         genboson_e_branch;
   Float_t         genboson_et_branch;
   Float_t         genboson_ez_branch;
   Float_t         genboson_p_branch;
   Float_t         genboson_pt_branch;
   Float_t         genboson_pz_branch;
   Float_t         x1quark;
   Float_t         x2quark;
   Int_t           fsrposiphotonsize;
   Int_t           fsrelecphotonsize;
   Float_t         energyfsrelec[20];   //[fsrelecphotonsize]
   Float_t         etfsrelec[20];   //[fsrelecphotonsize]
   Float_t         etafsrelec[20];   //[fsrelecphotonsize]
   Float_t         phifsrelec[20];   //[fsrelecphotonsize]
   Float_t         energyfsrposi[20];   //[fsrposiphotonsize]
   Float_t         etfsrposi[20];   //[fsrposiphotonsize]
   Float_t         etafsrposi[20];   //[fsrposiphotonsize]
   Float_t         phifsrposi[20];   //[fsrposiphotonsize]
   Float_t         scelecenergy;
   Float_t         sceleceta;
   Float_t         scelecphi;
   Float_t         scelecgsfmatched;
   Float_t         scposienergy;
   Float_t         scposieta;
   Float_t         scposiphi;
   Float_t         scposigsfmatched;
   Bool_t          genelechassc;
   Bool_t          genposihassc;
   Int_t           gsf_size;
   Float_t         gsf_theta[100];   //[gsf_size]
   Int_t           gsf_isEB[100];   //[gsf_size]
   Int_t           gsf_isEE[100];   //[gsf_size]
   Float_t         gsf_deltaEtaATcalo[100];   //[gsf_size]
   Float_t         gsf_deltaPhiATcalo[100];   //[gsf_size]
   Float_t         gsf_ecalEnergy[100];   //[gsf_size]
   Float_t         gsf_eOVERp[100];   //[gsf_size]
   Float_t         gsf_dxy[100];   //[gsf_size]
   Float_t         gsf_vz[100];   //[gsf_size]
   Int_t           gsf_nHits[100];   //[gsf_size]
   Float_t         gsf_fBrem[100];   //[gsf_size]
   Float_t         gsf_e1x5[100];   //[gsf_size]
   Float_t         gsf_e2x5[100];   //[gsf_size]
   Float_t         gsf_e5x5[100];   //[gsf_size]
   Float_t         gsf_eMax[100];   //[gsf_size]
   Float_t         gsf_SwissCross[100];   //[gsf_size]
   Float_t         gsf_e1x3[100];   //[gsf_size]
   Float_t         gsf_e3x1[100];   //[gsf_size]
   Float_t         gsf_e2x2[100];   //[gsf_size]
   Float_t         gsf_e3x2[100];   //[gsf_size]
   Float_t         gsf_e3x3[100];   //[gsf_size]
   Float_t         gsf_e4x4[100];   //[gsf_size]
   Float_t         gsf_e2x5Right[100];   //[gsf_size]
   Float_t         gsf_e2x5Left[100];   //[gsf_size]
   Float_t         gsf_e2x5Top[100];   //[gsf_size]
   Float_t         gsf_e2x5Bottom[100];   //[gsf_size]
   Float_t         gsf_e2x5Max[100];   //[gsf_size]
   Float_t         gsf_eLeft[100];   //[gsf_size]
   Float_t         gsf_eRight[100];   //[gsf_size]
   Float_t         gsf_eTop[100];   //[gsf_size]
   Float_t         gsf_eBottom[100];   //[gsf_size]
   Float_t         gsf_e2nd[100];   //[gsf_size]
   Float_t         gsf_p[100];   //[gsf_size]
   Float_t         gsf_e[100];   //[gsf_size]
   Float_t         gsf_pt[100];   //[gsf_size]
   Float_t         gsf_class[100];   //[gsf_size]
   Float_t         gsf_e2x5overe5x5[100];   //[gsf_size]
   Float_t         gsf_e1x5overe5x5[100];   //[gsf_size]
   Float_t         gsf_eta[100];   //[gsf_size]
   Float_t         gsf_phi[100];   //[gsf_size]
   Float_t         gsf_px[100];   //[gsf_size]
   Float_t         gsf_py[100];   //[gsf_size]
   Float_t         gsf_pz[100];   //[gsf_size]
   Float_t         gsf_deltaeta[100];   //[gsf_size]
   Float_t         gsf_deltaphi[100];   //[gsf_size]
   Float_t         gsf_hovere[100];   //[gsf_size]
   Float_t         gsf_trackiso[100];   //[gsf_size]
   Float_t         gsf_ecaliso[100];   //[gsf_size]
   Float_t         gsf_hcaliso1[100];   //[gsf_size]
   Float_t         gsf_hcaliso2[100];   //[gsf_size]
   Int_t           gsf_charge[100];   //[gsf_size]
   Float_t         gsf_sigmaetaeta[100];   //[gsf_size]
   Float_t         gsf_sigmaIetaIeta[100];   //[gsf_size]
   Float_t         gsf_isecaldriven[100];   //[gsf_size]
   Float_t         gsf_istrackerdriven[100];   //[gsf_size]
   Float_t         gsfsc_e[100];   //[gsf_size]
   Float_t         gsfsc_pt[100];   //[gsf_size]
   Float_t         gsfsc_eta[100];   //[gsf_size]
   Float_t         gsfsc_phi[100];   //[gsf_size]
   Float_t         gsfsc_px[100];   //[gsf_size]
   Float_t         gsfsc_py[100];   //[gsf_size]
   Float_t         gsfsc_pz[100];   //[gsf_size]
   Float_t         gsf_gsfet[100];   //[gsf_size]
   Int_t           scindexforgsf[100];   //[gsf_size]
   Int_t           gsfindexforgenelec;
   Int_t           gsfindexforgenposi;
   Int_t           scindexforgenelec;
   Int_t           scindexforgenposi;
   Bool_t          gsfpass_ET[100];   //[gsf_size]
   Bool_t          gsfpass_PT[100];   //[gsf_size]
   Bool_t          gsfpass_DETETA[100];   //[gsf_size]
   Bool_t          gsfpass_CRACK[100];   //[gsf_size]
   Bool_t          gsfpass_DETAIN[100];   //[gsf_size]
   Bool_t          gsfpass_DPHIIN[100];   //[gsf_size]
   Bool_t          gsfpass_HADEM[100];   //[gsf_size]
   Bool_t          gsfpass_SIGMAIETAIETA[100];   //[gsf_size]
   Bool_t          gsfpass_E2X5OVER5X5[100];   //[gsf_size]
   Bool_t          gsfpass_ISOLEMHADDEPTH1[100];   //[gsf_size]
   Bool_t          gsfpass_ISOLHADDEPTH2[100];   //[gsf_size]
   Bool_t          gsfpass_ISOLPTTRKS[100];   //[gsf_size]
   Bool_t          gsfpass_ECALDRIVEN[100];   //[gsf_size]
   Bool_t          gsfpass_INVALID[100];   //[gsf_size]
   Bool_t          gsfpass_HEEP[100];   //[gsf_size]
   Bool_t          gsfpass_ID[100];   //[gsf_size]
   Bool_t          gsfpass_ISO[100];   //[gsf_size]
   Int_t           scpixcharge[100];   //[gsf_size]
   Int_t           ctfcharge[100];   //[gsf_size]
   Int_t           gsfcharge[100];   //[gsf_size]
   Bool_t          gsfctfscpixconsistent[100];   //[gsf_size]
   Bool_t          gsfscpixconsistent[100];   //[gsf_size]
   Bool_t          gsfctfconsistent[100];   //[gsf_size]
   Int_t           gsftracksize;
   Float_t         gsftracketa[100];   //[gsftracksize]
   Float_t         gsftrackphi[100];   //[gsftracksize]
   Float_t         gsftrackp[100];   //[gsftracksize]
   Float_t         gsftrackpt[100];   //[gsftracksize]
   Float_t         gsftrackpx[100];   //[gsftracksize]
   Float_t         gsftrackpy[100];   //[gsftracksize]
   Float_t         gsftrackpz[100];   //[gsftracksize]

   // List of branches
   TBranch        *b_L1trigger_size;   //!
   TBranch        *b_L1trigger_bool;   //!
   TBranch        *b_PhysDecl_bool;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_nJetsAKT_pt15;   //!
   TBranch        *b_nJetsIC5_pt15;   //!
   TBranch        *b_calomet;   //!
   TBranch        *b_met;   //!
   TBranch        *b_jetIC5_size;   //!
   TBranch        *b_jetIC5_pt;   //!
   TBranch        *b_jetIC5_eta;   //!
   TBranch        *b_jetIC5_phi;   //!
   TBranch        *b_jetIC5_em;   //!
   TBranch        *b_muon_size;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_ptError;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_etaError;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaError;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_phiError;   //!
   TBranch        *b_muon_outerPt;   //!
   TBranch        *b_muon_outerEta;   //!
   TBranch        *b_muon_outerPhi;   //!
   TBranch        *b_muon_outerTheta;   //!
   TBranch        *b_muon_px;   //!
   TBranch        *b_muon_py;   //!
   TBranch        *b_muon_pz;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_nhitstotal;   //!
   TBranch        *b_muon_nhitstrack;   //!
   TBranch        *b_muon_nhitsmuons;   //!
   TBranch        *b_muon_nlosthits;   //!
   TBranch        *b_muon_chi2;   //!
   TBranch        *b_muon_ndof;   //!
   TBranch        *b_muon_normChi2;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_muon_d0Error;   //!
   TBranch        *b_muon_dz_cmsCenter;   //!
   TBranch        *b_muon_dz_beamSpot;   //!
   TBranch        *b_muon_dz_firstPVtx;   //!
   TBranch        *b_muon_dzError;   //!
   TBranch        *b_muon_dxy_cmsCenter;   //!
   TBranch        *b_muon_dxy_beamSpot;   //!
   TBranch        *b_muon_dxy_firstPVtx;   //!
   TBranch        *b_muon_dxyError;   //!
   TBranch        *b_muon_innerPosx;   //!
   TBranch        *b_muon_innerPosy;   //!
   TBranch        *b_muon_innerPosz;   //!
   TBranch        *b_muon_trackIso03;   //!
   TBranch        *b_muon_trackIso05;   //!
   TBranch        *b_muon_trackIso03_ptInVeto;   //!
   TBranch        *b_muon_trackIso05_ptInVeto;   //!
   TBranch        *b_runnumber;   //!
   TBranch        *b_eventnumber;   //!
   TBranch        *b_eventcounter;   //!
   TBranch        *b_processid;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_alphaqcd;   //!
   TBranch        *b_alphaqed;   //!
   TBranch        *b_qscale;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_sigmaZ;   //!
   TBranch        *b_sigmaZ0Error;   //!
   TBranch        *b_sq;   //!
   TBranch        *b_bsposx;   //!
   TBranch        *b_bsposy;   //!
   TBranch        *b_bsposz;   //!
   TBranch        *b_pvsize;   //!
   TBranch        *b_pvx;   //!
   TBranch        *b_pvy;   //!
   TBranch        *b_pvz;   //!
   TBranch        *b_scsize;   //!
   TBranch        *b_scenergy;   //!
   TBranch        *b_sceta;   //!
   TBranch        *b_scetacorr;   //!
   TBranch        *b_sctheta;   //!
   TBranch        *b_scthetacorr;   //!
   TBranch        *b_scet;   //!
   TBranch        *b_scphi;   //!
   TBranch        *b_scpx;   //!
   TBranch        *b_scpy;   //!
   TBranch        *b_scpz;   //!
   TBranch        *b_scx;   //!
   TBranch        *b_scy;   //!
   TBranch        *b_scz;   //!
   TBranch        *b_scgsfmatched;   //!
   TBranch        *b_genelec_e_branch;   //!
   TBranch        *b_genelec_eta_branch;   //!
   TBranch        *b_genelec_phi_branch;   //!
   TBranch        *b_genelec_et_branch;   //!
   TBranch        *b_genposi_e_branch;   //!
   TBranch        *b_genposi_eta_branch;   //!
   TBranch        *b_genposi_phi_branch;   //!
   TBranch        *b_genposi_et_branch;   //!
   TBranch        *b_genelec_hassc_branch;   //!
   TBranch        *b_genposi_hassc_branch;   //!
   TBranch        *b_unstablegenelec_e_branch;   //!
   TBranch        *b_unstablegenelec_eta_branch;   //!
   TBranch        *b_unstablegenelec_phi_branch;   //!
   TBranch        *b_unstablegenelec_et_branch;   //!
   TBranch        *b_unstablegenposi_e_branch;   //!
   TBranch        *b_unstablegenposi_eta_branch;   //!
   TBranch        *b_unstablegenposi_phi_branch;   //!
   TBranch        *b_unstablegenposi_et_branch;   //!
   TBranch        *b_genboson_m_branch;   //!
   TBranch        *b_genboson_eta_branch;   //!
   TBranch        *b_genboson_phi_branch;   //!
   TBranch        *b_genboson_e_branch;   //!
   TBranch        *b_genboson_et_branch;   //!
   TBranch        *b_genboson_ez_branch;   //!
   TBranch        *b_genboson_p_branch;   //!
   TBranch        *b_genboson_pt_branch;   //!
   TBranch        *b_genboson_pz_branch;   //!
   TBranch        *b_x1quark;   //!
   TBranch        *b_x2quark;   //!
   TBranch        *b_fsrposiphotonsize;   //!
   TBranch        *b_fsrelecphotonsize;   //!
   TBranch        *b_energyfsrelec;   //!
   TBranch        *b_etfsrelec;   //!
   TBranch        *b_etafsrelec;   //!
   TBranch        *b_phifsrelec;   //!
   TBranch        *b_energyfsrposi;   //!
   TBranch        *b_etfsrposi;   //!
   TBranch        *b_etafsrposi;   //!
   TBranch        *b_phifsrposi;   //!
   TBranch        *b_scelecenergy;   //!
   TBranch        *b_sceleceta;   //!
   TBranch        *b_scelecphi;   //!
   TBranch        *b_scelecgsfmatched;   //!
   TBranch        *b_scposienergy;   //!
   TBranch        *b_scposieta;   //!
   TBranch        *b_scposiphi;   //!
   TBranch        *b_scposigsfmatched;   //!
   TBranch        *b_genelechassc;   //!
   TBranch        *b_genposihassc;   //!
   TBranch        *b_gsf_size;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_isEB;   //!
   TBranch        *b_gsf_isEE;   //!
   TBranch        *b_gsf_deltaEtaATcalo;   //!
   TBranch        *b_gsf_deltaPhiATcalo;   //!
   TBranch        *b_gsf_ecalEnergy;   //!
   TBranch        *b_gsf_eOVERp;   //!
   TBranch        *b_gsf_dxy;   //!
   TBranch        *b_gsf_vz;   //!
   TBranch        *b_gsf_nHits;   //!
   TBranch        *b_gsf_fBrem;   //!
   TBranch        *b_gsf_e1x5;   //!
   TBranch        *b_gsf_e2x5;   //!
   TBranch        *b_gsf_e5x5;   //!
   TBranch        *b_gsf_eMax;   //!
   TBranch        *b_gsf_SwissCross;   //!
   TBranch        *b_gsf_e1x3;   //!
   TBranch        *b_gsf_e3x1;   //!
   TBranch        *b_gsf_e2x2;   //!
   TBranch        *b_gsf_e3x2;   //!
   TBranch        *b_gsf_e3x3;   //!
   TBranch        *b_gsf_e4x4;   //!
   TBranch        *b_gsf_e2x5Right;   //!
   TBranch        *b_gsf_e2x5Left;   //!
   TBranch        *b_gsf_e2x5Top;   //!
   TBranch        *b_gsf_e2x5Bottom;   //!
   TBranch        *b_gsf_e2x5Max;   //!
   TBranch        *b_gsf_eLeft;   //!
   TBranch        *b_gsf_eRight;   //!
   TBranch        *b_gsf_eTop;   //!
   TBranch        *b_gsf_eBottom;   //!
   TBranch        *b_gsf_e2nd;   //!
   TBranch        *b_gsf_p;   //!
   TBranch        *b_gsf_e;   //!
   TBranch        *b_gsf_pt;   //!
   TBranch        *b_gsf_class;   //!
   TBranch        *b_gsf_e2x5overe5x5;   //!
   TBranch        *b_gsf_e1x5overe5x5;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_px;   //!
   TBranch        *b_gsf_py;   //!
   TBranch        *b_gsf_pz;   //!
   TBranch        *b_gsf_deltaeta;   //!
   TBranch        *b_gsf_deltaphi;   //!
   TBranch        *b_gsf_hovere;   //!
   TBranch        *b_gsf_trackiso;   //!
   TBranch        *b_gsf_ecaliso;   //!
   TBranch        *b_gsf_hcaliso1;   //!
   TBranch        *b_gsf_hcaliso2;   //!
   TBranch        *b_gsf_charge;   //!
   TBranch        *b_gsf_sigmaetaeta;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_isecaldriven;   //!
   TBranch        *b_gsf_istrackerdriven;   //!
   TBranch        *b_gsfsc_e;   //!
   TBranch        *b_gsfsc_pt;   //!
   TBranch        *b_gsfsc_eta;   //!
   TBranch        *b_gsfsc_phi;   //!
   TBranch        *b_gsfsc_px;   //!
   TBranch        *b_gsfsc_py;   //!
   TBranch        *b_gsfsc_pz;   //!
   TBranch        *b_gsf_gsfet;   //!
   TBranch        *b_scindexforgsf;   //!
   TBranch        *b_gsfindexforgenelec;   //!
   TBranch        *b_gsfindexforgenposi;   //!
   TBranch        *b_scindexforgenelec;   //!
   TBranch        *b_scindexforgenposi;   //!
   TBranch        *b_gsfpass_ET;   //!
   TBranch        *b_gsfpass_PT;   //!
   TBranch        *b_gsfpass_DETETA;   //!
   TBranch        *b_gsfpass_CRACK;   //!
   TBranch        *b_gsfpass_DETAIN;   //!
   TBranch        *b_gsfpass_DPHIIN;   //!
   TBranch        *b_gsfpass_HADEM;   //!
   TBranch        *b_gsfpass_SIGMAIETAIETA;   //!
   TBranch        *b_gsfpass_E2X5OVER5X5;   //!
   TBranch        *b_gsfpass_ISOLEMHADDEPTH1;   //!
   TBranch        *b_gsfpass_ISOLHADDEPTH2;   //!
   TBranch        *b_gsfpass_ISOLPTTRKS;   //!
   TBranch        *b_gsfpass_ECALDRIVEN;   //!
   TBranch        *b_gsfpass_INVALID;   //!
   TBranch        *b_gsfpass_HEEP;   //!
   TBranch        *b_gsfpass_ID;   //!
   TBranch        *b_gsfpass_ISO;   //!
   TBranch        *b_scpixcharge;   //!
   TBranch        *b_ctfcharge;   //!
   TBranch        *b_gsfcharge;   //!
   TBranch        *b_gsfctfscpixconsistent;   //!
   TBranch        *b_gsfscpixconsistent;   //!
   TBranch        *b_gsfctfconsistent;   //!
   TBranch        *b_gsftracksize;   //!
   TBranch        *b_gsftracketa;   //!
   TBranch        *b_gsftrackphi;   //!
   TBranch        *b_gsftrackp;   //!
   TBranch        *b_gsftrackpt;   //!
   TBranch        *b_gsftrackpx;   //!
   TBranch        *b_gsftrackpy;   //!
   TBranch        *b_gsftrackpz;   //!

   FakeRate(TTree *tree=0);
   virtual ~FakeRate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
    virtual Bool_t   Notify();
  
 
};

#endif

#ifdef fakeRate_cxx
FakeRate::FakeRate(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/ZeeV6_OneEle.root");
      if (!f) {
         f = new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/ZeeV6_OneEle.root");
         f->cd("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/ZeeV6_OneEle.root:/gsfcheckerjob");
      }
      tree = (TTree*)gDirectory->Get("tree");

   }
   Init(tree);
}

FakeRate::~FakeRate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FakeRate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FakeRate::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FakeRate::Init(TTree *tree)
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

   fChain->SetBranchAddress("L1trigger_size", &L1trigger_size, &b_L1trigger_size);
   fChain->SetBranchAddress("L1trigger_bool", L1trigger_bool, &b_L1trigger_bool);
   fChain->SetBranchAddress("PhysDecl_bool", &PhysDecl_bool, &b_PhysDecl_bool);
   fChain->SetBranchAddress("HLTriggers", HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("nJetsAKT_pt15", &nJetsAKT_pt15, &b_nJetsAKT_pt15);
   fChain->SetBranchAddress("nJetsIC5_pt15", &nJetsIC5_pt15, &b_nJetsIC5_pt15);
   fChain->SetBranchAddress("calomet", &calomet, &b_calomet);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("jetIC5_size", &jetIC5_size, &b_jetIC5_size);
   fChain->SetBranchAddress("jetIC5_pt", jetIC5_pt, &b_jetIC5_pt);
   fChain->SetBranchAddress("jetIC5_eta", jetIC5_eta, &b_jetIC5_eta);
   fChain->SetBranchAddress("jetIC5_phi", jetIC5_phi, &b_jetIC5_phi);
   fChain->SetBranchAddress("jetIC5_em", jetIC5_em, &b_jetIC5_em);
   fChain->SetBranchAddress("muon_size", &muon_size, &b_muon_size);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_ptError", muon_ptError, &b_muon_ptError);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_etaError", muon_etaError, &b_muon_etaError);
   fChain->SetBranchAddress("muon_theta", muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaError", muon_thetaError, &b_muon_thetaError);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_phiError", muon_phiError, &b_muon_phiError);
   fChain->SetBranchAddress("muon_outerPt", muon_outerPt, &b_muon_outerPt);
   fChain->SetBranchAddress("muon_outerEta", muon_outerEta, &b_muon_outerEta);
   fChain->SetBranchAddress("muon_outerPhi", muon_outerPhi, &b_muon_outerPhi);
   fChain->SetBranchAddress("muon_outerTheta", muon_outerTheta, &b_muon_outerTheta);
   fChain->SetBranchAddress("muon_px", muon_px, &b_muon_px);
   fChain->SetBranchAddress("muon_py", muon_py, &b_muon_py);
   fChain->SetBranchAddress("muon_pz", muon_pz, &b_muon_pz);
   fChain->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_nhitstotal", muon_nhitstotal, &b_muon_nhitstotal);
   fChain->SetBranchAddress("muon_nhitstrack", muon_nhitstrack, &b_muon_nhitstrack);
   fChain->SetBranchAddress("muon_nhitsmuons", muon_nhitsmuons, &b_muon_nhitsmuons);
   fChain->SetBranchAddress("muon_nlosthits", muon_nlosthits, &b_muon_nlosthits);
   fChain->SetBranchAddress("muon_chi2", muon_chi2, &b_muon_chi2);
   fChain->SetBranchAddress("muon_ndof", muon_ndof, &b_muon_ndof);
   fChain->SetBranchAddress("muon_normChi2", muon_normChi2, &b_muon_normChi2);
   fChain->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("muon_d0Error", muon_d0Error, &b_muon_d0Error);
   fChain->SetBranchAddress("muon_dz_cmsCenter", muon_dz_cmsCenter, &b_muon_dz_cmsCenter);
   fChain->SetBranchAddress("muon_dz_beamSpot", muon_dz_beamSpot, &b_muon_dz_beamSpot);
   fChain->SetBranchAddress("muon_dz_firstPVtx", muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
   fChain->SetBranchAddress("muon_dzError", muon_dzError, &b_muon_dzError);
   fChain->SetBranchAddress("muon_dxy_cmsCenter", muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
   fChain->SetBranchAddress("muon_dxy_beamSpot", muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
   fChain->SetBranchAddress("muon_dxy_firstPVtx", muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
   fChain->SetBranchAddress("muon_dxyError", muon_dxyError, &b_muon_dxyError);
   fChain->SetBranchAddress("muon_innerPosx", muon_innerPosx, &b_muon_innerPosx);
   fChain->SetBranchAddress("muon_innerPosy", muon_innerPosy, &b_muon_innerPosy);
   fChain->SetBranchAddress("muon_innerPosz", muon_innerPosz, &b_muon_innerPosz);
   fChain->SetBranchAddress("muon_trackIso03", muon_trackIso03, &b_muon_trackIso03);
   fChain->SetBranchAddress("muon_trackIso05", muon_trackIso05, &b_muon_trackIso05);
   fChain->SetBranchAddress("muon_trackIso03_ptInVeto", muon_trackIso03_ptInVeto, &b_muon_trackIso03_ptInVeto);
   fChain->SetBranchAddress("muon_trackIso05_ptInVeto", muon_trackIso05_ptInVeto, &b_muon_trackIso05_ptInVeto);
   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("eventnumber", &eventnumber, &b_eventnumber);
   fChain->SetBranchAddress("eventcounter", &eventcounter, &b_eventcounter);
   fChain->SetBranchAddress("processid", &processid, &b_processid);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("alphaqcd", &alphaqcd, &b_alphaqcd);
   fChain->SetBranchAddress("alphaqed", &alphaqed, &b_alphaqed);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("sigmaZ", &sigmaZ, &b_sigmaZ);
   fChain->SetBranchAddress("sigmaZ0Error", &sigmaZ0Error, &b_sigmaZ0Error);
   fChain->SetBranchAddress("sq", &sq, &b_sq);
   fChain->SetBranchAddress("bsposx", &bsposx, &b_bsposx);
   fChain->SetBranchAddress("bsposy", &bsposy, &b_bsposy);
   fChain->SetBranchAddress("bsposz", &bsposz, &b_bsposz);
   fChain->SetBranchAddress("pvsize", &pvsize, &b_pvsize);
   fChain->SetBranchAddress("pvx", pvx, &b_pvx);
   fChain->SetBranchAddress("pvy", pvy, &b_pvy);
   fChain->SetBranchAddress("pvz", pvz, &b_pvz);
   fChain->SetBranchAddress("scsize", &scsize, &b_scsize);
   fChain->SetBranchAddress("scenergy", scenergy, &b_scenergy);
   fChain->SetBranchAddress("sceta", sceta, &b_sceta);
   fChain->SetBranchAddress("scetacorr", scetacorr, &b_scetacorr);
   fChain->SetBranchAddress("sctheta", sctheta, &b_sctheta);
   fChain->SetBranchAddress("scthetacorr", scthetacorr, &b_scthetacorr);
   fChain->SetBranchAddress("scet", scet, &b_scet);
   fChain->SetBranchAddress("scphi", scphi, &b_scphi);
   fChain->SetBranchAddress("scpx", scpx, &b_scpx);
   fChain->SetBranchAddress("scpy", scpy, &b_scpy);
   fChain->SetBranchAddress("scpz", scpz, &b_scpz);
   fChain->SetBranchAddress("scx", scx, &b_scx);
   fChain->SetBranchAddress("scy", scy, &b_scy);
   fChain->SetBranchAddress("scz", scz, &b_scz);
   fChain->SetBranchAddress("scgsfmatched", scgsfmatched, &b_scgsfmatched);
   fChain->SetBranchAddress("genelec_e_branch", &genelec_e_branch, &b_genelec_e_branch);
   fChain->SetBranchAddress("genelec_eta_branch", &genelec_eta_branch, &b_genelec_eta_branch);
   fChain->SetBranchAddress("genelec_phi_branch", &genelec_phi_branch, &b_genelec_phi_branch);
   fChain->SetBranchAddress("genelec_et_branch", &genelec_et_branch, &b_genelec_et_branch);
   fChain->SetBranchAddress("genposi_e_branch", &genposi_e_branch, &b_genposi_e_branch);
   fChain->SetBranchAddress("genposi_eta_branch", &genposi_eta_branch, &b_genposi_eta_branch);
   fChain->SetBranchAddress("genposi_phi_branch", &genposi_phi_branch, &b_genposi_phi_branch);
   fChain->SetBranchAddress("genposi_et_branch", &genposi_et_branch, &b_genposi_et_branch);
   fChain->SetBranchAddress("genelec_hassc_branch", &genelec_hassc_branch, &b_genelec_hassc_branch);
   fChain->SetBranchAddress("genposi_hassc_branch", &genposi_hassc_branch, &b_genposi_hassc_branch);
   fChain->SetBranchAddress("unstablegenelec_e_branch", &unstablegenelec_e_branch, &b_unstablegenelec_e_branch);
   fChain->SetBranchAddress("unstablegenelec_eta_branch", &unstablegenelec_eta_branch, &b_unstablegenelec_eta_branch);
   fChain->SetBranchAddress("unstablegenelec_phi_branch", &unstablegenelec_phi_branch, &b_unstablegenelec_phi_branch);
   fChain->SetBranchAddress("unstablegenelec_et_branch", &unstablegenelec_et_branch, &b_unstablegenelec_et_branch);
   fChain->SetBranchAddress("unstablegenposi_e_branch", &unstablegenposi_e_branch, &b_unstablegenposi_e_branch);
   fChain->SetBranchAddress("unstablegenposi_eta_branch", &unstablegenposi_eta_branch, &b_unstablegenposi_eta_branch);
   fChain->SetBranchAddress("unstablegenposi_phi_branch", &unstablegenposi_phi_branch, &b_unstablegenposi_phi_branch);
   fChain->SetBranchAddress("unstablegenposi_et_branch", &unstablegenposi_et_branch, &b_unstablegenposi_et_branch);
   fChain->SetBranchAddress("genboson_m_branch", &genboson_m_branch, &b_genboson_m_branch);
   fChain->SetBranchAddress("genboson_eta_branch", &genboson_eta_branch, &b_genboson_eta_branch);
   fChain->SetBranchAddress("genboson_phi_branch", &genboson_phi_branch, &b_genboson_phi_branch);
   fChain->SetBranchAddress("genboson_e_branch", &genboson_e_branch, &b_genboson_e_branch);
   fChain->SetBranchAddress("genboson_et_branch", &genboson_et_branch, &b_genboson_et_branch);
   fChain->SetBranchAddress("genboson_ez_branch", &genboson_ez_branch, &b_genboson_ez_branch);
   fChain->SetBranchAddress("genboson_p_branch", &genboson_p_branch, &b_genboson_p_branch);
   fChain->SetBranchAddress("genboson_pt_branch", &genboson_pt_branch, &b_genboson_pt_branch);
   fChain->SetBranchAddress("genboson_pz_branch", &genboson_pz_branch, &b_genboson_pz_branch);
   fChain->SetBranchAddress("x1quark", &x1quark, &b_x1quark);
   fChain->SetBranchAddress("x2quark", &x2quark, &b_x2quark);
   fChain->SetBranchAddress("fsrposiphotonsize", &fsrposiphotonsize, &b_fsrposiphotonsize);
   fChain->SetBranchAddress("fsrelecphotonsize", &fsrelecphotonsize, &b_fsrelecphotonsize);
   fChain->SetBranchAddress("energyfsrelec", energyfsrelec, &b_energyfsrelec);
   fChain->SetBranchAddress("etfsrelec", etfsrelec, &b_etfsrelec);
   fChain->SetBranchAddress("etafsrelec", etafsrelec, &b_etafsrelec);
   fChain->SetBranchAddress("phifsrelec", phifsrelec, &b_phifsrelec);
   fChain->SetBranchAddress("energyfsrposi", energyfsrposi, &b_energyfsrposi);
   fChain->SetBranchAddress("etfsrposi", etfsrposi, &b_etfsrposi);
   fChain->SetBranchAddress("etafsrposi", etafsrposi, &b_etafsrposi);
   fChain->SetBranchAddress("phifsrposi", phifsrposi, &b_phifsrposi);
   fChain->SetBranchAddress("scelecenergy", &scelecenergy, &b_scelecenergy);
   fChain->SetBranchAddress("sceleceta", &sceleceta, &b_sceleceta);
   fChain->SetBranchAddress("scelecphi", &scelecphi, &b_scelecphi);
   fChain->SetBranchAddress("scelecgsfmatched", &scelecgsfmatched, &b_scelecgsfmatched);
   fChain->SetBranchAddress("scposienergy", &scposienergy, &b_scposienergy);
   fChain->SetBranchAddress("scposieta", &scposieta, &b_scposieta);
   fChain->SetBranchAddress("scposiphi", &scposiphi, &b_scposiphi);
   fChain->SetBranchAddress("scposigsfmatched", &scposigsfmatched, &b_scposigsfmatched);
   fChain->SetBranchAddress("genelechassc", &genelechassc, &b_genelechassc);
   fChain->SetBranchAddress("genposihassc", &genposihassc, &b_genposihassc);
   fChain->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
   fChain->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_isEB", gsf_isEB, &b_gsf_isEB);
   fChain->SetBranchAddress("gsf_isEE", gsf_isEE, &b_gsf_isEE);
   fChain->SetBranchAddress("gsf_deltaEtaATcalo", gsf_deltaEtaATcalo, &b_gsf_deltaEtaATcalo);
   fChain->SetBranchAddress("gsf_deltaPhiATcalo", gsf_deltaPhiATcalo, &b_gsf_deltaPhiATcalo);
   fChain->SetBranchAddress("gsf_ecalEnergy", gsf_ecalEnergy, &b_gsf_ecalEnergy);
   fChain->SetBranchAddress("gsf_eOVERp", gsf_eOVERp, &b_gsf_eOVERp);
   fChain->SetBranchAddress("gsf_dxy", gsf_dxy, &b_gsf_dxy);
   fChain->SetBranchAddress("gsf_vz", gsf_vz, &b_gsf_vz);
   fChain->SetBranchAddress("gsf_nHits", gsf_nHits, &b_gsf_nHits);
   fChain->SetBranchAddress("gsf_fBrem", gsf_fBrem, &b_gsf_fBrem);
   fChain->SetBranchAddress("gsf_e1x5", gsf_e1x5, &b_gsf_e1x5);
   fChain->SetBranchAddress("gsf_e2x5", gsf_e2x5, &b_gsf_e2x5);
   fChain->SetBranchAddress("gsf_e5x5", gsf_e5x5, &b_gsf_e5x5);
   fChain->SetBranchAddress("gsf_eMax", gsf_eMax, &b_gsf_eMax);
   fChain->SetBranchAddress("gsf_SwissCross", gsf_SwissCross, &b_gsf_SwissCross);
   fChain->SetBranchAddress("gsf_e1x3", gsf_e1x3, &b_gsf_e1x3);
   fChain->SetBranchAddress("gsf_e3x1", gsf_e3x1, &b_gsf_e3x1);
   fChain->SetBranchAddress("gsf_e2x2", gsf_e2x2, &b_gsf_e2x2);
   fChain->SetBranchAddress("gsf_e3x2", gsf_e3x2, &b_gsf_e3x2);
   fChain->SetBranchAddress("gsf_e3x3", gsf_e3x3, &b_gsf_e3x3);
   fChain->SetBranchAddress("gsf_e4x4", gsf_e4x4, &b_gsf_e4x4);
   fChain->SetBranchAddress("gsf_e2x5Right", gsf_e2x5Right, &b_gsf_e2x5Right);
   fChain->SetBranchAddress("gsf_e2x5Left", gsf_e2x5Left, &b_gsf_e2x5Left);
   fChain->SetBranchAddress("gsf_e2x5Top", gsf_e2x5Top, &b_gsf_e2x5Top);
   fChain->SetBranchAddress("gsf_e2x5Bottom", gsf_e2x5Bottom, &b_gsf_e2x5Bottom);
   fChain->SetBranchAddress("gsf_e2x5Max", gsf_e2x5Max, &b_gsf_e2x5Max);
   fChain->SetBranchAddress("gsf_eLeft", gsf_eLeft, &b_gsf_eLeft);
   fChain->SetBranchAddress("gsf_eRight", gsf_eRight, &b_gsf_eRight);
   fChain->SetBranchAddress("gsf_eTop", gsf_eTop, &b_gsf_eTop);
   fChain->SetBranchAddress("gsf_eBottom", gsf_eBottom, &b_gsf_eBottom);
   fChain->SetBranchAddress("gsf_e2nd", gsf_e2nd, &b_gsf_e2nd);
   fChain->SetBranchAddress("gsf_p", gsf_p, &b_gsf_p);
   fChain->SetBranchAddress("gsf_e", gsf_e, &b_gsf_e);
   fChain->SetBranchAddress("gsf_pt", gsf_pt, &b_gsf_pt);
   fChain->SetBranchAddress("gsf_class", gsf_class, &b_gsf_class);
   fChain->SetBranchAddress("gsf_e2x5overe5x5", gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
   fChain->SetBranchAddress("gsf_e1x5overe5x5", gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
   fChain->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_px", gsf_px, &b_gsf_px);
   fChain->SetBranchAddress("gsf_py", gsf_py, &b_gsf_py);
   fChain->SetBranchAddress("gsf_pz", gsf_pz, &b_gsf_pz);
   fChain->SetBranchAddress("gsf_deltaeta", gsf_deltaeta, &b_gsf_deltaeta);
   fChain->SetBranchAddress("gsf_deltaphi", gsf_deltaphi, &b_gsf_deltaphi);
   fChain->SetBranchAddress("gsf_hovere", gsf_hovere, &b_gsf_hovere);
   fChain->SetBranchAddress("gsf_trackiso", gsf_trackiso, &b_gsf_trackiso);
   fChain->SetBranchAddress("gsf_ecaliso", gsf_ecaliso, &b_gsf_ecaliso);
   fChain->SetBranchAddress("gsf_hcaliso1", gsf_hcaliso1, &b_gsf_hcaliso1);
   fChain->SetBranchAddress("gsf_hcaliso2", gsf_hcaliso2, &b_gsf_hcaliso2);
   fChain->SetBranchAddress("gsf_charge", gsf_charge, &b_gsf_charge);
   fChain->SetBranchAddress("gsf_sigmaetaeta", gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_isecaldriven", gsf_isecaldriven, &b_gsf_isecaldriven);
   fChain->SetBranchAddress("gsf_istrackerdriven", gsf_istrackerdriven, &b_gsf_istrackerdriven);
   fChain->SetBranchAddress("gsfsc_e", gsfsc_e, &b_gsfsc_e);
   fChain->SetBranchAddress("gsfsc_pt", gsfsc_pt, &b_gsfsc_pt);
   fChain->SetBranchAddress("gsfsc_eta", gsfsc_eta, &b_gsfsc_eta);
   fChain->SetBranchAddress("gsfsc_phi", gsfsc_phi, &b_gsfsc_phi);
   fChain->SetBranchAddress("gsfsc_px", gsfsc_px, &b_gsfsc_px);
   fChain->SetBranchAddress("gsfsc_py", gsfsc_py, &b_gsfsc_py);
   fChain->SetBranchAddress("gsfsc_pz", gsfsc_pz, &b_gsfsc_pz);
   fChain->SetBranchAddress("gsf_gsfet", gsf_gsfet, &b_gsf_gsfet);
   fChain->SetBranchAddress("scindexforgsf", scindexforgsf, &b_scindexforgsf);
   fChain->SetBranchAddress("gsfindexforgenelec", &gsfindexforgenelec, &b_gsfindexforgenelec);
   fChain->SetBranchAddress("gsfindexforgenposi", &gsfindexforgenposi, &b_gsfindexforgenposi);
   fChain->SetBranchAddress("scindexforgenelec", &scindexforgenelec, &b_scindexforgenelec);
   fChain->SetBranchAddress("scindexforgenposi", &scindexforgenposi, &b_scindexforgenposi);
   fChain->SetBranchAddress("gsfpass_ET", gsfpass_ET, &b_gsfpass_ET);
   fChain->SetBranchAddress("gsfpass_PT", gsfpass_PT, &b_gsfpass_PT);
   fChain->SetBranchAddress("gsfpass_DETETA", gsfpass_DETETA, &b_gsfpass_DETETA);
   fChain->SetBranchAddress("gsfpass_CRACK", gsfpass_CRACK, &b_gsfpass_CRACK);
   fChain->SetBranchAddress("gsfpass_DETAIN", gsfpass_DETAIN, &b_gsfpass_DETAIN);
   fChain->SetBranchAddress("gsfpass_DPHIIN", gsfpass_DPHIIN, &b_gsfpass_DPHIIN);
   fChain->SetBranchAddress("gsfpass_HADEM", gsfpass_HADEM, &b_gsfpass_HADEM);
   fChain->SetBranchAddress("gsfpass_SIGMAIETAIETA", gsfpass_SIGMAIETAIETA, &b_gsfpass_SIGMAIETAIETA);
   fChain->SetBranchAddress("gsfpass_E2X5OVER5X5", gsfpass_E2X5OVER5X5, &b_gsfpass_E2X5OVER5X5);
   fChain->SetBranchAddress("gsfpass_ISOLEMHADDEPTH1", gsfpass_ISOLEMHADDEPTH1, &b_gsfpass_ISOLEMHADDEPTH1);
   fChain->SetBranchAddress("gsfpass_ISOLHADDEPTH2", gsfpass_ISOLHADDEPTH2, &b_gsfpass_ISOLHADDEPTH2);
   fChain->SetBranchAddress("gsfpass_ISOLPTTRKS", gsfpass_ISOLPTTRKS, &b_gsfpass_ISOLPTTRKS);
   fChain->SetBranchAddress("gsfpass_ECALDRIVEN", gsfpass_ECALDRIVEN, &b_gsfpass_ECALDRIVEN);
   fChain->SetBranchAddress("gsfpass_INVALID", gsfpass_INVALID, &b_gsfpass_INVALID);
   fChain->SetBranchAddress("gsfpass_HEEP", gsfpass_HEEP, &b_gsfpass_HEEP);
   fChain->SetBranchAddress("gsfpass_ID", gsfpass_ID, &b_gsfpass_ID);
   fChain->SetBranchAddress("gsfpass_ISO", gsfpass_ISO, &b_gsfpass_ISO);
   fChain->SetBranchAddress("scpixcharge", scpixcharge, &b_scpixcharge);
   fChain->SetBranchAddress("ctfcharge", ctfcharge, &b_ctfcharge);
   fChain->SetBranchAddress("gsfcharge", gsfcharge, &b_gsfcharge);
   fChain->SetBranchAddress("gsfctfscpixconsistent", gsfctfscpixconsistent, &b_gsfctfscpixconsistent);
   fChain->SetBranchAddress("gsfscpixconsistent", gsfscpixconsistent, &b_gsfscpixconsistent);
   fChain->SetBranchAddress("gsfctfconsistent", gsfctfconsistent, &b_gsfctfconsistent);
   fChain->SetBranchAddress("gsftracksize", &gsftracksize, &b_gsftracksize);
   fChain->SetBranchAddress("gsftracketa", gsftracketa, &b_gsftracketa);
   fChain->SetBranchAddress("gsftrackphi", gsftrackphi, &b_gsftrackphi);
   fChain->SetBranchAddress("gsftrackp", gsftrackp, &b_gsftrackp);
   fChain->SetBranchAddress("gsftrackpt", gsftrackpt, &b_gsftrackpt);
   fChain->SetBranchAddress("gsftrackpx", gsftrackpx, &b_gsftrackpx);
   fChain->SetBranchAddress("gsftrackpy", gsftrackpy, &b_gsftrackpy);
   fChain->SetBranchAddress("gsftrackpz", gsftrackpz, &b_gsftrackpz);
   Notify();
}

Bool_t FakeRate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

/* void FakeRate::Show(Long64_t entry) */
/* { */
/* // Print contents of entry. */
/* // If entry is not specified, print current entry */
/*    if (!fChain) return; */
/*    fChain->Show(entry); */
/* } */
Int_t FakeRate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


#endif // #ifdef fakeRate_cxx
