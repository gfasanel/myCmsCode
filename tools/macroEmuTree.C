#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <string>
#include <TLorentzVector.h>
#include <TDCacheFile.h>
using namespace std;


void macroEmuTree() {
  bool noCuts = 0;
  string inputline;
  string outputline;
  string blankline;
  //ifstream myfile ("listofdata.txt");
  ifstream myfile ("listofsamples.txt");
  //ifstream myfile ("listofsigsamples.txt");
  if (myfile.is_open())
    {
      while ( myfile.good() )
	{
	  getline (myfile,inputline);
	  getline (myfile,outputline);
	  getline (myfile,blankline);

	  //gSystem->Load("$ROOTSYS/test/libEvent");
	  cout << "File : " << inputline << endl;
	  //Get old file, old tree and set top branch address
	  TString inputfilepath =inputline; 
	  TString outputfilepath = outputline; 
	  TDCacheFile *oldfile = new TDCacheFile(inputfilepath);
	  //TFile *oldfile = new TFile(inputfilepath);
	  TTree *oldtree = (TTree*)oldfile->Get("gsfcheckerjob/tree");
	  Long64_t nentries = oldtree->GetEntries();

          // Declaration of leaf types
          UInt_t          runnumber;
          UInt_t          eventnumber;
          UInt_t          luminosityBlock;
          Int_t           HLT_Mu22_Photon22_CaloIdL;
          Int_t           HLT_Mu40_eta2p1;
          Int_t           HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
          Int_t           HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
          Int_t           prescale_HLT_Mu22_Photon22_CaloIdL;
          Int_t           prescale_HLT_Mu40_eta2p1;
          Int_t           prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
          Int_t           prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
          Bool_t           hltL1sMu16Eta2p1;
          Bool_t           hltL1sL1Mu3p5EG12;
          Bool_t           hltL1sL1SingleEG12;
          Bool_t           hltMu22Photon22CaloIdLHEFilter;
          Bool_t           hltL1Mu3p5EG12L3Filtered22;
          Bool_t           hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q;
          Float_t         rho;
          Float_t         pfmet;
          Int_t           pvsize;
          Int_t           JetColl_size;
          Float_t         Jet_pt[100];   //[JetColl_size]
          Float_t         Jet_eta[100];   //[JetColl_size]
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
          Float_t         genelemom_mass[100];
          Int_t           genelemom_pdgid[100];
          Float_t         genele_pt[100];
          Float_t         genele_eta[100];
          Float_t         genele_phi[100];
          Float_t         genmu_pt[100];
          Float_t         genmu_eta[100];
          Float_t         genmu_phi[100];
          Float_t         hardGenEle_pt[100];
          Float_t         hardGenEle_eta[100];
          Float_t         hardGenEle_phi[100];
          Float_t         hardGenMu_pt[100];
          Float_t         hardGenMu_eta[100];
          Float_t         hardGenMu_phi[100];
          Float_t         genPart_pt[100];
          Float_t         genPart_mass[100];
          Int_t           genPart_pdgid[100];
          Int_t           genPart_status[100];
          Float_t         genPair_mass;
          Float_t         emu_mass;
          Float_t         res_mass;
          Int_t           trueNVtx;
       
          // List of branches
          TBranch        *b_runnumber;   //!
          TBranch        *b_eventnumber;   //!
          TBranch        *b_luminosityBlock;   //!
          TBranch        *b_HLT_Mu22_Photon22_CaloIdL;   //!
          TBranch        *b_HLT_Mu40_eta2p1;   //!
          TBranch        *b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
          TBranch        *b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
          TBranch        *b_prescale_HLT_Mu22_Photon22_CaloIdL;   //!
          TBranch        *b_prescale_HLT_Mu40_eta2p1;   //!
          TBranch        *b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
          TBranch        *b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
          TBranch        *b_hltL1sMu16Eta2p1;
          TBranch        *b_hltL1sL1Mu3p5EG12;
          TBranch        *b_hltL1sL1SingleEG12;
          TBranch        *b_hltMu22Photon22CaloIdLHEFilter;
          TBranch        *b_hltL1Mu3p5EG12L3Filtered22;
          TBranch        *b_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q;
          TBranch        *b_rho;   //!
          TBranch        *b_pfmet;   //!
          TBranch        *b_pvsize;   //!
          TBranch        *b_JetColl_size;   //!
          TBranch        *b_Jet_pt;   //!
          TBranch        *b_Jet_eta;   //!
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
          TBranch        *b_genelemom_mass;   //!
          TBranch        *b_genelemom_pdgid;   //!
          TBranch        *b_genele_pt;   //!
          TBranch        *b_genele_eta;   //!
          TBranch        *b_genele_phi;   //!
          TBranch        *b_genmu_pt;   //!
          TBranch        *b_genmu_eta;   //!
          TBranch        *b_genmu_phi;   //!
          TBranch        *b_hardGenEle_pt;   //!
          TBranch        *b_hardGenEle_eta;   //!
          TBranch        *b_hardGenEle_phi;   //!
          TBranch        *b_hardGenMu_pt;   //!
          TBranch        *b_hardGenMu_eta;   //!
          TBranch        *b_hardGenMu_phi;   //!
          TBranch        *b_genPart_pt;   //!
          TBranch        *b_genPart_mass;   //!
          TBranch        *b_genPart_pdgid;   //!
          TBranch        *b_genPart_status;   //!
          TBranch        *b_genPair_mass;   //!
          TBranch        *b_emu_mass;   //!
          TBranch        *b_res_mass;   //!
          TBranch        *b_trueNVtx;   //!
       
          oldtree->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
          oldtree->SetBranchAddress("eventnumber", &eventnumber, &b_eventnumber);
          oldtree->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
          oldtree->SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &HLT_Mu22_Photon22_CaloIdL, &b_HLT_Mu22_Photon22_CaloIdL);
          oldtree->SetBranchAddress("HLT_Mu40_eta2p1", &HLT_Mu40_eta2p1, &b_HLT_Mu40_eta2p1);
          oldtree->SetBranchAddress("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
          oldtree->SetBranchAddress("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
          oldtree->SetBranchAddress("prescale_HLT_Mu22_Photon22_CaloIdL", &prescale_HLT_Mu22_Photon22_CaloIdL, &b_prescale_HLT_Mu22_Photon22_CaloIdL);
          oldtree->SetBranchAddress("prescale_HLT_Mu40_eta2p1", &prescale_HLT_Mu40_eta2p1, &b_prescale_HLT_Mu40_eta2p1);
          oldtree->SetBranchAddress("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
          oldtree->SetBranchAddress("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
          oldtree->SetBranchAddress("hltL1sMu16Eta2p1", &hltL1sMu16Eta2p1, &b_hltL1sMu16Eta2p1);
          oldtree->SetBranchAddress("hltL1sL1Mu3p5EG12", &hltL1sL1Mu3p5EG12, &b_hltL1sL1Mu3p5EG12);
          oldtree->SetBranchAddress("hltL1sL1SingleEG12", &hltL1sL1SingleEG12, &b_hltL1sL1SingleEG12);
          oldtree->SetBranchAddress("hltMu22Photon22CaloIdLHEFilter", &hltMu22Photon22CaloIdLHEFilter, &b_hltMu22Photon22CaloIdLHEFilter);
          oldtree->SetBranchAddress("hltL1Mu3p5EG12L3Filtered22", &hltL1Mu3p5EG12L3Filtered22, &b_hltL1Mu3p5EG12L3Filtered22);
          oldtree->SetBranchAddress("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q", &hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q, &b_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q);
          oldtree->SetBranchAddress("rho", &rho, &b_rho);
          oldtree->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
          oldtree->SetBranchAddress("pvsize", &pvsize, &b_pvsize);
          oldtree->SetBranchAddress("JetColl_size", &JetColl_size, &b_JetColl_size);
          oldtree->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
          oldtree->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
          oldtree->SetBranchAddress("muon_size", &muon_size, &b_muon_size);
          oldtree->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
          oldtree->SetBranchAddress("muon_ptError", muon_ptError, &b_muon_ptError);
          oldtree->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
          oldtree->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
          oldtree->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
          oldtree->SetBranchAddress("muon_nhitspixel", muon_nhitspixel, &b_muon_nhitspixel);
          oldtree->SetBranchAddress("muon_nhitstrack", muon_nhitstrack, &b_muon_nhitstrack);
          oldtree->SetBranchAddress("muon_nhitsmuons", muon_nhitsmuons, &b_muon_nhitsmuons);
          oldtree->SetBranchAddress("muon_nlayerswithhits", muon_nlayerswithhits, &b_muon_nlayerswithhits);
          oldtree->SetBranchAddress("muon_nSegmentMatch", muon_nSegmentMatch, &b_muon_nSegmentMatch);
          oldtree->SetBranchAddress("muon_isTrackerMuon", muon_isTrackerMuon, &b_muon_isTrackerMuon);
          oldtree->SetBranchAddress("muon_normChi2", muon_normChi2, &b_muon_normChi2);
          oldtree->SetBranchAddress("muon_dz_beamSpot", muon_dz_beamSpot, &b_muon_dz_beamSpot);
          oldtree->SetBranchAddress("muon_dz_firstPVtx", muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
          oldtree->SetBranchAddress("muon_dxy_cmsCenter", muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
          oldtree->SetBranchAddress("muon_dxy_beamSpot", muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
          oldtree->SetBranchAddress("muon_dxy_firstPVtx", muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
          oldtree->SetBranchAddress("muon_trackIso03", muon_trackIso03, &b_muon_trackIso03);
          oldtree->SetBranchAddress("muon_emIso03", muon_emIso03, &b_muon_emIso03);
          oldtree->SetBranchAddress("muon_hadIso03", muon_hadIso03, &b_muon_hadIso03);
          oldtree->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
          oldtree->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
          oldtree->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
          oldtree->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
          oldtree->SetBranchAddress("gsf_charge", gsf_charge, &b_gsf_charge);
          oldtree->SetBranchAddress("gsf_sigmaetaeta", gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
          oldtree->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
          oldtree->SetBranchAddress("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
          oldtree->SetBranchAddress("gsf_dz_beamSpot", gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
          oldtree->SetBranchAddress("gsf_dz_firstPVtx", gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
          oldtree->SetBranchAddress("gsf_nLostInnerHits", gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
          oldtree->SetBranchAddress("gsf_deltaeta", gsf_deltaeta, &b_gsf_deltaeta);
          oldtree->SetBranchAddress("gsf_deltaphi", gsf_deltaphi, &b_gsf_deltaphi);
          oldtree->SetBranchAddress("gsf_hovere", gsf_hovere, &b_gsf_hovere);
          oldtree->SetBranchAddress("gsf_trackiso", gsf_trackiso, &b_gsf_trackiso);
          oldtree->SetBranchAddress("gsf_ecaliso", gsf_ecaliso, &b_gsf_ecaliso);
          oldtree->SetBranchAddress("gsf_hcaliso1", gsf_hcaliso1, &b_gsf_hcaliso1);
          oldtree->SetBranchAddress("gsf_hcaliso2", gsf_hcaliso2, &b_gsf_hcaliso2);
          oldtree->SetBranchAddress("gsf_isecaldriven", gsf_isecaldriven, &b_gsf_isecaldriven);
          oldtree->SetBranchAddress("gsfsc_e", gsfsc_e, &b_gsfsc_e);
          oldtree->SetBranchAddress("gsfsc_eta", gsfsc_eta, &b_gsfsc_eta);
          oldtree->SetBranchAddress("gsfsc_phi", gsfsc_phi, &b_gsfsc_phi);
          oldtree->SetBranchAddress("gsf_e2x5overe5x5", gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
          oldtree->SetBranchAddress("gsf_e1x5overe5x5", gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
          oldtree->SetBranchAddress("gsf_gsfet", gsf_gsfet, &b_gsf_gsfet);
          oldtree->SetBranchAddress("genelemom_mass", genelemom_mass, &b_genelemom_mass);
          oldtree->SetBranchAddress("genelemom_pdgid", genelemom_pdgid, &b_genelemom_pdgid);
          oldtree->SetBranchAddress("genele_pt", genele_pt, &b_genele_pt);
          oldtree->SetBranchAddress("genele_eta", genele_eta, &b_genele_eta);
          oldtree->SetBranchAddress("genele_phi", genele_phi, &b_genele_phi);
          oldtree->SetBranchAddress("genmu_pt", genmu_pt, &b_genmu_pt);
          oldtree->SetBranchAddress("genmu_eta", genmu_eta, &b_genmu_eta);
          oldtree->SetBranchAddress("genmu_phi", genmu_phi, &b_genmu_phi);
          oldtree->SetBranchAddress("hardGenEle_pt", hardGenEle_pt, &b_hardGenEle_pt);
          oldtree->SetBranchAddress("hardGenEle_eta", hardGenEle_eta, &b_hardGenEle_eta);
          oldtree->SetBranchAddress("hardGenEle_phi", hardGenEle_phi, &b_hardGenEle_phi);
          oldtree->SetBranchAddress("hardGenMu_pt", hardGenMu_pt, &b_hardGenMu_pt);
          oldtree->SetBranchAddress("hardGenMu_eta", hardGenMu_eta, &b_hardGenMu_eta);
          oldtree->SetBranchAddress("hardGenMu_phi", hardGenMu_phi, &b_hardGenMu_phi);
          oldtree->SetBranchAddress("genPart_pt", genPart_pt, &b_genPart_pt);
          oldtree->SetBranchAddress("genPart_mass", genPart_mass, &b_genPart_mass);
          oldtree->SetBranchAddress("genPart_pdgid", genPart_pdgid, &b_genPart_pdgid);
          oldtree->SetBranchAddress("genPart_status", genPart_status, &b_genPart_status);
          oldtree->SetBranchAddress("genPair_mass", &genPair_mass, &b_genPair_mass);
          oldtree->SetBranchAddress("emu_mass", &emu_mass, &b_emu_mass);
          oldtree->SetBranchAddress("res_mass", &res_mass, &b_res_mass);
          oldtree->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);
       
          // enable only used branches
          oldtree->SetBranchStatus("*", 0);
          oldtree->SetBranchStatus("runnumber", 1);
          oldtree->SetBranchStatus("eventnumber", 1);
          oldtree->SetBranchStatus("luminosityBlock", 1);
          oldtree->SetBranchStatus("HLT_Mu22_Photon22_CaloIdL", 1);
          oldtree->SetBranchStatus("HLT_Mu40_eta2p1", 1);
          oldtree->SetBranchStatus("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
          oldtree->SetBranchStatus("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
          oldtree->SetBranchStatus("prescale_HLT_Mu22_Photon22_CaloIdL", 1);
          oldtree->SetBranchStatus("prescale_HLT_Mu40_eta2p1", 1);
          oldtree->SetBranchStatus("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
          oldtree->SetBranchStatus("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
          oldtree->SetBranchStatus("hltL1sMu16Eta2p1", 1);
          oldtree->SetBranchStatus("hltL1sL1Mu3p5EG12", 1);
          oldtree->SetBranchStatus("hltL1sL1SingleEG12", 1);
          oldtree->SetBranchStatus("hltMu22Photon22CaloIdLHEFilter", 1);
          oldtree->SetBranchStatus("hltL1Mu3p5EG12L3Filtered22", 1);
          oldtree->SetBranchStatus("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q", 1);
          oldtree->SetBranchStatus("rho", 1);
          oldtree->SetBranchStatus("pfmet", 1);
          oldtree->SetBranchStatus("pvsize", 1);
          oldtree->SetBranchStatus("JetColl_size", 1);
          oldtree->SetBranchStatus("Jet_pt", 1);
          oldtree->SetBranchStatus("Jet_eta", 1);
          oldtree->SetBranchStatus("muon_size", 1);
          oldtree->SetBranchStatus("muon_pt", 1);
          oldtree->SetBranchStatus("muon_ptError", 1);
          oldtree->SetBranchStatus("muon_eta", 1);
          oldtree->SetBranchStatus("muon_phi", 1);
          oldtree->SetBranchStatus("muon_charge", 1);
          oldtree->SetBranchStatus("muon_nhitspixel", 1);
          oldtree->SetBranchStatus("muon_nhitstrack", 1);
          oldtree->SetBranchStatus("muon_nhitsmuons", 1);
          oldtree->SetBranchStatus("muon_nlayerswithhits", 1);
          oldtree->SetBranchStatus("muon_nSegmentMatch", 1);
          oldtree->SetBranchStatus("muon_isTrackerMuon", 1);
          oldtree->SetBranchStatus("muon_normChi2", 1);
          oldtree->SetBranchStatus("muon_dz_beamSpot", 1);
          oldtree->SetBranchStatus("muon_dz_firstPVtx", 1);
          oldtree->SetBranchStatus("muon_dxy_cmsCenter", 1);
          oldtree->SetBranchStatus("muon_dxy_beamSpot", 1);
          oldtree->SetBranchStatus("muon_dxy_firstPVtx", 1);
          oldtree->SetBranchStatus("muon_trackIso03", 1);
          oldtree->SetBranchStatus("muon_emIso03", 1);
          oldtree->SetBranchStatus("muon_hadIso03", 1);
          oldtree->SetBranchStatus("gsf_size", 1);
          oldtree->SetBranchStatus("gsf_eta", 1);
          oldtree->SetBranchStatus("gsf_phi", 1);
          oldtree->SetBranchStatus("gsf_theta", 1);
          oldtree->SetBranchStatus("gsf_charge", 1);
          oldtree->SetBranchStatus("gsf_sigmaetaeta", 1);
          oldtree->SetBranchStatus("gsf_sigmaIetaIeta", 1);
          oldtree->SetBranchStatus("gsf_dxy_firstPVtx", 1);
          oldtree->SetBranchStatus("gsf_dz_beamSpot", 1);
          oldtree->SetBranchStatus("gsf_dz_firstPVtx", 1);
          oldtree->SetBranchStatus("gsf_nLostInnerHits", 1);
          oldtree->SetBranchStatus("gsf_deltaeta", 1);
          oldtree->SetBranchStatus("gsf_deltaphi", 1);
          oldtree->SetBranchStatus("gsf_hovere", 1);
          oldtree->SetBranchStatus("gsf_trackiso", 1);
          oldtree->SetBranchStatus("gsf_ecaliso", 1);
          oldtree->SetBranchStatus("gsf_hcaliso1", 1);
          oldtree->SetBranchStatus("gsf_hcaliso2", 1);
          oldtree->SetBranchStatus("gsf_isecaldriven", 1);
          oldtree->SetBranchStatus("gsfsc_e", 1);
          oldtree->SetBranchStatus("gsfsc_eta", 1);
          oldtree->SetBranchStatus("gsfsc_phi", 1);
          oldtree->SetBranchStatus("gsf_e2x5overe5x5", 1);
          oldtree->SetBranchStatus("gsf_e1x5overe5x5", 1);
          oldtree->SetBranchStatus("gsf_gsfet", 1);
          oldtree->SetBranchStatus("genelemom_mass", 1);
          oldtree->SetBranchStatus("genelemom_pdgid", 1);
          oldtree->SetBranchStatus("genele_pt", 1);
          oldtree->SetBranchStatus("genele_eta", 1);
          oldtree->SetBranchStatus("genele_phi", 1);
          oldtree->SetBranchStatus("genmu_pt", 1);
          oldtree->SetBranchStatus("genmu_eta", 1);
          oldtree->SetBranchStatus("genmu_phi", 1);
          oldtree->SetBranchStatus("hardGenEle_pt", 1);
          oldtree->SetBranchStatus("hardGenEle_eta", 1);
          oldtree->SetBranchStatus("hardGenEle_phi", 1);
          oldtree->SetBranchStatus("hardGenMu_pt", 1);
          oldtree->SetBranchStatus("hardGenMu_eta", 1);
          oldtree->SetBranchStatus("hardGenMu_phi", 1);
          oldtree->SetBranchStatus("genPart_pt", 1);
          oldtree->SetBranchStatus("genPart_mass", 1);
          oldtree->SetBranchStatus("genPart_pdgid", 1);
          oldtree->SetBranchStatus("genPart_status", 1);
          oldtree->SetBranchStatus("genPair_mass", 1);
          oldtree->SetBranchStatus("emu_mass", 1);
          oldtree->SetBranchStatus("res_mass", 1);
          oldtree->SetBranchStatus("trueNVtx", 1);

  
	  //Create a new file + a clone of old tree header. Do not copy events
	  TFile *newfile = new TFile(outputfilepath,"recreate");
	  newfile->mkdir("gsfcheckerjob");
	  newfile->cd("gsfcheckerjob");
	  TTree *newtree = oldtree->CloneTree(0);
	  
	  
	  for (Long64_t i = 0; i < nentries; ++i) {
	    if (i %100000 == 0) cout << "entry nb: " <<i << endl;
	    oldtree->GetEntry(i);
            if (noCuts) {
	      newtree->Fill();
            } else {
	      if(gsf_size < 1) continue;
	      if(muon_size < 1) continue;
	      bool passEleSelection = false;
	      bool passMuSelection = false;
	      for (int it = 0; it < gsf_size; ++it) {
	        if (gsf_gsfet[it] < 35.) continue;
	          passEleSelection = true;
	      }
	      for (int muIt = 0; muIt < muon_size; ++muIt) {
	        if (muon_pt[muIt] < 35.) continue;
	          passMuSelection = true;
	      }
	      if (passEleSelection && passMuSelection) newtree->Fill();
            }
	  }


	  newtree->Print();
	  newfile->Write();
	  delete oldfile;
	  delete newfile;
	}
      myfile.close();
    } 
}
