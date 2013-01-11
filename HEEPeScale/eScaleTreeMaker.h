
#ifndef invariantMass_h
#define invariantMass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"

class EScaleTreeMaker {
public :
   // constants
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TFile* _outF;
  
   double CalcInvariantMass (const int&, const int&, bool&);
   bool PassHEEP(const int &n);
   int Trigger(int &prescale);
   int GenEleDRMatch(const int&, double &);
   void CorrectEnergy();

   // Declaration of leaf types
   UInt_t          runnumber;
   UInt_t          eventnumber;
   UInt_t          luminosityBlock;
   UInt_t          eventcounter;
   Int_t           processid;
   Float_t         pthat;
   Float_t         qscale;
   Float_t         weight;
   Int_t           hltCount;
   Int_t           PhysDecl_bool;
   Int_t           nWasRun_;
   Int_t           nAccept_;
   Int_t           nErrors_;
   Int_t           HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           trueNVtx;
   Float_t         rho;
   Int_t           pvsize;
   Bool_t          pv_isValid[80];   //[pvsize]
   Int_t           gsf_size;
   Float_t         gsf_eta[150];   //[gsf_size]
   Float_t         gsf_phi[150];   //[gsf_size]
   Float_t         gsf_theta[150];   //[gsf_size]
   Float_t         gsf_sigmaIetaIeta[150];   //[gsf_size]
   Float_t         gsf_ecalEnergy[150];   //[gsf_size]
   Float_t         gsf_eOVERp[150];   //[gsf_size]
   Float_t         gsf_dxy_firstPVtx[150];   //[gsf_size]
   Int_t           gsf_nLostInnerHits[150];   //[gsf_size]
   Float_t         gsf_deltaeta[150];   //[gsf_size]
   Float_t         gsf_deltaphi[150];   //[gsf_size]
   Float_t         gsf_hovere[150];   //[gsf_size]
   Float_t         gsf_trackiso[150];   //[gsf_size]
   Float_t         gsf_ecaliso[150];   //[gsf_size]
   Float_t         gsf_hcaliso1[150];   //[gsf_size]
   Float_t         gsf_hcaliso2[150];   //[gsf_size]
   Bool_t          gsf_isecaldriven[150];   //[gsf_size]
   Float_t         gsfsc_e[150];   //[gsf_size]
   Float_t         gsfsc_eta[150];   //[gsf_size]
   Float_t         gsf_e2x5overe5x5[150];   //[gsf_size]
   Float_t         gsf_e1x5overe5x5[150];   //[gsf_size]
   Float_t         gsf_gsfet[150];   //[gsf_size]
   Float_t         genele_e[150];   //[genparticles_size]
   Float_t         genele_pt[150];   //[genparticles_size]
   Float_t         genele_eta[150];   //[genparticles_size]
   Float_t         genele_phi[150];   //[genparticles_size]
   Float_t         genelemom_mass[150];   //[genparticles_size]
   Float_t         unstableGenEle_e[150];   //[genparticles_size]

   // List of branches
   TBranch        *b_runnumber;   //!
   TBranch        *b_eventnumber;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_eventcounter;   //!
   TBranch        *b_processid;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_qscale;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_hltCount;   //!
   TBranch        *b_PhysDecl_bool;   //!
   TBranch        *b_nWasRun_;   //!
   TBranch        *b_nAccept_;   //!
   TBranch        *b_nErrors_;   //!
   TBranch        *b_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_trueNVtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pvsize;   //!
   TBranch        *b_pv_isValid;   //!
   TBranch        *b_gsf_size;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_ecalEnergy;   //!
   TBranch        *b_gsf_eOVERp;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
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
   TBranch        *b_gsf_e2x5overe5x5;   //!
   TBranch        *b_gsf_e1x5overe5x5;   //!
   TBranch        *b_gsf_gsfet;   //!
   TBranch        *b_genele_e;   //!
   TBranch        *b_genele_pt;   //!
   TBranch        *b_genele_eta;   //!
   TBranch        *b_genele_phi;   //!
   TBranch        *b_genelemom_mass;   //!
   TBranch        *b_unstableGenEle_e;   //!

   EScaleTreeMaker(TTree *tree=0);
   virtual ~EScaleTreeMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
  
 
};

#endif

#ifdef invariantMass_cxx
EScaleTreeMaker::EScaleTreeMaker(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //if (tree == 0) {
   //   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/ZeeV6_OneEle.root");
   //   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/user/treis/data2012/Photon_Run2012A-PromptReco-v1_AOD_Cert_190456-190688_8TeV_PromptReco_Collisions12_JSON_gct1_24.root");
   //   if (!f) {
   //      f = new TFile("/user/treis/data2012/Photon_Run2012A-PromptReco-v1_AOD_Cert_190456-190688_8TeV_PromptReco_Collisions12_JSON_gct1_24.root");
   //      f->cd("/user/treis/data2012/Photon_Run2012A-PromptReco-v1_AOD_Cert_190456-190688_8TeV_PromptReco_Collisions12_JSON_gct1_24.root:/gsfcheckerjob");
   //   }
   //   tree = (TTree*)gDirectory->Get("tree");

   //}
   //Init(tree);
}

EScaleTreeMaker::~EScaleTreeMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EScaleTreeMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EScaleTreeMaker::LoadTree(Long64_t entry)
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

void EScaleTreeMaker::Init(TTree *tree)
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
   fChain->SetBranchAddress("eventcounter", &eventcounter, &b_eventcounter);
   fChain->SetBranchAddress("processid", &processid, &b_processid);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("hltCount", &hltCount, &b_hltCount);
   fChain->SetBranchAddress("PhysDecl_bool", &PhysDecl_bool, &b_PhysDecl_bool);
   fChain->SetBranchAddress("nWasRun_", &nWasRun_, &b_nWasRun_);
   fChain->SetBranchAddress("nAccept_", &nAccept_, &b_nAccept_);
   fChain->SetBranchAddress("nErrors_", &nErrors_, &b_nErrors_);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pvsize", &pvsize, &b_pvsize);
   fChain->SetBranchAddress("pv_isValid", pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
   fChain->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_ecalEnergy", gsf_ecalEnergy, &b_gsf_ecalEnergy);
   fChain->SetBranchAddress("gsf_eOVERp", gsf_eOVERp, &b_gsf_eOVERp);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
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
   fChain->SetBranchAddress("gsf_e2x5overe5x5", gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
   fChain->SetBranchAddress("gsf_e1x5overe5x5", gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
   fChain->SetBranchAddress("gsf_gsfet", gsf_gsfet, &b_gsf_gsfet);
   fChain->SetBranchAddress("genele_e", genele_e, &b_genele_e);
   fChain->SetBranchAddress("genele_pt", genele_pt, &b_genele_pt);
   fChain->SetBranchAddress("genele_eta", genele_eta, &b_genele_eta);
   fChain->SetBranchAddress("genele_phi", genele_phi, &b_genele_phi);
   fChain->SetBranchAddress("genelemom_mass", genelemom_mass, &b_genelemom_mass);
   fChain->SetBranchAddress("unstableGenEle_e", unstableGenEle_e, &b_unstableGenEle_e);

   Notify();
}

Bool_t EScaleTreeMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

/* void EScaleTreeMaker::Show(Long64_t entry) */
/* { */
/* // Print contents of entry. */
/* // If entry is not specified, print current entry */
/*    if (!fChain) return; */
/*    fChain->Show(entry); */
/* } */
Int_t EScaleTreeMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


#endif // #ifdef invariantMass_cxx
