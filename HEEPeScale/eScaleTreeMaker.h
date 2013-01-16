
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
   Int_t           HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   Int_t           trueNVtx;
   Float_t         rho;
   Int_t           pvsize;
   Int_t           gsf_size;
   Float_t         gsf_eta[150];   //[gsf_size]
   Float_t         gsf_phi[150];   //[gsf_size]
   Float_t         gsf_theta[150];   //[gsf_size]
   Float_t         gsf_sigmaIetaIeta[150];   //[gsf_size]
   Float_t         gsf_dxy_firstPVtx[150];   //[gsf_size]
   Int_t           gsf_nLostInnerHits[150];   //[gsf_size]
   Float_t         gsf_deltaeta[150];   //[gsf_size]
   Float_t         gsf_deltaphi[150];   //[gsf_size]
   Float_t         gsf_hovere[150];   //[gsf_size]
   Float_t         gsf_trackiso[150];   //[gsf_size]
   Float_t         gsf_ecaliso[150];   //[gsf_size]
   Float_t         gsf_hcaliso1[150];   //[gsf_size]
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
   TBranch        *b_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;   //!
   TBranch        *b_trueNVtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pvsize;   //!
   TBranch        *b_gsf_size;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
   TBranch        *b_gsf_nLostInnerHits;   //!
   TBranch        *b_gsf_deltaeta;   //!
   TBranch        *b_gsf_deltaphi;   //!
   TBranch        *b_gsf_hovere;   //!
   TBranch        *b_gsf_trackiso;   //!
   TBranch        *b_gsf_ecaliso;   //!
   TBranch        *b_gsf_hcaliso1;   //!
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
   fChain->SetBranchAddress("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   fChain->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pvsize", &pvsize, &b_pvsize);
   fChain->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
   fChain->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   fChain->SetBranchAddress("gsf_nLostInnerHits", gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   fChain->SetBranchAddress("gsf_deltaeta", gsf_deltaeta, &b_gsf_deltaeta);
   fChain->SetBranchAddress("gsf_deltaphi", gsf_deltaphi, &b_gsf_deltaphi);
   fChain->SetBranchAddress("gsf_hovere", gsf_hovere, &b_gsf_hovere);
   fChain->SetBranchAddress("gsf_trackiso", gsf_trackiso, &b_gsf_trackiso);
   fChain->SetBranchAddress("gsf_ecaliso", gsf_ecaliso, &b_gsf_ecaliso);
   fChain->SetBranchAddress("gsf_hcaliso1", gsf_hcaliso1, &b_gsf_hcaliso1);
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

   // activate only used branches
   fChain->SetBranchStatus("*", 0);
   fChain->SetBranchStatus("runnumber", 1);
   fChain->SetBranchStatus("eventnumber", 1);
   fChain->SetBranchStatus("luminosityBlock", 1);
   fChain->SetBranchStatus("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
   fChain->SetBranchStatus("prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
   fChain->SetBranchStatus("trueNVtx", 1);
   fChain->SetBranchStatus("rho", 1);
   fChain->SetBranchStatus("pvsize", 1);
   fChain->SetBranchStatus("gsf_size", 1);
   fChain->SetBranchStatus("gsf_eta", 1);
   fChain->SetBranchStatus("gsf_phi", 1);
   fChain->SetBranchStatus("gsf_theta", 1);
   fChain->SetBranchStatus("gsf_sigmaIetaIeta", 1);
   fChain->SetBranchStatus("gsf_dxy_firstPVtx", 1);
   fChain->SetBranchStatus("gsf_nLostInnerHits", 1);
   fChain->SetBranchStatus("gsf_deltaeta", 1);
   fChain->SetBranchStatus("gsf_deltaphi", 1);
   fChain->SetBranchStatus("gsf_hovere", 1);
   fChain->SetBranchStatus("gsf_trackiso", 1);
   fChain->SetBranchStatus("gsf_ecaliso", 1);
   fChain->SetBranchStatus("gsf_hcaliso1", 1);
   fChain->SetBranchStatus("gsf_isecaldriven", 1);
   fChain->SetBranchStatus("gsfsc_e", 1);
   fChain->SetBranchStatus("gsfsc_eta", 1);
   fChain->SetBranchStatus("gsf_e2x5overe5x5", 1);
   fChain->SetBranchStatus("gsf_e1x5overe5x5", 1);
   fChain->SetBranchStatus("gsf_gsfet", 1);
   fChain->SetBranchStatus("genele_e", 1);
   fChain->SetBranchStatus("genele_pt", 1);
   fChain->SetBranchStatus("genele_eta", 1);
   fChain->SetBranchStatus("genele_phi", 1);
   fChain->SetBranchStatus("genelemom_mass", 1);
   fChain->SetBranchStatus("unstableGenEle_e", 1);

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
