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


void macroEeSpecTree() {
  string inputline;
  string outputline;
  string blankline;
  ifstream myfile ("listofsamples_eeSpectrum.txt");
  if (myfile.is_open())
    {
      while ( myfile.good() )
	{
	  getline(myfile, inputline);
	  getline(myfile, outputline);
	  getline(myfile, blankline);

	  //gSystem->Load("$ROOTSYS/test/libEvent");
	  cout << "File: " << inputline << endl;
	  //Get old file, old tree and set top branch address
	  TString inputfilepath = inputline; 
	  TString outputfilepath = outputline; 
	  TDCacheFile *oldfile = new TDCacheFile(inputfilepath);
	  //TFile *oldfile = new TFile(inputfilepath);
	  TTree *oldtree = (TTree*)oldfile->Get("gsfcheckerjob/tree");
	  Long64_t nentries = oldtree->GetEntries();

          cout << "total entries: " << nentries << endl;

          // Declaration of leaf types
          UInt_t runnumber;
          UInt_t eventnumber;
          UInt_t luminosityBlock;
          Int_t HLT_DoubleEle33_CaloIdL;
          Int_t HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
          Int_t HLT_DoubleEle33_CaloIdT;
          Int_t prescale_HLT_DoubleEle33_CaloIdL;
          Int_t prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
          Int_t prescale_HLT_DoubleEle33_CaloIdT;
          Float_t rho;
          Float_t calomet;
          Int_t gsf_size;
          Float_t gsf_eta[150];   //[gsf_size]
          Float_t gsf_phi[150];   //[gsf_size]
          Float_t gsf_theta[150];   //[gsf_size]
          Float_t gsf_sigmaIetaIeta[150];   //[gsf_size]
          Float_t gsf_dxy_firstPVtx[150];   //[gsf_size]
          Int_t gsf_nLostInnerHits[150];   //[gsf_size]
          Float_t gsf_deltaeta[150];   //[gsf_size]
          Float_t gsf_deltaphi[150];   //[gsf_size]
          Float_t gsf_hovere[150];   //[gsf_size]
          Float_t gsf_trackiso[150];   //[gsf_size]
          Float_t gsf_ecaliso[150];   //[gsf_size]
          Float_t gsf_hcaliso1[150];   //[gsf_size]
          Float_t gsf_hcaliso2[150];   //[gsf_size]
          Bool_t gsf_isecaldriven[150];   //[gsf_size]
          Float_t gsfsc_e[150];   //[gsf_size]
          Float_t gsfsc_eta[150];   //[gsf_size]
          Float_t gsf_e2x5overe5x5[150];   //[gsf_size]
          Float_t gsf_e1x5overe5x5[150];   //[gsf_size]
          Float_t gsf_gsfet[150];   //[gsf_size]
          Float_t genelemom_mass[150];   //[genparticles_size]
          Int_t trueNVtx;


          // List of branches
          TBranch *b_runnumber;
          TBranch *b_eventnumber;
          TBranch *b_luminosityBlock;
          TBranch *b_HLT_DoubleEle33_CaloIdL;
          TBranch *b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
          TBranch *b_HLT_DoubleEle33_CaloIdT;
          TBranch *b_prescale_HLT_DoubleEle33_CaloIdL;
          TBranch *b_prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
          TBranch *b_prescale_HLT_DoubleEle33_CaloIdT;
          TBranch *b_rho;
          TBranch *b_calomet;
          TBranch *b_gsf_size;
          TBranch *b_gsf_eta;   
          TBranch *b_gsf_phi;   
          TBranch *b_gsf_theta;   
          TBranch *b_gsf_sigmaIetaIeta;   
          TBranch *b_gsf_dxy_firstPVtx;   
          TBranch *b_gsf_nLostInnerHits;   
          TBranch *b_gsf_deltaeta;   
          TBranch *b_gsf_deltaphi;   
          TBranch *b_gsf_hovere;   
          TBranch *b_gsf_trackiso;   
          TBranch *b_gsf_ecaliso;   
          TBranch *b_gsf_hcaliso1;   
          TBranch *b_gsf_hcaliso2;   
          TBranch *b_gsf_isecaldriven;   
          TBranch *b_gsfsc_e;   
          TBranch *b_gsfsc_eta;   
          TBranch *b_gsf_e2x5overe5x5;   
          TBranch *b_gsf_e1x5overe5x5;   
          TBranch *b_gsf_gsfet;   
          TBranch *b_genelemom_mass;   
          TBranch *b_trueNVtx;

          oldtree->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
          oldtree->SetBranchAddress("eventnumber", &eventnumber, &b_eventnumber);
          oldtree->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
          oldtree->SetBranchAddress("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL, &b_HLT_DoubleEle33_CaloIdL);
          oldtree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
          oldtree->SetBranchAddress("HLT_DoubleEle33_CaloIdT", &HLT_DoubleEle33_CaloIdT, &b_HLT_DoubleEle33_CaloIdT);
          oldtree->SetBranchAddress("prescale_HLT_DoubleEle33_CaloIdL", &prescale_HLT_DoubleEle33_CaloIdL, &b_prescale_HLT_DoubleEle33_CaloIdL);
          oldtree->SetBranchAddress("prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b_prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
          oldtree->SetBranchAddress("prescale_HLT_DoubleEle33_CaloIdT", &prescale_HLT_DoubleEle33_CaloIdT, &b_prescale_HLT_DoubleEle33_CaloIdT);
          oldtree->SetBranchAddress("rho", &rho, &b_rho);
          oldtree->SetBranchAddress("calomet", &calomet, &b_calomet);
          oldtree->SetBranchAddress("gsf_size", &gsf_size, &b_gsf_size);
          oldtree->SetBranchAddress("gsf_eta", gsf_eta, &b_gsf_eta);
          oldtree->SetBranchAddress("gsf_phi", gsf_phi, &b_gsf_phi);
          oldtree->SetBranchAddress("gsf_theta", gsf_theta, &b_gsf_theta);
          oldtree->SetBranchAddress("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
          oldtree->SetBranchAddress("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
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
          oldtree->SetBranchAddress("gsf_e2x5overe5x5", gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
          oldtree->SetBranchAddress("gsf_e1x5overe5x5", gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
          oldtree->SetBranchAddress("gsf_gsfet", gsf_gsfet, &b_gsf_gsfet);
          oldtree->SetBranchAddress("genelemom_mass", &genelemom_mass, &b_genelemom_mass);
          oldtree->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);

   // select only used branches
          oldtree->SetBranchStatus("*", 0);
          oldtree->SetBranchStatus("runnumber", 1);
          oldtree->SetBranchStatus("eventnumber", 1);
          oldtree->SetBranchStatus("luminosityBlock", 1);
          oldtree->SetBranchStatus("eventcounter", 1);
          oldtree->SetBranchStatus("HLT_DoubleEle33_CaloIdL", 1);
          oldtree->SetBranchStatus("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", 1);
          oldtree->SetBranchStatus("HLT_DoubleEle33_CaloIdT", 1);
          oldtree->SetBranchStatus("prescale_HLT_DoubleEle33_CaloIdL", 1);
          oldtree->SetBranchStatus("prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", 1);
          oldtree->SetBranchStatus("prescale_HLT_DoubleEle33_CaloIdT", 1);
          oldtree->SetBranchStatus("rho", 1);
          oldtree->SetBranchStatus("calomet", 1);
          oldtree->SetBranchStatus("gsf_size", 1);
          oldtree->SetBranchStatus("gsf_eta", 1);
          oldtree->SetBranchStatus("gsf_phi", 1);
          oldtree->SetBranchStatus("gsf_theta", 1);
          oldtree->SetBranchStatus("gsf_sigmaIetaIeta", 1);
          oldtree->SetBranchStatus("gsf_dxy_firstPVtx", 1);
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
          oldtree->SetBranchStatus("gsf_e2x5overe5x5", 1);
          oldtree->SetBranchStatus("gsf_e1x5overe5x5", 1);
          oldtree->SetBranchStatus("gsf_gsfet", 1);
          oldtree->SetBranchStatus("genelemom_mass", 1);
          oldtree->SetBranchStatus("trueNVtx", 1);

	  //Create a new file + a clone of old tree header. Do not copy events
	  TFile *newfile = new TFile(outputfilepath,"recreate");
	  newfile->mkdir("gsfcheckerjob");
	  newfile->cd("gsfcheckerjob");
	  TTree *newtree = oldtree->CloneTree(0);
	  
	  for (Long64_t i = 0; i < nentries; ++i) {
	    if (i %100000 == 0) cout << "entry nb: " <<i << endl;
	    oldtree->GetEntry(i);
	    if(gsf_size < 2) continue;
	    bool passEleSelection1 = false;
	    bool passEleSelection2 = false;
	    for (int it = 0; it < gsf_size; ++it) {
	      if (gsf_gsfet[it] < 35.) continue;
		if (!passEleSelection1) passEleSelection1 = true;
                else if (!passEleSelection2) passEleSelection2 = true;
	    }
	    if (passEleSelection1 && passEleSelection2) newtree->Fill();
     
	  }

	  newtree->Print();
	  newfile->Write();
	  delete oldfile;
	  delete newfile;
	}
      myfile.close();
    } 
}
