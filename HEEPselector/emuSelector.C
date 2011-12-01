#define emuSelector_cxx
#include "emuSelector.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>

void EmuSelector::Loop()
{
  // parameters //////////////////////////////////////////////////////////////
  const float eleCut = 35;
  const float muCut = 35;

  string fileNamePrefix = "emuSkim";   // lumi and .root well be added

  // input file
  //TFile * input = new TFile("/user/treis/mcsamples/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_12.root","read");
  //TFile * input = new TFile("/user/treis/mcsamples/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_12.root","read");
  //TFile * input = new TFile("/user/treis/mcsamples/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_12.root","read");
  TFile * input = new TFile("/user/treis/mcsamples/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_12.root","read");
  ////////////////////////////////////////////////////////////////////////////

  stringstream ssOutfile;
  ssOutfile << fileNamePrefix << ".root";
  TFile *outFile = new TFile(ssOutfile.str().c_str(), "recreate");
  outFile->mkdir("gsfcheckerjob");
  outFile->cd("gsfcheckerjob");

  TTree *thetree = (TTree*)input->Get("gsfcheckerjob/tree");
  Init(thetree);
  Long64_t nentries = (*thetree).GetEntries();
  cout << nentries << " entries" << endl;

  TTree *theNewTree = thetree->CloneTree(0);

  unsigned int counter = 0;

  /////////////////////////////////////////////////////////////////////////
  // loop over events
  /////////////////////////////////////////////////////////////////////////
  //for (Long64_t jentry=0; jentry < 10000; ++jentry) {
  for (Long64_t jentry=0; jentry < nentries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    bool eleGood = false;
    bool muGood = false;

    thetree->GetEntry(jentry);
    if (jentry % 50000 == 0 ) cout << "entry " << jentry << endl;

    // at least one gsf electron and one muon above the threshold
    if (gsf_size < 1 || muon_size < 1) continue;

    // loop over electrons
    for (int n = 0; n < gsf_size && !eleGood; ++n) {
      if (gsf_gsfet[n] < eleCut) continue;
      eleGood = true;
    }
    if (!eleGood) continue;

    // loop over muons
    for (int m = 0; m < muon_size && !muGood; ++m) {
      if (muon_pt[m] < muCut) continue;
      muGood = true;
    }
    if (!muGood) continue;

    theNewTree->Fill();
    counter++;

  } // end of loop over events

  theNewTree->AutoSave();
  outFile->Close();

  cout << "New tree filled with " << counter << " e-mu entries." << endl;

} //end of method

