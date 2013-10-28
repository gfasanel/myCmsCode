#define invariantMass_cxx
#include "eScaleTreeMaker.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>

#include "TLorentzVector.h"
#include "TString.h"
#include "TStopwatch.h"

void EScaleTreeMaker::Loop()
{
  TStopwatch timer;
  timer.Start();
  // parameters //////////////////////////////////////////////////////////////
  const int ptcut = 0.;
  const float massCut = 0.;
  bool etShiftData = false;           // match data Z peak with MC peak

  double drMax = 1.;

  TString outFile("eScaleEvents");
  //TString outFile("eScaleEventsTest");

  bool usePUInfo = true;
  TString puFile = "file:////user/treis/data2013/pileup/pileupTrue_DoubleElectron_Run2012ABCDReReco22Jan2013.root";
  float bar_et = 25.;
  float end_et = 25.;
  ////////////////////////////////////////////////////////////////////////////

  TH1::SetDefaultSumw2(kTRUE);

  /////////////////////////////////////////////////////////////////////////
  // input files
  /////////////////////////////////////////////////////////////////////////
  vector<pair<TFile *, double> > input;
  vector<TString> inFileTag;
  vector<float> minPtDY;
  float lumi = 19712.;
  TFile *inData = TFile::Open("file:////user/treis/data2013/DoubleElectron_Run2012A+B+C+D-22Jan2013-v1_AOD_eScaleSkim_19712pb-1.root");
  input.push_back(make_pair(inData, 1 / lumi));
  inFileTag.push_back("Data");
  const unsigned int DATA = 0;

  TFile *inDY20 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1+2_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY20, 1915. / (3297045. + 42705454.)));
  inFileTag.push_back("DY20");
  minPtDY.push_back(20);
  const unsigned int DY20 = 1;

  TFile *inDY120 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY120, 11.89 / 99987. * 1915. / 1871.));
  inFileTag.push_back("DY120");
  minPtDY.push_back(120);
  const unsigned int DY120 = 2;

  TFile *inDY200 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY200, 1.483 / 99991. * 1915. / 1871.));
  inFileTag.push_back("DY200");
  minPtDY.push_back(200);
  const unsigned int DY200 = 3;

  TFile *inDY400 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-400_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY400, 0.1085 / 99991. * 1915. / 1871.));
  inFileTag.push_back("DY400");
  minPtDY.push_back(400);
  const unsigned int DY400 = 4;

  TFile *inDY500 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY500, 0.04409 / 99986. * 1915. / 1871.));
  inFileTag.push_back("DY500");
  minPtDY.push_back(500);
  const unsigned int DY500 = 5;

  TFile *inDY700 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-700_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY700, 0.01025 / 99990. * 1915. / 1871.));
  inFileTag.push_back("DY700");
  minPtDY.push_back(700);
  const unsigned int DY700 = 6;

  TFile *inDY800 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY800, 0.005491 / 99990. * 1915. / 1871.));
  inFileTag.push_back("DY800");
  minPtDY.push_back(800);
  const unsigned int DY800 = 7;

  TFile *inDY1000 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY1000, 0.001796 / 99992. * 1915. / 1871.));
  inFileTag.push_back("DY1000");
  minPtDY.push_back(1000);
  const unsigned int DY1000 = 8;

  TFile *inDY1500 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY1500, 1.705E-4 / 99999. * 1915. / 1871.));
  inFileTag.push_back("DY1500");
  minPtDY.push_back(1500);
  const unsigned int DY1500 = 9;

  TFile *inDY2000 = TFile::Open("file:////user/treis/mcsamples/mc2013/DYToEE_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inDY2000, 2.208E-5 / 99993. * 1915. / 1871.));
  inFileTag.push_back("DY2000");
  minPtDY.push_back(2000);
  const unsigned int DY2000 = 10;

  TFile *inZpPsi750 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-750_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi750, 0.1328 / 25200.));
  inFileTag.push_back("Zp750");
  const unsigned int ZP750 = 11;

  TFile *inZpPsi1000 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-1000_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi1000, 0.03933 / 25200.));
  inFileTag.push_back("Zp1000");
  const unsigned int ZP1000 = 12;

  TFile *inZpPsi1250 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-1250_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi1250, 0.01196 / 25200.));
  inFileTag.push_back("Zp1250");
  const unsigned int ZP1250 = 13;

  TFile *inZpPsi1500 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-1500_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi1500, 0.00437 / 25200.));
  inFileTag.push_back("Zp1500");
  const unsigned int ZP1500 = 14;

  TFile *inZpPsi1750 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-1750_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi1750, 0.00168 / 25280.));
  inFileTag.push_back("Zp1750");
  const unsigned int ZP1750 = 15;

  TFile *inZpPsi2000 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-2000_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi2000, 7.029E-4 / 25280.));
  inFileTag.push_back("Zp2000");
  const unsigned int ZP2000 = 16;

  TFile *inZpPsi2250 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-2250_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_25280ev.root");
  input.push_back(make_pair(inZpPsi2250, 2.895E-4 / 25280.));
  inFileTag.push_back("Zp2250");
  const unsigned int ZP2250 = 17;

  TFile *inZpPsi3000 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimePSIToEE_M-3000_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpPsi3000, 2.666E-5 / 25280.));
  inFileTag.push_back("Zp3000");
  const unsigned int ZP3000 = 18;

  TFile *inZpSsm2250 = TFile::Open("file:////user/treis/mcsamples/mc2013/ZprimeSSMToEE_M-2250_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_eScaleSkim.root");
  input.push_back(make_pair(inZpSsm2250, 0.001302 / 25280.));
  inFileTag.push_back("ZpSsm2250");
  const unsigned int ZPSSM2250 = 19;
  ////////////////////////////////////////////////////////////////////////////
  
  minPtDY.push_back(100000); // set ridiculously high to catch the tail of the last DY sample

  stringstream ssGoodHeepFileName;
  ssGoodHeepFileName << outFile << lumi << "pb-1.root";
  TFile *goodEvFile = new TFile(ssGoodHeepFileName.str().c_str(), "recreate");
  goodEvFile->Cd("");
  int evtRegion = -1;
  int evtRegionDrMatched = -1;
  float hHMass = -10000.;
  float hHMassDrMatched = -10000.;
  vector<TTree *> eleTrees;

  // data pileup histogram
  TFile *dataPuInfile = TFile::Open(puFile);
  dataPuInfile->Cd("");
  TH1F *puData = (TH1F *)gDirectory->Get("pileup");
  puData->SetDirectory(0);
  puData->SetName("puData");
  dataPuInfile->Close();
  TH1F *puDataNorm = (TH1F *)puData->DrawNormalized()->Clone("puDataNorm"); 

  /////////////////////////////////////////////////////////////////////////
  // loop over files
  /////////////////////////////////////////////////////////////////////////
  for (unsigned int p = 0; p < input.size(); ++p) {
  //for (unsigned int p = 0; p < 1; ++p) {
    string infile(input[p].first->GetName());
    cout << "file " << infile << endl;
    input[p].first->Cd("");
    TTree *thetree = (TTree*)input[p].first->Get("gsfcheckerjob/tree");
    Init(thetree);
    Long64_t nentries = (*thetree).GetEntries();
    cout << nentries << " entries" << endl;

    // get the histogram with the true number of vertices
    TH1F *puMc = new TH1F("puMc" + inFileTag[p], "puMc" + inFileTag[p], 100, 0., 100.);
    thetree->Draw("trueNVtx>>puMc" + inFileTag[p]);
    TH1F *puMcNorm = new TH1F("dummy" + inFileTag[p], "dummy" + inFileTag[p], 100, 0., 100.);
    if (p > DATA) puMcNorm = (TH1F *)puMc->DrawNormalized();
    puMcNorm->SetName("puMc" + inFileTag[p] + "Norm");
    puMcNorm->SetTitle("puMc" + inFileTag[p] + "Norm");
    // calculate the pu weights
    TH1F *puWeights = new TH1F("puWeight" + inFileTag[p], "puWeight" + inFileTag[p], 100, 0., 100.);
    if (p > DATA) puWeights->Divide(puDataNorm, puMcNorm);

    float mcWeight = input[p].second;
    float puWeight = 1.;
    float trueMass = -10.;
    float ele1Et = -10.;
    float ele1ScE = -10.;
    float ele1Eta = -10.;
    float ele2Et = -10.;
    float ele2ScE = -10.;
    float ele2Eta = -10.;
    float ele1EtDrMatched = -10.;
    float ele1ScEDrMatched = -10.;
    float ele1EtaDrMatched = -10.;
    float ele2EtDrMatched = -10.;
    float ele2ScEDrMatched = -10.;
    float ele2EtaDrMatched = -10.;
    float genEle1E = -10.;
    float genEle2E = -10.;
    float unstGenEle1E = -10.;
    float unstGenEle2E = -10.;
    // Make a new tree from this sample
    eleTrees.push_back(new TTree("ele" + inFileTag[p] + "Tree", "ele" + inFileTag[p] + "Tree"));
    eleTrees.back()->SetDirectory(0);
    eleTrees.back()->Branch("runnr", &runnumber, "runnr/i");
    eleTrees.back()->Branch("eventnr", &eventnumber, "eventnr/i");
    eleTrees.back()->Branch("lumiSec", &luminosityBlock, "lumiSec/i");
    eleTrees.back()->Branch("pvsize", &pvsize, "pvsize/I");
    eleTrees.back()->Branch("mcWeight", &mcWeight, "mcWeight/F");
    eleTrees.back()->Branch("puWeight", &puWeight, "puWeight/F");
    eleTrees.back()->Branch("trueMass", &trueMass, "trueMass/F");
    eleTrees.back()->Branch("evtRegion", &evtRegion, "evtRegion/I");
    eleTrees.back()->Branch("mass", &hHMass, "mass/F");
    eleTrees.back()->Branch("ele1Et", &ele1Et, "ele1Et/F");
    eleTrees.back()->Branch("ele1ScE", &ele1ScE, "ele1ScE/F");
    eleTrees.back()->Branch("ele1Eta", &ele1Eta, "ele1Eta/F");
    eleTrees.back()->Branch("ele2Et", &ele2Et, "ele2Et/F");
    eleTrees.back()->Branch("ele2ScE", &ele2ScE, "ele2ScE/F");
    eleTrees.back()->Branch("ele2Eta", &ele2Eta, "ele2Eta/F");
    eleTrees.back()->Branch("massDrMatched", &hHMassDrMatched, "massDrMatched/F");
    eleTrees.back()->Branch("evtRegionDrMatched", &evtRegionDrMatched, "evtRegionDrMatched/I");
    eleTrees.back()->Branch("ele1EtDrMatched", &ele1EtDrMatched, "ele1EtDrMatched/F");
    eleTrees.back()->Branch("ele1ScEDrMatched", &ele1ScEDrMatched, "ele1ScEDrMatched/F");
    eleTrees.back()->Branch("ele1EtaDrMatched", &ele1EtaDrMatched, "ele1EtaDrMatched/F");
    eleTrees.back()->Branch("ele2EtDrMatched", &ele2EtDrMatched, "ele2EtDrMatched/F");
    eleTrees.back()->Branch("ele2ScEDrMatched", &ele2ScEDrMatched, "ele2ScEDrMatched/F");
    eleTrees.back()->Branch("ele2EtaDrMatched", &ele2EtaDrMatched, "ele2EtaDrMatched/F");
    eleTrees.back()->Branch("genEle1E", &genEle1E, "genEle1E/F");
    eleTrees.back()->Branch("genEle2E", &genEle2E, "genEle2E/F");
    eleTrees.back()->Branch("unstGenEle1E", &unstGenEle1E, "unstGenEle1E/F");
    eleTrees.back()->Branch("unstGenEle2E", &unstGenEle2E, "unstGenEle2E/F");

    bool etShift = false;
    if (etShiftData && p == DATA)
      etShift = true;

    unsigned int oneHeepEleEv = 0;
    unsigned int twoHeepEleEv = 0;
    unsigned int threeHeepEleEv = 0;
    unsigned int fourHeepEleEv = 0;
    unsigned int fiveHeepEleEv = 0;

    unsigned int evCounter = 0;
    //Long64_t nbytes = 0, nb = 0;
    /////////////////////////////////////////////////////////////////////////
    // loop over events
    /////////////////////////////////////////////////////////////////////////
    //for (Long64_t jentry=0; (jentry < 100 && jentry < nentries); ++jentry) { // for tests
    for (Long64_t jentry=0; jentry < nentries; ++jentry) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      thetree->GetEntry(jentry);
      if (jentry % 50000 == 0 ) cout << "entry " << jentry << endl;

      if (p > DATA) trueMass = genelemom_mass[0];

      // preset to default values
      evtRegion = -1;
      hHMass = 0.;
      ele1Et = -10.;
      ele2Et = -10.;
      ele1ScE = -10.;
      ele2ScE = -10.;
      ele1Eta = -10.;
      ele2Eta = -10.;
      evtRegionDrMatched = -1;
      hHMassDrMatched = 0.;
      ele1EtDrMatched = -10.;
      ele2EtDrMatched = -10.;
      ele1ScEDrMatched = -10.;
      ele2ScEDrMatched = -10.;
      ele1EtaDrMatched = -10.;
      ele2EtaDrMatched = -10.;
      genEle1E = -10.;
      genEle2E = -10.;
      unstGenEle1E = -10.;
      unstGenEle2E = -10.;

      // first correct the energy
      //CorrectEnergy(); 

      // trigger fired?
      int prescale = 0;
      //if (p == DATA && Trigger(prescale) < 1) continue;
      if (p == DATA && Trigger(prescale) < 1) continue;
    
      // at least two gsf electrons
      if (gsf_size < 2) continue;

      // set the PU weight
      if (usePUInfo && p > DATA) puWeight = puWeights->GetBinContent(puWeights->FindBin(trueNVtx));

      //if (p >= DY20 && n_pvValid < 1) continue;
      ////////////////////////////////////////////////////////////////////////
      // find the two highest pt GSF and HEEP electrons
      ////////////////////////////////////////////////////////////////////////
      unsigned int heepCounter = 0;
      int iHeep1 = -1;
      int iHeep2 = -1;
      int iHeepDrMatched1 = -1;
      int iHeepDrMatched2 = -1;
      float highestHeepEt = ptcut;
      for (int n = 0; n < gsf_size; ++n) {
        bool passHeep = PassHEEP(n);
        if (passHeep) ++heepCounter;
        if (highestHeepEt < gsf_gsfet[n] && passHeep) {
          iHeep1 = n;
          highestHeepEt = gsf_gsfet[n];
        }
      }
      highestHeepEt = ptcut;
      for (int m = 0; m < gsf_size; ++m) {
        if (highestHeepEt < gsf_gsfet[m] && PassHEEP(m) && m != iHeep1) {
          iHeep2 = m;
          highestHeepEt = gsf_gsfet[m];
        }
      }

      int bestMatchE1 = -1;
      int bestMatchE2 = -1;
      double bestDrE1 = 1000.;
      double bestDrE2 = 1000.;
      if (p > 0) {
        for (int n = 0; n < gsf_size; ++n) {
          if (!PassHEEP(n)) continue;
          double dr = 1000.;
          int index = GenEleDRMatch(n, dr);
          if (index == 0 && dr < bestDrE1) {
            bestMatchE1 = n;
            bestDrE1 = dr;
          }
          if (index == 1 && dr < bestDrE2) {
            bestMatchE2 = n;
            bestDrE2 = dr;
          }
        }
        //if (bestMatchE1 == bestMatchE2 || bestMatchE1 > 1 || bestMatchE2 > 1 || bestMatchE1 == -1 || bestMatchE2 == -1 || bestDrE1 > 0.1 || bestDrE2 > 0.1) std::cout << "genele0 best matched by gsfele" << bestMatchE1  << " with dR=" << bestDrE1 << ". genele1 best matched by gsfele" << bestMatchE2  << " with dR=" << bestDrE2 << std::endl;
      }

      if (iHeep1 < 0 || iHeep2 < 0) continue;

      bool fillTree = false;
      // fill the HEEP-HEEP cases
      float heepHeepMass = CalcInvariantMass(iHeep1, iHeep2, etShift);
      if (heepHeepMass > 0.) {
        if (heepCounter > 4) ++fiveHeepEleEv;
        else if (heepCounter > 3) ++fourHeepEleEv;
        else if (heepCounter > 2) ++threeHeepEleEv;
        else if (heepCounter > 1) ++twoHeepEleEv;
        else if (heepCounter > 0) ++oneHeepEleEv;
      }
      if (heepHeepMass > massCut) {
        if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iHeep2]) < 2.5) {
          if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iHeep2]) > 1.56) {
            evtRegion = 2;
            hHMass = heepHeepMass;
            ele1Et = gsf_gsfet[iHeep1];
            ele2Et = gsf_gsfet[iHeep2];
            ele1ScE = gsfsc_e[iHeep1];
            ele2ScE = gsfsc_e[iHeep2];
            ele1Eta = gsfsc_eta[iHeep1];
            ele2Eta = gsfsc_eta[iHeep2];
            genEle1E = genele_e[0];
            genEle2E = genele_e[1];
            unstGenEle1E = unstableGenEle_e[0];
            unstGenEle2E = unstableGenEle_e[1];
          } else {
            if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iHeep2]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iHeep2]) < 1.442)) {
              evtRegion = 1;
              hHMass = heepHeepMass;
              ele1Et = gsf_gsfet[iHeep1];
              ele2Et = gsf_gsfet[iHeep2];
              ele1ScE = gsfsc_e[iHeep1];
              ele2ScE = gsfsc_e[iHeep2];
              ele1Eta = gsfsc_eta[iHeep1];
              ele2Eta = gsfsc_eta[iHeep2];
              genEle1E = genele_e[0];
              genEle2E = genele_e[1];
              unstGenEle1E = unstableGenEle_e[0];
              unstGenEle2E = unstableGenEle_e[1];
            } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iHeep2]) < 1.442) {
              evtRegion = 0;
              hHMass = heepHeepMass;
              ele1Et = gsf_gsfet[iHeep1];
              ele2Et = gsf_gsfet[iHeep2];
              ele1ScE = gsfsc_e[iHeep1];
              ele2ScE = gsfsc_e[iHeep2];
              ele1Eta = gsfsc_eta[iHeep1];
              ele2Eta = gsfsc_eta[iHeep2];
              genEle1E = genele_e[0];
              genEle2E = genele_e[1];
              unstGenEle1E = unstableGenEle_e[0];
              unstGenEle2E = unstableGenEle_e[1];
            }
          }
          fillTree = true;
        }
      }

      if (p > 0 && bestMatchE1 >= 0 && bestMatchE2 >= 0 && bestDrE1 <= drMax && bestDrE2 <= drMax) {
        // fill the HEEP-HEEP cases for deltaR matched electrons
        float heepHeepMassDrMatched = CalcInvariantMass(bestMatchE1, bestMatchE2, etShift);
        if (heepHeepMassDrMatched > massCut) {
          if (fabs(gsfsc_eta[bestMatchE1]) < 2.5 && fabs(gsfsc_eta[bestMatchE2]) < 2.5) {
            if (fabs(gsfsc_eta[bestMatchE1]) > 1.56 && fabs(gsfsc_eta[bestMatchE2]) > 1.56) {
              evtRegionDrMatched = 2;
              hHMassDrMatched = heepHeepMassDrMatched;
              ele1EtDrMatched = gsf_gsfet[bestMatchE1];
              ele2EtDrMatched = gsf_gsfet[bestMatchE2];
              ele1ScEDrMatched = gsfsc_e[bestMatchE1];
              ele2ScEDrMatched = gsfsc_e[bestMatchE2];
              ele1EtaDrMatched = gsfsc_eta[bestMatchE1];
              ele2EtaDrMatched = gsfsc_eta[bestMatchE2];
            } else {
              if ((fabs(gsfsc_eta[bestMatchE1]) < 1.442 && fabs(gsfsc_eta[bestMatchE2]) > 1.56) || (fabs(gsfsc_eta[bestMatchE1]) > 1.56 && fabs(gsfsc_eta[bestMatchE2]) < 1.442)) {
                evtRegionDrMatched = 1;
                hHMassDrMatched = heepHeepMassDrMatched;
                ele1EtDrMatched = gsf_gsfet[bestMatchE1];
                ele2EtDrMatched = gsf_gsfet[bestMatchE2];
                ele1ScEDrMatched = gsfsc_e[bestMatchE1];
                ele2ScEDrMatched = gsfsc_e[bestMatchE2];
                ele1EtaDrMatched = gsfsc_eta[bestMatchE1];
                ele2EtaDrMatched = gsfsc_eta[bestMatchE2];
              } else if (fabs(gsfsc_eta[bestMatchE1]) < 1.442 && fabs(gsfsc_eta[bestMatchE2]) < 1.442) {
                evtRegionDrMatched = 0;
                hHMassDrMatched = heepHeepMassDrMatched;
                ele1EtDrMatched = gsf_gsfet[bestMatchE1];
                ele2EtDrMatched = gsf_gsfet[bestMatchE2];
                ele1ScEDrMatched = gsfsc_e[bestMatchE1];
                ele2ScEDrMatched = gsfsc_e[bestMatchE2];
                ele1EtaDrMatched = gsfsc_eta[bestMatchE1];
                ele2EtaDrMatched = gsfsc_eta[bestMatchE2];
              }
            }
            fillTree = true;
          }
        }
      }
      if (fillTree) {
        eleTrees.back()->Fill();
        ++evCounter;
      }
    } // end of loop over events

    cout << "Number of selected events: " << evCounter << endl;
    cout << "Number of 1 / 2 / 3 / 4 / 5 HEEP events: " << oneHeepEleEv << " / " << twoHeepEleEv << " / " << threeHeepEleEv << " / " << fourHeepEleEv << " / " << fiveHeepEleEv << " / " << endl;

    // write root file with good HEEP-HEEP event data
    goodEvFile->cd();
    eleTrees.back()->Write();
    puMc->Write();
    puMcNorm->Write();
    puWeights->Write();
  } // end of loop over input files
   puData->Write();
   puDataNorm->Write();

   timer.Stop();
   timer.Print();
} //end of method


double
EScaleTreeMaker::CalcInvariantMass (const int& iEle1, const int& iEle2, bool &etShift)
{
  float dataEtShiftFactorEB = 1.;
  float dataEtShiftFactorEE = 1.;

  if (etShift) {
    dataEtShiftFactorEB = 1.0036;
    dataEtShiftFactorEE = 1.0256;
  }

  TLorentzVector ele1;
  TLorentzVector ele2;

  Float_t et1 = gsf_gsfet[iEle1];
  Float_t et2 = gsf_gsfet[iEle2];
  // correct energy in data
  (fabs(gsfsc_eta[iEle1]) > 1.56) ? et1 *= dataEtShiftFactorEE : et1 *= dataEtShiftFactorEB;
  (fabs(gsfsc_eta[iEle2]) > 1.56) ? et2 *= dataEtShiftFactorEE : et2 *= dataEtShiftFactorEB;

  ele1.SetPtEtaPhiE(et1, gsf_eta[iEle1], gsf_phi[iEle1], (et1 * cosh(gsf_eta[iEle1])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[iEle2], gsf_phi[iEle2], (et2 * cosh(gsf_eta[iEle2])));
  
  return (ele1+ele2).Mag();
}

bool
EScaleTreeMaker::PassHEEP(const int &n)
{
  const int selection = 1;

  // HEEP v4.0
  // barrel
  float bar_et = 25.;
  float bar_hoE = 0.05;
  float bar_DEta = 0.005;
  float bar_DPhi = 0.06;
  float bar_e2x5e5x5 = 0.94;
  float bar_e1x5e5x5 = 0.83;
  float bar_isoEcalHcal1_1 = 2.;
  float bar_isoEcalHcal1_2 = 0.03;
  float bar_isoEcalHcalRho = 0.28;
  float bar_isoTrack = 5.;
  float bar_dxy = 0.02;  // only for HEEP v4.1
  int bar_missInnerHits = 0;

  // endcap
  float end_et = 25.;
  float end_hoE = 0.05;
  float end_DEta = 0.007 ;
  float end_DPhi = 0.06;
  float end_sigmaietaieta = 0.03;
  float end_isoEcalHcal1_1_1 = 2.5;
  float end_isoEcalHcal1_1_2 = 1.;
  float end_isoEcalHcal1_2 = 0.03;
  float end_isoEcalHcalRho = 0.28;
  float end_isoTrack = 5.;
  float end_dxy = 0.05;   // only for HEEP v4.1
  int end_missInnerHits = 0;

  // HEEP v4.1
  if (selection == 1) {
    bar_missInnerHits = 1;
    end_missInnerHits = 1;
  }

  // HEEP
  // barrel
  if (fabs(gsfsc_eta[n]) < 1.442
      && gsf_gsfet[n] > bar_et
      && gsf_isecaldriven[n]
      && gsf_hovere[n] < bar_hoE
      && fabs(gsf_deltaeta[n]) < bar_DEta
      && fabs(gsf_deltaphi[n]) < bar_DPhi
      && (gsf_e2x5overe5x5[n] > bar_e2x5e5x5 || gsf_e1x5overe5x5[n] > bar_e1x5e5x5)
      && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * gsf_gsfet[n] + bar_isoEcalHcalRho * rho)
      && gsf_trackiso[n] < bar_isoTrack
      && gsf_nLostInnerHits[n] <= bar_missInnerHits
     ) if ((selection == 1 && fabs(gsf_dxy_firstPVtx[n]) <= bar_dxy) || selection != 1) return true;
  
  // endcap
  if ((fabs(gsfsc_eta[n]) > 1.56 && fabs(gsfsc_eta[n]) < 2.5)
      && gsf_gsfet[n] > end_et
      && gsf_isecaldriven[n]
      && gsf_hovere[n] < end_hoE
      && fabs(gsf_deltaeta[n]) < end_DEta
      && fabs(gsf_deltaphi[n]) < end_DPhi
      && gsf_sigmaIetaIeta[n] < end_sigmaietaieta
      && ((gsf_gsfet[n] < 50. && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (end_isoEcalHcal1_1_1 + end_isoEcalHcalRho * rho))
         || (gsf_gsfet[n] >= 50. && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (end_isoEcalHcal1_1_2 + end_isoEcalHcal1_2 * gsf_gsfet[n] + end_isoEcalHcalRho * rho)))
      && gsf_trackiso[n] < end_isoTrack
      && gsf_nLostInnerHits[n] <= end_missInnerHits
     ) if ((selection == 1 && fabs(gsf_dxy_firstPVtx[n]) <= end_dxy) || selection != 1) return true;
  return false;
}

int 
EScaleTreeMaker::Trigger(int &prescale)
{
  // trigger selection 2012
  int triggerBit = 0;
  if (prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 1 && HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 1) {
    prescale = prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
    triggerBit = HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  } 
  return triggerBit;

  cout << "Prescale alert! No unprescaled trigger found. Dropping this event." << endl;
  prescale = 0;
  return 0;
}

int
EScaleTreeMaker::GenEleDRMatch(const int& iEle, double &drMin)
{
  TLorentzVector gsfEle;
  TLorentzVector genEle;

  Float_t et = gsf_gsfet[iEle];

  gsfEle.SetPtEtaPhiE(et, gsf_eta[iEle], gsf_phi[iEle], (et * cosh(gsf_eta[iEle])));
  int iMinDr = -1;
  int i = 0;
  for (; i < 2; ++i) {
    genEle.SetPtEtaPhiE(genele_pt[i], genele_eta[i], genele_phi[i], genele_e[i]);
    double dr = gsfEle.DeltaR(genEle);
    if (fabs(dr) > drMin) continue;
    drMin = fabs(dr);
    iMinDr = i;
  }
  
  return iMinDr;
}

void
EScaleTreeMaker::CorrectEnergy() 
{
  for (int n = 0; n < gsf_size; ++n)
     gsf_gsfet[n] = gsfsc_e[n] * sin(gsf_theta[n]);
}

