#define invariantMass_cxx
#include "invariantMass.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>

#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TKey.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"

void InvariantMass::Loop()
{
  // parameters //////////////////////////////////////////////////////////////
  const int ptcut = 35;
  const float massCut = 50;
  const float lumiCorrFactor = 1.07895; // from 60GeV to 120GeV with 2177pb-1
                                     // not used if calcLumiCorrFactor=true
  bool calcLumiCorrFactor = true;    // first run over events once to
                                     // calculate correction and the second 
                                     // time with the factor applied
  float normToZPeakFrom = 60;        // lower mass for normalization
  float normToZPeakTo = 120;         // upper mass for normalization

  string fileNamePrefix = "invMassHistos";   // lumi and .root well be added
  //string fileNamePrefix = "test";

  const int massMin = 50;            // minimum Mass
  const int massMax = 1550;          // maximum Mass

  // define binning for histograms (variable or constant): 
  // [first]=max. Mass (in GeV) up to which [second] (in GeV) binning is used
  vector<pair<float, float> > binning;
  // VARIABLE BINNING
//  binning.push_back(make_pair(100, 1));
//  binning.push_back(make_pair(500, 10));
//  binning.push_back(make_pair(1500, 50));
  // CONSTANT BINNING
  binning.push_back(make_pair(massMax, 1));

  const int ratioRebin = 50;      // rebinning for FR ratio histograms
                                  // should be chosen compatible with binning
  ////////////////////////////////////////////////////////////////////////////


  vector<float> bins;
  bins.push_back(massMin);
  for (vector<pair<float,float> >::iterator it = binning.begin(); it < binning.end(); ++it) {
    while (bins.back() < it->first)
      bins.push_back(bins.back() + it->second);
  }
  if (bins.back() < massMax)
    bins.push_back(massMax);
  int nBins = bins.size() - 1;
  Float_t binArray[nBins + 1];
  for (int i = 0; i <= nBins; ++i)
    binArray[i] = bins.at(i);

  TH1::SetDefaultSumw2(kTRUE);

  vector<vector<TH1F *> > histosGsfGsfMassFR;
  vector<vector<TH1F *> > histosGsfGsfMassNoHeepFR;
  vector<vector<TH1F *> > histosHeepGsfMassNoHeepFR;
  vector<vector<TH1F *> > histosSectionDYCombinedHH;
  vector<vector<TH1F *> > histosSectionDYCombinedHGNH;
  vector<vector<TH1F *> > histosSectionDYCombinedGGNH;
  vector<vector<TH1F *> > histosSectionDYCombinedHGNHFR;
  vector<vector<TH1F *> > histosSectionDYCombinedGGNHFR;

  stringstream sStream;
  sStream << "# of events / " << ((massMax - massMin) / nBins) << "GeV/c^{2}";

  /////////////////////////////////////////////////////////////////////////
  // input files
  /////////////////////////////////////////////////////////////////////////
  vector<pair<TFile *, double> > input;
  vector<string> datasetTitle;
  vector<float> minPtDY;
  vector<float> minPtG;
  //float lumi = 204;
  //input.push_back(make_pair(new TFile("/user/treis/data2011/gsfcheckertree203_6pb-1.root","read"), 1 / lumi));
  float lumi = 2928;
  input.push_back(make_pair(new TFile("/user/lathomas/CMSSW_4_2_7/src/2928pb-1.root","read"), 1 / lumi));
  datasetTitle.push_back("Data");
  const unsigned int DATA = 0;

  input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-20_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2_GEN-SIM-RECO_HEEPSkim2ElePt30_gct1_6.root","read"), 7.54408E-4));
  datasetTitle.push_back("DY > 20GeV");
  minPtDY.push_back(20);
  const unsigned int DY20 = 1;

  input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-120_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_GEN-SIM-RECO_HEEPSkim2ElePt30_gct1_6.root","read"), 1.80678E-4));
  datasetTitle.push_back("DY > 120GeV");
  minPtDY.push_back(120);
  const unsigned int DY120 = 2;

  input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-200_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_GEN-SIM-RECO_HEEPSkim2ElePt30_gct1_6.root","read"), 2.18298E-5));
  datasetTitle.push_back("DY > 200GeV");
  minPtDY.push_back(200);
  const unsigned int DY200 = 3;

  input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-500_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_GEN-SIM-RECO_HEEPSkim2ElePt30_gct1_6.root","read"), 6.17792E-7));
  datasetTitle.push_back("DY > 500GeV");
  minPtDY.push_back(500);
  const unsigned int DY500 = 4;

  input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-800_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_GEN-SIM-RECO_HEEPSkim2ElePt30_gct1_6.root","read"), 7.47055E-8));
  datasetTitle.push_back("DY > 800GeV");
  minPtDY.push_back(800);
  const unsigned int DY800 = 5;

  input.push_back(make_pair(new TFile("/user/treis/mcsamples/TT_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_AODSIM_HEEPSkim2ElePt30.root","read"), 1.44545E-4));
  datasetTitle.push_back("ttbar");
  const unsigned int TTBAR = 6;

  input.push_back(make_pair(new TFile("/user/lathomas/CMSSW_4_2_3/src/UserCode/HEEPSkims/test/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola-Summer11-PU_S4_START42_V11-v1-AODSIM_HEEPSkim2ElePt30/res/TOTAL.root","read"), 0.000607645));
  datasetTitle.push_back("W+jets");
  const unsigned int WJETS = 7;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 7.6306E-3));
  datasetTitle.push_back("Gamma 30to50GeV");
  minPtG.push_back(30);
  const unsigned int G30T50 = 8;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 1.3365E-3));
  datasetTitle.push_back("Gamma 50to80GeV");
  minPtG.push_back(50);
  const unsigned int G50T80 = 9;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 2.1850E-4));
  datasetTitle.push_back("Gamma 80to120GeV");
  minPtG.push_back(80);
  const unsigned int G80T120 = 10;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 4.0307E-5));
  datasetTitle.push_back("Gamma 120to170GeV");
  minPtG.push_back(120);
  const unsigned int G120T170 = 11;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 1.0942E-5));
  datasetTitle.push_back("Gamma 170to300GeV");
  minPtG.push_back(170);
  const unsigned int G170T300 = 12;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-300to470_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 7.1887E-7));
  datasetTitle.push_back("Gamma 300to470GeV");
  minPtG.push_back(300);
  const unsigned int G300T470 = 13;
  
  input.push_back(make_pair(new TFile("/user/treis/mcsamples/G_Pt-470to800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_2gsfpt30skim_gct1_6.root","read"), 6.3386E-8));
  datasetTitle.push_back("Gamma 470to800GeV");
  minPtG.push_back(470);
  const unsigned int G470T800 = 14;
  ////////////////////////////////////////////////////////////////////////////
  
  minPtDY.push_back(100000); // set ridiculously high to catch the tail of the last DY sample

  vector<TString> acroSuffix;
  acroSuffix.push_back("");
  acroSuffix.push_back("BB");
  acroSuffix.push_back("BE");
  acroSuffix.push_back("EE");
  acroSuffix.push_back("EES");
  acroSuffix.push_back("EEO");

  vector<TString> suffix;
  suffix.push_back("");
  suffix.push_back(" EB-EB");
  suffix.push_back(" EB-EE");
  suffix.push_back(" EE-EE");
  suffix.push_back(" same EE-EE");
  suffix.push_back(" opposite EE-EE");

  vector<string> folders; // folder names in output root file

  stringstream ssOutfile;
  ssOutfile << fileNamePrefix << lumi << "pb-1.root";
  TFile *outFile = new TFile(ssOutfile.str().c_str(), "recreate");
  /////////////////////////////////////////////////////////////////////////
  // loop over files
  /////////////////////////////////////////////////////////////////////////
  for (unsigned int p = 0; p < input.size(); ++p) {
  //for (unsigned int p = 0; p < 1; ++p) {
    string infile(input[p].first->GetName());
    folders.push_back(infile.substr(infile.rfind("/") + 1, infile.rfind(".root") - infile.rfind("/") - 1));
    cout << "file " << infile << endl;
    input[p].first->cd();
    TTree *thetree = (TTree*)input[p].first->Get("gsfcheckerjob/tree");
    Init(thetree);
    Long64_t nentries = (*thetree).GetEntries();
    cout << nentries << " entries" << endl;

    // normalize bg
    if (p == 1 && !calcLumiCorrFactor) lumi *= lumiCorrFactor;

    /////////////////////////////////////////////////////////////////////////
    // set up histograms
    /////////////////////////////////////////////////////////////////////////
    vector<TH1F *> histoHeepHeepMass;
    vector<TH1F *> histoGsfGsfMass;
    vector<TH1F *> histoGsfGsfMassNoHeep;
    vector<TH1F *> histoHeepGsfMass;
    vector<TH1F *> histoHeepGsfMassNoHeep;
    vector<TH1F *> histoSectionDYCombinedHH;
    vector<TH1F *> histoSectionDYCombinedHGNH;
    vector<TH1F *> histoSectionDYCombinedGGNH;
    vector<TH1F *> histoGsfGsfMassFR;
    vector<TH1F *> histoGsfGsfMassNoHeepFR;
    vector<TH1F *> histoHeepGsfMassFR;
    vector<TH1F *> histoHeepGsfMassNoHeepFR;
    vector<TH1F *> histoSectionDYCombinedHGNHFR;
    vector<TH1F *> histoSectionDYCombinedGGNHFR;
    vector<TH1F *> histoHeepGsfMassNoHeepNoHOverEFR;
    vector<TH1F *> histoTAndPHeepGsfMassNoHeepFR;
    for (unsigned int i = BBBE; i <= EEO; ++i) {
      // setup HEEP-HEEP histograms
      histoHeepHeepMass.push_back(new TH1F("histoHeepHeepMass" + acroSuffix.at(i), "Dielectron invariant mass" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoHeepHeepMass.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup GSF-GSF histograms
      histoGsfGsfMass.push_back(new TH1F("histoGsfGsfMass" + acroSuffix.at(i), "GSF invariant mass" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoGsfGsfMass.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup GSF-GSF non HEEP histograms
      histoGsfGsfMassNoHeep.push_back(new TH1F("histoGsfGsfMassNoHeep" + acroSuffix.at(i), "GSF(non HEEP) invariant mass" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoGsfGsfMassNoHeep.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup HEEP-GSF histograms
      histoHeepGsfMass.push_back(new TH1F("histoHeepGsfMass" + acroSuffix.at(i), "Heep-GSF invariant mass" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoHeepGsfMass.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup HEEP-GSF(non HEEP) histograms
      histoHeepGsfMassNoHeep.push_back(new TH1F("histoHeepGsfMassNoHeep" + acroSuffix.at(i), "GSF(non HEEP) invariant mass" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoHeepGsfMassNoHeep.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup histogram sections for combined DY histograms
      if (p >= DY20 && p <= DY800) {
        histoSectionDYCombinedHH.push_back(new TH1F("histoSectionDYCombinedHH" + acroSuffix.at(i), "Section for DY combined HEEP-HEEP" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
        histoSectionDYCombinedHH.back()->GetYaxis()->SetTitle(sStream.str().c_str());

        histoSectionDYCombinedHGNH.push_back(new TH1F("histoSectionDYCombinedHGNH" + acroSuffix.at(i), "Section for DY combined HEEP-GSF non HEEP" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
        histoSectionDYCombinedHGNH.back()->GetYaxis()->SetTitle(sStream.str().c_str());


        histoSectionDYCombinedGGNH.push_back(new TH1F("histoSectionDYCombinedGGNH" + acroSuffix.at(i), "Section for DY combined GSF-GSF non HEEP" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
        histoSectionDYCombinedGGNH.back()->GetYaxis()->SetTitle(sStream.str().c_str());
      }

      // set up fake rate histograms
      if (i > EE) continue;

      histoGsfGsfMassFR.push_back(new TH1F("histoGsfGsfMassFR" + acroSuffix.at(i), "GSF invariant mass with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoGsfGsfMassFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoGsfGsfMassNoHeepFR.push_back(new TH1F("histoGsfGsfMassNoHeepFR" + acroSuffix.at(i), "GSF(non HEEP) invariant mass with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoGsfGsfMassNoHeepFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoHeepGsfMassFR.push_back(new TH1F("histoHeepGsfMassFR" + acroSuffix.at(i), "Heep-GSF invariant mass with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoHeepGsfMassFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoHeepGsfMassNoHeepFR.push_back(new TH1F("histoHeepGsfMassNoHeepFR" + acroSuffix.at(i), "GSF(non HEEP) invariant mass with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoHeepGsfMassNoHeepFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      if (p >= DY20 && p <= DY800) {
        histoSectionDYCombinedHGNHFR.push_back(new TH1F("histoSectionDYCombinedHGNHFR" + acroSuffix.at(i), "Section for DY combined GSF-GSF non HEEP with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
        histoSectionDYCombinedHGNHFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

        histoSectionDYCombinedGGNHFR.push_back(new TH1F("histoSectionDYCombinedGGNHFR" + acroSuffix.at(i), "Section for DY combined GSF-GSF non HEEP with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
        histoSectionDYCombinedGGNHFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());
      }

      if (i > BE) continue;

      // setup HEEP-GSF(non HEEP) histograms for tag and probe studies for Laurent
      histoHeepGsfMassNoHeepNoHOverEFR.push_back(new TH1F("histoHeepGsfMassNoHeepNoHOverEFR" + acroSuffix.at(i), "HEEP-GSF (non HEEP) with fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoHeepGsfMassNoHeepNoHOverEFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoTAndPHeepGsfMassNoHeepFR.push_back(new TH1F("histoTAndPHeepGsfMassNoHeepFR" + acroSuffix.at(i), "Invariant mass for tag and probe with additional cuts" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
      histoTAndPHeepGsfMassNoHeepFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());
    }

    sStream.str("");
  
    /////////////////////////////////////////////////////////////////////////
    // setup Pz histograms
    /////////////////////////////////////////////////////////////////////////
    vector<TH1F *> histoPZ;
    histoPZ.push_back(new TH1F("histoPZ", "Pz of the initial particle;Pz (GeV/c)", nBins, -massMax, massMax));
    histoPZ.push_back(new TH1F("histoPZSameEE", "Pz of the initial particle for same endcap events;Pz (GeV/c)", nBins, -massMax, massMax));
    histoPZ.push_back(new TH1F("histoPZOppositeEE", "Pz of the initial particle for opposite endcap events;Pz (GeV/c)", nBins, -massMax, massMax));
  
    sStream << "# of events / " << (2 * massMax / nBins) << "GeV/c";
    vector<TH1F *>::iterator iter;
    for (iter = histoPZ.begin(); iter < histoPZ.end(); ++iter) {
      (*iter)->GetYaxis()->SetTitle(sStream.str().c_str());
    }
    sStream.str("");

    //Long64_t nbytes = 0, nb = 0;
    /////////////////////////////////////////////////////////////////////////
    // loop over events
    /////////////////////////////////////////////////////////////////////////
    //for (Long64_t jentry=0; (jentry < 5000 && jentry < nentries); ++jentry) { // for tests
    for (Long64_t jentry=0; jentry < nentries; ++jentry) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      thetree->GetEntry(jentry);
      if (jentry % 50000 == 0 ) cout << "entry " << jentry << endl;
      // at least two gsf electrons
      if (gsf_size < 2) continue;

      ////////////////////////////////////////////////////////////////////////
      // find the two highest pt GSF and HEEP electrons
      ////////////////////////////////////////////////////////////////////////
      int iGsf1 = -1;
      int iGsf2 = -1;
      int iHeep1 = -1;
      int iHeep2 = -1;
      int iGsfNoHeep1 = -1;
      int iGsfNoHeep2 = -1;
      float highestEt = ptcut;
      float highestHeepEt = ptcut;
      float highestNoHeepEt = ptcut;
      for (int n = 0; n < gsf_size; ++n) {
        if (highestEt < gsf_gsfet[n]) {
          iGsf1 = n;
          highestEt = gsf_gsfet[n];
        }
        if (highestHeepEt < gsf_gsfet[n] && gsfpass_HEEP[n] && gsf_nLostInnerHits[n] == 0) {
          iHeep1 = n;
          highestHeepEt = gsf_gsfet[n];
        }
        if (highestNoHeepEt < gsf_gsfet[n] && (!gsfpass_HEEP[n] || gsf_nLostInnerHits[n] > 0)) {
          iGsfNoHeep1 = n;
          highestNoHeepEt = gsf_gsfet[n];
        }
      }
      highestEt = ptcut;
      highestHeepEt = ptcut;
      highestNoHeepEt = ptcut;
      for (int m = 0; m < gsf_size; ++m) {
        if (highestEt < gsf_gsfet[m] && m != iGsf1) {
          iGsf2 = m;
          highestEt = gsf_gsfet[m];
        }
        if (highestHeepEt < gsf_gsfet[m] && gsfpass_HEEP[m] && gsf_nLostInnerHits[m] == 0 && m != iHeep1) {
          iHeep2 = m;
          highestHeepEt = gsf_gsfet[m];
        }
        if (highestNoHeepEt < gsf_gsfet[m] && (!gsfpass_HEEP[m] || gsf_nLostInnerHits[m] > 0) && m != iGsfNoHeep1) {
          iGsfNoHeep2 = m;
          highestNoHeepEt = gsf_gsfet[m];
        }
      }

      if (iGsf2 < 0 || iGsf1 < 0) continue; // not enough good GSF electrons found 

      ////////////////////////////////////////////////////////////////////////
      // fill the histograms
      ////////////////////////////////////////////////////////////////////////
      // fill the GSF-GSF cases
      double gsfGsfMass = CalcInvariantMass(iGsf1, iGsf2);
      if (gsfGsfMass > massCut) {
        if (fabs(gsfsc_eta[iGsf1]) < 2.5 && fabs(gsfsc_eta[iGsf2]) < 2.5) {
          double fakeRate1 = FakeRate(gsf_gsfet[iGsf1], gsfsc_eta[iGsf1]);
          double fakeRate2 = FakeRate(gsf_gsfet[iGsf2], gsfsc_eta[iGsf2]);
          if (fabs(gsfsc_eta[iGsf1]) > 1.56 && fabs(gsfsc_eta[iGsf2]) > 1.56) {
            histoGsfGsfMass.at(EE)->Fill(gsfGsfMass, input.at(p).second * lumi);
            // fill the pz of the initial particle for EE-EE events
            double gsfGsfPZ = CalcPz(iGsf1, iGsf2);
            if (gsfsc_eta[iGsf1] * gsfsc_eta[iGsf2] > 0) {
              histoGsfGsfMass.at(EES)->Fill(gsfGsfMass, input.at(p).second * lumi);
              histoPZ.at(1)->Fill(gsfGsfPZ, input.at(p).second * lumi);
            } else {
              histoGsfGsfMass.at(EEO)->Fill(gsfGsfMass, input.at(p).second * lumi);
              histoPZ.at(2)->Fill(gsfGsfPZ, input.at(p).second * lumi);
            }
            if (PassFRPreSel(iGsf1, gsfsc_eta[iGsf1]) && PassFRPreSel(iGsf2, gsfsc_eta[iGsf2]))
              histoGsfGsfMassFR.at(EE)->Fill(gsfGsfMass, input.at(p).second * lumi * fakeRate1 * fakeRate2);
          } else if ((fabs(gsfsc_eta[iGsf1]) < 1.442 && fabs(gsfsc_eta[iGsf2]) > 1.56) || (fabs(gsfsc_eta[iGsf1]) > 1.56 && fabs(gsfsc_eta[iGsf2]) < 1.442)) {
            histoGsfGsfMass.at(BE)->Fill(gsfGsfMass, input.at(p).second * lumi);
            if (PassFRPreSel(iGsf1, gsfsc_eta[iGsf1]) && PassFRPreSel(iGsf2, gsfsc_eta[iGsf2]))
              histoGsfGsfMassFR.at(BE)->Fill(gsfGsfMass, input.at(p).second * lumi * fakeRate1 * fakeRate2);
          } else if (fabs(gsfsc_eta[iGsf1]) < 1.442 && fabs(gsfsc_eta[iGsf2]) < 1.442) {
            histoGsfGsfMass.at(BB)->Fill(gsfGsfMass, input.at(p).second * lumi);
            if (PassFRPreSel(iGsf1, gsfsc_eta[iGsf1]) && PassFRPreSel(iGsf2, gsfsc_eta[iGsf2]))
              histoGsfGsfMassFR.at(BB)->Fill(gsfGsfMass, input.at(p).second * lumi * fakeRate1 * fakeRate2);
          }
        }
      }

      if (iGsfNoHeep2 >= 0 && iGsfNoHeep1 >= 0 && iHeep1 < 0 && iHeep2 < 0) {
        // fill the GSF-GSF cases that do not pass the HEEP selection
        double gsfGsfMassNoHeep = CalcInvariantMass(iGsfNoHeep1, iGsfNoHeep2);
        if (gsfGsfMassNoHeep > massCut) {
          if (fabs(gsfsc_eta[iGsfNoHeep1]) < 2.5 && fabs(gsfsc_eta[iGsfNoHeep2]) < 2.5) {
            double fakeRate1 = FakeRate(gsf_gsfet[iGsfNoHeep1], gsfsc_eta[iGsfNoHeep1]);
            double fakeRate2 = FakeRate(gsf_gsfet[iGsfNoHeep2], gsfsc_eta[iGsfNoHeep2]);
            if (fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep2]) > 1.56) {
              histoGsfGsfMassNoHeep.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
              if (p >= DY20 && p <= DY800) {
                if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                  histoSectionDYCombinedGGNH.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
              }
              if (gsfsc_eta[iGsfNoHeep1] * gsfsc_eta[iGsfNoHeep2] > 0) {
                histoGsfGsfMassNoHeep.at(EES)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedGGNH.at(EES)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
                }
              } else {
                histoGsfGsfMassNoHeep.at(EEO)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedGGNH.at(EEO)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
                }
              }
              if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1]) && PassFRPreSel(iGsfNoHeep2, gsfsc_eta[iGsfNoHeep2])) {
                histoGsfGsfMassNoHeepFR.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedGGNHFR.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                }
              }
            } else if ((fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep2]) > 1.56) || (fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep2]) < 1.442)) {
              histoGsfGsfMassNoHeep.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
              if (p >= DY20 && p <= DY800) {
                if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                  histoSectionDYCombinedGGNH.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
              }
              if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1]) && PassFRPreSel(iGsfNoHeep2, gsfsc_eta[iGsfNoHeep2])) {
                histoGsfGsfMassNoHeepFR.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedGGNHFR.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                }
              }
            } else if (fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep2]) < 1.442) {
              histoGsfGsfMassNoHeep.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
              if (p >= DY20 && p <= DY800) {
                if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                  histoSectionDYCombinedGGNH.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi);
              }
              if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1]) && PassFRPreSel(iGsfNoHeep2, gsfsc_eta[iGsfNoHeep2])) {
                histoGsfGsfMassNoHeepFR.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedGGNHFR.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                }
              }
            }
          }
        }
      }

      if (iHeep1 >= 0) {
        // fill the HEEP-GSF and GSF-HEEP cases
        int iGsf = (iHeep1 == iGsf1) ? iGsf2 : iGsf1;
        float heepGsfMass = CalcInvariantMass(iHeep1, iGsf);
        if (heepGsfMass > massCut) {
          if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iGsf]) < 2.5) {
            double fakeRate = FakeRate(gsf_gsfet[iGsf], gsfsc_eta[iGsf]);
            if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsf]) > 1.56) {
              histoHeepGsfMass.at(EE)->Fill(heepGsfMass, input.at(p).second * lumi);
              if (gsfsc_eta[iHeep1] * gsfsc_eta[iGsf] > 0) {
                histoHeepGsfMass.at(EES)->Fill(heepGsfMass, input.at(p).second * lumi);
              } else {
                histoHeepGsfMass.at(EEO)->Fill(heepGsfMass, input.at(p).second * lumi);
              }
              if (PassFRPreSel(iGsf, gsfsc_eta[iGsf])) {
                histoHeepGsfMassFR.at(EE)->Fill(heepGsfMass, input.at(p).second * lumi * fakeRate);
              }
            } else if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsf]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsf]) < 1.442)) {
              histoHeepGsfMass.at(BE)->Fill(heepGsfMass, input.at(p).second * lumi);
              if (PassFRPreSel(iGsf, gsfsc_eta[iGsf])) {
                histoHeepGsfMassFR.at(BE)->Fill(heepGsfMass, input.at(p).second * lumi * fakeRate);
              }
            } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsf]) < 1.442) {
              histoHeepGsfMass.at(BB)->Fill(heepGsfMass, input.at(p).second * lumi);
              if (PassFRPreSel(iGsf, gsfsc_eta[iGsf])) {
                histoHeepGsfMassFR.at(BB)->Fill(heepGsfMass, input.at(p).second * lumi * fakeRate);
              }
            }
          }
        }

        // fill the HEEP-GSF and GSF-HEEP cases where the GSF does not pass the HEEP selection
        if (iGsfNoHeep1 >= 0 && iHeep2 < 0) {
          float heepGsfMassNoHeep = CalcInvariantMass(iHeep1, iGsfNoHeep1);
          if (heepGsfMassNoHeep > massCut) {
            if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iGsfNoHeep1]) < 2.5) {
              double fakeRate = FakeRate(gsf_gsfet[iGsfNoHeep1], gsfsc_eta[iGsfNoHeep1]);
              if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56) {
                histoHeepGsfMassNoHeep.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedHGNH.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                }
                if (gsfsc_eta[iHeep1] * gsfsc_eta[iGsfNoHeep1] > 0) {
                  histoHeepGsfMassNoHeep.at(EES)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                  if (p >= DY20 && p <= DY800) {
                    if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                      histoSectionDYCombinedHGNH.at(EES)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                  }
                } else {
                  histoHeepGsfMassNoHeep.at(EEO)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                  if (p >= DY20 && p <= DY800) {
                    if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                      histoSectionDYCombinedHGNH.at(EEO)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                  }
                }
                if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1])) {
                  histoHeepGsfMassNoHeepFR.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                  if (p >= DY20 && p <= DY800) {
                    if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                      histoSectionDYCombinedHGNHFR.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                  }
                }
              } else if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442)) {
                histoHeepGsfMassNoHeep.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedHGNH.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                }
                if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1])) {
                  histoHeepGsfMassNoHeepFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                  if (p >= DY20 && p <= DY800) {
                    if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                      histoSectionDYCombinedHGNHFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                  }
                }
                // special cuts for Laurent
                histoHeepGsfMassNoHeepNoHOverEFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                if (gsf_trackiso[iHeep1] < 1 && gsf_hovere[iHeep1] < 0.02 && fabs(gsf_deltaphi[iHeep1]) < 0.01 && calomet < 40) {
                  histoTAndPHeepGsfMassNoHeepFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                }
              } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442) {
                histoHeepGsfMassNoHeep.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedHGNH.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi);
                }
                if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1])) {
                  histoHeepGsfMassNoHeepFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                  if (p >= DY20 && p <= DY800) {
                    if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                      histoSectionDYCombinedHGNHFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                  }
                }
                // special cuts for Laurent
                histoHeepGsfMassNoHeepNoHOverEFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                if (gsf_trackiso[iHeep1] < 1 && gsf_hovere[iHeep1] < 0.02 && fabs(gsf_deltaphi[iHeep1]) < 0.01 && calomet < 40) {
                  histoTAndPHeepGsfMassNoHeepFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * fakeRate / (1 - fakeRate));
                }
              }
            }
          }
        } 
      }

      if (iHeep1 >= 0 && iHeep2 >= 0) {
        // fill the HEEP-HEEP cases
        float heepHeepMass = CalcInvariantMass(iHeep1, iHeep2);
        if (heepHeepMass > massCut) {
          if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iHeep2]) < 2.5) {
            if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iHeep2]) > 1.56) {
              histoHeepHeepMass.at(EE)->Fill(heepHeepMass, input.at(p).second * lumi);
              // fill the correct sector of the DY sample
              if (p >= DY20 && p <= DY800) {
                if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                  histoSectionDYCombinedHH.at(EE)->Fill(heepHeepMass, input.at(p).second * lumi);
              }
              if (gsfsc_eta[iHeep1] * gsfsc_eta[iHeep2] > 0) {
                histoHeepHeepMass.at(EES)->Fill(heepHeepMass, input.at(p).second * lumi);
                // fill the correct sector of the DY sample
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedHH.at(EES)->Fill(heepHeepMass, input.at(p).second * lumi);
                }
              } else {
                histoHeepHeepMass.at(EEO)->Fill(heepHeepMass, input.at(p).second * lumi);
                // fill the correct sector of the DY sample
                if (p >= DY20 && p <= DY800) {
                  if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                    histoSectionDYCombinedHH.at(EEO)->Fill(heepHeepMass, input.at(p).second * lumi);
                }
              }
            } else if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iHeep2]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iHeep2]) < 1.442)) {
              histoHeepHeepMass.at(BE)->Fill(heepHeepMass, input.at(p).second * lumi);
              // fill the correct sector of the DY sample
              if (p >= DY20 && p <= DY800) {
                if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                  histoSectionDYCombinedHH.at(BE)->Fill(heepHeepMass, input.at(p).second * lumi);
              }
            } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iHeep2]) < 1.442) {
              histoHeepHeepMass.at(BB)->Fill(heepHeepMass, input.at(p).second * lumi);
              // fill the correct sector of the DY sample
              if (p >= DY20 && p <= DY800) {
                if (genboson_m_branch >= minPtDY.at(p - DY20) && genboson_m_branch < minPtDY.at(p))
                  histoSectionDYCombinedHH.at(BB)->Fill(heepHeepMass, input.at(p).second * lumi);
              }
            }
          }
        }
      }
    } // end of loop over events

    outFile->cd();
    outFile->mkdir(folders.back().c_str());
    outFile->cd(folders.back().c_str());

    (*histoHeepHeepMass.at(BBBE)) = (*histoHeepHeepMass.at(BB)) + (*histoHeepHeepMass.at(BE));
    histoHeepHeepMass.at(BBBE)->SetName("histoHeepHeepMass");
    histoHeepHeepMass.at(BBBE)->SetTitle("Dielectron invariant mass");
    (*histoGsfGsfMass.at(BBBE)) = (*histoGsfGsfMass.at(BB)) + (*histoGsfGsfMass.at(BE));
    histoGsfGsfMass.at(BBBE)->SetName("histoGsfGsfMass");
    histoGsfGsfMass.at(BBBE)->SetTitle("GSF invariant mass");
    (*histoHeepGsfMass.at(BBBE)) = (*histoHeepGsfMass.at(BB)) + (*histoHeepGsfMass.at(BE));
    histoHeepGsfMass.at(BBBE)->SetName("histoHeepGsfMass");
    histoHeepGsfMass.at(BBBE)->SetTitle("Heep-GSF invariant mass");
    (*histoGsfGsfMassFR.at(BBBE)) = (*histoGsfGsfMassFR.at(BB)) + (*histoGsfGsfMassFR.at(BE));
    histoGsfGsfMassFR.at(BBBE)->SetName("histoGsfGsfMassFR");
    histoGsfGsfMassFR.at(BBBE)->SetTitle("GSF invariant mass with fake rate");
    (*histoHeepGsfMassFR.at(BBBE)) = (*histoHeepGsfMassFR.at(BB)) + (*histoHeepGsfMassFR.at(BE));
    histoHeepGsfMassFR.at(BBBE)->SetName("histoHeepGsfMassFR");
    histoHeepGsfMassFR.at(BBBE)->SetTitle("HEEP-GSF invariant mass with fake rate");
  
    (*histoGsfGsfMassNoHeep.at(BBBE)) = (*histoGsfGsfMassNoHeep.at(BB)) + (*histoGsfGsfMassNoHeep.at(BE));
    histoGsfGsfMassNoHeep.at(BBBE)->SetName("histoGsfGsfMassNoHeep");
    histoGsfGsfMassNoHeep.at(BBBE)->SetTitle("GSF(non HEEP) invariant mass");
    (*histoHeepGsfMassNoHeep.at(BBBE)) = (*histoHeepGsfMassNoHeep.at(BB)) + (*histoHeepGsfMassNoHeep.at(BE));
    histoHeepGsfMassNoHeep.at(BBBE)->SetName("histoHeepGsfMassNoHeep");
    histoHeepGsfMassNoHeep.at(BBBE)->SetTitle("GSF(non HEEP) invariant mass");
    (*histoGsfGsfMassNoHeepFR.at(BBBE)) = (*histoGsfGsfMassNoHeepFR.at(BB)) + (*histoGsfGsfMassNoHeepFR.at(BE));
    histoGsfGsfMassNoHeepFR.at(BBBE)->SetName("histoGsfGsfMassNoHeepFR");
    histoGsfGsfMassNoHeepFR.at(BBBE)->SetTitle("GSF(non HEEP) invariant mass with fake rate");
    (*histoHeepGsfMassNoHeepFR.at(BBBE)) = (*histoHeepGsfMassNoHeepFR.at(BB)) + (*histoHeepGsfMassNoHeepFR.at(BE));
    histoHeepGsfMassNoHeepFR.at(BBBE)->SetName("histoHeepGsfMassNoHeepFR");
    histoHeepGsfMassNoHeepFR.at(BBBE)->SetTitle("GSF(non HEEP) invariant mass with fake rate");
  
    if (p >= DY20 && p <= DY800) {
      (*histoSectionDYCombinedHH.at(BBBE)) = (*histoSectionDYCombinedHH.at(BB)) + (*histoSectionDYCombinedHH.at(BE));
      histoSectionDYCombinedHH.at(BBBE)->SetName("histoSectionDYCombinedHH");
      histoSectionDYCombinedHH.at(BBBE)->SetTitle("Section for DY combined");
      (*histoSectionDYCombinedHGNH.at(BBBE)) = (*histoSectionDYCombinedHGNH.at(BB)) + (*histoSectionDYCombinedHGNH.at(BE));
      histoSectionDYCombinedHGNH.at(BBBE)->SetName("histoSectionDYCombinedHGNH");
      histoSectionDYCombinedHGNH.at(BBBE)->SetTitle("Section for DY combined HEEP-GSF non HEEP");
      (*histoSectionDYCombinedGGNH.at(BBBE)) = (*histoSectionDYCombinedGGNH.at(BB)) + (*histoSectionDYCombinedGGNH.at(BE));
      histoSectionDYCombinedGGNH.at(BBBE)->SetName("histoSectionDYCombinedGGNH");
      histoSectionDYCombinedGGNH.at(BBBE)->SetTitle("Section for DY combined GSF-GSF non HEEP");
      (*histoSectionDYCombinedHGNHFR.at(BBBE)) = (*histoSectionDYCombinedHGNHFR.at(BB)) + (*histoSectionDYCombinedHGNHFR.at(BE));
      histoSectionDYCombinedHGNHFR.at(BBBE)->SetName("histoSectionDYCombinedHGNHFR");
      histoSectionDYCombinedHGNHFR.at(BBBE)->SetTitle("Section for DY combined GSF-GSF non HEEP with fake rate");
      (*histoSectionDYCombinedGGNHFR.at(BBBE)) = (*histoSectionDYCombinedGGNHFR.at(BB)) + (*histoSectionDYCombinedGGNHFR.at(BE));
      histoSectionDYCombinedGGNHFR.at(BBBE)->SetName("histoSectionDYCombinedGGNHFR");
      histoSectionDYCombinedGGNHFR.at(BBBE)->SetTitle("Section for DY combined GSF-GSF non HEEP with fake rate");
    }

    // ratio between GSF-GSF and HEEP-GSF fake rates 
    vector<TH1F *> histoRatioMassFR;
    for (unsigned int i = 0; i < histoGsfGsfMassFR.size(); ++i) {
      TH1F *numHisto = (TH1F *)histoGsfGsfMassFR.at(i)->Clone("numHistoUncorr" + acroSuffix.at(i));
      TH1F *denomHisto = (TH1F *)histoHeepGsfMassFR.at(i)->Clone("denomHistoUncorr" + acroSuffix.at(i));
      numHisto->Rebin(ratioRebin);
      denomHisto->Rebin(ratioRebin);
      TH1F *ratioHisto = (TH1F *)numHisto->Clone("ratioHistoUncorr" + acroSuffix.at(i)); // just to have a histogram with the same binning to write the ratio to
      ratioHisto->Divide(numHisto, denomHisto);
      histoRatioMassFR.push_back((TH1F *)ratioHisto->Clone("histoRatioMassFR" + acroSuffix.at(i)));
      histoRatioMassFR.back()->SetTitle("Ratio between GSF-GSF and HEEP-GSF" + suffix.at(i) + ";M_{ee} (GeV/c^{2});GSF-GSF / HEEP-GSF");
      histoRatioMassFR.back()->Write();
    }
    vector<TH1F *> histoRatioMassNoHeepFR;
    for (unsigned int i = 0; i < histoGsfGsfMassNoHeepFR.size(); ++i) {
      TH1F *numHisto = (TH1F *)histoGsfGsfMassNoHeepFR.at(i)->Clone("numHistoNoHeep" + acroSuffix.at(i));
      TH1F *denomHisto = (TH1F *)histoHeepGsfMassNoHeepFR.at(i)->Clone("denomHistoNoHeep" + acroSuffix.at(i));
      numHisto->Rebin(ratioRebin);
      denomHisto->Rebin(ratioRebin);
      TH1F *ratioHisto = (TH1F *)numHisto->Clone("ratioHistoNoHeep" + acroSuffix.at(i)); // just to have a histogram with the same binning to write the ratio to
      ratioHisto->Divide(numHisto, denomHisto);
      histoRatioMassNoHeepFR.push_back((TH1F *)ratioHisto->Clone("histoRatioMassNoHeepFR" + acroSuffix.at(i)));
      histoRatioMassNoHeepFR.back()->SetTitle("Ratio between GSF-GSF non HEEP and HEEP-GSF non HEEP" + suffix.at(i) + ";M_{ee} (GeV/c^{2});GSF-GSF / HEEP-GSF");
      histoRatioMassNoHeepFR.back()->Write();
    }

    (*histoPZ.at(0)) = (*histoPZ.at(1)) + (*histoPZ.at(2));
    histoPZ.at(0)->SetName("histoPZ");
    histoPZ.at(0)->SetTitle("Pz of the initial particle");

    //cout << "Nb of HEEP HEEP M>50 " << histoHeepHeepMass.at(BBBE)->GetEntries() << endl; 
    //cout << "Nb of HEEP GSF M>50 " << histoHeepGsfMass.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of GSF GSF M>50 " << histoGsfGsfMass.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of HEEP GSF with FR M>50 " << histoHeepGsfMassFR.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of GSF GSF with FR M>50 " << histoGsfGsfMassFR.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of HEEP GSF(non HEEP) M>50 " << histoHeepGsfMassNoHeep.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of GSF GSF non HEEP M>50 " << histoGsfGsfMassNoHeep.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of HEEP GSF(non HEEP) with FR M>50 " << histoHeepGsfMassNoHeepFR.at(BBBE)->GetEntries() << endl;
    //cout << "Nb of GSF GSF non HEEP with FR M>50 " << histoGsfGsfMassNoHeepFR.at(BBBE)->GetEntries() << endl;

    histosGsfGsfMassFR.push_back(histoGsfGsfMassFR);
    histosGsfGsfMassNoHeepFR.push_back(histoGsfGsfMassNoHeepFR);
    histosHeepGsfMassNoHeepFR.push_back(histoHeepGsfMassNoHeepFR);
    if (p >= DY20 && p <= DY800) {
      histosSectionDYCombinedHH.push_back(histoSectionDYCombinedHH);
      //histosSectionDYCombinedHGNH.push_back(histoSectionDYCombinedHGNH);
      //histosSectionDYCombinedGGNH.push_back(histoSectionDYCombinedGGNH);
      histosSectionDYCombinedHGNHFR.push_back(histoSectionDYCombinedHGNHFR);
      //histosSectionDYCombinedGGNHFR.push_back(histoSectionDYCombinedGGNHFR);
    }

    for (unsigned int i = BBBE; i <= EEO; ++i) {
      histoHeepHeepMass.at(i)->Write();
      histoGsfGsfMass.at(i)->Write();
      histoGsfGsfMassNoHeep.at(i)->Write();
      histoHeepGsfMass.at(i)->Write();
      histoHeepGsfMassNoHeep.at(i)->Write();
      if (p >= DY20 && p <= DY800) {
        histoSectionDYCombinedHH.at(i)->Write();
        histoSectionDYCombinedHGNH.at(i)->Write();
        histoSectionDYCombinedGGNH.at(i)->Write();
      }
      if (i > EE) continue;
      histoGsfGsfMassFR.at(i)->Write();
      histoGsfGsfMassNoHeepFR.at(i)->Write();
      histoHeepGsfMassFR.at(i)->Write();
      histoHeepGsfMassNoHeepFR.at(i)->Write();
      if (p >= DY20 && p <= DY800) {
        histoSectionDYCombinedHGNHFR.at(i)->Write();
        histoSectionDYCombinedGGNHFR.at(i)->Write();
      }
      if (i > BE) continue;
      histoHeepGsfMassNoHeepNoHOverEFR.at(i)->Write();
      histoTAndPHeepGsfMassNoHeepFR.at(i)->Write();
    }
    for (vector<TH1F *>::iterator it = histoPZ.begin(); it < histoPZ.end(); ++it)
      (*it)->Write();
  } // end of loop over input files

  outFile->cd();
  outFile->mkdir("combinations");  // WORKAROUND - should be in cd() but there seems to be a bug in root 4.27 that gives "R__unzip: error in header" in the rescaling process when histograms are in the base directory
  outFile->cd("combinations");

  vector<TH1F *> histoDYCombined;
  sStream << "# of events / pb^{-1} " << ((massMax - massMin) / nBins) << "GeV/c^{2}"; 
  for (unsigned int i = BBBE; i <= EEO; ++i) {
    histoDYCombined.push_back(new TH1F("histoDYCombined" + acroSuffix.at(i), "Combined Mass for DY samples" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
    histoDYCombined.back()->GetYaxis()->SetTitle(sStream.str().c_str());
    for (vector<vector<TH1F *> >::iterator it = histosSectionDYCombinedHH.begin(); it < histosSectionDYCombinedHH.end(); ++it) {
      it->at(i)->GetYaxis()->SetTitle(sStream.str().c_str());
      histoDYCombined.back()->Add(it->at(i));
    }
    histoDYCombined.back()->Write();
  }
  sStream.str("");

  // combined DY and gamma MC HEEP GSF fake rate
  vector<TH1F *> histosHGNHFRDYCombined;
  vector<TH1F *> histosHGNHFRGammaCombined;
  sStream << "# of events * fake rate/ pb^{-1} " << ((massMax - massMin) / nBins) << "GeV/c^{2}"; 
  // recalculate fake rate with corrections
  vector<TH1F *> histosHeepGsfCorr;
  vector<TH1F *> histosGsfGsfCorr;
  for (unsigned int i = BBBE; i <= EE; ++i) {
    histosHGNHFRDYCombined.push_back(new TH1F("histoHGNHFRDYCombined" + acroSuffix.at(i), "Combined HEEP-GSF(no HEEP) with fake rate from DY samples" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray));
    histosHGNHFRDYCombined.back()->GetYaxis()->SetTitle(sStream.str().c_str());
    for (vector<vector<TH1F *> >::iterator it = histosSectionDYCombinedHGNHFR.begin(); it < histosSectionDYCombinedHGNHFR.end(); ++it) {
      it->at(i)->GetYaxis()->SetTitle(sStream.str().c_str());
      histosHGNHFRDYCombined.back()->Add(it->at(i));
    }

    histosHeepGsfCorr.push_back(new TH1F((*histosHeepGsfMassNoHeepFR.at(DATA).at(i)) - (*histosHGNHFRDYCombined.back())));
    histosHeepGsfCorr.back()->SetName("histoHeepGsfCorr" + acroSuffix.at(i));
    histosHeepGsfCorr.back()->SetTitle("Corrected HEEP-GSF(no HEEP) fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})");

    //generate combined gamma HEEP-GSF non HEEP fake rate histogram
    TH1F *histoHGNHFRGammaCombined = new TH1F("histoHGNHFRGammaCombined" + acroSuffix.at(i), "Combined HEEP-GSF(no HEEP) with fake rate from gamma samples" + suffix.at(i) + ";M_{ee} (GeV/c^{2})", nBins, binArray);
    for (unsigned int j = G30T50; j <= G470T800; ++j) {
      histoHGNHFRGammaCombined->Add(histosHeepGsfMassNoHeepFR.at(j).at(i));
    }
    histosHGNHFRGammaCombined.push_back(new TH1F(*histoHGNHFRGammaCombined));

    histosGsfGsfCorr.push_back(new TH1F((*histosGsfGsfMassNoHeepFR.at(DATA).at(i)) + (*histosHeepGsfMassNoHeepFR.at(WJETS).at(i)) + (*histoHGNHFRGammaCombined)));
    histosGsfGsfCorr.back()->SetName("histoGsfGsfCorr" + acroSuffix.at(i));
    histosGsfGsfCorr.back()->SetTitle("Corrected GSF-GSF(no HEEP) fake rate" + suffix.at(i) + ";M_{ee} (GeV/c^{2})");

    histosHGNHFRDYCombined.back()->Write();
    histosHGNHFRGammaCombined.back()->Write();
    histosHeepGsfCorr.back()->Write();
    histosGsfGsfCorr.back()->Write();
  }
  sStream.str("");

  vector<TH1F *> histoRatioMassCorrFR;
  for (unsigned int i = 0; i < histosGsfGsfCorr.size(); ++i) {
    TH1F *numHisto = (TH1F *)histosGsfGsfCorr.at(i)->Clone("numHistoCorr" + acroSuffix.at(i));
    TH1F *denomHisto = (TH1F *)histosHeepGsfCorr.at(i)->Clone("denomHistoCorr" + acroSuffix.at(i));
    numHisto->Rebin(ratioRebin);
    denomHisto->Rebin(ratioRebin);
    TH1F *ratioHisto = (TH1F *)numHisto->Clone("ratioHistoCorr" + acroSuffix.at(i)); // just to have a histogram with the same binning to write the ratio to
    ratioHisto->Divide(numHisto, denomHisto);
    histoRatioMassCorrFR.push_back((TH1F *)ratioHisto->Clone("histoRatioMassCorrFR" + acroSuffix.at(i)));

    histoRatioMassCorrFR.back()->SetTitle("Ratio between corrected GSF-GSF and corrected HEEP-GSF" + suffix.at(i) + ";M_{ee} (GeV/c^{2});GSF-GSF / HEEP-GSF");
    histoRatioMassCorrFR.back()->Write();
  }

  outFile->Close();

  if (calcLumiCorrFactor)
    NormalizeToZPeak(normToZPeakFrom, normToZPeakTo, ssOutfile.str().c_str(), folders, ratioRebin);
  else if (lumiCorrFactor != 1.) {
    cout << endl << "-------------------------------------------------------------------------" << endl;
    cout << "Applied correction factor for MC: " << lumiCorrFactor << endl;
    cout << "-------------------------------------------------------------------------" << endl << endl;
  }

} //end of method

int
InvariantMass::NormalizeToZPeak(const float& lowMass, const float& highMass, const char* inputFile, vector<string> folders, int ratioRebin)
{
  float factor = 1.;
  const unsigned int DATA = 0;
  const unsigned int TTBAR = 6;
  const unsigned int WJETS = 7;

  vector<TString> acroSuffix;
  acroSuffix.push_back("");
  acroSuffix.push_back("BB");
  acroSuffix.push_back("BE");
  acroSuffix.push_back("EE");

  vector<TString> noHeepStr;
  noHeepStr.push_back("");
  noHeepStr.push_back("NoHeep");

  if (highMass == lowMass) {
    cout << endl << "Normalization to Z peak range is zero. Will apply no correction factor." << endl << endl;
    return 1;
  }

  folders.push_back("combinations");

  TFile *file = new TFile(inputFile, "update");

  // get the necessary histograms
  file->cd(folders.at(DATA).c_str());
  TH1F *dataHisto = (TH1F *)gDirectory->Get("histoHeepHeepMass");
  file->cd(folders.at(TTBAR).c_str());
  TH1F *ttbarBgHisto = (TH1F *)gDirectory->Get("histoHeepHeepMass");
  file->cd(folders.at(WJETS).c_str());
  TH1F *wjetBgHisto = (TH1F *)gDirectory->Get("histoHeepHeepMass");
  file->cd("combinations");
  TH1F *dyBgHisto = (TH1F *)gDirectory->Get("histoDYCombined");
  TH1F *qcdBgHisto = (TH1F *)gDirectory->Get("histoGsfGsfCorr");

  double intData = 0;
  double intBg = 0;

  // integrate over range
  intData = dataHisto->Integral(dataHisto->FindBin(lowMass), dataHisto->FindBin(highMass));
  intBg = dyBgHisto->Integral(dyBgHisto->FindBin(lowMass), dyBgHisto->FindBin(highMass));
  intBg += ttbarBgHisto->Integral(ttbarBgHisto->FindBin(lowMass), ttbarBgHisto->FindBin(highMass));
  intBg += wjetBgHisto->Integral(wjetBgHisto->FindBin(lowMass), wjetBgHisto->FindBin(highMass));
  intBg += qcdBgHisto->Integral(qcdBgHisto->FindBin(lowMass), qcdBgHisto->FindBin(highMass));

  cout << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "# of events (weighted) in the Z peak from " << lowMass << "GeV to " << highMass << "GeV" << endl;
  cout << " Data     " << intData << endl;
  cout << " DY       " << dyBgHisto->Integral(dyBgHisto->FindBin(lowMass), dyBgHisto->FindBin(highMass)) << endl;
  cout << " ttbar    " << ttbarBgHisto->Integral(ttbarBgHisto->FindBin(lowMass), ttbarBgHisto->FindBin(highMass)) << endl;
  cout << " W+jets   " << wjetBgHisto->Integral(wjetBgHisto->FindBin(lowMass), wjetBgHisto->FindBin(highMass)) << endl;
  cout << " QCD      " << qcdBgHisto->Integral(qcdBgHisto->FindBin(lowMass), qcdBgHisto->FindBin(highMass)) << endl << endl;

  // calculate the noramlization factor
  factor = intData / intBg;

  cout << "Applied correction factor for MC: " << factor << endl;
  cout << "-------------------------------------------------------------------------" << endl << endl;

  // scale histograms in file with normalization factor
  for (unsigned int i = 1; i < folders.size(); ++i) {
    file->cd(folders.at(i).c_str());
    //cout << "directory: /" << folders.at(i) << endl;
    // loop over all TH1F in the directory
    vector<TString> histoNames;
    for (int j = 0; j < gDirectory->GetListOfKeys()->GetSize(); ++j) {
      TKey *key = (TKey *)gDirectory->GetListOfKeys()->At(j);
      if (!key->ReadObj()->IsA()->InheritsFrom("TH1F")) continue;

      TH1F *histo = (TH1F *)key->ReadObj();
      string histoName = histo->GetName();
      if (histoName.find("histoRatioMass") != string::npos || histoName.rfind("GsfCorr") != string::npos) continue;
      histoNames.push_back(histoName.c_str());
      //cout << "histo: " << histoName << endl;
    }
    for (vector<TString>::iterator it = histoNames.begin(); it < histoNames.end(); ++it) {
      TH1F *histo = (TH1F *)gDirectory->Get(it->Data());
      histo->Scale(factor);
      histo->Write(0, TObject::kOverwrite);
    }

    // redo the ratio histograms and corrected histograms
    if (folders.at(i) == "combinations") {
      file->cd("combinations");
      for (unsigned int p = BBBE; p < EES; ++p) {
        file->cd(folders.at(0).c_str());
        TH1F *histoHeepGsfMassNoHeepFR = (TH1F *)gDirectory->Get("histoHeepGsfMassNoHeepFR" + acroSuffix[p]);
        TH1F *histoGsfGsfMassNoHeepFR = (TH1F *)gDirectory->Get("histoGsfGsfMassNoHeepFR" + acroSuffix[p]);
        file->cd(folders.at(WJETS).c_str());
        TH1F *histoWjetsFR = (TH1F *)gDirectory->Get("histoHeepGsfMassNoHeepFR" + acroSuffix[p] + ";1"); // ";1" do get the version from the file not the rebinned one in the memory
        file->cd("combinations");
        TH1F *histoDYCombFR = (TH1F *)gDirectory->Get("histoHGNHFRDYCombined" + acroSuffix[p]);
        TH1F *histoGammaCombFR = (TH1F *)gDirectory->Get("histoHGNHFRGammaCombined" + acroSuffix[p]);

        TH1F *numHisto = (TH1F *)gDirectory->Get("histoGsfGsfCorr" + acroSuffix[p]);
        TH1F *denomHisto = (TH1F *)gDirectory->Get("histoHeepGsfCorr" + acroSuffix[p]);
        TH1F *ratioHisto = (TH1F *)gDirectory->Get("histoRatioMassCorrFR" + acroSuffix[p]);

        TString name = numHisto->GetName();
        TString title = numHisto->GetTitle();
        (*numHisto) = (*histoGsfGsfMassNoHeepFR) + (*histoWjetsFR) + (*histoGammaCombFR);
        numHisto->SetName(name);
        numHisto->SetTitle(title);

        name = denomHisto->GetName();
        title = denomHisto->GetTitle();
        (*denomHisto) = (*histoHeepGsfMassNoHeepFR) - (*histoDYCombFR); 
        denomHisto->SetName(name);
        denomHisto->SetTitle(title);

        numHisto->Write(0, TObject::kOverwrite);
        denomHisto->Write(0, TObject::kOverwrite);

        numHisto->Rebin(ratioRebin);
        denomHisto->Rebin(ratioRebin);

        ratioHisto->Divide(numHisto, denomHisto);
        ratioHisto->Write(0, TObject::kOverwrite);
      } 
    } else {
      file->cd(folders.at(i).c_str());
      for (unsigned int k = 0; k < 2; ++k) {
        for (unsigned int p = BBBE; p < EES; ++p) {
          TH1F *numHisto = (TH1F *)gDirectory->Get("histoGsfGsfMass" + noHeepStr[k] + "FR" + acroSuffix[p]);
          TH1F *denomHisto = (TH1F *)gDirectory->Get("histoHeepGsfMass" + noHeepStr[k] + "FR" + acroSuffix[p]);
          TH1F *ratioHisto = (TH1F *)gDirectory->Get("histoRatioMass" + noHeepStr[k] + "FR" + acroSuffix[p]);

          numHisto->Rebin(ratioRebin);
          denomHisto->Rebin(ratioRebin);

          ratioHisto->Divide(numHisto, denomHisto);
          ratioHisto->Write(0, TObject::kOverwrite);
        }
      }
    }    
  }

  file->Close();
  return 0;
}

double
InvariantMass::CalcInvariantMass (const int& iEle1, const int& iEle2)
{
  const float dataEtShiftFactorEB = 1.0036;
  const float dataEtShiftFactorEE = 1.0256;
  //const float dataEtShiftFactorEB = 1.; // fix this factor
  //const float dataEtShiftFactorEE = 1.; // fix this factor

  TLorentzVector ele1;
  TLorentzVector ele2;

  Float_t et1 = gsf_gsfet[iEle1];
  Float_t et2 = gsf_gsfet[iEle2];
  // correct energy in data
  (fabs(gsfsc_eta[iEle1]) > 1.56 && runnumber > 1) ? et1 *= dataEtShiftFactorEE : et1 *= dataEtShiftFactorEB;
  (fabs(gsfsc_eta[iEle2]) > 1.56 && runnumber > 1) ? et2 *= dataEtShiftFactorEE : et2 *= dataEtShiftFactorEB;

  ele1.SetPtEtaPhiE(et1, gsf_eta[iEle1], gsf_phi[iEle1], (et1 * cosh(gsf_eta[iEle1])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[iEle2], gsf_phi[iEle2], (et2 * cosh(gsf_eta[iEle2])));
  
  return (ele1+ele2).Mag();
}

double
InvariantMass::CalcPz (const int& iEle1, const int& iEle2)
{
  const float dataEtShiftFactorEB = 1.0036;
  const float dataEtShiftFactorEE = 1.0256;
  //const float dataEtShiftFactorEB = 1.; // fix this factor
  //const float dataEtShiftFactorEE = 1.; // fix this factor

  TLorentzVector ele1;
  TLorentzVector ele2;

  Float_t et1 = gsf_gsfet[iEle1];
  Float_t et2 = gsf_gsfet[iEle2];
  // correct energy in data
  (fabs(gsfsc_eta[iEle1]) > 1.56 && runnumber > 1) ? et1 *= dataEtShiftFactorEE : et1 *= dataEtShiftFactorEB;
  (fabs(gsfsc_eta[iEle2]) > 1.56 && runnumber > 1) ? et2 *= dataEtShiftFactorEE : et2 *= dataEtShiftFactorEB;

  ele1.SetPtEtaPhiE(et1, gsf_eta[iEle1], gsf_phi[iEle1], (et1 * cosh(gsf_eta[iEle1])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[iEle2], gsf_phi[iEle2], (et2 * cosh(gsf_eta[iEle2])));
  
  return (ele1+ele2).Pz();
}

double
InvariantMass::FakeRate (const float& et, const float& eta)
{
  // constants 13/07/2011
  vector<pair<double, double> > fakeRateFitLowPtEB;
  fakeRateFitLowPtEB.push_back(make_pair(0.024, 0.00315));
  fakeRateFitLowPtEB.push_back(make_pair(-0.0001, 0.0000581));
  const pair<double, double> fakeRateFitHighPtEB(0.009, 0.);

  const pair<double, double> fakeRateFitLowPtEELow(0.081, -0.0006);
  const pair<double, double> fakeRateFitHighPtEELow(0.036, 0.);

  const pair<double, double> fakeRateFitLowPtEEHigh(0.114, -0.0006);
  const pair<double, double> fakeRateFitHighPtEEHigh(0.067, 0.);

  const float fitSelectPtEB = 150;
  const float fitSelectPtEELow = 75;
  const float fitSelectPtEEHigh = 78.3;
  const float etaDivideEE = 2.0;

  double fakeRate = 0;

  if (fabs(eta) < 2.5) {
    if (fabs(eta) > etaDivideEE) {
      // EE high eta
      if (et > fitSelectPtEEHigh)
        fakeRate = fakeRateFitHighPtEEHigh.first + fakeRateFitHighPtEEHigh.second * et;
      else
        fakeRate = fakeRateFitLowPtEEHigh.first + fakeRateFitLowPtEEHigh.second * et;
    }
    else if (fabs(eta) > 1.56) {
      // EE low eta
      if (et > fitSelectPtEELow)
        fakeRate = fakeRateFitHighPtEELow.first + fakeRateFitHighPtEELow.second * et;
      else
        fakeRate = fakeRateFitLowPtEELow.first + fakeRateFitLowPtEELow.second * et;
    }
    else if (fabs(eta) < 1.442) {
      // EB
      if (et > fitSelectPtEB)
        fakeRate = fakeRateFitHighPtEB.first + fakeRateFitHighPtEB.second * et;
      else {
        for (vector<pair<double, double> >::iterator it = fakeRateFitLowPtEB.begin(); it < fakeRateFitLowPtEB.end(); ++it)
          fakeRate += it->first * pow(et, distance(fakeRateFitLowPtEB.begin(), it));
      }
    }
  }

  return fakeRate;
}

bool
InvariantMass::PassFRPreSel(const int &n, const float &eta)
{
  // fake rate preselection
  const float hOverECutEB = 0.15;
  const float hOverECutEE = 0.1;
  const float sieieCutEB = 0.013;
  const float sieieCutEE = 0.034;
  const int missHitsCutEB = 0;
  const int missHitsCutEE = 0;

  if (fabs(eta) < 1.442) {
    if (gsf_hovere[n] > hOverECutEB) return false;
    if (gsf_nLostInnerHits[n] > missHitsCutEB) return false;
    if (gsf_sigmaIetaIeta[n] > sieieCutEB) return false;
  }
  else if (fabs(eta) > 1.56) {
    if (gsf_hovere[n] > hOverECutEE) return false;
    if (gsf_nLostInnerHits[n] > missHitsCutEE) return false;
    if (gsf_sigmaIetaIeta[n] > sieieCutEE) return false;
  }

  return true;
}

