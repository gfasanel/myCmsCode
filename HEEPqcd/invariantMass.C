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
  const float massCut = 0.;
  const float lumiCorrFactor = 1.;   // not used if calcLumiCorrFactor=true
  bool calcLumiCorrFactor = true;    // first run over events once to
                                     // calculate correction and the second 
                                     // time with the factor applied
  float normToZPeakFrom = 60;        // lower mass for normalization
  //float normToZPeakTo = 60;         // upper mass for normalization
  float normToZPeakTo = 120;         // upper mass for normalization
  bool etShiftData = false;           // match data Z peak with MC peak

  string fileNamePrefix = "invMassHistos";   // lumi and .root well be added
  //string fileNamePrefix = "test";

  const int massMin = 50;            // minimum Mass
  const int massMax = 2050;          // maximum Mass

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

  bool usePUInfo = true;
  TString puFile = "file:////user/treis/data2012/pileup/pileup_runA+B+C-ReReco+D-Prompt_puJSON-190389-208686_Photon+DoublePhotonHighPt.root";
  float bar_et = 35.;
  float end_et = 35.;
  
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
  sStream << "# of events / " << ((massMax - massMin) / nBins) << "GeV";

  /////////////////////////////////////////////////////////////////////////
  // input files
  /////////////////////////////////////////////////////////////////////////
  vector<pair<TFile *, double> > input;
  vector<TString> inFileTag;
  vector<float> minPtDY;
  vector<float> minPtG;
  float lumi = 19619.;
  TFile *inData = TFile::Open("file:////user/treis/data2012/Photon_Run2012A+DoublePhotonHighPt_Run2012B+C+D_13Jul2012+06Aug2012+24Aug2012+11Dec2012+PromptReco-Cv2+Dv1_Cert_190456-208686_gct1_45+46_19619pb-1.root");
  input.push_back(make_pair(inData, 1 / lumi));
  inFileTag.push_back("Data");
  const unsigned int DATA = 0;

  TFile *inDY20 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3297045ev.root");
  input.push_back(make_pair(inDY20, 1915. / 3297045.));
  inFileTag.push_back("DY20");
  minPtDY.push_back(20);
  const unsigned int DY20 = 1;

  TFile *inDY120 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-120_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99987ev.root");
  input.push_back(make_pair(inDY120, 11.89 / 99987. * 1915./1871.));
  inFileTag.push_back("DY120");
  minPtDY.push_back(120);
  const unsigned int DY120 = 2;

  TFile *inDY200 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99991ev.root");
  input.push_back(make_pair(inDY200, 1.483 / 99991. * 1915./1871.));
  inFileTag.push_back("DY200");
  minPtDY.push_back(200);
  const unsigned int DY200 = 3;

  TFile *inDY400 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-400_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99991ev.root");
  input.push_back(make_pair(inDY400, 0.1085 / 99991. * 1915./1871.));
  inFileTag.push_back("DY400");
  minPtDY.push_back(400);
  const unsigned int DY400 = 4;

  TFile *inDY500 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99986ev.root");
  input.push_back(make_pair(inDY500, 0.04409 / 99986. * 1915./1871.));
  inFileTag.push_back("DY500");
  minPtDY.push_back(500);
  const unsigned int DY500 = 5;

  TFile *inDY700 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-700_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99990ev.root");
  input.push_back(make_pair(inDY700, 0.01025 / 99990. * 1915./1871.));
  inFileTag.push_back("DY700");
  minPtDY.push_back(700);
  const unsigned int DY700 = 6;

  TFile *inDY800 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99990ev.root");
  input.push_back(make_pair(inDY800, 0.005491 / 99990. * 1915./1871.));
  inFileTag.push_back("DY800");
  minPtDY.push_back(800);
  const unsigned int DY800 = 7;

  TFile *inDY1000 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99992ev.root");
  input.push_back(make_pair(inDY1000, 0.001796 / 99992. * 1915./1871.));
  inFileTag.push_back("DY1000");
  minPtDY.push_back(1000);
  const unsigned int DY1000 = 8;

  TFile *inDY1500 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99999ev.root");
  input.push_back(make_pair(inDY1500, 1.705E-4 / 99999. * 1915./1871.));
  inFileTag.push_back("DY1500");
  minPtDY.push_back(1500);
  const unsigned int DY1500 = 9;

  TFile *inDY2000 = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99993ev.root");
  input.push_back(make_pair(inDY2000, 2.208E-5 / 99993. * 1915./1871.));
  inFileTag.push_back("DY2000");
  minPtDY.push_back(2000);
  const unsigned int DY2000 = 10;

  TFile *inTTBar = TFile::Open("file:////user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
  //input.push_back(make_pair(inTTBar, 225.197 / 28150723.));
  input.push_back(make_pair(inTTBar, 234 / 28150723.));
  inFileTag.push_back("ttbar");
  const unsigned int TTBAR = 11;

  TFile *inDYToTT = TFile::Open("file:////user/treis/mcsamples/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3295238ev.root");
  input.push_back(make_pair(inDYToTT, 1915.1 / 3295238.));
  inFileTag.push_back("DYToTauTau");
  const unsigned int DYTT = 12;

  TFile *inWW = TFile::Open("file:////user/treis/mcsamples/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_10000431ev.root");
  input.push_back(make_pair(inWW, 54.838 / 10000431.));
  inFileTag.push_back("WW");
  const unsigned int WW = 13;

  TFile *inWZ = TFile::Open("file:////user/treis/mcsamples/WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_10000283ev.root");
  input.push_back(make_pair(inWZ, 33.21 / 10000283.));
  inFileTag.push_back("WZ");
  const unsigned int WZ = 14;

  TFile *inZZ = TFile::Open("file:////user/treis/mcsamples/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_9799908ev.root");
  input.push_back(make_pair(inZZ, 17.654 / 9799908.));
  inFileTag.push_back("ZZ");
  const unsigned int ZZ = 15;

  TFile *inTW = TFile::Open("file:////user/treis/mcsamples/T+Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_991118ev.root");
  input.push_back(make_pair(inTW, 22.2 / 991118.));
  inFileTag.push_back("tW");
  const unsigned int TW = 16;

  TFile *inWJets = TFile::Open("file:////user/treis/mcsamples/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_76102995ev.root");
  input.push_back(make_pair(inWJets, 36257.2 / 76102995.));
  inFileTag.push_back("W+jets");
  const unsigned int WJETS = 17;
  
  TFile *inG30T50 = TFile::Open("file:////user/treis/mcsamples/G_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1985352ev.root");
  input.push_back(make_pair(inG30T50, 19931.62 / 1985352. * 1.3));
  inFileTag.push_back("G30to50");
  minPtG.push_back(30);
  const unsigned int G30T50 = 18;
  
  TFile *inG50T80 = TFile::Open("file:////user/treis/mcsamples/G_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1995062ev.root");
  input.push_back(make_pair(inG50T80, 3322.309 / 1995062. * 1.3));
  inFileTag.push_back("G50to80");
  minPtG.push_back(50);
  const unsigned int G50T80 = 19;
  
  TFile *inG80T120 = TFile::Open("file:////user/treis/mcsamples/G_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1976687ev.root");
  input.push_back(make_pair(inG80T120, 558.2865 / 1976687. * 1.3));
  inFileTag.push_back("G80to120");
  minPtG.push_back(80);
  const unsigned int G80T120 = 20;
  
  TFile *inG120T170 = TFile::Open("file:////user/treis/mcsamples/G_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_2000043ev.root");
  input.push_back(make_pair(inG120T170, 108.0068 / 2000043. * 1.3));
  inFileTag.push_back("G120to170");
  minPtG.push_back(120);
  const unsigned int G120T170 = 21;
  
  TFile *inG170T300 = TFile::Open("file:////user/treis/mcsamples/G_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_2000069ev.root");
  input.push_back(make_pair(inG170T300, 30.12207 / 2000069. * 1.3));
  inFileTag.push_back("G170to300");
  minPtG.push_back(170);
  const unsigned int G170T300 = 22;
  
  TFile *inG300T470 = TFile::Open("file:////user/treis/mcsamples/G_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_2000130ev.root");
  input.push_back(make_pair(inG300T470, 2.138632 / 2000130. * 1.3));
  inFileTag.push_back("G300to470");
  minPtG.push_back(300);
  const unsigned int G300T470 = 23;
  
  TFile *inG470T800 = TFile::Open("file:////user/treis/mcsamples/G_Pt-470to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1975231ev.root");
  input.push_back(make_pair(inG470T800, 0.2119244 / 1975231. * 1.3));
  inFileTag.push_back("G470to800");
  minPtG.push_back(470);
  const unsigned int G470T800 = 24;
  
  TFile *inG800T1400 = TFile::Open("file:////user/treis/mcsamples/G_Pt-800to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1973504ev.root");
  input.push_back(make_pair(inG800T1400, 0.007077847 / 1973504. * 1.3));
  inFileTag.push_back("G800to1400");
  minPtG.push_back(800);
  const unsigned int G800T1400 = 25;
  
  TFile *inG1400T1800 = TFile::Open("file:////user/treis/mcsamples/G_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1984890ev.root");
  input.push_back(make_pair(inG1400T1800, 4.510327E-5 / 1984890. * 1.3));
  inFileTag.push_back("G1400to1800");
  minPtG.push_back(1400);
  const unsigned int G1400T1800 = 26;
  
  TFile *inG1800 = TFile::Open("file:////user/treis/mcsamples/G_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_1939122ev.root");
  input.push_back(make_pair(inG1800, 1.867141E-6 / 1939122. * 1.3));
  inFileTag.push_back("G1800");
  minPtG.push_back(1800);
  const unsigned int G1800 = 27;
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

  // data pileup histogram
  TFile *dataPuInfile = TFile::Open(puFile);
  dataPuInfile->Cd("");
  TH1F *puData = (TH1F *)gDirectory->Get("pileup");
  puData->SetDirectory(0);
  puData->SetName("puData");
  dataPuInfile->Close();
  TH1F *puDataNorm = (TH1F *)puData->DrawNormalized()->Clone("puDataNorm"); 

  stringstream ssGoodHeepFileName;
  ssGoodHeepFileName << "goodHeepEvents" << lumi << "pb-1.root";
  TFile *goodEvFile = new TFile(ssGoodHeepFileName.str().c_str(), "recreate");
  goodEvFile->cd();
  TTree *eleDataTree = new TTree("eleDataTree", "eleDataTree");
  int evtRegion = -1;
  float evtWeight = 1.;
  float hHMass = 0.;
  eleDataTree->Branch("runnr", &runnumber, "runnr/I");
  eleDataTree->Branch("eventnr", &eventnumber, "eventnr/i");
  eleDataTree->Branch("lumiSec", &luminosityBlock, "lumiSec/I");
  eleDataTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
  eleDataTree->Branch("weight", &evtWeight, "weight/F");
  eleDataTree->Branch("mass", &hHMass, "mass/F");

  vector<string> folders; // folder names in output root file

  // percentage of data before run 191718 (no online ECAL laser correction)
  float ratioNoOnlineLaser = 0.;

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

    bool etShift = false;
    if (etShiftData && p == DATA)
      etShift = true;

    // normalize bg
    if (p == 1 && !calcLumiCorrFactor) lumi *= lumiCorrFactor;

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
      histoHeepHeepMass.push_back(new TH1F("histoHeepHeepMass" + acroSuffix.at(i), "Dielectron invariant mass" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoHeepHeepMass.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup GSF-GSF histograms
      histoGsfGsfMass.push_back(new TH1F("histoGsfGsfMass" + acroSuffix.at(i), "GSF invariant mass" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoGsfGsfMass.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup GSF-GSF non HEEP histograms
      histoGsfGsfMassNoHeep.push_back(new TH1F("histoGsfGsfMassNoHeep" + acroSuffix.at(i), "GSF(non HEEP) invariant mass" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoGsfGsfMassNoHeep.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup HEEP-GSF histograms
      histoHeepGsfMass.push_back(new TH1F("histoHeepGsfMass" + acroSuffix.at(i), "Heep-GSF invariant mass" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoHeepGsfMass.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup HEEP-GSF(non HEEP) histograms
      histoHeepGsfMassNoHeep.push_back(new TH1F("histoHeepGsfMassNoHeep" + acroSuffix.at(i), "GSF(non HEEP) invariant mass" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoHeepGsfMassNoHeep.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      // setup histogram sections for combined DY histograms
      if (p >= DY20 && p <= DY2000) {
        histoSectionDYCombinedHH.push_back(new TH1F("histoSectionDYCombinedHH" + acroSuffix.at(i), "Section for DY combined HEEP-HEEP" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
        histoSectionDYCombinedHH.back()->GetYaxis()->SetTitle(sStream.str().c_str());

        histoSectionDYCombinedHGNH.push_back(new TH1F("histoSectionDYCombinedHGNH" + acroSuffix.at(i), "Section for DY combined HEEP-GSF non HEEP" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
        histoSectionDYCombinedHGNH.back()->GetYaxis()->SetTitle(sStream.str().c_str());


        histoSectionDYCombinedGGNH.push_back(new TH1F("histoSectionDYCombinedGGNH" + acroSuffix.at(i), "Section for DY combined GSF-GSF non HEEP" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
        histoSectionDYCombinedGGNH.back()->GetYaxis()->SetTitle(sStream.str().c_str());
      }

      // set up fake rate histograms
      if (i > EE) continue;

      histoGsfGsfMassFR.push_back(new TH1F("histoGsfGsfMassFR" + acroSuffix.at(i), "GSF invariant mass with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoGsfGsfMassFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoGsfGsfMassNoHeepFR.push_back(new TH1F("histoGsfGsfMassNoHeepFR" + acroSuffix.at(i), "GSF(non HEEP) invariant mass with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoGsfGsfMassNoHeepFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoHeepGsfMassFR.push_back(new TH1F("histoHeepGsfMassFR" + acroSuffix.at(i), "Heep-GSF invariant mass with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoHeepGsfMassFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoHeepGsfMassNoHeepFR.push_back(new TH1F("histoHeepGsfMassNoHeepFR" + acroSuffix.at(i), "GSF(non HEEP) invariant mass with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoHeepGsfMassNoHeepFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      if (p >= DY20 && p <= DY2000) {
        histoSectionDYCombinedHGNHFR.push_back(new TH1F("histoSectionDYCombinedHGNHFR" + acroSuffix.at(i), "Section for DY combined GSF-GSF non HEEP with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
        histoSectionDYCombinedHGNHFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

        histoSectionDYCombinedGGNHFR.push_back(new TH1F("histoSectionDYCombinedGGNHFR" + acroSuffix.at(i), "Section for DY combined GSF-GSF non HEEP with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
        histoSectionDYCombinedGGNHFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());
      }

      if (i > BE) continue;

      // setup HEEP-GSF(non HEEP) histograms for tag and probe studies for Laurent
      histoHeepGsfMassNoHeepNoHOverEFR.push_back(new TH1F("histoHeepGsfMassNoHeepNoHOverEFR" + acroSuffix.at(i), "HEEP-GSF (non HEEP) with fake rate" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
      histoHeepGsfMassNoHeepNoHOverEFR.back()->GetYaxis()->SetTitle(sStream.str().c_str());

      histoTAndPHeepGsfMassNoHeepFR.push_back(new TH1F("histoTAndPHeepGsfMassNoHeepFR" + acroSuffix.at(i), "Invariant mass for tag and probe with additional cuts" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
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

    unsigned int heepCounter = 0;
    unsigned int heepCounterNoOnlineLaser = 0;
    //Long64_t nbytes = 0, nb = 0;
    /////////////////////////////////////////////////////////////////////////
    // loop over events
    /////////////////////////////////////////////////////////////////////////
    //for (Long64_t jentry=0; (jentry < 10000 && jentry < nentries); ++jentry) { // for tests
    for (Long64_t jentry=0; jentry < nentries; ++jentry) {
      hHMass = 0.;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      thetree->GetEntry(jentry);
      if (jentry % 50000 == 0 ) cout << "entry " << jentry << endl;

      // use trigger turn on without online laser correction for a percentage of MC events
      bool useSecondTurnOn = false;
      if (p > DATA && jentry < (ratioNoOnlineLaser * nentries)) useSecondTurnOn = true;

      // first correct the energy
      //CorrectEnergy(); 

      // trigger fired?
      int prescale = 0;
      if (p == DATA && Trigger(prescale) < 1) continue;

      // at least two gsf electrons
      if (gsf_size < 2) continue;

      // set the PU weight
      float puWeight = 1.;
      if (usePUInfo && p > DATA) puWeight = puWeights->GetBinContent(puWeights->FindBin(trueNVtx));

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
        if (highestHeepEt < gsf_gsfet[n] && PassHEEP(n)) {
          iHeep1 = n;
          highestHeepEt = gsf_gsfet[n];
        }
        if (highestNoHeepEt < gsf_gsfet[n] && !PassHEEP(n)) {
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
        if (highestHeepEt < gsf_gsfet[m] && PassHEEP(m) && m != iHeep1) {
          iHeep2 = m;
          highestHeepEt = gsf_gsfet[m];
        }
        if (highestNoHeepEt < gsf_gsfet[m] && !PassHEEP(m) && m != iGsfNoHeep1) {
          iGsfNoHeep2 = m;
          highestNoHeepEt = gsf_gsfet[m];
        }
      }

      double trgTurnOnWeight = 1.;

      if (iGsf2 < 0 || iGsf1 < 0) continue; // not enough good GSF electrons found 

      ////////////////////////////////////////////////////////////////////////
      // fill the histograms
      ////////////////////////////////////////////////////////////////////////
      // fill the GSF-GSF cases
      double gsfGsfMass = CalcInvariantMass(iGsf1, iGsf2, etShift);
      if (p > DATA) trgTurnOnWeight = TriggerTurnOn(iGsf1, useSecondTurnOn) * TriggerTurnOn(iGsf2, useSecondTurnOn);
      if (gsfGsfMass > massCut) {
        if (fabs(gsfsc_eta[iGsf1]) < 2.5 && fabs(gsfsc_eta[iGsf2]) < 2.5) {
          double fakeRate1 = FakeRate(gsf_gsfet[iGsf1], gsfsc_eta[iGsf1]);
          double fakeRate2 = FakeRate(gsf_gsfet[iGsf2], gsfsc_eta[iGsf2]);
          if (fabs(gsfsc_eta[iGsf1]) > 1.56 && fabs(gsfsc_eta[iGsf2]) > 1.56) {
            histoGsfGsfMass.at(EE)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
            // fill the pz of the initial particle for EE-EE events
            double gsfGsfPZ = CalcPz(iGsf1, iGsf2, etShift);
            if (gsfsc_eta[iGsf1] * gsfsc_eta[iGsf2] > 0) {
              histoGsfGsfMass.at(EES)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              histoPZ.at(1)->Fill(gsfGsfPZ, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
            } else {
              histoGsfGsfMass.at(EEO)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              histoPZ.at(2)->Fill(gsfGsfPZ, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
            }
            if (PassFRPreSel(iGsf1, gsfsc_eta[iGsf1]) && PassFRPreSel(iGsf2, gsfsc_eta[iGsf2]))
              histoGsfGsfMassFR.at(EE)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 * fakeRate2);
          } else if ((fabs(gsfsc_eta[iGsf1]) < 1.442 && fabs(gsfsc_eta[iGsf2]) > 1.56) || (fabs(gsfsc_eta[iGsf1]) > 1.56 && fabs(gsfsc_eta[iGsf2]) < 1.442)) {
            histoGsfGsfMass.at(BE)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
            if (PassFRPreSel(iGsf1, gsfsc_eta[iGsf1]) && PassFRPreSel(iGsf2, gsfsc_eta[iGsf2]))
              histoGsfGsfMassFR.at(BE)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 * fakeRate2);
          } else if (fabs(gsfsc_eta[iGsf1]) < 1.442 && fabs(gsfsc_eta[iGsf2]) < 1.442) {
            histoGsfGsfMass.at(BB)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
            if (PassFRPreSel(iGsf1, gsfsc_eta[iGsf1]) && PassFRPreSel(iGsf2, gsfsc_eta[iGsf2]))
              histoGsfGsfMassFR.at(BB)->Fill(gsfGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 * fakeRate2);
          }
        }
      }

      if (iGsfNoHeep2 >= 0 && iGsfNoHeep1 >= 0 && iHeep1 < 0 && iHeep2 < 0) {
        // fill the GSF-GSF cases that do not pass the HEEP selection
        double gsfGsfMassNoHeep = CalcInvariantMass(iGsfNoHeep1, iGsfNoHeep2, etShift);
        if (p > DATA) trgTurnOnWeight = TriggerTurnOn(iGsfNoHeep1, useSecondTurnOn) * TriggerTurnOn(iGsfNoHeep2, useSecondTurnOn);
        if (gsfGsfMassNoHeep > massCut) {
          if (fabs(gsfsc_eta[iGsfNoHeep1]) < 2.5 && fabs(gsfsc_eta[iGsfNoHeep2]) < 2.5) {
            double fakeRate1 = FakeRate(gsf_gsfet[iGsfNoHeep1], gsfsc_eta[iGsfNoHeep1]);
            double fakeRate2 = FakeRate(gsf_gsfet[iGsfNoHeep2], gsfsc_eta[iGsfNoHeep2]);
            if (fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep2]) > 1.56) {
              histoGsfGsfMassNoHeep.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              if (p >= DY20 && p <= DY2000) {
                if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                  histoSectionDYCombinedGGNH.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              }
              if (gsfsc_eta[iGsfNoHeep1] * gsfsc_eta[iGsfNoHeep2] > 0) {
                histoGsfGsfMassNoHeep.at(EES)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedGGNH.at(EES)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
              } else {
                histoGsfGsfMassNoHeep.at(EEO)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedGGNH.at(EEO)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
              }
              if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1]) && PassFRPreSel(iGsfNoHeep2, gsfsc_eta[iGsfNoHeep2])) {
                histoGsfGsfMassNoHeepFR.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedGGNHFR.at(EE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                }
              }
            } else if ((fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep2]) > 1.56) || (fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep2]) < 1.442)) {
              histoGsfGsfMassNoHeep.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              if (p >= DY20 && p <= DY2000) {
                if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                  histoSectionDYCombinedGGNH.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              }
              if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1]) && PassFRPreSel(iGsfNoHeep2, gsfsc_eta[iGsfNoHeep2])) {
                histoGsfGsfMassNoHeepFR.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedGGNHFR.at(BE)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                }
              }
            } else if (fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep2]) < 1.442) {
              histoGsfGsfMassNoHeep.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              if (p >= DY20 && p <= DY2000) {
                if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                  histoSectionDYCombinedGGNH.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              }
              if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1]) && PassFRPreSel(iGsfNoHeep2, gsfsc_eta[iGsfNoHeep2])) {
                histoGsfGsfMassNoHeepFR.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedGGNHFR.at(BB)->Fill(gsfGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate1 / (1 - fakeRate1) * fakeRate2 / (1 - fakeRate2));
                }
              }
            }
          }
        }
      }

      if (iHeep1 >= 0) {
        // fill the HEEP-GSF and GSF-HEEP cases
        int iGsf = (iHeep1 == iGsf1) ? iGsf2 : iGsf1;
        float heepGsfMass = CalcInvariantMass(iHeep1, iGsf, etShift);
        if (p > DATA) trgTurnOnWeight = TriggerTurnOn(iHeep1, useSecondTurnOn) * TriggerTurnOn(iGsf, useSecondTurnOn);
        if (heepGsfMass > massCut) {
          if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iGsf]) < 2.5) {
            double fakeRate = FakeRate(gsf_gsfet[iGsf], gsfsc_eta[iGsf]);
            if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsf]) > 1.56) {
              histoHeepGsfMass.at(EE)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              if (gsfsc_eta[iHeep1] * gsfsc_eta[iGsf] > 0) {
                histoHeepGsfMass.at(EES)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              } else {
                histoHeepGsfMass.at(EEO)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              }
              if (PassFRPreSel(iGsf, gsfsc_eta[iGsf])) {
                histoHeepGsfMassFR.at(EE)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate);
              }
            } else if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsf]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsf]) < 1.442)) {
              histoHeepGsfMass.at(BE)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              if (PassFRPreSel(iGsf, gsfsc_eta[iGsf])) {
                histoHeepGsfMassFR.at(BE)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate);
              }
            } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsf]) < 1.442) {
              histoHeepGsfMass.at(BB)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              if (PassFRPreSel(iGsf, gsfsc_eta[iGsf])) {
                histoHeepGsfMassFR.at(BB)->Fill(heepGsfMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate);
              }
            }
          }
        }

        // fill the HEEP-GSF and GSF-HEEP cases where the GSF does not pass the HEEP selection
        if (iGsfNoHeep1 >= 0 && iHeep2 < 0) {
          float heepGsfMassNoHeep = CalcInvariantMass(iHeep1, iGsfNoHeep1, etShift);
          if (p > DATA) trgTurnOnWeight = TriggerTurnOn(iHeep1, useSecondTurnOn) * TriggerTurnOn(iGsfNoHeep1, useSecondTurnOn);
          if (heepGsfMassNoHeep > massCut) {
            if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iGsfNoHeep1]) < 2.5) {
              double fakeRate = FakeRate(gsf_gsfet[iGsfNoHeep1], gsfsc_eta[iGsfNoHeep1]);
              if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56) {
                histoHeepGsfMassNoHeep.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHGNH.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
                if (gsfsc_eta[iHeep1] * gsfsc_eta[iGsfNoHeep1] > 0) {
                  histoHeepGsfMassNoHeep.at(EES)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                  if (p >= DY20 && p <= DY2000) {
                    if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                      histoSectionDYCombinedHGNH.at(EES)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                  }
                } else {
                  histoHeepGsfMassNoHeep.at(EEO)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                  if (p >= DY20 && p <= DY2000) {
                    if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                      histoSectionDYCombinedHGNH.at(EEO)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                  }
                }
                if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1])) {
                  histoHeepGsfMassNoHeepFR.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                  if (p >= DY20 && p <= DY2000) {
                    if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                      histoSectionDYCombinedHGNHFR.at(EE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                  }
                }
              } else if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep1]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442)) {
                histoHeepGsfMassNoHeep.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHGNH.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
                if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1])) {
                  histoHeepGsfMassNoHeepFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                  if (p >= DY20 && p <= DY2000) {
                    if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                      histoSectionDYCombinedHGNHFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                  }
                }
                // special cuts for Laurent
                histoHeepGsfMassNoHeepNoHOverEFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                if (gsf_trackiso[iHeep1] < 1 && gsf_hovere[iHeep1] < 0.02 && fabs(gsf_deltaphi[iHeep1]) < 0.01 && calomet < 40) {
                  histoTAndPHeepGsfMassNoHeepFR.at(BE)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                }
              } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iGsfNoHeep1]) < 1.442) {
                histoHeepGsfMassNoHeep.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHGNH.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
                if (PassFRPreSel(iGsfNoHeep1, gsfsc_eta[iGsfNoHeep1])) {
                  histoHeepGsfMassNoHeepFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                  if (p >= DY20 && p <= DY2000) {
                    if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                      histoSectionDYCombinedHGNHFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                  }
                }
                // special cuts for Laurent
                histoHeepGsfMassNoHeepNoHOverEFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                if (gsf_trackiso[iHeep1] < 1 && gsf_hovere[iHeep1] < 0.02 && fabs(gsf_deltaphi[iHeep1]) < 0.01 && calomet < 40) {
                  histoTAndPHeepGsfMassNoHeepFR.at(BB)->Fill(heepGsfMassNoHeep, input.at(p).second * lumi * puWeight * trgTurnOnWeight * fakeRate / (1 - fakeRate));
                }
              }
            }
          }
        } 
      }

      if (iHeep1 >= 0 && iHeep2 >= 0) {
        // fill the HEEP-HEEP cases
        float heepHeepMass = CalcInvariantMass(iHeep1, iHeep2, etShift);
        if (p > DATA) trgTurnOnWeight = TriggerTurnOn(iHeep1, useSecondTurnOn) * TriggerTurnOn(iHeep2, useSecondTurnOn);
        if (heepHeepMass > massCut) {
          if (fabs(gsfsc_eta[iHeep1]) < 2.5 && fabs(gsfsc_eta[iHeep2]) < 2.5) {
            if (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iHeep2]) > 1.56) {
              histoHeepHeepMass.at(EE)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              // fill the correct sector of the DY sample
              if (p >= DY20 && p <= DY2000) {
                if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                  histoSectionDYCombinedHH.at(EE)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
              }
              if (gsfsc_eta[iHeep1] * gsfsc_eta[iHeep2] > 0) {
                histoHeepHeepMass.at(EES)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                // fill the correct sector of the DY sample
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHH.at(EES)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
              } else {
                histoHeepHeepMass.at(EEO)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                // fill the correct sector of the DY sample
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHH.at(EEO)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
              }
            }
            else {
              evtRegion = -1;
              if ((fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iHeep2]) > 1.56) || (fabs(gsfsc_eta[iHeep1]) > 1.56 && fabs(gsfsc_eta[iHeep2]) < 1.442)) {
                evtRegion = 1;
                ++heepCounter;
                if (p == DATA) {
                  if (runnumber < 191718) ++heepCounterNoOnlineLaser;
                  if (heepHeepMass > 800.) cout << "High mass EB-EE event (M>800GeV):";
                }
                histoHeepHeepMass.at(BE)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                // fill the correct sector of the DY sample
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHH.at(BE)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
              } else if (fabs(gsfsc_eta[iHeep1]) < 1.442 && fabs(gsfsc_eta[iHeep2]) < 1.442) {
                evtRegion = 0;
                ++heepCounter;
                if (p == DATA) {
                  if (runnumber < 191718) ++heepCounterNoOnlineLaser;
                  if (heepHeepMass > 800.) cout << "High mass EB-EB event (M>800GeV):";
                }
                histoHeepHeepMass.at(BB)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                // fill the correct sector of the DY sample
                if (p >= DY20 && p <= DY2000) {
                  if (genelemom_mass[0] >= minPtDY.at(p - DY20) && genelemom_mass[0] < minPtDY.at(p))
                    histoSectionDYCombinedHH.at(BB)->Fill(heepHeepMass, input.at(p).second * lumi * puWeight * trgTurnOnWeight);
                }
              }
              if (p == DATA) {
                if (heepHeepMass > 800.)
                  cout << " event=" << runnumber << ":" << luminosityBlock << ":" << eventnumber 
                       << " , mass=" << heepHeepMass << "GeV/c^2" 
                       << " , ele1 et=" << gsf_gsfet[iHeep1] << "GeV" 
                       << " , ele1 eta = " << gsf_eta[iHeep1] 
                       << " , ele1 phi = " << gsf_phi[iHeep1] 
                       << " , ele1 trackIso = " << gsf_trackiso[iHeep1] 
                       << " , ele1 ecalIso = " << gsf_ecaliso[iHeep1] 
                       << " , ele1 hcalIso1 = " << gsf_hcaliso1[iHeep1] 
                       << " , ele1 hcalIso2 = " << gsf_hcaliso2[iHeep1] 
                       << " , ele2 et=" << gsf_gsfet[iHeep2] << "GeV" 
                       << " , ele2 eta = " << gsf_eta[iHeep2] 
                       << " , ele2 phi = " << gsf_phi[iHeep2] 
                       << " , ele2 trackIso = " << gsf_trackiso[iHeep2] 
                       << " , ele2 ecalIso = " << gsf_ecaliso[iHeep2] 
                       << " , ele2 hcalIso1 = " << gsf_hcaliso1[iHeep2] 
                       << " , ele2 hcalIso2 = " << gsf_hcaliso2[iHeep2] 
                       << endl;
                // fill tree with good events
                hHMass = heepHeepMass;
                eleDataTree->Fill();
              }
            }
          }
        }
      }
    } // end of loop over events

    cout << "Number of HEEP-HEEP events: " << heepCounter << endl;
    if (p == DATA) {
      ratioNoOnlineLaser = (float) heepCounterNoOnlineLaser / (float) heepCounter;
      cout << "Percentage of HEEP-HEEP events with no online ECAL laser correction: " << ratioNoOnlineLaser << endl;
    }

    // write root file with good HEEP-HEEP event data
    if (p == DATA) {
      goodEvFile->cd();
      eleDataTree->Write();
    }

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
  
    if (p >= DY20 && p <= DY2000) {
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
      histoRatioMassFR.back()->SetTitle("Ratio between GSF-GSF and HEEP-GSF" + suffix.at(i) + ";m_(ee) (GeV);GSF-GSF / HEEP-GSF");
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
      histoRatioMassNoHeepFR.back()->SetTitle("Ratio between GSF-GSF non HEEP and HEEP-GSF non HEEP" + suffix.at(i) + ";m_(ee) (GeV);GSF-GSF / HEEP-GSF");
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
    if (p >= DY20 && p <= DY2000) {
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
      if (p >= DY20 && p <= DY2000) {
        histoSectionDYCombinedHH.at(i)->Write();
        histoSectionDYCombinedHGNH.at(i)->Write();
        histoSectionDYCombinedGGNH.at(i)->Write();
      }
      if (i > EE) continue;
      histoGsfGsfMassFR.at(i)->Write();
      histoGsfGsfMassNoHeepFR.at(i)->Write();
      histoHeepGsfMassFR.at(i)->Write();
      histoHeepGsfMassNoHeepFR.at(i)->Write();
      if (p >= DY20 && p <= DY2000) {
        histoSectionDYCombinedHGNHFR.at(i)->Write();
        histoSectionDYCombinedGGNHFR.at(i)->Write();
      }
      if (i > BE) continue;
      histoHeepGsfMassNoHeepNoHOverEFR.at(i)->Write();
      histoTAndPHeepGsfMassNoHeepFR.at(i)->Write();
    }
    for (vector<TH1F *>::iterator it = histoPZ.begin(); it < histoPZ.end(); ++it)
      (*it)->Write();

    if (p > DATA) {
      puMc->Write();
      puMcNorm->Write();
      puWeights->Write();
    }
  } // end of loop over input files
  outFile->cd();
  puData->Write();
  puDataNorm->Write();

  outFile->mkdir("combinations");  // WORKAROUND - should be in cd() but there seems to be a bug in root 4.27 that gives "R__unzip: error in header" in the rescaling process when histograms are in the base directory
  outFile->cd("combinations");

  vector<TH1F *> histoDYCombined;
  sStream << "# of events / pb^{-1} " << ((massMax - massMin) / nBins) << "GeV"; 
  for (unsigned int i = BBBE; i <= EEO; ++i) {
    histoDYCombined.push_back(new TH1F("histoDYCombined" + acroSuffix.at(i), "Combined Mass for DY samples" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
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
  sStream << "# of events * fake rate/ pb^{-1} " << ((massMax - massMin) / nBins) << "GeV"; 
  // recalculate fake rate with corrections
  vector<TH1F *> histosHeepGsfCorr;
  vector<TH1F *> histosGsfGsfCorr;
  for (unsigned int i = BBBE; i <= EE; ++i) {
    histosHGNHFRDYCombined.push_back(new TH1F("histoHGNHFRDYCombined" + acroSuffix.at(i), "Combined HEEP-GSF(no HEEP) with fake rate from DY samples" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray));
    histosHGNHFRDYCombined.back()->GetYaxis()->SetTitle(sStream.str().c_str());
    for (vector<vector<TH1F *> >::iterator it = histosSectionDYCombinedHGNHFR.begin(); it < histosSectionDYCombinedHGNHFR.end(); ++it) {
      it->at(i)->GetYaxis()->SetTitle(sStream.str().c_str());
      histosHGNHFRDYCombined.back()->Add(it->at(i));
    }

    histosHeepGsfCorr.push_back(new TH1F((*histosHeepGsfMassNoHeepFR.at(DATA).at(i)) - (*histosHGNHFRDYCombined.back())));
    histosHeepGsfCorr.back()->SetName("histoHeepGsfCorr" + acroSuffix.at(i));
    histosHeepGsfCorr.back()->SetTitle("Corrected HEEP-GSF(no HEEP) fake rate" + suffix.at(i) + ";m_(ee) (GeV)");

    //generate combined gamma HEEP-GSF non HEEP fake rate histogram
    TH1F *histoHGNHFRGammaCombined = new TH1F("histoHGNHFRGammaCombined" + acroSuffix.at(i), "Combined HEEP-GSF(no HEEP) with fake rate from gamma samples" + suffix.at(i) + ";m_(ee) (GeV)", nBins, binArray);
    for (unsigned int j = G30T50; j <= G1800; ++j) {
      histoHGNHFRGammaCombined->Add(histosHeepGsfMassNoHeepFR.at(j).at(i));
    }
    histosHGNHFRGammaCombined.push_back(new TH1F(*histoHGNHFRGammaCombined));

    histosGsfGsfCorr.push_back(new TH1F((*histosGsfGsfMassNoHeepFR.at(DATA).at(i)) + (*histosHeepGsfMassNoHeepFR.at(WJETS).at(i)) + (*histoHGNHFRGammaCombined)));
    histosGsfGsfCorr.back()->SetName("histoGsfGsfCorr" + acroSuffix.at(i));
    histosGsfGsfCorr.back()->SetTitle("Corrected GSF-GSF(no HEEP) fake rate" + suffix.at(i) + ";m_(ee) (GeV)");

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

    histoRatioMassCorrFR.back()->SetTitle("Ratio between corrected GSF-GSF and corrected HEEP-GSF" + suffix.at(i) + ";m_(ee) (GeV);GSF-GSF / HEEP-GSF");
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
InvariantMass::NormalizeToZPeak(const float& lowMass, float highMass, const char* inputFile, vector<string> folders, int ratioRebin)
{
  float factor = 1.;
  const unsigned int DATA = 0;
  const unsigned int TTBAR = 11;
  const unsigned int DYTT = 12;
  const unsigned int WW = 13;
  const unsigned int WZ = 14;
  const unsigned int ZZ = 15;
  const unsigned int TW = 16;
  const unsigned int WJETS = 17;

  vector<TString> acroSuffix;
  acroSuffix.push_back("");
  acroSuffix.push_back("BB");
  acroSuffix.push_back("BE");
  acroSuffix.push_back("EE");

  vector<TString> noHeepStr;
  noHeepStr.push_back("");
  noHeepStr.push_back("NoHeep");

  highMass -= 0.001; // to exclude the bin with highMass if it starts there
  if (highMass <= lowMass) {
    cout << endl << "Normalization to Z peak range is zero or negative. Will apply no correction factor." << endl << endl;
    return 1;
  }

  folders.push_back("combinations");

  TFile *file = new TFile(inputFile, "update");

  // get the necessary histograms
  file->cd(folders.at(DATA).c_str());
  TH1F *dataHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *dataHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(TTBAR).c_str());
  TH1F *ttbarBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *ttbarBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(DYTT).c_str());
  TH1F *dyttBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *dyttBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(WW).c_str());
  TH1F *wwBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *wwBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(WZ).c_str());
  TH1F *wzBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *wzBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(ZZ).c_str());
  TH1F *zzBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *zzBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(TW).c_str());
  TH1F *twBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *twBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd(folders.at(WJETS).c_str());
  TH1F *wjetBgHistoBB = (TH1F *)gDirectory->Get("histoHeepHeepMassBB");
  TH1F *wjetBgHistoBE = (TH1F *)gDirectory->Get("histoHeepHeepMassBE");
  file->cd("combinations");
  TH1F *dyBgHistoBB = (TH1F *)gDirectory->Get("histoDYCombinedBB");
  TH1F *dyBgHistoBE = (TH1F *)gDirectory->Get("histoDYCombinedBE");
  TH1F *qcdBgHistoBB = (TH1F *)gDirectory->Get("histoGsfGsfCorrBB");
  TH1F *qcdBgHistoBE = (TH1F *)gDirectory->Get("histoGsfGsfCorrBE");

  double intData = 0;
  double intDataBB = 0;
  double intDataBE = 0;
  double intBg = 0;
  double intBgBB = 0;
  double intBgBE = 0;
  double intDyBB = 0;
  double intDyBE = 0;
  double intTtbarBB = 0;
  double intTtbarBE = 0;
  double intDyttBB = 0;
  double intDyttBE = 0;
  double intWwBB = 0;
  double intWwBE = 0;
  double intWzBB = 0;
  double intWzBE = 0;
  double intZzBB = 0;
  double intZzBE = 0;
  double intTwBB = 0;
  double intTwBE = 0;
  double intWjetBB = 0;
  double intWjetBE = 0;
  double intQcdBB = 0;
  double intQcdBE = 0;

  // integrate over range
  intDataBB = dataHistoBB->Integral(dataHistoBB->FindBin(lowMass), dataHistoBB->FindBin(highMass));
  intDyBB = dyBgHistoBB->Integral(dyBgHistoBB->FindBin(lowMass), dyBgHistoBB->FindBin(highMass));
  intTtbarBB = ttbarBgHistoBB->Integral(ttbarBgHistoBB->FindBin(lowMass), ttbarBgHistoBB->FindBin(highMass));
  intDyttBB = dyttBgHistoBB->Integral(dyttBgHistoBB->FindBin(lowMass), dyttBgHistoBB->FindBin(highMass));
  intWwBB = wwBgHistoBB->Integral(wwBgHistoBB->FindBin(lowMass), wwBgHistoBB->FindBin(highMass));
  intWzBB = wzBgHistoBB->Integral(wzBgHistoBB->FindBin(lowMass), wzBgHistoBB->FindBin(highMass));
  intZzBB = zzBgHistoBB->Integral(zzBgHistoBB->FindBin(lowMass), zzBgHistoBB->FindBin(highMass));
  intTwBB = twBgHistoBB->Integral(twBgHistoBB->FindBin(lowMass), twBgHistoBB->FindBin(highMass));
  intWjetBB = wjetBgHistoBB->Integral(wjetBgHistoBB->FindBin(lowMass), wjetBgHistoBB->FindBin(highMass));
  intQcdBB = qcdBgHistoBB->Integral(qcdBgHistoBB->FindBin(lowMass), qcdBgHistoBB->FindBin(highMass));

  intDataBE = dataHistoBE->Integral(dataHistoBE->FindBin(lowMass), dataHistoBE->FindBin(highMass));
  intDyBE = dyBgHistoBE->Integral(dyBgHistoBE->FindBin(lowMass), dyBgHistoBE->FindBin(highMass));
  intTtbarBE = ttbarBgHistoBE->Integral(ttbarBgHistoBE->FindBin(lowMass), ttbarBgHistoBE->FindBin(highMass));
  intDyttBE = dyttBgHistoBE->Integral(dyttBgHistoBE->FindBin(lowMass), dyttBgHistoBE->FindBin(highMass));
  intWwBE = wwBgHistoBE->Integral(wwBgHistoBE->FindBin(lowMass), wwBgHistoBE->FindBin(highMass));
  intWzBE = wzBgHistoBE->Integral(wzBgHistoBE->FindBin(lowMass), wzBgHistoBE->FindBin(highMass));
  intZzBE = zzBgHistoBE->Integral(zzBgHistoBE->FindBin(lowMass), zzBgHistoBE->FindBin(highMass));
  intTwBE = twBgHistoBE->Integral(twBgHistoBE->FindBin(lowMass), twBgHistoBE->FindBin(highMass));
  intWjetBE = wjetBgHistoBE->Integral(wjetBgHistoBE->FindBin(lowMass), wjetBgHistoBE->FindBin(highMass));
  intQcdBE = qcdBgHistoBE->Integral(qcdBgHistoBE->FindBin(lowMass), qcdBgHistoBE->FindBin(highMass));

  intData = intDataBB + intDataBE;
  intBgBB = intDyBB + intTtbarBB + intDyttBB + intWwBB + intWzBB + intZzBB + intTwBB + intWjetBB + intQcdBB;
  intBgBE = intDyBE + intTtbarBE + intDyttBE + intWwBE + intWzBE + intZzBE + intTwBE + intWjetBE + intQcdBE;
  intBg = intBgBB + intBgBE;

  // calculate the noramlization factor
  factor = intData / intBg;

  cout << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "# of events (weighted) in the Z peak from " << lowMass << "GeV to " << highMass + 0.001 << "GeV" << endl;
  cout << " Data       " << intData << endl;
  cout << " DY         " << intDyBB + intDyBE << endl;
  cout << " ttbar      " << intTtbarBB + intTtbarBE  << endl;
  cout << " DY->TauTau " << intDyttBB + intDyttBE  << endl;
  cout << " WW         " << intWwBB + intWwBE  << endl;
  cout << " WZ         " << intWzBB + intWzBE  << endl;
  cout << " ZZ         " << intZzBB + intZzBE  << endl;
  cout << " TW         " << intTwBB + intTwBE  << endl;
  cout << " W+jets     " << intWjetBB + intWjetBE  << endl;
  cout << " QCD        " << intQcdBB + intQcdBE  << endl << endl;;
  cout << "Applied correction factor for MC: " << factor << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "# of EB-EB events (weighted) in the Z peak from " << lowMass << "GeV to " << highMass + 0.001 << "GeV" << endl;
  cout << " Data       " << intDataBB << endl;
  cout << " DY         " << intDyBB << endl;
  cout << " ttbar      " << intTtbarBB << endl;
  cout << " DY->TauTau " << intDyttBB << endl;
  cout << " WW         " << intWwBB << endl;
  cout << " WZ         " << intWzBB << endl;
  cout << " ZZ         " << intZzBB << endl;
  cout << " TW         " << intTwBB << endl;
  cout << " W+jets     " << intWjetBB << endl;
  cout << " QCD        " << intQcdBB << endl << endl;
  cout << "EB-EB correction factor for MC: " << intDataBB / intBgBB << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "# of EB-EE events (weighted) in the Z peak from " << lowMass << "GeV to " << highMass + 0.001 << "GeV" << endl;
  cout << " Data       " << intDataBE << endl;
  cout << " DY         " << intDyBE << endl;
  cout << " ttbar      " << intTtbarBE << endl;
  cout << " DY->TauTau " << intDyttBE << endl;
  cout << " WW         " << intWwBE << endl;
  cout << " WZ         " << intWzBE << endl;
  cout << " ZZ         " << intZzBE << endl;
  cout << " TW         " << intTwBE << endl;
  cout << " W+jets     " << intWjetBE << endl;
  cout << " QCD        " << intQcdBE << endl << endl;
  cout << "EB-EE correction factor for MC: " << intDataBE / intBgBE << endl;
  cout << endl;
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
InvariantMass::CalcInvariantMass (const int& iEle1, const int& iEle2, bool &etShift)
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

double
InvariantMass::CalcPz (const int& iEle1, const int& iEle2, bool &etShift)
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
  
  return (ele1+ele2).Pz();
}

double
InvariantMass::FakeRate (const float& et, const float& eta)
{
  // fake rate 19.3/fb
  if (fabs(eta) < 1.5) {
    if (et < 189.3) return 0.0179 - 0.000056 * et;
    return 0.0073;
  } else if (abs(eta) < 2.0) {
    if (et < 96.6) return exp(-2.31 - 0.011 * et);
    else if (et < 178.0) return 0.040 - 0.000059 * et;
    return 0.0295;
  } else if (abs(eta) < 2.5) {
    if(et < 115.4) return 0.099 - 0.00035 * et;
    else return 0.0586;
  } else return 0; 
}

bool
InvariantMass::PassFRPreSel(const int &n, const float &eta)
{
  // fake rate preselection
  const float hOverECutEB = 0.15;
  const float hOverECutEE = 0.1;
  const float sieieCutEB = 0.013;
  const float sieieCutEE = 0.034;
  const float dxyCutEB = 0.02;
  const float dxyCutEE = 0.05;
  const int missHitsCutEB = 1;
  const int missHitsCutEE = 1;

  if (fabs(eta) < 1.442) {
    if (gsf_hovere[n] > hOverECutEB) return false;
    if (gsf_nLostInnerHits[n] > missHitsCutEB) return false;
    if (gsf_sigmaIetaIeta[n] > sieieCutEB) return false;
    if (fabs(gsf_dxy_firstPVtx[n]) > dxyCutEB) return false;
  }
  else if (fabs(eta) > 1.56) {
    if (gsf_hovere[n] > hOverECutEE) return false;
    if (gsf_nLostInnerHits[n] > missHitsCutEE) return false;
    if (gsf_sigmaIetaIeta[n] > sieieCutEE) return false;
    if (fabs(gsf_dxy_firstPVtx[n]) > dxyCutEE) return false;
  }

  return true;
}

bool
InvariantMass::PassHEEP(const int &n)
{
  //int selection = 0;  // HEEP v4.0
  int selection = 1;  // HEEP v4.1

// HEEP selection v3.2
//  // barrel
//  float bar_et = 35.;
//  float bar_hoE = 0.05;
//  float bar_DEta = 0.005;
//  float bar_DPhi = 0.06;
//  float bar_e2x5e5x5 = 0.94;
//  float bar_e1x5e5x5 = 0.83;
//  float bar_isoEcalHcal1_1 = 2.;
//  float bar_isoEcalHcal1_2 = 0.03;
//  float bar_isoTrack = 5.;
//  int bar_missInnerHits = 0;
//
//  // endcap
//  float end_et = 40.;
//  float end_hoE = 0.05;
//  float end_DEta = 0.007 ;
//  float end_DPhi = 0.06;
//  float end_e2x5e5x5 = 0.;
//  float end_e1x5e5x5 = 0.;
//  float end_sigmaietaieta = 0.03;
//  float end_isoEcalHcal1_1 = 2.5;
//  float end_isoEcalHcal1_2 = 0.03;
//  float end_isoTrack = 5.;
//  int end_missInnerHits = 0;

  // HEEP v4.0
  // barrel
  float bar_et = 35.;
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
  float end_et = 35.;
  float end_hoE = 0.05;
  float end_DEta = 0.007 ;
  float end_DPhi = 0.06;
  float end_sigmaietaieta = 0.03;
  float end_isoEcalHcal1_1_1 = 2.5;
  float end_isoEcalHcal1_1_2 = 1.;
  float end_isoEcalHcal1_2 = 0.03;
  float end_isoEcalHcalRho = 0.28;
  float end_isoTrack = 5.;
  float end_dxy = 0.05;  // only for HEEP v4.1
  int end_missInnerHits = 0;

  if (selection == 1) {
    bar_missInnerHits = 1;
    end_missInnerHits = 1;
  }

//  HEEP v3.2
//  // barrel
//  if (fabs(gsfsc_eta[n]) < 1.442
//      && gsf_gsfet[n] > bar_et
//      && gsf_isecaldriven[n]
//      && gsf_hovere[n] < bar_hoE
//      && fabs(gsf_deltaeta[n]) < bar_DEta
//      && fabs(gsf_deltaphi[n]) < bar_DPhi
//      && (gsf_e2x5overe5x5[n] > bar_e2x5e5x5 || gsf_e1x5overe5x5[n] > bar_e1x5e5x5)
//      && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * gsf_gsfet[n])
//      && gsf_trackiso[n] < bar_isoTrack
//      && gsf_nLostInnerHits[n] <= bar_missInnerHits
//     ) return true;
//
//  // endcap
//  else if ((fabs(gsfsc_eta[n]) > 1.56 && fabs(gsfsc_eta[n]) < 2.5)
//      && gsf_gsfet[n] > end_et
//      && gsf_isecaldriven[n]
//      && gsf_hovere[n] < end_hoE
//      && fabs(gsf_deltaeta[n]) < end_DEta
//      && fabs(gsf_deltaphi[n]) < end_DPhi
//      && (gsf_e2x5overe5x5[n] > end_e2x5e5x5 || gsf_e1x5overe5x5[n] > end_e1x5e5x5)
//      && gsf_sigmaIetaIeta[n] < end_sigmaietaieta
//      && ((gsf_gsfet[n] < 50. && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < end_isoEcalHcal1_1)
//          ||
//          (gsf_gsfet[n] >= 50. && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (end_isoEcalHcal1_1 + end_isoEcalHcal1_2 * (gsf_gsfet[n] - 50.))))
//      && gsf_trackiso[n] < end_isoTrack
//      && gsf_nLostInnerHits[n] <= end_missInnerHits
//     ) return true;

  // HEEP v4.0
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
InvariantMass::Trigger(int &prescale)
{
  // trigger selection 2012
  int triggerBit = 0;
  if (0 && prescale_HLT_DoubleEle33_CaloIdL == 1 && HLT_DoubleEle33_CaloIdL == 1) {  // switched off
    prescale = prescale_HLT_DoubleEle33_CaloIdL;
    triggerBit = HLT_DoubleEle33_CaloIdL;
  } else if (prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL == 1 && HLT_DoubleEle33_CaloIdL_GsfTrkIdVL == 1) {
    prescale = prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
    triggerBit = HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  } else if (0 && prescale_HLT_DoubleEle33_CaloIdT == 1 && HLT_DoubleEle33_CaloIdT == 1) { // switched off
    prescale = prescale_HLT_DoubleEle33_CaloIdT;
    triggerBit = HLT_DoubleEle33_CaloIdT;
  }
  return triggerBit;

  cout << "Prescale alert! No unprescaled trigger found. Dropping this event." << endl;
  prescale = 0;
  return 0;
}

double 
InvariantMass::TriggerTurnOn(const int &n, const bool &useSecond)
{
  float et = gsf_gsfet[n] * cosh(gsf_eta[n]) / cosh(gsfsc_eta[n]);
  float a0 = (fabs(gsfsc_eta[n]) < 1.5) ? 0.996 : 0.9948;
  float a1 = (fabs(gsfsc_eta[n]) < 1.5) ? 34.76 : 32.74;
  float a2 = (fabs(gsfsc_eta[n]) < 1.5) ? 0.85 : 2.36;

  if (useSecond) {
    a0 = (fabs(gsfsc_eta[n]) < 1.5) ? 0.996 : 0.979;
    a1 = (fabs(gsfsc_eta[n]) < 1.5) ? 34.76 : 37.2;
    a2 = (fabs(gsfsc_eta[n]) < 1.5) ? 0.85 : 2.3;
  }

  return 0.5 * a0 * (1 + erf((et - a1) / (sqrt(2) * a2)));
}

void
InvariantMass::CorrectEnergy() 
{
  for (int n = 0; n < gsf_size; ++n)
     gsf_gsfet[n] = gsfsc_e[n] * sin(gsf_theta[n]);
}

