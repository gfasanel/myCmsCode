#include "TFile.h"
#include "TDCacheFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TPaveLabel.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "TFile.h"
#endif // __CINT__

#define DATA 0
#define TTBAR 1
#define ZTT 2
#define WW 3
#define WZ 4
#define TW 5
#define WJET 6
#define ZMM 7
#define ZEE 8
#define ZZ 9
//#define QCD 10
#define QCD 9

#define ALL 0
#define LS 1
#define OS 2

using namespace std;

float CalcSystErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked = false);
float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked = false);
int Trigger(TFile *inFile, unsigned int &entry, int &prescale, unsigned int *trig, const int &selector = 0);

pair<unsigned int, unsigned int> runs_HLT_Mu15_Photon20_CaloIdL(99999999, 0);
pair<unsigned int, unsigned int> runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL(99999999, 0);
pair<unsigned int, unsigned int> runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL(99999999, 0);
//RUN ID
unsigned int c_runnumber;
unsigned int c_eventnumber;
unsigned int c_luminosityBlock;
//trigger
int c_HLT_Mu15;
int c_HLT_Mu30;

void emuSpectrum()
{
   // parameters /////////////////////////////////////////////////////////////
   float LumiFactor = 4684.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 4800.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 5035.; //Lumi in pb-1   -new LUMI FROM GOLDEN JSON
   //float LumiFactor = 2511.; //Lumi in pb-1   -LUMI run2011B FROM GOLDEN JSON
   //float LumiFactor = 2173.; //Lumi in pb-1   -LUMI run2011A FROM GOLDEN JSON
   //float LumiFactor = 4699.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 3534.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 3190.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 2179.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 1932.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 702.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON

   // DATA file
   TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-180252_7TeV_Collisions11_JSON_gct1_13_4684pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-180252_7TeV_Collisions11_JSON_4699pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-178078_7TeV_Collisions11_JSON_3534pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-177515_7TeV_Collisions11_JSON_3190pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6-Cert_160404-173692_7TeV_Collisions11_JSON_2179pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco+05Aug2011+PromptReco-AOD-Cert_160404-173244_7TeV_1932pb-1.root";
   //TString dataFile = "/user_mnt/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/total_MuEG-160404-163869-ReReco10May-GoldenJSON-27May-191pb_+_MuEG-165088-166861-PromptV4-GoldenJSON-17Jun-511pb__702pb.root";

   //string outfileName = "testEmuSpec";
   //string outfileName = "testEmuSpecSummerFallMix";
   //string outfileName = "testEmuSpecPureSummer11";
   string outfileName = "test";

   unsigned nPVtxMax = 25; // max number of primary vertices

   // scale factors
   // TODO errors and distinction EB - EE
   float Elec_trigger = 0.94;
   float Elec_ScaleFactor = 1.008 * 0.978;
   float Muon_ScaleFactor = 0.985;
   float Lumi_ScaleFactor = 1.0;

   // systematical errors
   vector<float> systErrMC;
   systErrMC.push_back(0.15);  //ttbar
   systErrMC.push_back(0.);    //z->tt  //FIXME
   systErrMC.push_back(0.);    //WW     //FIXME
   systErrMC.push_back(0.);    //WZ     //FIXME
   systErrMC.push_back(0.15);  //tW     //FIXME
   systErrMC.push_back(0.);    //WJets  //FIXME
   systErrMC.push_back(0.);    //Z->mm  //FIXME
   systErrMC.push_back(0.);    //Z->ee  //FIXME
   systErrMC.push_back(0.);    //ZZ     //FIXME
   float systErrLumi = 0.045;
   float systErrEff = 0.02;

   bool useQCDShape = true;
   bool calcQCDScaleFactor = true;
   float QCD_ScaleFactor = 1.10; // overwritten if calcQCDScaleFactor == true

   bool usePUInfo = true;
   bool generatePUFile = false;

   // selection cuts /////////////////////////////////////////////////////////
   float minInvMass = 0.;

   //MUON selection
   float muon_et = 35.;
   float muon_etaMax = 2.4;
   int muon_nHitsMinGlobal = 0;
   int muon_nHitsMinPixel = 1;
   int muon_nHitsMinMuon = 1;
   int muon_nLayersMin = 9; 
   float muon_impactParamMax = 0.2;   // in cm
   int muon_nSegMatchMin = 2;
   float muon_relIsoCutMax = 0.1;

   //HEEP selection v3.2
   //BARREL
   float bar_et = 35.;
   float bar_hoE = 0.05;
   float bar_DEta = 0.005;
   float bar_DPhi = 0.06;
   float bar_e2x5e5x5 = 0.94;
   float bar_e1x5e5x5 = 0.83;
   float bar_isoEcalHcal1_1 = 2.;
   float bar_isoEcalHcal1_2 = 0.03;
   float bar_isoTrack = 5.;
   int bar_missInnerHits = 0;

   //ENDCAP
   float end_et = 40.;
   float end_hoE = 0.05;
   float end_DEta = 0.007 ;
   float end_DPhi = 0.06;
   float end_e2x5e5x5 = 0.;
   float end_e1x5e5x5 = 0.;
   float end_sigmaietaieta = 0.03;
   float end_isoEcalHcal1_1 = 2.5;
   float end_isoEcalHcal1_2 = 0.03;
   float end_isoTrack = 5.;
   int end_missInnerHits = 0;
   ///////////////////////////////////////////////////////////////////////////

   TH1::SetDefaultSumw2(kTRUE);

   float MCemuScaleFactor = Elec_ScaleFactor * Muon_ScaleFactor * Lumi_ScaleFactor; 
   float MCScaleFactor = Elec_ScaleFactor * Muon_ScaleFactor * Lumi_ScaleFactor; 

   // calculate rate of syst errors
   float systErrLuEff = sqrt(systErrLumi*systErrLumi + systErrEff*systErrEff);
   vector<float> systErrMCLuEff;
   for (unsigned int i = TTBAR - 1; i < ZEE; ++i) 
      systErrMCLuEff.push_back(sqrt(systErrMC[i]*systErrMC[i] + systErrLuEff*systErrLuEff));

   ///////////////////////////////////////////////////////////////////////////
   // INPUT FILES
   vector<pair<TFile *, float> > input;

   // DATA
   input.push_back(make_pair(new TFile(dataFile, "open"), 1.)); //DATA        0 black

   // MC
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13_emuSkim.root", "open"), 4.441E-5)); //TTbar       1 red         (3701947 event - xsect NNLO 164.4pb) -- 22.02.2012
   input.push_back(make_pair(new TFile("/user/treis/mcsamples/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13_emuSkim.root", "open"), 4.40309E-5)); //TTbar       1 red         (3701947 event - xsect NNLO 163pb) -- 20.11.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v2_AODSIM_gct1_13_emuSkim.root", "open"), 0.000044031)); //TTbar       1 red         (3701947 event - xsect NNLO 163pb) -- 17.10.2011

   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13.root", "open"), 8.35612E-5)); //Ztautau     2 green       (19937479 event  - xsect 1666pb) -- 24.11.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_18_emuSkim.root", "open"), 1.5537E-4)); //Z+Jets to LL     2 green       (19617630 event  - xsect 3048pb AN-11-390) -- 27.2.2012  instead of all the othe r DY samples
   input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_AODSIM_gct1_13.root", "open"), 8.19666E-4)); //Ztautau     2 green       (2032536 event  - xsect 1666pb) -- 4.12.2011

   input.push_back(make_pair(new TFile("/user/treis/mcsamples/WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13.root", "open"), 2.20727E-5)); //WW          3 dark blue   (210667 event - xsect NLO 4.65pb (AN-11-259) ) -- 20.11.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_gct1_13.root", "open"), 2.20727E-5)); //WW          3 dark blue   (210667 event - xsect NLO 4.65pb (AN-11-259) ) -- 4.12.2011

   input.push_back(make_pair(new TFile("/user/treis/mcsamples/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13_emuSkim.root", "open"), 5.41102E-7)); //WZ          4 yellow      (1097759 event - xsect NLO 0.594pb (AN-11-259) ) -- 20.11.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_gct1_13.root", "open"), 0.000002901)); //WZ          4 yellow      (204725 event - xsect NLO 0.594pb (AN-11-259) ) -- 17.10.2011

   input.push_back(make_pair(new TFile("/user/treis/mcsamples/T-and-Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13.root","open"), 1.38228E-5)); //tW          5 pink        (814390+323401 event - xsect NNLO 7.87pb+7.87pb (note tW) ) -- 20.11.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/T-and-Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_gct1_13.root","open"), 0.000009690)); //tW          5 pink        (814390+809984 event - xsect NNLO 7.87pb+7.87pb (note tW) ) -- 14 NOV 2011

   input.push_back(make_pair(new TFile("/user/treis/mcsamples/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13.root", "open"), 3.84951E-4)); //W+jet       6 dark green  (81345381 event - xsect NNLO 31314pb) -- 20.11.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_gct1_13.root", "open"), 0.000384917)); //W+jet       6 dark green  (81352581 event - xsect NNLO 31314pb) -- 17.10.2011

   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13_emuSkim.root","open"), 5.60121E-5)); //Zmumu       7 light blue  (29743564 event  - xsect NNLO 1666pb ) -- 24 NOV 2011)
   input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2_AODSIM_gct1_13_emuSkim.root","open"), 0.00077549)); //Zmumu       7 light blue  (2148325 event  - xsect NNLO 1666pb ) -- 14 NOV 2011)

   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13_emuSkim.root","open"), 5.64799E-5)); //Zee         8 cyan        (29497207 event - xsect NNLO 1666pb ) -- 24 NOV 2011
   input.push_back(make_pair(new TFile("/user/treis/mcsamples/DYToEE_M-20_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2_AODSIM_gct1_13_emuSkim.root","open"), 0.00073630)); //Zee         8 cyan        (2262653 event - xsect NNLO 1666pb ) -- 14 NOV 2011

   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/ZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13.root","open"), 1.4394E-6)); //ZZ          9 violett     (4098843 event - xsect 5.9pb (AN-11-472) ) -- 22.02.2012
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/ZZTo4e_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1_AODSIM_gct1_13.root","open"), 0.000000164)); //ZZ          9 violett     (499929 event - xsect NLO 0.082pb (AN-11-259) ) -- 17.10.2011
   //input.push_back(make_pair(new TFile("/user/treis/mcsamples/ZZTo4e_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_gct1_13.root","open"), 0.000000164)); //ZZ          9 violett     (499929 event - xsect NLO 0.082pb (AN-11-259) ) -- 17.10.2011

   int nbFile = input.size();
   ///////////////////////////////////////////////////////////////////////////

   // strings for histogram names
   vector<TString> suffix;
   suffix.push_back("data");
   suffix.push_back("ttbar");
   suffix.push_back("ztautau");
   suffix.push_back("ww");
   suffix.push_back("wz");
   suffix.push_back("tw");
   suffix.push_back("wjets");
   suffix.push_back("zmumu");
   suffix.push_back("zee");
   //suffix.push_back("zz");
   suffix.push_back("qcd");

   TString sign[3] = {"", "LS_", "OS_"};

   stringstream ssGoodHeepFileName;
   ssGoodHeepFileName << "goodEmuEvents" << LumiFactor << "pb-1.root";
   TFile *goodEvFile = new TFile(ssGoodHeepFileName.str().c_str(), "recreate");
   goodEvFile->cd();
   float emuInvMass = 0.;
   TTree *emuDataTree = new TTree("emuDataTree", "emuDataTree");
   emuDataTree->Branch("runnr", &c_runnumber, "runnr/i");
   emuDataTree->Branch("eventnr", &c_eventnumber, "eventnr/i");
   emuDataTree->Branch("lumiSec", &c_luminosityBlock, "lumiSec/i");
   emuDataTree->Branch("mass", &emuInvMass, "mass/F");
   TTree *emuLSDataTree = new TTree("emuLSDataTree", "emuLSDataTree");
   emuLSDataTree->Branch("runnr", &c_runnumber, "runnr/i");
   emuLSDataTree->Branch("eventnr", &c_eventnumber, "eventnr/i");
   emuLSDataTree->Branch("lumiSec", &c_luminosityBlock, "lumiSec/i");
   emuLSDataTree->Branch("mass", &emuInvMass, "mass/F");
   TTree *emuOSDataTree = new TTree("emuOSDataTree", "emuOSDataTree");
   emuOSDataTree->Branch("runnr", &c_runnumber, "runnr/i");
   emuOSDataTree->Branch("eventnr", &c_eventnumber, "eventnr/i");
   emuOSDataTree->Branch("lumiSec", &c_luminosityBlock, "lumiSec/i");
   emuOSDataTree->Branch("mass", &emuInvMass, "mass/F");
   emuOSDataTree->Branch("HLT_Mu15", &c_HLT_Mu15, "HLT_Mu15/I");
   emuOSDataTree->Branch("HLT_Mu30", &c_HLT_Mu30, "HLT_Mu30/I");

   // counting variables
   int nb_plus_plus = 0;
   int nb_plus_minus = 0;
   int nb_minus_plus = 0;
   int nb_minus_minus = 0;

   // histogram containers
   vector<vector<TH1F *> > emuMass;
   vector<vector<TH1F *> > emuMassEB;
   vector<vector<TH1F *> > emuMassEE;
   vector<vector<TH1F *> > emu_mass_accVSgood;

   vector<vector<TH1F *> > met;
   vector<vector<TH1F *> > nVtx;
   vector<vector<TH1F *> > dPhi;
   vector<vector<TH1F *> > elePt;
   vector<vector<TH1F *> > eleEta;
   vector<vector<TH1F *> > elePhi;
   vector<vector<TH1F *> > eleId1;
   vector<vector<TH1F *> > eleId2;
   vector<vector<TH1F *> > eleId3;
   vector<vector<TH1F *> > eleIso1;
   vector<vector<TH1F *> > eleIso2;
   vector<vector<TH1F *> > eleIso3;

   vector<vector<TH1F *> > muIsoCombRel;
   vector<vector<TH1F *> > muPtEleOPtMu;
   vector<vector<TH1F *> > muPtPlusOPtMinus;
   vector<vector<TH1F *> > muPt;
   vector<vector<TH1F *> > muEta;
   vector<vector<TH1F *> > muPhi;
   vector<vector<TH1F *> > muId1;
   vector<vector<TH1F *> > muId2;
   vector<vector<TH1F *> > muId3;
   vector<vector<TH1F *> > muIso1;
   vector<vector<TH1F *> > muIso2;
   vector<vector<TH1F *> > muIso3;

   vector<vector<TH1F *> > numOfJets;
   vector<vector<TH1F *> > numOfJetsPt15;

   //vector<TH1F *> emu_plus_plus;
   //vector<TH1F *> emu_plus_minus;
   //vector<TH1F *> emu_minus_plus;
   //vector<TH1F *> emu_minus_minus;

   //
   TH1F *emu_dilepton = new TH1F("emu_dilepton", "emu_dilepton", 150, 0., 1500.);
   TH1F *emu_ewk = new TH1F("emu_ewk", "emu_ewk", 150, 0., 1500.);
   TH1F *emu_jet = new TH1F("emu_jet", "emu_jet", 150, 0., 1500.);

   vector<TH1F *> emuLoose_nValidPv;
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_data_nValidPv", "emuLoose_data_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_ttbar_nValidPv", "emuLoose_ttbar_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_ztautau_nValidPv", "emuLoose_ztautau_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_ww_nValidPv", "emuLoose_ww_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_wz_nValidPv", "emuLoose_wz_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_tw_nValidPv", "emuLoose_tw_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_wjets_nValidPv", "emuLoose_wjets_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_zmumu_nValidPv", "emuLoose_zmumu_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_nValidPv.push_back(new TH1F("emuLoose_zee_nValidPv", "emuLoose_zee_nValidPv", nPVtxMax, 0., nPVtxMax));
   //emuLoose_nValidPv.push_back(new TH1F("emuLoose_zz_nValidPv", "emuLoose_zz_nValidPv", nPVtxMax, 0., nPVtxMax));

   vector<TH1F *> emuLoose_dataOverX_nValidPv;
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverTtbar_nValidPv", "emuLoose_dataOverTtbar_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverZtautau_nValidPv", "emuLoose_dataOverZtautau_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverWw_nValidPv", "emuLoose_dataOverWw_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverWz_nValidPv", "emuLoose_dataOverWz_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverTw_nValidPv", "emuLoose_dataOverTw_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverWjets_nValidPv", "emuLoose_dataOverWjets_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverZmumu_nValidPv", "emuLoose_dataOverZmumu_nValidPv", nPVtxMax, 0., nPVtxMax));
   emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverZee_nValidPv", "emuLoose_dataOverZee_nValidPv", nPVtxMax, 0., nPVtxMax));
   //emuLoose_dataOverX_nValidPv.push_back(new TH1F("emuLoose_dataOverZz_nValidPv", "emuLoose_dataOverZz_nValidPv", nPVtxMax, 0., nPVtxMax));

   //GLOBAL
   int c_nJetsAKT_pt15;
   float c_calomet;
   float c_met;
   float c_pthat;
   float c_bsposx;
   float c_bsposy;
   float c_bsposz;

   //JETS AKT
   int c_jetAKT_size;
   float c_jetAKT_pt[100];
   float c_jetAKT_eta[100];
   float c_jetAKT_phi[100];
   float c_jetAKT_em[100];

   //PRIM VTX
   int c_pvsize;
   int c_pvz[20];
   bool c_pv_isValid[20];
   float c_pv_ndof[20];
   int c_pv_nTracks[20];
   float c_pv_normChi2[20];
   int c_pv_totTrackSize[20];

   //GSF
   int c_gsf_size;
   float c_gsf_gsfet[20];
   float c_gsf_px[20];
   float c_gsf_py[20];
   float c_gsf_pz[20];
   float c_gsf_pt[20];
   float c_gsf_eta[20];
   float c_gsf_theta[20];
   float c_gsf_phi[20];
   float c_gsf_dz[20];
   int c_gsf_isecaldriven[20];
   int c_gsf_istrackerdriven[20];
   int c_gsf_isEB[20];
   int c_gsf_isEE[20];
   int c_gsf_charge[20];
   float c_gsf_deltaeta[20];
   float c_gsf_deltaphi[20];
   float c_gsf_e1x5overe5x5[20];
   float c_gsf_e2x5overe5x5[20];
   float c_gsf_sigmaetaeta[20];
   float c_gsf_sigmaIetaIeta[20];
   float c_gsf_hovere[20];
   float c_gsf_eOVERp[20];
   float c_gsf_vz[20];
   int c_gsf_nHits[20];
   int c_gsf_nLostInnerHits[20];
   float c_gsf_fBrem[20];
   float c_gsf_ecaliso[20];
   float c_gsf_hcaliso1[20];
   float c_gsf_hcaliso2[20];
   float c_gsf_trackiso[20];

   float c_gsfsc_px[20];
   float c_gsfsc_py[20];
   float c_gsfsc_pt[20];
   float c_gsfsc_eta[20];
   float c_gsfsc_phi[20];

   bool c_gsfpass_ID[20];
   bool c_gsfpass_ISO[20];
   bool c_gsfpass_HEEP[20];

   //MUONS
   int c_muon_size;
   float c_muon_pt[20];
   float c_muon_eta[20];
   float c_muon_phi[20];
   float c_muon_theta[20];
   float c_muon_ptError[20];
   float c_muon_etaError[20];
   float c_muon_phiError[20];
   float c_muon_thetaError[20];
   float c_muon_outerPt[20];
   float c_muon_outerEta[20];
   float c_muon_outerPhi[20];
   float c_muon_outerTheta[20];
   float c_muon_px[20];
   float c_muon_py[20];
   float c_muon_pz[20];
   int c_muon_charge[20];
   int c_muon_nhitstrack[20];
   int c_muon_nhitspixel[20];
   int c_muon_nhitstotal[20];
   int c_muon_nhitsmuons[20];
   int c_muon_nlayerswithhits[20];
   int c_muon_nSegmentMatch[20];
   bool c_muon_isTrackerMuon[20];
   float c_muon_chi2[20];
   int c_muon_ndof[20];
   float c_muon_normChi2[20];
   float c_muon_d0[20];
   float c_muon_d0Error[20];
   float c_muon_dzError[20];
   float c_muon_dxyError[20];
   float c_muon_dz_cmsCenter[20];
   float c_muon_dz_beamSpot[20];
   float c_muon_dz_firstPVtx[20];
   float c_muon_dxy_cmsCenter[20];
   float c_muon_dxy_beamSpot[20];
   float c_muon_dxy_firstPVtx[20];
   float c_muon_trackIso03[20];
   float c_muon_trackIso05[20];
   float c_muon_emIso03[20];
   float c_muon_emIso05[20];
   float c_muon_hadIso03[20];
   float c_muon_hadIso05[20];
   float c_muon_trackIso03_ptInVeto[20];
   float c_muon_trackIso05_ptInVeto[20];
   float c_muon_emIso03_ptInVeto[20];
   float c_muon_emIso05_ptInVeto[20];
   float c_muon_hadIso03_ptInVeto[20];
   float c_muon_hadIso05_ptInVeto[20];
   float c_muon_innerPosx[20];
   float c_muon_innerPosy[20];
   float c_muon_innerPosz[20];

   //RUN ID
   TBranch        *b_runnumber;
   TBranch        *b_eventnumber;
   TBranch        *b_luminosityBlock;

   //GLOBAL
   TBranch        *b_nJetsAKT_pt15;
   TBranch        *b_calomet;
   TBranch        *b_met;
   TBranch        *b_pthat;
   TBranch        *b_bsposx;
   TBranch        *b_bsposy;
   TBranch        *b_bsposz;

   //JETS AKT
   TBranch        *b_jetAKT_size;
   TBranch        *b_jetAKT_pt;
   TBranch        *b_jetAKT_eta;
   TBranch        *b_jetAKT_phi;
   TBranch        *b_jetAKT_em;

   //PRIM VTX
   TBranch        *b_pvsize;
   TBranch        *b_pvz;
   TBranch        *b_pv_isValid;
   TBranch        *b_pv_ndof;
   TBranch        *b_pv_nTracks;
   TBranch        *b_pv_normChi2;
   TBranch        *b_pv_totTrackSize;

   //TRIGGER
   TBranch        *b_HLT_Mu15;
   TBranch        *b_HLT_Mu30;

   //GSF
   TBranch        *b_gsf_size;
   TBranch        *b_gsf_gsfet;
   TBranch        *b_gsf_px;
   TBranch        *b_gsf_py;
   TBranch        *b_gsf_pz;
   TBranch        *b_gsf_pt;
   TBranch        *b_gsf_eta;
   TBranch        *b_gsf_theta;
   TBranch        *b_gsf_phi;
   TBranch        *b_gsf_dz;
   TBranch        *b_gsf_isecaldriven;
   TBranch        *b_gsf_istrackerdriven;
   TBranch        *b_gsf_isEB;
   TBranch        *b_gsf_isEE;
   TBranch        *b_gsf_charge;
   TBranch        *b_gsf_deltaeta;
   TBranch        *b_gsf_deltaphi;
   TBranch        *b_gsf_e1x5overe5x5;
   TBranch        *b_gsf_e2x5overe5x5;
   TBranch        *b_gsf_sigmaetaeta;
   TBranch        *b_gsf_sigmaIetaIeta;
   TBranch        *b_gsf_hovere;
   TBranch        *b_gsf_eOVERp;
   TBranch        *b_gsf_vz;
   TBranch        *b_gsf_nHits;
   TBranch        *b_gsf_nLostInnerHits;
   TBranch        *b_gsf_fBrem;
   TBranch        *b_gsf_ecaliso;
   TBranch        *b_gsf_hcaliso1;
   TBranch        *b_gsf_hcaliso2;
   TBranch        *b_gsf_trackiso;

   TBranch        *b_gsfsc_px;
   TBranch        *b_gsfsc_py;
   TBranch        *b_gsfsc_pt;
   TBranch        *b_gsfsc_eta;
   TBranch        *b_gsfsc_phi;

   TBranch        *b_gsfpass_ID;
   TBranch        *b_gsfpass_ISO;
   TBranch        *b_gsfpass_HEEP;

   //MUONS
   TBranch        *b_muon_size;
   TBranch        *b_muon_pt;
   TBranch        *b_muon_eta;
   TBranch        *b_muon_phi;
   TBranch        *b_muon_theta;
   TBranch        *b_muon_ptError;
   TBranch        *b_muon_etaError;
   TBranch        *b_muon_phiError;
   TBranch        *b_muon_thetaError;
   TBranch        *b_muon_outerPt;
   TBranch        *b_muon_outerEta;
   TBranch        *b_muon_outerPhi;
   TBranch        *b_muon_outerTheta;
   TBranch        *b_muon_px;
   TBranch        *b_muon_py;
   TBranch        *b_muon_pz;
   TBranch        *b_muon_charge;
   TBranch        *b_muon_nhitstrack;
   TBranch        *b_muon_nhitspixel;
   TBranch        *b_muon_nhitstotal;
   TBranch        *b_muon_nhitsmuons;
   TBranch        *b_muon_nlayerswithhits;
   TBranch        *b_muon_nSegmentMatch;
   TBranch        *b_muon_isTrackerMuon;
   TBranch        *b_muon_chi2;
   TBranch        *b_muon_ndof;
   TBranch        *b_muon_normChi2;
   TBranch        *b_muon_d0;
   TBranch        *b_muon_dz_cmsCenter;
   TBranch        *b_muon_dz_beamSpot;
   TBranch        *b_muon_dz_firstPVtx;
   TBranch        *b_muon_dxy_cmsCenter;
   TBranch        *b_muon_dxy_beamSpot;
   TBranch        *b_muon_dxy_firstPVtx;
   TBranch        *b_muon_d0Error;
   TBranch        *b_muon_dzError;
   TBranch        *b_muon_dxyError;
   TBranch        *b_muon_trackIso03;
   TBranch        *b_muon_trackIso05;
   TBranch        *b_muon_emIso03;
   TBranch        *b_muon_emIso05;
   TBranch        *b_muon_hadIso03;
   TBranch        *b_muon_hadIso05;
   TBranch        *b_muon_trackIso03_ptInVeto;
   TBranch        *b_muon_trackIso05_ptInVeto;
   TBranch        *b_muon_emIso03_ptInVeto;
   TBranch        *b_muon_emIso05_ptInVeto;
   TBranch        *b_muon_hadIso03_ptInVeto;
   TBranch        *b_muon_hadIso05_ptInVeto;
   TBranch        *b_muon_innerPosx;
   TBranch        *b_muon_innerPosy;
   TBranch        *b_muon_innerPosz;

   unsigned int dataTrig[3] = {0, 0, 0};
   unsigned int dataEntries = 0;

   // primary vertex information /////////////////////////////////////////////
   vector<vector<float> > PVX_ScalingFactors;
   if (usePUInfo) {
      stringstream ssPUInfile;
      ssPUInfile << "emu_PUinfo" << nPVtxMax << "PVtx_eleBar" << bar_et << "_eleEnd" << end_et << "_mu" << muon_et << "_" << LumiFactor << "pb-1.root";
      // test if the PU file exists already. If yes -> use it; if no -> create it and fill it
      ifstream iFile(ssPUInfile.str().c_str());
      if (iFile) {
         iFile.close();
         TFile *inputPV = new TFile(ssPUInfile.str().c_str(), "open");
         inputPV->cd();

         vector<TH1F *> copy_emuLoose_dataOverX_nValidPv;
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverTtbar_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverZtautau_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverWw_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverWz_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverTw_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverWjets_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverZmumu_nValidPv"));
         copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverZee_nValidPv"));
         //copy_emuLoose_dataOverX_nValidPv.push_back((TH1F*)inputPV->Get("emuLoose_dataOverZz_nValidPv"));

         for (unsigned int p = 0; p < copy_emuLoose_dataOverX_nValidPv.size(); ++p) {
            vector<float> PVX_ScalingFactor;
            for (unsigned int i = 0; i < nPVtxMax; ++i) {
               //PVX_ScalingFactor.push_back(copy_emuLoose_dataOverX_nValidPv.at(p)->GetBinContent(i + 1) * nPVtxMax);
               PVX_ScalingFactor.push_back(copy_emuLoose_dataOverX_nValidPv.at(p)->GetBinContent(i + 1));
            }
            PVX_ScalingFactors.push_back(PVX_ScalingFactor);
         }
         cout << "Pile up information from file: " << ssPUInfile.str() << endl;
      } else {
         cout << "No file with pile up information for this parameter set found. Will create it." << endl;
         generatePUFile = true;
      }
   } else {
      PVX_ScalingFactors.clear();
      for (unsigned int p = 0; p < emuLoose_dataOverX_nValidPv.size(); ++p) {
         vector<float> PVX_ScalingFactor;
         for (unsigned int i = 0 ; i < nPVtxMax ; ++i) {
            PVX_ScalingFactor.push_back(1.);
         }
         PVX_ScalingFactors.push_back(PVX_ScalingFactor);
      }
   }
   if (usePUInfo && generatePUFile) {
      // 1st loop to get VERTEX information
      for (int p = 0; p < nbFile; ++p) {
         cout << "accessing file " << p + 1 << " for PU information: " << input[p].first->GetName() << endl;
         input[p].first->cd();

         // Get the TREE and connect the necessary variables
         TTree *thetree;
         thetree = (TTree*)(input[p].first)->Get("gsfcheckerjob/tree");

         //RUN ID
         thetree->SetBranchAddress("runnumber", &c_runnumber, &b_runnumber);
         //PRIM VTX
         thetree->SetBranchAddress("pvsize", &c_pvsize, &b_pvsize);
         thetree->SetBranchAddress("pv_ndof", &c_pv_ndof, &b_pv_ndof);
         thetree->SetBranchAddress("pv_nTracks", &c_pv_nTracks, &b_pv_nTracks);
         //GSF
         thetree->SetBranchAddress("gsf_size", &c_gsf_size, &b_gsf_size);
         thetree->SetBranchAddress("gsf_gsfet", &c_gsf_gsfet, &b_gsf_gsfet);
         thetree->SetBranchAddress("gsfsc_eta", &c_gsfsc_eta, &b_gsfsc_eta);
         thetree->SetBranchAddress("gsfsc_phi", &c_gsfsc_phi, &b_gsfsc_phi);
          //MUONS
         thetree->SetBranchAddress("muon_size", &c_muon_size, &b_muon_size);
         thetree->SetBranchAddress("muon_pt", &c_muon_pt, &b_muon_pt);
         thetree->SetBranchAddress("muon_eta", &c_muon_eta, &b_muon_eta);
         thetree->SetBranchAddress("muon_phi", &c_muon_eta, &b_muon_eta);

         Long64_t nentries = (*thetree).GetEntries();
         if (p == DATA) dataEntries = nentries;
         cout << nentries << " events" << endl;
         unsigned int trig[3] = {0, 0, 0}; // which trigger was used how often
         //LOOP OVER EVENTS
         for (unsigned int i = 0; i < nentries; ++i) {
            if (i % 50000 == 0) cout << i << endl;
            thetree->GetEntry(i);

            int prescale = 0;
            if (p == DATA && Trigger(input[p].first, i, prescale, dataTrig) < 1) continue;
            if (p != DATA && p != ZTT && p != ZMM && p != ZEE) {
            //if (p != DATA) {
               if (i < nentries * dataTrig[0] / dataEntries) {
                  if (Trigger(input[p].first, i, prescale, trig, 1) < 1) continue;
               } else if (i < nentries * (dataTrig[0] + dataTrig[1]) / dataEntries) {
                  if (Trigger(input[p].first, i, prescale, trig, 2) < 1) continue;
               } else {
                  if (Trigger(input[p].first, i, prescale, trig, 3) < 1) continue;
               }
            }

            //PRIMARY VTX COUNTING
            unsigned int n_pvValid = 0;
            for (int j = 0; j < c_pvsize; ++j) {
               if (c_pv_ndof[j] > 3 && c_pv_nTracks[j] > 3)
                  n_pvValid++;
            }
            if (n_pvValid < 1) continue;

            //CREATE VTX PONDERATION PLOT
            float gsfPtMaxB = 0.;
            float gsfPtMaxE = 0.;
            int gsfECAL = 0;
            float muPtMax = 0.;
            //LOOP OVER ELES
            for (int j = 0; j < c_gsf_size; ++j) {
               //CLEANING : FAKE ELES FROM MUONS
               bool fakeEle = false;
               for (int k = 0; k < c_muon_size; ++k) {
                  //if (c_muon_pt[k] < muon_et) continue;
                  float DeltaR = sqrt((c_gsf_eta[j] - c_muon_eta[k]) * (c_gsf_eta[j] - c_muon_eta[k]) + (c_gsf_phi[j] - c_muon_phi[k]) * (c_gsf_phi[j] - c_muon_phi[k]));
                  if (DeltaR < 0.1) {
                     fakeEle = true;
                     break;
                  }
               }
               if (fakeEle) continue;

               if ((fabs(c_gsfsc_eta[j]) < 1.442)  //BARREL
                   &&
                   c_gsf_gsfet[j] > gsfPtMaxB) {
                  gsfPtMaxB = c_gsf_gsfet[j];
                  gsfECAL = 1;
               }
               if ((fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5)  //ENDCAP
                   &&
                   c_gsf_gsfet[j] > gsfPtMaxE) {
                  gsfPtMaxE = c_gsf_gsfet[j];
                  gsfECAL = -1;
               }
            }
            //LOOP OVER MUS
            for (int j = 0; j < c_muon_size; ++j)
               if (fabs(c_muon_eta[j]) < muon_etaMax && c_muon_pt[j] > muPtMax) muPtMax = c_muon_pt[j];

            if (((gsfECAL == 1 && gsfPtMaxB > bar_et)  //BARREL
                 ||
                 (gsfECAL == -1 && gsfPtMaxE > end_et)) //ENDCAP
                && muPtMax > muon_et) {

               emuLoose_nValidPv.at(p)->Fill(n_pvValid);
            }
         } // end loop over events

         // calculate PU reweighting factors
         if (p == DATA) emuLoose_nValidPv.at(DATA)->Scale(1. / emuLoose_nValidPv.at(DATA)->Integral());
         if (p > DATA) {
            emuLoose_nValidPv.at(p)->Scale(1. / emuLoose_nValidPv.at(p)->Integral());
            emuLoose_dataOverX_nValidPv.at(p - 1)->Divide(emuLoose_nValidPv.at(DATA), emuLoose_nValidPv.at(p));
   
            cout << "Normalization factor PU reweighting: " << emuLoose_dataOverX_nValidPv.at(p - 1)->Integral() / nPVtxMax << endl;
         }
      } // end 1st loop over files

      stringstream ssPUOutfile;
      ssPUOutfile << "emu_PUinfo" << nPVtxMax << "PVtx_eleBar" << bar_et << "_eleEnd" << end_et << "_mu" << muon_et << "_" << LumiFactor << "pb-1.root";
      TFile *outputPV = new TFile(ssPUOutfile.str().c_str(), "recreate");
      outputPV->cd();

      emuLoose_nValidPv.at(DATA)->Write();
      for (unsigned int p = 0; p < emuLoose_dataOverX_nValidPv.size(); ++p) {
         emuLoose_nValidPv.at(p + 1)->Write();
         emuLoose_dataOverX_nValidPv.at(p)->Write();
      }
      outputPV->Close();

      PVX_ScalingFactors.clear();
      for (unsigned int p = 0; p < emuLoose_dataOverX_nValidPv.size(); ++p) {
         vector<float> PVX_ScalingFactor;
         for (unsigned int i = 0; i < nPVtxMax; ++i) {
            PVX_ScalingFactor.push_back(emuLoose_dataOverX_nValidPv.at(p)->GetBinContent(i + 1));
         }
         PVX_ScalingFactors.push_back(PVX_ScalingFactor);
      }
   }
   ///////////////////////////////////////////////////////////////////////////

   runs_HLT_Mu15_Photon20_CaloIdL.first = 99999999;
   runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.first = 99999999;
   runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.first = 99999999;
   dataTrig[0] = 0;
   dataTrig[1] = 0;
   dataTrig[2] = 0;
   dataEntries = 0;
   // 2nd loop for analysis 
   // GETTING FILES
   for (int p = 0; p < nbFile; ++p) {
      // correction for trigger
      if (p == ZTT || p == ZMM || p == ZEE) MCemuScaleFactor = MCScaleFactor * Elec_trigger;
      //if (0) MCemuScaleFactor = MCScaleFactor * Elec_trigger;
      else MCemuScaleFactor = MCScaleFactor;

      cout << "accessing file " << p + 1 << ": " << input[p].first->GetName() << endl;
      if (p > DATA) cout << "MC correction factor: " << MCemuScaleFactor << endl;
      input[p].first->cd();

      // Get the TREE
      TTree *thetree;
      thetree = (TTree*)(input[p].first)->Get("gsfcheckerjob/tree");

      //RUN ID
      thetree->SetBranchAddress("runnumber", &c_runnumber, &b_runnumber);
      thetree->SetBranchAddress("eventnumber", &c_eventnumber, &b_eventnumber);
      thetree->SetBranchAddress("luminosityBlock", &c_luminosityBlock, &b_luminosityBlock);

      //GLOBAL
      thetree->SetBranchAddress("nJetsAKT_pt15",&c_nJetsAKT_pt15,&b_nJetsAKT_pt15);
      thetree->SetBranchAddress("calomet", &c_calomet, &b_calomet);
      thetree->SetBranchAddress("met", &c_met, &b_met);
      thetree->SetBranchAddress("pthat", &c_pthat, &b_pthat);
      thetree->SetBranchAddress("bsposx", &c_bsposx, &b_bsposx);
      thetree->SetBranchAddress("bsposy", &c_bsposy, &b_bsposy);
      thetree->SetBranchAddress("bsposz", &c_bsposz, &b_bsposz);

      //JETS AKT
      thetree->SetBranchAddress("jetAKT_size",&c_jetAKT_size,&b_jetAKT_size);
      thetree->SetBranchAddress("jetAKT_pt",&c_jetAKT_pt,&b_jetAKT_pt);
      thetree->SetBranchAddress("jetAKT_eta",&c_jetAKT_eta,&b_jetAKT_eta);
      thetree->SetBranchAddress("jetAKT_phi",&c_jetAKT_phi,&b_jetAKT_phi);
      thetree->SetBranchAddress("jetAKT_em",&c_jetAKT_em,&b_jetAKT_em);

      //PRIM VTX
      thetree->SetBranchAddress("pvsize", &c_pvsize, &b_pvsize);
      thetree->SetBranchAddress("pvz", &c_pvz, &b_pvz);
      thetree->SetBranchAddress("pv_isValid", &c_pv_isValid, &b_pv_isValid);
      thetree->SetBranchAddress("pv_ndof", &c_pv_ndof, &b_pv_ndof);
      thetree->SetBranchAddress("pv_nTracks", &c_pv_nTracks, &b_pv_nTracks);
      thetree->SetBranchAddress("pv_normChi2", &c_pv_normChi2, &b_pv_normChi2);
      thetree->SetBranchAddress("pv_totTrackSize", &c_pv_totTrackSize, &b_pv_totTrackSize);

      //TRIGGER
      thetree->SetBranchAddress("HLT_Mu15", &c_HLT_Mu15, &b_HLT_Mu15);
      thetree->SetBranchAddress("HLT_Mu30", &c_HLT_Mu30, &b_HLT_Mu30);

      //GSF
      thetree->SetBranchAddress("gsf_size", &c_gsf_size, &b_gsf_size);
      thetree->SetBranchAddress("gsf_gsfet", &c_gsf_gsfet, &b_gsf_gsfet);
      thetree->SetBranchAddress("gsf_px", &c_gsf_px, &b_gsf_px);
      thetree->SetBranchAddress("gsf_py", &c_gsf_py, &b_gsf_py);
      thetree->SetBranchAddress("gsf_pz", &c_gsf_pz, &b_gsf_pz);
      thetree->SetBranchAddress("gsf_pt", &c_gsf_pt, &b_gsf_pt);
      thetree->SetBranchAddress("gsf_eta", &c_gsf_eta, &b_gsf_eta);
      thetree->SetBranchAddress("gsf_theta", &c_gsf_theta, &b_gsf_theta);
      thetree->SetBranchAddress("gsf_phi", &c_gsf_phi, &b_gsf_phi);
      thetree->SetBranchAddress("gsf_dz", &c_gsf_dz, &b_gsf_dz);
      thetree->SetBranchAddress("gsf_isecaldriven", &c_gsf_isecaldriven, &b_gsf_isecaldriven);
      thetree->SetBranchAddress("gsf_istrackerdriven", &c_gsf_istrackerdriven, &b_gsf_istrackerdriven);
      thetree->SetBranchAddress("gsf_isEB", &c_gsf_isEB, &b_gsf_isEB);
      thetree->SetBranchAddress("gsf_isEE", &c_gsf_isEE, &b_gsf_isEE);
      thetree->SetBranchAddress("gsf_charge", &c_gsf_charge, &b_gsf_charge);
      thetree->SetBranchAddress("gsf_deltaeta", &c_gsf_deltaeta, &b_gsf_deltaeta);
      thetree->SetBranchAddress("gsf_deltaphi", &c_gsf_deltaphi, &b_gsf_deltaphi);
      thetree->SetBranchAddress("gsf_e1x5overe5x5", &c_gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
      thetree->SetBranchAddress("gsf_e2x5overe5x5", &c_gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
      thetree->SetBranchAddress("gsf_sigmaetaeta", &c_gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
      thetree->SetBranchAddress("gsf_sigmaIetaIeta", &c_gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
      thetree->SetBranchAddress("gsf_hovere", &c_gsf_hovere, &b_gsf_hovere);
      thetree->SetBranchAddress("gsf_eOVERp", &c_gsf_eOVERp, &b_gsf_eOVERp);
      thetree->SetBranchAddress("gsf_vz", &c_gsf_vz, &b_gsf_vz);
      thetree->SetBranchAddress("gsf_nHits", &c_gsf_nHits, &b_gsf_nHits);
      thetree->SetBranchAddress("gsf_nLostInnerHits", &c_gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
      thetree->SetBranchAddress("gsf_fBrem", &c_gsf_fBrem, &b_gsf_fBrem);
      thetree->SetBranchAddress("gsf_ecaliso", &c_gsf_ecaliso, &b_gsf_ecaliso);
      thetree->SetBranchAddress("gsf_hcaliso1", &c_gsf_hcaliso1, &b_gsf_hcaliso1);
      thetree->SetBranchAddress("gsf_hcaliso2", &c_gsf_hcaliso2, &b_gsf_hcaliso2);
      thetree->SetBranchAddress("gsf_trackiso", &c_gsf_trackiso, &b_gsf_trackiso);
      thetree->SetBranchAddress("gsfsc_px", &c_gsfsc_px, &b_gsfsc_px);
      thetree->SetBranchAddress("gsfsc_py", &c_gsfsc_py, &b_gsfsc_py);
      thetree->SetBranchAddress("gsfsc_pt", &c_gsfsc_pt, &b_gsfsc_pt);
      thetree->SetBranchAddress("gsfsc_eta", &c_gsfsc_eta, &b_gsfsc_eta);
      thetree->SetBranchAddress("gsfsc_phi", &c_gsfsc_phi, &b_gsfsc_phi);
      thetree->SetBranchAddress("gsfpass_ID", &c_gsfpass_ID, &b_gsfpass_ID);
      thetree->SetBranchAddress("gsfpass_ISO", &c_gsfpass_ISO, &b_gsfpass_ISO);
      thetree->SetBranchAddress("gsfpass_HEEP", &c_gsfpass_HEEP, &b_gsfpass_HEEP);

      //MUONS
      thetree->SetBranchAddress("muon_size", &c_muon_size, &b_muon_size);
      thetree->SetBranchAddress("muon_pt", &c_muon_pt, &b_muon_pt);
      thetree->SetBranchAddress("muon_eta", &c_muon_eta, &b_muon_eta);
      thetree->SetBranchAddress("muon_phi", &c_muon_phi, &b_muon_phi);
      thetree->SetBranchAddress("muon_theta", &c_muon_theta, &b_muon_theta);
      thetree->SetBranchAddress("muon_ptError", &c_muon_ptError, &b_muon_ptError);
      thetree->SetBranchAddress("muon_etaError", &c_muon_etaError, &b_muon_etaError);
      thetree->SetBranchAddress("muon_phiError", &c_muon_phiError, &b_muon_phiError);
      thetree->SetBranchAddress("muon_thetaError", &c_muon_thetaError, &b_muon_thetaError);
      thetree->SetBranchAddress("muon_outerPt", &c_muon_outerPt, &b_muon_outerPt);
      thetree->SetBranchAddress("muon_outerEta", &c_muon_outerEta, &b_muon_outerEta);
      thetree->SetBranchAddress("muon_outerPhi", &c_muon_outerPhi, &b_muon_outerPhi);
      thetree->SetBranchAddress("muon_outerTheta", &c_muon_outerTheta, &b_muon_outerTheta);
      thetree->SetBranchAddress("muon_px", &c_muon_px, &b_muon_px);
      thetree->SetBranchAddress("muon_py", &c_muon_py, &b_muon_py);
      thetree->SetBranchAddress("muon_pz", &c_muon_pz, &b_muon_pz);
      thetree->SetBranchAddress("muon_charge", &c_muon_charge, &b_muon_charge);
      thetree->SetBranchAddress("muon_nhitstrack", &c_muon_nhitstrack, &b_muon_nhitstrack);
      thetree->SetBranchAddress("muon_nhitspixel", &c_muon_nhitspixel, &b_muon_nhitspixel);
      thetree->SetBranchAddress("muon_nhitstotal", &c_muon_nhitstotal, &b_muon_nhitstotal);
      thetree->SetBranchAddress("muon_nhitsmuons", &c_muon_nhitsmuons, &b_muon_nhitsmuons);
      thetree->SetBranchAddress("muon_nlayerswithhits", &c_muon_nlayerswithhits, &b_muon_nlayerswithhits);
      thetree->SetBranchAddress("muon_nSegmentMatch", &c_muon_nSegmentMatch, &b_muon_nSegmentMatch);
      thetree->SetBranchAddress("muon_isTrackerMuon", &c_muon_isTrackerMuon, &b_muon_isTrackerMuon);
      thetree->SetBranchAddress("muon_chi2", &c_muon_chi2, &b_muon_chi2);
      thetree->SetBranchAddress("muon_ndof", &c_muon_ndof, &b_muon_ndof);
      thetree->SetBranchAddress("muon_normChi2", &c_muon_normChi2, &b_muon_normChi2);
      thetree->SetBranchAddress("muon_d0", &c_muon_d0, &b_muon_d0);
      thetree->SetBranchAddress("muon_dz_cmsCenter", &c_muon_dz_cmsCenter, &b_muon_dz_cmsCenter);
      thetree->SetBranchAddress("muon_dz_beamSpot", &c_muon_dz_beamSpot, &b_muon_dz_beamSpot);
      thetree->SetBranchAddress("muon_dz_firstPVtx", &c_muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
      thetree->SetBranchAddress("muon_dxy_cmsCenter", &c_muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
      thetree->SetBranchAddress("muon_dxy_beamSpot", &c_muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
      thetree->SetBranchAddress("muon_dxy_firstPVtx", &c_muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
      thetree->SetBranchAddress("muon_d0Error", &c_muon_d0Error, &b_muon_d0Error);
      thetree->SetBranchAddress("muon_dzError", &c_muon_dzError, &b_muon_dzError);
      thetree->SetBranchAddress("muon_dxyError", &c_muon_dxyError, &b_muon_dxyError);
      thetree->SetBranchAddress("muon_trackIso03", &c_muon_trackIso03, &b_muon_trackIso03);
      thetree->SetBranchAddress("muon_trackIso05", &c_muon_trackIso05, &b_muon_trackIso05);
      thetree->SetBranchAddress("muon_emIso03", &c_muon_emIso03, &b_muon_emIso03);
      thetree->SetBranchAddress("muon_emIso05", &c_muon_emIso05, &b_muon_emIso05);
      thetree->SetBranchAddress("muon_hadIso03", &c_muon_hadIso03, &b_muon_hadIso03);
      thetree->SetBranchAddress("muon_hadIso05", &c_muon_hadIso05, &b_muon_hadIso05);
      thetree->SetBranchAddress("muon_trackIso03_ptInVeto", &c_muon_trackIso03_ptInVeto, &b_muon_trackIso03_ptInVeto);
      thetree->SetBranchAddress("muon_trackIso05_ptInVeto", &c_muon_trackIso05_ptInVeto, &b_muon_trackIso05_ptInVeto);
      thetree->SetBranchAddress("muon_emIso03_ptInVeto", &c_muon_emIso03_ptInVeto, &b_muon_emIso03_ptInVeto);
      thetree->SetBranchAddress("muon_emIso05_ptInVeto", &c_muon_emIso05_ptInVeto, &b_muon_emIso05_ptInVeto);
      thetree->SetBranchAddress("muon_hadIso03_ptInVeto", &c_muon_hadIso03_ptInVeto, &b_muon_hadIso03_ptInVeto);
      thetree->SetBranchAddress("muon_hadIso05_ptInVeto", &c_muon_hadIso05_ptInVeto, &b_muon_hadIso05_ptInVeto);
      thetree->SetBranchAddress("muon_innerPosx", &c_muon_innerPosx, &b_muon_innerPosx);
      thetree->SetBranchAddress("muon_innerPosy", &c_muon_innerPosy, &b_muon_innerPosy);
      thetree->SetBranchAddress("muon_innerPosz", &c_muon_innerPosz, &b_muon_innerPosz);

      // set up histograms
      vector<TH1F *> helper;
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMass_" + sign[k] + suffix[p], "emuMass_" + sign[k] + suffix[p], 150, 0., 1500.));
      emuMass.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMassEB_" + sign[k] + suffix[p], "emuMassEB_" + sign[k] + suffix[p], 150, 0., 1500.));
      emuMassEB.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMassEE_" + sign[k] + suffix[p], "emuMassEE_" + sign[k] + suffix[p], 150, 0., 1500.));
      emuMassEE.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emu_mass_accVSgood_" + sign[k] + suffix[p], "emu_mass_accVSgood_" + sign[k] + suffix[p], 150, 0., 1500.));
      emu_mass_accVSgood.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("MET_" + sign[k] + suffix[p], "MET_" + sign[k] + suffix[p], 50, 0., 500.));
      met.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("nVtx_" + sign[k] + suffix[p], "nVtx_" + sign[k] + suffix[p], nPVtxMax, 0., nPVtxMax));
      nVtx.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("dPhi_" + sign[k] + suffix[p], "dPhi_" + sign[k] + suffix[p], 64, 0., 3.2));
      dPhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("elePt_" + sign[k] + suffix[p], "elePt_" + sign[k] + suffix[p], 50, 0., 500.));
      elePt.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleEta_" + sign[k] + suffix[p], "eleEta_" + sign[k] + suffix[p], 25, -2.5, 2.5));
      eleEta.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("elePhi_" + sign[k] + suffix[p], "elePhi_" + sign[k] + suffix[p], 25, -3.14, 3.14));
      elePhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleId1_" + sign[k] + suffix[p], "eleId1_" + sign[k] + suffix[p], 50, -0.008, 0.008));
      eleId1.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleId2_" + sign[k] + suffix[p], "eleId2_" + sign[k] + suffix[p], 50, -0.1, 0.1));
      eleId2.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleId3_" + sign[k] + suffix[p], "eleId3_" + sign[k] + suffix[p], 50, 0., 0.04));
      eleId3.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleIso1_" + sign[k] + suffix[p], "eleIso1_" + sign[k] + suffix[p], 50, 0., 10.));
      eleIso1.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleIso2_" + sign[k] + suffix[p], "eleIso2_" + sign[k] + suffix[p], 50, 0., 10.));
      eleIso2.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleIso3_" + sign[k] + suffix[p], "eleIso3_" + sign[k] + suffix[p], 50, 0., 20.));
      eleIso3.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muIsoCombRel_" + sign[k] + suffix[p], "muIsoCombRel_" + sign[k] + suffix[p], 50, 0., 0.2));
      muIsoCombRel.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPtEleOPtMu_" + sign[k] + suffix[p], "muPtEleOPtMu_" + sign[k] + suffix[p], 50, 0., 4.0));
      muPtEleOPtMu.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPtPlusOPtMinus_" + sign[k] + suffix[p], "muPtPlusOPtMinus_" + sign[k] + suffix[p], 50, 0., 4.0));
      muPtPlusOPtMinus.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPt_" + sign[k] + suffix[p], "muPt_" + sign[k] + suffix[p], 50, 0., 500.));
      muPt.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muEta_" + sign[k] + suffix[p], "muEta_" + sign[k] + suffix[p], 25, -2.5, 2.5));
      muEta.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPhi_" + sign[k] + suffix[p], "muPhi_" + sign[k] + suffix[p], 25, -3.14, 3.14));
      muPhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muId1_" + sign[k] + suffix[p], "muId1_" + sign[k] + suffix[p], 50, 0., 10.));
      muId1.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muId2_" + sign[k] + suffix[p], "muId2_" + sign[k] + suffix[p], 32, 0., 32));
      muId2.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muId3_" + sign[k] + suffix[p], "muId3_" + sign[k] + suffix[p], 52, 0., 52));
      muId3.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muIso1_" + sign[k] + suffix[p], "muIso1_" + sign[k] + suffix[p], 50, 0., 10.));
      muIso1.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muIso2_" + sign[k] + suffix[p], "muIso2_" + sign[k] + suffix[p], 50, 0., 10.));
      muIso2.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muIso3_" + sign[k] + suffix[p], "muIso3_" + sign[k] + suffix[p], 50, 0., 20.));
      muIso3.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("numOfJets_" + sign[k] + suffix[p], "numOfJets_" + sign[k] + suffix[p], 20, 0., 20.));
      numOfJets.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("numOfJetsPt15_" + sign[k] + suffix[p], "numOfJetsPt15_" + sign[k] + suffix[p], 20, 0., 20.));
      numOfJetsPt15.push_back(helper);
      helper.clear();

      Long64_t nentries = (*thetree).GetEntries();
      cout << nentries << " events" << endl;
      if (p == DATA) {
         dataEntries = nentries;
         cout << "-----------------------------------------------------------------------------------------------------------" << endl;
         cout << "M_emu > 600GeV/c^2       |  run   | lumi |    event   |   M_emu  |"
              << " muon pt |"
              << " muon eta |"
              //<< "  muon phi  |"
              //<< "  muon trackIso3  |"
              //<< "  muon emIso03  |"
              //<< "  muon hadIso03  |"
              //<< "  muon VETOtrackIso03  |"
              //<< "  muon VETOemIso03  |"
              //<< "  muon VETOhadIso03  |"
              //<< "  muon normChi2  |"
              //<< "  muon dxy_beamSpot  |"
              //<< "  muon nhitstrack  |"
              //<< "  muon nhitsmuons  |"
              << " ele pt  |"
              << "  ele eta |"
              //<< "  ele phi  |"
              //<< "  ele trackIso  |"
              //<< "  ele ecalIso  |"
              //<< "  ele hcalIso1  |"
              //<< "  ele hcalIso2  |"
              << endl;;
         cout << "-----------------------------------------------------------------------------------------------------------" << endl;

      }
      unsigned int trig[3] = {0, 0, 0};
      //LOOP OVER EVENTS
      //for (unsigned int i = 0; i < 10000; ++i) {
      for (unsigned int i = 0; i < nentries; ++i) {
         if (i % 50000 == 0) cout << "Processing event " << i << endl;
         thetree->GetEntry(i);

         int prescale = 0;
         if (p == DATA && Trigger(input[p].first, i, prescale, dataTrig) < 1) continue;
         if (p != DATA && p != ZTT && p != ZMM && p != ZEE) {
         //if (p != DATA) {
            if (i < nentries * dataTrig[0] / dataEntries) {
               if (Trigger(input[p].first, i, prescale, trig, 1) < 1) continue;
            } else if (i < nentries * (dataTrig[0] + dataTrig[1]) / dataEntries) {
               if (Trigger(input[p].first, i, prescale, trig, 2) < 1) continue;
            } else {
               if (Trigger(input[p].first, i, prescale, trig, 3) < 1) continue;
            }
         }

         vector<int> GSF_passHEEP;
         vector<int> GSF_passACC;

         vector<int> MU_passGOOD;
         vector<int> MU_passACC;

         //PRIMARY VTX COUNTING
         unsigned int n_pvValid = 0;
         for (int j = 0; j < c_pvsize; ++j) {
            if (c_pv_ndof[j] > 3 && c_pv_nTracks[j] > 3)
               n_pvValid++;
         }
         if (n_pvValid < 1) continue;

         //FILL THE VTX WEIGTH
         float npv_weight = 1.;
         if (p >= TTBAR) {
           if (n_pvValid < nPVtxMax) npv_weight = PVX_ScalingFactors.at(p - 1).at(n_pvValid);
           else cout << "Event has more than " << nPVtxMax << " vertices. Number of valid primary vertices: " << n_pvValid << endl;
         }

         //LOOP OVER ELES
         for (int j = 0; j < c_gsf_size; ++j) {
            //CLEANING : FAKE ELES FROM MUONS
            bool fakeEle = false;
            for (int k = 0; k < c_muon_size; ++k) {
               //if (c_muon_pt[k] < muon_et) continue;
               float DeltaR = sqrt((c_gsf_eta[j] - c_muon_eta[k]) * (c_gsf_eta[j] - c_muon_eta[k]) + (c_gsf_phi[j] - c_muon_phi[k]) * (c_gsf_phi[j] - c_muon_phi[k]));
               if (DeltaR < 0.1) {
                  fakeEle = true;
                  break;
               }
            }
            if (fakeEle) continue;

            //BARREL HEEP
            if (fabs(c_gsfsc_eta[j]) < 1.442
                && c_gsf_gsfet[j] > bar_et
                && c_gsf_isecaldriven[j]
                && c_gsf_hovere[j] < bar_hoE
                && fabs(c_gsf_deltaeta[j]) < bar_DEta
                && fabs(c_gsf_deltaphi[j]) < bar_DPhi
                && (c_gsf_e2x5overe5x5[j] > bar_e2x5e5x5 || c_gsf_e1x5overe5x5[j] > bar_e1x5e5x5)
                && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * c_gsf_gsfet[j])
                && c_gsf_trackiso[j] < bar_isoTrack
                && c_gsf_nLostInnerHits[j] <= bar_missInnerHits
               ) GSF_passHEEP.push_back(j);

            //ENDCAP HEEP
            if ((fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5)
                && c_gsf_gsfet[j] > end_et
                && c_gsf_isecaldriven[j]
                && c_gsf_hovere[j] < end_hoE
                && fabs(c_gsf_deltaeta[j]) < end_DEta
                && fabs(c_gsf_deltaphi[j]) < end_DPhi
                && (c_gsf_e2x5overe5x5[j] > end_e2x5e5x5 || c_gsf_e1x5overe5x5[j] > end_e1x5e5x5)
                && c_gsf_sigmaIetaIeta[j] < end_sigmaietaieta
                && ((c_gsf_gsfet[j] < 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < end_isoEcalHcal1_1)
                    ||
                    (c_gsf_gsfet[j] >= 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (end_isoEcalHcal1_1 + end_isoEcalHcal1_2 * (c_gsf_gsfet[j] - 50.))))
                && c_gsf_trackiso[j] < end_isoTrack
                && c_gsf_nLostInnerHits[j] <= end_missInnerHits
               ) GSF_passHEEP.push_back(j);

            //GSF IN ACCEPTANCE
            if (c_gsf_gsfet[j] > bar_et
                && (fabs(c_gsfsc_eta[j]) < 1.442 || (fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5))
               ) GSF_passACC.push_back(j);
         }

         //LOOP OVER MUS
         for (int j = 0; j < c_muon_size; ++j) {
            //MU PASS GOOD
            if (c_muon_pt[j] > muon_et
                && fabs(c_muon_eta[j]) < muon_etaMax
                && c_muon_nhitstrack[j] >= muon_nHitsMinGlobal
                && c_muon_nhitspixel[j] >= muon_nHitsMinPixel
                && c_muon_nhitsmuons[j] >= muon_nHitsMinMuon
                && c_muon_nlayerswithhits[j] >= muon_nLayersMin
                && fabs(c_muon_dxy_beamSpot[j]) < muon_impactParamMax
                && c_muon_nSegmentMatch[j] >= muon_nSegMatchMin
                && c_muon_isTrackerMuon[j]
                && c_muon_trackIso03[j] / c_muon_pt[j] < muon_relIsoCutMax
               ) MU_passGOOD.push_back(j);

            //MU PASS ACC
            if (c_muon_pt[j] > muon_et
                && fabs(c_muon_eta[j]) < muon_etaMax
               ) MU_passACC.push_back(j);
         }

         //search highest pt passing muon if there are more than one
         unsigned int MU_leadingPassGOOD = 0;
         if (MU_passGOOD.size() > 1) {
            for (unsigned int j = 0; j < MU_passGOOD.size(); ++j)
               if (c_muon_pt[MU_passGOOD[j]] > c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]) MU_leadingPassGOOD = j;
         }

         //HEEP - MU_GOOD
         if (GSF_passHEEP.size() > 0 && MU_passGOOD.size() > 0) {
            // remove duplicates
            // TODO implement deltaR <0.1 correctly
            //if (fabs(c_gsf_phi[GSF_passHEEP[0]] - c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]]) < 0.1) continue;

            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(c_gsf_gsfet[GSF_passHEEP[0]], c_gsf_eta[GSF_passHEEP[0]], c_gsf_phi[GSF_passHEEP[0]], 0.000511);
            mu1.SetPtEtaPhiM(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], 0.10566);

            double invMass = (ele1 + mu1).M();

            // shift MC spectrum
            //if (p > DATA) invMass -= 5.;

            //MASS CUT
            if (invMass < minInvMass) continue;

            float CombRelIso = (c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]]) / c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];

            // fill histograms for all, LS and OS
            for (unsigned int k = 0; k< 3; ++k) {
               if (k == 1 && c_gsf_charge[GSF_passHEEP[0]] * c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) continue;
               if (k == 2 && c_gsf_charge[GSF_passHEEP[0]] * c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) continue;

               emuMass.at(p).at(k)->Fill(invMass, npv_weight);
               if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 1.442) emuMassEB.at(p).at(k)->Fill(invMass, npv_weight);
               if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 2.5) emuMassEE.at(p).at(k)->Fill(invMass, npv_weight);

               // fill test histograms
               met.at(p).at(k)->Fill(c_calomet, npv_weight);
               nVtx.at(p).at(k)->Fill(n_pvValid, npv_weight);
               if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) < 3.14) dPhi.at(p).at(k)->Fill(fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), npv_weight);
               if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) > 3.14) dPhi.at(p).at(k)->Fill(6.28 - fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), npv_weight);
               elePt.at(p).at(k)->Fill(c_gsf_gsfet[GSF_passHEEP[0]], npv_weight);
               eleEta.at(p).at(k)->Fill(c_gsf_eta[GSF_passHEEP[0]], npv_weight);
               elePhi.at(p).at(k)->Fill(c_gsf_phi[GSF_passHEEP[0]], npv_weight);
               eleId1.at(p).at(k)->Fill(c_gsf_deltaeta[GSF_passHEEP[0]], npv_weight);
               eleId2.at(p).at(k)->Fill(c_gsf_deltaphi[GSF_passHEEP[0]], npv_weight);
               eleId3.at(p).at(k)->Fill(c_gsf_sigmaIetaIeta[GSF_passHEEP[0]], npv_weight);
               eleIso1.at(p).at(k)->Fill(c_gsf_ecaliso[GSF_passHEEP[0]], npv_weight);
               eleIso2.at(p).at(k)->Fill(c_gsf_hcaliso1[GSF_passHEEP[0]] + c_gsf_hcaliso2[GSF_passHEEP[0]], npv_weight);
               eleIso3.at(p).at(k)->Fill(c_gsf_trackiso[GSF_passHEEP[0]], npv_weight);

               muIsoCombRel.at(p).at(k)->Fill(CombRelIso, npv_weight);
               muPtEleOPtMu.at(p).at(k)->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) muPtPlusOPtMinus.at(p).at(k)->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) muPtPlusOPtMinus.at(p).at(k)->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]/c_gsf_gsfet[GSF_passHEEP[0]], npv_weight);
               muPt.at(p).at(k)->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muEta.at(p).at(k)->Fill(c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muPhi.at(p).at(k)->Fill(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muId1.at(p).at(k)->Fill(c_muon_normChi2[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muId2.at(p).at(k)->Fill(c_muon_nhitstrack[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muId3.at(p).at(k)->Fill(c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muIso1.at(p).at(k)->Fill(c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muIso2.at(p).at(k)->Fill(c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);
               muIso3.at(p).at(k)->Fill(c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]], npv_weight);

               numOfJets.at(p).at(k)->Fill(c_jetAKT_size, npv_weight);
               numOfJetsPt15.at(p).at(k)->Fill(c_nJetsAKT_pt15, npv_weight);

               if (p == DATA) {
               //if (p == DATA && c_HLT_Mu15) {
               //if (p == DATA && c_HLT_Mu30) {
                  // fill tree with good events
                  emuInvMass = invMass;
                  if (k == ALL) emuDataTree->Fill();
                  if (k == LS) emuLSDataTree->Fill();
                  if (k == OS) emuOSDataTree->Fill();
               }

            }

            //LIKE SIGN OPPOSITE SIGN
            if (p == DATA) {
               if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) nb_plus_plus++;
               if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) nb_plus_minus++;
               if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) nb_minus_plus++;
               if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) nb_minus_minus++;

               if (invMass > 600. || invMass < 1.) {
                  if (invMass > 600.) cout << "M_emu > 600GeV/c^2 event | " << setw(6) << c_runnumber << " | " << setw(4) << c_luminosityBlock << " | " << setw(10) << c_eventnumber << " | "; 
                  if (invMass < 1.) cout << "M_emu < 1GeV/c^2 event   | " << setw(6) << c_runnumber << " | " << setw(4) << c_luminosityBlock << " | " << setw(10) << c_eventnumber << " | "; 
                       cout << setw(8) << invMass << " | " << 
                       setw(7) << c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                       setw(8) << c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_trackIso03_ptInVeto[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_emIso03_ptInVeto[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_hadIso03_ptInVeto[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_normChi2[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_dxy_beamSpot[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_nhitstrack[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //    setw(8) << c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]] << " | " <<
                   //                                          " << " | " <<
                       setw(7) << c_gsf_gsfet[GSF_passHEEP[0]] << " | " << 
                       setw(8) << c_gsfsc_eta[GSF_passHEEP[0]] << " | " <<
                   //    setw(8) << c_gsfsc_phi[GSF_passHEEP[0]] << " | " <<
                   //    setw(8) << c_gsf_trackiso[GSF_passHEEP[0]] << " | " <<
                   //    setw(8) << c_gsf_ecaliso[GSF_passHEEP[0]] << " | " <<
                   //    setw(8) << c_gsf_hcaliso1[GSF_passHEEP[0]] << " | " <<
                   //    setw(8) << c_gsf_hcaliso2[GSF_passHEEP[0]] << " | " <<
                       endl;
               }
            }
         }

         //GSF ACC - MU_GOOD
         if (GSF_passACC.size() > 0 && MU_passGOOD.size() > 0) {

            //REMOVE DUPLICATES
            //if (fabs(c_gsf_phi[GSF_passACC[0]] - c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]]) < 0.1) continue;

            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(c_gsf_gsfet[GSF_passACC[0]], c_gsf_eta[GSF_passACC[0]], c_gsf_phi[GSF_passACC[0]], 0.000511);
            mu1.SetPtEtaPhiM(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], 0.10566);

            double invMass = (ele1 + mu1).M();

            for (unsigned int k = 0; k < 3; ++k) {
               if (k == 1 && c_gsf_charge[GSF_passACC[0]] * c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) continue;
               if (k == 2 && c_gsf_charge[GSF_passACC[0]] * c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) continue;
               emu_mass_accVSgood.at(p).at(k)->Fill(invMass, npv_weight);
            }
         }

      } // END LOOP OVER ENTRIES

      if (p == DATA) {
         cout << "HLT_Mu15_Photon20_CaloIdL: " << dataTrig[0] << " , HLT_Mu8_Ele17_CaloIdT_CaloIsoVL: " << dataTrig[1] << " , HLT_Mu17_Ele8_CaloIdT_CaloIsoVL: " << dataTrig[2] << " , sum: " << dataTrig[0]+dataTrig[1]+dataTrig[2] << endl;
         cout << "Runrange HLT_Mu15_Photon20_CaloIdL: " << runs_HLT_Mu15_Photon20_CaloIdL.first << " - " << runs_HLT_Mu15_Photon20_CaloIdL.second << endl;
         cout << "Runrange HLT_Mu8_Ele17_CaloIdT_CaloIsoVL: " << runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.first << " - " << runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.second << endl;
         cout << "Runrange HLT_Mu17_Ele8_CaloIdT_CaloIsoVL: " << runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.first << " - " << runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.second << endl;
      }
      else cout << "HLT_Mu15_Photon20_CaloIdL: " << trig[0] << " , HLT_Mu8_Ele17_CaloIdT_CaloIsoVL: " << trig[1] << " , HLT_Mu17_Ele8_CaloIdT_CaloIsoVL: " << trig[2] << " , sum: " << trig[0]+trig[1]+trig[2] << endl;

      if (p == DATA) {
         // write root file with good emu event data
         goodEvFile->cd();
         emuDataTree->Write();
         emuLSDataTree->Write();
         emuOSDataTree->Write();
      }

      if (p > DATA) {
         for (unsigned int k = 0; k < 3; ++k) {
            //SCALE MC
            emuMass.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            emuMassEB.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            emuMassEE.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            emu_mass_accVSgood.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);

            met.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            nVtx.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            dPhi.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            elePt.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleEta.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            elePhi.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleId1.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleId2.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleId3.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleIso1.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleIso2.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            eleIso3.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);

            muIsoCombRel.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muPtEleOPtMu.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muPtPlusOPtMinus.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muPt.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muEta.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muPhi.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muId1.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muId2.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muId3.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muIso1.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muIso2.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            muIso3.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);

            numOfJets.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);
            numOfJetsPt15.at(p).at(k)->Scale(input[p].second * LumiFactor * MCemuScaleFactor);

            // set systematic errors
            //for (unsigned int i = 0; i < emuMass.at(p).at(k)->GetNbinsX() + 2; ++i) emuMass.at(p).at(k)->SetBinError(i, emuMass.at(p).at(k)->GetBinContent(i) * systErrMCLuEff[p - 1]);
            //cout << "Systematic error for file " << p << ": " << systErrMCLuEff[p - 1] << endl;
         }
      }

   }//END FILE LOOP

   // define groups of MC samples
   vector<bool> ttLikeSamples(5, true);
   vector<bool> contamSamples(5, false);
   contamSamples.push_back(true); // WJets
   contamSamples.push_back(true); // Zmm
   contamSamples.push_back(true); // Zee
   vector<bool> allSamples(QCD-1, true);

   //PRINT INTEGRAL

   cout << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   cout << "HEEP - TIGHT MU        Lumi        = " << LumiFactor << "pb-1" << endl;
   cout << "                       e pT EB     > " << bar_et << "GeV/c" << endl;
   cout << "                       e pT EE     > " << end_et << "GeV/c" << endl;
   cout << "                       mu pT       > " << muon_et << "GeV/c" << endl;
   cout << "                       mu |eta|    < " << muon_etaMax << endl;
   cout << endl;
   cout << "Systematic errors" << endl;
   cout << " Luminosity:  " << systErrLumi * 100 << "%" << endl;
   cout << " Efficiency:  " << systErrEff * 100 << "%" << endl;
   cout << " ttbar:      " << systErrMC[TTBAR-1] * 100 << "%" << endl;
   cout << " Z->tautau:   " << systErrMC[ZTT-1] * 100 << "%" << endl;
   cout << " WW:          " << systErrMC[WW-1] * 100 << "%" << endl;
   cout << " WZ:          " << systErrMC[WZ-1] * 100 << "%" << endl;
   cout << " tW, tbarW:  " << systErrMC[TW-1] * 100 << "%" << endl;
   cout << " Z->mumu:     " << systErrMC[ZMM-1] * 100 << "%" << endl;
   cout << " Z->ee:       " << systErrMC[ZEE-1] * 100 << "%" << endl;
   cout << " ZZ:          " << systErrMC[ZZ-1] * 100 << "%" << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("M_emu         |       >%5.0fGeV/c^2          |        > 120GeV/c^2          |        > 200GeV/c^2          |\n", minInvMass);
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
          emuMass.at(DATA).at(ALL)->Integral(), sqrt(emuMass.at(DATA).at(ALL)->Integral()),         
          emuMass.at(DATA).at(ALL)->Integral(13,151), sqrt(emuMass.at(DATA).at(ALL)->Integral(13,151)),
          emuMass.at(DATA).at(ALL)->Integral(21,151), sqrt(emuMass.at(DATA).at(ALL)->Integral(21,151)));
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("nb ttbar      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(TTBAR).at(ALL)->Integral(), emuMass.at(TTBAR).at(ALL)->Integral() * systErrMCLuEff[TTBAR-1], 
          emuMass.at(TTBAR).at(ALL)->Integral(13,151), emuMass.at(TTBAR).at(ALL)->Integral(13,151) * systErrMCLuEff[TTBAR-1], 
          emuMass.at(TTBAR).at(ALL)->Integral(21,151), emuMass.at(TTBAR).at(ALL)->Integral(21,151) * systErrMCLuEff[TTBAR-1]);
   printf("nb Ztautau    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(ZTT).at(ALL)->Integral(), emuMass.at(ZTT).at(ALL)->Integral() * systErrMCLuEff[ZTT-1], 
          emuMass.at(ZTT).at(ALL)->Integral(13,151), emuMass.at(ZTT).at(ALL)->Integral(13,151) * systErrMCLuEff[ZTT-1], 
          emuMass.at(ZTT).at(ALL)->Integral(21,151), emuMass.at(ZTT).at(ALL)->Integral(21,151) * systErrMCLuEff[ZTT-1]);
   printf("nb WW         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(WW).at(ALL)->Integral(), emuMass.at(WW).at(ALL)->Integral() * systErrMCLuEff[WW-1], 
          emuMass.at(WW).at(ALL)->Integral(13,151), emuMass.at(WW).at(ALL)->Integral(13,151) * systErrMCLuEff[WW-1], 
          emuMass.at(WW).at(ALL)->Integral(21,151), emuMass.at(WW).at(ALL)->Integral(21,151) * systErrMCLuEff[WW-1]);
   printf("nb WZ         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(WZ).at(ALL)->Integral(), emuMass.at(WZ).at(ALL)->Integral() * systErrMCLuEff[WZ-1], 
          emuMass.at(WZ).at(ALL)->Integral(13,151), emuMass.at(WZ).at(ALL)->Integral(13,151) * systErrMCLuEff[WZ-1], 
          emuMass.at(WZ).at(ALL)->Integral(21,151), emuMass.at(WZ).at(ALL)->Integral(21,151) * systErrMCLuEff[WZ-1]);
   printf("nb tW         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(TW).at(ALL)->Integral(), emuMass.at(TW).at(ALL)->Integral() * systErrMCLuEff[TW-1], 
          emuMass.at(TW).at(ALL)->Integral(13,151), emuMass.at(TW).at(ALL)->Integral(13,151) * systErrMCLuEff[TW-1], 
          emuMass.at(TW).at(ALL)->Integral(21,151), emuMass.at(TW).at(ALL)->Integral(21,151) * systErrMCLuEff[TW-1]);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   printf("nb WJets      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(WJET).at(ALL)->Integral(), emuMass.at(WJET).at(ALL)->Integral() * systErrMCLuEff[WJET-1], 
          emuMass.at(WJET).at(ALL)->Integral(13,151), emuMass.at(WJET).at(ALL)->Integral(13,151) * systErrMCLuEff[WJET-1], 
          emuMass.at(WJET).at(ALL)->Integral(21,151), emuMass.at(WJET).at(ALL)->Integral(21,151) * systErrMCLuEff[WJET-1]);
   printf("nb Zmumu      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(ZMM).at(ALL)->Integral(), emuMass.at(ZMM).at(ALL)->Integral() * systErrMCLuEff[ZMM-1], 
          emuMass.at(ZMM).at(ALL)->Integral(13,151), emuMass.at(ZMM).at(ALL)->Integral(13,151) * systErrMCLuEff[ZMM-1], 
          emuMass.at(ZMM).at(ALL)->Integral(21,151), emuMass.at(ZMM).at(ALL)->Integral(21,151) * systErrMCLuEff[ZMM-1]);
   printf("nb Zee        | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(ZEE).at(ALL)->Integral(), emuMass.at(ZEE).at(ALL)->Integral() * systErrMCLuEff[ZEE-1], 
          emuMass.at(ZEE).at(ALL)->Integral(13,151), emuMass.at(ZEE).at(ALL)->Integral(13,151) * systErrMCLuEff[ZEE-1], 
          emuMass.at(ZEE).at(ALL)->Integral(21,151), emuMass.at(ZEE).at(ALL)->Integral(21,151) * systErrMCLuEff[ZEE-1]);
   //printf("nb ZZ         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
   //       emuMass.at(ZZ).at(ALL)->Integral(), emuMass.at(ZZ).at(ALL)->Integral() * systErrMCLuEff[ZZ-1], 
   //       emuMass.at(ZZ).at(ALL)->Integral(13,151), emuMass.at(ZZ).at(ALL)->Integral(13,151) * systErrMCLuEff[ZZ-1], 
   //       emuMass.at(ZZ).at(ALL)->Integral(21,151), emuMass.at(ZZ).at(ALL)->Integral(21,151) * systErrMCLuEff[ZZ-1]);
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;

   //printf("TOT ttlike    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
   //       emuMass.at(TTBAR).at(ALL)->Integral() + emuMass.at(ZTT).at(ALL)->Integral() + emuMass.at(WW).at(ALL)->Integral() + emuMass.at(WZ).at(ALL)->Integral() + emuMass.at(TW).at(ALL)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, ttLikeSamples, ALL, 7, 151),  
   //       emuMass.at(TTBAR).at(ALL)->Integral(13,151) + emuMass.at(ZTT).at(ALL)->Integral(13,151) + emuMass.at(WW).at(ALL)->Integral(13,151) + emuMass.at(WZ).at(ALL)->Integral(13,151) + emuMass.at(TW).at(ALL)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, ttLikeSamples, ALL, 13, 151), 
   //       emuMass.at(TTBAR).at(ALL)->Integral(21,151) + emuMass.at(ZTT).at(ALL)->Integral(21,151) + emuMass.at(WW).at(ALL)->Integral(21,151) + emuMass.at(WZ).at(ALL)->Integral(21,151) + emuMass.at(TW).at(ALL)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, ttLikeSamples, ALL, 21, 151));
   //printf("TOT contam    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
   //       emuMass.at(WJET).at(ALL)->Integral() + emuMass.at(ZMM).at(ALL)->Integral() + emuMass.at(ZEE).at(ALL)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, contamSamples, ALL, 7, 151), 
   //       emuMass.at(WJET).at(ALL)->Integral(13,151) + emuMass.at(ZMM).at(ALL)->Integral(13,151) + emuMass.at(ZEE).at(ALL)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, contamSamples, ALL, 13, 151), 
   //       emuMass.at(WJET).at(ALL)->Integral(21,151) + emuMass.at(ZMM).at(ALL)->Integral(21,151) + emuMass.at(ZEE).at(ALL)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, contamSamples, ALL, 21, 151));
   //cout << "-------------------------------------------------------------------------------------------------" << endl;

   //SUM MC
   for (int p = 1; p < nbFile; ++p) {
      for (unsigned int k = 0; k < 3; ++k) {
         for (int q = p + 1; q < nbFile; ++q) {
            emuMass.at(p).at(k)->Add(emuMass.at(q).at(k));
            emuMassEB.at(p).at(k)->Add(emuMassEB.at(q).at(k));
            emuMassEE.at(p).at(k)->Add(emuMassEE.at(q).at(k));
            emu_mass_accVSgood.at(p).at(k)->Add(emu_mass_accVSgood.at(q).at(k));

            met.at(p).at(k)->Add(met.at(q).at(k));
            nVtx.at(p).at(k)->Add(nVtx.at(q).at(k));
            dPhi.at(p).at(k)->Add(dPhi.at(q).at(k));
            elePt.at(p).at(k)->Add(elePt.at(q).at(k));
            eleEta.at(p).at(k)->Add(eleEta.at(q).at(k));
            elePhi.at(p).at(k)->Add(elePhi.at(q).at(k));
            eleId1.at(p).at(k)->Add(eleId1.at(q).at(k));
            eleId2.at(p).at(k)->Add(eleId2.at(q).at(k));
            eleId3.at(p).at(k)->Add(eleId3.at(q).at(k));
            eleIso1.at(p).at(k)->Add(eleIso1.at(q).at(k));
            eleIso2.at(p).at(k)->Add(eleIso2.at(q).at(k));
            eleIso3.at(p).at(k)->Add(eleIso3.at(q).at(k));
 
            muIsoCombRel.at(p).at(k)->Add(muIsoCombRel.at(q).at(k));
            muPtEleOPtMu.at(p).at(k)->Add(muPtEleOPtMu.at(q).at(k));
            muPtPlusOPtMinus.at(p).at(k)->Add(muPtPlusOPtMinus.at(q).at(k));
            muPt.at(p).at(k)->Add(muPt.at(q).at(k));
            muEta.at(p).at(k)->Add(muEta.at(q).at(k));
            muPhi.at(p).at(k)->Add(muPhi.at(q).at(k));
            muId1.at(p).at(k)->Add(muId1.at(q).at(k));
            muId2.at(p).at(k)->Add(muId2.at(q).at(k));
            muId3.at(p).at(k)->Add(muId3.at(q).at(k));
            muIso1.at(p).at(k)->Add(muIso1.at(q).at(k));
            muIso2.at(p).at(k)->Add(muIso2.at(q).at(k));
            muIso3.at(p).at(k)->Add(muIso3.at(q).at(k));

            numOfJets.at(p).at(k)->Add(numOfJets.at(q).at(k));
            numOfJetsPt15.at(p).at(k)->Add(numOfJetsPt15.at(q).at(k));
         }
      }
   }
   printf("TOT ttlike    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(TTBAR).at(ALL)->Integral() - emuMass.at(WJET).at(ALL)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, ttLikeSamples, ALL, 1, 151, true),  
          emuMass.at(TTBAR).at(ALL)->Integral(13,151) - emuMass.at(WJET).at(ALL)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, ttLikeSamples, ALL, 13, 151, true), 
          emuMass.at(TTBAR).at(ALL)->Integral(21,151) - emuMass.at(WJET).at(ALL)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, ttLikeSamples, ALL, 21, 151, true));
   printf("TOT contam    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(WJET).at(ALL)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, contamSamples, ALL, 1, 151, true), 
          emuMass.at(WJET).at(ALL)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, contamSamples, ALL, 13, 151, true), 
          emuMass.at(WJET).at(ALL)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, contamSamples, ALL, 21, 151, true));
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;

   printf("TOT MC        | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMass.at(1).at(ALL)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, allSamples, ALL, 1, 151, true),
          emuMass.at(1).at(ALL)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, allSamples, ALL, 13, 151, true),
          emuMass.at(1).at(ALL)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, allSamples, ALL, 21, 151, true));
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   cout << endl << endl << endl;


   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("nb LS DATA    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n", 
          emuMass.at(DATA).at(LS)->Integral(), sqrt(emuMass.at(DATA).at(LS)->Integral()), 
          emuMass.at(DATA).at(LS)->Integral(13,151), sqrt(emuMass.at(DATA).at(LS)->Integral(13,151)), 
          emuMass.at(DATA).at(LS)->Integral(21,151), sqrt(emuMass.at(DATA).at(LS)->Integral(21,151))); 
   printf("nb LS MC      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(1).at(LS)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, allSamples, LS, 1, 151, true), 
          emuMass.at(1).at(LS)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, allSamples, LS, 13, 151, true), 
          emuMass.at(1).at(LS)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, allSamples, LS, 21, 151, true)); 
   cout << endl;
   printf("nb OS DATA    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
          emuMass.at(DATA).at(OS)->Integral(), sqrt(emuMass.at(DATA).at(OS)->Integral()), 
          emuMass.at(DATA).at(OS)->Integral(13,151), sqrt(emuMass.at(DATA).at(OS)->Integral(13,151)), 
          emuMass.at(DATA).at(OS)->Integral(21,151), sqrt(emuMass.at(DATA).at(OS)->Integral(21,151))); 
   printf("nb OS MC      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
          emuMass.at(1).at(OS)->Integral(), CalcSystErr(emuMass, systErrMCLuEff, allSamples, OS, 1, 151, true), 
          emuMass.at(1).at(OS)->Integral(13,151), CalcSystErr(emuMass, systErrMCLuEff, allSamples, OS, 13, 151, true), 
          emuMass.at(1).at(OS)->Integral(21,151), CalcSystErr(emuMass, systErrMCLuEff, allSamples, OS, 21, 151, true)); 
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   cout << "nb e+mu+     = " << nb_plus_plus << endl;
   cout << "nb e+mu-     = " << nb_plus_minus << endl;
   cout << "nb e-mu+     = " << nb_minus_plus << endl;
   cout << "nb e-mu-     = " << nb_minus_minus << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;

   vector<TH1F *> helper;
   for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMass_" + sign[k] + suffix[QCD], "emuMass_" + sign[k] + suffix[QCD], 150, 0., 1500.));
   emuMass.push_back(helper);

   if (useQCDShape) {
      TH1F *qcdBg = new TH1F("qcdBg", "qcdBg", 150, 0., 1500.);
      qcdBg->Add(emuMass.at(DATA).at(LS));
      qcdBg->Add(emuMass.at(1).at(LS), -1.);
      for (int i = 0; i < qcdBg->GetNbinsX() + 2; ++i) {
         if (qcdBg->GetBinContent(i) < 0) qcdBg->SetBinContent(i, 0.);
      }

      // add QCD histo to histogram stack
      for (int p = 1; p < nbFile + 1; ++p) {
         for (unsigned int k = 0; k < 3; ++k) {
            float factor = 1.;
            if (k == 0) factor = 2.;
            emuMass.at(p).at(k)->Add(qcdBg, factor);
         }
      }
   }
   else {
      // FIXME calculate qcd scale factor
      float QCDScaleFactor = 1 + 2 * (emuMass.at(DATA).at(LS)->Integral() - emuMass.at(1).at(LS)->Integral()) / emuMass.at(1).at(ALL)->Integral();
      cout << endl << endl << "calculated QCD scale factor from data: " << QCDScaleFactor << endl;
      // add QCD histo to histogram stack
      for (int p = 1; p < nbFile + 1; ++p) {
         for (unsigned int k = 0; k < 3; ++k) {
            TH1F *qcdBg = new TH1F("qcdBg" + sign[k] + suffix[p], "qcdBg" + sign[k] + suffix[p], 150, 0., 1500.);
            qcdBg->Add(emuMass.at(DATA).at(k));
            float factor = 1.;
            // apply QCD scale factor
            if (calcQCDScaleFactor) factor = QCDScaleFactor;
            else {
               cout << "applied preset QCD scale factor:       " << QCD_ScaleFactor << endl;
               factor = QCD_ScaleFactor;
            }
            emuMass.at(p).at(k)->Add(qcdBg, factor - 1.);
         }
      }
   }

   systErrMC.push_back(2 * sqrt(emuMass.at(DATA).at(LS)->Integral() + pow(CalcSystErr(emuMass, systErrMCLuEff, allSamples, LS, 1, 151, true), 2)) / emuMass.at(QCD).at(ALL)->Integral());
   systErrMCLuEff.push_back(systErrMC[QCD]);

   vector<bool> onlyQCD(QCD-1, false);
   onlyQCD.push_back(true);
   contamSamples.push_back(true);
   allSamples.push_back(true);

   cout << endl;
   cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << "QCD events from LS spectrum:" << endl;
   printf("nb QCD        | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) |\n",
          emuMass.at(QCD).at(ALL)->Integral(), CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 1, 151, true), 100 * CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 7, 151, true) / emuMass.at(QCD).at(ALL)->Integral(), 
          emuMass.at(QCD).at(ALL)->Integral(13,151), CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 13, 151, true), 100 * CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 13, 151, true) / emuMass.at(QCD).at(ALL)->Integral(13,151),
          emuMass.at(QCD).at(ALL)->Integral(21,151), CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 21, 151, true), 100 * CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 21, 151, true) / emuMass.at(QCD).at(ALL)->Integral(21,151));
   printf("%% of total MC |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |\n", 
          100 * emuMass.at(QCD).at(ALL)->Integral() / (emuMass.at(1).at(ALL)->Integral() - emuMass.at(QCD).at(ALL)->Integral()), 100 * CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 1, 151, true) / (emuMass.at(1).at(ALL)->Integral() - emuMass.at(QCD).at(ALL)->Integral()), 
          100 * emuMass.at(QCD).at(ALL)->Integral(13,151) / (emuMass.at(1).at(ALL)->Integral(13,151) - emuMass.at(QCD).at(ALL)->Integral(13,151)), 100 * CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 13, 151, true) / (emuMass.at(1).at(ALL)->Integral(13,151) - emuMass.at(QCD).at(ALL)->Integral(13,151)),
          100 * emuMass.at(QCD).at(ALL)->Integral(21,151) / (emuMass.at(1).at(ALL)->Integral(21,151) - emuMass.at(QCD).at(ALL)->Integral(21,151)), 100 * CalcSystErrWithQCD(emuMass, systErrMCLuEff, onlyQCD, ALL, 21, 151, true) / (emuMass.at(1).at(ALL)->Integral(21,151) - emuMass.at(QCD).at(ALL)->Integral(21,151)));
   cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;

   cout << endl;
   cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   if (useQCDShape) cout << "AFTER ADDING QCD CONTRIBUTION:" << endl;
   else cout << "AFTER CORRECTION WITH QCD FACTOR:" << endl;
   cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   printf("M_emu         |      >%6.0fGeV/c^2          |        > 120GeV/c^2          |         > 200GeV/c^2         |         > 500GeV/c^2         |\n", minInvMass);
   cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n", 
           emuMass.at(DATA).at(ALL)->Integral(), sqrt(emuMass.at(DATA).at(ALL)->Integral()), 
           emuMass.at(DATA).at(ALL)->Integral(13, 151), sqrt((emuMass.at(DATA).at(ALL))->Integral(13, 151)), 
           emuMass.at(DATA).at(ALL)->Integral(21, 151), sqrt(emuMass.at(DATA).at(ALL)->Integral(21, 151)), 
           emuMass.at(DATA).at(ALL)->Integral(51, 151), sqrt(emuMass.at(DATA).at(ALL)->Integral(51, 151)));
   printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
           emuMass.at(1).at(ALL)->Integral(), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 1, 151, true), 
           emuMass.at(1).at(ALL)->Integral(13, 151), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 13, 151, true), 
           emuMass.at(1).at(ALL)->Integral(21, 151), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 21, 151, true), 
           emuMass.at(1).at(ALL)->Integral(51, 151), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 51, 151, true));
   cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl << endl;

   cout << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("M_emu         |    %6.0f - 120GeV/c^2       |      120 - 200GeV/c^2        |       200 - 400GeV/c^2       |\n", minInvMass);
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
           emuMass.at(DATA).at(ALL)->Integral(1, 12), sqrt(emuMass.at(DATA).at(ALL)->Integral(1, 12)), 
           emuMass.at(DATA).at(ALL)->Integral(13, 20), sqrt((emuMass.at(DATA).at(ALL))->Integral(13, 20)), 
           emuMass.at(DATA).at(ALL)->Integral(21, 40), sqrt(emuMass.at(DATA).at(ALL)->Integral(21, 40)));
   printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
           emuMass.at(1).at(ALL)->Integral(1, 12), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 1, 12, true), 
           emuMass.at(1).at(ALL)->Integral(13, 20), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 13, 20, true), 
           emuMass.at(1).at(ALL)->Integral(21, 40), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, ALL, 21, 40, true)); 
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   printf("nb data OS    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
           emuMass.at(DATA).at(OS)->Integral(1, 12), sqrt(emuMass.at(DATA).at(OS)->Integral(1, 12)), 
           emuMass.at(DATA).at(OS)->Integral(13, 20), sqrt((emuMass.at(DATA).at(OS))->Integral(13, 20)), 
           emuMass.at(DATA).at(OS)->Integral(21, 40), sqrt(emuMass.at(DATA).at(OS)->Integral(21, 40)));
   printf("nb MC OS      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
           emuMass.at(1).at(OS)->Integral(1, 12), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, OS, 1, 12, true), 
           emuMass.at(1).at(OS)->Integral(13, 20), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, OS, 13, 20, true), 
           emuMass.at(1).at(OS)->Integral(21, 40), CalcSystErrWithQCD(emuMass, systErrMCLuEff, allSamples, OS, 21, 40, true)); 
   cout << "-----------------------------------------------------------------------------------------------------------" << endl << endl;

   emu_dilepton->Add(emuMass.at(1).at(ALL));
   emu_ewk->Add(emuMass.at(6).at(ALL));
   //emu_jet->Add(emuMass.at(8).at(ALL));

   //WRITING
   stringstream ssOutfile;
   ssOutfile << outfileName << LumiFactor << "pb-1.root";
   TFile *output = new TFile(ssOutfile.str().c_str(), "recreate");

   output->cd();

   emu_dilepton->Write();
   emu_ewk->Write();
   emu_jet->Write();

   //emuLoose_data_nValidPv->Write();
   //emuLoose_ttbar_nValidPv->Write();
   //emuLoose_ztautau_nValidPv->Write();

   //emuLoose_dataOverTtbar_nValidPv->Write();
   //emuLoose_dataOverZtautau_nValidPv->Write();

   for (int p = 0; p < nbFile; ++p) {
      for (unsigned int k = 0; k < 3; ++k) {
         emuMass.at(p).at(k)->Write();
         emuMassEB.at(p).at(k)->Write();
         emuMassEE.at(p).at(k)->Write();
         emu_mass_accVSgood.at(p).at(k)->Write();

         met.at(p).at(k)->Write();
         nVtx.at(p).at(k)->Write();
         dPhi.at(p).at(k)->Write();
         elePt.at(p).at(k)->Write();
         eleEta.at(p).at(k)->Write();
         elePhi.at(p).at(k)->Write();
         eleId1.at(p).at(k)->Write();
         eleId2.at(p).at(k)->Write();
         eleId3.at(p).at(k)->Write();
         eleIso1.at(p).at(k)->Write();
         eleIso2.at(p).at(k)->Write();
         eleIso3.at(p).at(k)->Write();

         muIsoCombRel.at(p).at(k)->Write();
         muPtEleOPtMu.at(p).at(k)->Write();
         muPtPlusOPtMinus.at(p).at(k)->Write();
         muPt.at(p).at(k)->Write();
         muEta.at(p).at(k)->Write();
         muPhi.at(p).at(k)->Write();
         muId1.at(p).at(k)->Write();
         muId2.at(p).at(k)->Write();
         muId3.at(p).at(k)->Write();
         muIso1.at(p).at(k)->Write();
         muIso2.at(p).at(k)->Write();
         muIso3.at(p).at(k)->Write();

         numOfJets.at(p).at(k)->Write();
         numOfJetsPt15.at(p).at(k)->Write();
      }
   }
   for (unsigned int k = 0; k < 3; ++k) emuMass.at(QCD).at(k)->Write();

   output->Close();
}

float CalcSystErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked)
{
   float err2 = 0.;

   // for stacked histograms subtract all contributions below the current one 
   if (stacked) {
      for (unsigned int i = 0; i < samples.size() && i < histos.size()-1; ++i) {
         if (samples[i]) {
            float numEv;
            if (i < histos.size()-2) numEv = histos.at(i+1).at(region)->Integral(lowerBin, upperBin) - histos.at(i+2).at(region)->Integral(lowerBin, upperBin);
            else numEv = histos.at(i+1).at(region)->Integral(lowerBin, upperBin);
            err2 += pow(numEv, 2) * errors[i] * errors[i];
         }
      }
   } else {
      for (unsigned int i = 0; i < samples.size(); ++i) {
         if (samples[i])
            err2 += histos.at(i+1).at(region)->Integral(lowerBin, upperBin) * histos.at(i+1).at(region)->Integral(lowerBin, upperBin) * errors[i] * errors[i];
      }
   }
 
   return sqrt(err2);
}

float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked)
{
   float qcdErr = 0.;
   vector<bool> allButQCD(QCD-1, true);

   qcdErr = 2 * sqrt(histos.at(DATA).at(LS)->Integral(lowerBin, upperBin) + pow(CalcSystErr(histos, errors, allButQCD, LS, lowerBin, upperBin, stacked), 2));
   errors.back() = qcdErr / histos.at(QCD).at(region)->Integral(lowerBin, upperBin);

   return CalcSystErr(histos, errors, samples, region, lowerBin, upperBin, stacked);
}

int Trigger(TFile *inFile, unsigned int &entry, int &prescale, unsigned int *trig, const int &selector)
{
   int c_HLT_Mu15_Photon20_CaloIdL;
   int c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
   int c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;
   int c_prescale_HLT_Mu15_Photon20_CaloIdL;
   int c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
   int c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;

   TBranch *b_HLT_Mu15_Photon20_CaloIdL;
   TBranch *b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
   TBranch *b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;
   TBranch *b_prescale_HLT_Mu15_Photon20_CaloIdL;
   TBranch *b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
   TBranch *b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;

   // Get the tree and connect the necessary variables
   inFile->cd();
   TTree *thetree;
   thetree = (TTree*)(inFile)->Get("gsfcheckerjob/tree");

   thetree->SetBranchAddress("HLT_Mu15_Photon20_CaloIdL", &c_HLT_Mu15_Photon20_CaloIdL, &b_HLT_Mu15_Photon20_CaloIdL);
   thetree->SetBranchAddress("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL", &c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL, &b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL);
   thetree->SetBranchAddress("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL", &c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL, &b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL);
   thetree->SetBranchAddress("prescale_HLT_Mu15_Photon20_CaloIdL", &c_prescale_HLT_Mu15_Photon20_CaloIdL, &b_prescale_HLT_Mu15_Photon20_CaloIdL);
   thetree->SetBranchAddress("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL", &c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL, &b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL);
   thetree->SetBranchAddress("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL", &c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL, &b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL);

   thetree->GetEntry(entry);

   if (selector == 1) {
      prescale = c_prescale_HLT_Mu15_Photon20_CaloIdL;
      if (c_HLT_Mu15_Photon20_CaloIdL >= 0) trig[0]++;
      return c_HLT_Mu15_Photon20_CaloIdL;
   } else if (selector == 2) {
      prescale = c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
      if (c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL >= 0) trig[1]++;
      return c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
   } else if (selector == 3) {
      prescale = c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;
      if (c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL >= 0) trig[2]++;
      return c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;
   } else {
      // select an unprescaled trigger
      if (c_prescale_HLT_Mu15_Photon20_CaloIdL == 1) {
         prescale = c_prescale_HLT_Mu15_Photon20_CaloIdL;
         trig[0]++;
         if (c_runnumber < runs_HLT_Mu15_Photon20_CaloIdL.first) runs_HLT_Mu15_Photon20_CaloIdL.first = c_runnumber;
         if (c_runnumber > runs_HLT_Mu15_Photon20_CaloIdL.second) runs_HLT_Mu15_Photon20_CaloIdL.second = c_runnumber;
         return c_HLT_Mu15_Photon20_CaloIdL;
      } else if (c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL == 1) {
         prescale = c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
         trig[1]++;
         if (c_runnumber < runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.first) runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.first = c_runnumber;
         if (c_runnumber > runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.second) runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL.second = c_runnumber;
         return c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
      } else if (c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL == 1) {
         prescale = c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;
         trig[2]++;
         if (c_runnumber < runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.first) runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.first = c_runnumber;
         if (c_runnumber > runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.second) runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL.second = c_runnumber;
         return c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;
      }
   }

   cout << "Prescale alert! No unprescaled trigger found." << endl;
   prescale = 0;
   return -1;
}
