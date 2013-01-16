#include "TFile.h"
#include "TDCacheFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TPaveLabel.h"
#include "TParameter.h"
#include "THashList.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TStopwatch.h"

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
#define TTBAR700TO1000 2
#define TTBAR1000UP 3
#define ZTT 4
#define WW 5
#define WZ 6
#define TW 7
#define WJET 8
#define ZMM 9
#define ZEE 10
#define ZZ 11

#define ALL 0
#define LS 1
#define OS 2

using namespace std;

int Trigger(TFile *inFile, unsigned int &entry, int &prescale, unsigned int *trig, const int &selector = 0);
float CorrectEnergy(TFile *inFile, int &pos);

pair<unsigned int, unsigned int> runs_HLT_Mu22_Photon22_CaloIdL(99999999, 0);
pair<unsigned int, unsigned int> runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL(99999999, 0);
pair<unsigned int, unsigned int> runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL(99999999, 0);
//RUN ID
unsigned int c_runnumber;
unsigned int c_eventnumber;
unsigned int c_luminosityBlock;

void emuSpectrum()
{
   TStopwatch timer;
   timer.Start();
   // parameters /////////////////////////////////////////////////////////////
   float LumiFactor = 19619.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   TParameter<float> lumi("lumi", LumiFactor);

   // DATA file
   //TString dataFile = "dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/data2012/MuEG_Run2012A+B+C_13Jul2012+06Aug2012+PromptReco-v1+v2_gct1_37_12201pb-1.root";
   TString dataFile = "file:////user/treis/data2012/MuEG_Run2012A+B+C+D_13Jul2012+06Aug2012+24Aug2012+11Dec2012+PromptReco-Cv2+Dv1_Cert_190456-208686_gct1_46_19619pb-1.root";
   // pile up histogram
   TString puFile = "file:////user/treis/data2012/pileup/pileup_runA+B+C-ReReco+D-Prompt_puJSON-190389-208686_MuEG.root";

   string outfileName = "testEmuSpec";
   //string outfileName = "test";

   //unsigned selection = 0;   // HEEP v3.2
   //unsigned selection = 1;   // HEEP v4.0strong
   //unsigned selection = 2;   // HEEP v4.0
   unsigned selection = 3;   // HEEP v4.1

   // scale factors
   // muon factors  https://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=214870
   TParameter<float> trgEff("trgEff", 0.99 * 0.891); // ele L1 eff times Mu40 eff measured by Z' to mumu
   TParameter<float> trgEffLowEta("trgEffLowEta", 0.99 * 0.941); // ele L1 eff times Mu40 eff measured by Z' to mumu for |eta|<0.9
   TParameter<float> trgEffMidEta("trgEffMidEta", 0.99 * 0.844); // ele L1 eff times Mu40 eff measured by Z' to mumu for 0.9<|eta|<1.2
   TParameter<float> trgEffHighEta("trgEffHighEta", 0.99 * 0.827); // ele L1 eff times Mu40 eff measured by Z' to mumu for 1.2<|eta|
   TParameter<float> trgDataMcScaleFactor("trgDataMcScaleFactor", 0.981); // scale factor between data and mc measured by Z' to mumu for Mu40
   TParameter<float> trgDataMcScaleFactorLowEta("trgDataMcScaleFactorLowEta", 0.981); // scale factor between data and mc measured by Z' to mumu for Mu40 for |eta|<0.9
   TParameter<float> trgDataMcScaleFactorMidEta("trgDataMcScaleFactorMidEta", 0.962); // scale factor between data and mc measured by Z' to mumu for Mu40 for 0.9<|eta|<1.2
   TParameter<float> trgDataMcScaleFactorHighEta("trgDataMcScaleFactorHighEta", 0.990); // scale factor between data and mc measured by Z' to mumu for Mu40 for 1.2<|eta|
   // epsilon_cand from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgCommissioningAndPhysicsDeliverables#Electron_reconstruction_effi_AN1
   TParameter<float> eleScaleFactorEB("eleScaleFactorEB", 0.995 * 1.007); // data/MC scale for epsilon_cand (>50GeV) * epsilon_id (at Z peak)
   TParameter<float> eleScaleFactorEE("eleScaleFactorEE", 0.994 * 1.004); // data/MC scale for epsilon_cand (>50GeV) * epsilon_id (at Z peak)
   TParameter<float> muScaleFactor("muScaleFactor", 0.982);
   TParameter<float> muScaleFactorLowEta("muScaleFactorLowEta", 0.979); // for |eta|<0.9
   TParameter<float> muScaleFactorMidEta("muScaleFactorMidEta", 0.967); // for 0.9<|eta|<1.2
   TParameter<float> muScaleFactorHighEta("muScaleFactorHighEta", 0.993); // for 1.2<|eta|
   TParameter<float> lumiScaleFactor("lumiScaleFactor", 0.980);  // not used // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1
   TParameter<float> lumiScaleFactorEB("lumiScaleFactorEB", 0.993);  // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1
   TParameter<float> lumiScaleFactorEE("lumiScaleFactorEE", 0.948);  // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1

   // global systematic errors
   TParameter<float> systErrLumi("systErrLumi", 0.022);
   TParameter<float> systErrEff("systErrEff", 0.007); // muon err (0.002) & ele err (0.007)

   bool usePUInfo = true;
   bool lowMassPuOnly = false;
   float puMassCut = 120.;
   // selection cuts /////////////////////////////////////////////////////////
   float minInvMass = 0.;

   //MUON selection
   float muon_pt = 35.;
   float muon_dptOverPt = 0.3;
   float muon_etaMax = 2.4;
   int muon_nHitsMinGlobal = 0;
   int muon_nHitsMinPixel = 1;
   int muon_nHitsMinMuon = 1;
   int muon_nLayersMin = 6; 
   float muon_impactParamMaxXY = 0.2;   // in cm
   //float muon_impactParamMaxZ = 0.5;   // in cm , not used
   int muon_nSegMatchMin = 2;
   float muon_relIsoCutMax = 0.1;

   //BARREL
   float bar_et = 0.;
   float bar_hoE = 0.;
   float bar_DEta = 0.;
   float bar_DPhi = 0.;
   float bar_e2x5e5x5 = 0.;
   float bar_e1x5e5x5 = 0.;
   float bar_isoEcalHcal1_1 = 0.;
   float bar_isoEcalHcal1_2 = 0.;
   float bar_isoEcalHcalRho = 0.;
   float bar_isoTrack = 0.;
   float bar_dxy = 0.;
   int bar_missInnerHits = 0;

   //ENDCAP
   float end_et = 0.;
   float end_hoE = 0.;
   float end_DEta = 0.;
   float end_DPhi = 0.;
   float end_sigmaietaieta = 0.;
   float end_isoEcalHcal1_1 = 0.;
   float end_isoEcalHcal1_1_1 = 0.;
   float end_isoEcalHcal1_1_2 = 0.;
   float end_isoEcalHcal1_2 = 0.;
   float end_isoEcalHcalRho = 0.;
   float end_isoTrack = 0.;
   float end_dxy = 0.;
   int end_missInnerHits = 0;

   switch (selection) {
      case 0: {  // HEEP v3.2
         outfileName += "HEEP32_";
         //BARREL
         bar_et = 35.;
         bar_hoE = 0.05;
         bar_DEta = 0.005;
         bar_DPhi = 0.06;
         bar_e2x5e5x5 = 0.94;
         bar_e1x5e5x5 = 0.83;
         bar_isoEcalHcal1_1 = 2.;
         bar_isoEcalHcal1_2 = 0.03;
         bar_isoTrack = 5.;
         bar_missInnerHits = 0;
     
         //ENDCAP
         end_et = 35.;
         end_hoE = 0.05;
         end_DEta = 0.007 ;
         end_DPhi = 0.06;
         end_sigmaietaieta = 0.03;
         end_isoEcalHcal1_1 = 2.5;
         end_isoEcalHcal1_2 = 0.03;
         end_isoTrack = 5.;
         end_missInnerHits = 0;
         break;
      }
      case 1: {   // HEEP v4.0 strong
         outfileName += "HEEP4strong";
         //BARREL
         bar_et = 35.;
         bar_hoE = 0.05;
         bar_DEta = 0.005;
         bar_DPhi = 0.06;
         bar_e2x5e5x5 = 0.94;
         bar_e1x5e5x5 = 0.83;
         bar_isoEcalHcal1_1 = 2.;
         bar_isoEcalHcal1_2 = 0.03;
         bar_isoEcalHcalRho = 0.28;
         bar_isoTrack = 5.;
         bar_missInnerHits = 0;
     
         //ENDCAP
         end_et = 35.;
         end_hoE = 0.05;
         end_DEta = 0.007 ;
         end_DPhi = 0.06;
         end_sigmaietaieta = 0.03;
         end_isoEcalHcal1_1_1 = 2.5;
         end_isoEcalHcal1_1_2 = 1.;
         end_isoEcalHcal1_2 = 0.03;
         end_isoEcalHcalRho = 0.28;
         end_isoTrack = 5.;
         end_missInnerHits = 0;
         break;
 
      }
      case 2: {  // HEEP v4.0
         outfileName += "HEEP4_";
         //BARREL
         bar_et = 35.;
         bar_hoE = 0.05;
         bar_DEta = 0.005;
         bar_DPhi = 0.06;
         bar_e2x5e5x5 = 0.94;
         bar_e1x5e5x5 = 0.83;
         bar_isoEcalHcal1_1 = 2.;
         bar_isoEcalHcal1_2 = 0.03;
         bar_isoEcalHcalRho = 0.28;
         bar_isoTrack = 5.;
         bar_missInnerHits = 0;
     
         //ENDCAP
         end_et = 35.;
         end_hoE = 0.05;
         end_DEta = 0.007 ;
         end_DPhi = 0.06;
         end_sigmaietaieta = 0.03;
         end_isoEcalHcal1_1_1 = 2.5;
         end_isoEcalHcal1_1_2 = 1.;
         end_isoEcalHcal1_2 = 0.03;
         end_isoEcalHcalRho = 0.28;
         end_isoTrack = 5.;
         end_missInnerHits = 0;
         break;
      }
      case 3: {  // HEEP v4.1
         outfileName += "HEEP41_";
         //BARREL
         bar_et = 35.;
         bar_hoE = 0.05;
         bar_DEta = 0.005;
         bar_DPhi = 0.06;
         bar_e2x5e5x5 = 0.94;
         bar_e1x5e5x5 = 0.83;
         bar_isoEcalHcal1_1 = 2.;
         bar_isoEcalHcal1_2 = 0.03;
         bar_isoEcalHcalRho = 0.28;
         bar_isoTrack = 5.;
         bar_dxy = 0.02;
         bar_missInnerHits = 1;
     
         //ENDCAP
         end_et = 35.;
         end_hoE = 0.05;
         end_DEta = 0.007 ;
         end_DPhi = 0.06;
         end_sigmaietaieta = 0.03;
         end_isoEcalHcal1_1_1 = 2.5;
         end_isoEcalHcal1_1_2 = 1.;
         end_isoEcalHcal1_2 = 0.03;
         end_isoEcalHcalRho = 0.28;
         end_isoTrack = 5.;
         end_dxy = 0.05;
         end_missInnerHits = 1;
         break;
      }
    }
   ///////////////////////////////////////////////////////////////////////////

   TH1::SetDefaultSumw2(kTRUE);

   ///////////////////////////////////////////////////////////////////////////
   // INPUT FILES
   vector<pair<TFile *, float> > input;
   vector<bool> storeGenMTtbar;
   THashList systErrMCs;
   THashList mcWeights;

   // DATA
   TFile *inData = TFile::Open(dataFile);
   input.push_back(make_pair(inData, 1.)); //DATA        0 black
   storeGenMTtbar.push_back(0);

   // MC
   TFile *inTTbar = TFile::Open("file:////user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
   //input.push_back(make_pair(inTTbar, 225.197 / 28150723.)); // NLO
   input.push_back(make_pair(inTTbar, 234. / 28150723.));  // approx NNLO
   //TFile *inTTbar = TFile::Open("file:////user/treis/mcsamples/TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1_AODSIM_gct1_41_6736135ev.root");
   //input.push_back(make_pair(inTTbar, 225.197 / 6736135.)); 
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar", 0.067));
   storeGenMTtbar.push_back(1);

   TFile *inTTbar700to1000 = TFile::Open("file:////user/treis/mcsamples/TT_Mtt-700to1000_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3082812ev.root");
   input.push_back(make_pair(inTTbar700to1000, 15.614 / 3082812. * 234./211.));  // ttbar  mtt 700to1000
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar700to1000", 0.15));
   storeGenMTtbar.push_back(1);

   TFile *inTTbar1000up = TFile::Open("file:////user/treis/mcsamples/TT_Mtt-1000toInf_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_1249111ev.root");
   input.push_back(make_pair(inTTbar1000up, 2.954 / 1249111. * 234./211.));  // ttbar  mtt>1000
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar1000up", 0.15));
   storeGenMTtbar.push_back(1);

   TFile *inZtt = TFile::Open("file:////user/treis/mcsamples/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3295238ev.root");
   input.push_back(make_pair(inZtt, 1915.1 / 3295238.)); //Ztautau
   systErrMCs.Add(new TParameter<float>("systErrMcDyTauTau", 0.054));
   storeGenMTtbar.push_back(0);

   TFile *inWW = TFile::Open("file:////user/treis/mcsamples/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_10000431ev.root");
   input.push_back(make_pair(inWW, 54.838 / 10000431.)); //WW
   systErrMCs.Add(new TParameter<float>("systErrMcWW", 0.035));
   storeGenMTtbar.push_back(0);

   TFile *inWZ = TFile::Open("file:////user/treis/mcsamples/WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_10000283ev.root");
   input.push_back(make_pair(inWZ, 33.21 / 10000283.)); //WZ
   systErrMCs.Add(new TParameter<float>("systErrMcWZ", 0.038));
   storeGenMTtbar.push_back(0);

   TFile *inTW = TFile::Open("file:////user/treis/mcsamples/T+Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_991118ev.root");
   input.push_back(make_pair(inTW, 22.2 / 991118.)); //tW
   systErrMCs.Add(new TParameter<float>("systErrMcTW", 0.069));
   storeGenMTtbar.push_back(0);

   TFile *inWJet = TFile::Open("file:////user/treis/mcsamples/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_76102995ev.root");
   input.push_back(make_pair(inWJet, 36257.2 / 76102995.)); //W+jet
   systErrMCs.Add(new TParameter<float>("systErrMcWJets", 0.05));
   storeGenMTtbar.push_back(0);

   TFile *inZmm = TFile::Open("file:////user/treis/mcsamples/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3293740ev.root");
   input.push_back(make_pair(inZmm, 1915. / 3293740.)); //Zmumu
   systErrMCs.Add(new TParameter<float>("systErrMcDyMuMu", 0.054));
   storeGenMTtbar.push_back(0);

   TFile *inZee = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3297045ev.root");
   input.push_back(make_pair(inZee, 1915. / 3297045.)); //Zee
   systErrMCs.Add(new TParameter<float>("systErrMcDyEE", 0.054));
   storeGenMTtbar.push_back(0);

   TFile *inZZ = TFile::Open("file:////user/treis/mcsamples/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_9799908ev.root");
   input.push_back(make_pair(inZZ, 17.654 / 9799908.)); //ZZ
   systErrMCs.Add(new TParameter<float>("systErrMcZZ", 0.025));
   storeGenMTtbar.push_back(0);

   TFile *inTTbar22l = TFile::Open("file:////user/treis/mcsamples/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_16365457ev.root");
   input.push_back(make_pair(inTTbar22l, 13.43 / 16365457. * 234./(13.43+53.2+53.4))); //TT to 2l
   systErrMCs.Add(new TParameter<float>("systErrMcTtJets2l", 0.067));
   storeGenMTtbar.push_back(1);

   TFile *inTTbar21l = TFile::Open("file:////user/treis/mcsamples/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_gct1_46_25424818ev.root");
   input.push_back(make_pair(inTTbar21l, 53.2 / 25424818. * 234./(13.43+53.2+53.4))); //TT to 1l1jet
   systErrMCs.Add(new TParameter<float>("systErrMcTtJets1l1jet", 0.067));
   storeGenMTtbar.push_back(1);

   TFile *inTTW = TFile::Open("file:////user/treis/mcsamples/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_196046ev.root");
   input.push_back(make_pair(inTTW, 0.2149 / 196046.)); //TTW
   systErrMCs.Add(new TParameter<float>("systErrMcTtW", 0.));
   storeGenMTtbar.push_back(1);

   TFile *inTTWW = TFile::Open("file:////user/treis/mcsamples/TTWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_217820ev.root");
   input.push_back(make_pair(inTTWW, 0.002037 / 217820.)); //TTWW
   systErrMCs.Add(new TParameter<float>("systErrMcTtWW", 0.));
   storeGenMTtbar.push_back(1);

   //TFile *inWJet2 = TFile::Open("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/mcsamples2012/WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_35_22447541ev.root");
   //input.push_back(make_pair(inWJet2, 1.91068E-5)); //W+jet
   //systErrMCs.Add(new TParameter<float>("systErrMcWJets70t0100", 0.05));
   //storeGenMTtbar.push_back(0);
   //TFile *inWJet3 = TFile::Open("file:////user/treis/mcsamples/WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1_AODSIM_gct1_33_12766180ev.root");
   //input.push_back(make_pair(inWJet3, 1.79302E-5)); //W+jet
   //systErrMCs.Add(new TParameter<float>("systErrMcWJets100up", 0.05));
   //storeGenMTtbar.push_back(0);
   int nbFile = input.size();
   ///////////////////////////////////////////////////////////////////////////

   // strings for histogram names
   vector<TString> suffix;
   suffix.push_back("data");
   suffix.push_back("ttbar");
   suffix.push_back("ttbar700to1000");
   suffix.push_back("ttbar1000up");
   suffix.push_back("ztautau");
   suffix.push_back("ww");
   suffix.push_back("wz");
   suffix.push_back("tw");
   suffix.push_back("wjets");
   suffix.push_back("zmumu");
   suffix.push_back("zee");
   suffix.push_back("zz");
   suffix.push_back("ttbarto2l");
   suffix.push_back("ttbarto1l1jet");
   suffix.push_back("ttbarw");
   suffix.push_back("ttbarww");

   // counting variables
   int nb_plus_plus = 0;
   int nb_plus_minus = 0;
   int nb_minus_plus = 0;
   int nb_minus_minus = 0;

   //vector<TH1F *> emu_plus_plus;
   //vector<TH1F *> emu_plus_minus;
   //vector<TH1F *> emu_minus_plus;
   //vector<TH1F *> emu_minus_minus;

   // data pileup histogram
   TFile *dataPuInfile = TFile::Open(puFile);
   dataPuInfile->Cd("");
   TH1F *puData = (TH1F *)gDirectory->Get("pileup");
   puData->SetDirectory(0);
   puData->SetName("puData");
   dataPuInfile->Close();
   TH1F *puDataNorm = (TH1F *)puData->DrawNormalized()->Clone("puDataNorm");

   //GLOBAL
   float c_pfmet;
   int c_jetColl_size;
   float c_jetPt[100];

   //PRIM VTX
   int c_pvsize;

   float c_rho;
   int c_trueNVtx;

   //GSF
   int c_gsf_size;
   float c_gsf_gsfet[100];
   float c_gsf_eta[100];
   float c_gsf_theta[100];
   float c_gsf_phi[100];
   float c_gsf_dxy_firstPVtx[100];
   float c_gsf_dz_beamSpot[100];
   float c_gsf_dz_firstPVtx[100];
   bool c_gsf_isecaldriven[100];
   int c_gsf_charge[100];
   float c_gsf_deltaeta[100];
   float c_gsf_deltaphi[100];
   float c_gsf_e1x5overe5x5[100];
   float c_gsf_e2x5overe5x5[100];
   float c_gsf_sigmaetaeta[100];
   float c_gsf_sigmaIetaIeta[100];
   float c_gsf_hovere[100];
   int c_gsf_nLostInnerHits[100];
   float c_gsf_ecaliso[100];
   float c_gsf_hcaliso1[100];
   float c_gsf_hcaliso2[100];
   float c_gsf_trackiso[100];

   float c_gsfsc_e[100];
   float c_gsfsc_eta[100];
   float c_gsfsc_phi[100];

   //MUONS
   int c_muon_size;
   float c_muon_pt[100];
   float c_muon_eta[100];
   float c_muon_phi[100];
   float c_muon_ptError[100];
   int c_muon_charge[100];
   int c_muon_nhitstrack[100];
   int c_muon_nhitspixel[100];
   int c_muon_nhitsmuons[100];
   int c_muon_nlayerswithhits[100];
   int c_muon_nSegmentMatch[100];
   bool c_muon_isTrackerMuon[100];
   float c_muon_normChi2[100];
   float c_muon_dz_beamSpot[100];
   float c_muon_dz_firstPVtx[100];
   float c_muon_dxy_cmsCenter[100];
   float c_muon_dxy_beamSpot[100];
   float c_muon_dxy_firstPVtx[100];
   float c_muon_trackIso03[100];
   float c_muon_emIso03[100];
   float c_muon_hadIso03[100];

   float c_genPair_mass;

   //RUN ID
   TBranch *b_runnumber;
   TBranch *b_eventnumber;
   TBranch *b_luminosityBlock;

   //GLOBAL
   TBranch *b_pfmet;
   TBranch *b_jetColl_size;
   TBranch *b_jetPt;

   //PRIM VTX
   TBranch *b_pvsize;

   TBranch *b_rho;
   TBranch *b_trueNVtx;

   //GSF
   TBranch *b_gsf_size;
   TBranch *b_gsf_gsfet;
   TBranch *b_gsf_eta;
   TBranch *b_gsf_theta;
   TBranch *b_gsf_phi;
   TBranch *b_gsf_dxy_firstPVtx;
   TBranch *b_gsf_dz_beamSpot;
   TBranch *b_gsf_dz_firstPVtx;
   TBranch *b_gsf_isecaldriven;
   TBranch *b_gsf_charge;
   TBranch *b_gsf_deltaeta;
   TBranch *b_gsf_deltaphi;
   TBranch *b_gsf_e1x5overe5x5;
   TBranch *b_gsf_e2x5overe5x5;
   TBranch *b_gsf_sigmaetaeta;
   TBranch *b_gsf_sigmaIetaIeta;
   TBranch *b_gsf_hovere;
   TBranch *b_gsf_nLostInnerHits;
   TBranch *b_gsf_ecaliso;
   TBranch *b_gsf_hcaliso1;
   TBranch *b_gsf_hcaliso2;
   TBranch *b_gsf_trackiso;

   TBranch *b_gsfsc_e;
   TBranch *b_gsfsc_eta;
   TBranch *b_gsfsc_phi;

   //MUONS
   TBranch *b_muon_size;
   TBranch *b_muon_pt;
   TBranch *b_muon_eta;
   TBranch *b_muon_phi;
   TBranch *b_muon_ptError;
   TBranch *b_muon_charge;
   TBranch *b_muon_nhitstrack;
   TBranch *b_muon_nhitspixel;
   TBranch *b_muon_nhitsmuons;
   TBranch *b_muon_nlayerswithhits;
   TBranch *b_muon_nSegmentMatch;
   TBranch *b_muon_isTrackerMuon;
   TBranch *b_muon_normChi2;
   TBranch *b_muon_dz_beamSpot;
   TBranch *b_muon_dz_firstPVtx;
   TBranch *b_muon_dxy_cmsCenter;
   TBranch *b_muon_dxy_beamSpot;
   TBranch *b_muon_dxy_firstPVtx;
   TBranch *b_muon_trackIso03;
   TBranch *b_muon_emIso03;
   TBranch *b_muon_hadIso03;

   TBranch *b_genPair_mass;

   unsigned int dataTrig[3] = {0, 0, 0};
   unsigned int dataEntries = 0;

   // output file
   stringstream ssOutfile;
   ssOutfile << outfileName << LumiFactor << "pb-1.root";
   TFile *output = new TFile(ssOutfile.str().c_str(), "recreate");

   // write parameters
   lumi.Write();
   trgEff.Write();
   trgEffLowEta.Write();
   trgEffMidEta.Write();
   trgEffHighEta.Write();
   trgDataMcScaleFactor.Write();
   trgDataMcScaleFactorLowEta.Write();
   trgDataMcScaleFactorMidEta.Write();
   trgDataMcScaleFactorHighEta.Write();
   lumiScaleFactor.Write();
   lumiScaleFactorEB.Write();
   lumiScaleFactorEE.Write();
   eleScaleFactorEB.Write();
   eleScaleFactorEE.Write();
   muScaleFactor.Write();
   muScaleFactorLowEta.Write();
   muScaleFactorMidEta.Write();
   muScaleFactorHighEta.Write();
   systErrLumi.Write();
   systErrEff.Write();
   systErrMCs.Write("systErrMCs", TObject::kSingleKey);

   runs_HLT_Mu22_Photon22_CaloIdL.first = 99999999;
   runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first = 99999999;
   runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first = 99999999;
   dataTrig[0] = 0;
   dataTrig[1] = 0;
   dataTrig[2] = 0;
   dataEntries = 0;
   // GETTING FILES
   for (int p = 0; p < nbFile; ++p) {
      cout << "accessing file " << p + 1 << ": " << input[p].first->GetName() << endl;
      input[p].first->Cd("");

      // Get the TREE
      TTree *thetree;
      thetree = (TTree*)(input[p].first)->Get("gsfcheckerjob/tree");

      //RUN ID
      thetree->SetBranchAddress("runnumber", &c_runnumber, &b_runnumber);
      thetree->SetBranchAddress("eventnumber", &c_eventnumber, &b_eventnumber);
      thetree->SetBranchAddress("luminosityBlock", &c_luminosityBlock, &b_luminosityBlock);

      //GLOBAL
      thetree->SetBranchAddress("pfmet", &c_pfmet, &b_pfmet);
      thetree->SetBranchAddress("JetColl_size", &c_jetColl_size, &b_jetColl_size);
      thetree->SetBranchAddress("Jet_pt", &c_jetPt, &b_jetPt);

      //PRIM VTX
      thetree->SetBranchAddress("pvsize", &c_pvsize, &b_pvsize);

      thetree->SetBranchAddress("rho", &c_rho, &b_rho);
      thetree->SetBranchAddress("trueNVtx", &c_trueNVtx, &b_trueNVtx);

      //GSF
      thetree->SetBranchAddress("gsf_size", &c_gsf_size, &b_gsf_size);
      thetree->SetBranchAddress("gsf_gsfet", &c_gsf_gsfet, &b_gsf_gsfet);
      thetree->SetBranchAddress("gsf_eta", &c_gsf_eta, &b_gsf_eta);
      thetree->SetBranchAddress("gsf_theta", &c_gsf_theta, &b_gsf_theta);
      thetree->SetBranchAddress("gsf_phi", &c_gsf_phi, &b_gsf_phi);
      thetree->SetBranchAddress("gsf_dxy_firstPVtx", &c_gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
      thetree->SetBranchAddress("gsf_dz_beamSpot", &c_gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
      thetree->SetBranchAddress("gsf_dz_firstPVtx", &c_gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
      thetree->SetBranchAddress("gsf_isecaldriven", &c_gsf_isecaldriven, &b_gsf_isecaldriven);
      thetree->SetBranchAddress("gsf_charge", &c_gsf_charge, &b_gsf_charge);
      thetree->SetBranchAddress("gsf_deltaeta", &c_gsf_deltaeta, &b_gsf_deltaeta);
      thetree->SetBranchAddress("gsf_deltaphi", &c_gsf_deltaphi, &b_gsf_deltaphi);
      thetree->SetBranchAddress("gsf_e1x5overe5x5", &c_gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
      thetree->SetBranchAddress("gsf_e2x5overe5x5", &c_gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
      thetree->SetBranchAddress("gsf_sigmaetaeta", &c_gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
      thetree->SetBranchAddress("gsf_sigmaIetaIeta", &c_gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
      thetree->SetBranchAddress("gsf_hovere", &c_gsf_hovere, &b_gsf_hovere);
      thetree->SetBranchAddress("gsf_nLostInnerHits", &c_gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
      thetree->SetBranchAddress("gsf_ecaliso", &c_gsf_ecaliso, &b_gsf_ecaliso);
      thetree->SetBranchAddress("gsf_hcaliso1", &c_gsf_hcaliso1, &b_gsf_hcaliso1);
      thetree->SetBranchAddress("gsf_hcaliso2", &c_gsf_hcaliso2, &b_gsf_hcaliso2);
      thetree->SetBranchAddress("gsf_trackiso", &c_gsf_trackiso, &b_gsf_trackiso);
      thetree->SetBranchAddress("gsfsc_e", &c_gsfsc_e, &b_gsfsc_e);
      thetree->SetBranchAddress("gsfsc_eta", &c_gsfsc_eta, &b_gsfsc_eta);
      thetree->SetBranchAddress("gsfsc_phi", &c_gsfsc_phi, &b_gsfsc_phi);

      //MUONS
      thetree->SetBranchAddress("muon_size", &c_muon_size, &b_muon_size);
      thetree->SetBranchAddress("muon_pt", &c_muon_pt, &b_muon_pt);
      thetree->SetBranchAddress("muon_eta", &c_muon_eta, &b_muon_eta);
      thetree->SetBranchAddress("muon_phi", &c_muon_phi, &b_muon_phi);
      thetree->SetBranchAddress("muon_ptError", &c_muon_ptError, &b_muon_ptError);
      thetree->SetBranchAddress("muon_charge", &c_muon_charge, &b_muon_charge);
      thetree->SetBranchAddress("muon_nhitstrack", &c_muon_nhitstrack, &b_muon_nhitstrack);
      thetree->SetBranchAddress("muon_nhitspixel", &c_muon_nhitspixel, &b_muon_nhitspixel);
      thetree->SetBranchAddress("muon_nhitsmuons", &c_muon_nhitsmuons, &b_muon_nhitsmuons);
      thetree->SetBranchAddress("muon_nlayerswithhits", &c_muon_nlayerswithhits, &b_muon_nlayerswithhits);
      thetree->SetBranchAddress("muon_nSegmentMatch", &c_muon_nSegmentMatch, &b_muon_nSegmentMatch);
      thetree->SetBranchAddress("muon_isTrackerMuon", &c_muon_isTrackerMuon, &b_muon_isTrackerMuon);
      thetree->SetBranchAddress("muon_normChi2", &c_muon_normChi2, &b_muon_normChi2);
      thetree->SetBranchAddress("muon_dz_beamSpot", &c_muon_dz_beamSpot, &b_muon_dz_beamSpot);
      thetree->SetBranchAddress("muon_dz_firstPVtx", &c_muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
      thetree->SetBranchAddress("muon_dxy_cmsCenter", &c_muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
      thetree->SetBranchAddress("muon_dxy_beamSpot", &c_muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
      thetree->SetBranchAddress("muon_dxy_firstPVtx", &c_muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
      thetree->SetBranchAddress("muon_trackIso03", &c_muon_trackIso03, &b_muon_trackIso03);
      thetree->SetBranchAddress("muon_emIso03", &c_muon_emIso03, &b_muon_emIso03);
      thetree->SetBranchAddress("muon_hadIso03", &c_muon_hadIso03, &b_muon_hadIso03);
      
      if (storeGenMTtbar[p]) thetree->SetBranchAddress("genPair_mass", &c_genPair_mass, &b_genPair_mass);

      // enable only used branches
      thetree->SetBranchStatus("*", 0);
      thetree->SetBranchStatus("runnumber", 1);
      thetree->SetBranchStatus("eventnumber", 1);
      thetree->SetBranchStatus("luminosityBlock", 1);
      thetree->SetBranchStatus("HLT_Mu22_Photon22_CaloIdL", 1);
      thetree->SetBranchStatus("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
      thetree->SetBranchStatus("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
      thetree->SetBranchStatus("prescale_HLT_Mu22_Photon22_CaloIdL", 1);
      thetree->SetBranchStatus("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
      thetree->SetBranchStatus("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1);
      thetree->SetBranchStatus("pfmet", 1);
      thetree->SetBranchStatus("JetColl_size", 1);
      thetree->SetBranchStatus("Jet_pt", 1);
      thetree->SetBranchStatus("pvsize", 1);
      thetree->SetBranchStatus("rho", 1);
      thetree->SetBranchStatus("trueNVtx", 1);
      thetree->SetBranchStatus("gsf_size", 1);
      thetree->SetBranchStatus("gsf_gsfet", 1);
      thetree->SetBranchStatus("gsf_eta", 1);
      thetree->SetBranchStatus("gsf_theta", 1);
      thetree->SetBranchStatus("gsf_phi", 1);
      thetree->SetBranchStatus("gsf_dxy_firstPVtx", 1);
      thetree->SetBranchStatus("gsf_dz_beamSpot", 1);
      thetree->SetBranchStatus("gsf_dz_firstPVtx", 1);
      thetree->SetBranchStatus("gsf_isecaldriven", 1);
      thetree->SetBranchStatus("gsf_charge", 1);
      thetree->SetBranchStatus("gsf_deltaeta", 1);
      thetree->SetBranchStatus("gsf_deltaphi", 1);
      thetree->SetBranchStatus("gsf_e1x5overe5x5", 1);
      thetree->SetBranchStatus("gsf_e2x5overe5x5", 1);
      thetree->SetBranchStatus("gsf_sigmaetaeta", 1);
      thetree->SetBranchStatus("gsf_sigmaIetaIeta", 1);
      thetree->SetBranchStatus("gsf_hovere", 1);
      thetree->SetBranchStatus("gsf_nLostInnerHits", 1);
      thetree->SetBranchStatus("gsf_ecaliso", 1);
      thetree->SetBranchStatus("gsf_hcaliso1", 1);
      thetree->SetBranchStatus("gsf_hcaliso2", 1);
      thetree->SetBranchStatus("gsf_trackiso", 1);
      thetree->SetBranchStatus("gsfsc_e", 1);
      thetree->SetBranchStatus("gsfsc_eta", 1);
      thetree->SetBranchStatus("gsfsc_phi", 1);
      thetree->SetBranchStatus("muon_size", 1);  
      thetree->SetBranchStatus("muon_pt", 1);
      thetree->SetBranchStatus("muon_eta", 1);
      thetree->SetBranchStatus("muon_phi", 1);
      thetree->SetBranchStatus("muon_ptError", 1); 
      thetree->SetBranchStatus("muon_charge", 1);
      thetree->SetBranchStatus("muon_nhitstrack", 1);
      thetree->SetBranchStatus("muon_nhitspixel", 1);
      thetree->SetBranchStatus("muon_nhitsmuons", 1);
      thetree->SetBranchStatus("muon_nlayerswithhits", 1);
      thetree->SetBranchStatus("muon_nSegmentMatch", 1);
      thetree->SetBranchStatus("muon_isTrackerMuon", 1);
      thetree->SetBranchStatus("muon_normChi2", 1);
      thetree->SetBranchStatus("muon_dz_beamSpot", 1);
      thetree->SetBranchStatus("muon_dz_firstPVtx", 1);
      thetree->SetBranchStatus("muon_dxy_beamSpot", 1);
      thetree->SetBranchStatus("muon_dxy_firstPVtx", 1);
      thetree->SetBranchStatus("muon_trackIso03", 1);
      thetree->SetBranchStatus("muon_emIso03", 1);
      thetree->SetBranchStatus("muon_hadIso03", 1);
      if (storeGenMTtbar[p]) thetree->SetBranchStatus("genPair_mass", 1);      // get the histogram with the true number of vertices

      TH1F *puMc = new TH1F("puMc" + suffix[p], "puMc" + suffix[p], 100, 0., 100.);
      thetree->Draw("trueNVtx>>puMc" + suffix[p]);
      TH1F *puMcNorm = new TH1F("dummy" + suffix[p], "dummy" + suffix[p], 100, 0., 100.);
      if (p > DATA) puMcNorm = (TH1F *)puMc->DrawNormalized();
      puMcNorm->SetName("puMc" + suffix[p] + "Norm");
      puMcNorm->SetTitle("puMc" + suffix[p] + "Norm");
      // calculate the pu weights
      TH1F *puWeights = new TH1F("puWeight" + suffix[p], "puWeight" + suffix[p], 100, 0., 100.);
      if (p > DATA) puWeights->Divide(puDataNorm, puMcNorm);

      // tree with event data
      bool passTrg = false;
      float emuInvMass = 0.;
      int evtRegion;
      int eCharge;
      int muCharge;
      float totWeight = 1.;
      float puWeight = 1.;
      TTree *emuTree = new TTree("emuTree_" + suffix[p], "emuTree_" + suffix[p]);
      emuTree->Branch("runnr", &c_runnumber, "runnr/i");
      emuTree->Branch("eventnr", &c_eventnumber, "eventnr/i");
      emuTree->Branch("lumiSec", &c_luminosityBlock, "lumiSec/i");
      emuTree->Branch("passTrg", &passTrg, "passTrg/O");
      emuTree->Branch("mass", &emuInvMass, "mass/F");
      emuTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
      emuTree->Branch("eCharge", &eCharge, "eCharge/I");
      emuTree->Branch("muCharge", &muCharge, "muCharge/I");
      emuTree->Branch("weight", &totWeight, "weight/F");
      emuTree->Branch("puWeight", &puWeight, "puWeight/F");
      if (storeGenMTtbar[p]) emuTree->Branch("genMTtbar", &c_genPair_mass, "genMTtbar/F");
      // control variables
      float nVtx = 0.;
      float dZFstPVtx = 0.;
      float dXYFstPVtx = 0.;
      float dPhi = 0.;
      float dEta = 0.;
      float eleEt = 0.;
      float eleEta = 0.;
      float elePhi = 0.;
      float eleDEta = 0.;
      float eleDPhi = 0.;
      float eleHOE = 0.;
      float eleSigmaIEIE = 0.;
      float eleEcalIso = 0.;
      float eleHcalIso12 = 0.;
      float eleTrkIso = 0.;
      float eleLostHits = 0.;
      float muIsoCombRel = 0.;
      float muEtEleOPtMu = 0.;
      float muPtPlusOPtMinus = 0.;
      float muPt = 0.;
      float muEta = 0.;
      float muPhi = 0.;
      float muHitLayers = 0.;
      float muPxlHits = 0.;
      float muMuHits = 0.;
      float muDZFstPVtx = 0.;
      float muDXYFstPVtx = 0.;
      float muNSeg = 0.;
      float muTrkIso03 = 0.;
      float numOfJets = 0.;
      float numOfJetsPt20 = 0.;
      float numOfJetsPt30 = 0.;
      emuTree->Branch("pfMet", &c_pfmet, "pfMet/F");
      emuTree->Branch("nVtx", &nVtx, "nVtx/F");
      emuTree->Branch("dZFstPVtx", &dZFstPVtx, "dZFstPVtx/F");
      emuTree->Branch("dXYFstPVtx", &dXYFstPVtx, "dXYFstPVtx/F");
      emuTree->Branch("rho", &c_rho, "rho/F");
      emuTree->Branch("dPhi", &dPhi, "dPhi/F");
      emuTree->Branch("dEta", &dEta, "dEta/F");
      emuTree->Branch("eleEt", &eleEt, "eleEt/F");
      emuTree->Branch("eleEta", &eleEta, "eleEta/F");
      emuTree->Branch("elePhi", &elePhi, "elePhi/F");
      emuTree->Branch("eleDEta", &eleDEta, "eleDEta/F");
      emuTree->Branch("eleDPhi", &eleDPhi, "eleDPhi/F");
      emuTree->Branch("eleHOE", &eleHOE, "eleHOE/F");
      emuTree->Branch("eleSigmaIEIE", &eleSigmaIEIE, "eleSigmaIEIE/F");
      emuTree->Branch("eleEcalIso", &eleEcalIso, "eleEcalIso/F");
      emuTree->Branch("eleHcalIso12", &eleHcalIso12, "eleHcalIso12/F");
      emuTree->Branch("eleTrkIso", &eleTrkIso, "eleTrkIso/F");
      emuTree->Branch("eleLostHits", &eleLostHits, "eleLostHits/F");
      emuTree->Branch("muIsoCombRel", &muIsoCombRel, "muIsoCombRel/F");
      emuTree->Branch("muEtEleOPtMu", &muEtEleOPtMu, "muEtEleOPtMu/F");
      emuTree->Branch("muPtPlusOPtMinus", &muPtPlusOPtMinus, "muPtPlusOPtMinus/F");
      emuTree->Branch("muPt", &muPt, "muPt/F");
      emuTree->Branch("muEta", &muEta, "muEta/F");
      emuTree->Branch("muPhi", &muPhi, "muPhi/F");
      emuTree->Branch("muHitLayers", &muHitLayers, "muHitLayers/F");
      emuTree->Branch("muPxlHits", &muPxlHits, "muPxlHits/F");
      emuTree->Branch("muMuHits", &muMuHits, "muMuHits/F");
      emuTree->Branch("muDZFstPVtx", &muDZFstPVtx, "muDZFstPVtx/F");
      emuTree->Branch("muDXYFstPVtx", &muDXYFstPVtx, "muDXYFstPVtx/F");
      emuTree->Branch("muNSeg", &muNSeg, "muNSeg/F");
      emuTree->Branch("muTrkIso03", &muTrkIso03, "muTrkIso03/F");
      emuTree->Branch("numOfJets", &numOfJets, "numOfJets/F");
      emuTree->Branch("numOfJetsPt20", &numOfJetsPt20, "numOfJetsPt20/F");
      emuTree->Branch("numOfJetsPt30", &numOfJetsPt30, "numOfJetsPt30/F");

      Long64_t nentries = (*thetree).GetEntries();
      cout << nentries << " events" << endl;
      if (p == DATA) {
         cout << "-----------------------------------------------------------------------------------------------------------" << endl;
         cout << "M_emu > 600GeV/c^2       |     run:lumi:event    |   M_emu  |"
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
      unsigned int evCounter = 0;
      //LOOP OVER EVENTS
      //for (unsigned int i = 0; i < 10000; ++i) {
      for (unsigned int i = 0; i < nentries; ++i) {
         if (i % 50000 == 0) cout << "Processing event " << i << endl;
         thetree->GetEntry(i);

         // at least one gsf electron and one muon above the threshold
         if (c_gsf_size < 1 || c_muon_size < 1) continue;

         int prescale = 0;
         passTrg = true;
         if (p == DATA && Trigger(input[p].first, i, prescale, dataTrig) < 1) passTrg = false;
         if (p != DATA) {
            if (i < nentries * dataTrig[0] / dataEntries) {
               if (Trigger(input[p].first, i, prescale, trig, 1) < 1) passTrg = false;
            } else if (i < nentries * (dataTrig[0] + dataTrig[1]) / dataEntries) {
               if (Trigger(input[p].first, i, prescale, trig, 2) < 1) passTrg = false;
            } else {
               if (Trigger(input[p].first, i, prescale, trig, 3) < 1) passTrg = false;
            }
         }

         vector<int> GSF_passHEEP;
         vector<int> GSF_passACC;

         vector<int> MU_passGOOD;
         vector<int> MU_passACC;

         // set the PU weight
         if (usePUInfo && p > DATA) puWeight = puWeights->GetBinContent(puWeights->FindBin(c_trueNVtx));
         float weight = puWeight;

         //LOOP OVER ELES
         for (int j = 0; j < c_gsf_size; ++j) {
            // correct the energy
            //if (p == 1) c_gsf_gsfet[j] = CorrectEnergy(input[p].first, j);
            //CLEANING : FAKE ELES FROM MUONS
            bool fakeEle = false;
            for (int k = 0; k < c_muon_size; ++k) {
               //if (c_muon_pt[k] < muon_pt) continue;
               float DeltaR = sqrt((c_gsf_eta[j] - c_muon_eta[k]) * (c_gsf_eta[j] - c_muon_eta[k]) + (c_gsf_phi[j] - c_muon_phi[k]) * (c_gsf_phi[j] - c_muon_phi[k]));
               if (DeltaR < 0.1) {
                  fakeEle = true;
                  break;
               }
            }
            if (fakeEle) continue;

            if (selection == 0) {
               // HEEP v3.2 /////////////////////////////////////////////////////////////////////////////
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
                   && c_gsf_sigmaIetaIeta[j] < end_sigmaietaieta
                   && ((c_gsf_gsfet[j] < 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < end_isoEcalHcal1_1)
                       ||
                       (c_gsf_gsfet[j] >= 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (end_isoEcalHcal1_1 + end_isoEcalHcal1_2 * (c_gsf_gsfet[j] - 50.))))
                   && c_gsf_trackiso[j] < end_isoTrack
                   && c_gsf_nLostInnerHits[j] <= end_missInnerHits
                  ) GSF_passHEEP.push_back(j);
            } else {
               //BARREL HEEP
               if (fabs(c_gsfsc_eta[j]) < 1.442
                   && c_gsf_gsfet[j] > bar_et
                   && c_gsf_isecaldriven[j]
                   && c_gsf_hovere[j] < bar_hoE
                   && fabs(c_gsf_deltaeta[j]) < bar_DEta
                   && fabs(c_gsf_deltaphi[j]) < bar_DPhi
                   && ((c_gsf_e2x5overe5x5[j] > bar_e2x5e5x5) || (c_gsf_e1x5overe5x5[j] > bar_e1x5e5x5))
                   && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * c_gsf_gsfet[j] + bar_isoEcalHcalRho * c_rho)
                   && c_gsf_trackiso[j] < bar_isoTrack
                   && c_gsf_nLostInnerHits[j] <= bar_missInnerHits
                  ) if ((selection == 3 && fabs(c_gsf_dxy_firstPVtx[j]) <= bar_dxy) || selection != 3) GSF_passHEEP.push_back(j);
   
               //ENDCAP HEEP
               if ((fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5)
                   && c_gsf_gsfet[j] > end_et
                   && c_gsf_isecaldriven[j]
                   && c_gsf_hovere[j] < end_hoE
                   && fabs(c_gsf_deltaeta[j]) < end_DEta
                   && fabs(c_gsf_deltaphi[j]) < end_DPhi
                   && c_gsf_sigmaIetaIeta[j] < end_sigmaietaieta
                   && ((c_gsf_gsfet[j] < 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (end_isoEcalHcal1_1_1 + end_isoEcalHcalRho * c_rho))
                      || (c_gsf_gsfet[j] >= 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (end_isoEcalHcal1_1_2 + end_isoEcalHcal1_2 * c_gsf_gsfet[j] + end_isoEcalHcalRho * c_rho)))
                   && c_gsf_trackiso[j] < end_isoTrack
                   && c_gsf_nLostInnerHits[j] <= end_missInnerHits
                  ) if ((selection == 3 && fabs(c_gsf_dxy_firstPVtx[j]) <= end_dxy) || selection != 3) GSF_passHEEP.push_back(j);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////

            //GSF IN ACCEPTANCE
            if (c_gsf_gsfet[j] > bar_et
                && (fabs(c_gsfsc_eta[j]) < 1.442 || (fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5))
               ) GSF_passACC.push_back(j);
         }

         //LOOP OVER MUS
         for (int j = 0; j < c_muon_size; ++j) {
            //MU PASS GOOD
            if (c_muon_pt[j] > muon_pt
                && c_muon_ptError[j] / c_muon_pt[j] < muon_dptOverPt
                && fabs(c_muon_eta[j]) < muon_etaMax
                && c_muon_nhitstrack[j] >= muon_nHitsMinGlobal
                && c_muon_nhitspixel[j] >= muon_nHitsMinPixel
                && c_muon_nhitsmuons[j] >= muon_nHitsMinMuon
                && c_muon_nlayerswithhits[j] >= muon_nLayersMin
                && fabs(c_muon_dxy_firstPVtx[j]) < muon_impactParamMaxXY
                //&& fabs(c_muon_dz_firstPVtx[j]) < muon_impactParamMaxZ
                && c_muon_nSegmentMatch[j] >= muon_nSegMatchMin
                && c_muon_isTrackerMuon[j]
                && c_muon_trackIso03[j] / c_muon_pt[j] < muon_relIsoCutMax
               ) MU_passGOOD.push_back(j);

            //MU PASS ACC
            if (c_muon_pt[j] > muon_pt
                && fabs(c_muon_eta[j]) < muon_etaMax
               ) MU_passACC.push_back(j);
         }

         // veto when there are more than one good candidates
         if (GSF_passHEEP.size() > 1) continue;
         if (MU_passGOOD.size() > 1) continue;

         //search highest pt passing muon if there are more than one
         unsigned int MU_leadingPassGOOD = 0;
         if (MU_passGOOD.size() > 1) {
            for (unsigned int j = 0; j < MU_passGOOD.size(); ++j)
               if (c_muon_pt[MU_passGOOD[j]] > c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]) MU_leadingPassGOOD = j;
         }

         //HEEP - MU_GOOD
         if (GSF_passHEEP.size() > 0 && MU_passGOOD.size() > 0) {
            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(c_gsf_gsfet[GSF_passHEEP[0]], c_gsf_eta[GSF_passHEEP[0]], c_gsf_phi[GSF_passHEEP[0]], 0.000511);
            mu1.SetPtEtaPhiM(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], 0.10566);

            double invMass = (ele1 + mu1).M();

            // vertex matching of tracks by dz
            //if (fabs(c_gsf_dz_beamSpot[GSF_passHEEP[0]] - c_muon_dz_beamSpot[MU_passGOOD[MU_leadingPassGOOD]]) > 0.05) continue;

            //MASS CUT
            if (invMass < minInvMass) continue;

            ++evCounter;

            float CombRelIso = (c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]]) / c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];

            // set correction factors according to detector region
            float Elec_ScaleFactor = 1.;
            float Lumi_ScaleFactor = 1.;
            if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 1.442) {
              Elec_ScaleFactor = eleScaleFactorEB.GetVal();
              Lumi_ScaleFactor = lumiScaleFactorEB.GetVal();
            }
            else if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 2.5) {
              Elec_ScaleFactor = eleScaleFactorEE.GetVal();
              Lumi_ScaleFactor = lumiScaleFactorEE.GetVal();
            }
            float MCemuScaleFactor = Elec_ScaleFactor * muScaleFactor.GetVal() * Lumi_ScaleFactor * trgDataMcScaleFactor.GetVal(); 
            float MCScaleFactor = Elec_ScaleFactor * muScaleFactor.GetVal() * Lumi_ScaleFactor * trgDataMcScaleFactor.GetVal(); 

            // correction for trigger
            //if (p == ZMM) MCemuScaleFactor = MCScaleFactor * trgEff.GetVal();
            //if (1) MCemuScaleFactor = MCScaleFactor * trgEff.GetVal();
            //else MCemuScaleFactor = MCScaleFactor;
            MCemuScaleFactor = MCScaleFactor;

            if (lowMassPuOnly && invMass > puMassCut) weight = 1.;
            if (p > DATA) weight *= MCemuScaleFactor;

            int jetsPt20 = 0;
            int jetsPt30 = 0;
            for (int jetIt = 0; jetIt < c_jetColl_size; ++ jetIt) {
               if (c_jetPt[jetIt] > 20.) ++jetsPt20;
               if (c_jetPt[jetIt] > 30.) ++jetsPt30;
            }
            //if (c_jetColl_size > 10) continue;

            // fill the data tree
            emuInvMass = invMass;
            eCharge = c_gsf_charge[GSF_passHEEP[0]];
            muCharge = c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]];
            totWeight = weight;
            if (p > DATA) totWeight *= input[p].second * LumiFactor;
            if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 1.442) evtRegion = 0;
            else if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 2.5) evtRegion = 1;
            else evtRegion = -1;
            // fill control variables
            nVtx = c_pvsize;
            dZFstPVtx = fabs(c_gsf_dz_firstPVtx[GSF_passHEEP[0]] - c_muon_dz_firstPVtx[MU_passGOOD[MU_leadingPassGOOD]]);
            dXYFstPVtx = fabs(c_gsf_dxy_firstPVtx[GSF_passHEEP[0]] - c_muon_dxy_firstPVtx[MU_passGOOD[MU_leadingPassGOOD]]);
            if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) < 3.14) dPhi = fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]);
            if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) > 3.14) dPhi = 6.28 - fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]);
            dEta = fabs(c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_eta[GSF_passHEEP[0]]);
            eleEt = c_gsf_gsfet[GSF_passHEEP[0]];
            eleEta = c_gsf_eta[GSF_passHEEP[0]];
            elePhi = c_gsf_phi[GSF_passHEEP[0]];
            eleDEta = c_gsf_deltaeta[GSF_passHEEP[0]];
            eleDPhi = c_gsf_deltaphi[GSF_passHEEP[0]];
            eleHOE = c_gsf_hovere[GSF_passHEEP[0]];
            eleSigmaIEIE = c_gsf_sigmaIetaIeta[GSF_passHEEP[0]];
            eleEcalIso = c_gsf_ecaliso[GSF_passHEEP[0]];
            eleHcalIso12 = c_gsf_hcaliso1[GSF_passHEEP[0]] + c_gsf_hcaliso2[GSF_passHEEP[0]];
            eleTrkIso = c_gsf_trackiso[GSF_passHEEP[0]];
            eleLostHits = c_gsf_nLostInnerHits[GSF_passHEEP[0]];
            muIsoCombRel = CombRelIso;
            muEtEleOPtMu = c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];
            if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) muPtPlusOPtMinus = c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];
            if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) muPtPlusOPtMinus = c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]/c_gsf_gsfet[GSF_passHEEP[0]];
            muPt = c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];
            muEta = c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]];
            muPhi = c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]];
            muHitLayers = c_muon_nlayerswithhits[MU_passGOOD[MU_leadingPassGOOD]];
            muPxlHits = c_muon_nhitspixel[MU_passGOOD[MU_leadingPassGOOD]];
            muMuHits = c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]];
            muDZFstPVtx = c_muon_dz_firstPVtx[MU_passGOOD[MU_leadingPassGOOD]];
            muDXYFstPVtx = c_muon_dxy_firstPVtx[MU_passGOOD[MU_leadingPassGOOD]];
            muNSeg = c_muon_nSegmentMatch[MU_passGOOD[MU_leadingPassGOOD]];
            muTrkIso03 = c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]];
            numOfJets = c_jetColl_size;
            numOfJetsPt20 = jetsPt20;
            numOfJetsPt30 = jetsPt30;

            emuTree->Fill();

            //LIKE SIGN OPPOSITE SIGN
            if (p == DATA) {
               if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) ++nb_plus_plus;
               if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) ++nb_plus_minus;
               if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) ++nb_minus_plus;
               if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) ++nb_minus_minus;

               if (invMass > 600. || invMass < 1.) {
                  if (invMass > 600.) cout << "M_emu > 600GeV/c^2 event | " << c_runnumber << ":" << c_luminosityBlock << ":" << c_eventnumber << " | "; 
                  if (invMass < 1.) cout << "M_emu < 1GeV/c^2 event   | " << c_runnumber << ":" << c_luminosityBlock << ":" << c_eventnumber << " | "; 
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
      } // END LOOP OVER ENTRIES

      cout << "Number of selected events: " << evCounter << endl;

      if (p == DATA) {
         dataEntries = dataTrig[0] + dataTrig[1] + dataTrig[2];
         cout << "HLT_Mu22_Photon22_CaloIdL: " << dataTrig[0] 
              << " , HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << dataTrig[1] 
              << " , HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << dataTrig[2] 
              << " , sum: " << dataTrig[0]+dataTrig[1]+dataTrig[2] << endl;
         cout << "Runrange HLT_Mu22_Photon22_CaloIdL: " << runs_HLT_Mu22_Photon22_CaloIdL.first 
              << " - " << runs_HLT_Mu22_Photon22_CaloIdL.second << endl;
         cout << "Runrange HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " 
              << runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first << " - " 
              << runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second << endl;
         cout << "Runrange HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " 
              << runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first << " - " 
              << runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second << endl;
      }
      else {
         cout << "HLT_Mu22_Photon22_CaloIdL: " << trig[0] 
                << " , HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << trig[1] 
                << " , HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << trig[2] 
                << " , sum: " << trig[0]+trig[1]+trig[2] << endl;
      }
      output->Cd("");

      TParameter<float> *mcWeight = new TParameter<float>((const char *)suffix[p], input[p].second);
      mcWeights.Add(mcWeight);

      emuTree->Write();
      puMc->Write();
      puMcNorm->Write();
      puWeights->Write();
   }//END FILE LOOP
   mcWeights.Write("mcWeights", TObject::kSingleKey);
   puData->Write();
   puDataNorm->Write();

   //PRINT INTEGRAL

   cout << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   cout << "HEEP - TIGHT MU        Lumi        = " << LumiFactor << "pb-1" << endl;
   cout << "                       e pT EB     > " << bar_et << "GeV/c" << endl;
   cout << "                       e pT EE     > " << end_et << "GeV/c" << endl;
   cout << "                       mu pT       > " << muon_pt << "GeV/c" << endl;
   cout << "                       mu |eta|    < " << muon_etaMax << endl;
   cout << endl;
   cout << "Systematic errors" << endl;
   cout << "(xsec errors from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV)" << endl;
   cout << " Luminosity:  " << systErrLumi.GetVal() * 100 << "%" << endl;
   cout << " Efficiency:  " << systErrEff.GetVal() * 100 << "%" << endl;
   cout << " ttbar:      " <<  ((TParameter<float> *)systErrMCs.FindObject("systErrMcTtbar"))->GetVal() * 100 << "%" << endl;
   cout << " Z->tautau:   " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcDyTauTau"))->GetVal() * 100 << "%" << endl;
   cout << " WW:          " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcWW"))->GetVal() * 100 << "%" << endl;
   cout << " WZ:          " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcWZ"))->GetVal() * 100 << "%" << endl;
   cout << " tW, tbarW:   " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcTW"))->GetVal() * 100 << "%" << endl;
   cout << " W+Jets:      " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcWJets"))->GetVal() * 100 << "%" << endl;
   cout << " Z->mumu:     " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcDyMuMu"))->GetVal() * 100 << "%" << endl;
   cout << " Z->ee:       " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcDyEE"))->GetVal() * 100 << "%" << endl;
   cout << " ZZ:          " << ((TParameter<float> *)systErrMCs.FindObject("systErrMcZZ"))->GetVal() * 100 << "%" << endl;

   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   cout << "nb e+mu+     = " << nb_plus_plus << endl;
   cout << "nb e+mu-     = " << nb_plus_minus << endl;
   cout << "nb e-mu+     = " << nb_minus_plus << endl;
   cout << "nb e-mu-     = " << nb_minus_minus << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;

   output->Close();
   timer.Stop();
   timer.Print();
}

int Trigger(TFile *inFile, unsigned int &entry, int &prescale, unsigned int *trig, const int &selector)
{
   int c_HLT_Mu22_Photon22_CaloIdL;
   int c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   int c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   int c_prescale_HLT_Mu22_Photon22_CaloIdL;
   int c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   int c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;

   TBranch *b_HLT_Mu22_Photon22_CaloIdL;
   TBranch *b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   TBranch *b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   TBranch *b_prescale_HLT_Mu22_Photon22_CaloIdL;
   TBranch *b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   TBranch *b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;

   // Get the tree and connect the necessary variables
   inFile->Cd("");
   TTree *thetree;
   thetree = (TTree*)(inFile)->Get("gsfcheckerjob/tree");

   thetree->SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &c_HLT_Mu22_Photon22_CaloIdL, &b_HLT_Mu22_Photon22_CaloIdL);
   thetree->SetBranchAddress("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   thetree->SetBranchAddress("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   thetree->SetBranchAddress("prescale_HLT_Mu22_Photon22_CaloIdL", &c_prescale_HLT_Mu22_Photon22_CaloIdL, &b_prescale_HLT_Mu22_Photon22_CaloIdL);
   thetree->SetBranchAddress("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
   thetree->SetBranchAddress("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, &b_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);

   thetree->GetEntry(entry);

   // switch off trigger
   //if (selector < 2) trig[0]++;
   //return 1;

   if (selector == 1) {
      prescale = c_prescale_HLT_Mu22_Photon22_CaloIdL;
      if (c_HLT_Mu22_Photon22_CaloIdL >= 0) trig[0]++;
      return c_HLT_Mu22_Photon22_CaloIdL;
   } else if (selector == 2) {
      prescale = c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      if (c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL >= 0) trig[1]++;
      return c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   } else if (selector == 3) {
      prescale = c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      if (c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL >= 0) trig[2]++;
      return c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   } else {
      // select an unprescaled trigger
      if (c_prescale_HLT_Mu22_Photon22_CaloIdL == 1) {
         prescale = c_prescale_HLT_Mu22_Photon22_CaloIdL;
         trig[0]++;
         if (c_runnumber < runs_HLT_Mu22_Photon22_CaloIdL.first) runs_HLT_Mu22_Photon22_CaloIdL.first = c_runnumber;
         if (c_runnumber > runs_HLT_Mu22_Photon22_CaloIdL.second) runs_HLT_Mu22_Photon22_CaloIdL.second = c_runnumber;
         return c_HLT_Mu22_Photon22_CaloIdL;
      } else if (c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 1) {
         prescale = c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
         trig[1]++;
         if (c_runnumber < runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first) runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first = c_runnumber;
         if (c_runnumber > runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second) runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second = c_runnumber;
         return c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      } else if (c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 1) {
         prescale = c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
         trig[2]++;
         if (c_runnumber < runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first) runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first = c_runnumber;
         if (c_runnumber > runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second) runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second = c_runnumber;
         return c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      }
   }

   if (c_prescale_HLT_Mu22_Photon22_CaloIdL == 0 && c_prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 0 && c_prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 0) {
      //cout << "Prescale alert! All selecting triggers have prescale 0. Event: " << c_runnumber << ":" << c_luminosityBlock << ":" << c_eventnumber << " Triggers: " << c_HLT_Mu22_Photon22_CaloIdL << c_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL << c_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL << endl;
   } else {
      //cout << "Prescale alert! No unprescaled trigger found. Event: " << c_runnumber << ":" << c_luminosityBlock << ":" << c_eventnumber << endl;
   }
   prescale = 0;
   return -1;
}

float CorrectEnergy(TFile *inFile, int &pos) {
   float c_gsfsc_e[100];
   float c_gsf_theta[100];

   TBranch *b_gsfsc_e;
   TBranch *b_gsf_theta;

   // Get the tree and connect the necessary variables
   inFile->Cd("");
   TTree *thetree;
   thetree = (TTree*)(inFile)->Get("gsfcheckerjob/tree");

   thetree->SetBranchAddress("gsfsc_e", &c_gsfsc_e, &b_gsfsc_e);
   thetree->SetBranchAddress("gsf_theta", &c_gsf_theta, &b_gsf_theta);

   return c_gsfsc_e[pos] * sin(c_gsf_theta[pos]);
}
