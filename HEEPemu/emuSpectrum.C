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

   int nPVtxMax = 60; // maximum number of primary vertices

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

   int nBins = 75;
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
   float muon_impactParamMaxZ = 0.5;   // in cm , not used
   int muon_nSegMatchMin = 2;
   float muon_relIsoCutMax = 0.1;

   //BARREL
   float bar_et = 0.;
   float bar_hoE = 0.;
   float bar_coshFactor = 0.;
   float bar_coshOffset = 0.;
   float bar_etaFactor = 0.;
   float bar_DEta = 0.;
   float bar_DPhi = 0.;
   float bar_e2x5e5x5 = 0.;
   float bar_e2x5e5x5Rho = 0.;
   float bar_e1x5e5x5 = 0.;
   float bar_e1x5e5x5Rho = 0.;
   float bar_isoEcalHcal1_1 = 0.;
   float bar_isoEcalHcal1_2 = 0.;
   float bar_isoEcalHcalRho = 0.;
   float bar_isoTrack = 0.;
   float bar_dxy = 0.;
   int bar_missInnerHits = 0;

   //ENDCAP
   float end_et = 0.;
   float end_hoE = 0.;
   float end_coshFactor = 0.;
   float end_coshOffset = 0.;
   float end_etaFactor = 0.;
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
         bar_coshFactor = 0.13;
         bar_coshOffset = 0.085;
         bar_etaFactor = 0.41;
         bar_DEta = 0.005;
         bar_DPhi = 0.06;
         bar_e2x5e5x5 = 0.94;
         bar_e2x5e5x5Rho = 0.0054;
         bar_e1x5e5x5 = 0.83;
         bar_e1x5e5x5Rho = 0.0045;
         bar_isoEcalHcal1_1 = 2.;
         bar_isoEcalHcal1_2 = 0.03;
         bar_isoEcalHcalRho = 0.28;
         bar_isoTrack = 5.;
         bar_missInnerHits = 0;
     
         //ENDCAP
         end_et = 35.;
         end_hoE = 0.05;
         end_coshFactor = 0.13;
         end_coshOffset = 0.085;
         end_etaFactor = 0.41;
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

   TString sign[3] = {"", "LS_", "OS_"};

   // counting variables
   int nb_plus_plus = 0;
   int nb_plus_minus = 0;
   int nb_minus_plus = 0;
   int nb_minus_minus = 0;

   // histogram containers
   vector<vector<TH1F *> > h_emuMass;
   vector<vector<TH1F *> > h_emuMassEB;
   vector<vector<TH1F *> > h_emuMassEE;
   vector<vector<TH1F *> > h_emu_mass_accVSgood;

   vector<vector<TH1F *> > h_pfMet;
   vector<vector<TH1F *> > h_nVtx;
   vector<vector<TH1F *> > h_dDz;
   vector<vector<TH1F *> > h_dDzBeamSpot;
   vector<vector<TH1F *> > h_dDzFirstPVtx;
   vector<vector<TH1F *> > h_rho;
   vector<vector<TH1F *> > h_nJets;
   vector<vector<TH1F *> > h_nJetsPt20;
   vector<vector<TH1F *> > h_nJetsPt30;
   vector<vector<TH1F *> > h_dPhi;
   vector<vector<TH1F *> > h_eleEt;
   vector<vector<TH1F *> > h_eleEta;
   vector<vector<TH1F *> > h_elePhi;
   vector<vector<TH1F *> > h_eleDEta;
   vector<vector<TH1F *> > h_eleDPhi;
   vector<vector<TH1F *> > h_eleHOE;
   vector<vector<TH1F *> > h_eleSigmaIEIE;
   vector<vector<TH1F *> > h_eleEcalIso;
   vector<vector<TH1F *> > h_eleHcalIso12;
   vector<vector<TH1F *> > h_eleTrkIso;
   vector<vector<TH1F *> > h_eleLostHits;

   vector<vector<TH1F *> > h_muIsoCombRel;
   vector<vector<TH1F *> > h_muEtEleOPtMu;
   vector<vector<TH1F *> > h_muPtPlusOPtMinus;
   vector<vector<TH1F *> > h_muPt;
   vector<vector<TH1F *> > h_muEta;
   vector<vector<TH1F *> > h_muPhi;
   vector<vector<TH1F *> > h_muHitLayers;
   vector<vector<TH1F *> > h_muPxlHits;
   vector<vector<TH1F *> > h_muMuHits;
   vector<vector<TH1F *> > h_muDxyBeamSpot;
   vector<vector<TH1F *> > h_muNSeg;
   vector<vector<TH1F *> > h_muTrkIso03;

   //vector<TH1F *> emu_plus_plus;
   //vector<TH1F *> emu_plus_minus;
   //vector<TH1F *> emu_minus_plus;
   //vector<TH1F *> emu_minus_minus;

   //
   TH1F *emu_dilepton = new TH1F("emu_dilepton", "emu_dilepton", nBins, 0., 1500.);
   TH1F *emu_ewk = new TH1F("emu_ewk", "emu_ewk", nBins, 0., 1500.);
   TH1F *emu_jet = new TH1F("emu_jet", "emu_jet", nBins, 0., 1500.);

   // data pileup histogram
   TFile *dataPuInfile = TFile::Open(puFile);
   dataPuInfile->Cd("");
   TH1F *puData = (TH1F *)gDirectory->Get("pileup");
   puData->SetDirectory(0);
   puData->SetName("puData");
   dataPuInfile->Close();
   TH1F *puDataNorm = (TH1F *)puData->DrawNormalized()->Clone("puDataNorm");

   //GLOBAL
   float c_calomet;
   float c_met;
   float c_pfmet;
   float c_pthat;
   float c_bsposx;
   float c_bsposy;
   float c_bsposz;
   int c_jetColl_size;
   float c_jetPt[100];

   //PRIM VTX
   int c_pvsize;
   int c_pvz[100];
   bool c_pv_isValid[100];
   int c_pv_ndof[100];
   int c_pv_nTracks[100];
   float c_pv_normChi2[100];
   int c_pv_totTrackSize[100];

   float c_rho;
   int c_trueNVtx;

   //GSF
   int c_gsf_size;
   float c_gsf_gsfet[100];
   float c_gsf_px[100];
   float c_gsf_py[100];
   float c_gsf_pz[100];
   float c_gsf_pt[100];
   float c_gsf_eta[100];
   float c_gsf_theta[100];
   float c_gsf_phi[100];
   float c_gsf_dxy[100];
   float c_gsf_dxy_beamSpot[100];
   float c_gsf_dxy_firstPVtx[100];
   float c_gsf_dz[100];
   float c_gsf_dz_beamSpot[100];
   float c_gsf_dz_firstPVtx[100];
   bool c_gsf_isecaldriven[100];
   int c_gsf_charge[100];
   float c_gsf_deltaeta[100];
   float c_gsf_deltaphi[100];
   float c_gsf_e5x5[100];
   float c_gsf_e1x5overe5x5[100];
   float c_gsf_e2x5overe5x5[100];
   float c_gsf_sigmaetaeta[100];
   float c_gsf_sigmaIetaIeta[100];
   float c_gsf_hovere[100];
   int c_gsf_nHits[100];
   int c_gsf_nLostInnerHits[100];
   float c_gsf_ecaliso[100];
   float c_gsf_hcaliso1[100];
   float c_gsf_hcaliso2[100];
   float c_gsf_trackiso[100];

   float c_gsfsc_e[100];
   float c_gsfsc_eta[100];
   float c_gsfsc_phi[100];

   bool c_gsfpass_HEEP[100];

   //MUONS
   int c_muon_size;
   float c_muon_pt[100];
   float c_muon_eta[100];
   float c_muon_phi[100];
   float c_muon_theta[100];
   float c_muon_ptError[100];
   float c_muon_etaError[100];
   float c_muon_phiError[100];
   float c_muon_thetaError[100];
   float c_muon_outerPt[100];
   float c_muon_outerEta[100];
   float c_muon_outerPhi[100];
   float c_muon_outerTheta[100];
   float c_muon_px[100];
   float c_muon_py[100];
   float c_muon_pz[100];
   int c_muon_charge[100];
   int c_muon_nhitstrack[100];
   int c_muon_nhitspixel[100];
   int c_muon_nhitstotal[100];
   int c_muon_nhitsmuons[100];
   int c_muon_nlayerswithhits[100];
   int c_muon_nSegmentMatch[100];
   bool c_muon_isTrackerMuon[100];
   float c_muon_normChi2[100];
   float c_muon_dzError[100];
   float c_muon_dxyError[100];
   float c_muon_dz_cmsCenter[100];
   float c_muon_dz_beamSpot[100];
   float c_muon_dz_firstPVtx[100];
   float c_muon_dxy_cmsCenter[100];
   float c_muon_dxy_beamSpot[100];
   float c_muon_dxy_firstPVtx[100];
   float c_muon_trackIso03[100];
   float c_muon_emIso03[100];
   float c_muon_hadIso03[100];
   float c_muon_trackIso03_ptInVeto[100];
   float c_muon_emIso03_ptInVeto[100];
   float c_muon_hadIso03_ptInVeto[100];

   float c_genPair_mass;

   //RUN ID
   TBranch *b_runnumber;
   TBranch *b_eventnumber;
   TBranch *b_luminosityBlock;

   //GLOBAL
   TBranch *b_calomet;
   TBranch *b_met;
   TBranch *b_pfmet;
   TBranch *b_pthat;
   TBranch *b_bsposx;
   TBranch *b_bsposy;
   TBranch *b_bsposz;
   TBranch *b_jetColl_size;
   TBranch *b_jetPt;

   //PRIM VTX
   TBranch *b_pvsize;
   TBranch *b_pvz;
   TBranch *b_pv_isValid;
   TBranch *b_pv_ndof;
   TBranch *b_pv_nTracks;
   TBranch *b_pv_normChi2;
   TBranch *b_pv_totTrackSize;

   TBranch *b_rho;
   TBranch *b_trueNVtx;

   //GSF
   TBranch *b_gsf_size;
   TBranch *b_gsf_gsfet;
   TBranch *b_gsf_px;
   TBranch *b_gsf_py;
   TBranch *b_gsf_pz;
   TBranch *b_gsf_pt;
   TBranch *b_gsf_eta;
   TBranch *b_gsf_theta;
   TBranch *b_gsf_phi;
   TBranch *b_gsf_dxy;
   TBranch *b_gsf_dxy_beamSpot;
   TBranch *b_gsf_dxy_firstPVtx;
   TBranch *b_gsf_dz;
   TBranch *b_gsf_dz_beamSpot;
   TBranch *b_gsf_dz_firstPVtx;
   TBranch *b_gsf_isecaldriven;
   TBranch *b_gsf_charge;
   TBranch *b_gsf_deltaeta;
   TBranch *b_gsf_deltaphi;
   TBranch *b_gsf_e5x5;
   TBranch *b_gsf_e1x5overe5x5;
   TBranch *b_gsf_e2x5overe5x5;
   TBranch *b_gsf_sigmaetaeta;
   TBranch *b_gsf_sigmaIetaIeta;
   TBranch *b_gsf_hovere;
   TBranch *b_gsf_nHits;
   TBranch *b_gsf_nLostInnerHits;
   TBranch *b_gsf_ecaliso;
   TBranch *b_gsf_hcaliso1;
   TBranch *b_gsf_hcaliso2;
   TBranch *b_gsf_trackiso;

   TBranch *b_gsfsc_e;
   TBranch *b_gsfsc_eta;
   TBranch *b_gsfsc_phi;

   TBranch *b_gsfpass_HEEP;

   //MUONS
   TBranch *b_muon_size;
   TBranch *b_muon_pt;
   TBranch *b_muon_eta;
   TBranch *b_muon_phi;
   TBranch *b_muon_theta;
   TBranch *b_muon_ptError;
   TBranch *b_muon_etaError;
   TBranch *b_muon_phiError;
   TBranch *b_muon_thetaError;
   TBranch *b_muon_outerPt;
   TBranch *b_muon_outerEta;
   TBranch *b_muon_outerPhi;
   TBranch *b_muon_outerTheta;
   TBranch *b_muon_px;
   TBranch *b_muon_py;
   TBranch *b_muon_pz;
   TBranch *b_muon_charge;
   TBranch *b_muon_nhitstrack;
   TBranch *b_muon_nhitspixel;
   TBranch *b_muon_nhitstotal;
   TBranch *b_muon_nhitsmuons;
   TBranch *b_muon_nlayerswithhits;
   TBranch *b_muon_nSegmentMatch;
   TBranch *b_muon_isTrackerMuon;
   TBranch *b_muon_normChi2;
   TBranch *b_muon_dz_cmsCenter;
   TBranch *b_muon_dz_beamSpot;
   TBranch *b_muon_dz_firstPVtx;
   TBranch *b_muon_dxy_cmsCenter;
   TBranch *b_muon_dxy_beamSpot;
   TBranch *b_muon_dxy_firstPVtx;
   TBranch *b_muon_dzError;
   TBranch *b_muon_dxyError;
   TBranch *b_muon_trackIso03;
   TBranch *b_muon_emIso03;
   TBranch *b_muon_hadIso03;
   TBranch *b_muon_trackIso03_ptInVeto;
   TBranch *b_muon_emIso03_ptInVeto;
   TBranch *b_muon_hadIso03_ptInVeto;

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
      thetree->SetBranchAddress("calomet", &c_calomet, &b_calomet);
      thetree->SetBranchAddress("met", &c_met, &b_met);
      thetree->SetBranchAddress("pfmet", &c_pfmet, &b_pfmet);
      thetree->SetBranchAddress("pthat", &c_pthat, &b_pthat);
      thetree->SetBranchAddress("bsposx", &c_bsposx, &b_bsposx);
      thetree->SetBranchAddress("bsposy", &c_bsposy, &b_bsposy);
      thetree->SetBranchAddress("bsposz", &c_bsposz, &b_bsposz);
      thetree->SetBranchAddress("JetColl_size", &c_jetColl_size, &b_jetColl_size);
      thetree->SetBranchAddress("Jet_pt", &c_jetPt, &b_jetPt);

      //PRIM VTX
      thetree->SetBranchAddress("pvsize", &c_pvsize, &b_pvsize);
      thetree->SetBranchAddress("pvz", &c_pvz, &b_pvz);
      thetree->SetBranchAddress("pv_isValid", &c_pv_isValid, &b_pv_isValid);
      thetree->SetBranchAddress("pv_ndof", &c_pv_ndof, &b_pv_ndof);
      thetree->SetBranchAddress("pv_nTracks", &c_pv_nTracks, &b_pv_nTracks);
      thetree->SetBranchAddress("pv_normChi2", &c_pv_normChi2, &b_pv_normChi2);
      thetree->SetBranchAddress("pv_totTrackSize", &c_pv_totTrackSize, &b_pv_totTrackSize);

      thetree->SetBranchAddress("rho", &c_rho, &b_rho);
      thetree->SetBranchAddress("trueNVtx", &c_trueNVtx, &b_trueNVtx);

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
      thetree->SetBranchAddress("gsf_dxy", &c_gsf_dxy, &b_gsf_dxy);
      thetree->SetBranchAddress("gsf_dxy_beamSpot", &c_gsf_dxy_beamSpot, &b_gsf_dxy_beamSpot);
      thetree->SetBranchAddress("gsf_dxy_firstPVtx", &c_gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
      thetree->SetBranchAddress("gsf_dz", &c_gsf_dz, &b_gsf_dz);
      thetree->SetBranchAddress("gsf_dz_beamSpot", &c_gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
      thetree->SetBranchAddress("gsf_dz_firstPVtx", &c_gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
      thetree->SetBranchAddress("gsf_isecaldriven", &c_gsf_isecaldriven, &b_gsf_isecaldriven);
      thetree->SetBranchAddress("gsf_charge", &c_gsf_charge, &b_gsf_charge);
      thetree->SetBranchAddress("gsf_deltaeta", &c_gsf_deltaeta, &b_gsf_deltaeta);
      thetree->SetBranchAddress("gsf_deltaphi", &c_gsf_deltaphi, &b_gsf_deltaphi);
      thetree->SetBranchAddress("gsf_e5x5", &c_gsf_e5x5, &b_gsf_e5x5);
      thetree->SetBranchAddress("gsf_e1x5overe5x5", &c_gsf_e1x5overe5x5, &b_gsf_e1x5overe5x5);
      thetree->SetBranchAddress("gsf_e2x5overe5x5", &c_gsf_e2x5overe5x5, &b_gsf_e2x5overe5x5);
      thetree->SetBranchAddress("gsf_sigmaetaeta", &c_gsf_sigmaetaeta, &b_gsf_sigmaetaeta);
      thetree->SetBranchAddress("gsf_sigmaIetaIeta", &c_gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
      thetree->SetBranchAddress("gsf_hovere", &c_gsf_hovere, &b_gsf_hovere);
      thetree->SetBranchAddress("gsf_nHits", &c_gsf_nHits, &b_gsf_nHits);
      thetree->SetBranchAddress("gsf_nLostInnerHits", &c_gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
      thetree->SetBranchAddress("gsf_ecaliso", &c_gsf_ecaliso, &b_gsf_ecaliso);
      thetree->SetBranchAddress("gsf_hcaliso1", &c_gsf_hcaliso1, &b_gsf_hcaliso1);
      thetree->SetBranchAddress("gsf_hcaliso2", &c_gsf_hcaliso2, &b_gsf_hcaliso2);
      thetree->SetBranchAddress("gsf_trackiso", &c_gsf_trackiso, &b_gsf_trackiso);
      thetree->SetBranchAddress("gsfsc_e", &c_gsfsc_e, &b_gsfsc_e);
      thetree->SetBranchAddress("gsfsc_eta", &c_gsfsc_eta, &b_gsfsc_eta);
      thetree->SetBranchAddress("gsfsc_phi", &c_gsfsc_phi, &b_gsfsc_phi);
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
      thetree->SetBranchAddress("muon_normChi2", &c_muon_normChi2, &b_muon_normChi2);
      thetree->SetBranchAddress("muon_dz_cmsCenter", &c_muon_dz_cmsCenter, &b_muon_dz_cmsCenter);
      thetree->SetBranchAddress("muon_dz_beamSpot", &c_muon_dz_beamSpot, &b_muon_dz_beamSpot);
      thetree->SetBranchAddress("muon_dz_firstPVtx", &c_muon_dz_firstPVtx, &b_muon_dz_firstPVtx);
      thetree->SetBranchAddress("muon_dxy_cmsCenter", &c_muon_dxy_cmsCenter, &b_muon_dxy_cmsCenter);
      thetree->SetBranchAddress("muon_dxy_beamSpot", &c_muon_dxy_beamSpot, &b_muon_dxy_beamSpot);
      thetree->SetBranchAddress("muon_dxy_firstPVtx", &c_muon_dxy_firstPVtx, &b_muon_dxy_firstPVtx);
      thetree->SetBranchAddress("muon_dzError", &c_muon_dzError, &b_muon_dzError);
      thetree->SetBranchAddress("muon_dxyError", &c_muon_dxyError, &b_muon_dxyError);
      thetree->SetBranchAddress("muon_trackIso03", &c_muon_trackIso03, &b_muon_trackIso03);
      thetree->SetBranchAddress("muon_emIso03", &c_muon_emIso03, &b_muon_emIso03);
      thetree->SetBranchAddress("muon_hadIso03", &c_muon_hadIso03, &b_muon_hadIso03);
      thetree->SetBranchAddress("muon_trackIso03_ptInVeto", &c_muon_trackIso03_ptInVeto, &b_muon_trackIso03_ptInVeto);
      thetree->SetBranchAddress("muon_emIso03_ptInVeto", &c_muon_emIso03_ptInVeto, &b_muon_emIso03_ptInVeto);
      thetree->SetBranchAddress("muon_hadIso03_ptInVeto", &c_muon_hadIso03_ptInVeto, &b_muon_hadIso03_ptInVeto);
      
      if (storeGenMTtbar[p]) thetree->SetBranchAddress("genPair_mass", &c_genPair_mass, &b_genPair_mass);

      // get the histogram with the true number of vertices
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

      // set up histograms
      vector<TH1F *> helper;
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMass_" + sign[k] + suffix[p], "emuMass_" + sign[k] + suffix[p], nBins, 0., 1500.));
      h_emuMass.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMassEB_" + sign[k] + suffix[p], "emuMassEB_" + sign[k] + suffix[p], nBins, 0., 1500.));
      h_emuMassEB.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emuMassEE_" + sign[k] + suffix[p], "emuMassEE_" + sign[k] + suffix[p], nBins, 0., 1500.));
      h_emuMassEE.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("emu_mass_accVSgood_" + sign[k] + suffix[p], "emu_mass_accVSgood_" + sign[k] + suffix[p], nBins, 0., 1500.));
      h_emu_mass_accVSgood.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("PFMET_" + sign[k] + suffix[p], "PFMET_" + sign[k] + suffix[p], 50, 0., 500.));
      h_pfMet.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("nVtx_" + sign[k] + suffix[p], "nVtx_" + sign[k] + suffix[p], nPVtxMax, 0., nPVtxMax));
      h_nVtx.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("dDz_" + sign[k] + suffix[p], "dDz_" + sign[k] + suffix[p], 50, 0., 2.));
      h_dDz.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("dDzBeamSpot_" + sign[k] + suffix[p], "dDzBeamSpot_" + sign[k] + suffix[p], 100, 0., 0.2));
      h_dDzBeamSpot.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("dDzFirstPVtx_" + sign[k] + suffix[p], "dDzFirstPVtx_" + sign[k] + suffix[p], 100, 0., 0.2));
      h_dDzFirstPVtx.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("rho_" + sign[k] + suffix[p], "rho_" + sign[k] + suffix[p], 50, 0., 50.));
      h_rho.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("dPhi_" + sign[k] + suffix[p], "dPhi_" + sign[k] + suffix[p], 64, 0., 3.2));
      h_dPhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleEt_" + sign[k] + suffix[p], "eleEt_" + sign[k] + suffix[p], 50, 0., 500.));
      h_eleEt.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleEta_" + sign[k] + suffix[p], "eleEta_" + sign[k] + suffix[p], 25, -2.5, 2.5));
      h_eleEta.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("elePhi_" + sign[k] + suffix[p], "elePhi_" + sign[k] + suffix[p], 25, -3.14, 3.14));
      h_elePhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleDEta_" + sign[k] + suffix[p], "eleDEta_" + sign[k] + suffix[p], 50, -0.008, 0.008));
      h_eleDEta.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleDPhi_" + sign[k] + suffix[p], "eleDPhi_" + sign[k] + suffix[p], 50, -0.06, 0.06));
      h_eleDPhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleHOE_" + sign[k] + suffix[p], "eleHOE_" + sign[k] + suffix[p], 100, 0., 0.06));
      h_eleHOE.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleSigmaIEIE_" + sign[k] + suffix[p], "eleSigmaIEIE_" + sign[k] + suffix[p], 50, 0., 0.04));
      h_eleSigmaIEIE.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleEcalIso_" + sign[k] + suffix[p], "eleEcalIso_" + sign[k] + suffix[p], 45, -10., 35.));
      h_eleEcalIso.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleHcalIso12_" + sign[k] + suffix[p], "eleHcalIso12_" + sign[k] + suffix[p], 20, 0., 20.));
      h_eleHcalIso12.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleTrkIso_" + sign[k] + suffix[p], "eleTrkIso_" + sign[k] + suffix[p], 50, 0., 5.));
      h_eleTrkIso.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("eleLostHits_" + sign[k] + suffix[p], "eleLostHits_" + sign[k] + suffix[p], 3, 0., 3.));
      h_eleLostHits.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muIsoCombRel_" + sign[k] + suffix[p], "muIsoCombRel_" + sign[k] + suffix[p], 50, 0., 0.2));
      h_muIsoCombRel.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muEtEleOPtMu_" + sign[k] + suffix[p], "muEtEleOPtMu_" + sign[k] + suffix[p], 50, 0., 4.0));
      h_muEtEleOPtMu.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPtPlusOPtMinus_" + sign[k] + suffix[p], "muPtPlusOPtMinus_" + sign[k] + suffix[p], 50, 0., 4.0));
      h_muPtPlusOPtMinus.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPt_" + sign[k] + suffix[p], "muPt_" + sign[k] + suffix[p], 50, 0., 500.));
      h_muPt.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muEta_" + sign[k] + suffix[p], "muEta_" + sign[k] + suffix[p], 25, -2.5, 2.5));
      h_muEta.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPhi_" + sign[k] + suffix[p], "muPhi_" + sign[k] + suffix[p], 25, -3.14, 3.14));
      h_muPhi.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muHitLayers_" + sign[k] + suffix[p], "muHitLayers_" + sign[k] + suffix[p], 20, 0., 20.));
      h_muHitLayers.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muPxlHits_" + sign[k] + suffix[p], "muPxlHits_" + sign[k] + suffix[p], 10, 0., 10));
      h_muPxlHits.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muMuHits_" + sign[k] + suffix[p], "muMuHits_" + sign[k] + suffix[p], 12, 0., 12));
      h_muMuHits.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muDxyBeamSpot_" + sign[k] + suffix[p], "muDxyBeamSpot_" + sign[k] + suffix[p], 200, -0.2, 0.2));
      h_muDxyBeamSpot.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muNSeg_" + sign[k] + suffix[p], "muNSeg_" + sign[k] + suffix[p], 10, 0., 10.));
      h_muNSeg.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("muTrkIso03_" + sign[k] + suffix[p], "muTrkIso03_" + sign[k] + suffix[p], 50, 0., 20.));
      h_muTrkIso03.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("numOfJets_" + sign[k] + suffix[p], "numOfJets_" + sign[k] + suffix[p], 30, 0., 30.));
      h_nJets.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("numOfJetsPt20_" + sign[k] + suffix[p], "numOfJetsPt20_" + sign[k] + suffix[p], 30, 0., 30.));
      h_nJetsPt20.push_back(helper);
      helper.clear();
      for (unsigned int k = 0; k < 3; ++k) helper.push_back(new TH1F("numOfJetsPt30_" + sign[k] + suffix[p], "numOfJetsPt30_" + sign[k] + suffix[p], 30, 0., 30.));
      h_nJetsPt30.push_back(helper);

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
                   && (c_gsf_hovere[j] * c_gsf_gsfet[j]) < (bar_hoE * c_gsf_gsfet[j] + (bar_coshFactor * cosh(bar_etaFactor * c_gsfsc_eta[j]) - bar_coshOffset) * c_rho)
                   && fabs(c_gsf_deltaeta[j]) < bar_DEta
                   && fabs(c_gsf_deltaphi[j]) < bar_DPhi
                   && ((c_gsf_e2x5overe5x5[j] * c_gsf_e5x5[j] / cosh(c_gsfsc_eta[j])) > (bar_e2x5e5x5 * c_gsf_e5x5[j] / cosh(c_gsfsc_eta[j]) - bar_e2x5e5x5Rho * c_rho) 
                      || (c_gsf_e1x5overe5x5[j] * c_gsf_e5x5[j] / cosh(c_gsfsc_eta[j])) > (bar_e1x5e5x5 * c_gsf_e5x5[j] / cosh(c_gsfsc_eta[j]) - bar_e1x5e5x5Rho * c_rho))
                   && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * c_gsf_gsfet[j] + bar_isoEcalHcalRho * c_rho)
                   && c_gsf_trackiso[j] < bar_isoTrack
                   && c_gsf_nLostInnerHits[j] <= bar_missInnerHits
                  ) if ((selection == 3 && fabs(c_gsf_dxy_firstPVtx[j]) <= bar_dxy) || selection != 3) GSF_passHEEP.push_back(j);
   
               //ENDCAP HEEP
               if ((fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5)
                   && c_gsf_gsfet[j] > end_et
                   && c_gsf_isecaldriven[j]
                   && (c_gsf_hovere[j] * c_gsf_gsfet[j]) < (end_hoE * c_gsf_gsfet[j] + (end_coshFactor * cosh(end_etaFactor * c_gsfsc_eta[j]) - end_coshOffset) * c_rho)
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

            // fill histograms for all, LS and OS
            for (unsigned int k = 0; k< 3; ++k) {
               if (k == 1 && c_gsf_charge[GSF_passHEEP[0]] * c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) continue;
               if (k == 2 && c_gsf_charge[GSF_passHEEP[0]] * c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) continue;

               h_emuMass.at(p).at(k)->Fill(invMass, weight);
               if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 1.442) h_emuMassEB.at(p).at(k)->Fill(invMass, weight);
               if (fabs(c_gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(c_gsfsc_eta[GSF_passHEEP[0]]) < 2.5) h_emuMassEE.at(p).at(k)->Fill(invMass, weight);

               if (!passTrg) continue;

               // fill test histograms
               h_pfMet.at(p).at(k)->Fill(c_pfmet, weight);
               h_nVtx.at(p).at(k)->Fill(c_pvsize, weight);
               h_dDz.at(p).at(k)->Fill(fabs(c_gsf_dz[GSF_passHEEP[0]] - c_muon_dz_cmsCenter[MU_passGOOD[MU_leadingPassGOOD]]), weight);
               h_dDzBeamSpot.at(p).at(k)->Fill(fabs(c_gsf_dz_beamSpot[GSF_passHEEP[0]] - c_muon_dz_beamSpot[MU_passGOOD[MU_leadingPassGOOD]]), weight); // change this
               h_dDzFirstPVtx.at(p).at(k)->Fill(fabs(c_gsf_dz_firstPVtx[GSF_passHEEP[0]] - c_muon_dz_firstPVtx[MU_passGOOD[MU_leadingPassGOOD]]), weight); // change this
               h_rho.at(p).at(k)->Fill(c_rho, weight);
               h_nJets.at(p).at(k)->Fill(c_jetColl_size, weight);
               h_nJetsPt20.at(p).at(k)->Fill(jetsPt20, weight);
               h_nJetsPt30.at(p).at(k)->Fill(jetsPt30, weight);
               if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) < 3.14) h_dPhi.at(p).at(k)->Fill(fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), weight);
               if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) > 3.14) h_dPhi.at(p).at(k)->Fill(6.28 - fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), weight);
               h_eleEt.at(p).at(k)->Fill(c_gsf_gsfet[GSF_passHEEP[0]], weight);
               h_eleEta.at(p).at(k)->Fill(c_gsf_eta[GSF_passHEEP[0]], weight);
               h_elePhi.at(p).at(k)->Fill(c_gsf_phi[GSF_passHEEP[0]], weight);
               h_eleDEta.at(p).at(k)->Fill(c_gsf_deltaeta[GSF_passHEEP[0]], weight);
               h_eleDPhi.at(p).at(k)->Fill(c_gsf_deltaphi[GSF_passHEEP[0]], weight);
               h_eleHOE.at(p).at(k)->Fill(c_gsf_hovere[GSF_passHEEP[0]], weight);
               h_eleSigmaIEIE.at(p).at(k)->Fill(c_gsf_sigmaIetaIeta[GSF_passHEEP[0]], weight);
               h_eleEcalIso.at(p).at(k)->Fill(c_gsf_ecaliso[GSF_passHEEP[0]], weight);
               h_eleHcalIso12.at(p).at(k)->Fill(c_gsf_hcaliso1[GSF_passHEEP[0]] + c_gsf_hcaliso2[GSF_passHEEP[0]], weight);
               h_eleTrkIso.at(p).at(k)->Fill(c_gsf_trackiso[GSF_passHEEP[0]], weight);
               h_eleLostHits.at(p).at(k)->Fill(c_gsf_nLostInnerHits[GSF_passHEEP[0]], weight);

               h_muIsoCombRel.at(p).at(k)->Fill(CombRelIso, weight);
               h_muEtEleOPtMu.at(p).at(k)->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], weight);
               if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) h_muPtPlusOPtMinus.at(p).at(k)->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], weight);
               if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) h_muPtPlusOPtMinus.at(p).at(k)->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]/c_gsf_gsfet[GSF_passHEEP[0]], weight);
               h_muPt.at(p).at(k)->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muEta.at(p).at(k)->Fill(c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muPhi.at(p).at(k)->Fill(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muHitLayers.at(p).at(k)->Fill(c_muon_nlayerswithhits[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muPxlHits.at(p).at(k)->Fill(c_muon_nhitspixel[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muMuHits.at(p).at(k)->Fill(c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muDxyBeamSpot.at(p).at(k)->Fill(c_muon_dxy_beamSpot[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muNSeg.at(p).at(k)->Fill(c_muon_nSegmentMatch[MU_passGOOD[MU_leadingPassGOOD]], weight);
               h_muTrkIso03.at(p).at(k)->Fill(c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]], weight);
            }

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
               h_emu_mass_accVSgood.at(p).at(k)->Fill(invMass, weight);
            }
         }

      } // END LOOP OVER ENTRIES

      cout << "Number of selected events: " << evCounter << endl;

      if (p == DATA) {
         dataEntries = dataTrig[0] + dataTrig[1] + dataTrig[2];
         cout << "HLT_Mu22_Photon22_CaloIdL: " << dataTrig[0] << " , HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << dataTrig[1] << " , HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << dataTrig[2] << " , sum: " << dataTrig[0]+dataTrig[1]+dataTrig[2] << endl;
         cout << "Runrange HLT_Mu22_Photon22_CaloIdL: " << runs_HLT_Mu22_Photon22_CaloIdL.first << " - " << runs_HLT_Mu22_Photon22_CaloIdL.second << endl;
         cout << "Runrange HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first << " - " << runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second << endl;
         cout << "Runrange HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first << " - " << runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second << endl;
      }
      else cout << "HLT_Mu22_Photon22_CaloIdL: " << trig[0] << " , HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << trig[1] << " , HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL: " << trig[2] << " , sum: " << trig[0]+trig[1]+trig[2] << endl;

      if (p > DATA) {
         for (unsigned int k = 0; k < 3; ++k) {
            //SCALE MC
            h_emuMass.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_emuMassEB.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_emuMassEE.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_emu_mass_accVSgood.at(p).at(k)->Scale(input[p].second * LumiFactor);

            h_pfMet.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_nVtx.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_dDz.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_dDzBeamSpot.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_dDzFirstPVtx.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_rho.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_nJets.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_nJetsPt20.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_nJetsPt30.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_dPhi.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleEt.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleEta.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_elePhi.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleDEta.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleDPhi.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleHOE.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleSigmaIEIE.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleEcalIso.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleHcalIso12.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleTrkIso.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_eleLostHits.at(p).at(k)->Scale(input[p].second * LumiFactor);

            h_muIsoCombRel.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muEtEleOPtMu.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muPtPlusOPtMinus.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muPt.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muEta.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muPhi.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muHitLayers.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muPxlHits.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muMuHits.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muDxyBeamSpot.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muNSeg.at(p).at(k)->Scale(input[p].second * LumiFactor);
            h_muTrkIso03.at(p).at(k)->Scale(input[p].second * LumiFactor);
         }
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

   //SUM MC
   for (int p = 1; p < nbFile; ++p) {
      for (unsigned int k = 0; k < 3; ++k) {
         for (int q = p + 1; q < nbFile; ++q) {
            h_emuMass.at(p).at(k)->Add(h_emuMass.at(q).at(k));
            h_emuMassEB.at(p).at(k)->Add(h_emuMassEB.at(q).at(k));
            h_emuMassEE.at(p).at(k)->Add(h_emuMassEE.at(q).at(k));
            h_emu_mass_accVSgood.at(p).at(k)->Add(h_emu_mass_accVSgood.at(q).at(k));

            h_pfMet.at(p).at(k)->Add(h_pfMet.at(q).at(k));
            h_nVtx.at(p).at(k)->Add(h_nVtx.at(q).at(k));
            h_dDz.at(p).at(k)->Add(h_dDz.at(q).at(k));
            h_dDzBeamSpot.at(p).at(k)->Add(h_dDzBeamSpot.at(q).at(k));
            h_dDzFirstPVtx.at(p).at(k)->Add(h_dDzFirstPVtx.at(q).at(k));
            h_rho.at(p).at(k)->Add(h_rho.at(q).at(k));
            h_nJets.at(p).at(k)->Add(h_nJets.at(q).at(k));
            h_nJetsPt20.at(p).at(k)->Add(h_nJetsPt20.at(q).at(k));
            h_nJetsPt30.at(p).at(k)->Add(h_nJetsPt30.at(q).at(k));
            h_dPhi.at(p).at(k)->Add(h_dPhi.at(q).at(k));
            h_eleEt.at(p).at(k)->Add(h_eleEt.at(q).at(k));
            h_eleEta.at(p).at(k)->Add(h_eleEta.at(q).at(k));
            h_elePhi.at(p).at(k)->Add(h_elePhi.at(q).at(k));
            h_eleDEta.at(p).at(k)->Add(h_eleDEta.at(q).at(k));
            h_eleDPhi.at(p).at(k)->Add(h_eleDPhi.at(q).at(k));
            h_eleHOE.at(p).at(k)->Add(h_eleHOE.at(q).at(k));
            h_eleSigmaIEIE.at(p).at(k)->Add(h_eleSigmaIEIE.at(q).at(k));
            h_eleEcalIso.at(p).at(k)->Add(h_eleEcalIso.at(q).at(k));
            h_eleHcalIso12.at(p).at(k)->Add(h_eleHcalIso12.at(q).at(k));
            h_eleTrkIso.at(p).at(k)->Add(h_eleTrkIso.at(q).at(k));
            h_eleLostHits.at(p).at(k)->Add(h_eleLostHits.at(q).at(k));
 
            h_muIsoCombRel.at(p).at(k)->Add(h_muIsoCombRel.at(q).at(k));
            h_muEtEleOPtMu.at(p).at(k)->Add(h_muEtEleOPtMu.at(q).at(k));
            h_muPtPlusOPtMinus.at(p).at(k)->Add(h_muPtPlusOPtMinus.at(q).at(k));
            h_muPt.at(p).at(k)->Add(h_muPt.at(q).at(k));
            h_muEta.at(p).at(k)->Add(h_muEta.at(q).at(k));
            h_muPhi.at(p).at(k)->Add(h_muPhi.at(q).at(k));
            h_muHitLayers.at(p).at(k)->Add(h_muHitLayers.at(q).at(k));
            h_muPxlHits.at(p).at(k)->Add(h_muPxlHits.at(q).at(k));
            h_muMuHits.at(p).at(k)->Add(h_muMuHits.at(q).at(k));
            h_muDxyBeamSpot.at(p).at(k)->Add(h_muDxyBeamSpot.at(q).at(k));
            h_muNSeg.at(p).at(k)->Add(h_muNSeg.at(q).at(k));
            h_muTrkIso03.at(p).at(k)->Add(h_muTrkIso03.at(q).at(k));
         }
      }
   }
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   cout << "nb e+mu+     = " << nb_plus_plus << endl;
   cout << "nb e+mu-     = " << nb_plus_minus << endl;
   cout << "nb e-mu+     = " << nb_minus_plus << endl;
   cout << "nb e-mu-     = " << nb_minus_minus << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;

   emu_dilepton->Add(h_emuMass.at(1).at(ALL));
   //emu_ewk->Add(h_emuMass.at(6).at(ALL));
   //emu_jet->Add(h_emuMass.at(8).at(ALL));

   //WRITING
   output->Cd("");

   emu_dilepton->Write();
   //emu_ewk->Write();
   //emu_jet->Write();

   //emuLoose_data_nValidPv->Write();
   //emuLoose_ttbar_nValidPv->Write();
   //emuLoose_ztautau_nValidPv->Write();

   //emuLoose_dataOverTtbar_nValidPv->Write();
   //emuLoose_dataOverZtautau_nValidPv->Write();

   for (int p = 0; p < nbFile; ++p) {
      for (unsigned int k = 0; k < 3; ++k) {
         h_emuMass.at(p).at(k)->Write();
         h_emuMassEB.at(p).at(k)->Write();
         h_emuMassEE.at(p).at(k)->Write();
         h_emu_mass_accVSgood.at(p).at(k)->Write();

         h_pfMet.at(p).at(k)->Write();
         h_nVtx.at(p).at(k)->Write();
         h_dDz.at(p).at(k)->Write();
         h_dDzBeamSpot.at(p).at(k)->Write();
         h_dDzFirstPVtx.at(p).at(k)->Write();
         h_rho.at(p).at(k)->Write();
         h_nJets.at(p).at(k)->Write();
         h_nJetsPt20.at(p).at(k)->Write();
         h_nJetsPt30.at(p).at(k)->Write();
         h_dPhi.at(p).at(k)->Write();
         h_eleEt.at(p).at(k)->Write();
         h_eleEta.at(p).at(k)->Write();
         h_elePhi.at(p).at(k)->Write();
         h_eleDEta.at(p).at(k)->Write();
         h_eleDPhi.at(p).at(k)->Write();
         h_eleHOE.at(p).at(k)->Write();
         h_eleSigmaIEIE.at(p).at(k)->Write();
         h_eleEcalIso.at(p).at(k)->Write();
         h_eleHcalIso12.at(p).at(k)->Write();
         h_eleTrkIso.at(p).at(k)->Write();
         h_eleLostHits.at(p).at(k)->Write();

         h_muIsoCombRel.at(p).at(k)->Write();
         h_muEtEleOPtMu.at(p).at(k)->Write();
         h_muPtPlusOPtMinus.at(p).at(k)->Write();
         h_muPt.at(p).at(k)->Write();
         h_muEta.at(p).at(k)->Write();
         h_muPhi.at(p).at(k)->Write();
         h_muHitLayers.at(p).at(k)->Write();
         h_muPxlHits.at(p).at(k)->Write();
         h_muMuHits.at(p).at(k)->Write();
         h_muDxyBeamSpot.at(p).at(k)->Write();
         h_muNSeg.at(p).at(k)->Write();
         h_muTrkIso03.at(p).at(k)->Write();
      }
   }

   output->Close();
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
