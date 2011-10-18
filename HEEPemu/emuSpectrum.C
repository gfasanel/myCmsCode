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
#include <algorithm>

using namespace std;

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "TFile.h"

#endif // __CINT__

void emuSpectrum()
{
   bool bool_accessPUFile = true;
   //bool bool_accessPUFile=false;

   float LumiFactor = 3190.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 2179.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 1932.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 702.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON

   // max number of primary vertices
   unsigned nPVtxMax = 20;

   //
   //float Elec_trigger = 0.97;
   float Elec_trigger = 0.94;
   //float Elec_trigger = 0.98;
   //float Elec_trigger = 1.0;
   float Elec_ScaleFactor = 1.008 * 0.978;
   float Muon_ScaleFactor = 0.985;
   float Lumi_ScaleFactor = 1.0;

   //float MCemuScaleFactor = Elec_trigger * Elec_ScaleFactor * Muon_ScaleFactor;
   float MCemuScaleFactor = Elec_trigger * Elec_ScaleFactor * Muon_ScaleFactor * Lumi_ScaleFactor; //2011

   // 2*(LS_data - LS_mc)/totMC
   float QCDScaleFactor = 1.10;

   TH1::SetDefaultSumw2(kTRUE);

   //
   //ACCESSING PILE-UP FILE
   vector<float> PVX_ScalingFactor_ttba;
   vector<float> PVX_ScalingFactor_ztau;
   if (bool_accessPUFile == false) {
      for (unsigned int i = 0 ; i < nPVtxMax ; i++) {
         PVX_ScalingFactor_ttba.push_back(1.);
         PVX_ScalingFactor_ztau.push_back(1.);
      }
   }
   if (bool_accessPUFile) {
      TFile *inputPV = new TFile("EMu_3190pb-1_veto35GeV_PVXinfo.root", "open");
      //TFile *inputPV = new TFile("EMu_2179pb-1_veto35GeV_PVXinfo.root", "open");
      //TFile *inputPV = new TFile("EMu_1932pb-1_veto35GeV_PVXinfo.root", "open");
      //TFile *inputPV = new TFile("EMu_702pb-1_veto35GeV_PVXinfo.root", "open");
      inputPV->cd();

      TH1F *copy_emuloose_dataoverttba_nvalidpv;
      TH1F *copy_emuloose_dataoverztau_nvalidpv;
      //TH1F * copy_emuloose_dataoverttba_nvalidpv = new TH1F("copy_emuloose_dataoverttba_nvalidpv","copy_emuloose_dataoverttba_nvalidpv",nPVtxMax,0.,nPVtxMax);
      //TH1F * copy_emuloose_dataoverztau_nvalidpv = new TH1F("copy_emuloose_dataoverztau_nvalidpv","copy_emuloose_dataoverztau_nvalidpv",nPVtxMax,0.,nPVtxMax);

      copy_emuloose_dataoverttba_nvalidpv = (TH1F*)inputPV->Get("emuloose_dataoverttba_nvalidpv");
      copy_emuloose_dataoverztau_nvalidpv = (TH1F*)inputPV->Get("emuloose_dataoverztau_nvalidpv");

      for (unsigned int i = 0 ; i < nPVtxMax ; i++) {
         PVX_ScalingFactor_ttba.push_back(copy_emuloose_dataoverttba_nvalidpv->GetBinContent(i + 1) * nPVtxMax);
         PVX_ScalingFactor_ztau.push_back(copy_emuloose_dataoverztau_nvalidpv->GetBinContent(i + 1) * nPVtxMax);
      }
   }

   //
   // INPUT FILES
   vector<TFile *> input;
   vector<float> weight;

   //DATA
   input.push_back(new TFile("/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-177515_7TeV_Collisions11_JSON_3190pb-1.root", "open"));
   //input.push_back(new TFile("/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6-Cert_160404-173692_7TeV_Collisions11_JSON_2179pb-1.root", "open"));
   //input.push_back(new TFile("/user/treis/data2011/MuEG-Run2011A-May10ReReco+05Aug2011+PromptReco-AOD-Cert_160404-173244_7TeV_1932pb-1.root", "open"));
   //input.push_back(new TFile("/user_mnt/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/total_MuEG-160404-163869-ReReco10May-GoldenJSON-27May-191pb_+_MuEG-165088-166861-PromptV4-GoldenJSON-17Jun-511pb__702pb.root","open"));

   //MC
   //input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2-AODSIM_TreeEMuSkim25_RUN2/res/total_tree.root", "open"));
   input.push_back(new TFile("/user/treis/mcsamples/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v2_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root", "open"));
   input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2-AODSIM_TreeEMuSkim35_RUN1/res/total_missing1And13.root", "open"));
   //input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_WW_TuneZ2_7TeV_pythia6_tauola-Summer11-PU_S4_START42_V11-v1-AODSIM_TreeEMuSkim35_RUN2/res/total_tree.root", "open"));
   input.push_back(new TFile("/user/treis/mcsamples/WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root", "open"));
   //input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_WZ_TuneZ2_7TeV_pythia6_tauola-Summer11-PU_S4_START42_V11-v1-AODSIM_TreeEMuSkim35_RUN1/res/total_tree.root", "open"));
   input.push_back(new TFile("/user/treis/mcsamples/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root", "open"));
   input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/ForTreeProdCMSSW_4_2_3/src/UserCode/HEEPSkims/test/T_plus_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola-Summer11-PU_S4_START42_V11-v1-AODSIM_TreeEMuSkim35_total_tree.root","open")); 
   //input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola-Summer11-PU_S4_START42_V11-v1-AODSIM_TreeEMuSkim35_RUN1/res/total_tree.root", "open"));
   input.push_back(new TFile("/user/treis/mcsamples/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root", "open"));
   input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/ForTreeProdCMSSW_4_2_3/src/UserCode/HEEPSkims/test/MC-2011-v3_DYToMuMu_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v1-AODSIM_TreeEMuSkim35_S3V11v1_RUN2/res/total_tree.root","open"));
   input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/ForTreeProdCMSSW_4_2_3/src/UserCode/HEEPSkims/test/MC-2011-v3_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2-AODSIM_TreeEMuSkim35_RUN4/res/total_tree.root","open"));
   //input.push_back(new TFile("/user/treis/mcsamples/ZZTo4e_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root","open"));

   int nbFile = input.size();

   //WEIGHTS FOR 1pb-1
   weight.push_back(1.);                                  //DATA        0 black

   //weight.push_back(0.000149593);                         //TTbar       1 red         (1089625 event - xsect NNLO 163pb) -- 21 JUIN 2011
   weight.push_back(0.000044031);                         //TTbar       1 red         (3701947 event - xsect NNLO 163pb) -- 17.10.2011
   weight.push_back(0.000936761);                         //Ztautau     2 green       (2032536 event * 14/16 - xsect 1666pb)
   //weight.push_back(0.000010175);                         //WW          3 dark blue   (4225916 event - xsect NNLO 43.pb (note muons) ) -- 28 JUIL 2011
   weight.push_back(0.000020727);                         //WW          3 dark blue   (210667 event - xsect NLO 4.65pb (AN-11-259) ) -- 17.10.2011
   //weight.push_back(0.000004220);                         //WZ          4 yellow      (4265243 event - xsect NNLO 18.pb (note muons) ) -- 28 JUIL 2011
   weight.push_back(0.000002901);                         //WZ          4 yellow      (204725 event - xsect NLO 0.594pb (AN-11-259) ) -- 17.10.2011
   weight.push_back(0.000009664);                         //tW          5 pink        (814390/809984 event - xsect NNLO 7.87pb (note tW) ) -- 13 OCT 2011
   //weight.push_back(0.00063259);                          //W+jet       6 dark green  (49501047 event - xsect NNLO 31314pb) -- 23 JUIN 2011
   weight.push_back(0.000384917);                         //W+jet       6 dark green  (81352581 event - xsect NNLO 31314pb) -- 17.10.2011
   weight.push_back(0.00078075);                          //Zmumu       7 light blue  (2133856 event  - xsect NNLO 1666pb ) -- 14 OCT 2011) 
   weight.push_back(0.00073630*(26./25.));                //Zee         8 cyan        (2262653 event - xsect NNLO 1666pb ) -- 13 OCT 2011 
   //weight.push_back(0.000000164);                         //ZZ          9 violett     (499929 event - xsect NLO 0.082pb (AN-11-259) ) -- 17.10.2011 

   //MUON SELECTION
   float muon_et = 35.;
   float muon_etaMax = 2.4;
   int muon_nHitsMinGlobal = 11;
   int muon_nHitsMinPixel = 1;
   int muon_nHitsMinMuon = 1;
   float muon_impactParamMax = 0.2;   // in cm
   int muon_nSegMatchMin = 2;
   float muon_relIsoCutMax = 0.1;

   //HEEP CUTS

   //BARREL
   float bar_et = 35.;
   float bar_hoE = 0.05;
   float bar_DEta = 0.005;
   float bar_DPhi = 0.09;
   float bar_e2x5e5x5 = 0.94;
   float bar_e1x5e5x5 = 0.83;
   float bar_isoEcalHcal1_1 = 2.;
   float bar_isoEcalHcal1_2 = 0.03;
   float bar_isoTrack = 7.5;
   int bar_missInnerHits = 0;

   //ENDCAP
   float end_et = 40.;
   float end_hoE = 0.05;
   float end_DEta = 0.007 ;
   float end_DPhi = 0.09;
   float end_e2x5e5x5 = 0.;
   float end_e1x5e5x5 = 0.;
   float end_sigmaietaieta = 0.03;
   float end_isoEcalHcal1_1 = 2.5;
   float end_isoEcalHcal1_2 = 0.03;
   float end_isoHcal2 = 0.5;
   float end_isoTrack = 15.;
   int end_missInnerHits = 0;

   int nb_plus_plus = 0;
   int nb_plus_minus = 0;
   int nb_minus_plus = 0;
   int nb_minus_minus = 0;

   float sumNPV[5];
   float sumNEvt[5];

   sumNPV[0] = 0;
   sumNPV[1] = 0;
   sumNPV[2] = 0;
   sumNPV[3] = 0;
   sumNPV[4] = 0;
   sumNEvt[0] = 0;
   sumNEvt[1] = 0;
   sumNEvt[2] = 0;
   sumNEvt[3] = 0;
   sumNEvt[4] = 0;

   //
   vector<TH1F *> emu_mass;
   vector<TH1F *> emu_mass_passHEEP;
   vector<TH1F *> emu_mass_accVSgood;

   vector<TH1F *> emu_plus_plus;
   vector<TH1F *> emu_plus_minus;
   vector<TH1F *> emu_minus_plus;
   vector<TH1F *> emu_minus_minus;
   vector<TH1F *> emu_LS;
   vector<TH1F *> emu_OS;

   vector<TH1F *> test1;
   vector<TH1F *> test2;
   vector<TH1F *> test3;
   vector<TH1F *> test4;
   vector<TH1F *> test5;
   vector<TH1F *> test6;
   vector<TH1F *> test7;
   vector<TH1F *> test8;
   vector<TH1F *> test9;
   vector<TH1F *> test10;
   vector<TH1F *> test11;
   vector<TH1F *> test12;

   vector<TH1F *> test21;
   vector<TH1F *> test22;

   //
   TH1F *emu_data = new TH1F("emu_data", "emu_data", 100, 0., 1000.);
   TH1F *emu_dilepton = new TH1F("emu_dilepton", "emu_dilepton", 100, 0., 1000.);
   TH1F *emu_ewk = new TH1F("emu_ewk", "emu_ewk", 100, 0., 1000.);
   TH1F *emu_jet = new TH1F("emu_jet", "emu_jet", 100, 0., 1000.);

   TH1F *emu_ttbar = new TH1F("emu_ttbar", "emu_ttbar", 100, 0., 1000.);
   TH1F *emu_ztautau = new TH1F("emu_ztautau", "emu_ztautau", 100, 0., 1000.);
   TH1F *emu_ww = new TH1F("emu_ww", "emu_ww", 100, 0., 1000.);
   TH1F *emu_wz = new TH1F("emu_wz", "emu_wz", 100, 0., 1000.);
   TH1F *emu_tw = new TH1F("emu_tw", "emu_tw", 100, 0., 1000.);
   TH1F *emu_wjet = new TH1F("emu_wjet", "emu_wjet", 100, 0., 1000.);
   TH1F *emu_zmumu = new TH1F("emu_zmumu", "emu_zmumu", 100, 0., 1000.);
   TH1F *emu_zee = new TH1F("emu_zee", "emu_zee", 100, 0., 1000.);
   //TH1F *emu_zz = new TH1F("emu_zz", "emu_zz", 100, 0., 1000.);
   TH1F *emu_qcd = new TH1F("emu_qcd", "emu_qcd", 100, 0., 1000.);

   //LS
   TH1F *OS_data = new TH1F("OS_data", "OS_data", 100, 0., 1000.);
   TH1F *OS_ttbar = new TH1F("OS_ttbar", "OS_ttbar", 100, 0., 1000.);
   TH1F *OS_ztautau = new TH1F("OS_ztautau", "OS_ztautau", 100, 0., 1000.);
   TH1F *OS_ww = new TH1F("OS_ww", "OS_ww", 100, 0., 1000.);
   TH1F *OS_wz = new TH1F("OS_wz", "OS_wz", 100, 0., 1000.);
   TH1F *OS_tw = new TH1F("OS_tw", "OS_tw", 100, 0., 1000.);
   TH1F *OS_wjet = new TH1F("OS_wjet", "OS_wjet", 100, 0., 1000.);
   TH1F *OS_zmum = new TH1F("OS_zmum", "OS_zmum", 100, 0., 1000.);
   TH1F *OS_zele = new TH1F("OS_zele", "OS_zele", 100, 0., 1000.);
   //TH1F *OS_zz = new TH1F("OS_zz", "OS_zz", 100, 0., 1000.);

   TH1F *LS_data = new TH1F("LS_data", "LS_data", 100, 0., 1000.);
   TH1F *LS_ttbar = new TH1F("LS_ttbar", "LS_ttbar", 100, 0., 1000.);
   TH1F *LS_ztautau = new TH1F("LS_ztautau", "LS_ztautau", 100, 0., 1000.);
   TH1F *LS_ww = new TH1F("LS_ww", "LS_ww", 100, 0., 1000.);
   TH1F *LS_wz = new TH1F("LS_wz", "LS_wz", 100, 0., 1000.);
   TH1F *LS_tw = new TH1F("LS_tw", "LS_tw", 100, 0., 1000.);
   TH1F *LS_wjet = new TH1F("LS_wjet", "LS_wjet", 100, 0., 1000.);
   TH1F *LS_zmum = new TH1F("LS_zmum", "LS_zmum", 100, 0., 1000.);
   TH1F *LS_zele = new TH1F("LS_zele", "LS_zele", 100, 0., 1000.);
   //TH1F *LS_zz = new TH1F("LS_zz", "LS_zz", 100, 0., 1000.);

   TH1F *emuloose_data_nvalidpv = new TH1F("emuloose_data_nvalidpv", "emuloose_data_nvalidpv", nPVtxMax, 0., nPVtxMax);
   TH1F *emuloose_ttba_nvalidpv = new TH1F("emuloose_ttba_nvalidpv", "emuloose_ttba_nvalidpv", nPVtxMax, 0., nPVtxMax);
   TH1F *emuloose_ztau_nvalidpv = new TH1F("emuloose_ztau_nvalidpv", "emuloose_ztau_nvalidpv", nPVtxMax, 0., nPVtxMax);

   TH1F *emuloose_dataoverttba_nvalidpv = new TH1F("emuloose_dataoverttba_nvalidpv", "emuloose_dataoverttba_nvalidpv", nPVtxMax, 0., nPVtxMax);
   TH1F *emuloose_dataoverztau_nvalidpv = new TH1F("emuloose_dataoverztau_nvalidpv", "emuloose_dataoverztau_nvalidpv", nPVtxMax, 0., nPVtxMax);

   TH1F *test = new TH1F("test", "test", nPVtxMax, 0., nPVtxMax);

   //
   //ELES

//   TH1F *test1_data = new TH1F("test1_data","test1_data",50,0.,400.);
//   TH1F *test1_ttba = new TH1F("test1_ttba","test1_ttba",50,0.,400.);
//   TH1F *test1_ztau = new TH1F("test1_ztau","test1_ztau",50,0.,400.);
//   TH1F *test1_wwtw = new TH1F("test1_wwtw","test1_wwtw",50,0.,400.);
//   TH1F *test1_wjet = new TH1F("test1_wjet","test1_wjet",50,0.,400.);
//   TH1F *test1_zmum = new TH1F("test1_zmum","test1_zmum",50,0.,400.); //MET
//   TH1F *test1_zele = new TH1F("test1_zele","test1_zele",50,0.,400.); //MET

   TH1F *test1_data = new TH1F("test1_data", "test1_data", nPVtxMax, 0., nPVtxMax);
   TH1F *test1_ttba = new TH1F("test1_ttba", "test1_ttba", nPVtxMax, 0., nPVtxMax);
   TH1F *test1_ztau = new TH1F("test1_ztau", "test1_ztau", nPVtxMax, 0., nPVtxMax);
   TH1F *test1_wwtw = new TH1F("test1_wwtw", "test1_wwtw", nPVtxMax, 0., nPVtxMax);
   TH1F *test1_wjet = new TH1F("test1_wjet", "test1_wjet", nPVtxMax, 0., nPVtxMax);
   TH1F *test1_zmum = new TH1F("test1_zmum", "test1_zmum", nPVtxMax, 0., nPVtxMax); //MET
   TH1F *test1_zele = new TH1F("test1_zele", "test1_zele", nPVtxMax, 0., nPVtxMax); //MET

   TH1F *test2_data = new TH1F("test2_data", "test2_data", nPVtxMax, 0., nPVtxMax);
   TH1F *test2_ttba = new TH1F("test2_ttba", "test2_ttba", nPVtxMax, 0., nPVtxMax);
   TH1F *test2_ztau = new TH1F("test2_ztau", "test2_ztau", nPVtxMax, 0., nPVtxMax);
   TH1F *test2_wwtw = new TH1F("test2_wwtw", "test2_wwtw", nPVtxMax, 0., nPVtxMax);
   TH1F *test2_wjet = new TH1F("test2_wjet", "test2_wjet", nPVtxMax, 0., nPVtxMax);
   TH1F *test2_zmum = new TH1F("test2_zmum", "test2_zmum", nPVtxMax, 0., nPVtxMax); //NVtx
   TH1F *test2_zele = new TH1F("test2_zele", "test2_zele", nPVtxMax, 0., nPVtxMax); //NVtx

   TH1F *test3_data = new TH1F("test3_data", "test3_data", 50, 0., 3.14);
   TH1F *test3_ttba = new TH1F("test3_ttba", "test3_ttba", 50, 0., 3.14);
   TH1F *test3_ztau = new TH1F("test3_ztau", "test3_ztau", 50, 0., 3.14);
   TH1F *test3_wwtw = new TH1F("test3_wwtw", "test3_wwtw", 50, 0., 3.14);
   TH1F *test3_wjet = new TH1F("test3_wjet", "test3_wjet", 50, 0., 3.14);
   TH1F *test3_zmum = new TH1F("test3_zmum", "test3_zmum", 50, 0., 3.14); // Dphi
   TH1F *test3_zele = new TH1F("test3_zele", "test3_zele", 50, 0., 3.14); // Dphi

   TH1F *test4_data = new TH1F("test4_data", "test4_data", 50, 0., 500.);
   TH1F *test4_ttba = new TH1F("test4_ttba", "test4_ttba", 50, 0., 500.);
   TH1F *test4_ztau = new TH1F("test4_ztau", "test4_ztau", 50, 0., 500.);
   TH1F *test4_wwtw = new TH1F("test4_wwtw", "test4_wwtw", 50, 0., 500.);
   TH1F *test4_wjet = new TH1F("test4_wjet", "test4_wjet", 50, 0., 500.);
   TH1F *test4_zmum = new TH1F("test4_zmum", "test4_zmum", 50, 0., 500.); //pt
   TH1F *test4_zele = new TH1F("test4_zele", "test4_zele", 50, 0., 500.); //pt

   TH1F *test5_data = new TH1F("test5_data", "test5_data", 25, -2.5, 2.5);
   TH1F *test5_ttba = new TH1F("test5_ttba", "test5_ttba", 25, -2.5, 2.5);
   TH1F *test5_ztau = new TH1F("test5_ztau", "test5_ztau", 25, -2.5, 2.5);
   TH1F *test5_wwtw = new TH1F("test5_wwtw", "test5_wwtw", 25, -2.5, 2.5);
   TH1F *test5_wjet = new TH1F("test5_wjet", "test5_wjet", 25, -2.5, 2.5);
   TH1F *test5_zmum = new TH1F("test5_zmum", "test5_zmum", 25, -2.5, 2.5); //eta
   TH1F *test5_zele = new TH1F("test5_zele", "test5_zele", 25, -2.5, 2.5); //eta

   TH1F *test6_data = new TH1F("test6_data", "test6_data", 25, -3.14, 3.14);
   TH1F *test6_ttba = new TH1F("test6_ttba", "test6_ttba", 25, -3.14, 3.14);
   TH1F *test6_ztau = new TH1F("test6_ztau", "test6_ztau", 25, -3.14, 3.14);
   TH1F *test6_wwtw = new TH1F("test6_wwtw", "test6_wwtw", 25, -3.14, 3.14);
   TH1F *test6_wjet = new TH1F("test6_wjet", "test6_wjet", 25, -3.14, 3.14);
   TH1F *test6_zmum = new TH1F("test6_zmum", "test6_zmum", 25, -3.14, 3.14); //phi
   TH1F *test6_zele = new TH1F("test6_zele", "test6_zele", 25, -3.14, 3.14); //phi

   TH1F *test7_data = new TH1F("test7_data", "test7_data", 50, -0.008, 0.008);
   TH1F *test7_ttba = new TH1F("test7_ttba", "test7_ttba", 50, -0.008, 0.008);
   TH1F *test7_ztau = new TH1F("test7_ztau", "test7_ztau", 50, -0.008, 0.008);
   TH1F *test7_wwtw = new TH1F("test7_wwtw", "test7_wwtw", 50, -0.008, 0.008);
   TH1F *test7_wjet = new TH1F("test7_wjet", "test7_wjet", 50, -0.008, 0.008);
   TH1F *test7_zmum = new TH1F("test7_zmum", "test7_zmum", 50, -0.008, 0.008); //id1 ele
   TH1F *test7_zele = new TH1F("test7_zele", "test7_zele", 50, -0.008, 0.008); //id1 ele

   TH1F *test8_data = new TH1F("test8_data", "test8_data", 50, -0.1, 0.1);
   TH1F *test8_ttba = new TH1F("test8_ttba", "test8_ttba", 50, -0.1, 0.1);
   TH1F *test8_ztau = new TH1F("test8_ztau", "test8_ztau", 50, -0.1, 0.1);
   TH1F *test8_wwtw = new TH1F("test8_wwtw", "test8_wwtw", 50, -0.1, 0.1);
   TH1F *test8_wjet = new TH1F("test8_wjet", "test8_wjet", 50, -0.1, 0.1);
   TH1F *test8_zmum = new TH1F("test8_zmum", "test8_zmum", 50, -0.1, 0.1); //id2 ele
   TH1F *test8_zele = new TH1F("test8_zele", "test8_zele", 50, -0.1, 0.1); //id2 ele

   TH1F *test9_data = new TH1F("test9_data", "test9_data", 50, 0., 0.04);
   TH1F *test9_ttba = new TH1F("test9_ttba", "test9_ttba", 50, 0., 0.04);
   TH1F *test9_ztau = new TH1F("test9_ztau", "test9_ztau", 50, 0., 0.04);
   TH1F *test9_wwtw = new TH1F("test9_wwtw", "test9_wwtw", 50, 0., 0.04);
   TH1F *test9_wjet = new TH1F("test9_wjet", "test9_wjet", 50, 0., 0.04);
   TH1F *test9_zmum = new TH1F("test9_zmum", "test9_zmum", 50, 0., 0.04); //id3 ele
   TH1F *test9_zele = new TH1F("test9_zele", "test9_zele", 50, 0., 0.04); //id3 ele

   TH1F *test10_data = new TH1F("test10_data", "test10_data", 50, 0., 10.);
   TH1F *test10_ttba = new TH1F("test10_ttba", "test10_ttba", 50, 0., 10.);
   TH1F *test10_ztau = new TH1F("test10_ztau", "test10_ztau", 50, 0., 10.);
   TH1F *test10_wwtw = new TH1F("test10_wwtw", "test10_wwtw", 50, 0., 10.);
   TH1F *test10_wjet = new TH1F("test10_wjet", "test10_wjet", 50, 0., 10.);
   TH1F *test10_zmum = new TH1F("test10_zmum", "test10_zmum", 50, 0., 10.); //iso1 ele
   TH1F *test10_zele = new TH1F("test10_zele", "test10_zele", 50, 0., 10.); //iso1 ele

   TH1F *test11_data = new TH1F("test11_data", "test11_data", 50, 0., 10.);
   TH1F *test11_ttba = new TH1F("test11_ttba", "test11_ttba", 50, 0., 10.);
   TH1F *test11_ztau = new TH1F("test11_ztau", "test11_ztau", 50, 0., 10.);
   TH1F *test11_wwtw = new TH1F("test11_wwtw", "test11_wwtw", 50, 0., 10.);
   TH1F *test11_wjet = new TH1F("test11_wjet", "test11_wjet", 50, 0., 10.);
   TH1F *test11_zmum = new TH1F("test11_zmum", "test11_zmum", 50, 0., 10.); //iso2 ele
   TH1F *test11_zele = new TH1F("test11_zele", "test11_zele", 50, 0., 10.); //iso2 ele

   TH1F *test12_data = new TH1F("test12_data", "test12_data", 50, 0., 20.);
   TH1F *test12_ttba = new TH1F("test12_ttba", "test12_ttba", 50, 0., 20.);
   TH1F *test12_ztau = new TH1F("test12_ztau", "test12_ztau", 50, 0., 20.);
   TH1F *test12_wwtw = new TH1F("test12_wwtw", "test12_wwtw", 50, 0., 20.);
   TH1F *test12_wjet = new TH1F("test12_wjet", "test12_wjet", 50, 0., 20.);
   TH1F *test12_zmum = new TH1F("test12_zmum", "test12_zmum", 50, 0., 20.); //iso3 ele
   TH1F *test12_zele = new TH1F("test12_zele", "test12_zele", 50, 0., 20.); //iso3 ele

   //

//   //MUONS

//   TH1F *test1_data = new TH1F("test1_data", "test1_data", 50, 0., 0.2);
//   TH1F *test1_ttba = new TH1F("test1_ttba", "test1_ttba", 50, 0., 0.2);
//   TH1F *test1_ztau = new TH1F("test1_ztau", "test1_ztau", 50, 0., 0.2);
//   TH1F *test1_wwtw = new TH1F("test1_wwtw", "test1_wwtw", 50, 0., 0.2);
//   TH1F *test1_wjet = new TH1F("test1_wjet", "test1_wjet", 50, 0., 0.2);
//   TH1F *test1_zmum = new TH1F("test1_zmum", "test1_zmum", 50, 0., 0.2); // isol rel combinnée mu
//   TH1F *test1_zele = new TH1F("test1_zele", "test1_zele", 50, 0., 0.2); // isol rel combinnée mu

//   TH1F *test2_data = new TH1F("test2_data", "test2_data", 50, 0., 4.0);
//   TH1F *test2_ttba = new TH1F("test2_ttba", "test2_ttba", 50, 0., 4.0);
//   TH1F *test2_ztau = new TH1F("test2_ztau", "test2_ztau", 50, 0., 4.0);
//   TH1F *test2_wwtw = new TH1F("test2_wwtw", "test2_wwtw", 50, 0., 4.0);
//   TH1F *test2_wjet = new TH1F("test2_wjet", "test2_wjet", 50, 0., 4.0);
//   TH1F *test2_zmum = new TH1F("test2_zmum", "test2_zmum", 50, 0., 4.0); //
//   TH1F *test2_zele = new TH1F("test2_zele", "test2_zele", 50, 0., 4.0); // pt ele /pt mu

//   TH1F *test3_data = new TH1F("test3_data", "test3_data", 50, 0., 4.0);
//   TH1F *test3_ttba = new TH1F("test3_ttba", "test3_ttba", 50, 0., 4.0);
//   TH1F *test3_ztau = new TH1F("test3_ztau", "test3_ztau", 50, 0., 4.0);
//   TH1F *test3_wwtw = new TH1F("test3_wwtw", "test3_wwtw", 50, 0., 4.0);
//   TH1F *test3_wjet = new TH1F("test3_wjet", "test3_wjet", 50, 0., 4.0);
//   TH1F *test3_zmum = new TH1F("test3_zmum", "test3_zmum", 50, 0., 4.0); //
//   TH1F *test3_zele = new TH1F("test3_zele", "test3_zele", 50, 0., 4.0); // pt + /pt -

//   TH1F *test4_data = new TH1F("test4_data", "test4_data", 50, 0., 500.);
//   TH1F *test4_ttba = new TH1F("test4_ttba", "test4_ttba", 50, 0., 500.);
//   TH1F *test4_ztau = new TH1F("test4_ztau", "test4_ztau", 50, 0., 500.);
//   TH1F *test4_wwtw = new TH1F("test4_wwtw", "test4_wwtw", 50, 0., 500.);
//   TH1F *test4_wjet = new TH1F("test4_wjet", "test4_wjet", 50, 0., 500.);
//   TH1F *test4_zmum = new TH1F("test4_zmum", "test4_zmum", 50, 0., 500.); //pt
//   TH1F *test4_zele = new TH1F("test4_zele", "test4_zele", 50, 0., 500.); //pt

//   TH1F *test5_data = new TH1F("test5_data", "test5_data", 25, -2.5, 2.5);
//   TH1F *test5_ttba = new TH1F("test5_ttba", "test5_ttba", 25, -2.5, 2.5);
//   TH1F *test5_ztau = new TH1F("test5_ztau", "test5_ztau", 25, -2.5, 2.5);
//   TH1F *test5_wwtw = new TH1F("test5_wwtw", "test5_wwtw", 25, -2.5, 2.5);
//   TH1F *test5_wjet = new TH1F("test5_wjet", "test5_wjet", 25, -2.5, 2.5);
//   TH1F *test5_zmum = new TH1F("test5_zmum", "test5_zmum", 25, -2.5, 2.5); //eta
//   TH1F *test5_zele = new TH1F("test5_zele", "test5_zele", 25, -2.5, 2.5); //eta

//   TH1F *test6_data = new TH1F("test6_data", "test6_data", 25, -3.14, 3.14);
//   TH1F *test6_ttba = new TH1F("test6_ttba", "test6_ttba", 25, -3.14, 3.14);
//   TH1F *test6_ztau = new TH1F("test6_ztau", "test6_ztau", 25, -3.14, 3.14);
//   TH1F *test6_wwtw = new TH1F("test6_wwtw", "test6_wwtw", 25, -3.14, 3.14);
//   TH1F *test6_wjet = new TH1F("test6_wjet", "test6_wjet", 25, -3.14, 3.14);
//   TH1F *test6_zmum = new TH1F("test6_zmum", "test6_zmum", 25, -3.14, 3.14); //phi
//   TH1F *test6_zele = new TH1F("test6_zele", "test6_zele", 25, -3.14, 3.14); //phi

//   TH1F *test7_data = new TH1F("test7_data", "test7_data", 50, 0., 10.);
//   TH1F *test7_ttba = new TH1F("test7_ttba", "test7_ttba", 50, 0., 10.);
//   TH1F *test7_ztau = new TH1F("test7_ztau", "test7_ztau", 50, 0., 10.);
//   TH1F *test7_wwtw = new TH1F("test7_wwtw", "test7_wwtw", 50, 0., 10.);
//   TH1F *test7_wjet = new TH1F("test7_wjet", "test7_wjet", 50, 0., 10.);
//   TH1F *test7_zmum = new TH1F("test7_zmum", "test7_zmum", 50, 0., 10.); //id1 mu
//   TH1F *test7_zele = new TH1F("test7_zele", "test7_zele", 50, 0., 10.); //id1 mu

//   TH1F *test8_data = new TH1F("test8_data", "test8_data", 32, 0., 32);
//   TH1F *test8_ttba = new TH1F("test8_ttba", "test8_ttba", 32, 0., 32);
//   TH1F *test8_ztau = new TH1F("test8_ztau", "test8_ztau", 32, 0., 32);
//   TH1F *test8_wwtw = new TH1F("test8_wwtw", "test8_wwtw", 32, 0., 32);
//   TH1F *test8_wjet = new TH1F("test8_wjet", "test8_wjet", 32, 0., 32);
//   TH1F *test8_zmum = new TH1F("test8_zmum", "test8_zmum", 32, 0., 32); //id2 mu
//   TH1F *test8_zele = new TH1F("test8_zele", "test8_zele", 32, 0., 32); //id2 mu

//   TH1F *test9_data = new TH1F("test9_data", "test9_data", 52, 0., 52.);
//   TH1F *test9_ttba = new TH1F("test9_ttba", "test9_ttba", 52, 0., 52.);
//   TH1F *test9_ztau = new TH1F("test9_ztau", "test9_ztau", 52, 0., 52.);
//   TH1F *test9_wwtw = new TH1F("test9_wwtw", "test9_wwtw", 52, 0., 52.);
//   TH1F *test9_wjet = new TH1F("test9_wjet", "test9_wjet", 52, 0., 52.);
//   TH1F *test9_zmum = new TH1F("test9_zmum", "test9_zmum", 52, 0., 52.); //id3 mu
//   TH1F *test9_zele = new TH1F("test9_zele", "test9_zele", 52, 0., 52.); //id3 mu

//   TH1F *test10_data = new TH1F("test10_data", "test10_data", 50, 0., 10.);
//   TH1F *test10_ttba = new TH1F("test10_ttba", "test10_ttba", 50, 0., 10.);
//   TH1F *test10_ztau = new TH1F("test10_ztau", "test10_ztau", 50, 0., 10.);
//   TH1F *test10_wwtw = new TH1F("test10_wwtw", "test10_wwtw", 50, 0., 10.);
//   TH1F *test10_wjet = new TH1F("test10_wjet", "test10_wjet", 50, 0., 10.);
//   TH1F *test10_zmum = new TH1F("test10_zmum", "test10_zmum", 50, 0., 10.); //iso1 mu
//   TH1F *test10_zele = new TH1F("test10_zele", "test10_zele", 50, 0., 10.); //iso1 mu

//   TH1F *test11_data = new TH1F("test11_data", "test11_data", 50, 0., 10.);
//   TH1F *test11_ttba = new TH1F("test11_ttba", "test11_ttba", 50, 0., 10.);
//   TH1F *test11_ztau = new TH1F("test11_ztau", "test11_ztau", 50, 0., 10.);
//   TH1F *test11_wwtw = new TH1F("test11_wwtw", "test11_wwtw", 50, 0., 10.);
//   TH1F *test11_wjet = new TH1F("test11_wjet", "test11_wjet", 50, 0., 10.);
//   TH1F *test11_zmum = new TH1F("test11_zmum", "test11_zmum", 50, 0., 10.); //iso2 mu
//   TH1F *test11_zele = new TH1F("test11_zele", "test11_zele", 50, 0., 10.); //iso2 mu

//   TH1F *test12_data = new TH1F("test12_data", "test12_data", 50, 0., 20.);
//   TH1F *test12_ttba = new TH1F("test12_ttba", "test12_ttba", 50, 0., 20.);
//   TH1F *test12_ztau = new TH1F("test12_ztau", "test12_ztau", 50, 0., 20.);
//   TH1F *test12_wwtw = new TH1F("test12_wwtw", "test12_wwtw", 50, 0., 20.);
//   TH1F *test12_wjet = new TH1F("test12_wjet", "test12_wjet", 50, 0., 20.);
//   TH1F *test12_zmum = new TH1F("test12_zmum", "test12_zmum", 50, 0., 20.); //iso3 mu
//   TH1F *test12_zele = new TH1F("test12_zele", "test12_zele", 50, 0., 20.); //iso3 mu

   //
   //DECLARE HISTOS FOR THE SAMPLES
   for (int i = 0; i < nbFile; ++i) {

      //EE
      //emu_mass.push_back(new TH1F("emu_mass", "emu_mass", 200, 0., 1000.));
      emu_mass.push_back(new TH1F("emu_mass", "emu_mass", 100, 0., 1000.));
      emu_mass_passHEEP.push_back(new TH1F("emu_mass_passHEEP", "emu_mass_passHEEP", 100, 0., 1000.));
      emu_mass_accVSgood.push_back(new TH1F("emu_mass_accVSgood", "emu_mass_accVSgood", 100, 0., 1000.));

//     emu_plus_plus.push_back(new TH1F("emu_plus_plus", "emu_plus_plus", 100, 0., 1000.));
//     emu_plus_minus.push_back(new TH1F("emu_plus_minus", "emu_plus_minus", 100, 0., 1000.));
//     emu_minus_plus.push_back(new TH1F("emu_minus_plus", "emu_minus_plus", 100, 0., 1000.));
//     emu_minus_minus.push_back(new TH1F("emu_minus_minus", "emu_minus_minus", 100, 0., 1000.));
      emu_LS.push_back(new TH1F("emu_LS", "emu_LS", 100, 0., 1000.));
      emu_OS.push_back(new TH1F("emu_OS", "emu_OS", 100, 0., 1000.));

      //ele
      //test1.push_back(new TH1F("test1", "test1", 50, 0., 400.));
      test1.push_back(new TH1F("test1", "test1", nPVtxMax, 0., nPVtxMax));
      test2.push_back(new TH1F("test2", "test2", nPVtxMax, 0, nPVtxMax));
      test3.push_back(new TH1F("test3", "test3", 50, 0., 3.14));

      test4.push_back(new TH1F("test4", "test4", 50, 0., 500.));
      test5.push_back(new TH1F("test5", "test5", 25, -2.5, 2.5));
      test6.push_back(new TH1F("test6", "test6", 25, -3.14, 3.14));

      test7.push_back(new TH1F("test7", "test7", 50, -0.008, 0.008));
      test8.push_back(new TH1F("test8", "test8", 50, -0.1, 0.1));
      test9.push_back(new TH1F("test9", "test9", 50, 0., 0.04));
      test10.push_back(new TH1F("test10", "test10", 50, 0., 10.));
      test11.push_back(new TH1F("test11", "test11", 50, 0., 10.));
      test12.push_back(new TH1F("test12", "test12", 50, 0., 20.));

      // //muons
      // test1.push_back(new TH1F("test1", "test1", 50, 0., 0.2));
//     test2.push_back(new TH1F("test2", "test2", 50, 0., 4.));
//     test3.push_back(new TH1F("test3", "test3", 50, 0., 4.));

//     test4.push_back(new TH1F("test4", "test4", 50, 0., 500.));
//     test5.push_back(new TH1F("test5", "test5", 25, -2.5, 2.5 ));
//     test6.push_back(new TH1F("test6", "test6", 25, -3.14, 3.14 ));

//     test7.push_back(new TH1F("test7", "test7", 50, 0., 10.));
//     test8.push_back(new TH1F("test8", "test8", 32, 0., 32.));
//     test9.push_back(new TH1F("test9", "test9", 52, 0., 52.));
//     test10.push_back(new TH1F("test10", "test10", 50, 0., 10.));
//     test11.push_back(new TH1F("test11", "test11", 50, 0., 10.));
//     test12.push_back(new TH1F("test12", "test12", 50, 0., 20.));

   }

   TH1F* emu_total = new TH1F("emu_total", "emu_total", 100, 0., 1000.);

   //
   // GETTING FILES
   for (int p = 0; p < nbFile; ++p) {

      cout << "accessing file nr " << p + 1 << endl;

      // Get the TREE
      (input[p])->cd();

      TTree *thetree;

      //if (p!=0) thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");
      //if (p==0) thetree = (TTree*)(input[p])->Get("tree");

      thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");

      //RUN ID

      int c_runnumber;
      int c_eventnumber;

      int c_HLT_Mu15_Photon20_CaloIdL;

      //GLOBAL
      int c_nJetsAKT_pt15;
      int c_nJetsIC5_pt15;
      float c_calomet;
      float c_met;
      float c_mass;
      float c_pthat;
      float c_bsposx;
      float c_bsposy;
      float c_bsposz;

      //JETS IC5
      int c_jetIC5_size;
      float c_jetIC5_pt[100];
      float c_jetIC5_eta[100];
      float c_jetIC5_phi[100];
      float c_jetIC5_em[100];

      //PRIM VTX
      int c_pvsize;
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
      float c_gsf_isecaldriven[20];
      float c_gsf_istrackerdriven[20];
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
      float c_gsf_SwissCross[20];

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

      //
      //RUN ID
      TBranch        *b_runnumber;
      TBranch        *b_eventnumber;

      TBranch        *b_HLT_Mu15_Photon20_CaloIdL;

      //GLOBAL
      TBranch        *b_nJetsAKT_pt15;
      TBranch        *b_nJetsIC5_pt15;
      TBranch        *b_calomet;
      TBranch        *b_met;
      TBranch        *b_mass;
      TBranch        *b_pthat;
      TBranch        *b_bsposx;
      TBranch        *b_bsposy;
      TBranch        *b_bsposz;

      //JETS IC5
      TBranch        *b_jetIC5_size;
      TBranch        *b_jetIC5_pt;
      TBranch        *b_jetIC5_eta;
      TBranch        *b_jetIC5_phi;
      TBranch        *b_jetIC5_em;

      //PRIM VTX
      TBranch        *b_pvsize;
      TBranch        *b_pv_isValid;
      TBranch        *b_pv_ndof;
      TBranch        *b_pv_nTracks;
      TBranch        *b_pv_normChi2;
      TBranch        *b_pv_totTrackSize;

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
      TBranch        *b_gsf_SwissCross;

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

      //
//    //RUN ID
      thetree->SetBranchAddress("runnumber", &c_runnumber, &b_runnumber);
      thetree->SetBranchAddress("eventnumber", &c_eventnumber, &b_eventnumber);

      //HLT TRIGGER BITS
      thetree->SetBranchAddress("HLT_Mu15_Photon20_CaloIdL", &c_HLT_Mu15_Photon20_CaloIdL, &b_HLT_Mu15_Photon20_CaloIdL);

//    //GLOBAL
//     thetree->SetBranchAddress("nJetsAKT_pt15",&c_nJetsAKT_pt15,&b_nJetsAKT_pt15);
//     thetree->SetBranchAddress("nJetsIC5_pt15",&c_nJetsIC5_pt15,&b_nJetsIC5_pt15);
      thetree->SetBranchAddress("calomet", &c_calomet, &b_calomet);
      thetree->SetBranchAddress("met", &c_met, &b_met);
//     thetree->SetBranchAddress("mass",&c_mass,&b_mass);
      thetree->SetBranchAddress("pthat", &c_pthat, &b_pthat);
      thetree->SetBranchAddress("bsposx", &c_bsposx, &b_bsposx);
      thetree->SetBranchAddress("bsposy", &c_bsposy, &b_bsposy);
      thetree->SetBranchAddress("bsposz", &c_bsposz, &b_bsposz);

      //  //JETS IC5
      // thetree->SetBranchAddress("jetIC5_size",&c_jetIC5_size,&b_jetIC5_size);
//     thetree->SetBranchAddress("jetIC5_pt",&c_jetIC5_pt,&b_jetIC5_pt);
//     thetree->SetBranchAddress("jetIC5_eta",&c_jetIC5_eta,&b_jetIC5_eta);
//     thetree->SetBranchAddress("jetIC5_phi",&c_jetIC5_phi,&b_jetIC5_phi);
//     thetree->SetBranchAddress("jetIC5_em",&c_jetIC5_em,&b_jetIC5_em);

      //PRIM VTX
      thetree->SetBranchAddress("pvsize", &c_pvsize, &b_pvsize);
      thetree->SetBranchAddress("pv_isValid", &c_pv_isValid, &b_pv_isValid);
      thetree->SetBranchAddress("pv_ndof", &c_pv_ndof, &b_pv_ndof);
      thetree->SetBranchAddress("pv_nTracks", &c_pv_nTracks, &b_pv_nTracks);
      thetree->SetBranchAddress("pv_normChi2", &c_pv_normChi2, &b_pv_normChi2);
      thetree->SetBranchAddress("pv_totTrackSize", &c_pv_totTrackSize, &b_pv_totTrackSize);

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
      thetree->SetBranchAddress("gsf_SwissCross", &c_gsf_SwissCross, &b_gsf_SwissCross);
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
      //

      Long64_t nentries = (*thetree).GetEntries();

      //LOOP OVER EVENTS
      for (unsigned int i = 0; i < nentries; ++i) {
         thetree->GetEntry(i);

         if (p == 0 && c_HLT_Mu15_Photon20_CaloIdL == 0) continue; //ASK MuPhoton trigger bit

         vector<int> GSF_passHEEP;
         vector<int> GSF_passACC;

         vector<int> MU_passGOOD;
         vector<int> MU_passACC;

         //QCD removing
//       if (p==8 && c_pthat > 30) continue; //pt15
//       if (p==9 && c_pthat > 80) continue; //pt30
//       if (p==10 && c_pthat > 170) continue; //pt80

         //PRIMARY VTX COUNTING
         //int n_pvValid = 1;
         unsigned int n_pvValid = 0;
         for (int j = 0; j < c_pvsize; ++j) {
            if (c_pv_ndof[j] > 3 && c_pv_nTracks[j] > 3.)
               n_pvValid++;
         }
         if (p == 0 && n_pvValid < 1) continue;
         if ((p == 1 || p == 2) && n_pvValid < 1) continue; //ttbar

         //FILL THE VTX WEIGTH
         float npv_weight = 1.;

         for (unsigned int v = 0; v < nPVtxMax; ++v) {//LOOP OVER NPV
            if (p == 1 && n_pvValid == v) npv_weight = PVX_ScalingFactor_ttba[v]; //ttbar
            //if (p==2 && n_pvValid == v) npv_weight = PVX_ScalingFactor_ztau[v]; //ztau
         }

         //FOR VTX MC PONDERATION

         //CREATE VTX PONDERATION PLOT
         float gsfPtMaxB = 0.;
         float gsfPtMaxE = 0.;
         int gsfECAL = 0;
         float muPtMax = 0.;
         //LOOP OVER ELES
         for (int j = 0; j < c_gsf_size; ++j) {
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
            if (c_muon_pt[j] > muPtMax) muPtMax = c_muon_pt[j];

         if (((gsfECAL == 1 && gsfPtMaxB > 35.)  //BARREL
              ||
              (gsfECAL == -1 && gsfPtMaxE > 40.)) //ENDCAP
             && muPtMax > 35.) {
            if (p == 0) emuloose_data_nvalidpv->Fill(n_pvValid);
            else if (p == 1) emuloose_ttba_nvalidpv->Fill(n_pvValid);
            else if (p == 2) emuloose_ztau_nvalidpv->Fill(n_pvValid);

            //TEST
            //if (p==1 || p==0) (test1[p])->Fill(n_pvValid,npv_weight);
            if (p == 1) test->Fill(n_pvValid, npv_weight);
         }

         //LOOP OVER ELES
         for (int j = 0; j < c_gsf_size; ++j) {

            //CLEANING : FAKE ELES FROM MUONS
            bool fakeEle = false;
            for (int k = 0; k < c_muon_size; ++k) {
               if (c_muon_pt[k] < muon_et) continue;
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
                && fabs(c_gsf_isecaldriven[j] - 1.) < 0.01 //!!!
                && c_gsf_hovere[j] < bar_hoE
                && fabs(c_gsf_deltaeta[j]) < bar_DEta
                && fabs(c_gsf_deltaphi[j]) < bar_DPhi
                && (c_gsf_e2x5overe5x5[j] > bar_e2x5e5x5 || c_gsf_e1x5overe5x5[j] > bar_e1x5e5x5)
                && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * c_gsf_gsfet[j])
                && c_gsf_trackiso[j] < bar_isoTrack
               ) GSF_passHEEP.push_back(j);

            //ENDCAP HEEP
            if ((fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5)
                && c_gsf_gsfet[j] > end_et
                && fabs(c_gsf_isecaldriven[j] - 1.) < 0.01 //!!!
                && c_gsf_hovere[j] < end_hoE
                && fabs(c_gsf_deltaeta[j]) < end_DEta
                && fabs(c_gsf_deltaphi[j]) < end_DPhi
                && (c_gsf_e2x5overe5x5[j] > end_e2x5e5x5 || c_gsf_e1x5overe5x5[j] > end_e1x5e5x5)
                && c_gsf_sigmaIetaIeta[j] < end_sigmaietaieta
                && ((c_gsf_gsfet[j] < 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < end_isoEcalHcal1_1)
                    ||
                    (c_gsf_gsfet[j] >= 50. && (c_gsf_ecaliso[j] + c_gsf_hcaliso1[j]) < (end_isoEcalHcal1_1 + end_isoEcalHcal1_2 * (c_gsf_gsfet[j] - 50.))))
                && c_gsf_hcaliso2[j] < end_isoHcal2
                && c_gsf_trackiso[j] < end_isoTrack
               ) GSF_passHEEP.push_back(j);

            //GSF IN ACCEPTANCE
            if (c_gsf_gsfet[j] > bar_et
                && (fabs(c_gsfsc_eta[j]) < 1.442 || (fabs(c_gsfsc_eta[j]) > 1.56 && fabs(c_gsfsc_eta[j]) < 2.5))
                && (fabs(c_gsfsc_eta[j]) > 1.56 || c_gsf_SwissCross[j] < 0.95)
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
                && c_muon_dxy_beamSpot[j] < muon_impactParamMax
                && c_muon_nSegmentMatch[j] >= muon_nSegMatchMin
                && c_muon_isTrackerMuon[j]
                && c_muon_trackIso03[j]/ c_muon_pt[j] < muon_relIsoCutMax
               ) MU_passGOOD.push_back(j);

            //MU PASS ACC
            if (c_muon_pt[j] > muon_et
                && fabs(c_muon_eta[j]) < 2.4
               ) MU_passACC.push_back(j);
         }

         //search highest pt passing muon if there are more than one
         unsigned int MU_leadingPassGOOD = 0;
         if (MU_passGOOD.size() > 1) {
            for (unsigned int j = 0; j < MU_passGOOD.size(); ++j)
               if (c_muon_pt[MU_passGOOD[j]] > c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]) MU_leadingPassGOOD = j;
         }

         //HEEP - MU_GOOD
         if (GSF_passHEEP.size() > 0 && MU_passGOOD.size() > 0
             //&& c_calomet > 50.
             //&& c_calomet < 50.
             //&& c_met > 40.
            ) {

            //REMOVE DUPLICATES
            if (fabs(c_gsf_phi[GSF_passHEEP[0]] - c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]]) < 0.1) continue;

            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(c_gsf_gsfet[GSF_passHEEP[0]], c_gsf_eta[GSF_passHEEP[0]], c_gsf_phi[GSF_passHEEP[0]], 0.000511);
            mu1.SetPtEtaPhiM(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], 0.10566);

            double invMass = (ele1 + mu1).M();

            //MASS CUT
            if (invMass < 60.) continue;
            //if (invMass < 120.) continue;

            //ASKING PV
            //
            //if (p==0 && n_pvValid < 2) continue; // ASKING >1 PV

            //LS et SS CUT !!
            //if (c_gsf_charge[GSF_passHEEP[0]]*c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) continue; // ASK SS
            //if (c_gsf_charge[GSF_passHEEP[0]]*c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) continue; // ASK OS

            //
            (emu_mass[p])->Fill(invMass, npv_weight);

            //
            float CombRelIso = (c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]]) / c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];
            //

            //ELES

            (test1[p])->Fill(c_calomet, npv_weight);
            //(test2[p])->Fill(n_pvValid,npv_weight);
            if (p == 1)(test2[p])->Fill(n_pvValid, npv_weight);
            ////(test2[p])->Fill(c_pvsize,npv_weight);
            if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) < 3.14)(test3[p])->Fill(fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), npv_weight);
            if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) > 3.14)(test3[p])->Fill(6.28 - fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), npv_weight);

            (test4[p])->Fill(c_gsf_gsfet[GSF_passHEEP[0]], npv_weight);
            (test5[p])->Fill(c_gsf_eta[GSF_passHEEP[0]], npv_weight);
            (test6[p])->Fill(c_gsf_phi[GSF_passHEEP[0]], npv_weight);

            (test7[p])->Fill(c_gsf_deltaeta[GSF_passHEEP[0]], npv_weight);
            (test8[p])->Fill(c_gsf_deltaphi[GSF_passHEEP[0]], npv_weight);
            (test9[p])->Fill(c_gsf_sigmaIetaIeta[GSF_passHEEP[0]], npv_weight);

            (test10[p])->Fill(c_gsf_ecaliso[GSF_passHEEP[0]], npv_weight);
            (test11[p])->Fill(c_gsf_hcaliso1[GSF_passHEEP[0]] + c_gsf_hcaliso2[GSF_passHEEP[0]], npv_weight);
            (test12[p])->Fill(c_gsf_trackiso[GSF_passHEEP[0]], npv_weight);

            // //MUONS

            // (test1[p])->Fill(CombRelIso,npv_weight);
//    (test2[p])->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) (test3[p])->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) (test3[p])->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]/c_gsf_gsfet[GSF_passHEEP[0]],npv_weight);

//    (test4[p])->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    (test5[p])->Fill(c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    (test6[p])->Fill(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);

//    (test7[p])->Fill(c_muon_normChi2[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    (test8[p])->Fill(c_muon_nhitstrack[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    (test9[p])->Fill(c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);

//    (test10[p])->Fill(c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    (test11[p])->Fill(c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    (test12[p])->Fill(c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
//    //

            //LIKE SIGN

            if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) {
               //(emu_plus_plus[p])->Fill(invMass);
               (emu_LS[p])->Fill(invMass, npv_weight);
               if (p == 0) nb_plus_plus++;
            }

            if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) {
               //(emu_plus_minus[p])->Fill(invMass);
               (emu_OS[p])->Fill(invMass, npv_weight);
               if (p == 0) nb_plus_minus++;
            }
            if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) {
               //(emu_minus_plus[p])->Fill(invMass);
               (emu_OS[p])->Fill(invMass, npv_weight);
               if (p == 0) nb_minus_plus++;
            }
            if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) {
               //(emu_minus_minus[p])->Fill(invMass);
               (emu_LS[p])->Fill(invMass, npv_weight);
               if (p == 0) nb_minus_minus++;
            }

            if (p == 0 && invMass > 500.) cout << "HEEP-TightMU event in DATA (M>500) : run number = " << c_runnumber << " , event number = " << c_eventnumber << " , mass = " << invMass <<
               " , muon_pt = " << c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_eta = " << c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_phi = " << c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_trackIso03 = " << c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_emIso03 = " << c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_hadIso03 = " << c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_VETOtrackIso03 = " << c_muon_trackIso03_ptInVeto[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_VETOemIso03 = " << c_muon_emIso03_ptInVeto[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_VETOhadIso03 = " << c_muon_hadIso03_ptInVeto[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_normChi2 = " << c_muon_normChi2[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_dxy_beamSpot = " << c_muon_dxy_beamSpot[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_nhitstrack = " << c_muon_nhitstrack[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    " , muon_nhitsmuons = " << c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]] <<
           //    "                                                          " <<
               " , gsf_et = " << c_gsf_gsfet[GSF_passHEEP[0]] <<
           //    " , gsf_eta = " << c_gsf_eta[GSF_passHEEP[0]] <<
           //    " , gsf_phi = " << c_gsf_phi[GSF_passHEEP[0]] <<
           //    " , gsf_trackIso = " << c_gsf_trackiso[GSF_passHEEP[0]] <<
           //    " , gsf_ecalIso = " << c_gsf_ecaliso[GSF_passHEEP[0]] <<
           //    " , gsf_hcalIso1 = " << c_gsf_hcaliso1[GSF_passHEEP[0]] <<
           //    " , gsf_hcalIso2 = " << c_gsf_hcaliso2[GSF_passHEEP[0]] <<
               endl;

         }

         //GSF ACC - MU_GOOD
         if (GSF_passACC.size() > 0 && MU_passGOOD.size() > 0) {

            //REMOVE DUPLICATES
            if (fabs(c_gsf_phi[GSF_passACC[0]] - c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]]) < 0.1) continue;

            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(c_gsf_gsfet[GSF_passACC[0]], c_gsf_eta[GSF_passACC[0]], c_gsf_phi[GSF_passACC[0]], 0.000511);
            mu1.SetPtEtaPhiM(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], 0.10566);

            double invMass = (ele1 + mu1).M();
            (emu_mass_accVSgood[p])->Fill(invMass, npv_weight);
         }

      } // END LOOP OVER ENTRIES

   }//END FILE LOOP

   //PRIM VTX PRINTING

   // cout << "----PRIM VTX PART ------------------------- " << endl;
//   if (sumNEvt[0] != 0) cout << "nb event in 1st sample = " << sumNEvt[0] <<  " , average nb of PV =  " << sumNPV[0]/sumNEvt[0] << endl;
//   if (sumNEvt[1] != 0) cout << "nb event in 2nd sample = " << sumNEvt[1] <<  " , average nb of PV =  " << sumNPV[1]/sumNEvt[1] << endl;
//   if (sumNEvt[2] != 0) cout << "nb event in 3nd sample = " << sumNEvt[2] <<  " , average nb of PV =  " << sumNPV[2]/sumNEvt[2] << endl;
//   if (sumNEvt[3] != 0) cout << "nb event in 4nd sample = " << sumNEvt[3] <<  " , average nb of PV =  " << sumNPV[3]/sumNEvt[3] << endl;
//   if (sumNEvt[4] != 0) cout << "nb event in 5nd sample = " << sumNEvt[4] <<  " , average nb of PV =  " << sumNPV[4]/sumNEvt[4] << endl;
//   cout << "----END PRIM VTX PART --------------------- " << endl;
//   cout << "" << endl;

   //FOR PV SCALING FACTORS
   //CREATING PILE-UP FILE
   //emuloose_data_nvalidpv->Scale( 1./emuloose_data_nvalidpv->Integral() );
   //emuloose_ttba_nvalidpv->Scale( 1./emuloose_ttba_nvalidpv->Integral() );
   //emuloose_ztau_nvalidpv->Scale( 1./emuloose_ztau_nvalidpv->Integral() );

   emuloose_dataoverttba_nvalidpv->Divide(emuloose_data_nvalidpv, emuloose_ttba_nvalidpv);
   emuloose_dataoverttba_nvalidpv->Scale(1. / emuloose_dataoverttba_nvalidpv->Integral());

   TH1F * ttbarScaled = new TH1F("ttbarScaled", "ttbarScaled", nPVtxMax, 0., nPVtxMax);
   ttbarScaled->Multiply(emuloose_ttba_nvalidpv, emuloose_dataoverttba_nvalidpv);
   double ttbarEvents =  emuloose_ttba_nvalidpv->Integral();
   double scaledTtbarEvents =  ttbarScaled->Integral();

   cout << "++++++++++++ ttbarEvents " << ttbarEvents << endl;
   cout << "++++++++++++ scaledTtbarEvents " << scaledTtbarEvents << endl;
   cout << "++++++++++++ Correction factor " << ttbarEvents / scaledTtbarEvents / nPVtxMax << endl;

   emuloose_dataoverttba_nvalidpv->Scale(ttbarEvents / scaledTtbarEvents / nPVtxMax);
   // emuloose_dataoverztau_nvalidpv->Divide(emuloose_data_nvalidpv,emuloose_ztau_nvalidpv);
//   emuloose_dataoverztau_nvalidpv->Scale( 1./emuloose_dataoverztau_nvalidpv->Integral() );

   //
   //PLOTTING

   //POUR AVOIR DU STYLE
   gStyle->SetOptStat(011111110);
   gStyle->SetCanvasColor(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetHistFillColor(0);

   gStyle->SetTitleFillColor(0);
   gStyle->SetTitleFontSize(0.07);
   gStyle->SetStatColor(0);
   gStyle->SetStatFontSize(0.05);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.05);

   gStyle->SetOptTitle(0);

   //SCALE MC
   for (int p = 1; p < nbFile; p++) {
      (emu_mass[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (emu_mass_accVSgood[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);

//    (emu_plus_plus[p])->Scale(weight[p]*LumiFactor);
//     (emu_plus_minus[p])->Scale(weight[p]*LumiFactor);
//     (emu_minus_plus[p])->Scale(weight[p]*LumiFactor);
//     (emu_minus_minus[p])->Scale(weight[p]*LumiFactor);
      (emu_LS[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (emu_OS[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);

      (test1[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test2[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test3[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test4[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test5[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test6[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test7[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test8[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test9[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test10[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test11[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
      (test12[p])->Scale(weight[p]*LumiFactor * MCemuScaleFactor);
   }

   //PRINT INTEGRAL

   cout << "------------------------------------" << endl;
   cout << "HEEP - TIGHT MU        Lumi = " << LumiFactor << endl;
   cout << "                       pT   = " << bar_et << endl;
   cout << "                       M > 60 GeV" << endl;
   cout << "nb data      = " << (emu_mass[0])->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "nb TTbar     = " << (emu_mass[1])->Integral() << endl;
   cout << "nb Ztautau   = " << (emu_mass[2])->Integral() << endl;
   cout << "nb WW        = " << (emu_mass[3])->Integral() << endl;
   cout << "nb WZ        = " << (emu_mass[4])->Integral() << endl;
   cout << "nb tW        = " << (emu_mass[5])->Integral() << endl;
   cout << "- - - - - - - - - - - - - - - - - - " << endl;
   cout << "nb WJets     = " << (emu_mass[6])->Integral() << endl;
   cout << "nb Zmumu     = " << (emu_mass[7])->Integral() << endl;
   cout << "nb Zee       = " << (emu_mass[8])->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "TOT ttlike   = " << (emu_mass[1])->Integral() + (emu_mass[2])->Integral() + (emu_mass[3])->Integral() + (emu_mass[4])->Integral() + (emu_mass[5])->Integral() << endl;
   cout << "TOT contam   = " << (emu_mass[6])->Integral() + (emu_mass[7])->Integral() + (emu_mass[8])->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "TOT MC       = " << (emu_mass[1])->Integral() + (emu_mass[2])->Integral() + (emu_mass[3])->Integral() + (emu_mass[4])->Integral() + (emu_mass[5])->Integral() + (emu_mass[6])->Integral() + (emu_mass[7])->Integral() + (emu_mass[8])->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "                                    " << endl;


   //SUM MC
   for (int p = 1; p < nbFile; p++) {
      for (int q = p + 1; q < nbFile; q++) {
         (emu_mass[p])->Add(emu_mass[q]);
         (emu_mass_accVSgood[p])->Add(emu_mass_accVSgood[q]);

//       (emu_plus_plus[p])->Add(emu_plus_plus[q]);
//       (emu_plus_minus[p])->Add(emu_plus_minus[q]);
//       (emu_minus_plus[p])->Add(emu_minus_plus[q]);
//       (emu_minus_minus[p])->Add(emu_minus_minus[q]);
         (emu_LS[p])->Add(emu_LS[q]);
         (emu_OS[p])->Add(emu_OS[q]);

         (test1[p])->Add(test1[q]);
         (test2[p])->Add(test2[q]);
         (test3[p])->Add(test3[q]);
         (test4[p])->Add(test4[q]);
         (test5[p])->Add(test5[q]);
         (test6[p])->Add(test6[q]);
         (test7[p])->Add(test7[q]);
         (test8[p])->Add(test8[q]);
         (test9[p])->Add(test9[q]);
         (test10[p])->Add(test10[q]);
         (test11[p])->Add(test11[q]);
         (test12[p])->Add(test12[q]);
      }
   }

   cout << "                                    " << endl;
   cout << "                                    " << endl;
   cout << "------------------------------------" << endl;
   cout << "  Like Sign                         " << endl;
   cout << "                                    " << endl;
   cout << "nb LS DATA   = " << (emu_LS[0])->Integral() << " +- " << sqrt((emu_LS[0])->Integral()) << endl;
   cout << "nb LS MC     = " << (emu_LS[1])->Integral() << endl;
   cout << "                                    " << endl;
   cout << "nb OS DATA   = " << (emu_OS[0])->Integral() << " +- " << sqrt((emu_OS[0])->Integral()) << endl;
   cout << "nb OS MC     = " << (emu_OS[1])->Integral() << endl;
   cout << "- - - - - - - - - - - - - - - - - - " << endl;
   cout << "nb e+mu+     = " << nb_plus_plus << endl;
   cout << "nb e+mu-     = " << nb_plus_minus << endl;
   cout << "nb e-mu+     = " << nb_minus_plus << endl;
   cout << "nb e-mu-     = " << nb_minus_minus << endl;
   cout << "------------------------------------" << endl;

   cout << "--For ttlike-- (NO QCD CORR)" << endl;
   cout << "INTEGRAL M> 60" << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(7, 100)  - (emu_mass[6])->Integral(7, 100)   << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>120" << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(13, 100) - (emu_mass[6])->Integral(13, 100)  << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>200" << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(21, 100) - (emu_mass[6])->Integral(21, 100)  << endl;
   cout << "--------------" << endl;

   cout << "--For ttlike-- (NO QCD CORR) + W+Jet + Zmumu" << endl;
   cout << "INTEGRAL M> 60" << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(7, 100)  << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>120" << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(13, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>200" << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(21, 100) << endl;
   cout << "--------------" << endl;

   // SCALE BY QCD CORR FACTOR
   (emu_mass[1])->Scale(QCDScaleFactor);

   cout << "" << endl;
   cout << "" << endl;
   cout << "APRES CORRECTION PR QCD CONTRIBUTION :" << endl;
   cout << "" << endl;
   cout << "TOTAL INTEGRAL EMU" << endl;
   cout << "nb data     = " << (emu_mass[0])->Integral() << " +- " << sqrt((emu_mass[0])->Integral()) << endl;
   cout << "nb MC       = " << (emu_mass[1])->Integral() << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL EMU M>60 " << endl;
   cout << "nb data     = " << (emu_mass[0])->Integral(7, 100) << " +- " << sqrt((emu_mass[0])->Integral(7, 100)) << endl;
   cout << "nb MC       = " << (emu_mass[1])->Integral(7, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL EMU M>120" << endl;
   cout << "nb data     = " << (emu_mass[0])->Integral(13, 100) << " +- " << sqrt((emu_mass[0])->Integral(13, 100)) << endl;
   cout << "nb MC       = " << (emu_mass[1])->Integral(13, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL EMU M>200" << endl;
   cout << "nb data     = " << (emu_mass[0])->Integral(21, 100) << " +- " << sqrt((emu_mass[0])->Integral(21, 100)) << endl;
   cout << "nb MC       = " << (emu_mass[1])->Integral(21, 100) << endl;
   cout << "--------------" << endl;

   cout << "--For ttlike-- QCD CORR + W+Jet + Zmumu" << endl;
   cout << "60  < M < 120 " << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(7, 12)  << endl;
   cout << "--------------" << endl;
   cout << "120 < M < 200 " << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(13, 20) << endl;
   cout << "--------------" << endl;
   cout << "M > 200     " << endl;
   cout << "nb tt       = " << (emu_mass[1])->Integral(21, 100) << endl;
   cout << "--------------" << endl;


   // INTEGRATING FROM THE RIGHT SIDE!

//  //SCALE MC
//   for (int p = 0; p<nbFile; p++) {
//     //LOOPING ON BINS
//     for (int i = 1; i<51; i++) {
//       (emu_mass[p])->SetBinContent(i, (emu_mass[p])->Integral(i,50) );
//     }
//   }

   emu_total->Add(emu_mass[0]);

//   //E-E SPECTRUM
//   TCanvas *c0 = new TCanvas("c0", "c0", 100, 100, 1200, 725);
//
//   gStyle->SetOptStat(1111);
//
//   c0->cd(1); //
//
//   for (int p = 0; p < nbFile; p++) {
//
//      (emu_mass[p])->SetFillColor(p + 1);
//
//      if (p == 1) emu_mass[p]->Draw("HIST");
//      if (p > 1) emu_mass[p]->Draw("HISTsames");
//
//   }
//
//   (emu_mass[0])->SetFillColor(0);
//   (emu_mass[0])->SetLineWidth(2);
//   (emu_mass[0])->Draw("sames");
//   //(emu_mass[0])->Draw();
//
//   //----
//
//   //E-E SPECTRUM PASS HEEP
//   TCanvas *c1 = new TCanvas("c1", "c1", 100, 100, 1200, 725);
//
//   gStyle->SetOptStat(0000000000);
//   //gStyle->SetOptStat(1111);
//
//   c1->cd(1); //
//
//   for (int p = 0; p < nbFile; p++) {
//
//      (emu_mass_accVSgood[p])->SetFillColor(p + 1);
//
//      if (p == 1) {
//         (emu_mass_accVSgood[p])->GetXaxis()->SetTitle("M_{e #mu}");
//         (emu_mass_accVSgood[p])->Draw("HIST");
//      }
//      if (p > 1)(emu_mass_accVSgood[p])->Draw("HISTsames");
//
//   }
//
//   (emu_mass_accVSgood[0])->SetFillColor(0);
//   (emu_mass_accVSgood[0])->SetLineWidth(2);
//   (emu_mass_accVSgood[0])->Draw("sames");
//   //(emu_mass_passHEEP[0])->Draw();
//
//   TPaveLabel label0(0.300, 0.89, 0.740, 0.99, "7 TeV, #int L dt = 6.965 pb^{-1}", "brNDC");
//   label0.SetFillColor(0);
//   label0.SetFillStyle(0);
//   label0.SetBorderSize(0);
//   label0.SetTextSize(0.40);
//   label0.Draw();
//
//   //-------------------------------------

   //c0->Print("emu_mass.eps");


   emu_data->Add((emu_mass[0]));
   emu_dilepton->Add((emu_mass[1]));
   emu_ewk->Add((emu_mass[6]));
   //emu_jet->Add( (emu_mass[8]) );


   emu_ttbar->Add((emu_mass[1]));
   emu_ztautau->Add((emu_mass[2]));
   emu_ww->Add((emu_mass[3]));
   emu_wz->Add((emu_mass[4]));
   emu_tw->Add((emu_mass[5]));
   emu_wjet->Add((emu_mass[6]));
   emu_zmumu->Add((emu_mass[7]));
   emu_zee->Add((emu_mass[8]));
   //emu_zz->Add((emu_mass[9]));

   //LS
   OS_data->Add((emu_OS[0]));
   OS_ttbar->Add((emu_OS[1]));
   OS_ztautau->Add((emu_OS[2]));
   OS_ww->Add((emu_OS[3]));
   OS_wz->Add((emu_OS[4]));
   OS_tw->Add((emu_OS[5]));
   OS_wjet->Add((emu_OS[6]));
   OS_zmum->Add((emu_OS[7]));
   OS_zele->Add((emu_OS[8]));
   //OS_zz->Add((emu_OS[9]));

   //LS
   LS_data->Add((emu_LS[0]));
   LS_ttbar->Add((emu_LS[1]));
   LS_ztautau->Add((emu_LS[2]));
   LS_ww->Add((emu_LS[3]));
   LS_wz->Add((emu_LS[4]));
   LS_tw->Add((emu_LS[5]));
   LS_wjet->Add((emu_LS[6]));
   LS_zmum->Add((emu_LS[7]));
   LS_zele->Add((emu_LS[8]));
   //LS_zz->Add((emu_LS[9]));

   //tests

   test1_data->Add((test1[0]));
   test1_ttba->Add((test1[1]));
   test1_ztau->Add((test1[2]));
   test1_wwtw->Add((test1[3]));
   test1_wjet->Add((test1[6]));
   test1_zmum->Add((test1[7]));
   test1_zele->Add((test1[8]));

   test2_data->Add((test2[0]));
   test2_ttba->Add((test2[1]));
   test2_ztau->Add((test2[2]));
   test2_wwtw->Add((test2[3]));
   test2_wjet->Add((test2[6]));
   test2_zmum->Add((test2[7]));
   test2_zele->Add((test2[8]));

   test3_data->Add((test3[0]));
   test3_ttba->Add((test3[1]));
   test3_ztau->Add((test3[2]));
   test3_wwtw->Add((test3[3]));
   test3_wjet->Add((test3[6]));
   test3_zmum->Add((test3[7]));
   test3_zele->Add((test3[8]));

   test4_data->Add((test4[0]));
   test4_ttba->Add((test4[1]));
   test4_ztau->Add((test4[2]));
   test4_wwtw->Add((test4[3]));
   test4_wjet->Add((test4[6]));
   test4_zmum->Add((test4[7]));
   test4_zele->Add((test4[8]));

   test5_data->Add((test5[0]));
   test5_ttba->Add((test5[1]));
   test5_ztau->Add((test5[2]));
   test5_wwtw->Add((test5[3]));
   test5_wjet->Add((test5[6]));
   test5_zmum->Add((test5[7]));
   test5_zele->Add((test5[8]));

   test6_data->Add((test6[0]));
   test6_ttba->Add((test6[1]));
   test6_ztau->Add((test6[2]));
   test6_wwtw->Add((test6[3]));
   test6_wjet->Add((test6[6]));
   test6_zmum->Add((test6[7]));
   test6_zele->Add((test6[8]));

   test7_data->Add((test7[0]));
   test7_ttba->Add((test7[1]));
   test7_ztau->Add((test7[2]));
   test7_wwtw->Add((test7[3]));
   test7_wjet->Add((test7[6]));
   test7_zmum->Add((test7[7]));
   test7_zele->Add((test7[8]));

   test8_data->Add((test8[0]));
   test8_ttba->Add((test8[1]));
   test8_ztau->Add((test8[2]));
   test8_wwtw->Add((test8[3]));
   test8_wjet->Add((test8[6]));
   test8_zmum->Add((test8[7]));
   test8_zele->Add((test8[8]));

   test9_data->Add((test9[0]));
   test9_ttba->Add((test9[1]));
   test9_ztau->Add((test9[2]));
   test9_wwtw->Add((test9[3]));
   test9_wjet->Add((test9[6]));
   test9_zmum->Add((test9[7]));
   test9_zele->Add((test9[8]));

   test10_data->Add((test10[0]));
   test10_ttba->Add((test10[1]));
   test10_ztau->Add((test10[2]));
   test10_wwtw->Add((test10[3]));
   test10_wjet->Add((test10[6]));
   test10_zmum->Add((test10[7]));
   test10_zele->Add((test10[8]));

   test11_data->Add((test11[0]));
   test11_ttba->Add((test11[1]));
   test11_ztau->Add((test11[2]));
   test11_wwtw->Add((test11[3]));
   test11_wjet->Add((test11[6]));
   test11_zmum->Add((test11[7]));
   test11_zele->Add((test11[8]));

   test12_data->Add((test12[0]));
   test12_ttba->Add((test12[1]));
   test12_ztau->Add((test12[2]));
   test12_wwtw->Add((test12[3]));
   test12_wjet->Add((test12[6]));
   test12_zmum->Add((test12[7]));
   test12_zele->Add((test12[8]));


   //WRITING

   TFile output("testEmuSpec.root", "recreate");
   //TFile output("emu_mass_L6965nb_pT25_13October.root","recreate");

   output.cd();


   emu_data->Write();
   emu_dilepton->Write();
   emu_ewk->Write();
   emu_jet->Write();

   emu_ttbar->Write();
   emu_ztautau->Write();
   emu_ww->Write();
   emu_wz->Write();
   emu_tw->Write();
   emu_wjet->Write();
   emu_zmumu->Write();
   emu_zee->Write();
   //emu_zz->Write();
   emu_qcd->Write();

   emu_total->Write();

   OS_data->Write();
   OS_ttbar->Write();
   OS_ztautau->Write();
   OS_ww->Write();
   OS_wz->Write();
   OS_tw->Write();
   OS_wjet->Write();
   OS_zmum->Write();
   OS_zele->Write();
   //OS_zz->Write();

   LS_data->Write();
   LS_ttbar->Write();
   LS_ztautau->Write();
   LS_ww->Write();
   LS_wz->Write();
   LS_tw->Write();
   LS_wjet->Write();
   LS_zmum->Write();
   LS_zele->Write();
   //LS_zz->Write();

   emuloose_data_nvalidpv->Write();
   emuloose_ttba_nvalidpv->Write();
   emuloose_ztau_nvalidpv->Write();

   emuloose_dataoverttba_nvalidpv->Write();
   emuloose_dataoverztau_nvalidpv->Write();

   test->Write();

   test1_data->Write();
   test1_ttba->Write();
   test1_ztau->Write();
   test1_wwtw->Write();
   test1_wjet->Write();
   test1_zmum->Write();
   test1_zele->Write();

   test2_data->Write();
   test2_ttba->Write();
   test2_ztau->Write();
   test2_wwtw->Write();
   test2_wjet->Write();
   test2_zmum->Write();
   test2_zele->Write();

   test3_data->Write();
   test3_ttba->Write();
   test3_ztau->Write();
   test3_wwtw->Write();
   test3_wjet->Write();
   test3_zmum->Write();
   test3_zele->Write();

   test4_data->Write();
   test4_ttba->Write();
   test4_ztau->Write();
   test4_wwtw->Write();
   test4_wjet->Write();
   test4_zmum->Write();
   test4_zele->Write();

   test5_data->Write();
   test5_ttba->Write();
   test5_ztau->Write();
   test5_wwtw->Write();
   test5_wjet->Write();
   test5_zmum->Write();
   test5_zele->Write();

   test6_data->Write();
   test6_ttba->Write();
   test6_ztau->Write();
   test6_wwtw->Write();
   test6_wjet->Write();
   test6_zmum->Write();
   test6_zele->Write();

   test7_data->Write();
   test7_ttba->Write();
   test7_ztau->Write();
   test7_wwtw->Write();
   test7_wjet->Write();
   test7_zmum->Write();
   test7_zele->Write();

   test8_data->Write();
   test8_ttba->Write();
   test8_ztau->Write();
   test8_wwtw->Write();
   test8_wjet->Write();
   test8_zmum->Write();
   test8_zele->Write();

   test9_data->Write();
   test9_ttba->Write();
   test9_ztau->Write();
   test9_wwtw->Write();
   test9_wjet->Write();
   test9_zmum->Write();
   test9_zele->Write();

   test10_data->Write();
   test10_ttba->Write();
   test10_ztau->Write();
   test10_wwtw->Write();
   test10_wjet->Write();
   test10_zmum->Write();
   test10_zele->Write();

   test11_data->Write();
   test11_ttba->Write();
   test11_ztau->Write();
   test11_wwtw->Write();
   test11_wjet->Write();
   test11_zmum->Write();
   test11_zele->Write();

   test12_data->Write();
   test12_ttba->Write();
   test12_ztau->Write();
   test12_wwtw->Write();
   test12_wjet->Write();
   test12_zmum->Write();
   test12_zele->Write();

   output.Close();

//   //LIKE SIGN
//   TCanvas *c2 = new TCanvas("c2", "c2", 100, 100, 1200, 725);
//
//   gStyle->SetOptStat(1111);
//
//   c2->cd(1);
//
//   for (int p = 0; p < nbFile; p++) {
//
//      (emu_LS[p])->SetFillColor(p + 1);
//
//      if (p == 1)(emu_LS[p])->Draw("HIST");
//      if (p > 1)(emu_LS[p])->Draw("HISTsames");
//
//   }
//
//   (emu_LS[0])->SetFillColor(0);
//   (emu_LS[0])->SetLineWidth(2);
//   (emu_LS[0])->Draw("sames");
//
//
//   TCanvas *c3 = new TCanvas("c3", "c3", 100, 100, 1200, 725);
//
//   gStyle->SetOptStat(1111);
//
//   c3->cd(1);
//
//   for (int p = 0; p < nbFile; p++) {
//
//      (emu_OS[p])->SetFillColor(p + 1);
//
//      if (p == 1)(emu_OS[p])->Draw("HIST");
//      if (p > 1)(emu_OS[p])->Draw("HISTsames");
//
//   }
//
//   (emu_OS[0])->SetFillColor(0);
//   (emu_OS[0])->SetLineWidth(2);
//   (emu_OS[0])->Draw("sames");
//
}
