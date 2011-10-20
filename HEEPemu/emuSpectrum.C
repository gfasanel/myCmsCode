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

using namespace std;

void emuSpectrum()
{
   // parameters /////////////////////////////////////////////////////////////
   float LumiFactor = 3534.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 3190.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 2179.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 1932.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   //float LumiFactor = 702.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON

   // DATA file
   TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-178078_7TeV_Collisions11_JSON_3534pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-Cert_160404-177515_7TeV_Collisions11_JSON_3190pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco-v1+05Aug2011-v1+PromptReco-v4+PromptReco-v6-Cert_160404-173692_7TeV_Collisions11_JSON_2179pb-1.root";
   //TString dataFile = "/user/treis/data2011/MuEG-Run2011A-May10ReReco+05Aug2011+PromptReco-AOD-Cert_160404-173244_7TeV_1932pb-1.root";
   //TString dataFile = "/user_mnt/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/total_MuEG-160404-163869-ReReco10May-GoldenJSON-27May-191pb_+_MuEG-165088-166861-PromptV4-GoldenJSON-17Jun-511pb__702pb.root";

   unsigned nPVtxMax = 20; // max number of primary vertices

   // scale factors
   float Elec_trigger = 0.94;
   float Elec_ScaleFactor = 1.008 * 0.978;
   float Muon_ScaleFactor = 0.985;
   float Lumi_ScaleFactor = 1.0;

   bool calcQCDScaleFactor = false;
   float QCD_ScaleFactor = 1.10; // overwritten if calcQCDScaleFactor == true

   bool usePUInfo = false;
   bool generatePUFile = false;

   // selection cuts /////////////////////////////////////////////////////////
   float minInvMass = 60.;

   //MUON selection
   float muon_et = 35.;
   float muon_etaMax = 2.4;
   int muon_nHitsMinGlobal = 11;
   int muon_nHitsMinPixel = 1;
   int muon_nHitsMinMuon = 1;
   float muon_impactParamMax = 0.2;   // in cm
   int muon_nSegMatchMin = 2;
   float muon_relIsoCutMax = 0.1;

   //HEEP selection
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
   ///////////////////////////////////////////////////////////////////////////

   TH1::SetDefaultSumw2(kTRUE);

   float MCemuScaleFactor = Elec_trigger * Elec_ScaleFactor * Muon_ScaleFactor * Lumi_ScaleFactor; 

   ///////////////////////////////////////////////////////////////////////////
   // INPUT FILES
   vector<TFile *> input;
   vector<float> weight;

   // DATA
   input.push_back(new TFile(dataFile, "open"));

   // MC
   //input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2-AODSIM_TreeEMuSkim25_RUN2/res/total_tree.root", "open"));
   input.push_back(new TFile("/user/treis/mcsamples/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v2_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root", "open"));
   //input.push_back(new TFile("/user/vdero/ProdTreeSpring2011/CMSSW_4_2_1_patch2/src/UserCode/HEEPSkims/test/MC-2011-v2_DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2-AODSIM_TreeEMuSkim35_RUN1/res/total_missing1And13.root", "open"));
   input.push_back(new TFile("/user/treis/mcsamples/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_AODSIM_HEEPSkim1Ele1MuPt35_gct1_6.root", "open"));
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
   //weight.push_back(0.000936761);                         //Ztautau     2 green       (2032536 event * 14/16 - xsect 1666pb)
   weight.push_back(0.000821491);                         //Ztautau     2 green       (2028020 event  - xsect 1666pb)
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

   // counting variables
   int nb_plus_plus = 0;
   int nb_plus_minus = 0;
   int nb_minus_plus = 0;
   int nb_minus_minus = 0;

   // histogram containers
   vector<TH1F *> emuMass;
   vector<TH1F *> emuMass_LS;
   vector<TH1F *> emuMass_OS;
   vector<TH1F *> emu_mass_accVSgood;

   //vector<TH1F *> emu_plus_plus;
   //vector<TH1F *> emu_plus_minus;
   //vector<TH1F *> emu_minus_plus;
   //vector<TH1F *> emu_minus_minus;

   vector<TH1F *> eleMET;
   vector<TH1F *> eleNVtx;
   vector<TH1F *> eleDphi;
   vector<TH1F *> elePt;
   vector<TH1F *> eleEta;
   vector<TH1F *> elePhi;
   vector<TH1F *> eleId1;
   vector<TH1F *> eleId2;
   vector<TH1F *> eleId3;
   vector<TH1F *> eleIso1;
   vector<TH1F *> eleIso2;
   vector<TH1F *> eleIso3;

   vector<TH1F *> muIsoCombRel;
   vector<TH1F *> muPtEleOPtMu;
   vector<TH1F *> muPtPlusOPtMinus;
   vector<TH1F *> muPt;
   vector<TH1F *> muEta;
   vector<TH1F *> muPhi;
   vector<TH1F *> muId1;
   vector<TH1F *> muId2;
   vector<TH1F *> muId3;
   vector<TH1F *> muIso1;
   vector<TH1F *> muIso2;
   vector<TH1F *> muIso3;

   //
   TH1F *emu_dilepton = new TH1F("emu_dilepton", "emu_dilepton", 100, 0., 1000.);
   TH1F *emu_ewk = new TH1F("emu_ewk", "emu_ewk", 100, 0., 1000.);
   TH1F *emu_jet = new TH1F("emu_jet", "emu_jet", 100, 0., 1000.);

   TH1F *emuLoose_data_nValidPv = new TH1F("emuLoose_data_nValidPv", "emuLoose_data_nValidPv", nPVtxMax, 0., nPVtxMax);
   TH1F *emuLoose_ttbar_nValidPv = new TH1F("emuLoose_ttbar_nValidPv", "emuLoose_ttbar_nValidPv", nPVtxMax, 0., nPVtxMax);
   TH1F *emuLoose_ztautau_nValidPv = new TH1F("emuLoose_ztautau_nValidPv", "emuLoose_ztautau_nValidPv", nPVtxMax, 0., nPVtxMax);

   TH1F *emuLoose_dataOverTtbar_nValidPv = new TH1F("emuLoose_dataOverTtbar_nValidPv", "emuLoose_dataOverTtbar_nValidPv", nPVtxMax, 0., nPVtxMax);
   TH1F *emuLoose_dataOverZtautau_nValidPv = new TH1F("emuLoose_dataOverZtautau_nValidPv", "emuLoose_dataOverZtautau_nValidPv", nPVtxMax, 0., nPVtxMax);

   //RUN ID
   int c_runnumber;
   int c_eventnumber;

   //TRIGGER
   int c_HLT_Mu15_Photon20_CaloIdL;

   //GLOBAL
   //int c_nJetsAKT_pt15;
   //int c_nJetsIC5_pt15;
   float c_calomet;
   float c_met;
   //float c_mass;
   float c_pthat;
   float c_bsposx;
   float c_bsposy;
   float c_bsposz;

   //JETS IC5
   //int c_jetIC5_size;
   //float c_jetIC5_pt[100];
   //float c_jetIC5_eta[100];
   //float c_jetIC5_phi[100];
   //float c_jetIC5_em[100];

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

   //RUN ID
   TBranch        *b_runnumber;
   TBranch        *b_eventnumber;

   //TRIGGER
   TBranch        *b_HLT_Mu15_Photon20_CaloIdL;

   //GLOBAL
   //TBranch        *b_nJetsAKT_pt15;
   //TBranch        *b_nJetsIC5_pt15;
   TBranch        *b_calomet;
   TBranch        *b_met;
   //TBranch        *b_mass;
   TBranch        *b_pthat;
   TBranch        *b_bsposx;
   TBranch        *b_bsposy;
   TBranch        *b_bsposz;

   //JETS IC5
   //TBranch        *b_jetIC5_size;
   //TBranch        *b_jetIC5_pt;
   //TBranch        *b_jetIC5_eta;
   //TBranch        *b_jetIC5_phi;
   //TBranch        *b_jetIC5_em;

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

   // primary vertex information /////////////////////////////////////////////
   vector<float> PVX_ScalingFactor_ttbar;
   vector<float> PVX_ScalingFactor_ztautau;
   if (usePUInfo) {
      stringstream ssPUInfile;
      ssPUInfile << "emu_PUinfo_eleBar" << bar_et << "_eleEnd" << end_et << "_mu" << muon_et << "_" << LumiFactor << "pb-1.root";
      // test if the PU file exists already. If yes -> use it; if no -> create it and fill it
      ifstream iFile(ssPUInfile.str().c_str());
      if (iFile) {
         iFile.close();
         TFile *inputPV = new TFile(ssPUInfile.str().c_str(), "open");
         inputPV->cd();

         TH1F *copy_emuLoose_dataOverTtbar_nValidPv;
         TH1F *copy_emuLoose_dataOverZtautau_nValidPv;

         copy_emuLoose_dataOverTtbar_nValidPv = (TH1F*)inputPV->Get("emuLoose_dataOverTtbar_nValidPv");
         copy_emuLoose_dataOverZtautau_nValidPv = (TH1F*)inputPV->Get("emuLoose_dataOverZtautau_nValidPv");

         for (unsigned int i = 0 ; i < nPVtxMax ; ++i) {
            PVX_ScalingFactor_ttbar.push_back(copy_emuLoose_dataOverTtbar_nValidPv->GetBinContent(i + 1) * nPVtxMax);
            PVX_ScalingFactor_ztautau.push_back(copy_emuLoose_dataOverZtautau_nValidPv->GetBinContent(i + 1) * nPVtxMax);
         }
         cout << "Pile up information from file: " << ssPUInfile.str() << endl;
      } else {
         cout << "Not file with pile up information for this parameter set found. Will create it." << endl;
         generatePUFile = true;
      }
   } else {
      PVX_ScalingFactor_ttbar.clear();
      PVX_ScalingFactor_ztautau.clear();
      for (unsigned int i = 0 ; i < nPVtxMax ; ++i) {
         PVX_ScalingFactor_ttbar.push_back(1.);
         PVX_ScalingFactor_ztautau.push_back(1.);
      }
   }
   if (usePUInfo && generatePUFile) {
      // 1st loop to get VERTEX information
      //for (int p = 0; p < nbFile; ++p) {
      for (int p = 0; p <= ZTT; ++p) { // only data, ttbar and ztautau
         cout << "accessing file " << p + 1 << " for PU information: " << input[p]->GetName() << endl;
         input[p]->cd();

         // Get the TREE and connect the necessary variables
         TTree *thetree;
         thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");

         //HLT TRIGGER BITS
         thetree->SetBranchAddress("HLT_Mu15_Photon20_CaloIdL", &c_HLT_Mu15_Photon20_CaloIdL, &b_HLT_Mu15_Photon20_CaloIdL);
         //PRIM VTX
         thetree->SetBranchAddress("pvsize", &c_pvsize, &b_pvsize);
         thetree->SetBranchAddress("pv_ndof", &c_pv_ndof, &b_pv_ndof);
         thetree->SetBranchAddress("pv_nTracks", &c_pv_nTracks, &b_pv_nTracks);
         //GSF
         thetree->SetBranchAddress("gsf_size", &c_gsf_size, &b_gsf_size);
         thetree->SetBranchAddress("gsf_gsfet", &c_gsf_gsfet, &b_gsf_gsfet);
         thetree->SetBranchAddress("gsfsc_eta", &c_gsfsc_eta, &b_gsfsc_eta);
         //MUONS
         thetree->SetBranchAddress("muon_size", &c_muon_size, &b_muon_size);
         thetree->SetBranchAddress("muon_pt", &c_muon_pt, &b_muon_pt);

         Long64_t nentries = (*thetree).GetEntries();
         cout << nentries << " events" << endl;
         //LOOP OVER EVENTS
         for (unsigned int i = 0; i < nentries; ++i) {
            if (i % 50000 == 0) cout << i << endl;
            thetree->GetEntry(i);

            if (p == DATA && c_HLT_Mu15_Photon20_CaloIdL == 0) continue; // ask MuPhoton trigger bit

            //PRIMARY VTX COUNTING
            unsigned int n_pvValid = 0;
            for (int j = 0; j < c_pvsize; ++j) {
               if (c_pv_ndof[j] > 3 && c_pv_nTracks[j] > 3)
                  n_pvValid++;
            }
            if (p == DATA && n_pvValid < 1) continue;
            if ((p == TTBAR || p == ZTT) && n_pvValid < 1) continue; //ttbar

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

            if (((gsfECAL == 1 && gsfPtMaxB > bar_et)  //BARREL
                 ||
                 (gsfECAL == -1 && gsfPtMaxE > end_et)) //ENDCAP
                && muPtMax > muon_et) {
               if (p == DATA) emuLoose_data_nValidPv->Fill(n_pvValid);
               else if (p == TTBAR) emuLoose_ttbar_nValidPv->Fill(n_pvValid);
               else if (p == ZTT) emuLoose_ztautau_nValidPv->Fill(n_pvValid);
            }
         } // end loop over events
      } // end 1st loop over files

      emuLoose_dataOverTtbar_nValidPv->Divide(emuLoose_data_nValidPv, emuLoose_ttbar_nValidPv);
      emuLoose_dataOverTtbar_nValidPv->Scale(1. / emuLoose_dataOverTtbar_nValidPv->Integral());

      TH1F * ttbarScaled = new TH1F("ttbarScaled", "ttbarScaled", nPVtxMax, 0., nPVtxMax);
      ttbarScaled->Multiply(emuLoose_ttbar_nValidPv, emuLoose_dataOverTtbar_nValidPv);
      double ttbarEvents =  emuLoose_ttbar_nValidPv->Integral();
      double scaledTtbarEvents =  ttbarScaled->Integral();

      cout << "PU normalization factor for ttbar: " << ttbarEvents / scaledTtbarEvents / nPVtxMax << endl;

      emuLoose_dataOverTtbar_nValidPv->Scale(ttbarEvents / scaledTtbarEvents / nPVtxMax);

      if (generatePUFile) {
         stringstream ssPUOutfile;
         ssPUOutfile << "emu_PUinfo_eleBar" << bar_et << "_eleEnd" << end_et << "_mu" << muon_et << "_" << LumiFactor << "pb-1.root";
         TFile *outputPV = new TFile(ssPUOutfile.str().c_str(), "recreate");
         outputPV->cd();

         emuLoose_data_nValidPv->Write();
         emuLoose_ttbar_nValidPv->Write();
         emuLoose_ztautau_nValidPv->Write();
         emuLoose_dataOverTtbar_nValidPv->Write();
         emuLoose_dataOverZtautau_nValidPv->Write();
      }

      for (unsigned int i = 0 ; i < nPVtxMax ; ++i) {
         PVX_ScalingFactor_ttbar.push_back(emuLoose_dataOverTtbar_nValidPv->GetBinContent(i + 1) * nPVtxMax);
         PVX_ScalingFactor_ztautau.push_back(emuLoose_dataOverZtautau_nValidPv->GetBinContent(i + 1) * nPVtxMax);
      }
   }
   ///////////////////////////////////////////////////////////////////////////

   // 2nd loop for analysis 
   // GETTING FILES
   for (int p = 0; p < nbFile; ++p) {
      cout << "accessing file " << p + 1 << ": " << input[p]->GetName() << endl;
      input[p]->cd();

      // Get the TREE
      TTree *thetree;
      thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");

      //RUN ID
      thetree->SetBranchAddress("runnumber", &c_runnumber, &b_runnumber);
      thetree->SetBranchAddress("eventnumber", &c_eventnumber, &b_eventnumber);

      //HLT TRIGGER BITS
      thetree->SetBranchAddress("HLT_Mu15_Photon20_CaloIdL", &c_HLT_Mu15_Photon20_CaloIdL, &b_HLT_Mu15_Photon20_CaloIdL);

      //GLOBAL
      //thetree->SetBranchAddress("nJetsAKT_pt15",&c_nJetsAKT_pt15,&b_nJetsAKT_pt15);
      //thetree->SetBranchAddress("nJetsIC5_pt15",&c_nJetsIC5_pt15,&b_nJetsIC5_pt15);
      thetree->SetBranchAddress("calomet", &c_calomet, &b_calomet);
      thetree->SetBranchAddress("met", &c_met, &b_met);
      //thetree->SetBranchAddress("mass",&c_mass,&b_mass);
      thetree->SetBranchAddress("pthat", &c_pthat, &b_pthat);
      thetree->SetBranchAddress("bsposx", &c_bsposx, &b_bsposx);
      thetree->SetBranchAddress("bsposy", &c_bsposy, &b_bsposy);
      thetree->SetBranchAddress("bsposz", &c_bsposz, &b_bsposz);

      //JETS IC5
      //thetree->SetBranchAddress("jetIC5_size",&c_jetIC5_size,&b_jetIC5_size);
      //thetree->SetBranchAddress("jetIC5_pt",&c_jetIC5_pt,&b_jetIC5_pt);
      //thetree->SetBranchAddress("jetIC5_eta",&c_jetIC5_eta,&b_jetIC5_eta);
      //thetree->SetBranchAddress("jetIC5_phi",&c_jetIC5_phi,&b_jetIC5_phi);
      //thetree->SetBranchAddress("jetIC5_em",&c_jetIC5_em,&b_jetIC5_em);

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

      // set up invariant mass histograms
      emuMass.push_back(new TH1F("emuMass_" + suffix[p], "emuMass_" + suffix[p], 100, 0., 1000.));
      emuMass_LS.push_back(new TH1F("emuMass_LS_" + suffix[p], "emuMass_LS_" + suffix[p], 100, 0., 1000.));
      emuMass_OS.push_back(new TH1F("emuMass_OS_" + suffix[p], "emuMass_OS_" + suffix[p], 100, 0., 1000.));
      emu_mass_accVSgood.push_back(new TH1F("emu_mass_accVSgood", "emu_mass_accVSgood", 100, 0., 1000.));

      //emu_plus_plus.push_back(new TH1F("emu_plus_plus", "emu_plus_plus", 100, 0., 1000.));
      //emu_plus_minus.push_back(new TH1F("emu_plus_minus", "emu_plus_minus", 100, 0., 1000.));
      //emu_minus_plus.push_back(new TH1F("emu_minus_plus", "emu_minus_plus", 100, 0., 1000.));
      //emu_minus_minus.push_back(new TH1F("emu_minus_minus", "emu_minus_minus", 100, 0., 1000.));

      // set up test histograms
      eleMET.push_back(new TH1F("eleMET_" + suffix[p], "eleMET_" + suffix[p], nPVtxMax, 0., nPVtxMax));
      eleNVtx.push_back(new TH1F("eleNVtx_" + suffix[p], "eleNVtx_" + suffix[p], nPVtxMax, 0., nPVtxMax));
      eleDphi.push_back(new TH1F("eleDphi_" + suffix[p], "eleDphi_" + suffix[p], 50, 0., 3.14));
      elePt.push_back(new TH1F("elePt_" + suffix[p], "elePt_" + suffix[p], 50, 0., 500.));
      eleEta.push_back(new TH1F("eleEta_" + suffix[p], "eleEta_" + suffix[p], 25, -2.5, 2.5));
      elePhi.push_back(new TH1F("elePhi_" + suffix[p], "elePhi_" + suffix[p], 25, -3.14, 3.14));
      eleId1.push_back(new TH1F("eleId1_" + suffix[p], "eleId1_" + suffix[p], 50, -0.008, 0.008));
      eleId2.push_back(new TH1F("eleId2_" + suffix[p], "eleId2_" + suffix[p], 50, -0.1, 0.1));
      eleId3.push_back(new TH1F("eleId3_" + suffix[p], "eleId3_" + suffix[p], 50, 0., 0.04));
      eleIso1.push_back(new TH1F("eleIso1_" + suffix[p], "eleIso1_" + suffix[p], 50, 0., 10.));
      eleIso2.push_back(new TH1F("eleIso2_" + suffix[p], "eleIso2_" + suffix[p], 50, 0., 10.));
      eleIso3.push_back(new TH1F("eleIso3_" + suffix[p], "eleIso3_" + suffix[p], 50, 0., 20.));
      
      muIsoCombRel.push_back(new TH1F("muIsoCombRel_" + suffix[p], "muIsoCombRel_" + suffix[p], 50, 0., 0.2));
      muPtEleOPtMu.push_back(new TH1F("muPtEleOPtMu_" + suffix[p], "muPtEleOPtMu_" + suffix[p], 50, 0., 4.0));
      muPtPlusOPtMinus.push_back(new TH1F("muPtPlusOPtMinus_" + suffix[p], "muPtPlusOPtMinus_" + suffix[p], 50, 0., 4.0));
      muPt.push_back(new TH1F("muPt_" + suffix[p], "muPt_" + suffix[p], 50, 0., 500.));
      muEta.push_back(new TH1F("muEta_" + suffix[p], "muEta_" + suffix[p], 25, -2.5, 2.5));
      muPhi.push_back(new TH1F("muPhi_" + suffix[p], "muPhi_" + suffix[p], 25, -3.14, 3.14));
      muId1.push_back(new TH1F("muId1_" + suffix[p], "muId1_" + suffix[p], 50, 0., 10.));
      muId2.push_back(new TH1F("muId2_" + suffix[p], "muId2_" + suffix[p], 32, 0., 32));
      muId3.push_back(new TH1F("muId3_" + suffix[p], "muId3_" + suffix[p], 52, 0., 52));
      muIso1.push_back(new TH1F("muIso1_" + suffix[p], "muIso1_" + suffix[p], 50, 0., 10.));
      muIso2.push_back(new TH1F("muIso2_" + suffix[p], "muIso2_" + suffix[p], 50, 0., 10.));
      muIso3.push_back(new TH1F("muIso3_" + suffix[p], "muIso3_" + suffix[p], 50, 0., 20.));

      Long64_t nentries = (*thetree).GetEntries();
      cout << nentries << " events" << endl;
      //LOOP OVER EVENTS
      for (unsigned int i = 0; i < nentries; ++i) {
         if (i % 50000 == 0) cout << i << endl;
         thetree->GetEntry(i);

         if (p == DATA && c_HLT_Mu15_Photon20_CaloIdL == 0) continue; // ask MuPhoton trigger bit

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
         if (p == DATA && n_pvValid < 1) continue;
         if ((p == TTBAR || p == ZTT) && n_pvValid < 1) continue; //ttbar

         //FILL THE VTX WEIGTH
         float npv_weight = 1.;
         if (p == TTBAR && n_pvValid < nPVtxMax) npv_weight = PVX_ScalingFactor_ttbar[n_pvValid]; //ttbar
         //if (p == ZTT && n_pvValid == v) npv_weight = PVX_ScalingFactor_ztautau[v]; //ztau

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
                && c_gsf_nLostInnerHits[j] <= bar_missInnerHits
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
                && c_muon_dxy_beamSpot[j] < muon_impactParamMax
                && c_muon_nSegmentMatch[j] >= muon_nSegMatchMin
                && c_muon_isTrackerMuon[j]
                && c_muon_trackIso03[j]/ c_muon_pt[j] < muon_relIsoCutMax
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
            if (fabs(c_gsf_phi[GSF_passHEEP[0]] - c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]]) < 0.1) continue;

            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(c_gsf_gsfet[GSF_passHEEP[0]], c_gsf_eta[GSF_passHEEP[0]], c_gsf_phi[GSF_passHEEP[0]], 0.000511);
            mu1.SetPtEtaPhiM(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]], c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]], c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]], 0.10566);

            double invMass = (ele1 + mu1).M();

            //MASS CUT
            if (invMass < minInvMass) continue;

            emuMass[p]->Fill(invMass, npv_weight);

            float CombRelIso = (c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]] + c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]]) / c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]];

            // fill test histograms
            eleMET[p]->Fill(c_calomet, npv_weight);
            eleNVtx[p]->Fill(n_pvValid, npv_weight);
            if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) < 3.14) eleDphi[p]->Fill(fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), npv_weight);
            if (fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]) > 3.14) eleDphi[p]->Fill(6.28 - fabs(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]] - c_gsf_phi[GSF_passHEEP[0]]), npv_weight);
            elePt[p]->Fill(c_gsf_gsfet[GSF_passHEEP[0]], npv_weight);
            eleEta[p]->Fill(c_gsf_eta[GSF_passHEEP[0]], npv_weight);
            elePhi[p]->Fill(c_gsf_phi[GSF_passHEEP[0]], npv_weight);
            eleId1[p]->Fill(c_gsf_deltaeta[GSF_passHEEP[0]], npv_weight);
            eleId2[p]->Fill(c_gsf_deltaphi[GSF_passHEEP[0]], npv_weight);
            eleId3[p]->Fill(c_gsf_sigmaIetaIeta[GSF_passHEEP[0]], npv_weight);
            eleIso1[p]->Fill(c_gsf_ecaliso[GSF_passHEEP[0]], npv_weight);
            eleIso2[p]->Fill(c_gsf_hcaliso1[GSF_passHEEP[0]] + c_gsf_hcaliso2[GSF_passHEEP[0]], npv_weight);
            eleIso3[p]->Fill(c_gsf_trackiso[GSF_passHEEP[0]], npv_weight);

            muIsoCombRel[p]->Fill(CombRelIso,npv_weight);
            muPtEleOPtMu[p]->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) muPtPlusOPtMinus[p]->Fill(c_gsf_gsfet[GSF_passHEEP[0]]/c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) muPtPlusOPtMinus[p]->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]]/c_gsf_gsfet[GSF_passHEEP[0]],npv_weight);
            muPt[p]->Fill(c_muon_pt[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muEta[p]->Fill(c_muon_eta[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muPhi[p]->Fill(c_muon_phi[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muId1[p]->Fill(c_muon_normChi2[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muId2[p]->Fill(c_muon_nhitstrack[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muId3[p]->Fill(c_muon_nhitsmuons[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muIso1[p]->Fill(c_muon_emIso03[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muIso2[p]->Fill(c_muon_hadIso03[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);
            muIso3[p]->Fill(c_muon_trackIso03[MU_passGOOD[MU_leadingPassGOOD]],npv_weight);

            //LIKE SIGN OPPOSITE SIGN
            if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) {
               //(emu_plus_plus[p])->Fill(invMass);
               emuMass_LS[p]->Fill(invMass, npv_weight);
               if (p == DATA) nb_plus_plus++;
            }
            if (c_gsf_charge[GSF_passHEEP[0]] > 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) {
               //(emu_plus_minus[p])->Fill(invMass);
               emuMass_OS[p]->Fill(invMass, npv_weight);
               if (p == DATA) nb_plus_minus++;
            }
            if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] > 0) {
               //(emu_minus_plus[p])->Fill(invMass);
               emuMass_OS[p]->Fill(invMass, npv_weight);
               if (p == DATA) nb_minus_plus++;
            }
            if (c_gsf_charge[GSF_passHEEP[0]] < 0 && c_muon_charge[MU_passGOOD[MU_leadingPassGOOD]] < 0) {
               //(emu_minus_minus[p])->Fill(invMass);
               emuMass_LS[p]->Fill(invMass, npv_weight);
               if (p == DATA) nb_minus_minus++;
            }

            if (p == DATA && invMass > 500.) cout << "HEEP-TightMU event in DATA (M>500) : run number = " << c_runnumber << " , event number = " << c_eventnumber << " , mass = " << invMass <<
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

   //SCALE MC
   for (int p = 1; p < nbFile; p++) {
      //emu_plus_plus[p]->Scale(weight[p] * LumiFactor);
      //emu_plus_minus[p]->Scale(weight[p] * LumiFactor);
      //emu_minus_plus[p]->Scale(weight[p] * LumiFactor);
      //emu_minus_minus[p]->Scale(weight[p] * LumiFactor);

      emuMass[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      emuMass_LS[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      emuMass_OS[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      emu_mass_accVSgood[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);

      eleMET[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleNVtx[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleDphi[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      elePt[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleEta[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      elePhi[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleId1[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleId2[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleId3[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleIso1[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleIso2[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      eleIso3[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);

      muIsoCombRel[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muPtEleOPtMu[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muPtPlusOPtMinus[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muPt[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muEta[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muPhi[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muId1[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muId2[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muId3[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muIso1[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muIso2[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
      muIso3[p]->Scale(weight[p] * LumiFactor * MCemuScaleFactor);
   }

   //PRINT INTEGRAL

   cout << "------------------------------------" << endl;
   cout << "HEEP - TIGHT MU        Lumi = " << LumiFactor << "pb-1" << endl;
   cout << "                       pT   = " << bar_et << "GeV/c" << endl;
   cout << "                       M > " << minInvMass << "GeV/c^2" << endl;
   cout << "nb data      = " << emuMass[DATA]->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "nb TTbar     = " << emuMass[TTBAR]->Integral() << endl;
   cout << "nb Ztautau   = " << emuMass[ZTT]->Integral() << endl;
   cout << "nb WW        = " << emuMass[WW]->Integral() << endl;
   cout << "nb WZ        = " << emuMass[WZ]->Integral() << endl;
   cout << "nb tW        = " << emuMass[TW]->Integral() << endl;
   cout << "- - - - - - - - - - - - - - - - - - " << endl;
   cout << "nb WJets     = " << emuMass[WJET]->Integral() << endl;
   cout << "nb Zmumu     = " << emuMass[ZMM]->Integral() << endl;
   cout << "nb Zee       = " << emuMass[ZEE]->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "TOT ttlike   = " << emuMass[TTBAR]->Integral() + emuMass[ZTT]->Integral() + emuMass[WW]->Integral() + emuMass[WZ]->Integral() + emuMass[TW]->Integral() << endl;
   cout << "TOT contam   = " << emuMass[WJET]->Integral() + emuMass[ZMM]->Integral() + emuMass[ZEE]->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << "TOT MC       = " << emuMass[TTBAR]->Integral() + emuMass[ZTT]->Integral() + emuMass[WW]->Integral() + emuMass[WZ]->Integral() + emuMass[TW]->Integral() + emuMass[WJET]->Integral() + emuMass[ZMM]->Integral() + emuMass[ZEE]->Integral() << endl;
   cout << "------------------------------------" << endl;
   cout << endl;

   //SUM MC
   for (int p = 1; p < nbFile; ++p) {
      for (int q = p + 1; q < nbFile; ++q) {
         //emu_plus_plus[p]->Add(emu_plus_plus[q]);
         //emu_plus_minus[p]->Add(emu_plus_minus[q]);
         //emu_minus_plus[p]->Add(emu_minus_plus[q]);
         //emu_minus_minus[p]->Add(emu_minus_minus[q]);

         emuMass[p]->Add(emuMass[q]);
         emuMass_LS[p]->Add(emuMass_LS[q]);
         emuMass_OS[p]->Add(emuMass_OS[q]);
         emu_mass_accVSgood[p]->Add(emu_mass_accVSgood[q]);

         eleMET[p]->Add(eleMET[q]);
         eleNVtx[p]->Add(eleNVtx[q]);
         eleDphi[p]->Add(eleDphi[q]);
         elePt[p]->Add(elePt[q]);
         eleEta[p]->Add(eleEta[q]);
         elePhi[p]->Add(elePhi[q]);
         eleId1[p]->Add(eleId1[q]);
         eleId2[p]->Add(eleId2[q]);
         eleId3[p]->Add(eleId3[q]);
         eleIso1[p]->Add(eleIso1[q]);
         eleIso2[p]->Add(eleIso2[q]);
         eleIso3[p]->Add(eleIso3[q]);
 
         muIsoCombRel[p]->Add(muIsoCombRel[q]);
         muPtEleOPtMu[p]->Add(muPtEleOPtMu[q]);
         muPtPlusOPtMinus[p]->Add(muPtPlusOPtMinus[q]);
         muPt[p]->Add(muPt[q]);
         muEta[p]->Add(muEta[q]);
         muPhi[p]->Add(muPhi[q]);
         muId1[p]->Add(muId1[q]);
         muId2[p]->Add(muId2[q]);
         muId3[p]->Add(muId3[q]);
         muIso1[p]->Add(muIso1[q]);
         muIso2[p]->Add(muIso2[q]);
         muIso3[p]->Add(muIso3[q]);
      }
   }

   cout << endl;
   cout << endl;
   cout << "------------------------------------" << endl;
   cout << "  Like Sign                         " << endl;
   cout << endl;
   cout << "nb LS DATA   = " << emuMass_LS[DATA]->Integral() << " +- " << sqrt(emuMass_LS[DATA]->Integral()) << endl;
   cout << "nb LS MC     = " << emuMass_LS[1]->Integral() << endl;
   cout << endl;
   cout << "nb OS DATA   = " << emuMass_OS[DATA]->Integral() << " +- " << sqrt(emuMass_OS[DATA]->Integral()) << endl;
   cout << "nb OS MC     = " << emuMass_OS[1]->Integral() << endl;
   cout << "- - - - - - - - - - - - - - - - - - " << endl;
   cout << "nb e+mu+     = " << nb_plus_plus << endl;
   cout << "nb e+mu-     = " << nb_plus_minus << endl;
   cout << "nb e-mu+     = " << nb_minus_plus << endl;
   cout << "nb e-mu-     = " << nb_minus_minus << endl;
   cout << "------------------------------------" << endl;
   cout << "--For ttlike-- (NO QCD CORR)" << endl;
   cout << "INTEGRAL M> 60" << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(7, 100)  - emuMass[6]->Integral(7, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>120" << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(13, 100) - emuMass[6]->Integral(13, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>200" << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(21, 100) - emuMass[6]->Integral(21, 100) << endl;
   cout << "--------------" << endl;
   cout << "--For ttlike-- (NO QCD CORR) + W+Jet + Zmumu" << endl;
   cout << "INTEGRAL M> 60" << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(7, 100)  << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>120" << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(13, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL M>200" << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(21, 100) << endl;
   cout << "--------------" << endl;

   // calculate qcd scale factor
   float QCDScaleFactor = 1 + 2 * (emuMass_LS[DATA]->Integral() - emuMass_LS[1]->Integral()) / emuMass[1]->Integral();
   cout << endl << endl << "calculated QCD scale factor from data: " << QCDScaleFactor << endl;
   // apply QCD scale factor
   if (calcQCDScaleFactor) emuMass[1]->Scale(QCDScaleFactor);
   else {
      cout << "applied preset QCD scale factor:       " << QCD_ScaleFactor << endl;
      emuMass[1]->Scale(QCD_ScaleFactor);
   }

   cout << endl;
   cout << "AFTER CORRECTION WITH QCD FACTOR:" << endl;
   cout << endl;
   cout << "TOTAL INTEGRAL EMU" << endl;
   cout << "nb data     = " << emuMass[DATA]->Integral() << " +- " << sqrt(emuMass[DATA]->Integral()) << endl;
   cout << "nb MC       = " << emuMass[1]->Integral() << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL EMU M>60 " << endl;
   cout << "nb data     = " << emuMass[DATA]->Integral(7, 100) << " +- " << sqrt(emuMass[DATA]->Integral(7, 100)) << endl;
   cout << "nb MC       = " << emuMass[1]->Integral(7, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL EMU M>120" << endl;
   cout << "nb data     = " << emuMass[DATA]->Integral(13, 100) << " +- " << sqrt((emuMass[DATA])->Integral(13, 100)) << endl;
   cout << "nb MC       = " << emuMass[1]->Integral(13, 100) << endl;
   cout << "--------------" << endl;
   cout << "INTEGRAL EMU M>200" << endl;
   cout << "nb data     = " << emuMass[DATA]->Integral(21, 100) << " +- " << sqrt(emuMass[DATA]->Integral(21, 100)) << endl;
   cout << "nb MC       = " << emuMass[1]->Integral(21, 100) << endl;
   cout << "--------------" << endl;

   cout << "--For ttlike-- QCD CORR + W+Jet + Zmumu" << endl;
   cout << "60  < M < 120 " << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(7, 12)  << endl;
   cout << "--------------" << endl;
   cout << "120 < M < 200 " << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(13, 20) << endl;
   cout << "--------------" << endl;
   cout << "M > 200     " << endl;
   cout << "nb tt       = " << emuMass[1]->Integral(21, 100) << endl;
   cout << "--------------" << endl;

   emu_dilepton->Add(emuMass[1]);
   emu_ewk->Add(emuMass[6]);
   //emu_jet->Add(emuMass[8]);

   //WRITING
   stringstream ssOutfile;
   ssOutfile << "testEmuSpec" << LumiFactor << "pb-1.root";
   TFile *output = new TFile(ssOutfile.str().c_str(), "recreate");

   output->cd();

   emu_dilepton->Write();
   emu_ewk->Write();
   emu_jet->Write();

   emuLoose_data_nValidPv->Write();
   emuLoose_ttbar_nValidPv->Write();
   emuLoose_ztautau_nValidPv->Write();

   emuLoose_dataOverTtbar_nValidPv->Write();
   emuLoose_dataOverZtautau_nValidPv->Write();

   for (int p = 0; p < nbFile; ++p) {
      emuMass[p]->Write();
      emuMass_LS[p]->Write();
      emuMass_OS[p]->Write();

      eleMET[p]->Write();
      eleNVtx[p]->Write();
      eleDphi[p]->Write();
      elePt[p]->Write();
      eleEta[p]->Write();
      elePhi[p]->Write();
      eleId1[p]->Write();
      eleId2[p]->Write();
      eleId3[p]->Write();
      eleIso1[p]->Write();
      eleIso2[p]->Write();
      eleIso3[p]->Write();

      muIsoCombRel[p]->Write();
      muPtEleOPtMu[p]->Write();
      muPtPlusOPtMinus[p]->Write();
      muPt[p]->Write();
      muEta[p]->Write();
      muPhi[p]->Write();
      muId1[p]->Write();
      muId2[p]->Write();
      muId3[p]->Write();
      muIso1[p]->Write();
      muIso2[p]->Write();
      muIso3[p]->Write();
   }

   output->Close();
}
