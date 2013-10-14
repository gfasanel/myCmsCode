#define emuSpectrum_cxx
#include "emuSpectrum.h"
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void EmuSpectrum::Loop()
{
   TStopwatch timer;
   timer.Start();
   // parameters /////////////////////////////////////////////////////////////
   float LumiFactor = 19712.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   TParameter<float> lumi("lumi", LumiFactor);

   // DATA file
   TString dataFile = "file:////user/treis/data2013/MuEG_Run2012A+B+C+D-ReReco22Jan2013_1e1muSkim_19780pb-1.root";
   // pile up histogram
   TString puFile = "file:////user/lathomas/data2012/ReReco22Jan2013/pileupHistoMuEG_Run2012ABCDReReco22Jan2013.root";

   string outfileName = "emuSpec";
   //string outfileName = "test";

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
   TParameter<float> eleScaleFactorEB("eleScaleFactorEB", 0.995 * 1.001); // data/MC scale for epsilon_cand (>50GeV) * epsilon_id (at Z peak, pt>35 GeV)
   TParameter<float> eleScaleFactorEE("eleScaleFactorEE", 0.994 * 1.000); // data/MC scale for epsilon_cand (>50GeV) * epsilon_id (at Z peak, pt>35 GeV)
   // muon scale factors from https://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=214870
   TParameter<float> muScaleFactor("muScaleFactor", 0.994);
   TParameter<float> muScaleFactorLowEta("muScaleFactorLowEta", 0.993); // for |eta|<0.9
   TParameter<float> muScaleFactorMidEta("muScaleFactorMidEta", 0.991); // for 0.9<|eta|<1.2
   TParameter<float> muScaleFactorHighEta("muScaleFactorHighEta", 0.998); // for 1.2<|eta|<2.1
   TParameter<float> lumiScaleFactor("lumiScaleFactor", 0.980);  // not used // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1
   TParameter<float> lumiScaleFactorEB("lumiScaleFactorEB", 0.993);  // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1
   TParameter<float> lumiScaleFactorEE("lumiScaleFactorEE", 0.948);  // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1

   // global systematic errors
   TParameter<float> systErrLumi("systErrLumi", 0.026);
   TParameter<float> systErrEff("systErrEff", 0.007); // muon err (0.002) & ele err (0.007)

   bool usePUInfo = true;
   bool lowMassPuOnly = false;
   float puMassCut = 120.;
   // selection cuts /////////////////////////////////////////////////////////
   float minInvMass = 0.;

   TH1::SetDefaultSumw2(kTRUE);

   ///////////////////////////////////////////////////////////////////////////
   // INPUT FILES
   vector<pair<TFile *, float> > input;
   vector<bool> storeGenMTtbar;
   THashList systErrMCs;
   THashList mcWeights;
   THashList eleVetoRatios;

   // DATA
   TFile *inData = TFile::Open(dataFile);
   input.push_back(make_pair(inData, 1.)); //DATA
   storeGenMTtbar.push_back(0);

   // MC
   TFile *inTTbar = TFile::Open("file:////user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
   //input.push_back(make_pair(inTTbar, 225.197 / 28150723.)); // NLO
   input.push_back(make_pair(inTTbar, 245.8 / 28150723.));  // NNLO http://arxiv.org/pdf/1303.6254.pdf
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar", 0.034));
   storeGenMTtbar.push_back(1);

   TFile *inTTbar700to1000 = TFile::Open("file:////user/treis/mcsamples/TT_Mtt-700to1000_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_3082812ev.root");
   input.push_back(make_pair(inTTbar700to1000, 15.614 / 3082812. * 245.8./211.));  // ttbar  mtt 700to1000
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar700to1000", 0.15));
   storeGenMTtbar.push_back(1);

   TFile *inTTbar1000up = TFile::Open("file:////user/treis/mcsamples/TT_Mtt-1000toInf_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_46_1249111ev.root");
   input.push_back(make_pair(inTTbar1000up, 2.954 / 1249111. * 245.8/211.));  // ttbar  mtt>1000
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

   TFile *inTTbarMg = TFile::Open("file:////user/treis/mcsamples/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_8288533ev.root");
   input.push_back(make_pair(inTTbarMg, 234./ 8288533.)); //TTjets from MadGraph
   systErrMCs.Add(new TParameter<float>("systErrMcTtJets", 0.067));
   storeGenMTtbar.push_back(1);

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

   TFile *inSig1 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   input.push_back(make_pair(inSig1, 0.000239 / 10000.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig1", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig2 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-750_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig2, 8.956e-5 / 9999.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig2", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig3 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig3, 3.867e-5 / 9999.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig3", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig4 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig4, 1.778e-5 / 9999.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig4", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig5 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   input.push_back(make_pair(inSig5, 8.501e-6 / 10000.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig5", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig6 = TFile::Open("file://///user/treis/mcsamples/ZprimeToEMu_M-1750_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   input.push_back(make_pair(inSig6, 4.124e-6 / 10000.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig6", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig7 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig7, 2.014e-6 / 9999.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig7", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig8 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   input.push_back(make_pair(inSig8, 4.724e-7 / 10000.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig8", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig9 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig9, 1.059e-7 / 9999.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig9", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig10 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9997ev.root");
   input.push_back(make_pair(inSig10, 2.198e-8 / 9997.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig10", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig11 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9993ev.root");
   input.push_back(make_pair(inSig11, 4.157e-9 / 9993.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig11", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig12 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9959ev.root");
   input.push_back(make_pair(inSig12, 1.359e-10 / 9959.)); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig12", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig13 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig13, 0.000239 / 9999.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig13", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig14 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   input.push_back(make_pair(inSig14, 8.956e-5 / 10000.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig14", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig15 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9998ev.root");
   input.push_back(make_pair(inSig15, 3.862e-5 / 9998.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig15", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig16 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9998ev.root");
   input.push_back(make_pair(inSig16, 1.781e-5 / 9998.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig16", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig17 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9997ev.root");
   input.push_back(make_pair(inSig17, 8.503e-6 / 9997.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig17", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig18 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9997ev.root");
   input.push_back(make_pair(inSig18, 4.127e-6 / 9997.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig18", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig19 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig19, 2.014e-6 / 9999.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig19", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig20 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   input.push_back(make_pair(inSig20, 4.735e-7 / 9999.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig20", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig21 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   input.push_back(make_pair(inSig21, 1.059e-7 / 10000.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig21", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig22 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9898ev.root");
   input.push_back(make_pair(inSig22, 2.197e-8 / 9898.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig22", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig23 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9998ev.root");
   input.push_back(make_pair(inSig23, 4.155e-9 / 9998.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig23", 0.));
   storeGenMTtbar.push_back(0);

   TFile *inSig24 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9966ev.root");
   input.push_back(make_pair(inSig24, 1.293e-10 / 9966.)); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig24", 0.));
   storeGenMTtbar.push_back(0);

   int nbFile = input.size();
   ///////////////////////////////////////////////////////////////////////////

   // counting variables
   int nb_plus_plus = 0;
   int nb_plus_minus = 0;
   int nb_minus_plus = 0;
   int nb_minus_minus = 0;

   // data pileup histogram
   TFile *dataPuInfile = TFile::Open(puFile);
   dataPuInfile->Cd("");
   TH1F *puData = (TH1F *)gDirectory->Get("pileup");
   puData->SetDirectory(0);
   puData->SetName("puData");
   dataPuInfile->Close();
   TH1F *puDataNorm = (TH1F *)puData->DrawNormalized()->Clone("puDataNorm");

   stringstream ssGoodEmuFileName;
   ssGoodEmuFileName << "goodEmuEvents" << LumiFactor << "pb-1.root";
   TFile *goodEvFile = new TFile(ssGoodEmuFileName.str().c_str(), "recreate");
   goodEvFile->cd();
   TTree *eleDataTree = new TTree("eleDataTree", "eleDataTree");
   int evRegion = -1;
   float emuMass = 0.;
   eleDataTree->Branch("runnr", &runnumber, "runnr/I");
   eleDataTree->Branch("eventnr", &eventnumber, "eventnr/i");
   eleDataTree->Branch("lumiSec", &luminosityBlock, "lumiSec/I");
   eleDataTree->Branch("evtRegion", &evRegion, "evtRegion/I");
   eleDataTree->Branch("mass", &emuMass, "mass/F");

   unsigned int dataTrig[3] = {0, 0, 0};
   unsigned int dataEntries = 0;

   // output file
   stringstream ssOutfile;
   ssOutfile << outfileName << "_" << LumiFactor << "pb-1.root";
   TFile *output = new TFile(ssOutfile.str().c_str(), "recreate");
   output->cd();

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
   ///////////////////////////////////////////////////////////////////////////
   //LOOP OVER FILES 
   ///////////////////////////////////////////////////////////////////////////
   for (int p = 0; p < nbFile; ++p) {
      cout << "accessing file " << p + 1 << ": " << input[p].first->GetName() << endl;
      input[p].first->Cd("");
      TTree *thetree = (TTree*)input[p].first->Get("gsfcheckerjob/tree");
      Init(thetree);
      Long64_t nentries = fChain->GetEntriesFast();
      cout << nentries << " events" << endl;

      // get the histogram with the true number of vertices
      TH1F *puMc = new TH1F("puMc" + suffix[p], "puMc" + suffix[p], 100, 0., 100.);
      thetree->Draw("trueNVtx>>puMc" + suffix[p]);
      TH1F *puMcNorm = new TH1F("dummy" + suffix[p], "dummy" + suffix[p], 100, 0., 100.);
      if (p > 0) puMcNorm = (TH1F *)puMc->DrawNormalized();
      puMcNorm->SetName("puMc" + suffix[p] + "Norm");
      puMcNorm->SetTitle("puMc" + suffix[p] + "Norm");
      // calculate the pu weights
      TH1F *puWeights = new TH1F("puWeight" + suffix[p], "puWeight" + suffix[p], 100, 0., 100.);
      if (p > 0) puWeights->Divide(puDataNorm, puMcNorm);

      // tree with event data
      bool passTrg = false;
      bool passHeep = false;
      float emuInvMass = 0.;
      float trueMass = 0.;
      int evtRegion;
      int eCharge;
      int muCharge;
      float totWeight = 1.;
      float puWeight = 1.;
      TTree *emuTree = new TTree("emuTree_" + suffix[p], "emuTree_" + suffix[p]);
      emuTree->Branch("runnr", &runnumber, "runnr/i");
      emuTree->Branch("eventnr", &eventnumber, "eventnr/i");
      emuTree->Branch("lumiSec", &luminosityBlock, "lumiSec/i");
      emuTree->Branch("passTrg", &passTrg, "passTrg/O");
      emuTree->Branch("mass", &emuInvMass, "mass/F");
      if (p > 0) emuTree->Branch("trueMass", &trueMass, "trueMass/F");
      emuTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
      emuTree->Branch("eCharge", &eCharge, "eCharge/I");
      emuTree->Branch("muCharge", &muCharge, "muCharge/I");
      emuTree->Branch("weight", &totWeight, "weight/F");
      emuTree->Branch("puWeight", &puWeight, "puWeight/F");
      if (storeGenMTtbar[p]) emuTree->Branch("genMTtbar", &genPair_mass, "genMTtbar/F");
      // control variables
      float nVtx = 0.;
      float eDzMinusMuDz = 0.;
      float eDxyMinusMuDxy = 0.;
      float dPhi = 0.;
      float dEta = 0.;
      float eleEt = 0.;
      float eleEta = 0.;
      float elePhi = 0.;
      float eleDEta = 0.;
      float eleDPhi = 0.;
      float eleHOE = 0.;
      float eleE1x5overE5x5 = 0.;
      float eleE2x5overE5x5 = 0.;
      float eleSigmaIEIE = 0.;
      float eleEcalIso = 0.;
      float eleHcalIso1 = 0.;
      float eleHcalIso2 = 0.;
      float eleHeepIso = -1.;
      float eleTrkIso = 0.;
      float eleLostHits = 0.;
      float eleDZFstPVtx = 0.;
      float eleDXYFstPVtx = 0.;
      float etEleOPtMu = 0.;
      float lepPtPlusOPtMinus = 0.;
      float eChTimesMuCh = 0;
      float muPt = 0.;
      float muPtErr = 0.;
      float muEta = 0.;
      float muPhi = 0.;
      float muHitLayers = 0.;
      float muTrkHits = 0.;
      float muPxlHits = 0.;
      float muMuHits = 0.;
      float muDZFstPVtx = 0.;
      float muDXYFstPVtx = 0.;
      float muNSeg = 0.;
      float muIsoCombRel = 0.;
      float muTrkIso03 = 0.;
      float numOfJets = 0.;
      float numOfJetsPt20 = 0.;
      float numOfJetsPt30 = 0.;
      emuTree->Branch("pfMet", &pfmet, "pfMet/F");
      emuTree->Branch("nVtx", &nVtx, "nVtx/F");
      emuTree->Branch("eDzMinusMuDz", &eDzMinusMuDz, "eDzMinusMuDz/F");
      emuTree->Branch("eDxyMinusMuDxy", &eDxyMinusMuDxy, "eDxyMinusMuDxy/F");
      emuTree->Branch("rho", &rho, "rho/F");
      emuTree->Branch("dPhi", &dPhi, "dPhi/F");
      emuTree->Branch("dEta", &dEta, "dEta/F");
      emuTree->Branch("eleEt", &eleEt, "eleEt/F");
      emuTree->Branch("eleEta", &eleEta, "eleEta/F");
      emuTree->Branch("elePhi", &elePhi, "elePhi/F");
      emuTree->Branch("eleDEta", &eleDEta, "eleDEta/F");
      emuTree->Branch("eleDPhi", &eleDPhi, "eleDPhi/F");
      emuTree->Branch("eleHOE", &eleHOE, "eleHOE/F");
      emuTree->Branch("eleE1x5overE5x5", &eleE1x5overE5x5, "eleE1x5overE5x5/F");
      emuTree->Branch("eleE2x5overE5x5", &eleE2x5overE5x5, "eleE2x5overE5x5/F");
      emuTree->Branch("eleSigmaIEIE", &eleSigmaIEIE, "eleSigmaIEIE/F");
      emuTree->Branch("eleEcalIso", &eleEcalIso, "eleEcalIso/F");
      emuTree->Branch("eleHcalIso1", &eleHcalIso1, "eleHcalIso1/F");
      emuTree->Branch("eleHcalIso2", &eleHcalIso2, "eleHcalIso2/F");
      emuTree->Branch("eleHeepIso", &eleHeepIso, "eleHeepIso/F");
      emuTree->Branch("eleTrkIso", &eleTrkIso, "eleTrkIso/F");
      emuTree->Branch("eleLostHits", &eleLostHits, "eleLostHits/F");
      emuTree->Branch("eleDZFstPVtx", &eleDZFstPVtx, "eleDZFstPVtx/F");
      emuTree->Branch("eleDXYFstPVtx", &eleDXYFstPVtx, "eleDXYFstPVtx/F");
      emuTree->Branch("etEleOPtMu", &etEleOPtMu, "etEleOPtMu/F");
      emuTree->Branch("lepPtPlusOPtMinus", &lepPtPlusOPtMinus, "lepPtPlusOPtMinus/F");
      emuTree->Branch("eChTimesMuCh", &eChTimesMuCh, "eChTimesMuCh/F");
      emuTree->Branch("muPt", &muPt, "muPt/F");
      emuTree->Branch("muPtErr", &muPtErr, "muPtErr/F");
      emuTree->Branch("muEta", &muEta, "muEta/F");
      emuTree->Branch("muPhi", &muPhi, "muPhi/F");
      emuTree->Branch("muHitLayers", &muHitLayers, "muHitLayers/F");
      emuTree->Branch("muTrkHits", &muTrkHits, "muTrkHits/F");
      emuTree->Branch("muPxlHits", &muPxlHits, "muPxlHits/F");
      emuTree->Branch("muMuHits", &muMuHits, "muMuHits/F");
      emuTree->Branch("muDZFstPVtx", &muDZFstPVtx, "muDZFstPVtx/F");
      emuTree->Branch("muDXYFstPVtx", &muDXYFstPVtx, "muDXYFstPVtx/F");
      emuTree->Branch("muNSeg", &muNSeg, "muNSeg/F");
      emuTree->Branch("muTrkIso03", &muTrkIso03, "muTrkIso03/F");
      emuTree->Branch("muIsoCombRel", &muIsoCombRel, "muIsoCombRel/F");
      emuTree->Branch("numOfJets", &numOfJets, "numOfJets/F");
      emuTree->Branch("numOfJetsPt20", &numOfJetsPt20, "numOfJetsPt20/F");
      emuTree->Branch("numOfJetsPt30", &numOfJetsPt30, "numOfJetsPt30/F");

      TTree *frEmuTree = new TTree("frEmuTree_" + suffix[p], "frEmuTree_" + suffix[p]);
      float fakeRate = 0.;
      {
         frEmuTree->Branch("passTrg", &passTrg, "passTrg/O");
         frEmuTree->Branch("passHeep", &passHeep, "passHeep/O");
         frEmuTree->Branch("mass", &emuInvMass, "mass/F");
         frEmuTree->Branch("trueMass", &trueMass, "trueMass/F");
         frEmuTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
         frEmuTree->Branch("eCharge", &eCharge, "eCharge/I");
         frEmuTree->Branch("muCharge", &muCharge, "muCharge/I");
         frEmuTree->Branch("fakeRate", &fakeRate, "fakeRate/F");
         frEmuTree->Branch("puWeight", &puWeight, "puWeight/F");
         if (storeGenMTtbar[p]) frEmuTree->Branch("genMTtbar", &genPair_mass, "genMTtbar/F");
         frEmuTree->Branch("pfMet", &pfmet, "pfMet/F");
         frEmuTree->Branch("nVtx", &nVtx, "nVtx/F");
         frEmuTree->Branch("eDzMinusMuDz", &eDzMinusMuDz, "eDzMinusMuDz/F");
         frEmuTree->Branch("eDxyMinusMuDxy", &eDxyMinusMuDxy, "eDxyMinusMuDxy/F");
         frEmuTree->Branch("rho", &rho, "rho/F");
         frEmuTree->Branch("dPhi", &dPhi, "dPhi/F");
         frEmuTree->Branch("dEta", &dEta, "dEta/F");
         frEmuTree->Branch("eleEt", &eleEt, "eleEt/F");
         frEmuTree->Branch("eleEta", &eleEta, "eleEta/F");
         frEmuTree->Branch("elePhi", &elePhi, "elePhi/F");
         frEmuTree->Branch("eleDEta", &eleDEta, "eleDEta/F");
         frEmuTree->Branch("eleDPhi", &eleDPhi, "eleDPhi/F");
         frEmuTree->Branch("eleHOE", &eleHOE, "eleHOE/F");
         frEmuTree->Branch("eleE1x5overE5x5", &eleE1x5overE5x5, "eleE1x5overE5x5/F");
         frEmuTree->Branch("eleE2x5overE5x5", &eleE2x5overE5x5, "eleE2x5overE5x5/F");
         frEmuTree->Branch("eleSigmaIEIE", &eleSigmaIEIE, "eleSigmaIEIE/F");
         frEmuTree->Branch("eleEcalIso", &eleEcalIso, "eleEcalIso/F");
         frEmuTree->Branch("eleHcalIso1", &eleHcalIso1, "eleHcalIso1/F");
         frEmuTree->Branch("eleHcalIso2", &eleHcalIso2, "eleHcalIso2/F");
         frEmuTree->Branch("eleHeepIso", &eleHeepIso, "eleHeepIso/F");
         frEmuTree->Branch("eleTrkIso", &eleTrkIso, "eleTrkIso/F");
         frEmuTree->Branch("eleLostHits", &eleLostHits, "eleLostHits/F");
         frEmuTree->Branch("eleDZFstPVtx", &eleDZFstPVtx, "eleDZFstPVtx/F");
         frEmuTree->Branch("eleDXYFstPVtx", &eleDXYFstPVtx, "eleDXYFstPVtx/F");
         frEmuTree->Branch("etEleOPtMu", &etEleOPtMu, "etEleOPtMu/F");
         frEmuTree->Branch("lepPtPlusOPtMinus", &lepPtPlusOPtMinus, "lepPtPlusOPtMinus/F");
         frEmuTree->Branch("eChTimesMuCh", &eChTimesMuCh, "eChTimesMuCh/F");
         frEmuTree->Branch("muPt", &muPt, "muPt/F");
         frEmuTree->Branch("muPtErr", &muPtErr, "muPtErr/F");
         frEmuTree->Branch("muEta", &muEta, "muEta/F");
         frEmuTree->Branch("muPhi", &muPhi, "muPhi/F");
         frEmuTree->Branch("muHitLayers", &muHitLayers, "muHitLayers/F");
         frEmuTree->Branch("muTrkHits", &muTrkHits, "muTrkHits/F");
         frEmuTree->Branch("muPxlHits", &muPxlHits, "muPxlHits/F");
         frEmuTree->Branch("muMuHits", &muMuHits, "muMuHits/F");
         frEmuTree->Branch("muDZFstPVtx", &muDZFstPVtx, "muDZFstPVtx/F");
         frEmuTree->Branch("muDXYFstPVtx", &muDXYFstPVtx, "muDXYFstPVtx/F");
         frEmuTree->Branch("muNSeg", &muNSeg, "muNSeg/F");
         frEmuTree->Branch("muTrkIso03", &muTrkIso03, "muTrkIso03/F");
         frEmuTree->Branch("muIsoCombRel", &muIsoCombRel, "muIsoCombRel/F");
         frEmuTree->Branch("numOfJets", &numOfJets, "numOfJets/F");
         frEmuTree->Branch("numOfJetsPt20", &numOfJetsPt20, "numOfJetsPt20/F");
         frEmuTree->Branch("numOfJetsPt30", &numOfJetsPt30, "numOfJetsPt30/F");

         if (p == 0) {
            cout << "-----------------------------------------------------------------------------------------------------------" << endl;
            cout << "M_emu > 600 GeV/c^2       |     run:lumi:event    |   M_emu  |"
                 << " muon pt |"
                 << " muon eta |"
                 //<< "  muon phi  |"
                 //<< "  muon trackIso3  |"
                 //<< "  muon emIso03  |"
                 //<< "  muon hadIso03  |"
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
      }
      unsigned int trig[3] = {0, 0, 0};
      unsigned int evCounter = 0;
      unsigned int goodHeepCounter = 0;
      unsigned int goodMuCounter = 0;
      unsigned int moreHeepCounter = 0;
      unsigned int moreMuCounter = 0;
      unsigned int moreHeepMuCounter = 0;
      /////////////////////////////////////////////////////////////////////////////////////////////
      //LOOP OVER EVENTS
      /////////////////////////////////////////////////////////////////////////////////////////////
      Long64_t nbytes = 0, nb = 0;
      //nentries = 10000;
      for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
         emuMass = 0.;
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         // if (Cut(ientry) < 0) continue;
         if (jentry % 50000 == 0) cout << "Processing event " << jentry << endl;
         thetree->GetEntry(jentry);

         // at least one gsf electron and one muon above the threshold
         if (gsf_size < 1 || muon_size < 1) continue;

         // get trigger info
         int prescale = 0;
         passTrg = true;
         if (p == 0 && Trigger(prescale, dataTrig) < 1) passTrg = false;
         if (p != 0) {
            if (jentry < nentries * dataTrig[0] / dataEntries) {
               if (Trigger(prescale, trig, 1) < 1) passTrg = false;
            } else if (jentry < nentries * (dataTrig[0] + dataTrig[1]) / dataEntries) {
               if (Trigger(prescale, trig, 2) < 1) passTrg = false;
            } else {
               if (Trigger(prescale, trig, 3) < 1) passTrg = false;
            }
         }

         // set the PU weight
         if (usePUInfo && p > 0) puWeight = puWeights->GetBinContent(puWeights->FindBin(trueNVtx));
         float weight = puWeight;

         vector<int> GSF_passHEEP;
         vector<int> GSF_passFrPre;
         vector<int> MU_passGOOD;
         /////////////////////////////////////////////////////////////////////////////////////////////
         //loop over electrons
         for (int j = 0; j < gsf_size; ++j) {
            //cleaning: fake electrons from muons skim electrons if there is a muon within dR<0.1
            bool fakeEle = false;
            for (int k = 0; k < muon_size; ++k) {
               float DeltaR = deltaR(gsf_eta[j], gsf_phi[j], muon_eta[k], muon_phi[k]);
               if (DeltaR < 0.1) {
                  fakeEle = true;
                  break;
               }
            }
            if (fakeEle) continue;

            if (PassHEEP(j)) GSF_passHEEP.push_back(j);
            if ((fabs(gsfsc_eta[j]) < 1.442 && gsf_gsfet[j] > bar_et) || (fabs(gsfsc_eta[j]) > 1.56 && fabs(gsfsc_eta[j]) < 2.5 && gsf_gsfet[j] > end_et)) {
               if (PassFRPreSel(j)) GSF_passFrPre.push_back(j);
            }
         }

         //loop over muons
         for (int j = 0; j < muon_size; ++j) {
            if (PassHighPtMu(j)) MU_passGOOD.push_back(j);
         }

         if (GSF_passHEEP.size() > 1) ++moreHeepCounter;

         // veto when there are more than one good muon
         if (MU_passGOOD.size() > 1) {
            //cout << "Muon masses: ";
            //vector<double> muInvMasses;
            //for (unsigned int mu1It = 0; mu1It < MU_passGOOD.size(); ++mu1It) {
            //   for (unsigned int mu2It = MU_passGOOD.size()-1; mu2It > mu1It; --mu2It) {
            //      TLorentzVector mu1;
            //      TLorentzVector mu2;

            //      mu1.SetPtEtaPhiM(muon_pt[MU_passGOOD[mu1It]], muon_eta[MU_passGOOD[mu1It]], muon_phi[MU_passGOOD[mu1It]], 0.10566);
            //      mu2.SetPtEtaPhiM(muon_pt[MU_passGOOD[mu2It]], muon_eta[MU_passGOOD[mu2It]], muon_phi[MU_passGOOD[mu2It]], 0.10566);

            //      muInvMasses.push_back((mu1 + mu2).M());
            //      cout << muInvMasses.back() << " ";
            //   }
            //}
            ////cout << MU_passGOOD.size() << " high pt muons found" << endl;
            ++moreMuCounter;
            if (GSF_passHEEP.size() > 1) ++moreHeepMuCounter;
            continue;
         }
         //if (GSF_passHEEP.size() > 1 || MU_passGOOD.size() > 1) cout << endl;

         //GSF-nonHEEP ele + GOOD muon
         if (GSF_passFrPre.size() > 0 && MU_passGOOD.size() > 0) {
            TLorentzVector ele1;
            TLorentzVector mu1;
   
            ele1.SetPtEtaPhiM(gsf_gsfet[GSF_passFrPre[0]], gsf_eta[GSF_passFrPre[0]], gsf_phi[GSF_passFrPre[0]], 0.000511);
            mu1.SetPtEtaPhiM(muon_pt[MU_passGOOD[0]], muon_eta[MU_passGOOD[0]], muon_phi[MU_passGOOD[0]], 0.10566);
   
            double invMass = (ele1 + mu1).M();
   
            //MASS CUT
            if (invMass > minInvMass) {
               fakeRate = FakeRate(gsf_gsfet[GSF_passFrPre[0]], gsfsc_eta[GSF_passFrPre[0]]);
   
               float CombRelIso = (muon_emIso03[MU_passGOOD[0]] + muon_hadIso03[MU_passGOOD[0]] + muon_trackIso03[MU_passGOOD[0]]) / muon_pt[MU_passGOOD[0]];
   
               int jetsPt20 = 0;
               int jetsPt30 = 0;
               for (int jetIt = 0; jetIt < JetColl_size; ++jetIt) {
                  if (Jet_pt[jetIt] > 20.) ++jetsPt20;
                  if (Jet_pt[jetIt] > 30.) ++jetsPt30;
               }
   
               // fill the data tree
               passHeep = PassHEEP(GSF_passFrPre[0]);
               emuInvMass = invMass;
               if (p > 0) trueMass = genelemom_mass[0];
               eCharge = gsf_charge[GSF_passFrPre[0]];
               muCharge = muon_charge[MU_passGOOD[0]];
               totWeight = weight;
               if (p > 0) totWeight *= input[p].second * LumiFactor;
               if (fabs(gsfsc_eta[GSF_passFrPre[0]]) < 1.442) evtRegion = 0;
               else if (fabs(gsfsc_eta[GSF_passFrPre[0]]) > 1.56 && fabs(gsfsc_eta[GSF_passFrPre[0]]) < 2.5) evtRegion = 1;
               else evtRegion = -1;
               // fill control variables
               nVtx = pvsize;
               eDzMinusMuDz = fabs(gsf_dz_firstPVtx[GSF_passFrPre[0]] - muon_dz_firstPVtx[MU_passGOOD[0]]);
               eDxyMinusMuDxy = fabs(gsf_dxy_firstPVtx[GSF_passFrPre[0]] - muon_dxy_firstPVtx[MU_passGOOD[0]]);
               if (fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passFrPre[0]]) < 3.14) dPhi = fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passFrPre[0]]);
               if (fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passFrPre[0]]) > 3.14) dPhi = 6.28 - fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passFrPre[0]]);
               dEta = fabs(muon_eta[MU_passGOOD[0]] - gsf_eta[GSF_passFrPre[0]]);
               eleEt = gsf_gsfet[GSF_passFrPre[0]];
               eleEta = gsf_eta[GSF_passFrPre[0]];
               elePhi = gsf_phi[GSF_passFrPre[0]];
               eleDEta = gsf_deltaeta[GSF_passFrPre[0]];
               eleDPhi = gsf_deltaphi[GSF_passFrPre[0]];
               eleHOE = gsf_hovere[GSF_passFrPre[0]];
               eleE1x5overE5x5 = gsf_e1x5overe5x5[GSF_passFrPre[0]];
               eleE2x5overE5x5 = gsf_e2x5overe5x5[GSF_passFrPre[0]];
               eleSigmaIEIE = gsf_sigmaIetaIeta[GSF_passFrPre[0]];
               eleEcalIso = gsf_ecaliso[GSF_passFrPre[0]];
               eleHcalIso1 = gsf_hcaliso1[GSF_passFrPre[0]];
               eleHcalIso2 = gsf_hcaliso2[GSF_passFrPre[0]];
               if (evtRegion == 0) eleHeepIso = -1. * (gsf_ecaliso[GSF_passFrPre[0]] + gsf_hcaliso1[GSF_passFrPre[0]] - 2. - rho * 0.28 - 0.03 * gsf_gsfet[GSF_passFrPre[0]]);
               else if (evtRegion == 1 && gsf_gsfet[GSF_passFrPre[0]] < 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passFrPre[0]] + gsf_hcaliso1[GSF_passFrPre[0]] - 2.5 - rho * 0.28);
               else if (evtRegion == 1 && gsf_gsfet[GSF_passFrPre[0]] >= 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passFrPre[0]] + gsf_hcaliso1[GSF_passFrPre[0]] - 2.5 - rho * 0.28 - 0.03 * (gsf_gsfet[GSF_passFrPre[0]] - 50.));
               eleTrkIso = gsf_trackiso[GSF_passFrPre[0]];
               eleLostHits = gsf_nLostInnerHits[GSF_passFrPre[0]];
               eleDXYFstPVtx = gsf_dxy_firstPVtx[GSF_passFrPre[0]];
               eleDZFstPVtx = gsf_dz_firstPVtx[GSF_passFrPre[0]];
               etEleOPtMu = gsf_gsfet[GSF_passFrPre[0]]/muon_pt[MU_passGOOD[0]];
               if (gsf_charge[GSF_passFrPre[0]] > 0 && muon_charge[MU_passGOOD[0]] < 0) lepPtPlusOPtMinus = gsf_gsfet[GSF_passFrPre[0]]/muon_pt[MU_passGOOD[0]];
               else if (gsf_charge[GSF_passFrPre[0]] < 0 && muon_charge[MU_passGOOD[0]] > 0) lepPtPlusOPtMinus = muon_pt[MU_passGOOD[0]]/gsf_gsfet[GSF_passFrPre[0]];
               else lepPtPlusOPtMinus = 0.;
               eChTimesMuCh = eCharge * muCharge;
               muPt = muon_pt[MU_passGOOD[0]];
               muPtErr = muon_ptError[MU_passGOOD[0]];
               muEta = muon_eta[MU_passGOOD[0]];
               muPhi = muon_phi[MU_passGOOD[0]];
               muHitLayers = muon_nlayerswithhits[MU_passGOOD[0]];
               muTrkHits = muon_nhitstrack[MU_passGOOD[0]];
               muPxlHits = muon_nhitspixel[MU_passGOOD[0]];
               muMuHits = muon_nhitsmuons[MU_passGOOD[0]];
               muDZFstPVtx = muon_dz_firstPVtx[MU_passGOOD[0]];
               muDXYFstPVtx = muon_dxy_firstPVtx[MU_passGOOD[0]];
               muNSeg = muon_nSegmentMatch[MU_passGOOD[0]];
               muTrkIso03 = muon_trackIso03[MU_passGOOD[0]];
               muIsoCombRel = CombRelIso;
               numOfJets = JetColl_size;
               numOfJetsPt20 = jetsPt20;
               numOfJetsPt30 = jetsPt30;
   
               frEmuTree->Fill();
            }
         }

         //HEEP ele + GOOD muon
         if (GSF_passHEEP.size() > 0 && MU_passGOOD.size() > 0) {
            TLorentzVector ele1;
            TLorentzVector mu1;

            ele1.SetPtEtaPhiM(gsf_gsfet[GSF_passHEEP[0]], gsf_eta[GSF_passHEEP[0]], gsf_phi[GSF_passHEEP[0]], 0.000511);
            mu1.SetPtEtaPhiM(muon_pt[MU_passGOOD[0]], muon_eta[MU_passGOOD[0]], muon_phi[MU_passGOOD[0]], 0.10566);

            double invMass = (ele1 + mu1).M();

            // vertex matching of tracks by dz
            //if (fabs(gsf_dz_beamSpot[GSF_passHEEP[0]] - muon_dz_beamSpot[MU_passGOOD[0]]) > 0.05) continue;

            //MASS CUT
            if (invMass < minInvMass) continue;
            ++goodHeepCounter;
            ++goodMuCounter;

            if (GSF_passHEEP.size() > 1) {
               //cout << "Electron masses: ";
               //vector<double> eleInvMasses;
               //for (unsigned int heep1It = 0; heep1It < GSF_passHEEP.size(); ++heep1It) {
               //   for (unsigned int heep2It = GSF_passHEEP.size()-1; heep2It > heep1It; --heep2It) {
               //      TLorentzVector ele1;
               //      TLorentzVector ele2;

               //      ele1.SetPtEtaPhiM(gsf_gsfet[GSF_passHEEP[heep1It]], gsf_eta[GSF_passHEEP[heep1It]], gsf_phi[GSF_passHEEP[heep1It]], 0.000511);
               //      ele2.SetPtEtaPhiM(gsf_gsfet[GSF_passHEEP[heep2It]], gsf_eta[GSF_passHEEP[heep2It]], gsf_phi[GSF_passHEEP[heep2It]], 0.000511);

               //      eleInvMasses.push_back((ele1 + ele2).M());
               //      cout << eleInvMasses.back() << " ";
               //   }
               //}
               ////cout << GSF_passHEEP.size() << " HEEP electrons found" << endl;
               //++moreHeepCounter;
               continue;
            }

            ++evCounter;

            float CombRelIso = (muon_emIso03[MU_passGOOD[0]] + muon_hadIso03[MU_passGOOD[0]] + muon_trackIso03[MU_passGOOD[0]]) / muon_pt[MU_passGOOD[0]];

            // set correction factors according to detector region
            float Elec_ScaleFactor = 1.;
            float Lumi_ScaleFactor = 1.;
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.442) {
              Elec_ScaleFactor = eleScaleFactorEB.GetVal();
              Lumi_ScaleFactor = lumiScaleFactorEB.GetVal();
            }
            else if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(gsfsc_eta[GSF_passHEEP[0]]) < 2.5) {
              Elec_ScaleFactor = eleScaleFactorEE.GetVal();
              Lumi_ScaleFactor = lumiScaleFactorEE.GetVal();
            }
            if (lowMassPuOnly && invMass > puMassCut) weight = 1.;
            if (p > 0) weight *= Elec_ScaleFactor * muScaleFactor.GetVal() * Lumi_ScaleFactor * trgDataMcScaleFactor.GetVal();

            int jetsPt20 = 0;
            int jetsPt30 = 0;
            for (int jetIt = 0; jetIt < JetColl_size; ++ jetIt) {
               if (Jet_pt[jetIt] > 20.) ++jetsPt20;
               if (Jet_pt[jetIt] > 30.) ++jetsPt30;
            }

            // fill the data tree
            emuInvMass = invMass;
            if (p > 0) trueMass = genelemom_mass[0];
            eCharge = gsf_charge[GSF_passHEEP[0]];
            muCharge = muon_charge[MU_passGOOD[0]];
            totWeight = weight;
            if (p > 0) totWeight *= input[p].second * LumiFactor;
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.442) evtRegion = 0;
            else if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(gsfsc_eta[GSF_passHEEP[0]]) < 2.5) evtRegion = 1;
            else evtRegion = -1;
            // fill control variables
            nVtx = pvsize;
            eDzMinusMuDz = fabs(gsf_dz_firstPVtx[GSF_passHEEP[0]] - muon_dz_firstPVtx[MU_passGOOD[0]]);
            eDxyMinusMuDxy = fabs(gsf_dxy_firstPVtx[GSF_passHEEP[0]] - muon_dxy_firstPVtx[MU_passGOOD[0]]);
            if (fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]) < 3.14) dPhi = fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]);
            if (fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]) > 3.14) dPhi = 6.28 - fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]);
            dEta = fabs(muon_eta[MU_passGOOD[0]] - gsf_eta[GSF_passHEEP[0]]);
            eleEt = gsf_gsfet[GSF_passHEEP[0]];
            eleEta = gsf_eta[GSF_passHEEP[0]];
            elePhi = gsf_phi[GSF_passHEEP[0]];
            eleDEta = gsf_deltaeta[GSF_passHEEP[0]];
            eleDPhi = gsf_deltaphi[GSF_passHEEP[0]];
            eleHOE = gsf_hovere[GSF_passHEEP[0]];
            eleE1x5overE5x5 = gsf_e1x5overe5x5[GSF_passHEEP[0]];
            eleE2x5overE5x5 = gsf_e2x5overe5x5[GSF_passHEEP[0]];
            eleSigmaIEIE = gsf_sigmaIetaIeta[GSF_passHEEP[0]];
            eleEcalIso = gsf_ecaliso[GSF_passHEEP[0]];
            eleHcalIso1 = gsf_hcaliso1[GSF_passHEEP[0]];
            eleHcalIso2 = gsf_hcaliso2[GSF_passHEEP[0]];
            if (evtRegion == 0) eleHeepIso = -1. * (gsf_ecaliso[GSF_passHEEP[0]] + gsf_hcaliso1[GSF_passHEEP[0]] - 2. - rho * 0.28 - 0.03 * gsf_gsfet[GSF_passHEEP[0]]);
            else if (evtRegion == 1 && gsf_gsfet[GSF_passHEEP[0]] < 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passHEEP[0]] + gsf_hcaliso1[GSF_passHEEP[0]] - 2.5 - rho * 0.28);
            else if (evtRegion == 1 && gsf_gsfet[GSF_passHEEP[0]] >= 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passHEEP[0]] + gsf_hcaliso1[GSF_passHEEP[0]] - 2.5 - rho * 0.28 - 0.03 * (gsf_gsfet[GSF_passHEEP[0]] - 50.));
            eleTrkIso = gsf_trackiso[GSF_passHEEP[0]];
            eleLostHits = gsf_nLostInnerHits[GSF_passHEEP[0]];
            eleDXYFstPVtx = gsf_dxy_firstPVtx[GSF_passHEEP[0]];
            eleDZFstPVtx = gsf_dz_firstPVtx[GSF_passHEEP[0]];
            etEleOPtMu = gsf_gsfet[GSF_passHEEP[0]]/muon_pt[MU_passGOOD[0]];
            if (gsf_charge[GSF_passHEEP[0]] > 0 && muon_charge[MU_passGOOD[0]] < 0) lepPtPlusOPtMinus = gsf_gsfet[GSF_passHEEP[0]]/muon_pt[MU_passGOOD[0]];
            else if (gsf_charge[GSF_passHEEP[0]] < 0 && muon_charge[MU_passGOOD[0]] > 0) lepPtPlusOPtMinus = muon_pt[MU_passGOOD[0]]/gsf_gsfet[GSF_passHEEP[0]];
            else lepPtPlusOPtMinus = 0.;
            eChTimesMuCh = eCharge * muCharge;
            muPt = muon_pt[MU_passGOOD[0]];
            muPtErr = muon_ptError[MU_passGOOD[0]];
            muEta = muon_eta[MU_passGOOD[0]];
            muPhi = muon_phi[MU_passGOOD[0]];
            muHitLayers = muon_nlayerswithhits[MU_passGOOD[0]];
            muTrkHits = muon_nhitstrack[MU_passGOOD[0]];
            muPxlHits = muon_nhitspixel[MU_passGOOD[0]];
            muMuHits = muon_nhitsmuons[MU_passGOOD[0]];
            muDZFstPVtx = muon_dz_firstPVtx[MU_passGOOD[0]];
            muDXYFstPVtx = muon_dxy_firstPVtx[MU_passGOOD[0]];
            muNSeg = muon_nSegmentMatch[MU_passGOOD[0]];
            muTrkIso03 = muon_trackIso03[MU_passGOOD[0]];
            muIsoCombRel = CombRelIso;
            numOfJets = JetColl_size;
            numOfJetsPt20 = jetsPt20;
            numOfJetsPt30 = jetsPt30;

            emuTree->Fill();

            // fill the good emu events tree
            evRegion = evtRegion;
            emuMass = invMass;
            if (passTrg) eleDataTree->Fill();
            

            //same sign and opposite sign
            if (p == 0) {
               if (gsf_charge[GSF_passHEEP[0]] > 0 && muon_charge[MU_passGOOD[0]] > 0) ++nb_plus_plus;
               if (gsf_charge[GSF_passHEEP[0]] > 0 && muon_charge[MU_passGOOD[0]] < 0) ++nb_plus_minus;
               if (gsf_charge[GSF_passHEEP[0]] < 0 && muon_charge[MU_passGOOD[0]] > 0) ++nb_minus_plus;
               if (gsf_charge[GSF_passHEEP[0]] < 0 && muon_charge[MU_passGOOD[0]] < 0) ++nb_minus_minus;

               if (invMass > 600. || invMass < 1.) {
                  if (invMass > 600.) cout << "M_emu > 600 GeV/c^2 event | " << runnumber << ":" << luminosityBlock << ":" << eventnumber << " | "; 
                  if (invMass < 1.) cout << "M_emu < 1 GeV/c^2 event   | " << runnumber << ":" << luminosityBlock << ":" << eventnumber << " | "; 
                       cout << setw(8) << invMass << " | " << 
                       setw(7) << muon_pt[MU_passGOOD[0]] << " | " <<
                       setw(8) << muon_eta[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_phi[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_trackIso03[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_emIso03[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_hadIso03[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_normChi2[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_dxy_beamSpot[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_nhitstrack[MU_passGOOD[0]] << " | " <<
                       //setw(8) << muon_nhitsmuons[MU_passGOOD[0]] << " | " <<
                       setw(7) << gsf_gsfet[GSF_passHEEP[0]] << " | " << 
                       setw(8) << gsfsc_eta[GSF_passHEEP[0]] << " | " <<
                       //setw(8) << gsfsc_phi[GSF_passHEEP[0]] << " | " <<
                       //setw(8) << gsf_trackiso[GSF_passHEEP[0]] << " | " <<
                       //setw(8) << gsf_ecaliso[GSF_passHEEP[0]] << " | " <<
                       //setw(8) << gsf_hcaliso1[GSF_passHEEP[0]] << " | " <<
                       //setw(8) << gsf_hcaliso2[GSF_passHEEP[0]] << " | " <<
                       endl;
               }
            }
         }
        ///////////////////////////////////////////////////////////////////////
      } //END LOOP OVER EVENTS
        //////////////////////////////////////////////////////////////////////

      cout << "Number of selected events: " << evCounter << endl;
      cout << "moreHeep/goodHeep: " << moreHeepCounter << " / " << goodHeepCounter << " = " << (float)moreHeepCounter/(float)goodHeepCounter << endl;
      cout << "moreMu/goodMu: " << moreMuCounter << " / " << goodMuCounter << " = " << (float)moreMuCounter/(float)goodMuCounter << endl;
      cout << "moreHeepMu: " << moreHeepMuCounter << endl;

      if (p == 0) {
         //write root file with good emu events
         goodEvFile->cd();
         eleDataTree->Write();

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
      output->cd();

      TParameter<float> *mcWeight = new TParameter<float>((const char *)suffix[p], input[p].second);
      mcWeights.Add(mcWeight);

      TParameter<float> *eleVetoRatio = new TParameter<float>((const char *)suffix[p], (float)moreHeepCounter/(float)goodHeepCounter);
      eleVetoRatios.Add(eleVetoRatio);

      emuTree->Write();
      frEmuTree->Write();
      puMc->Write();
      puMcNorm->Write();
      puWeights->Write();
     //////////////////////////////////////////////////////////////////////////
   } //END LOOP OVER FILES
     //////////////////////////////////////////////////////////////////////////
   mcWeights.Write("mcWeights", TObject::kSingleKey);
   eleVetoRatios.Write("eleVetoRatios", TObject::kSingleKey);
   puData->Write();
   puDataNorm->Write();

   //PRINT INTEGRAL

   cout << endl;
   cout << "-----------------------------------------------------------------------------------------------------------" << endl;
   cout << "HEEP - TIGHT MU        Lumi        = " << LumiFactor << "pb-1" << endl;
   cout << "                       e pT EB     > " << bar_et << "GeV/c" << endl;
   cout << "                       e pT EE     > " << end_et << "GeV/c" << endl;
   cout << "                       mu pT       > " << muon_pt_min << "GeV/c" << endl;
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

bool EmuSpectrum::PassHEEP(const int &n)
{
   // HEEPv4.1
   // barrel
   if (fabs(gsfsc_eta[n]) < 1.442
       && gsf_gsfet[n] > bar_et
       && gsf_isecaldriven[n]
       && gsf_hovere[n] < bar_hoE
       && fabs(gsf_deltaeta[n]) < bar_DEta
       && fabs(gsf_deltaphi[n]) < bar_DPhi
       && ((gsf_e2x5overe5x5[n] > bar_e2x5e5x5) || (gsf_e1x5overe5x5[n] > bar_e1x5e5x5))
       && (gsf_ecaliso[n] + gsf_hcaliso1[n]) < (bar_isoEcalHcal1_1 + bar_isoEcalHcal1_2 * gsf_gsfet[n] + bar_isoEcalHcalRho * rho)
       && gsf_trackiso[n] < bar_isoTrack
       && gsf_nLostInnerHits[n] <= bar_missInnerHits
       && fabs(gsf_dxy_firstPVtx[n]) <= bar_dxy
      ) return true;

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
       && fabs(gsf_dxy_firstPVtx[n]) <= end_dxy
      ) return true;
   return false;
}

bool EmuSpectrum::PassHighPtMu(const int &n)
{
   if (muon_pt[n] > muon_pt_min
       && muon_ptError[n] / muon_pt[n] < muon_dptOverPt
       && fabs(muon_eta[n]) < muon_etaMax
       && muon_nhitstrack[n] >= muon_nHitsMinGlobal
       && muon_nhitspixel[n] >= muon_nHitsMinPixel
       && muon_nhitsmuons[n] >= muon_nHitsMinMuon
       && muon_nlayerswithhits[n] >= muon_nLayersMin
       && fabs(muon_dxy_firstPVtx[n]) < muon_impactParamMaxXY
       //&& fabs(muon_dz_firstPVtx[n]) < muon_impactParamMaxZ
       && muon_nSegmentMatch[n] >= muon_nSegMatchMin
       && muon_isTrackerMuon[n]
       && muon_trackIso03[n] / muon_pt[n] < muon_relIsoCutMax
      ) return true;

   return false;
}

bool
EmuSpectrum::PassFRPreSel(const int &n)
{
  if (fabs(gsfsc_eta[n]) < 1.442) {
    if (gsf_hovere[n] > frps_hOverECutEB) return false;
    if (gsf_nLostInnerHits[n] > frps_missHitsCutEB) return false;
    if (gsf_sigmaIetaIeta[n] > frps_sieieCutEB) return false;
    if (fabs(gsf_dxy_firstPVtx[n]) > frps_dxyCutEB) return false;
  }
  else if (fabs(gsfsc_eta[n]) > 1.56) {
    if (gsf_hovere[n] > frps_hOverECutEE) return false;
    if (gsf_nLostInnerHits[n] > frps_missHitsCutEE) return false;
    if (gsf_sigmaIetaIeta[n] > frps_sieieCutEE) return false;
    if (fabs(gsf_dxy_firstPVtx[n]) > frps_dxyCutEE) return false;
  }
  
  return true;
}

double
EmuSpectrum::FakeRate (const float& et, const float& eta)
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

int EmuSpectrum::Trigger(int &prescale, unsigned int *trig, const int &selector)
{
   // switch off trigger
   //if (selector < 2) trig[0]++;
   //return 1;

   if (selector == 1) {
      prescale = prescale_HLT_Mu22_Photon22_CaloIdL;
      if (HLT_Mu22_Photon22_CaloIdL >= 0) trig[0]++;
      return HLT_Mu22_Photon22_CaloIdL;
   } else if (selector == 2) {
      prescale = prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      if (HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL >= 0) trig[1]++;
      return HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   } else if (selector == 3) {
      prescale = prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      if (HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL >= 0) trig[2]++;
      return HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
   } else {
      // select an unprescaled trigger
      if (prescale_HLT_Mu22_Photon22_CaloIdL == 1) {
         prescale = prescale_HLT_Mu22_Photon22_CaloIdL;
         trig[0]++;
         if (runnumber < runs_HLT_Mu22_Photon22_CaloIdL.first) runs_HLT_Mu22_Photon22_CaloIdL.first = runnumber;
         if (runnumber > runs_HLT_Mu22_Photon22_CaloIdL.second) runs_HLT_Mu22_Photon22_CaloIdL.second = runnumber;
         return HLT_Mu22_Photon22_CaloIdL;
      } else if (prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 1) {
         prescale = prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
         trig[1]++;
         if (runnumber < runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first) runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first = runnumber;
         if (runnumber > runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second) runs_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second = runnumber;
         return HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      } else if (prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 1) {
         prescale = prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
         trig[2]++;
         if (runnumber < runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first) runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.first = runnumber;
         if (runnumber > runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second) runs_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL.second = runnumber;
         return HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      }
   }

   if (prescale_HLT_Mu22_Photon22_CaloIdL == 0 && prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 0 && prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL == 0) {
      //cout << "Prescale alert! All selecting triggers have prescale 0. Event: " << runnumber << ":" << luminosityBlock << ":" << eventnumber << " Triggers: " << HLT_Mu22_Photon22_CaloIdL << HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL << HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL << endl;
   } else {
      //cout << "Prescale alert! No unprescaled trigger found. Event: " << runnumber << ":" << luminosityBlock << ":" << eventnumber << endl;
   }
   prescale = 0;
   return -1;
}

