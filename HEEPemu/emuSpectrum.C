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
   ////-----------------------------------------------------------------------
   //// MuEG dataset triggered by Mu22_photon22 trigger
   //float LumiFactor = 19703.; //Lumi in pb-1
   //// DATA file
   //TString dataFile = "file:////user/treis/data2013/MuEG_Run2012A+B+C+D-ReReco22Jan2013_1e1muSkim_19703pb-1.root";
   //// pile up histogram
   //TString puFile = "file:////user/treis/data2013/pileup/pileupTrue_MuEG_Run2012ABCDReReco22Jan2013.root";
   //string outfileName = "emuSpec_MuGammaTrg";
   //float eleL1Eff = 0.99;
   //int trgSelector = 0;
   ////-----------------------------------------------------------------------
   
   //-----------------------------------------------------------------------
   // SingleMu dataset triggered by Mu40_eta2p1 trigger
   float LumiFactor = 19706.; //Lumi in pb-1
   // DATA file
   TString dataFile = "file:////user/treis/data2013/SingleMu_Run2012ABCD-22Jan2013-v1_AOD_ele20mu20_19706pb-1.root";
   // pile up histogram
   TString puFile = "file:////user/treis/data2013/pileup/pileupTrue_SingleMu_Run2012ABCD-22Jan2013.root";
   //string outfileName = "emuSpec_singleMuTrg_altdiboson";
   string outfileName = "emuSpec_test";
   muon_pt_min = 45.;
   muon_etaMax = 2.1;
   float eleL1Eff = 1.; // there is no electron L1 for the trigger used and so just set it to 1.
   int trgSelector = 4;
   mResV7C1->SetParameters(1.358e-2, 5.474e-1, 6.146e+3);
   mResV7C2->SetParameters(1.224e-2, 6.18e-2, 5.229e+2);
   //-----------------------------------------------------------------------

   TParameter<float> lumi("lumi", LumiFactor);
   // scale factors
   // muon factors: mu high_pt id trk iso https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=257000  
   TParameter<float> trgL1Eff("trgL1Eff", eleL1Eff); // ele L1 eff
   // epsilon_cand from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgCommissioningAndPhysicsDeliverables#Electron_reconstruction_effi_AN1
   TParameter<float> eps_cand_sf_0p8("eps_cand_sf_0p8", 0.990); // data/MC scale for epsilon_cand (>50GeV) |eta|<0.8
   TParameter<float> eps_cand_sf_err_0p8("eps_cand_sf_err_0p8", 0.004); // data/MC scale error (stat. + syst.) for epsilon_cand (>50GeV) |eta|<0.8
   TParameter<float> eps_cand_sf_0p8to1p4442("eps_cand_sf_0p8to1p4442", 0.991); // data/MC scale for epsilon_cand (>50GeV) 0.8<|eta|<1.4442
   TParameter<float> eps_cand_sf_err_0p8to1p4442("eps_cand_sf_err_0p8to1p4442", 0.004); // data/MC scale error (stat. + syst.) for epsilon_cand (>50GeV) 0.8<|eta|<1.4442
   TParameter<float> eps_cand_sf_1p566to2p0("eps_cand_sf_1p566to2p0", 0.990); // data/MC scale for epsilon_cand (>50GeV) 1.566<|eta|<2.0
   TParameter<float> eps_cand_sf_err_1p566to2p0("eps_cand_sf_err_1p566to2p0", 0.005); // data/MC scale error (stat. + syst.) for epsilon_cand (>50GeV) 1.566<|eta|<2.0
   TParameter<float> eps_cand_sf_2p0to2p5("eps_cand_sf_2p0to2p5", 0.998); // data/MC scale for epsilon_cand (>50GeV) 2.0<|eta|<2.5
   TParameter<float> eps_cand_sf_err_2p0to2p5("eps_cand_sf_err_2p0to2p5", 0.006); // data/MC scale error (stat. + syst.) for epsilon_cand (>50GeV) 2.0<|eta|<2.5
   // HEEP eff scale factors from AN-13-359 Table 5 reReco
   TParameter<float> eps_heep_sf_eb_pt35("eps_heep_sf_eb_pt35", 0.997); // HEEP eff scale factor
   TParameter<float> eps_heep_sf_err_eb_pt35("eps_heep_sf_err_eb_pt35", 0.007); // HEEP eff scale factor error
   TParameter<float> eps_heep_sf_eb_pt100("eps_heep_sf_eb_pt100", 0.985); // HEEP eff scale factor
   TParameter<float> eps_heep_sf_err_eb_pt100("eps_heep_sf_err_eb_pt100", 0.014); // HEEP eff scale factor error
   TParameter<float> eps_heep_sf_ee_pt35("eps_heep_sf_eb_pt35", 0.979); // HEEP eff scale factor
   TParameter<float> eps_heep_sf_err_ee_pt35("eps_heep_sf_err_eb_pt35", 0.006); // HEEP eff scale factor error
   TParameter<float> eps_heep_sf_ee_pt100("eps_heep_sf_eb_pt100", 0.981); // HEEP eff scale factor
   TParameter<float> eps_heep_sf_err_ee_pt100("eps_heep_sf_err_eb_pt100", 0.007); // HEEP eff scale factor error
   // muon scale factors from https://indico.cern.ch/getFile.py/access?contribId=1&resId=2&materialId=slides&confId=257630
   TParameter<float> lumiScaleFactorEB("lumiScaleFactorEB", 0.997);  // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1
   TParameter<float> lumiScaleFactorEE("lumiScaleFactorEE", 0.934);  // powheg - from normalization to the Z peak of the Z->ee spectrum HEEP v4.1

   // global systematic errors
   TParameter<float> systErrLumi("systErrLumi", 0.026);
   TParameter<float> systErrEff("systErrEff", 0.010); // muon err (0.0057) & ele err (0.0086)

   float muPtSmearFactor = 0.3; // smear 1/muon_pt with 30% of the resolution for shape uncertainties
   float jetLeptonVetoDR = 0.5;
   float bDiscrWP_M = 0.679; // b jet discriminator for working point M (1% mistag) of CSV tagger
   bool usePUInfo = true;
   bool lowMassPuOnly = false;
   float puMassCut = 120.;
   // selection cuts /////////////////////////////////////////////////////////
   float minInvMass = 0.;

   TH1::SetDefaultSumw2(kTRUE);

   ///////////////////////////////////////////////////////////////////////////
   // INPUT FILES
   vector<pair<TFile *, float> > input;
   vector<float> nGenEvtsV;
   vector<bool> storeGenMTtbar;
   vector<bool> storeEmuMass; // store emu_mass variable directly from the tree
   vector<bool> storeShapes; // store +/- 1 sigma trees to build shape histogram
   THashList systErrMCs;
   THashList mcWeights;
   THashList nGenEvents;
   THashList eleVetoRatios;

   // DATA
   TFile *inData = TFile::Open(dataFile);
   input.push_back(make_pair(inData, 1.)); //DATA
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   // MC
   TFile *inTTbar = TFile::Open("file:////user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_emuSkim_28150723ev.root");
   //input.push_back(make_pair(inTTbar, 225.197 / 28150723.)); // NLO
   //input.push_back(make_pair(inTTbar, 234. / 28150723.));  // approx NNLO
   nGenEvtsV.push_back(28150723.);
   input.push_back(make_pair(inTTbar, 245.8 / nGenEvtsV.back()));  // NNLO http://arxiv.org/pdf/1303.6254.pdf
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inTTbar700to1000 = TFile::Open("file:////user/treis/mcsamples/TT_Mtt-700to1000_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_3082812ev.root");
   nGenEvtsV.push_back(3082812.);
   input.push_back(make_pair(inTTbar700to1000, 15.614 / nGenEvtsV.back() * 245.8/211.));  // ttbar  mtt 700to1000
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar700to1000", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inTTbar1000up = TFile::Open("file:////user/treis/mcsamples/TT_Mtt-1000toInf_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_1249111ev.root");
   nGenEvtsV.push_back(1249111.);
   input.push_back(make_pair(inTTbar1000up, 2.954 / nGenEvtsV.back() * 245.8/211.));  // ttbar  mtt>1000
   systErrMCs.Add(new TParameter<float>("systErrMcTtbar1000up", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inTTbarPriv600up = TFile::Open("file:////user/treis/mcsamples/ttbar_tail_600_inf_35066ev.root");
   nGenEvtsV.push_back(35066.);
   input.push_back(make_pair(inTTbarPriv600up, 0.00431 * 1.166 / nGenEvtsV.back()));  // ttbar private production from Andreas with M^gen_emu=600+ scaled to Czakon et al. NNLO
   systErrMCs.Add(new TParameter<float>("systErrMcTtbarPriv600up", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inZtt = TFile::Open("file:////user/treis/mcsamples/DYToTauTau_M-20_CT10_TuneZ2star_v1+v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_52085800ev.root");
   nGenEvtsV.push_back(52085800.);
   input.push_back(make_pair(inZtt, 1915.1 / nGenEvtsV.back())); //Ztautau
   systErrMCs.Add(new TParameter<float>("systErrMcDyTauTau", 0.01));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inZjetsll = TFile::Open("file:////user/treis/mcsamples/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_30459503ev.root");
   nGenEvtsV.push_back(30459503.);
   input.push_back(make_pair(inZjetsll, 3531.9 / nGenEvtsV.back())); //ZJetsLL
   systErrMCs.Add(new TParameter<float>("systErrMcDyJetsLL", 0.01));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inWW = TFile::Open("file:////user/treis/mcsamples/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_10000431ev.root");
   nGenEvtsV.push_back(10000431.);
   input.push_back(make_pair(inWW, 54.838 / nGenEvtsV.back())); //WW
   systErrMCs.Add(new TParameter<float>("systErrMcWW", 0.04));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inWWpow = TFile::Open("file:////user/treis/mcsamples/WWJetTo2L2Nu_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_999864ev.root");
   nGenEvtsV.push_back(999864.);
   input.push_back(make_pair(inWWpow, 5.67 / nGenEvtsV.back())); //WW
   systErrMCs.Add(new TParameter<float>("systErrMcWWpow", 0.04));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inWWeMinusMuPlusPriv600up = TFile::Open("file:////user/treis/mcsamples/WW_emu_tail_Muplus_Eminus_25625ev.root");
   nGenEvtsV.push_back(25625.);
   input.push_back(make_pair(inWWeMinusMuPlusPriv600up, 1.531e-3 * 1.071 / nGenEvtsV.back())); //WW to e- mu+ from private production Andreas with SF to MCFM
   systErrMCs.Add(new TParameter<float>("systErrMcWWeMinusMuPlusPriv600up", 0.04));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inWWePlusMuMinusPriv600up = TFile::Open("file:////user/treis/mcsamples/WW_emu_tail_Eplus_Muminus_26537ev.root");
   nGenEvtsV.push_back(26537.);
   input.push_back(make_pair(inWWePlusMuMinusPriv600up, 1.586e-3 * 1.071 / nGenEvtsV.back())); //WW to e- mu+ from private production Andreas with SF to MCFM
   systErrMCs.Add(new TParameter<float>("systErrMcWWePlusMuMinusPriv600up", 0.04));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inWZ = TFile::Open("file:////user/treis/mcsamples/WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_10000283ev.root");
   nGenEvtsV.push_back(10000283.);
   input.push_back(make_pair(inWZ, 33.21 / nGenEvtsV.back())); //WZ
   systErrMCs.Add(new TParameter<float>("systErrMcWZ", 0.04));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inWZmg = TFile::Open("file:////user/treis/mcsamples/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_2017979ev.root");
   nGenEvtsV.push_back(2017979.);
   input.push_back(make_pair(inWZmg, 1.09 / nGenEvtsV.back())); //WZ
   systErrMCs.Add(new TParameter<float>("systErrMcWZmg", 0.04));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inTW = TFile::Open("file:////user/treis/mcsamples/T+Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_991118ev.root");
   nGenEvtsV.push_back(991118.);
   input.push_back(make_pair(inTW, 22.2 / nGenEvtsV.back())); //tW
   systErrMCs.Add(new TParameter<float>("systErrMcTW", 0.03));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inTWpow = TFile::Open("file:////user/treis/mcsamples/T+TBarToDilepton_tW-channel-DR_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_5947992ev.root");
   nGenEvtsV.push_back(5947992.);
   input.push_back(make_pair(inTWpow, 2.34 / nGenEvtsV.back())); //tW
   systErrMCs.Add(new TParameter<float>("systErrMcTWpow", 0.03));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inWJet = TFile::Open("file:////user/treis/mcsamples/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_emuSkim_76102995ev.root");
   nGenEvtsV.push_back(76102995.);
   input.push_back(make_pair(inWJet, 36257.2 / nGenEvtsV.back())); //W+jet
   systErrMCs.Add(new TParameter<float>("systErrMcWJets", 0.05));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inZmm = TFile::Open("file:////user/treis/mcsamples/DYToMuMu_M-20_CT10_TuneZ2star_v1+v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_52113126ev.root");
   nGenEvtsV.push_back(52113126.);
   input.push_back(make_pair(inZmm, 1915. / nGenEvtsV.back())); //Zmumu
   systErrMCs.Add(new TParameter<float>("systErrMcDyMuMu", 0.01));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inZee = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-20_CT10_TuneZ2star_v1+v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_46002499ev.root");
   nGenEvtsV.push_back(46002499.);
   input.push_back(make_pair(inZee, 1915. / nGenEvtsV.back())); //Zee
   systErrMCs.Add(new TParameter<float>("systErrMcDyEE", 0.01));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inZZ = TFile::Open("file:////user/treis/mcsamples/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_9799908ev.root");
   nGenEvtsV.push_back(9799908.);
   input.push_back(make_pair(inZZ, 17.654 / nGenEvtsV.back())); //ZZ
   systErrMCs.Add(new TParameter<float>("systErrMcZZ", 0.03));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inZZmg = TFile::Open("file:////user/treis/mcsamples/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_4807893ev.root");
   nGenEvtsV.push_back(4807893.);
   input.push_back(make_pair(inZZmg, 0.18 / nGenEvtsV.back())); //ZZ
   systErrMCs.Add(new TParameter<float>("systErrMcZZmg", 0.03));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inTTbarMg = TFile::Open("file:////user/treis/mcsamples/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_6923750ev.root");
   nGenEvtsV.push_back(6923750.);
   input.push_back(make_pair(inTTbarMg, 245.8/ nGenEvtsV.back())); //TTjets from MadGraph
   systErrMCs.Add(new TParameter<float>("systErrMcTtJets", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inTTbar22l = TFile::Open("file:////user/treis/mcsamples/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_emuSkim_12119013ev.root");
   nGenEvtsV.push_back(12119013.);
   input.push_back(make_pair(inTTbar22l, 13.43 / nGenEvtsV.back() * 245.8/(13.43+53.2+53.4))); //TT to 2l
   systErrMCs.Add(new TParameter<float>("systErrMcTtJets2l", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inTTbar21l = TFile::Open("file:////user/treis/mcsamples/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_11229902ev.root");
   nGenEvtsV.push_back(11229902.);
   input.push_back(make_pair(inTTbar21l, 53.2 / nGenEvtsV.back() * 245.8/(13.43+53.2+53.4))); //TT to 1l1jet
   systErrMCs.Add(new TParameter<float>("systErrMcTtJets1l1jet", 0.05));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(1);
   storeShapes.push_back(1);

   TFile *inTTW = TFile::Open("file:////user/treis/mcsamples/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_196046ev.root");
   nGenEvtsV.push_back(196046.);
   input.push_back(make_pair(inTTW, 0.232 / nGenEvtsV.back())); //TTW
   systErrMCs.Add(new TParameter<float>("systErrMcTtW", 0.29));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   TFile *inTTWW = TFile::Open("file:////user/treis/mcsamples/TTWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_emuSkim_217820ev.root");
   nGenEvtsV.push_back(217820.);
   input.push_back(make_pair(inTTWW, 0.002037 / nGenEvtsV.back())); //TTWW
   systErrMCs.Add(new TParameter<float>("systErrMcTtWW", 0.));
   storeGenMTtbar.push_back(1);
   storeEmuMass.push_back(0);
   storeShapes.push_back(1);

   //TFile *inWJet2 = TFile::Open("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/mcsamples2012/WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_35_22447541ev.root");
   //nGenEvtsV.push_back(22447541.);
   //input.push_back(make_pair(inWJet2, 1.91068E-5)); //W+jet
   //systErrMCs.Add(new TParameter<float>("systErrMcWJets70t0100", 0.05));
   //storeGenMTtbar.push_back(0);
   //storeEmuMass.push_back(0);
   //storeShapes.push_back(1);

   //TFile *inWJet3 = TFile::Open("file:////user/treis/mcsamples/WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1_AODSIM_gct1_33_12766180ev.root");
   //nGenEvtsV.push_back(12766180.);
   //input.push_back(make_pair(inWJet3, 1.79302E-5)); //W+jet
   //systErrMCs.Add(new TParameter<float>("systErrMcWJets100up", 0.05));
   //storeGenMTtbar.push_back(0);
   //storeEmuMass.push_back(0);
   //storeShapes.push_back(1);

   TFile *inSig0 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-250_TuneZ2star_8TeV_madgraph_v1_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig0, 0.000947 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig0", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig1 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-500_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9400ev.root");
   nGenEvtsV.push_back(9400.);
   input.push_back(make_pair(inSig1, 0.000239 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig1", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig2 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-750_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig2, 8.956e-5 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig2", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig3 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig3, 3.867e-5 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig3", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig4 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig4, 1.778e-5 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig4", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig5 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig5, 8.501e-6 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig5", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig6 = TFile::Open("file://///user/treis/mcsamples/ZprimeToEMu_M-1750_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig6, 4.124e-6 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig6", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig7 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig7, 2.014e-6 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig7", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig8 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig8, 4.724e-7 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig8", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig9 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig9, 1.059e-7 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig9", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig10 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9997ev.root");
   nGenEvtsV.push_back(9997.);
   input.push_back(make_pair(inSig10, 2.198e-8 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig10", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig11 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9993ev.root");
   nGenEvtsV.push_back(9993.);
   input.push_back(make_pair(inSig11, 4.157e-9 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig11", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig12 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9959ev.root");
   nGenEvtsV.push_back(9959.);
   input.push_back(make_pair(inSig12, 1.359e-10 / nGenEvtsV.back())); //Signal
   systErrMCs.Add(new TParameter<float>("systErrMcSig12", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig13 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-250_noAccCuts_TuneZ2star_8TeV_madgraph_v1_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig13, 0.000953 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig13", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig13a = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig13a, 0.000239 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig13a", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig14 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig14, 8.9745e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig14", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig15 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig15, 3.862e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig15", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig16 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig16, 1.781e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig16", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig17 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9997ev.root");
   nGenEvtsV.push_back(9997.);
   input.push_back(make_pair(inSig17, 8.503e-6 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig17", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig18 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9997ev.root");
   nGenEvtsV.push_back(9997.);
   input.push_back(make_pair(inSig18, 4.127e-6 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig18", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig19 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig19, 2.014e-6 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig19", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig20 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_noAccCuts_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig20, 4.735e-7 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig20", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig21 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig21, 1.060e-7 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig21", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig22 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_noAccCuts_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig22, 2.197e-8 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig22", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig23 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig23, 4.155e-9 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig23", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig24 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_USER_emuSkim_9966ev.root");
   nGenEvtsV.push_back(9966.);
   input.push_back(make_pair(inSig24, 1.293e-10 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig24", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig25 = TFile::Open("file:////user/treis/mcsamples/ZprimeLFVToEMu_M-1000_TuneZ2star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V19E-v1_AODSIM_emuSkim_9996ev.root");
   nGenEvtsV.push_back(9996.);
   input.push_back(make_pair(inSig25, 3.87e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts, central production
   systErrMCs.Add(new TParameter<float>("systErrMcSig25", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig26 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_Mz-1000_Ma-20000_TuneZ2star_8TeV_madgraph_v1_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig26, 1.018e-5 / nGenEvtsV.back())); //Signal with different masses for Z' and a'
   systErrMCs.Add(new TParameter<float>("systErrMcSig26", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig26a = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_Ma-1000_Mz-20000_TuneZ2star_8TeV_madgraph_v1_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig26a, 2.765e-05 / nGenEvtsV.back())); //Signal with different masses for Z' and a'
   systErrMCs.Add(new TParameter<float>("systErrMcSig26a", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig27 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-250_noAccCuts_TuneZ2star_8TeV_madgraph_v1_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig27, 0.000953 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig27", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig28 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig28, 0.000239 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig28", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig29 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig29, 8.9745e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig29", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig30 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig30, 3.862e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig30", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig31 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig31, 1.781e-5 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig31", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig32 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9997ev.root");
   nGenEvtsV.push_back(9997.);
   input.push_back(make_pair(inSig32, 8.503e-6 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig32", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig33 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-1750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9997ev.root");
   nGenEvtsV.push_back(9997.);
   input.push_back(make_pair(inSig33, 4.127e-6 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig33", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig34 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9999ev.root");
   nGenEvtsV.push_back(9999.);
   input.push_back(make_pair(inSig34, 2.014e-6 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig34", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig35 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_noAccCuts_TuneZ2star_8TeV_madgraph_v2_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig35, 4.735e-7 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig35", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig36 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_10000ev.root");
   nGenEvtsV.push_back(10000.);
   input.push_back(make_pair(inSig36, 1.060e-7 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig36", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig37 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_noAccCuts_TuneZ2star_8TeV_madgraph_v2_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig37, 2.197e-8 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig37", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig38 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9998ev.root");
   nGenEvtsV.push_back(9998.);
   input.push_back(make_pair(inSig38, 4.155e-9 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig38", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

   TFile *inSig39 = TFile::Open("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-EXOMCRECO_Private14_DR53X_PU_S10_START53_V7C2-v1_USER_emuSkim_9966ev.root");
   nGenEvtsV.push_back(9966.);
   input.push_back(make_pair(inSig39, 1.293e-10 / nGenEvtsV.back())); //Signal with no acceptance cuts
   systErrMCs.Add(new TParameter<float>("systErrMcSig39", 0.));
   storeGenMTtbar.push_back(0);
   storeEmuMass.push_back(1);
   storeShapes.push_back(0);

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
   lumiScaleFactorEB.Write();
   lumiScaleFactorEE.Write();
   eps_cand_sf_0p8.Write();
   eps_cand_sf_err_0p8.Write();
   eps_cand_sf_0p8to1p4442.Write();
   eps_cand_sf_err_0p8to1p4442.Write();
   eps_cand_sf_1p566to2p0.Write();
   eps_cand_sf_err_1p566to2p0.Write();
   eps_cand_sf_2p0to2p5.Write();
   eps_cand_sf_err_2p0to2p5.Write();
   eps_heep_sf_eb_pt35.Write();
   eps_heep_sf_err_eb_pt35.Write();
   eps_heep_sf_eb_pt100.Write();
   eps_heep_sf_err_eb_pt100.Write();
   eps_heep_sf_ee_pt35.Write();
   eps_heep_sf_err_ee_pt35.Write();
   eps_heep_sf_ee_pt100.Write();
   eps_heep_sf_err_ee_pt100.Write();
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

      
      ///////////////////////////////////////////////////////////////////////////
      //LOOP OVER SHAPE UNCERTAINTIES 
      ///////////////////////////////////////////////////////////////////////////
      for (unsigned int shUnc = 0; shUnc < shapeUncNames.size(); ++shUnc) {
         TString shapeUncName = shapeUncNames[shUnc];
         if (shUnc > 0) {
            if (p == 0) output->mkdir(shapeUncName);
            output->cd(shapeUncName);
            shapeUncName.Prepend("_");
            if (!storeShapes[p]) 
               continue;
         }
         else output->cd();

         // tree with event data
         bool passTrg = false;
         bool passAltTrg = false;
         bool passHeep = false;
         float emuInvMass = 0.;
         float trueMass = 0.;
         float genEmuMass = 0.;
         int evtRegion;
         int eCharge;
         int muCharge;
         float totWeight = 1.;
         float puWeight = 1.;
         float trgEff = 1.;
         float trgEffSf = 1.;
         float eleEffSf = 1.;
         float eleEffSfErr = 0.;
         float muEffSf = 1.;
         float topRewSf = 1.;
         TTree *emuTree = new TTree("emuTree_" + suffix[p] + shapeUncName, "emuTree_" + suffix[p] + shapeUncName);
         emuTree->Branch("runnr", &runnumber, "runnr/i");
         emuTree->Branch("eventnr", &eventnumber, "eventnr/i");
         emuTree->Branch("lumiSec", &luminosityBlock, "lumiSec/i");
         emuTree->Branch("passTrg", &passTrg, "passTrg/O");
         emuTree->Branch("passAltTrg", &passAltTrg, "passAltTrg/O");
         emuTree->Branch("mass", &emuInvMass, "mass/F");
         if (p > 0) {
            emuTree->Branch("trueMass", &trueMass, "trueMass/F");
            emuTree->Branch("genEmuMass", &genEmuMass, "genEmuMass/F");
         }
         emuTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
         emuTree->Branch("eCharge", &eCharge, "eCharge/I");
         emuTree->Branch("muCharge", &muCharge, "muCharge/I");
         emuTree->Branch("weight", &totWeight, "weight/F");
         emuTree->Branch("puWeight", &puWeight, "puWeight/F");
         emuTree->Branch("trgEff", &trgEff, "trgEff/F");
         emuTree->Branch("trgEffSf", &trgEffSf, "trgEffSf/F");
         emuTree->Branch("eleEffSf", &eleEffSf, "eleEffSf/F");
         emuTree->Branch("eleEffSfErr", &eleEffSfErr, "eleEffSfErr/F");
         emuTree->Branch("muEffSf", &muEffSf, "muEffSf/F");
         if (storeGenMTtbar[p]) {
           emuTree->Branch("topRewSf", &topRewSf, "topRewSf/F");
           emuTree->Branch("genMTtbar", &genPair_mass, "genMTtbar/F");
         }
         if (storeEmuMass[p]) {
           emuTree->Branch("emu_mass", &emu_mass, "emu_mass/F");
         }
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
         float numOfBJetsPt30 = 0.;
         float numOfBJetsMVAPt30 = 0.;
         float jetPt = 0.;
         float genElePt = 0.;
         float genEleEta = 0.;
         float genMuPt = 0.;
         float genMuEta = 0.;
         int genMuCharge = 0.;
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
         emuTree->Branch("numOfBJetsPt30", &numOfBJetsPt30, "numOfBJetsPt30/F");
         emuTree->Branch("numOfBJetsMVAPt30", &numOfBJetsMVAPt30, "numOfBJetsMVAPt30/F");
         emuTree->Branch("jetPt", &jetPt, "jetPt/F");
         emuTree->Branch("genElePt", &genElePt, "genElePt/F");
         emuTree->Branch("genEleEta", &genEleEta, "genEleEta/F");
         emuTree->Branch("genMuPt", &genMuPt, "genMuPt/F");
         emuTree->Branch("genMuEta", &genMuEta, "genMuEta/F");
         emuTree->Branch("genMuCharge", &genMuCharge, "genMuCharge/I");

         TTree *frEmuTree = new TTree("frEmuTree_" + suffix[p] + shapeUncName, "frEmuTree_" + suffix[p] + shapeUncName);
         float fakeRate = 0.;
         {
            frEmuTree->Branch("passTrg", &passTrg, "passTrg/O");
            frEmuTree->Branch("passAltTrg", &passAltTrg, "passAltTrg/O");
            frEmuTree->Branch("passHeep", &passHeep, "passHeep/O");
            frEmuTree->Branch("mass", &emuInvMass, "mass/F");
            frEmuTree->Branch("trueMass", &trueMass, "trueMass/F");
            frEmuTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
            frEmuTree->Branch("eCharge", &eCharge, "eCharge/I");
            frEmuTree->Branch("muCharge", &muCharge, "muCharge/I");
            frEmuTree->Branch("fakeRate", &fakeRate, "fakeRate/F");
            frEmuTree->Branch("puWeight", &puWeight, "puWeight/F");
            if (storeGenMTtbar[p]) {
              frEmuTree->Branch("genMTtbar", &genPair_mass, "genMTtbar/F");
              frEmuTree->Branch("topRewSf", &topRewSf, "topRewSf/F");
            }
            if (storeEmuMass[p]) {
              frEmuTree->Branch("emu_mass", &emu_mass, "emu_mass/F");
            }
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
            frEmuTree->Branch("numOfBJetsPt30", &numOfBJetsPt30, "numOfBJetsPt30/F");
            frEmuTree->Branch("numOfBJetsMVAPt30", &numOfBJetsMVAPt30, "numOfBJetsMVAPt30/F");
            frEmuTree->Branch("jetPt", &jetPt, "jetPt/F");

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
                    << "  triggered |"
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
            topRewSf = 1.;

            // modify eleScale, muScale, muonRes by +/-1 sigma for shape uncertainties
            switch (shUnc) {
               case 1:
                  ScaleEle(true);
                  break;
               case 2:
                  ScaleEle(false);
                  break;
               case 3:
                  ScaleMu(true);
                  break;
               case 4:
                  ScaleMu(false);
                  break;
               case 5:
                  ResMu(true);
                  break;
               case 6:
                  ResMu(false);
                  break;
               case 7:
                  SmearOneOverMuPt(muPtSmearFactor);
                  break;
            }
 
            // at least one gsf electron and one muon above the threshold
            if (gsf_size < 1 || muon_size < 1) continue;

            // get trigger info
            int prescale = 0;
            passTrg = true;
            if (p == 0 && Trigger(prescale, dataTrig, trgSelector) < 1) passTrg = false;
            if (p != 0) {
               if (trgSelector == 4) {    // single muon trigger
                  if (Trigger(prescale, trig, 4) < 1) passTrg = false;
               } else { 
                  if (jentry < nentries * dataTrig[0] / dataEntries) {
                     if (Trigger(prescale, trig, 1) < 1) passTrg = false;
                  } else if (jentry < nentries * (dataTrig[0] + dataTrig[1]) / dataEntries) {
                     if (Trigger(prescale, trig, 2) < 1) passTrg = false;
                  } else {
                     if (Trigger(prescale, trig, 3) < 1) passTrg = false;
                  }
               }
            }
            int altPrescale = 0;
            passAltTrg = true;
            unsigned int altTrig[3] = {0, 0, 0};
            if (trgSelector == 0) {
               if (Trigger(altPrescale, altTrig, 4) < 1) passAltTrg = false;
            } else if (trgSelector == 4) {
               if (Trigger(altPrescale, altTrig, 0) < 1) passAltTrg = false;
            } else {
               passAltTrg = false;
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
               //cleaning: fake electrons from muons skim electrons if there is a muon with pt>5 GeV within dR<0.1
               bool fakeEle = false;
               for (int k = 0; k < muon_size; ++k) {
                  if (muon_pt[k] < 5.) continue;
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
               //continue;
            }
            //if (GSF_passHEEP.size() > 1 || MU_passGOOD.size() > 1) cout << endl;

            // top pt reweighting factor
            if (storeGenMTtbar[p]) {
               float topPt = 0.;
               float atopPt = 0.;
               for (int j = 0; j < genPart_size; ++j) {
                  if (genPart_pdgid[j] == 6) topPt = genPart_pt[j];
                  else if (genPart_pdgid[j] == -6) atopPt = genPart_pt[j];
               }
               topRewSf = TopScaleFactor(topPt, atopPt);
            }

            //GSF-nonHEEP ele + GOOD muon
            if (GSF_passFrPre.size() > 0 && MU_passGOOD.size() > 0) {
               TLorentzVector ele1;
               TLorentzVector mu1;
               double invMass = 0.;

               // find the e-mu pair with the maximum invariant mass
               double maxInvMass = 0.;
               unsigned int eleInd = 0;
               unsigned int muInd = 0;
               for (unsigned int eleIt = 0; eleIt < GSF_passFrPre.size(); ++eleIt) {
                  for (unsigned int muIt = 0; muIt < MU_passGOOD.size(); ++muIt) {
                     ele1.SetPtEtaPhiM(gsf_gsfet[GSF_passFrPre[eleIt]], gsf_eta[GSF_passFrPre[eleIt]], gsf_phi[GSF_passFrPre[eleIt]], 0.000511);
                     mu1.SetPtEtaPhiM(muon_pt[MU_passGOOD[muIt]], muon_eta[MU_passGOOD[muIt]], muon_phi[MU_passGOOD[muIt]], 0.10566);

                     invMass = (ele1 + mu1).M();
                     if (invMass > maxInvMass) {
                        maxInvMass = invMass;
                        eleInd = eleIt;
                        muInd = muIt;
                     }
                  }
               }
               invMass = maxInvMass;

               //MASS CUT
               if (invMass > minInvMass) {
                  fakeRate = FakeRate(gsf_gsfet[GSF_passFrPre[eleInd]], gsfsc_eta[GSF_passFrPre[eleInd]]);
   
                  float CombRelIso = (muon_emIso03[MU_passGOOD[muInd]] + muon_hadIso03[MU_passGOOD[muInd]] + muon_trackIso03[MU_passGOOD[muInd]]) / muon_pt[MU_passGOOD[muInd]];
   
                  int jetsPt20 = 0;
                  int jetsPt30 = 0;
                  int bJetsPt30 = 0;
                  int bJetsMVAPt30 = 0;
                  for (int jetIt = 0; jetIt < JetColl_size; ++jetIt) {
                     // lepton veto
                     bool vetoed = false;
                     for (unsigned int eleIt = 0; eleIt < GSF_passFrPre.size(); ++eleIt) {
                        float dR = deltaR(Jet_eta[jetIt], Jet_phi[jetIt], gsf_eta[GSF_passFrPre[eleIt]], gsf_phi[GSF_passFrPre[eleIt]]);
                        if (dR < jetLeptonVetoDR) vetoed = true;
                     }
                     if (vetoed) continue;
                     for (unsigned int muIt = 0; muIt < MU_passGOOD.size(); ++muIt) {
                        float dR = deltaR(Jet_eta[jetIt], Jet_phi[jetIt], muon_eta[MU_passGOOD[muIt]], muon_phi[MU_passGOOD[muIt]]);
                        if (dR < jetLeptonVetoDR) vetoed = true;
                     }
                     if (vetoed) continue;
                     if (fabs(Jet_eta[jetIt]) < 2.4) {
                        if (Jet_pt[jetIt] > 20.) ++jetsPt20;
                        if (Jet_pt[jetIt] > 30.) {
                           ++jetsPt30;
                           if (cSecVertBTags[jetIt] > bDiscrWP_M) ++bJetsPt30;
                           if (cSecVertMVABTags[jetIt] > bDiscrWP_M) ++bJetsMVAPt30;
                        }
                     }
                  }
   
                  // fill the data tree
                  passHeep = PassHEEP(GSF_passFrPre[eleInd]);
                  emuInvMass = invMass;
                  if (p > 0) trueMass = genelemom_mass[0];
                  eCharge = gsf_charge[GSF_passFrPre[eleInd]];
                  muCharge = muon_charge[MU_passGOOD[muInd]];
                  totWeight = weight;
                  if (p > 0) totWeight *= input[p].second * LumiFactor;
                  if (fabs(gsfsc_eta[GSF_passFrPre[eleInd]]) < 1.442) evtRegion = 0;
                  else if (fabs(gsfsc_eta[GSF_passFrPre[eleInd]]) > 1.56 && fabs(gsfsc_eta[GSF_passFrPre[eleInd]]) < 2.5) evtRegion = 1;
                  else evtRegion = -1;
                  // fill control variables
                  nVtx = pvsize;
                  eDzMinusMuDz = fabs(gsf_dz_firstPVtx[GSF_passFrPre[eleInd]] - muon_dz_firstPVtx[MU_passGOOD[muInd]]);
                  eDxyMinusMuDxy = fabs(gsf_dxy_firstPVtx[GSF_passFrPre[eleInd]] - muon_dxy_firstPVtx[MU_passGOOD[muInd]]);
                  if (fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passFrPre[eleInd]]) < 3.14) dPhi = fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passFrPre[eleInd]]);
                  if (fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passFrPre[eleInd]]) > 3.14) dPhi = 6.28 - fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passFrPre[eleInd]]);
                  dEta = fabs(muon_eta[MU_passGOOD[muInd]] - gsf_eta[GSF_passFrPre[eleInd]]);
                  eleEt = gsf_gsfet[GSF_passFrPre[eleInd]];
                  eleEta = gsf_eta[GSF_passFrPre[eleInd]];
                  elePhi = gsf_phi[GSF_passFrPre[eleInd]];
                  eleDEta = gsf_deltaeta[GSF_passFrPre[eleInd]];
                  eleDPhi = gsf_deltaphi[GSF_passFrPre[eleInd]];
                  eleHOE = gsf_hovere[GSF_passFrPre[eleInd]];
                  eleE1x5overE5x5 = gsf_e1x5overe5x5[GSF_passFrPre[eleInd]];
                  eleE2x5overE5x5 = gsf_e2x5overe5x5[GSF_passFrPre[eleInd]];
                  eleSigmaIEIE = gsf_sigmaIetaIeta[GSF_passFrPre[eleInd]];
                  eleEcalIso = gsf_ecaliso[GSF_passFrPre[eleInd]];
                  eleHcalIso1 = gsf_hcaliso1[GSF_passFrPre[eleInd]];
                  eleHcalIso2 = gsf_hcaliso2[GSF_passFrPre[eleInd]];
                  if (evtRegion == 0) eleHeepIso = -1. * (gsf_ecaliso[GSF_passFrPre[eleInd]] + gsf_hcaliso1[GSF_passFrPre[eleInd]] - 2. - rho * 0.28 - 0.03 * gsf_gsfet[GSF_passFrPre[eleInd]]);
                  else if (evtRegion == 1 && gsf_gsfet[GSF_passFrPre[eleInd]] < 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passFrPre[eleInd]] + gsf_hcaliso1[GSF_passFrPre[eleInd]] - 2.5 - rho * 0.28);
                  else if (evtRegion == 1 && gsf_gsfet[GSF_passFrPre[eleInd]] >= 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passFrPre[eleInd]] + gsf_hcaliso1[GSF_passFrPre[eleInd]] - 2.5 - rho * 0.28 - 0.03 * (gsf_gsfet[GSF_passFrPre[eleInd]] - 50.));
                  eleTrkIso = gsf_trackiso[GSF_passFrPre[eleInd]];
                  eleLostHits = gsf_nLostInnerHits[GSF_passFrPre[eleInd]];
                  eleDXYFstPVtx = gsf_dxy_firstPVtx[GSF_passFrPre[eleInd]];
                  eleDZFstPVtx = gsf_dz_firstPVtx[GSF_passFrPre[eleInd]];
                  etEleOPtMu = gsf_gsfet[GSF_passFrPre[eleInd]]/muon_pt[MU_passGOOD[muInd]];
                  if (gsf_charge[GSF_passFrPre[eleInd]] > 0 && muon_charge[MU_passGOOD[muInd]] < 0) lepPtPlusOPtMinus = gsf_gsfet[GSF_passFrPre[eleInd]]/muon_pt[MU_passGOOD[muInd]];
                  else if (gsf_charge[GSF_passFrPre[eleInd]] < 0 && muon_charge[MU_passGOOD[muInd]] > 0) lepPtPlusOPtMinus = muon_pt[MU_passGOOD[muInd]]/gsf_gsfet[GSF_passFrPre[eleInd]];
                  else lepPtPlusOPtMinus = 0.;
                  eChTimesMuCh = eCharge * muCharge;
                  muPt = muon_pt[MU_passGOOD[muInd]];
                  muPtErr = muon_ptError[MU_passGOOD[muInd]];
                  muEta = muon_eta[MU_passGOOD[muInd]];
                  muPhi = muon_phi[MU_passGOOD[muInd]];
                  muHitLayers = muon_nlayerswithhits[MU_passGOOD[muInd]];
                  muTrkHits = muon_nhitstrack[MU_passGOOD[muInd]];
                  muPxlHits = muon_nhitspixel[MU_passGOOD[muInd]];
                  muMuHits = muon_nhitsmuons[MU_passGOOD[muInd]];
                  muDZFstPVtx = muon_dz_firstPVtx[MU_passGOOD[muInd]];
                  muDXYFstPVtx = muon_dxy_firstPVtx[MU_passGOOD[muInd]];
                  muNSeg = muon_nSegmentMatch[MU_passGOOD[muInd]];
                  muTrkIso03 = muon_trackIso03[MU_passGOOD[muInd]];
                  muIsoCombRel = CombRelIso;
                  numOfJets = JetColl_size;
                  numOfJetsPt20 = jetsPt20;
                  numOfJetsPt30 = jetsPt30;
                  numOfBJetsPt30 = bJetsPt30;
                  numOfBJetsMVAPt30 = bJetsMVAPt30;

                  frEmuTree->Fill();
               }
            }

            //HEEP ele + GOOD muon
            if (GSF_passHEEP.size() > 0 && MU_passGOOD.size() > 0) {
               TLorentzVector ele1;
               TLorentzVector mu1;
               double invMass = 0.;

               // find the e-mu pair with the maximum invariant mass
               double maxInvMass = 0.;
               unsigned int eleInd = 0;
               unsigned int muInd = 0;
               for (unsigned int eleIt = 0; eleIt < GSF_passHEEP.size(); ++eleIt) {
                  for (unsigned int muIt = 0; muIt < MU_passGOOD.size(); ++muIt) {
                     ele1.SetPtEtaPhiM(gsf_gsfet[GSF_passHEEP[eleIt]], gsf_eta[GSF_passHEEP[eleIt]], gsf_phi[GSF_passHEEP[eleIt]], 0.000511);
                     mu1.SetPtEtaPhiM(muon_pt[MU_passGOOD[muIt]], muon_eta[MU_passGOOD[muIt]], muon_phi[MU_passGOOD[muIt]], 0.10566);

                     invMass = (ele1 + mu1).M();
                     if (invMass > maxInvMass) {
                        maxInvMass = invMass;
                        eleInd = eleIt;
                        muInd = muIt;
                     }
                  }
               }
               invMass = maxInvMass;

               // vertex matching of tracks by dz
               //if (fabs(gsf_dz_beamSpot[GSF_passHEEP[eleInd]] - muon_dz_beamSpot[MU_passGOOD[muInd]]) > 0.05) continue;

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

                  //continue;
               }

               ++evCounter;

               float CombRelIso = (muon_emIso03[MU_passGOOD[muInd]] + muon_hadIso03[MU_passGOOD[muInd]] + muon_trackIso03[MU_passGOOD[muInd]]) / muon_pt[MU_passGOOD[muInd]];

               // set correction factors according to detector region
               float heepEffSf = 1.;
               float heepEffSfErr = 0.;
               float Lumi_ScaleFactor = 1.;
               if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.442) {
                 Lumi_ScaleFactor = lumiScaleFactorEB.GetVal();
                 if (gsf_gsfet[GSF_passHEEP[eleInd]] > 35.) {
                   heepEffSf = eps_heep_sf_eb_pt35.GetVal();
                   heepEffSfErr = eps_heep_sf_err_eb_pt35.GetVal();
                 }
                 if (gsf_gsfet[GSF_passHEEP[eleInd]] > 100.) {
                   heepEffSf = eps_heep_sf_eb_pt100.GetVal();
                   heepEffSfErr = eps_heep_sf_err_eb_pt100.GetVal();
                 }
               }
               else if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.56 && fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 2.5) {
                 Lumi_ScaleFactor = lumiScaleFactorEE.GetVal();
                 if (gsf_gsfet[GSF_passHEEP[eleInd]] > 35.) {
                   heepEffSf = eps_heep_sf_ee_pt35.GetVal();
                   heepEffSfErr = eps_heep_sf_err_ee_pt35.GetVal();
                 }
                 if (gsf_gsfet[GSF_passHEEP[eleInd]] > 100.) {
                   heepEffSf = eps_heep_sf_ee_pt100.GetVal();
                   heepEffSfErr = eps_heep_sf_err_ee_pt100.GetVal();
                 }
               }
               if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 0.8) {
                 eleEffSf = heepEffSf * eps_cand_sf_0p8.GetVal();
                 eleEffSfErr = sqrt(heepEffSfErr*heepEffSfErr + eps_cand_sf_err_0p8.GetVal()*eps_cand_sf_err_0p8.GetVal());
               }
               else if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.442) {
                 eleEffSf = heepEffSf * eps_cand_sf_0p8to1p4442.GetVal();
                 eleEffSfErr = sqrt(heepEffSfErr*heepEffSfErr + eps_cand_sf_err_0p8to1p4442.GetVal()*eps_cand_sf_err_0p8to1p4442.GetVal());
               }
               else if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.56 && fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 2.0) {
                 eleEffSf = heepEffSf * eps_cand_sf_1p566to2p0.GetVal();
                 eleEffSfErr = sqrt(heepEffSfErr*heepEffSfErr + eps_cand_sf_err_1p566to2p0.GetVal()*eps_cand_sf_err_1p566to2p0.GetVal());
               }
               else if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 2.0 && fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 2.5) {
                 eleEffSf = heepEffSf * eps_cand_sf_2p0to2p5.GetVal();
                 eleEffSfErr = sqrt(heepEffSfErr*heepEffSfErr + eps_cand_sf_err_2p0to2p5.GetVal()*eps_cand_sf_err_2p0to2p5.GetVal());
               }
               // muon scale factors for id and trg plus trg efficiencies
               WeightMuonRecoIsoTrigger(muon_pt[MU_passGOOD[muInd]], muon_eta[MU_passGOOD[muInd]], muEffSf, trgEffSf, trgEff, suffix[p], trgL1Eff.GetVal());

               if (lowMassPuOnly && invMass > puMassCut) weight = 1.;
               if (p > 0) weight *= eleEffSf * muEffSf * Lumi_ScaleFactor * trgEffSf;

               int jetsPt20 = 0;
               int jetsPt30 = 0;
               int bJetsPt30 = 0;
               int bJetsMVAPt30 = 0;
               for (int jetIt = 0; jetIt < JetColl_size; ++ jetIt) {
                  // lepton veto
                  bool vetoed = false;
                  for (unsigned int eleIt = 0; eleIt < GSF_passHEEP.size(); ++eleIt) {
                     float dR = deltaR(Jet_eta[jetIt], Jet_phi[jetIt], gsf_eta[GSF_passHEEP[eleIt]], gsf_phi[GSF_passHEEP[eleIt]]);
                     if (dR < jetLeptonVetoDR) vetoed = true;
                  }
                  if (vetoed) continue;
                  for (unsigned int muIt = 0; muIt < MU_passGOOD.size(); ++muIt) {
                     float dR = deltaR(Jet_eta[jetIt], Jet_phi[jetIt], muon_eta[MU_passGOOD[muIt]], muon_phi[MU_passGOOD[muIt]]);
                     if (dR < jetLeptonVetoDR) vetoed = true;
                  }
                  if (vetoed) continue;
                  if (fabs(Jet_eta[jetIt]) < 2.4) {
                     if (Jet_pt[jetIt] > 20.) ++jetsPt20;
                     if (Jet_pt[jetIt] > 30.) {
                        ++jetsPt30;
                        if (cSecVertBTags[jetIt] > bDiscrWP_M) ++bJetsPt30;
                        if (cSecVertMVABTags[jetIt] > bDiscrWP_M) ++bJetsMVAPt30;
                     }
                  }
               }

               TLorentzVector genEle1;
               TLorentzVector genMu1;
               genEle1.SetPtEtaPhiM(genele_pt[0], genele_eta[0], genele_phi[0], 0.000511);
               genMu1.SetPtEtaPhiM(genmu_pt[0], genmu_eta[0], genmu_phi[0], 0.10566);
   
               // fill the data tree
               emuInvMass = invMass;
               if (p > 0) {
                 trueMass = genelemom_mass[0];
                 genEmuMass = (genEle1 + genMu1).M();
                 genElePt = genele_pt[0];
                 genEleEta = genele_eta[0];
                 genMuPt = genmu_pt[0];
                 genMuEta = genmu_eta[0];
                 genMuCharge = genmu_charge[0];
               }
               eCharge = gsf_charge[GSF_passHEEP[eleInd]];
               muCharge = muon_charge[MU_passGOOD[muInd]];
               totWeight = weight;
               if (p > 0) totWeight *= input[p].second * LumiFactor;
               if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.442) evtRegion = 0;
               else if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.56 && fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 2.5) evtRegion = 1;
               else evtRegion = -1;
               // fill control variables
               nVtx = pvsize;
               eDzMinusMuDz = fabs(gsf_dz_firstPVtx[GSF_passHEEP[eleInd]] - muon_dz_firstPVtx[MU_passGOOD[muInd]]);
               eDxyMinusMuDxy = fabs(gsf_dxy_firstPVtx[GSF_passHEEP[eleInd]] - muon_dxy_firstPVtx[MU_passGOOD[muInd]]);
               if (fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passHEEP[eleInd]]) < 3.14) dPhi = fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passHEEP[eleInd]]);
               if (fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passHEEP[eleInd]]) > 3.14) dPhi = 6.28 - fabs(muon_phi[MU_passGOOD[muInd]] - gsf_phi[GSF_passHEEP[eleInd]]);
               dEta = fabs(muon_eta[MU_passGOOD[muInd]] - gsf_eta[GSF_passHEEP[eleInd]]);
               eleEt = gsf_gsfet[GSF_passHEEP[eleInd]];
               eleEta = gsf_eta[GSF_passHEEP[eleInd]];
               elePhi = gsf_phi[GSF_passHEEP[eleInd]];
               eleDEta = gsf_deltaeta[GSF_passHEEP[eleInd]];
               eleDPhi = gsf_deltaphi[GSF_passHEEP[eleInd]];
               eleHOE = gsf_hovere[GSF_passHEEP[eleInd]];
               eleE1x5overE5x5 = gsf_e1x5overe5x5[GSF_passHEEP[eleInd]];
               eleE2x5overE5x5 = gsf_e2x5overe5x5[GSF_passHEEP[eleInd]];
               eleSigmaIEIE = gsf_sigmaIetaIeta[GSF_passHEEP[eleInd]];
               eleEcalIso = gsf_ecaliso[GSF_passHEEP[eleInd]];
               eleHcalIso1 = gsf_hcaliso1[GSF_passHEEP[eleInd]];
               eleHcalIso2 = gsf_hcaliso2[GSF_passHEEP[eleInd]];
               if (evtRegion == 0) eleHeepIso = -1. * (gsf_ecaliso[GSF_passHEEP[eleInd]] + gsf_hcaliso1[GSF_passHEEP[eleInd]] - 2. - rho * 0.28 - 0.03 * gsf_gsfet[GSF_passHEEP[eleInd]]);
               else if (evtRegion == 1 && gsf_gsfet[GSF_passHEEP[eleInd]] < 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passHEEP[eleInd]] + gsf_hcaliso1[GSF_passHEEP[eleInd]] - 2.5 - rho * 0.28);
               else if (evtRegion == 1 && gsf_gsfet[GSF_passHEEP[eleInd]] >= 50.) eleHeepIso = -1. * (gsf_ecaliso[GSF_passHEEP[eleInd]] + gsf_hcaliso1[GSF_passHEEP[eleInd]] - 2.5 - rho * 0.28 - 0.03 * (gsf_gsfet[GSF_passHEEP[eleInd]] - 50.));
               eleTrkIso = gsf_trackiso[GSF_passHEEP[eleInd]];
               eleLostHits = gsf_nLostInnerHits[GSF_passHEEP[eleInd]];
               eleDXYFstPVtx = gsf_dxy_firstPVtx[GSF_passHEEP[eleInd]];
               eleDZFstPVtx = gsf_dz_firstPVtx[GSF_passHEEP[eleInd]];
               etEleOPtMu = gsf_gsfet[GSF_passHEEP[eleInd]]/muon_pt[MU_passGOOD[muInd]];
               if (gsf_charge[GSF_passHEEP[eleInd]] > 0 && muon_charge[MU_passGOOD[muInd]] < 0) lepPtPlusOPtMinus = gsf_gsfet[GSF_passHEEP[eleInd]]/muon_pt[MU_passGOOD[muInd]];
               else if (gsf_charge[GSF_passHEEP[eleInd]] < 0 && muon_charge[MU_passGOOD[muInd]] > 0) lepPtPlusOPtMinus = muon_pt[MU_passGOOD[muInd]]/gsf_gsfet[GSF_passHEEP[eleInd]];
               else lepPtPlusOPtMinus = 0.;
               eChTimesMuCh = eCharge * muCharge;
               muPt = muon_pt[MU_passGOOD[muInd]];
               muPtErr = muon_ptError[MU_passGOOD[muInd]];
               muEta = muon_eta[MU_passGOOD[muInd]];
               muPhi = muon_phi[MU_passGOOD[muInd]];
               muHitLayers = muon_nlayerswithhits[MU_passGOOD[muInd]];
               muTrkHits = muon_nhitstrack[MU_passGOOD[muInd]];
               muPxlHits = muon_nhitspixel[MU_passGOOD[muInd]];
               muMuHits = muon_nhitsmuons[MU_passGOOD[muInd]];
               muDZFstPVtx = muon_dz_firstPVtx[MU_passGOOD[muInd]];
               muDXYFstPVtx = muon_dxy_firstPVtx[MU_passGOOD[muInd]];
               muNSeg = muon_nSegmentMatch[MU_passGOOD[muInd]];
               muTrkIso03 = muon_trackIso03[MU_passGOOD[muInd]];
               muIsoCombRel = CombRelIso;
               numOfJets = JetColl_size;
               numOfJetsPt20 = jetsPt20;
               numOfJetsPt30 = jetsPt30;
               numOfBJetsPt30 = bJetsPt30;
               numOfBJetsMVAPt30 = bJetsMVAPt30;

               emuTree->Fill();

               evRegion = evtRegion;
               emuMass = invMass;
               if (p == 0) {
                  // fill the good emu events tree
                  if (shUnc == 0 && passTrg) eleDataTree->Fill();
                  //same sign and opposite sign
                  if (gsf_charge[GSF_passHEEP[eleInd]] > 0 && muon_charge[MU_passGOOD[muInd]] > 0) ++nb_plus_plus;
                  if (gsf_charge[GSF_passHEEP[eleInd]] > 0 && muon_charge[MU_passGOOD[muInd]] < 0) ++nb_plus_minus;
                  if (gsf_charge[GSF_passHEEP[eleInd]] < 0 && muon_charge[MU_passGOOD[muInd]] > 0) ++nb_minus_plus;
                  if (gsf_charge[GSF_passHEEP[eleInd]] < 0 && muon_charge[MU_passGOOD[muInd]] < 0) ++nb_minus_minus;

                  if (invMass > 600. || invMass < 1.) {
                     if (invMass > 600.) cout << "M_emu > 600 GeV/c^2 event | " << runnumber << ":" << luminosityBlock << ":" << eventnumber << " | "; 
                     if (invMass < 1.) cout << "M_emu < 1 GeV/c^2 event   | " << runnumber << ":" << luminosityBlock << ":" << eventnumber << " | "; 
                          cout << setw(8) << invMass << " | " << 
                          setw(7) << muon_pt[MU_passGOOD[muInd]] << " | " <<
                          setw(8) << muon_eta[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_phi[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_trackIso03[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_emIso03[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_hadIso03[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_normChi2[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_dxy_beamSpot[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_nhitstrack[MU_passGOOD[muInd]] << " | " <<
                          //setw(8) << muon_nhitsmuons[MU_passGOOD[muInd]] << " | " <<
                          setw(7) << gsf_gsfet[GSF_passHEEP[eleInd]] << " | " << 
                          setw(8) << gsfsc_eta[GSF_passHEEP[eleInd]] << " | " <<
                          //setw(8) << gsfsc_phi[GSF_passHEEP[eleInd]] << " | " <<
                          //setw(8) << gsf_trackiso[GSF_passHEEP[eleInd]] << " | " <<
                          //setw(8) << gsf_ecaliso[GSF_passHEEP[eleInd]] << " | " <<
                          //setw(8) << gsf_hcaliso1[GSF_passHEEP[eleInd]] << " | " <<
                          //setw(8) << gsf_hcaliso2[GSF_passHEEP[eleInd]] << " | " <<
                          setw(4) << passTrg << " | " <<
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
            if (shUnc == 0) {
               goodEvFile->cd();
               eleDataTree->Write();
            }

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
         if (shUnc > 0) output->cd(shapeUncNames[shUnc]);
         else output->cd();

         TParameter<float> *nGenEvts = new TParameter<float>((const char *)(suffix[p] + shapeUncName), nGenEvtsV[p-1]);
         nGenEvents.Add(nGenEvts);
         TParameter<float> *mcWeight = new TParameter<float>((const char *)(suffix[p] + shapeUncName), input[p].second);
         mcWeights.Add(mcWeight);

         TParameter<float> *eleVetoRatio = new TParameter<float>((const char *)(suffix[p] + shapeUncName), (float)moreHeepCounter/(float)goodHeepCounter);
         eleVetoRatios.Add(eleVetoRatio);

         emuTree->Write();
         frEmuTree->Write();
         if (shUnc > 0) output->cd("..");
        //////////////////////////////////////////////////////////////////////////
      } //END LOOP OVER SHAPE UNCERTAINTIES
        //////////////////////////////////////////////////////////////////////////
      output->cd();
      puMc->Write();
      puMcNorm->Write();
      puWeights->Write();
     //////////////////////////////////////////////////////////////////////////
   } //END LOOP OVER FILES
     //////////////////////////////////////////////////////////////////////////
   nGenEvents.Write("nGenEvents", TObject::kSingleKey);
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
       && fabs(muon_dz_firstPVtx[n]) < muon_impactParamMaxZ
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
  // fake rate full 2012 rereco 19.7/fb from AN-13-359 table 10
  if (fabs(eta) < 1.5) {
    if (et >= 35. && et < 98.) return 0.0226 - 0.000153*et;
    else if (et >= 98. && et < 191.9) return 0.0115 - 3.98e-5*et;
    else if (et >= 191.9) return 0.00382;
    else return 0.;
  }
  else if (fabs(eta) > 1.5) {
    if (et >= 35. && et < 89.9) return 0.0823 - 0.000522*et + (fabs(eta)-1.9)*0.065;
    else if (et >= 89.9 && et < 166.4) return 0.0403 - 5.45e-5*et + (fabs(eta)-1.9)*0.065;
    else if (et >= 166.4) return 0.0290 + 1.32e-5*et + (fabs(eta)-1.9)*0.065;
    else return 0.;
  }
  else return 0.; 

  //// fake rate 19.3/fb
  //if (fabs(eta) < 1.5) {
  //  if (et < 189.3) return 0.0179 - 0.000056 * et;
  //  return 0.0073;
  //} else if (abs(eta) < 2.0) {
  //  if (et < 96.6) return exp(-2.31 - 0.011 * et);
  //  else if (et < 178.0) return 0.040 - 0.000059 * et;
  //  return 0.0295;
  //} else if (abs(eta) < 2.5) {
  //  if(et < 115.4) return 0.099 - 0.00035 * et;
  //  else return 0.0586;
  //} else return 0;
}

void
EmuSpectrum::ScaleEle (const bool up)
{
   int sign = -1;
   if (up) sign = 1;
   for (int i = 0; i < gsf_size; ++i) {
      if (fabs(gsfsc_eta[i]) < 1.442) gsf_gsfet[i] += sign*gsf_gsfet[i]*0.006;
      else if (fabs(gsfsc_eta[i]) > 1.56) gsf_gsfet[i] += sign*gsf_gsfet[i]*0.015;
   }
}

void
EmuSpectrum::ScaleMu (const bool up)
{
   int sign = -1;
   if (up) sign = 1;
   for (int i = 0; i < muon_size; ++i) {
      muon_pt[i] += sign*muon_pt[i]*0.00005*muon_pt[i]; // 5%/TeV
   }
}

void
EmuSpectrum::ResMu (const bool up)
{
   int sign = -1;
   if (up) sign = 1;
   // first find the reconstructed muons matching best to the generated ones
   for (int i = 0; i < genMu_size; ++i) {
      float bestDr = 0.5;
      int bestJ = -1;
      for (int j = 0; j < muon_size; ++j) {
         float dR = deltaR(genmu_eta[i], genmu_phi[i], muon_eta[j], muon_phi[j]);
         if (dR < bestDr) {
            bestDr = dR;
            bestJ = j;
         }
      }
      // if a match was found for the generated muon
      if (bestJ > -1) {
         // smear the generated muons pt according to the resolution +/- resulution uncertainty
         // take random pt from a Gaussian distribution centered at the generated muons pt 
         // with width from the resolution and assign it as the new muon_pt
         double v7c1Res = mResV7C1->Eval(genmu_pt[i]);
         muon_pt[bestJ] = randGen.Gaus(genmu_pt[i], genmu_pt[i]*(v7c1Res + sign*fabs(v7c1Res-mResV7C2->Eval(genmu_pt[i]))));
      }
   }
}

void
EmuSpectrum::SmearOneOverMuPt (float smearFactor)
{
   double v7c1Res = 0.;
   for (int i = 0; i < muon_size; ++i) {
      // smear 1/muon_pt with the smearFactor*resolution at the muon_pt
      // take random 1/pt from a Gaussian distribution centered at the 1/muon_pt 
      v7c1Res = mResV7C1->Eval(muon_pt[i]);
      muon_pt[i] = 1./randGen.Gaus(1./muon_pt[i], 1./muon_pt[i] * smearFactor*v7c1Res);
   }
}

double
EmuSpectrum::TopScaleFactor (const float& top_genpt, const float& atop_genpt)
{
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  // constants for 8 TeV dilepton
  float a = 0.159;
  float b = -0.00141;
  return sqrt(exp(a+b*top_genpt) * exp(a+b*atop_genpt)); 
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
   } else if (selector == 4) {
      prescale = prescale_HLT_Mu40_eta2p1;
      if (HLT_Mu40_eta2p1 >= 0) trig[0]++;
      return HLT_Mu40_eta2p1;
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

void 
EmuSpectrum::WeightMuonRecoIsoTrigger(float MuonPt, float MuonEta, float &weight_muon_reco, float &weight_trigger, float &eff_trigger, TString &type, float l1_eff)
{
   // muon scale factors from https://indico.cern.ch/getFile.py/access?contribId=1&resId=2&materialId=slides&confId=257630
   if(fabs(MuonEta)<0.9) {
      if (MuonPt>35. && MuonPt<40.) weight_muon_reco=0.994067;
      else if (MuonPt>40. && MuonPt<50.) weight_muon_reco=0.993037;
      else if (MuonPt>50. && MuonPt<60.) weight_muon_reco=0.991306;
      else if (MuonPt>60. && MuonPt<90.) weight_muon_reco=0.989358;
      else if (MuonPt>90. && MuonPt<140.) weight_muon_reco=1.0029;
      else if (MuonPt>140. && MuonPt<300.) weight_muon_reco=1.01755;
      else if (MuonPt>300.) weight_muon_reco=1.0;
   }
   if(fabs(MuonEta)>0.9 && fabs(MuonEta)<1.2) {
      if (MuonPt>35. && MuonPt<40.) weight_muon_reco=0.993689;
      else if (MuonPt>40. && MuonPt<50.) weight_muon_reco=0.993637;
      else if (MuonPt>50. && MuonPt<60.) weight_muon_reco=0.994911;
      else if (MuonPt>60. && MuonPt<90.) weight_muon_reco=0.990111;
      else if (MuonPt>90. && MuonPt<140.) weight_muon_reco=1.0091;
      else if (MuonPt>140. && MuonPt<300.) weight_muon_reco=1.00899;
      else if (MuonPt>300.) weight_muon_reco=1.0;
   }
   if(fabs(MuonEta)>1.2 && fabs(MuonEta)<2.1) {
      if (MuonPt>35. && MuonPt<40.) weight_muon_reco=0.996594;
      else if (MuonPt>40. && MuonPt<50.) weight_muon_reco=0.997472;
      else if (MuonPt>50. && MuonPt<60.) weight_muon_reco=0.996475;
      else if (MuonPt>60. && MuonPt<90.) weight_muon_reco=0.991681;
      else if (MuonPt>90. && MuonPt<140.) weight_muon_reco=1.01997;
      else if (MuonPt>140. && MuonPt<300.) weight_muon_reco=0.983776;
      else if (MuonPt>300.) weight_muon_reco=1.0;
   }
   if(fabs(MuonEta)>2.1 && fabs(MuonEta)<2.4) {
      if (MuonPt>35. && MuonPt<40.) weight_muon_reco=0.99407;
      else if (MuonPt>40. && MuonPt<50.) weight_muon_reco=0.996346;
      else if (MuonPt>50. && MuonPt<60.) weight_muon_reco=0.991808;
      else if (MuonPt>60. && MuonPt<90.) weight_muon_reco=0.986235;
      else if (MuonPt>90. && MuonPt<140.) weight_muon_reco=1.04176;
      else if (MuonPt>140. && MuonPt<300.) weight_muon_reco=0.789917;
      else if (MuonPt>300.) weight_muon_reco=1.0;
   }

   // muon factors: mu-high_pt-id-trk_iso https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=257000  
   if(type.BeginsWith("sig")) {
      if (MuonEta>-2.1 && MuonEta<-1.6) {weight_trigger=0.974366; eff_trigger=0.743592;}
      else if (MuonEta>-1.6 && MuonEta<-1.2) {weight_trigger=0.988066; eff_trigger=0.838491;}
      else if (MuonEta>-1.2 && MuonEta<-0.9) {weight_trigger=0.955999; eff_trigger=0.833293;}
      else if (MuonEta>-0.9 && MuonEta<-0.6) {weight_trigger=0.977317; eff_trigger=0.928455;}
      else if (MuonEta>-0.6 && MuonEta<-0.3) {weight_trigger=0.987245; eff_trigger=0.955694;}
      else if (MuonEta>-0.3 && MuonEta<-0.2) {weight_trigger=0.926888; eff_trigger=0.816532;}
      else if (MuonEta>-0.2 && MuonEta<0.2) {weight_trigger=0.982564; eff_trigger=0.945548;}
      else if (MuonEta>0.2 && MuonEta<0.3) {weight_trigger=0.945051; eff_trigger=0.826780;}
      else if (MuonEta>0.3 && MuonEta<0.6) {weight_trigger=0.981904; eff_trigger=0.952209;}
      else if (MuonEta>0.6 && MuonEta<0.9) {weight_trigger=0.98001; eff_trigger=0.931446;}
      else if (MuonEta>0.9 && MuonEta<1.2) {weight_trigger=0.955033; eff_trigger=0.829441;}
      else if (MuonEta>1.2 && MuonEta<1.6) {weight_trigger=0.967649; eff_trigger=0.807104;}
      else if (MuonEta>1.6 && MuonEta<2.1) {weight_trigger=1.00618; eff_trigger=0.814047;}
   } else {
      if(fabs(MuonEta)<0.9) {
         if (MuonPt<50.) {weight_trigger=0.977615; eff_trigger=0.929521;}
         else if (MuonPt>50. && MuonPt<60.) {weight_trigger=0.976165; eff_trigger=0.928760;}
         else if (MuonPt>60. && MuonPt<90.) {weight_trigger=0.974155; eff_trigger=0.924693;}
         else if (MuonPt>90. && MuonPt<140.) {weight_trigger=0.976909; eff_trigger=0.921874;}
         else if (MuonPt>140. && MuonPt<500.) {weight_trigger=0.990006; eff_trigger=0.929174;}
         else if (MuonPt>500.) {weight_trigger=0.990006; eff_trigger=0.929174;}
      }
      else if(fabs(MuonEta)>0.9 && fabs(MuonEta)<1.2) {
         if (MuonPt<50.) {weight_trigger=0.956918; eff_trigger=0.831272;}
         else if (MuonPt>50. && MuonPt<60.) {weight_trigger=0.954277; eff_trigger=0.832350;}
         else if (MuonPt>60. && MuonPt<90.) {weight_trigger=0.946978; eff_trigger=0.825169;}
         else if (MuonPt>90. && MuonPt<140.) {weight_trigger=0.952839; eff_trigger=0.825886;}
         else if (MuonPt>140. && MuonPt<500.) {weight_trigger=0.961069; eff_trigger=0.841223;}
         else if (MuonPt>500.) {weight_trigger=0.961069; eff_trigger=0.841223;}
      }
      else if(fabs(MuonEta)>1.2 && fabs(MuonEta)<2.1) {
         if (MuonPt<50.) {weight_trigger=0.987869; eff_trigger=0.803493;}
         else if (MuonPt>50. && MuonPt<60.) {weight_trigger=0.981194; eff_trigger=0.802700;}
         else if (MuonPt>60. && MuonPt<90.) {weight_trigger=0.972754; eff_trigger=0.796577;}
         else if (MuonPt>90. && MuonPt<140.) {weight_trigger=0.984463; eff_trigger=0.807423;}
         else if (MuonPt>140. && MuonPt<500.) {weight_trigger=0.986544; eff_trigger=0.783592;}
         else if (MuonPt>500.) {weight_trigger=0.986544; eff_trigger=0.783592;}
      }
   }
   if (fabs(MuonEta)>=2.1) {weight_trigger=1.0; eff_trigger=0.78;} // assumption

   eff_trigger *= l1_eff;
}

