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
   float LumiFactor = 19619.; //Lumi in pb-1   -LUMI FROM GOLDEN JSON
   TParameter<float> lumi("lumi", LumiFactor);

   // DATA file
   TString dataFile = "file:////user/treis/data2012/MuEG_Run2012A+B+C+D_13Jul2012+06Aug2012+24Aug2012+11Dec2012+PromptReco-Cv2+Dv1_Cert_190456-208686_gct1_46_19619pb-1.root";
   // pile up histogram
   TString puFile = "file:////user/treis/data2012/pileup/pileup_runA+B+C-ReReco+D-Prompt_puJSON-190389-208686_MuEG.root";

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
   TParameter<float> systErrLumi("systErrLumi", 0.022);
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

   // DATA
   TFile *inData = TFile::Open(dataFile);
   input.push_back(make_pair(inData, 1.)); //DATA
   storeGenMTtbar.push_back(0);

   // MC
   TFile *inTTbar = TFile::Open("file:////user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
   //input.push_back(make_pair(inTTbar, 225.197 / 28150723.)); // NLO
   input.push_back(make_pair(inTTbar, 234. / 28150723.));  // approx NNLO
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

   unsigned int dataTrig[3] = {0, 0, 0};
   unsigned int dataEntries = 0;

   // output file
   stringstream ssOutfile;
   ssOutfile << outfileName << " " << LumiFactor << "pb-1.root";
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
      float emuInvMass = 0.;
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
      emuTree->Branch("evtRegion", &evtRegion, "evtRegion/I");
      emuTree->Branch("eCharge", &eCharge, "eCharge/I");
      emuTree->Branch("muCharge", &muCharge, "muCharge/I");
      emuTree->Branch("weight", &totWeight, "weight/F");
      emuTree->Branch("puWeight", &puWeight, "puWeight/F");
      if (storeGenMTtbar[p]) emuTree->Branch("genMTtbar", &genPair_mass, "genMTtbar/F");
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
      emuTree->Branch("pfMet", &pfmet, "pfMet/F");
      emuTree->Branch("nVtx", &nVtx, "nVtx/F");
      emuTree->Branch("dZFstPVtx", &dZFstPVtx, "dZFstPVtx/F");
      emuTree->Branch("dXYFstPVtx", &dXYFstPVtx, "dXYFstPVtx/F");
      emuTree->Branch("rho", &rho, "rho/F");
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
      unsigned int trig[3] = {0, 0, 0};
      unsigned int evCounter = 0;
      /////////////////////////////////////////////////////////////////////////////////////////////
      //LOOP OVER EVENTS
      /////////////////////////////////////////////////////////////////////////////////////////////
      Long64_t nbytes = 0, nb = 0;
      //nentries = 10000;
      for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
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
         vector<int> MU_passGOOD;
         /////////////////////////////////////////////////////////////////////////////////////////////
         //loop over electrons
         for (int j = 0; j < gsf_size; ++j) {
            //cleaning: fake electrons from muons
            bool fakeEle = false;
            for (int k = 0; k < muon_size; ++k) {
               float DeltaR = sqrt((gsf_eta[j] - muon_eta[k]) * (gsf_eta[j] - muon_eta[k]) + (gsf_phi[j] - muon_phi[k]) * (gsf_phi[j] - muon_phi[k]));
               if (DeltaR < 0.1) {
                  fakeEle = true;
                  break;
               }
            }
            if (fakeEle) continue;

            if (PassHEEP(j)) GSF_passHEEP.push_back(j);
         }

         //loop over muons
         for (int j = 0; j < muon_size; ++j) {
            if (PassHighPtMu(j)) MU_passGOOD.push_back(j);
         }

         // veto when there are more than one good candidates
         if (GSF_passHEEP.size() > 1) continue;
         if (MU_passGOOD.size() > 1) continue;

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
            eCharge = gsf_charge[GSF_passHEEP[0]];
            muCharge = muon_charge[MU_passGOOD[0]];
            totWeight = weight;
            if (p > 0) totWeight *= input[p].second * LumiFactor;
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.442) evtRegion = 0;
            else if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.56 && fabs(gsfsc_eta[GSF_passHEEP[0]]) < 2.5) evtRegion = 1;
            else evtRegion = -1;
            // fill control variables
            nVtx = pvsize;
            dZFstPVtx = fabs(gsf_dz_firstPVtx[GSF_passHEEP[0]] - muon_dz_firstPVtx[MU_passGOOD[0]]);
            dXYFstPVtx = fabs(gsf_dxy_firstPVtx[GSF_passHEEP[0]] - muon_dxy_firstPVtx[MU_passGOOD[0]]);
            if (fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]) < 3.14) dPhi = fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]);
            if (fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]) > 3.14) dPhi = 6.28 - fabs(muon_phi[MU_passGOOD[0]] - gsf_phi[GSF_passHEEP[0]]);
            dEta = fabs(muon_eta[MU_passGOOD[0]] - gsf_eta[GSF_passHEEP[0]]);
            eleEt = gsf_gsfet[GSF_passHEEP[0]];
            eleEta = gsf_eta[GSF_passHEEP[0]];
            elePhi = gsf_phi[GSF_passHEEP[0]];
            eleDEta = gsf_deltaeta[GSF_passHEEP[0]];
            eleDPhi = gsf_deltaphi[GSF_passHEEP[0]];
            eleHOE = gsf_hovere[GSF_passHEEP[0]];
            eleSigmaIEIE = gsf_sigmaIetaIeta[GSF_passHEEP[0]];
            eleEcalIso = gsf_ecaliso[GSF_passHEEP[0]];
            eleHcalIso12 = gsf_hcaliso1[GSF_passHEEP[0]] + gsf_hcaliso2[GSF_passHEEP[0]];
            eleTrkIso = gsf_trackiso[GSF_passHEEP[0]];
            eleLostHits = gsf_nLostInnerHits[GSF_passHEEP[0]];
            muIsoCombRel = CombRelIso;
            muEtEleOPtMu = gsf_gsfet[GSF_passHEEP[0]]/muon_pt[MU_passGOOD[0]];
            if (gsf_charge[GSF_passHEEP[0]] > 0 && muon_charge[MU_passGOOD[0]] < 0) muPtPlusOPtMinus = gsf_gsfet[GSF_passHEEP[0]]/muon_pt[MU_passGOOD[0]];
            if (gsf_charge[GSF_passHEEP[0]] < 0 && muon_charge[MU_passGOOD[0]] > 0) muPtPlusOPtMinus = muon_pt[MU_passGOOD[0]]/gsf_gsfet[GSF_passHEEP[0]];
            muPt = muon_pt[MU_passGOOD[0]];
            muEta = muon_eta[MU_passGOOD[0]];
            muPhi = muon_phi[MU_passGOOD[0]];
            muHitLayers = muon_nlayerswithhits[MU_passGOOD[0]];
            muPxlHits = muon_nhitspixel[MU_passGOOD[0]];
            muMuHits = muon_nhitsmuons[MU_passGOOD[0]];
            muDZFstPVtx = muon_dz_firstPVtx[MU_passGOOD[0]];
            muDXYFstPVtx = muon_dxy_firstPVtx[MU_passGOOD[0]];
            muNSeg = muon_nSegmentMatch[MU_passGOOD[0]];
            muTrkIso03 = muon_trackIso03[MU_passGOOD[0]];
            numOfJets = JetColl_size;
            numOfJetsPt20 = jetsPt20;
            numOfJetsPt30 = jetsPt30;

            emuTree->Fill();

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

      if (p == 0) {
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
     //////////////////////////////////////////////////////////////////////////
   } //END LOOP OVER FILES
     //////////////////////////////////////////////////////////////////////////
   mcWeights.Write("mcWeights", TObject::kSingleKey);
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

