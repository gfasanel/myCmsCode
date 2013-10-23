#define accTimesEff_cxx
#include "accTimesEff.h"
#include <TH1.h>
#include <TF1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

void AccTimesEff::Loop()
{
   TStopwatch timer;
   timer.Start();
   // parameters /////////////////////////////////////////////////////////////
   vector<TString> files;
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeLFVToEMu_M-1000_TuneZ2star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V19E-v1_AODSIM_9996ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9998ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9998ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9997ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1750_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9997ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_noAccCuts_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_9998ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_noAccCuts_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_9998ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9998ev.root");
   files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_noAccCuts_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9966ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-500_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_9400ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-750_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1000_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_9999ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1250_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-1750_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-2000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-2500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_10000ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-3000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9999ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-3500_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9997ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-4000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9993ev.root");
   //files.push_back("file:////user/treis/mcsamples/ZprimeToEMu_M-5000_TuneZ2star_8TeV_madgraph_treis-Summer12_DR53X_PU_S10_START53_V7C1-v1_9959ev.root");
   string outfileName = "accTimesEffHistos";

   // output file formats
   const bool saveSpec = 0;
   const bool saveAsPdf = 0;
   const bool saveAsPng = 1;
   const bool saveAsRoot = 1;
   TString plotDir = "./plots/";

   unsigned int triggerInd = 1;  // 0: HLT_Mu22_Photon22_CaloIdL; 1: HLT_Mu40_eta2p1

   int font = 42; //62
   // selection cuts /////////////////////////////////////////////////////////
   float elePtCut = 35.;
   float muPtCut = 35.;
   float muEtaCut = 2.4;
   float minInvMass = 0.;

   TH1::SetDefaultSumw2(kTRUE);

   ///////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////
   bool triggerBit = 0;
   TString triggerName = "";
   if (triggerInd == 0) {
      triggerName = "HLT_Mu22_Photon22_CaloIdL";
   } else {
      triggerName = "HLT_Mu40_eta2p1";
      muon_pt_min = 45.;
      muPtCut = muon_pt_min;
      muon_etaMax = 2.1;
      muEtaCut = muon_etaMax;
   }

   TH1F *hGenEvts = new TH1F("hGenEvts", "hGenEvts", 145, 0., 5010.);
   hGenEvts->GetXaxis()->SetTitle("M_{Z'}^{truth}");
   hGenEvts->GetXaxis()->SetTitleFont(font);
   hGenEvts->GetXaxis()->SetLabelFont(font);
   hGenEvts->GetYaxis()->SetTitle("Events");
   hGenEvts->GetYaxis()->SetTitleFont(font);
   hGenEvts->GetYaxis()->SetLabelFont(font);
   hGenEvts->GetYaxis()->SetTitleOffset(1.2);
   hGenEvts->SetLineColor(kBlack);
   hGenEvts->SetLineWidth(2);
   hGenEvts->SetMarkerStyle(20);
   hGenEvts->SetMarkerColor(kBlack);
   TH1F *hGenEvts_ele_pt = (TH1F*)hGenEvts->Clone("hGenEvts_ele_pt");
   hGenEvts_ele_pt->SetBins(50, 0., 3000.);
   hGenEvts_ele_pt->GetXaxis()->SetTitle("e p_{T}^{truth}");
   TH1F *hGenEvts_ele_eta = (TH1F*)hGenEvts->Clone("hGenEvts_ele_eta");
   hGenEvts_ele_eta->SetBins(30, -3., 3.);
   hGenEvts_ele_eta->GetXaxis()->SetTitle("e #eta^{truth}");
   TH1F *hGenEvts_mu_pt = (TH1F*)hGenEvts_ele_pt->Clone("hGenEvts_mu_pt");
   hGenEvts_mu_pt->GetXaxis()->SetTitle("#mu p_{T}^{truth}");
   TH1F *hGenEvts_mu_eta = (TH1F*)hGenEvts_ele_eta->Clone("hGenEvts_mu_eta");
   hGenEvts_mu_eta->GetXaxis()->SetTitle("#mu #eta^{truth}");

   TH1F *hPassPtCuts = (TH1F*)hGenEvts->Clone("hPassPtCuts");
   TH1F *hPassEtaCuts = (TH1F*)hGenEvts->Clone("hPassEtaCuts");
   TH1F *hPassPtEtaCuts = (TH1F*)hGenEvts->Clone("hPassPtEtaCuts");
   TH1F *hPassPtEff = (TH1F*)hGenEvts->Clone("hPassPtEff");
   TH1F *hPassEtaEff = (TH1F*)hGenEvts->Clone("hPassEtaEff");
   TH1F *hPassPtEtaEff = (TH1F*)hGenEvts->Clone("hPassPtEtaEff");
   TH1F *hGenEvtsEleInAcc = (TH1F*)hGenEvts->Clone("hGenEvtsEleInAcc");
   TH1F *hGenEvtsEleInAccEB = (TH1F*)hGenEvts->Clone("hGenEvtsEleInAccEB");
   TH1F *hGenEvtsEleInAccEE = (TH1F*)hGenEvts->Clone("hGenEvtsEleInAccEE");
   TH1F *hGenEvtsMuInAcc = (TH1F*)hGenEvts->Clone("hGenEvtsMuInAcc");
   TH1F *hGenEvtsInAcc = (TH1F*)hGenEvts->Clone("hGenEvtsInAcc");
   TH1F *hTrgEvts = (TH1F*)hGenEvts->Clone("hTrgEvts");
   hTrgEvts->SetTitle("hTrgEvts");
   TH1F *hTrgEvts_mu22ph22 = (TH1F*)hGenEvts->Clone("hTrgEvts_mu22ph22");
   hTrgEvts_mu22ph22->SetTitle("hTrgEvts_mu22ph22");
   TH1F *hTrgEvts_mu22ph22_ele_pt = (TH1F*)hGenEvts_ele_pt->Clone("hTrgEvts_mu22ph22_ele_pt");
   hTrgEvts_mu22ph22_ele_pt->SetTitle("hTrgEvts_mu22ph22_ele_pt");
   TH1F *hTrgEvts_mu22ph22_ele_eta = (TH1F*)hGenEvts_ele_eta->Clone("hTrgEvts_mu22ph22_ele_eta");
   hTrgEvts_mu22ph22_ele_eta->SetTitle("hTrgEvts_mu22ph22_ele_eta");
   TH1F *hTrgEvts_mu22ph22_mu_pt = (TH1F*)hGenEvts_mu_pt->Clone("hTrgEvts_mu22ph22_mu_pt");
   hTrgEvts_mu22ph22_mu_pt->SetTitle("hTrgEvts_mu22ph22_mu_pt");
   TH1F *hTrgEvts_mu22ph22_mu_eta = (TH1F*)hGenEvts_mu_eta->Clone("hTrgEvts_mu22ph22_mu_eta");
   hTrgEvts_mu22ph22_mu_eta->SetTitle("hTrgEvts_mu22ph22_mu_eta");
   TH1F *hTrgEvts_mu40 = (TH1F*)hGenEvts->Clone("hTrgEvts_mu40");
   hTrgEvts_mu40->SetTitle("hTrgEvts_mu40");
   TH1F *hTrgEvts_mu40_ele_pt = (TH1F*)hGenEvts_ele_pt->Clone("hTrgEvts_mu40_ele_pt");
   hTrgEvts_mu40_ele_pt->SetTitle("hTrgEvts_mu40_ele_pt");
   TH1F *hTrgEvts_mu40_ele_eta = (TH1F*)hGenEvts_ele_eta->Clone("hTrgEvts_mu40_ele_eta");
   hTrgEvts_mu40_ele_eta->SetTitle("hTrgEvts_mu40_ele_eta");
   TH1F *hTrgEvts_mu40_mu_pt = (TH1F*)hGenEvts_mu_pt->Clone("hTrgEvts_mu40_mu_pt");
   hTrgEvts_mu40_mu_pt->SetTitle("hTrgEvts_mu40_mu_pt");
   TH1F *hTrgEvts_mu40_mu_eta = (TH1F*)hGenEvts_mu_eta->Clone("hTrgEvts_mu40_mu_eta");
   hTrgEvts_mu40_mu_eta->SetTitle("hTrgEvts_mu40_mu_eta");
   TH1F *hFltrMatchEvts_mu22ph22_phLeg = (TH1F*)hGenEvts->Clone("hFltrMatchEvts_mu22ph22_phLeg");
   hFltrMatchEvts_mu22ph22_phLeg->SetTitle("hFltrMatchEvts_mu22ph22_phLeg");
   TH1F *hFltrMatchEvts_mu22ph22_muLeg = (TH1F*)hGenEvts->Clone("hFltrMatchEvts_mu22ph22_muLeg");
   hFltrMatchEvts_mu22ph22_muLeg->SetTitle("hFltrMatchEvts_mu22ph22_muLeg");
   TH1F *hFltrMatchEvts_mu40 = (TH1F*)hGenEvts->Clone("hFltrMatchEvts_mu40");
   hFltrMatchEvts_mu40->SetTitle("hFltrMatchEvts_mu40");
   TH1F *hRecoEvts = (TH1F*)hGenEvts->Clone("hRecoEvts");
   hRecoEvts->SetTitle("hRecoEvts");
   TH1F *hRecoEvtsEB = (TH1F*)hRecoEvts->Clone("hRecoEvtsEB");
   TH1F *hRecoEvtsEE = (TH1F*)hRecoEvts->Clone("hRecoEvtsEE");
   TH1F *hRecoEleEvts = (TH1F*)hRecoEvts->Clone("hRecoEleEvts");
   TH1F *hRecoEleEvtsEB = (TH1F*)hRecoEvts->Clone("hRecoEleEvtsEB");
   TH1F *hRecoEleEvtsEE = (TH1F*)hRecoEvts->Clone("hRecoEleEvtsEE");
   TH1F *hRecoMuEvts = (TH1F*)hRecoEvts->Clone("hRecoMuEvts");
   TH1F *hRecoNoTrgEvts = (TH1F*)hGenEvts->Clone("hRecoNoTrgEvts");
   hRecoNoTrgEvts->SetTitle("hRecoNoTrgEvts");
   TH1F *hRecoNoTrgEvtsEB = (TH1F*)hRecoNoTrgEvts->Clone("hRecoNoTrgEvtsEB");
   TH1F *hRecoNoTrgEvtsEE = (TH1F*)hRecoNoTrgEvts->Clone("hRecoNoTrgEvtsEE");
   TH1F *hRecoNoTrgEleEvts = (TH1F*)hRecoNoTrgEvts->Clone("hRecoNoTrgEleEvts");
   TH1F *hRecoNoTrgEleEvtsEB = (TH1F*)hRecoNoTrgEvts->Clone("hRecoNoTrgEleEvtsEB");
   TH1F *hRecoNoTrgEleEvtsEE = (TH1F*)hRecoNoTrgEvts->Clone("hRecoNoTrgEleEvtsEE");
   TH1F *hRecoNoTrgMuEvts = (TH1F*)hRecoNoTrgEvts->Clone("hRecoNoTrgMuEvts");
   TH1F* hAcc;
   TH1F* hAccEle;
   TH1F* hAccEleEB;
   TH1F* hAccEleEE;
   TH1F* hAccMu;
   TH1F* hTrgEff;
   TH1F* hTrgEff_mu22ph22;
   TH1F* hTrgEff_mu22ph22_ele_pt;
   TH1F* hTrgEff_mu22ph22_ele_eta;
   TH1F* hTrgEff_mu22ph22_mu_pt;
   TH1F* hTrgEff_mu22ph22_mu_eta;
   TH1F* hTrgEff_mu40;
   TH1F* hTrgEff_mu40_ele_pt;
   TH1F* hTrgEff_mu40_ele_eta;
   TH1F* hTrgEff_mu40_mu_pt;
   TH1F* hTrgEff_mu40_mu_eta;
   TH1F* hAccTimesEff;
   TH1F* hAccTimesEffEB;
   TH1F* hAccTimesEffEE;
   TH1F* hAccTimesEffEle;
   TH1F* hAccTimesEffEleEB;
   TH1F* hAccTimesEffEleEE;
   TH1F* hAccTimesEffMu;
   TH1F* hAccTimesEffNoTrg;
   //TH1F* hAccTimesEffNoTrgEB;
   //TH1F* hAccTimesEffNoTrgEE;
   TH1F* hAccTimesEffNoTrgEle;
   TH1F* hAccTimesEffNoTrgEleEB;
   TH1F* hAccTimesEffNoTrgEleEE;
   TH1F* hAccTimesEffNoTrgMu;
   TH1F* hEffAftTrg;
   TH1F* hEffAftTrgEle;
   TH1F* hEffAftTrgEleEB;
   TH1F* hEffAftTrgEleEE;
   TH1F* hEffAftTrgMu;
   TH1F* hTrgRecoVsReco;
   //TH1F* hTrgRecoVsRecoEB;
   //TH1F* hTrgRecoVsRecoEE;
   TH1F* hTrgRecoVsRecoEle;
   TH1F* hTrgRecoVsRecoEleEB;
   TH1F* hTrgRecoVsRecoEleEE;
   TH1F* hTrgRecoVsRecoMu;
   TH1F* hFltrMatchVsReco_mu22ph22_phLeg;
   TH1F* hFltrMatchVsReco_mu22ph22_muLeg;
   TH1F* hFltrMatchVsReco_mu40;

   // output file
   stringstream ssOutfile;
   ssOutfile << outfileName << ".root";
   TFile *output = new TFile(ssOutfile.str().c_str(), "recreate");

   ///////////////////////////////////////////////////////////////////////////
   //LOOP OVER FILES 
   ///////////////////////////////////////////////////////////////////////////
   for (unsigned int p = 0; p < files.size(); ++p) {
      TFile* input = new TFile(files[p]);
      TTree *thetree = (TTree*)input->Get("gsfcheckerjob/tree");
      Init(thetree);
      Long64_t nentries = fChain->GetEntriesFast();
      cout << nentries << " events" << endl;

      unsigned int evCounter = 0;
      /////////////////////////////////////////////////////////////////////////////////////////////
      //LOOP OVER EVENTS
      /////////////////////////////////////////////////////////////////////////////////////////////
      //nentries = 10000;
      for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         fChain->GetEntry(jentry);
         // if (Cut(ientry) < 0) continue;
         //if (jentry % 50000 == 0) cout << "Processing event " << jentry << endl;
         thetree->GetEntry(jentry);

         //hard ele + hard muon
         TLorentzVector hardEle1;
         TLorentzVector hardMu1;
         hardEle1.SetPtEtaPhiM(hardGenEle_pt[0], hardGenEle_eta[0], hardGenEle_phi[0], 0.000511);
         hardMu1.SetPtEtaPhiM(hardGenMu_pt[0], hardGenMu_eta[0], hardGenMu_phi[0], 0.10566);
         double hardInvMass = (hardEle1 + hardMu1).M();
         //gen ele + gen muon
         TLorentzVector genEle1;
         TLorentzVector genMu1;
         genEle1.SetPtEtaPhiM(genele_pt[0], genele_eta[0], genele_phi[0], 0.000511);
         genMu1.SetPtEtaPhiM(genmu_pt[0], genmu_eta[0], genmu_phi[0], 0.10566);
         double genInvMass = (genEle1 + genMu1).M();

         //cout << genelemom_mass[0] << "     " << genInvMass << "     " << hardInvMass << endl;
         if (fabs(genelemom_pdgid[0]) < 7) continue;
         // fill the gen histograms
         hGenEvts->Fill(hardInvMass);
         hGenEvts_ele_pt->Fill(genele_pt[0]);
         hGenEvts_ele_eta->Fill(genele_eta[0]);
         hGenEvts_mu_pt->Fill(genmu_pt[0]);
         hGenEvts_mu_eta->Fill(genmu_eta[0]);

         // fill the acc histograms
         bool passEleEta = 0;
         bool passMuEta = 0;
         if (fabs(hardGenEle_eta[0]) < 1.442 || (fabs(hardGenEle_eta[0]) > 1.56 && fabs(hardGenEle_eta[0]) < 2.5)) passEleEta = 1;
         if (fabs(hardGenMu_eta[0]) < muEtaCut) passMuEta = 1;
         if (passEleEta && passMuEta) hPassEtaCuts->Fill(hardInvMass);
         if (hardGenEle_pt[0] > elePtCut) {
            if (fabs(hardGenEle_eta[0]) < 1.442) {
               hGenEvtsEleInAcc->Fill(hardInvMass);
               hGenEvtsEleInAccEB->Fill(hardInvMass);
            }
            else if (fabs(hardGenEle_eta[0]) > 1.56 && fabs(hardGenEle_eta[0]) < 2.5) {
               hGenEvtsEleInAcc->Fill(hardInvMass);
               hGenEvtsEleInAccEE->Fill(hardInvMass);
            }
            if (hardGenMu_pt[0] > muPtCut) {
               hPassPtCuts->Fill(hardInvMass);
               if (passEleEta && passMuEta) hPassPtEtaCuts->Fill(hardInvMass);
            }
         }
         if (hardGenMu_pt[0] > muPtCut && passMuEta) {
            hGenEvtsMuInAcc->Fill(hardInvMass);
            if (passEleEta) {
               if (hardGenEle_pt[0] > elePtCut) hGenEvtsInAcc->Fill(hardInvMass);
            }
         }

         // trigger?
         if (triggerInd == 0) triggerBit = HLT_Mu22_Photon22_CaloIdL;
         else triggerBit = HLT_Mu40_eta2p1;

         if (triggerBit) hTrgEvts->Fill(hardInvMass);
         if (HLT_Mu22_Photon22_CaloIdL) {
            hTrgEvts_mu22ph22->Fill(hardInvMass);
            hTrgEvts_mu22ph22_ele_pt->Fill(genele_pt[0]);
            hTrgEvts_mu22ph22_ele_eta->Fill(genele_eta[0]);
            hTrgEvts_mu22ph22_mu_pt->Fill(genmu_pt[0]);
            hTrgEvts_mu22ph22_mu_eta->Fill(genmu_eta[0]);
         }
         if (HLT_Mu40_eta2p1) {
            hTrgEvts_mu40->Fill(hardInvMass);
            hTrgEvts_mu40_ele_pt->Fill(genele_pt[0]);
            hTrgEvts_mu40_ele_eta->Fill(genele_eta[0]);
            hTrgEvts_mu40_mu_pt->Fill(genmu_pt[0]);
            hTrgEvts_mu40_mu_eta->Fill(genmu_eta[0]);
         }

         // at least one gsf electron and one muon above the threshold
         if (gsf_size < 1 || muon_size < 1) continue;

         vector<int> GSF_passHEEP;
         vector<int> MU_passGOOD;
         /////////////////////////////////////////////////////////////////////////////////////////////
         //loop over electrons
         for (int j = 0; j < gsf_size; ++j) {
            //cleaning: fake electrons from muons
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
         }

         //loop over muons
         for (int j = 0; j < muon_size; ++j) {
            if (PassHighPtMu(j)) MU_passGOOD.push_back(j);
         }

         if (GSF_passHEEP.size() == 1) {
            if (triggerBit) {
               hRecoEleEvts->Fill(hardInvMass);
               if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.5) hRecoEleEvtsEB->Fill(hardInvMass);
               if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.5) hRecoEleEvtsEE->Fill(hardInvMass);
            }
            hRecoNoTrgEleEvts->Fill(hardInvMass);
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.5) hRecoNoTrgEleEvtsEB->Fill(hardInvMass);
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.5) hRecoNoTrgEleEvtsEE->Fill(hardInvMass);
         } 
         if (MU_passGOOD.size() == 1) {
            if (triggerBit) hRecoMuEvts->Fill(hardInvMass);
            hRecoNoTrgMuEvts->Fill(hardInvMass);
         }

         // veto when there are more than one good candidates
         if (GSF_passHEEP.size() != 1 || MU_passGOOD.size() != 1) continue;

         //HEEP ele + GOOD muon
         TLorentzVector ele1;
         TLorentzVector mu1;

         ele1.SetPtEtaPhiM(gsf_gsfet[GSF_passHEEP[0]], gsf_eta[GSF_passHEEP[0]], gsf_phi[GSF_passHEEP[0]], 0.000511);
         mu1.SetPtEtaPhiM(muon_pt[MU_passGOOD[0]], muon_eta[MU_passGOOD[0]], muon_phi[MU_passGOOD[0]], 0.10566);

         double invMass = (ele1 + mu1).M();

         //MASS CUT
         if (invMass < minInvMass) continue;

         if (muMatch_hltL1Mu3p5EG12L3Filtered22[MU_passGOOD[0]]) hFltrMatchEvts_mu22ph22_muLeg->Fill(hardInvMass);
         if (gsfmatch_hltMu22Photon22CaloIdLHEFilter[GSF_passHEEP[0]]) hFltrMatchEvts_mu22ph22_phLeg->Fill(hardInvMass);
         if (muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q[MU_passGOOD[0]]) hFltrMatchEvts_mu40->Fill(hardInvMass);

         if (triggerBit) {
            hRecoEvts->Fill(hardInvMass);
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.5) hRecoEvtsEB->Fill(hardInvMass);
            if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.5) hRecoEvtsEE->Fill(hardInvMass);
         }
         hRecoNoTrgEvts->Fill(hardInvMass);
         if (fabs(gsfsc_eta[GSF_passHEEP[0]]) < 1.5) hRecoNoTrgEvtsEB->Fill(hardInvMass);
         if (fabs(gsfsc_eta[GSF_passHEEP[0]]) > 1.5) hRecoNoTrgEvtsEE->Fill(hardInvMass);
         ++evCounter;
        ///////////////////////////////////////////////////////////////////////
      } //END LOOP OVER EVENTS
        //////////////////////////////////////////////////////////////////////

     //////////////////////////////////////////////////////////////////////////
   } //END LOOP OVER FILES
     //////////////////////////////////////////////////////////////////////////
   hAcc = (TH1F*)hGenEvtsInAcc->Clone("hAcc");
   hAcc->Divide(hGenEvts);
   hAccEle = (TH1F*)hGenEvtsEleInAcc->Clone("hAccEle");
   hAccEle->Divide(hGenEvts);
   hAccEleEB = (TH1F*)hGenEvtsEleInAccEB->Clone("hAccEleEB");
   hAccEleEB->Divide(hGenEvts);
   hAccEleEE = (TH1F*)hGenEvtsEleInAccEE->Clone("hAccEleEE");
   hAccEleEE->Divide(hGenEvts);
   hAccMu = (TH1F*)hGenEvtsMuInAcc->Clone("hAccMu");
   hAccMu->Divide(hGenEvts);
   hAccTimesEff = (TH1F*)hRecoEvts->Clone("hAccTimesEff");
   hAccTimesEff->Divide(hGenEvts);
   hAccTimesEffEB = (TH1F*)hRecoEvtsEB->Clone("hAccTimesEffEB");
   hAccTimesEffEB->Divide(hGenEvts);
   hAccTimesEffEE = (TH1F*)hRecoEvtsEE->Clone("hAccTimesEffEE");
   hAccTimesEffEE->Divide(hGenEvts);
   hAccTimesEffEle = (TH1F*)hRecoEleEvts->Clone("hAccTimesEffEle");
   hAccTimesEffEle->Divide(hGenEvts);
   hAccTimesEffEleEB = (TH1F*)hRecoEleEvtsEB->Clone("hAccTimesEffEleEB");
   hAccTimesEffEleEB->Divide(hGenEvts);
   hAccTimesEffEleEE = (TH1F*)hRecoEleEvtsEE->Clone("hAccTimesEffEleEE");
   hAccTimesEffEleEE->Divide(hGenEvts);
   hAccTimesEffMu = (TH1F*)hRecoMuEvts->Clone("hAccTimesEffMu");
   hAccTimesEffMu->Divide(hGenEvts);
   hAccTimesEffNoTrg = (TH1F*)hRecoNoTrgEvts->Clone("hAccTimesEffNoTrg");
   hAccTimesEffNoTrg->Divide(hGenEvts);
   //hAccTimesEffNoTrgEB = (TH1F*)hRecoNoTrgEvtsEB->Clone("hAccTimesEffNoTrgEB");
   //hAccTimesEffNoTrgEB->Divide(hGenEvts);
   //hAccTimesEffNoTrgEE = (TH1F*)hRecoNoTrgEvtsEE->Clone("hAccTimesEffNoTrgEE");
   //hAccTimesEffNoTrgEE->Divide(hGenEvts);
   hAccTimesEffNoTrgEle = (TH1F*)hRecoNoTrgEleEvts->Clone("hAccTimesEffNoTrgEle");
   hAccTimesEffNoTrgEle->Divide(hGenEvts);
   hAccTimesEffNoTrgEleEB = (TH1F*)hRecoNoTrgEleEvtsEB->Clone("hAccTimesEffNoTrgEleEB");
   hAccTimesEffNoTrgEleEB->Divide(hGenEvts);
   hAccTimesEffNoTrgEleEE = (TH1F*)hRecoNoTrgEleEvtsEE->Clone("hAccTimesEffNoTrgEleEE");
   hAccTimesEffNoTrgEleEE->Divide(hGenEvts);
   hAccTimesEffNoTrgMu = (TH1F*)hRecoNoTrgMuEvts->Clone("hAccTimesEffNoTrgMu");
   hAccTimesEffNoTrgMu->Divide(hGenEvts);
   hEffAftTrg = (TH1F*)hRecoEvts->Clone("hEffAftTrg");
   hEffAftTrg->Divide(hTrgEvts);
   hEffAftTrgEle = (TH1F*)hRecoEleEvts->Clone("hEffAftTrgEle");
   hEffAftTrgEle->Divide(hTrgEvts);
   hEffAftTrgEleEB = (TH1F*)hRecoEleEvtsEB->Clone("hEffAftTrgEleEB");
   hEffAftTrgEleEB->Divide(hTrgEvts);
   hEffAftTrgEleEE = (TH1F*)hRecoEleEvtsEE->Clone("hEffAftTrgEleEE");
   hEffAftTrgEleEE->Divide(hTrgEvts);
   hEffAftTrgMu = (TH1F*)hRecoMuEvts->Clone("hEffAftTrgMu");
   hEffAftTrgMu->Divide(hTrgEvts);
   hTrgRecoVsReco = (TH1F*)hRecoEvts->Clone("hTrgRecoVsReco");
   hTrgRecoVsReco->Divide(hRecoNoTrgEvts);
   //hTrgRecoVsRecoEB = (TH1F*)hRecoEvtsEB->Clone("hTrgRecoVsRecoEB");
   //hTrgRecoVsRecoEB->Divide(hRecoNoTrgEvtsEB);
   //hTrgRecoVsRecoEE = (TH1F*)hRecoEvtsEE->Clone("hTrgRecoVsRecoEE");
   //hTrgRecoVsRecoEE->Divide(hRecoNoTrgEvtsEE);
   hTrgRecoVsRecoEle = (TH1F*)hRecoEleEvts->Clone("hTrgRecoVsRecoEle");
   hTrgRecoVsRecoEle->Divide(hRecoNoTrgEleEvts);
   hTrgRecoVsRecoEleEB = (TH1F*)hRecoEleEvtsEB->Clone("hTrgRecoVsRecoEleEB");
   hTrgRecoVsRecoEleEB->Divide(hRecoNoTrgEleEvtsEB);
   hTrgRecoVsRecoEleEE = (TH1F*)hRecoEleEvtsEE->Clone("hTrgRecoVsRecoEleEE");
   hTrgRecoVsRecoEleEE->Divide(hRecoNoTrgEleEvtsEE);
   hTrgRecoVsRecoMu = (TH1F*)hRecoMuEvts->Clone("hTrgRecoVsRecoMu");
   hTrgRecoVsRecoMu->Divide(hRecoNoTrgMuEvts);
   hFltrMatchVsReco_mu22ph22_phLeg = (TH1F*)hFltrMatchEvts_mu22ph22_phLeg->Clone("hFltrMatchVsReco_mu22ph22_phLeg");
   hFltrMatchVsReco_mu22ph22_phLeg->Divide(hRecoNoTrgEvts);
   hFltrMatchVsReco_mu22ph22_muLeg = (TH1F*)hFltrMatchEvts_mu22ph22_muLeg->Clone("hFltrMatchVsReco_mu22ph22_muLeg");
   hFltrMatchVsReco_mu22ph22_muLeg->Divide(hRecoNoTrgEvts);
   hFltrMatchVsReco_mu40 = (TH1F*)hFltrMatchEvts_mu40->Clone("hFltrMatchVsReco_mu40");
   hFltrMatchVsReco_mu40->Divide(hRecoNoTrgEvts);
   hTrgEff = (TH1F*)hTrgEvts->Clone("hTrgEff");
   hTrgEff->Divide(hGenEvts);
   hTrgEff_mu22ph22 = (TH1F*)hTrgEvts_mu22ph22->Clone("hTrgEff_mu22ph22");
   hTrgEff_mu22ph22->Divide(hGenEvts);
   hTrgEff_mu22ph22_ele_pt = (TH1F*)hTrgEvts_mu22ph22_ele_pt->Clone("hTrgEff_mu22ph22_ele_pt");
   hTrgEff_mu22ph22_ele_pt->Divide(hGenEvts_ele_pt);
   hTrgEff_mu22ph22_ele_eta = (TH1F*)hTrgEvts_mu22ph22_ele_eta->Clone("hTrgEff_mu22ph22_ele_eta");
   hTrgEff_mu22ph22_ele_eta->Divide(hGenEvts_ele_eta);
   hTrgEff_mu22ph22_mu_pt = (TH1F*)hTrgEvts_mu22ph22_mu_pt->Clone("hTrgEff_mu22ph22_mu_pt");
   hTrgEff_mu22ph22_mu_pt->Divide(hGenEvts_mu_pt);
   hTrgEff_mu22ph22_mu_eta = (TH1F*)hTrgEvts_mu22ph22_mu_eta->Clone("hTrgEff_mu22ph22_mu_eta");
   hTrgEff_mu22ph22_mu_eta->Divide(hGenEvts_mu_eta);
   hTrgEff_mu40 = (TH1F*)hTrgEvts_mu40->Clone("hTrgEff_mu40");
   hTrgEff_mu40->Divide(hGenEvts);
   hTrgEff_mu40_ele_pt = (TH1F*)hTrgEvts_mu40_ele_pt->Clone("hTrgEff_mu40_ele_pt");
   hTrgEff_mu40_ele_pt->Divide(hGenEvts_ele_pt);
   hTrgEff_mu40_ele_eta = (TH1F*)hTrgEvts_mu40_ele_eta->Clone("hTrgEff_mu40_ele_eta");
   hTrgEff_mu40_ele_eta->Divide(hGenEvts_ele_eta);
   hTrgEff_mu40_mu_pt = (TH1F*)hTrgEvts_mu40_mu_pt->Clone("hTrgEff_mu40_mu_pt");
   hTrgEff_mu40_mu_pt->Divide(hGenEvts_mu_pt);
   hTrgEff_mu40_mu_eta = (TH1F*)hTrgEvts_mu40_mu_eta->Clone("hTrgEff_mu40_mu_eta");
   hTrgEff_mu40_mu_eta->Divide(hGenEvts_mu_eta);
   hPassPtEff = (TH1F*)hPassPtCuts->Clone("hPassPtEff");
   hPassPtEff->Divide(hGenEvts);
   hPassEtaEff = (TH1F*)hPassEtaCuts->Clone("hPassEtaEff");
   hPassEtaEff->Divide(hGenEvts);
   hPassPtEtaEff = (TH1F*)hPassPtEtaCuts->Clone("hPassPtEtaEff");
   hPassPtEtaEff->Divide(hGenEvts);

   // plot
   TCanvas *accTimesEffPlot = new TCanvas("accTimesEffPlot", "acc x eff", 100, 100, 600, 600);
   TPad *accTimesEffPad = new TPad("accTimesEffPad", "acc x eff pad", 0., 0., 1., 1.);
   accTimesEffPad->SetBottomMargin(0.12);
   accTimesEffPad->SetBorderMode(0);
   accTimesEffPad->SetBorderSize(2);
   accTimesEffPad->SetFrameBorderMode(0);
   accTimesEffPad->SetFillColor(0);
   accTimesEffPad->SetFrameFillColor(0);
   accTimesEffPad->SetLeftMargin(0.11);
   accTimesEffPad->SetRightMargin(0.09);
   accTimesEffPad->SetTopMargin(0.08);
   accTimesEffPad->SetTickx(1);
   accTimesEffPad->SetTicky(1);
   accTimesEffPad->Draw();
   accTimesEffPad->cd();

   gStyle->SetTitleFont(font);
   gStyle->SetLabelFont(font);
   gStyle->SetLegendFont(font);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetOptFit(1111);
   gStyle->SetTitleXOffset(1.);
   gStyle->SetTitleYOffset(1.3);
   gPad->SetTicks(1, 1);
   gPad->SetGrid(1, 1);

   TH1F* hAccTimesEff2 = (TH1F*)hAccTimesEff->Clone("hAccTimesEff2");

   TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]/ (x + [2]) + [3]*x", 10., 5010.);
   //TF1 *fitFuncEB = new TF1("fitFuncEB", "[0] + [1]/ (x + [2])", 10., 5010.);
   TF1 *fitFuncEB = new TF1("fitFuncEB", "[0] + [1]/ (x + [2]) + [3]*x", 10., 5010.);
   TF1 *fitFuncEE = new TF1("fitFuncEE", "[0] + [1]/ (x*x + [2])", 10., 5010.);
   fitFunc->SetLineColor(kBlue);
   fitFuncEB->SetLineColor(kBlue);
   fitFuncEE->SetLineColor(kBlue);
   hAccTimesEff->Fit("fitFunc", "", "", 480., 5010.);
   hAccTimesEffEB->Fit("fitFuncEB", "", "", 480., 5010.);
   hAccTimesEffEE->Fit("fitFuncEE", "", "", 480., 5010.);
   cout << "Chi^2 / NDF: " << fitFunc->GetChisquare() << " / " << fitFunc->GetNDF() << ", prob: " << fitFunc->GetProb() << endl;
   cout << "Chi^2 / NDF EB: " << fitFuncEB->GetChisquare() << " / " << fitFuncEB->GetNDF() << ", prob: " << fitFuncEB->GetProb() << endl;
   cout << "Chi^2 / NDF EE: " << fitFuncEE->GetChisquare() << " / " << fitFuncEE->GetNDF() << ", prob: " << fitFuncEE->GetProb() << endl;

   hAccTimesEff->GetYaxis()->SetTitle("acc x eff");
   hAccTimesEff->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEff->Draw();
   TLatex *tex = new TLatex(0.22, 0.21, "P(M|p0,p1,p2,p3) = p0 + #frac{p1}{M+p2} + p3*M");
   tex->SetNDC();
   tex->SetTextFont(font);
   tex->SetLineWidth(2);
   tex->SetTextSize(0.03);
   tex->Draw();
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.17, 0.85, "trg + electron + muon");
   tex->DrawLatex(0.17, 0.80, "trg: " + triggerName);

   TCanvas *accTimesEffPlotEB = new TCanvas("accTimesEffPlotEB", "acc x eff, barrel electron + muon", 100, 100, 600, 600);
   TPad *accTimesEffPadEB = (TPad*)accTimesEffPad->Clone("accTimesEffPadEB");
   accTimesEffPadEB->Draw(); 
   accTimesEffPadEB->cd();
   hAccTimesEffEB->GetYaxis()->SetTitle("acc x eff");
   hAccTimesEffEB->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEffEB->Draw();
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   //tex->DrawLatex(0.46, 0.21, "P(M|p0,p1,p2) = p0 + #frac{p1}{M+p2}");
   tex->DrawLatex(0.22, 0.21, "P(M|p0,p1,p2,p3) = p0 + #frac{p1}{M+p2} + p3*M");
   tex->DrawLatex(0.17, 0.85, "trg + barrel electron + muon");
   tex->DrawLatex(0.17, 0.80, "trg: " + triggerName);

   TCanvas *accTimesEffPlotEE = new TCanvas("accTimesEffPlotEE", "acc x eff, endcap electron + muon", 100, 100, 600, 600);
   TPad *accTimesEffPadEE = (TPad*)accTimesEffPad->Clone("accTimesEffPadEE");
   accTimesEffPadEE->Draw(); 
   accTimesEffPadEE->cd();
   hAccTimesEffEE->GetYaxis()->SetTitle("acc x eff");
   hAccTimesEffEE->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEffEE->Draw();
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.45, 0.38, "P(M|p0,p1,p2) = p0 + #frac{p1}{M^{2}+p2}");
   tex->DrawLatex(0.17, 0.85, "trg + endcap electron + muon");
   tex->DrawLatex(0.17, 0.80, "trg: " + triggerName);

   TCanvas *accTimesEffObjPlot = new TCanvas("accTimesEffObjPlot", "acc x eff, objects", 100, 100, 600, 600);
   TPad *accTimesEffObjPad = (TPad*)accTimesEffPad->Clone("accTimesEffObjPad");
   accTimesEffObjPad->Draw(); 
   accTimesEffObjPad->cd();
   hAccTimesEffEle->GetYaxis()->SetTitle("acc x eff");
   hAccTimesEffEle->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEffEle->SetMarkerStyle(kFullSquare);
   hAccTimesEffEle->SetMarkerColor(kViolet);
   hAccTimesEffEle->SetLineColor(kViolet);
   hAccTimesEffEle->Draw();
   hAccTimesEffEleEB->SetMarkerStyle(kFullTriangleUp);
   hAccTimesEffEleEB->SetMarkerColor(kRed);
   hAccTimesEffEleEB->SetLineColor(kRed);
   hAccTimesEffEleEB->Draw("same");
   hAccTimesEffEleEE->SetMarkerStyle(kFullTriangleDown);
   hAccTimesEffEleEE->SetMarkerColor(kBlue);
   hAccTimesEffEleEE->SetLineColor(kBlue);
   hAccTimesEffEleEE->Draw("same");
   hAccTimesEffMu->SetMarkerStyle(34);
   hAccTimesEffMu->SetMarkerColor(kGreen+1);
   hAccTimesEffMu->SetLineColor(kGreen+1);
   hAccTimesEffMu->Draw("same");
   hAccTimesEff2->Draw("same");
   TLegend* legend = new TLegend(0.592, 0.279, 0.881, 0.467);
   legend->SetTextFont(font);
   legend->SetTextSize(0.03);
   legend->SetBorderSize(0);
   legend->SetLineColor(1);
   legend->SetLineStyle(1);
   legend->SetLineWidth(1);
   legend->SetFillColor(19);
   legend->SetFillStyle(0);
   legend->AddEntry(hAccTimesEff, "total acc x eff");
   legend->AddEntry(hAccTimesEffMu, "muons");
   legend->AddEntry(hAccTimesEffEle, "all electrons");
   legend->AddEntry(hAccTimesEffEleEB, "barrel electrons");
   legend->AddEntry(hAccTimesEffEleEE, "endcap electrons");
   legend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.17, 0.85, "trg: " + triggerName);

   // acc x eff with no trg applied
   TCanvas *accTimesEffNoTrgObjPlot = new TCanvas("accTimesEffNoTrgObjPlot", "acc x eff, no trigger, objects", 100, 100, 600, 600);
   TPad *accTimesEffNoTrgObjPad = (TPad*)accTimesEffPad->Clone("accTimesEffNoTrgObjPad");
   accTimesEffNoTrgObjPad->Draw(); 
   accTimesEffNoTrgObjPad->cd();
   hAccTimesEffNoTrgEle->GetYaxis()->SetTitle("acc x eff");
   hAccTimesEffNoTrgEle->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEffNoTrgEle->SetMarkerStyle(kFullSquare);
   hAccTimesEffNoTrgEle->SetMarkerColor(kViolet);
   hAccTimesEffNoTrgEle->SetLineColor(kViolet);
   hAccTimesEffNoTrgEle->Draw();
   hAccTimesEffNoTrgEleEB->SetMarkerStyle(kFullTriangleUp);
   hAccTimesEffNoTrgEleEB->SetMarkerColor(kRed);
   hAccTimesEffNoTrgEleEB->SetLineColor(kRed);
   hAccTimesEffNoTrgEleEB->Draw("same");
   hAccTimesEffNoTrgEleEE->SetMarkerStyle(kFullTriangleDown);
   hAccTimesEffNoTrgEleEE->SetMarkerColor(kBlue);
   hAccTimesEffNoTrgEleEE->SetLineColor(kBlue);
   hAccTimesEffNoTrgEleEE->Draw("same");
   hAccTimesEffNoTrgMu->SetMarkerStyle(34);
   hAccTimesEffNoTrgMu->SetMarkerColor(kGreen+1);
   hAccTimesEffNoTrgMu->SetLineColor(kGreen+1);
   hAccTimesEffNoTrgMu->Draw("same");
   //hPassPtEff->SetMarkerStyle(24);
   //hPassPtEff->SetMarkerColor(kGreen+2);
   //hPassPtEff->SetLineColor(kGreen+2);
   //hPassPtEff->Draw("same");
   //hPassEtaEff->SetMarkerStyle(25);
   //hPassEtaEff->SetMarkerColor(kGreen+3);
   //hPassEtaEff->SetLineColor(kGreen+3);
   //hPassEtaEff->Draw("same");
   //hPassPtEtaEff->SetMarkerStyle(26);
   //hPassPtEtaEff->SetMarkerColor(kGreen+4);
   //hPassPtEtaEff->SetLineColor(kGreen+4);
   //hPassPtEtaEff->Draw("same");
   hAccTimesEffNoTrg->Draw("same");
   TLegend* accXeffNoTrg = (TLegend*)legend->Clone("effAftTrgLegend");
   accXeffNoTrg->Clear();
   //accXeffNoTrg->AddEntry(hPassPtEff, "passing p_{T} cuts");
   //accXeffNoTrg->AddEntry(hPassEtaEff, "passing #eta cuts");
   //accXeffNoTrg->AddEntry(hPassPtEtaEff, "passing p_{T} and #eta cuts");
   accXeffNoTrg->AddEntry(hAccTimesEffNoTrg, "total acc x eff");
   accXeffNoTrg->AddEntry(hAccTimesEffNoTrgMu, "muons");
   accXeffNoTrg->AddEntry(hAccTimesEffNoTrgEle, "all electrons");
   accXeffNoTrg->AddEntry(hAccTimesEffNoTrgEleEB, "barrel electrons");
   accXeffNoTrg->AddEntry(hAccTimesEffNoTrgEleEE, "endcap electrons");
   accXeffNoTrg->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.17, 0.85, "trg: none");
   tex->DrawLatex(0.14, 0.15, "#mu_{Id} for " + triggerName);

   // efficiency on triggered events
   TCanvas *effAftTrgPlot = new TCanvas("effAftTrgPlot", "efficiency after trigger", 100, 100, 600, 600);
   TPad *effAftTrgPad = (TPad*)accTimesEffPad->Clone("effAftTrgPad");
   effAftTrgPad->Draw(); 
   effAftTrgPad->cd();
   hEffAftTrgEle->GetYaxis()->SetTitle("eff for triggered events");
   hEffAftTrgEle->GetYaxis()->SetRangeUser(0., 1.);
   hEffAftTrgEle->SetMarkerStyle(kFullSquare);
   hEffAftTrgEle->SetMarkerColor(kViolet);
   hEffAftTrgEle->SetLineColor(kViolet);
   hEffAftTrgEle->Draw();
   hEffAftTrgEleEB->SetMarkerStyle(kFullTriangleUp);
   hEffAftTrgEleEB->SetMarkerColor(kRed);
   hEffAftTrgEleEB->SetLineColor(kRed);
   hEffAftTrgEleEB->Draw("same");
   hEffAftTrgEleEE->SetMarkerStyle(kFullTriangleDown);
   hEffAftTrgEleEE->SetMarkerColor(kBlue);
   hEffAftTrgEleEE->SetLineColor(kBlue);
   hEffAftTrgEleEE->Draw("same");
   hEffAftTrgMu->SetMarkerStyle(34);
   hEffAftTrgMu->SetMarkerColor(kGreen+1);
   hEffAftTrgMu->SetLineColor(kGreen+1);
   hEffAftTrgMu->Draw("same");
   hEffAftTrg->Draw("same");
   TLegend* effAftTrgLegend = (TLegend*)legend->Clone("effAftTrgLegend");
   effAftTrgLegend->Clear();
   effAftTrgLegend->AddEntry(hEffAftTrg, "total eff after trigger");
   effAftTrgLegend->AddEntry(hEffAftTrgMu, "muons");
   effAftTrgLegend->AddEntry(hEffAftTrgEle, "all electrons");
   effAftTrgLegend->AddEntry(hEffAftTrgEleEB, "barrel electrons");
   effAftTrgLegend->AddEntry(hEffAftTrgEleEE, "endcap electrons");
   effAftTrgLegend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.14, 0.15, "trg: " + triggerName);

   // acceptance
   TCanvas *accPlot = new TCanvas("accPlot", "acc", 100, 100, 600, 600);
   TPad *accPad = (TPad*)accTimesEffPad->Clone("accPad");
   accPad->Draw(); 
   accPad->cd();
   hAccEle->GetYaxis()->SetTitle("acc");
   hAccEle->GetYaxis()->SetRangeUser(0., 1.);
   hAccEle->SetMarkerStyle(kFullSquare);
   hAccEle->SetMarkerColor(kViolet);
   hAccEle->SetLineColor(kViolet);
   hAccEle->Draw();
   hAccEleEB->SetMarkerStyle(kFullTriangleUp);
   hAccEleEB->SetMarkerColor(kRed);
   hAccEleEB->SetLineColor(kRed);
   hAccEleEB->Draw("same");
   hAccEleEE->SetMarkerStyle(kFullTriangleDown);
   hAccEleEE->SetMarkerColor(kBlue);
   hAccEleEE->SetLineColor(kBlue);
   hAccEleEE->Draw("same");
   hAccMu->SetMarkerStyle(34);
   hAccMu->SetMarkerColor(kGreen+1);
   hAccMu->SetLineColor(kGreen+1);
   hAccMu->Draw("same");
   hAcc->Draw("same");
   TLegend* accLegend = (TLegend*)legend->Clone("accLegend");
   accLegend->Clear();
   accLegend->AddEntry(hAcc, "total acceptance");
   accLegend->AddEntry(hAccMu, "muons");
   accLegend->AddEntry(hAccEle, "all electrons");
   accLegend->AddEntry(hAccEleEB, "barrel electrons");
   accLegend->AddEntry(hAccEleEE, "endcap electrons");
   accLegend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.14, 0.15, "#mu_{Id} for " + triggerName);

   // reco with trg vs. reco
   TCanvas *trgRecoVsRecoPlot = new TCanvas("trgRecoVsRecoPlot", "reco with trg vs. reco", 100, 100, 600, 600);
   TPad *trgRecoVsRecoPad = (TPad*)accTimesEffPad->Clone("trgRecoVsRecoPad");
   trgRecoVsRecoPad->Draw(); 
   trgRecoVsRecoPad->cd();
   hTrgRecoVsRecoEle->GetYaxis()->SetTitle("trg eff for selected objects");
   hTrgRecoVsRecoEle->GetYaxis()->SetRangeUser(0., 1.);
   hTrgRecoVsRecoEle->SetMarkerStyle(kFullSquare);
   hTrgRecoVsRecoEle->SetMarkerColor(kViolet);
   hTrgRecoVsRecoEle->SetLineColor(kViolet);
   hTrgRecoVsRecoEle->Draw();
   //hTrgRecoVsRecoEleEB->SetMarkerStyle(kFullTriangleUp);
   //hTrgRecoVsRecoEleEB->SetMarkerColor(kRed);
   //hTrgRecoVsRecoEleEB->SetLineColor(kRed);
   //hTrgRecoVsRecoEleEB->Draw("same");
   //hTrgRecoVsRecoEleEE->SetMarkerStyle(kFullTriangleDown);
   //hTrgRecoVsRecoEleEE->SetMarkerColor(kBlue);
   //hTrgRecoVsRecoEleEE->SetLineColor(kBlue);
   //hTrgRecoVsRecoEleEE->Draw("same");
   hTrgRecoVsRecoMu->SetMarkerStyle(34);
   hTrgRecoVsRecoMu->SetMarkerColor(kGreen+1);
   hTrgRecoVsRecoMu->SetLineColor(kGreen+1);
   hTrgRecoVsRecoMu->Draw("same");
   hTrgRecoVsReco->Draw("same");
   TLegend* trgRecoVsRecoLegend = (TLegend*)legend->Clone("trgRecoVsRecoLegend");
   trgRecoVsRecoLegend->Clear();
   trgRecoVsRecoLegend->AddEntry(hTrgRecoVsReco, "total");
   trgRecoVsRecoLegend->AddEntry(hTrgRecoVsRecoMu, "muons");
   trgRecoVsRecoLegend->AddEntry(hTrgRecoVsRecoEle, "electrons");
   //trgRecoVsRecoLegend->AddEntry(hTrgRecoVsRecoEleEB, "barrel electrons");
   //trgRecoVsRecoLegend->AddEntry(hTrgRecoVsRecoEleEE, "endcap electrons");
   trgRecoVsRecoLegend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.14, 0.15, "trg: " + triggerName);

   TCanvas *trgEffPlot = new TCanvas("trgEffPlot", "trigger eff", 100, 100, 600, 600);
   TPad *trgEffPad = (TPad*)accTimesEffPad->Clone("trgEffPad");
   trgEffPad->Draw(); 
   trgEffPad->cd();
   hTrgEff_mu22ph22->GetYaxis()->SetTitle("trg eff");
   hTrgEff_mu22ph22->GetYaxis()->SetRangeUser(0., 1.);
   hTrgEff_mu22ph22->SetMarkerStyle(21);
   hTrgEff_mu22ph22->SetMarkerColor(kMagenta);
   hTrgEff_mu22ph22->SetLineColor(kMagenta);
   hTrgEff_mu22ph22->Draw();
   hTrgEff_mu40->SetMarkerStyle(20);
   hTrgEff_mu40->SetMarkerColor(kCyan+1);
   hTrgEff_mu40->SetLineColor(kCyan+1);
   hTrgEff_mu40->Draw("same");
   TLegend* trgLegend = (TLegend*)legend->Clone("trgLegend");
   trgLegend->Clear();
   trgLegend->SetX1(0.38);
   trgLegend->SetX2(0.67);
   trgLegend->AddEntry(hTrgEff_mu22ph22, "eff HLT_Mu22_Photon22_CaloIdL");
   trgLegend->AddEntry(hTrgEff_mu40, "eff HLT_Mu40_eta2p1");
   trgLegend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");

   TCanvas *trgEffPlot_ele_pt = new TCanvas("trgEffPlot_ele_pt", "trigger eff wrt ele pt", 100, 100, 600, 600);
   TPad *trgEffPad_ele_pt = (TPad*)accTimesEffPad->Clone("trgEffPad_ele_pt");
   trgEffPad_ele_pt->Draw(); 
   trgEffPad_ele_pt->cd();
   hTrgEff_mu22ph22_ele_pt->GetYaxis()->SetTitle("trg eff");
   hTrgEff_mu22ph22_ele_pt->GetYaxis()->SetRangeUser(0., 1.);
   hTrgEff_mu22ph22_ele_pt->SetMarkerStyle(21);
   hTrgEff_mu22ph22_ele_pt->SetMarkerColor(kMagenta);
   hTrgEff_mu22ph22_ele_pt->SetLineColor(kMagenta);
   hTrgEff_mu22ph22_ele_pt->Draw();
   hTrgEff_mu40_ele_pt->SetMarkerStyle(20);
   hTrgEff_mu40_ele_pt->SetMarkerColor(kCyan+1);
   hTrgEff_mu40_ele_pt->SetLineColor(kCyan+1);
   hTrgEff_mu40_ele_pt->Draw("same");
   TLegend* trgLegend_ele_pt = (TLegend*)legend->Clone("trgLegend_ele_pt");
   trgLegend_ele_pt->Clear();
   trgLegend_ele_pt->SetX1(0.38);
   trgLegend_ele_pt->SetX2(0.67);
   trgLegend_ele_pt->AddEntry(hTrgEff_mu22ph22_ele_pt, "eff HLT_Mu22_Photon22_CaloIdL");
   trgLegend_ele_pt->AddEntry(hTrgEff_mu40_ele_pt, "eff HLT_Mu40_eta2p1");
   trgLegend_ele_pt->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");

   TCanvas *trgEffPlot_ele_eta = new TCanvas("trgEffPlot_ele_eta", "trigger eff wrt ele eta", 100, 100, 600, 600);
   TPad *trgEffPad_ele_eta = (TPad*)accTimesEffPad->Clone("trgEffPad_ele_eta");
   trgEffPad_ele_eta->Draw(); 
   trgEffPad_ele_eta->cd();
   hTrgEff_mu22ph22_ele_eta->GetYaxis()->SetTitle("trg eff");
   hTrgEff_mu22ph22_ele_eta->GetYaxis()->SetRangeUser(0., 1.);
   hTrgEff_mu22ph22_ele_eta->SetMarkerStyle(21);
   hTrgEff_mu22ph22_ele_eta->SetMarkerColor(kMagenta);
   hTrgEff_mu22ph22_ele_eta->SetLineColor(kMagenta);
   hTrgEff_mu22ph22_ele_eta->Draw();
   hTrgEff_mu40_ele_eta->SetMarkerStyle(20);
   hTrgEff_mu40_ele_eta->SetMarkerColor(kCyan+1);
   hTrgEff_mu40_ele_eta->SetLineColor(kCyan+1);
   hTrgEff_mu40_ele_eta->Draw("same");
   TLegend* trgLegend_ele_eta = (TLegend*)legend->Clone("trgLegend_ele_eta");
   trgLegend_ele_eta->Clear();
   trgLegend_ele_eta->SetX1(0.38);
   trgLegend_ele_eta->SetX2(0.67);
   trgLegend_ele_eta->AddEntry(hTrgEff_mu22ph22_ele_eta, "eff HLT_Mu22_Photon22_CaloIdL");
   trgLegend_ele_eta->AddEntry(hTrgEff_mu40_ele_eta, "eff HLT_Mu40_eta2p1");
   trgLegend_ele_eta->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");

   TCanvas *trgEffPlot_mu_pt = new TCanvas("trgEffPlot_mu_pt", "trigger eff wrt muon pt", 100, 100, 600, 600);
   TPad *trgEffPad_mu_pt = (TPad*)accTimesEffPad->Clone("trgEffPad_mu_pt");
   trgEffPad_mu_pt->Draw(); 
   trgEffPad_mu_pt->cd();
   hTrgEff_mu22ph22_mu_pt->GetYaxis()->SetTitle("trg eff");
   hTrgEff_mu22ph22_mu_pt->GetYaxis()->SetRangeUser(0., 1.);
   hTrgEff_mu22ph22_mu_pt->SetMarkerStyle(21);
   hTrgEff_mu22ph22_mu_pt->SetMarkerColor(kMagenta);
   hTrgEff_mu22ph22_mu_pt->SetLineColor(kMagenta);
   hTrgEff_mu22ph22_mu_pt->Draw();
   hTrgEff_mu40_mu_pt->SetMarkerStyle(20);
   hTrgEff_mu40_mu_pt->SetMarkerColor(kCyan+1);
   hTrgEff_mu40_mu_pt->SetLineColor(kCyan+1);
   hTrgEff_mu40_mu_pt->Draw("same");
   TLegend* trgLegend_mu_pt = (TLegend*)legend->Clone("trgLegend_mu_pt");
   trgLegend_mu_pt->Clear();
   trgLegend_mu_pt->SetX1(0.38);
   trgLegend_mu_pt->SetX2(0.67);
   trgLegend_mu_pt->AddEntry(hTrgEff_mu22ph22_mu_pt, "eff HLT_Mu22_Photon22_CaloIdL");
   trgLegend_mu_pt->AddEntry(hTrgEff_mu40_mu_pt, "eff HLT_Mu40_eta2p1");
   trgLegend_mu_pt->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");

   TCanvas *trgEffPlot_mu_eta = new TCanvas("trgEffPlot_mu_eta", "trigger eff wrt muon eta", 100, 100, 600, 600);
   TPad *trgEffPad_mu_eta = (TPad*)accTimesEffPad->Clone("trgEffPad_mu_eta");
   trgEffPad_mu_eta->Draw(); 
   trgEffPad_mu_eta->cd();
   hTrgEff_mu22ph22_mu_eta->GetYaxis()->SetTitle("trg eff");
   hTrgEff_mu22ph22_mu_eta->GetYaxis()->SetRangeUser(0., 1.);
   hTrgEff_mu22ph22_mu_eta->SetMarkerStyle(21);
   hTrgEff_mu22ph22_mu_eta->SetMarkerColor(kMagenta);
   hTrgEff_mu22ph22_mu_eta->SetLineColor(kMagenta);
   hTrgEff_mu22ph22_mu_eta->Draw();
   hTrgEff_mu40_mu_eta->SetMarkerStyle(20);
   hTrgEff_mu40_mu_eta->SetMarkerColor(kCyan+1);
   hTrgEff_mu40_mu_eta->SetLineColor(kCyan+1);
   hTrgEff_mu40_mu_eta->Draw("same");
   TLegend* trgLegend_mu_eta = (TLegend*)legend->Clone("trgLegend_mu_eta");
   trgLegend_mu_eta->Clear();
   trgLegend_mu_eta->SetX1(0.38);
   trgLegend_mu_eta->SetX2(0.67);
   trgLegend_mu_eta->AddEntry(hTrgEff_mu22ph22_mu_eta, "eff HLT_Mu22_Photon22_CaloIdL");
   trgLegend_mu_eta->AddEntry(hTrgEff_mu40_mu_eta, "eff HLT_Mu40_eta2p1");
   trgLegend_mu_eta->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");

   TCanvas *trgFltrEffPlot = new TCanvas("trgFltrEffPlot", "trigger filter eff wrt selected event", 100, 100, 600, 600);
   TPad *trgFltrEffPad = (TPad*)accTimesEffPad->Clone("trgFltrEffPad");
   trgFltrEffPad->Draw(); 
   trgFltrEffPad->cd();
   TLegend* trgFltrLegend = (TLegend*)legend->Clone("trgFltrLegend");
   trgFltrLegend->Clear();
   trgFltrLegend->SetX1(0.33);
   trgFltrLegend->SetX2(0.62);
   if (triggerInd == 0) {
      hFltrMatchVsReco_mu22ph22_muLeg->GetYaxis()->SetTitle("eff");
      hFltrMatchVsReco_mu22ph22_muLeg->GetYaxis()->SetRangeUser(0.75, 1.);
      hFltrMatchVsReco_mu22ph22_muLeg->SetMarkerStyle(kFullTriangleUp);
      hFltrMatchVsReco_mu22ph22_muLeg->SetMarkerColor(kGreen+1);
      hFltrMatchVsReco_mu22ph22_muLeg->SetLineColor(kGreen+1);
      hFltrMatchVsReco_mu22ph22_muLeg->Draw();
      hFltrMatchVsReco_mu22ph22_phLeg->SetMarkerStyle(kFullTriangleDown);
      hFltrMatchVsReco_mu22ph22_phLeg->SetMarkerColor(kBlue);
      hFltrMatchVsReco_mu22ph22_phLeg->SetLineColor(kBlue);
      hFltrMatchVsReco_mu22ph22_phLeg->Draw("same");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_mu22ph22_muLeg, "eff #mu-leg HLT_Mu22_Photon22_CaloIdL");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_mu22ph22_phLeg, "eff #gamma-leg HLT_Mu22_Photon22_CaloIdL");
   } else {
      hFltrMatchVsReco_mu40->GetYaxis()->SetTitle("eff");
      hFltrMatchVsReco_mu40->GetYaxis()->SetRangeUser(0.75, 1.);
      hFltrMatchVsReco_mu40->SetMarkerStyle(21);
      hFltrMatchVsReco_mu40->SetMarkerColor(kRed);
      hFltrMatchVsReco_mu40->SetLineColor(kRed);
      hFltrMatchVsReco_mu40->Draw();
      trgFltrLegend->AddEntry(hFltrMatchVsReco_mu40, "eff #mu-leg HLT_Mu40_eta2p1");
   }
   trgFltrLegend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");

   // safe in various file formats
   if (saveSpec) {
     if (saveAsPdf) {
        trgEffPlot->Print(plotDir + trgEffPlot->GetName() + ".pdf", "pdf");
        trgEffPlot_ele_pt->Print(plotDir + trgEffPlot_ele_pt->GetName() + ".pdf", "pdf");
        trgEffPlot_ele_eta->Print(plotDir + trgEffPlot_ele_eta->GetName() + ".pdf", "pdf");
        trgEffPlot_mu_pt->Print(plotDir + trgEffPlot_mu_pt->GetName() + ".pdf", "pdf");
        trgEffPlot_mu_eta->Print(plotDir + trgEffPlot_mu_eta->GetName() + ".pdf", "pdf");
        trgFltrEffPlot->Print(plotDir + trgFltrEffPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
        accTimesEffPlot->Print(plotDir + accTimesEffPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
        accTimesEffPlotEB->Print(plotDir + accTimesEffPlotEB->GetName() + "_" + triggerName + ".pdf", "pdf");
        accTimesEffPlotEE->Print(plotDir + accTimesEffPlotEE->GetName() + "_" + triggerName + ".pdf", "pdf");
        accTimesEffObjPlot->Print(plotDir + accTimesEffObjPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
        accTimesEffNoTrgObjPlot->Print(plotDir + accTimesEffNoTrgObjPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
        effAftTrgPlot->Print(plotDir + effAftTrgPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
        trgRecoVsRecoPlot->Print(plotDir + trgRecoVsRecoPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
        accPlot->Print(plotDir + accPlot->GetName() + "_" + triggerName + ".pdf", "pdf");
     }
     if (saveAsPng) {
        trgEffPlot->Print(plotDir + trgEffPlot->GetName() + ".png", "png");
        trgEffPlot_ele_pt->Print(plotDir + trgEffPlot_ele_pt->GetName() + ".png", "png");
        trgEffPlot_ele_eta->Print(plotDir + trgEffPlot_ele_eta->GetName() + ".png", "png");
        trgEffPlot_mu_pt->Print(plotDir + trgEffPlot_mu_pt->GetName() + ".png", "png");
        trgEffPlot_mu_eta->Print(plotDir + trgEffPlot_mu_eta->GetName() + ".png", "png");
        trgFltrEffPlot->Print(plotDir + trgFltrEffPlot->GetName() + "_" + triggerName + ".png", "png");
        accTimesEffPlot->Print(plotDir + accTimesEffPlot->GetName() + "_" + triggerName + ".png", "png");
        accTimesEffPlotEB->Print(plotDir + accTimesEffPlotEB->GetName() + "_" + triggerName + ".png", "png");
        accTimesEffPlotEE->Print(plotDir + accTimesEffPlotEE->GetName() + "_" + triggerName + ".png", "png");
        accTimesEffObjPlot->Print(plotDir + accTimesEffObjPlot->GetName() + "_" + triggerName + ".png", "png");
        accTimesEffNoTrgObjPlot->Print(plotDir + accTimesEffNoTrgObjPlot->GetName() + "_" + triggerName + ".png", "png");
        effAftTrgPlot->Print(plotDir + effAftTrgPlot->GetName() + "_" + triggerName + ".png", "png");
        trgRecoVsRecoPlot->Print(plotDir + trgRecoVsRecoPlot->GetName() + "_" + triggerName + ".png", "png");
        accPlot->Print(plotDir + accPlot->GetName() + "_" + triggerName + ".png", "png");
     }
     if (saveAsRoot) {
        trgEffPlot->Print(plotDir + trgEffPlot->GetName() + ".root", "root");
        trgEffPlot_ele_pt->Print(plotDir + trgEffPlot_ele_pt->GetName() + ".root", "root");
        trgEffPlot_ele_eta->Print(plotDir + trgEffPlot_ele_eta->GetName() + ".root", "root");
        trgEffPlot_mu_pt->Print(plotDir + trgEffPlot_mu_pt->GetName() + ".root", "root");
        trgEffPlot_mu_eta->Print(plotDir + trgEffPlot_mu_eta->GetName() + ".root", "root");
        trgFltrEffPlot->Print(plotDir + trgFltrEffPlot->GetName() + "_" + triggerName + ".root", "root");
        accTimesEffPlot->Print(plotDir + accTimesEffPlot->GetName() + "_" + triggerName + ".root", "root");
        accTimesEffPlotEB->Print(plotDir + accTimesEffPlotEB->GetName() + "_" + triggerName + ".root", "root");
        accTimesEffPlotEE->Print(plotDir + accTimesEffPlotEE->GetName() + "_" + triggerName + ".root", "root");
        accTimesEffObjPlot->Print(plotDir + accTimesEffObjPlot->GetName() + "_" + triggerName + ".root", "root");
        accTimesEffNoTrgObjPlot->Print(plotDir + accTimesEffNoTrgObjPlot->GetName() + "_" + triggerName + ".root", "root");
        effAftTrgPlot->Print(plotDir + effAftTrgPlot->GetName() + "_" + triggerName + ".root", "root");
        trgRecoVsRecoPlot->Print(plotDir + trgRecoVsRecoPlot->GetName() + "_" + triggerName + ".root", "root");
        accPlot->Print(plotDir + accPlot->GetName() + "_" + triggerName + ".root", "root");
     }
   }

   // write histos to file
   output->cd();
   hGenEvts->Write();
   hGenEvts_ele_pt->Write();
   hGenEvts_ele_eta->Write();
   hGenEvts_mu_pt->Write();
   hGenEvts_mu_eta->Write();
   hGenEvtsEleInAcc->Write();
   hGenEvtsEleInAccEB->Write();
   hGenEvtsEleInAccEE->Write();
   hGenEvtsMuInAcc->Write();
   hGenEvtsInAcc->Write();
   hTrgEvts->Write();
   hTrgEvts_mu22ph22->Write();
   hTrgEvts_mu22ph22_ele_pt->Write();
   hTrgEvts_mu22ph22_ele_eta->Write();
   hTrgEvts_mu22ph22_mu_pt->Write();
   hTrgEvts_mu22ph22_mu_eta->Write();
   hTrgEvts_mu40->Write();
   hTrgEvts_mu40_ele_pt->Write();
   hTrgEvts_mu40_ele_eta->Write();
   hTrgEvts_mu40_mu_pt->Write();
   hTrgEvts_mu40_mu_eta->Write();
   hRecoEvts->Write();
   hRecoEvtsEB->Write();
   hRecoEvtsEE->Write();
   hRecoEleEvts->Write();
   hRecoEleEvtsEB->Write();
   hRecoEleEvtsEE->Write();
   hRecoMuEvts->Write();
   hRecoNoTrgEvts->Write();
   hRecoNoTrgEvtsEB->Write();
   hRecoNoTrgEvtsEE->Write();
   hRecoNoTrgEleEvts->Write();
   hRecoNoTrgEleEvtsEB->Write();
   hRecoNoTrgEleEvtsEE->Write();
   hRecoNoTrgMuEvts->Write();
   hFltrMatchEvts_mu22ph22_phLeg->Write();
   hFltrMatchEvts_mu22ph22_muLeg->Write();
   hFltrMatchEvts_mu40->Write();
   hFltrMatchVsReco_mu22ph22_phLeg->Write();
   hFltrMatchVsReco_mu22ph22_muLeg->Write();
   hFltrMatchVsReco_mu40->Write();
   hAccEle->Write();
   hAccEleEB->Write();
   hAccEleEE->Write();
   hAccMu->Write();
   hTrgEff->Write();
   hTrgEff_mu22ph22->Write();
   hTrgEff_mu22ph22_ele_pt->Write();
   hTrgEff_mu22ph22_ele_eta->Write();
   hTrgEff_mu22ph22_mu_pt->Write();
   hTrgEff_mu22ph22_mu_eta->Write();
   hTrgEff_mu40->Write();
   hTrgEff_mu40_ele_pt->Write();
   hTrgEff_mu40_ele_eta->Write();
   hTrgEff_mu40_mu_pt->Write();
   hTrgEff_mu40_mu_eta->Write();
   hAccTimesEff->Write();
   hAccTimesEffEB->Write();
   hAccTimesEffEE->Write();
   hAccTimesEffEle->Write();
   hAccTimesEffEleEB->Write();
   hAccTimesEffEleEE->Write();
   hAccTimesEffMu->Write();
   hEffAftTrg->Write();
   hEffAftTrgEle->Write();
   hEffAftTrgEleEB->Write();
   hEffAftTrgEleEE->Write();
   hEffAftTrgMu->Write();
   hPassPtCuts->Write();
   hPassEtaCuts->Write();
   hPassPtEtaCuts->Write();
   hPassPtEff->Write();
   hPassEtaEff->Write();
   hPassPtEtaEff->Write();
   fitFunc->Write();
   fitFuncEB->Write();
   fitFuncEE->Write();

   output->Close();
   timer.Stop();
   timer.Print();
}

bool AccTimesEff::PassHEEP(const int &n)
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

bool AccTimesEff::PassHighPtMu(const int &n)
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


