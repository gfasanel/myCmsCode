#define accTimesEffRpv_cxx
#include "accTimesEff_rpv.h"
#include <TH1.h>
#include <TF1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

void AccTimesEffRpv::Loop()
{
   TStopwatch timer;
   timer.Start();
   // parameters /////////////////////////////////////////////////////////////
   vector<TString> files;
   files.push_back("file:////user/treis/mcsamples/RPVresonantToEMu_M-scan_TuneZ2star_8TeV-calchep-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_emuSkim_319984ev.root");
   string outfileName = "accTimesEffRpvHistos";

   // output file formats
   const bool saveSpec = 0;
   const bool saveAsPdf = 0;
   const bool saveAsPng = 1;
   const bool saveAsRoot = 0;
   TString plotDir = "./plots/";
   TString fileNameExtra = "";
   //fileNameExtra = "_test";

   unsigned int triggerInd = 1;  // 0: HLT_Mu22_Photon22_CaloIdL; 1: HLT_Mu40_eta2p1
   unsigned int eleDetRegion = 0;  // electron in subdetector. 0: EB or EE; 1: EB; 2: EE
   const bool useScaleFactors = 1;

   int font = 42; //62
   // selection cuts /////////////////////////////////////////////////////////
   float elePtCut = 35.;
   float muPtCut = 35.;
   float muEtaCut = 2.4;
   float minInvMass = 0.;

   // scale factors
   // epsilon_cand from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgCommissioningAndPhysicsDeliverables#Electron_reconstruction_effi_AN1
   float eps_cand_sf_0p8 = 0.990; // data/MC scale for epsilon_cand (>50GeV) |eta|<0.8
   float eps_cand_sf_0p8to1p4442 = 0.991; // data/MC scale for epsilon_cand (>50GeV) 0.8<|eta|<1.4442
   float eps_cand_sf_1p566to2p0 = 0.990; // data/MC scale for epsilon_cand (>50GeV) 1.566<|eta|<2.0
   float eps_cand_sf_2p0to2p5 = 0.998; // data/MC scale for epsilon_cand (>50GeV) 2.0<|eta|<2.5
   // HEEP eff scale factors from AN-13-359 Table 5 reReco
   float eps_heep_sf_eb_pt35 = 0.997; // HEEP eff scale factor
   float eps_heep_sf_eb_pt100 = 0.985; // HEEP eff scale factor
   float eps_heep_sf_ee_pt35 = 0.979; // HEEP eff scale factor
   float eps_heep_sf_ee_pt100 = 0.981; // HEEP eff scale factor
 
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
   TString eleDetRegSuffix = "";
   if (eleDetRegion == 1) eleDetRegSuffix = "_EB";
   else if (eleDetRegion == 2) eleDetRegSuffix = "_EE";

   TH1F *hGenEvts = new TH1F("hGenEvts", "hGenEvts", 27, -50., 2650.);
   //TH1F *hGenEvts = new TH1F("hGenEvts", "hGenEvts", 125, 0., 3510.);
   hGenEvts->GetXaxis()->SetTitle("M_{Z'}^{truth} (GeV)");
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
   TH1F *hRecoEvtsLepInAcc = (TH1F*)hGenEvts->Clone("hRecoEvtsLepInAcc");
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
   TH1F *hFltrMatchEvts_l1SingleEG12 = (TH1F*)hGenEvts->Clone("hFltrMatchEvts_l1SingleEG12");
   hFltrMatchEvts_l1SingleEG12->SetTitle("hFltrMatchEvts_l1SingleEG12");
   TH1F *hFltrMatchEvts_l1Mu3p5EG12 = (TH1F*)hGenEvts->Clone("hFltrMatchEvts_l1Mu3p5EG12");
   hFltrMatchEvts_l1Mu3p5EG12->SetTitle("hFltrMatchEvts_l1Mu3p5EG12");
   TH1F *hFltrMatchEvts_l1Mu16Eta2p1 = (TH1F*)hGenEvts->Clone("hFltrMatchEvts_l1Mu16Eta2p1");
   hFltrMatchEvts_l1Mu16Eta2p1->SetTitle("hFltrMatchEvts_l1Mu16Eta2p1");
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
   TH1F* hRecoLepAcc;
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
   TH1F* hFltrMatchVsReco_l1SingleEG12;
   TH1F* hFltrMatchVsReco_l1Mu3p5EG12;
   TH1F* hFltrMatchVsReco_l1Mu16Eta2p1;
   TH1F* hFltrMatchVsReco_mu22ph22_phLeg;
   TH1F* hFltrMatchVsReco_mu22ph22_muLeg;
   TH1F* hFltrMatchVsPrev_mu22ph22_phLeg;
   TH1F* hFltrMatchVsPrev_mu22ph22_muLeg;
   TH1F* hFltrMatchVsReco_mu40;
   TH1F* hFltrMatchVsPrev_mu40;

   // output file
   stringstream ssOutfile;
   ssOutfile << outfileName << "_" <<  triggerName << eleDetRegSuffix << fileNameExtra << ".root";
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

         if (emu_mass > 2050.) continue;

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
         //if (fabs(genelemom_pdgid[0]) < 7) continue;
         // fill the gen histograms
         hGenEvts->Fill(emu_mass);
         hGenEvts_ele_pt->Fill(genele_pt[0]);
         hGenEvts_ele_eta->Fill(genele_eta[0]);
         hGenEvts_mu_pt->Fill(genmu_pt[0]);
         hGenEvts_mu_eta->Fill(genmu_eta[0]);

         // fill the acc histograms
         bool passEleEta = 0;
         bool passMuEta = 0;
         if (fabs(hardGenEle_eta[0]) < 1.442 || (fabs(hardGenEle_eta[0]) > 1.56 && fabs(hardGenEle_eta[0]) < 2.5)) passEleEta = 1;
         if (fabs(hardGenMu_eta[0]) < muEtaCut) passMuEta = 1;
         if (passEleEta && passMuEta) hPassEtaCuts->Fill(emu_mass);
         if (hardGenEle_pt[0] > elePtCut) {
            if (fabs(hardGenEle_eta[0]) < 1.442) {
               hGenEvtsEleInAcc->Fill(emu_mass);
               hGenEvtsEleInAccEB->Fill(emu_mass);
            }
            else if (fabs(hardGenEle_eta[0]) > 1.56 && fabs(hardGenEle_eta[0]) < 2.5) {
               hGenEvtsEleInAcc->Fill(emu_mass);
               hGenEvtsEleInAccEE->Fill(emu_mass);
            }
            if (hardGenMu_pt[0] > muPtCut) {
               hPassPtCuts->Fill(emu_mass);
               if (passEleEta && passMuEta) hPassPtEtaCuts->Fill(emu_mass);
            }
         }
         if (hardGenMu_pt[0] > muPtCut && passMuEta) {
            hGenEvtsMuInAcc->Fill(emu_mass);
            if (passEleEta) {
               if (hardGenEle_pt[0] > elePtCut) hGenEvtsInAcc->Fill(emu_mass);
            }
         }

         // dummy scale factors
         float eleSf = 1.;
         float muSf = 1.;
         float trgSf = 1.;
         float trgEff = 1.;
         TString sigString = "sig";

         if (useScaleFactors) {
            // set trg scale factor according to generated muon
            WeightMuonRecoIsoTrigger(genmu_pt[0], genmu_eta[0], muSf, trgSf, trgEff, sigString);
         }

         // trigger?
         if (triggerInd == 0) triggerBit = HLT_Mu22_Photon22_CaloIdL;
         else triggerBit = HLT_Mu40_eta2p1;

         if (triggerBit) hTrgEvts->Fill(emu_mass);
         if (HLT_Mu22_Photon22_CaloIdL) {
            hTrgEvts_mu22ph22->Fill(emu_mass, trgSf);
            hTrgEvts_mu22ph22_ele_pt->Fill(genele_pt[0], trgSf);
            hTrgEvts_mu22ph22_ele_eta->Fill(genele_eta[0], trgSf);
            hTrgEvts_mu22ph22_mu_pt->Fill(genmu_pt[0], trgSf);
            hTrgEvts_mu22ph22_mu_eta->Fill(genmu_eta[0], trgSf);
         }
         if (HLT_Mu40_eta2p1) {
            hTrgEvts_mu40->Fill(emu_mass, trgSf);
            hTrgEvts_mu40_ele_pt->Fill(genele_pt[0], trgSf);
            hTrgEvts_mu40_ele_eta->Fill(genele_eta[0], trgSf);
            hTrgEvts_mu40_mu_pt->Fill(genmu_pt[0], trgSf);
            hTrgEvts_mu40_mu_eta->Fill(genmu_eta[0], trgSf);
         }

         // at least one gsf electron and one muon above the threshold
         //if (gsf_size < 1 || muon_size < 1) continue;

         vector<int> GSF_passHEEP;
         vector<int> MU_passGOOD;
         bool recoElePassAcc = false;
         bool recoMuPassAcc = false;
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
            if (fabs(gsf_eta[j]) < 1.442 || (fabs(gsf_eta[j]) > 1.56 && fabs(gsf_eta[j]) < 2.5)) 
               if (gsf_gsfet[j] > elePtCut) recoElePassAcc = true;

            if (fakeEle) continue;

            if (PassHEEP(j)) GSF_passHEEP.push_back(j);
         }

         //loop over muons
         for (int j = 0; j < muon_size; ++j) {
            if (fabs(muon_eta[j]) < muEtaCut && muon_pt[j] > muPtCut) recoMuPassAcc = true;

            if (PassHighPtMu(j)) MU_passGOOD.push_back(j);
         }

         // reco leptons in acceptance
         if (recoElePassAcc && recoMuPassAcc) hRecoEvtsLepInAcc->Fill(emu_mass);

         // find the e-mu pair with the maximum invariant mass if it exists
         unsigned int eleInd = 0;
         unsigned int muInd = 0;
         double invMass = 0.;
         if (GSF_passHEEP.size() > 0 && MU_passGOOD.size() > 0) {
            TLorentzVector ele1;
            TLorentzVector mu1;

            // find the e-mu pair with the maximum invariant mass
            double maxInvMass = 0.;
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
         }

         if (useScaleFactors && MU_passGOOD.size() > 0) {
            // set trg and muon scale factor according to selected muon
            WeightMuonRecoIsoTrigger(muon_pt[MU_passGOOD[muInd]], muon_eta[MU_passGOOD[muInd]], muSf, trgSf, trgEff, sigString);
         }
         float sf = trgSf*muSf;
         float leptSf = muSf;

         if (GSF_passHEEP.size() > 0) {
            if (useScaleFactors) {
               // set ele scale factor according to selected electron
               if (fabs(GSF_passHEEP[eleInd]) < 0.8) eleSf = eps_cand_sf_0p8;
               else if (fabs(GSF_passHEEP[eleInd]) < 1.4442) eleSf = eps_cand_sf_0p8to1p4442;
               else if (fabs(GSF_passHEEP[eleInd]) > 1.566 && fabs(GSF_passHEEP[eleInd]) < 2.) eleSf = eps_cand_sf_1p566to2p0;
               else if (fabs(GSF_passHEEP[eleInd]) > 2. && fabs(GSF_passHEEP[eleInd]) < 2.5) eleSf = eps_cand_sf_2p0to2p5;
               if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.5) {
                  if (gsfsc_eta[GSF_passHEEP[eleInd]] < 100.) eleSf *= eps_heep_sf_eb_pt35;
                  else eleSf *= eps_heep_sf_eb_pt100;
               } else {
                  if (gsfsc_eta[GSF_passHEEP[eleInd]] < 100.) eleSf *= eps_heep_sf_ee_pt35;
                  else eleSf *= eps_heep_sf_ee_pt100;
               }
               sf *= eleSf;
               leptSf *= eleSf;
            }

            if (triggerBit) {
               hRecoEleEvts->Fill(emu_mass, sf);
               if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.5) hRecoEleEvtsEB->Fill(emu_mass, sf);
               if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.5) hRecoEleEvtsEE->Fill(emu_mass, sf);
            }
            hRecoNoTrgEleEvts->Fill(emu_mass, leptSf);
            if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.5) hRecoNoTrgEleEvtsEB->Fill(emu_mass, leptSf);
            if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.5) hRecoNoTrgEleEvtsEE->Fill(emu_mass, leptSf);
         } 
         if (MU_passGOOD.size() > 0) {
            if (triggerBit) hRecoMuEvts->Fill(emu_mass, sf);
            hRecoNoTrgMuEvts->Fill(emu_mass, leptSf);
         }

         // veto when there are more than one good candidates
         if (GSF_passHEEP.size() < 1 || MU_passGOOD.size() < 1) continue;

         // detector region
         if (eleDetRegion == 1 && fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.5) continue;  // EB
         else if (eleDetRegion == 2 && fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.5) continue;  // EE

         //MASS CUT
         if (invMass < minInvMass) continue;

         if (hltL1sMu16Eta2p1) hFltrMatchEvts_l1Mu16Eta2p1->Fill(emu_mass, sf);
         if (hltL1Mu3p5EG12L3Filtered22) hFltrMatchEvts_mu22ph22_muLeg->Fill(emu_mass, sf);
         if (hltL1sL1SingleEG12) hFltrMatchEvts_l1SingleEG12->Fill(emu_mass, sf);
         if (hltL1sL1Mu3p5EG12) hFltrMatchEvts_l1Mu3p5EG12->Fill(emu_mass, sf);
         if (hltMu22Photon22CaloIdLHEFilter) hFltrMatchEvts_mu22ph22_phLeg->Fill(emu_mass, sf);
         if (hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q) hFltrMatchEvts_mu40->Fill(emu_mass, sf);

         if (triggerBit) {
            hRecoEvts->Fill(emu_mass, sf);
            if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.5) hRecoEvtsEB->Fill(emu_mass, sf);
            if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.5) hRecoEvtsEE->Fill(emu_mass, sf);
         }
         hRecoNoTrgEvts->Fill(emu_mass, leptSf);
         if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) < 1.5) hRecoNoTrgEvtsEB->Fill(emu_mass, leptSf);
         if (fabs(gsfsc_eta[GSF_passHEEP[eleInd]]) > 1.5) hRecoNoTrgEvtsEE->Fill(emu_mass, leptSf);
         ++evCounter;
        ///////////////////////////////////////////////////////////////////////
      } //END LOOP OVER EVENTS
        //////////////////////////////////////////////////////////////////////

     //////////////////////////////////////////////////////////////////////////
   } //END LOOP OVER FILES
     //////////////////////////////////////////////////////////////////////////
   hAcc = (TH1F*)hGenEvtsInAcc->Clone("hAcc");
   hAcc->Divide(hAcc, hGenEvts, 1, 1, "B");
   hAccEle = (TH1F*)hGenEvtsEleInAcc->Clone("hAccEle");
   hAccEle->Divide(hAccEle, hGenEvts, 1, 1, "B");
   hAccEleEB = (TH1F*)hGenEvtsEleInAccEB->Clone("hAccEleEB");
   hAccEleEB->Divide(hAccEleEB, hGenEvts, 1, 1, "B");
   hAccEleEE = (TH1F*)hGenEvtsEleInAccEE->Clone("hAccEleEE");
   hAccEleEE->Divide(hAccEleEE, hGenEvts, 1, 1, "B");
   hAccMu = (TH1F*)hGenEvtsMuInAcc->Clone("hAccMu");
   hAccMu->Divide(hAccMu, hGenEvts, 1, 1, "B");
   hRecoLepAcc = (TH1F*)hRecoEvtsLepInAcc->Clone("hRecoLepAcc");
   hRecoLepAcc->Divide(hRecoLepAcc, hGenEvts, 1, 1, "B");
   hAccTimesEff = (TH1F*)hRecoEvts->Clone("hAccTimesEff");
   hAccTimesEff->Divide(hAccTimesEff, hGenEvts, 1, 1, "B");
   hAccTimesEffEB = (TH1F*)hRecoEvtsEB->Clone("hAccTimesEffEB");
   hAccTimesEffEB->Divide(hAccTimesEffEB, hGenEvts, 1, 1, "B");
   hAccTimesEffEE = (TH1F*)hRecoEvtsEE->Clone("hAccTimesEffEE");
   hAccTimesEffEE->Divide(hAccTimesEffEE, hGenEvts, 1, 1, "B");
   hAccTimesEffEle = (TH1F*)hRecoEleEvts->Clone("hAccTimesEffEle");
   hAccTimesEffEle->Divide(hAccTimesEffEle, hGenEvts, 1, 1, "B");
   hAccTimesEffEleEB = (TH1F*)hRecoEleEvtsEB->Clone("hAccTimesEffEleEB");
   hAccTimesEffEleEB->Divide(hAccTimesEffEleEB, hGenEvts, 1, 1, "B");
   hAccTimesEffEleEE = (TH1F*)hRecoEleEvtsEE->Clone("hAccTimesEffEleEE");
   hAccTimesEffEleEE->Divide(hAccTimesEffEleEE, hGenEvts, 1, 1, "B");
   hAccTimesEffMu = (TH1F*)hRecoMuEvts->Clone("hAccTimesEffMu");
   hAccTimesEffMu->Divide(hAccTimesEffMu, hGenEvts, 1, 1, "B");
   hAccTimesEffNoTrg = (TH1F*)hRecoNoTrgEvts->Clone("hAccTimesEffNoTrg");
   hAccTimesEffNoTrg->Divide(hAccTimesEffNoTrg, hGenEvts, 1, 1, "B");
   //hAccTimesEffNoTrgEB = (TH1F*)hRecoNoTrgEvtsEB->Clone("hAccTimesEffNoTrgEB");
   //hAccTimesEffNoTrgEB->Divide(hAccTimesEffNoTrgEB, hGenEvts, 1, 1, "B");
   //hAccTimesEffNoTrgEE = (TH1F*)hRecoNoTrgEvtsEE->Clone("hAccTimesEffNoTrgEE");
   //hAccTimesEffNoTrgEE->Divide(hAccTimesEffNoTrgEE, hGenEvts, 1, 1, "B");
   hAccTimesEffNoTrgEle = (TH1F*)hRecoNoTrgEleEvts->Clone("hAccTimesEffNoTrgEle");
   hAccTimesEffNoTrgEle->Divide(hAccTimesEffNoTrgEle, hGenEvts, 1, 1, "B");
   hAccTimesEffNoTrgEleEB = (TH1F*)hRecoNoTrgEleEvtsEB->Clone("hAccTimesEffNoTrgEleEB");
   hAccTimesEffNoTrgEleEB->Divide(hAccTimesEffNoTrgEleEB, hGenEvts, 1, 1, "B");
   hAccTimesEffNoTrgEleEE = (TH1F*)hRecoNoTrgEleEvtsEE->Clone("hAccTimesEffNoTrgEleEE");
   hAccTimesEffNoTrgEleEE->Divide(hAccTimesEffNoTrgEleEE, hGenEvts, 1, 1, "B");
   hAccTimesEffNoTrgMu = (TH1F*)hRecoNoTrgMuEvts->Clone("hAccTimesEffNoTrgMu");
   hAccTimesEffNoTrgMu->Divide(hAccTimesEffNoTrgMu, hGenEvts, 1, 1, "B");
   hEffAftTrg = (TH1F*)hRecoEvts->Clone("hEffAftTrg");
   hEffAftTrg->Divide(hEffAftTrg, hTrgEvts, 1, 1, "B");
   hEffAftTrgEle = (TH1F*)hRecoEleEvts->Clone("hEffAftTrgEle");
   hEffAftTrgEle->Divide(hEffAftTrgEle, hTrgEvts, 1, 1, "B");
   hEffAftTrgEleEB = (TH1F*)hRecoEleEvtsEB->Clone("hEffAftTrgEleEB");
   hEffAftTrgEleEB->Divide(hEffAftTrgEleEB, hTrgEvts, 1, 1, "B");
   hEffAftTrgEleEE = (TH1F*)hRecoEleEvtsEE->Clone("hEffAftTrgEleEE");
   hEffAftTrgEleEE->Divide(hEffAftTrgEleEE, hTrgEvts, 1, 1, "B");
   hEffAftTrgMu = (TH1F*)hRecoMuEvts->Clone("hEffAftTrgMu");
   hEffAftTrgMu->Divide(hEffAftTrgMu, hTrgEvts, 1, 1, "B");
   hTrgRecoVsReco = (TH1F*)hRecoEvts->Clone("hTrgRecoVsReco");
   hTrgRecoVsReco->Divide(hTrgRecoVsReco, hRecoNoTrgEvts, 1, 1, "B");
   //hTrgRecoVsRecoEB = (TH1F*)hRecoEvtsEB->Clone("hTrgRecoVsRecoEB");
   //hTrgRecoVsRecoEB->Divide(hTrgRecoVsRecoEB, hRecoNoTrgEvtsEB, 1, 1, "B");
   //hTrgRecoVsRecoEE = (TH1F*)hRecoEvtsEE->Clone("hTrgRecoVsRecoEE");
   //hTrgRecoVsRecoEE->Divide(hTrgRecoVsRecoEE, hRecoNoTrgEvtsEE, 1, 1, "B");
   hTrgRecoVsRecoEle = (TH1F*)hRecoEleEvts->Clone("hTrgRecoVsRecoEle");
   hTrgRecoVsRecoEle->Divide(hTrgRecoVsRecoEle, hRecoNoTrgEleEvts, 1, 1, "B");
   hTrgRecoVsRecoEleEB = (TH1F*)hRecoEleEvtsEB->Clone("hTrgRecoVsRecoEleEB");
   hTrgRecoVsRecoEleEB->Divide(hTrgRecoVsRecoEleEB, hRecoNoTrgEleEvtsEB, 1, 1, "B");
   hTrgRecoVsRecoEleEE = (TH1F*)hRecoEleEvtsEE->Clone("hTrgRecoVsRecoEleEE");
   hTrgRecoVsRecoEleEE->Divide(hTrgRecoVsRecoEleEE, hRecoNoTrgEleEvtsEE, 1, 1, "B");
   hTrgRecoVsRecoMu = (TH1F*)hRecoMuEvts->Clone("hTrgRecoVsRecoMu");
   hTrgRecoVsRecoMu->Divide(hTrgRecoVsRecoMu, hRecoNoTrgMuEvts, 1, 1, "B");
   hFltrMatchVsReco_l1SingleEG12 = (TH1F*)hFltrMatchEvts_l1SingleEG12->Clone("hFltrMatchVsReco_l1SingleEG12");
   hFltrMatchVsReco_l1SingleEG12->Divide(hFltrMatchVsReco_l1SingleEG12, hRecoNoTrgEvts, 1, 1, "B");
   hFltrMatchVsReco_l1Mu3p5EG12 = (TH1F*)hFltrMatchEvts_l1Mu3p5EG12->Clone("hFltrMatchVsReco_l1Mu3p5EG12");
   hFltrMatchVsReco_l1Mu3p5EG12->Divide(hFltrMatchVsReco_l1Mu3p5EG12, hRecoNoTrgEvts, 1, 1, "B");
   hFltrMatchVsReco_l1Mu16Eta2p1 = (TH1F*)hFltrMatchEvts_l1Mu16Eta2p1->Clone("hFltrMatchVsReco_l1Mu16Eta2p1");
   hFltrMatchVsReco_l1Mu16Eta2p1->Divide(hFltrMatchVsReco_l1Mu16Eta2p1, hRecoNoTrgEvts, 1, 1, "B");
   hFltrMatchVsReco_mu22ph22_phLeg = (TH1F*)hFltrMatchEvts_mu22ph22_phLeg->Clone("hFltrMatchVsReco_mu22ph22_phLeg");
   hFltrMatchVsReco_mu22ph22_phLeg->Divide(hFltrMatchVsReco_mu22ph22_phLeg, hRecoNoTrgEvts, 1, 1, "B");
   hFltrMatchVsReco_mu22ph22_muLeg = (TH1F*)hFltrMatchEvts_mu22ph22_muLeg->Clone("hFltrMatchVsReco_mu22ph22_muLeg");
   hFltrMatchVsReco_mu22ph22_muLeg->Divide(hFltrMatchVsReco_mu22ph22_muLeg, hRecoNoTrgEvts, 1, 1, "B");
   hFltrMatchVsPrev_mu22ph22_phLeg = (TH1F*)hFltrMatchEvts_mu22ph22_phLeg->Clone("hFltrMatchVsPrev_mu22ph22_phLeg");
   hFltrMatchVsPrev_mu22ph22_phLeg->Divide(hFltrMatchVsPrev_mu22ph22_phLeg, hFltrMatchEvts_mu22ph22_muLeg, 1, 1, "B");
   hFltrMatchVsPrev_mu22ph22_muLeg = (TH1F*)hFltrMatchEvts_mu22ph22_muLeg->Clone("hFltrMatchVsPrev_mu22ph22_muLeg");
   hFltrMatchVsPrev_mu22ph22_muLeg->Divide(hFltrMatchVsPrev_mu22ph22_muLeg, hFltrMatchEvts_l1Mu3p5EG12, 1, 1, "B");
   hFltrMatchVsReco_mu40 = (TH1F*)hFltrMatchEvts_mu40->Clone("hFltrMatchVsReco_mu40");
   hFltrMatchVsReco_mu40->Divide(hFltrMatchVsReco_mu40, hRecoNoTrgEvts, 1, 1, "B");
   hFltrMatchVsPrev_mu40 = (TH1F*)hFltrMatchEvts_mu40->Clone("hFltrMatchVsPrev_mu40");
   hFltrMatchVsPrev_mu40->Divide(hFltrMatchVsPrev_mu40, hFltrMatchEvts_l1Mu16Eta2p1, 1, 1, "B");
   hTrgEff = (TH1F*)hTrgEvts->Clone("hTrgEff");
   hTrgEff->Divide(hTrgEff, hGenEvts, 1, 1, "B");
   hTrgEff_mu22ph22 = (TH1F*)hTrgEvts_mu22ph22->Clone("hTrgEff_mu22ph22");
   hTrgEff_mu22ph22->Divide(hTrgEff_mu22ph22, hGenEvts, 1, 1, "B");
   hTrgEff_mu22ph22_ele_pt = (TH1F*)hTrgEvts_mu22ph22_ele_pt->Clone("hTrgEff_mu22ph22_ele_pt");
   hTrgEff_mu22ph22_ele_pt->Divide(hTrgEff_mu22ph22_ele_pt, hGenEvts_ele_pt, 1, 1, "B");
   hTrgEff_mu22ph22_ele_eta = (TH1F*)hTrgEvts_mu22ph22_ele_eta->Clone("hTrgEff_mu22ph22_ele_eta");
   hTrgEff_mu22ph22_ele_eta->Divide(hTrgEff_mu22ph22_ele_eta, hGenEvts_ele_eta, 1, 1, "B");
   hTrgEff_mu22ph22_mu_pt = (TH1F*)hTrgEvts_mu22ph22_mu_pt->Clone("hTrgEff_mu22ph22_mu_pt");
   hTrgEff_mu22ph22_mu_pt->Divide(hTrgEff_mu22ph22_mu_pt, hGenEvts_mu_pt, 1, 1, "B");
   hTrgEff_mu22ph22_mu_eta = (TH1F*)hTrgEvts_mu22ph22_mu_eta->Clone("hTrgEff_mu22ph22_mu_eta");
   hTrgEff_mu22ph22_mu_eta->Divide(hTrgEff_mu22ph22_mu_eta, hGenEvts_mu_eta, 1, 1, "B");
   hTrgEff_mu40 = (TH1F*)hTrgEvts_mu40->Clone("hTrgEff_mu40");
   hTrgEff_mu40->Divide(hTrgEff_mu40, hGenEvts, 1, 1, "B");
   hTrgEff_mu40_ele_pt = (TH1F*)hTrgEvts_mu40_ele_pt->Clone("hTrgEff_mu40_ele_pt");
   hTrgEff_mu40_ele_pt->Divide(hTrgEff_mu40_ele_pt, hGenEvts_ele_pt, 1, 1, "B");
   hTrgEff_mu40_ele_eta = (TH1F*)hTrgEvts_mu40_ele_eta->Clone("hTrgEff_mu40_ele_eta");
   hTrgEff_mu40_ele_eta->Divide(hTrgEff_mu40_ele_eta, hGenEvts_ele_eta, 1, 1, "B");
   hTrgEff_mu40_mu_pt = (TH1F*)hTrgEvts_mu40_mu_pt->Clone("hTrgEff_mu40_mu_pt");
   hTrgEff_mu40_mu_pt->Divide(hTrgEff_mu40_mu_pt, hGenEvts_mu_pt, 1, 1, "B");
   hTrgEff_mu40_mu_eta = (TH1F*)hTrgEvts_mu40_mu_eta->Clone("hTrgEff_mu40_mu_eta");
   hTrgEff_mu40_mu_eta->Divide(hTrgEff_mu40_mu_eta, hGenEvts_mu_eta, 1, 1, "B");
   hPassPtEff = (TH1F*)hPassPtCuts->Clone("hPassPtEff");
   hPassPtEff->Divide(hPassPtEff, hGenEvts, 1, 1, "B");
   hPassEtaEff = (TH1F*)hPassEtaCuts->Clone("hPassEtaEff");
   hPassEtaEff->Divide(hPassEtaEff, hGenEvts, 1, 1, "B");
   hPassPtEtaEff = (TH1F*)hPassPtEtaCuts->Clone("hPassPtEtaEff");
   hPassPtEtaEff->Divide(hPassPtEtaEff, hGenEvts, 1, 1, "B");

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

   TH1F* hAccTimesEffNote = (TH1F*)hAccTimesEff->Clone("hAccTimesEffNote");
   TH1F* hAccTimesEff2 = (TH1F*)hAccTimesEff->Clone("hAccTimesEff2");

   TF1 *andreasFunc = new TF1("andreasFunc", "[0] + [1]/ (x + [2]) + [3]*x", 200., 5010.);
   andreasFunc->SetParameters(0.78, -74.4, 46.7, -3.5e-5);
   andreasFunc->SetLineColor(kGreen);
   TF1 *otherDatasetFunc = new TF1("otherDatasetFunc", "[0] + [1]/ (x + [2]) + [3]*x", 200., 5010.);
   if (triggerInd == 0)
      otherDatasetFunc->SetParameters(7.41590e-1, -1.41313e2, 1.63499e2, -2.87806e-5);
   else
      otherDatasetFunc->SetParameters(7.50730e-1, -1.15163e2, 1.21807e2, -2.58317e-5);
   otherDatasetFunc->SetLineColor(kRed);
   TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]/ (x + [2]) + [3]*x", 200., 5010.);
   TF1 *fitFunc2 = new TF1("fitFunc2", "[0] + [1]/ (x + [2]) + [3]/ (x*x + [4])", 200., 5010.);
   //TF1 *fitFuncEB = new TF1("fitFuncEB", "[0] + [1]/ (x + [2])", 10., 5010.);
   TF1 *fitFuncEB = new TF1("fitFuncEB", "[0] + [1]/ (x + [2]) + [3]*x", 10., 5010.);
   //TF1 *fitFuncEE = new TF1("fitFuncEE", "[0] + [1]/ (x*x + [2])", 10., 5010.);
   TF1 *fitFuncEE = new TF1("fitFuncEE", "[0] + [1]/ (x + [2]) + [3]/ (x*x + [4])", 10., 5010.);
   fitFunc->SetLineColor(kBlue);
   fitFunc2->SetLineColor(kRed);
   fitFuncEB->SetLineColor(kBlue);
   fitFuncEE->SetLineColor(kBlue);
   hAccTimesEff->Fit("fitFunc", "", "", 50., 2600.);
   //hAccTimesEff->Fit("fitFunc2", "", "", 50., 2600.);
   hAccTimesEffEB->Fit("fitFuncEB", "", "", 50., 2600.);
   hAccTimesEffEE->Fit("fitFuncEE", "", "", 50., 2600.);
   cout << "Chi^2 / NDF F1: " << fitFunc->GetChisquare() << " / " << fitFunc->GetNDF() << ", prob: " << fitFunc->GetProb() << endl;
   cout << "Chi^2 / NDF F2: " << fitFunc2->GetChisquare() << " / " << fitFunc2->GetNDF() << ", prob: " << fitFunc2->GetProb() << endl;
   cout << "Chi^2 / NDF EB: " << fitFuncEB->GetChisquare() << " / " << fitFuncEB->GetNDF() << ", prob: " << fitFuncEB->GetProb() << endl;
   cout << "Chi^2 / NDF EE: " << fitFuncEE->GetChisquare() << " / " << fitFuncEE->GetNDF() << ", prob: " << fitFuncEE->GetProb() << endl;

   hAccTimesEff->GetYaxis()->SetTitle("acc x eff");
   hAccTimesEff->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEff->Draw();
   fitFunc->Draw("same");
   //fitFunc2->Draw("same");
   andreasFunc->Draw("same");
   TLatex *tex = new TLatex(0.15, 0.25, "P(M|p0,p1,p2,p3) = p0 + #frac{p1}{M+p2} + p3*M");
   tex->SetNDC();
   tex->SetTextFont(font);
   tex->SetLineWidth(2);
   tex->SetTextSize(0.03);
   tex->Draw();
   //tex->DrawLatex(0.15, 0.17, "Andreas (green): A#dot#epsilon(M) = 0.61 + #frac{280.1}{M + 2008.7} - #frac{30537.3}{M^{2} + 75925.2}");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   tex->DrawLatex(0.17, 0.85, "trg + electron + muon");
   tex->DrawLatex(0.17, 0.80, "trg: " + triggerName);

   // plot
   TCanvas *accTimesEffPlotNote = new TCanvas("accTimesEffPlotNote", "Ax#epsilon", 100, 100, 600, 600);
   TPad *accTimesEffPadNote = (TPad*)accTimesEffPad->Clone("accTimesEffPadNote");
   accTimesEffPadNote->Draw();
   accTimesEffPadNote->cd();

   gStyle->SetOptFit(0);
   gPad->SetTicks(1, 1);
   gPad->SetGrid(0, 0);

   TH1F* hAccNote = (TH1F*)hAcc->Clone("hAccNote");
   TH1F* hAccTimesEffNoTrgNote = (TH1F*)hAccTimesEffNoTrg->Clone("hAccTimesEffNoTrgNote");
   TF1* fitFuncNote = (TF1 *)fitFunc->Clone("fitFuncNote");

   cout << "Note plot:" << endl;
   hAccTimesEffNote->GetYaxis()->SetTitle("Ax#epsilon");
   hAccTimesEffNote->GetYaxis()->SetRangeUser(0., 1.);
   hAccTimesEffNote->GetXaxis()->SetRangeUser(0., 3500.);
   hAccTimesEffNote->SetLineColor(kWhite);
   hAccTimesEffNote->Draw();
   hAccTimesEffNote->Fit("fitFuncNote", "", "", 50., 2600.);
   hAccNote->SetMarkerColor(kBlue);
   hAccNote->Draw("histpsame");
   hAccTimesEffNoTrgNote->SetMarkerColor(kRed);
   hAccTimesEffNoTrgNote->Draw("histpsame");
   //andreasFunc->Draw("same");
   TLegend* legendNote = new TLegend(0.223, 0.120, 0.910, 0.521);
   legendNote->SetTextFont(font);
   legendNote->SetTextSize(0.035);
   legendNote->SetBorderSize(1);
   legendNote->SetLineColor(1);
   legendNote->SetLineStyle(1);
   legendNote->SetLineWidth(1);
   legendNote->SetFillColor(0);
   //legendNote->SetFillStyle(0);
   legendNote->AddEntry(hAccNote, "gen leptons within acceptance", "p");
   legendNote->AddEntry(hAccTimesEffNoTrgNote, "reco leptons within acceptance", "p");
   legendNote->AddEntry(hAccTimesEffNote, "full selection with trigger", "p");
   legendNote->AddEntry(fitFuncNote, "#splitline{fit 200 GeV < M_{Z'}^{truth} < 3.6 TeV}{A+B/(M+C)+DxM}", "l");
   legendNote->Draw("same");
   //tex->DrawLatex(0.15, 0.17, "Andreas (green): A#dot#epsilon(M) = 0.61 + #frac{280.1}{M + 2008.7} - #frac{30537.3}{M^{2} + 75925.2}");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   //tex->DrawLatex(0.17, 0.85, "trg + electron + muon");
   //tex->DrawLatex(0.17, 0.80, "trg: " + triggerName);

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
   //tex->DrawLatex(0.45, 0.38, "P(M|p0,p1,p2) = p0 + #frac{p1}{M^{2}+p2}");
   tex->DrawLatex(0.35, 0.55, "P(M|p0,p1,p2,p3,p4) = p0 + #frac{p1}{M+p2} + #frac{p3}{M^{2}+p4}");
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
      hFltrMatchVsReco_l1SingleEG12->SetMarkerStyle(kOpenCircle);
      hFltrMatchVsReco_l1SingleEG12->SetMarkerColor(kCyan);
      hFltrMatchVsReco_l1SingleEG12->SetLineColor(kCyan);
      hFltrMatchVsReco_l1SingleEG12->Draw("same");
      hFltrMatchVsReco_l1Mu3p5EG12->SetMarkerStyle(kOpenSquare);
      hFltrMatchVsReco_l1Mu3p5EG12->SetMarkerColor(kMagenta);
      hFltrMatchVsReco_l1Mu3p5EG12->SetLineColor(kMagenta);
      hFltrMatchVsReco_l1Mu3p5EG12->Draw("same");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_mu22ph22_muLeg, "eff #mu-leg HLT_Mu22_Photon22_CaloIdL");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_mu22ph22_phLeg, "eff #gamma-leg HLT_Mu22_Photon22_CaloIdL");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_l1SingleEG12, "eff L1 L1sL1SingleEG12");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_l1Mu3p5EG12, "eff L1 HLT_Mu22_Photon22_CaloIdL");
   } else {
      hFltrMatchVsReco_mu40->GetYaxis()->SetTitle("eff");
      hFltrMatchVsReco_mu40->GetYaxis()->SetRangeUser(0.75, 1.);
      hFltrMatchVsReco_mu40->SetMarkerStyle(21);
      hFltrMatchVsReco_mu40->SetMarkerColor(kRed);
      hFltrMatchVsReco_mu40->SetLineColor(kRed);
      hFltrMatchVsReco_mu40->Draw();
      hFltrMatchVsReco_l1Mu16Eta2p1->SetMarkerStyle(kOpenSquare);
      hFltrMatchVsReco_l1Mu16Eta2p1->SetMarkerColor(kMagenta);
      hFltrMatchVsReco_l1Mu16Eta2p1->SetLineColor(kMagenta);
      hFltrMatchVsReco_l1Mu16Eta2p1->Draw("same");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_mu40, "eff #mu-leg HLT_Mu40_eta2p1");
      trgFltrLegend->AddEntry(hFltrMatchVsReco_l1Mu16Eta2p1, "eff L1 HLT_Mu40_eta2p1");
   }
   trgFltrLegend->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   if (eleDetRegion == 1) tex->DrawLatex(0.60, 0.935, "barrel electron");
   else if (eleDetRegion == 2) tex->DrawLatex(0.60, 0.935, "endcap electron");

   TCanvas *trgFltrEffPlotByFilter = new TCanvas("trgFltrEffPlotByFilter", "trigger filter eff wrt previous filter", 100, 100, 600, 600);
   TPad *trgFltrEffPadByFilter = (TPad*)accTimesEffPad->Clone("trgFltrEffPadByFilter");
   trgFltrEffPadByFilter->Draw(); 
   trgFltrEffPadByFilter->cd();
   TLegend* trgFltrLegendByFilter = (TLegend*)legend->Clone("trgFltrLegendByFilter");
   trgFltrLegendByFilter->Clear();
   trgFltrLegendByFilter->SetX1(0.33);
   trgFltrLegendByFilter->SetX2(0.62);
   if (triggerInd == 0) {
      hFltrMatchVsPrev_mu22ph22_muLeg->GetYaxis()->SetTitle("filter eff");
      hFltrMatchVsPrev_mu22ph22_muLeg->GetYaxis()->SetRangeUser(0.75, 1.);
      hFltrMatchVsPrev_mu22ph22_muLeg->SetMarkerStyle(kFullTriangleUp);
      hFltrMatchVsPrev_mu22ph22_muLeg->SetMarkerColor(kGreen+1);
      hFltrMatchVsPrev_mu22ph22_muLeg->SetLineColor(kGreen+1);
      hFltrMatchVsPrev_mu22ph22_muLeg->Draw();
      hFltrMatchVsPrev_mu22ph22_phLeg->SetMarkerStyle(kFullTriangleDown);
      hFltrMatchVsPrev_mu22ph22_phLeg->SetMarkerColor(kBlue);
      hFltrMatchVsPrev_mu22ph22_phLeg->SetLineColor(kBlue);
      hFltrMatchVsPrev_mu22ph22_phLeg->Draw("same");
      //hFltrMatchVsReco_l1SingleEG12->SetMarkerStyle(kOpenCircle);
      //hFltrMatchVsReco_l1SingleEG12->SetMarkerColor(kCyan);
      //hFltrMatchVsReco_l1SingleEG12->SetLineColor(kCyan);
      //hFltrMatchVsReco_l1SingleEG12->Draw("same");
      hFltrMatchVsReco_l1Mu3p5EG12->SetMarkerStyle(kOpenSquare);
      hFltrMatchVsReco_l1Mu3p5EG12->SetMarkerColor(kMagenta);
      hFltrMatchVsReco_l1Mu3p5EG12->SetLineColor(kMagenta);
      hFltrMatchVsReco_l1Mu3p5EG12->Draw("same");
      trgFltrLegendByFilter->AddEntry(hFltrMatchVsPrev_mu22ph22_muLeg, "eff #mu-leg HLT_Mu22_Photon22_CaloIdL");
      trgFltrLegendByFilter->AddEntry(hFltrMatchVsPrev_mu22ph22_phLeg, "eff #gamma-leg HLT_Mu22_Photon22_CaloIdL");
      //trgFltrLegendByFilter->AddEntry(hFltrMatchVsReco_l1SingleEG12, "eff L1 L1sL1SingleEG12");
      trgFltrLegendByFilter->AddEntry(hFltrMatchVsReco_l1Mu3p5EG12, "eff L1 HLT_Mu22_Photon22_CaloIdL");
   } else {
      hFltrMatchVsPrev_mu40->GetYaxis()->SetTitle("filter eff");
      hFltrMatchVsPrev_mu40->GetYaxis()->SetRangeUser(0.75, 1.);
      hFltrMatchVsPrev_mu40->SetMarkerStyle(21);
      hFltrMatchVsPrev_mu40->SetMarkerColor(kRed);
      hFltrMatchVsPrev_mu40->SetLineColor(kRed);
      hFltrMatchVsPrev_mu40->Draw();
      hFltrMatchVsReco_l1Mu16Eta2p1->SetMarkerStyle(kOpenSquare);
      hFltrMatchVsReco_l1Mu16Eta2p1->SetMarkerColor(kMagenta);
      hFltrMatchVsReco_l1Mu16Eta2p1->SetLineColor(kMagenta);
      hFltrMatchVsReco_l1Mu16Eta2p1->Draw("same");
      trgFltrLegendByFilter->AddEntry(hFltrMatchVsReco_mu40, "eff #mu-leg HLT_Mu40_eta2p1");
      trgFltrLegendByFilter->AddEntry(hFltrMatchVsReco_l1Mu16Eta2p1, "eff L1 HLT_Mu40_eta2p1");
   }
   trgFltrLegendByFilter->Draw("same");
   tex->DrawLatex(0.109, 0.935, "CMS Simulation, 8 TeV");
   if (eleDetRegion == 1) tex->DrawLatex(0.60, 0.935, "barrel electron");
   else if (eleDetRegion == 2) tex->DrawLatex(0.60, 0.935, "endcap electron");

   // safe in various file formats
   if (saveSpec) {
     if (saveAsPdf) {
        if (eleDetRegion == 0) {
           trgEffPlot->Print(plotDir + trgEffPlot->GetName() + fileNameExtra + ".pdf", "pdf");
           trgEffPlot_ele_pt->Print(plotDir + trgEffPlot_ele_pt->GetName() + fileNameExtra + ".pdf", "pdf");
           trgEffPlot_ele_eta->Print(plotDir + trgEffPlot_ele_eta->GetName() + fileNameExtra + ".pdf", "pdf");
           trgEffPlot_mu_pt->Print(plotDir + trgEffPlot_mu_pt->GetName() + fileNameExtra + ".pdf", "pdf");
           trgEffPlot_mu_eta->Print(plotDir + trgEffPlot_mu_eta->GetName() + fileNameExtra + ".pdf", "pdf");
           accTimesEffPlotNote->Print(plotDir + accTimesEffPlotNote->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           accTimesEffPlot->Print(plotDir + accTimesEffPlot->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           accTimesEffPlotEB->Print(plotDir + accTimesEffPlotEB->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           accTimesEffPlotEE->Print(plotDir + accTimesEffPlotEE->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           accTimesEffObjPlot->Print(plotDir + accTimesEffObjPlot->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           accTimesEffNoTrgObjPlot->Print(plotDir + accTimesEffNoTrgObjPlot->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           effAftTrgPlot->Print(plotDir + effAftTrgPlot->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           trgRecoVsRecoPlot->Print(plotDir + trgRecoVsRecoPlot->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
           accPlot->Print(plotDir + accPlot->GetName() + "_" + triggerName + fileNameExtra + ".pdf", "pdf");
        }
        trgFltrEffPlot->Print(plotDir + trgFltrEffPlot->GetName() + "_" + triggerName + eleDetRegSuffix + fileNameExtra + ".pdf", "pdf");
        trgFltrEffPlotByFilter->Print(plotDir + trgFltrEffPlotByFilter->GetName() + "_" + triggerName + eleDetRegSuffix + fileNameExtra + ".pdf", "pdf");
     }
     if (saveAsPng) {
        if (eleDetRegion == 0) {
           trgEffPlot->Print(plotDir + trgEffPlot->GetName() + fileNameExtra + ".png", "png");
           trgEffPlot_ele_pt->Print(plotDir + trgEffPlot_ele_pt->GetName() + fileNameExtra + ".png", "png");
           trgEffPlot_ele_eta->Print(plotDir + trgEffPlot_ele_eta->GetName() + fileNameExtra + ".png", "png");
           trgEffPlot_mu_pt->Print(plotDir + trgEffPlot_mu_pt->GetName() + fileNameExtra + ".png", "png");
           trgEffPlot_mu_eta->Print(plotDir + trgEffPlot_mu_eta->GetName() + fileNameExtra + ".png", "png");
           accTimesEffPlotNote->Print(plotDir + accTimesEffPlotNote->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           accTimesEffPlot->Print(plotDir + accTimesEffPlot->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           accTimesEffPlotEB->Print(plotDir + accTimesEffPlotEB->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           accTimesEffPlotEE->Print(plotDir + accTimesEffPlotEE->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           accTimesEffObjPlot->Print(plotDir + accTimesEffObjPlot->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           accTimesEffNoTrgObjPlot->Print(plotDir + accTimesEffNoTrgObjPlot->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           effAftTrgPlot->Print(plotDir + effAftTrgPlot->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           trgRecoVsRecoPlot->Print(plotDir + trgRecoVsRecoPlot->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
           accPlot->Print(plotDir + accPlot->GetName() + "_" + triggerName + fileNameExtra + ".png", "png");
        }
        trgFltrEffPlot->Print(plotDir + trgFltrEffPlot->GetName() + "_" + triggerName + eleDetRegSuffix + fileNameExtra + ".png", "png");
        trgFltrEffPlotByFilter->Print(plotDir + trgFltrEffPlotByFilter->GetName() + "_" + triggerName + eleDetRegSuffix + fileNameExtra + ".png", "png");
     }
     if (saveAsRoot) {
        if (eleDetRegion == 0) {
           trgEffPlot->Print(plotDir + trgEffPlot->GetName() + fileNameExtra + ".root", "root");
           trgEffPlot_ele_pt->Print(plotDir + trgEffPlot_ele_pt->GetName() + fileNameExtra + ".root", "root");
           trgEffPlot_ele_eta->Print(plotDir + trgEffPlot_ele_eta->GetName() + fileNameExtra + ".root", "root");
           trgEffPlot_mu_pt->Print(plotDir + trgEffPlot_mu_pt->GetName() + fileNameExtra + ".root", "root");
           trgEffPlot_mu_eta->Print(plotDir + trgEffPlot_mu_eta->GetName() + fileNameExtra + ".root", "root");
           accTimesEffPlotNote->Print(plotDir + accTimesEffPlotNote->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           accTimesEffPlot->Print(plotDir + accTimesEffPlot->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           accTimesEffPlotEB->Print(plotDir + accTimesEffPlotEB->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           accTimesEffPlotEE->Print(plotDir + accTimesEffPlotEE->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           accTimesEffObjPlot->Print(plotDir + accTimesEffObjPlot->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           accTimesEffNoTrgObjPlot->Print(plotDir + accTimesEffNoTrgObjPlot->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           effAftTrgPlot->Print(plotDir + effAftTrgPlot->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           trgRecoVsRecoPlot->Print(plotDir + trgRecoVsRecoPlot->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
           accPlot->Print(plotDir + accPlot->GetName() + "_" + triggerName + fileNameExtra + ".root", "root");
        }
        trgFltrEffPlot->Print(plotDir + trgFltrEffPlot->GetName() + "_" + triggerName + eleDetRegSuffix + fileNameExtra + ".root", "root");
        trgFltrEffPlotByFilter->Print(plotDir + trgFltrEffPlotByFilter->GetName() + "_" + triggerName + eleDetRegSuffix + fileNameExtra + ".root", "root");
     }
   }

   // write to file
   output->cd();
   trgEffPlot->Write();
   trgEffPlot_ele_pt->Write();
   trgEffPlot_ele_eta->Write();
   trgEffPlot_mu_pt->Write();
   trgEffPlot_mu_eta->Write();
   trgFltrEffPlot->Write();
   trgFltrEffPlotByFilter->Write();
   accTimesEffPlotNote->Write();
   accTimesEffPlot->Write();
   accTimesEffPlotEB->Write();
   accTimesEffPlotEE->Write();
   accTimesEffObjPlot->Write();
   accTimesEffNoTrgObjPlot->Write();
   effAftTrgPlot->Write();
   trgRecoVsRecoPlot->Write();
   accPlot->Write();
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
   hFltrMatchEvts_l1SingleEG12->Write();
   hFltrMatchEvts_l1Mu3p5EG12->Write();
   hFltrMatchEvts_l1Mu16Eta2p1->Write();
   hFltrMatchEvts_mu22ph22_phLeg->Write();
   hFltrMatchEvts_mu22ph22_muLeg->Write();
   hFltrMatchEvts_mu40->Write();
   hFltrMatchVsReco_l1SingleEG12->Write();
   hFltrMatchVsReco_l1Mu3p5EG12->Write();
   hFltrMatchVsReco_l1Mu16Eta2p1->Write();
   hFltrMatchVsReco_mu22ph22_phLeg->Write();
   hFltrMatchVsReco_mu22ph22_muLeg->Write();
   hFltrMatchVsPrev_mu22ph22_phLeg->Write();
   hFltrMatchVsPrev_mu22ph22_muLeg->Write();
   hFltrMatchVsReco_mu40->Write();
   hFltrMatchVsPrev_mu40->Write();
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

bool AccTimesEffRpv::PassHEEP(const int &n)
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

bool AccTimesEffRpv::PassHighPtMu(const int &n)
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

void 
AccTimesEffRpv::WeightMuonRecoIsoTrigger(float MuonPt, float MuonEta, float &weight_muon_reco, float &weight_trigger, float &eff_trigger, TString &type, float l1_eff)
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

