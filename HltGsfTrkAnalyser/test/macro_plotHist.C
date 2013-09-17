#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>

#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TH1F.h"
#include "TFile.h"
#include <TCanvas.h>
#include <TStyle.h>

void macro_plotHist()
{
   TH1::AddDirectory(kFALSE);

   const bool compToReco = 0;
   const bool matched = 0; // when comparing to reco
   int font = 42;
   const bool ratioPlot = 1;
   const bool saveAsPng = 1;
   TString plotsDir = "./plots_20130917/";
   //TString fileNameExtra = "_+nc2-nc6";
   TString fileNameExtra = ""; 

   vector<TString> files;
   if (compToReco) {
      files.push_back("./recoVsNc6_pt10_dr0p5_passHltDEle33.root");
      files.push_back("./recoVsNc5_pt10_dr0p5_passHltDEle33.root");
      files.push_back("./recoVsNc4_pt10_dr0p5_passHltDEle33.root");
      files.push_back("./recoVsNc3_pt10_dr0p5_passHltDEle33.root");
      files.push_back("./recoVsNc2_pt10_dr0p5_passHltDEle33.root");
   } else {
      files.push_back("./nc6VsNc5_pt10_sharedHits3_passHltDEle33.root");
      files.push_back("./nc6VsNc4_pt10_sharedHits3_passHltDEle33.root");
      files.push_back("./nc6VsNc3_pt10_sharedHits3_passHltDEle33.root");
      files.push_back("./nc6VsNc2_pt10_sharedHits3_passHltDEle33.root");
   }
   //files.push_back("./nc6VsNc5_pt30_dr0p5_passHltDEle33.root");
   //files.push_back("./nc6VsNc4_pt30_dr0p5_passHltDEle33.root");
   //files.push_back("./nc6VsNc3_pt30_dr0p5_passHltDEle33.root");
   //files.push_back("./nc6VsNc2_pt30_dr0p5_passHltDEle33.root");
   //files.push_back("./nc6VsNc5_pt30_dr0p5_passHltDEle33TrkNC5_noPassHltDEle33TrkNC6.root");
   //files.push_back("./nc6VsNc4_pt30_dr0p5_passHltDEle33TrkNC4_noPassHltDEle33TrkNC6.root");
   //files.push_back("./nc6VsNc3_pt30_dr0p5_passHltDEle33TrkNC3_noPassHltDEle33TrkNC6.root");
   //files.push_back("./nc6VsNc2_pt30_dr0p5_passHltDEle33TrkNC2_noPassHltDEle33TrkNC6.root");
   //files.push_back("./nc6VsNc5_pt30_dr0p5_passHltDEle33TrkNC6_noPassHltDEle33TrkNC5.root");
   //files.push_back("./nc6VsNc4_pt30_dr0p5_passHltDEle33TrkNC6_noPassHltDEle33TrkNC4.root");
   //files.push_back("./nc6VsNc3_pt30_dr0p5_passHltDEle33TrkNC6_noPassHltDEle33TrkNC3.root");
   //files.push_back("./nc6VsNc2_pt30_dr0p5_passHltDEle33TrkNC6_noPassHltDEle33TrkNC2.root");
   TString dir("HltGsfTrkAna");

   vector<TString> legEntries;
   if (compToReco) legEntries.push_back("RECO");
   legEntries.push_back("NC6");
   legEntries.push_back("NC5");
   legEntries.push_back("NC4");
   legEntries.push_back("NC3");
   legEntries.push_back("NC2");

   TString hPrefixR = "hR_";
   TString hPrefixT = "hT_";
   if (matched) {
      hPrefixR = "hR_matched_";
      hPrefixT = "hT_matched_";
   }
   vector<TString> histoNames;
   if (!matched) histoNames.push_back("size");
   histoNames.push_back("p");
   histoNames.push_back("pt");
   histoNames.push_back("eta");
   histoNames.push_back("phi");
   histoNames.push_back("ptError");
   histoNames.push_back("dPtOverPt");
   histoNames.push_back("etaError");
   histoNames.push_back("phiError");
   histoNames.push_back("outerP");
   histoNames.push_back("outerPt");
   histoNames.push_back("outerEta");
   histoNames.push_back("outerPhi");
   histoNames.push_back("outerPOverP");
   histoNames.push_back("outerPtOverPt");
   histoNames.push_back("chi2");
   histoNames.push_back("ndof");
   histoNames.push_back("normalizedChi2");
   histoNames.push_back("charge");
   histoNames.push_back("dxy");
   histoNames.push_back("dz");
   histoNames.push_back("dxyError");
   histoNames.push_back("dzError");
   histoNames.push_back("numberOfValidHits");
   histoNames.push_back("numberOfLostHits");
   histoNames.push_back("validFraction");
   if (!matched) {
      histoNames.push_back("dEta_mapSize");
      histoNames.push_back("dEta_mapSizeMatch");
      histoNames.push_back("dEta_dEta");
      histoNames.push_back("dEta_nPassCut");
      histoNames.push_back("dEta_recoEcalCandsSize");
      histoNames.push_back("dEta_recoEcalCandsSizeMatch");
      histoNames.push_back("dEta_dEta_trgdEcalCands");
      histoNames.push_back("dEta_nPassCut_trgdEcalCands");
      histoNames.push_back("dPhi_mapSize");
      histoNames.push_back("dPhi_mapSizeMatch");
      histoNames.push_back("dPhi_dPhi");
      histoNames.push_back("dPhi_nPassCut");
      histoNames.push_back("dPhi_recoEcalCandsSize");
      histoNames.push_back("dPhi_recoEcalCandsSizeMatch");
      histoNames.push_back("dPhi_dPhi_trgdEcalCands");
      histoNames.push_back("dPhi_nPassCut_trgdEcalCands");
   }

   int colours[] = {kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};

   TFile inFile(files[0], "read");
   inFile.cd(dir);
   for (unsigned int i = 0; i < histoNames.size(); ++i) {
   //for (unsigned int i = 0; i < 3; ++i) {
      TH1F* hR = (TH1F *)gDirectory->Get(hPrefixR + histoNames[i]);
      vector<TH1F*> hT;
      vector<TH1F*> hRatio;
      cout << hPrefixR << histoNames[i] << " number of entries: " << hR->Integral();
      for (unsigned int j = 0; j < files.size(); ++j) {
         TFile inFileT(files[j], "read");
         inFileT.cd(dir);
         hT.push_back((TH1F *)gDirectory->Get(hPrefixT + histoNames[i]));
         hRatio.push_back((TH1F *)gDirectory->Get(hPrefixT + histoNames[i]));
         hRatio.back()->Add(hR, -1.);
         hRatio.back()->Divide(hR);
      }
      cout << endl;

      //float normFactor = hT[0]->Integral() / hR->Integral();
      //normFactor = 1.;
      //hR->Scale(normFactor);
      //cout << "Normalisation factor: " << normFactor << endl;

      TCanvas *c0;
      TPad *specPad;
      float fullPadTxtHeight;
      if (ratioPlot) {
         c0 = new TCanvas("c0_" + histoNames[i], "c0_" + histoNames[i], 100, 100, 800, 700);
         specPad = new TPad("specPad" + histoNames[i], "emu Spectrum" + histoNames[i], 0., 0., 1., 1.);
         fullPadTxtHeight = specPad->GetHNDC();
         specPad->SetPad(0., 0.33, 1., 1.);
         specPad->SetBottomMargin(0.06);
      } else {
         c0 = new TCanvas("c0_" + histoNames[i], "c0_" + histoNames[i], 100, 100, 800, 600);
         specPad = new TPad("specPad" + histoNames[i], "emu Spectrum" + histoNames[i], 0., 0., 1., 1.);
         fullPadTxtHeight = specPad->GetHNDC();
         hR->GetXaxis()->SetTitle(hR->GetTitle()); 
         hR->GetXaxis()->SetTitleOffset(1.1);
      }
      specPad->SetBorderMode(0);
      specPad->SetFrameBorderMode(0);
      specPad->SetFillColor(0);
      specPad->SetFrameFillColor(0);
      specPad->SetRightMargin(0.123);
      specPad->SetLogy();
      specPad->Draw();
      specPad->cd();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      float fontScaleTop = fullPadTxtHeight / specPad->GetHNDC();

      hR->GetXaxis()->SetTitleSize(0.03 * fontScaleTop);
      hR->GetXaxis()->SetLabelSize(0.03 * fontScaleTop);
      hR->GetYaxis()->SetTitleSize(0.03 * fontScaleTop);
      hR->GetYaxis()->SetTitleOffset(1.1 / fontScaleTop);
      hR->GetYaxis()->SetLabelSize(0.03 * fontScaleTop);
      hR->SetLineColor(colours[0]);
      hR->SetLineWidth(2);
      hR->SetFillColor(colours[0]);
      hR->Draw("h");

      TLegend *legend = new TLegend(0.88, 0.60, 0.99, 0.90);
      legend->SetTextSize(0.03 * fontScaleTop);
      legend->SetTextFont(font);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      //legend->AddEntry(hR, Form("Ref x%.2f", normFactor), "f");
      legend->AddEntry(hR, legEntries[0], "f");

      for (unsigned int j = 0; j < files.size(); ++j) {
         hT[j]->SetLineColor(colours[j+1]);
         hT[j]->SetLineWidth(2);
         hT[j]->Draw("hsame");
         legend->AddEntry(hT[j], legEntries[j+1], "l");
      }

      legend->Draw("same");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.035);
      if (compToReco) tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref: electronGsfTracks - RECO");
      else tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref Trigger: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7 (NC6)");

      float fontScaleBot = 1.;
      TPad *pullPad = new TPad("pullPad" + histoNames[i], "(data - bg) / bg" + histoNames[i], 0., 0., 1., 1.);
      if (ratioPlot) {
         c0->cd();
         pullPad->SetPad(0., 0., 1., 0.33);
         pullPad->SetBottomMargin(0.22);
         pullPad->SetTopMargin(0.);
         fontScaleBot = fullPadTxtHeight / pullPad->GetHNDC();
         pullPad->Draw();
         pullPad->cd();
         pullPad->SetBorderMode(0);
         pullPad->SetBorderSize(2);
         pullPad->SetFrameBorderMode(0);
         pullPad->SetFillColor(0);
         pullPad->SetFrameFillColor(0);
         pullPad->SetRightMargin(0.123);
         pullPad->SetTickx(1);
         pullPad->SetTicky(1);
         //specPad->SetFillColor(0);
         //specPad->SetFrameFillColor(0);
         specPad->SetBottomMargin(0.06);

         for (unsigned int j = 0; j < files.size(); ++j) {
            hRatio[j]->SetLineColor(colours[j+1]);
            hRatio[j]->SetLineWidth(2);
            hRatio[j]->GetXaxis()->SetTitle(hR->GetTitle());
            hRatio[j]->GetXaxis()->SetTitleFont(font);
            hRatio[j]->GetXaxis()->SetTitleSize(0.03 * fontScaleBot);
            hRatio[j]->GetXaxis()->SetTitleOffset(1.1);
            hRatio[j]->GetXaxis()->SetLabelFont(font);
            hRatio[j]->GetXaxis()->SetLabelSize(0.03 * fontScaleBot);
            hRatio[j]->GetYaxis()->SetTitle("(test-ref)/ref");
            hRatio[j]->GetYaxis()->SetTitleFont(font);
            hRatio[j]->GetYaxis()->SetTitleSize(0.03 * fontScaleBot);
            hRatio[j]->GetYaxis()->SetTitleOffset(1.2 / fontScaleBot);
            hRatio[j]->GetYaxis()->SetLabelFont(font);
            hRatio[j]->GetYaxis()->SetLabelSize(0.03 * fontScaleBot);

            if (j == 0) hRatio[j]->Draw("h");
            else hRatio[j]->Draw("hsame");
         }
      }

      if (saveAsPng) {
         if (compToReco) {
            if (matched) c0->Print(plotsDir + "h_vsReco_drMatched0p5_" + histoNames[i] + fileNameExtra + ".png", "png");
            else c0->Print(plotsDir + "h_vsReco_" + histoNames[i] + fileNameExtra + ".png", "png");
         }
         else c0->Print(plotsDir + "h_" + histoNames[i] + fileNameExtra + ".png", "png");
      }
   }
}
