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

void macro_plotPullHist()
{
   TH1::AddDirectory(kFALSE);

   const bool compToReco = 1;
   const bool underflow = 1;
   const bool overflow = 1;
   int font = 42;
   const bool saveAsPng = 0;
   TString plotsDir = "./plots_20130930/";

   vector<TString> files;
   if (compToReco) {
      files.push_back("./recoVsNc6_pt10_dr0p5_passHltDEle33_narrowPull.root");
      files.push_back("./recoVsNc5_pt10_dr0p5_passHltDEle33_narrowPull.root");
      files.push_back("./recoVsNc4_pt10_dr0p5_passHltDEle33_narrowPull.root");
      files.push_back("./recoVsNc3_pt10_dr0p5_passHltDEle33_narrowPull.root");
      files.push_back("./recoVsNc2_pt10_dr0p5_passHltDEle33_narrowPull.root");
      //files.push_back("./recoVsNc6_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc5_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc4_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc3_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc2_pt10_dr0p5_passHltDEle25_narrowerPull.root");
   } else {
      files.push_back("./nc6VsNc5_pt10_sharedHits3_passHltDEle33.root");
      files.push_back("./nc6VsNc4_pt10_sharedHits3_passHltDEle33.root");
      files.push_back("./nc6VsNc3_pt10_sharedHits3_passHltDEle33.root");
      files.push_back("./nc6VsNc2_pt10_sharedHits3_passHltDEle33.root");
      //files.push_back("./nc6VsNc5_pt10_sharedHits3_passHltDEle25_narrow.root");
      //files.push_back("./nc6VsNc4_pt10_sharedHits3_passHltDEle25_narrow.root");
      //files.push_back("./nc6VsNc3_pt10_sharedHits3_passHltDEle25_narrow.root");
      //files.push_back("./nc6VsNc2_pt10_sharedHits3_passHltDEle25_narrow.root");
   }
   TString dir("HltGsfTrkAna");

   vector<TString> legEntries;
   if (compToReco) legEntries.push_back("NC6");
   legEntries.push_back("NC5");
   legEntries.push_back("NC4");
   legEntries.push_back("NC3");
   legEntries.push_back("NC2");

   vector<TString> histoNames;
   histoNames.push_back("hMatch_pull_p");
   histoNames.push_back("hMatch_pull_pt");
   histoNames.push_back("hMatch_pull_eta");
   histoNames.push_back("hMatch_pull_ptError");
   histoNames.push_back("hMatch_pull_dPtOverPt");
   histoNames.push_back("hMatch_pull_etaError");
   histoNames.push_back("hMatch_pull_phiError");
   histoNames.push_back("hMatch_pull_outerP");
   histoNames.push_back("hMatch_pull_outerPt");
   histoNames.push_back("hMatch_pull_outerEta");
   histoNames.push_back("hMatch_pull_outerPOverP");
   histoNames.push_back("hMatch_pull_outerPtOverPt");
   histoNames.push_back("hMatch_pull_chi2");
   histoNames.push_back("hMatch_pull_ndof");
   histoNames.push_back("hMatch_pull_normalizedChi2");
   histoNames.push_back("hMatch_pull_dxy");
   histoNames.push_back("hMatch_pull_dz");
   histoNames.push_back("hMatch_pull_dxyError");
   histoNames.push_back("hMatch_pull_dzError");
   histoNames.push_back("hMatch_pull_numberOfValidHits");
   histoNames.push_back("hMatch_pull_numberOfLostHits");
   histoNames.push_back("hMatch_pull_validFraction");

   int colours[] = {kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};

   TFile inFile(files[0], "read");
   inFile.cd(dir);
   for (unsigned int i = 0; i < histoNames.size(); ++i) {
   //for (unsigned int i = 0; i < 3; ++i) {
      TCanvas *c0 = new TCanvas("c0_" + histoNames[i], "c0_" + histoNames[i], 100, 100, 800, 600);
      c0->cd();
      c0->SetBorderMode(0);
      c0->SetFrameBorderMode(0);
      c0->SetFillColor(0);
      c0->SetFrameFillColor(0);
      c0->SetRightMargin(0.123);
      c0->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      TLegend *legend = new TLegend(0.88, 0.68, 0.99, 0.90);
      legend->SetTextSize(0.03);
      legend->SetTextFont(font);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      vector<TH1F*> h;
      cout << histoNames[i] << " number of matched entries: ";
      for (unsigned int j = 0; j < files.size(); ++j) {
         TFile inFileT(files[j], "read");
         inFileT.cd(dir);
         h.push_back((TH1F *)gDirectory->Get(histoNames[i]));
         cout << "  " << h.back()->Integral();
         inFile.Close();

         if (underflow) {
            h.back()->SetBinContent(1, h.back()->GetBinContent(1) + h.back()->GetBinContent(0));
            //h.back()->SetBinError(1, sqrt(h.back()->GetBinError(1)*h.back()->GetBinError(1) + h.back()->GetBinError(0)*h.back()->GetBinError(0)));
         }
         if (overflow) {
            int nBin = h.back()->GetNbinsX();
            h.back()->SetBinContent(nBin, h.back()->GetBinContent(nBin) + h.back()->GetBinContent(nBin + 1));
            //h.back()->SetBinError(nBin, sqrt(h.back()->GetBinContent(nBin)*h.back()->GetBinContent(nBin) + h.back()->GetBinContent(nBin + 1)*h.back()->GetBinContent(nBin + 1)));
         }

         h.back()->GetXaxis()->SetTitle(h.back()->GetTitle()); 
         h.back()->GetXaxis()->SetTitleOffset(1.1);
         h.back()->SetLineColor(colours[j]);
         h.back()->SetLineWidth(2);
         if (j == 0) h.back()->Draw("h");
         else h.back()->Draw("hsame");
         legend->AddEntry(h.back(), legEntries[j], "l");
      }

      legend->Draw("same");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.03);
      if (compToReco) tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref: electronGsfTracks - RECO");
      else tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref Trigger: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7 (NC6)");
      //else tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref Trigger: HLT_Ele25_CaloIdVT_GsfTrkIdT_v2 (NC6)");

      if (saveAsPng) {
         if (compToReco) c0->Print(plotsDir + "hPull_vsReco_" + histoNames[i].Remove(0, 7) + ".png", "png");
         else c0->Print(plotsDir + "hPull_" + histoNames[i].Remove(0, 7) + ".png", "png");
      }
      cout << endl;
   }
}
