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

void macro_resolutions()
{
   TH1::AddDirectory(kFALSE);

   const bool compToReco = 1;
   int font = 42;
   const bool saveAsPng = 0;
   const bool saveFitsAsPng = 1;
   TString plotsDir = "./plots_20130930/";

   vector<TString> files;
   if (compToReco) {
      files.push_back("./recoVsNc6_pt10_dr0p5_passHltDEle33_narrowerPull.root");
      files.push_back("./recoVsNc5_pt10_dr0p5_passHltDEle33_narrowerPull.root");
      files.push_back("./recoVsNc4_pt10_dr0p5_passHltDEle33_narrowerPull.root");
      files.push_back("./recoVsNc3_pt10_dr0p5_passHltDEle33_narrowerPull.root");
      files.push_back("./recoVsNc2_pt10_dr0p5_passHltDEle33_narrowerPull.root");
      //files.push_back("./recoVsNc6_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc5_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc4_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc3_pt10_dr0p5_passHltDEle25_narrowerPull.root");
      //files.push_back("./recoVsNc2_pt10_dr0p5_passHltDEle25_narrowerPull.root");
   } else {
      files.push_back("./nc6VsNc5_pt10_sharedHits3_passHltDEle33_narrowDiff.root");
      files.push_back("./nc6VsNc4_pt10_sharedHits3_passHltDEle33_narrowDiff.root");
      files.push_back("./nc6VsNc3_pt10_sharedHits3_passHltDEle33_narrowDiff.root");
      files.push_back("./nc6VsNc2_pt10_sharedHits3_passHltDEle33_narrowDiff.root");
   }
   TString dir("HltGsfTrkAna");

   vector<TString> legEntries;
   if (compToReco) legEntries.push_back("NC6");
   legEntries.push_back("NC5");
   legEntries.push_back("NC4");
   legEntries.push_back("NC3");
   legEntries.push_back("NC2");

   vector<TString> axisLabels;
   if (compToReco) axisLabels.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7");
   axisLabels.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC5");
   axisLabels.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC4");
   axisLabels.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC3");
   axisLabels.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC2");
   //if (compToReco) axisLabels.push_back("HLT_Ele25_CaloIdVT_GsfTrkIdT_v2");
   //axisLabels.push_back("HLT_Ele25_CaloIdVT_GsfTrkIdT_nC5");
   //axisLabels.push_back("HLT_Ele25_CaloIdVT_GsfTrkIdT_nC4");
   //axisLabels.push_back("HLT_Ele25_CaloIdVT_GsfTrkIdT_nC3");
   //axisLabels.push_back("HLT_Ele25_CaloIdVT_GsfTrkIdT_nC2");

   vector<TString> histoNames;
   //histoNames.push_back("hMatch_diff_p");
   //histoNames.push_back("hMatch_diff_pt");
   //histoNames.push_back("hMatch_diff_eta");
   //histoNames.push_back("hMatch_diff_phi");
   //histoNames.push_back("hMatch_diff_ptError");
   //histoNames.push_back("hMatch_diff_dPtOverPt");
   //histoNames.push_back("hMatch_diff_etaError");
   //histoNames.push_back("hMatch_diff_phiError");
   //histoNames.push_back("hMatch_diff_outerP");
   //histoNames.push_back("hMatch_diff_outerPt");
   //histoNames.push_back("hMatch_diff_outerEta");
   //histoNames.push_back("hMatch_diff_outerPhi");
   //histoNames.push_back("hMatch_diff_outerPOverP");
   //histoNames.push_back("hMatch_diff_outerPtOverPt");
   //histoNames.push_back("hMatch_diff_dxy");
   //histoNames.push_back("hMatch_diff_dz");
   //histoNames.push_back("hMatch_diff_dxyError");
   //histoNames.push_back("hMatch_diff_dzError");
   //histoNames.push_back("hMatch_pull_p");
   histoNames.push_back("hMatch_pull_pt");
   //histoNames.push_back("hMatch_pull_eta");
   //histoNames.push_back("hMatch_pull_ptError");
   histoNames.push_back("hMatch_pull_dPtOverPt");
   //histoNames.push_back("hMatch_pull_etaError");
   //histoNames.push_back("hMatch_pull_phiError");
   //histoNames.push_back("hMatch_pull_outerP");
   //histoNames.push_back("hMatch_pull_outerPt");
   //histoNames.push_back("hMatch_pull_outerEta");
   //histoNames.push_back("hMatch_pull_outerPOverP");
   histoNames.push_back("hMatch_pull_outerPtOverPt");
   //histoNames.push_back("hMatch_pull_dxy");
   //histoNames.push_back("hMatch_pull_dz");

   int colours[] = {kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   TFile inFile(files[0], "read");
   inFile.cd(dir);
   for (unsigned int i = 0; i < histoNames.size(); ++i) {
   //for (unsigned int i = 0; i < 3; ++i) {
      TString plotVar = histoNames[i];
      plotVar.Remove(0, 12);

      TCanvas *c0 = new TCanvas("c0_" + histoNames[i], "c0_" + histoNames[i], 100, 100, 800, 600);
      c0->cd();
      c0->SetBorderMode(0);
      c0->SetFrameBorderMode(0);
      c0->SetFillColor(0);
      c0->SetFrameFillColor(0);
      c0->SetLeftMargin(0.162);
      c0->SetRightMargin(0.304);
      c0->SetBottomMargin(0.215);
      TCanvas *c1 = new TCanvas("c1_" + histoNames[i], "c1_" + histoNames[i], 100, 100, 800, 600);
      c1->cd();
      c1->SetBorderMode(0);
      c1->SetFrameBorderMode(0);
      c1->SetFillColor(0);
      c1->SetFrameFillColor(0);
      c1->SetLeftMargin(0.162);
      c1->SetRightMargin(0.304);
      c1->SetBottomMargin(0.215);
      TCanvas *c2 = new TCanvas("c2_" + histoNames[i], "c2_" + histoNames[i], 100, 100, 800, 600);
      c2->cd();
      c2->SetBorderMode(0);
      c2->SetFrameBorderMode(0);
      c2->SetFillColor(0);
      c2->SetFrameFillColor(0);
      c2->SetLeftMargin(0.162);
      c2->SetRightMargin(0.304);
      c2->SetBottomMargin(0.215);

      TH1F* h_means = new TH1F("h_means", "mean of Gaussian fit", files.size(), 0, files.size());
      TH1F* h_sigmas = new TH1F("h_sigmas", "sigma of Gaussian fit", files.size(), 0, files.size());
      TH1F* h_chi2s = new TH1F("h_chi2s", "#chi^2 of Gaussian fit", files.size(), 0, files.size());

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.03);

      h_means->SetLineColor(colours[0]);
      h_means->SetLineWidth(2);
      h_sigmas->SetLineColor(colours[0]);
      h_sigmas->SetLineWidth(2);
      h_chi2s->SetLineColor(colours[0]);
      h_chi2s->SetLineWidth(2);
      h_means->GetYaxis()->SetTitle("mean_{Gauss}");
      h_means->GetYaxis()->SetTitleOffset(1.7);
      h_sigmas->GetYaxis()->SetTitle("#sigma_{Gauss}");
      h_sigmas->GetYaxis()->SetTitleOffset(1.7);
      h_chi2s->GetYaxis()->SetTitle("#chi^{2}");
      h_chi2s->GetYaxis()->SetTitleOffset(1.2);

      vector<TH1F*> h;
      for (unsigned int j = 0; j < files.size(); ++j) {
      //for (unsigned int j = 0; j < 3; ++j) {
         TFile inFileT(files[j], "read");
         inFileT.cd(dir);
         h.push_back((TH1F *)gDirectory->Get(histoNames[i]));
         inFile.Close();
         TCanvas *c0f = new TCanvas("c0f_" + histoNames[i] + legEntries[j], "c0_" + histoNames[i] + legEntries[j], 100, 100, 800, 600); 
         c0f->cd();

         h.back()->GetXaxis()->SetTitle(h.back()->GetTitle()); 
         h.back()->Draw("h");

         TLegend *legend = new TLegend(0.88, 0.68, 0.99, 0.90);
         legend->SetTextSize(0.03);
         legend->SetTextFont(font);
         legend->SetBorderSize(0);
         legend->SetFillStyle(0);
         legend->AddEntry(h.back(), legEntries[j], "l");

         TF1 *gauss = new TF1("gauss", "gaus(0)", h.back()->GetXaxis()->GetXmin(), h.back()->GetXaxis()->GetXmax());
         gauss->SetParameters(1000., h.back()->GetMean(), h.back()->GetRMS());
         float factor = 0.7;
         h.back()->Fit("gauss", "", "", h.back()->GetMean() - factor*h.back()->GetRMS(), h.back()->GetMean() + factor*h.back()->GetRMS());
         legend->Draw("same");
         if (compToReco) tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref: electronGsfTracks - RECO");
         else tex->DrawLatex(0.10, 0.91, "CMS Preliminary,   Ref Trigger: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7 (NC6)");

         h_means->SetBinContent(j+1, gauss->GetParameter(1));
         h_means->SetBinError(j+1, gauss->GetParError(1));
         h_sigmas->SetBinContent(j+1, fabs(gauss->GetParameter(2)));
         h_sigmas->SetBinError(j+1, fabs(gauss->GetParError(2)));
         h_chi2s->SetBinContent(j+1, gauss->GetChisquare());
         h_chi2s->SetBinError(j+1, 0.);
         cout << "chi2: " << gauss->GetChisquare() << " , NDF: " << gauss->GetNDF() << endl;

         const char* axisLabel = axisLabels[j];
         h_means->GetXaxis()->SetBinLabel(j+1, axisLabel);
         h_sigmas->GetXaxis()->SetBinLabel(j+1, axisLabel);
         h_chi2s->GetXaxis()->SetBinLabel(j+1, axisLabel);

         if (saveAsPng && saveFitsAsPng) {
            if (compToReco) {
               c0f->Print(plotsDir + "hFit_vsReco_" + plotVar + "_" + legEntries[j] + ".png", "png");
            } else {
               c0f->Print(plotsDir + "hFit_" + plotVar + "_" + legEntries[j] + ".png", "png");
            }
         }
      }

      c0->cd();
      h_means->Draw("e0");
      if (compToReco) tex->DrawLatex(0.24, 0.91, "CMS Preliminary,   Ref: electronGsfTracks - RECO");
      else tex->DrawLatex(0.24, 0.91, "CMS Preliminary,   Ref Trigger: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7 (NC6)");
      tex->DrawLatex(0.71, 0.87, plotVar);
      c1->cd();
      h_sigmas->Draw("e0");
      if (compToReco) tex->DrawLatex(0.24, 0.91, "CMS Preliminary,   Ref: electronGsfTracks - RECO");
      else tex->DrawLatex(0.24, 0.91, "CMS Preliminary,   Ref Trigger: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7 (NC6)");
      tex->DrawLatex(0.71, 0.87, plotVar);
      c2->cd();
      h_chi2s->Draw("e0");
      if (compToReco) tex->DrawLatex(0.24, 0.91, "CMS Preliminary,   Ref: electronGsfTracks - RECO");
      else tex->DrawLatex(0.24, 0.91, "CMS Preliminary,   Ref Trigger: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7 (NC6)");
      tex->DrawLatex(0.71, 0.87, plotVar);

      if (saveAsPng) {
         if (compToReco) {
            c0->Print(plotsDir + "hFit_mean_vsReco_" + histoNames[i].Remove(0, 7) + ".png", "png");
            c1->Print(plotsDir + "hFit_sigma_vsReco_" + histoNames[i] + ".png", "png");
            c2->Print(plotsDir + "hFit_chi2_vsReco_" + histoNames[i] + ".png", "png");
         } else {
            c0->Print(plotsDir + "hFit_mean_" + histoNames[i].Remove(0, 7) + ".png", "png");
            c1->Print(plotsDir + "hFit_sigma_" + histoNames[i] + ".png", "png");
            c2->Print(plotsDir + "hFit_chi2_" + histoNames[i] + ".png", "png");
         }
      }
   }
}
