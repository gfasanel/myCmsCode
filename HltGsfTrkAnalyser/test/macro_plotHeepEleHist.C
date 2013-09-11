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

void macro_plotHeepEleHist()
{
   TH1::AddDirectory(kFALSE);

   int font = 42;
   const bool saveAsPng = 0;
   TString plotsDir = "../plots_20130726/";
   TString fileNameExtra = ""; 

   vector<TString> files;
   //files.push_back("../hltHeepEleHistos.root");
   files.push_back("../histoele.root");
   TString dir("HltGsfEleAna");

   vector<TString> h1Names;
   h1Names.push_back("h_total");
   h1Names.push_back("h_1PassHeep_total");
   h1Names.push_back("h_2PassHeep_total");
   h1Names.push_back("h_nPassHeep_total");
   h1Names.push_back("h_nPassHeep_trgd");
   h1Names.push_back("h_nPassHeep_notTrgd");
   h1Names.push_back("h_nPassHeep_nPassEt_total");
   h1Names.push_back("h_nPassHeep_nPassEt_trgd");
   h1Names.push_back("h_nPassHeep_nPassEt_notTrgd");
   h1Names.push_back("h_nHeep");
   h1Names.push_back("h_trgd");
   h1Names.push_back("h_notTrgd");
   h1Names.push_back("hRatio_nPassHeep_trgdVsTotal");
   h1Names.push_back("hRatio_nPassHeep_notTrgdVsTotal");
   h1Names.push_back("hRatio_nPassHeep_nPassEt_trgdVsTotal");
   h1Names.push_back("hRatio_nPassHeep_nPassEt_notTrgdVsTotal");
   vector<TString> h2Names;
   h2Names.push_back("h2Ratio_nPassHeep_trgd");
   h2Names.push_back("h2Ratio_nPassHeep_notTrgd");
   h2Names.push_back("h2Ratio_nPassHeep_nPassEt_trgd");
   h2Names.push_back("h2Ratio_nPassHeep_nPassEt_notTrgd");
   h2Names.push_back("h2Ratio_trgd");
   h2Names.push_back("h2Ratio_notTrgd");

   vector<pair<double, double> > zPlotRanges;
   zPlotRanges.push_back(make_pair(0.99, 1.01));
   zPlotRanges.push_back(make_pair(0.8, 1.2));
   zPlotRanges.push_back(make_pair(0.995, 1.005));
   zPlotRanges.push_back(make_pair(0.75, 1.25));
   zPlotRanges.push_back(make_pair(0.96, 1.04));
   zPlotRanges.push_back(make_pair(0.6, 1.4));

   int colours[] = {kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};

   TFile inFile(files[0], "read");
   inFile.cd(dir);
   for (unsigned int i = 0; i < h1Names.size(); ++i) {
   //for (unsigned int i = 0; i < 1; ++i) {
      TH1F* h = (TH1F*)gDirectory->Get(h1Names[i]);

      TCanvas *c0;
      TPad *h1Pad;
      c0 = new TCanvas("c0_" + h1Names[i], "c0_" + h1Names[i], 100, 100, 800, 600);
      h1Pad = new TPad("h1Pad" + h1Names[i], "emu Spectrum" + h1Names[i], 0., 0., 1., 1.);
      h1Pad->SetBorderMode(0);
      h1Pad->SetFrameBorderMode(0);
      h1Pad->SetFillColor(0);
      h1Pad->SetFrameFillColor(0);
      h1Pad->SetRightMargin(0.16);
      h1Pad->SetBottomMargin(0.18);
      //h1Pad->SetLogy();
      h1Pad->Draw();
      h1Pad->cd();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      h->SetLineColor(colours[0]);
      h->SetLineWidth(2);
      h->Draw("e");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.035);
      tex->DrawLatex(0.03, 0.96, h->GetTitle());
      tex->SetTextSize(0.03);
      tex->DrawLatex(0.10, 0.91, "CMS Preliminary");

      if (saveAsPng) {
         c0->Print(plotsDir + h1Names[i] + fileNameExtra + ".png", "png");
      }
   }

   for (unsigned int i = 0; i < h2Names.size(); ++i) {
   //for (unsigned int i = 0; i < 1; ++i) {
      TH2F* h2 = (TH2F*)gDirectory->Get(h2Names[i]);

      TCanvas *c0;
      TPad *h2Pad;
      c0 = new TCanvas("c0_" + h2Names[i], "c0_" + h2Names[i], 100, 100, 800, 600);
      h2Pad = new TPad("h2Pad" + h2Names[i], h2Names[i], 0., 0., 1., 1.);
      h2Pad->SetBorderMode(0);
      h2Pad->SetFrameBorderMode(0);
      h2Pad->SetFillColor(0);
      h2Pad->SetFrameFillColor(0);
      h2Pad->SetLeftMargin(0.32);
      h2Pad->SetRightMargin(0.18);
      h2Pad->SetBottomMargin(0.18);
      h2Pad->SetGridx();
      h2Pad->SetGridy();
      //h2Pad->SetLogy();
      h2Pad->Draw();
      h2Pad->cd();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      h2->SetAxisRange(zPlotRanges[i].first, zPlotRanges[i].second, "Z");
      h2->DrawClone("COLZ");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.035);
      tex->DrawLatex(0.03, 0.96, h2->GetTitle());
      tex->SetTextSize(0.03);
      tex->DrawLatex(0.32, 0.91, "CMS Preliminary");
      tex->DrawLatex(0.687, 0.91, "Trigger X / Trigger Y");

      if (saveAsPng) {
         c0->Print(plotsDir + h2Names[i] + fileNameExtra + ".png", "png");
      }

      //TCanvas *c1;
      //TPad *h2Pad3D;
      //c1 = new TCanvas("c1_" + h2Names[i], "c1_" + h2Names[i], 100, 100, 800, 600);
      //h2Pad3D = new TPad("h2Pad3D" + h2Names[i], h2Names[i], 0., 0., 1., 1.);
      //h2Pad3D->SetBorderMode(0);
      //h2Pad3D->SetFrameBorderMode(0);
      //h2Pad3D->SetFillColor(0);
      //h2Pad3D->SetFrameFillColor(0);
      //h2Pad->SetRightMargin(0.25);
      ////h2Pad3D->SetLogy();
      //h2Pad3D->Draw();
      //h2Pad3D->cd();
      //gStyle->SetOptStat(0);
      //gStyle->SetOptTitle(0);
      //gStyle->SetPadTickX(1);
      //gStyle->SetPadTickY(1);

      //for (unsigned int i = 1; i <= h2->GetNbinsX(); ++i) {
      //   h2->GetXaxis()->SetBinLabel(i, Form("%i", i));
      //   h2->GetYaxis()->SetBinLabel(i, Form("%i", i));
      //}

      //h2->Draw("LEGO2");

      //TLatex *tex = new TLatex();
      //tex->SetNDC();
      //tex->SetTextFont(font);
      //tex->SetLineWidth(2);
      //tex->SetTextSize(0.035);
      //tex->DrawLatex(0.03, 0.96, h2->GetTitle());
      //tex->SetTextSize(0.03);
      //tex->DrawLatex(0.03, 0.91, "CMS Preliminary");
      //tex->DrawLatex(0.03, 0.86, "Trigger X / Trigger Y");

      //if (saveAsPng) {
      //   c1->Print(plotsDir + h2Names[i] + "_3d_" + fileNameExtra + ".png", "png");
      //}
   }
}
