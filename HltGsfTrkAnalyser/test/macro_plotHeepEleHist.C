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
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include <TCanvas.h>
#include <TStyle.h>

void macro_plotHeepEleHist()
{
   TH1::AddDirectory(kFALSE);

   const bool plotHistos = 0;
   const bool plotEff = 1;
   const bool plot2DRatios = 0;

   int font = 42;
   const bool saveAsPng = 0;
   TString plotsDir = "./plots_20140318/";
   TString fileNameExtra = "";
   //TString fileNameExtra = "_realEle";
   //TString fileNameExtra = "_qcd";

   vector<TString> files;
   //files.push_back("../hltHeepEleHistos.root");
   files.push_back("./histoele_dataZskim.root");
   files.push_back("./histoele_dy_8tev50ns.root");
   files.push_back("./histoele_qcd_8tev50ns.root");
   //files.push_back("./histoele_dataZskim_lowEtMenu.root");
   //files.push_back("./histoele_dy_8tev50ns_lowEtMenu.root");
   //files.push_back("./histoele_qcd_8tev50ns_lowEtMenu.root");
   TString dir("HltGsfEleAna");
   vector<TString> fileTitles;
   fileTitles.push_back("_dataZskim");
   fileTitles.push_back("_dy");
   fileTitles.push_back("_qcd");
   vector<TString> legTitles;
   legTitles.push_back("#splitline{Data}{(Z skim)}");
   legTitles.push_back("DY");
   legTitles.push_back("QCD");

   vector<TString> h1Names;
   h1Names.push_back("h_nHeep");
   h1Names.push_back("h_total");
   h1Names.push_back("h_refTrgd");
   h1Names.push_back("h_trgd");
   h1Names.push_back("h_notTrgdRefTrgd");
   h1Names.push_back("h_1PassHeep_refTrgd");
   h1Names.push_back("h_2PassHeep_refTrgd");
   h1Names.push_back("h_nPassHeep_refTrgd");
   h1Names.push_back("h_nPassHeep_trgd");
   h1Names.push_back("h_nPassHeep_notTrgdRefTrgd");
   h1Names.push_back("h_nPassHeep_nPassEt_refTrgd");
   h1Names.push_back("h_nPassHeep_nPassEt_trgd");
   h1Names.push_back("h_nPassHeep_nPassEt_notTrgdRefTrgd");
 
   // for efficiency plots
   vector<TString> gEffHistNamesNum;
   gEffHistNamesNum.push_back("h_trgd");
   gEffHistNamesNum.push_back("h_notTrgdRefTrgd");
   gEffHistNamesNum.push_back("h_nPassHeep_trgd");
   gEffHistNamesNum.push_back("h_nPassHeep_notTrgdRefTrgd");
   gEffHistNamesNum.push_back("h_nPassHeep_nPassEt_trgd");
   gEffHistNamesNum.push_back("h_nPassHeep_nPassEt_notTrgdRefTrgd");
   vector<TString> gEffHistNamesDen;
   gEffHistNamesDen.push_back("h_refTrgd");
   gEffHistNamesDen.push_back("h_refTrgd");
   gEffHistNamesDen.push_back("h_nPassHeep_refTrgd");
   gEffHistNamesDen.push_back("h_nPassHeep_refTrgd");
   gEffHistNamesDen.push_back("h_nPassHeep_nPassEt_refTrgd");
   gEffHistNamesDen.push_back("h_nPassHeep_nPassEt_refTrgd");
   vector<TString> gEffTitles;
   gEffTitles.push_back("# triggered / # reference triggered events");
   gEffTitles.push_back("# not triggered but reference triggered / # reference triggered events");
   gEffTitles.push_back("# triggered / # reference triggered n-HEEP events");
   gEffTitles.push_back("# not triggered but reference triggered / # reference triggered n-HEEP events");
   gEffTitles.push_back("# triggered / # reference triggered n-HEEP events passing Et cuts");
   gEffTitles.push_back("# not triggered but reference triggered / # reference triggered n-HEEP events passing Et cuts");
   //vector<TString> gEffHistNamesDen;
   //gEffHistNamesDen.push_back("h_refTrgd");
   //gEffHistNamesDen.push_back("h_refTrgd");
   //gEffHistNamesDen.push_back("h_refTrgd");
   //gEffHistNamesDen.push_back("h_refTrgd");
   //gEffHistNamesDen.push_back("h_refTrgd");
   //gEffHistNamesDen.push_back("h_refTrgd");
   //vector<TString> gEffTitles;
   //gEffTitles.push_back("# triggered / # reference triggered events");
   //gEffTitles.push_back("# not triggered but reference triggered / # reference triggered events");
   //gEffTitles.push_back("# triggered n-HEEP events / # reference triggered events");
   //gEffTitles.push_back("# not triggered but reference triggered n-HEEP events / # reference triggered events");
   //gEffTitles.push_back("# triggered n-HEEP events passing Et cuts / # reference triggered events");
   //gEffTitles.push_back("# not triggered but reference triggered n-HEEP events passing Et cuts / # reference triggered events");

   // 2D efficiency matrix 
   vector<TString> h2Names;
   h2Names.push_back("h2Ratio_trgd");
   h2Names.push_back("h2Ratio_notTrgdRefTrgd");
   h2Names.push_back("h2Ratio_nPassHeep_trgd");
   h2Names.push_back("h2Ratio_nPassHeep_notTrgdRefTrgd");
   h2Names.push_back("h2Ratio_nPassHeep_nPassEt_trgd");
   h2Names.push_back("h2Ratio_nPassHeep_nPassEt_notTrgdRefTrgd");

   vector<pair<double, double> > zPlotRanges;
   // data
   zPlotRanges.push_back(make_pair(0.994, 1.006));
   zPlotRanges.push_back(make_pair(0.92, 1.08));
   zPlotRanges.push_back(make_pair(0.996, 1.004));
   zPlotRanges.push_back(make_pair(0.85, 1.15));
   zPlotRanges.push_back(make_pair(0.998, 1.002));
   zPlotRanges.push_back(make_pair(0.7, 1.3));
   // dy
   zPlotRanges.push_back(make_pair(0.996, 1.004));
   zPlotRanges.push_back(make_pair(0.88, 1.12));
   zPlotRanges.push_back(make_pair(0.997, 1.003));
   zPlotRanges.push_back(make_pair(0.85, 1.15));
   zPlotRanges.push_back(make_pair(0.997, 1.003));
   zPlotRanges.push_back(make_pair(0.4, 2.2));
   // qcd
   zPlotRanges.push_back(make_pair(0.98, 1.02));
   zPlotRanges.push_back(make_pair(0.95, 1.05));
   zPlotRanges.push_back(make_pair(0.97, 1.03));
   zPlotRanges.push_back(make_pair(0.92, 1.08));
   zPlotRanges.push_back(make_pair(0.97, 1.03));
   zPlotRanges.push_back(make_pair(0.8, 1.2));
   //// data
   //zPlotRanges.push_back(make_pair(0.9985, 1.0015));
   //zPlotRanges.push_back(make_pair(0.95, 1.05));
   //zPlotRanges.push_back(make_pair(0.99985, 1.00015));
   //zPlotRanges.push_back(make_pair(0.98, 1.02));
   //zPlotRanges.push_back(make_pair(0.99975, 1.00025));
   //zPlotRanges.push_back(make_pair(0.95, 1.05));
   //// dy
   //zPlotRanges.push_back(make_pair(0.9996, 1.0004));
   //zPlotRanges.push_back(make_pair(0.99, 1.01));
   //zPlotRanges.push_back(make_pair(0.9996, 1.0004));
   //zPlotRanges.push_back(make_pair(0.96, 1.04));
   //zPlotRanges.push_back(make_pair(0.9996, 1.0004));
   //zPlotRanges.push_back(make_pair(0.94, 1.06));
   //// qcd
   //zPlotRanges.push_back(make_pair(0.992, 1.008));
   //zPlotRanges.push_back(make_pair(0.996, 1.004));
   //zPlotRanges.push_back(make_pair(0.97, 1.03));
   //zPlotRanges.push_back(make_pair(0.92, 1.08));
   //zPlotRanges.push_back(make_pair(0.97, 1.03));
   //zPlotRanges.push_back(make_pair(0.8, 1.2));

   int colours[] = {kBlack, kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};
   //int colours[] = {kRed, kGreen, kMagenta, kOrange, kCyan};

   TFile* inFile;

   if (plotHistos) {
      for (unsigned int i = 0; i < h1Names.size(); ++i) {
      //for (unsigned int i = 1; i < 2; ++i) {
         TCanvas *c0;
         TPad *h1Pad;
         c0 = new TCanvas("c0_" + h1Names[i], "c0_" + h1Names[i], 100, 100, 800, 600);
         h1Pad = new TPad("h1Pad" + h1Names[i], "Spectrum" + h1Names[i], 0., 0., 1., 1.);
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
         TLegend *legend = new TLegend(0.85, 0.68, 0.99, 0.90);
         legend->SetTextSize(0.03);
         legend->SetTextFont(font);
         legend->SetBorderSize(0);
         legend->SetFillStyle(0);
   
         vector<TH1F*> h;
         double minimum = 1.e99;
         double maximum = 0.;
         for (unsigned int fInd = 0; fInd < files.size(); ++fInd) {
            inFile = new TFile(files[fInd], "read");
            inFile->cd(dir);
            h.push_back((TH1F*)gDirectory->Get(h1Names[i]));
            inFile->Close();
            h.back()->SetLineColor(colours[fInd]);
            h.back()->SetLineWidth(2);
            if (h.back()->GetMinimum() < minimum) minimum = 0.8*h.back()->GetMinimum();
            if (h.back()->GetMaximum() > maximum) maximum = 1.2*h.back()->GetMaximum();
         }
         for (unsigned int fInd = 0; fInd < files.size(); ++fInd) {
            h.at(fInd)->SetMinimum(minimum);
            h.at(fInd)->SetMaximum(maximum);
            if (fInd == 0) h.at(fInd)->Draw("e");
            else h.at(fInd)->Draw("esame");
            legend->AddEntry(h.at(fInd), legTitles[fInd], "l");
         }
         TLine* l = new TLine(7, minimum, 7, maximum);
         if (i > 0) l->Draw("same");
         legend->Draw("same");
   
         TLatex *tex = new TLatex();
         tex->SetNDC();
         tex->SetTextFont(font);
         tex->SetLineWidth(2);
         tex->SetTextSize(0.035);
         tex->DrawLatex(0.02, 0.96, h.front()->GetTitle());
         tex->SetTextSize(0.03);
         tex->DrawLatex(0.10, 0.91, "CMS Preliminary");
   
         if (saveAsPng) {
            c0->Print(plotsDir + h1Names[i] + fileNameExtra + ".png", "png");
         }
      }
   }

   if (plotEff) {
      for (unsigned int i = 0; i < gEffHistNamesNum.size(); ++i) {
      //for (unsigned int i = 5; i < 6; ++i) {
      //for (unsigned int i = 0; i < 1; ++i) {
         TCanvas *c1;
         TPad *gEffPad;
         c1 = new TCanvas("c1_" + gEffHistNamesNum[i], "c1_" + gEffHistNamesNum[i], 100, 100, 800, 600);
         gEffPad = new TPad("gEffPad" + gEffHistNamesNum[i], "gEffPad" + gEffHistNamesNum[i], 0., 0., 1., 1.);
         gEffPad->SetBorderMode(0);
         gEffPad->SetFrameBorderMode(0);
         gEffPad->SetFillColor(0);
         gEffPad->SetFrameFillColor(0);
         gEffPad->SetRightMargin(0.16);
         gEffPad->SetBottomMargin(0.18);
         gEffPad->Draw();
         gEffPad->cd();
         gStyle->SetOptStat(0);
         gStyle->SetOptTitle(0);
         gStyle->SetPadTickX(1);
         gStyle->SetPadTickY(1);
         gStyle->SetEndErrorSize(0);
         TLegend *legend = new TLegend(0.85, 0.68, 0.99, 0.90);
         legend->SetTextSize(0.03);
         legend->SetTextFont(font);
         legend->SetBorderSize(0);
         legend->SetFillStyle(0);
   
         vector<TH1F*> h;
         vector<TH1F*> hTot;
         vector<TH1F*> hDiv;
         double minimum = 1.;
         double maximum = 0.;
         for (unsigned int fInd = 0; fInd < files.size(); ++fInd) {
            inFile = new TFile(files[fInd], "read");
            inFile->cd(dir);
            h.push_back((TH1F*)gDirectory->Get(gEffHistNamesNum[i]));
            hTot.push_back((TH1F*)gDirectory->Get(gEffHistNamesDen[i]));
            // necessary to plot x axis with labels
            hDiv.push_back((TH1F*)h.back()->Clone("hDiv"));
            hDiv.back()->Divide(hTot.back());
            TGraphAsymmErrors* gEff = new TGraphAsymmErrors(h.back(), hTot.back(), "cp");
            inFile->Close();
   
            hDiv.back()->SetLineColor(gEffPad->GetFrameFillColor());
            for (int bin = 0; bin < hDiv.back()->GetNbinsX(); ++bin) {
               if (hDiv.back()->GetMinimum() < minimum) minimum = hDiv.back()->GetMinimum();
               if (minimum > hDiv.back()->GetBinContent(bin+1) - 1.0*gEff->GetErrorYlow(bin)) minimum = hDiv.back()->GetBinContent(bin+1) - 1.0*gEff->GetErrorYlow(bin);
               if (hDiv.back()->GetMaximum() > maximum) maximum = hDiv.back()->GetMaximum();
               if (maximum < hDiv.back()->GetBinContent(bin+1) + 1.4*gEff->GetErrorYhigh(bin)) maximum = hDiv.back()->GetBinContent(bin+1) + 1.4*gEff->GetErrorYhigh(bin);
            }
         }
         TString gEffName = "gEff_" + gEffHistNamesNum[i] + "_vs_" + gEffHistNamesDen[i];
         gEffName.ReplaceAll("h_", "");
         for (unsigned int fInd = 0; fInd < files.size(); ++fInd) {
            TGraphAsymmErrors* gEff = new TGraphAsymmErrors(h.at(fInd), hTot.at(fInd), "cp");
            hDiv.at(fInd)->SetMinimum(minimum);
            hDiv.at(fInd)->SetMaximum(maximum);
            if (fInd == 0) hDiv.at(fInd)->Draw("hist");
            else hDiv.at(fInd)->Draw("histsame");
            gEff->SetLineColor(colours[fInd]);
            gEff->SetLineWidth(2);
            gEff->SetName(gEffName);
            gEff->SetTitle(gEffTitles[i]);
            gEff->DrawClone("p");
            legend->AddEntry(gEff, legTitles[fInd], "l");
         }
         TLine* l = new TLine(7, minimum, 7, maximum);
         l->Draw("same");
         legend->Draw("same");
   
         TLatex *tex = new TLatex();
         tex->SetNDC();
         tex->SetTextFont(font);
         tex->SetLineWidth(2);
         tex->SetTextSize(0.035);
         tex->DrawLatex(0.02, 0.96, gEffTitles[i]);
         tex->SetTextSize(0.03);
         tex->DrawLatex(0.10, 0.91, "CMS Preliminary");
         tex->DrawLatex(0.62, 0.91, "Errors: Clopper-Pearson");

         if (saveAsPng) {
            c1->Print(plotsDir + gEffName + fileNameExtra + ".png", "png");
         }
      }
   }

   if (plot2DRatios) {
      for (unsigned int fInd = 0; fInd < files.size(); ++fInd) {
         inFile = new TFile(files[fInd], "read");
         inFile->cd(dir);
         for (unsigned int i = 0; i < h2Names.size(); ++i) {
         //for (unsigned int i = 0; i < 1; ++i) {
            TH2F* h2 = (TH2F*)gDirectory->Get(h2Names[i]);
   
            TCanvas *c2;
            TPad *h2Pad;
            c2 = new TCanvas("c2_" + h2Names[i] + fileTitles[fInd], "c2_" + h2Names[i] + fileTitles[fInd], 100, 100, 800, 600);
            h2Pad = new TPad("h2Pad" + h2Names[i] + fileTitles[fInd], h2Names[i] + fileTitles[fInd], 0., 0., 1., 1.);
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
   
            h2->SetAxisRange(zPlotRanges[i + fInd*h2Names.size()].first, zPlotRanges[i + fInd*h2Names.size()].second, "Z");
            h2->DrawClone("COLZ");
   
            TLatex *tex = new TLatex();
            tex->SetNDC();
            tex->SetTextFont(font);
            tex->SetLineWidth(2);
            tex->SetTextSize(0.035);
            tex->DrawLatex(0.02, 0.96, h2->GetTitle());
            tex->SetTextSize(0.03);
            tex->DrawLatex(0.32, 0.91, "CMS Preliminary");
            tex->DrawLatex(0.687, 0.91, "Trigger X / Trigger Y");
            tex->DrawLatex(0.23, 0.1, legTitles[fInd]);
   
            if (saveAsPng) {
               c2->Print(plotsDir + h2Names[i] + fileNameExtra + fileTitles[fInd] + ".png", "png");
            }
   
            //TCanvas *c3;
            //TPad *h2Pad3D;
            //c3 = new TCanvas("c3_" + h2Names[i], "c3_" + h2Names[i], 100, 100, 800, 600);
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
            //   c3->Print(plotsDir + h2Names[i] + "_3d_" + fileNameExtra + ".png", "png");
            //}
         }
         inFile->Close();
      }
   }
}
