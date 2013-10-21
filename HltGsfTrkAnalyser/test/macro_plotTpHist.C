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

void macro_plotTpHist()
{
   TH1::AddDirectory(kFALSE);

   const bool plotHistos = 1;
   const bool plotEtHistos = 0;
   const bool plotEff = 0;
   const bool plotEtEff = 0;

   int font = 42;
   const bool saveAsPng = 0;
   TString plotsDir = "./plots_20130926/";
   //TString fileNameExtra = "_realEle";
   TString fileNameExtra = "";

   vector<TString> files;
   files.push_back("./histoele_dataZskim.root");
   files.push_back("./histoele_dy_8tev50ns.root");
   files.push_back("./histoele_qcd_8tev50ns.root");
   TString dir("HltGsfEleAna");
   vector<TString> fileTitles;
   fileTitles.push_back("_dataZskim");
   fileTitles.push_back("_dy");
   fileTitles.push_back("_qcd");
   vector<TString> legTitles;
   legTitles.push_back("#splitline{Data}{(Z skim)}");
   legTitles.push_back("DY");
   legTitles.push_back("QCD");

   vector<TString> hNames;
   hNames.push_back("hTp_tags");
   hNames.push_back("hTp_probes");
   hNames.push_back("hTp_allProbesGsf");
   hNames.push_back("hTp_passProbesGsf");
   hNames.push_back("hTp_allProbesGsf_et33");
   hNames.push_back("hTp_passProbesGsf_et33");
   hNames.push_back("hTp_allProbesGsf_et35");
   hNames.push_back("hTp_passProbesGsf_et35");
   hNames.push_back("hTp_allProbesGsf_et80");
   hNames.push_back("hTp_passProbesGsf_et80");
   hNames.push_back("hTp_allProbesGsf_et100");
   hNames.push_back("hTp_passProbesGsf_et100");
   hNames.push_back("hTp_allProbesGsf_withEtCut");
   hNames.push_back("hTp_passProbesGsf_withEtCut");
   hNames.push_back("hTp_allProbesHeep");
   hNames.push_back("hTp_passProbesHeep");
   hNames.push_back("hTp_allProbesHeep_et33");
   hNames.push_back("hTp_passProbesHeep_et33");
   hNames.push_back("hTp_allProbesHeep_et35");
   hNames.push_back("hTp_passProbesHeep_et35");
   hNames.push_back("hTp_allProbesHeep_et80");
   hNames.push_back("hTp_passProbesHeep_et80");
   hNames.push_back("hTp_allProbesHeep_et100");
   hNames.push_back("hTp_passProbesHeep_et100");
   hNames.push_back("hTp_allProbesHeep_withEtCut");
   hNames.push_back("hTp_passProbesHeep_withEtCut");
   //hNames.push_back("hTpRatio_passProbesGsf_vs_allProbesGsf");
   //hNames.push_back("hTpRatio_passProbesGsf_et33_vs_allProbesGsf_et33");
   //hNames.push_back("hTpRatio_passProbesGsf_et35_vs_allProbesGsf_et35");
   //hNames.push_back("hTpRatio_passProbesGsf_et80_vs_allProbesGsf_et80");
   //hNames.push_back("hTpRatio_passProbesGsf_et100_vs_allProbesGsf_et100");
   //hNames.push_back("hTpRatio_passProbesGsf_withEtCut_vs_allProbesGsf_withEtCut");
   //hNames.push_back("hTpRatio_passProbesHeep_vs_allProbesHeep");
   //hNames.push_back("hTpRatio_passProbesHeep_et33_vs_allProbesHeep_et33");
   //hNames.push_back("hTpRatio_passProbesHeep_et35_vs_allProbesHeep_et35");
   //hNames.push_back("hTpRatio_passProbesHeep_et80_vs_allProbesHeep_et80");
   //hNames.push_back("hTpRatio_passProbesHeep_et100_vs_allProbesHeep_et100");
   //hNames.push_back("hTpRatio_passProbesHeep_withEtCut_vs_allProbesHeep_withEtCut");

   vector<TString> gEffHistNamesPrefNum;
   gEffHistNamesPrefNum.push_back("hTp_passProbesGsf_et");
   gEffHistNamesPrefNum.push_back("hTp_passProbesHeep_et");
   vector<TString> gEffHistNamesPrefDen;
   gEffHistNamesPrefDen.push_back("hTp_allProbesGsf_et");
   gEffHistNamesPrefDen.push_back("hTp_allProbesHeep_et");

   // for efficiency plots
   vector<TString> gEffTitles;

   int colours[] = {kBlack, kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   TFile* inFile;

   if (plotHistos) {
      for (unsigned int i = 0; i < hNames.size(); ++i) {
      //for (unsigned int i = 1; i < 2; ++i) {
         TCanvas *c0;
         TPad *h1Pad;
         c0 = new TCanvas("c0_" + hNames[i], "c0_" + hNames[i], 100, 100, 800, 600);
         h1Pad = new TPad("h1Pad" + hNames[i], "Spectrum" + hNames[i], 0., 0., 1., 1.);
         h1Pad->SetBorderMode(0);
         h1Pad->SetFrameBorderMode(0);
         h1Pad->SetFillColor(0);
         h1Pad->SetFrameFillColor(0);
         h1Pad->SetRightMargin(0.16);
         h1Pad->SetBottomMargin(0.18);
         //h1Pad->SetLogy();
         h1Pad->Draw();
         h1Pad->cd();
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
            h.push_back((TH1F*)gDirectory->Get(hNames[i]));
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
         h.front()->Draw("sameaxis");
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
            c0->Print(plotsDir + hNames[i] + fileNameExtra + ".png", "png");
         }
      }
   }

   // find out how many paths we are plotting
   inFile = new TFile(files[0], "read");
   inFile->cd(dir);
   TH1F* histHelper = (TH1F*)gDirectory->Get("h_total");
   inFile->Close();

   vector<TString> pathNames;
   for (int bin = 1; bin <= histHelper->GetNbinsX(); ++bin) {
      pathNames.push_back(histHelper->GetXaxis()->GetBinLabel(bin));
      pathNames.back().ReplaceAll(" [ref]", "");
   }

   if (plotEtHistos) {
      unsigned int probeInd = 0;
      unsigned int totCtr = 0;
      while (totCtr < gEffHistNamesPrefNum.size() + gEffHistNamesPrefDen.size()) {
         for (unsigned int i = 0; i < pathNames.size(); ++i) {
         //for (unsigned int i = 0; i < 1; ++i) {
            TString histoName = gEffHistNamesPrefNum[probeInd] + "_" + pathNames[i];
            if (totCtr % 2 == 1) histoName = gEffHistNamesPrefDen[probeInd] + "_" + pathNames[i];

            TCanvas *c0;
            TPad *h1Pad;
            c0 = new TCanvas("c0_" + histoName, "c0_" + histoName, 100, 100, 800, 600);
            h1Pad = new TPad("h1Pad" + histoName, "Spectrum" + histoName, 0., 0., 1., 1.);
            h1Pad->SetBorderMode(0);
            h1Pad->SetFrameBorderMode(0);
            h1Pad->SetFillColor(0);
            h1Pad->SetFrameFillColor(0);
            h1Pad->SetRightMargin(0.16);
            h1Pad->SetBottomMargin(0.18);
            //h1Pad->SetLogy();
            h1Pad->Draw();
            h1Pad->cd();
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
               h.push_back((TH1F*)gDirectory->Get(histoName));
               inFile->Close();
               h.back()->SetLineColor(colours[fInd]);
               h.back()->SetLineWidth(2);
               if (h.back()->GetMinimum() < minimum) minimum = 0.7*h.back()->GetMinimum();
               if (h.back()->GetMaximum() > maximum) maximum = 1.3*h.back()->GetMaximum();
            }
            for (unsigned int fInd = 0; fInd < files.size(); ++fInd) {
               h.at(fInd)->SetMinimum(minimum);
               h.at(fInd)->SetMaximum(maximum);
               if (fInd == 0) h.at(fInd)->Draw("e");
               else h.at(fInd)->Draw("esame");
               legend->AddEntry(h.at(fInd), legTitles[fInd], "l");
            }
            h.front()->Draw("sameaxis");
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
            tex->DrawLatex(0.32, 0.91, pathNames[i]);
   
            if (saveAsPng) {
               c0->Print(plotsDir + histoName + fileNameExtra + ".png", "png");
            }
         }
         if (totCtr % 2 == 1) ++probeInd;
         ++totCtr;
      }
   }

   if (plotEff) {
      for (unsigned int i = 0; i < hNames.size(); ++i) {
         TString hNumName;
         TString hDenName;
         if (hNames[i].BeginsWith("hTp_all")) {
            hDenName = hNames[i];
            hNumName = hNames[i+1];
            ++i;
         }

         if (hNumName.Length() == 0 || hDenName.Length() == 0) continue;

         TCanvas *c1;
         TPad *gEffPad;
         c1 = new TCanvas("c1_" + hNumName + "_vs_" + hDenName, "c1_" + hNumName + "_vs_" + hDenName, 100, 100, 800, 600);
         gEffPad = new TPad("gEffPad" + hNumName + "_vs_" + hDenName, "gEffPad" + hNumName + "_vs_" + hDenName, 0., 0., 1., 1.);
         gEffPad->SetBorderMode(0);
         gEffPad->SetFrameBorderMode(0);
         gEffPad->SetFillColor(0);
         gEffPad->SetFrameFillColor(0);
         gEffPad->SetRightMargin(0.16);
         gEffPad->SetBottomMargin(0.11);
         gEffPad->Draw();
         gEffPad->cd();
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
            h.push_back((TH1F*)gDirectory->Get(hNumName));
            hTot.push_back((TH1F*)gDirectory->Get(hDenName));
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
         TString gEffName = "gEff_" + hNumName + "_vs_" + hDenName;
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
            //gEff->SetTitle(gEffTitles[i]);
            gEff->DrawClone("p");
            legend->AddEntry(gEff, legTitles[fInd], "l");
         }
         hDiv.front()->Draw("sameaxis");
         TLine* l = new TLine(7, minimum, 7, maximum);
         l->Draw("same");
         legend->Draw("same");
   
         TLatex *tex = new TLatex();
         tex->SetNDC();
         tex->SetTextFont(font);
         tex->SetLineWidth(2);
         tex->SetTextSize(0.035);
         //tex->DrawLatex(0.02, 0.96, gEffTitles[i]);
         tex->SetTextSize(0.03);
         tex->DrawLatex(0.10, 0.91, "CMS Preliminary");
         tex->DrawLatex(0.62, 0.91, "Errors: Clopper-Pearson");

         if (saveAsPng) {
            c1->Print(plotsDir + gEffName + fileNameExtra + ".png", "png");
         }
      }
   }

   if (plotEtEff) {
      for (unsigned int probeInd = 0; probeInd < gEffHistNamesPrefNum.size(); ++probeInd) {
         for (unsigned int i = 0; i < pathNames.size(); ++i) {
         //for (unsigned int i = 0; i < 1; ++i) {
            TString hNumName = gEffHistNamesPrefNum[probeInd] + "_" + pathNames[i];
            TString hDenName = gEffHistNamesPrefDen[probeInd] + "_" + pathNames[i];

            TCanvas *c1;
            TPad *gEffPad;
            c1 = new TCanvas("c1_" + hNumName + "_vs_" + hDenName, "c1_" + hNumName + "_vs_" + hDenName, 100, 100, 800, 600);
            gEffPad = new TPad("gEffPad" + hNumName + "_vs_" + hDenName, "gEffPad" + hNumName + "_vs_" + hDenName, 0., 0., 1., 1.);
            gEffPad->SetBorderMode(0);
            gEffPad->SetFrameBorderMode(0);
            gEffPad->SetFillColor(0);
            gEffPad->SetFrameFillColor(0);
            gEffPad->SetRightMargin(0.16);
            gEffPad->SetBottomMargin(0.18);
            gEffPad->Draw();
            gEffPad->cd();
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
               h.push_back((TH1F*)gDirectory->Get(hNumName));
               hTot.push_back((TH1F*)gDirectory->Get(hDenName));
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
            TString gEffName = "gEff_" + hNumName + "_vs_" + hDenName;
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
               //gEff->SetTitle(gEffTitles[i]);
               gEff->DrawClone("p");
               legend->AddEntry(gEff, legTitles[fInd], "l");
            }
            hDiv.front()->Draw("sameaxis");
            TLine* l = new TLine(7, minimum, 7, maximum);
            l->Draw("same");
            legend->Draw("same");
   
            TLatex *tex = new TLatex();
            tex->SetNDC();
            tex->SetTextFont(font);
            tex->SetLineWidth(2);
            tex->SetTextSize(0.035);
            //tex->DrawLatex(0.02, 0.96, gEffTitles[i]);
            tex->SetTextSize(0.03);
            tex->DrawLatex(0.10, 0.91, "CMS Preliminary");
            tex->DrawLatex(0.32, 0.91, pathNames[i]);
            tex->DrawLatex(0.62, 0.91, "Errors: Clopper-Pearson");

            if (saveAsPng) {
               c1->Print(plotsDir + gEffName + fileNameExtra + ".png", "png");
            }
         }
      }
   }
}
