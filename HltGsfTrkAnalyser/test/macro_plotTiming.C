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

void macro_plotTiming()
{
   TH1::AddDirectory(kFALSE);

   const bool plotDoubleEleTrg = 1;

   int font = 42;
   const bool saveAsPng = 0;
   TString plotsDir = "plots_20130917/";

   TString fileDir("/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/");
   vector<TString> files;
   files.push_back("./DQM_V0001_R000204113__HLT__FastTimerService__All_vocms110.root");
   TString dir("DQMData/Run 204113/HLT/Run summary/TimerService/");
   TString subDir("Paths/");

   vector<TString> trgNames;
   if (plotDoubleEleTrg) {
      trgNames.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7");
      trgNames.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC5");
      trgNames.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC4");
      trgNames.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC3");
      trgNames.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC2");
   } else {
      trgNames.push_back("HLT_Ele80_CaloIdVT_GsfTrkIdT_v2");
      trgNames.push_back("HLT_Ele80_CaloIdVT_GsfTrkIdT_nC5");
      trgNames.push_back("HLT_Ele80_CaloIdVT_GsfTrkIdT_nC4");
      trgNames.push_back("HLT_Ele80_CaloIdVT_GsfTrkIdT_nC3");
      trgNames.push_back("HLT_Ele80_CaloIdVT_GsfTrkIdT_nC2");
   }

   vector<TString> moduleNames;
   if (plotDoubleEleTrg) {
      moduleNames.push_back("hltActivityElectronGsfTracks");
      moduleNames.push_back("hltActivityElectronGsfTracksNC5");
      moduleNames.push_back("hltActivityElectronGsfTracksNC4");
      moduleNames.push_back("hltActivityElectronGsfTracksNC3");
      moduleNames.push_back("hltActivityElectronGsfTracksNC2");
   } else {
      moduleNames.push_back("hltL1SeededElectronGsfTracks");
      moduleNames.push_back("hltL1SeededElectronGsfTracksNC5");
      moduleNames.push_back("hltL1SeededElectronGsfTracksNC4");
      moduleNames.push_back("hltL1SeededElectronGsfTracksNC3");
      moduleNames.push_back("hltL1SeededElectronGsfTracksNC2");
   }

   vector<TString> legEntries;
   legEntries.push_back("NC6");
   legEntries.push_back("NC5");
   legEntries.push_back("NC4");
   legEntries.push_back("NC3");
   legEntries.push_back("NC2");

   vector<TString> histoNames;
   histoNames.push_back("_total");

   int colours[] = {kAzure -5, kRed, kGreen, kMagenta, kOrange, kCyan};

   TFile inFile(fileDir + files[0], "read");
   inFile.cd(dir + subDir);
   for (unsigned int i = 0; i < histoNames.size(); ++i) {
   //for (unsigned int i = 0; i < 3; ++i) {
      TCanvas *c0 = new TCanvas("c0_" + histoNames[i], "c0_" + histoNames[i], 100, 100, 800, 600);
      c0->cd();
      c0->SetBorderMode(0);
      c0->SetFrameBorderMode(0);
      c0->SetFillColor(0);
      c0->SetFrameFillColor(0);
      //c0->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
      TLegend *legend = new TLegend(0.63, 0.63, 0.75, 0.85);
      legend->SetTextSize(0.03);
      legend->SetTextFont(font);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      vector<TH1F*> h;
      for (unsigned int j = 0; j < trgNames.size(); ++j) {
         inFile.cd(dir + subDir);
         h.push_back((TH1F *)gDirectory->Get(trgNames[j] + histoNames[i]));

         h.back()->SetLineColor(colours[j]);
         if (j == 0) h.back()->SetFillColor(colours[j]);
         h.back()->SetLineWidth(2);
         h.back()->Rebin(10);
         if (j == 0) h.back()->Draw("h");
         else h.back()->Draw("hsame");
         if (j == 0) legend->AddEntry(h.back(), legEntries[j] + Form(", Mean: %.1f ms", h.back()->GetMean()), "f");
         else legend->AddEntry(h.back(), legEntries[j] + Form(", Mean: %.1f ms", h.back()->GetMean()), "l");
      }

      legend->Draw("same");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.03);
      if (plotDoubleEleTrg) tex->DrawLatex(0.10, 0.91, "Paths: HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_*");
      else tex->DrawLatex(0.10, 0.91, "Paths: HLT_Ele80_CaloIdVT_GsfTrkIdT_*");
      tex->DrawLatex(0.75, 0.91, "CMS Preliminary");

      if (plotDoubleEleTrg) {
         if (saveAsPng) c0->Print(plotsDir + "hTrgTime_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7" + histoNames[i] + ".png", "png");
      } else {
         if (saveAsPng) c0->Print(plotsDir + "hTrgTime_HLT_Ele80_CaloIdVT_GsfTrkIdT" + histoNames[i] + ".png", "png");
      }
   }

   {
      TCanvas *c1 = new TCanvas("c1_", "c1_", 100, 100, 800, 600);
      c1->cd();
      c1->SetBorderMode(0);
      c1->SetFrameBorderMode(0);
      c1->SetFillColor(0);
      c1->SetFrameFillColor(0);
      //c1->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      TLegend *legend = new TLegend(0.63, 0.63, 0.75, 0.85);
      legend->SetTextSize(0.03);
      legend->SetTextFont(font);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);

      int rebin = 5;
      double max = 0.;
      vector<TH1F*> hMod;
      for (unsigned int j = 0; j < moduleNames.size(); ++j) {
         inFile.cd(dir + "Modules/");
         hMod.push_back((TH1F *)gDirectory->Get(moduleNames[j]));
         double hMax = hMod.back()->GetMaximum();
         if (max < hMax) max = hMax;
         hMod.back()->GetXaxis()->SetTitle(hMod.back()->GetYaxis()->GetTitle());
         hMod.back()->GetYaxis()->SetTitle("");
         hMod.back()->SetLineColor(colours[j]);
         if (j == 0) hMod.back()->SetFillColor(colours[j]);
         hMod.back()->SetLineWidth(2);
         hMod.back()->Rebin(rebin);
      }
      hMod.at(0)->SetMaximum(max * rebin * 0.75);
      for (unsigned int j = 0; j < moduleNames.size(); ++j) {
         if (j == 0) {
            hMod[j]->Draw("h");
            legend->AddEntry(hMod[j], legEntries[j] + Form(", Mean: %.1f ms", hMod[j]->GetMean()), "f");
         } else {
            hMod[j]->Draw("hsame");
            legend->AddEntry(hMod[j], legEntries[j] + Form(", Mean: %.1f ms", hMod[j]->GetMean()), "l");
         }
      }
      hMod[0]->Draw("sameaxis");

      legend->Draw("same");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.03);
      if (plotDoubleEleTrg) tex->DrawLatex(0.10, 0.91, "Modules: hltActivityElectronGsfTracks*");
      else tex->DrawLatex(0.10, 0.91, "Modules: hltL1SeededElectronGsfTracks*");
      tex->DrawLatex(0.75, 0.91, "CMS Preliminary");

      if (saveAsPng) c1->Print(plotsDir + "hModTime_" + moduleNames[0] + ".png", "png");
   }
   {
      TCanvas *c2 = new TCanvas("c2_", "c2_", 100, 100, 800, 600);
      c2->cd();
      c2->SetBorderMode(0);
      c2->SetFrameBorderMode(0);
      c2->SetFillColor(0);
      c2->SetFrameFillColor(0);
      //c2->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      inFile.cd(dir);
      TH1F* hAll;
      hAll = (TH1F *)gDirectory->Get("all_paths");
      hAll->SetLineColor(colours[0]);
      hAll->SetFillColor(colours[0]);
      hAll->SetLineWidth(2);
      //hAll->Rebin(4);
      hAll->Draw("h");

      TLegend *legend = new TLegend(0.63, 0.80, 0.75, 0.85);
      legend->SetTextSize(0.03);
      legend->SetTextFont(font);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->AddEntry(hAll, "total", "f");
      legend->Draw("same");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.03);
      tex->DrawLatex(0.10, 0.91, "Menu time");
      tex->DrawLatex(0.75, 0.91, "CMS Preliminary");

      if (saveAsPng) c2->Print(plotsDir + "hMenuTime.png", "png");
   }
   {
      TCanvas *c3 = new TCanvas("c3_", "c3_", 100, 100, 800, 600);
      c3->cd();
      c3->SetBorderMode(0);
      c3->SetFrameBorderMode(0);
      c3->SetFillColor(0);
      c3->SetFrameFillColor(0);
      c3->SetBottomMargin(0.17);
      //c3->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      inFile.cd(dir);
      TProfile* pAll;
      pAll = (TProfile*)gDirectory->Get("paths_total_time");
      pAll->SetLineColor(colours[0]);
      pAll->SetLineWidth(2);
      pAll->Draw();

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      tex->SetTextSize(0.03);
      tex->DrawLatex(0.75, 0.91, "CMS Preliminary");

      if (saveAsPng) c3->Print(plotsDir + "hAllPathsTotalTime.png", "png");
   }
}
