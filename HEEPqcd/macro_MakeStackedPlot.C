#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>

#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TColor.h"
#include "TH1F.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include <TCanvas.h>
#include <TStyle.h>

void macro_MakeStackedPlot()
{
  // parameters //////////////////////////////////////////////////////////////
  // luminosity of data
  const float lumi = 14818.;

  string inputFilePrefix = "invMassHistos";  // lumi and .root will be added
  
  // which plot type (cumulated or not) should be plotted?
  bool plotType[2];
  plotType[0] = true;  // Invariant mass plots
  plotType[1] = true;  // Integrated invariant mass plots

  // which histograms should be plotted?
  bool plotHisto[4];
  plotHisto[0] = true;  // EB-EB + EB-EE
  plotHisto[1] = true;  // EB-EB
  plotHisto[2] = true;  // EB-EE
  plotHisto[3] = true;  // EE-EE

  // plot style
  int ttbarColour = TColor::GetColor("#ff6666");
  int zttColour = TColor::GetColor("#ff4d4d");
  int wwColour = TColor::GetColor("#ff3333");
  int wzColour = TColor::GetColor("#ff0f0f");
  int twColour = TColor::GetColor("#eb0000");
  int wjetColour=  TColor::GetColor("#66b3ff");
  int zmmColour=  TColor::GetColor("#80bfff");
  int zeeColour=  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66");

  //int ttbarColour = kRed;
  //int zttColour = kRed-7;
  //int wwColour = kRed-4;
  //int wzColour = kRed+1;
  //int twColour = kRed+2;
  //int wjetColour = kCyan-2;
  //int zmmColour = kCyan-1;
  //int zeeColour = kCyan;
  //int jetBfgColour = kYellow;

  int font = 42; //62
  ////////////////////////////////////////////////////////////////////////////

//  const int rebin = 10;
  float massMin = 50;
  float massMax = 2050;
  int nBins = 200;
  vector<pair<float, float> > binning;
  // VARIABLE BINNING
//  binning.push_back(make_pair(100, 1));
//  binning.push_back(make_pair(500, 10));
//  binning.push_back(make_pair(1500, 50));
  // CONSTANT BINNING
  binning.push_back(make_pair(massMax, 10));

  vector<float> bins;
  bins.push_back(massMin);
  for (vector<pair<float,float> >::iterator it = binning.begin(); it < binning.end(); ++it) {
    while (bins.back() < it->first)
      bins.push_back(bins.back() + it->second);
  }
  if (bins.back() < massMax)
    bins.push_back(massMax);
  nBins = bins.size() - 1;
  Double_t binArray[nBins + 1];
  for (int i = 0; i <= nBins; ++i)
    binArray[i] = (Double_t)bins.at(i);

  vector<TString> acroSuffix;
  acroSuffix.push_back("");
  acroSuffix.push_back("BB");
  acroSuffix.push_back("BE");
  acroSuffix.push_back("EE");

  vector<TString> suffix;
  suffix.push_back("");
  suffix.push_back(" EB-EB");
  suffix.push_back(" EB-EE");
  suffix.push_back(" EE-EE");

  stringstream sStream;
  sStream << "invMassHistos" << lumi << "pb-1.root";
  //sStream << inputFilePrefix << lumi << "pb-1.root";
  TFile input(sStream.str().c_str(), "read");
  input.cd();

  cout << endl << "Input file: " << sStream.str() << endl;
  //vector<THStack *> histoStacks;

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  for (unsigned int j = 0; j < 2; ++j) {
    if (!plotType[j]) continue;
    for (unsigned int p = 0; p < 4; ++p) {
      if (!plotHisto[p]) continue;

      // get the histograms 
      input.cd("combinations");
      TH1F *histoBgQCD = (TH1F *)gDirectory->Get("histoGsfGsfCorr" + acroSuffix[p]);
      TH1F *histoBgDY = (TH1F *)gDirectory->Get("histoDYCombined" + acroSuffix[p]);
      input.cd("WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_45_76102995ev");
      TH1F *histoBgWjet = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      input.cd("TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_45_28150723ev");
      TH1F *histoBgTtbar = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      input.cd("DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_3295238ev");
      TH1F *histoBgDYTT = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      input.cd("WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_10000431ev");
      TH1F *histoBgWW = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      input.cd("WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_10000283ev");
      TH1F *histoBgWZ = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      input.cd("T+Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_991118ev");
      TH1F *histoBgTW = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      sStream.str("");
      sStream << "Photon_Run2012A+DoublePhotonHighPt_Run2012B+C+D_13Jul2012+06Aug2012+PromptReco-v1+v2_Cert_190456-206098_gct1_41_14818pb-1";
      //sStream << lumi << "pb-1";
      input.cd(sStream.str().c_str());
      TH1F *histoData = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);

      if (j == 1) {
        acroSuffix[p] = "Cumulative" + acroSuffix[p];
        suffix[p] = ", Cumulative Plot"  + suffix[p];
      }

      //cout << histoData->Integral(histoData->FindBin(60), histoData->FindBin(120)) << endl;
      //cout << histoBgDY->Integral(histoBgDY->FindBin(60), histoBgDY->FindBin(120)) << endl;
      //cout << histoBgQCD->Integral(histoBgQCD->FindBin(60), histoBgQCD->FindBin(120)) << endl;
      //cout << histoBgWjet->Integral(histoBgWjet->FindBin(60), histoBgWjet->FindBin(120)) << endl;
      //cout << histoBgTtbar->Integral(histoBgTtbar->FindBin(60), histoBgTtbar->FindBin(120)) << endl << endl;

      TCanvas *c0 = new TCanvas("invMassStack" + acroSuffix[p], "Invariant Mass" + suffix[p], 100, 100, 800, 600);
      c0->cd();
      c0->SetBorderMode(0);
      c0->SetFrameBorderMode(0);
      c0->SetFillColor(0);
      c0->SetFrameFillColor(0);
      c0->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);

      THStack *histoStack = new THStack("bgStack" + acroSuffix[p], "Invariant Mass" + suffix[p]);
      //histoStacks.push_back(new THStack("bgStack" + acroSuffix[p], "Invariant Mass" + suffix[p]));

      // rebin
      TH1F *histoBgQCDRebinned = (TH1F *)histoBgQCD->Rebin(nBins, "histoBgQCDRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgDYRebinned = (TH1F *)histoBgDY->Rebin(nBins, "histoBgDYRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgWjetRebinned = (TH1F *)histoBgWjet->Rebin(nBins, "histoBgWjetRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgTtbarRebinned = (TH1F *)histoBgTtbar->Rebin(nBins, "histoBgTtbarRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgDYTTRebinned = (TH1F *)histoBgDYTT->Rebin(nBins, "histoBgDYTTRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgWWRebinned = (TH1F *)histoBgWW->Rebin(nBins, "histoBgWWRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgWZRebinned = (TH1F *)histoBgWZ->Rebin(nBins, "histoBgWZRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgTWRebinned = (TH1F *)histoBgTW->Rebin(nBins, "histoBgTWRebinned" + acroSuffix[p], binArray);
      TH1F *histoDataRebinned = (TH1F *)histoData->Rebin(nBins, "histoDataRebinned" + acroSuffix[p], binArray);

      histoBgQCDRebinned->SetFillColor(jetBkgColour);
      histoBgQCDRebinned->SetLineColor(kBlack);
      histoBgQCDRebinned->SetLineWidth(2);
      histoBgWjetRebinned->SetFillColor(wjetColour);
      histoBgWjetRebinned->SetLineColor(kBlack);
      histoBgWjetRebinned->SetLineWidth(2);
      histoBgTtbarRebinned->SetFillColor(ttbarColour);
      histoBgTtbarRebinned->SetLineColor(kBlack);
      histoBgTtbarRebinned->SetLineWidth(2);
      histoBgDYTTRebinned->SetFillColor(zttColour);
      histoBgDYTTRebinned->SetLineColor(kBlack);
      histoBgDYTTRebinned->SetLineWidth(2);
      histoBgWWRebinned->SetFillColor(wwColour);
      histoBgWWRebinned->SetLineColor(kBlack);
      histoBgWWRebinned->SetLineWidth(2);
      histoBgWZRebinned->SetFillColor(wzColour);
      histoBgWZRebinned->SetLineColor(kBlack);
      histoBgWZRebinned->SetLineWidth(2);
      histoBgTWRebinned->SetFillColor(twColour);
      histoBgTWRebinned->SetLineColor(kBlack);
      histoBgTWRebinned->SetLineWidth(2);
      histoBgDYRebinned->SetFillColor(zeeColour);
      histoBgDYRebinned->SetLineColor(kBlack);
      histoBgDYRebinned->SetLineWidth(2);

      histoDataRebinned->SetLineColor(kBlack);
      histoDataRebinned->SetFillColor(0);
      histoDataRebinned->SetMarkerColor(kBlack);
      histoDataRebinned->SetMarkerStyle(20);
      histoDataRebinned->SetMarkerSize(1.1);

      // integrate (from high mass) if flag is set
      if (j == 1) {
        double error;
        for (int i = histoDataRebinned->GetNbinsX(); i > 0; i--) {
          histoDataRebinned->SetBinContent(i, histoDataRebinned->IntegralAndError(i, i + 1, error));
          histoDataRebinned->SetBinError(i, error);
        }
        for (int i = histoBgQCDRebinned->GetNbinsX(); i > 0; i--) {
          histoBgQCDRebinned->SetBinContent(i, histoBgQCDRebinned->IntegralAndError(i, i + 1, error));
          histoBgQCDRebinned->SetBinError(i, error);
        }
        for (int i = histoBgWjetRebinned->GetNbinsX(); i > 0; i--) {
          histoBgWjetRebinned->SetBinContent(i, histoBgWjetRebinned->IntegralAndError(i, i + 1, error));
          histoBgWjetRebinned->SetBinError(i, error);
        }
        for (int i = histoBgTtbarRebinned->GetNbinsX(); i > 0; i--) {
          histoBgTtbarRebinned->SetBinContent(i, histoBgTtbarRebinned->IntegralAndError(i, i + 1, error));
          histoBgTtbarRebinned->SetBinError(i, error);
        }
        for (int i = histoBgDYTTRebinned->GetNbinsX(); i > 0; i--) {
          histoBgDYTTRebinned->SetBinContent(i, histoBgDYTTRebinned->IntegralAndError(i, i + 1, error));
          histoBgDYTTRebinned->SetBinError(i, error);
        }
        for (int i = histoBgWWRebinned->GetNbinsX(); i > 0; i--) {
          histoBgWWRebinned->SetBinContent(i, histoBgWWRebinned->IntegralAndError(i, i + 1, error));
          histoBgWWRebinned->SetBinError(i, error);
        }
        for (int i = histoBgWZRebinned->GetNbinsX(); i > 0; i--) {
          histoBgWZRebinned->SetBinContent(i, histoBgWZRebinned->IntegralAndError(i, i + 1, error));
          histoBgWZRebinned->SetBinError(i, error);
        }
        for (int i = histoBgTWRebinned->GetNbinsX(); i > 0; i--) {
          histoBgTWRebinned->SetBinContent(i, histoBgTWRebinned->IntegralAndError(i, i + 1, error));
          histoBgTWRebinned->SetBinError(i, error);
        }
        for (int i = histoBgDYRebinned->GetNbinsX(); i > 0; i--) {
          histoBgDYRebinned->SetBinContent(i, histoBgDYRebinned->IntegralAndError(i, i + 1, error));
          histoBgDYRebinned->SetBinError(i, error);
        }
      }
      histoStack->Add(histoBgQCDRebinned);
      histoStack->Add(histoBgWjetRebinned);
      histoStack->Add(histoBgTWRebinned);
      histoStack->Add(histoBgWZRebinned);
      histoStack->Add(histoBgWWRebinned);
      histoStack->Add(histoBgDYTTRebinned);
      histoStack->Add(histoBgTtbarRebinned);
      histoStack->Add(histoBgDYRebinned);
       
      histoStack->Draw("hist");
      histoStack->GetXaxis()->SetTitle(histoDataRebinned->GetXaxis()->GetTitle());
      if (j == 1) {
        histoStack->GetYaxis()->SetTitle("Events >= m_{ee}");
      } else {
        sStream.str("");
        if (binning.size() > 1)
          sStream << "Events / GeV";
        else
          sStream << "Events / " << binning.begin()->second << "GeV";
        histoStack->GetYaxis()->SetTitle(sStream.str().c_str());
      }
    
      sStream.str("");
      sStream << "#sqrt{s} = 8TeV,  #int L dt = " << lumi << "pb^{-1}";
      TPaveLabel *label0 = new TPaveLabel(0.6, 0.7, 0.9, 0.8, sStream.str().c_str(), "brNDC");
      label0->SetFillColor(0);
      label0->SetFillStyle(0);
      label0->SetBorderSize(0);
      label0->SetTextSize(0.30);
      label0->SetTextFont(font);
      label0->Draw("sames");
      TPaveLabel *label1 = new TPaveLabel(0.6, 0.8, 0.9, 0.9, "CMS preliminary", "brNDC");
      label1->SetFillColor(0);
      label1->SetFillStyle(0);
      label1->SetBorderSize(0);
      label1->SetTextSize(0.40);
      label1->SetTextFont(font);
      label1->Draw("sames");
 
      TLegend *legend = new TLegend(0.38, 0.6, 0.53, 0.9);
      legend->SetTextSize(0.03);
      legend->SetTextFont(font);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
    
      legend->AddEntry(histoDataRebinned, "DATA", "lep");
      legend->AddEntry(histoBgDYRebinned, "Z/#gamma* #rightarrow ee", "f");
      legend->AddEntry(histoBgTtbarRebinned, "t#bar{t}", "f");
      legend->AddEntry(histoBgDYTTRebinned, "Z/#gamma* #rightarrow #tau#tau", "f");
      legend->AddEntry(histoBgWWRebinned, "WW", "f");
      legend->AddEntry(histoBgWZRebinned, "WZ", "f");
      legend->AddEntry(histoBgTWRebinned, "tW & #bar{t}W", "f");
      legend->AddEntry(histoBgWjetRebinned, "W + jets", "f");
      legend->AddEntry(histoBgQCDRebinned, "QCD", "f");
      legend->Draw("sames");
 
      histoDataRebinned->Draw("sames");

    } // end loop over eta regions
  } // end loop over plotTypes
  input.Close();
}
