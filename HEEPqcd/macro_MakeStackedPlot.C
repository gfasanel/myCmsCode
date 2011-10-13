#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>

#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
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
  const float lumi = 2928;

  string inputFilePrefix = "invMassHistos";  // lumi and .root will be added
  //string inputFilePrefix = "test";

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
  ////////////////////////////////////////////////////////////////////////////

//  const int rebin = 10;
  float massMin = 50;
  float massMax = 1550;
  int nBins = 150;
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
  //sStream << "invMassHistos" << lumi << "pb-1.root";
  //sStream << "invMassHistos" << lumi << "pb-1_save.root";
  sStream << inputFilePrefix << lumi << "pb-1.root";
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
      input.cd("TOTAL");
      TH1F *histoBgWjet = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      input.cd("TT_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_AODSIM_HEEPSkim2ElePt30");
      TH1F *histoBgTtbar = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);
      sStream.str("");
      //sStream << "Photon-Run2011A-May10ReReco-v1+05Aug2011-v1+Run2011A-PromptReco-v4+PromptReco-v6-AOD-Cert_160404-173664_7TeV_Collisions11_JSON_" << lumi << "pb-1";
      //sStream << "Photon-Run2011A-May10ReReco-v1+05Aug2011-v1+Run2011A-PromptReco-v4+PromptReco-v6-AOD-Cert_160404-173692_7TeV_Collisions11_JSON_" << lumi << "pb-1";
      sStream << lumi << "pb-1";
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
      gStyle->SetPadTickY(1);

      THStack *histoStack = new THStack("bgStack" + acroSuffix[p], "Invariant Mass" + suffix[p]);
      //histoStacks.push_back(new THStack("bgStack" + acroSuffix[p], "Invariant Mass" + suffix[p]));

      // rebin
      TH1F *histoBgQCDRebinned = (TH1F *)histoBgQCD->Rebin(nBins, "histoBgQCDRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgDYRebinned = (TH1F *)histoBgDY->Rebin(nBins, "histoBgDYRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgWjetRebinned = (TH1F *)histoBgWjet->Rebin(nBins, "histoBgWjetRebinned" + acroSuffix[p], binArray);
      TH1F *histoBgTtbarRebinned = (TH1F *)histoBgTtbar->Rebin(nBins, "histoBgTtbarRebinned" + acroSuffix[p], binArray);
      TH1F *histoDataRebinned = (TH1F *)histoData->Rebin(nBins, "histoDataRebinned" + acroSuffix[p], binArray);

      histoBgQCDRebinned->SetFillColor(5);
      histoBgWjetRebinned->SetFillColor(4);
      histoBgTtbarRebinned->SetFillColor(3);
      histoBgDYRebinned->SetFillColor(2);

      histoDataRebinned->SetLineColor(1);
      histoDataRebinned->SetFillColor(0);
      histoDataRebinned->SetMarkerColor(1);
      histoDataRebinned->SetMarkerStyle(20);

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
        for (int i = histoBgDYRebinned->GetNbinsX(); i > 0; i--) {
          histoBgDYRebinned->SetBinContent(i, histoBgDYRebinned->IntegralAndError(i, i + 1, error));
          histoBgDYRebinned->SetBinError(i, error);
        }
      }
      histoStack->Add(histoBgQCDRebinned);
      histoStack->Add(histoBgWjetRebinned);
      histoStack->Add(histoBgTtbarRebinned);
      histoStack->Add(histoBgDYRebinned);
       
      histoStack->Draw("hist");
      histoStack->GetXaxis()->SetTitle(histoDataRebinned->GetXaxis()->GetTitle());
      if (j == 1) {
        histoStack->GetYaxis()->SetTitle("# of events >= M_{ee}");
      } else {
        sStream.str("");
        if (binning.size() > 1)
          sStream << "# of events / bin";
        else
          sStream << "# of events / " << binning.begin()->second << "GeV/c^{2}";
        histoStack->GetYaxis()->SetTitle(sStream.str().c_str());
      }
    
      sStream.str("");
      sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
      TPaveLabel *label0 = new TPaveLabel(0.6, 0.7, 0.9, 0.8, sStream.str().c_str(), "brNDC");
      label0->SetFillColor(0);
      label0->SetFillStyle(0);
      label0->SetBorderSize(0);
      label0->SetTextSize(0.30);
      label0->Draw("sames");
      TPaveLabel *label1 = new TPaveLabel(0.6, 0.8, 0.9, 0.9, "CMS preliminary", "brNDC");
      label1->SetFillColor(0);
      label1->SetFillStyle(0);
      label1->SetBorderSize(0);
      label1->SetTextSize(0.40);
      label1->Draw("sames");
 
      TLegend *legend = new TLegend(0.38, 0.6, 0.53, 0.9);
      legend->SetTextSize(0.03);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
    
      legend->AddEntry(histoDataRebinned, "Data", "lep");
      legend->AddEntry(histoBgDYRebinned, "DY", "f");
      legend->AddEntry(histoBgTtbarRebinned, "ttbar", "f");
      legend->AddEntry(histoBgWjetRebinned, "W + jets", "f");
      legend->AddEntry(histoBgQCDRebinned, "QCD bg", "f");
      legend->Draw("sames");
 
      histoDataRebinned->Draw("sames");

    } // end loop over eta regions
  } // end loop over plotTypes
  input.Close();
}
