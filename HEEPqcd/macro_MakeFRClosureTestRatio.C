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

void macro_MakeFRClosureTestRatio()
{
  // parameters //////////////////////////////////////////////////////////////
  // luminosity of data
  const float lumi = 2928;

  string inputFilePrefix = "invMassHistos";  // lumi and .root will be added
  //string inputFilePrefix = "test";

  // which closure test ratios should be plotted?
  bool plotClosureTest[4];
  plotClosureTest[0] = false;  // no corrections
  plotClosureTest[1] = false;  // GSF electrons not passing HEEP
  plotClosureTest[2] = false;  // all above + HEEP-GSF corrected with DY contribuion
  plotClosureTest[3] = true;  // all above + GSF-GSF corrected with W+jet and gamma+jet contribution

  // which histograms should be plotted?
  bool plotHisto[4];
  plotHisto[0] = true;  // EB-EB + EB-EE
  plotHisto[1] = true;  // EB-EB
  plotHisto[2] = true;  // EB-EE
  plotHisto[3] = true;  // EE-EE
  ////////////////////////////////////////////////////////////////////////////

//  const int rebin = 10;
  float massMin = 50;
  float massMax = 1500;
  int nBins = 1450;
  vector<pair<float, float> > binning;
  // VARIABLE BINNING
  binning.push_back(make_pair(100, 10));
  binning.push_back(make_pair(500, 25));
  binning.push_back(make_pair(1500, 250));
  // CONSTANT BINNING
//  binning.push_back(make_pair(massMax, 50));

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

  stringstream sStream;
  //sStream << "Photon-Run2011A-May10ReReco-v1+05Aug2011-v1+Run2011A-PromptReco-v4+PromptReco-v6-AOD-Cert_160404-173692_7TeV_Collisions11_JSON_" << lumi << "pb-1";
  sStream << lumi << "pb-1";
  TString folderDataHisto = sStream.str().c_str();
  // select histograms dynamically depending on state of correction //////////
  vector<TString> folderHeepGsfHisto;
  folderHeepGsfHisto.push_back(sStream.str().c_str());
  folderHeepGsfHisto.push_back(sStream.str().c_str());
  folderHeepGsfHisto.push_back("combinations");
  folderHeepGsfHisto.push_back("combinations");

  vector<TString> folderGsfGsfHisto;
  folderGsfGsfHisto.push_back(sStream.str().c_str());
  folderGsfGsfHisto.push_back(sStream.str().c_str());
  folderGsfGsfHisto.push_back(sStream.str().c_str());
  folderGsfGsfHisto.push_back("combinations");

  vector<TString> heepGsfHisto;
  heepGsfHisto.push_back("histoHeepGsfMassFR");
  heepGsfHisto.push_back("histoHeepGsfMassNoHeepFR");
  heepGsfHisto.push_back("histoHeepGsfCorr");
  heepGsfHisto.push_back("histoHeepGsfCorr");

  vector<TString> gsfGsfHisto;
  gsfGsfHisto.push_back("histoGsfGsfMassFR");
  gsfGsfHisto.push_back("histoGsfGsfMassNoHeepFR");
  gsfGsfHisto.push_back("histoGsfGsfMassNoHeepFR");
  gsfGsfHisto.push_back("histoGsfGsfCorr");
  ////////////////////////////////////////////////////////////////////////////

  vector<TString> canvasName;
  canvasName.push_back("ratioFR");
  canvasName.push_back("ratioFRNoHeep");
  canvasName.push_back("ratioFRCorrDY");
  canvasName.push_back("ratioFRCorrFull");

  vector<TString> canvasTitle;
  canvasTitle.push_back("Fake rate ratio GSF-GSF / HEEP-GSF - uncorrected");
  canvasTitle.push_back("Fake rate ratio GSF-GSF / HEEP-GSF - non HEEP");
  canvasTitle.push_back("Fake rate ratio GSF-GSF / HEEP-GSF - DY corrected");
  canvasTitle.push_back("Fake rate ratio GSF-GSF / HEEP-GSF - full corrected");

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

  sStream.str("");
  sStream << inputFilePrefix << lumi << "pb-1.root";
  TFile input(sStream.str().c_str(), "read");
  input.cd();

  cout << endl << "Input file: " << sStream.str() << endl;

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);
  
  for (unsigned int j = 0; j < 4; ++j) {
    if (!plotClosureTest[j]) continue;
    for (unsigned int p = 0; p < 4; ++p) {
      if (!plotHisto[p]) continue;
  
      TCanvas *c0 = new TCanvas(canvasName[j] + acroSuffix[p], canvasTitle[j] + suffix[p], 100, 100, 800, 600);
      c0->cd();
      c0->SetBorderMode(0);
      c0->SetFrameBorderMode(0);
      c0->SetFillColor(0);
      c0->SetFrameFillColor(0);
      gStyle->SetOptStat(0);
      gStyle->SetPadTickY(1);
  
      // get the histograms 
      input.cd(folderGsfGsfHisto[j]);
      TH1F *numHisto = (TH1F *)gDirectory->Get(gsfGsfHisto[j] + acroSuffix[p]);
      input.cd(folderHeepGsfHisto[j]);
      TH1F *denomHisto = (TH1F *)gDirectory->Get(heepGsfHisto[j] + acroSuffix[p]);
  
      TH1F *numHistoRebinned = (TH1F *)numHisto->Rebin(nBins, "numHistoRebinned" + acroSuffix[p], binArray);
      TH1F *denomHistoRebinned = (TH1F *)denomHisto->Rebin(nBins, "denomHistoRebinned" + acroSuffix[p], binArray);
  
      TH1F *ratioHisto = new TH1F("histoRatioCorr" + acroSuffix[p], canvasTitle[j] + suffix[p], nBins, binArray);
      ratioHisto->Divide(numHistoRebinned, denomHistoRebinned);
  
      ratioHisto->SetLineColor(4);
      ratioHisto->SetMarkerColor(4);
      ratioHisto->SetMarkerStyle(20);
  
      ratioHisto->Draw();
      ratioHisto->Fit("pol0", "+", "lep", 120, 500);
  
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
  
      //TLegend *legend = new TLegend(0.38, 0.6, 0.53, 0.9);
      //legend->SetTextSize(0.03);
      //legend->SetBorderSize(0);
      //legend->SetFillStyle(0);
      //legend->AddEntry(ratioHisto, "Ratio GSF-GSF / HEEP-GSF (non HEEP)", "lep");
      //legend->Draw("sames");

    } // end loop over eta ranges
  } // end loop over corrections
  input.Close();
}
