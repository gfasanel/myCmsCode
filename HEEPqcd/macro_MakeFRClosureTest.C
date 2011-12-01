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

void macro_MakeFRClosureTest()
{
  // parameters //////////////////////////////////////////////////////////////
  // luminosity of data
  const float lumi = 4679;

  string inputFilePrefix = "invMassHistos";  // lumi and .root will be added
  //string inputFilePrefix = "test";

  // plot data?
  const bool plotData = true;
  //const bool plotData = false;

  // which closure tests should be plotted?
  bool plotClosureTest[4];
  plotClosureTest[0] = true;  // no corrections
  plotClosureTest[1] = true;  // GSF electrons not passing HEEP
  plotClosureTest[2] = true;  // all above + HEEP-GSF corrected with DY contribuion
  plotClosureTest[3] = true;  // all above + GSF-GSF corrected with W+jet and gamma+jet contribution

  // which histograms should be plotted?
  bool plotHisto[4];
  plotHisto[0] = true;  // EB-EB + EB-EE
  plotHisto[1] = true;  // EB-EB
  plotHisto[2] = true;  // EB-EE
  plotHisto[3] = true;  // EE-EE
  ////////////////////////////////////////////////////////////////////////////

  const double yAxisMin = 0.001;
//  const int rebin = 10;
  float massMin = 50;
  float massMax = 1550;
  int nBins = 1500;
  vector<pair<float, float> > binning;
  // VARIABLE BINNING
//  binning.push_back(make_pair(100, 1));
//  binning.push_back(make_pair(500, 10));
//  binning.push_back(make_pair(massMax, 50));
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

  stringstream sStream;
  sStream << "Photon-Run2011A-May10ReReco-v1+05Aug2011-v1+Run2011A-PromptReco-v4+PromptReco-v6+Run2011B-PromptReco-v1-AOD-Cert_160404-180252_7TeV_Collisions11_JSON_" << lumi << "pb-1";
  //sStream << lumi << "pb-1";
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
  canvasName.push_back("closureFR");
  canvasName.push_back("closureFRNoHeep");
  canvasName.push_back("closureFRCorrDY");
  canvasName.push_back("closureFRCorrFull");

  vector<TString> canvasTitle;
  canvasTitle.push_back("Fake rate closure test uncorrected");
  canvasTitle.push_back("Fake rate closure test GSF (non HEEP)");
  canvasTitle.push_back("Fake rate closure test DY corrected");
  canvasTitle.push_back("Fake rate closure test fully corrected");

  vector<TString> legendHeepGsf;
  legendHeepGsf.push_back("HEEP-GSF uncorrected");
  legendHeepGsf.push_back("HEEP-GSF(non HEEP)");
  legendHeepGsf.push_back("HEEP-GSF DY corrected");
  legendHeepGsf.push_back("HEEP-GSF DY corrected");

  vector<TString> legendGsfGsf;
  legendGsfGsf.push_back("GSF-GSF uncorrected");
  legendGsfGsf.push_back("GSF-GSF (non HEEP)");
  legendGsfGsf.push_back("GSF-GSF (non HEEP)");
  legendGsfGsf.push_back("GSF-GSF W+jet + \\gamma+jet corrected");

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
      c0->SetLogy();
      gStyle->SetOptStat(0);
      gStyle->SetPadTickY(1);
  
      // get the histograms 
      input.cd(folderGsfGsfHisto[j]);
      TH1F *histoGsfGsf = (TH1F *)gDirectory->Get(gsfGsfHisto[j] + acroSuffix[p]);
      input.cd(folderHeepGsfHisto[j]);
      TH1F *histoHeepGsf = (TH1F *)gDirectory->Get(heepGsfHisto[j] + acroSuffix[p]);
      input.cd(folderDataHisto);
      TH1F *histoData = (TH1F *)gDirectory->Get("histoHeepHeepMass" + acroSuffix[p]);

      TH1F *histoGsfGsfRebinned = (TH1F *)histoGsfGsf->Rebin(nBins, gsfGsfHisto[j] + "Rebinned" + acroSuffix[p], binArray);
      TH1F *histoHeepGsfRebinned = (TH1F *)histoHeepGsf->Rebin(nBins, heepGsfHisto[j] + "Rebinned" + acroSuffix[p], binArray);
      TH1F *histoDataRebinned = (TH1F *)histoData->Rebin(nBins, "histoHeepHeepMassRebinned" + acroSuffix[p], binArray);
  
      histoGsfGsfRebinned->SetLineColor(4);
      histoGsfGsfRebinned->SetMarkerColor(4);
      histoGsfGsfRebinned->SetMarkerStyle(20);
      histoHeepGsfRebinned->SetLineColor(2);
      histoHeepGsfRebinned->SetMarkerColor(2);
      histoHeepGsfRebinned->SetMarkerStyle(21);
      histoDataRebinned->SetLineColor(1);
      histoDataRebinned->SetMarkerColor(1);
      histoDataRebinned->SetMarkerStyle(20);
 
      sStream.str("");
      if (binning.size() > 1)
        sStream << "# of events / bin";
      else
        sStream << "# of events / " << binning.begin()->second << "GeV/c^{2}";

      if (plotData) { 
        histoDataRebinned->SetMinimum(yAxisMin);
        histoDataRebinned->SetTitle(canvasTitle[j] + suffix[p]);
        histoDataRebinned->GetYaxis()->SetTitle(sStream.str().c_str());

        histoDataRebinned->Draw();
        histoHeepGsfRebinned->Draw("sames");
      } else {
        histoHeepGsfRebinned->SetTitle(canvasTitle[j] + suffix[p]);
        histoHeepGsfRebinned->GetYaxis()->SetTitle(sStream.str().c_str());
        histoHeepGsfRebinned->Draw();
      }
      histoGsfGsfRebinned->Draw("sames");
  
      sStream.str("");
      sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
      TPaveLabel *label0 = new TPaveLabel(0.6, 0.79, 0.9, 0.89, sStream.str().c_str(), "brNDC");
      label0->SetFillColor(0);
      label0->SetFillStyle(0);
      label0->SetBorderSize(0);
      label0->SetTextSize(0.30);
      label0->Draw("sames");
      TPaveLabel *label1 = new TPaveLabel(0.7, 0.89, 0.91, 0.98, "CMS preliminary", "brNDC");
      label1->SetFillColor(0);
      label1->SetFillStyle(0);
      label1->SetBorderSize(0);
      label1->SetTextSize(0.40);
      label1->Draw("sames");
  
      TLegend *legend = new TLegend(0.38, 0.6, 0.53, 0.9);
      if (!plotData) legend->SetY2(0.8);
      legend->SetTextSize(0.03);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      if (plotData) legend->AddEntry(histoDataRebinned, "Data", "lep");
      legend->AddEntry(histoHeepGsfRebinned, legendHeepGsf[j], "lep");
      legend->AddEntry(histoGsfGsfRebinned, legendGsfGsf[j], "lep");
      legend->Draw("sames"); 
    } // end loop over eta ranges
  } // end loop over corrections
  input.Close();
}
