#include <stdio>
#include <string>
#include <sstream>
#include <vector>
#include <sstream>
#include "TH1F.h"
#include "TCanvas.h"

#define DATA 0
#define TTBAR 1
#define ZTT 2
#define WW 3
#define WZ 4
#define TW 5
#define WJET 6
#define ZMM 7
#define ZEE 8
#define ZZ 9
//#define QCD 10
#define QCD 9

#define ALL 0
#define LS 1
#define OS 2

#define EBEE 0
#define EB 1
#define EE 2

void macro_MakeDataOverBkgPlot()
{
  // parameters //////////////////////////////////////////////////////////////
  TFile input("testEmuSpec4684pb-1.root", "open");
  //TFile input("testEmuSpecSummerFallMix4684pb-1.root", "open");
  //TFile input("testEmuSpecPureSummer114684pb-1.root", "open");

  const float lumi = 4684.;

  bool plotSign[3];
  plotSign[0] = true;  // all
  plotSign[1] = true;  // LS like sign
  plotSign[2] = true;  // OS opposite sign

  // systematical errors
  vector<float> systErrMC;
  systErrMC.push_back(0.15);  //ttbar
  systErrMC.push_back(0.);    //z->tt  //FIXME
  systErrMC.push_back(0.);    //WW     //FIXME
  systErrMC.push_back(0.);    //WZ     //FIXME
  systErrMC.push_back(0.15);  //tW     //FIXME
  systErrMC.push_back(0.);    //WJets  //FIXME
  systErrMC.push_back(0.);    //Z->mm  //FIXME
  systErrMC.push_back(0.);    //Z->ee  //FIXME
  systErrMC.push_back(0.);    //ZZ     //FIXME
  float systErrLumi = 0.036;
  float systErrEff = 0.02;
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString histoSign[3] = {"", "LS_", "OS_"};
  TString histoTitleSign[3] = {"", " e^{#pm}#mu^{#pm}", " e^{+}#mu^{-}/e^{-}#mu^{+}"};
  TString xAxisTitle[3] = {"m(e#mu)", " m(e^{#pm}#mu^{#pm})", " m(e^{+}#mu^{-}/m(e^{-}#mu^{+})"};

  vector<TH1F *> emuMass_data;
  vector<TH1F *> emuMass_bg;

  int nBins = 42;
  Double_t binArray[43] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 560, 620, 680, 740, 800, 1000, 1500};

  // loop over full spectrum, LS and OS
  for (unsigned int k = 0; k < 3; ++k) {
    input.cd();

    // get the histograms
    emuMass_data.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "data"));
    emuMass_bg.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "ttbar"));

    if (!plotSign[k]) continue;

    TH1F *dataHistoRebinned = (TH1F *)emuMass_data.back()->Rebin(nBins, "dataHistoRebinned", binArray);
    TH1F *bgHistoRebinned = (TH1F *)emuMass_bg.back()->Rebin(nBins, "bgHistoRebinned", binArray);

    // set systematic error for ttbar MC
    for (int n = 0; n < emuMass_bg.back()->GetNbinsX() + 2; ++n) emuMass_bg.back()->SetBinError(n, emuMass_bg.back()->GetBinContent(n) * sqrt(systErrLumi*systErrLumi + systErrMC[0]*systErrMC[0]));

    // calculate (data-bkg)/bkg
    dataHistoRebinned->Add(bgHistoRebinned, -1);
    dataHistoRebinned->Divide(bgHistoRebinned);

    TCanvas *emuPlot = new TCanvas("dataOverBkg" + histoSign[k], "(data-bkg)/bkg" + histoTitleSign[k], 100, 100, 800, 600);
    emuPlot->SetBorderMode(0);
    emuPlot->SetFrameBorderMode(0);
    emuPlot->SetFillColor(0);
    emuPlot->SetFrameFillColor(0);
    emuPlot->cd();
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleXOffset(1.);
    gStyle->SetTitleYOffset(1.3);

    dataHistoRebinned->GetXaxis()->SetTitle(xAxisTitle[k] + " [GeV]");
    dataHistoRebinned->GetXaxis()->SetTitleSize(0.04);
    dataHistoRebinned->GetXaxis()->SetLabelSize(0.04);
    dataHistoRebinned->GetXaxis()->SetRangeUser(10., 800.);
    //dataHistoRebinned->GetXaxis()->SetRangeUser(0., 400.);
    //dataHistoRebinned->GetXaxis()->SetRangeUser(0., 900.);
    
    dataHistoRebinned->GetYaxis()->SetTitle("(data-bkg)/bkg");
    dataHistoRebinned->GetYaxis()->SetTitleSize(0.04);
    dataHistoRebinned->GetYaxis()->SetTitleOffset(1.2);
    
    // plot spectrum
    dataHistoRebinned->SetLineColor(kBlue);
    dataHistoRebinned->SetMarkerColor(kBlue);
    dataHistoRebinned->SetMarkerStyle(20);
    dataHistoRebinned->Draw("e");

    TF1 *f0 = new TF1("f0" + histoSign[k], "[0]");
    dataHistoRebinned->Fit("f0" + histoSign[k], "", "", 10., 800.);

    // redraw axis
    dataHistoRebinned->Draw("sameaxis");

    stringstream sStream;
    sStream.str("");
    sStream << "CMS preliminary    #sqrt{s} = 7 TeV    #int L dt = 4.7 fb^{-1}";
    //sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
    TPaveLabel labelLumi(0.13, 0.79, 0.69, 0.88, sStream.str().c_str(), "brNDC");
    labelLumi.SetFillColor(0);
    labelLumi.SetFillStyle(0);
    labelLumi.SetBorderSize(0);
    labelLumi.SetTextSize(0.40);
    labelLumi.DrawClone("sames");
  } // end loop over full, LS and OS
}

