#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include "TColor.h"
#include "TH1F.h"
#include "THStack.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "THashList.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Math/QuantFuncMathCore.h"

#include "makeHistoFromBranch.C"

#define DATA 0
#define TTBAR 1
#define ZTT 2
#define WW 3
#define WZ 4
#define ZZ 5
#define TW 6
#define ZMM 7
#define ZEE 8
#define WJET 9

#define ALL 0
#define ALLCUM 1
#define SS 1
#define SSCUM 3
#define OS 2
#define OSCUM 5

float CalcBgSum(vector<vector<TH1F *> > &histos, vector<bool> &samples, int region, int lowerBin, int upperBin = -1);
float CalcSystErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin = -1);
float CalcAllErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples,  int region, int lowerBin, int upperBin = -1);
float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin = -1, bool calcQcdErr = 0);
float CalcSSQcdErr(vector<vector<TH1F *> > &histos, vector<float> &errors, int lowerBin, int upperBin = -1);
float CalcAllErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int sample, int region, int lowerBin, int upperBin = -1, bool calcQcdErr = 0);

void macro_MakeQcdClosureTest_SSvsOS()
{
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("./emuSpec_19703pb-1.root", "open");
  TFile input("./emuSpec_singleMuTrg_altdiboson_19706pb-1.root", "open");
  input.cd();

  TParameter<float> *lumi = (TParameter<float> *)input.Get("lumi");

  const int nBins = 75;
  const int qcdEst = 1; // estimation method of QCD contribution. none(0), from SS spectrum(1), from fake rate(2)
  const bool topReweighting = 0;

  int eRegion = 2; // electron region EB(0), EE(1), EB+EE(2)

  bool plotType[2];
  plotType[0] = 1;  // emu spectrum
  plotType[1] = 1;  // cumulative emu spectrum

  const bool plotRatio = 1; // plot frMethod/ssMethod
  const bool plotRatioBelowSpec = 1; // plot frMethod/ssMethod below spectrum
  const bool logPlotX = 0;
  const bool logPlotY = 1;
  const bool prelim = 1;
  const bool overflowBin = 1;

  float xRangeMin = 60.;
  float xRangeMax = 1200.;
  //float xRangeMin = 0.;
  //float xRangeMax = 1500.;
  float yRangeMin[6] = {0.002, 0.002, 0.002, 0.4, 0.4, 0.4};
  float yRangeMax[6] = {30, 10, 30, 3000, 2000, 3000};
  float yRangeMinRatio[3] = {0.3, 0.3, 0.3};
  float yRangeMaxRatio[3] = {1.7, 1.7, 1.7};
  float fitMin = xRangeMin;
  float fitMax = xRangeMax; // set to highest bin with a data point
  float xRangeMinRatio = fitMin;
  float xRangeMaxRatio = fitMax;

  // output file formats
  const bool saveRatio = 0;
  const bool saveSpec = 0;
  const bool saveCumSpec = 0;
  const bool saveAsPdf = 1;
  const bool saveAsPng = 1;
  const bool saveAsRoot = 0;
  const char *fileNameExtra = "";
  //const char *fileNameExtra = "madgraphTTbar_";
  const char *plotDir = "./plottemp/";

  int font = 42; //62
  ////////////////////////////////////////////////////////////////////////////

  // systematic errors
  float systErrLumi = ((TParameter<float> *)input.Get("systErrLumi"))->GetVal();
  systErrLumi = 0.; // since we normalize to the Z peak
  float systErrEff = ((TParameter<float> *)input.Get("systErrEff"))->GetVal(); // muon err & ele err
  THashList *systErrMCs = (THashList *)input.Get("systErrMCs");
  vector<float> systErrMC;
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcTtbar"))->GetVal());  // NNLO ttbar
  //systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcTtbar700to1000"))->GetVal());  // NLO ttbar700to1000
  //systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcTtbar1000up"))->GetVal());  // NLO ttbar1000up
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcDyTauTau"))->GetVal()); //z->tt
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcWW"))->GetVal()); //WW
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcWZ"))->GetVal()); //WZ
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcZZ"))->GetVal()); //ZZ
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcTW"))->GetVal()); //tW
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcDyMuMu"))->GetVal()); //Z->mm
  systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcDyEE"))->GetVal()); //Z->ee
  if (qcdEst == 2) systErrMC.push_back(0.4); // qcd error
  else systErrMC.push_back(((TParameter<float> *)systErrMCs->FindObject("systErrMcWJets"))->GetVal());  //WJets

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString xAxisTitle[1] = {"m(e#mu)"};
  TString nameSuffix[2] = {"", "cumul"};
  TString titleSuffix[2] = {"", " - Cumulative"};

  vector<TH1F *> emuMass_qcdFromFakeSS;
  vector<TH1F *> emuMass_qcdFromFakeOS;
  vector<TH1F *> dataOverBgHist;

  // define the binning
  vector<float> binning;
  if (logPlotX) {
    //for (float bin = 0.; bin < 100.; bin += 5.)
    //  binning.push_back(bin);
    for (float bin = 0.; bin < 200.; bin += 10.)
      binning.push_back(bin);
    for (float bin = 200.; bin < 400.; bin += 20.)
      binning.push_back(bin);
    for (float bin = 400.; bin < 500.; bin += 25.)
      binning.push_back(bin);
    for (float bin = 500.; bin <= 620.; bin += 40.)
      binning.push_back(bin);
      binning.push_back(670.);
      binning.push_back(720.);
      binning.push_back(780.);
      binning.push_back(840.);
      binning.push_back(920.);
      binning.push_back(1000.);
      binning.push_back(1100.);
      binning.push_back(1220.);
      binning.push_back(1380.);
      binning.push_back(1500.);
  } else {
    //for (float bin = 0.; bin <= 1500.; bin += 20.)
    //  binning.push_back(bin);
    for (float bin = 0.; bin < 200.; bin += 20.)
      binning.push_back(bin);
    for (float bin = 200.; bin < 400.; bin += 40.)
      binning.push_back(bin);
    for (float bin = 400.; bin < 700.; bin += 50.)
      binning.push_back(bin);
    for (float bin = 700.; bin < 1000.; bin += 75.)
      binning.push_back(bin);
    for (float bin = 1000.; bin < 1200.; bin += 100.)
      binning.push_back(bin);
    for (float bin = 1200.; bin <= 1500.; bin += 150.)
      binning.push_back(bin);
  }

  THashList *mcWeights = (THashList *)input.Get("mcWeights");
  THashList *nGenEvents = (THashList *)input.Get("nGenEvents");

  // vectors for processes with multiple samples 
  vector<const char *> cutVarsEmpty;
  vector<float> lowCutsEmpty;
  vector<float> highCutsEmpty;
  vector<float> mcWeightsForCutRangesEmpty;

  // loop to get normal and cumulated spectrum
  for (unsigned int j = 0; j < 2; ++j) {
    input.cd();

    bool normToBin = true;
    if (j > 0) normToBin = false;

    // qcd contribution
    emuMass_qcdFromFakeSS.push_back((MakeHistoFromBranch(&input, "frEmuTree_data", "", "mass", 1, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x300)));
    emuMass_qcdFromFakeOS.push_back((MakeHistoFromBranch(&input, "frEmuTree_data", "", "mass", -1, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x300)));
    // normalize to bin width
    if (j < 1) {
      for (int i = 1; i < emuMass_qcdFromFakeSS.back()->GetNbinsX() + 1; ++i) {
        emuMass_qcdFromFakeSS.back()->SetBinContent(i, emuMass_qcdFromFakeSS.back()->GetBinContent(i) / emuMass_qcdFromFakeSS.back()->GetBinWidth(i));
        emuMass_qcdFromFakeSS.back()->SetBinError(i, emuMass_qcdFromFakeSS.back()->GetBinError(i) / emuMass_qcdFromFakeSS.back()->GetBinWidth(i));
        emuMass_qcdFromFakeOS.back()->SetBinContent(i, emuMass_qcdFromFakeOS.back()->GetBinContent(i) / emuMass_qcdFromFakeOS.back()->GetBinWidth(i));
        emuMass_qcdFromFakeOS.back()->SetBinError(i, emuMass_qcdFromFakeOS.back()->GetBinError(i) / emuMass_qcdFromFakeOS.back()->GetBinWidth(i));
      }
    }

    dataOverBgHist.push_back((TH1F*)emuMass_qcdFromFakeSS.back()->Clone("dataOverBgHist"));

    // add overflow in last bin
    if (j == 0 && overflowBin) {
      dataOverBgHist.back()->SetBinContent(dataOverBgHist.back()->GetNbinsX(), dataOverBgHist.back()->GetBinContent(dataOverBgHist.back()->GetNbinsX()) + dataOverBgHist.back()->GetBinContent(dataOverBgHist.back()->GetNbinsX() + 1));
      emuMass_qcdFromFakeSS.back()->SetBinContent(emuMass_qcdFromFakeSS.back()->GetNbinsX(), emuMass_qcdFromFakeSS.back()->GetBinContent(emuMass_qcdFromFakeSS.back()->GetNbinsX()) + emuMass_qcdFromFakeSS.back()->GetBinContent(emuMass_qcdFromFakeSS.back()->GetNbinsX() + 1));
      emuMass_qcdFromFakeOS.back()->SetBinContent(emuMass_qcdFromFakeOS.back()->GetNbinsX(), emuMass_qcdFromFakeOS.back()->GetBinContent(emuMass_qcdFromFakeOS.back()->GetNbinsX()) + emuMass_qcdFromFakeOS.back()->GetBinContent(emuMass_qcdFromFakeOS.back()->GetNbinsX() + 1));
    }

    // integrate from the right side
    if (j == 1) { 
      // loop over bins
      double error;
      for (int i = 1; i < nBins + 1; ++i) {
        emuMass_qcdFromFakeSS.back()->SetBinContent(i, emuMass_qcdFromFakeSS.back()->IntegralAndError(i, nBins, error));
        emuMass_qcdFromFakeSS.back()->SetBinError(i, error);
        emuMass_qcdFromFakeOS.back()->SetBinContent(i, emuMass_qcdFromFakeOS.back()->IntegralAndError(i, nBins, error));
        emuMass_qcdFromFakeOS.back()->SetBinError(i, error);
      }     
    }

    if (!plotType[j]) continue;

    TCanvas *emuPlot;
    TPad *specPad;
    if (plotRatioBelowSpec && j == 0) {
      emuPlot = new TCanvas("emuPlot" + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 100, 100, 900, 720);
      specPad = new TPad("specPad" + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 0., 0.33, 1., 1.);
      specPad->SetBottomMargin(0.06);
    } else {
      emuPlot = new TCanvas("emuPlot" + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 100, 100, 900, 600);
      specPad = new TPad("specPad" + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 0., 0., 1., 1.);
      specPad->SetBottomMargin(0.12);
    }
    specPad->SetBorderMode(0);
    specPad->SetBorderSize(2);
    specPad->SetFrameBorderMode(0);
    specPad->SetFillColor(0);
    specPad->SetFrameFillColor(0);
    if (logPlotX) specPad->SetLogx();
    if (logPlotY) specPad->SetLogy();
    specPad->SetLeftMargin(0.11);
    specPad->SetRightMargin(0.09);
    specPad->SetTopMargin(0.08);
    specPad->SetTickx(1);
    specPad->SetTicky(1);
    specPad->Draw();
    specPad->cd();
 
    gStyle->SetTitleFont(font);
    gStyle->SetLabelFont(font);
    gStyle->SetLegendFont(font);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTitleXOffset(1.);
    gStyle->SetTitleYOffset(1.3);
    gPad->SetTicks(1, 1);

    // plot spectrum
    emuMass_qcdFromFakeOS.back()->SetLineColor(kOrange);
    emuMass_qcdFromFakeOS.back()->SetFillColor(TColor::GetColor("#ffff66"));
    emuMass_qcdFromFakeOS.back()->SetLineWidth(2);
    emuMass_qcdFromFakeOS.back()->Draw("hist");
    emuMass_qcdFromFakeSS.back()->SetLineColor(kRed);
    emuMass_qcdFromFakeSS.back()->SetLineWidth(2);
    emuMass_qcdFromFakeSS.back()->Draw("esame");

    if (plotRatioBelowSpec && j == 0) {
      emuMass_qcdFromFakeOS.back()->GetXaxis()->SetTitle("");
    } else {
      emuMass_qcdFromFakeOS.back()->GetXaxis()->SetTitle(xAxisTitle[0] + " [GeV]");
    }
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetTitleFont(font);
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetTitleSize(0.047);
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetTitleOffset(0.9);
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetLabelFont(font);
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetLabelSize(0.05);
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetMoreLogLabels();
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetNoExponent();
    emuMass_qcdFromFakeOS.back()->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax); 
    //emuMass_qcdFromFakeOS.back()->GetXaxis()->SetLimits(xRangeMin, xRangeMax); 
    if (j == 1) emuMass_qcdFromFakeOS.back()->GetYaxis()->SetTitle("Events #geq " + xAxisTitle[0]);
    else emuMass_qcdFromFakeOS.back()->GetYaxis()->SetTitle("Events / GeV");
    emuMass_qcdFromFakeOS.back()->GetYaxis()->SetTitleFont(font);
    emuMass_qcdFromFakeOS.back()->GetYaxis()->SetTitleSize(0.047);
    emuMass_qcdFromFakeOS.back()->GetYaxis()->SetTitleOffset(1.1);
    emuMass_qcdFromFakeOS.back()->GetYaxis()->SetLabelFont(font);
    emuMass_qcdFromFakeOS.back()->GetYaxis()->SetLabelSize(0.05);
    emuMass_qcdFromFakeOS.back()->SetMinimum(yRangeMin[j * 3]); 
    emuMass_qcdFromFakeOS.back()->SetMaximum(yRangeMax[j * 3]); 

    // redraw axis
    emuMass_qcdFromFakeOS.back()->Draw("sameaxis");

    // legend and labels
    TLegend legend(0.7, 0.646, 0.901, 0.885);
    legend.SetTextFont(font);
    legend.SetTextSize(0.03);
    legend.SetBorderSize(0);
    legend.SetLineColor(1);
    legend.SetLineStyle(1);
    legend.SetLineWidth(1);
    legend.SetFillColor(19);
    legend.SetFillStyle(0);
    legend.AddEntry(emuMass_qcdFromFakeSS.back(), "jets (Fake Rate) SS" ,"le");
    legend.AddEntry(emuMass_qcdFromFakeOS.back(), "jets (Fake Rate) OS" ,"f");
    legend.DrawClone("sames");
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(font);
    tex->SetLineWidth(2);
    tex->SetTextSize(0.042);
    if (prelim) tex->DrawLatex(0.325, 0.853, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
    else tex->DrawLatex(0.405, 0.853, "CMS, 8 TeV, 19.7 fb^{-1}");
    if (eRegion == 0) tex->DrawLatex(0.325, 0.775, "e in barrel");
    if (eRegion == 1) tex->DrawLatex(0.325, 0.775, "e in endcap");

    // denominator histogram for frMethod/ssMethod plot
    TH1F *denHist = (TH1F *)emuMass_qcdFromFakeOS.back()->Clone("denHist");

    stringstream sStream;
    if (!plotRatioBelowSpec || j > 0) {
      sStream << plotDir << "qcdClosureTestSpec";
      if (eRegion == 0) sStream << "EB_";
      if (eRegion == 1) sStream << "EE_";
      sStream << fileNameExtra << nameSuffix[j];
      if (j > 0) sStream << "_";
      if (!logPlotY) sStream << "lin_";
      sStream << lumi->GetVal() << "pb-1";
      TString saveFileName = sStream.str();
      if ((j == 0 && saveSpec) || (j > 0 && saveCumSpec)) {
        if (saveAsPdf) emuPlot->Print(saveFileName + ".pdf", "pdf");
        if (saveAsPng) emuPlot->Print(saveFileName + ".png", "png");
        if (saveAsRoot) emuPlot->Print(saveFileName + ".root", "root");
      }
    }

    if (j > 0 || !plotRatio) continue;
    // plot a frMethod/ssMethod histogram
    dataOverBgHist.back()->Divide(denHist);

    float fontScaleBot = 1.;
    TPad *pullPad = new TPad("pullPad" + nameSuffix[j], "Ratio" + titleSuffix[j], 0., 0., 1., 1.);
    TCanvas *dataOverBgPlot;
    if (plotRatioBelowSpec) {
      emuPlot->cd();
      pullPad->SetPad(0., 0., 1., 0.33);
      pullPad->SetBottomMargin(0.22);
      pullPad->SetTopMargin(0.);
      fontScaleBot = specPad->GetHNDC() / pullPad->GetHNDC();
    } else {
      dataOverBgPlot = new TCanvas("dataOverBgPlot" + nameSuffix[j], "Ratio" + titleSuffix[j], 100, 100, 900, 600);
      pullPad->SetBottomMargin(0.12);
      pullPad->SetTopMargin(0.08);
    }
    pullPad->Draw();
    pullPad->cd();
    pullPad->SetBorderMode(0);
    pullPad->SetBorderSize(2);
    pullPad->SetFrameBorderMode(0);
    pullPad->SetFillColor(0);
    pullPad->SetFrameFillColor(0);
    if (logPlotX) pullPad->SetLogx();
    pullPad->SetLeftMargin(0.11);
    pullPad->SetRightMargin(0.09);
    pullPad->SetTickx(1);
    pullPad->SetTicky(1);
    pullPad->SetGridy(1);

    dataOverBgHist.back()->SetLineWidth(2);
    dataOverBgHist.back()->SetLineColor(kRed);
    //dataOverBgHist.back()->SetMarkerStyle(20);
    //dataOverBgHist.back()->SetMarkerSize(1.1);

    dataOverBgHist.back()->GetXaxis()->SetTitle(xAxisTitle[0] + " [GeV]");
    dataOverBgHist.back()->GetXaxis()->SetTitleFont(font);
    dataOverBgHist.back()->GetXaxis()->SetTitleSize(0.047 * fontScaleBot);
    dataOverBgHist.back()->GetXaxis()->SetTitleOffset(0.9);
    dataOverBgHist.back()->GetXaxis()->SetLabelFont(font);
    dataOverBgHist.back()->GetXaxis()->SetLabelSize(0.05 * fontScaleBot);
    dataOverBgHist.back()->GetXaxis()->SetMoreLogLabels();
    dataOverBgHist.back()->GetXaxis()->SetNoExponent();
    dataOverBgHist.back()->GetXaxis()->SetRangeUser(xRangeMinRatio, xRangeMaxRatio - 0.1);
    dataOverBgHist.back()->GetYaxis()->SetTitle("Ratio ");
    dataOverBgHist.back()->GetYaxis()->SetTitleFont(font);
    dataOverBgHist.back()->GetYaxis()->SetTitleSize(0.047 * fontScaleBot);
    dataOverBgHist.back()->GetYaxis()->SetTitleOffset(1.1 / fontScaleBot);
    dataOverBgHist.back()->GetYaxis()->SetLabelFont(font);
    dataOverBgHist.back()->GetYaxis()->SetLabelSize(0.05 * fontScaleBot);
    dataOverBgHist.back()->GetYaxis()->SetRangeUser(yRangeMinRatio[0], yRangeMaxRatio[0]);

    dataOverBgHist.back()->Draw();

    //TF1 *f0 = new TF1("f0", "[0]");
    //dataOverBgHist.back()->Fit("f0", "", "", fitMin, fitMax);

    if (!plotRatioBelowSpec) {
      tex->SetTextSize(0.042 * fontScaleBot);
      if (prelim) tex->DrawLatex(0.150, 0.853, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
      else tex->DrawLatex(0.150, 0.853, "CMS, 8 TeV, 19.7 fb^{-1}");
      tex->SetTextSize(0.042 * fontScaleBot);
      //tex->DrawLatex(0.150, 0.764, Form("#chi^{2} / ndf: %.2f / %i", f0->GetChisquare(), f0->GetNDF()));
      //tex->DrawLatex(0.150, 0.720, Form("p0: %.4f #pm %0.4f", f0->GetParameter(0), f0->GetParError(0)));
      if (eRegion == 0) tex->DrawLatex(0.708, 0.853, "e in barrel");
      if (eRegion == 1) tex->DrawLatex(0.708, 0.853, "e in endcap");
    } else {
      tex->SetTextSize(0.042 * fontScaleBot);
      //tex->DrawLatex(0.150, 0.875, Form("#chi^{2} / ndf: %.2f / %i", f0->GetChisquare(), f0->GetNDF()));
      //tex->DrawLatex(0.150, 0.775, Form("p0: %.4f #pm %0.4f", f0->GetParameter(0), f0->GetParError(0)));
    }

    // safe in various file formats
    if (saveRatio && !plotRatioBelowSpec) {
      sStream.str("");
      sStream << plotDir << "qcdClosureTestRatio";
      if (eRegion == 0) sStream << "EB_";
      if (eRegion == 1) sStream << "EE_";
      sStream << fileNameExtra << lumi->GetVal() << "pb-1";
      TString saveFileName = sStream.str();
      if (saveAsPdf) dataOverBgPlot->Print(saveFileName + ".pdf", "pdf");
      if (saveAsPng) dataOverBgPlot->Print(saveFileName + ".png", "png");
      if (saveAsRoot) dataOverBgPlot->Print(saveFileName + ".root", "root");
    }
    if ((saveSpec || saveRatio) && plotRatioBelowSpec) {
      sStream.str("");
      sStream << plotDir << "qcdClosureTestSpec";
      if (eRegion == 0) sStream << "EB_";
      if (eRegion == 1) sStream << "EE_";
      sStream << fileNameExtra << nameSuffix[j];
      if (j > 0) sStream << "_";
      if (!logPlotY) sStream << "lin_";
      sStream << lumi->GetVal() << "pb-1";
      TString saveFileName = sStream.str();
      if ((j == 0 && saveSpec) || (j > 0 && saveCumSpec)) {
        if (saveAsPdf) emuPlot->Print(saveFileName + ".pdf", "pdf");
        if (saveAsPng) emuPlot->Print(saveFileName + ".png", "png");
        if (saveAsRoot) emuPlot->Print(saveFileName + ".root", "root");
      }
    }
  } // end loop over normal or cumulated
}

// calculate the sum of several cumulated background histograms
float CalcBgSum(vector<vector<TH1F *> > &histos, vector<bool> &samples, int region, int lowerBin, int upperBin)
{
   float sumLow = 0.;
   float sumUp = 0.;

   for (unsigned int i = 0; i < samples.size(); ++i ) {
      if (samples[i]) {
         sumLow += histos.at(i+1).at(region)->GetBinContent(lowerBin);
         if (upperBin > -1) sumUp += histos.at(i+1).at(region)->GetBinContent(upperBin);
      }
   }
   return sumLow - sumUp;
}

float CalcSystErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin)
{
   float err2 = 0.;

   for (unsigned int i = 0; i < samples.size(); ++i) {
      if (samples[i]) {
         float contLow = histos.at(i+1).at(region)->GetBinContent(lowerBin);
         float contUp = 0.;
         if (upperBin > -1) contUp = histos.at(i+1).at(region)->GetBinContent(upperBin);
         err2 += (contLow - contUp) * (contLow - contUp) * errors[i] * errors[i];
      }
   }
   return sqrt(err2);
}

float CalcAllErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin)
{
   float numEv;

   for (unsigned int i = 0; i < samples.size(); ++i) {
      if (samples[i])
         numEv += histos.at(i+1).at(region)->Integral(lowerBin, upperBin);
   }

   float systErr = CalcSystErr(histos, errors, samples, region, lowerBin, upperBin);

   //return sqrt(numEv + systErr * systErr);
   return sqrt(systErr * systErr);
}

float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool calcQcdErr)
{
   if (!calcQcdErr) return CalcSystErr(histos, errors, samples, region, lowerBin, upperBin);

   errors.back() = CalcSSQcdErr(histos, errors, lowerBin, upperBin);

   return CalcSystErr(histos, errors, samples, region, lowerBin, upperBin);
}

float CalcSSQcdErr(vector<vector<TH1F *> > &histos, vector<float> &errors, int lowerBin, int upperBin)
{
   float qcdErr = 0.;
   vector<bool> allButQCD(histos.size() - 2, true);
   allButQCD.push_back(false);

   float allErrButQCD = CalcSystErr(histos, errors, allButQCD, SSCUM, lowerBin, upperBin);

   float qcdCont = histos.back().at(SSCUM)->GetBinContent(lowerBin) - histos.back().at(SSCUM)->GetBinContent(upperBin);
   qcdErr = allErrButQCD;
   return qcdErr / qcdCont;
}

float CalcAllErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int sample, int region, int lowerBin, int upperBin, bool calcQcdErr)
{
   float statErr = sqrt(histos.at(sample).at(region)->Integral(lowerBin, upperBin));
   float systErr = CalcSystErrWithQCD(histos, errors, samples, region, lowerBin, upperBin, calcQcdErr);

   //return sqrt(statErr * statErr + systErr * systErr);
   return sqrt(systErr * systErr);
}

