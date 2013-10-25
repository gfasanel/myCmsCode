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
TGraphAsymmErrors * makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBars);
TH1F * MakeHistoFromBranch(TFile *input, const char *treeName, const char *brName, int signs, int &region, const char *cutVariable, float cutLow, float cutHigh, vector<float> &binning, unsigned int flags, bool normToBinWidth = false, float userScale = 1.);

void macro_MakeEMuInvMassPlot()
{
  // parameters //////////////////////////////////////////////////////////////
  TFile input("./emuSpec_19703pb-1.root", "open");
  //TFile input("test_19619pb-1.root", "open");
  input.cd();

  TParameter<float> *lumi = (TParameter<float> *)input.Get("lumi");

  const int nBins = 75;
  const int qcdEst = 1; // estimation method of QCD contribution. none(0), from SS spectrum(1), from fake rate(2)

  int eRegion = 2; // electron region EB(0), EE(1), EB+EE(2)

  bool plotSign[3];
  plotSign[0] = 1;  // all
  plotSign[1] = 1;  // SS same sign
  plotSign[2] = 1;  // OS opposite sign

  bool plotType[2];
  plotType[0] = 1;  // emu spectrum
  plotType[1] = 1;  // cumulative emu spectrum

  const bool plotPull = 1; // plot (data-bkg)/bkg
  const bool plotPullBelowSpec = 1; // plot (data-bkg)/bkg below spectrum
  const bool logPlotX = 0;
  const bool logPlotY = 1;
  const bool prelim = 1;
  const bool groupedPlot = 1;
  const bool overflowBin = 0;

  float xRangeMin = 60.;
  float xRangeMax = 1200.;
  //float xRangeMin = 0.;
  //float xRangeMax = 1500.;
  float yRangeMin[6] = {0.002, 0.002, 0.002, 0.2, 0.2, 0.2};
  float yRangeMax[6] = {250, 40, 250, 30000, 4000, 30000};
  float yRangeMinRatio[3] = {-0.7, -0.7, -0.7};
  float yRangeMaxRatio[3] = {0.7, 0.7, 0.7};
  float fitMin = xRangeMin;
  float fitMax = 1200.; // set to highest bin with a data point
  float xRangeMinRatio = xRangeMin;
  float xRangeMaxRatio = xRangeMax;

  // output file formats
  const bool savePull = 0;
  const bool saveSpec = 0;
  const bool saveCumSpec = 0;
  const bool saveAsPdf = 1;
  const bool saveAsPng = 1;
  const bool saveAsRoot = 1;
  const char *fileNameExtra = "";
  //const char *fileNameExtra = "madgraphTTbar_";
  const char *plotDir = "./plottemp/";

  // plot style
  int ttbarColour = TColor::GetColor("#ff6666");
  int zttColour = TColor::GetColor("#ff4d4d");
  int wwColour = TColor::GetColor("#ff3333");
  int wzColour = TColor::GetColor("#ff0f0f");
  int zzColour = TColor::GetColor("#eb0000");
  int twColour = TColor::GetColor("#db0000");
  int zmmColour=  TColor::GetColor("#66b3ff");
  int zeeColour=  TColor::GetColor("#99ccff");
  int wjetColour=  TColor::GetColor("#ffd324");
  int jetBkgColour = TColor::GetColor("#ffff66"); 

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

  TString histoSign[3] = {"", "SS_", "OS_"};
  TString xAxisTitle[3] = {"m(e#mu)", "m(e^{#pm}#mu^{#pm})", "m(e^{#pm}#mu^{#mp})"};
  TString nameSuffix[2] = {"", "cumul"};
  TString titleSuffix[2] = {"", " - Cumulative"};

  vector<TH1F *> emuMass_data;
  vector<TH1F *> emuMass_ttbar;
  vector<TH1F *> emuMass_ztautau;
  vector<TH1F *> emuMass_ww;
  vector<TH1F *> emuMass_wz;
  vector<TH1F *> emuMass_zz;
  vector<TH1F *> emuMass_tw;
  vector<TH1F *> emuMass_zmumu;
  vector<TH1F *> emuMass_zee;
  vector<TH1F *> emuMass_wjets;
  vector<TH1F *> emuMass_qcd;
  vector<TH1F *> emuMass_ttLike;
  vector<TH1F *> emuMass_fakePairs;

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
  TParameter<float> *mcWeight = (TParameter<float> *)mcWeights->FindObject("ttbar");
  TParameter<float> *mcWeight700to1000 = (TParameter<float> *)mcWeights->FindObject("ttbar700to1000");
  TParameter<float> *mcWeight1000up = (TParameter<float> *)mcWeights->FindObject("ttbar1000up");

  float totMcWeight = 1.;
  // determine qcd contribution
  TH1F *qcdContrib;
  if (qcdEst == 1) {
    TH1F *ssData = MakeHistoFromBranch(&input, "emuTree_data", "mass", SS, eRegion, "", 0., 0., binning, 0x100);
    TH1F *ssBg = MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", SS, eRegion, "genMTtbar", 0., 700., binning, 0x1DF);
    totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight700to1000->GetVal());
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", SS, eRegion, "genMTtbar", 700., 1000., binning, 0x19F), totMcWeight);
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", "mass", SS, eRegion, "genMTtbar", 700., 1000., binning, 0x19F), totMcWeight);
    totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight1000up->GetVal());
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", SS, eRegion, "genMTtbar", 1000., 1000000000., binning, 0x19F), totMcWeight);
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", "mass", SS, eRegion, "genMTtbar", 1000., 1000000000., binning, 0x19F), totMcWeight);
    //TH1F *ssBg = MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF);
    //TH1F *ssBg = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF);
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ztautau", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ww", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wz", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zz", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_tw", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zmumu", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zee", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wjets", "mass", SS, eRegion, "", 0., 0., binning, 0x1DF));
    qcdContrib = (TH1F *)ssData->Clone("qcdContrib_SS");
    qcdContrib->Add(ssBg, -1);
    for (int i = 0; i < qcdContrib->GetNbinsX() + 2; ++i) {
      if (qcdContrib->GetBinContent(i) < 0) qcdContrib->SetBinContent(i, 0.);
    }
    cout << "expected SS QCD events: " << ssData->Integral() - ssBg->Integral() << endl;
    cout << "derived SS QCD events: " << qcdContrib->Integral() << endl;
    cout << "scale factor: " << (ssData->Integral() - ssBg->Integral()) / qcdContrib->Integral()<< endl;
    qcdContrib->Scale((ssData->Integral() - ssBg->Integral()) / qcdContrib->Integral());
  } 

  // loop over full spectrum, SS and OS
  for (int k = 0; k < 3; ++k) {
    // loop to get normal and cumulated spectrum
    for (unsigned int j = 0; j < 2; ++j) {
      input.cd();

      bool normToBin = true;
      if (j > 0) normToBin = false;

      if (k == 2) k = -1;
      // make the histograms
      TH1F *dataOverBgHist = (TH1F *)MakeHistoFromBranch(&input, "emuTree_data", "mass", k, eRegion, "", 0., 0., binning, 0x100, normToBin);
      emuMass_data.push_back(MakeHistoFromBranch(&input, "emuTree_data", "mass", k, eRegion, "", 0., 0., binning, 0x100));

      TH1F *ttbarComb = MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", k, eRegion, "genMTtbar", 0., 700., binning, 0x1DF, normToBin);
      totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight700to1000->GetVal());
      ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", k, eRegion, "genMTtbar", 700., 1000., binning, 0x19F, normToBin), totMcWeight);
      ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", "mass", k, eRegion, "genMTtbar", 700., 1000., binning, 0x19F, normToBin), totMcWeight);
      totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight1000up->GetVal());
      ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", k, eRegion, "genMTtbar", 1000., 1000000000., binning, 0x19F, normToBin), totMcWeight);
      ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", "mass", k, eRegion, "genMTtbar", 1000., 1000000000., binning, 0x19F, normToBin), totMcWeight);
      emuMass_ttbar.push_back(ttbarComb);
      //emuMass_ttbar.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      //emuMass_ttbar.push_back(MakeHistoFromBranch(&input, "emuTree_ttbarto2l", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_ztautau.push_back(MakeHistoFromBranch(&input, "emuTree_ztautau", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_ww.push_back(MakeHistoFromBranch(&input, "emuTree_ww", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_wz.push_back(MakeHistoFromBranch(&input, "emuTree_wz", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_zz.push_back(MakeHistoFromBranch(&input, "emuTree_zz", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_tw.push_back(MakeHistoFromBranch(&input, "emuTree_tw", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_zmumu.push_back(MakeHistoFromBranch(&input, "emuTree_zmumu", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      emuMass_zee.push_back(MakeHistoFromBranch(&input, "emuTree_zee", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      if (qcdEst != 2) emuMass_wjets.push_back(MakeHistoFromBranch(&input, "emuTree_wjets", "mass", k, eRegion, "", 0., 0., binning, 0x1DF, normToBin));
      if (k == -1) k = 2;
      emuMass_data.back()->SetName("emuMass_" + histoSign[k] + "data" + nameSuffix[j]);
      emuMass_ttbar.back()->SetName("emuMass_" + histoSign[k] + "ttbar" + nameSuffix[j]);
      emuMass_ztautau.back()->SetName("emuMass_" + histoSign[k] + "ztautau" + nameSuffix[j]);
      emuMass_ww.back()->SetName("emuMass_" + histoSign[k] + "ww" + nameSuffix[j]);
      emuMass_wz.back()->SetName("emuMass_" + histoSign[k] + "wz" + nameSuffix[j]);
      emuMass_zz.back()->SetName("emuMass_" + histoSign[k] + "zz" + nameSuffix[j]);
      emuMass_tw.back()->SetName("emuMass_" + histoSign[k] + "tw" + nameSuffix[j]);
      emuMass_zmumu.back()->SetName("emuMass_" + histoSign[k] + "zmumu" + nameSuffix[j]);
      emuMass_zee.back()->SetName("emuMass_" + histoSign[k] + "zee" + nameSuffix[j]);
      if (qcdEst != 2) emuMass_wjets.back()->SetName("emuMass_" + histoSign[k] + "wjets" + nameSuffix[j]);

      // qcd contribution
      if (qcdEst > 0) {
        if (k == 2) k = -1;
        if (qcdEst == 2) {
          //qcdContrib = MakeHistoFromBranch(&input, "frEmuTree_data", "mass", ALL, eRegion, "", 0., 0., binning, 0x300);
          //if (k != ALL) emuMass_qcd.back()->Scale(0.5);
          qcdContrib = MakeHistoFromBranch(&input, "frEmuTree_data", "mass", k, eRegion, "", 0., 0., binning, 0x300);
        }
        emuMass_qcd.push_back((TH1F *)qcdContrib->Clone("emuMass_" + histoSign[k] + "qcd"));
        if (qcdEst == 1 && k == ALL) emuMass_qcd.back()->Scale(2.);
        // normalize to bin width
        if (j < 1) {
          for (int i = 1; i < emuMass_qcd.back()->GetNbinsX() + 1; ++i) {
            emuMass_qcd.back()->SetBinContent(i, emuMass_qcd.back()->GetBinContent(i) / emuMass_qcd.back()->GetBinWidth(i));
            emuMass_qcd.back()->SetBinError(i, emuMass_qcd.back()->GetBinError(i) / emuMass_qcd.back()->GetBinWidth(i));
          }
        }
        if (k == -1) k = 2;
      }

      // add overflow in last bin
      if (j == 0 && overflowBin) {
        dataOverBgHist->SetBinContent(dataOverBgHist->GetNbinsX(), dataOverBgHist->GetBinContent(dataOverBgHist->GetNbinsX()) + dataOverBgHist->GetBinContent(dataOverBgHist->GetNbinsX() + 1));
        emuMass_data.back()->SetBinContent(emuMass_data.back()->GetNbinsX(), emuMass_data.back()->GetBinContent(emuMass_data.back()->GetNbinsX()) + emuMass_data.back()->GetBinContent(emuMass_data.back()->GetNbinsX() + 1));
        emuMass_ttbar.back()->SetBinContent(emuMass_ttbar.back()->GetNbinsX(), emuMass_ttbar.back()->GetBinContent(emuMass_ttbar.back()->GetNbinsX()) + emuMass_ttbar.back()->GetBinContent(emuMass_ttbar.back()->GetNbinsX() + 1));
        emuMass_ztautau.back()->SetBinContent(emuMass_ztautau.back()->GetNbinsX(), emuMass_ztautau.back()->GetBinContent(emuMass_ztautau.back()->GetNbinsX()) + emuMass_ztautau.back()->GetBinContent(emuMass_ztautau.back()->GetNbinsX() + 1));
        emuMass_ww.back()->SetBinContent(emuMass_ww.back()->GetNbinsX(), emuMass_ww.back()->GetBinContent(emuMass_ww.back()->GetNbinsX()) + emuMass_ww.back()->GetBinContent(emuMass_ww.back()->GetNbinsX() + 1));
        emuMass_wz.back()->SetBinContent(emuMass_wz.back()->GetNbinsX(), emuMass_wz.back()->GetBinContent(emuMass_wz.back()->GetNbinsX()) + emuMass_wz.back()->GetBinContent(emuMass_wz.back()->GetNbinsX() + 1));
        emuMass_zz.back()->SetBinContent(emuMass_zz.back()->GetNbinsX(), emuMass_zz.back()->GetBinContent(emuMass_zz.back()->GetNbinsX()) + emuMass_zz.back()->GetBinContent(emuMass_zz.back()->GetNbinsX() + 1));
        emuMass_tw.back()->SetBinContent(emuMass_tw.back()->GetNbinsX(), emuMass_tw.back()->GetBinContent(emuMass_tw.back()->GetNbinsX()) + emuMass_tw.back()->GetBinContent(emuMass_tw.back()->GetNbinsX() + 1));
        emuMass_zmumu.back()->SetBinContent(emuMass_zmumu.back()->GetNbinsX(), emuMass_zmumu.back()->GetBinContent(emuMass_zmumu.back()->GetNbinsX()) + emuMass_zmumu.back()->GetBinContent(emuMass_zmumu.back()->GetNbinsX() + 1));
        emuMass_zee.back()->SetBinContent(emuMass_zee.back()->GetNbinsX(), emuMass_zee.back()->GetBinContent(emuMass_zee.back()->GetNbinsX()) + emuMass_zee.back()->GetBinContent(emuMass_zee.back()->GetNbinsX() + 1));
        if (qcdEst != 2) emuMass_wjets.back()->SetBinContent(emuMass_wjets.back()->GetNbinsX(), emuMass_wjets.back()->GetBinContent(emuMass_wjets.back()->GetNbinsX()) + emuMass_wjets.back()->GetBinContent(emuMass_wjets.back()->GetNbinsX() + 1));
        if (qcdEst > 0) emuMass_qcd.back()->SetBinContent(emuMass_qcd.back()->GetNbinsX(), emuMass_qcd.back()->GetBinContent(emuMass_qcd.back()->GetNbinsX()) + emuMass_qcd.back()->GetBinContent(emuMass_qcd.back()->GetNbinsX() + 1));
      }

      // make grouped histograms
      emuMass_ttLike.push_back((TH1F *)emuMass_ttbar.back()->Clone("emuMass_" + histoSign[k] + "ttLike" + nameSuffix[j]));
      emuMass_ttLike.back()->Add(emuMass_ztautau.back());
      emuMass_ttLike.back()->Add(emuMass_ww.back());
      emuMass_ttLike.back()->Add(emuMass_wz.back());
      emuMass_ttLike.back()->Add(emuMass_zz.back());
      emuMass_ttLike.back()->Add(emuMass_tw.back());
      emuMass_fakePairs.push_back((TH1F *)emuMass_zmumu.back()->Clone("emuMass_" + histoSign[k] + "fakePairs" + nameSuffix[j]));
      emuMass_fakePairs.back()->Add(emuMass_zee.back());
      if (qcdEst != 2) emuMass_fakePairs.back()->Add(emuMass_wjets.back());
      if (qcdEst > 0) emuMass_fakePairs.back()->Add(emuMass_qcd.back());

      // integrate from the right side
      if (j == 1) { 
        // loop over bins
        double error;
        for (int i = 1; i < nBins + 1; ++i) {
          emuMass_data.back()->SetBinContent(i, emuMass_data.back()->IntegralAndError(i, nBins, error));
          emuMass_data.back()->SetBinError(i, error);
          emuMass_ttbar.back()->SetBinContent(i, emuMass_ttbar.back()->IntegralAndError(i, nBins, error));
          emuMass_ttbar.back()->SetBinError(i, error);
          emuMass_ztautau.back()->SetBinContent(i, emuMass_ztautau.back()->IntegralAndError(i, nBins, error));
          emuMass_ztautau.back()->SetBinError(i, error);
          emuMass_ww.back()->SetBinContent(i, emuMass_ww.back()->IntegralAndError(i, nBins, error));
          emuMass_ww.back()->SetBinError(i, error);
          emuMass_wz.back()->SetBinContent(i, emuMass_wz.back()->IntegralAndError(i, nBins, error));
          emuMass_wz.back()->SetBinError(i, error);
          emuMass_zz.back()->SetBinContent(i, emuMass_zz.back()->IntegralAndError(i, nBins, error));
          emuMass_zz.back()->SetBinError(i, error);
          emuMass_tw.back()->SetBinContent(i, emuMass_tw.back()->IntegralAndError(i, nBins, error));
          emuMass_tw.back()->SetBinError(i, error);
          emuMass_zmumu.back()->SetBinContent(i, emuMass_zmumu.back()->IntegralAndError(i, nBins, error));
          emuMass_zmumu.back()->SetBinError(i, error);
          emuMass_zee.back()->SetBinContent(i, emuMass_zee.back()->IntegralAndError(i, nBins, error));
          emuMass_zee.back()->SetBinError(i, error);
          if (qcdEst != 2) {
            emuMass_wjets.back()->SetBinContent(i, emuMass_wjets.back()->IntegralAndError(i, nBins, error));
            emuMass_wjets.back()->SetBinError(i, error);
          }
          if (qcdEst > 0) {
            emuMass_qcd.back()->SetBinContent(i, emuMass_qcd.back()->IntegralAndError(i, nBins, error));
            emuMass_qcd.back()->SetBinError(i, error);
          }

          // grouped histograms
          emuMass_ttLike.back()->SetBinContent(i, emuMass_ttLike.back()->IntegralAndError(i, nBins, error));
          emuMass_ttLike.back()->SetBinError(i, error);
          emuMass_fakePairs.back()->SetBinContent(i, emuMass_fakePairs.back()->IntegralAndError(i, nBins, error));
          emuMass_fakePairs.back()->SetBinError(i, error);
        }     
      }

      if (!plotSign[k]) continue;
      if (!plotType[j]) continue;

      TCanvas *emuPlot;
      TPad *specPad;
      if (plotPullBelowSpec && j == 0) {
        emuPlot = new TCanvas("emuPlot" + histoSign[k] + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 100, 100, 900, 900);
        specPad = new TPad("specPad" + histoSign[k] + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 0., 0.33, 1., 1.);
        specPad->SetBottomMargin(0.06);
      } else {
        emuPlot = new TCanvas("emuPlot" + histoSign[k] + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 100, 100, 900, 600);
        specPad = new TPad("specPad" + histoSign[k] + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 0., 0., 1., 1.);
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
      if (!groupedPlot) {
        gStyle->SetTitleXOffset(1.);
        gStyle->SetTitleYOffset(1.3);
      }
      gPad->SetTicks(1, 1);

      // make a histogram stack with the bg 
      THStack *bgStack = new THStack("bgStack" + histoSign[k] + nameSuffix[j], "Invariant Mass" + titleSuffix[j]);
      if (qcdEst > 0) bgStack->Add(emuMass_qcd.back());
      if (qcdEst != 2) bgStack->Add(emuMass_wjets.back());
      bgStack->Add(emuMass_zee.back());
      bgStack->Add(emuMass_zmumu.back());
      bgStack->Add(emuMass_tw.back());
      bgStack->Add(emuMass_zz.back());
      bgStack->Add(emuMass_wz.back());
      bgStack->Add(emuMass_ww.back());
      bgStack->Add(emuMass_ztautau.back());
      bgStack->Add(emuMass_ttbar.back());

      THStack *bgStackGrouped = new THStack("bgStackGrouped" + histoSign[k] + nameSuffix[j], "Invariant Mass" + titleSuffix[j]);
      bgStackGrouped->Add(emuMass_fakePairs.back());
      bgStackGrouped->Add(emuMass_ttLike.back());

      // bkg histogram for (data-bkg)/bkg plot
      TH1F *bgHist = (TH1F *)emuMass_ttbar.back()->Clone();
      if (qcdEst > 0) bgHist->Add(emuMass_qcd.back());
      if (qcdEst != 2) bgHist->Add(emuMass_wjets.back());
      bgHist->Add(emuMass_zee.back());
      bgHist->Add(emuMass_zmumu.back());
      bgHist->Add(emuMass_tw.back());
      bgHist->Add(emuMass_zz.back());
      bgHist->Add(emuMass_wz.back());
      bgHist->Add(emuMass_ww.back());
      bgHist->Add(emuMass_ztautau.back());

      // plot spectrum
      if (groupedPlot) {
        emuMass_ttLike.back()->SetFillColor(ttbarColour);
        emuMass_ttLike.back()->SetLineColor(kBlack);
        emuMass_ttLike.back()->SetLineWidth(2);
        //emuMass_ttLike.back()->DrawClone("HIST");
        emuMass_fakePairs.back()->SetFillColor(jetBkgColour);
        emuMass_fakePairs.back()->SetMarkerColor(jetBkgColour);
        emuMass_fakePairs.back()->SetLineColor(kBlack);
        emuMass_fakePairs.back()->SetLineWidth(2);
        //emuMass_fakePairs.back()->Draw("HISTsames");
        bgStackGrouped->Draw("hist");

        if (plotPullBelowSpec && j == 0) {
          bgStackGrouped->GetXaxis()->SetTitle("");
        } else {
          bgStackGrouped->GetXaxis()->SetTitle(xAxisTitle[k] + " [GeV]");
        }
        bgStackGrouped->GetXaxis()->SetTitleFont(font);
        bgStackGrouped->GetXaxis()->SetTitleSize(0.047);
        bgStackGrouped->GetXaxis()->SetTitleOffset(0.9);
        bgStackGrouped->GetXaxis()->SetLabelFont(font);
        bgStackGrouped->GetXaxis()->SetLabelSize(0.05);
        bgStackGrouped->GetXaxis()->SetMoreLogLabels();
        bgStackGrouped->GetXaxis()->SetNoExponent();
        bgStackGrouped->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax); 
        //bgStackGrouped->GetXaxis()->SetLimits(xRangeMin, xRangeMax); 
        if (j == 1) bgStackGrouped->GetYaxis()->SetTitle("Events #geq " + xAxisTitle[k]);
        else bgStackGrouped->GetYaxis()->SetTitle("Events / GeV");
        bgStackGrouped->GetYaxis()->SetTitleFont(font);
        bgStackGrouped->GetYaxis()->SetTitleSize(0.047);
        bgStackGrouped->GetYaxis()->SetTitleOffset(1.1);
        bgStackGrouped->GetYaxis()->SetLabelFont(font);
        bgStackGrouped->GetYaxis()->SetLabelSize(0.05);
        bgStackGrouped->SetMinimum(yRangeMin[k + j * 3]); 
        bgStackGrouped->SetMaximum(yRangeMax[k + j * 3]); 
      } else {
        emuMass_ttbar.back()->SetFillColor(ttbarColour);
        emuMass_ttbar.back()->SetLineColor(kBlack);
        emuMass_ttbar.back()->SetLineWidth(2);
        //emuMass_ttbar.back()->DrawClone("HIST");
        emuMass_ztautau.back()->SetFillColor(zttColour);
        emuMass_ztautau.back()->SetMarkerColor(zttColour);
        emuMass_ztautau.back()->SetLineColor(kBlack);
        emuMass_ztautau.back()->SetLineWidth(2);
        //emuMass_ztautau.back()->Draw("HISTsames");
        emuMass_ww.back()->SetFillColor(wwColour);
        emuMass_ww.back()->SetMarkerColor(wwColour);
        emuMass_ww.back()->SetLineColor(kBlack);
        emuMass_ww.back()->SetLineWidth(2);
        //emuMass_ww.back()->Draw("HISTsames");
        emuMass_wz.back()->SetFillColor(wzColour);
        emuMass_wz.back()->SetMarkerColor(wzColour);
        emuMass_wz.back()->SetLineColor(kBlack);
        emuMass_wz.back()->SetLineWidth(2);
        //emuMass_wz.back()->Draw("HISTsames");
        emuMass_zz.back()->SetFillColor(zzColour);
        emuMass_zz.back()->SetMarkerColor(zzColour);
        emuMass_zz.back()->SetLineColor(kBlack);
        emuMass_zz.back()->SetLineWidth(2);
        //emuMass_zz.back()->Draw("HISTsames");
        emuMass_tw.back()->SetFillColor(twColour);
        emuMass_tw.back()->SetMarkerColor(twColour);
        emuMass_tw.back()->SetLineColor(kBlack);
        emuMass_tw.back()->SetLineWidth(2);
        //emuMass_tw.back()->Draw("HISTsames");
        emuMass_zmumu.back()->SetFillColor(zmmColour);
        emuMass_zmumu.back()->SetMarkerColor(zmmColour);
        emuMass_zmumu.back()->SetLineColor(kBlack);
        emuMass_zmumu.back()->SetLineWidth(2);
        //emuMass_zmumu.back()->Draw("HISTsames");
        emuMass_zee.back()->SetFillColor(zeeColour);
        emuMass_zee.back()->SetMarkerColor(zeeColour);
        emuMass_zee.back()->SetLineColor(kBlack);
        emuMass_zee.back()->SetLineWidth(2);
        //emuMass_zee.back()->Draw("HISTsames");
        if (qcdEst != 2) {
          emuMass_wjets.back()->SetFillColor(wjetColour);
          emuMass_wjets.back()->SetMarkerColor(wjetColour);
          emuMass_wjets.back()->SetLineColor(kBlack);
          emuMass_wjets.back()->SetLineWidth(2);
          //emuMass_wjets.back()->Draw("HISTsames");
        }
        if (qcdEst > 0) {
          emuMass_qcd.back()->SetFillColor(jetBkgColour);
          emuMass_qcd.back()->SetMarkerColor(jetBkgColour);
          emuMass_qcd.back()->SetLineColor(kBlack);
          emuMass_qcd.back()->SetLineWidth(2);
          //emuMass_qcd.back()->Draw("HISTsames");
        }
        bgStack->Draw("hist");

        if (plotPullBelowSpec && j == 0) {
          bgStack->GetXaxis()->SetTitle("");
        } else {
          bgStack->GetXaxis()->SetTitle(xAxisTitle[k] + " [GeV]");
        }
        bgStack->GetXaxis()->SetTitleFont(font);
        bgStack->GetXaxis()->SetTitleSize(0.047);
        bgStack->GetXaxis()->SetTitleOffset(0.9);
        bgStack->GetXaxis()->SetLabelFont(font);
        bgStack->GetXaxis()->SetLabelSize(0.05);
        bgStack->GetXaxis()->SetMoreLogLabels();
        bgStack->GetXaxis()->SetNoExponent();
        bgStack->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax); 
        //bgStack->GetXaxis()->SetLimits(xRangeMin, xRangeMax); 
        if (j == 1) bgStack->GetYaxis()->SetTitle("Events #geq " + xAxisTitle[k]);
        else bgStack->GetYaxis()->SetTitle("Events / GeV");
        bgStack->GetYaxis()->SetTitleFont(font);
        bgStack->GetYaxis()->SetTitleSize(0.047);
        bgStack->GetYaxis()->SetTitleOffset(1.1);
        bgStack->GetYaxis()->SetLabelFont(font);
        bgStack->GetYaxis()->SetLabelSize(0.05);
        bgStack->SetMinimum(yRangeMin[k + j * 3]); 
        bgStack->SetMaximum(yRangeMax[k + j * 3]); 
      }

      emuMass_data.back()->SetLineWidth(1);
      emuMass_data.back()->SetLineColor(kBlack);
      emuMass_data.back()->SetMarkerStyle(20);
      emuMass_data.back()->SetMarkerSize(1.1);
      //emuMass_data.back()->Draw("e1 sames");

      float normalize = 1.;
      if (j > 0) normalize = -1.;
      TGraphAsymmErrors* dataHist = makeDataGraph(emuMass_data.back(), normalize, false);      
      dataHist->SetMarkerStyle(20);
      dataHist->SetMarkerSize(1.1);
      dataHist->Draw("PZ");

      // redraw axis
      emuMass_ttbar.back()->Draw("sameaxis");

      // legend and labels
      TLegend legend(0.731, 0.446, 0.921, 0.885);
      legend.SetTextFont(font);
      legend.SetTextSize(0.03);
      legend.SetBorderSize(0);
      legend.SetLineColor(1);
      legend.SetLineStyle(1);
      legend.SetLineWidth(1);
      legend.SetFillColor(19);
      legend.SetFillStyle(0);
      if (groupedPlot) {
        if (prelim) {
          legend.SetX1(0.460);
          legend.SetX2(0.720);
        } else {
          legend.SetX1(0.511);
          legend.SetX2(0.771);
        }
        legend.SetY1(0.643);
        legend.SetY2(0.830);
        legend.SetTextSize(0.047);
      }
      legend.AddEntry(emuMass_data.back(), "DATA");
      if (groupedPlot) {
        legend.AddEntry(emuMass_ttLike.back(), "t#bar{t}, tW, VV, #tau#tau" ,"F");
        legend.AddEntry(emuMass_fakePairs.back(), "#mu#mu, ee, W+jets, jets (data)" ,"F");
      } else {
        legend.AddEntry(emuMass_ttbar.back(), "t#bar{t}" ,"F");
        legend.AddEntry(emuMass_ztautau.back(), "#gamma/Z#rightarrow#tau#tau" ,"F");
        legend.AddEntry(emuMass_ww.back(), "WW" ,"F");
        legend.AddEntry(emuMass_wz.back(), "WZ" ,"F");
        legend.AddEntry(emuMass_zz.back(), "ZZ" ,"F");
        legend.AddEntry(emuMass_tw.back(), "tW" ,"F");
        legend.AddEntry(emuMass_zmumu.back(), "#gamma/Z#rightarrow#mu#mu" ,"F");
        legend.AddEntry(emuMass_zee.back(), "#gamma/Z#rightarrowee" ,"F");
        if (qcdEst != 2) legend.AddEntry(emuMass_wjets.back(), "W+jets" ,"F");
        if (qcdEst > 0) legend.AddEntry(emuMass_qcd.back(), "jets (data)" ,"F");
      }
      legend.DrawClone("sames");
      
      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(font);
      tex->SetLineWidth(2);
      if (groupedPlot) {
        tex->SetTextSize(0.047);
        if (prelim) tex->DrawLatex(0.467, 0.846, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
        else tex->DrawLatex(0.518, 0.846, "CMS, 8 TeV, 19.7 fb^{-1}");
        if (eRegion == 0) tex->DrawLatex(0.275, 0.846, "e in barrel");
        if (eRegion == 1) tex->DrawLatex(0.275, 0.846, "e in endcap");
      } else {
        tex->SetTextSize(0.042);
        if (prelim) tex->DrawLatex(0.325, 0.853, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
        else tex->DrawLatex(0.405, 0.853, "CMS, 8 TeV, 19.7 fb^{-1}");
        if (eRegion == 0) tex->DrawLatex(0.325, 0.775, "e in barrel");
        if (eRegion == 1) tex->DrawLatex(0.325, 0.775, "e in endcap");
      }

      // safe in various file formats
      stringstream sStream;
      if (!plotPullBelowSpec || j > 0) {
        sStream << plotDir << "emuSpec";
        if (k == 0) sStream << "_";
        sStream << histoSign[k];
        if (eRegion == 0) sStream << "EB_";
        if (eRegion == 1) sStream << "EE_";
        sStream << fileNameExtra << nameSuffix[j];
        if (j > 0) sStream << "_";
        if (groupedPlot) sStream << "grouped_";
        if (!logPlotY) sStream << "lin_";
        sStream << lumi->GetVal() << "pb-1";
        TString saveFileName = sStream.str();
        if ((j == 0 && saveSpec) || (j > 0 && saveCumSpec)) {
          if (saveAsPdf) emuPlot->Print(saveFileName + ".pdf", "pdf");
          if (saveAsPng) emuPlot->Print(saveFileName + ".png", "png");
          if (saveAsRoot) emuPlot->Print(saveFileName + ".root", "root");
        }
      }

      if (j > 0 || !plotPull) continue;
      // plot a (data - bg)/bg histogram
      dataOverBgHist->Add(bgHist, -1.);
      dataOverBgHist->Divide(bgHist);

      float fontScaleBot = 1.;
      TPad *pullPad = new TPad("pullPad" + histoSign[k] + nameSuffix[j], "(data - bg) / bg" + titleSuffix[j], 0., 0., 1., 1.);
      TCanvas *dataOverBgPlot;
      if (plotPullBelowSpec) {
        emuPlot->cd();
        pullPad->SetPad(0., 0., 1., 0.33);
        pullPad->SetBottomMargin(0.22);
        pullPad->SetTopMargin(0.);
        fontScaleBot = specPad->GetHNDC() / pullPad->GetHNDC();
      } else {
        dataOverBgPlot = new TCanvas("dataOverBgPlot" + histoSign[k] + nameSuffix[j], "(data - bg) / bg" + titleSuffix[j], 100, 100, 900, 600);
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

      dataOverBgHist->SetLineWidth(1);
      dataOverBgHist->SetLineColor(kBlack);
      dataOverBgHist->SetMarkerStyle(20);
      dataOverBgHist->SetMarkerSize(1.1);

      dataOverBgHist->GetXaxis()->SetTitle(xAxisTitle[k] + " [GeV]");
      dataOverBgHist->GetXaxis()->SetTitleFont(font);
      dataOverBgHist->GetXaxis()->SetTitleSize(0.047 * fontScaleBot);
      dataOverBgHist->GetXaxis()->SetTitleOffset(0.9);
      dataOverBgHist->GetXaxis()->SetLabelFont(font);
      dataOverBgHist->GetXaxis()->SetLabelSize(0.05 * fontScaleBot);
      dataOverBgHist->GetXaxis()->SetMoreLogLabels();
      dataOverBgHist->GetXaxis()->SetNoExponent();
      dataOverBgHist->GetXaxis()->SetRangeUser(xRangeMinRatio, xRangeMaxRatio);
      dataOverBgHist->GetYaxis()->SetTitle("(data-bkg)/bkg");
      dataOverBgHist->GetYaxis()->SetTitleFont(font);
      dataOverBgHist->GetYaxis()->SetTitleSize(0.047 * fontScaleBot);
      dataOverBgHist->GetYaxis()->SetTitleOffset(1.1 / fontScaleBot);
      dataOverBgHist->GetYaxis()->SetLabelFont(font);
      dataOverBgHist->GetYaxis()->SetLabelSize(0.05 * fontScaleBot);
      dataOverBgHist->GetYaxis()->SetRangeUser(yRangeMinRatio[k], yRangeMaxRatio[k]);

      dataOverBgHist->Draw();

      TF1 *f0 = new TF1("f0" + histoSign[k], "[0]");
      dataOverBgHist->Fit("f0" + histoSign[k], "", "", fitMin, fitMax);
      //dataOverBgHist->Draw("sameaxis");

      if (!plotPullBelowSpec) {
        if (groupedPlot) tex->SetTextSize(0.047 * fontScaleBot);
        else tex->SetTextSize(0.042 * fontScaleBot);
        if (prelim) tex->DrawLatex(0.150, 0.853, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
        else tex->DrawLatex(0.150, 0.853, "CMS, 8 TeV, 19.7 fb^{-1}");
        tex->SetTextSize(0.042 * fontScaleBot);
        tex->DrawLatex(0.150, 0.764, Form("#chi^{2} / ndf: %.2f / %i", f0->GetChisquare(), f0->GetNDF()));
        tex->DrawLatex(0.150, 0.720, Form("p0: %.4f #pm %0.4f", f0->GetParameter(0), f0->GetParError(0)));
        if (eRegion == 0) tex->DrawLatex(0.708, 0.853, "e in barrel");
        if (eRegion == 1) tex->DrawLatex(0.708, 0.853, "e in endcap");
      } else {
        tex->SetTextSize(0.042 * fontScaleBot);
        tex->DrawLatex(0.150, 0.875, Form("#chi^{2} / ndf: %.2f / %i", f0->GetChisquare(), f0->GetNDF()));
        tex->DrawLatex(0.150, 0.775, Form("p0: %.4f #pm %0.4f", f0->GetParameter(0), f0->GetParError(0)));
      }

      // safe in various file formats
      if (savePull && !plotPullBelowSpec) {
        sStream.str("");
        sStream << plotDir << "emuDataOverBkg";
        if (k == 0) sStream << "_";
        sStream << histoSign[k];
        if (eRegion == 0) sStream << "EB_";
        if (eRegion == 1) sStream << "EE_";
        sStream << fileNameExtra << lumi->GetVal() << "pb-1";
        TString saveFileName = sStream.str();
        if (saveAsPdf) dataOverBgPlot->Print(saveFileName + ".pdf", "pdf");
        if (saveAsPng) dataOverBgPlot->Print(saveFileName + ".png", "png");
        if (saveAsRoot) dataOverBgPlot->Print(saveFileName + ".root", "root");
      }
      if ((saveSpec || savePull) && plotPullBelowSpec) {
        sStream.str("");
        sStream << plotDir << "emuSpec";
        if (k == 0) sStream << "_";
        sStream << histoSign[k];
        if (eRegion == 0) sStream << "EB_";
        if (eRegion == 1) sStream << "EE_";
        sStream << fileNameExtra << nameSuffix[j];
        if (j > 0) sStream << "_";
        if (groupedPlot) sStream << "grouped_";
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
  } // end loop over full, SS and OS

  // generate one object containing everything
  vector<vector<TH1F *> > emuMasses;
  emuMasses.push_back(emuMass_data);
  emuMasses.push_back(emuMass_ttbar);
  emuMasses.push_back(emuMass_ztautau);
  emuMasses.push_back(emuMass_ww);
  emuMasses.push_back(emuMass_wz);
  emuMasses.push_back(emuMass_zz);
  emuMasses.push_back(emuMass_tw);
  emuMasses.push_back(emuMass_zmumu);
  emuMasses.push_back(emuMass_zee);
  if (qcdEst != 2) emuMasses.push_back(emuMass_wjets);
  if (qcdEst > 0) emuMasses.push_back(emuMass_qcd);

  // define groups of MC samples
  vector<bool> ttLikeSamples(6, true);
  vector<bool> contamSamples(6, false);
  contamSamples.push_back(true); // Zmm
  contamSamples.push_back(true); // Zee
  contamSamples.push_back(true); // WJets or QCD
  vector<bool> contamSamplesNoQcd(contamSamples);
  vector<bool> allSamples(9, true);
  vector<bool> onlyQCD(emuMasses.size() - 1, false);
  if (qcdEst > 0) {
    onlyQCD.back() = true;
    if (qcdEst != 2) {
      allSamples.push_back(true);
      contamSamples.push_back(true);
      contamSamplesNoQcd.push_back(false);
      systErrMC.push_back(0.); // QCD error will be calculated later
    } else {
      contamSamplesNoQcd.back() = false;
    }
  }
  vector<bool> allSamplesNoQcd(allSamples);
  if (qcdEst > 0) allSamplesNoQcd.back() = false;
  unsigned int qcdInd = onlyQCD.size();
  unsigned int qcdErrInd = qcdInd - 1;

  // calculate rate of syst errors
  float systErrLuEff = sqrt(systErrLumi*systErrLumi + systErrEff*systErrEff);
  vector<float> systErrMCLuEff;
  for (unsigned int it = 0; it < systErrMC.size(); ++it)
     systErrMCLuEff.push_back(sqrt(systErrMC[it]*systErrMC[it] + systErrLuEff*systErrLuEff));

  bool calcQcdErr = false;
  if (qcdEst == 1) calcQcdErr = true;

  //cout << "qcdInd " << qcdInd << ", emuMasses.size() " << emuMasses.size() << ", systErrMC.size() " << systErrMC.size() 
  //     << ", systErrMCLuEff.size() " << systErrMCLuEff.size() << ", allSamples.size() " << allSamples.size() 
  //     << ", allSamplesNoQcd.size() " << allSamplesNoQcd.size() << ", contamSamplesNoQcd.size() " << contamSamplesNoQcd.size() 
  //     << ", contamSamples.size() " << contamSamples.size() << ", onlyQCD.size() " << onlyQCD.size() << endl;
  //for (unsigned int sIt = 0; sIt < emuMasses.size() - 1; ++sIt) {
  //   cout << "allSamples " << allSamples[sIt] << ", allSamplesNoQcd " << allSamplesNoQcd[sIt] 
  //        << ", contamSamples " << contamSamples[sIt] << ", contamSamplesNoQcd " << contamSamplesNoQcd[sIt] 
  //        << ", onlyQCD " << onlyQCD[sIt] << ", systErrMC " << systErrMC[sIt] << ", systErrMCLuEff " << systErrMCLuEff[sIt] << endl;
  //}

  // define special bins corresponding to specific masses
  int bin60 = emuMass_data.at(ALL)->FindBin(60.);
  int bin120 = emuMass_data.at(ALL)->FindBin(120.);
  int bin200 = emuMass_data.at(ALL)->FindBin(200.); 
  int bin400 = emuMass_data.at(ALL)->FindBin(400.); 
  int bin500 = emuMass_data.at(ALL)->FindBin(500.); 

  vector<const char *> sampleNames;
  sampleNames.push_back("data   ");
  sampleNames.push_back("ttbar  ");
  sampleNames.push_back("Ztautau");
  sampleNames.push_back("WW     ");
  sampleNames.push_back("WZ     ");
  sampleNames.push_back("ZZ     ");
  sampleNames.push_back("tW     ");
  sampleNames.push_back("Zmumu  ");
  sampleNames.push_back("Zee    ");
  if (qcdEst != 2) sampleNames.push_back("WJets  ");
  if (qcdEst > 0) sampleNames.push_back("QCD    ");

  // write numbers
  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "HEEP - TIGHT MU        Lumi        = " << lumi->GetVal() << "pb-1" << endl;
  //cout << "                       e pT EB     > " << bar_et << "GeV/c" << endl;
  //cout << "                       e pT EE     > " << end_et << "GeV/c" << endl;
  //cout << "                       mu pT       > " << muon_et << "GeV/c" << endl;
  //cout << "                       mu |eta|    < " << muon_etaMax << endl;
  cout << endl;
  cout << "Systematic errors" << endl;
  cout << " Luminosity:  " << systErrLumi * 100 << "%" << endl;
  cout << " Efficiency:  " << systErrEff * 100 << "%" << endl;
  cout << " ttbar:       " << systErrMC[TTBAR-1] * 100 << "%" << endl;
  cout << " Z->tautau:   " << systErrMC[ZTT-1] * 100 << "%" << endl;
  cout << " WW:          " << systErrMC[WW-1] * 100 << "%" << endl;
  cout << " WZ:          " << systErrMC[WZ-1] * 100 << "%" << endl;
  cout << " ZZ:          " << systErrMC[ZZ-1] * 100 << "%" << endl;
  cout << " tW, tbarW:   " << systErrMC[TW-1] * 100 << "%" << endl;
  cout << " Z->mumu:     " << systErrMC[ZMM-1] * 100 << "%" << endl;
  cout << " Z->ee:       " << systErrMC[ZEE-1] * 100 << "%" << endl;
  if (qcdEst != 2) cout << " W+Jets:      " << systErrMC[WJET-1] * 100 << "%" << endl;
  else cout << " QCD:        " << systErrMC.back() * 100 << "%" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  for (unsigned int signIt = 1; signIt < 6; signIt += 2) {
    if (signIt == 3) cout << "-SS--------------------------------------------------------------------------------------------------------" << endl;
    if (signIt == 5) cout << "-OS--------------------------------------------------------------------------------------------------------" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "M_emu         |         >  60GeV/c^2          |        > 120GeV/c^2          |        > 200GeV/c^2         |        > 400GeV/c^2          |" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  
    printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
           emuMass_data.at(signIt)->GetBinContent(bin60), sqrt(emuMass_data.at(signIt)->GetBinContent(bin60)),         
           emuMass_data.at(signIt)->GetBinContent(bin120), sqrt(emuMass_data.at(signIt)->GetBinContent(bin120)),
           emuMass_data.at(signIt)->GetBinContent(bin200), sqrt(emuMass_data.at(signIt)->GetBinContent(bin200)),
           emuMass_data.at(signIt)->GetBinContent(bin400), sqrt(emuMass_data.at(signIt)->GetBinContent(bin400)));
    cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
    for (unsigned int sampleIt = 1; sampleIt < sampleNames.size(); ++sampleIt) {
      if (sampleIt == 7) cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      if (qcdEst == 1 && sampleIt == sampleNames.size() - 1) {
        printf("nb %7s    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", sampleNames[sampleIt],
               emuMass_qcd.at(signIt)->GetBinContent(bin60), emuMass_qcd.at(signIt)->GetBinContent(bin60) * CalcSSQcdErr(emuMasses, systErrMCLuEff, bin60), 
               emuMass_qcd.at(signIt)->GetBinContent(bin120), emuMass_qcd.at(signIt)->GetBinContent(bin120) * CalcSSQcdErr(emuMasses, systErrMCLuEff, bin120), 
               emuMass_qcd.at(signIt)->GetBinContent(bin200), emuMass_qcd.at(signIt)->GetBinContent(bin200) * CalcSSQcdErr(emuMasses, systErrMCLuEff, bin200),
               emuMass_qcd.at(signIt)->GetBinContent(bin400), emuMass_qcd.at(signIt)->GetBinContent(bin400) * CalcSSQcdErr(emuMasses, systErrMCLuEff, bin400));
      } else {
        printf("nb %7s    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", sampleNames[sampleIt],
               emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin60), emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin60) * systErrMCLuEff[sampleIt-1], 
               emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin120), emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin120) * systErrMCLuEff[sampleIt-1], 
               emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin200), emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin200) * systErrMCLuEff[sampleIt-1],
               emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin400), emuMasses.at(sampleIt).at(signIt)->GetBinContent(bin400) * systErrMCLuEff[sampleIt-1]);
      }
    }
    cout << endl;
    cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
    printf("TOT ttlike    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
           CalcBgSum(emuMasses, ttLikeSamples, signIt, bin60), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, signIt, bin60),
           CalcBgSum(emuMasses, ttLikeSamples, signIt, bin120), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, signIt, bin120),
           CalcBgSum(emuMasses, ttLikeSamples, signIt, bin200), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, signIt, bin200),
           CalcBgSum(emuMasses, ttLikeSamples, signIt, bin400), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, signIt, bin400));
    printf("TOT contam    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
           CalcBgSum(emuMasses, contamSamples, signIt, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, contamSamples, signIt, bin60, -1, calcQcdErr),
           CalcBgSum(emuMasses, contamSamples, signIt, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, contamSamples, signIt, bin120, -1, calcQcdErr),
           CalcBgSum(emuMasses, contamSamples, signIt, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, contamSamples, signIt, bin200, -1, calcQcdErr),
           CalcBgSum(emuMasses, contamSamples, signIt, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, contamSamples, signIt, bin400, -1, calcQcdErr));
    cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  
    printf("TOT Bkg       | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
           CalcBgSum(emuMasses, allSamples, signIt, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin60, -1, calcQcdErr),
           CalcBgSum(emuMasses, allSamples, signIt, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin120, -1, calcQcdErr),
           CalcBgSum(emuMasses, allSamples, signIt, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin200, -1, calcQcdErr),
           CalcBgSum(emuMasses, allSamples, signIt, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin400, -1, calcQcdErr));
    cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << endl << endl;
  }
  cout << endl;

  cout << "--Without adding QCD contribution:--------------------------------------------------------------------------------------------------------" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |         > 60GeV/c^2          |        > 120GeV/c^2          |         > 200GeV/c^2         |         > 400GeV/c^2         |" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  for (unsigned int signIt = 1; signIt < 6; signIt += 2) {
    if (signIt == 3) cout << "-SS-------------------------------------------------------------------------------------------------------------------------------------" << endl;
    if (signIt == 5) cout << "-OS-------------------------------------------------------------------------------------------------------------------------------------" << endl;
    printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
           emuMasses.at(DATA).at(signIt)->GetBinContent(bin60), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin60)),
           emuMasses.at(DATA).at(signIt)->GetBinContent(bin120), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin120)),
           emuMasses.at(DATA).at(signIt)->GetBinContent(bin200), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin200)),
           emuMasses.at(DATA).at(signIt)->GetBinContent(bin400), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin400)));
    printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
           CalcBgSum(emuMasses, allSamplesNoQcd, signIt, bin60), CalcSystErr(emuMasses, systErrMCLuEff, allSamplesNoQcd, signIt, bin60),
           CalcBgSum(emuMasses, allSamplesNoQcd, signIt, bin120), CalcSystErr(emuMasses, systErrMCLuEff, allSamplesNoQcd, signIt, bin120),
           CalcBgSum(emuMasses, allSamplesNoQcd, signIt, bin200), CalcSystErr(emuMasses, systErrMCLuEff, allSamplesNoQcd, signIt, bin200),
           CalcBgSum(emuMasses, allSamplesNoQcd, signIt, bin400), CalcSystErr(emuMasses, systErrMCLuEff, allSamplesNoQcd, signIt, bin400));
  }
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;

  if (qcdEst == 1) {
    //systErrMC.back() = 2 * sqrt(emuMasses.at(DATA).at(SS)->Integral() + pow(CalcSystErr(emuMasses, systErrMCLuEff, allSamplesNoQcd, SS, 1), 2)) / emuMasses.at(qcdInd).at(ALL)->Integral();
    //systErrMC.back() = CalcSystErr(emuMasses, systErrMCLuEff, allSamplesNoQcd, SSCUM, 1) / emuMass_qcd.at(SSCUM)->GetBinContent(1);
    //systErrMCLuEff.back() = systErrMC[qcdErrInd];

    cout << endl;
      cout << "---QCD events from SS spectrum:----------------------------------------------------------------------------------------------------------------------------------" << endl;
      cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
      printf("nb QCD SS+OS  | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) |\n",
             emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin60, -1, calcQcdErr), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin60, -1, calcQcdErr) / emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin60),
             emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin120, -1, calcQcdErr), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin120, -1, calcQcdErr) / emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin120),
             emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin200, -1, calcQcdErr), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin200, -1, calcQcdErr) / emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin200),
             emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin400, -1, calcQcdErr), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin400, -1, calcQcdErr) / emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin400));
      printf("%% of total MC |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |\n",
             100 * emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin60) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin60), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin60, -1, calcQcdErr) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin60),
             100 * emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin120) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin120), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin120, -1, calcQcdErr) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin120),
             100 * emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin200) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin200), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin200, -1, calcQcdErr) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin200),
             100 * emuMasses.at(qcdInd).at(ALLCUM)->GetBinContent(bin400) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin400), 
             100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALLCUM, bin400, -1, calcQcdErr) / CalcBgSum(emuMasses, allSamplesNoQcd, ALLCUM, bin400));
      cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  }

  // top up bg contribution with qcd
  if (qcdEst > 0) {
    cout << endl;
    cout << "--After adding QCD contribution:----------------------------------------------------------------------------------------------------------" << endl;
    cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "M_emu         |         > 60GeV/c^2          |        > 120GeV/c^2          |         > 200GeV/c^2         |         > 400GeV/c^2         |" << endl;
    cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    for (unsigned int signIt = 1; signIt < 6; signIt += 2) {
      if (signIt == 3) cout << "-SS-------------------------------------------------------------------------------------------------------------------------------------" << endl;
      if (signIt == 5) cout << "-OS-------------------------------------------------------------------------------------------------------------------------------------" << endl;
      printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
              emuMasses.at(DATA).at(signIt)->GetBinContent(bin60), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin60)),
              emuMasses.at(DATA).at(signIt)->GetBinContent(bin120), sqrt((emuMasses.at(DATA).at(signIt))->GetBinContent(bin120)),
              emuMasses.at(DATA).at(signIt)->GetBinContent(bin200), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin200)),
              emuMasses.at(DATA).at(signIt)->GetBinContent(bin400), sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin400)));
      printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
              CalcBgSum(emuMasses, allSamples, signIt, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin60, -1, calcQcdErr),
              CalcBgSum(emuMasses, allSamples, signIt, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin120, -1, calcQcdErr),
              CalcBgSum(emuMasses, allSamples, signIt, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin200, -1, calcQcdErr),
              CalcBgSum(emuMasses, allSamples, signIt, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin400, -1, calcQcdErr));
    }
    cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  }

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |        60 - 120GeV/c^2       |      120 - 200GeV/c^2        |       200 - 400GeV/c^2       |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  for (unsigned int signIt = 1; signIt < 6; signIt += 2) {
    printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
            emuMasses.at(DATA).at(signIt)->GetBinContent(bin60) - emuMasses.at(DATA).at(signIt)->GetBinContent(bin120), 
            sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin60) - emuMasses.at(DATA).at(signIt)->GetBinContent(bin120)),
            emuMasses.at(DATA).at(signIt)->GetBinContent(bin120) - emuMasses.at(DATA).at(signIt)->GetBinContent(bin200), 
            sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin120) - emuMasses.at(DATA).at(signIt)->GetBinContent(bin200)),
            emuMasses.at(DATA).at(signIt)->GetBinContent(bin200) - emuMasses.at(DATA).at(signIt)->GetBinContent(bin400), 
            sqrt(emuMasses.at(DATA).at(signIt)->GetBinContent(bin200) - emuMasses.at(DATA).at(signIt)->GetBinContent(bin400)));
    printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
            CalcBgSum(emuMasses, allSamples, signIt, bin60, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin60, bin120, calcQcdErr),
            CalcBgSum(emuMasses, allSamples, signIt, bin120, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin120, bin200, calcQcdErr),
            CalcBgSum(emuMasses, allSamples, signIt, bin200, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, signIt, bin200, bin400, calcQcdErr));
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  }

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "|Event yield table                                                                                        |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "\\begin{table}[tbh]" << endl;
  cout << "\\centering" << endl;
  cout << "\\caption{Number of $e\\mu$ events with different charge combinations from data and Monte Carlo simulation. The listed errors are the systematic errors}" << endl;
  cout << "\\begin{tabular}{|c|c|c|c|c|c|c|}" << endl;
  cout << "\\hline" << endl;
  cout << " & \\multicolumn{6}{c|}{number of events} \\\\" << endl;
  cout << "$m_{e\\mu}$ &  \\multicolumn{2}{c|}{same-sign}  & \\multicolumn{2}{c|}{opposite-sign} & \\multicolumn{2}{c|}{combined}  \\\\" << endl;
  cout << " &  data & MC & data & MC & data & MC \\\\" << endl;
  cout << "\\hline" << endl;
  printf("$>$ 60~$\\gevsq$   & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin60), CalcBgSum(emuMasses, allSamples, SSCUM, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin60, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin60), CalcBgSum(emuMasses, allSamples, OSCUM, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin60, -1, calcQcdErr), 
          emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin60), CalcBgSum(emuMasses, allSamples, ALLCUM, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin60, -1, calcQcdErr));
  printf("$>$ 120~$\\gevsq$  & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin120), CalcBgSum(emuMasses, allSamples, SSCUM, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin120, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin120), CalcBgSum(emuMasses, allSamples, OSCUM, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin120, -1, calcQcdErr), 
          emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin120), CalcBgSum(emuMasses, allSamples, ALLCUM, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin120, -1, calcQcdErr));
  printf("$>$ 200~$\\gevsq$  & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin200), CalcBgSum(emuMasses, allSamples, SSCUM, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin200, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin200), CalcBgSum(emuMasses, allSamples, OSCUM, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin200, -1, calcQcdErr), 
          emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin200), CalcBgSum(emuMasses, allSamples, ALLCUM, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin200, -1, calcQcdErr));
  printf("$>$ 400~$\\gevsq$  & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin400), CalcBgSum(emuMasses, allSamples, SSCUM, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin400, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin400), CalcBgSum(emuMasses, allSamples, OSCUM, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin400, -1, calcQcdErr), 
          emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin400), CalcBgSum(emuMasses, allSamples, ALLCUM, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin400, -1, calcQcdErr));
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\label{tab:emu_event_yield}" << endl;
  cout << "\\end{table}" << endl;

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "|Event yield table with combined = SS + OS                                                                |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "\\begin{table}[tbh]" << endl;
  cout << "\\centering" << endl;
  cout << "\\caption{Number of $e\\mu$ events with different charge combinations from data and Monte Carlo simulation. The listed errors are the systematic errors}" << endl;
  cout << "\\begin{tabular}{|c|c|c|c|c|c|c|}" << endl;
  cout << "\\hline" << endl;
  cout << " & \\multicolumn{6}{c|}{number of events} \\\\" << endl;
  cout << "$m_{e\\mu}$ &  \\multicolumn{2}{c|}{same-sign}  & \\multicolumn{2}{c|}{opposite-sign} & \\multicolumn{2}{c|}{combined}  \\\\" << endl;
  cout << " &  data & MC & data & MC & data & MC \\\\" << endl;
  cout << "\\hline" << endl;
  printf("$>$ 60~$\\gevsq$   & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin60), CalcBgSum(emuMasses, allSamples, SSCUM, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin60, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin60), CalcBgSum(emuMasses, allSamples, OSCUM, bin60), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin60, -1, calcQcdErr), 
          floor(0.5 + emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin60)) + floor(0.5 + emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin60)), 
          floor(0.5 + CalcBgSum(emuMasses, allSamples, SSCUM, bin60)) + floor(0.5 + CalcBgSum(emuMasses, allSamples, OSCUM, bin60)), 
          CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin60, -1, calcQcdErr)); 
  printf("$>$ 120~$\\gevsq$  & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin120), CalcBgSum(emuMasses, allSamples, SSCUM, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin120, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin120), CalcBgSum(emuMasses, allSamples, OSCUM, bin120), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin120, -1, calcQcdErr), 
          floor(0.5 + emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin120)) + floor(0.5 + emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin120)), 
          floor(0.5 + CalcBgSum(emuMasses, allSamples, SSCUM, bin120)) + floor(0.5 + CalcBgSum(emuMasses, allSamples, OSCUM, bin120)), 
          CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin120, -1, calcQcdErr)); 
  printf("$>$ 200~$\\gevsq$  & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin200), CalcBgSum(emuMasses, allSamples, SSCUM, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin200, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin200), CalcBgSum(emuMasses, allSamples, OSCUM, bin200), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin200, -1, calcQcdErr), 
          floor(0.5 + emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin200)) + floor(0.5 + emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin200)), 
          floor(0.5 + CalcBgSum(emuMasses, allSamples, SSCUM, bin200)) + floor(0.5 + CalcBgSum(emuMasses, allSamples, OSCUM, bin200)), 
          CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin200, -1, calcQcdErr));
  printf("$>$ 400~$\\gevsq$  & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f & %.0f & %.0f $\\pm$ %.0f \\\\\n", 
          emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin400), CalcBgSum(emuMasses, allSamples, SSCUM, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, SSCUM, bin400, -1, calcQcdErr), 
          emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin400), CalcBgSum(emuMasses, allSamples, OSCUM, bin400), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OSCUM, bin400, -1, calcQcdErr), 
          floor(0.5 + emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin400)) + floor(0.5 + emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin400)), 
          floor(0.5 + CalcBgSum(emuMasses, allSamples, SSCUM, bin400)) + floor(0.5 + CalcBgSum(emuMasses, allSamples, OSCUM, bin400)), 
          CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALLCUM, bin400, -1, calcQcdErr));
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\label{tab:emu_event_yield}" << endl;
  cout << "\\end{table}" << endl;

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "|Event yield table per sample                                                                             |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "\\begin{table}[b]" << endl;
  cout << "\\centering" << endl;
  cout << "\\caption{Number of $e\\mu$ events with invariant mass in different regions and with different charge combinations.}" << endl;
  cout << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" << endl;
  cout << "\\hline\\hline" << endl;
  cout << "Source & \\multicolumn{9}{c|}{number of events} \\\\" << endl;
  cout << " &  \\multicolumn{3}{c|}{[$120-200$]~GeV/c$^2$}  & \\multicolumn{3}{c|}{[$200-400$]~GeV/c$^2$} & \\multicolumn{3}{|c|}{$>$ 400~GeV/c$^2$}  \\\\" << endl;
  cout << " &  OS & SS & Combined & OS & SS & Combined & OS & SS & Combined \\\\\\hline" << endl;

  printf("CMS data  & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f \\\\\n", 
         emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin120) - emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin200), emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin120) - emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin200), emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin120) - emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin200), 
         emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin200) - emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin400), emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin200) - emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin400), emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin200) - emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin400), 
         emuMasses.at(DATA).at(OSCUM)->GetBinContent(bin400), emuMasses.at(DATA).at(SSCUM)->GetBinContent(bin400), emuMasses.at(DATA).at(ALLCUM)->GetBinContent(bin400));
  printf("Total Bkg & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         CalcBgSum(emuMasses, allSamples, OSCUM, bin120, bin200), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, OSCUM, bin120, bin200, calcQcdErr), 
         CalcBgSum(emuMasses, allSamples, SSCUM, bin120, bin200), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, SSCUM, bin120, bin200, calcQcdErr), 
         CalcBgSum(emuMasses, allSamples, ALLCUM, bin120, bin200), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, ALLCUM, bin120, bin200, calcQcdErr),
         CalcBgSum(emuMasses, allSamples, OSCUM, bin200, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, OSCUM, bin200, bin400, calcQcdErr), 
         CalcBgSum(emuMasses, allSamples, SSCUM, bin200, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, SSCUM, bin200, bin400, calcQcdErr), 
         CalcBgSum(emuMasses, allSamples, ALLCUM, bin200, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, ALLCUM, bin200, bin400, calcQcdErr),
         CalcBgSum(emuMasses, allSamples, OSCUM, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, OSCUM, bin400, -1, calcQcdErr), 
         CalcBgSum(emuMasses, allSamples, SSCUM, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, SSCUM, bin400, -1, calcQcdErr), 
         CalcBgSum(emuMasses, allSamples, ALLCUM, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, ALLCUM, bin400, -1, calcQcdErr)); 
  cout << "\\hline\\hline" << endl;
  printf("\\ttbar +  \\ttbar-like & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         CalcBgSum(emuMasses, ttLikeSamples, OSCUM, bin120, bin200), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, OSCUM, bin120, bin200),
         CalcBgSum(emuMasses, ttLikeSamples, SSCUM, bin120, bin200), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, SSCUM, bin120, bin200),
         CalcBgSum(emuMasses, ttLikeSamples, ALLCUM, bin120, bin200), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALLCUM, bin120, bin200),
         CalcBgSum(emuMasses, ttLikeSamples, OSCUM, bin200, bin400), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, OSCUM, bin200, bin400),
         CalcBgSum(emuMasses, ttLikeSamples, SSCUM, bin200, bin400), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, SSCUM, bin200, bin400),
         CalcBgSum(emuMasses, ttLikeSamples, ALLCUM, bin200, bin400), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALLCUM, bin200, bin400),
         CalcBgSum(emuMasses, ttLikeSamples, OSCUM, bin400), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, OSCUM, bin400),
         CalcBgSum(emuMasses, ttLikeSamples, SSCUM, bin400), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, SSCUM, bin400),
         CalcBgSum(emuMasses, ttLikeSamples, ALLCUM, bin400), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALLCUM, bin400));

  printf("contaminations        & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         CalcBgSum(emuMasses, contamSamplesNoQcd, OSCUM, bin120, bin200), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, OSCUM, bin120, bin200),
         CalcBgSum(emuMasses, contamSamplesNoQcd, SSCUM, bin120, bin200), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, SSCUM, bin120, bin200),
         CalcBgSum(emuMasses, contamSamplesNoQcd, ALLCUM, bin120, bin200), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, ALLCUM, bin120, bin200),
         CalcBgSum(emuMasses, contamSamplesNoQcd, OSCUM, bin200, bin400), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, OSCUM, bin200, bin400),
         CalcBgSum(emuMasses, contamSamplesNoQcd, SSCUM, bin200, bin400), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, SSCUM, bin200, bin400),
         CalcBgSum(emuMasses, contamSamplesNoQcd, ALLCUM, bin200, bin400), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, ALLCUM, bin200, bin400),
         CalcBgSum(emuMasses, contamSamplesNoQcd, OSCUM, bin400), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, OSCUM, bin400),
         CalcBgSum(emuMasses, contamSamplesNoQcd, SSCUM, bin400), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, SSCUM, bin400),
         CalcBgSum(emuMasses, contamSamplesNoQcd, ALLCUM, bin400), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, ALLCUM, bin400));

  if (qcdEst > 0) {
    printf("multi-jet             & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
           CalcBgSum(emuMasses, onlyQCD, OSCUM, bin120, bin200), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, OSCUM, bin120, bin200, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, SSCUM, bin120, bin200), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, SSCUM, bin120, bin200, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, ALLCUM, bin120, bin200), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, ALLCUM, bin120, bin200, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, OSCUM, bin200, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, OSCUM, bin200, bin400, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, SSCUM, bin200, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, SSCUM, bin200, bin400, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, ALLCUM, bin200, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, ALLCUM, bin200, bin400, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, OSCUM, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, OSCUM, bin400, -1, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, SSCUM, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, SSCUM, bin400, -1, calcQcdErr),
           CalcBgSum(emuMasses, onlyQCD, ALLCUM, bin400), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, qcdInd, ALLCUM, bin400, -1, calcQcdErr));
  }

  cout << "\\hline\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\label{tab:emu_event_yield}" << endl;
  cout << "\\end{table}" << endl;
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

TGraphAsymmErrors* makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBars)
{
  std::vector<double> xPoint,yPoint,xErrLow,xErrHigh,yErrLow,yErrHigh;
  for(int binNr=1;binNr<=dataHist->GetNbinsX();binNr++){
    double nrData = dataHist->GetBinContent(binNr);

    float scale = 1;
    if(normToBinWidth>0) scale= normToBinWidth/dataHist->GetBinWidth(binNr);

    const double alpha = (1 - 0.6827)/2;
    const double beta  = (1 - 0.6827)/2;

    double dataLowBound=0;
    double dataHighBound=0;
    if(nrData!=0){
      dataLowBound = 0.5*ROOT::Math::chisquared_quantile_c(1-alpha, 2*nrData);
      dataHighBound = 0.5*ROOT::Math::chisquared_quantile_c(beta, 2*(nrData+1));
    }
    double binCentre= dataHist->GetBinLowEdge(binNr)+0.5*dataHist->GetBinWidth(binNr);
    xPoint.push_back(binCentre);
    if(xErrBars){
      xErrLow.push_back(dataHist->GetBinWidth(binNr)*0.5);
      xErrHigh.push_back(dataHist->GetBinWidth(binNr)*0.5);
    }else{
      xErrHigh.push_back(0);
      xErrLow.push_back(0);
    }
    yPoint.push_back(nrData*scale);
    yErrLow.push_back((nrData-dataLowBound)*scale);
    yErrHigh.push_back((dataHighBound-nrData)*scale);
  }
  TGraphAsymmErrors* resultGraph = new TGraphAsymmErrors(xPoint.size(),&xPoint[0],&yPoint[0],&xErrLow[0],&xErrHigh[0],&yErrLow[0],&yErrHigh[0]);
   return resultGraph;
}

// flags: [apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
TH1F *
MakeHistoFromBranch(TFile *input, const char *treeName, const char *brName, int signs, int &region, const char *cutVariable, float cutLow, float cutHigh, vector<float> &binning, unsigned int flags, bool normToBinWidth, float userScale)
{
  input->cd();
  // prepare the histogram
  TH1F *histo = new TH1F("dummy", "dummy", 150, 0., 1500.);
  histo->Sumw2();
  histo->SetName(brName);
  histo->SetTitle(brName);
  float *bins = &binning[0];
  histo->GetXaxis()->Set(binning.size() - 1, bins);

  if (flags & 1<<7) userScale *= ((TParameter<float> *)input->Get("lumi"))->GetVal();
  if (flags & 1<<6) {
    THashList *mcWeights = (THashList *)input->Get("mcWeights");
    unsigned int charOffset = 8;
    if (flags & 1<<9) charOffset += 2;
    userScale *= ((TParameter<float> *)mcWeights->FindObject(treeName + charOffset))->GetVal();
  }
  float lumiScaleFactorEB = ((TParameter<float> *)input->Get("lumiScaleFactorEB"))->GetVal();
  float lumiScaleFactorEE = ((TParameter<float> *)input->Get("lumiScaleFactorEE"))->GetVal();

  //cout << "Scalefactor = " << userScale << endl;

  // get the tree
  TTree *tree;
  tree = (TTree *)input->Get(treeName);

  // get branches
  float var;
  bool passTrg;
  bool passHeep;
  float puWeight = 1.;
  int eCharge;
  int muCharge;
  int evtRegion;
  float cutVar = 0.;
  float fakeRate = 0.;
  float trgEff = 1.;
  float trgEffSf = 1.;
  float eleEffSf = 1.;
  float muEffSf = 1.;
  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus(brName,1);
  tree->SetBranchAddress(brName, &var);
  if (flags & 1<<8) {
    tree->SetBranchStatus("passTrg",1);
    tree->SetBranchAddress("passTrg", &passTrg);
  }
  if (signs != 0) {
    tree->SetBranchStatus("eCharge",1);
    tree->SetBranchStatus("muCharge",1);
    tree->SetBranchAddress("eCharge", &eCharge);
    tree->SetBranchAddress("muCharge", &muCharge);
  }
  if (region < 2) {
    tree->SetBranchStatus("evtRegion",1);
    tree->SetBranchAddress("evtRegion", &evtRegion);
  }
  if (flags & 1) {
    tree->SetBranchStatus("puWeight",1);
    tree->SetBranchAddress("puWeight", &puWeight);
  }
  if (cutVariable[0] != '\0') {
    tree->SetBranchStatus(cutVariable,1);
    tree->SetBranchAddress(cutVariable, &cutVar);
  }
  if (flags & 1<<5) {
    tree->SetBranchStatus("trgEff",1);
    tree->SetBranchAddress("trgEff", &trgEff);
  }
  if (flags & 1<<4) {
    tree->SetBranchStatus("trgEffSf",1);
    tree->SetBranchAddress("trgEffSf", &trgEffSf);
  }
  if (flags & 1<<2) {
    tree->SetBranchStatus("eleEffSf",1);
    tree->SetBranchAddress("eleEffSf", &eleEffSf);
  }
  if (flags & 1<<1) {
    tree->SetBranchStatus("muEffSf",1);
    tree->SetBranchAddress("muEffSf", &muEffSf);
  }
  if (flags & 1<<9) {
    tree->SetBranchStatus("passHeep",1);
    tree->SetBranchStatus("fakeRate",1);
    tree->SetBranchAddress("passHeep", &passHeep);
    tree->SetBranchAddress("fakeRate", &fakeRate);
  }

  Long64_t nEntries = (*tree).GetEntries();
  for (unsigned int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // trigger fired?
    if ((flags & 1<<8) && passTrg == false) continue;

    // select electron region
    if (evtRegion == 0 && region == 1) continue;
    if (evtRegion == 1 && region == 0) continue;

    float scaleFactor = userScale;
    if (flags & 1<<5) scaleFactor *= trgEff;
    if (flags & 1<<4) scaleFactor *= trgEffSf;
    if (flags & 1<<2) scaleFactor *= eleEffSf;
    if (flags & 1<<1) scaleFactor *= muEffSf;
    // set lumi according to detector region
    if (evtRegion == 0 && flags & 1<<3) scaleFactor *= lumiScaleFactorEB;
    if (evtRegion == 1 && flags & 1<<3) scaleFactor *= lumiScaleFactorEE;

    // PU reweight
    if (flags & 1) scaleFactor *= puWeight;

    // get only the desired charge combination. Scheme emu -3 to +3: -+, +-, OS, ALL, SS, ++, --
    if (signs < 0 && (eCharge * muCharge) > 0) continue; // OS
    if (signs > 0 && (eCharge * muCharge) < 0) continue; // SS
    if (abs(signs) == 3 && eCharge > 0) continue; // e-mu+ or e-mu-
    if (abs(signs) == 2 && eCharge < 0) continue; // e+mu- or e+mu+

    // user defined cut
    if (cutVariable[0] != '\0')
      if (cutVar < cutLow || cutVar >= cutHigh) continue;
    
    if (flags & 1<<9) {
      if (!passHeep) scaleFactor *= fakeRate / (1 - fakeRate);
      else continue;
    }

    if (normToBinWidth) scaleFactor /= histo->GetBinWidth(histo->FindBin(var));
    histo->Fill(var, scaleFactor);
  }

  //cout << "integral: " << histo->Integral() << "       overflow: " << histo->GetBinContent(histo->GetNbinsX() + 1) << endl;
  return histo;
}

