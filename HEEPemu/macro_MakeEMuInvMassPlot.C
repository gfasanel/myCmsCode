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
TGraphAsymmErrors * makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBars);

void macro_MakeEMuInvMassPlot()
{
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("./emuSpec_19703pb-1.root", "open");
  TFile input("./emuSpec_singleMuTrg_19706pb-1.root", "open");
  input.cd();

  TParameter<float> *lumi = (TParameter<float> *)input.Get("lumi");

  //const bool usePu = 1;
  const bool topReweighting = 0;
  const int qcdEst = 2; // estimation method of QCD contribution. none(0), from SS spectrum(1), from fake rate(2)

  const bool ttbar_sample_type = 0; // 0=powheg, 1=madgraph
  int eRegion = 2; // electron region EB(0), EE(1), EB+EE(2)

  bool plotSign[5];
  plotSign[0] = 1;  // all
  plotSign[1] = 0;  // SS same sign
  plotSign[2] = 0;  // OS opposite sign
  plotSign[3] = 0;  // e-mu+
  plotSign[4] = 0;  // e+mu-

  bool plotType[2];
  plotType[0] = 1;  // emu spectrum
  plotType[1] = 0;  // cumulative emu spectrum

  float plotSignal[4] = {100., 100., 100., 100.}; // signal scale factors. 0 for off
  const bool plotPull = 1; // plot (data-bkg)/bkg
  const bool plotPullBelowSpec = 1; // plot (data-bkg)/bkg below spectrum
  const bool varBinning = 1;
  const bool logPlotX = 0;
  const bool logPlotY = 1;
  const bool prelim = 1;
  const bool groupedPlot = 0;
  const bool overflowBin = 0;

  float xRangeMin = 60.;
  //float xRangeMax = 1200.;
  //float xRangeMin = 0.;
  float xRangeMax = 1500.;
  float yRangeMin[10] = {0.002, 0.002, 0.002,  0.002,  0.0002, 0.2, 0.2, 0.2, 0.2, 0.2};
  float yRangeMax[10] = {250, 40, 250, 250, 250, 30000, 4000, 30000, 30000, 30000};
  float yRangeMinRatio[5] = {-0.7, -0.7, -0.7, -0.7, -0.7};
  float yRangeMaxRatio[5] = {0.7, 0.7, 0.7, 0.7, 0.7};
  float fitMin = xRangeMin;
  float fitMax = 1500.; // set to highest bin with a data point
  float xRangeMinRatio = xRangeMin;
  float xRangeMaxRatio = xRangeMax;

  // output file formats
  const bool savePull = 0;
  const bool saveSpec = 0;
  const bool saveCumSpec = 0;
  const bool saveAsPdf = 0;
  const bool saveAsPng = 1;
  const bool saveAsRoot = 0;
  const char *fileNameExtra = "";
  //const char *fileNameExtra = "madgraphTTbar_";
  const char *plotDir = "../plots/20130523_withSignal/";

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

  TString histoSign[5] = {"", "SS_", "OS_", "e-mu+_", "e+mu-_"};
  TString xAxisTitle[5] = {"m(e#mu)", "m(e^{#pm}#mu^{#pm})", "m(e^{#pm}#mu^{#mp})", "m(e^{-}#mu^{+})", "m(e^{+}#mu^{-})"};
  TString nameSuffix[2] = {"", "cumul"};
  TString titleSuffix[2] = {"", " - Cumulative"};

  vector<TH1F *> emuMass_data;
  vector<TH1F *> emuMass_ttbar;
  //std::vector<TH1F *> emuMass_ttbar700to1000;
  //std::vector<TH1F *> emuMass_ttbar1000up;
  //std::vector<TH1F *> emuMass_ttbarPriv600up;
  vector<TH1F *> emuMass_ztautau;
  vector<TH1F *> emuMass_ww;
  //vector<TH1F *> emuMass_wwPriv600upEminusMuPlus;
  //vector<TH1F *> emuMass_wwPriv600upEplusMuMinus;
  vector<TH1F *> emuMass_wz;
  vector<TH1F *> emuMass_zz;
  vector<TH1F *> emuMass_tw;
  vector<TH1F *> emuMass_zmumu;
  vector<TH1F *> emuMass_zee;
  vector<TH1F *> emuMass_wjets;
  vector<TH1F *> emuMass_qcd;
  vector<TH1F *> emuMass_ttLike;
  vector<TH1F *> emuMass_fakePairs;
  vector<TH1F *> emuMass_sig1;
  vector<TH1F *> emuMass_sig2;
  vector<TH1F *> emuMass_sig3;
  vector<TH1F *> emuMass_sig4;
  vector<TH1F *> dataOverBgHist;

  // define the binning
  vector<float> binning;
  if (varBinning) {
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
  } else {
    for (float bin = 0.; bin <= 1500.; bin += 20.)
      binning.push_back(bin);
  }
  int nBins = binning.size();

  THashList *mcWeights = (THashList *)input.Get("mcWeights");
  THashList *nGenEvents = (THashList *)input.Get("nGenEvents");
  TParameter<float> *ttbarMcWeight = (TParameter<float> *)mcWeights->FindObject("ttbar");
  TParameter<float> *ttbarMcWeight700to1000 = (TParameter<float> *)mcWeights->FindObject("ttbar700to1000");
  TParameter<float> *ttbarMcWeight1000up = (TParameter<float> *)mcWeights->FindObject("ttbar1000up");
  TParameter<float> *ttbarMcWeightTo2l = (TParameter<float> *)mcWeights->FindObject("ttbarto2l");
  TParameter<float> *ttbarMcWeightTo1l1jet = (TParameter<float> *)mcWeights->FindObject("ttbarto1l1jet");
  TParameter<float> *ttbarMcWeight600up = (TParameter<float> *)mcWeights->FindObject("ttbarPriv600up");
  TParameter<float> *wwMcWeight = (TParameter<float> *)mcWeights->FindObject("ww");
  TParameter<float> *wwMcWeight600upEminusMuPlus = (TParameter<float> *)mcWeights->FindObject("wwPriv600upEminusMuPlus");
  TParameter<float> *wwMcWeight600upEplusMuMinus = (TParameter<float> *)mcWeights->FindObject("wwPriv600upEplusMuMinus");
  TParameter<float> *wwNGen600upEminusMuPlus = (TParameter<float> *)nGenEvents->FindObject("wwPriv600upEminusMuPlus");
  TParameter<float> *wwNGen600upEplusMuMinus = (TParameter<float> *)nGenEvents->FindObject("wwPriv600upEplusMuMinus");

  // vectors for processes with multiple samples 
  vector<const char *> cutVarsTtbar;
  vector<float> lowCutsTtbar;
  vector<float> highCutsTtbar;
  vector<float> mcWeightsForCutRangesTtbar;
  // additional ones for madgraph samples
  vector<float> mcWeightsForCutRangesTtbarTo2l;
  vector<float> mcWeightsForCutRangesTtbarTo1l1jet;
  vector<float> lowCutsTtbarPriv600up;
  vector<float> highCutsTtbarPriv600up;
  vector<float> mcWeightsForCutRangesTtbarPriv600up;
  if (ttbar_sample_type) {
    // for madgraph we do not know the number of events generated above emu_mass=600 GeV so we cut there and use the private sample above
    cutVarsTtbar.push_back("emu_mass");
    lowCutsTtbar.push_back(0.);
    highCutsTtbar.push_back(600.);
    highCutsTtbarPriv600up.push_back(1.e9);
    mcWeightsForCutRangesTtbarTo2l.push_back(ttbarMcWeightTo2l->GetVal());
    mcWeightsForCutRangesTtbarTo1l1jet.push_back(ttbarMcWeightTo1l1jet->GetVal());
    mcWeightsForCutRangesTtbarPriv600up.push_back(ttbarMcWeight600up->GetVal());
  } else {
    cutVarsTtbar.push_back("");
    lowCutsTtbar.push_back(0.);
    highCutsTtbar.push_back(1.e9);
    mcWeightsForCutRangesTtbar.push_back(ttbarMcWeight->GetVal());
    cutVarsTtbar.push_back("genMTtbar");
    lowCutsTtbar.push_back(700.);
    highCutsTtbar.push_back(1000.);
    mcWeightsForCutRangesTtbar.push_back(ttbarMcWeight700to1000->GetVal());
    cutVarsTtbar.push_back("genMTtbar");
    lowCutsTtbar.push_back(1000.);
    highCutsTtbar.push_back(1.e9);
    mcWeightsForCutRangesTtbar.push_back(ttbarMcWeight1000up->GetVal());
    cutVarsTtbar.push_back("emu_mass");
    lowCutsTtbar.push_back(600.);
    highCutsTtbar.push_back(1.e9);
    mcWeightsForCutRangesTtbar.push_back(ttbarMcWeight600up->GetVal());
  }
  vector<const char *> cutVarsWw;
  vector<float> lowCutsWw;
  vector<float> highCutsWw;
  vector<float> mcWeightsForCutRangesWw;
  cutVarsWw.push_back("");
  lowCutsWw.push_back(0.);
  highCutsWw.push_back(1.e9);
  mcWeightsForCutRangesWw.push_back(wwMcWeight->GetVal());
  cutVarsWw.push_back("emu_mass");
  lowCutsWw.push_back(600.);
  highCutsWw.push_back(1.e9);
  mcWeightsForCutRangesWw.push_back((wwMcWeight600upEminusMuPlus->GetVal()*wwNGen600upEminusMuPlus->GetVal() + wwMcWeight600upEplusMuMinus->GetVal()*wwNGen600upEplusMuMinus->GetVal())/(wwNGen600upEminusMuPlus->GetVal()+wwNGen600upEplusMuMinus->GetVal()));
  vector<const char *> cutVarsEmpty;
  vector<float> lowCutsEmpty;
  vector<float> highCutsEmpty;
  vector<float> mcWeightsForCutRangesEmpty;

  // determine qcd contribution
  TH1F *qcdContrib;
  if (qcdEst == 1) {
    TH1F *ssData = MakeHistoFromBranch(&input, "emuTree_data", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x100);
    unsigned int flags = 0x1DF; // flags: [apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
    unsigned int topFlags = flags;
    if (topReweighting) topFlags += (1<<10);
    TH1F *ssBg;
    if (ttbar_sample_type) {
      ssBg = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo2l, binning, topFlags);
      ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbarto1l1jet", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo1l1jet, binning, topFlags));
      ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbarPriv600up, mcWeightsForCutRangesTtbarPriv600up, binning, topFlags));
    } else {
      ssBg = MakeHistoFromBranch(&input, "emuTree_ttbar", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags);
      ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
      ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
      ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", "", "mass", SS, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
    }
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ztautau", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ww", "", "mass", SS, eRegion, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEminusMuPlus", "", "mass", SS, eRegion, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEplusMuMinus", "", "mass", SS, eRegion, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wz", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zz", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_tw", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zmumu", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zee", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wjets", "", "mass", SS, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
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

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(font);
  tex->SetLineWidth(2);

  // loop over full spectrum, SS and OS
  for (int k = 0; k < 5; ++k) {
    // loop to get normal and cumulated spectrum
    for (unsigned int j = 0; j < 2; ++j) {
      input.cd();

      bool normToBin = true;
      if (j > 0) normToBin = false;

      if (k == 2) k = -1;
      else if (k == 3) k = -3;
      else if (k == 4) k = -2;
      // make the histograms
      dataOverBgHist.push_back(MakeHistoFromBranch(&input, "emuTree_data", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x100, normToBin));
      emuMass_data.push_back(MakeHistoFromBranch(&input, "emuTree_data", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x100));

      unsigned int flags = 0x1DF; // flags: [apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
      unsigned int topFlags = flags;
      if (topReweighting) topFlags += (1<<10);
      TH1F *ttbarComb;
      if (ttbar_sample_type) {
        ttbarComb = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo2l, binning, topFlags, normToBin);
        ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbarto1l1jet", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo1l1jet, binning, topFlags, normToBin));
        ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbarPriv600up, mcWeightsForCutRangesTtbarPriv600up, binning, topFlags, normToBin));
      } else {
        ttbarComb = MakeHistoFromBranch(&input, "emuTree_ttbar", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags, normToBin);
        ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags, normToBin));
        ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags, normToBin));
        ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", "", "mass", k, eRegion, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags, normToBin));
      }
      emuMass_ttbar.push_back(ttbarComb);
      emuMass_ztautau.push_back(MakeHistoFromBranch(&input, "emuTree_ztautau", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      TH1F *wwComb = MakeHistoFromBranch(&input, "emuTree_ww", "", "mass", k, eRegion, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags, normToBin);
      wwComb->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEminusMuPlus", "", "mass", k, eRegion, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags, normToBin));
      wwComb->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEplusMuMinus", "", "mass", k, eRegion, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags, normToBin));
      emuMass_ww.push_back(wwComb);
      emuMass_wz.push_back(MakeHistoFromBranch(&input, "emuTree_wz", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      emuMass_zz.push_back(MakeHistoFromBranch(&input, "emuTree_zz", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      emuMass_tw.push_back(MakeHistoFromBranch(&input, "emuTree_tw", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      emuMass_zmumu.push_back(MakeHistoFromBranch(&input, "emuTree_zmumu", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      emuMass_zee.push_back(MakeHistoFromBranch(&input, "emuTree_zee", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      if (plotSignal[0] > 0.) emuMass_sig1.push_back(MakeHistoFromBranch(&input, "emuTree_sig500", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      if (plotSignal[1] > 0.) emuMass_sig2.push_back(MakeHistoFromBranch(&input, "emuTree_sig750", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      if (plotSignal[2] > 0.) emuMass_sig3.push_back(MakeHistoFromBranch(&input, "emuTree_sig1000", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      if (plotSignal[3] > 0.) emuMass_sig4.push_back(MakeHistoFromBranch(&input, "emuTree_sig1250", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      if (qcdEst != 2) emuMass_wjets.push_back(MakeHistoFromBranch(&input, "emuTree_wjets", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags, normToBin));
      if (k == -1) k = 2;
      else if (k == -2) k = 4;
      else if (k == -3) k = 3;
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
      if (plotSignal[0] > 0.) {
        emuMass_sig1.back()->SetName("emuMass_" + histoSign[k] + "sig1" + nameSuffix[j]);
        emuMass_sig1.back()->Scale(plotSignal[0]);
      }
      if (plotSignal[1] > 0.) {
        emuMass_sig2.back()->SetName("emuMass_" + histoSign[k] + "sig2" + nameSuffix[j]);
        emuMass_sig2.back()->Scale(plotSignal[1]);
      }
      if (plotSignal[2] > 0.) {
        emuMass_sig3.back()->SetName("emuMass_" + histoSign[k] + "sig3" + nameSuffix[j]);
        emuMass_sig3.back()->Scale(plotSignal[2]);
      }
      if (plotSignal[3] > 0.) {
        emuMass_sig4.back()->SetName("emuMass_" + histoSign[k] + "sig4" + nameSuffix[j]);
        emuMass_sig4.back()->Scale(plotSignal[3]);
      }

      // qcd contribution
      if (qcdEst > 0) {
        if (k == 2) k = -1;
        else if (k == 3) k = -3;
        else if (k == 4) k = -2;
        if (qcdEst == 2) {
          //qcdContrib = MakeHistoFromBranch(&input, "frEmuTree_data", "", "mass", ALL, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x300);
          //if (k != ALL) emuMass_qcd.back()->Scale(0.5);
          qcdContrib = MakeHistoFromBranch(&input, "frEmuTree_data", "", "mass", k, eRegion, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x300);
        }
        if (k == -1) k = 2;
        else if (k == -2) k = 4;
        else if (k == -3) k = 3;
        emuMass_qcd.push_back((TH1F *)qcdContrib->Clone("emuMass_" + histoSign[k] + "qcd"));
        if (qcdEst == 1 && k == ALL) emuMass_qcd.back()->Scale(2.);
        // normalize to bin width
        if (j < 1) {
          for (int i = 1; i < emuMass_qcd.back()->GetNbinsX() + 1; ++i) {
            emuMass_qcd.back()->SetBinContent(i, emuMass_qcd.back()->GetBinContent(i) / emuMass_qcd.back()->GetBinWidth(i));
            emuMass_qcd.back()->SetBinError(i, emuMass_qcd.back()->GetBinError(i) / emuMass_qcd.back()->GetBinWidth(i));
          }
        }
      }

      // add overflow in last bin
      if (j == 0 && overflowBin) {
        dataOverBgHist.back()->SetBinContent(dataOverBgHist.back()->GetNbinsX(), dataOverBgHist.back()->GetBinContent(dataOverBgHist.back()->GetNbinsX()) + dataOverBgHist.back()->GetBinContent(dataOverBgHist.back()->GetNbinsX() + 1));
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
      if (plotSignal[0] > 0.) emuMass_sig1.back()->SetBinContent(emuMass_sig1.back()->GetNbinsX(), emuMass_sig1.back()->GetBinContent(emuMass_sig1.back()->GetNbinsX()) + emuMass_sig1.back()->GetBinContent(emuMass_sig1.back()->GetNbinsX() + 1));
      if (plotSignal[1] > 0.) emuMass_sig2.back()->SetBinContent(emuMass_sig2.back()->GetNbinsX(), emuMass_sig2.back()->GetBinContent(emuMass_sig2.back()->GetNbinsX()) + emuMass_sig2.back()->GetBinContent(emuMass_sig2.back()->GetNbinsX() + 1));
      if (plotSignal[2] > 0.) emuMass_sig3.back()->SetBinContent(emuMass_sig3.back()->GetNbinsX(), emuMass_sig3.back()->GetBinContent(emuMass_sig3.back()->GetNbinsX()) + emuMass_sig3.back()->GetBinContent(emuMass_sig3.back()->GetNbinsX() + 1));
      if (plotSignal[3] > 0.) emuMass_sig4.back()->SetBinContent(emuMass_sig4.back()->GetNbinsX(), emuMass_sig4.back()->GetBinContent(emuMass_sig4.back()->GetNbinsX()) + emuMass_sig4.back()->GetBinContent(emuMass_sig4.back()->GetNbinsX() + 1));
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
          if (plotSignal[0] > 0.) {
            emuMass_sig1.back()->SetBinContent(i, emuMass_sig1.back()->IntegralAndError(i, nBins, error));
            emuMass_sig1.back()->SetBinError(i, error);
          }
          if (plotSignal[1] > 0.) {
            emuMass_sig2.back()->SetBinContent(i, emuMass_sig2.back()->IntegralAndError(i, nBins, error));
            emuMass_sig2.back()->SetBinError(i, error);
          }
          if (plotSignal[2] > 0.) {
            emuMass_sig3.back()->SetBinContent(i, emuMass_sig3.back()->IntegralAndError(i, nBins, error));
            emuMass_sig3.back()->SetBinError(i, error);
          }
          if (plotSignal[3] > 0.) {
            emuMass_sig4.back()->SetBinContent(i, emuMass_sig4.back()->IntegralAndError(i, nBins, error));
            emuMass_sig4.back()->SetBinError(i, error);
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
        emuPlot = new TCanvas("emuPlot" + histoSign[k] + nameSuffix[j], "emu Spectrum" + titleSuffix[j], 100, 100, 900, 720);
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
      
      if (plotSignal[0] > 0.) {
        emuMass_sig1.back()->SetLineColor(kGreen);
        emuMass_sig1.back()->SetLineWidth(2);
      }
      if (plotSignal[1] > 0.) {
        emuMass_sig2.back()->SetLineColor(kBlue);
        emuMass_sig2.back()->SetLineWidth(2);
      }
      if (plotSignal[2] > 0.) {
        emuMass_sig3.back()->SetLineColor(kCyan);
        emuMass_sig3.back()->SetLineWidth(2);
      }
      if (plotSignal[3] > 0.) {
        emuMass_sig4.back()->SetLineColor(kMagenta);
        emuMass_sig4.back()->SetLineWidth(2);
      }

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
        if (plotSignal[0] > 0.) emuMass_sig1.back()->Draw("histsame");
        if (plotSignal[1] > 0.) emuMass_sig2.back()->Draw("histsame");
        if (plotSignal[2] > 0.) emuMass_sig3.back()->Draw("histsame");
        if (plotSignal[3] > 0.) emuMass_sig4.back()->Draw("histsame");

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
        //bgStackGrouped->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax); 
        bgStackGrouped->GetXaxis()->SetLimits(xRangeMin, xRangeMax); 
        if (j == 1) bgStackGrouped->GetYaxis()->SetTitle("Events #geq " + xAxisTitle[k]);
        else bgStackGrouped->GetYaxis()->SetTitle("Events / GeV");
        bgStackGrouped->GetYaxis()->SetTitleFont(font);
        bgStackGrouped->GetYaxis()->SetTitleSize(0.047);
        bgStackGrouped->GetYaxis()->SetTitleOffset(1.1);
        bgStackGrouped->GetYaxis()->SetLabelFont(font);
        bgStackGrouped->GetYaxis()->SetLabelSize(0.05);
        bgStackGrouped->SetMinimum(yRangeMin[k + j * 5]); 
        bgStackGrouped->SetMaximum(yRangeMax[k + j * 5]); 
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
        if (plotSignal[0] > 0.) emuMass_sig1.back()->Draw("histsame");
        if (plotSignal[1] > 0.) emuMass_sig2.back()->Draw("histsame");
        if (plotSignal[2] > 0.) emuMass_sig3.back()->Draw("histsame");
        if (plotSignal[3] > 0.) emuMass_sig4.back()->Draw("histsame");

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
        //bgStack->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax); 
        bgStack->GetXaxis()->SetLimits(xRangeMin, xRangeMax); 
        if (j == 1) bgStack->GetYaxis()->SetTitle("Events #geq " + xAxisTitle[k]);
        else bgStack->GetYaxis()->SetTitle("Events / GeV");
        bgStack->GetYaxis()->SetTitleFont(font);
        bgStack->GetYaxis()->SetTitleSize(0.047);
        bgStack->GetYaxis()->SetTitleOffset(1.1);
        bgStack->GetYaxis()->SetLabelFont(font);
        bgStack->GetYaxis()->SetLabelSize(0.05);
        bgStack->SetMinimum(yRangeMin[k + j * 5]); 
        bgStack->SetMaximum(yRangeMax[k + j * 5]); 
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
        legend.AddEntry(emuMass_ttLike.back(), "t#bar{t} + other prompt leptons" ,"F");
        legend.AddEntry(emuMass_fakePairs.back(), "fake e#mu pairs" ,"F");
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
        if (plotSignal[0] > 0.) legend.AddEntry(emuMass_sig1.back(), Form("Z'/#gamma'500 (x%.0f)", plotSignal[0]) ,"l");
        if (plotSignal[1] > 0.) legend.AddEntry(emuMass_sig2.back(), Form("Z'/#gamma'750 (x%.0f)", plotSignal[1]) ,"l");
        if (plotSignal[2] > 0.) legend.AddEntry(emuMass_sig3.back(), Form("Z'/#gamma'1000 (x%.0f)", plotSignal[2]) ,"l");
        if (plotSignal[3] > 0.) legend.AddEntry(emuMass_sig4.back(), Form("Z'/#gamma'1250 (x%.0f)", plotSignal[3]) ,"l");
      }
      legend.DrawClone("sames");
      
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
        if (!varBinning) sStream << "constBin_";
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
      dataOverBgHist.back()->Add(bgHist, -1.);
      dataOverBgHist.back()->Divide(bgHist);

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

      dataOverBgHist.back()->SetLineWidth(1);
      dataOverBgHist.back()->SetLineColor(kBlack);
      dataOverBgHist.back()->SetMarkerStyle(20);
      dataOverBgHist.back()->SetMarkerSize(1.1);

      dataOverBgHist.back()->GetXaxis()->SetTitle(xAxisTitle[k] + " [GeV]");
      dataOverBgHist.back()->GetXaxis()->SetTitleFont(font);
      dataOverBgHist.back()->GetXaxis()->SetTitleSize(0.047 * fontScaleBot);
      dataOverBgHist.back()->GetXaxis()->SetTitleOffset(0.9);
      dataOverBgHist.back()->GetXaxis()->SetLabelFont(font);
      dataOverBgHist.back()->GetXaxis()->SetLabelSize(0.05 * fontScaleBot);
      dataOverBgHist.back()->GetXaxis()->SetMoreLogLabels();
      dataOverBgHist.back()->GetXaxis()->SetNoExponent();
      dataOverBgHist.back()->GetXaxis()->SetRangeUser(xRangeMinRatio, xRangeMaxRatio - 0.1);
      dataOverBgHist.back()->GetYaxis()->SetTitle("(data-bkg)/bkg");
      dataOverBgHist.back()->GetYaxis()->SetTitleFont(font);
      dataOverBgHist.back()->GetYaxis()->SetTitleSize(0.047 * fontScaleBot);
      dataOverBgHist.back()->GetYaxis()->SetTitleOffset(1.1 / fontScaleBot);
      dataOverBgHist.back()->GetYaxis()->SetLabelFont(font);
      dataOverBgHist.back()->GetYaxis()->SetLabelSize(0.05 * fontScaleBot);
      dataOverBgHist.back()->GetYaxis()->SetRangeUser(yRangeMinRatio[k], yRangeMaxRatio[k]);

      dataOverBgHist.back()->Draw();

      TF1 *f0 = new TF1("f0" + histoSign[k], "[0]");
      dataOverBgHist.back()->Fit("f0" + histoSign[k], "", "", fitMin, fitMax);

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
        if (!varBinning) sStream << "constBin_";
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
        if (!varBinning) sStream << "constBin_";
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

  ////////////////////////////////////////////////////////////////////////////
  // background parametrisation histograms
  gStyle->SetOptFit(1111);
  TCanvas* bgParamPlot;
  TPad* bgParamPad;
  if (plotPull) {
    bgParamPlot = new TCanvas("bgParamPlot", "Bkg Parametrisation", 100, 100, 600, 720);
    bgParamPad = new TPad("bgParamPad", "Bkg Parametrisation", 0., 0.33, 1., 1.);
    bgParamPad->SetBottomMargin(0.06);
  } else {
    bgParamPlot = new TCanvas("bgParamPlot", "Bkg Parametrisation", 100, 100, 600, 600);
    bgParamPad = new TPad("bgParamPad", "Bkg Parametrisation", 0., 0., 1., 1.);
  }
  bgParamPad->SetBorderMode(0);
  bgParamPad->SetBorderSize(2);
  bgParamPad->SetFrameBorderMode(0);
  bgParamPad->SetFillColor(0);
  bgParamPad->SetFrameFillColor(0);
  if (logPlotX) bgParamPad->SetLogx();
  bgParamPad->SetLogy();
  bgParamPad->SetLeftMargin(0.11);
  bgParamPad->SetRightMargin(0.09);
  bgParamPad->SetTopMargin(0.08);
  bgParamPad->SetTicks(1, 1);
  bgParamPad->Draw();
  TPad* bgParamDiffPad = (TPad*)bgParamPad->Clone("bgParamDiffPad");
  if (plotPull) {
    bgParamDiffPad->SetPad(0., 0., 1., 0.33);
    bgParamDiffPad->SetLogy(0);
    bgParamDiffPad->SetGridy();
    bgParamDiffPad->SetTopMargin(0.05);
    bgParamDiffPad->SetBottomMargin(0.22);
    bgParamDiffPad->Draw();
  }
  float fontScaleLow = bgParamPad->GetHNDC() / bgParamDiffPad->GetHNDC();
  bgParamPad->cd();
  TH1F* emuMass_allBkg = (TH1F*)emuMass_ttbar[0]->Clone("");
  if (qcdEst > 0) emuMass_allBkg->Add(emuMass_qcd[0]);
  if (qcdEst != 2) emuMass_allBkg->Add(emuMass_wjets[0]);
  emuMass_allBkg->Add(emuMass_zee[0]);
  emuMass_allBkg->Add(emuMass_zmumu[0]);
  emuMass_allBkg->Add(emuMass_tw[0]);
  emuMass_allBkg->Add(emuMass_zz[0]);
  emuMass_allBkg->Add(emuMass_wz[0]);
  emuMass_allBkg->Add(emuMass_ww[0]);
  emuMass_allBkg->Add(emuMass_ztautau[0]);
  if (!plotPull) emuMass_allBkg->GetXaxis()->SetTitle(xAxisTitle[0] + " [GeV]");
  emuMass_allBkg->GetYaxis()->SetTitle("Events / GeV");
  emuMass_allBkg->GetYaxis()->SetTitleOffset(1.2);
  emuMass_allBkg->Draw();
  TFile* outFile = new TFile("bgFit.root", "recreate");
  //TF1 *bgParamFunc = new TF1("bgParamFunc", "x^(-1*[3]) * exp([0] + [1]*x + [2]*x*x)", 150., 1500.);
  //bgParamFunc->SetParLimits(0, 0., 1000.);
  //bgParamFunc->SetParLimits(1, -10., 0.);
  //bgParamFunc->SetParLimits(2, 0., 1.);
  //bgParamFunc->SetParLimits(3, 0., 5.);
  //bgParamFunc->SetParNames("a", "b", "c", "#kappa");
  TF1 *bgParamFunc = new TF1("bgParamFunc", "1/[1]*(1+([2]*(x-[0]))/([1]))**(-1/[2]-1)", 0., 6000.);
  bgParamFunc->SetParLimits(0, 100., 1.e5);
  bgParamFunc->SetParLimits(1, 10., 1000.);
  bgParamFunc->SetParLimits(2, 0.01, 1.);
  bgParamFunc->SetParNames("m_{min}", "#alpha", "#beta");
  //emuMass_allBkg->Rebin(2);
  emuMass_allBkg->Fit("bgParamFunc", "", "", 150., 1500.);
  bgParamFunc->Write();
  outFile->Close();
  cout << "Chi^2 / NDF: " << bgParamFunc->GetChisquare() << " / " << bgParamFunc->GetNDF() << ", prob: " << bgParamFunc->GetProb() << endl;
  tex->SetTextSize(0.035);
  if (prelim) tex->DrawLatex(0.109, 0.935, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
  else tex->DrawLatex(0.109, 0.935, "CMS, 8 TeV, 19.7 fb^{-1}");
  if (eRegion == 0) tex->DrawLatex(0.325, 0.845, "e in barrel");
  if (eRegion == 1) tex->DrawLatex(0.325, 0.845, "e in endcap");
  //tex->DrawLatex(0.159, 0.180, "P(m|a,b,c,#kappa) = m^{-#kappa} e^{a + b m + c m^{2}}");
  tex->DrawLatex(0.159, 0.180, "P(m|m_{min},#alpha,#beta) = #frac{1}{#alpha}#left(1+#frac{#beta(m-m_{min})}{#alpha}#right)^{-#frac{1}{#beta}-1}");
  if (plotPull) {
    bgParamDiffPad->cd();
    TH1F* bgParamDiff = (TH1F*)emuMass_allBkg->Clone("bgParamDiff");
    bgParamDiff->GetXaxis()->SetTitle(xAxisTitle[0] + " [GeV]");
    bgParamDiff->GetXaxis()->SetTitleSize(emuMass_allBkg->GetXaxis()->GetTitleSize() * fontScaleLow);
    bgParamDiff->GetXaxis()->SetLabelSize(emuMass_allBkg->GetXaxis()->GetTitleSize() * fontScaleLow);
    bgParamDiff->GetYaxis()->SetTitle("(bgk-fit)/fit");
    bgParamDiff->GetYaxis()->SetTitleSize(emuMass_allBkg->GetYaxis()->GetTitleSize() * fontScaleLow);
    bgParamDiff->GetYaxis()->SetTitleOffset(emuMass_allBkg->GetYaxis()->GetTitleOffset() / fontScaleLow);
    bgParamDiff->GetYaxis()->SetLabelSize(emuMass_allBkg->GetYaxis()->GetTitleSize() * fontScaleLow);
    bgParamDiff->Eval(bgParamFunc, "R");
    bgParamDiff->Add(emuMass_allBkg, -1.);
    bgParamDiff->Divide(bgParamFunc, -1.);
    bgParamDiff->Draw();
    bgParamDiff->GetYaxis()->SetRangeUser(-1., 2.5);
  }

  ////////////////////////////////////////////////////////////////////////////
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
  cout << "\\label{tab:emu_event_yield}" << endl;
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
  cout << "\\end{table}" << endl;

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "|Event yield table with combined = SS + OS                                                                |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "\\begin{table}[tbh]" << endl;
  cout << "\\centering" << endl;
  cout << "\\caption{Number of $e\\mu$ events with different charge combinations from data and Monte Carlo simulation. The listed errors are the systematic errors}" << endl;
  cout << "\\label{tab:emu_event_yield}" << endl;
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
  cout << "\\end{table}" << endl;

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "|Event yield table per sample                                                                             |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "\\begin{table}[htb]" << endl;
  cout << "\\centering" << endl;
  cout << "\\caption{Number of $e\\mu$ events with invariant mass in different regions and with different charge combinations.}" << endl;
  cout << "\\label{tab:emu_event_yield_by_source}" << endl;
  cout << "\\resizebox{16cm}{!}{" << endl;
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
  cout << "\\end{tabular}}" << endl;
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

