//===============================================================
// RooFit Macro to perform unbinned likelihood fit 
//================================================================
#ifndef HighMassRes_h
#define HighMassRes_h

#include <sstream>
#include <iostream>
#include <utility>
#include <vector>

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooFormula.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "THashTable.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include "tdrstyle.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooDCBShape.h"
#endif

using namespace RooFit;

class HighMassRes {
public:
  HighMassRes(const char *_inFile = "../forest/emuSpec_19780pb-1.root", const int _lumi = 19780);
  inline virtual ~HighMassRes() {};

  void RunCryBall();

protected:
  const char *inFile;
  int lumi;
  const unsigned int fitModelType; // (0, 1) = (CB, double CB)
  const unsigned int fitModelType2; // (0, 1) = (CB, double CB)
  double cutoff_cb;

  const char *plotOpt;
  const int nBins;
  bool plotReg[7];

  const int font;
  const bool plotPull;

  const bool saveFitsAsPdf;
  const bool saveFitsAsPng;
  const bool saveFitsAsRoot;
  const bool saveResAsPdf;
  const bool saveResAsPng;
  const bool saveResAsRoot;
  const char *fileNameExtra;
  const char *plotDir;

  const bool useRootTermForFit;

  std::vector<std::pair<Float_t, Float_t> > sigmaExtras;

  float fitRangeMin;
  float fitRangeMax;
  std::vector<std::pair<Int_t, const char *> > zPrimeGenMasses;

  const int fitColor;
  const int fitColorRes;

  std::vector<TString> regTxt;
  std::vector<TString> regFileNameSuffix;

  TString biasName;
  TString sigmaName;
  TString cutLName;
  TString cutRName;
  TString powerLName;
  TString powerRName;

  // histos for the Z'PSI signal samples
  std::vector<TH1F *> sigmaHistosZpPsi;
  //std::vector<TH1F *> dmHistosZpPsi;
  //std::vector<TH1F *> acbHistosZpPsi;
  //std::vector<TH1F *> ncbHistosZpPsi;

  RooWorkspace * cryBall(
                   TTree* tree,
                   const int regions = 3,
                   double minMass = 500,
                   double maxMass = 600,
                   double trueMassMin = 500,
                   double trueMassMax = 600,
                   double cutoff_cb = 1.0,
                   const char* plotOpt = "NEU",
                   const int nBins = 80,
                   const unsigned int fitModelType = 1); // (0, 1) = (CB, double CB)
};

#endif

HighMassRes::HighMassRes(const char *_inFile, const int _lumi) :
                   fitModelType(0), // (0, 1) = (CB, double CB)
                   fitModelType2(0), // (0, 1) = (CB, double CB)
                   cutoff_cb(2.),
                   plotOpt("NEU"),
                   nBins(80),
                   font(42),
                   plotPull(true),
                   saveFitsAsPdf(0),
                   saveFitsAsPng(0),
                   saveFitsAsRoot(0),
                   saveResAsPdf(0),
                   saveResAsPng(0),
                   saveResAsRoot(0),
                   fileNameExtra(""),
                   plotDir("./"),
                   useRootTermForFit(0),
                   fitColor(kRed),
                   fitColorRes(kBlue),
                   biasName("#Deltam_{CB}"),
                   sigmaName("#sigma_{CB}"),
                   cutLName("a_{CB}"),
                   cutRName("a_{CB}"),
                   powerLName("n_{CB}"),
                   powerRName("n_{CB}")
{
  inFile = _inFile;
  lumi = _lumi;

  fitRangeMin = 0.;
  fitRangeMax = 5000.;

  plotReg[0] = 0; // EB-EB
  plotReg[1] = 0; // EB-EE
  plotReg[2] = 1; // EB-EB + EB-EE

  // powheg START53 19.6fb-1 HEEP v4.1 double Crystal Ball
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));

  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (500, "emuTree_sigNoAccCuts500"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (750, "emuTree_sigNoAccCuts750"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1000, "emuTree_sigNoAccCuts1000"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1250, "emuTree_sigNoAccCuts1250"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1500, "emuTree_sigNoAccCuts1500"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1750, "emuTree_sigNoAccCuts1750"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (2000, "emuTree_sigNoAccCuts2000"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (2500, "emuTree_sigNoAccCuts2500"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (3000, "emuTree_sigNoAccCuts3000"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (3500, "emuTree_sigNoAccCuts3500"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (4000, "emuTree_sigNoAccCuts4000"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (5000, "emuTree_sigNoAccCuts5000"));

  regTxt.push_back("EB-EB");
  regTxt.push_back("EB-EE");
  regTxt.push_back("EB-EB + EB-EE");

  regFileNameSuffix.push_back("BB");
  regFileNameSuffix.push_back("BE");
  regFileNameSuffix.push_back("BBBE");

  if (fitModelType > 0) {
    biasName = "#Deltam_{DCB}";
    sigmaName = "#sigma_{DCB}";
    cutLName = "al_{DCB}";
    cutRName = "ar_{DCB}";
    powerLName = "nl_{DCB}";
    powerRName = "nr_{DCB}";
  }

  // histos for the Z'PSI signal samples
  sigmaHistosZpPsi.reserve(4);
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBB", "sigma of fitted function for Signal MC EB-EB", 338, 0., 5010.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBE", "sigma of fitted function for Signal MC EB-EE", 338, 0., 5010.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 0., 5010.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiEE", "sigma of fitted function for Signal MC EE-EE", 338, 0., 5010.));

  //dmHistosZpPsi.reserve(4);
  //dmHistosZpPsi.push_back(new TH1F("dmZpPsiBB", "dm of fitted function EB-EB", nBinsHisto, binArray));
  //dmHistosZpPsi.push_back(new TH1F("dmZpPsiBE", "dm of fitted function EB-EE", nBinsHisto, binArray));
  //dmHistosZpPsi.push_back(new TH1F("dmZpPsiBBBE", "dm of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  //dmHistosZpPsi.push_back(new TH1F("dmZpPsiEE", "dm of fitted function EE-EE", nBinsHisto, binArray));

  //acbHistosZpPsi.reserve(4);
  //acbHistosZpPsi.push_back(new TH1F("acbZpPsiBB", "acb of fitted function EB-EB", nBinsHisto, binArray));
  //acbHistosZpPsi.push_back(new TH1F("acbZpPsiBE", "acb of fitted function EB-EE", nBinsHisto, binArray));
  //acbHistosZpPsi.push_back(new TH1F("acbZpPsiBBBE", "acb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  //acbHistosZpPsi.push_back(new TH1F("acbZpPsiEE", "acb of fitted function EE-EE", nBinsHisto, binArray));

  //ncbHistosZpPsi.reserve(4);
  //ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBB", "ncb of fitted function EB-EB", nBinsHisto, binArray));
  //ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBE", "ncb of fitted function EB-EE", nBinsHisto, binArray));
  //ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBBBE", "ncb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  //ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiEE", "ncb of fitted function EE-EE", nBinsHisto, binArray));

  TH1::SetDefaultSumw2(kTRUE);
}


