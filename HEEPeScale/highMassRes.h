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
  HighMassRes(const char *_inFile = "plots_28dec2012/eScaleEvents19616pb-1.root", const int _lumi = 19616);
  inline virtual ~HighMassRes() {};

  void RunCryBall();
  void CompareCryBall();

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

  std::vector<std::pair<Float_t, const char *> > dyFitRanges;
  std::vector<std::pair<Int_t, const char *> > zPrimeGenMasses;
  unsigned int ssmStart;

  const int fitColorDy;
  const int fitColorZpPsi;
  const int fitColorZpSsm;
  const int fitColorRes;
  const int fit2ColorDy;
  const int fit2ColorZpPsi;
  const int fit2ColorZpSsm;
  const int fitColorRes2;

  std::vector<TString> regTxt;
  std::vector<TString> regFileNameSuffix;

  TString biasName;
  TString sigmaName;
  TString cutLName;
  TString cutRName;
  TString powerLName;
  TString powerRName;

  // histos for the DY samples
  std::vector<TH1F *> sigmaHistos;
  std::vector<TH1F *> dmHistos;
  std::vector<TH1F *> acbHistos;
  std::vector<TH1F *> ncbHistos;
  // histos for the Z'PSI signal samples
  std::vector<TH1F *> sigmaHistosZpPsi;
  std::vector<TH1F *> dmHistosZpPsi;
  std::vector<TH1F *> acbHistosZpPsi;
  std::vector<TH1F *> ncbHistosZpPsi;
  // histos for the Z'SSM signal samples
  std::vector<TH1F *> sigmaHistosZpSsm;
  std::vector<TH1F *> dmHistosZpSsm;
  std::vector<TH1F *> acbHistosZpSsm;
  std::vector<TH1F *> ncbHistosZpSsm;

  std::vector<TH1F *> sigmaDCBHistos;
  std::vector<TH1F *> sigmaDCBHistosZpPsi;
  std::vector<TH1F *> sigmaDCBHistosZpSsm;

  RooWorkspace * cryBall(
                   const char *treeName = "eleDY500Tree",
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
                   fitModelType(1), // (0, 1) = (CB, double CB)
                   fitModelType2(1), // (0, 1) = (CB, double CB)
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
                   plotDir("./plots_28dec2012/"),
                   useRootTermForFit(0),
                   fitColorDy(kRed),
                   fitColorZpPsi(kGreen+1),
                   fitColorZpSsm(kMagenta+1),
                   fitColorRes(kBlue),
                   fit2ColorDy(kCyan+1),
                   fit2ColorZpPsi(kOrange+7),
                   fit2ColorZpSsm(kOrange),
                   fitColorRes2(kOrange-3),
                   biasName("#Deltam_{CB}"),
                   sigmaName("#sigma_{CB}"),
                   cutLName("a_{CB}"),
                   cutRName("a_{CB}"),
                   powerLName("n_{CB}"),
                   powerRName("n_{CB}")
{
  inFile = _inFile;
  lumi = _lumi;

  plotReg[0] = 1; // EB-EB
  plotReg[1] = 1; // EB-EE
  plotReg[2] = 1; // EB-EB + EB-EE
  plotReg[3] = 0; // EE-EE
  plotReg[4] = 0; // EB-EB + EE-EE
  plotReg[5] = 0; // EB-EE + EE-EE
  plotReg[6] = 0; // EB-EB + EB-EE + EE-EE

  // powheg START53 14.8fb-1 HEEP v4.1 Crystal Ball
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.82, 0.01));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.19, 0.02));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.93, 0.01));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.68, 0.03));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  // powheg START53 14.8fb-1 HEEP v4.1 double Crystal Ball
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.89, 0.02));
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.34, 0.02));
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.18, 0.02));
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.92, 0.04));
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  //sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  // powheg START53 19.6fb-1 HEEP v4.1 double Crystal Ball
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.92, 0.02));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.43, 0.02));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.25, 0.02));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.99, 0.04));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));

  dyFitRanges.push_back(std::pair<Float_t, const char *> (100., "eleDY20Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (120., "eleDY120Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (200., "eleDY200Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (300., "eleDY200Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (400., "eleDY400Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (500., "eleDY500Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (600., "eleDY500Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (700., "eleDY700Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (800., "eleDY800Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (1000., "eleDY1000Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (1250., "eleDY1000Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (1500., "eleDY1500Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (1750., "eleDY1500Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (2000., "eleDY2000Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (2250., "eleDY2000Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (2500., "eleDY2000Tree"));
  dyFitRanges.push_back(std::pair<Float_t, const char *> (3500., "eleDY2000Tree"));

  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (750, "eleZp750Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1000, "eleZp1000Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1250, "eleZp1250Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1500, "eleZp1500Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (1750, "eleZp1750Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (2000, "eleZp2000Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (2250, "eleZp2250Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (3000, "eleZp3000Tree"));
  zPrimeGenMasses.push_back(std::pair<Int_t, const char *> (2250, "eleZpSsm2250Tree"));
  ssmStart = 8;

  regTxt.push_back("EB-EB");
  regTxt.push_back("EB-EE");
  regTxt.push_back("EB-EB + EB-EE");
  regTxt.push_back("EE-EE");
  regTxt.push_back("EB-EB + EE-EE");
  regTxt.push_back("EB-EE + EE-EE");
  regTxt.push_back("EB-EB + EB-EE + EE-EE");

  regFileNameSuffix.push_back("BB");
  regFileNameSuffix.push_back("BE");
  regFileNameSuffix.push_back("BBBE");
  regFileNameSuffix.push_back("EE");
  regFileNameSuffix.push_back("BBEE");
  regFileNameSuffix.push_back("BEEE");
  regFileNameSuffix.push_back("BBBEEE");

  if (fitModelType > 0) {
    biasName = "#Deltam_{DCB}";
    sigmaName = "#sigma_{DCB}";
    cutLName = "al_{DCB}";
    cutRName = "ar_{DCB}";
    powerLName = "nl_{DCB}";
    powerRName = "nr_{DCB}";
  }

  // set up histograms for the high mass resolution parametrisation
  Int_t nBinsHisto = dyFitRanges.size() - 1;
  Float_t binArray[16];
  for (Int_t it = 0; it <= nBinsHisto; ++it)
    binArray[it] = dyFitRanges[it].first;
  // histos for the DY samples
  sigmaHistos.push_back(new TH1F("sigmaBB", "sigma of fitted function EB-EB", nBinsHisto, binArray));
  sigmaHistos.push_back(new TH1F("sigmaBE", "sigma of fitted function EB-EE", nBinsHisto, binArray));
  sigmaHistos.push_back(new TH1F("sigmaBBBE", "sigma of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  sigmaHistos.push_back(new TH1F("sigmaEE", "sigma of fitted function EE-EE", nBinsHisto, binArray));

  dmHistos.push_back(new TH1F("dmBB", "dm of fitted function EB-EB", nBinsHisto, binArray));
  dmHistos.push_back(new TH1F("dmBE", "dm of fitted function EB-EE", nBinsHisto, binArray));
  dmHistos.push_back(new TH1F("dmBBBE", "dm of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  dmHistos.push_back(new TH1F("dmEE", "dm of fitted function EE-EE", nBinsHisto, binArray));

  acbHistos.push_back(new TH1F("acbBB", "acb of fitted function EB-EB", nBinsHisto, binArray));
  acbHistos.push_back(new TH1F("acbBE", "acb of fitted function EB-EE", nBinsHisto, binArray));
  acbHistos.push_back(new TH1F("acbBBBE", "acb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  acbHistos.push_back(new TH1F("acbEE", "acb of fitted function EE-EE", nBinsHisto, binArray));
  
  ncbHistos.push_back(new TH1F("ncbBB", "ncb of fitted function EB-EB", nBinsHisto, binArray));
  ncbHistos.push_back(new TH1F("ncbBE", "ncb of fitted function EB-EE", nBinsHisto, binArray));
  ncbHistos.push_back(new TH1F("ncbBBBE", "ncb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  ncbHistos.push_back(new TH1F("ncbEE", "ncb of fitted function EE-EE", nBinsHisto, binArray));
  // histos for the Z'PSI signal samples
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));

  dmHistosZpPsi.push_back(new TH1F("dmZpPsiBB", "dm of fitted function EB-EB", nBinsHisto, binArray));
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiBE", "dm of fitted function EB-EE", nBinsHisto, binArray));
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiBBBE", "dm of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiEE", "dm of fitted function EE-EE", nBinsHisto, binArray));

  acbHistosZpPsi.push_back(new TH1F("acbZpPsiBB", "acb of fitted function EB-EB", nBinsHisto, binArray));
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiBE", "acb of fitted function EB-EE", nBinsHisto, binArray));
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiBBBE", "acb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiEE", "acb of fitted function EE-EE", nBinsHisto, binArray));

  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBB", "ncb of fitted function EB-EB", nBinsHisto, binArray));
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBE", "ncb of fitted function EB-EE", nBinsHisto, binArray));
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBBBE", "ncb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiEE", "ncb of fitted function EE-EE", nBinsHisto, binArray));
  // histos for the Z'SSM signal samples
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));

  dmHistosZpSsm.push_back(new TH1F("dmZpSsmBB", "dm of fitted function EB-EB", nBinsHisto, binArray));
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmBE", "dm of fitted function EB-EE", nBinsHisto, binArray));
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmBBBE", "dm of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmEE", "dm of fitted function EE-EE", nBinsHisto, binArray));

  acbHistosZpSsm.push_back(new TH1F("acbZpSsmBB", "acb of fitted function EB-EB", nBinsHisto, binArray));
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmBE", "acb of fitted function EB-EE", nBinsHisto, binArray));
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmBBBE", "acb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmEE", "acb of fitted function EE-EE", nBinsHisto, binArray));

  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmBB", "ncb of fitted function EB-EB", nBinsHisto, binArray));
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmBE", "ncb of fitted function EB-EE", nBinsHisto, binArray));
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmBBBE", "ncb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmEE", "ncb of fitted function EE-EE", nBinsHisto, binArray));

  ///////
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBBB", "sigma of fitted function EB-EB", nBinsHisto, binArray));
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBBE", "sigma of fitted function EB-EE", nBinsHisto, binArray));
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBBBBE", "sigma of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBEE", "sigma of fitted function EE-EE", nBinsHisto, binArray));

  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));

  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));


  TH1::SetDefaultSumw2(kTRUE);
}


