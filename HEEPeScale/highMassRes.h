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
  HighMassRes(const char *_inFile = "./eScaleEvents19712pb-1.root", const int _lumi = 19712);
  inline virtual ~HighMassRes() {};

  void RunCryBall();
  void PlotRes();
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
  bool mcOnly;
  int eleEbReg;

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

  std::vector<const char *> dyTreeNames;
  std::vector<Float_t> dyFitRanges;
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
  std::vector<TH1F *> dmHistosRel;
  std::vector<TH1F *> acbHistos;
  std::vector<TH1F *> ncbHistos;
  std::vector<TH1F *> ardcbHistos;
  std::vector<TH1F *> nrdcbHistos;
  // histos for the Z'PSI signal samples
  std::vector<TH1F *> sigmaHistosZpPsi;
  std::vector<TH1F *> dmHistosZpPsi;
  std::vector<TH1F *> dmHistosRelZpPsi;
  std::vector<TH1F *> acbHistosZpPsi;
  std::vector<TH1F *> ncbHistosZpPsi;
  std::vector<TH1F *> ardcbHistosZpPsi;
  std::vector<TH1F *> nrdcbHistosZpPsi;
  // histos for the Z'SSM signal samples
  std::vector<TH1F *> sigmaHistosZpSsm;
  std::vector<TH1F *> dmHistosZpSsm;
  std::vector<TH1F *> dmHistosRelZpSsm;
  std::vector<TH1F *> acbHistosZpSsm;
  std::vector<TH1F *> ncbHistosZpSsm;
  std::vector<TH1F *> ardcbHistosZpSsm;
  std::vector<TH1F *> nrdcbHistosZpSsm;

  std::vector<TH1F *> sigmaDCBHistos;
  std::vector<TH1F *> sigmaDCBHistosZpPsi;
  std::vector<TH1F *> sigmaDCBHistosZpSsm;

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
                   fitModelType(1), // (0, 1) = (CB, double CB)
                   fitModelType2(1), // (0, 1) = (CB, double CB)
                   cutoff_cb(2.),
                   plotOpt("NEU"),
                   nBins(80),
                   font(42),
                   plotPull(true),
                   saveFitsAsPdf(1),
                   saveFitsAsPng(1),
                   saveFitsAsRoot(1),
                   saveResAsPdf(1),
                   saveResAsPng(1),
                   saveResAsRoot(1),
                   fileNameExtra(""),
                   plotDir("./plots_20141014_res_splitRangeFit/"),
                   //plotDir("./plots_20141014_bb_noSigmaExtra/"),
                   //plotDir("./plots_20141014_bb_eleEtaSmaller0p7/"),
                   //plotDir("./plots_20141014_bb_eleEtaLarger0p7/"),
                   //plotDir("./plots_20141014_bb_oneEleEtaSmaller0p7/"),
                   //plotDir("./plots_20141014_bb_oneEleEtaLarger0p7/"),
                   //plotDir("./plots_20141014_bb_oneEleEtaLargerOneSmaller0p7/"),
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
                   cutLName("#alpha_{CB}"),
                   cutRName("#alpha_{CB}"),
                   powerLName("n_{CB}"),
                   powerRName("n_{CB}")
{
  inFile = _inFile;
  lumi = _lumi;

  mcOnly = 0; // no use of sigma_extra
  eleEbReg = 0; // electron eta regions. 0: full, 1: 2e<0.7, 2: 2e>=0.7, 3: 1e<0.7, 4: 1e>=0.7, 5: 1e<0.7 && 1e>=0.7

  plotReg[0] = 1; // EB-EB
  plotReg[1] = 1; // EB-EE
  plotReg[2] = 0; // EB-EB + EB-EE
  plotReg[3] = 1; // EE-EE
  plotReg[4] = 0; // EB-EB + EE-EE
  plotReg[5] = 0; // EB-EE + EE-EE
  plotReg[6] = 0; // EB-EB + EB-EE + EE-EE

  // powheg START53 19.7fb-1 HEEP v4.1 double Crystal Ball
  sigmaExtras.reserve(7);
  if (mcOnly) {
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  } else {
    if (eleEbReg == 1) sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.70, 0.01)); // for |eta_e1| < 0.7 && |eta_e2| < 0.7
    else if (eleEbReg == 2) sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.96, 0.01)); // for |eta_e1| >= 0.7 && |eta_e2| >= 0.7
    else if (eleEbReg == 3) sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.69, 0.01)); // for |eta_e1| < 0.7 || |eta_e2| < 0.7
    else if (eleEbReg == 4) sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.84, 0.01)); // for |eta_e1| >= 0.7 || |eta_e2| >= 0.7
    else if (eleEbReg == 5) sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.84, 0.01)); // for (|eta_e1| < 0.7 && |eta_e2| >= 0.7) || (|eta_e1| >= 0.7 && |eta_e2| < 0.7)
    else sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.73, 0.01));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0.93, 0.01));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (1.2, 0.01));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
    sigmaExtras.push_back(std::pair<Float_t, Float_t> (0., 0.));
  }

  dyFitRanges.reserve(16);
  dyFitRanges.push_back(100.);
  dyFitRanges.push_back(120.);
  dyFitRanges.push_back(200.);
  dyFitRanges.push_back(300.);
  dyFitRanges.push_back(400.);
  dyFitRanges.push_back(500.);
  dyFitRanges.push_back(600.);
  dyFitRanges.push_back(700.);
  dyFitRanges.push_back(800.);
  dyFitRanges.push_back(1000.);
  dyFitRanges.push_back(1250.);
  dyFitRanges.push_back(1500.);
  dyFitRanges.push_back(1750.);
  dyFitRanges.push_back(2000.);
  dyFitRanges.push_back(2250.);
  dyFitRanges.push_back(2500.);
  dyFitRanges.push_back(3500.);

  zPrimeGenMasses.reserve(9);
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

  dyTreeNames.reserve(10);
  dyTreeNames.push_back("eleDY20Tree");
  dyTreeNames.push_back("eleDY120Tree");
  dyTreeNames.push_back("eleDY200Tree");
  dyTreeNames.push_back("eleDY400Tree");
  dyTreeNames.push_back("eleDY500Tree");
  dyTreeNames.push_back("eleDY700Tree");
  dyTreeNames.push_back("eleDY800Tree");
  dyTreeNames.push_back("eleDY1000Tree");
  dyTreeNames.push_back("eleDY1500Tree");
  dyTreeNames.push_back("eleDY2000Tree");

  regTxt.reserve(7);
  regTxt.push_back("EB-EB");
  regTxt.push_back("EB-EE");
  regTxt.push_back("EB-EB + EB-EE");
  regTxt.push_back("EE-EE");
  regTxt.push_back("EB-EB + EE-EE");
  regTxt.push_back("EB-EE + EE-EE");
  regTxt.push_back("EB-EB + EB-EE + EE-EE");

  regFileNameSuffix.reserve(7);
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
    cutLName = "#alpha_{L}^{DCB}";
    cutRName = "#alpha_{R}^{DCB}";
    powerLName = "n_{L}^{DCB}";
    powerRName = "n_{R}^{DCB}";
  }

  // set up histograms for the high mass resolution parametrisation
  Int_t nBinsHisto = dyFitRanges.size() - 1;
  Float_t binArray[16];
  for (Int_t it = 0; it <= nBinsHisto; ++it)
    binArray[it] = dyFitRanges[it];
  // histos for the DY samples
  sigmaHistos.reserve(4);
  sigmaHistos.push_back(new TH1F("sigmaBB", "sigma of fitted function EB-EB", nBinsHisto, binArray));
  sigmaHistos.push_back(new TH1F("sigmaBE", "sigma of fitted function EB-EE", nBinsHisto, binArray));
  sigmaHistos.push_back(new TH1F("sigmaBBBE", "sigma of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  sigmaHistos.push_back(new TH1F("sigmaEE", "sigma of fitted function EE-EE", nBinsHisto, binArray));

  dmHistos.reserve(4);
  dmHistos.push_back(new TH1F("dmBB", "dm of fitted function EB-EB", nBinsHisto, binArray));
  dmHistos.push_back(new TH1F("dmBE", "dm of fitted function EB-EE", nBinsHisto, binArray));
  dmHistos.push_back(new TH1F("dmBBBE", "dm of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  dmHistos.push_back(new TH1F("dmEE", "dm of fitted function EE-EE", nBinsHisto, binArray));

  dmHistosRel.reserve(4);
  dmHistosRel.push_back(new TH1F("dmRelBB", "relative dm of fitted function EB-EB", nBinsHisto, binArray));
  dmHistosRel.push_back(new TH1F("dmRelBE", "relative dm of fitted function EB-EE", nBinsHisto, binArray));
  dmHistosRel.push_back(new TH1F("dmRelBBBE", "relative dm of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  dmHistosRel.push_back(new TH1F("dmRelEE", "relative dm of fitted function EE-EE", nBinsHisto, binArray));

  acbHistos.reserve(4);
  acbHistos.push_back(new TH1F("acbBB", "acb of fitted function EB-EB", nBinsHisto, binArray));
  acbHistos.push_back(new TH1F("acbBE", "acb of fitted function EB-EE", nBinsHisto, binArray));
  acbHistos.push_back(new TH1F("acbBBBE", "acb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  acbHistos.push_back(new TH1F("acbEE", "acb of fitted function EE-EE", nBinsHisto, binArray));
  
  ncbHistos.reserve(4);
  ncbHistos.push_back(new TH1F("ncbBB", "ncb of fitted function EB-EB", nBinsHisto, binArray));
  ncbHistos.push_back(new TH1F("ncbBE", "ncb of fitted function EB-EE", nBinsHisto, binArray));
  ncbHistos.push_back(new TH1F("ncbBBBE", "ncb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  ncbHistos.push_back(new TH1F("ncbEE", "ncb of fitted function EE-EE", nBinsHisto, binArray));

  ardcbHistos.reserve(4);
  ardcbHistos.push_back(new TH1F("ardcbBB", "ardcb of fitted function EB-EB", nBinsHisto, binArray));
  ardcbHistos.push_back(new TH1F("ardcbBE", "ardcb of fitted function EB-EE", nBinsHisto, binArray));
  ardcbHistos.push_back(new TH1F("ardcbBBBE", "ardcb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  ardcbHistos.push_back(new TH1F("ardcbEE", "ardcb of fitted function EE-EE", nBinsHisto, binArray));
  
  nrdcbHistos.reserve(4);
  nrdcbHistos.push_back(new TH1F("nrdcbBB", "nrdcb of fitted function EB-EB", nBinsHisto, binArray));
  nrdcbHistos.push_back(new TH1F("nrdcbBE", "nrdcb of fitted function EB-EE", nBinsHisto, binArray));
  nrdcbHistos.push_back(new TH1F("nrdcbBBBE", "nrdcb of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  nrdcbHistos.push_back(new TH1F("nrdcbEE", "nrdcb of fitted function EE-EE", nBinsHisto, binArray));
  // histos for the Z'PSI signal samples
  sigmaHistosZpPsi.reserve(4);
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaHistosZpPsi.push_back(new TH1F("sigmaZpPsiEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));

  dmHistosZpPsi.reserve(4);
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiBB", "dm of fitted function EB-EB", 338, 120., 3500.));
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiBE", "dm of fitted function EB-EE", 338, 120., 3500.));
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiBBBE", "dm of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  dmHistosZpPsi.push_back(new TH1F("dmZpPsiEE", "dm of fitted function EE-EE", 338, 120., 3500.));

  dmHistosRelZpPsi.reserve(4);
  dmHistosRelZpPsi.push_back(new TH1F("dmRelZpPsiBB", "relative dm of fitted function EB-EB", 338, 120., 3500.));
  dmHistosRelZpPsi.push_back(new TH1F("dmRelZpPsiBE", "relative dm of fitted function EB-EE", 338, 120., 3500.));
  dmHistosRelZpPsi.push_back(new TH1F("dmRelZpPsiBBBE", "relative dm of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  dmHistosRelZpPsi.push_back(new TH1F("dmRelZpPsiEE", "relative dm of fitted function EE-EE", 338, 120., 3500.));

  acbHistosZpPsi.reserve(4);
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiBB", "acb of fitted function EB-EB", 338, 120., 3500.));
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiBE", "acb of fitted function EB-EE", 338, 120., 3500.));
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiBBBE", "acb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  acbHistosZpPsi.push_back(new TH1F("acbZpPsiEE", "acb of fitted function EE-EE", 338, 120., 3500.));

  ncbHistosZpPsi.reserve(4);
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBB", "ncb of fitted function EB-EB", 338, 120., 3500.));
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBE", "ncb of fitted function EB-EE", 338, 120., 3500.));
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiBBBE", "ncb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  ncbHistosZpPsi.push_back(new TH1F("ncbZpPsiEE", "ncb of fitted function EE-EE", 338, 120., 3500.));

  ardcbHistosZpPsi.reserve(4);
  ardcbHistosZpPsi.push_back(new TH1F("ardcbZpPsiBB", "ardcb of fitted function EB-EB", 338, 120., 3500.));
  ardcbHistosZpPsi.push_back(new TH1F("ardcbZpPsiBE", "ardcb of fitted function EB-EE", 338, 120., 3500.));
  ardcbHistosZpPsi.push_back(new TH1F("ardcbZpPsiBBBE", "ardcb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  ardcbHistosZpPsi.push_back(new TH1F("ardcbZpPsiEE", "ardcb of fitted function EE-EE", 338, 120., 3500.));

  nrdcbHistosZpPsi.reserve(4);
  nrdcbHistosZpPsi.push_back(new TH1F("nrdcbZpPsiBB", "nrdcb of fitted function EB-EB", 338, 120., 3500.));
  nrdcbHistosZpPsi.push_back(new TH1F("nrdcbZpPsiBE", "nrdcb of fitted function EB-EE", 338, 120., 3500.));
  nrdcbHistosZpPsi.push_back(new TH1F("nrdcbZpPsiBBBE", "nrdcb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  nrdcbHistosZpPsi.push_back(new TH1F("nrdcbZpPsiEE", "nrdcb of fitted function EE-EE", 338, 120., 3500.));
  // histos for the Z'SSM signal samples
  sigmaHistosZpSsm.reserve(4);
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaHistosZpSsm.push_back(new TH1F("sigmaZpSsmEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));

  dmHistosZpSsm.reserve(4);
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmBB", "dm of fitted function EB-EB", 338, 120., 3500.));
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmBE", "dm of fitted function EB-EE", 338, 120., 3500.));
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmBBBE", "dm of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  dmHistosZpSsm.push_back(new TH1F("dmZpSsmEE", "dm of fitted function EE-EE", 338, 120., 3500.));

  dmHistosRelZpSsm.reserve(4);
  dmHistosRelZpSsm.push_back(new TH1F("dmRelZpSsmBB", "relative dm of fitted function EB-EB", 338, 120., 3500.));
  dmHistosRelZpSsm.push_back(new TH1F("dmRelZpSsmBE", "relative dm of fitted function EB-EE", 338, 120., 3500.));
  dmHistosRelZpSsm.push_back(new TH1F("dmRelZpSsmBBBE", "relative dm of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  dmHistosRelZpSsm.push_back(new TH1F("dmRelZpSsmEE", "relative dm of fitted function EE-EE", 338, 120., 3500.));

  acbHistosZpSsm.reserve(4);
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmBB", "acb of fitted function EB-EB", 338, 120., 3500.));
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmBE", "acb of fitted function EB-EE", 338, 120., 3500.));
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmBBBE", "acb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  acbHistosZpSsm.push_back(new TH1F("acbZpSsmEE", "acb of fitted function EE-EE", 338, 120., 3500.));

  ncbHistosZpSsm.reserve(4);
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmBB", "ncb of fitted function EB-EB", 338, 120., 3500.));
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmBE", "ncb of fitted function EB-EE", 338, 120., 3500.));
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmBBBE", "ncb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  ncbHistosZpSsm.push_back(new TH1F("ncbZpSsmEE", "ncb of fitted function EE-EE", 338, 120., 3500.));

  ardcbHistosZpSsm.reserve(4);
  ardcbHistosZpSsm.push_back(new TH1F("ardcbZpSsmBB", "ardcb of fitted function EB-EB", 338, 120., 3500.));
  ardcbHistosZpSsm.push_back(new TH1F("ardcbZpSsmBE", "ardcb of fitted function EB-EE", 338, 120., 3500.));
  ardcbHistosZpSsm.push_back(new TH1F("ardcbZpSsmBBBE", "ardcb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  ardcbHistosZpSsm.push_back(new TH1F("ardcbZpSsmEE", "ardcb of fitted function EE-EE", 338, 120., 3500.));

  nrdcbHistosZpSsm.reserve(4);
  nrdcbHistosZpSsm.push_back(new TH1F("nrdcbZpSsmBB", "nrdcb of fitted function EB-EB", 338, 120., 3500.));
  nrdcbHistosZpSsm.push_back(new TH1F("nrdcbZpSsmBE", "nrdcb of fitted function EB-EE", 338, 120., 3500.));
  nrdcbHistosZpSsm.push_back(new TH1F("nrdcbZpSsmBBBE", "nrdcb of fitted function EB-EB + EB-EE", 338, 120., 3500.));
  nrdcbHistosZpSsm.push_back(new TH1F("nrdcbZpSsmEE", "nrdcb of fitted function EE-EE", 338, 120., 3500.));

  ///////
  sigmaDCBHistos.reserve(4);
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBBB", "sigma of fitted function EB-EB", nBinsHisto, binArray));
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBBE", "sigma of fitted function EB-EE", nBinsHisto, binArray));
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBBBBE", "sigma of fitted function EB-EB + EB-EE", nBinsHisto, binArray));
  sigmaDCBHistos.push_back(new TH1F("sigmaDCBEE", "sigma of fitted function EE-EE", nBinsHisto, binArray));

  sigmaDCBHistosZpPsi.reserve(4);
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpPsi.push_back(new TH1F("sigmaDCBZpPsiEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));

  sigmaDCBHistosZpSsm.reserve(4);
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmBB", "sigma of fitted function for Signal MC EB-EB", 338, 120., 3500.));
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmBE", "sigma of fitted function for Signal MC EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmBBBE", "sigma of fitted function for Signal MC EB-EB + EB-EE", 338, 120., 3500.));
  sigmaDCBHistosZpSsm.push_back(new TH1F("sigmaDCBZpSsmEE", "sigma of fitted function for Signal MC EE-EE", 338, 120., 3500.));


  TH1::SetDefaultSumw2(kTRUE);
}


