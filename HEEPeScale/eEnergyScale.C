//===============================================================
// RooFit Macro to perform unbinned likelihood fit to the Z peak
//================================================================

#include <sstream>
#include <iostream>
#include <utility>
#include <vector>

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include "tdrstyle.C"


#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooDCBShape.h"
#endif

using namespace RooFit;

RooWorkspace * cryBallBW(
                   const char *filename = "eScaleEvents920pb-1.root",
                   const char *treeName = "eleDataTree",
                   const int regions = 3,
                   double minMass = 60,
                   double maxMass = 120,
                   double mean_bw = 91.1876,
                   double gamma_bw = 2.4952,
                   double cutoff_cb = 1.0,
                   const char* plotOpt = "NEU",
                   const int nbins = 40,
                   const unsigned int fitModelType = 1,  // (0, 1) = (CB, double CB)
                   const unsigned int nCpu = 1)
{
  setTDRStyle();

  // make a workspace to return
  RooWorkspace *w = new RooWorkspace("w","workspace");

  RooRealVar mass("mass", "m(ee)", minMass, maxMass, "GeV");
  RooRealVar evtRegion("evtRegion", "Detector region", 0., 2.);
  RooRealVar ele1Eta("ele1Eta", "ele1Eta", -3., 3.);
  RooRealVar ele2Eta("ele2Eta", "ele2Eta", -3., 3.);

  // get invariant mass histo from file
  TFile input(filename, "read");
  input.cd();
  TTree *tree = (TTree *)gDirectory->Get(treeName);

  // Read dataset from tree
  RooDataSet *dataAll = new RooDataSet("dataAll", "complete dataset", RooArgSet(mass, evtRegion, ele1Eta, ele2Eta), Import(*tree));
  dataAll->Print();

  input.Close(); 

  RooDataSet *data = new RooDataSet("data", "pruned dataset", RooArgSet(mass));

  // prune dataset according to the desired regions. EE(4)/BE(2)/BB(1)
  //if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 0 && abs(ele1Eta) < 0.7 && abs(ele2Eta) < 0.7")); // both electrons |eta| < 0.7
  //if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 0 && abs(ele1Eta) >= 0.7 && abs(ele2Eta) >= 0.7")); // both electrons |eta| >= 0.7
  //if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 0 && (abs(ele1Eta) < 0.7 || abs(ele2Eta) < 0.7)")); // one electron |eta| < 0.7
  //if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 0 && (abs(ele1Eta) >= 0.7 || abs(ele2Eta) >= 0.7)")); // one electron |eta| >= 0.7
  //if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 0 && ((abs(ele1Eta) >= 0.7 && abs(ele2Eta) < 0.7) || (abs(ele2Eta) >= 0.7 && abs(ele1Eta) < 0.7))")); // one electron |eta| < 0.7 and one |eta| >= 0.7
  if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 0"));
  if ((regions & 2) == 2) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 1"));
  if ((regions & 4) == 4) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(mass), "evtRegion == 2"));
  if (data->numEntries() == 0) { // return if no falid region selected
    std::cout << "No event region selected. No data is imported." << std::endl;
    return w;
  }
  data->Print("v");

  // Build p.d.f.
  ////////////////////////////////////////////////
  //             Parameters                     //
  ////////////////////////////////////////////////

  //  Signal p.d.f. parameters
  //  Parameters for a Gaussian and a Crystal Ball Lineshape
  RooRealVar gBias ("#Deltam_{Gauss}", "Gaussian Bias", -.2, -20, 20,"GeV");
  RooRealVar gSigma("#sigma_{Gauss}","Gaussian Width", 2., 0.02, 20.,"GeV");
  //  Parameters for a Gaussian and a Crystal Ball Lineshape
  RooRealVar cbBias ("#Deltam_{CB}", "CB Bias", -.2, -20, 20,"GeV");
  RooRealVar cbSigma("#sigma_{CB}","CB Width", 2., 0.02, 20.,"GeV");
  RooRealVar cbCut  ("a_{CB}","CB Cut", 1., 0.1, 10.);
  RooRealVar cbPower("n_{CB}","CB Power", 2., 0.2, 50.);
  cbCut.setVal(cutoff_cb);
  //  cbCut.setConstant(kTRUE);  // her you can fix the cbCut_off

  RooRealVar dCBBias ("#Deltam_{DCB}", "Double CB Bias", -.2, -20, 20, "GeV");
  RooRealVar dCBSigma ("#sigma_{DCB}", "Double CB Width", 2., 0.02, 20., "GeV");
  RooRealVar dCBCutL ("al_{DCB}", "Double CB Cut left", 1., 0.1, 50.);
  RooRealVar dCBCutR ("ar_{DCB}", "Double CB Cut right", 1., 0.1, 50.);
  RooRealVar dCBPowerL ("nl_{DCB}", "Double CB Power left", 2., 0.2, 50.);
  RooRealVar dCBPowerR ("nr_{DCB}", "Double CB Power right", 2., 0.2, 50.);
  dCBCutL.setVal(cutoff_cb);
  dCBCutR.setVal(cutoff_cb);

  //  Parameters for Breit-Wigner
  RooRealVar bwMean("m_{Z}","BW Mean", 91.1876, "GeV");
  bwMean.setVal(mean_bw);
  RooRealVar bwWidth("#Gamma_{Z}", "BW Width", 2.4952, "GeV");
  //  RooRealVar bwWidth("#Gamma_{Z}", "BW Width", gamma_bw, "GeV");
  bwWidth.setVal(gamma_bw);

  // Keep Breit-Wigner parameters fixed to the PDG values
  bwMean.setConstant(kTRUE);
  bwWidth.setConstant(kTRUE);

  //  Background p.d.f. parameters
  // Parameters for exponential
  RooRealVar expRate("#lambda_{exp}", "Exponential Rate", -0.064, -1, 1);
  RooRealVar  c0("c_{0}", "c0", 1., 0., 50.);

  //  expRate.setConstant(kTRUE);

  // fraction of signal
  //  RooRealVar  frac("frac", "Signal Fraction", 0.1,0.,0.3.);
  RooRealVar  nsig("N_{S}", "#signal events", 524, 0.1, 10000000000.);
  nsig.setVal(data->numEntries());
  RooRealVar  nbkg("N_{B}", "#background events", 43, 1., 10000000.);
  
  //  nbkg.setVal(0);
  //  nbkg.setConstant(kTRUE);

  ////////////////////////////////////////////////
  //               P.D.F.s                      //
  ////////////////////////////////////////////////
 
  // Di-Electron mass signal p.d.f.
  RooBreitWigner bw("bw", "bw", mass, bwMean, bwWidth);
  RooGaussian gauss("gauss", "A Gaussian Lineshape", mass, gBias, gSigma);
  RooCBShape cball("cball", "A Crystal Ball Lineshape", mass, cbBias, cbSigma, cbCut, cbPower);
#ifdef __CINT__
  gROOT->ProcessLineSync(".x RooDCBShape.cxx+") ;
#endif  
  RooDCBShape dCBall("dCBall", "A double Crystal Ball lineshape", mass, dCBBias, dCBSigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR);
 
  // mass.setBins(100000, "fft");
  RooFFTConvPdf BWxGauss("BWxGauss","Breit - Wigner (X) Gaussian", mass, bw, gauss);
  RooFFTConvPdf BWxCB("BWxCB","Breit - Wigner (X) Crystal Ball", mass, bw, cball);
  RooFFTConvPdf BWxDCB("BWxDCB","Breit - Wigner (X) double Crystal Ball", mass, bw, dCBall);
 
  // Di-Electron mass background  p.d.f.
  RooExponential bg("bg","bkgd exp", mass, expRate);
  // RooPolynomial bg("bg","bkg Polynomial",mass,c0);
 
  // Di-Electron mass model p.d.f.
  //RooAddPdf model("model", "signal + background mass model", RooArgList(BWxDCB, bg), RooArgList(nsig, nbkg)); // Signal + Bkg
  RooArgList pdfs("pdfs");
  RooArgSet plotParam("plotParam");
  if (fitModelType == 0) {
    pdfs.add(BWxCB);
    plotParam.add(RooArgSet(cbBias, cbSigma));
  }
  else if (fitModelType == 1) {
    pdfs.add(BWxDCB);
    plotParam.add(RooArgSet(dCBBias, dCBSigma));
  }
  else {
    pdfs.add(BWxGauss);
    plotParam.add(RooArgSet(gBias, gSigma));
  }
  pdfs.Print("v");
  plotParam.Print("v");
  RooAddPdf model("model", "signal", pdfs, RooArgList(nsig)); // Signal Only
  
  
  TStopwatch t;
  t.Start();
  //model.fitTo(*data, FitOptions("mh"), NumCPU(2), Optimize(kFALSE), Timer(kTRUE), Range(70,110));
  //model.fitTo(*data, FitOptions("mh"), NumCPU(2), Optimize(kFALSE), Timer(kTRUE));
  model.fitTo(*data, NumCPU(nCpu), Timer(kTRUE), Save());
  t.Print();

  //fit->Print();
    
  //w->import(mass);
  w->import(model);
  w->import(bg);
  w->import(*data);  
  w->defineSet("plotParam", plotParam, kTRUE);
  //w->Print("v"); 
  return w;
}

void RunCryBall()
{
  // parameters //////////////////////////////////////////////////////////////
  const char *inFile = "eScaleEvents19712pb-1.root";
  const int lumi = 19712;
  double minMass = 60;
  double maxMass = 120;
  double mean_bw = 91.1876;
  double gamma_bw = 2.4952;
  double cutoff_cb = 1.0;
  const unsigned int fitModelType = 1; // (0, 1, 2) = (CB, double CB, Gaussian)

  const char* plotOpt = "NEU";
  const int nBins = 40;
  bool plotReg[4];
  plotReg[0] = 1; // EB-EB
  plotReg[1] = 0; // EB-EE
  plotReg[2] = 0; // EB-EB + EB-EE
  plotReg[3] = 0; // EE-EE

  bool logY = 0; // log y axis
  bool plotMC = 1;
  bool normalizeMC = 1;
  bool printLatexTable = 1;
  bool printHtmlTable = 1;

  const bool saveAsPdf = 1;
  const bool saveAsPng = 1;
  const bool saveAsRoot = 1;
  const char *fileNameExtra = "";
  //const char *fileNameExtra = "_runCv2";
  const char *plotDir = "./plots_20141009_bb/";
  ////////////////////////////////////////////////////////////////////////////

  vector<TString> treeNames;
  treeNames.push_back("eleDataTree");
  treeNames.push_back("eleDY20Tree");

  vector<TString> regTxt;
  regTxt.push_back("EB-EB");
  regTxt.push_back("EB-EE");
  regTxt.push_back("EB-EB & EB-EE");
  regTxt.push_back("EE-EE");

  std::vector<TString> regFileNameSuffix;
  regFileNameSuffix.push_back("BB");
  regFileNameSuffix.push_back("BE");
  regFileNameSuffix.push_back("BBBE");
  regFileNameSuffix.push_back("EE");
  regFileNameSuffix.push_back("BBEE");
  regFileNameSuffix.push_back("BEEE");
  regFileNameSuffix.push_back("BBBEEE");

  vector<pair<float, float> > logYplotRanges;
  logYplotRanges.push_back(make_pair(1000., 5000000.));
  logYplotRanges.push_back(make_pair(1000., 2000000.));
  logYplotRanges.push_back(make_pair(3000., 6000000.));
  logYplotRanges.push_back(make_pair(300., 700000.));

  // fit results
  vector<vector<Double_t> > fitResData;
  vector<vector<Double_t> > fitErrData;
  vector<vector<Double_t> > fitResDy20;
  vector<vector<Double_t> > fitErrDy20;

  for (unsigned int reg = 0; reg < 4; ++reg) {
    if (!plotReg[reg]) continue;
    unsigned int nCpu = 2;
    if (reg == 2) nCpu = 1;
    // fit
    RooWorkspace *dataWkSpc = cryBallBW(inFile, treeNames[0].Data(), reg + 1, minMass, maxMass, mean_bw, gamma_bw, cutoff_cb, plotOpt, nBins, fitModelType, nCpu);
    dataWkSpc->Print("v");
    RooRealVar *mass = dataWkSpc->var("mass");
    RooDataSet *data = (RooDataSet *)dataWkSpc->data("data");
    RooAddPdf *model = (RooAddPdf *)dataWkSpc->pdf("model");
    RooExponential *bkgPdf = (RooExponential *)dataWkSpc->pdf("bg");

    RooWorkspace *dy20WkSpc = cryBallBW(inFile, treeNames[1].Data(), reg + 1, minMass, maxMass, mean_bw, gamma_bw, cutoff_cb, plotOpt, nBins, fitModelType, nCpu);
    RooDataSet *dataMc = (RooDataSet *)dy20WkSpc->data("data");
    RooAddPdf *modelMc = (RooAddPdf *)dy20WkSpc->pdf("model");

    // plot the data
    TCanvas* c = new TCanvas("c" + regTxt[reg], "Unbinned Invariant Mass Fit " + regTxt[reg], 0, 0, 800, 600);
    if (logY) gPad->SetLogy();
    // plot the fit results
    RooPlot* plot = mass->frame(Range(minMass,maxMass),Bins(nBins));

    // calculate normalization factor
    Double_t normFactor = (Double_t)data->numEntries() / (Double_t)dataMc->numEntries();
    if (!normalizeMC) normFactor = 1.;

    vector<TString> delMName(1, "#Deltam_{CB}");
    vector<TString> sigmaName(1, "#sigma_{CB}");
    vector<TString> aLName(1, "a_{CB}");
    vector<TString> aRName(1, "a_{CB}");
    vector<TString> powerLName(1, "n_{CB}");
    vector<TString> powerRName(1, "n_{CB}");
    delMName.push_back("#Deltam_{DCB}");
    sigmaName.push_back("#sigma_{DCB}");
    aLName.push_back("al_{DCB}");
    aRName.push_back("ar_{DCB}");
    powerLName.push_back("nl_{DCB}");
    powerRName.push_back("nr_{DCB}");
    delMName.push_back("#Deltam_{Gauss}");
    sigmaName.push_back("#sigma_{Gauss}");

    // get the fit results
    RooRealVar *bias = dataWkSpc->var(delMName[fitModelType].Data());
    RooRealVar *sigma = dataWkSpc->var(sigmaName[fitModelType].Data());
    RooRealVar *cutL;
    RooRealVar *cutR;
    RooRealVar *powerL;
    RooRealVar *powerR;
    if (fitModelType == 0 || fitModelType == 1) {
      cutL = dataWkSpc->var(aLName[fitModelType].Data());
      cutR = dataWkSpc->var(aRName[fitModelType].Data());
      powerL = dataWkSpc->var(powerLName[fitModelType].Data());
      powerR = dataWkSpc->var(powerRName[fitModelType].Data());
    }
    RooRealVar *nsig = dataWkSpc->var("N_{S}");
    //RooRealVar *nbkg = dataWkSpc->var("N_{B}");
    //RooRealVar *lambdaExp = dataWkSpc->var("#lambda_{exp}");

    RooRealVar *biasMC = dy20WkSpc->var(delMName[fitModelType].Data());
    RooRealVar *sigmaMC = dy20WkSpc->var(sigmaName[fitModelType].Data());
    RooRealVar *cutLMC;
    RooRealVar *cutRMC;
    RooRealVar *powerLMC;
    RooRealVar *powerRMC;
    if (fitModelType == 0 || fitModelType == 1) {
      cutLMC = dy20WkSpc->var(aLName[fitModelType].Data());
      cutRMC = dy20WkSpc->var(aRName[fitModelType].Data());
      powerLMC = dy20WkSpc->var(powerLName[fitModelType].Data());
      powerRMC = dy20WkSpc->var(powerRName[fitModelType].Data());
    }
    RooRealVar *nsigMC = dy20WkSpc->var("N_{S}");
    //RooRealVar *nbkgMC = dy20WkSpc->var("N_{B}");
    //RooRealVar *lambdaExpMC = dy20WkSpc->var("#lambda_{exp}");

    // save fit results
    vector<Double_t> helper;
    helper.push_back(bias->getVal());
    helper.push_back(sigma->getVal());
    if (fitModelType == 0 || fitModelType == 1) {
      helper.push_back(cutL->getVal());
      if (fitModelType > 0) helper.push_back(cutR->getVal());
      helper.push_back(powerL->getVal());
      if (fitModelType > 0) helper.push_back(powerR->getVal());
    }
    helper.push_back(nsig->getVal());
    //helper.push_back(nbkg->getVal());
    //helper.push_back(lambdaExp->getVal());
    fitResData.push_back(helper);
    helper.clear();
    helper.push_back(bias->getError());
    helper.push_back(sigma->getError());
    if (fitModelType == 0 || fitModelType == 1) {
      helper.push_back(cutL->getError());
      if (fitModelType > 0) helper.push_back(cutR->getError());
      helper.push_back(powerL->getError());
      if (fitModelType > 0) helper.push_back(powerR->getError());
    }
    helper.push_back(nsig->getError());
    //helper.push_back(nbkg->getError());
    //helper.push_back(lambdaExp->getError());
    fitErrData.push_back(helper);
    helper.clear();
    helper.push_back(biasMC->getVal());
    helper.push_back(sigmaMC->getVal());
    if (fitModelType == 0 || fitModelType == 1) {
      helper.push_back(cutLMC->getVal());
      if (fitModelType > 0) helper.push_back(cutRMC->getVal());
      helper.push_back(powerLMC->getVal());
      if (fitModelType > 0) helper.push_back(powerRMC->getVal());
    }
    helper.push_back(nsigMC->getVal());
    //helper.push_back(nbkgMC->getVal());
    //helper.push_back(lambdaExpMC->getVal());
    fitResDy20.push_back(helper);
    helper.clear();
    helper.push_back(biasMC->getError());
    helper.push_back(sigmaMC->getError());
    if (fitModelType == 0 || fitModelType == 1) {
      helper.push_back(cutLMC->getError());
      if (fitModelType > 0) helper.push_back(cutRMC->getError());
      helper.push_back(powerLMC->getError());
      if (fitModelType > 0) helper.push_back(powerRMC->getError());
    }
    helper.push_back(nsigMC->getError());
    //helper.push_back(nbkgMC->getError());
    //helper.push_back(lambdaExpMC->getError());
    fitErrDy20.push_back(helper);
    helper.clear();

    // plot
    data->plotOn(plot, Name("dataData"));
    model->plotOn(plot, Name("modelData"));
    //bkgPdf->plotOn(plot, Name("bkgData"), Normalization(nbkg->getVal()/(nbkg->getVal()+nsig->getVal()), RooAbsReal::Relative));
    gStyle->SetTextColor(kBlue);
    model->paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(*dataWkSpc->set("plotParam")), Layout(0.12, 0.5, 0.95));
    if (plotMC) dataMc->plotOn(plot, Name("dataDy20"), Rescale(normFactor), MarkerColor(kRed));
    modelMc->plotOn(plot, Name("modelDy20"), Normalization(normFactor, RooAbsReal::Relative), LineStyle(kDashed), LineColor(kRed));
    gStyle->SetTextColor(kRed);  
    modelMc->paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(*dy20WkSpc->set("plotParam")), Layout(0.12, 0.5, 0.832168));
    gStyle->SetTextColor(kBlack);  

    if (logY) plot->GetYaxis()->SetRangeUser(logYplotRanges[reg].first, logYplotRanges[reg].second);

    //plot->drawAfter("bkgData", "modelData");
    plot->drawAfter("modelData", "modelDy20");
    plot->drawAfter("modelDy20", "dataDy20");
    plot->drawAfter("dataDy20", "dataData");

    int ndf = 5+2*fitModelType;
    if (fitModelType == 2) ndf = 3;
    cout << "Chi^2 data: " << plot->chiSquare("modelData", "dataData", ndf) << endl;
    cout << "Chi^2 DY20: " << plot->chiSquare("modelDy20", "dataDy20", ndf) << endl;

    // draw a legend
    TLegend *legend = new TLegend(0.75, 0.65, 0.91, 0.85);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    TLegendEntry *entry=legend->AddEntry("data", "DATA", "lp");
    entry->SetLineColor(kBlack);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1);
    entry=legend->AddEntry("datafit", "fit to data", "l");
    entry->SetLineColor(kBlue);
    entry->SetLineWidth(3);
    entry=legend->AddEntry("mc", "DY MC", "lp");
    entry->SetLineColor(kBlack);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(kRed);
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1);
    entry=legend->AddEntry("mcfit", "fit to DY MC", "l");
    entry->SetLineColor(kRed);
    entry->SetLineStyle(kDashed);
    entry->SetLineWidth(3);

    plot->addObject(legend);

    plot->Draw();
    //plot->Print("v");
   
    //TLatex *tex = new TLatex(0.63, 0.85, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    //tex->Draw();
    tex->SetTextSize(0.035);
    tex->DrawLatex(0.52, 0.88, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.12, 0.65, (const char *)regTxt[reg]);
    //tex->SetTextColor(kBlue);
    //tex->DrawLatex(0.2, 0.575, Form("#chi^{2} = %.2f",plot->chiSquare("modelData", "dataData")));
    //tex->SetTextColor(kRed);
    //tex->DrawLatex(0.2, 0.525, Form("#chi^{2} MC = %.2f",plot->chiSquare("modelDy20", "dataDy20")));
    //tex->SetTextColor(kBlack);
    //   tex->DrawLatex(0.6, 0.575, Form("CB Mean = %.1f",bias.getVal()));
    //   tex->DrawLatex(0.6, 0.525, Form("CB #sigma = %.1f",sigma.getVal()));

    // safe in various file formats
    stringstream sStream;
    sStream << plotDir << "zPeakFit60-120" << regFileNameSuffix[reg];
    if (fitModelType == 0) sStream << "_cb";
    else if (fitModelType == 1) sStream << "_dcb";
    else sStream << "_gauss";
    sStream << fileNameExtra;
    if (logY) sStream << "_log";
    sStream << "_" << lumi << "pb-1";
    TString saveFileName = sStream.str();
    if (saveAsPdf) c->Print(saveFileName + ".pdf", "pdf");
    if (saveAsPng) c->Print(saveFileName + ".png", "png");
    if (saveAsRoot) c->Print(saveFileName + ".root", "root");
  }

  // print the fit parameters
  unsigned int i = 0;
  cout << "===============================================================================" << endl;
  cout << "| Fit parameters                                                              |" << endl;
  for (unsigned int reg = 0; reg < 4; ++reg) {
    if (!plotReg[reg]) continue;
    cout << "===============================================================================" << endl;
    cout << "| " << regTxt[reg] << "   |         data         |           MC" << endl;
    cout << "===============================================================================" << endl;
    Float_t sigmaCbExtra = 0.;
    if (fitModelType == 0) {
      cout << "deltaM_CB: " << fitResData[i][0] << " +/- " << fitErrData[i][0] << "    |    " << fitResDy20[i][0] << " +/- " << fitErrDy20[i][0] << endl;
      cout << "sigma_CB:  " << fitResData[i][1] << " +/- " << fitErrData[i][1] << "    |    " << fitResDy20[i][1] << " +/- " << fitErrDy20[i][1] << endl;
      cout << "a_CB:      " << fitResData[i][2] << " +/- " << fitErrData[i][2] << "    |    " << fitResDy20[i][2] << " +/- " << fitErrDy20[i][2] << endl;
      cout << "n_CB:      " << fitResData[i][3] << " +/- " << fitErrData[i][3] << "    |    " << fitResDy20[i][3] << " +/- " << fitErrDy20[i][3] << endl;
      cout << "N_S:       " << fitResData[i][4] << " +/- " << fitErrData[i][4] << "    |    " << fitResDy20[i][4] << " +/- " << fitErrDy20[i][4] << endl;
      //cout << "N_B:       " << fitResData[i][5] << " +/- " << fitErrData[i][5] << "    |    " << fitResDy20[i][5] << " +/- " << fitErrDy20[i][5] << endl;
      //cout << "lambda_exp:" << fitResData[i][6] << " +/- " << fitErrData[i][6] << "    |    " << fitResDy20[i][6] << " +/- " << fitErrDy20[i][6] << endl << endl;
      cout << "(deltaM_CB(data) - deltaM_CB(MC)) / M(Z), %: " << (fitResData[i][0] - fitResDy20[i][0]) * 100. / mean_bw << " +/- " << sqrt(fitErrData[i][0] * fitErrData[i][0] + fitErrDy20[i][0] * fitErrDy20[i][0]) * 100. / mean_bw << endl << endl;
      if (fitResData[i][1] > fitResDy20[i][1]) sigmaCbExtra = sqrt(fitResData[i][1] * fitResData[i][1] - fitResDy20[i][1] * fitResDy20[i][1]);
      cout << "sigma_CB(extra), %: " << sigmaCbExtra * 100. / mean_bw << " +/- " << sqrt(fitErrData[i][1] * fitErrData[i][1] + fitErrDy20[i][1] * fitErrDy20[i][1]) * 100. / mean_bw << endl;
    }
    else if (fitModelType == 1) {
      cout << "deltaM_DCB: " << fitResData[i][0] << " +/- " << fitErrData[i][0] << "    |    " << fitResDy20[i][0] << " +/- " << fitErrDy20[i][0] << endl;
      cout << "sigma_DCB:  " << fitResData[i][1] << " +/- " << fitErrData[i][1] << "    |    " << fitResDy20[i][1] << " +/- " << fitErrDy20[i][1] << endl;
      cout << "al_DCB:      " << fitResData[i][2] << " +/- " << fitErrData[i][2] << "    |    " << fitResDy20[i][2] << " +/- " << fitErrDy20[i][2] << endl;
      cout << "ar_DCB:      " << fitResData[i][3] << " +/- " << fitErrData[i][3] << "    |    " << fitResDy20[i][3] << " +/- " << fitErrDy20[i][3] << endl;
      cout << "nl_DCB:      " << fitResData[i][4] << " +/- " << fitErrData[i][4] << "    |    " << fitResDy20[i][4] << " +/- " << fitErrDy20[i][4] << endl;
      cout << "nr_DCB:      " << fitResData[i][5] << " +/- " << fitErrData[i][5] << "    |    " << fitResDy20[i][5] << " +/- " << fitErrDy20[i][5] << endl;
      cout << "N_S:       " << fitResData[i][6] << " +/- " << fitErrData[i][6] << "    |    " << fitResDy20[i][6] << " +/- " << fitErrDy20[i][6] << endl;
      //cout << "N_B:       " << fitResData[i][7] << " +/- " << fitErrData[i][7] << "    |    " << fitResDy20[i][7] << " +/- " << fitErrDy20[i][7] << endl << endl;
      //cout << "lambda_exp:" << fitResData[i][8] << " +/- " << fitErrData[i][8] << "    |    " << fitResDy20[i][8] << " +/- " << fitErrDy20[i][8] << endl << endl;
      cout << "(deltaM_DCB(data) - deltaM_DCB(MC)) / M(Z), %: " << (fitResData[i][0] - fitResDy20[i][0]) * 100. / mean_bw << " +/- " << sqrt(fitErrData[i][0] * fitErrData[i][0] + fitErrDy20[i][0] * fitErrDy20[i][0]) * 100. / mean_bw << endl << endl;
      if (fitResData[i][1] > fitResDy20[i][1]) sigmaCbExtra = sqrt(fitResData[i][1] * fitResData[i][1] - fitResDy20[i][1] * fitResDy20[i][1]);
      cout << "sigma_CB(extra), %: " << sigmaCbExtra * 100. / mean_bw << " +/- " << sqrt(fitErrData[i][1] * fitErrData[i][1] + fitErrDy20[i][1] * fitErrDy20[i][1]) * 100. / mean_bw << endl;
    }
    else {
      cout << "deltaM_Gauss: " << fitResData[i][0] << " +/- " << fitErrData[i][0] << "    |    " << fitResDy20[i][0] << " +/- " << fitErrDy20[i][0] << endl;
      cout << "sigma_Gauss:  " << fitResData[i][1] << " +/- " << fitErrData[i][1] << "    |    " << fitResDy20[i][1] << " +/- " << fitErrDy20[i][1] << endl;
      cout << "(deltaM_Gauss(data) - deltaM_Gauss(MC)) / M(Z), %: " << (fitResData[i][0] - fitResDy20[i][0]) * 100. / mean_bw << " +/- " << sqrt(fitErrData[i][0] * fitErrData[i][0] + fitErrDy20[i][0] * fitErrDy20[i][0]) * 100. / mean_bw << endl << endl;
      if (fitResData[i][1] > fitResDy20[i][1]) sigmaCbExtra = sqrt(fitResData[i][1] * fitResData[i][1] - fitResDy20[i][1] * fitResDy20[i][1]);
      cout << "sigma_Gauss(extra), %: " << sigmaCbExtra * 100. / mean_bw << " +/- " << sqrt(fitErrData[i][1] * fitErrData[i][1] + fitErrDy20[i][1] * fitErrDy20[i][1]) * 100. / mean_bw << endl;
    }
    ++i;
  }
  cout << "===============================================================================" << endl << endl;

  // print the results as a HTML table
  if (printHtmlTable) {
    cout << "== HTML tables ================================================================" << endl;
    cout << "<br>Electron energy scale summary for different ECAL regions as measured at Z peak. Comparison between the data and MC." << endl;
    cout << "<table border=1>" << endl;
    if (fitModelType == 0) cout << "<tr align=center><td>region</td><td>&Delta;m<sub>CB</sub>(data), GeV/c<sup>2</sup></td><td>&Delta;m<sub>CB</sub>(MC), GeV/c<sup>2</sup></td><td>(&Delta;m<sub>CB</sub>(Data)-&Delta;m<sub>CB</sub>(MC))/m(Z<sup>0</sup>), %</td></tr>" << endl;
    else if (fitModelType == 1) cout << "<tr align=center><td>region</td><td>&Delta;m<sub>DCB</sub>(data), GeV/c<sup>2</sup></td><td>&Delta;m<sub>DCB</sub>(MC), GeV/c<sup>2</sup></td><td>(&Delta;m<sub>DCB</sub>(Data)-&Delta;m<sub>DCB</sub>(MC))/m(Z<sup>0</sup>), %</td></tr>" << endl;
    else cout << "<tr align=center><td>region</td><td>&Delta;m<sub>Gauss</sub>(data), GeV/c<sup>2</sup></td><td>&Delta;m<sub>Gauss</sub>(MC), GeV/c<sup>2</sup></td><td>(&Delta;m<sub>Gauss</sub>(Data)-&Delta;m<sub>Gaus</sub>(MC))/m(Z<sup>0</sup>), %</td></tr>" << endl;
    i = 0;
    for (unsigned int reg = 0; reg < 4; ++reg) {
      if (!plotReg[reg]) continue;
      Float_t relShift = (fitResData[i][0] - fitResDy20[i][0]) * 100. / mean_bw;
      Float_t relShiftErr = sqrt(fitErrData[i][0] * fitErrData[i][0] + fitErrDy20[i][0] * fitErrDy20[i][0]) * 100. / mean_bw;
      cout << "<tr align=right><td align=center>" << regTxt[reg];
      printf("</td><td> %.2f &pm; %.2f </td><td> %.2f &pm; %.2f </td><td> %.2f &pm; %.2f </td></tr>\n", fitResData[i][0], fitErrData[i][0], fitResDy20[i][0], fitErrDy20[i][0], relShift, relShiftErr);
      ++i;
    }
    cout << "</table>" << endl;

    cout << endl;
    cout << "<br>Di-electron mass resolution for different ECAL regions as measured at the Z peak for the data and MC." << endl;
    cout << "<table border=1>" << endl;
    if (fitModelType == 0) cout << "<tr align=center><td>region</td><td>&sigma;<sub>CB</sub>(data), %</td><td>&sigma;<sub>CB</sub>(MC), %</td><td>&sigma;<sub>CB</sub>(extra), %</td></tr>" << endl;
    else if (fitModelType == 1) cout << "<tr align=center><td>region</td><td>&sigma;<sub>DCB</sub>(data), %</td><td>&sigma;<sub>DCB</sub>(MC), %</td><td>&sigma;<sub>DCB</sub>(extra), %</td></tr>" << endl;
    else cout << "<tr align=center><td>region</td><td>&sigma;<sub>Gauss</sub>(data), %</td><td>&sigma;<sub>Gauss</sub>(MC), %</td><td>&sigma;<sub>Gauss</sub>(extra), %</td></tr>" << endl;
    i = 0;
    for (unsigned int reg = 0; reg < 4; ++reg) {
      if (!plotReg[reg]) continue;
      cout << "<tr align=right><td align=center>" << regTxt[reg];
      Float_t relFac = 100. / mean_bw;
      Float_t sigmaCbExtra = 0.;
      if (fitResData[i][1] > fitResDy20[i][1]) sigmaCbExtra = sqrt(fitResData[i][1] * fitResData[i][1] - fitResDy20[i][1] * fitResDy20[i][1]) * relFac;
      Float_t sigmaCbExtraErr = sqrt(fitErrData[i][1] * fitErrData[i][1] + fitErrDy20[i][1] * fitErrDy20[i][1]) * relFac;
      printf("</td><td>%.2f &pm; %.2f </td><td>%.2f &pm; %.2f </td><td>%.2f &pm; %.2f </td></tr>\n", fitResData[i][1] * relFac, fitErrData[i][1] * relFac, fitResDy20[i][1] * relFac, fitErrDy20[i][1] * relFac, sigmaCbExtra, sigmaCbExtraErr);
      ++i;
    }
    cout << "</table>" << endl;
    cout << "===============================================================================" << endl << endl;
  }

  // print the results as a Latex table
  if (printLatexTable) {
    cout << "== LaTex tables ===============================================================" << endl;
    cout << "\\begin{table}[t]" << endl;
    cout << "\\begin{center}" << endl;
    cout << "\\begin{tabular}{|c|c|c|c|}" << endl;
    cout << "\\hline" << endl;
    if (fitModelType == 0) cout << "region & $\\Delta{m}_{CB}(data), GeV/c^2 $ & $\\Delta{m}_{CB}(MC), GeV/c^2 $ & $\\frac{\\Delta{m}_{CB}(Data)-\\Delta{m}_{CB}(MC)}{m(Z^{0})}, \\% $ \\\\" << endl;
    else if (fitModelType == 1) cout << "region & $\\Delta{m}_{DCB}(data), GeV/c^2 $ & $\\Delta{m}_{DCB}(MC), GeV/c^2 $ & $\\frac{\\Delta{m}_{DCB}(Data)-\\Delta{m}_{DCB}(MC)}{m(Z^{0})}, \\% $ \\\\" << endl;
    else cout << "region & $\\Delta{m}_{Gauss}(data), GeV/c^2 $ & $\\Delta{m}_{Gauss}(MC), GeV/c^2 $ & $\\frac{\\Delta{m}_{Gauss}(Data)-\\Delta{m}_{Gauss}(MC)}{m(Z^{0})}, \\% $ \\\\" << endl;
    cout << "\\hline" << endl;
    i = 0;
    for (unsigned int reg = 0; reg < 4; ++reg) {
      if (!plotReg[reg]) continue;
      Float_t relShift = (fitResData[i][0] - fitResDy20[i][0]) * 100. / mean_bw;
      Float_t relShiftErr = sqrt(fitErrData[i][0] * fitErrData[i][0] + fitErrDy20[i][0] * fitErrDy20[i][0]) * 100. / mean_bw;
      cout << regTxt[reg].ReplaceAll("&", "\\&");
      printf("          & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n", fitResData[i][0], fitErrData[i][0], fitResDy20[i][0], fitErrDy20[i][0], relShift, relShiftErr);
      ++i;
    }
    cout << "\\hline" << endl;
    cout << "\\end{tabular}" << endl;
    cout << "\\caption{Electron energy scale summary for different ECAL regions as measured at the Z peak. Comparison between the data and MC.}" << endl;
    cout << "\\label{tab:energyScale}" << endl;
    cout << "\\end{center}" << endl;
    cout << "\\end{table}" << endl;

    cout << endl;
    cout << "\\begin{table}[t]" << endl;
    cout << "\\begin{center}" << endl;
    cout << "\\begin{tabular}{|c|c|c|c|}" << endl;
    cout << "\\hline" << endl;
    if (fitModelType == 0) cout << "region & $\\frac{\\sigma_{CB}(data)}{m(Z^{0})}, \\%$ & $\\frac{\\sigma_{CB}(MC)}{m(Z^{0})}, \\%$ & $\\sigma_{CB}(extra), \\%$ \\\\" << endl;
    else if (fitModelType == 1) cout << "region & $\\frac{\\sigma_{DCB}(data)}{m(Z^{0})}, \\%$ & $\\frac{\\sigma_{DCB}(MC)}{m(Z^{0})}, \\%$ & $\\sigma_{DCB}(extra), \\%$ \\\\" << endl;
    else cout << "region & $\\frac{\\sigma_{Gauss}(data)}{m(Z^{0})}, \\%$ & $\\frac{\\sigma_{Gauss}(MC)}{m(Z^{0})}, \\%$ & $\\sigma_{Gauss}(extra), \\%$ \\\\" << endl;
    cout << "\\hline" << endl;
    i = 0;
    for (unsigned int reg = 0; reg < 4; ++reg) {
      if (!plotReg[reg]) continue;
      cout << regTxt[reg];
      Float_t relFac = 100. / mean_bw;
      Float_t sigmaCbExtra = 0.;
      if (fitResData[i][1] > fitResDy20[i][1]) sigmaCbExtra = sqrt(fitResData[i][1] * fitResData[i][1] - fitResDy20[i][1] * fitResDy20[i][1]) * relFac;
      Float_t sigmaCbExtraErr = sqrt(fitErrData[i][1] * fitErrData[i][1] + fitErrDy20[i][1] * fitErrDy20[i][1]) * relFac;
      printf("          & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n", fitResData[i][1] * relFac, fitErrData[i][1] * relFac, fitResDy20[i][1] * relFac, fitErrDy20[i][1] * relFac, sigmaCbExtra, sigmaCbExtraErr);
      ++i;
    }
    cout << "\\hline" << endl;
    cout << "\\end{tabular}" << endl;
    cout << "\\caption{Di-electron mass resolution for different ECAL regions as measured at the Z peak for the data and MC. The resolutions for data and MC are normalized to the Z mass.}" << endl;
    cout << "\\label{tab:resolution}" << endl;
    cout << "\\end{center}" << endl;
    cout << "\\end{table}" << endl;
    cout << "===============================================================================" << endl << endl;
  }
}

