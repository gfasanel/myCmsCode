#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TF1.h"
#include "THStack.h"
#include "TColor.h"
#include "TString.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THashList.h"

#include "makeHistoFromBranch.C"

// class for plotting of control variable plots
class ContVarPlot {
  public:
    TString fName;
    TString fTitle;
    TString fXaxisTitle;
    float fXmin;
    float fXmax;
    int fNbins; // variable binning if fNbins < 1
    bool fLogPlot;
    bool fUnderFlow;
    bool fOverFlow;

    ContVarPlot(const char* name, const char* title, const char* xAxisTitle, float xmin, float xmax, int nBins, bool logPlot=0, bool underFlow=0, bool overFlow=0) : fName(name), fTitle(title), fXaxisTitle(xAxisTitle), fXmin(xmin), fXmax(xmax), fNbins(nBins), fLogPlot(logPlot), fUnderFlow(underFlow), fOverFlow(overFlow) { } 
    ~ContVarPlot() { }
};

void macro_MakeTopReweightShapes(unsigned int var = 0, int sig = 0, unsigned int reg = 0)
{ 
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("./emuSpec_MuGammaTrg_19703pb-1.root");
  TFile input("./emuSpec_singleMuTrg_19706pb-1.root");
  input.cd();

  TParameter<float> *lumi = (TParameter<float> *)input.Get("lumi");

  const bool logPlotX = 0; // only for variable binning

  // plot style
  int ttbarColour = kGreen;
  int ttbarRewUpColour = kBlue;
  int ttbarRewDownColour = kRed;
  int ttbarRewColour = kMagenta;

  // output file formats
  const bool saveShapeHistos = 1;
  const bool saveSpec = 0;
  const bool saveAsPdf = 0;
  const bool saveAsPng = 1;
  const bool saveAsRoot = 0;
  const char *fileNameExtra = "";
  //const char *fileNameExtra = "madgraphTTbar_";
  const char *plotDir = "../plots/";

  int font = 42; //62
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  std::vector<ContVarPlot> testPlots;
  testPlots.reserve(40);
  // flags: logPlot | underflow in first bin | overflow in last bin
  testPlots.push_back(ContVarPlot("mass", "e#mu invariant mass", "m(e#mu) [GeV]", 0., 3000., 3000, 1, 0, 0));

  TString sign[7] = {"_e-mu+", "_e+mu-", "_OS", "", "_SS", "_++", "_--"};
  TString nameSign[7] = {" e-mu+", " e+mu-", " OS", "", " SS", " ++", " --"};
  TString region[3] = {"", "_EB", "_EE"};
  TString nameReg[3] = {"", " EB", " EE"};

  // plot a list with possible test histograms
  if (var > testPlots.size()) var = 0;
  if (var == 0) {
    cout << "Use macro_MakeMultiSampleClosureTest.C(x, y, z) with x {1-" << testPlots.size() << "} being the \n";
    cout << "number of the test histogram to plot, y {-3 - +3} selects the charge with \n";
    cout << "the scheme e-mu: -+, +-, OS, ALL, SS, ++, -- and z {0-2} selects the \n";
    cout << "detector EB+EE, EB or EE events for the electron." << endl;
    cout << "-----------------------------------" << endl;
    for (unsigned int i = 0; i < testPlots.size(); ++i) {
      cout << testPlots[i].fTitle;
      unsigned int j = 35 - testPlots[i].fTitle.Sizeof();
      if (i > 8) --j;
      for (; j > 0; --j)
        cout << " ";
      cout << i + 1 << endl;
    }
    cout << "-----------------------------------" << endl;
    input.Close();
    return;
  }
  --var;
  // sanity checks for sign and region input
  if (abs(sig) > 3) sig = 0;
  if (reg > 2) reg = 0;
  // the makeHistoFromBranch function uses a different scheme for barrel and endcap selection
  unsigned int histoReg = 2;
  if (reg == 1) histoReg = 0;
  else if (reg == 2) histoReg = 1;

  std::vector<TH1F *> emuTest_ttbar;
  std::vector<TH1F *> emuTest_ttbar700to1000;
  std::vector<TH1F *> emuTest_ttbar1000up;
  std::vector<TH1F *> emuTest_ttbarPriv600up;
  std::vector<TH1F *> emuTest_ttbar_rew;
  std::vector<TH1F *> emuTest_ttbar700to1000_rew;
  std::vector<TH1F *> emuTest_ttbar1000up_rew;
  std::vector<TH1F *> emuTest_ttbarPriv600up_rew;
  std::vector<TH1F *> emuTest_ttbar_tot;
  std::vector<TH1F *> emuTest_ttbar_tot_rew;
  std::vector<TH1F *> mass_ttbar_topPtReweightUp;
  std::vector<TH1F *> mass_ttbar_topPtReweightDown;

  // configure plot style
  TString testVar = testPlots[var].fName;
  bool logPlot = testPlots[var].fLogPlot;

  // define the binning
  std::vector<float> binning;
  bool normToBin = false;
  float lowBin = testPlots[var].fXmin;
  float highBin = testPlots[var].fXmax;
  int nBins = testPlots[var].fNbins;
  if (testPlots[var].fNbins < 1) {
    normToBin = true;
    if (logPlotX) {
      //for (float bin = 0.; bin < 100.; bin += 5.)
      //  binning.push_back(bin);
      for (float bin = 50.; bin < 200.; bin += 10.)
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
      binning.push_back(highBin);
    } else {
      for (float bin = lowBin; bin < 200.; bin += 20.)
        binning.push_back(bin);
      for (float bin = 200.; bin < 400.; bin += 40.)
        binning.push_back(bin);
      for (float bin = 400.; bin < 700.; bin += 50.)
        binning.push_back(bin);
      for (float bin = 700.; bin < 1000.; bin += 75.)
        binning.push_back(bin);
      for (float bin = 1000.; bin < 1100.; bin += 100.)
        binning.push_back(bin);
      for (float bin = 1100.; bin <= highBin; bin += 200.)
        binning.push_back(bin);
    }
  } else {
    for (float bin = lowBin; bin <= highBin; bin += (highBin - lowBin) / nBins)
      binning.push_back(bin);
  }
  nBins = binning.size();

  input.cd();

  THashList *mcWeights = (THashList *)input.Get("mcWeights");
  TParameter<float> *ttbarMcWeight = (TParameter<float> *)mcWeights->FindObject("ttbar");
  TParameter<float> *ttbarMcWeight700to1000 = (TParameter<float> *)mcWeights->FindObject("ttbar700to1000");
  TParameter<float> *ttbarMcWeight1000up = (TParameter<float> *)mcWeights->FindObject("ttbar1000up");
  TParameter<float> *ttbarMcWeight600up = (TParameter<float> *)mcWeights->FindObject("ttbarPriv600up");

  unsigned int topFlags = 0x1DF; // flags: [apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
  // get the histograms
  vector<const char *> cutVars;
  vector<float> lowCuts;
  vector<float> highCuts;
  vector<float> mcWeightsForCutRanges;
 
  // get combined ttbar histograms with MC weights for different samples mixed on an event-by-event basis
  cutVars.push_back("");
  lowCuts.push_back(0.);
  highCuts.push_back(1.e9);
  mcWeightsForCutRanges.push_back(ttbarMcWeight->GetVal());
  cutVars.push_back("genMTtbar");
  lowCuts.push_back(700.);
  highCuts.push_back(1000.);
  mcWeightsForCutRanges.push_back(ttbarMcWeight700to1000->GetVal());
  cutVars.push_back("genMTtbar");
  lowCuts.push_back(1000.);
  highCuts.push_back(1.e9);
  mcWeightsForCutRanges.push_back(ttbarMcWeight1000up->GetVal());
  cutVars.push_back("emu_mass");
  lowCuts.push_back(600.);
  highCuts.push_back(1.e9);
  mcWeightsForCutRanges.push_back(ttbarMcWeight600up->GetVal());
  emuTest_ttbar.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  emuTest_ttbar700to1000.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  emuTest_ttbar1000up.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  emuTest_ttbarPriv600up.push_back(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  topFlags += (1<<10);
  emuTest_ttbar_rew.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  emuTest_ttbar700to1000_rew.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  emuTest_ttbar1000up_rew.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  emuTest_ttbarPriv600up_rew.push_back(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  sig+=3;

  emuTest_ttbar_tot.push_back((TH1F *)emuTest_ttbar.back()->Clone("emuTest_ttbar_tot"));
  emuTest_ttbar_tot.back()->Add(emuTest_ttbar700to1000.back());
  emuTest_ttbar_tot.back()->Add(emuTest_ttbar1000up.back());
  emuTest_ttbar_tot.back()->Add(emuTest_ttbarPriv600up.back());
  emuTest_ttbar_tot_rew.push_back((TH1F *)emuTest_ttbar_rew.back()->Clone("emuTest_ttbar_tot_rew"));
  emuTest_ttbar_tot_rew.back()->Add(emuTest_ttbar700to1000_rew.back());
  emuTest_ttbar_tot_rew.back()->Add(emuTest_ttbar1000up_rew.back());
  emuTest_ttbar_tot_rew.back()->Add(emuTest_ttbarPriv600up_rew.back());

  // calculate up and down shapes
  mass_ttbar_topPtReweightDown.push_back((TH1F *)emuTest_ttbar_tot_rew.back()->Clone("mass_ttbar_topPtReweightDown"));
  mass_ttbar_topPtReweightUp.push_back((TH1F *)emuTest_ttbar_tot.back()->Clone("mass_ttbar_topPtReweightUp"));
  // decided to take Up fluctuation spectrum same as nominal
  //mass_ttbar_topPtReweightUp.back()->Add(emuTest_ttbar_tot.back());
  //mass_ttbar_topPtReweightUp.back()->Add(emuTest_ttbar_tot_rew.back(), -1.);

  TCanvas *emuPlot;
  TPad *specPad;
  emuPlot = new TCanvas("emuPlot" + testVar + sign[sig] + region[reg], "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg], 100, 100, 900, 600);
  specPad = new TPad("specPad" + testVar + sign[sig] + region[reg], "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg], 0., 0., 1., 1.);
  specPad->SetBottomMargin(0.12);
  specPad->SetBorderMode(0);
  specPad->SetBorderSize(2);
  specPad->SetFrameBorderMode(0);
  specPad->SetFillColor(0);
  specPad->SetFrameFillColor(0);
  if (logPlot) specPad->SetLogy();
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

  //// plot spectrum
  emuTest_ttbar_tot.back()->SetLineColor(ttbarColour);
  emuTest_ttbar_tot_rew.back()->SetLineColor(ttbarRewColour);
  mass_ttbar_topPtReweightUp.back()->SetLineColor(ttbarRewUpColour);
  mass_ttbar_topPtReweightDown.back()->SetLineColor(ttbarRewDownColour);

  emuTest_ttbar_tot.back()->Draw("hist");
  emuTest_ttbar_tot.back()->GetXaxis()->SetTitle(testPlots[var].fXaxisTitle);
  emuTest_ttbar_tot.back()->GetXaxis()->SetTitleFont(font);
  emuTest_ttbar_tot.back()->GetXaxis()->SetTitleSize(0.047);
  emuTest_ttbar_tot.back()->GetXaxis()->SetLabelSize(0.05);
  emuTest_ttbar_tot.back()->GetXaxis()->SetLabelFont(font);
  emuTest_ttbar_tot.back()->GetXaxis()->SetMoreLogLabels();
  emuTest_ttbar_tot.back()->GetXaxis()->SetNoExponent();
  emuTest_ttbar_tot.back()->GetYaxis()->SetTitle("Events");
  emuTest_ttbar_tot.back()->GetYaxis()->SetTitleFont(font);
  emuTest_ttbar_tot.back()->GetYaxis()->SetTitleSize(0.047);
  emuTest_ttbar_tot.back()->GetYaxis()->SetTitleOffset(1.1);
  emuTest_ttbar_tot.back()->GetYaxis()->SetLabelFont(font);
  emuTest_ttbar_tot.back()->GetYaxis()->SetLabelSize(0.05);
  //emuTest_ttbar_tot_rew.back()->Draw("histsame");
  mass_ttbar_topPtReweightUp.back()->Draw("histsame");
  mass_ttbar_topPtReweightDown.back()->Draw("histsame");

  // legent and labels
  TLegend legend(0.5, 0.64, 0.75, 0.885);
  legend.SetTextFont(font);
  legend.SetTextSize(0.042);
  legend.SetBorderSize(0);
  legend.SetLineColor(1);
  legend.SetLineStyle(1);
  legend.SetLineWidth(1);
  legend.SetFillColor(19);
  legend.SetFillStyle(0);
  legend.AddEntry(emuTest_ttbar_tot.back(), "t#bar{t} no reweighting", "l");
  legend.AddEntry(mass_ttbar_topPtReweightDown.back(), "t#bar{t} reweighting (Down)", "l");
  //legend.AddEntry(mass_ttbar_topPtReweightUp.back(), "2x t#bar{t}_{no rew} - t#bar{t}_{rew} (Up)", "l");
  legend.SetBorderSize(0);
  legend.DrawClone("same");
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(font);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.042);
  tex->DrawLatex(0.596, 0.937, "CMS Simulation, 8 TeV");
  //if (m == 1) tex->DrawLatex(0.430, 0.849, "e in barrel");
  //if (m == 2) tex->DrawLatex(0.430, 0.849, "e in endcap");
  stringstream sStream;
  sStream << testPlots[var].fTitle << nameSign[sig].Data() << nameReg[reg].Data();
  tex->DrawLatex(0.109, 0.937, sStream.str().c_str());

  // safe in various file formats
  if (saveSpec) {
    sStream.str("");
    sStream << plotDir << "emuControlSpec_" << testVar << sign[sig] << region[reg] << "_" << fileNameExtra << lumi->GetVal() << "pb-1";
    TString saveFileName = sStream.str();
    if (saveAsPdf) emuPlot->Print(saveFileName + ".pdf", "pdf");
    if (saveAsPng) emuPlot->Print(saveFileName + ".png", "png");
    if (saveAsRoot) emuPlot->Print(saveFileName + ".root", "root");
  }

  if (saveShapeHistos) {
    TFile *shapeFile = new TFile("histograms.root", "update");
    mass_ttbar_topPtReweightUp.back()->Write();
    mass_ttbar_topPtReweightDown.back()->Write();
    shapeFile->Close();
  }

  cout << testVar.Data() << nameSign[sig].Data() << nameReg[reg].Data() << " plotted" << endl;
}

