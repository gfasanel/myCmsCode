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

void macro_MakeMultiSampleClosureTest(unsigned int var = 0, int sig = 0, unsigned int reg = 0)
{ 
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("./emuSpec_19703pb-1.root", "open");
  TFile input("./emuSpec_singleMuTrg_19706pb-1.root", "open");
  input.cd();

  TParameter<float> *lumi = (TParameter<float> *)input.Get("lumi");

  const bool topReweighting = 0;
  const bool plot_ttbar = 1; // plot ttbar or WW

  const bool plotPull = 1; // plot (ref - test) / test
  const bool pullGridY = 1; // grid lines on y axis of pull plot
  const bool logPlotX = 0; // only for variable binning

  float yRangeMinRatio = -0.7;
  float yRangeMaxRatio = 1.3;

  // plot style
  int ttbarColour = kGreen;
  int ttbar700to1000Colour = kBlue;
  int ttbar1000upColour = kRed;
  int ttbarPriv600upColour = kCyan;
  int wwColour = kGreen;
  int wwPriv600upEminusMuPlusColour = kMagenta;
  int wwPriv600upEplusMuMinusColour = kCyan;

  // output file formats
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
  testPlots.push_back(ContVarPlot("mass", "e#mu invariant mass", "m(e#mu) [GeV]", 0., 1500., 0, 1, 0, 0));
  testPlots.push_back(ContVarPlot("pfMet", "PF MET", "E^{T}_{miss} [GeV]", 0., 500., 50, 1, 1, 1));
  testPlots.push_back(ContVarPlot("nVtx", "Number of primary vertices", "# PV", 0., 50., 50, 0, 1, 1));
  testPlots.push_back(ContVarPlot("rho", "rho", "#rho", 0., 50., 50, 0, 1, 1));
  testPlots.push_back(ContVarPlot("numOfJetsPt20", "Number of jets > 20 GeV", "# jets_{p_{T}>20}", 0., 20., 20, 0, 1, 1));
  testPlots.push_back(ContVarPlot("numOfJetsPt30", "Number of jets > 30 GeV", "# jets_{p_{T}>30}", 0., 15., 15, 0, 1, 1));
  testPlots.push_back(ContVarPlot("eDxyMinusMuDxy", "|Dxy_e - Dxy_mu|", "|#Deltaxy_{e}-#Deltaxy_{#mu}|", 0., 0.06, 40, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eDzMinusMuDz", "|Dz_e - Dz_mu|", "|#Deltaz_{e}-#Deltaz_{#mu}|", 0., 0.5, 40, 1, 1, 1));
  testPlots.push_back(ContVarPlot("dEta", "|eta_e - eta_mu|", "|#eta_{e}-#eta_{#mu}|", 0., 5., 50, 0, 1, 1));
  testPlots.push_back(ContVarPlot("dPhi", "|phi_e - phi_mu|", "|#varphi_{e}-#varphi_{#mu}|", 0., 3.2, 64, 0, 1, 1));
  testPlots.push_back(ContVarPlot("eleEt", "Electron Et", "E_{T}^{e} [GeV]", 0., 500., 50, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleEta", "Electron eta", "#eta_{e}", -2.5, 2.5, 50, 0, 1, 1));
  testPlots.push_back(ContVarPlot("elePhi", "Electron phi", "#phi_{e}", -3.2, 3.2, 64, 0, 1, 1));
  testPlots.push_back(ContVarPlot("eleDEta", "Electron dEta", "#Delta#eta", -0.007, 0.007, 29, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleDPhi", "Electron dPhi", "#Delta#varphi", -0.06, 0.06, 50, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleHOE", "Electron H/E", "H/E", 0., 0.051, 51, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleE1x5overE5x5", "Electron E1x5/E5x5", "E1x5/E5x5", 0., 1.5, 30, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleE2x5overE5x5", "Electron E2x5/E5x5", "E2x5/E5x5", 0., 1.5, 30, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleSigmaIEIE", "Electron #sigma_i#etai#eta", "", 0., 0.04, 40, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleEcalIso", "Electron ECAL iso", "ECAL iso", 0., 15., 30, 0, 0, 0));
  testPlots.push_back(ContVarPlot("eleHcalIso1", "Electron HCAL iso1", "HCAL iso1", 0., 15., 30, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleHcalIso2", "Electron HCAL iso2", "HCAL iso2", 0., 15., 30, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleHeepIso", "Electron HEEP iso", "HEEP iso", 0., 20., 40, 1, 0, 1));
  testPlots.push_back(ContVarPlot("eleTrkIso", "Electron track iso", "trk_{iso}", 0., 5.1, 51, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleLostHits", "Electron number of lost hits", "", 0., 2., 2, 0, 0, 0));
  testPlots.push_back(ContVarPlot("eleDXYFstPVtx", "Electron Dxy at first PV", "#Deltaxy_{e}", -0.05, 0.05, 40, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleDZFstPVtx", "Electron Dz at first PV", "#Deltaz_{e}", -1., 1., 80, 1, 1, 1));
  testPlots.push_back(ContVarPlot("etEleOPtMu", "Et_e / pt_mu", "E^{e}_{T} / p^{#mu}_{T}", 0., 13., 26, 1, 1, 1));
  testPlots.push_back(ContVarPlot("lepPtPlusOPtMinus", "lepton pt+ / pt-", "p^{+}_{T}/p^{-}_{T}", 0., 10., 20, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eChTimesMuCh", "charge_e * charge_mu", "charge_{e}*charge_{#mu}", -1., 2., 3, 1, 0, 0));
  testPlots.push_back(ContVarPlot("muPt", "Muon pt", "p_{T}^{#mu} [GeV]", 0., 500., 50, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muPtErr", "Muon pt error", "p_{T}^{#mu} [GeV]", 0., 50., 50, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muEta", "Muon eta", "#eta_{#mu}", -2.5, 2.5, 50, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muPhi", "Muon phi", "#varphi_{#mu}", -3.2, 3.2, 64, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muHitLayers", "Muon layers with hits", "# layers hits", 4., 20., 16, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muTrkHits", "Muon tracker hits", "# tracker hits", 0., 40., 40, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muPxlHits", "Muon pixel hits", "# pixel hits", 0., 10., 10, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muMuHits", "Muon system hits", "# #mu system hits", 0., 55., 55, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muDXYFstPVtx", "Muon Dxy at first PV", "#Deltaxy_{#mu}", -0.2, 0.2, 80, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muDZFstPVtx", "Muon Dz at first PV", "#Deltaz_{#mu}", -0.1, 0.1, 40, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muNSeg", "Muon number of segments", "# segments", 0., 7., 7, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muTrkIso03", "Muon track iso 03", "#mu trk iso", 0., 20., 40, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muIsoCombRel", "Muon combined iso / pt", "#mu (iso_{em}+iso_{had}+iso_{trk})/p_{T}", 0., 3.5, 35, 1, 1, 1));

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
    cout << "Use PlotRange(y, z, x-start, x-stop) to plot a range of test histograms." << endl;
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

  std::vector<TH1F *> emuTest_ref;
  std::vector<TH1F *> emuTest_ttbar;
  std::vector<TH1F *> emuTest_ttbar700to1000;
  std::vector<TH1F *> emuTest_ttbar1000up;
  std::vector<TH1F *> emuTest_ttbarPriv600up;
  std::vector<TH1F *> emuTest_ww;
  std::vector<TH1F *> emuTest_wwPriv600upEminusMuPlus;
  std::vector<TH1F *> emuTest_wwPriv600upEplusMuMinus;

  input.cd();

  // configure plot style
  TString testVar = testPlots[var].fName;
  bool logPlot = testPlots[var].fLogPlot;
  const bool overflowBin = testPlots[var].fOverFlow;
  const bool underflowBin = testPlots[var].fUnderFlow; 

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


  THashList *mcWeights = (THashList *)input.Get("mcWeights");
  THashList *nGenEvents = (THashList *)input.Get("nGenEvents");
  TParameter<float> *ttbarMcWeight = (TParameter<float> *)mcWeights->FindObject("ttbar");
  TParameter<float> *ttbarMcWeight700to1000 = (TParameter<float> *)mcWeights->FindObject("ttbar700to1000");
  TParameter<float> *ttbarMcWeight1000up = (TParameter<float> *)mcWeights->FindObject("ttbar1000up");
  TParameter<float> *ttbarMcWeight600up = (TParameter<float> *)mcWeights->FindObject("ttbarPriv600up");
  TParameter<float> *wwMcWeight = (TParameter<float> *)mcWeights->FindObject("ww");
  TParameter<float> *wwMcWeight600upEminusMuPlus = (TParameter<float> *)mcWeights->FindObject("wwPriv600upEminusMuPlus");
  TParameter<float> *wwMcWeight600upEplusMuMinus = (TParameter<float> *)mcWeights->FindObject("wwPriv600upEplusMuMinus");
  TParameter<float> *wwNGen600upEminusMuPlus = (TParameter<float> *)nGenEvents->FindObject("wwPriv600upEminusMuPlus");
  TParameter<float> *wwNGen600upEplusMuMinus = (TParameter<float> *)nGenEvents->FindObject("wwPriv600upEplusMuMinus");

  unsigned int flags = 0x1DF; // flags: [apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
  unsigned int topFlags = flags;
  if (topReweighting) topFlags += (1<<10);
  // get the histograms
  vector<const char *> cutVars;
  vector<float> lowCuts;
  vector<float> highCuts;
  vector<float> mcWeightsForCutRanges;
  TH1F *refOverTestHist;
  if (plot_ttbar) {
    refOverTestHist = (TH1F *)MakeHistoFromBranch(&input, "emuTree_ttbar", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin);
    emuTest_ref.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, topFlags, normToBin));
  } else {
    refOverTestHist = (TH1F *)MakeHistoFromBranch(&input, "emuTree_ww", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, flags, normToBin);
    emuTest_ref.push_back(MakeHistoFromBranch(&input, "emuTree_ww", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, flags, normToBin));
  }
 
  // get combined ttbar histograms with MC weights for different samples mixed on an event-by-event basis
  if (plot_ttbar) {
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
  } else {
    // get combined WW histograms with different final states
    // total weight above 600 GeV (emu_mass): w_t = 1/(1/w_central + (N_p1+N_p2)/(w_p1*N_p1 + w_p2*N_p2))
    cutVars.clear();
    lowCuts.clear();
    highCuts.clear();
    mcWeightsForCutRanges.clear();
    cutVars.push_back("");
    lowCuts.push_back(0.);
    highCuts.push_back(1.e9);
    mcWeightsForCutRanges.push_back(wwMcWeight->GetVal());
    cutVars.push_back("emu_mass");
    lowCuts.push_back(600.);
    highCuts.push_back(1.e9);
    mcWeightsForCutRanges.push_back((wwMcWeight600upEminusMuPlus->GetVal()*wwNGen600upEminusMuPlus->GetVal() + wwMcWeight600upEplusMuMinus->GetVal()*wwNGen600upEplusMuMinus->GetVal())/(wwNGen600upEminusMuPlus->GetVal()+wwNGen600upEplusMuMinus->GetVal()));
    emuTest_ww.push_back(MakeHistoFromBranch(&input, "emuTree_ww", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, flags, normToBin));
    emuTest_wwPriv600upEminusMuPlus.push_back(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEminusMuPlus", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, flags, normToBin));
    emuTest_wwPriv600upEplusMuMinus.push_back(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEplusMuMinus", "", testVar, sig, histoReg, cutVars, lowCuts, highCuts, mcWeightsForCutRanges, binning, flags, normToBin));
  }

  sig+=3;

  // add overflow to last bin
  if (overflowBin) {
    refOverTestHist->SetBinContent(refOverTestHist->GetNbinsX(), refOverTestHist->GetBinContent(refOverTestHist->GetNbinsX()) + refOverTestHist->GetBinContent(refOverTestHist->GetNbinsX() + 1));
    emuTest_ref.back()->SetBinContent(emuTest_ref.back()->GetNbinsX(), emuTest_ref.back()->GetBinContent(emuTest_ref.back()->GetNbinsX()) + emuTest_ref.back()->GetBinContent(emuTest_ref.back()->GetNbinsX() + 1));
    refOverTestHist->SetBinError(refOverTestHist->GetNbinsX(), sqrt(refOverTestHist->GetBinContent(refOverTestHist->GetNbinsX())));
    emuTest_ref.back()->SetBinError(emuTest_ref.back()->GetNbinsX(), sqrt(emuTest_ref.back()->GetBinContent(emuTest_ref.back()->GetNbinsX())));

    if (plot_ttbar) {
      emuTest_ttbar.back()->SetBinContent(emuTest_ttbar.back()->GetNbinsX(), emuTest_ttbar.back()->GetBinContent(emuTest_ttbar.back()->GetNbinsX()) + emuTest_ttbar.back()->GetBinContent(emuTest_ttbar.back()->GetNbinsX() + 1));
      emuTest_ttbar700to1000.back()->SetBinContent(emuTest_ttbar700to1000.back()->GetNbinsX(), emuTest_ttbar700to1000.back()->GetBinContent(emuTest_ttbar700to1000.back()->GetNbinsX()) + emuTest_ttbar700to1000.back()->GetBinContent(emuTest_ttbar700to1000.back()->GetNbinsX() + 1));
      emuTest_ttbar1000up.back()->SetBinContent(emuTest_ttbar1000up.back()->GetNbinsX(), emuTest_ttbar1000up.back()->GetBinContent(emuTest_ttbar1000up.back()->GetNbinsX()) + emuTest_ttbar1000up.back()->GetBinContent(emuTest_ttbar1000up.back()->GetNbinsX() + 1));
      emuTest_ttbarPriv600up.back()->SetBinContent(emuTest_ttbarPriv600up.back()->GetNbinsX(), emuTest_ttbarPriv600up.back()->GetBinContent(emuTest_ttbarPriv600up.back()->GetNbinsX()) + emuTest_ttbarPriv600up.back()->GetBinContent(emuTest_ttbarPriv600up.back()->GetNbinsX() + 1));
      emuTest_ttbar.back()->SetBinError(emuTest_ttbar.back()->GetNbinsX(), sqrt(emuTest_ttbar.back()->GetBinContent(emuTest_ttbar.back()->GetNbinsX())));
      emuTest_ttbar700to1000.back()->SetBinError(emuTest_ttbar700to1000.back()->GetNbinsX(), sqrt(emuTest_ttbar700to1000.back()->GetBinContent(emuTest_ttbar700to1000.back()->GetNbinsX())));
      emuTest_ttbar1000up.back()->SetBinError(emuTest_ttbar1000up.back()->GetNbinsX(), sqrt(emuTest_ttbar1000up.back()->GetBinContent(emuTest_ttbar1000up.back()->GetNbinsX())));
      emuTest_ttbarPriv600up.back()->SetBinError(emuTest_ttbarPriv600up.back()->GetNbinsX(), sqrt(emuTest_ttbarPriv600up.back()->GetBinContent(emuTest_ttbarPriv600up.back()->GetNbinsX())));
    } else {
      emuTest_ww.back()->SetBinContent(emuTest_ww.back()->GetNbinsX(), emuTest_ww.back()->GetBinContent(emuTest_ww.back()->GetNbinsX()) + emuTest_ww.back()->GetBinContent(emuTest_ww.back()->GetNbinsX() + 1));
      emuTest_wwPriv600upEminusMuPlus.back()->SetBinContent(emuTest_wwPriv600upEminusMuPlus.back()->GetNbinsX(), emuTest_wwPriv600upEminusMuPlus.back()->GetBinContent(emuTest_wwPriv600upEminusMuPlus.back()->GetNbinsX()) + emuTest_wwPriv600upEminusMuPlus.back()->GetBinContent(emuTest_wwPriv600upEminusMuPlus.back()->GetNbinsX() + 1));
      emuTest_wwPriv600upEplusMuMinus.back()->SetBinContent(emuTest_wwPriv600upEplusMuMinus.back()->GetNbinsX(), emuTest_wwPriv600upEplusMuMinus.back()->GetBinContent(emuTest_wwPriv600upEplusMuMinus.back()->GetNbinsX()) + emuTest_wwPriv600upEplusMuMinus.back()->GetBinContent(emuTest_wwPriv600upEplusMuMinus.back()->GetNbinsX() + 1));
      emuTest_ww.back()->SetBinError(emuTest_ww.back()->GetNbinsX(), sqrt(emuTest_ww.back()->GetBinContent(emuTest_ww.back()->GetNbinsX())));
      emuTest_wwPriv600upEminusMuPlus.back()->SetBinError(emuTest_wwPriv600upEminusMuPlus.back()->GetNbinsX(), sqrt(emuTest_wwPriv600upEminusMuPlus.back()->GetBinContent(emuTest_wwPriv600upEminusMuPlus.back()->GetNbinsX())));
      emuTest_wwPriv600upEplusMuMinus.back()->SetBinError(emuTest_wwPriv600upEplusMuMinus.back()->GetNbinsX(), sqrt(emuTest_wwPriv600upEplusMuMinus.back()->GetBinContent(emuTest_wwPriv600upEplusMuMinus.back()->GetNbinsX())));
    }
  }
  // add underflow to first bin
  if (underflowBin) {
    refOverTestHist->SetBinContent(1, refOverTestHist->GetBinContent(1) + refOverTestHist->GetBinContent(0));
    emuTest_ref.back()->SetBinContent(1, emuTest_ref.back()->GetBinContent(1) + emuTest_ref.back()->GetBinContent(0));
    refOverTestHist->SetBinError(1, sqrt(refOverTestHist->GetBinContent(1)));
    emuTest_ref.back()->SetBinError(1, sqrt(emuTest_ref.back()->GetBinContent(1)));
    if (plot_ttbar) {
      emuTest_ttbar.back()->SetBinContent(1, emuTest_ttbar.back()->GetBinContent(1) + emuTest_ttbar.back()->GetBinContent(0));
      emuTest_ttbar700to1000.back()->SetBinContent(1, emuTest_ttbar700to1000.back()->GetBinContent(1) + emuTest_ttbar700to1000.back()->GetBinContent(0));
      emuTest_ttbar1000up.back()->SetBinContent(1, emuTest_ttbar1000up.back()->GetBinContent(1) + emuTest_ttbar1000up.back()->GetBinContent(0));
      emuTest_ttbarPriv600up.back()->SetBinContent(1, emuTest_ttbarPriv600up.back()->GetBinContent(1) + emuTest_ttbarPriv600up.back()->GetBinContent(0));
      emuTest_ttbar.back()->SetBinError(1, sqrt(emuTest_ttbar.back()->GetBinContent(1)));
      emuTest_ttbar700to1000.back()->SetBinError(1, sqrt(emuTest_ttbar700to1000.back()->GetBinContent(1)));
      emuTest_ttbar1000up.back()->SetBinError(1, sqrt(emuTest_ttbar1000up.back()->GetBinContent(1)));
      emuTest_ttbarPriv600up.back()->SetBinError(1, sqrt(emuTest_ttbarPriv600up.back()->GetBinContent(1)));
    } else {
      emuTest_ww.back()->SetBinContent(1, emuTest_ww.back()->GetBinContent(1) + emuTest_ww.back()->GetBinContent(0));
      emuTest_wwPriv600upEminusMuPlus.back()->SetBinContent(1, emuTest_wwPriv600upEminusMuPlus.back()->GetBinContent(1) + emuTest_wwPriv600upEminusMuPlus.back()->GetBinContent(0));
      emuTest_wwPriv600upEplusMuMinus.back()->SetBinContent(1, emuTest_wwPriv600upEplusMuMinus.back()->GetBinContent(1) + emuTest_wwPriv600upEplusMuMinus.back()->GetBinContent(0));
      emuTest_ww.back()->SetBinError(1, sqrt(emuTest_ww.back()->GetBinContent(1)));
      emuTest_wwPriv600upEminusMuPlus.back()->SetBinError(1, sqrt(emuTest_wwPriv600upEminusMuPlus.back()->GetBinContent(1)));
      emuTest_wwPriv600upEplusMuMinus.back()->SetBinError(1, sqrt(emuTest_wwPriv600upEplusMuMinus.back()->GetBinContent(1)));
    }
  }

  // make a histogram stack with the test 
  THStack *testStack = new THStack("testStack" + sign[sig] + region[reg], testVar + sign[sig] + nameReg[reg]);
  if (plot_ttbar) {
    testStack->Add(emuTest_ttbar.back());
    testStack->Add(emuTest_ttbar700to1000.back());
    testStack->Add(emuTest_ttbar1000up.back());
    testStack->Add(emuTest_ttbarPriv600up.back());
  } else {
    testStack->Add(emuTest_ww.back());
    testStack->Add(emuTest_wwPriv600upEminusMuPlus.back());
    testStack->Add(emuTest_wwPriv600upEplusMuMinus.back());
  }

  // test histogram for (ref - test) / test
  TH1F *testHist;
  if (plot_ttbar) {
    testHist = (TH1F *)emuTest_ttbar.back()->Clone("testHist_ttbar");
    testHist->Add(emuTest_ttbar700to1000.back());
    testHist->Add(emuTest_ttbar1000up.back());
    testHist->Add(emuTest_ttbarPriv600up.back());
  } else {
    testHist = (TH1F *)emuTest_ww.back()->Clone("testHist_ww");
    testHist->Add(emuTest_wwPriv600upEminusMuPlus.back());
    testHist->Add(emuTest_wwPriv600upEplusMuMinus.back());
  }

  TCanvas *emuPlot;
  TPad *specPad;
  if (plotPull) {
    emuPlot = new TCanvas("emuPlot" + testVar + sign[sig] + region[reg], "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg], 100, 100, 900, 740);
    specPad = new TPad("specPad" + testVar + + sign[sig] + region[reg], "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg], 0., 0.33, 1., 1.);
    specPad->SetBottomMargin(0.06);
  } else {
    emuPlot = new TCanvas("emuPlot" + testVar + sign[sig] + region[reg], "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg], 100, 100, 900, 600);
    specPad = new TPad("specPad" + testVar + sign[sig] + region[reg], "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg], 0., 0., 1., 1.);
    specPad->SetBottomMargin(0.12);
  }
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

  // make sure that data and bkg are visible on plot
  if (emuTest_ref.back()->GetMaximum() > testHist->GetMaximum()) {
    if (!logPlot) testStack->SetMaximum(emuTest_ref.back()->GetMaximum() * 1.1);
    else testStack->SetMaximum(emuTest_ref.back()->GetMaximum() * 1.3);
  }

  //// plot spectrum
  if (plot_ttbar) {
    emuTest_ttbar.back()->SetFillColor(ttbarColour);
    emuTest_ttbar.back()->SetLineColor(ttbarColour);
    //emuTest_ttbar.back()->SetLineWidth(2);
    emuTest_ttbar700to1000.back()->SetFillColor(ttbar700to1000Colour);
    emuTest_ttbar700to1000.back()->SetLineColor(ttbar700to1000Colour);
    //emuTest_ttbar700to1000.back()->SetLineWidth(2);
    emuTest_ttbar1000up.back()->SetFillColor(ttbar1000upColour);
    emuTest_ttbar1000up.back()->SetLineColor(ttbar1000upColour);
    //emuTest_ttbar1000up.back()->SetLineWidth(2);
    emuTest_ttbarPriv600up.back()->SetFillColor(ttbarPriv600upColour);
    emuTest_ttbarPriv600up.back()->SetLineColor(ttbarPriv600upColour);
    //emuTest_ttbarPriv600up.back()->SetLineWidth(2);
  } else {
    emuTest_ww.back()->SetFillColor(wwColour);
    emuTest_ww.back()->SetLineColor(wwColour);
    //emuTest_ww.back()->SetLineWidth(2);
    emuTest_wwPriv600upEminusMuPlus.back()->SetFillColor(wwPriv600upEminusMuPlusColour);
    emuTest_wwPriv600upEminusMuPlus.back()->SetLineColor(wwPriv600upEminusMuPlusColour);
    //emuTest_wwPriv600upEminusMuPlus.back()->SetLineWidth(2);
    emuTest_wwPriv600upEplusMuMinus.back()->SetFillColor(wwPriv600upEplusMuMinusColour);
    emuTest_wwPriv600upEplusMuMinus.back()->SetLineColor(wwPriv600upEplusMuMinusColour);
    //emuTest_wwPriv600upEplusMuMinus.back()->SetLineWidth(2);
  }

  testStack->Draw("hist");
  if (!plotPull) testStack->GetXaxis()->SetTitle(testPlots[var].fXaxisTitle);
  testStack->GetXaxis()->SetTitleFont(font);
  testStack->GetXaxis()->SetTitleSize(0.047);
  testStack->GetXaxis()->SetLabelSize(0.05);
  testStack->GetXaxis()->SetLabelFont(font);
  testStack->GetXaxis()->SetMoreLogLabels();
  testStack->GetXaxis()->SetNoExponent();
  testStack->GetYaxis()->SetTitle("Events");
  testStack->GetYaxis()->SetTitleFont(font);
  testStack->GetYaxis()->SetTitleSize(0.047);
  testStack->GetYaxis()->SetTitleOffset(1.1);
  testStack->GetYaxis()->SetLabelFont(font);
  testStack->GetYaxis()->SetLabelSize(0.05);
  if (!logPlot) testStack->SetMinimum(0.);
  else testStack->SetMinimum(0.5);
  testStack->Draw("hist");

  emuTest_ref.back()->SetLineWidth(2);
  emuTest_ref.back()->SetLineColor(kBlack);
  emuTest_ref.back()->SetMarkerStyle(8);
  emuTest_ref.back()->SetMarkerSize(0.8);
  emuTest_ref.back()->Draw("same");

  // legent and labels
  TLegend legend(0.56, 0.53, 0.75, 0.885);
  legend.SetTextFont(font);
  legend.SetTextSize(0.042);
  legend.SetBorderSize(0);
  legend.SetLineColor(1);
  legend.SetLineStyle(1);
  legend.SetLineWidth(1);
  legend.SetFillColor(19);
  legend.SetFillStyle(0);
  if (plot_ttbar) {
    legend.AddEntry(emuTest_ref.back(), "t#bar{t} (reference)");
    legend.AddEntry(emuTest_ttbar.back(), "t#bar{t}_{rew}", "F");
    legend.AddEntry(emuTest_ttbar700to1000.back(), "t#bar{t}_{rew} 700 GeV < M_{t#bar{t}} < 1000 GeV", "F");
    legend.AddEntry(emuTest_ttbar1000up.back(), "t#bar{t}_{rew} M_{t#bar{t}} > 1000 GeV", "F");
    legend.AddEntry(emuTest_ttbarPriv600up.back(), "t#bar{t}_{rew} (Priv.)", "F");
  } else {
    legend.AddEntry(emuTest_ref.back(), "WW (reference)");
    legend.AddEntry(emuTest_ww.back(), "WW_{rew}", "F");
    legend.AddEntry(emuTest_wwPriv600upEminusMuPlus.back(), "WW #rightarrow e^{-}#mu^{+}_{rew} (priv.)", "F");
    legend.AddEntry(emuTest_wwPriv600upEplusMuMinus.back(), "WW #rightarrow e^{+}#mu^{-}_{rew} (priv.)", "F");
  }
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

  // plot a (ref - test) / test histogram below the spectrum
  if (plotPull) {
    refOverTestHist->Add(testHist, -1.);
    refOverTestHist->Divide(testHist);

    float fontScaleBot = 1.;
    TPad *pullPad = new TPad("pullPad" + testVar + + sign[sig] + region[reg], "(ref - test) / test" + testVar + nameSign[sig] + nameReg[reg], 0., 0., 1., 0.33);
    emuPlot->cd();
    fontScaleBot = specPad->GetHNDC() / pullPad->GetHNDC();
    pullPad->Draw();
    pullPad->cd();
    pullPad->SetBorderMode(0);
    pullPad->SetBorderSize(2);
    pullPad->SetFrameBorderMode(0);
    pullPad->SetFillColor(0);
    pullPad->SetFrameFillColor(0);
    pullPad->SetTopMargin(0.);
    pullPad->SetBottomMargin(0.22);
    pullPad->SetLeftMargin(0.11);
    pullPad->SetRightMargin(0.09);
    pullPad->SetTickx(1);
    pullPad->SetTicky(1);
    if (pullGridY) pullPad->SetGridy();

    refOverTestHist->SetLineWidth(1);
    refOverTestHist->SetLineColor(kBlack);
    refOverTestHist->SetMarkerStyle(20);
    refOverTestHist->SetMarkerSize(1.1);

    refOverTestHist->GetXaxis()->SetTitle(testPlots[var].fXaxisTitle);
    refOverTestHist->GetXaxis()->SetTitleFont(font);
    refOverTestHist->GetXaxis()->SetTitleSize(0.047 * fontScaleBot);
    refOverTestHist->GetXaxis()->SetTitleOffset(1.);
    refOverTestHist->GetXaxis()->SetLabelFont(font);
    refOverTestHist->GetXaxis()->SetLabelSize(0.05 * fontScaleBot);
    refOverTestHist->GetXaxis()->SetMoreLogLabels();
    refOverTestHist->GetXaxis()->SetNoExponent();
    //refOverTestHist->GetXaxis()->SetRangeUser(xRangeMinRatio, xRangeMaxRatio);
    refOverTestHist->GetYaxis()->SetTitle("(ref-test)/test");
    refOverTestHist->GetYaxis()->SetTitleFont(font);
    refOverTestHist->GetYaxis()->SetTitleSize(0.047 * fontScaleBot);
    refOverTestHist->GetYaxis()->SetTitleOffset(1.1 / fontScaleBot);
    refOverTestHist->GetYaxis()->SetLabelFont(font);
    refOverTestHist->GetYaxis()->SetLabelSize(0.05 * fontScaleBot);
    refOverTestHist->GetYaxis()->SetRangeUser(yRangeMinRatio, yRangeMaxRatio);

    refOverTestHist->Draw();
  }

  // safe in various file formats
  if (saveSpec) {
    sStream.str("");
    sStream << plotDir << "emuControlSpec_" << testVar << sign[sig] << region[reg] << "_" << fileNameExtra << lumi->GetVal() << "pb-1";
    TString saveFileName = sStream.str();
    if (saveAsPdf) emuPlot->Print(saveFileName + ".pdf", "pdf");
    if (saveAsPng) emuPlot->Print(saveFileName + ".png", "png");
    if (saveAsRoot) emuPlot->Print(saveFileName + ".root", "root");
  }

  cout << testVar.Data() << nameSign[sig].Data() << nameReg[reg].Data() << " plotted" << endl;
}

// plot a range of control variables
void PlotRange(int sign = 0, unsigned int region = 0, unsigned int from = 1, unsigned int to = 43)
{
  if (to == 0 || to > 43) to = 43;
  if (from == 0 || from > to) from = 1;
  for (unsigned int i = from; i <= to; ++i)
    macro_MakeMultiSampleClosureTest(i, sign, region);
}

