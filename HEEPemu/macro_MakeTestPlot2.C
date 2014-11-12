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
#include "TGraphAsymmErrors.h"
#include "TColor.h"
#include "TString.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THashList.h"
#include "Math/QuantFuncMathCore.h"

#include "makeHistoFromBranch.C"

TGraphAsymmErrors * makeDataGraph(TH1* dataHist,float normToBinWidth,bool xErrBars);

// class for plotting of control variable plots
class ContVarPlot {
  public:
    TString fName;
    TString fTitle;
    TString fXaxisTitle;
    float fXmin;
    float fXmax;
    int fNbins;
    bool fPlotQcd;
    bool fLogPlot;
    bool fUnderFlow;
    bool fOverFlow;

    ContVarPlot(const char* name, const char* title, const char* xAxisTitle, float xmin, float xmax, int nBins, bool plotQcd=1, bool logPlot=0, bool underFlow=0, bool overFlow=0) : fName(name), fTitle(title), fXaxisTitle(xAxisTitle), fXmin(xmin), fXmax(xmax), fNbins(nBins), fPlotQcd(plotQcd), fLogPlot(logPlot), fUnderFlow(underFlow), fOverFlow(overFlow) { } 
    ~ContVarPlot() { }
};

// Plot e-mu+ minus e+mu- histogram
///////////////////////////////////////////////////////////////////////
void macro_MakeTestPlot2(unsigned int var = 0, int sig = 0, unsigned int reg = 0)
{ 
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("./emuSpec_MuGammaTrg_topxsect245p8_19703pb-1.root", "open");
  //TFile input("./emuSpec_MuGammaTrg_topxsect252p89_19703pb-1.root", "open");
  //TFile input("./emuSpec_singleMuTrg_topxsect245p8_19706pb-1.root", "open");
  TFile input("./emuSpec_singleMuTrg_altdiboson_19706pb-1.root", "open");
  //TFile input("./emuSpec_singleMuTrg_topxsect252p89_19706pb-1.root", "open");
  input.cd();

  TParameter<float> *lumi = (TParameter<float> *)input.Get("lumi");

  const bool topReweighting = 0;
  const bool ttbar_sample_type = 0; // 0=powheg, 1=madgraph
  const bool dy_sample_type = 1; // 0=dytoee/mumu/tautau, 1=dyjetstoll
  TString ww_base = "wwpow";
  TString wz_base = "wzmg";
  TString zz_base = "zzmg";
  TString tw_base = "twpow";
  const int qcdEst = 2; // estimation method of QCD contribution. none(0), from SS spectrum(1), from fake rate(2)
  const float plotSignal[4] = {10., 10., 10., 10.}; // signal scale factors. 0 for off
  const bool plotPull = 0; // plot (data-bkg)/bkg
  const bool plotShapeUnc = 0; // plot up/down shape uncertainty histograms
  const bool pullGridY = 1; // grid lines on y axis of pull plot
  const bool prelim = 1; // print Preliminary
  const bool writeHistosToFile = 1;

  float yRangeMinRatio = -0.7;
  float yRangeMaxRatio = 0.7;

  // plot style
  int ttbarColour = TColor::GetColor("#ff6666");
  int zttColour = TColor::GetColor("#ff4d4d");
  int zjetsllColour = TColor::GetColor("#99ccff");
  int wwColour = TColor::GetColor("#ff3333");
  int wzColour = TColor::GetColor("#ff0f0f");
  int zzColour = TColor::GetColor("#eb0000");
  int twColour = TColor::GetColor("#db0000");
  int zmmColour=  TColor::GetColor("#66b3ff");
  int zeeColour=  TColor::GetColor("#99ccff");
  int wjetColour=  TColor::GetColor("#ffd324");
  int jetBkgColour = TColor::GetColor("#ffff66");

  // output file formats
  const bool saveSpec = 0;
  const bool saveAsPdf = 1;
  const bool saveAsPng = 1;
  const bool saveAsRoot = 0;
  const char *fileNameExtra = "";
  //const char *fileNameExtra = "madgraphTTbar_";
  const char *plotDir = "./testplots/";
  TString outfileName = "test_plots";

  int font = 42; //62
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  std::vector<ContVarPlot> testPlots;
  testPlots.reserve(40);
  // flags: plotQCD | logPlot | underflow in first bin | overflow in last bin
  testPlots.push_back(ContVarPlot("mass", "e#mu invariant mass", "m(e#mu) [GeV]", 0., 1500., 75, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("pfMet", "PF MET", "E^{T}_{miss} [GeV]", 0., 500., 50, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("nVtx", "Number of primary vertices", "# PV", 0., 50., 50, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("rho", "rho", "#rho", 0., 50., 50, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("numOfJetsPt20", "Number of jets > 20 GeV", "# jets_{p_{T}>20}", 0., 20., 20, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("numOfJetsPt30", "Number of jets > 30 GeV", "# jets_{p_{T}>30}", 0., 15., 15, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("numOfBJetsPt30", "Number of b-jets > 30 GeV", "# b-jets_{p_{T}>30}", 0., 6., 6, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("numOfBJetsMVAPt30", "Number of b-jets_{MVA} > 30 GeV", "# b-jets_{p_{T}>30}", 0., 6., 6, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("eDxyMinusMuDxy", "|Dxy_e - Dxy_mu|", "|#Deltaxy_{e}-#Deltaxy_{#mu}|", 0., 0.06, 40, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eDzMinusMuDz", "|Dz_e - Dz_mu|", "|#Deltaz_{e}-#Deltaz_{#mu}|", 0., 0.5, 40, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("dEta", "|eta_e - eta_mu|", "|#eta_{e}-#eta_{#mu}|", 0., 5., 50, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("dPhi", "|phi_e - phi_mu|", "|#varphi_{e}-#varphi_{#mu}|", 0., 3.2, 64, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("eleEt", "Electron Et", "E_{T}^{e} [GeV]", 0., 500., 50, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleEta", "Electron eta", "#eta_{e}", -2.5, 2.5, 50, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("elePhi", "Electron phi", "#phi_{e}", -3.2, 3.2, 64, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("eleDEta", "Electron dEta", "#Delta#eta", -0.007, 0.007, 29, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleDPhi", "Electron dPhi", "#Delta#varphi", -0.06, 0.06, 50, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleHOE", "Electron H/E", "H/E", 0., 0.051, 51, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleE1x5overE5x5", "Electron E1x5/E5x5", "E1x5/E5x5", 0., 1.5, 30, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleE2x5overE5x5", "Electron E2x5/E5x5", "E2x5/E5x5", 0., 1.5, 30, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleSigmaIEIE", "Electron #sigma_i#etai#eta", "", 0., 0.04, 40, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleEcalIso", "Electron ECAL iso", "ECAL iso", 0., 15., 30, 1, 0, 0, 0));
  testPlots.push_back(ContVarPlot("eleHcalIso1", "Electron HCAL iso1", "HCAL iso1", 0., 15., 30, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleHcalIso2", "Electron HCAL iso2", "HCAL iso2", 0., 15., 30, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleHeepIso", "Electron HEEP iso", "HEEP iso", 0., 20., 40, 1, 1, 0, 1));
  testPlots.push_back(ContVarPlot("eleTrkIso", "Electron track iso", "trk_{iso}", 0., 5.1, 51, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("eleLostHits", "Electron number of lost hits", "", 0., 2., 2, 1, 0, 0, 0));
  testPlots.push_back(ContVarPlot("eleDXYFstPVtx", "Electron Dxy at first PV", "#Deltaxy_{e}", -0.05, 0.05, 40, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eleDZFstPVtx", "Electron Dz at first PV", "#Deltaz_{e}", -1., 1., 80, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("etEleOPtMu", "Et_e / pt_mu", "E^{e}_{T} / p^{#mu}_{T}", 0., 13., 26, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("lepPtPlusOPtMinus", "lepton pt+ / pt-", "p^{+}_{T}/p^{-}_{T}", 0., 10., 20, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("eChTimesMuCh", "charge_e * charge_mu", "charge_{e}*charge_{#mu}", -1., 2., 3, 1, 1, 0, 0));
  testPlots.push_back(ContVarPlot("muPt", "Muon pt", "p_{T}^{#mu} [GeV]", 0., 500., 50, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muPtErr", "Muon pt error", "p_{T}^{#mu} [GeV]", 0., 50., 50, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muEta", "Muon eta", "#eta_{#mu}", -2.5, 2.5, 50, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muPhi", "Muon phi", "#varphi_{#mu}", -3.2, 3.2, 64, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muHitLayers", "Muon layers with hits", "# layers hits", 4., 20., 16, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muTrkHits", "Muon tracker hits", "# tracker hits", 0., 40., 40, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muPxlHits", "Muon pixel hits", "# pixel hits", 0., 10., 10, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muMuHits", "Muon system hits", "# #mu system hits", 0., 55., 55, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muDXYFstPVtx", "Muon Dxy at first PV", "#Deltaxy_{#mu}", -0.2, 0.2, 80, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muDZFstPVtx", "Muon Dz at first PV", "#Deltaz_{#mu}", -0.1, 0.1, 40, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muNSeg", "Muon number of segments", "# segments", 0., 7., 7, 1, 0, 1, 1));
  testPlots.push_back(ContVarPlot("muTrkIso03", "Muon track iso 03", "#mu trk iso", 0., 20., 40, 1, 1, 1, 1));
  testPlots.push_back(ContVarPlot("muIsoCombRel", "Muon combined iso / pt", "#mu (iso_{em}+iso_{had}+iso_{trk})/p_{T}", 0., 3.5, 35, 1, 1, 1, 1));

  TString sign[8] = {"_e-mu+_-_e+mu-", "_e-mu+", "_e+mu-", "_OS", "", "_SS", "_++", "_--"};
  TString nameSign[8] = {" e^{-}#mu^{+} - e^{+}#mu^{-}", " e-mu+", " e+mu-", " OS", "", " SS", " ++", " --"};
  TString region[3] = {"", "_EB", "_EE"};
  TString nameReg[3] = {"", " EB", " EE"};
  TString shapeUncNames[9] = {"", "eleScaleUp", "eleScaleDown", "muScaleUp", "muScaleDown", "muonResUp", "muonResDown", "muonResSmearUp", "muonResSmearDown"};

  // plot a list with possible test histograms
  if (var > testPlots.size()) var = 0;
  if (var == 0) {
    cout << "Use macro_MakeTestPlot.C(x, y, z) with x {1-" << testPlots.size() << "} being the \n";
    cout << "number of the test histogram to plot, y {-4 - +3} selects the charge with \n";
    cout << "the scheme e-mu: (e-mu+)-(e+mu-), -+, +-, OS, ALL, SS, ++, -- and z {0-2} selects the \n";
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
  if (abs(sig) > 4) sig = 0;
  if (reg > 2) reg = 0;
  // the makeHistoFromBranch function uses a different scheme for barrel and endcap selection
  unsigned int histoReg = 2;
  if (reg == 1) histoReg = 0;
  else if (reg == 2) histoReg = 1;

  std::vector<TH1F *> emuTest_data;
  std::vector<TH1F *> emuTest_ttbar;
  //std::vector<TH1F *> emuTest_ttbar700to1000;
  //std::vector<TH1F *> emuTest_ttbar1000up;
  //std::vector<TH1F *> emuTest_ttbarPriv600up;
  std::vector<TH1F *> emuTest_ztautau;
  std::vector<TH1F *> emuTest_zjetsll;
  std::vector<TH1F *> emuTest_ww;
  //std::vector<TH1F *> emuTest_wwPriv600upEminusMuPlus;
  //std::vector<TH1F *> emuTest_wwPriv600upEplusMuMinus;
  std::vector<TH1F *> emuTest_wz;
  std::vector<TH1F *> emuTest_zz;
  std::vector<TH1F *> emuTest_tw;
  std::vector<TH1F *> emuTest_zmumu;
  std::vector<TH1F *> emuTest_zee;
  std::vector<TH1F *> emuTest_wjets;
  std::vector<TH1F *> emuTest_qcd;
  std::vector<TH1F *> emuTest_sig1;
  std::vector<TH1F *> emuTest_sig2;
  std::vector<TH1F *> emuTest_sig3;
  std::vector<TH1F *> emuTest_sig4;

  std::vector<TH1F *> emuTest_data2;
  std::vector<TH1F *> emuTest_ttbar2;
  //std::vector<TH1F *> emuTest_ttbar700to10002;
  //std::vector<TH1F *> emuTest_ttbar1000up2;
  //std::vector<TH1F *> emuTest_ttbarPriv600up2;
  std::vector<TH1F *> emuTest_ztautau2;
  std::vector<TH1F *> emuTest_zjetsll2;
  std::vector<TH1F *> emuTest_ww2;
  //std::vector<TH1F *> emuTest_wwPriv600upEminusMuPlus2;
  //std::vector<TH1F *> emuTest_wwPriv600upEplusMuMinus2;
  std::vector<TH1F *> emuTest_wz2;
  std::vector<TH1F *> emuTest_zz2;
  std::vector<TH1F *> emuTest_tw2;
  std::vector<TH1F *> emuTest_zmumu2;
  std::vector<TH1F *> emuTest_zee2;
  std::vector<TH1F *> emuTest_wjets2;
  std::vector<TH1F *> emuTest_qcd2;
  std::vector<TH1F *> emuTest_sig12;
  std::vector<TH1F *> emuTest_sig22;
  std::vector<TH1F *> emuTest_sig32;
  std::vector<TH1F *> emuTest_sig42;

  outfileName = plotDir+outfileName+".root";
  TFile *outfile = new TFile(outfileName, "update");

  input.cd();

  // configure plot style
  TString testVar = testPlots[var].fName;
  bool plotQcd = testPlots[var].fPlotQcd;
  bool logPlot = testPlots[var].fLogPlot;
  if (abs(sig) > 3) logPlot = false;
  const bool overflowBin = testPlots[var].fOverFlow;
  const bool underflowBin = testPlots[var].fUnderFlow; 

  // make the correct binning
  std::vector<float> binning;
  float lowBin = testPlots[var].fXmin;
  float highBin = testPlots[var].fXmax;
  int nBins = testPlots[var].fNbins;
  for (float bin = lowBin; bin <= highBin; bin += (highBin - lowBin) / nBins)
    binning.push_back(bin);

  THashList *mcWeights = (THashList *)input.Get("mcWeights");
  THashList *nGenEvents = (THashList *)input.Get("nGenEvents");
  TParameter<float> *ttbarMcWeight = (TParameter<float> *)mcWeights->FindObject("ttbar");
  TParameter<float> *ttbarMcWeight700to1000 = (TParameter<float> *)mcWeights->FindObject("ttbar700to1000");
  TParameter<float> *ttbarMcWeight1000up = (TParameter<float> *)mcWeights->FindObject("ttbar1000up");
  TParameter<float> *ttbarMcWeightTo2l = (TParameter<float> *)mcWeights->FindObject("ttbarto2l");
  TParameter<float> *ttbarMcWeightTo1l1jet = (TParameter<float> *)mcWeights->FindObject("ttbarto1l1jet");
  TParameter<float> *ttbarMcWeight600up = (TParameter<float> *)mcWeights->FindObject("ttbarPriv600up");
  TParameter<float> *wwMcWeight = (TParameter<float> *)mcWeights->FindObject((const char*)ww_base);
  TParameter<float> *wwMcWeight600upEminusMuPlus = (TParameter<float> *)mcWeights->FindObject("wwPriv600upEminusMuPlus");
  TParameter<float> *wwMcWeight600upEplusMuMinus = (TParameter<float> *)mcWeights->FindObject("wwPriv600upEplusMuMinus");
  TParameter<float> *wwNGen600upEminusMuPlus = (TParameter<float> *)nGenEvents->FindObject("wwPriv600upEminusMuPlus");
  TParameter<float> *wwNGen600upEplusMuMinus = (TParameter<float> *)nGenEvents->FindObject("wwPriv600upEplusMuMinus");

  TFile* outFile;
  if (writeHistosToFile) outFile = new TFile("histograms.root", "recreate");

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
    lowCutsTtbarPriv600up.push_back(600.);
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

  // loop over shape uncertainty up and downscalings
  for (unsigned int shUnc = 0; shUnc < 9; ++shUnc) {
     if (shUnc > 0 && !plotShapeUnc) continue;
     TString shapeUncName = shapeUncNames[shUnc];
     input.cd();
     input.cd(shapeUncNames[shUnc]);
     if (shUnc > 0) {
        shapeUncName.Prepend("_");
        sig-=4;
     }

     // determine qcd contribution
     TH1F *ssData;
     TH1F *ssBg;
     TH1F *qcdContrib;
     TH1F *qcdContrib2;
     if (plotQcd && qcdEst == 1) {
       ssData = MakeHistoFromBranch(&input, "emuTree_data", shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x100);
       unsigned int flags = 0x1DF; // flags: [apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
       unsigned int topFlags = flags;
       if (topReweighting) topFlags += (1<<10);
       if (ttbar_sample_type) {
         ssBg = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo2l, binning, topFlags);
         ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbarto1l1jet", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo1l1jet, binning, topFlags));
         ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbarPriv600up, highCutsTtbarPriv600up, mcWeightsForCutRangesTtbarPriv600up, binning, topFlags));
       } else {
         ssBg = MakeHistoFromBranch(&input, "emuTree_ttbar", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags);
         ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
         ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
         ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", shapeUncName, testVar, 1, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
       }
       ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ztautau", shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, (const char*)("emuTree_"+ww_base), shapeUncName, testVar, 1, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEminusMuPlus", shapeUncName, testVar, 1, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEplusMuMinus", shapeUncName, testVar, 1, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, (const char*)("emuTree_"+wz_base), shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, (const char*)("emuTree_"+zz_base), shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, (const char*)("emuTree_"+tw_base), shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zmumu", shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zee", shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wjets", shapeUncName, testVar, 1, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       qcdContrib = (TH1F *)ssData->Clone("qcdContrib_SS");
       qcdContrib->Add(ssBg, -1);
       for (int i = 0; i < qcdContrib->GetNbinsX() + 2; ++i) {
         if (qcdContrib->GetBinContent(i) < 0) qcdContrib->SetBinContent(i, 0.);
       }
       cout << "Expected SS QCD events: " << ssData->Integral() - ssBg->Integral();
       cout << "; Derived SS QCD events: " << qcdContrib->Integral();
       cout << "; scale factor: " << (ssData->Integral() - ssBg->Integral()) / qcdContrib->Integral()<< endl;
       qcdContrib->Scale((ssData->Integral() - ssBg->Integral()) / qcdContrib->Integral());

       emuTest_qcd.push_back((TH1F *)qcdContrib->Clone(testVar + sign[sig+4] + nameReg[reg] + "qcd"));
       if (sig == 0) emuTest_qcd.back()->Scale(2.);
       else if (abs(sig) > 1) emuTest_qcd.back()->Scale(0.5);
     }

     // get the histograms
     TH1F *dataOverBgHist = (TH1F *)MakeHistoFromBranch(&input, "emuTree_data", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x100);
     if (shUnc == 0) emuTest_data.push_back((TH1F *)dataOverBgHist->Clone("mass_data_obs"));
     else emuTest_data.push_back((TH1F *)dataOverBgHist->Clone());
     unsigned int flags = 0x1DF; // flags: [apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
     unsigned int topFlags = flags;
     if (topReweighting) topFlags += (1<<10);
     TH1F *ttbarComb;
     if (ttbar_sample_type) {
       ttbarComb = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo2l, binning, topFlags);
       ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbarto1l1jet", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo1l1jet, binning, topFlags));
       ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbarPriv600up, highCutsTtbarPriv600up, mcWeightsForCutRangesTtbarPriv600up, binning, topFlags));
     } else {
       ttbarComb = MakeHistoFromBranch(&input, "emuTree_ttbar", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags);
       ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
       ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
       ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", shapeUncName, testVar, -3, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
     }
     emuTest_ttbar.push_back(ttbarComb);
     if (dy_sample_type == 0 )emuTest_ztautau.push_back(MakeHistoFromBranch(&input, "emuTree_ztautau", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     TH1F *wwComb = MakeHistoFromBranch(&input, (const char*)("emuTree_"+ww_base), shapeUncName, testVar, -3, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags);
     wwComb->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEminusMuPlus", shapeUncName, testVar, -3, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
     wwComb->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEplusMuMinus", shapeUncName, testVar, -3, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
     emuTest_ww.push_back(wwComb);
     emuTest_wz.push_back(MakeHistoFromBranch(&input, (const char*)("emuTree_"+wz_base), shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     emuTest_zz.push_back(MakeHistoFromBranch(&input, (const char*)("emuTree_"+zz_base), shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     emuTest_tw.push_back(MakeHistoFromBranch(&input, (const char*)("emuTree_"+tw_base), shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     if (dy_sample_type == 0 ) {
       emuTest_zmumu.push_back(MakeHistoFromBranch(&input, "emuTree_zmumu", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_zee.push_back(MakeHistoFromBranch(&input, "emuTree_zee", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     } else {
       emuTest_zjetsll.push_back(MakeHistoFromBranch(&input, "emuTree_zjetsll", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     }
     if (qcdEst != 2) emuTest_wjets.push_back(MakeHistoFromBranch(&input, "emuTree_wjets", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     if (plotSignal[0] > 0.) {
       emuTest_sig1.push_back(MakeHistoFromBranch(&input, "emuTree_sig500", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig1.back()->Scale(plotSignal[0]);
     }
     if (plotSignal[1] > 0.) {
       emuTest_sig2.push_back(MakeHistoFromBranch(&input, "emuTree_sig750", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig2.back()->Scale(plotSignal[1]);
     }
     if (plotSignal[2] > 0.) {
       emuTest_sig3.push_back(MakeHistoFromBranch(&input, "emuTree_sig1000", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig3.back()->Scale(plotSignal[2]);
     }
     if (plotSignal[3] > 0.) {
       emuTest_sig4.push_back(MakeHistoFromBranch(&input, "emuTree_sig1250", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig4.back()->Scale(plotSignal[3]);
     }

     emuTest_data2.push_back((TH1F *)MakeHistoFromBranch(&input, "emuTree_data", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x100));
     TH1F *ttbarComb2;
     if (ttbar_sample_type) {
       ttbarComb2 = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo2l, binning, topFlags);
       ttbarComb2->Add(MakeHistoFromBranch(&input, "emuTree_ttbarto1l1jet", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbarTo1l1jet, binning, topFlags));
       ttbarComb2->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbarPriv600up, highCutsTtbarPriv600up, mcWeightsForCutRangesTtbarPriv600up, binning, topFlags));
     } else {
       ttbarComb2 = MakeHistoFromBranch(&input, "emuTree_ttbar", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags);
       ttbarComb2->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
       ttbarComb2->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
       ttbarComb2->Add(MakeHistoFromBranch(&input, "emuTree_ttbarPriv600up", shapeUncName, testVar, -2, histoReg, cutVarsTtbar, lowCutsTtbar, highCutsTtbar, mcWeightsForCutRangesTtbar, binning, topFlags));
     }
     emuTest_ttbar2.push_back(ttbarComb2);
     if (dy_sample_type == 0 )emuTest_ztautau2.push_back(MakeHistoFromBranch(&input, "emuTree_ztautau", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     TH1F *wwComb2 = MakeHistoFromBranch(&input, (const char*)("emuTree_"+ww_base), shapeUncName, testVar, -2, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags);
     wwComb2->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEminusMuPlus", shapeUncName, testVar, -2, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
     wwComb2->Add(MakeHistoFromBranch(&input, "emuTree_wwPriv600upEplusMuMinus", shapeUncName, testVar, -2, histoReg, cutVarsWw, lowCutsWw, highCutsWw, mcWeightsForCutRangesWw, binning, flags));
     emuTest_ww2.push_back(wwComb2);
     emuTest_wz2.push_back(MakeHistoFromBranch(&input, (const char*)("emuTree_"+wz_base), shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     emuTest_zz2.push_back(MakeHistoFromBranch(&input, (const char*)("emuTree_"+zz_base), shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     emuTest_tw2.push_back(MakeHistoFromBranch(&input, (const char*)("emuTree_"+tw_base), shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     if (dy_sample_type == 0 ) {
       emuTest_zmumu2.push_back(MakeHistoFromBranch(&input, "emuTree_zmumu", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_zee2.push_back(MakeHistoFromBranch(&input, "emuTree_zee", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     } else {
       emuTest_zjetsll2.push_back(MakeHistoFromBranch(&input, "emuTree_zjetsll", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     }
     if (qcdEst != 2) emuTest_wjets2.push_back(MakeHistoFromBranch(&input, "emuTree_wjets", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
     if (plotSignal[0] > 0.) {
       emuTest_sig12.push_back(MakeHistoFromBranch(&input, "emuTree_sig500", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig12.back()->Scale(plotSignal[0]);
     }
     if (plotSignal[1] > 0.) {
       emuTest_sig22.push_back(MakeHistoFromBranch(&input, "emuTree_sig750", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig22.back()->Scale(plotSignal[1]);
     }
     if (plotSignal[2] > 0.) {
       emuTest_sig32.push_back(MakeHistoFromBranch(&input, "emuTree_sig1000", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig32.back()->Scale(plotSignal[2]);
     }
     if (plotSignal[3] > 0.) {
       emuTest_sig42.push_back(MakeHistoFromBranch(&input, "emuTree_sig1250", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, flags));
       emuTest_sig42.back()->Scale(plotSignal[3]);
     }

     // qcd contribution
     if (plotQcd && qcdEst == 2) {
       qcdContrib = MakeHistoFromBranch(&input, "frEmuTree_data", shapeUncName, testVar, -3, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x300);
       emuTest_qcd.push_back((TH1F *)qcdContrib->Clone(testVar + sign[1] + nameReg[reg] + "qcd"));
       qcdContrib2 = MakeHistoFromBranch(&input, "frEmuTree_data", shapeUncName, testVar, -2, histoReg, cutVarsEmpty, lowCutsEmpty, highCutsEmpty, mcWeightsForCutRangesEmpty, binning, 0x300);
       emuTest_qcd2.push_back((TH1F *)qcdContrib2->Clone(testVar + sign[2] + nameReg[reg] + "qcd"));
       // normalize to bin width
       //for (int i = 1; i < emuTest_qcd.back()->GetNbinsX() + 1; ++i) {
       //  emuTest_qcd.back()->SetBinContent(i, emuTest_qcd.back()->GetBinContent(i) / emuTest_qcd.back()->GetBinWidth(i));
       //  emuTest_qcd.back()->SetBinError(i, emuTest_qcd.back()->GetBinError(i) / emuTest_qcd.back()->GetBinWidth(i));
       //}
     }

     // write histograms in a file
     TH1F* emuTest_allBkg;
     TH1F* emuTest_allBkg2;
     if (writeHistosToFile) {
       TCanvas *fitCanvas;
       outFile->cd();
       emuTest_allBkg = (TH1F*)emuTest_ttbar.back()->Clone("emuTree_allBkg" + shapeUncName);
       emuTest_allBkg2 = (TH1F*)emuTest_ttbar2.back()->Clone("emuTree_allBkg2" + shapeUncName);
       if (dy_sample_type == 0) {
         emuTest_allBkg->Add(emuTest_ztautau.back());
         emuTest_allBkg2->Add(emuTest_ztautau2.back());
       }
       emuTest_allBkg->Add(emuTest_ww.back());
       emuTest_allBkg2->Add(emuTest_ww2.back());
       emuTest_allBkg->Add(emuTest_wz.back());
       emuTest_allBkg2->Add(emuTest_wz2.back());
       emuTest_allBkg->Add(emuTest_zz.back());
       emuTest_allBkg2->Add(emuTest_zz2.back());
       emuTest_allBkg->Add(emuTest_tw.back());
       emuTest_allBkg2->Add(emuTest_tw2.back());
       if (dy_sample_type == 0) {
         emuTest_allBkg->Add(emuTest_zmumu.back());
         emuTest_allBkg2->Add(emuTest_zmumu2.back());
         emuTest_allBkg->Add(emuTest_zee.back());
         emuTest_allBkg2->Add(emuTest_zee2.back());
       } else {
         emuTest_allBkg->Add(emuTest_zjetsll.back());
         emuTest_allBkg2->Add(emuTest_zjetsll2.back());
       }
       if (qcdEst != 2) {
         emuTest_allBkg->Add(emuTest_wjets.back());
         emuTest_allBkg2->Add(emuTest_wjets2.back());
       }
       if (plotQcd) {
         emuTest_allBkg->Add(emuTest_qcd.back());
         emuTest_allBkg2->Add(emuTest_qcd2.back());
       }
       emuTest_data.back()->Write();
       emuTest_ttbar.back()->Write();
       if (dy_sample_type == 0) emuTest_ztautau.back()->Write();
       emuTest_ww.back()->Write();
       emuTest_wz.back()->Write();
       emuTest_zz.back()->Write();
       emuTest_tw.back()->Write();
       if (dy_sample_type == 0) {
         emuTest_zmumu.back()->Write();
         emuTest_zee.back()->Write();
       } else {
         emuTest_zjetsll.back()->Write();
       }
       if (qcdEst != 2) emuTest_wjets.back()->Write();
       if (plotQcd) emuTest_qcd.back()->Write(testVar + "_qcd" + shapeUncName);
       emuTest_allBkg->Write(testVar + "_allBkg" + shapeUncName);
       if (var == 0 && abs(sig) < 4) {
         fitCanvas = new TCanvas("fitCanvas" + shapeUncName, "fitCanvas" + shapeUncName, 100, 100, 700, 600);
         fitCanvas->SetLogy();
         TF1 *bgParamFunc = new TF1("bgParamFunc" + shapeUncName, "1/[1]*(1+([2]*(x-[0]))/([1]))**(-1/[2]-1)", 0., 6000.);
         bgParamFunc->SetParLimits(0, 100., 1.e5);
         bgParamFunc->SetParLimits(1, 10., 1000.);
         bgParamFunc->SetParLimits(2, 0.01, 1.);
         bgParamFunc->SetParNames("m_{min}", "#alpha", "#beta");
         fitCanvas->cd();
         emuTest_allBkg->Fit("bgParamFunc" + shapeUncName, "", "", 150., 1500.);
         std::cout << "Fit chi^2, ndf, chi^2/ndf: " << bgParamFunc->GetChisquare() << ", " << bgParamFunc->GetNDF() << ", " << bgParamFunc->GetChisquare()/bgParamFunc->GetNDF() << std::endl;
         outFile->cd();
         bgParamFunc->Write();
       }
     }
     input.cd();

     TH1F *emuDiff_data = (TH1F *)emuTest_data.back()->Clone("emuDiff_data");
     TH1F *emuDiff_allBkg = (TH1F *)emuTest_allBkg->Clone("emuDiff_allBkg");
     TH1F *emuDiff_sig1 = (TH1F *)emuTest_sig1.back()->Clone("emuDiff_sig1");
     TH1F *emuDiff_sig2 = (TH1F *)emuTest_sig2.back()->Clone("emuDiff_sig2");
     TH1F *emuDiff_sig3 = (TH1F *)emuTest_sig3.back()->Clone("emuDiff_sig3");
     TH1F *emuDiff_sig4 = (TH1F *)emuTest_sig4.back()->Clone("emuDiff_sig4");
     emuDiff_data->Add(emuTest_data2.back(), -1);
     emuDiff_allBkg->Add(emuTest_allBkg2, -1);
     emuDiff_sig1->Add(emuTest_sig12.back(), -1);
     emuDiff_sig2->Add(emuTest_sig22.back(), -1);
     emuDiff_sig3->Add(emuTest_sig32.back(), -1);
     emuDiff_sig4->Add(emuTest_sig42.back(), -1);
     //emuDiff_data->Divide(emuTest_data2.back());
     //emuDiff_allBkg->Divide(emuTest_allBkg2);
     //emuDiff_sig1->Divide(emuTest_sig12.back());
     //emuDiff_sig2->Divide(emuTest_sig22.back());
     //emuDiff_sig3->Divide(emuTest_sig32.back());
     //emuDiff_sig4->Divide(emuTest_sig42.back());

     sig+=4;

     // add overflow to last bin
     if (overflowBin) {
       dataOverBgHist->SetBinContent(dataOverBgHist->GetNbinsX(), dataOverBgHist->GetBinContent(dataOverBgHist->GetNbinsX()) + dataOverBgHist->GetBinContent(dataOverBgHist->GetNbinsX() + 1));
       emuTest_data.back()->SetBinContent(emuTest_data.back()->GetNbinsX(), emuTest_data.back()->GetBinContent(emuTest_data.back()->GetNbinsX()) + emuTest_data.back()->GetBinContent(emuTest_data.back()->GetNbinsX() + 1));
       emuTest_ttbar.back()->SetBinContent(emuTest_ttbar.back()->GetNbinsX(), emuTest_ttbar.back()->GetBinContent(emuTest_ttbar.back()->GetNbinsX()) + emuTest_ttbar.back()->GetBinContent(emuTest_ttbar.back()->GetNbinsX() + 1));
       if (dy_sample_type == 0) emuTest_ztautau.back()->SetBinContent(emuTest_ztautau.back()->GetNbinsX(), emuTest_ztautau.back()->GetBinContent(emuTest_ztautau.back()->GetNbinsX()) + emuTest_ztautau.back()->GetBinContent(emuTest_ztautau.back()->GetNbinsX() + 1));
       emuTest_ww.back()->SetBinContent(emuTest_ww.back()->GetNbinsX(), emuTest_ww.back()->GetBinContent(emuTest_ww.back()->GetNbinsX()) + emuTest_ww.back()->GetBinContent(emuTest_ww.back()->GetNbinsX() + 1));
       emuTest_wz.back()->SetBinContent(emuTest_wz.back()->GetNbinsX(), emuTest_wz.back()->GetBinContent(emuTest_wz.back()->GetNbinsX()) + emuTest_wz.back()->GetBinContent(emuTest_wz.back()->GetNbinsX() + 1));
       emuTest_zz.back()->SetBinContent(emuTest_zz.back()->GetNbinsX(), emuTest_zz.back()->GetBinContent(emuTest_zz.back()->GetNbinsX()) + emuTest_zz.back()->GetBinContent(emuTest_zz.back()->GetNbinsX() + 1));
       emuTest_tw.back()->SetBinContent(emuTest_tw.back()->GetNbinsX(), emuTest_tw.back()->GetBinContent(emuTest_tw.back()->GetNbinsX()) + emuTest_tw.back()->GetBinContent(emuTest_tw.back()->GetNbinsX() + 1));
       if (dy_sample_type == 0) {
         emuTest_zmumu.back()->SetBinContent(emuTest_zmumu.back()->GetNbinsX(), emuTest_zmumu.back()->GetBinContent(emuTest_zmumu.back()->GetNbinsX()) + emuTest_zmumu.back()->GetBinContent(emuTest_zmumu.back()->GetNbinsX() + 1));
         emuTest_zee.back()->SetBinContent(emuTest_zee.back()->GetNbinsX(), emuTest_zee.back()->GetBinContent(emuTest_zee.back()->GetNbinsX()) + emuTest_zee.back()->GetBinContent(emuTest_zee.back()->GetNbinsX() + 1));
       } else {
         emuTest_zjetsll.back()->SetBinContent(emuTest_zjetsll.back()->GetNbinsX(), emuTest_zjetsll.back()->GetBinContent(emuTest_zjetsll.back()->GetNbinsX()) + emuTest_zjetsll.back()->GetBinContent(emuTest_zjetsll.back()->GetNbinsX() + 1));
       }
       if (qcdEst != 2) emuTest_wjets.back()->SetBinContent(emuTest_wjets.back()->GetNbinsX(), emuTest_wjets.back()->GetBinContent(emuTest_wjets.back()->GetNbinsX()) + emuTest_wjets.back()->GetBinContent(emuTest_wjets.back()->GetNbinsX() + 1));
       if (plotQcd) emuTest_qcd.back()->SetBinContent(emuTest_qcd.back()->GetNbinsX(), emuTest_qcd.back()->GetBinContent(emuTest_qcd.back()->GetNbinsX()) + emuTest_qcd.back()->GetBinContent(emuTest_qcd.back()->GetNbinsX() + 1));
       if (plotSignal[0] > 0.) emuTest_sig1.back()->SetBinContent(emuTest_sig1.back()->GetNbinsX(), emuTest_sig1.back()->GetBinContent(emuTest_sig1.back()->GetNbinsX()) + emuTest_sig1.back()->GetBinContent(emuTest_sig1.back()->GetNbinsX() + 1));
       if (plotSignal[1] > 0.) emuTest_sig2.back()->SetBinContent(emuTest_sig2.back()->GetNbinsX(), emuTest_sig2.back()->GetBinContent(emuTest_sig2.back()->GetNbinsX()) + emuTest_sig2.back()->GetBinContent(emuTest_sig2.back()->GetNbinsX() + 1));
       if (plotSignal[2] > 0.) emuTest_sig3.back()->SetBinContent(emuTest_sig3.back()->GetNbinsX(), emuTest_sig3.back()->GetBinContent(emuTest_sig3.back()->GetNbinsX()) + emuTest_sig3.back()->GetBinContent(emuTest_sig3.back()->GetNbinsX() + 1));
       if (plotSignal[3] > 0.) emuTest_sig4.back()->SetBinContent(emuTest_sig4.back()->GetNbinsX(), emuTest_sig4.back()->GetBinContent(emuTest_sig4.back()->GetNbinsX()) + emuTest_sig4.back()->GetBinContent(emuTest_sig4.back()->GetNbinsX() + 1));

       dataOverBgHist->SetBinError(dataOverBgHist->GetNbinsX(), sqrt(dataOverBgHist->GetBinContent(dataOverBgHist->GetNbinsX())));
       emuTest_data.back()->SetBinError(emuTest_data.back()->GetNbinsX(), sqrt(emuTest_data.back()->GetBinContent(emuTest_data.back()->GetNbinsX())));
       emuTest_ttbar.back()->SetBinError(emuTest_ttbar.back()->GetNbinsX(), sqrt(emuTest_ttbar.back()->GetBinContent(emuTest_ttbar.back()->GetNbinsX())));
       if (dy_sample_type == 0) emuTest_ztautau.back()->SetBinError(emuTest_ztautau.back()->GetNbinsX(), sqrt(emuTest_ztautau.back()->GetBinContent(emuTest_ztautau.back()->GetNbinsX())));
       emuTest_ww.back()->SetBinError(emuTest_ww.back()->GetNbinsX(), sqrt(emuTest_ww.back()->GetBinContent(emuTest_ww.back()->GetNbinsX())));
       emuTest_wz.back()->SetBinError(emuTest_wz.back()->GetNbinsX(), sqrt(emuTest_wz.back()->GetBinContent(emuTest_wz.back()->GetNbinsX())));
       emuTest_zz.back()->SetBinError(emuTest_zz.back()->GetNbinsX(), sqrt(emuTest_zz.back()->GetBinContent(emuTest_zz.back()->GetNbinsX())));
       emuTest_tw.back()->SetBinError(emuTest_tw.back()->GetNbinsX(), sqrt(emuTest_tw.back()->GetBinContent(emuTest_tw.back()->GetNbinsX())));
       if (dy_sample_type == 0) {
         emuTest_zmumu.back()->SetBinError(emuTest_zmumu.back()->GetNbinsX(), sqrt(emuTest_zmumu.back()->GetBinContent(emuTest_zmumu.back()->GetNbinsX())));
         emuTest_zee.back()->SetBinError(emuTest_zee.back()->GetNbinsX(), sqrt(emuTest_zee.back()->GetBinContent(emuTest_zee.back()->GetNbinsX())));
       } else {
         emuTest_zjetsll.back()->SetBinError(emuTest_zjetsll.back()->GetNbinsX(), sqrt(emuTest_zjetsll.back()->GetBinContent(emuTest_zjetsll.back()->GetNbinsX())));
       }
       if (qcdEst != 2) emuTest_wjets.back()->SetBinError(emuTest_wjets.back()->GetNbinsX(), sqrt(emuTest_wjets.back()->GetBinContent(emuTest_wjets.back()->GetNbinsX())));
       if (plotQcd) emuTest_qcd.back()->SetBinError(emuTest_qcd.back()->GetNbinsX(), sqrt(emuTest_qcd.back()->GetBinContent(emuTest_qcd.back()->GetNbinsX())));
       if (plotSignal[0]) emuTest_sig1.back()->SetBinError(emuTest_sig1.back()->GetNbinsX(), sqrt(emuTest_sig1.back()->GetBinContent(emuTest_sig1.back()->GetNbinsX())));
       if (plotSignal[1]) emuTest_sig2.back()->SetBinError(emuTest_sig2.back()->GetNbinsX(), sqrt(emuTest_sig2.back()->GetBinContent(emuTest_sig2.back()->GetNbinsX())));
       if (plotSignal[2]) emuTest_sig3.back()->SetBinError(emuTest_sig3.back()->GetNbinsX(), sqrt(emuTest_sig3.back()->GetBinContent(emuTest_sig3.back()->GetNbinsX())));
       if (plotSignal[3]) emuTest_sig4.back()->SetBinError(emuTest_sig4.back()->GetNbinsX(), sqrt(emuTest_sig4.back()->GetBinContent(emuTest_sig4.back()->GetNbinsX())));
     }
     // add underflow to first bin
     if (underflowBin) {
       dataOverBgHist->SetBinContent(1, dataOverBgHist->GetBinContent(1) + dataOverBgHist->GetBinContent(0));
       emuTest_data.back()->SetBinContent(1, emuTest_data.back()->GetBinContent(1) + emuTest_data.back()->GetBinContent(0));
       emuTest_ttbar.back()->SetBinContent(1, emuTest_ttbar.back()->GetBinContent(1) + emuTest_ttbar.back()->GetBinContent(0));
       if (dy_sample_type == 0) emuTest_ztautau.back()->SetBinContent(1, emuTest_ztautau.back()->GetBinContent(1) + emuTest_ztautau.back()->GetBinContent(0));
       emuTest_ww.back()->SetBinContent(1, emuTest_ww.back()->GetBinContent(1) + emuTest_ww.back()->GetBinContent(0));
       emuTest_wz.back()->SetBinContent(1, emuTest_wz.back()->GetBinContent(1) + emuTest_wz.back()->GetBinContent(0));
       emuTest_zz.back()->SetBinContent(1, emuTest_zz.back()->GetBinContent(1) + emuTest_zz.back()->GetBinContent(0));
       emuTest_tw.back()->SetBinContent(1, emuTest_tw.back()->GetBinContent(1) + emuTest_tw.back()->GetBinContent(0));
       if (dy_sample_type == 0) {
         emuTest_zmumu.back()->SetBinContent(1, emuTest_zmumu.back()->GetBinContent(1) + emuTest_zmumu.back()->GetBinContent(0));
         emuTest_zee.back()->SetBinContent(1, emuTest_zee.back()->GetBinContent(1) + emuTest_zee.back()->GetBinContent(0));
       } else {
         emuTest_zjetsll.back()->SetBinContent(1, emuTest_zjetsll.back()->GetBinContent(1) + emuTest_zjetsll.back()->GetBinContent(0));
       }
       if (qcdEst != 2) emuTest_wjets.back()->SetBinContent(1, emuTest_wjets.back()->GetBinContent(1) + emuTest_wjets.back()->GetBinContent(0));
       if (plotQcd) emuTest_qcd.back()->SetBinContent(1, emuTest_qcd.back()->GetBinContent(1) + emuTest_qcd.back()->GetBinContent(0));
       if (plotSignal[0]) emuTest_sig1.back()->SetBinContent(1, emuTest_sig1.back()->GetBinContent(1) + emuTest_sig1.back()->GetBinContent(0));
       if (plotSignal[1]) emuTest_sig2.back()->SetBinContent(1, emuTest_sig2.back()->GetBinContent(1) + emuTest_sig2.back()->GetBinContent(0));
       if (plotSignal[2]) emuTest_sig3.back()->SetBinContent(1, emuTest_sig3.back()->GetBinContent(1) + emuTest_sig3.back()->GetBinContent(0));
       if (plotSignal[3]) emuTest_sig4.back()->SetBinContent(1, emuTest_sig4.back()->GetBinContent(1) + emuTest_sig4.back()->GetBinContent(0));

       dataOverBgHist->SetBinError(1, sqrt(dataOverBgHist->GetBinContent(1)));
       emuTest_data.back()->SetBinError(1, sqrt(emuTest_data.back()->GetBinContent(1)));
       emuTest_ttbar.back()->SetBinError(1, sqrt(emuTest_ttbar.back()->GetBinContent(1)));
       if (dy_sample_type == 0) emuTest_ztautau.back()->SetBinError(1, sqrt(emuTest_ztautau.back()->GetBinContent(1)));
       emuTest_ww.back()->SetBinError(1, sqrt(emuTest_ww.back()->GetBinContent(1)));
       emuTest_wz.back()->SetBinError(1, sqrt(emuTest_wz.back()->GetBinContent(1)));
       emuTest_zz.back()->SetBinError(1, sqrt(emuTest_zz.back()->GetBinContent(1)));
       emuTest_tw.back()->SetBinError(1, sqrt(emuTest_tw.back()->GetBinContent(1)));
       if (dy_sample_type == 0) {
         emuTest_zmumu.back()->SetBinError(1, sqrt(emuTest_zmumu.back()->GetBinContent(1)));
         emuTest_zee.back()->SetBinError(1, sqrt(emuTest_zee.back()->GetBinContent(1)));
       } else {
         emuTest_zjetsll.back()->SetBinError(1, sqrt(emuTest_zjetsll.back()->GetBinContent(1)));
       }
       if (qcdEst != 2) emuTest_wjets.back()->SetBinError(1, sqrt(emuTest_wjets.back()->GetBinContent(1)));
       if (plotQcd) emuTest_qcd.back()->SetBinError(1, sqrt(emuTest_qcd.back()->GetBinContent(1)));
       if (plotSignal[0]) emuTest_sig1.back()->SetBinError(1, sqrt(emuTest_sig1.back()->GetBinContent(1)));
       if (plotSignal[1]) emuTest_sig2.back()->SetBinError(1, sqrt(emuTest_sig2.back()->GetBinContent(1)));
       if (plotSignal[2]) emuTest_sig3.back()->SetBinError(1, sqrt(emuTest_sig3.back()->GetBinContent(1)));
       if (plotSignal[3]) emuTest_sig4.back()->SetBinError(1, sqrt(emuTest_sig4.back()->GetBinContent(1)));
     }

     // make a histogram stack with the bg 
     THStack *bgStack = new THStack("bgStack" + sign[sig] + region[reg], testVar + sign[sig] + nameReg[reg]);
     if (plotQcd) bgStack->Add(emuTest_qcd.back());
     if (qcdEst != 2) bgStack->Add(emuTest_wjets.back());
     if (dy_sample_type == 0) {
       bgStack->Add(emuTest_zee.back());
       bgStack->Add(emuTest_zmumu.back());
     } else {
       bgStack->Add(emuTest_zjetsll.back());
     }
     bgStack->Add(emuTest_tw.back());
     bgStack->Add(emuTest_zz.back());
     bgStack->Add(emuTest_wz.back());
     bgStack->Add(emuTest_ww.back());
     if (dy_sample_type == 0) bgStack->Add(emuTest_ztautau.back());
     bgStack->Add(emuTest_ttbar.back());

     // bkg histogram for (data-bkg)/bkg plot
     TH1F *bgHist;
     if (plotQcd) {
       bgHist = (TH1F *)emuTest_qcd.back()->Clone();
       if (qcdEst != 2) bgHist->Add(emuTest_wjets.back());
       if (dy_sample_type == 0) bgHist->Add(emuTest_zee.back());
       else bgHist->Add(emuTest_zjetsll.back());
     } else {
       if (emuTest_wjets.size() > 0) {
         bgHist = (TH1F *)emuTest_wjets.back()->Clone();
         if (dy_sample_type == 0) bgHist->Add(emuTest_zee.back());
         else bgHist->Add(emuTest_zjetsll.back());
       } else {
         if (dy_sample_type == 0) bgHist = (TH1F *)emuTest_zee.back()->Clone();
         else bgHist = (TH1F *)emuTest_zjetsll.back()->Clone();
       }    
     }
     if (dy_sample_type == 0) bgHist->Add(emuTest_zmumu.back());
     bgHist->Add(emuTest_tw.back());
     bgHist->Add(emuTest_zz.back());
     bgHist->Add(emuTest_wz.back());
     bgHist->Add(emuTest_ww.back());
     if (dy_sample_type == 0) bgHist->Add(emuTest_ztautau.back());
     bgHist->Add(emuTest_ttbar.back());

     if (plotSignal[0] > 0.) {
       emuTest_sig1.back()->SetLineColor(kGreen);
       emuTest_sig1.back()->SetLineWidth(2);
     }
     if (plotSignal[1] > 0.) {
       emuTest_sig2.back()->SetLineColor(kBlue);
       emuTest_sig2.back()->SetLineWidth(2);
     }
     if (plotSignal[2] > 0.) {
       emuTest_sig3.back()->SetLineColor(kCyan);
       emuTest_sig3.back()->SetLineWidth(2);
     }
     if (plotSignal[3] > 0.) {
       emuTest_sig4.back()->SetLineColor(kMagenta);
       emuTest_sig4.back()->SetLineWidth(2);
     }
     //cout << "underflow (data/bkg): " << emuTest_data.back()->GetBinContent(0) << "/" << bgHist->GetBinContent(0) << "  overflow: " << emuTest_data.back()->GetBinContent(emuTest_data.back()->GetNbinsX()+1) << "/" << bgHist->GetBinContent(bgHist->GetNbinsX()+1) << endl;

     TCanvas *emuPlot;
     TPad *specPad;
     if (plotPull) {
       emuPlot = new TCanvas("emuPlot" + testVar + sign[sig] + region[reg] + shapeUncName, "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg] + shapeUncName, 100, 100, 900, 740);
       specPad = new TPad("specPad" + testVar + + sign[sig] + region[reg] + shapeUncName, "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg] + shapeUncName, 0., 0.33, 1., 1.);
       specPad->SetBottomMargin(0.06);
     } else {
       emuPlot = new TCanvas("emuPlot" + testVar + sign[sig] + region[reg] + shapeUncName, "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg] + shapeUncName, 100, 100, 900, 600);
       specPad = new TPad("specPad" + testVar + sign[sig] + region[reg] + shapeUncName, "emu Spectrum" + testVar + nameSign[sig] + nameReg[reg] + shapeUncName, 0., 0., 1., 1.);
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
     if (emuTest_data.back()->GetMaximum() > bgHist->GetMaximum()) {
       if (!logPlot) bgStack->SetMaximum(emuTest_data.back()->GetMaximum() * 1.1);
       else bgStack->SetMaximum(emuTest_data.back()->GetMaximum() * 1.3);
     }

     //// plot spectrum
     emuTest_ttbar.back()->SetFillColor(ttbarColour);
     emuTest_ttbar.back()->SetLineColor(kBlack);
     emuTest_ttbar.back()->SetLineWidth(2);
     if (dy_sample_type == 0) {
       emuTest_ztautau.back()->SetFillColor(zttColour);
       emuTest_ztautau.back()->SetMarkerColor(zttColour);
       emuTest_ztautau.back()->SetLineColor(kBlack);
       emuTest_ztautau.back()->SetLineWidth(2);
     } else {
       emuTest_zjetsll.back()->SetFillColor(zjetsllColour);
       emuTest_zjetsll.back()->SetMarkerColor(zjetsllColour);
       emuTest_zjetsll.back()->SetLineColor(kBlack);
       emuTest_zjetsll.back()->SetLineWidth(2);
     }
     emuTest_ww.back()->SetFillColor(wwColour);
     emuTest_ww.back()->SetMarkerColor(wwColour);
     emuTest_ww.back()->SetLineColor(kBlack);
     emuTest_ww.back()->SetLineWidth(2);
     emuTest_wz.back()->SetFillColor(wzColour);
     emuTest_wz.back()->SetMarkerColor(wzColour);
     emuTest_wz.back()->SetLineColor(kBlack);
     emuTest_wz.back()->SetLineWidth(2);
     emuTest_zz.back()->SetFillColor(zzColour);
     emuTest_zz.back()->SetMarkerColor(zzColour);
     emuTest_zz.back()->SetLineColor(kBlack);
     emuTest_zz.back()->SetLineWidth(2);
     emuTest_tw.back()->SetFillColor(twColour);
     emuTest_tw.back()->SetMarkerColor(twColour);
     emuTest_tw.back()->SetLineColor(kBlack);
     emuTest_tw.back()->SetLineWidth(2);
     if (dy_sample_type == 0) {
       emuTest_zmumu.back()->SetFillColor(zmmColour);
       emuTest_zmumu.back()->SetMarkerColor(zmmColour);
       emuTest_zmumu.back()->SetLineColor(kBlack);
       emuTest_zmumu.back()->SetLineWidth(2);
       emuTest_zee.back()->SetFillColor(zeeColour);
       emuTest_zee.back()->SetMarkerColor(zeeColour);
       emuTest_zee.back()->SetLineColor(kBlack);
       emuTest_zee.back()->SetLineWidth(2);
     }
     if (qcdEst != 2) {
       emuTest_wjets.back()->SetFillColor(wjetColour);
       emuTest_wjets.back()->SetMarkerColor(wjetColour);
       emuTest_wjets.back()->SetLineColor(kBlack);
       emuTest_wjets.back()->SetLineWidth(2);
     }
     if (plotQcd) {
       emuTest_qcd.back()->SetFillColor(jetBkgColour);
       emuTest_qcd.back()->SetMarkerColor(jetBkgColour);
       emuTest_qcd.back()->SetLineColor(kBlack);
       emuTest_qcd.back()->SetLineWidth(2);
     }

     emuDiff_allBkg->SetFillColor(ttbarColour);
     emuDiff_allBkg->SetLineColor(kBlack);
     emuDiff_allBkg->SetLineWidth(2);
     emuDiff_allBkg->Draw("hist");
     if (!plotPull) emuDiff_allBkg->GetXaxis()->SetTitle(testPlots[var].fXaxisTitle);
     emuDiff_allBkg->GetXaxis()->SetTitleFont(font);
     emuDiff_allBkg->GetXaxis()->SetTitleSize(0.047);
     emuDiff_allBkg->GetXaxis()->SetLabelSize(0.05);
     emuDiff_allBkg->GetXaxis()->SetLabelFont(font);
     emuDiff_allBkg->GetXaxis()->SetMoreLogLabels();
     emuDiff_allBkg->GetXaxis()->SetNoExponent();
     emuDiff_allBkg->GetYaxis()->SetTitle("Events");
     emuDiff_allBkg->GetYaxis()->SetTitleFont(font);
     emuDiff_allBkg->GetYaxis()->SetTitleSize(0.047);
     emuDiff_allBkg->GetYaxis()->SetTitleOffset(1.1);
     emuDiff_allBkg->GetYaxis()->SetLabelFont(font);
     emuDiff_allBkg->GetYaxis()->SetLabelSize(0.05);
     if (!logPlot) emuDiff_allBkg->SetMinimum(0.);
     else emuDiff_allBkg->SetMinimum(0.5);
     emuDiff_allBkg->Draw("hist");
     TH1F *emuDiff_allBkg_staterr = (TH1F *)emuDiff_allBkg->Clone("emuDiff_allBkg_staterr");
     emuDiff_allBkg_staterr->SetFillColor(kBlue);
     emuDiff_allBkg_staterr->SetFillStyle(3002);
     emuDiff_allBkg_staterr->SetLineColor(0);
     emuDiff_allBkg_staterr->SetLineWidth(0);
     emuDiff_allBkg_staterr->Draw("E2same");

     emuDiff_data->SetLineWidth(2);
     emuDiff_data->SetLineColor(kBlack);
     emuDiff_data->SetMarkerStyle(8);
     emuDiff_data->SetMarkerSize(0.8);
     emuDiff_data->Draw("same");

     if (plotSignal[0] > 0.) {
       emuDiff_sig1->SetLineColor(kGreen);
       emuDiff_sig1->SetLineWidth(2);
     }
     if (plotSignal[1] > 0.) {
       emuDiff_sig2->SetLineColor(kBlue);
       emuDiff_sig2->SetLineWidth(2);
     }
     if (plotSignal[2] > 0.) {
       emuDiff_sig3->SetLineColor(kCyan);
       emuDiff_sig3->SetLineWidth(2);
     }
     if (plotSignal[3] > 0.) {
       emuDiff_sig4->SetLineColor(kMagenta);
       emuDiff_sig4->SetLineWidth(2);
     }
     if (plotSignal[0] > 0.) emuDiff_sig1->Draw("histsame");
     if (plotSignal[1] > 0.) emuDiff_sig2->Draw("histsame");
     if (plotSignal[2] > 0.) emuDiff_sig3->Draw("histsame");
     if (plotSignal[3] > 0.) emuDiff_sig4->Draw("histsame");

     //bgStack->Draw("hist");
     //if (!plotPull) bgStack->GetXaxis()->SetTitle(testPlots[var].fXaxisTitle);
     //bgStack->GetXaxis()->SetTitleFont(font);
     //bgStack->GetXaxis()->SetTitleSize(0.047);
     //bgStack->GetXaxis()->SetLabelSize(0.05);
     //bgStack->GetXaxis()->SetLabelFont(font);
     //bgStack->GetXaxis()->SetMoreLogLabels();
     //bgStack->GetXaxis()->SetNoExponent();
     //bgStack->GetYaxis()->SetTitle("Events");
     //bgStack->GetYaxis()->SetTitleFont(font);
     //bgStack->GetYaxis()->SetTitleSize(0.047);
     //bgStack->GetYaxis()->SetTitleOffset(1.1);
     //bgStack->GetYaxis()->SetLabelFont(font);
     //bgStack->GetYaxis()->SetLabelSize(0.05);
     //if (!logPlot) bgStack->SetMinimum(0.);
     //else bgStack->SetMinimum(0.5);
     //bgStack->Draw("hist");
     //bgHist->SetFillColor(kBlue);
     //bgHist->SetFillStyle(3002);
     //bgHist->Draw("E2same");

     //emuTest_data.back()->SetLineWidth(2);
     //emuTest_data.back()->SetLineColor(kBlack);
     //emuTest_data.back()->SetMarkerStyle(8);
     //emuTest_data.back()->SetMarkerSize(0.8);
     //emuTest_data.back()->Draw("same");
     ////TGraphAsymmErrors* dataHist = makeDataGraph(emuTest_data.back(), 1, false);      
     ////dataHist->SetMarkerStyle(20);
     ////dataHist->SetMarkerSize(1.1);
     ////dataHist->Draw("PZ");

     //if (plotSignal[0] > 0.) emuTest_sig1.back()->Draw("histsame");
     //if (plotSignal[1] > 0.) emuTest_sig2.back()->Draw("histsame");
     //if (plotSignal[2] > 0.) emuTest_sig3.back()->Draw("histsame");
     //if (plotSignal[3] > 0.) emuTest_sig4.back()->Draw("histsame");
     
     // legent and labels
     TLegend legend(0.700, 0.566, 0.891, 0.885);
     legend.SetTextFont(font);
     legend.SetTextSize(0.042);
     legend.SetBorderSize(0);
     legend.SetLineColor(1);
     legend.SetLineStyle(1);
     legend.SetLineWidth(1);
     legend.SetFillColor(19);
     legend.SetFillStyle(0);
     legend.AddEntry(emuTest_data.back(), "DATA");
     legend.AddEntry(emuTest_ttbar.back(), "Backgrounds", "F");
     //if (dy_sample_type == 0) legend.AddEntry(emuTest_ztautau.back(), "#gamma/Z#rightarrow#tau#tau", "F");
     //legend.AddEntry(emuTest_ww.back(), "WW", "F");
     //legend.AddEntry(emuTest_wz.back(), "WZ", "F");
     //legend.AddEntry(emuTest_zz.back(), "ZZ", "F");
     //legend.AddEntry(emuTest_tw.back(), "tW", "F");
     //if (dy_sample_type == 0) {
     //  legend.AddEntry(emuTest_zmumu.back(), "#gamma/Z#rightarrow#mu#mu", "F");
     //  legend.AddEntry(emuTest_zee.back(), "#gamma/Z #rightarrow ee", "F");
     //} else {
     //  legend.AddEntry(emuTest_zjetsll.back(), "#gamma/Z + jets #rightarrow ll", "F");
     //}
     //if (qcdEst != 2) legend.AddEntry(emuTest_wjets.back(), "W+jets", "F");
     //if (plotQcd) legend.AddEntry(emuTest_qcd.back(), "jets (data)", "F");
     legend.AddEntry(emuDiff_allBkg_staterr, "stat error", "F");
     if (plotSignal[0] > 1.) legend.AddEntry(emuTest_sig1.back(), Form("Z'/#gamma'500 (x%.0f)", plotSignal[0]) ,"l");
     else if (plotSignal[0] == 1.) legend.AddEntry(emuTest_sig1.back(), "Z'/#gamma'500" ,"l");
     if (plotSignal[1] > 1.) legend.AddEntry(emuTest_sig2.back(), Form("Z'/#gamma'750 (x%.0f)", plotSignal[1]) ,"l");
     else if (plotSignal[1] == 1.) legend.AddEntry(emuTest_sig2.back(), "Z'/#gamma'750" ,"l");
     if (plotSignal[2] > 1.) legend.AddEntry(emuTest_sig3.back(), Form("Z'/#gamma'1000 (x%.0f)", plotSignal[2]) ,"l");
     else if (plotSignal[2] == 1.) legend.AddEntry(emuTest_sig3.back(), "Z'/#gamma'1000" ,"l");
     if (plotSignal[3] > 1.) legend.AddEntry(emuTest_sig4.back(), Form("Z'/#gamma'1250 (x%.0f)", plotSignal[3]) ,"l");
     else if (plotSignal[3] == 1.) legend.AddEntry(emuTest_sig4.back(), "Z'/#gamma'1250" ,"l");
     legend.SetBorderSize(0);
     legend.DrawClone("same");
     
     TLatex *tex = new TLatex();
     tex->SetNDC();
     tex->SetTextFont(font);
     tex->SetLineWidth(2);
     tex->SetTextSize(0.042);
     if (prelim) tex->DrawLatex(0.596, 0.937, "CMS Preliminary, 8 TeV, 19.7 fb^{-1}");
     else tex->DrawLatex(0.596, 0.937, "CMS, 8 TeV, 19.7 fb^{-1}");
     //if (m == 1) tex->DrawLatex(0.430, 0.849, "e in barrel");
     //if (m == 2) tex->DrawLatex(0.430, 0.849, "e in endcap");
     stringstream sStream;
     sStream << testPlots[var].fTitle << nameSign[sig].Data() << nameReg[reg].Data();
     tex->DrawLatex(0.109, 0.937, sStream.str().c_str());

     // plot a (data - bg)/bg histogram below the spectrum
     if (plotPull) {
       dataOverBgHist->Add(bgHist, -1.);
       dataOverBgHist->Divide(bgHist);

       float fontScaleBot = 1.;
       TPad *pullPad = new TPad("pullPad" + testVar + + sign[sig] + region[reg] + shapeUncName, "(data - bg) / bg" + testVar + nameSign[sig] + nameReg[reg], 0., 0., 1., 0.33);
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

       dataOverBgHist->SetLineWidth(1);
       dataOverBgHist->SetLineColor(kBlack);
       dataOverBgHist->SetMarkerStyle(20);
       dataOverBgHist->SetMarkerSize(1.1);

       dataOverBgHist->GetXaxis()->SetTitle(testPlots[var].fXaxisTitle);
       dataOverBgHist->GetXaxis()->SetTitleFont(font);
       dataOverBgHist->GetXaxis()->SetTitleSize(0.047 * fontScaleBot);
       dataOverBgHist->GetXaxis()->SetTitleOffset(1.);
       dataOverBgHist->GetXaxis()->SetLabelFont(font);
       dataOverBgHist->GetXaxis()->SetLabelSize(0.05 * fontScaleBot);
       dataOverBgHist->GetXaxis()->SetMoreLogLabels();
       dataOverBgHist->GetXaxis()->SetNoExponent();
       //dataOverBgHist->GetXaxis()->SetRangeUser(xRangeMinRatio, xRangeMaxRatio);
       dataOverBgHist->GetYaxis()->SetTitle("(data-bkg)/bkg");
       dataOverBgHist->GetYaxis()->SetTitleFont(font);
       dataOverBgHist->GetYaxis()->SetTitleSize(0.047 * fontScaleBot);
       dataOverBgHist->GetYaxis()->SetTitleOffset(1.1 / fontScaleBot);
       dataOverBgHist->GetYaxis()->SetLabelFont(font);
       dataOverBgHist->GetYaxis()->SetLabelSize(0.05 * fontScaleBot);
       dataOverBgHist->GetYaxis()->SetRangeUser(yRangeMinRatio, yRangeMaxRatio);

       dataOverBgHist->Draw();
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
     outfile->cd();
     emuPlot->Write();

     cout << testVar.Data() << shapeUncName << nameSign[sig].Data() << nameReg[reg].Data() << " plotted" << endl;
     if (shUnc > 0) input.cd("..");
  }
  outFile->Close();
  outfile->Close();
  std::cout << "Written canvases to " << outfileName << std::endl;
}

// plot a range of control variables
void PlotRange(int sign = 0, unsigned int region = 0, unsigned int from = 1, unsigned int to = 45)
{
  if (to == 0 || to > 45) to = 45;
  if (from == 0 || from > to) from = 1;
  for (unsigned int i = from; i <= to; ++i)
    macro_MakeTestPlot2(i, sign, region);
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

