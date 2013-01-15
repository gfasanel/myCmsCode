#include <string>
#include <sstream>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"

TH1F * MakeHistoFromBranch(TFile *input, const char *treeName, const char *brName, int &signs, int &region, const char *cutVariable, float &cutLow, float &cutHigh, vector<float> &binning, unsigned int &flags, bool normToBinWidth = false, float userScale = 1.);

void macro_MakeTestPlot(unsigned int k = 0, int l = 0, unsigned int m = 0)
{ 
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("testEmuSpecHEEP4_3692pb-1.root", "open");
  TFile input("./testEmuSpecHEEP41_19619pb-1.root", "open");

  const float lumi = 19619.;

  bool plotType[2];
  plotType[0] = true;  // normal plot
  plotType[1] = false;  // cumulative plot

  bool plotQcd = true;
  bool logPlot = false;
  const unsigned int numVars = 33;

  // plot style
  int ttbarColour = TColor::GetColor("#ff6666");
  int zttColour = TColor::GetColor("#ff4d4d");
  int wwColour = TColor::GetColor("#ff3333");
  int wzColour = TColor::GetColor("#ff0f0f");
  int zzColour = TColor::GetColor("#eb0000");
  int twColour = TColor::GetColor("#db0000");
  int wjetColour=  TColor::GetColor("#66b3ff");
  int zmmColour=  TColor::GetColor("#80bfff");
  int zeeColour=  TColor::GetColor("#99ccff");
  int jetBkgColour = TColor::GetColor("#ffff66");
  
  int font = 42; //62
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString testVar[33] = {"mass",
                         "pfMet",
                         "nVtx",
                         "dXYFstPVtx",
                         "dZFstPVtx",
                         "rho",
                         "numOfJets",
                         "numOfJetsPt20",
                         "numOfJetsPt30",
                         "dPhi",
                         "eleEt",
                         "eleEta",
                         "elePhi",
                         "eleDEta",
                         "eleDPhi",
                         "eleHOE",
                         "eleSigmaIEIE",
                         "eleEcalIso",
                         "eleHcalIso12",
                         "eleTrkIso",
                         "eleLostHits",
                         "muIsoCombRel",
                         "muEtEleOPtMu",
                         "muPtPlusOPtMinus",
                         "muPt",
                         "muEta",
                         "muPhi",
                         "muHitLayers",
                         "muPxlHits",
                         "muMuHits",
                         "muDZFstPVtx",
                         "muNSeg",
                         "muTrkIso03"
                        };
  float lowBins[33] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2.5, -3.2, -0.008, -0.06, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2.5, -3.2, 0., 0., 0., -0.2, 0., 0.};
  float highBins[33] = {1500., 500., 80., 0.2, 20., 50., 20., 20., 15., 3.2, 500., 2.5, 3.2, 0.008, 0.06, 0.06, 0.04, 15., 15., 5., 3., 5., 20., 20., 500., 2.5, 3.2, 20., 10., 60., 0.2, 10., 20.};
  int nBinss[33] = {75, 50, 80, 40, 40, 50, 20, 20, 15, 64, 50, 50, 64, 32, 48, 100, 50, 30, 30, 50, 3, 50, 40, 40, 50, 50, 64, 20, 10, 60, 100, 10, 40};

  TString sign[7] = {"_e-mu+", "_e+mu-", "_OS", "", "_SS", "_++", "_--"};
  TString nameSign[7] = {" e-mu+", " e+mu-", " OS", "", " SS", " ++", " --"};
  TString reg[3] = {"", "_EB", "_EE"};
  TString nameReg[3] = {"", " EB", " EE"};

  // plot a list with possible test histograms
  if (k > numVars) k = 0;
  if (k == 0) {
    cout << "Use macro_MakeTestPlot.C(x, y, z) with x {1-" << numVars << "} being the \n";
    cout << "number of the test histogram to plot, y {-3 - +3} selects the charge with \n";
    cout << "the scheme e-mu: -+, +-, OS, ALL, SS, ++, -- and z {0-2} selects the \n";
    cout << "detector EB+EE, EB or EE events for the electron." << endl;
    cout << "-----------------------------------" << endl;
    for (unsigned int i = 0; i < numVars; ++i) {
      cout << testVar[i];
      unsigned int j = 30 - testVar[i].Sizeof();
      if (i > 8) --j;
      for (; j > 0; --j)
        cout << " ";
      cout << i + 1 << endl;
    }
    cout << "-----------------------------------" << endl;
    input.Close();
    return;
  }
  --k;
  // sanity checks for sign and region input
  if (abs(l) > 3) l = 0;
  if (m > 2) m = 0;

  std::vector<TH1F *> emuTest_data;
  std::vector<TH1F *> emuTest_ttbar;
  std::vector<TH1F *> emuTest_ztautau;
  std::vector<TH1F *> emuTest_ww;
  std::vector<TH1F *> emuTest_wz;
  std::vector<TH1F *> emuTest_tw;
  std::vector<TH1F *> emuTest_wjets;
  std::vector<TH1F *> emuTest_zmumu;
  std::vector<TH1F *> emuTest_zee;
  std::vector<TH1F *> emuTest_zz;
  std::vector<TH1F *> emuTest_qcd;

  input.cd();

  // make the correct binning
  std::vector<float> binning;
  float lowBin = lowBins[k];
  float highBin = highBins[k];
  int nBins = nBinss[k];
  for (float bin = lowBin; bin <= highBin; bin += (highBin - lowBin) / nBins)
    binning.push_back(bin);

  float totMcWeight = 1.;
  THashList *mcWeights = (THashList *)input.Get("mcWeights");
  TParameter<float> *mcWeight = (TParameter<float> *)mcWeights->FindObject("ttbar");
  TParameter<float> *mcWeight700to1000 = (TParameter<float> *)mcWeights->FindObject("ttbar700to1000");
  TParameter<float> *mcWeight1000up = (TParameter<float> *)mcWeights->FindObject("ttbar1000up");

  // determine qcd contribution
  TH1F *ssData;
  TH1F *ssBg;
  TH1F *qcdContrib;
  if (plotQcd) {
    ssData = MakeHistoFromBranch(&input, "emuTree_data", testVar[k], 1, m, "", 0., 0., binning, 0x100);
    ssBg = MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], 1, m, "genMTtbar", 0., 700., binning, 0x1DF);
    totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight700to1000->GetVal());
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], 1, m, "genMTtbar", 700., 1000., binning, 0x19F), totMcWeight);
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", testVar[k], 1, m, "genMTtbar", 700., 1000., binning, 0x19F), totMcWeight);
    totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight1000up->GetVal());
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], 1, m, "genMTtbar", 1000., 1000000000., binning, 0x19F), totMcWeight);
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", testVar[k], 1, m, "genMTtbar", 1000., 1000000000., binning, 0x19F), totMcWeight);
    //TH1F *ssBg = MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], 1, m, "", 0., 0., binning, 0x1DF);
    //TH1F *ssBg = MakeHistoFromBranch(&input, "emuTree_ttbarto2l", testVar[k], 1, m, "", 0., 0., binning, 0x1DF);
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ztautau", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_ww", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wz", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zz", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_tw", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_wjets", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zmumu", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    ssBg->Add(MakeHistoFromBranch(&input, "emuTree_zee", testVar[k], 1, m, "", 0., 0., binning, 0x1DF));
    qcdContrib = (TH1F *)ssData->Clone("qcdContrib_SS");
    qcdContrib->Add(ssBg, -1);
    for (int i = 0; i < qcdContrib->GetNbinsX() + 2; ++i) {
      if (qcdContrib->GetBinContent(i) < 0) qcdContrib->SetBinContent(i, 0.);
    }
    cout << "expected SS QCD events: " << ssData->Integral() - ssBg->Integral() << endl;
    cout << "derived SS QCD events: " << qcdContrib->Integral() << endl;
    cout << "scale factor: " << (ssData->Integral() - ssBg->Integral()) / qcdContrib->Integral()<< endl;
    qcdContrib->Scale((ssData->Integral() - ssBg->Integral()) / qcdContrib->Integral());

    emuTest_qcd.push_back((TH1F *)qcdContrib->Clone(testVar[k] + sign[l+3] + nameReg[m] + "qcd"));
    if (l == 0) emuTest_qcd.back()->Scale(2.);
    else if (abs(l) > 1) emuTest_qcd.back()->Scale(0.5);
  }

  // get the histograms
  emuTest_data.push_back(MakeHistoFromBranch(&input, "emuTree_data", testVar[k], l, m, "", 0., 0., binning, 0x100));
  TH1F *ttbarComb = MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], l, m, "genMTtbar", 0., 700., binning, 0x1DF);
  totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight700to1000->GetVal());
  ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], l, m, "genMTtbar", 700., 1000., binning, 0x19F), totMcWeight);
  ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar700to1000", testVar[k], l, m, "genMTtbar", 700., 1000., binning, 0x19F), totMcWeight);
  totMcWeight = 1. / (1 / mcWeight->GetVal() + 1 / mcWeight1000up->GetVal());
  ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], l, m, "genMTtbar", 1000., 1000000000., binning, 0x19F), totMcWeight);
  ttbarComb->Add(MakeHistoFromBranch(&input, "emuTree_ttbar1000up", testVar[k], l, m, "genMTtbar", 1000., 1000000000., binning, 0x19F), totMcWeight);
  emuTest_ttbar.push_back(ttbarComb);
  //emuTest_ttbar.push_back(MakeHistoFromBranch(&input, "emuTree_ttbar", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  //emuTest_ttbar.push_back(MakeHistoFromBranch(&input, "emuTree_ttbarto2l", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_ztautau.push_back(MakeHistoFromBranch(&input, "emuTree_ztautau", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_ww.push_back(MakeHistoFromBranch(&input, "emuTree_ww", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_wz.push_back(MakeHistoFromBranch(&input, "emuTree_wz", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_tw.push_back(MakeHistoFromBranch(&input, "emuTree_tw", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_wjets.push_back(MakeHistoFromBranch(&input, "emuTree_wjets", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_zmumu.push_back(MakeHistoFromBranch(&input, "emuTree_zmumu", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_zee.push_back(MakeHistoFromBranch(&input, "emuTree_zee", testVar[k], l, m, "", 0., 0., binning, 0x1DF));
  emuTest_zz.push_back(MakeHistoFromBranch(&input, "emuTree_zz", testVar[k], l, m, "", 0., 0., binning, 0x1DF));

  // make a histogram stack with the bg 
  THStack *bgStack = new THStack("bgStack" + sign[l+3] + reg[m], testVar[k] + sign[l+3] + nameReg[m]);
  if (plotQcd) bgStack->Add(emuTest_qcd.back());
  bgStack->Add(emuTest_zee.back());
  bgStack->Add(emuTest_zmumu.back());
  bgStack->Add(emuTest_wjets.back());
  bgStack->Add(emuTest_tw.back());
  bgStack->Add(emuTest_zz.back());
  bgStack->Add(emuTest_wz.back());
  bgStack->Add(emuTest_ww.back());
  bgStack->Add(emuTest_ztautau.back());
  bgStack->Add(emuTest_ttbar.back());

  TCanvas *emuPlot = new TCanvas("emuPlot" + testVar[k] + sign[l+3] + reg[m], "emu Spectrum " + testVar[k] + nameSign[l+3] + nameReg[m], 100, 100, 800, 600);
  emuPlot->cd();
  emuPlot->SetBorderMode(0);
  emuPlot->SetBorderSize(2);
  emuPlot->SetFrameBorderMode(0);
  emuPlot->SetFillColor(0);
  emuPlot->SetFrameFillColor(0);
  if (logPlot) emuPlot->SetLogy();
  emuPlot->SetLeftMargin(0.11);
  emuPlot->SetRightMargin(0.09);
  emuPlot->SetBottomMargin(0.12);
  emuPlot->SetTopMargin(0.08);
  emuPlot->SetTickx(1);
  emuPlot->SetTicky(1);

  gStyle->SetTitleFont(font);
  gStyle->SetLabelFont(font);
  gStyle->SetLegendFont(font);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleXOffset(1.);
  gStyle->SetTitleYOffset(1.3);
  gPad->SetTicks(1, 1);

  emuTest_ttbar.back()->GetXaxis()->SetTitle("");
  emuTest_ttbar.back()->GetXaxis()->SetTitleFont(font);
  emuTest_ttbar.back()->GetXaxis()->SetTitleSize(0.047);
  emuTest_ttbar.back()->GetXaxis()->SetLabelSize(0.047);
  emuTest_ttbar.back()->GetXaxis()->SetLabelFont(font);
  emuTest_ttbar.back()->GetXaxis()->SetMoreLogLabels();
  emuTest_ttbar.back()->GetXaxis()->SetNoExponent();
  //emuTest_ttbar.back()->GetXaxis()->SetRangeUser(60., 1000.);
  emuTest_ttbar.back()->GetYaxis()->SetTitle("# of " + testVar[k] + nameSign[l+3] + " events");
  emuTest_ttbar.back()->GetYaxis()->SetTitleFont(font);
  emuTest_ttbar.back()->GetYaxis()->SetTitleSize(0.047);
  emuTest_ttbar.back()->GetYaxis()->SetTitleOffset(1.2);
  emuTest_ttbar.back()->GetYaxis()->SetLabelFont(font);
  emuTest_ttbar.back()->GetYaxis()->SetLabelSize(0.047);

  if (emuTest_data.back()->GetMaximum() > emuTest_ttbar.back()->GetMaximum()) {
    if (!logPlot)  emuTest_ttbar.back()->SetMaximum(emuTest_data.back()->GetMaximum() * 1.1);
    else emuTest_ttbar.back()->SetMaximum(emuTest_data.back()->GetMaximum() * 1.3);
  }

  //// plot spectrum
  emuTest_ttbar.back()->SetFillColor(ttbarColour);
  emuTest_ttbar.back()->SetLineColor(kBlack);
  emuTest_ttbar.back()->SetLineWidth(2);
  emuTest_ztautau.back()->SetFillColor(zttColour);
  emuTest_ztautau.back()->SetMarkerColor(zttColour);
  emuTest_ztautau.back()->SetLineColor(kBlack);
  emuTest_ztautau.back()->SetLineWidth(2);
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
  emuTest_wjets.back()->SetFillColor(wjetColour);
  emuTest_wjets.back()->SetMarkerColor(wjetColour);
  emuTest_wjets.back()->SetLineColor(kBlack);
  emuTest_wjets.back()->SetLineWidth(2);
  emuTest_zmumu.back()->SetFillColor(zmmColour);
  emuTest_zmumu.back()->SetMarkerColor(zmmColour);
  emuTest_zmumu.back()->SetLineColor(kBlack);
  emuTest_zmumu.back()->SetLineWidth(2);
  emuTest_zee.back()->SetFillColor(zeeColour);
  emuTest_zee.back()->SetMarkerColor(zeeColour);
  emuTest_zee.back()->SetLineColor(kBlack);
  emuTest_zee.back()->SetLineWidth(2);
  if (plotQcd) {
    emuTest_qcd.back()->SetFillColor(jetBkgColour);
    emuTest_qcd.back()->SetMarkerColor(jetBkgColour);
    emuTest_qcd.back()->SetLineColor(kBlack);
    emuTest_qcd.back()->SetLineWidth(2);
  }

  bgStack->Draw("hist");

  emuTest_data.back()->SetLineWidth(2);
  emuTest_data.back()->SetLineColor(kBlack);
  emuTest_data.back()->SetMarkerStyle(8);
  emuTest_data.back()->SetMarkerSize(0.8);
  emuTest_data.back()->Draw("sames");
  
  // redraw axis
  emuTest_ttbar.back()->Draw("sameaxis");

  // legent and labels
  TLegend legend(0.65, 0.8, 0.88, 0.44);
  legend.SetTextSize(0.035);
  legend.SetTextFont(font);
  legend.SetFillColor(0);
  legend.AddEntry(emuTest_data.back(), "DATA");
  legend.AddEntry(emuTest_ttbar.back(), "t#bar{t}", "F");
  legend.AddEntry(emuTest_ztautau.back(), "#gamma/Z#rightarrow#tau#tau", "F");
  legend.AddEntry(emuTest_ww.back(), "WW", "F");
  legend.AddEntry(emuTest_wz.back(), "WZ", "F");
  legend.AddEntry(emuTest_zz.back(), "ZZ", "F");
  legend.AddEntry(emuTest_tw.back(), "tW", "F");
  legend.AddEntry(emuTest_wjets.back(), "W+jets", "F");
  legend.AddEntry(emuTest_zmumu.back(), "#gamma/Z#rightarrow#mu#mu", "F");
  legend.AddEntry(emuTest_zee.back(), "#gamma/Z #rightarrow ee", "F");
  if (plotQcd) legend.AddEntry(emuTest_qcd.back(), "QCD", "F");
  legend.SetBorderSize(0);
  legend.DrawClone("sames");
  
  TPaveLabel labelPrelim(0.7, 0.91, 0.9, 0.82, "CMS Preliminary", "brNDC");
  labelPrelim.SetFillColor(0);
  labelPrelim.SetFillStyle(0);
  labelPrelim.SetBorderSize(0);
  labelPrelim.SetTextSize(0.42);
  labelPrelim.SetTextFont(font);
  labelPrelim.DrawClone("sames");
  
  stringstream sStream;
  sStream.str("");
  sStream << "#sqrt{s} = 8TeV,  #int L dt = " << lumi << "pb^{-1}";
  TPaveLabel labelLumi(0.3, 0.88, 0.65, 0.78, sStream.str().c_str(), "brNDC");
  labelLumi.SetFillColor(0);
  labelLumi.SetFillStyle(0);
  labelLumi.SetBorderSize(0);
  labelLumi.SetTextSize(0.42);
  labelLumi.SetTextFont(font);
  labelLumi.DrawClone("sames");

  sStream.str("");
  sStream << testVar[k].Data() << nameSign[l+3].Data() << nameReg[m].Data();
  TPaveLabel labelVar(0.5, 0.90, 0.91, 0.99, sStream.str().c_str(), "brNDC");
  labelVar.SetFillColor(0);
  labelVar.SetFillStyle(0);
  labelVar.SetBorderSize(0);
  labelVar.SetTextSize(0.42);
  labelVar.SetTextFont(font);
  labelVar.DrawClone("sames");

  cout << testVar[k].Data() << nameSign[l+3].Data() << nameReg[m].Data() << " plotted" << endl;
}

void PlotRange(int sign = 0, unsigned int region = 0, unsigned int from = 1, unsigned int to = 33)
{
  if (to == 0 || to > 33) to = 33;
  if (from == 0 || from > to) from = 1;
  for (unsigned int i = from; i <= to; ++i)
    macro_MakeTestPlot(i, sign, region);
}

// flags: [pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF|ele SF | mu SF | PU reweight]
TH1F *
MakeHistoFromBranch(TFile *input, const char *treeName, const char *brName, int &signs, int &region, const char *cutVariable, float &cutLow, float &cutHigh, vector<float> &binning, unsigned int &flags, bool normToBinWidth, float userScale)
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
    userScale *= ((TParameter<float> *)mcWeights->FindObject(treeName + 8))->GetVal();
  }
  if (flags & 1<<5) userScale *= ((TParameter<float> *)input->Get("trgEff"))->GetVal();
  if (flags & 1<<4) userScale *= ((TParameter<float> *)input->Get("trgDataMcScaleFactor"))->GetVal();
  if (flags & 1<<1) userScale *= ((TParameter<float> *)input->Get("muScaleFactor"))->GetVal();
  float eleScaleFactorEB = ((TParameter<float> *)input->Get("eleScaleFactorEB"))->GetVal();
  float eleScaleFactorEE = ((TParameter<float> *)input->Get("eleScaleFactorEE"))->GetVal();
  float lumiScaleFactorEB = ((TParameter<float> *)input->Get("lumiScaleFactorEB"))->GetVal();
  float lumiScaleFactorEE = ((TParameter<float> *)input->Get("lumiScaleFactorEE"))->GetVal();

  //cout << "Scalefactor = " << userScale << endl;

  // get the tree
  TTree *tree;
  tree = (TTree *)input->Get(treeName);

  // get the branch
  float var;
  TBranch *bVar;
  tree->SetBranchAddress(brName, &var, &bVar);

  // get auxillary branches
  bool passTrg;
  float puWeight;
  int eCharge;
  int muCharge;
  int evtRegion;
  float cutVar = 0.;
  TBranch *bPassTrg;
  TBranch *bPuWeight;
  TBranch *bECharge;
  TBranch *bMuCharge;
  TBranch *bEvtRegion;
  TBranch *bCutVar;
  tree->SetBranchAddress("passTrg", &passTrg, &bPassTrg);
  tree->SetBranchAddress("puWeight", &puWeight, &bPuWeight);
  tree->SetBranchAddress("eCharge", &eCharge, &bECharge);
  tree->SetBranchAddress("muCharge", &muCharge, &bMuCharge);
  tree->SetBranchAddress("evtRegion", &evtRegion, &bEvtRegion);
  if (cutVariable[0] != '\0') tree->SetBranchAddress(cutVariable, &cutVar, &bCutVar);

  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus(brName,1);
  tree->SetBranchStatus("passTrg",1);
  tree->SetBranchStatus("puWeight",1);
  tree->SetBranchStatus("eCharge",1);
  tree->SetBranchStatus("muCharge",1);
  tree->SetBranchStatus("evtRegion",1);
  if (cutVariable[0] != '\0') tree->SetBranchStatus(cutVariable,1);

  Long64_t nEntries = (*tree).GetEntries();
  for (unsigned int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // trigger fired?
    if ((flags & 1<<8) && passTrg == false) continue;

    float scaleFactor = userScale;
    // set lumi and electron scalefactor according to detector region
    if (evtRegion == 0 && flags & 1<<2) scaleFactor *= eleScaleFactorEB;
    if (evtRegion == 1 && flags & 1<<2) scaleFactor *= eleScaleFactorEE;
    if (evtRegion == 0 && flags & 1<<3) scaleFactor *= lumiScaleFactorEB;
    if (evtRegion == 1 && flags & 1<<3) scaleFactor *= lumiScaleFactorEE;

    // PU reweight
    if (flags & 1) scaleFactor *= puWeight;

    // get only the desired charge combination. Scheme e-mu -3 to +3: -+, +-, OS, ALL, SS, ++, --
    if (signs < 0 && (eCharge * muCharge) > 0) continue; // OS
    if (signs > 0 && (eCharge * muCharge) < 0) continue; // SS
    if (abs(signs) == 3 && eCharge > 0) continue; // e-mu+ or e-mu-
    if (abs(signs) == 2 && eCharge < 0) continue; // e+mu- or e+mu+

    // get only desired detector region events. Scheme: EB+EE=0, EB=1, EE=2
    if (region == 1 && evtRegion == 1) continue;
    if (region == 2 && evtRegion == 0) continue;

    // user defined cut
    if (cutVariable[0] != '\0')
      if (cutVar < cutLow || cutVar >= cutHigh) continue;

    if (normToBinWidth) scaleFactor /= histo->GetBinWidth(histo->FindBin(var));
    histo->Fill(var, scaleFactor);
  }

  //cout << "integral: " << histo->Integral() << "       overflow: " << histo->GetBinContent(histo->GetNbinsX() + 1) << endl;
  return histo;
}

