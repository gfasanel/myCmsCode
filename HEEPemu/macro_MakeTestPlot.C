#include <string>
#include <sstream>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"

void macro_MakeTestPlot(unsigned int k = 0, unsigned int l = 0)
{ 
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("testEmuSpecHEEP4_3692pb-1.root", "open");
  TFile input("./plots_25dec2012/testEmuSpecHEEP41_19619pb-1.root", "open");

  const float lumi = 19619.;

  bool plotType[2];
  plotType[0] = true;  // normal plot
  plotType[1] = false;  // cumulative plot

  bool logPlot = false;
  const unsigned int numVars = 35;

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

  TString testVar[35] = {"emuMassEB",
                         "emuMassEE",
                         "PFMET",
                         "nVtx",
                         "dDz",
                         "dDzBeamSpot",
                         "dDzFirstPVtx",
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
                         "muDxyBeamSpot",
                         "muNSeg",
                         "muTrkIso03"
                        };
  TString histoTitleTestVar[35] = {"emuMassEB",
                                   "emuMassEE",
                                   "PFMET",
                                   "nVtx",
                                   "dDz",
                                   "dDzBeamSpot",
                                   "dDzFirstPVtx",
                                   "rho",
                                   "nJets",
                                   "nJetsPt20",
                                   "nJetsPt30",
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
                                   "muDxyBeamSpot",
                                   "muNSeg",
                                   "muTrkIso03"
                                  };
  TString sign[3] = {"", "_LS", "_OS"};
  TString nameSign[3] = {"", " LS", " OS"};

  // plot a list with possible test histograms
  if (k > numVars) k = 0;
  if (k == 0) {
    cout << "Use macro_MakeTestPlot.C(x, y) with x {1-" << numVars << "} being the number of the test histogram to plot and y {0-2} selects ALL, LS or OS." << endl;
    cout << "-----------------------------------" << endl;
    for (unsigned int i = 0; i < numVars; ++i) {
      cout << histoTitleTestVar[i];
      unsigned int j = 30 - histoTitleTestVar[i].Sizeof();
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
  //std::vector<TH1F *> emuTest_qcd;

  input.cd();

  // get the histograms
  emuTest_data.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_data"));
  emuTest_ttbar.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_ttbar"));
  emuTest_ztautau.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_ztautau"));
  emuTest_ww.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_ww"));
  emuTest_wz.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_wz"));
  emuTest_tw.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_tw"));
  emuTest_wjets.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_wjets"));
  emuTest_zmumu.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_zmumu"));
  emuTest_zee.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_zee"));
  emuTest_zz.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_zz"));
  //emuTest_qcd.push_back((TH1F *)gDirectory->Get(testVar[k] + sign[l] + "_qcd"));

  TCanvas *emuPlot = new TCanvas("emuPlot" + testVar[k] + sign[l], "emu Spectrum " + histoTitleTestVar[k] + nameSign[l], 100, 100, 800, 600);
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
  emuTest_ttbar.back()->GetYaxis()->SetTitle("# of " + histoTitleTestVar[k] + nameSign[l] + " events");
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
  emuTest_ttbar.back()->DrawClone("HIST");
  emuTest_ztautau.back()->SetFillColor(zttColour);
  emuTest_ztautau.back()->SetMarkerColor(zttColour);
  emuTest_ztautau.back()->SetLineColor(kBlack);
  emuTest_ztautau.back()->SetLineWidth(2);
  emuTest_ztautau.back()->Draw("HISTsames");
  emuTest_ww.back()->SetFillColor(wwColour);
  emuTest_ww.back()->SetMarkerColor(wwColour);
  emuTest_ww.back()->SetLineColor(kBlack);
  emuTest_ww.back()->SetLineWidth(2);
  emuTest_ww.back()->Draw("HISTsames");
  emuTest_wz.back()->SetFillColor(wzColour);
  emuTest_wz.back()->SetMarkerColor(wzColour);
  emuTest_wz.back()->SetLineColor(kBlack);
  emuTest_wz.back()->SetLineWidth(2);
  emuTest_wz.back()->Draw("HISTsames");
  emuTest_zz.back()->SetFillColor(zzColour);
  emuTest_zz.back()->SetMarkerColor(zzColour);
  emuTest_zz.back()->SetLineColor(kBlack);
  emuTest_zz.back()->SetLineWidth(2);
  emuTest_zz.back()->Draw("HISTsames");
  emuTest_tw.back()->SetFillColor(twColour);
  emuTest_tw.back()->SetMarkerColor(twColour);
  emuTest_tw.back()->SetLineColor(kBlack);
  emuTest_tw.back()->SetLineWidth(2);
  emuTest_tw.back()->Draw("HISTsames");
  emuTest_wjets.back()->SetFillColor(wjetColour);
  emuTest_wjets.back()->SetMarkerColor(wjetColour);
  emuTest_wjets.back()->SetLineColor(kBlack);
  emuTest_wjets.back()->SetLineWidth(2);
  emuTest_wjets.back()->Draw("HISTsames");
  emuTest_zmumu.back()->SetFillColor(zmmColour);
  emuTest_zmumu.back()->SetMarkerColor(zmmColour);
  emuTest_zmumu.back()->SetLineColor(kBlack);
  emuTest_zmumu.back()->SetLineWidth(2);
  emuTest_zmumu.back()->Draw("HISTsames");
  emuTest_zee.back()->SetFillColor(zeeColour);
  emuTest_zee.back()->SetMarkerColor(zeeColour);
  emuTest_zee.back()->SetLineColor(kBlack);
  emuTest_zee.back()->SetLineWidth(2);
  emuTest_zee.back()->Draw("HISTsames");
  //emuTest_qcd.back()->SetFillColor(jetBkgColour);
  //emuTest_qcd.back()->SetMarkerColor(jetBkgColour);
  //emuTest_qcd.back()->SetLineColor(kBlack);
  //emuTest_qcd.back()->SetLineWidth(2);
  //emuTest_qcd.back()->Draw("HISTsames");

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
  //legend.AddEntry(emuTest_qcd.back(), "QCD", "F");
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

  cout << histoTitleTestVar[k].Data() << nameSign[l].Data() << " plotted" << endl;
}

void PlotRange(unsigned int sign = 0, unsigned int from = 1, unsigned int to = 35)
{
  if (to == 0 || to > 35) to = 35;
  if (from == 0 || from > to) from = 1;
  for (unsigned int i = from; i <= to; ++i)
    macro_MakeTestPlot(i, sign);
}

