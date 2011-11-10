  #include <string>
  #include <sstream>
  #include <vector>
  #include "TH1F.h"
  #include "TCanvas.h"

  vector<TH1F *> emuTest_data;
  vector<TH1F *> emuTest_ttbar;
  vector<TH1F *> emuTest_ztautau;
  vector<TH1F *> emuTest_ww;
  vector<TH1F *> emuTest_wz;
  vector<TH1F *> emuTest_tw;
  vector<TH1F *> emuTest_wjets;
  vector<TH1F *> emuTest_zmumu;
  vector<TH1F *> emuTest_zee;
  //vector<TH1F *> emuTest_zz;
  //vector<TH1F *> emuTest_qcd;

int macro_MakeTestPlot(unsigned int k = 0)
{ 
  // parameters //////////////////////////////////////////////////////////////
  TFile input("testEmuSpec4699pb-1.root", "open");
  //TFile input("testEmuSpec3534pb-1.root", "open");
  //TFile input("testEmuSpec_withoutTriggerBit_3534pb-1.root", "open");
  //TFile input("./PUreweightTests20111020/testEmuSpec_PUrew20_3534pb-1.root", "open");

  const float lumi = 4699;

  bool plotType[2];
  plotType[0] = true;  // normal plot
  plotType[1] = false;  // cumulative plot

  bool logPlot = false;
  ////////////////////////////////////////////////////////////////////////////

  if (k > 26) k = 0;
  if (k == 0) {
    cout << "Use .x macro_MakeTestPlot.C(x) with x being the number of the test histogram to plot." << endl;
    cout << "-----------------------------------" << endl;
    cout << "MET                 1" << endl;
    cout << "nVtx                2" << endl;
    cout << "dPhi                3" << endl;
    cout << "elePt               4" << endl;
    cout << "eleEta              5" << endl;
    cout << "elePhi              6" << endl;
    cout << "eleId1              7" << endl;
    cout << "eleId2              8" << endl;
    cout << "eleId3              9" << endl;
    cout << "eleIso1            10" << endl;
    cout << "eleIso2            11" << endl;
    cout << "eleIso3            12" << endl;
    cout << "muIsoCombRel       13" << endl;
    cout << "muPtEleOPtMu       14" << endl;
    cout << "muPtPlusOPtMinus   15" << endl;
    cout << "muPt               16" << endl;
    cout << "muEta              17" << endl;
    cout << "muPhi              18" << endl;
    cout << "muId1              19" << endl;
    cout << "muId2              20" << endl;
    cout << "muId3              21" << endl;
    cout << "muIso1             22" << endl;
    cout << "muIso2             23" << endl;
    cout << "muIso3             24" << endl;
    cout << "numOfJets          25" << endl;
    cout << "numOfJetsPt15      26" << endl;
    cout << "-----------------------------------" << endl;
    return 0;
  }
  --k;

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString testVar[26] = {"MET",
                         "nVtx",
                         "dPhi",
                         "elePt",
                         "eleEta",
                         "elePhi",
                         "eleId1",
                         "eleId2",
                         "eleId3",
                         "eleIso1",
                         "eleIso2",
                         "eleIso3",
                         "muIsoCombRel",
                         "muPtEleOPtMu",
                         "muPtPlusOPtMinus",
                         "muPt",
                         "muEta",
                         "muPhi",
                         "muId1",
                         "muId2",
                         "muId3",
                         "muIso1",
                         "muIso2",
                         "muIso3",
                         "numOfJets",
                         "numOfJetsPt15"
                        };
  TString histoTitleTestVar[26] = {"MET",
                                   "nVtx",
                                   "dPhi",
                                   "elePt",
                                   "eleEta",
                                   "elePhi",
                                   "eleId1",
                                   "eleId2",
                                   "eleId3",
                                   "eleIso1",
                                   "eleIso2",
                                   "eleIso3",
                                   "muIsoCombRel",
                                   "muPtEleOPtMu",
                                   "muPtPlusOPtMinus",
                                   "muPt",
                                   "muEta",
                                   "muPhi",
                                   "muId1",
                                   "muId2",
                                   "muId3",
                                   "muIso1",
                                   "muIso2",
                                   "muIso3",
                                   "numOfJets",
                                   "numOfJetsPt15"
                                  };
  TString nameSuffix[2] = {"", "Cumul"};
  TString titleSuffix[2] = {"", " - Cumulative"};

  // loop to get normal and cumulated spectrum
  for (unsigned int j = 0; j < 2; ++j) {
    if (!plotType[j]) continue;
    input.cd();

    // get the histograms
    emuTest_data.push_back((TH1F *)gDirectory->Get(testVar[k] + "_data"));
    emuTest_ttbar.push_back((TH1F *)gDirectory->Get(testVar[k] + "_ttbar"));
    emuTest_ztautau.push_back((TH1F *)gDirectory->Get(testVar[k] + "_ztautau"));
    emuTest_ww.push_back((TH1F *)gDirectory->Get(testVar[k] + "_ww"));
    emuTest_wz.push_back((TH1F *)gDirectory->Get(testVar[k] + "_wz"));
    emuTest_tw.push_back((TH1F *)gDirectory->Get(testVar[k] + "_tw"));
    emuTest_wjets.push_back((TH1F *)gDirectory->Get(testVar[k] + "_wjets"));
    emuTest_zmumu.push_back((TH1F *)gDirectory->Get(testVar[k] + "_zmumu"));
    emuTest_zee.push_back((TH1F *)gDirectory->Get(testVar[k] + "_zee"));
    //emuTest_zz.push_back((TH1F *)gDirectory->Get(testVar[k] + "_zz"));
    //emuTest_qcd.push_back((TH1F *)gDirectory->Get(testVar[k] + "_qcd"));

    // set unique name
    emuTest_data.back()->SetName(emuTest_data.back()->GetName() + nameSuffix[j]);
    emuTest_ttbar.back()->SetName(emuTest_ttbar.back()->GetName() + nameSuffix[j]);
    emuTest_ztautau.back()->SetName(emuTest_ztautau.back()->GetName() + nameSuffix[j]);
    emuTest_ww.back()->SetName(emuTest_ww.back()->GetName() + nameSuffix[j]);
    emuTest_wz.back()->SetName(emuTest_wz.back()->GetName() + nameSuffix[j]);
    emuTest_tw.back()->SetName(emuTest_tw.back()->GetName() + nameSuffix[j]);
    emuTest_wjets.back()->SetName(emuTest_wjets.back()->GetName() + nameSuffix[j]);
    emuTest_zmumu.back()->SetName(emuTest_zmumu.back()->GetName() + nameSuffix[j]);
    emuTest_zee.back()->SetName(emuTest_zee.back()->GetName() + nameSuffix[j]);
    //emuTest_zz.back()->SetName(emuTest_zz.back()->GetName() + nameSuffix[j]);
    //emuTest_qcd.back()->SetName(emuTest_qcd.back()->GetName() + nameSuffix[j]);

    // integrate from the right side
    if (j == 1) { 
      // loop over bins
      double error;
      for (int i = 1; i < 101; ++i) {
        emuTest_data.back()->SetBinContent(i, emuTest_data.back()->IntegralAndError(i, 100, error));
        emuTest_data.back()->SetBinError(i, error);
        emuTest_ttbar.back()->SetBinContent(i, emuTest_ttbar.back()->IntegralAndError(i, 100, error));
        emuTest_ttbar.back()->SetBinError(i, error);
        emuTest_ztautau.back()->SetBinContent(i, emuTest_ztautau.back()->IntegralAndError(i, 100, error));
        emuTest_ztautau.back()->SetBinError(i, error);
        emuTest_ww.back()->SetBinContent(i, emuTest_ww.back()->IntegralAndError(i, 100, error));
        emuTest_ww.back()->SetBinError(i, error);
        emuTest_wz.back()->SetBinContent(i, emuTest_wz.back()->IntegralAndError(i, 100, error));
        emuTest_wz.back()->SetBinError(i, error);
        emuTest_tw.back()->SetBinContent(i, emuTest_tw.back()->IntegralAndError(i, 100, error));
        emuTest_tw.back()->SetBinError(i, error);
        emuTest_wjets.back()->SetBinContent(i, emuTest_wjets.back()->IntegralAndError(i, 100, error));
        emuTest_wjets.back()->SetBinError(i, error);
        emuTest_zmumu.back()->SetBinContent(i, emuTest_zmumu.back()->IntegralAndError(i, 100, error));
        emuTest_zmumu.back()->SetBinError(i, error);
        emuTest_zee.back()->SetBinContent(i, emuTest_zee.back()->IntegralAndError(i, 100, error));
        emuTest_zee.back()->SetBinError(i, error);
        //emuTest_zz.back()->SetBinContent(i, emuTest_zz.back()->IntegralAndError(i, 100, error));
        //emuTest_zz.back()->SetBinError(i, error);
        //emuTest_qcd.back()->SetBinContent(i, emuTest_qcd.back()->IntegralAndError(i, 100, error));
        //emuTest_qcd.back()->SetBinError(i, error);
      }     
    }

    TCanvas *emuPlot = new TCanvas("emuPlot" + testVar[k] + nameSuffix[j], "emu Spectrum " + histoTitleTestVar[k] + titleSuffix[j], 100, 100, 800, 600);
    emuPlot->SetBorderMode(0);
    emuPlot->SetFrameBorderMode(0);
    emuPlot->SetFillColor(0);
    emuPlot->SetFrameFillColor(0);
    if (logPlot) emuPlot->SetLogy();
    emuPlot->cd();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTitleXOffset(1.);
    gStyle->SetTitleYOffset(1.3);

    emuTest_ttbar.back()->GetXaxis()->SetTitle("");
    emuTest_ttbar.back()->GetXaxis()->SetTitleSize(0.04);
    emuTest_ttbar.back()->GetXaxis()->SetLabelSize(0.04);
    //emuTest_ttbar.back()->GetXaxis()->SetRangeUser(60., 1000.);

    if (j == 1) emuTest_ttbar.back()->GetYaxis()->SetTitle("cumulated # of " + histoTitleTestVar[k] + " events");
    else emuTest_ttbar.back()->GetYaxis()->SetTitle("# of " + histoTitleTestVar[k] + " events");
    emuTest_ttbar.back()->GetYaxis()->SetTitleSize(0.04);
    emuTest_ttbar.back()->GetYaxis()->SetTitleOffset(1.2);

    if (emuTest_data.back()->GetMaximum() > emuTest_ttbar.back()->GetMaximum()) {
      if (!logPlot)  emuTest_ttbar.back()->SetMaximum(emuTest_data.back()->GetMaximum() * 1.1);
      else emuTest_ttbar.back()->SetMaximum(emuTest_data.back()->GetMaximum() * 1.3);
    }

    // plot spectrum
    emuTest_ttbar.back()->SetFillColor(2);
    emuTest_ttbar.back()->Draw("HIST");
    emuTest_ztautau.back()->SetFillColor(6);
    emuTest_ztautau.back()->Draw("HISTsames");
    emuTest_ww.back()->SetFillColor(5);
    emuTest_ww.back()->Draw("HISTsames");
    emuTest_wz.back()->SetFillColor(8);
    emuTest_wz.back()->Draw("HISTsames");
    emuTest_tw.back()->SetFillColor(46);
    emuTest_tw.back()->Draw("HISTsames");
    emuTest_wjets.back()->SetFillColor(4);
    emuTest_wjets.back()->Draw("HISTsames");
    emuTest_zmumu.back()->SetFillColor(7);
    emuTest_zmumu.back()->Draw("HISTsames");
    emuTest_zee.back()->SetFillColor(3);
    emuTest_zee.back()->Draw("HISTsames");
    //emuTest_zz.back()->SetFillColor(9);
    //emuTest_zz.back()->Draw("HISTsames");
    //emuTest_qcd.back()->SetFillColor(9);
    //emuTest_qcd.back()->Draw("HISTsames");
    //
    emuTest_data.back()->SetLineWidth(2);
    emuTest_data.back()->SetMarkerStyle(8);
    emuTest_data.back()->SetMarkerSize(0.8);
    emuTest_data.back()->Draw("sames");
  
    // redraw axis
    emuTest_ttbar.back()->Draw("sameaxis");

    // legent and labels
    TLegend legend(0.7, 0.8, 0.88, 0.54);
    legend.SetTextSize(0.03);
    legend.SetFillColor(0);
  
    legend.AddEntry(emuTest_data.back(), "Data");
    legend.AddEntry(emuTest_ttbar.back(), "t #bar{t} (MC)");
    legend.AddEntry(emuTest_ztautau.back(), "Z #rightarrow #tau #tau (MC)");
    legend.AddEntry(emuTest_ww.back(), "WW, (MC)");
    legend.AddEntry(emuTest_wz.back(), "WZ, (MC)");
    legend.AddEntry(emuTest_tw.back(), "tW, (MC)");
    legend.AddEntry(emuTest_wjets.back(), "W+jets (MC)");
    legend.AddEntry(emuTest_zmumu.back(), "Z #rightarrow #mu #mu (MC)");
    legend.AddEntry(emuTest_zee.back(), "Z #rightarrow ee (MC)");
    //legend.AddEntry(emuTest_zz.back(), "ZZ, (MC)");
    //legend.AddEntry(emuTest_qcd.back(), "QCD");
  
    legend.SetBorderSize(0);
    legend.DrawClone("sames");
    
    TPaveLabel labelPrelim(0.7, 0.91, 0.9, 0.82, "CMS Preliminary", "brNDC");
    labelPrelim.SetFillColor(0);
    labelPrelim.SetFillStyle(0);
    labelPrelim.SetBorderSize(0);
    labelPrelim.SetTextSize(0.40);
    labelPrelim.DrawClone("sames");
  
    stringstream sStream;
    sStream.str("");
    sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
    TPaveLabel labelLumi(0.3, 0.88, 0.65, 0.78, sStream.str().c_str(), "brNDC");
    labelLumi.SetFillColor(0);
    labelLumi.SetFillStyle(0);
    labelLumi.SetBorderSize(0);
    labelLumi.SetTextSize(0.40);
    labelLumi.DrawClone("sames");
  } // end loop over normal or cumulated
  return 0;
}
