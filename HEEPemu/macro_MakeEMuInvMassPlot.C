#include <stdio>
#include <string>
#include <sstream>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"

#define DATA 0
#define TTBAR 1
#define ZTT 2
#define WW 3
#define WZ 4
#define TW 5
#define WJET 6
#define ZMM 7
#define ZEE 8
#define ZZ 9
//#define QCD 10
#define QCD 9

#define ALL 0
#define LS 1
#define OS 2

float CalcSystErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked = false);
float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked = false);

void macro_MakeEMuInvMassPlot()
{
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("test4684pb-1.root", "open");
  //TFile input("testEmuSpec4684pb-1.root", "open");
  //TFile input("testEmuSpecSummerFallMix4684pb-1.root", "open");
  TFile input("testEmuSpecSummerFallMix4980pb-1.root", "open");
  //TFile input("testEmuSpecPureSummer114684pb-1.root", "open");
  //TFile input("testEmuSpecNewDYNoPURew4684pb-1.root", "open");
  //TFile input("testEmuSpec4699pb-1.root", "open");
  //TFile input("testEmuSpecNewDY4699pb-1.root", "open");
  //TFile input("testEmuSpec3534pb-1.root", "open");
  //TFile input("testEmuSpec_withoutTriggerBit_3534pb-1.root", "open");
  //TFile input("./PUreweightTests20111020/testEmuSpec_PUrew20_3534pb-1.root", "open");

  const float lumi = 4980.;
  const float minInvMass = 60.;

  bool plotSign[3];
  plotSign[0] = true;  // all
  plotSign[1] = true;  // LS like sign
  plotSign[2] = true;  // OS opposite sign

  bool plotType[2];
  plotType[0] = true;  // emu spectrum
  plotType[1] = true;  // cumulative emu spectrum

  bool logPlot = true;
  bool prelim = true;

  // systematical errors
  vector<float> systErrMC;
  systErrMC.push_back(0.15);  //ttbar
  systErrMC.push_back(0.);    //z->tt  //FIXME
  systErrMC.push_back(0.);    //WW     //FIXME
  systErrMC.push_back(0.);    //WZ     //FIXME
  systErrMC.push_back(0.15);  //tW     //FIXME
  systErrMC.push_back(0.);    //WJets  //FIXME
  systErrMC.push_back(0.);    //Z->mm  //FIXME
  systErrMC.push_back(0.);    //Z->ee  //FIXME
  systErrMC.push_back(0.);    //ZZ     //FIXME
  float systErrLumi = 0.022;
  float systErrEff = 0.02;
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString histoSign[3] = {"", "LS_", "OS_"};
  TString histoTitleSign[3] = {"E", "e^{#pm}#mu^{#pm} e", "e^{+}#mu^{-}/e^{-}#mu^{+} e"};
  TString xAxisTitle[3] = {"m(e#mu)", "m(e^{#pm}#mu^{#pm})", "m(e^{+}#mu^{-})/m(e^{-}#mu^{+})"};
  TString nameSuffix[2] = {"", "Cumul"};
  TString titleSuffix[2] = {"", " - Cumulative"};

  vector<TH1F *> emuMass_data;
  vector<TH1F *> emuMass_ttbar;
  vector<TH1F *> emuMass_ztautau;
  vector<TH1F *> emuMass_ww;
  vector<TH1F *> emuMass_wz;
  vector<TH1F *> emuMass_tw;
  vector<TH1F *> emuMass_wjets;
  vector<TH1F *> emuMass_zmumu;
  vector<TH1F *> emuMass_zee;
  //vector<TH1F *> emuMass_zz;
  vector<TH1F *> emuMass_qcd;

  // loop over full spectrum, LS and OS
  for (unsigned int k = 0; k < 3; ++k) {
    // loop to get nowmal and cumulated spectrum
    for (unsigned int j = 0; j < 2; ++j) {
      if (!plotType[j]) continue;
      input.cd();

      // get the histograms
      emuMass_data.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "data"));
      emuMass_ttbar.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "ttbar"));
      emuMass_ztautau.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "ztautau"));
      emuMass_ww.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "ww"));
      emuMass_wz.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "wz"));
      emuMass_tw.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "tw"));
      emuMass_wjets.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "wjets"));
      emuMass_zmumu.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "zmumu"));
      emuMass_zee.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "zee"));
      //emuMass_zz.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "zz"));
      emuMass_qcd.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "qcd"));

      // set unique name
      emuMass_data.back()->SetName(emuMass_data.back()->GetName() + nameSuffix[j]);
      emuMass_ttbar.back()->SetName(emuMass_ttbar.back()->GetName() + nameSuffix[j]);
      emuMass_ztautau.back()->SetName(emuMass_ztautau.back()->GetName() + nameSuffix[j]);
      emuMass_ww.back()->SetName(emuMass_ww.back()->GetName() + nameSuffix[j]);
      emuMass_wz.back()->SetName(emuMass_wz.back()->GetName() + nameSuffix[j]);
      emuMass_tw.back()->SetName(emuMass_tw.back()->GetName() + nameSuffix[j]);
      emuMass_wjets.back()->SetName(emuMass_wjets.back()->GetName() + nameSuffix[j]);
      emuMass_zmumu.back()->SetName(emuMass_zmumu.back()->GetName() + nameSuffix[j]);
      emuMass_zee.back()->SetName(emuMass_zee.back()->GetName() + nameSuffix[j]);
      //emuMass_zz.back()->SetName(emuMass_zz.back()->GetName() + nameSuffix[j]);
      emuMass_qcd.back()->SetName(emuMass_qcd.back()->GetName() + nameSuffix[j]);

      if (!plotSign[k]) continue;

      // integrate from the right side
      if (j == 1) { 
        // loop over bins
        double error;
	for (int i = 1; i < 151; ++i) {
	  emuMass_data.back()->SetBinContent(i, emuMass_data.back()->IntegralAndError(i, 150, error));
	  emuMass_data.back()->SetBinError(i, error);
	  emuMass_ttbar.back()->SetBinContent(i, emuMass_ttbar.back()->IntegralAndError(i, 150, error));
	  emuMass_ttbar.back()->SetBinError(i, error);
	  emuMass_ztautau.back()->SetBinContent(i, emuMass_ztautau.back()->IntegralAndError(i, 150, error));
	  emuMass_ztautau.back()->SetBinError(i, error);
	  emuMass_ww.back()->SetBinContent(i, emuMass_ww.back()->IntegralAndError(i, 150, error));
	  emuMass_ww.back()->SetBinError(i, error);
	  emuMass_wz.back()->SetBinContent(i, emuMass_wz.back()->IntegralAndError(i, 150, error));
	  emuMass_wz.back()->SetBinError(i, error);
	  emuMass_tw.back()->SetBinContent(i, emuMass_tw.back()->IntegralAndError(i, 150, error));
	  emuMass_tw.back()->SetBinError(i, error);
	  emuMass_wjets.back()->SetBinContent(i, emuMass_wjets.back()->IntegralAndError(i, 150, error));
	  emuMass_wjets.back()->SetBinError(i, error);
	  emuMass_zmumu.back()->SetBinContent(i, emuMass_zmumu.back()->IntegralAndError(i, 150, error));
	  emuMass_zmumu.back()->SetBinError(i, error);
	  emuMass_zee.back()->SetBinContent(i, emuMass_zee.back()->IntegralAndError(i, 150, error));
	  emuMass_zee.back()->SetBinError(i, error);
	  //emuMass_zz.back()->SetBinContent(i, emuMass_zz.back()->IntegralAndError(i, 150, error));
	  //emuMass_zz.back()->SetBinError(i, error);
	  emuMass_qcd.back()->SetBinContent(i, emuMass_qcd.back()->IntegralAndError(i, 150, error));
	  emuMass_qcd.back()->SetBinError(i, error);
	}     
      }

      TCanvas *emuPlot = new TCanvas("emuPlot" + histoSign[k] + nameSuffix[j], "emu Spectrum" + histoTitleSign[k] + titleSuffix[j], 100, 100, 800, 600);
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
      gPad->SetTicks(1, 1);

      emuMass_ttbar.back()->GetXaxis()->SetTitle(xAxisTitle[k] + " [GeV]");
      emuMass_ttbar.back()->GetXaxis()->SetTitleSize(0.04);
      emuMass_ttbar.back()->GetXaxis()->SetLabelSize(0.04);
      //emuMass_ttbar.back()->GetXaxis()->SetRangeUser(60., 1110.);
    
      if (j == 1) emuMass_ttbar.back()->GetYaxis()->SetTitle(histoTitleSign[k] + "vents #geq " + xAxisTitle[k]);
      else emuMass_ttbar.back()->GetYaxis()->SetTitle(histoTitleSign[k] + "vents / 10 GeV");
      emuMass_ttbar.back()->GetYaxis()->SetTitleSize(0.04);
      emuMass_ttbar.back()->GetYaxis()->SetTitleOffset(1.2);
    
      if (emuMass_data.back()->GetMaximum() > emuMass_ttbar.back()->GetMaximum()) {
        if (!logPlot)  emuMass_ttbar.back()->SetMaximum(emuMass_data.back()->GetMaximum() * 1.1);
        else emuMass_ttbar.back()->SetMaximum(emuMass_data.back()->GetMaximum() * 1.3);
      }
    
      // plot spectrum
      emuMass_ttbar.back()->SetFillColor(2);
      emuMass_ttbar.back()->Draw("HIST");
      emuMass_ztautau.back()->SetFillColor(6);
      emuMass_ztautau.back()->Draw("HISTsames");
      emuMass_ww.back()->SetFillColor(5);
      emuMass_ww.back()->Draw("HISTsames");
      emuMass_wz.back()->SetFillColor(8);
      emuMass_wz.back()->Draw("HISTsames");
      emuMass_tw.back()->SetFillColor(46);
      emuMass_tw.back()->Draw("HISTsames");
      emuMass_wjets.back()->SetFillColor(4);
      emuMass_wjets.back()->Draw("HISTsames");
      emuMass_zmumu.back()->SetFillColor(7);
      emuMass_zmumu.back()->Draw("HISTsames");
      emuMass_zee.back()->SetFillColor(3);
      emuMass_zee.back()->Draw("HISTsames");
      //emuMass_zz.back()->SetFillColor(9);
      //emuMass_zz.back()->Draw("HISTsames");
      emuMass_qcd.back()->SetFillColor(9);
      emuMass_qcd.back()->Draw("HISTsames");
      //
      emuMass_data.back()->SetLineWidth(2);
      emuMass_data.back()->SetMarkerStyle(8);
      emuMass_data.back()->SetMarkerSize(0.8);
      emuMass_data.back()->Draw("sames");
    
      // redraw axis
      emuMass_ttbar.back()->Draw("sameaxis");

      // legent and labels
      TLegend legend(0.69, 0.37, 0.85, 0.76);
      legend.SetTextSize(0.03);
      legend.SetFillColor(0);
    
      legend.AddEntry(emuMass_data.back(), "DATA");
      legend.AddEntry(emuMass_ttbar.back(), "t #bar{t} (MC)" ,"F");
      legend.AddEntry(emuMass_ztautau.back(), "Z/#gamma* #rightarrow #tau#tau (MC)" ,"F");
      legend.AddEntry(emuMass_ww.back(), "WW (MC)" ,"F");
      legend.AddEntry(emuMass_wz.back(), "WZ (MC)" ,"F");
      legend.AddEntry(emuMass_tw.back(), "tW (MC)" ,"F");
      legend.AddEntry(emuMass_wjets.back(), "W+jets (MC)" ,"F");
      legend.AddEntry(emuMass_zmumu.back(), "Z/#gamma* #rightarrow #mu#mu (MC)" ,"F");
      legend.AddEntry(emuMass_zee.back(), "Z/#gamma* #rightarrow ee (MC)" ,"F");
      //legend.AddEntry(emuMass_zz.back(), "ZZ, (MC)" ,"F");
      legend.AddEntry(emuMass_qcd.back(), "QCD" ,"F");
    
      legend.SetBorderSize(0);
      legend.DrawClone("sames");
      
      stringstream sStream;
      sStream.str("");
      if (prelim) sStream << "CMS preliminary    #sqrt{s} = 7 TeV    #int L dt = 5.0 fb^{-1}";
      else sStream << "CMS     #sqrt{s} = 7 TeV    #int L dt = 5.0 fb^{-1}";
      //sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
      TPaveLabel labelLumi(0.32, 0.78, 0.88, 0.87, sStream.str().c_str(), "brNDC");
      labelLumi.SetFillColor(0);
      labelLumi.SetFillStyle(0);
      labelLumi.SetBorderSize(0);
      labelLumi.SetTextSize(0.40);
      labelLumi.DrawClone("sames");
    } // end loop over normal or cumulated
  } // end loop over full, LS and OS

  vector<vector<TH1F *> > emuMasses;

  emuMasses.push_back(emuMass_data);
  emuMasses.push_back(emuMass_ttbar);
  emuMasses.push_back(emuMass_ztautau);
  emuMasses.push_back(emuMass_ww);
  emuMasses.push_back(emuMass_wz);
  emuMasses.push_back(emuMass_tw);
  emuMasses.push_back(emuMass_wjets);
  emuMasses.push_back(emuMass_zmumu);
  emuMasses.push_back(emuMass_zee);
  //emuMasses.push_back(emuMass_zz);
  emuMasses.push_back(emuMass_qcd);
  if (plotType[1]) {
    for (vector<vector<TH1F *> >::iterator iter = emuMasses.begin(); iter < emuMasses.end(); ++iter) {
      iter->erase(iter->begin() + 5);
      iter->erase(iter->begin() + 3);
      iter->erase(iter->begin() + 1);
    }
  }

  float systErrLuEff = sqrt(systErrLumi*systErrLumi + systErrEff*systErrEff);

  // calculate rate of syst errors
  float systErrLuEff = sqrt(systErrLumi*systErrLumi + systErrEff*systErrEff);
  vector<float> systErrMCLuEff;
  for (unsigned int it = TTBAR - 1; it < ZEE; ++it) 
     systErrMCLuEff.push_back(sqrt(systErrMC[it]*systErrMC[it] + systErrLuEff*systErrLuEff));

  // define groups of MC samples
  vector<bool> ttLikeSamples(5, true);
  vector<bool> contamSamples(5, false);
  contamSamples.push_back(true); // WJets
  contamSamples.push_back(true); // Zmm
  contamSamples.push_back(true); // Zee
  vector<bool> allSamples(8, true);

  // write numbers
  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "HEEP - TIGHT MU        Lumi        = " << lumi << "pb-1" << endl;
  //cout << "                       e pT EB     > " << bar_et << "GeV/c" << endl;
  //cout << "                       e pT EE     > " << end_et << "GeV/c" << endl;
  //cout << "                       mu pT       > " << muon_et << "GeV/c" << endl;
  //cout << "                       mu |eta|    < " << muon_etaMax << endl;
  cout << endl;
  cout << "Systematic errors" << endl;
  cout << " Luminosity:  " << systErrLumi * 100 << "%" << endl;
  cout << " Efficiency:  " << systErrEff * 100 << "%" << endl;
  cout << " ttbar:      " << systErrMC[TTBAR-1] * 100 << "%" << endl;
  cout << " Z->tautau:   " << systErrMC[ZTT-1] * 100 << "%" << endl;
  cout << " WW:          " << systErrMC[WW-1] * 100 << "%" << endl;
  cout << " WZ:          " << systErrMC[WZ-1] * 100 << "%" << endl;
  cout << " tW, tbarW:  " << systErrMC[TW-1] * 100 << "%" << endl;
  cout << " Z->mumu:     " << systErrMC[ZMM-1] * 100 << "%" << endl;
  cout << " Z->ee:       " << systErrMC[ZEE-1] * 100 << "%" << endl;
  cout << " ZZ:          " << systErrMC[ZZ-1] * 100 << "%" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |         >  60GeV/c^2          |        > 120GeV/c^2          |        > 200GeV/c^2          |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
         emuMass_data.at(ALL)->Integral(7,151), sqrt(emuMass_data.at(ALL)->Integral(7,151)),         
         emuMass_data.at(ALL)->Integral(13,151), sqrt(emuMass_data.at(ALL)->Integral(13,151)),
         emuMass_data.at(ALL)->Integral(21,151), sqrt(emuMass_data.at(ALL)->Integral(21,151)));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb ttbar      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_ttbar.at(ALL)->Integral(7,151) - emuMass_ztautau.at(ALL)->Integral(7,151), (emuMass_ttbar.at(ALL)->Integral(7,151) - emuMass_ztautau.at(ALL)->Integral(7,151)) * systErrMCLuEff[TTBAR-1], 
         emuMass_ttbar.at(ALL)->Integral(13,151) - emuMass_ztautau.at(ALL)->Integral(13,151), (emuMass_ttbar.at(ALL)->Integral(13,151) - emuMass_ztautau.at(ALL)->Integral(13,151)) * systErrMCLuEff[TTBAR-1], 
         emuMass_ttbar.at(ALL)->Integral(21,151) - emuMass_ztautau.at(ALL)->Integral(21,151), (emuMass_ttbar.at(ALL)->Integral(21,151) - emuMass_ztautau.at(ALL)->Integral(21,151)) * systErrMCLuEff[TTBAR-1]);
  printf("nb Ztautau    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_ztautau.at(ALL)->Integral(7,151) - emuMass_ww.at(ALL)->Integral(7,151), (emuMass_ztautau.at(ALL)->Integral(7,151) - emuMass_ww.at(ALL)->Integral(7,151)) * systErrMCLuEff[ZTT-1], 
         emuMass_ztautau.at(ALL)->Integral(13,151) - emuMass_ww.at(ALL)->Integral(13,151), (emuMass_ztautau.at(ALL)->Integral(13,151) - emuMass_ww.at(ALL)->Integral(13,151)) * systErrMCLuEff[ZTT-1], 
         emuMass_ztautau.at(ALL)->Integral(21,151) - emuMass_ww.at(ALL)->Integral(21,151), (emuMass_ztautau.at(ALL)->Integral(21,151) - emuMass_ww.at(ALL)->Integral(21,151)) * systErrMCLuEff[ZTT-1]); 
  printf("nb WW         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_ww.at(ALL)->Integral(7,151) - emuMass_wz.at(ALL)->Integral(7,151), (emuMass_ww.at(ALL)->Integral(7,151) - emuMass_wz.at(ALL)->Integral(7,151)) * systErrMCLuEff[WW-1], 
         emuMass_ww.at(ALL)->Integral(13,151) - emuMass_wz.at(ALL)->Integral(13,151), (emuMass_ww.at(ALL)->Integral(13,151) - emuMass_wz.at(ALL)->Integral(13,151)) * systErrMCLuEff[WW-1], 
         emuMass_ww.at(ALL)->Integral(21,151) - emuMass_wz.at(ALL)->Integral(21,151), (emuMass_ww.at(ALL)->Integral(21,151) - emuMass_wz.at(ALL)->Integral(21,151)) * systErrMCLuEff[WW-1]);
  printf("nb WZ         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_wz.at(ALL)->Integral(7,151) - emuMass_tw.at(ALL)->Integral(7,151), (emuMass_wz.at(ALL)->Integral(7,151) - emuMass_tw.at(ALL)->Integral(7,151)) * systErrMCLuEff[WZ-1], 
         emuMass_wz.at(ALL)->Integral(13,151) - emuMass_tw.at(ALL)->Integral(13,151), (emuMass_wz.at(ALL)->Integral(13,151) - emuMass_tw.at(ALL)->Integral(13,151)) * systErrMCLuEff[WZ-1], 
         emuMass_wz.at(ALL)->Integral(21,151) - emuMass_tw.at(ALL)->Integral(21,151), (emuMass_wz.at(ALL)->Integral(21,151) - emuMass_tw.at(ALL)->Integral(21,151)) * systErrMCLuEff[WZ-1]);
  printf("nb tW         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_tw.at(ALL)->Integral(7,151) - emuMass_wjets.at(ALL)->Integral(7,151), (emuMass_tw.at(ALL)->Integral(7,151) - emuMass_wjets.at(ALL)->Integral(7,151)) * systErrMCLuEff[TW-1], 
         emuMass_tw.at(ALL)->Integral(13,151) - emuMass_wjets.at(ALL)->Integral(13,151), (emuMass_tw.at(ALL)->Integral(13,151) - emuMass_wjets.at(ALL)->Integral(13,151)) * systErrMCLuEff[TW-1], 
         emuMass_tw.at(ALL)->Integral(21,151) - emuMass_wjets.at(ALL)->Integral(21,151), (emuMass_tw.at(ALL)->Integral(21,151) - emuMass_wjets.at(ALL)->Integral(21,151)) * systErrMCLuEff[TW-1]);
  cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
  printf("nb WJets      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_wjets.at(ALL)->Integral(7,151) - emuMass_zmumu.at(ALL)->Integral(7,151), (emuMass_wjets.at(ALL)->Integral(7,151) - emuMass_zmumu.at(ALL)->Integral(7,151)) * systErrMCLuEff[WJET-1], 
         emuMass_wjets.at(ALL)->Integral(13,151) - emuMass_zmumu.at(ALL)->Integral(13,151), (emuMass_wjets.at(ALL)->Integral(13,151) - emuMass_zmumu.at(ALL)->Integral(13,151)) * systErrMCLuEff[WJET-1], 
         emuMass_wjets.at(ALL)->Integral(21,151) - emuMass_zmumu.at(ALL)->Integral(21,151), (emuMass_wjets.at(ALL)->Integral(21,151) - emuMass_zmumu.at(ALL)->Integral(21,151)) * systErrMCLuEff[WJET-1]);
  printf("nb Zmumu      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_zmumu.at(ALL)->Integral(7,151) - emuMass_zee.at(ALL)->Integral(7,151), (emuMass_zmumu.at(ALL)->Integral(7,151) - emuMass_zee.at(ALL)->Integral(7,151)) * systErrMCLuEff[ZMM-1], 
         emuMass_zmumu.at(ALL)->Integral(13,151) - emuMass_zee.at(ALL)->Integral(13,151), (emuMass_zmumu.at(ALL)->Integral(13,151) - emuMass_zee.at(ALL)->Integral(13,151)) * systErrMCLuEff[ZMM-1], 
         emuMass_zmumu.at(ALL)->Integral(21,151) - emuMass_zee.at(ALL)->Integral(21,151), (emuMass_zmumu.at(ALL)->Integral(21,151) - emuMass_zee.at(ALL)->Integral(21,151)) * systErrMCLuEff[ZMM-1]);
  printf("nb Zee        | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_zee.at(ALL)->Integral(7,151) - emuMass_qcd.at(ALL)->Integral(7,151), (emuMass_zee.at(ALL)->Integral(7,151) - emuMass_qcd.at(ALL)->Integral(7,151)) * systErrMCLuEff[ZEE-1], 
         emuMass_zee.at(ALL)->Integral(13,151) - emuMass_qcd.at(ALL)->Integral(13,151), (emuMass_zee.at(ALL)->Integral(13,151) - emuMass_qcd.at(ALL)->Integral(13,151)) * systErrMCLuEff[ZEE-1], 
         emuMass_zee.at(ALL)->Integral(21,151) - emuMass_qcd.at(ALL)->Integral(21,151), (emuMass_zee.at(ALL)->Integral(21,151) - emuMass_qcd.at(ALL)->Integral(21,151)) * systErrMCLuEff[ZEE-1]); 
  //printf("nb ZZ         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
  //       emuMass_zz.at(ALL)->Integral(7,151) - emuMass_qcd.at(ALL)->Integral(7,151), (emuMass_zz.at(ALL)->Integral(7,151) - emuMass_qcd.at(ALL)->Integral(7,151)) * systErrMCLuEff[ZZ-1], 
  //       emuMass_zz.at(ALL)->Integral(13,151) - emuMass_qcd.at(ALL)->Integral(13,151), (emuMass_zz.at(ALL)->Integral(13,151) - emuMass_qcd.at(ALL)->Integral(13,151)) * systErrMCLuEff[ZZ-1], 
  //       emuMass_zz.at(ALL)->Integral(21,151) - emuMass_qcd.at(ALL)->Integral(21,151), (emuMass_zz.at(ALL)->Integral(21,151) - emuMass_qcd.at(ALL)->Integral(21,151)) * systErrMCLuEff[ZZ-1]); 
  cout << endl;
  cout << "---Without QCD correction:---------------------------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("TOT ttlike    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMass_ttbar.at(ALL)->Integral(7,151) - emuMass_wjets.at(ALL)->Integral(7,151), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, 7, 151, true),
         emuMass_ttbar.at(ALL)->Integral(13,151) - emuMass_wjets.at(ALL)->Integral(13,151), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, 13, 151, true),
         emuMass_ttbar.at(ALL)->Integral(21,151) - emuMass_wjets.at(ALL)->Integral(21,151), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, 21, 151, true));
  printf("TOT contam    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMass_wjets.at(ALL)->Integral(7,151) - emuMasses.at(QCD).at(ALL)->Integral(7,151), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, 7, 151, true),
         emuMass_wjets.at(ALL)->Integral(13,151) - emuMasses.at(QCD).at(ALL)->Integral(13,151), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, 13, 151, true),
         emuMass_wjets.at(ALL)->Integral(21,151) - emuMasses.at(QCD).at(ALL)->Integral(21,151), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, 21, 151, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;

  printf("TOT MC        | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMasses.at(1).at(ALL)->Integral(7,151) - emuMasses.at(QCD).at(ALL)->Integral(7,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, 7, 151, true),
         emuMasses.at(1).at(ALL)->Integral(13,151) - emuMasses.at(QCD).at(ALL)->Integral(13,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, 13, 151, true),
         emuMasses.at(1).at(ALL)->Integral(21,151) - emuMasses.at(QCD).at(ALL)->Integral(21,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, 21, 151, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << endl << endl << endl;


  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb LS DATA    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
         emuMasses.at(DATA).at(LS)->Integral(7,151), sqrt(emuMasses.at(DATA).at(LS)->Integral(7,151)),
         emuMasses.at(DATA).at(LS)->Integral(13,151), sqrt(emuMasses.at(DATA).at(LS)->Integral(13,151)),
         emuMasses.at(DATA).at(LS)->Integral(21,151), sqrt(emuMasses.at(DATA).at(LS)->Integral(21,151)));
  printf("nb LS MC      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMasses.at(1).at(LS)->Integral(7,151) - emuMasses.back().at(LS)->Integral(7,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, 7, 151, true),
         emuMasses.at(1).at(LS)->Integral(13,151) - emuMasses.back().at(LS)->Integral(13,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, 13, 151, true),
         emuMasses.at(1).at(LS)->Integral(21,151) - emuMasses.back().at(LS)->Integral(21,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, 21, 151, true));
  cout << endl;
  printf("nb OS DATA    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
         emuMasses.at(DATA).at(OS)->Integral(7,151), sqrt(emuMasses.at(DATA).at(OS)->Integral(7,151)),
         emuMasses.at(DATA).at(OS)->Integral(13,151), sqrt(emuMasses.at(DATA).at(OS)->Integral(13,151)),
         emuMasses.at(DATA).at(OS)->Integral(21,151), sqrt(emuMasses.at(DATA).at(OS)->Integral(21,151)));
  printf("nb OS MC      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMasses.at(1).at(OS)->Integral(7,151) - emuMasses.back().at(OS)->Integral(7,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, 7, 151, true),
         emuMasses.at(1).at(OS)->Integral(13,151) - emuMasses.back().at(OS)->Integral(13,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, 13, 151, true),
         emuMasses.at(1).at(OS)->Integral(21,151) - emuMasses.back().at(OS)->Integral(21,151), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, 21, 151, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;

  systErrMC.push_back(2 * sqrt(emuMasses.at(DATA).at(LS)->Integral() + pow(CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, 1, 151, true), 2)) / emuMasses.at(QCD).at(ALL)->Integral());
  systErrMCLuEff.push_back(systErrMC[QCD]);

  vector<bool> onlyQCD(8, false);
  onlyQCD.push_back(true);
  contamSamples.push_back(true);
  allSamples.push_back(true);

  cout << endl;
  cout << "---QCD events from LS spectrum:-----------------------------------------------------------------------------------------------------" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb QCD        | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) |\n",
         emuMasses.at(QCD).at(ALL)->Integral(7,151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 1, 151, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 7, 151, true) / emuMasses.at(QCD).at(ALL)->Integral(7,151),
         emuMasses.at(QCD).at(ALL)->Integral(13,151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 13, 151, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 13, 151, true) / emuMasses.at(QCD).at(ALL)->Integral(13,151),
         emuMasses.at(QCD).at(ALL)->Integral(21,151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 21, 151, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 21, 151, true) / emuMasses.at(QCD).at(ALL)->Integral(21,151));
  printf("%% of total MC |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |\n",
         100 * emuMasses.at(QCD).at(ALL)->Integral(7,151) / (emuMasses.at(1).at(ALL)->Integral(7,151) - emuMasses.at(QCD).at(ALL)->Integral(7,151)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 7, 151, true) / (emuMasses.at(1).at(ALL)->Integral(7,151) - emuMasses.at(QCD).at(ALL)->Integral(7,151)),
         100 * emuMasses.at(QCD).at(ALL)->Integral(13,151) / (emuMasses.at(1).at(ALL)->Integral(13,151) - emuMasses.at(QCD).at(ALL)->Integral(13,151)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 13, 151, true) / (emuMasses.at(1).at(ALL)->Integral(13,151) - emuMasses.at(QCD).at(ALL)->Integral(13,151)),
         100 * emuMasses.at(QCD).at(ALL)->Integral(21,151) / (emuMasses.at(1).at(ALL)->Integral(21,151) - emuMasses.at(QCD).at(ALL)->Integral(21,151)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, 21, 151, true) / (emuMasses.at(1).at(ALL)->Integral(21,151) - emuMasses.at(QCD).at(ALL)->Integral(21,151)));
  cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;

  cout << endl;
  cout << "--After adding QCD contribution:----------------------------------------------------------------------------------------------------------" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |         > 60GeV/c^2          |        > 120GeV/c^2          |         > 200GeV/c^2         |         > 500GeV/c^2         |" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
          emuMasses.at(DATA).at(ALL)->Integral(7, 151), sqrt(emuMasses.at(DATA).at(ALL)->Integral(7, 151)),
          emuMasses.at(DATA).at(ALL)->Integral(13, 151), sqrt((emuMasses.at(DATA).at(ALL))->Integral(13, 151)),
          emuMasses.at(DATA).at(ALL)->Integral(21, 151), sqrt(emuMasses.at(DATA).at(ALL)->Integral(21, 151)),
          emuMasses.at(DATA).at(ALL)->Integral(51, 151), sqrt(emuMasses.at(DATA).at(ALL)->Integral(51, 151)));
  printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(ALL)->Integral(7, 151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 7, 151, true),
          emuMasses.at(1).at(ALL)->Integral(13, 151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 13, 151, true),
          emuMasses.at(1).at(ALL)->Integral(21, 151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 21, 151, true),
          emuMasses.at(1).at(ALL)->Integral(51, 151), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 51, 151, true));
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl << endl;

  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |        60 - 120GeV/c^2       |      120 - 200GeV/c^2        |       200 - 400GeV/c^2       |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
          emuMasses.at(DATA).at(ALL)->Integral(7, 12), sqrt(emuMasses.at(DATA).at(ALL)->Integral(7, 12)),
          emuMasses.at(DATA).at(ALL)->Integral(13, 20), sqrt((emuMasses.at(DATA).at(ALL))->Integral(13, 20)),
          emuMasses.at(DATA).at(ALL)->Integral(21, 40), sqrt(emuMasses.at(DATA).at(ALL)->Integral(21, 40)));
  printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(ALL)->Integral(7, 12), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 7, 12, true),
          emuMasses.at(1).at(ALL)->Integral(13, 20), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 13, 20, true),
          emuMasses.at(1).at(ALL)->Integral(21, 40), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, 21, 40, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data OS    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
          emuMasses.at(DATA).at(OS)->Integral(7, 12), sqrt(emuMasses.at(DATA).at(OS)->Integral(7, 12)),
          emuMasses.at(DATA).at(OS)->Integral(13, 20), sqrt((emuMasses.at(DATA).at(OS))->Integral(13, 20)),
          emuMasses.at(DATA).at(OS)->Integral(21, 40), sqrt(emuMasses.at(DATA).at(OS)->Integral(21, 40)));
  printf("nb MC OS      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(OS)->Integral(7, 12), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, 7, 12, true),
          emuMasses.at(1).at(OS)->Integral(13, 20), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, 13, 20, true),
          emuMasses.at(1).at(OS)->Integral(21, 40), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, 21, 40, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl << endl;
}

float CalcSystErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked)
{
   float err2 = 0.;

   // for stacked histograms subtract all contributions below the current one 
   if (stacked) {
      for (unsigned int i = 0; i < samples.size() && i < histos.size()-1; ++i) {
         if (samples[i]) {
            float numEv;
            if (i < histos.size()-2) numEv = histos.at(i+1).at(region)->Integral(lowerBin, upperBin) - histos.at(i+2).at(region)->Integral(lowerBin, upperBin);
            else numEv = histos.at(i+1).at(region)->Integral(lowerBin, upperBin);
            err2 += pow(numEv, 2) * errors[i] * errors[i];
         }
      }
   } else {
      for (unsigned int i = 0; i < samples.size(); ++i) {
         if (samples[i])
            err2 += histos.at(i+1).at(region)->Integral(lowerBin, upperBin) * histos.at(i+1).at(region)->Integral(lowerBin, upperBin) * errors[i] * errors[i];
      }
   }

   return sqrt(err2);
}

float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked)
{
   float qcdErr = 0.;
   vector<bool> allButQCD(8, true);

   qcdErr = 2 * sqrt(histos.at(DATA).at(LS)->Integral(lowerBin, upperBin) + pow(CalcSystErr(histos, errors, allButQCD, LS, lowerBin, upperBin, stacked), 2));
   errors.back() = qcdErr / histos.at(QCD).at(region)->Integral(lowerBin, upperBin);

   return CalcSystErr(histos, errors, samples, region, lowerBin, upperBin, stacked);
}

