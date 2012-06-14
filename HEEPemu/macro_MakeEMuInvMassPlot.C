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
float CalcAllErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples,  int region, int lowerBin, int upperBin, bool stacked = false);
float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked = false);
float CalcAllErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int sample, int region, int lowerBin, int upperBin, bool stacked = false);

void macro_MakeEMuInvMassPlot()
{
  // parameters //////////////////////////////////////////////////////////////
  //TFile input("testEmuSpecHEEP32_291pb-1.root", "open");
  //TFile input("testEmuSpecHEEP4strong_291pb-1.root", "open");
  //TFile input("testEmuSpecHEEP4light_816pb-1.root", "open");
  //TFile input("testEmuSpecHEEP4_920pb-1.root", "open");
  TFile input("testEmuSpecHEEP4_3692pb-1.root", "open");
  //TFile input("testEmuSpecHEEP4_840pb-1_noTrg.root", "open");
  //TFile input("testHEEP4_525pb-1.root", "open");

  const float lumi = 3692.;
  const float minInvMass = 0.;

  int nBins = 75;

  bool plotSign[3];
  plotSign[0] = true;  // all
  plotSign[1] = true;  // LS like sign
  plotSign[2] = true;  // OS opposite sign

  bool plotType[2];
  plotType[0] = true;  // emu spectrum
  plotType[1] = true;  // cumulative emu spectrum

  bool logPlot = true;
  bool prelim = true;
  bool groupedPlot = false;

  // systematical errors
  vector<float> systErrMC;
  systErrMC.push_back(0.15);  //ttbar
  systErrMC.push_back(0.054); //z->tt
  systErrMC.push_back(0.035); //WWE
  systErrMC.push_back(0.038); //WZ
  systErrMC.push_back(0.075); //tW
  systErrMC.push_back(0.05);  //WJets
  systErrMC.push_back(0.054); //Z->mm
  systErrMC.push_back(0.054); //Z->ee
  //  systErrMC.push_back(0.025); //ZZ
  float systErrLumi = 0.022;
  float systErrEff = 0.008; // muon err & ele err
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString histoSign[3] = {"", "LS_", "OS_"};
  TString histoTitleSign[3] = {"E", "E", "E"};
  TString xAxisTitle[3] = {"m(e#mu)", "m(e^{#pm}#mu^{#pm})", "m(e^{#pm}#mu^{#mp})"};
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
//  //vector<TH1F *> emuMass_zz;
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
//      //emuMass_zz.push_back((TH1F *)gDirectory->Get("emuMass_" + histoSign[k] + "zz"));
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
//      //emuMass_zz.back()->SetName(emuMass_zz.back()->GetName() + nameSuffix[j]);
      emuMass_qcd.back()->SetName(emuMass_qcd.back()->GetName() + nameSuffix[j]);

      if (!plotSign[k]) continue;

      // integrate from the right side
      if (j == 1) { 
        // loop over bins
        double error;
        for (int i = 1; i < nBins + 1; ++i) {
          emuMass_data.back()->SetBinContent(i, emuMass_data.back()->IntegralAndError(i, nBins, error));
          emuMass_data.back()->SetBinError(i, error);
          emuMass_ttbar.back()->SetBinContent(i, emuMass_ttbar.back()->IntegralAndError(i, nBins, error));
          emuMass_ttbar.back()->SetBinError(i, error);
          emuMass_ztautau.back()->SetBinContent(i, emuMass_ztautau.back()->IntegralAndError(i, nBins, error));
          emuMass_ztautau.back()->SetBinError(i, error);
          emuMass_ww.back()->SetBinContent(i, emuMass_ww.back()->IntegralAndError(i, nBins, error));
          emuMass_ww.back()->SetBinError(i, error);
          emuMass_wz.back()->SetBinContent(i, emuMass_wz.back()->IntegralAndError(i, nBins, error));
          emuMass_wz.back()->SetBinError(i, error);
          emuMass_tw.back()->SetBinContent(i, emuMass_tw.back()->IntegralAndError(i, nBins, error));
          emuMass_tw.back()->SetBinError(i, error);
          emuMass_wjets.back()->SetBinContent(i, emuMass_wjets.back()->IntegralAndError(i, nBins, error));
          emuMass_wjets.back()->SetBinError(i, error);
          emuMass_zmumu.back()->SetBinContent(i, emuMass_zmumu.back()->IntegralAndError(i, nBins, error));
          emuMass_zmumu.back()->SetBinError(i, error);
          emuMass_zee.back()->SetBinContent(i, emuMass_zee.back()->IntegralAndError(i, nBins, error));
          emuMass_zee.back()->SetBinError(i, error);
//          //emuMass_zz.back()->SetBinContent(i, emuMass_zz.back()->IntegralAndError(i, nBins, error));
//          //emuMass_zz.back()->SetBinError(i, error);
          emuMass_qcd.back()->SetBinContent(i, emuMass_qcd.back()->IntegralAndError(i, nBins, error));
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
      if (groupedPlot) {
        emuPlot->SetLeftMargin(0.11);
        emuPlot->SetRightMargin(0.09);
        emuPlot->SetBottomMargin(0.12);
        emuPlot->SetTopMargin(0.08);
      }
 
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      if (!groupedPlot) {
        gStyle->SetTitleXOffset(1.);
        gStyle->SetTitleYOffset(1.3);
      }
      gPad->SetTicks(1, 1);

      emuMass_ttbar.back()->GetXaxis()->SetTitle(xAxisTitle[k] + " (GeV)");
      emuMass_ttbar.back()->GetXaxis()->SetTitleSize(0.04);
      emuMass_ttbar.back()->GetXaxis()->SetLabelSize(0.04);
      if (groupedPlot) {
        emuMass_ttbar.back()->GetXaxis()->SetTitleSize(0.05);
        emuMass_ttbar.back()->GetXaxis()->SetLabelSize(0.05);
      }
      //emuMass_ttbar.back()->GetXaxis()->SetRangeUser(60., 1110.);
    
      if (j == 1) emuMass_ttbar.back()->GetYaxis()->SetTitle(histoTitleSign[k] + "vents #geq " + xAxisTitle[k]);
      else emuMass_ttbar.back()->GetYaxis()->SetTitle(histoTitleSign[k] + "vents / 20 GeV");
      emuMass_ttbar.back()->GetYaxis()->SetTitleSize(0.05);
      if (!groupedPlot) {
        emuMass_ttbar.back()->GetYaxis()->SetTitleSize(0.04);
        emuMass_ttbar.back()->GetYaxis()->SetTitleOffset(1.2);
      }
    
      if (emuMass_data.back()->GetMaximum() > emuMass_ttbar.back()->GetMaximum()) {
        if (!logPlot)  emuMass_ttbar.back()->SetMaximum(emuMass_data.back()->GetMaximum() * 1.1);
        else emuMass_ttbar.back()->SetMaximum(emuMass_data.back()->GetMaximum() * 1.3);
      }
    
      // plot spectrum
      //emuMass_ttbar.back()->SetFillColor(kRed-9);
      emuMass_ttbar.back()->SetFillColor(kRed);
      //emuMass_ttbar.back()->SetLineColor(kRed);
      //emuMass_ttbar.back()->SetMarkerColor(kRed);
      //emuMass_ttbar.back()->SetLineWidth(2);
      emuMass_ttbar.back()->DrawClone("HIST");
      //emuMass_ttbar.back()->Draw("e0sames");
      if (groupedPlot) {
        emuMass_wjets.back()->SetFillColor(kYellow);
        //emuMass_wjets.back()->SetLineColor(kYellow);
        emuMass_wjets.back()->SetMarkerColor(kYellow);
        emuMass_wjets.back()->Draw("HISTsames");
      } else {
        emuMass_ztautau.back()->SetFillColor(kRed-7);
        //emuMass_ztautau.back()->SetLineColor(kRed-7);
        emuMass_ztautau.back()->SetMarkerColor(kRed-7);
        emuMass_ztautau.back()->Draw("HISTsames");
        emuMass_ww.back()->SetFillColor(kRed-4);
        //emuMass_ww.back()->SetLineColor(kRed-4);
        emuMass_ww.back()->SetMarkerColor(kRed-4);
        emuMass_ww.back()->Draw("HISTsames");
        emuMass_wz.back()->SetFillColor(kRed+1);
        //emuMass_wz.back()->SetLineColor(kRed+1);
        emuMass_wz.back()->SetMarkerColor(kRed+1);
        emuMass_wz.back()->Draw("HISTsames");
        emuMass_tw.back()->SetFillColor(kRed+2);
        //emuMass_tw.back()->SetLineColor(kRed+2);
        emuMass_tw.back()->SetMarkerColor(kRed+2);
        emuMass_tw.back()->Draw("HISTsames");
        emuMass_wjets.back()->SetFillColor(kCyan);
        //emuMass_wjets.back()->SetLineColor(kCyan);
        emuMass_wjets.back()->SetMarkerColor(kCyan);
        emuMass_wjets.back()->Draw("HISTsames");
        emuMass_zmumu.back()->SetFillColor(kCyan+1);
        //emuMass_zmumu.back()->SetLineColor(kCyan+1);
        emuMass_zmumu.back()->SetMarkerColor(kCyan+1);
        emuMass_zmumu.back()->Draw("HISTsames");
        emuMass_zee.back()->SetFillColor(kCyan+2);
        //emuMass_zee.back()->SetLineColor(kCyan+2);
        emuMass_zee.back()->SetMarkerColor(kCyan+2);
        emuMass_zee.back()->Draw("HISTsames");
        //emuMass_zz.back()->Draw("HISTsames");
        emuMass_qcd.back()->SetFillColor(kYellow);
        //emuMass_qcd.back()->SetLineColor(kYellow);
        emuMass_qcd.back()->SetMarkerColor(kYellow);
        emuMass_qcd.back()->Draw("HISTsames");
      }
      
      emuMass_data.back()->SetLineWidth(2);
      emuMass_data.back()->SetMarkerStyle(8);
      emuMass_data.back()->SetMarkerSize(0.8);
      emuMass_data.back()->Draw("sames");

      // redraw axis
      emuMass_ttbar.back()->Draw("sameaxis");

      // legend and labels
      TLegend legend(0.69, 0.37, 0.85, 0.76);
      legend.SetTextSize(0.03);
      legend.SetFillColor(0);
      if (groupedPlot) {
        legend.SetX1(0.44);
        legend.SetY1(0.60);
        legend.SetX2(0.7);
        legend.SetY2(0.79);
        legend.SetTextSize(0.045);
      }
    
      legend.AddEntry(emuMass_data.back(), "DATA");
      if (groupedPlot) {
        legend.AddEntry(emuMass_ttbar.back(), "t#bar{t} & other prompt leptons" ,"F");
        legend.AddEntry(emuMass_wjets.back(), "fake e#mu pairs" ,"F");
      } else {
        legend.AddEntry(emuMass_ttbar.back(), "t#bar{t} (MC)" ,"F");
        legend.AddEntry(emuMass_ztautau.back(), "Z/#gamma* #rightarrow #tau#tau (MC)" ,"F");
        legend.AddEntry(emuMass_ww.back(), "WW (MC)" ,"F");
        legend.AddEntry(emuMass_wz.back(), "WZ (MC)" ,"F");
        legend.AddEntry(emuMass_tw.back(), "tW (MC)" ,"F");
        legend.AddEntry(emuMass_wjets.back(), "W+jets (MC)" ,"F");
        legend.AddEntry(emuMass_zmumu.back(), "Z/#gamma* #rightarrow #mu#mu (MC)" ,"F");
        legend.AddEntry(emuMass_zee.back(), "Z/#gamma* #rightarrow ee (MC)" ,"F");
        //legend.AddEntry(emuMass_zz.back(), "ZZ, (MC)" ,"F");
        legend.AddEntry(emuMass_qcd.back(), "QCD" ,"F");
      }
    
      legend.SetBorderSize(0);
      legend.DrawClone("sames");
      
      stringstream sStream;
      sStream.str("");
      //if (prelim) sStream << "CMS preliminary    #sqrt{s} = 8 TeV    #int L dt = " << lumi << " pb^{-1}";
      if (prelim) sStream << "CMS preliminary    #sqrt{s} = 8 TeV    #int L dt = 3.7 fb^{-1}";
      else sStream << "CMS     #sqrt{s} = 8 TeV    #int L dt = " << lumi << " pb^{-1}";
      //sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
      TPaveLabel labelLumi(0.32, 0.78, 0.88, 0.87, sStream.str().c_str(), "brNDC");
      labelLumi.SetFillColor(0);
      labelLumi.SetFillStyle(0);
      labelLumi.SetBorderSize(0);
      labelLumi.SetTextSize(0.4);
      if (groupedPlot) {
        labelLumi.SetX1(0.22);
        labelLumi.SetTextSize(0.5);
      }
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
//  //emuMasses.push_back(emuMass_zz);
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

  int bin60 = (int)60 / 1500. * nBins + 1;
  int bin120 = (int)120 / 1500. * nBins + 1;
  int bin200 = (int)200 / 1500. * nBins + 1;
  int bin400 = (int)400 / 1500. * nBins + 1;
  int bin500 = (int)500 / 1500. * nBins + 1;

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
  cout << " W+Jets:      " << systErrMC[WJET-1] * 100 << "%" << endl;
  cout << " Z->mumu:     " << systErrMC[ZMM-1] * 100 << "%" << endl;
  cout << " Z->ee:       " << systErrMC[ZEE-1] * 100 << "%" << endl;
//  cout << " ZZ:          " << systErrMC[ZZ-1] * 100 << "%" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |         >  60GeV/c^2          |        > 120GeV/c^2          |        > 200GeV/c^2         |        > 400GeV/c^2          |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n", 
         emuMass_data.at(ALL)->Integral(bin60,nBins + 1), sqrt(emuMass_data.at(ALL)->Integral(bin60,nBins + 1)),         
         emuMass_data.at(ALL)->Integral(bin120,nBins + 1), sqrt(emuMass_data.at(ALL)->Integral(bin120,nBins + 1)),
         emuMass_data.at(ALL)->Integral(bin200,nBins + 1), sqrt(emuMass_data.at(ALL)->Integral(bin200,nBins + 1)),
         emuMass_data.at(ALL)->Integral(bin400,nBins + 1), sqrt(emuMass_data.at(ALL)->Integral(bin400,nBins + 1)));
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb ttbar      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_ttbar.at(ALL)->Integral(bin60,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin60,nBins + 1), (emuMass_ttbar.at(ALL)->Integral(bin60,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[TTBAR-1], 
         emuMass_ttbar.at(ALL)->Integral(bin120,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin120,nBins + 1), (emuMass_ttbar.at(ALL)->Integral(bin120,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[TTBAR-1], 
         emuMass_ttbar.at(ALL)->Integral(bin200,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin200,nBins + 1), (emuMass_ttbar.at(ALL)->Integral(bin200,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[TTBAR-1],
         emuMass_ttbar.at(ALL)->Integral(bin400,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin400,nBins + 1), (emuMass_ttbar.at(ALL)->Integral(bin400,nBins + 1) - emuMass_ztautau.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[TTBAR-1]);
  printf("nb Ztautau    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_ztautau.at(ALL)->Integral(bin60,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin60,nBins + 1), (emuMass_ztautau.at(ALL)->Integral(bin60,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[ZTT-1], 
         emuMass_ztautau.at(ALL)->Integral(bin120,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin120,nBins + 1), (emuMass_ztautau.at(ALL)->Integral(bin120,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[ZTT-1], 
         emuMass_ztautau.at(ALL)->Integral(bin200,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin200,nBins + 1), (emuMass_ztautau.at(ALL)->Integral(bin200,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[ZTT-1], 
         emuMass_ztautau.at(ALL)->Integral(bin400,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin400,nBins + 1), (emuMass_ztautau.at(ALL)->Integral(bin400,nBins + 1) - emuMass_ww.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[ZTT-1]); 
  printf("nb WW         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_ww.at(ALL)->Integral(bin60,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin60,nBins + 1), (emuMass_ww.at(ALL)->Integral(bin60,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[WW-1], 
         emuMass_ww.at(ALL)->Integral(bin120,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin120,nBins + 1), (emuMass_ww.at(ALL)->Integral(bin120,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[WW-1], 
         emuMass_ww.at(ALL)->Integral(bin200,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin200,nBins + 1), (emuMass_ww.at(ALL)->Integral(bin200,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[WW-1],
         emuMass_ww.at(ALL)->Integral(bin400,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin400,nBins + 1), (emuMass_ww.at(ALL)->Integral(bin400,nBins + 1) - emuMass_wz.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[WW-1]);
  printf("nb WZ         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_wz.at(ALL)->Integral(bin60,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin60,nBins + 1), (emuMass_wz.at(ALL)->Integral(bin60,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[WZ-1], 
         emuMass_wz.at(ALL)->Integral(bin120,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin120,nBins + 1), (emuMass_wz.at(ALL)->Integral(bin120,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[WZ-1], 
         emuMass_wz.at(ALL)->Integral(bin200,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin200,nBins + 1), (emuMass_wz.at(ALL)->Integral(bin200,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[WZ-1],
         emuMass_wz.at(ALL)->Integral(bin400,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin400,nBins + 1), (emuMass_wz.at(ALL)->Integral(bin400,nBins + 1) - emuMass_tw.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[WZ-1]);
  printf("nb tW         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_tw.at(ALL)->Integral(bin60,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin60,nBins + 1), (emuMass_tw.at(ALL)->Integral(bin60,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[TW-1], 
         emuMass_tw.at(ALL)->Integral(bin120,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin120,nBins + 1), (emuMass_tw.at(ALL)->Integral(bin120,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[TW-1], 
         emuMass_tw.at(ALL)->Integral(bin200,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin200,nBins + 1), (emuMass_tw.at(ALL)->Integral(bin200,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[TW-1],
         emuMass_tw.at(ALL)->Integral(bin400,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin400,nBins + 1), (emuMass_tw.at(ALL)->Integral(bin400,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[TW-1]);
  cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
  printf("nb WJets      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_wjets.at(ALL)->Integral(bin60,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin60,nBins + 1), (emuMass_wjets.at(ALL)->Integral(bin60,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[WJET-1], 
         emuMass_wjets.at(ALL)->Integral(bin120,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin120,nBins + 1), (emuMass_wjets.at(ALL)->Integral(bin120,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[WJET-1], 
         emuMass_wjets.at(ALL)->Integral(bin200,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin200,nBins + 1), (emuMass_wjets.at(ALL)->Integral(bin200,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[WJET-1],
         emuMass_wjets.at(ALL)->Integral(bin400,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin400,nBins + 1), (emuMass_wjets.at(ALL)->Integral(bin400,nBins + 1) - emuMass_zmumu.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[WJET-1]);
  printf("nb Zmumu      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_zmumu.at(ALL)->Integral(bin60,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin60,nBins + 1), (emuMass_zmumu.at(ALL)->Integral(bin60,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[ZMM-1], 
         emuMass_zmumu.at(ALL)->Integral(bin120,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin120,nBins + 1), (emuMass_zmumu.at(ALL)->Integral(bin120,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[ZMM-1], 
         emuMass_zmumu.at(ALL)->Integral(bin200,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin200,nBins + 1), (emuMass_zmumu.at(ALL)->Integral(bin200,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[ZMM-1],
         emuMass_zmumu.at(ALL)->Integral(bin400,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin400,nBins + 1), (emuMass_zmumu.at(ALL)->Integral(bin400,nBins + 1) - emuMass_zee.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[ZMM-1]);
  printf("nb Zee        | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
         emuMass_zee.at(ALL)->Integral(bin60,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin60,nBins + 1), (emuMass_zee.at(ALL)->Integral(bin60,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[ZEE-1], 
         emuMass_zee.at(ALL)->Integral(bin120,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin120,nBins + 1), (emuMass_zee.at(ALL)->Integral(bin120,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[ZEE-1], 
         emuMass_zee.at(ALL)->Integral(bin200,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin200,nBins + 1), (emuMass_zee.at(ALL)->Integral(bin200,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[ZEE-1], 
         emuMass_zee.at(ALL)->Integral(bin400,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin400,nBins + 1), (emuMass_zee.at(ALL)->Integral(bin400,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[ZEE-1]); 
//  //printf("nb ZZ         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n", 
//  //       emuMass_zz.at(ALL)->Integral(bin60,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin60,nBins + 1), (emuMass_zz.at(ALL)->Integral(bin60,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin60,nBins + 1)) * systErrMCLuEff[ZZ-1], 
//  //       emuMass_zz.at(ALL)->Integral(bin120,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin120,nBins + 1), (emuMass_zz.at(ALL)->Integral(bin120,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin120,nBins + 1)) * systErrMCLuEff[ZZ-1], 
//  //       emuMass_zz.at(ALL)->Integral(bin200,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin200,nBins + 1), (emuMass_zz.at(ALL)->Integral(bin200,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin200,nBins + 1)) * systErrMCLuEff[ZZ-1], 
//  //       emuMass_zz.at(ALL)->Integral(bin400,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin400,nBins + 1), (emuMass_zz.at(ALL)->Integral(bin400,nBins + 1) - emuMass_qcd.at(ALL)->Integral(bin400,nBins + 1)) * systErrMCLuEff[ZZ-1]); 
  cout << endl;
  cout << "---Without QCD correction:--------------------------------------------------------------------------------------------------------------" << endl;
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("TOT ttlike    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMass_ttbar.at(ALL)->Integral(bin60,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin60,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin60, nBins + 1, true),
         emuMass_ttbar.at(ALL)->Integral(bin120,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin120,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin120, nBins + 1, true),
         emuMass_ttbar.at(ALL)->Integral(bin200,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin200,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin200, nBins + 1, true),
         emuMass_ttbar.at(ALL)->Integral(bin400,nBins + 1) - emuMass_wjets.at(ALL)->Integral(bin400,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin400, nBins + 1, true));
  printf("TOT contam    | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMass_wjets.at(ALL)->Integral(bin60,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin60, nBins + 1, true),
         emuMass_wjets.at(ALL)->Integral(bin120,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin120, nBins + 1, true),
         emuMass_wjets.at(ALL)->Integral(bin200,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin200, nBins + 1, true),
         emuMass_wjets.at(ALL)->Integral(bin400,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin400, nBins + 1, true));
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;

  printf("TOT MC        | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMasses.at(1).at(ALL)->Integral(bin60,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, bin60, nBins + 1, true),
         emuMasses.at(1).at(ALL)->Integral(bin120,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, bin120, nBins + 1, true),
         emuMasses.at(1).at(ALL)->Integral(bin200,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, bin200, nBins + 1, true),
         emuMasses.at(1).at(ALL)->Integral(bin400,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, ALL, bin400, nBins + 1, true));
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << endl << endl << endl;


  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb LS DATA    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        | %5.0f +- %-.3f (stat)        |\n",
         emuMasses.at(DATA).at(LS)->Integral(bin60,nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin60,nBins + 1)),
         emuMasses.at(DATA).at(LS)->Integral(bin120,nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin120,nBins + 1)),
         emuMasses.at(DATA).at(LS)->Integral(bin200,nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin200,nBins + 1)),
         emuMasses.at(DATA).at(LS)->Integral(bin400,nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin400,nBins + 1)));
  printf("nb LS MC      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMasses.at(1).at(LS)->Integral(bin60,nBins + 1) - emuMasses.back().at(LS)->Integral(bin60,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, bin60, nBins + 1, true),
         emuMasses.at(1).at(LS)->Integral(bin120,nBins + 1) - emuMasses.back().at(LS)->Integral(bin120,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, bin120, nBins + 1, true),
         emuMasses.at(1).at(LS)->Integral(bin200,nBins + 1) - emuMasses.back().at(LS)->Integral(bin200,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, bin200, nBins + 1, true),
         emuMasses.at(1).at(LS)->Integral(bin400,nBins + 1) - emuMasses.back().at(LS)->Integral(bin400,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, bin400, nBins + 1, true));
  cout << endl;
  printf("nb OS DATA    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
         emuMasses.at(DATA).at(OS)->Integral(bin60,nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin60,nBins + 1)),
         emuMasses.at(DATA).at(OS)->Integral(bin120,nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin120,nBins + 1)),
         emuMasses.at(DATA).at(OS)->Integral(bin200,nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin200,nBins + 1)),
         emuMasses.at(DATA).at(OS)->Integral(bin400,nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin400,nBins + 1)));
  printf("nb OS MC      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
         emuMasses.at(1).at(OS)->Integral(bin60,nBins + 1) - emuMasses.back().at(OS)->Integral(bin60,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, bin60, nBins + 1, true),
         emuMasses.at(1).at(OS)->Integral(bin120,nBins + 1) - emuMasses.back().at(OS)->Integral(bin120,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, bin120, nBins + 1, true),
         emuMasses.at(1).at(OS)->Integral(bin200,nBins + 1) - emuMasses.back().at(OS)->Integral(bin200,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, bin200, nBins + 1, true),
         emuMasses.at(1).at(OS)->Integral(bin400,nBins + 1) - emuMasses.back().at(OS)->Integral(bin400,nBins + 1), CalcSystErr(emuMasses, systErrMCLuEff, allSamples, OS, bin400, nBins + 1, true));
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;

  systErrMC.push_back(2 * sqrt(emuMasses.at(DATA).at(LS)->Integral() + pow(CalcSystErr(emuMasses, systErrMCLuEff, allSamples, LS, 1, nBins + 1, true), 2)) / emuMasses.at(QCD).at(ALL)->Integral());
  systErrMCLuEff.push_back(systErrMC[QCD]);

  vector<bool> onlyQCD(8, false);
  onlyQCD.push_back(true);
  contamSamples.push_back(true);
  allSamples.push_back(true);

  cout << endl;
  cout << "---QCD events from LS spectrum:----------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb QCD LS+OS  | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) | %9.3f +- %8.3f (%.1f%%) (syst) |\n",
         emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin60, nBins + 1, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin60, nBins + 1, true) / emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1),
         emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin120, nBins + 1, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin120, nBins + 1, true) / emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1),
         emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin200, nBins + 1, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin200, nBins + 1, true) / emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1),
         emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin400, nBins + 1, true), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin400, nBins + 1, true) / emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1));
  printf("%% of total MC |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |  %7.3f%% +- %7.3f%% (syst)         |\n",
         100 * emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1) / (emuMasses.at(1).at(ALL)->Integral(bin60,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin60, nBins + 1, true) / (emuMasses.at(1).at(ALL)->Integral(bin60,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin60,nBins + 1)),
         100 * emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1) / (emuMasses.at(1).at(ALL)->Integral(bin120,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin120, nBins + 1, true) / (emuMasses.at(1).at(ALL)->Integral(bin120,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin120,nBins + 1)),
         100 * emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1) / (emuMasses.at(1).at(ALL)->Integral(bin200,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin200, nBins + 1, true) / (emuMasses.at(1).at(ALL)->Integral(bin200,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin200,nBins + 1)),
         100 * emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1) / (emuMasses.at(1).at(ALL)->Integral(bin400,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1)), 100 * CalcSystErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, ALL, bin400, nBins + 1, true) / (emuMasses.at(1).at(ALL)->Integral(bin400,nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin400,nBins + 1)));
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

  cout << endl;
  cout << "--After adding QCD contribution:----------------------------------------------------------------------------------------------------------" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |         > 60GeV/c^2          |        > 120GeV/c^2          |         > 200GeV/c^2         |         > 400GeV/c^2         |" << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
          emuMasses.at(DATA).at(ALL)->Integral(bin60, nBins + 1), sqrt(emuMasses.at(DATA).at(ALL)->Integral(bin60, nBins + 1)),
          emuMasses.at(DATA).at(ALL)->Integral(bin120, nBins + 1), sqrt((emuMasses.at(DATA).at(ALL))->Integral(bin120, nBins + 1)),
          emuMasses.at(DATA).at(ALL)->Integral(bin200, nBins + 1), sqrt(emuMasses.at(DATA).at(ALL)->Integral(bin200, nBins + 1)),
          emuMasses.at(DATA).at(ALL)->Integral(bin400, nBins + 1), sqrt(emuMasses.at(DATA).at(ALL)->Integral(bin400, nBins + 1)));
  printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(ALL)->Integral(bin60, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin60, nBins + 1, true),
          emuMasses.at(1).at(ALL)->Integral(bin120, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin120, nBins + 1, true),
          emuMasses.at(1).at(ALL)->Integral(bin200, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin200, nBins + 1, true),
          emuMasses.at(1).at(ALL)->Integral(bin400, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin400, nBins + 1, true));
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data LS    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
          emuMasses.at(DATA).at(LS)->Integral(bin60, nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin60, nBins + 1)),
          emuMasses.at(DATA).at(LS)->Integral(bin120, nBins + 1), sqrt((emuMasses.at(DATA).at(LS))->Integral(bin120, nBins + 1)),
          emuMasses.at(DATA).at(LS)->Integral(bin200, nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin200, nBins + 1)),
          emuMasses.at(DATA).at(LS)->Integral(bin400, nBins + 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin400, nBins + 1)));
  printf("nb MC LS      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(LS)->Integral(bin60, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin60, nBins + 1, true),
          emuMasses.at(1).at(LS)->Integral(bin120, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin120, nBins + 1, true),
          emuMasses.at(1).at(LS)->Integral(bin200, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin200, nBins + 1, true),
          emuMasses.at(1).at(LS)->Integral(bin400, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin400, nBins + 1, true));
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data OS    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)        |\n",
          emuMasses.at(DATA).at(OS)->Integral(bin60, nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin60, nBins + 1)),
          emuMasses.at(DATA).at(OS)->Integral(bin120, nBins + 1), sqrt((emuMasses.at(DATA).at(OS))->Integral(bin120, nBins + 1)),
          emuMasses.at(DATA).at(OS)->Integral(bin200, nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin200, nBins + 1)),
          emuMasses.at(DATA).at(OS)->Integral(bin400, nBins + 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin400, nBins + 1)));
  printf("nb MC OS      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(OS)->Integral(bin60, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin60, nBins + 1, true),
          emuMasses.at(1).at(OS)->Integral(bin120, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin120, nBins + 1, true),
          emuMasses.at(1).at(OS)->Integral(bin200, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin200, nBins + 1, true),
          emuMasses.at(1).at(OS)->Integral(bin400, nBins + 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin400, nBins + 1, true));
  cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl << endl;



  cout << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  cout << "M_emu         |        60 - 120GeV/c^2       |      120 - 200GeV/c^2        |       200 - 400GeV/c^2       |" << endl;
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
          emuMasses.at(DATA).at(ALL)->Integral(bin60, bin120 - 1), sqrt(emuMasses.at(DATA).at(ALL)->Integral(bin60, bin120 - 1)),
          emuMasses.at(DATA).at(ALL)->Integral(bin120, bin200 - 1), sqrt((emuMasses.at(DATA).at(ALL))->Integral(bin120, bin200 - 1)),
          emuMasses.at(DATA).at(ALL)->Integral(bin200, bin400 - 1), sqrt(emuMasses.at(DATA).at(ALL)->Integral(bin200, bin400 - 1)));
  printf("nb MC         | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(ALL)->Integral(bin60, bin120 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin60, bin120 - 1, true),
          emuMasses.at(1).at(ALL)->Integral(bin120, bin200 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin120, bin200 - 1, true),
          emuMasses.at(1).at(ALL)->Integral(bin200, bin400 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, ALL, bin200, bin400 - 1, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data LS    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
          emuMasses.at(DATA).at(LS)->Integral(bin60, bin120 - 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin60, bin120 - 1)),
          emuMasses.at(DATA).at(LS)->Integral(bin120, bin200 - 1), sqrt((emuMasses.at(DATA).at(LS))->Integral(bin120, bin200 - 1)),
          emuMasses.at(DATA).at(LS)->Integral(bin200, bin400 - 1), sqrt(emuMasses.at(DATA).at(LS)->Integral(bin200, bin400 - 1)));
  printf("nb MC LS      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(LS)->Integral(bin60, bin120 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin60, bin120 - 1, true),
          emuMasses.at(1).at(LS)->Integral(bin120, bin200 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin120, bin200 - 1, true),
          emuMasses.at(1).at(LS)->Integral(bin200, bin400 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, LS, bin200, bin400 - 1, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl;
  printf("nb data OS    | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       | %5.0f +- %-.3f (stat)       |\n",
          emuMasses.at(DATA).at(OS)->Integral(bin60, bin120 - 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin60, bin120 - 1)),
          emuMasses.at(DATA).at(OS)->Integral(bin120, bin200 - 1), sqrt((emuMasses.at(DATA).at(OS))->Integral(bin120, bin200 - 1)),
          emuMasses.at(DATA).at(OS)->Integral(bin200, bin400 - 1), sqrt(emuMasses.at(DATA).at(OS)->Integral(bin200, bin400 - 1)));
  printf("nb MC OS      | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) | %9.3f +- %8.3f (syst) |\n",
          emuMasses.at(1).at(OS)->Integral(bin60, bin120 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin60, bin120 - 1, true),
          emuMasses.at(1).at(OS)->Integral(bin120, bin200 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin120, bin200 - 1, true),
          emuMasses.at(1).at(OS)->Integral(bin200, bin400 - 1), CalcSystErrWithQCD(emuMasses, systErrMCLuEff, allSamples, OS, bin200, bin400 - 1, true));
  cout << "-----------------------------------------------------------------------------------------------------------" << endl << endl;

  cout << endl;
  cout << "\\begin{table}[b]" << endl;
  cout << "\\centering" << endl;
  cout << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" << endl;
  cout << "\\hline\\hline" << endl;
  cout << "Source & \\multicolumn{9}{c|}{number of events} \\\\" << endl;
  cout << " &  \\multicolumn{3}{c|}{[$120-200$]~GeV/c$^2$}  & \\multicolumn{3}{c|}{[$200-400$]~GeV/c$^2$} & \\multicolumn{3}{|c|}{$>$ 400~GeV/c$^2$}  \\\\" << endl;
  cout << " &  OS & LS & Combined & OS & LS & Combined & OS & LS & Combined \\\\\\hline" << endl;

  //printf("          & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"    
  printf("CMS data  & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f \\\\\n", 
emuMasses.at(DATA).at(OS)->Integral(bin120, bin200 - 1), emuMasses.at(DATA).at(LS)->Integral(bin120, bin200 - 1), emuMasses.at(DATA).at(ALL)->Integral(bin120, bin200 - 1), 
emuMasses.at(DATA).at(OS)->Integral(bin200, bin400 - 1), emuMasses.at(DATA).at(LS)->Integral(bin200, bin400 - 1), emuMasses.at(DATA).at(ALL)->Integral(bin200, bin400 - 1), 
emuMasses.at(DATA).at(OS)->Integral(bin400, nBins + 1), emuMasses.at(DATA).at(LS)->Integral(bin400, nBins + 1), emuMasses.at(DATA).at(ALL)->Integral(bin400, nBins + 1));
  printf("Total Bkg & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         emuMasses.at(1).at(OS)->Integral(bin120, bin200 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, OS, bin120, bin200 - 1, true), 
         emuMasses.at(1).at(LS)->Integral(bin120, bin200 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, LS, bin120, bin200 - 1, true), 
         emuMasses.at(1).at(ALL)->Integral(bin120, bin200 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, ALL, bin120, bin200 - 1, true), 
         emuMasses.at(1).at(OS)->Integral(bin200, bin400 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, OS, bin200, bin400 - 1, true), 
         emuMasses.at(1).at(LS)->Integral(bin200, bin400 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, LS, bin200, bin400 - 1, true), 
         emuMasses.at(1).at(ALL)->Integral(bin200, bin400 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, ALL, bin200, bin400 - 1, true),
         emuMasses.at(1).at(OS)->Integral(bin400, nBins + 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, OS, bin400, nBins + 1, true), 
         emuMasses.at(1).at(LS)->Integral(bin400, nBins + 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, LS, bin400, nBins + 1, true), 
         emuMasses.at(1).at(ALL)->Integral(bin400, nBins + 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, allSamples, 1, ALL, bin400, nBins + 1, true)); 
  cout << "\\hline\\hline" << endl;
  printf("\\ttbar +  \\ttbar-like & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         emuMasses.at(TTBAR).at(OS)->Integral(bin120, bin200 - 1) - emuMasses.at(WJET).at(OS)->Integral(bin120, bin200 - 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, OS, bin120, bin200 - 1, true),
         emuMasses.at(TTBAR).at(LS)->Integral(bin120, bin200 - 1) - emuMasses.at(WJET).at(LS)->Integral(bin120, bin200 - 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, LS, bin120, bin200 - 1, true),
         emuMasses.at(TTBAR).at(ALL)->Integral(bin120, bin200 - 1) - emuMasses.at(WJET).at(ALL)->Integral(bin120, bin200 - 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin120, bin200 - 1, true),
         emuMasses.at(TTBAR).at(OS)->Integral(bin200, bin400 - 1) - emuMasses.at(WJET).at(OS)->Integral(bin200, bin400 - 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, OS, bin200, bin400 - 1, true),
         emuMasses.at(TTBAR).at(LS)->Integral(bin200, bin400 - 1) - emuMasses.at(WJET).at(LS)->Integral(bin200, bin400 - 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, LS, bin200, bin400 - 1, true),
         emuMasses.at(TTBAR).at(ALL)->Integral(bin200, bin400 - 1) - emuMasses.at(WJET).at(ALL)->Integral(bin200, bin400 - 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin200, bin400 - 1, true),
         emuMasses.at(TTBAR).at(OS)->Integral(bin400, nBins + 1) - emuMasses.at(WJET).at(OS)->Integral(bin400, nBins + 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, OS, bin400, nBins + 1, true),
         emuMasses.at(TTBAR).at(LS)->Integral(bin400, nBins + 1) - emuMasses.at(WJET).at(LS)->Integral(bin400, nBins + 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, LS, bin400, nBins + 1, true),
         emuMasses.at(TTBAR).at(ALL)->Integral(bin400, nBins + 1) - emuMasses.at(WJET).at(ALL)->Integral(bin400, nBins + 1), CalcAllErr(emuMasses, systErrMCLuEff, ttLikeSamples, ALL, bin400, nBins + 1, true));

  printf("contaminations        & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         emuMasses.at(WJET).at(OS)->Integral(bin120, bin200 - 1) - emuMasses.at(QCD).at(OS)->Integral(bin120, bin200 - 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, OS, bin120, bin200 - 1, true),
         emuMasses.at(WJET).at(LS)->Integral(bin120, bin200 - 1) - emuMasses.at(QCD).at(LS)->Integral(bin120, bin200 - 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, LS, bin120, bin200 - 1, true),
         emuMasses.at(WJET).at(ALL)->Integral(bin120, bin200 - 1) - emuMasses.at(QCD).at(ALL)->Integral(bin120, bin200 - 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin120, bin200 - 1, true),
         emuMasses.at(WJET).at(OS)->Integral(bin200, bin400 - 1) - emuMasses.at(QCD).at(OS)->Integral(bin200, bin400 - 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, OS, bin200, bin400 - 1, true),
         emuMasses.at(WJET).at(LS)->Integral(bin200, bin400 - 1) - emuMasses.at(QCD).at(LS)->Integral(bin200, bin400 - 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, LS, bin200, bin400 - 1, true),
         emuMasses.at(WJET).at(ALL)->Integral(bin200, bin400 - 1) - emuMasses.at(QCD).at(ALL)->Integral(bin200, bin400 - 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin200, bin400 - 1, true),
         emuMasses.at(WJET).at(OS)->Integral(bin400, nBins + 1) - emuMasses.at(QCD).at(OS)->Integral(bin400, nBins + 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, OS, bin400, nBins + 1, true),
         emuMasses.at(WJET).at(LS)->Integral(bin400, nBins + 1) - emuMasses.at(QCD).at(LS)->Integral(bin400, nBins + 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, LS, bin400, nBins + 1, true),
         emuMasses.at(WJET).at(ALL)->Integral(bin400, nBins + 1) - emuMasses.at(QCD).at(ALL)->Integral(bin400, nBins + 1), CalcAllErr(emuMasses, systErrMCLuEff, contamSamples, ALL, bin400, nBins + 1, true));

  printf("multi-jet             & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f & %.0f $\\pm$ %.0f \\\\\n", 
         emuMasses.at(QCD).at(OS)->Integral(bin120, bin200 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, OS, bin120, bin200 - 1, true),
         emuMasses.at(QCD).at(LS)->Integral(bin120, bin200 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, LS, bin120, bin200 - 1, true),
         emuMasses.at(QCD).at(ALL)->Integral(bin120, bin200 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, ALL, bin120, bin200 - 1, true),
         emuMasses.at(QCD).at(OS)->Integral(bin200, bin400 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, OS, bin200, bin400 - 1, true),
         emuMasses.at(QCD).at(LS)->Integral(bin200, bin400 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, LS, bin200, bin400 - 1, true),
         emuMasses.at(QCD).at(ALL)->Integral(bin200, bin400 - 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, ALL, bin200, bin400 - 1, true),
         emuMasses.at(QCD).at(OS)->Integral(bin400, nBins + 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, OS, bin400, nBins + 1, true),
         emuMasses.at(QCD).at(LS)->Integral(bin400, nBins + 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, LS, bin400, nBins + 1, true),
         emuMasses.at(QCD).at(ALL)->Integral(bin400, nBins + 1), CalcAllErrWithQCD(emuMasses, systErrMCLuEff, onlyQCD, QCD, ALL, bin400, nBins + 1, true));

  cout << "\\hline\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\caption{Number of $e\\mu$ events with invariant mass in different regions and with different charge combinations.}" << endl;
  cout << "\\label{tab:emu_event_yield}" << endl;
  cout << "\\end{table}" << endl;
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

float CalcAllErr(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked)
{
   float numEv;

   // for stacked histograms subtract all contributions below the current one 
   if (stacked) {
      for (unsigned int i = 0; i < samples.size() && i < histos.size()-1; ++i) {
         if (samples[i]) {
            if (i < histos.size()-2) numEv = histos.at(i+1).at(region)->Integral(lowerBin, upperBin) - histos.at(i+2).at(region)->Integral(lowerBin, upperBin);
            else numEv = histos.at(i+1).at(region)->Integral(lowerBin, upperBin);
         }
      }
   } else {
      for (unsigned int i = 0; i < samples.size(); ++i) {
         if (samples[i])
            numEv += histos.at(i+1).at(region)->Integral(lowerBin, upperBin);
      }
   }

   float systErr = CalcSystErr(histos, errors, samples, region, lowerBin, upperBin, stacked);

   //return sqrt(numEv + systErr * systErr);
   return sqrt(systErr * systErr);
}

float CalcSystErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int region, int lowerBin, int upperBin, bool stacked)
{
   float qcdErr = 0.;
   vector<bool> allButQCD(8, true);

   qcdErr = 2 * sqrt(histos.at(DATA).at(LS)->Integral(lowerBin, upperBin) + pow(CalcSystErr(histos, errors, allButQCD, LS, lowerBin, upperBin, stacked), 2));
   errors.back() = qcdErr / histos.at(QCD).at(region)->Integral(lowerBin, upperBin);

   return CalcSystErr(histos, errors, samples, region, lowerBin, upperBin, stacked);
}

float CalcAllErrWithQCD(vector<vector<TH1F *> > &histos, vector<float> &errors, vector<bool> &samples, int sample, int region, int lowerBin, int upperBin, bool stacked)
{
   float statErr = sqrt(histos.at(sample).at(region)->Integral(lowerBin, upperBin));
   float systErr = CalcSystErrWithQCD(histos, errors, samples, region, lowerBin, upperBin, stacked);

   //return sqrt(statErr * statErr + systErr * systErr);
   return sqrt(systErr * systErr);
}
