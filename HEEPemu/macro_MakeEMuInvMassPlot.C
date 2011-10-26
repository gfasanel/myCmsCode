{
  #include <string>
  #include <sstream>
  #include <vector>
  #include "TH1F.h"
  #include "TCanvas.h"
 
  // parameters //////////////////////////////////////////////////////////////
  TFile input("testEmuSpec3534pb-1.root", "open");
  //TFile input("testEmuSpec_withoutTriggerBit_3534pb-1.root", "open");
  //TFile input("./PUreweightTests20111020/testEmuSpec_PUrew20_3534pb-1.root", "open");

  const float lumi = 3534;

  bool plotSign[3];
  plotSign[0] = true;  // all
  plotSign[1] = true;  // LS like sign
  plotSign[2] = true;  // OS opposite sign

  bool plotType[2];
  plotType[0] = true;  // emu spectrum
  plotType[1] = true;  // cumulative emu spectrum

  bool logPlot = false;
  ////////////////////////////////////////////////////////////////////////////

  // to keep the histogram when the file is closed
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2(kTRUE);

  TString histoSign[3] = {"", "LS_", "OS_"};
  TString histoTitleSign[3] = {"", " LS", " OS"};
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
    if (!plotSign[k]) continue;
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

      // integrate from the right side
      if (j == 1) { 
        // loop over bins
        double error;
        for (int i = 1; i < 101; ++i) {
          emuMass_data.back()->SetBinContent(i, emuMass_data.back()->IntegralAndError(i, 100, error));
          emuMass_data.back()->SetBinError(i, error);
          emuMass_ttbar.back()->SetBinContent(i, emuMass_ttbar.back()->IntegralAndError(i, 100, error));
          emuMass_ttbar.back()->SetBinError(i, error);
          emuMass_ztautau.back()->SetBinContent(i, emuMass_ztautau.back()->IntegralAndError(i, 100, error));
          emuMass_ztautau.back()->SetBinError(i, error);
          emuMass_ww.back()->SetBinContent(i, emuMass_ww.back()->IntegralAndError(i, 100, error));
          emuMass_ww.back()->SetBinError(i, error);
          emuMass_wz.back()->SetBinContent(i, emuMass_wz.back()->IntegralAndError(i, 100, error));
          emuMass_wz.back()->SetBinError(i, error);
          emuMass_tw.back()->SetBinContent(i, emuMass_tw.back()->IntegralAndError(i, 100, error));
          emuMass_tw.back()->SetBinError(i, error);
          emuMass_wjets.back()->SetBinContent(i, emuMass_wjets.back()->IntegralAndError(i, 100, error));
          emuMass_wjets.back()->SetBinError(i, error);
          emuMass_zmumu.back()->SetBinContent(i, emuMass_zmumu.back()->IntegralAndError(i, 100, error));
          emuMass_zmumu.back()->SetBinError(i, error);
          emuMass_zee.back()->SetBinContent(i, emuMass_zee.back()->IntegralAndError(i, 100, error));
          emuMass_zee.back()->SetBinError(i, error);
          //emuMass_zz.back()->SetBinContent(i, emuMass_zz.back()->IntegralAndError(i, 100, error));
          //emuMass_zz.back()->SetBinError(i, error);
          emuMass_qcd.back()->SetBinContent(i, emuMass_qcd.back()->IntegralAndError(i, 100, error));
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

      emuMass_ttbar.back()->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
      emuMass_ttbar.back()->GetXaxis()->SetTitleSize(0.04);
      emuMass_ttbar.back()->GetXaxis()->SetLabelSize(0.04);
      emuMass_ttbar.back()->GetXaxis()->SetRangeUser(60., 1000.);
    
      if (j == 1) emuMass_ttbar.back()->GetYaxis()->SetTitle("# of" + histoTitleSign[k] + " events >= M_{e#mu}");
      else emuMass_ttbar.back()->GetYaxis()->SetTitle("# of" + histoTitleSign[k] + " events / 10 GeV/c^{2}");
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
      TLegend legend(0.7, 0.8, 0.88, 0.54);
      legend.SetTextSize(0.03);
      legend.SetFillColor(0);
    
      legend.AddEntry(emuMass_data.back(), "Data");
      legend.AddEntry(emuMass_ttbar.back(), "t #bar{t} (MC)");
      legend.AddEntry(emuMass_ztautau.back(), "Z #rightarrow #tau #tau (MC)");
      legend.AddEntry(emuMass_ww.back(), "WW, (MC)");
      legend.AddEntry(emuMass_wz.back(), "WZ, (MC)");
      legend.AddEntry(emuMass_tw.back(), "tW, (MC)");
      legend.AddEntry(emuMass_wjets.back(), "W+jets (MC)");
      legend.AddEntry(emuMass_zmumu.back(), "Z #rightarrow #mu #mu (MC)");
      legend.AddEntry(emuMass_zee.back(), "Z #rightarrow ee (MC)");
      //legend.AddEntry(emuMass_zz.back(), "ZZ, (MC)");
      legend.AddEntry(emuMass_qcd.back(), "QCD");
    
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
  } // end loop over full, LS and OS
}
