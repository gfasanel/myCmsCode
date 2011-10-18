{
  #include <string>
  #include <sstream>
  #include "TH1F.h"
 
  // parameters //////////////////////////////////////////////////////////////
  TFile input("testEmuSpec3534pb-1.root", "open");
  const float lumi = 3534;

  bool logPlot = false;
  bool integrate = false;
  ////////////////////////////////////////////////////////////////////////////

  input.cd();

  // INTEGRATING FROM THE RIGHT SIDE!
  if (integrate) { 
    // //LOOPING ON BINS
    double error;
    for (int i = 1; i < 101; ++i) {
      emuMass_LS_ttbar->SetBinContent(i, emuMass_LS_ttbar->IntegralAndError(i, 100, error));
      emuMass_LS_ttbar->SetBinError(i, error);
      emuMass_LS_ztautau->SetBinContent(i, emuMass_LS_ztautau->IntegralAndError(i, 100, error));
      emuMass_LS_ztautau->SetBinError(i, error);
      emuMass_LS_ww->SetBinContent(i, emuMass_LS_ww->IntegralAndError(i, 100, error));
      emuMass_LS_ww->SetBinError(i, error);
      emuMass_LS_wz->SetBinContent(i, emuMass_LS_wz->IntegralAndError(i, 100, error));
      emuMass_LS_wz->SetBinError(i, error);
      emuMass_LS_tw->SetBinContent(i, emuMass_LS_tw->IntegralAndError(i, 100, error));
      emuMass_LS_tw->SetBinError(i, error);
      emuMass_LS_wjets->SetBinContent(i, emuMass_LS_wjets->IntegralAndError(i, 100, error));
      emuMass_LS_wjets->SetBinError(i, error);
      emuMass_LS_zmumu->SetBinContent(i, emuMass_LS_zmumu->IntegralAndError(i, 100, error));
      emuMass_LS_zmumu->SetBinError(i, error);
      emuMass_LS_zee->SetBinContent(i, emuMass_LS_zee->IntegralAndError(i, 100, error));
      emuMass_LS_zee->SetBinError(i, error);
      //emuMass_LS_zz->SetBinContent(i, emuMass_LS_zz->IntegralAndError(i, 100, error));
      //emuMass_LS_zz->SetBinError(i, error);
      emuMass_LS_data->SetBinContent(i, emuMass_LS_data->IntegralAndError(i, 100, error));   
      emuMass_LS_data->SetBinError(i, error);

      emuMass_OS_ttbar->SetBinContent(i, emuMass_OS_ttbar->IntegralAndError(i, 100, error));
      emuMass_OS_ttbar->SetBinError(i, error);
      emuMass_OS_ztautau->SetBinContent(i, emuMass_OS_ztautau->IntegralAndError(i, 100, error));
      emuMass_OS_ztautau->SetBinError(i, error);
      emuMass_OS_ww->SetBinContent(i, emuMass_OS_ww->IntegralAndError(i, 100, error));
      emuMass_OS_ww->SetBinError(i, error);
      emuMass_OS_wz->SetBinContent(i, emuMass_OS_wz->IntegralAndError(i, 100, error));
      emuMass_OS_wz->SetBinError(i, error);
      emuMass_OS_tw->SetBinContent(i, emuMass_OS_tw->IntegralAndError(i, 100, error));
      emuMass_OS_tw->SetBinError(i, error);
      emuMass_OS_wjets->SetBinContent(i, emuMass_OS_wjets->IntegralAndError(i, 100, error));
      emuMass_OS_wjets->SetBinError(i, error);
      emuMass_OS_zmumu->SetBinContent(i, emuMass_OS_zmumu->IntegralAndError(i, 100, error));
      emuMass_OS_zmumu->SetBinError(i, error);
      emuMass_OS_zee->SetBinContent(i, emuMass_OS_zee->IntegralAndError(i, 100, error));
      emuMass_OS_zee->SetBinError(i, error);
      //emuMass_OS_zz->SetBinContent(i, emuMass_OS_zz->IntegralAndError(i, 100, error));
      //emuMass_OS_zz->SetBinError(i, error);
      emuMass_OS_data->SetBinContent(i, emuMass_OS_data->IntegralAndError(i, 100, error));   
      emuMass_OS_data->SetBinError(i, error);
    }     
  }
  
  TCanvas *cLS = new TCanvas("cLS", "emu Spectrum LS", 100, 100, 800, 600);
  cLS->SetBorderMode(0);
  cLS->SetFrameBorderMode(0);
  cLS->SetFillColor(0);
  cLS->SetFrameFillColor(0);
  if (logPlot) cLS->SetLogy();
  TCanvas *cOS = new TCanvas("cOS", "emu Spectrum OS", 100, 100, 800, 600);
  cOS->SetBorderMode(0);
  cOS->SetFrameBorderMode(0);
  cOS->SetFillColor(0);
  cOS->SetFrameFillColor(0);
  if (logPlot) cOS->SetLogy();
 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleXOffset(1.);
  gStyle->SetTitleYOffset(1.3);

  ////////////////////////////////////////////////////////////////////////////
  // LS---------
  cLS->cd();
  
  emuMass_LS_ttbar->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
  emuMass_LS_ttbar->GetXaxis()->SetTitleSize(0.04);
  emuMass_LS_ttbar->GetXaxis()->SetLabelSize(0.04);

  if (integrate) emuMass_LS_ttbar->GetYaxis()->SetTitle("# of LS events >= M_{e#mu}");
  else emuMass_LS_ttbar->GetYaxis()->SetTitle("# of LS events / 10 GeV/c^{2}");
  emuMass_LS_ttbar->GetYaxis()->SetTitleSize(0.04);
  emuMass_LS_ttbar->GetYaxis()->SetTitleOffset(1.2);

  //
  emuMass_LS_ttbar->SetFillColor(2);
  emuMass_LS_ttbar->Draw("HIST");

  emuMass_LS_ztautau->SetFillColor(6);
  emuMass_LS_ztautau->Draw("HISTsames");

  emuMass_LS_ww->SetFillColor(5);
  emuMass_LS_ww->Draw("HISTsames");

  emuMass_LS_wz->SetFillColor(8);
  emuMass_LS_wz->Draw("HISTsames");

  emuMass_LS_tw->SetFillColor(46);
  emuMass_LS_tw->Draw("HISTsames");

  emuMass_LS_wjets->SetFillColor(4);
  emuMass_LS_wjets->Draw("HISTsames");

  emuMass_LS_zmumu->SetFillColor(7);
  emuMass_LS_zmumu->Draw("HISTsames");

  emuMass_LS_zee->SetFillColor(3);
  emuMass_LS_zee->Draw("HISTsames");

  //emuMass_LS_zz->SetFillColor(9);
  //emuMass_LS_zz->Draw("HISTsames");

  //
  emuMass_LS_data->SetLineWidth(2);
  emuMass_LS_data->SetMarkerStyle(8);
  emuMass_LS_data->SetMarkerSize(0.8);
  emuMass_LS_data->Draw("sames");

  double dataMax = emuMass_LS_data->GetMaximum();
  double ttbarMax = emuMass_LS_ttbar->GetMaximum();
  if (dataMax > ttbarMax) emuMass_LS_ttbar->SetMaximum(dataMax * 1.1);
  //if (emuMass_LS_data->GetMaximum() > emuMass_LS_ttbar->GetMaximum()) emuMass_LS_ttbar->SetMaximum(emuMass_LS_data->GetMaximum() * 1.1);

  //Redraw axis
  emuMass_LS_ttbar->Draw("sameaxis");

  //
  TLegend legendLS(0.7, 0.8, 0.88, 0.54);
  legendLS.SetTextSize(0.03);
  legendLS.SetFillColor(0);

  legendLS.AddEntry(emuMass_LS_data, "Data");
  legendLS.AddEntry(emuMass_LS_ttbar, "t #bar{t} (MC)");
  legendLS.AddEntry(emuMass_LS_ztautau, "Z #rightarrow #tau #tau (MC)");
  legendLS.AddEntry(emuMass_LS_ww, "WW, (MC)");
  legendLS.AddEntry(emuMass_LS_wz, "WZ, (MC)");
  legendLS.AddEntry(emuMass_LS_tw, "tW, (MC)");
  legendLS.AddEntry(emuMass_LS_wjets, "W+jets (MC)");
  legendLS.AddEntry(emuMass_LS_zmumu, "Z #rightarrow #mu #mu (MC)");
  legendLS.AddEntry(emuMass_LS_zee, "Z #rightarrow ee (MC)");
  //legendLS.AddEntry(emuMass_LS_zz, "ZZ, (MC)");

  legendLS.SetBorderSize(0);
  legendLS.Draw("sames");
  
  TPaveLabel labelPrelimLS(0.7, 0.91, 0.9, 0.82, "CMS Preliminary", "brNDC");
  labelPrelimLS.SetFillColor(0);
  labelPrelimLS.SetFillStyle(0);
  labelPrelimLS.SetBorderSize(0);
  labelPrelimLS.SetTextSize(0.40);
  labelPrelimLS.Draw("sames");

  stringstream sStream;
  sStream.str("");
  sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
  TPaveLabel labelLumiLS(0.3, 0.88, 0.65, 0.78, sStream.str().c_str(), "brNDC");
  labelLumiLS.SetFillColor(0);
  labelLumiLS.SetFillStyle(0);
  labelLumiLS.SetBorderSize(0);
  labelLumiLS.SetTextSize(0.40);
  labelLumiLS.Draw("sames");

  ////////////////////////////////////////////////////////////////////////////
  // OS-------
  cOS->cd();

  emuMass_OS_ttbar->GetXaxis()->SetTitle("M_{e #mu} (GeV/c^{2})");
  emuMass_OS_ttbar->GetXaxis()->SetTitleSize(0.04);
  emuMass_OS_ttbar->GetXaxis()->SetLabelSize(0.04);

  if (integrate) emuMass_OS_ttbar->GetYaxis()->SetTitle("# of OS events >= M_{e#mu}");
  else emuMass_OS_ttbar->GetYaxis()->SetTitle("# OS events / 10 GeV/c^{2}");
  emuMass_OS_ttbar->GetYaxis()->SetTitleSize(0.04);
  emuMass_OS_ttbar->GetYaxis()->SetTitleOffset(1.1);
  
  emuMass_OS_ttbar->SetFillColor(2);
  emuMass_OS_ttbar->Draw("HIST");

  emuMass_OS_ztautau->SetFillColor(6);
  emuMass_OS_ztautau->Draw("HISTsames");

  emuMass_OS_ww->SetFillColor(5);
  emuMass_OS_ww->Draw("HISTsames");

  emuMass_OS_wz->SetFillColor(8);
  emuMass_OS_wz->Draw("HISTsames");

  emuMass_OS_tw->SetFillColor(46);
  emuMass_OS_tw->Draw("HISTsames");

  emuMass_OS_wjets->SetFillColor(4);
  emuMass_OS_wjets->Draw("HISTsames");

  emuMass_OS_zmumu->SetFillColor(7);
  emuMass_OS_zmumu->Draw("HISTsames");

  emuMass_OS_zee->SetFillColor(3);
  emuMass_OS_zee->Draw("HISTsames");

  //emuMass_OS_zz->SetFillColor(9);
  //emuMass_OS_zz->Draw("HISTsames");

  //
  emuMass_OS_data->SetLineWidth(2);
  emuMass_OS_data->SetMarkerStyle(8);
  emuMass_OS_data->SetMarkerSize(0.8);
  emuMass_OS_data->Draw("sames");

  dataMax = emuMass_OS_data->GetMaximum();
  ttbarMax = emuMass_OS_ttbar->GetMaximum();
  if (dataMax > ttbarMax) emuMass_OS_ttbar->SetMaximum(dataMax * 1.1);
  //if (emuMass_OS_data->GetMaximum() > emuMass_OS_ttbar->GetMaximum()) emuMass_OS_ttbar->SetMaximum(emuMass_OS_data->GetMaximum() * 0.1);

  //Redraw axis
  emuMass_OS_ttbar->Draw("sameaxis");

  //
  TLegend legendOS(0.7, 0.8, 0.88, 0.54);
  legendOS.SetTextSize(0.03);
  legendOS.SetFillColor(0);

  legendOS.AddEntry(emuMass_OS_data, "Data");
  legendOS.AddEntry(emuMass_OS_ttbar, "t #bar{t} (MC)");
  legendOS.AddEntry(emuMass_OS_ztautau, "Z #rightarrow #tau #tau (MC)");
  legendOS.AddEntry(emuMass_OS_ww, "WW, (MC)");
  legendOS.AddEntry(emuMass_OS_wz, "WZ, (MC)");
  legendOS.AddEntry(emuMass_OS_tw, "tW, (MC)");
  legendOS.AddEntry(emuMass_OS_wjets, "W+jets (MC)");
  legendOS.AddEntry(emuMass_OS_zmumu, "Z #rightarrow #mu #mu (MC)");
  legendOS.AddEntry(emuMass_OS_zee, "Z #rightarrow ee (MC)");
  //legendOS.AddEntry(emuMass_OS_zz, "ZZ, (MC)");

  legendOS.SetBorderSize(0);
  legendOS.Draw("sames");

  TPaveLabel labelPrelimOS(0.7, 0.91, 0.9, 0.82, "CMS Preliminary", "brNDC");
  labelPrelimOS.SetFillColor(0);
  labelPrelimOS.SetFillStyle(0);
  labelPrelimOS.SetBorderSize(0);
  labelPrelimOS.SetTextSize(0.40);
  labelPrelimOS.Draw("sames");

  TPaveLabel labelLumiOS(0.3, 0.88, 0.65, 0.78, sStream.str().c_str(), "brNDC");
  labelLumiOS.SetFillColor(0);
  labelLumiOS.SetFillStyle(0);
  labelLumiOS.SetBorderSize(0);
  labelLumiOS.SetTextSize(0.40);
  labelLumiOS.Draw("sames");

  //gPad->Print("Plot_EMuMethod_R.eps");
  
}
