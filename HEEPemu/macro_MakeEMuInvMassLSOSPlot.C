{
  #include <string>
  #include <sstream>

  #include "TH1F.h"
  //#include "plotstyle.C"
 
  // parameters //////////////////////////////////////////////////////////////
  // luminosity of data 
  TFile input("testEmuSpec.root", "open");
  const float lumi = 3190;

  bool logPlot = false;
  bool integrate = false;
  ////////////////////////////////////////////////////////////////////////////

  input.cd();

  // INTEGRATING FROM THE RIGHT SIDE!
  if (integrate) { 
    // //LOOPING ON BINS
    double error;
    for (int i = 1; i < 101; ++i) {
      LS_ttbar->SetBinContent(i, LS_ttbar->IntegralAndError(i, 100, error));
      LS_ttbar->SetBinError(i, error);
      LS_ztautau->SetBinContent(i, LS_ztautau->IntegralAndError(i, 100, error));
      LS_ztautau->SetBinError(i, error);
      LS_ww->SetBinContent(i, LS_ww->IntegralAndError(i, 100, error));
      LS_ww->SetBinError(i, error);
      LS_wz->SetBinContent(i, LS_wz->IntegralAndError(i, 100, error));
      LS_wz->SetBinError(i, error);
      LS_tw->SetBinContent(i, LS_tw->IntegralAndError(i, 100, error));
      LS_tw->SetBinError(i, error);
      LS_wjet->SetBinContent(i, LS_wjet->IntegralAndError(i, 100, error));
      LS_wjet->SetBinError(i, error);
      LS_zmum->SetBinContent(i, LS_zmum->IntegralAndError(i, 100, error));
      LS_zmum->SetBinError(i, error);
      LS_zele->SetBinContent(i, LS_zele->IntegralAndError(i, 100, error));
      LS_zele->SetBinError(i, error);
      LS_data->SetBinContent(i, LS_data->IntegralAndError(i, 100, error));   
      LS_data->SetBinError(i, error);

      OS_ttbar->SetBinContent(i, OS_ttbar->IntegralAndError(i, 100, error));
      OS_ttbar->SetBinError(i, error);
      OS_ztautau->SetBinContent(i, OS_ztautau->IntegralAndError(i, 100, error));
      OS_ztautau->SetBinError(i, error);
      OS_ww->SetBinContent(i, OS_ww->IntegralAndError(i, 100, error));
      OS_ww->SetBinError(i, error);
      OS_wz->SetBinContent(i, OS_wz->IntegralAndError(i, 100, error));
      OS_wz->SetBinError(i, error);
      OS_tw->SetBinContent(i, OS_tw->IntegralAndError(i, 100, error));
      OS_tw->SetBinError(i, error);
      OS_wjet->SetBinContent(i, OS_wjet->IntegralAndError(i, 100, error));
      OS_wjet->SetBinError(i, error);
      OS_zmum->SetBinContent(i, OS_zmum->IntegralAndError(i, 100, error));
      OS_zmum->SetBinError(i, error);
      OS_zele->SetBinContent(i, OS_zele->IntegralAndError(i, 100, error));
      OS_zele->SetBinError(i, error);
      OS_data->SetBinContent(i, OS_data->IntegralAndError(i, 100, error));   
      OS_data->SetBinError(i, error);
    }     
  }
  
  //setPlotStyle();
 
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
  
  LS_ttbar->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
  LS_ttbar->GetXaxis()->SetTitleSize(0.04);
  LS_ttbar->GetXaxis()->SetLabelSize(0.04);

  if (integrate) LS_ttbar->GetYaxis()->SetTitle("# of LS events >= M_{e#mu}");
  else LS_ttbar->GetYaxis()->SetTitle("# of LS events / 10 GeV/c^{2}");
  LS_ttbar->GetYaxis()->SetTitleSize(0.04);
  LS_ttbar->GetYaxis()->SetTitleOffset(1.2);

  //
  LS_ttbar->SetFillColor(2);
  LS_ttbar->Draw("HIST");

  LS_ztautau->SetFillColor(6);
  LS_ztautau->Draw("HISTsames");

  LS_ww->SetFillColor(5);
  LS_ww->Draw("HISTsames");

  LS_wz->SetFillColor(8);
  LS_wz->Draw("HISTsames");

  LS_tw->SetFillColor(46);
  LS_tw->Draw("HISTsames");

  LS_wjet->SetFillColor(4);
  LS_wjet->Draw("HISTsames");

  LS_zmum->SetFillColor(7);
  LS_zmum->Draw("HISTsames");

  LS_zele->SetFillColor(3);
  LS_zele->Draw("HISTsames");

  //
  LS_data->SetLineWidth(2);
  LS_data->SetMarkerStyle(8);
  LS_data->SetMarkerSize(0.8);
  LS_data->Draw("sames");

  double dataMax = LS_data->GetMaximum();
  double ttbarMax = LS_ttbar->GetMaximum();
  if (dataMax > ttbarMax) LS_ttbar->SetMaximum(dataMax * 1.1);
  //if (LS_data->GetMaximum() > LS_ttbar->GetMaximum()) LS_ttbar->GetYaxis()->SetRangeUser(0., LS_data->GetMaximum() * 1.1);

  //Redraw axis
  LS_ttbar->Draw("sameaxis");

  //
  TLegend legendLS(0.7, 0.8, 0.88, 0.54);
  legendLS.SetTextSize(0.03);
  legendLS.SetFillColor(0);

  legendLS.AddEntry(LS_data, "Data");
  legendLS.AddEntry(LS_ttbar, "t #bar{t} (MC)");
  legendLS.AddEntry(LS_ztautau, "Z #rightarrow #tau #tau (MC)");
  legendLS.AddEntry(LS_ww, "WW, (MC)");
  legendLS.AddEntry(LS_wz, "WZ, (MC)");
  legendLS.AddEntry(LS_tw, "tW, (MC)");
  legendLS.AddEntry(LS_wjet, "W+jet (MC)");
  legendLS.AddEntry(LS_zmum, "Z #rightarrow #mu #mu (MC)");
  legendLS.AddEntry(LS_zele, "Z #rightarrow ee (MC)");

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

  OS_ttbar->GetXaxis()->SetTitle("M_{e #mu} (GeV/c^{2})");
  OS_ttbar->GetXaxis()->SetTitleSize(0.04);
  OS_ttbar->GetXaxis()->SetLabelSize(0.04);

  if (integrate) OS_ttbar->GetYaxis()->SetTitle("# of OS events >= M_{e#mu}");
  else OS_ttbar->GetYaxis()->SetTitle("# OS events / 10 GeV/c^{2}");
  OS_ttbar->GetYaxis()->SetTitleSize(0.04);
  OS_ttbar->GetYaxis()->SetTitleOffset(1.1);
  
  OS_ttbar->SetFillColor(2);
  OS_ttbar->Draw("HIST");

  OS_ztautau->SetFillColor(6);
  OS_ztautau->Draw("HISTsames");

  OS_ww->SetFillColor(5);
  OS_ww->Draw("HISTsames");

  OS_wz->SetFillColor(8);
  OS_wz->Draw("HISTsames");

  OS_tw->SetFillColor(46);
  OS_tw->Draw("HISTsames");

  OS_wjet->SetFillColor(4);
  OS_wjet->Draw("HISTsames");

  OS_zmum->SetFillColor(7);
  OS_zmum->Draw("HISTsames");

  OS_zele->SetFillColor(3);
  OS_zele->Draw("HISTsames");

  //
  OS_data->SetLineWidth(2);
  OS_data->SetMarkerStyle(8);
  OS_data->SetMarkerSize(0.8);
  OS_data->Draw("sames");

  dataMax = OS_data->GetMaximum();
  ttbarMax = OS_ttbar->GetMaximum();
  if (dataMax > ttbarMax) OS_ttbar->SetMaximum(dataMax * 1.1);
  //if (OS_data->GetMaximum() > OS_ttbar->GetMaximum()) OS_ttbar->SetMaximum(OS_data->GetMaximum() * 0.1);

  //Redraw axis
  OS_ttbar->Draw("sameaxis");

  //
  TLegend legendOS(0.7, 0.8, 0.88, 0.54);
  legendOS.SetTextSize(0.03);
  legendOS.SetFillColor(0);

  legendOS.AddEntry(OS_data, "Data");
  legendOS.AddEntry(OS_ttbar, "t #bar{t} (MC)");
  legendOS.AddEntry(OS_ztautau, "Z #rightarrow #tau #tau (MC)");
  legendOS.AddEntry(OS_ww, "WW, (MC)");
  legendOS.AddEntry(OS_wz, "WZ, (MC)");
  legendOS.AddEntry(OS_tw, "tW, (MC)");
  legendOS.AddEntry(OS_wjet, "W+jet (MC)");
  legendOS.AddEntry(OS_zmum, "Z #rightarrow #mu #mu (MC)");
  legendOS.AddEntry(OS_zele, "Z #rightarrow ee (MC)");

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
