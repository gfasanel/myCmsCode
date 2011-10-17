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

  TFile input("testEmuSpec.root","open");
  //TFile input("testEmuSpec_PUcut20_3190pb-1.root","open");
  
  input.cd();

  // INTEGRATING FROM THE RIGHT SIDE!
   if (integrate) { 
    // //LOOPING ON BINS
    double error;
    for (int i = 1; i < 101; ++i) {
      //emu_ttbar->GetNbinsX();
      emu_ttbar->SetBinContent(i, emu_ttbar->IntegralAndError(i, 100, error));
      emu_ttbar->SetBinError(i, error);
      emu_ztautau->SetBinContent(i, emu_ztautau->IntegralAndError(i, 100, error));
      emu_ztautau->SetBinError(i, error);
      emu_ww->SetBinContent(i, emu_ww->IntegralAndError(i, 100, error));
      emu_ww->SetBinError(i, error);
      emu_wz->SetBinContent(i, emu_wz->IntegralAndError(i, 100, error));
      emu_wz->SetBinError(i, error);
      emu_tw->SetBinContent(i, emu_tw->IntegralAndError(i, 100, error));
      emu_tw->SetBinError(i, error);
      emu_wjet->SetBinContent(i, emu_wjet->IntegralAndError(i, 100, error));
      emu_wjet->SetBinError(i, error);
      emu_zmumu->SetBinContent(i, emu_zmumu->IntegralAndError(i, 100, error));
      emu_zmumu->SetBinError(i, error);
      emu_zee->SetBinContent(i, emu_zee->IntegralAndError(i, 100, error));
      emu_zee->SetBinError(i, error);
      emu_data->SetBinContent(i, emu_data->IntegralAndError(i, 100, error));   
      emu_data->SetBinError(i, error);
    }     
  }
 
  //setPlotStyle();

  TCanvas *cEmu = new TCanvas("cEmu", "emu Spectrum", 100, 100, 800, 600);
  cEmu->SetBorderMode(0);
  cEmu->SetFrameBorderMode(0);
  cEmu->SetFillColor(0);
  cEmu->SetFrameFillColor(0);
  if (logPlot) cEmu->SetLogy();
  cEmu->cd();
 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleXOffset(1.);
  gStyle->SetTitleYOffset(1.3);
  
  emu_ttbar->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
  emu_ttbar->GetXaxis()->SetTitleSize(0.04);
  emu_ttbar->GetXaxis()->SetLabelSize(0.04);

  if (integrate) emu_ttbar->GetYaxis()->SetTitle("# of events >= M_{e#mu}");
  else emu_ttbar->GetYaxis()->SetTitle("# of events / 10 GeV/c^{2}");
  emu_ttbar->GetYaxis()->SetTitleSize(0.04);
  emu_ttbar->GetYaxis()->SetTitleOffset(1.2);
  
  // //
//   emu_dilepton->SetFillColor(0);
//   emu_dilepton->Draw("HIST");
    
//   emu_ewk->SetFillColor(2);
//   emu_ewk->Draw("HISTsames");
    
//   emu_jet->SetFillColor(3);
//   emu_jet->Draw("HISTsames");

  emu_ttbar->SetFillColor(2);
  emu_ttbar->Draw("HIST");

  emu_ztautau->SetFillColor(6);
  emu_ztautau->Draw("HISTsames");

  emu_ww->SetFillColor(5);
  emu_ww->Draw("HISTsames");

  emu_wz->SetFillColor(8);
  emu_wz->Draw("HISTsames");

  emu_tw->SetFillColor(46);
  emu_tw->Draw("HISTsames");
  
  emu_wjet->SetFillColor(4);
  emu_wjet->Draw("HISTsames");

  emu_zmumu->SetFillColor(7);
  emu_zmumu->Draw("HISTsames");

  emu_zee->SetFillColor(3);
  emu_zee->Draw("HISTsames");

//   emu_qcd->SetFillColor(1);
//   emu_qcd->Draw("HISTsames");


  emu_data->SetLineWidth(2);
  emu_data->SetMarkerStyle(8);
  emu_data->SetMarkerSize(0.8);
  emu_data->Draw("sames");

  double dataMax = emu_data->GetMaximum();
  double ttbarMax = emu_ttbar->GetMaximum();
  if (dataMax > ttbarMax) emu_ttbar->SetMaximum(dataMax * 1.1);
  //if (emu_data->GetMaximum() > emu_ttbar->GetMaximum()) emu_ttbar->GetYaxis()->SetRangeUser(0., emu_data->GetMaximum() * 1.1);

  //Redraw axis
  emu_ttbar->Draw("sameaxis");

//   emu_LS->SetLineWidth(2);
//   emu_LS->SetLineColor(4);
//   emu_LS->SetMarkerStyle(23);
//   emu_LS->SetMarkerColor(4);
//   emu_LS->Draw("sames");
  
  //
  TLegend legend0(0.7, 0.8, 0.88, 0.54);
  legend0.SetTextSize(0.03);
  legend0.SetFillColor(0);
  
  legend0.AddEntry(emu_data, "Data");
//   legend0.AddEntry(emu_dilepton, "dilepton (MC)");
//   legend0.AddEntry(emu_ewk, "EWK (MC)");
//   legend0.AddEntry(emu_jet, "Jet BG (MC)");

   legend0.AddEntry(emu_ttbar, "t #bar{t} (MC)");
   legend0.AddEntry(emu_ztautau, "Z #rightarrow #tau #tau (MC)");
   legend0.AddEntry(emu_ww, "WW (MC)");
   legend0.AddEntry(emu_wz, "WZ (MC)");
   legend0.AddEntry(emu_tw, "tW (MC)");
   legend0.AddEntry(emu_wjet, "W+jets (MC)");
   legend0.AddEntry(emu_zmumu, "Z #rightarrow #mu #mu (MC)");
   legend0.AddEntry(emu_zee, "Z #rightarrow ee (MC)");
//   legend0.AddEntry(emu_qcd, "Jet BG (MC)");


  //legend0.AddEntry(emu_LS, "2 * Like Sign emu");
  
  legend0.SetBorderSize(0);
  legend0.Draw("sames");
  
  TPaveLabel label1(0.7, 0.91, 0.9, 0.82, "CMS Preliminary","brNDC");
  label1.SetFillColor(0);
  label1.SetFillStyle(0);
  label1.SetBorderSize(0);
  label1.SetTextSize(0.40);
  label1.Draw("sames");

  stringstream sStream;
  sStream.str("");
  sStream << "#sqrt{s} = 7TeV,  #int L dt = " << lumi << "pb^{-1}";
  TPaveLabel label0(0.3, 0.88, 0.65, 0.78, sStream.str().c_str(), "brNDC");
  label0.SetFillColor(0);
  label0.SetFillStyle(0);
  label0.SetBorderSize(0);
  label0.SetTextSize(0.40);
  label0.Draw("sames");


  
  //gPad->Print("Plot_EMuMethod_R.eps");
  
}
