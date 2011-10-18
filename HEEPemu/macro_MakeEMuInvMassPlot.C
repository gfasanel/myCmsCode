{
  #include <string>
  #include <sstream>
  #include "TH1F.h"
      
  // parameters //////////////////////////////////////////////////////////////
  TFile input("testEmuSpec3534pb-1.root", "open");
  //TFile input("testEmuSpec_PUcut20_3190pb-1.root","open");
  
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
      //emuMass_ttbar->GetNbinsX();
      emuMass_ttbar->SetBinContent(i, emuMass_ttbar->IntegralAndError(i, 100, error));
      emuMass_ttbar->SetBinError(i, error);
      emuMass_ztautau->SetBinContent(i, emuMass_ztautau->IntegralAndError(i, 100, error));
      emuMass_ztautau->SetBinError(i, error);
      emuMass_ww->SetBinContent(i, emuMass_ww->IntegralAndError(i, 100, error));
      emuMass_ww->SetBinError(i, error);
      emuMass_wz->SetBinContent(i, emuMass_wz->IntegralAndError(i, 100, error));
      emuMass_wz->SetBinError(i, error);
      emuMass_tw->SetBinContent(i, emuMass_tw->IntegralAndError(i, 100, error));
      emuMass_tw->SetBinError(i, error);
      emuMass_wjets->SetBinContent(i, emuMass_wjets->IntegralAndError(i, 100, error));
      emuMass_wjets->SetBinError(i, error);
      emuMass_zmumu->SetBinContent(i, emuMass_zmumu->IntegralAndError(i, 100, error));
      emuMass_zmumu->SetBinError(i, error);
      emuMass_zee->SetBinContent(i, emuMass_zee->IntegralAndError(i, 100, error));
      emuMass_zee->SetBinError(i, error);
      //emuMass_zz->SetBinContent(i, emuMass_zz->IntegralAndError(i, 100, error));
      //emuMass_zz->SetBinError(i, error);
      emuMass_data->SetBinContent(i, emuMass_data->IntegralAndError(i, 100, error));   
      emuMass_data->SetBinError(i, error);
    }     
  }
 
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
  
  emuMass_ttbar->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
  emuMass_ttbar->GetXaxis()->SetTitleSize(0.04);
  emuMass_ttbar->GetXaxis()->SetLabelSize(0.04);

  if (integrate) emuMass_ttbar->GetYaxis()->SetTitle("# of events >= M_{e#mu}");
  else emuMass_ttbar->GetYaxis()->SetTitle("# of events / 10 GeV/c^{2}");
  emuMass_ttbar->GetYaxis()->SetTitleSize(0.04);
  emuMass_ttbar->GetYaxis()->SetTitleOffset(1.2);
  
  // //
//   emu_dilepton->SetFillColor(0);
//   emu_dilepton->Draw("HIST");
    
//   emu_ewk->SetFillColor(2);
//   emu_ewk->Draw("HISTsames");
    
//   emu_jet->SetFillColor(3);
//   emu_jet->Draw("HISTsames");

  emuMass_ttbar->SetFillColor(2);
  emuMass_ttbar->Draw("HIST");

  emuMass_ztautau->SetFillColor(6);
  emuMass_ztautau->Draw("HISTsames");

  emuMass_ww->SetFillColor(5);
  emuMass_ww->Draw("HISTsames");

  emuMass_wz->SetFillColor(8);
  emuMass_wz->Draw("HISTsames");

  emuMass_tw->SetFillColor(46);
  emuMass_tw->Draw("HISTsames");
  
  emuMass_wjets->SetFillColor(4);
  emuMass_wjets->Draw("HISTsames");

  emuMass_zmumu->SetFillColor(7);
  emuMass_zmumu->Draw("HISTsames");

  emuMass_zee->SetFillColor(3);
  emuMass_zee->Draw("HISTsames");

  //emuMass_zz->SetFillColor(9);
  //emuMass_zz->Draw("HISTsames");

  emuMass_data->SetLineWidth(2);
  emuMass_data->SetMarkerStyle(8);
  emuMass_data->SetMarkerSize(0.8);
  emuMass_data->Draw("sames");

  double dataMax = emuMass_data->GetMaximum();
  double ttbarMax = emuMass_ttbar->GetMaximum();
  if (dataMax > ttbarMax) emuMass_ttbar->SetMaximum(dataMax * 1.1);
  //if (emuMass_data->GetMaximum() > emuMass_ttbar->GetMaximum()) emuMass_ttbar->SetMaximum(emuMass_data->GetMaximum() * 1.1);

  //Redraw axis
  emuMass_ttbar->Draw("sameaxis");

  //
  TLegend legend0(0.7, 0.8, 0.88, 0.54);
  legend0.SetTextSize(0.03);
  legend0.SetFillColor(0);
  
  legend0.AddEntry(emuMass_data, "Data");
//   legend0.AddEntry(emu_dilepton, "dilepton (MC)");
//   legend0.AddEntry(emu_ewk, "EWK (MC)");
//   legend0.AddEntry(emu_jet, "Jet BG (MC)");

  legend0.AddEntry(emuMass_ttbar, "t #bar{t} (MC)");
  legend0.AddEntry(emuMass_ztautau, "Z #rightarrow #tau #tau (MC)");
  legend0.AddEntry(emuMass_ww, "WW (MC)");
  legend0.AddEntry(emuMass_wz, "WZ (MC)");
  legend0.AddEntry(emuMass_tw, "tW (MC)");
  legend0.AddEntry(emuMass_wjets, "W+jets (MC)");
  legend0.AddEntry(emuMass_zmumu, "Z #rightarrow #mu #mu (MC)");
  legend0.AddEntry(emuMass_zee, "Z #rightarrow ee (MC)");
  //legend0.AddEntry(emuMass_zz, "ZZ (MC)");
  
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
  
}
