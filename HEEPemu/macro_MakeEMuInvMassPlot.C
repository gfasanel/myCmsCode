{

  //#include "plotstyle.C"

  //TFile input("EMu_1092pb_QCD10pc_Trig94pc_PVXinfo_EE40_electrons_PURew28jun__11July.root","open");
  TFile input("testEmuSpec.root","open");
  
  input.cd();

  // INTEGRATING FROM THE RIGHT SIDE!

 //  //LOOPING ON BINS
//   for (int i = 1; i<101; i++) {
//     emu_ttbar->SetBinContent(i, emu_ttbar->Integral(i,100) );
//     emu_ztautau->SetBinContent(i, emu_ztautau->Integral(i,100) );
//     emu_ww->SetBinContent(i, emu_ww->Integral(i,100) );
//     emu_wjet->SetBinContent(i, emu_wjet->Integral(i,100) );
//     emu_data->SetBinContent(i, emu_data->Integral(i,100) );   
//     //emu_data->SetBinError(i, sqrt(emu_data->GetBinContent(i) ) );   
//   }     
  
  //setPlotStyle();
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.3);
  
  emu_ttbar->GetXaxis()->SetTitle("M_{e #mu} (GeV/c^{2})");
  emu_ttbar->GetXaxis()->SetTitleSize(0.047);
  emu_ttbar->GetXaxis()->SetTitleOffset(0.9);
  emu_ttbar->GetXaxis()->SetLabelSize(0.040);

  emu_ttbar->GetYaxis()->SetTitle("Nb events / 10 GeV/c^{2}");
  emu_ttbar->GetYaxis()->SetTitleSize(0.047);
  emu_ttbar->GetYaxis()->SetTitleOffset(1.3);
  
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

//   emu_tw->SetFillColor(46);
//   emu_tw->Draw("HISTsames");
  
  emu_wjet->SetFillColor(4);
  emu_wjet->Draw("HISTsames");

//   emu_zee->SetFillColor(3);
//   emu_zee->Draw("HISTsames");

//   emu_zmumu->SetFillColor(7);
//   emu_zmumu->Draw("HISTsames");

//   emu_qcd->SetFillColor(1);
//   emu_qcd->Draw("HISTsames");


  emu_data->SetLineWidth(2);
  emu_data->SetMarkerStyle(8);
  emu_data->Draw("sames");

  //Redraw axis
  emu_ttbar->Draw("sameaxis");

//   emu_LS->SetLineWidth(2);
//   emu_LS->SetLineColor(4);
//   emu_LS->SetMarkerStyle(23);
//   emu_LS->SetMarkerColor(4);
//   emu_LS->Draw("sames");
  
  //
  TLegend legend0(0.48,0.60,0.68,0.88);
  legend0.SetTextSize(0.04);
  legend0.SetFillColor(0);
  
  legend0.AddEntry(emu_data, "Data");
//   legend0.AddEntry(emu_dilepton, "dilepton (MC)");
//   legend0.AddEntry(emu_ewk, "EWK (MC)");
//   legend0.AddEntry(emu_jet, "Jet BG (MC)");

  legend0.AddEntry(emu_ttbar, "t #bar{t} (MC)");
  legend0.AddEntry(emu_ztautau, "Z #rightarrow #tau #tau (MC)");
  //legend0.AddEntry(emu_ww, "WW, tW (MC)");
  legend0.AddEntry(emu_ww, "WW (MC)");
  legend0.AddEntry(emu_wjet, "W+jet + Z #rightarrow #mu #mu (MC)");

  // legend0.AddEntry(emu_ttbar, "t #bar{t} (MC)");
//   legend0.AddEntry(emu_ztautau, "Z #rightarrow #tau #tau (MC)");
//   legend0.AddEntry(emu_ww, "WW (MC)");
//   legend0.AddEntry(emu_tw, "tW (MC)");
//   legend0.AddEntry(emu_wjet, "W+jets (MC)");
//   legend0.AddEntry(emu_zee, "Z #rightarrow ee (MC)");
//   legend0.AddEntry(emu_zmumu, "Z #rightarrow #mu #mu (MC)");
//   legend0.AddEntry(emu_qcd, "Jet BG (MC)");


  //legend0.AddEntry(emu_LS, "2 * Like Sign emu");
  
  legend0.SetBorderSize(0);
  legend0.Draw("sames");
  
  TPaveLabel label0(0.300,0.89,0.740,0.99,"7 TeV, #int L dt = 3190pb^{-1}","brNDC");
  label0.SetFillColor(0);
  label0.SetFillStyle(0);
  label0.SetBorderSize(0);
  label0.SetTextSize(0.40);
  label0.Draw("sames");

  TPaveLabel label1(0.300,0.89,0.740,0.99,"CMS Preliminary","brNDC");
  label1.SetFillColor(0);
  label1.SetFillStyle(0);
  label1.SetBorderSize(0);
  label1.SetTextSize(0.40);
  label1.Draw("sames");


  
  //gPad->Print("Plot_EMuMethod_R.eps");
  
}
