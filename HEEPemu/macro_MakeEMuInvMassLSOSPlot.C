{

  //#include "plotstyle.C"
  
  //TFile input("EMu_1092pb_QCD10pc_Trig94pc_PVXinfo_EE40_electrons_PURew28jun__11July.root","open");
  TFile input("testEmuSpec.root","open");

  input.cd();

  // INTEGRATING FROM THE RIGHT SIDE!

  // //LOOPING ON BINS
//   for (int i = 1; i<101; i++) {
//     LS_ttbar->SetBinContent(i, LS_ttbar->Integral(i,100) );
//     LS_ztautau->SetBinContent(i, LS_ztautau->Integral(i,100) );
//     LS_ww->SetBinContent(i, LS_ww->Integral(i,100) );
//     LS_wjet->SetBinContent(i, LS_wjet->Integral(i,100) );
//     LS_data->SetBinContent(i, LS_data->Integral(i,100) );   
//   }     
  
  //setPlotStyle();
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.3);

  // LS---------
  
   LS_ttbar->GetXaxis()->SetTitle("M_{e #mu} (GeV/c^{2})");
   LS_ttbar->GetXaxis()->SetTitleSize(0.047);
   LS_ttbar->GetXaxis()->SetTitleOffset(0.9);
   LS_ttbar->GetXaxis()->SetLabelSize(0.040);

   LS_ttbar->GetYaxis()->SetTitle("Nb LS events  / 10 GeV/c^{2}");
   LS_ttbar->GetYaxis()->SetTitleSize(0.047);
   LS_ttbar->GetYaxis()->SetTitleOffset(1.3);

   //
   LS_ttbar->SetFillColor(2);
   LS_ttbar->Draw("HIST");

   LS_ztautau->SetFillColor(6);
   LS_ztautau->Draw("HISTsames");

   LS_ww->SetFillColor(5);
   LS_ww->Draw("HISTsames");

   LS_wjet->SetFillColor(4);
   LS_wjet->Draw("HISTsames");

   //LS_zmum->SetFillColor(3);
   //LS_zmum->Draw("HISTsames");

   //LS_zele->SetFillColor(7);
   //LS_zele->Draw("HISTsames");

   //
   LS_data->SetLineWidth(2);
   LS_data->Draw("sames");

   //Redraw axis
   LS_ttbar->Draw("sameaxis");

   //
   TLegend legend0(0.48,0.60,0.68,0.88);
   legend0.SetTextSize(0.04);
   legend0.SetFillColor(0);

   legend0.AddEntry(LS_data, "Data");
   legend0.AddEntry(LS_ttbar, "t #bar{t} (MC)");
   legend0.AddEntry(LS_ztautau, "Z #rightarrow #tau #tau (MC)");
   legend0.AddEntry(LS_ww, "WW, (MC)");
   legend0.AddEntry(LS_wjet, "W+jet + Z #rightarrow #mu #mu (MC)");

  // OS-------
//
//  OS_ttbar->GetXaxis()->SetTitle("M_{e #mu} (GeV/c^{2})");
//  OS_ttbar->GetXaxis()->SetTitleSize(0.047);
//  OS_ttbar->GetXaxis()->SetTitleOffset(0.9);
//  OS_ttbar->GetXaxis()->SetLabelSize(0.040);
//
//  OS_ttbar->GetYaxis()->SetTitle("Nb OS events / 10 GeV/c^{2}");
//  OS_ttbar->GetYaxis()->SetTitleSize(0.047);
//  OS_ttbar->GetYaxis()->SetTitleOffset(1.3);
//  
//  //
//  OS_ttbar->SetFillColor(2);
//  OS_ttbar->Draw("HIST");
//
//  OS_ztautau->SetFillColor(6);
//  OS_ztautau->Draw("HISTsames");
//
//  OS_ww->SetFillColor(5);
//  OS_ww->Draw("HISTsames");
//
//  OS_wjet->SetFillColor(4);
//  OS_wjet->Draw("HISTsames");
//
////   OS_zmum->SetFillColor(3);
////   OS_zmum->Draw("HISTsames");
//
////   OS_zele->SetFillColor(7);
////   OS_zele->Draw("HISTsames");
//
//  //
//  OS_data->SetLineWidth(2);
//  OS_data->Draw("sames");
//
//  //Redraw axis
//  OS_ttbar->Draw("sameaxis");
//
//  //
//  TLegend legend0(0.48,0.60,0.68,0.88);
//  legend0.SetTextSize(0.04);
//  legend0.SetFillColor(0);
//
//  legend0.AddEntry(OS_data, "Data");
//  legend0.AddEntry(OS_ttbar, "t #bar{t} (MC)");
//  legend0.AddEntry(OS_ztautau, "Z #rightarrow #tau #tau (MC)");
//  legend0.AddEntry(OS_ww, "WW, (MC)");
//  legend0.AddEntry(OS_wjet, "W+jet + Z #rightarrow #mu #mu (MC)");
  



  //---------------------------------------------------------
  legend0.SetBorderSize(0);
  legend0.Draw("sames");
  
  TPaveLabel label1(0.300,0.89,0.740,0.99,"CMS Preliminary","brNDC");
  label1.SetFillColor(0);
  label1.SetFillStyle(0);
  label1.SetBorderSize(0);
  label1.SetTextSize(0.40);
  label1.Draw("sames");

  TPaveLabel label0(0.300,0.89,0.740,0.99,"7 TeV, #int L dt = 3190pb^{-1}","brNDC");
  label0.SetFillColor(0);
  label0.SetFillStyle(0);
  label0.SetBorderSize(0);
  label0.SetTextSize(0.40);
  label0.Draw("sames");


  
  //gPad->Print("Plot_EMuMethod_R.eps");
  
}
