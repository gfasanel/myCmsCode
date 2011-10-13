
#define fakeRate_cxx
#include "fakeRate.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <stdio.h>
#include "TImage.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TPaveLabel.h"
//#include "/beo5/charafbis/HomeSafe03062010/CMSSW_2_2_6/src/HEEPRoutine/plostyle.C"
#include "/user/lathomas/TandP/CMSSW_3_8_4/src/DataFormats/Math/interface/deltaR.h"

bool saveplots = false;
const float valmin = 0.8;

void FakeRate::Loop()
{

  vector<TH1F *> denomEB;
  vector<TH1F *> denomEE;
  vector<TH1F *> numEB;
  vector<TH1F *> numEE;
  
  const int vtxmax = 15;

  // plot variables
  const int nbin = 100;
  const int ptmax = 1000;  

  bool requireTagEoverP=true;
  bool requireTaginBarrel=true; 
  bool requiretagheep = true; 
  float MassMin =50;
  float MassMax = 10000; 

  const int nbinspt = 20; 
  const int binptmin = 100;
  const int binptmax =500; 
  const int ptcut = 20;      // minimum pt of required gsf electrons
  const int maxNumGSFEle = 1;   // number of required gsf electrons
  const float maxSigmaIEtaIEtaEB = 0.013;
  const float maxSigmaIEtaIEtaEE = 0.034;
  const float maxHOverEEB = 0.15;
  const float maxHOverEEE = 0.1;
  const float maxPFMet = 20;

  vector<TFile *> input;
  //input.push_back(new TFile("/user/lathomas/TandP/CMSSW_3_8_4/MCFiles/MC2011highpileuppythia.root","open"));
  input.push_back(new TFile("/user/lathomas/TandP/CMSSW_3_8_4/DataFiles/data2011_101ppb.root","open"));

  const int nbFile = input.size();

   //input.push_back(new TFile("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/QCD_Pt_300to470_TuneZ2_7TeV_pythia6/QCD_Pt_300to470_TuneZ2_7TeV_pythia6-Spring11-PU_S1_START311_V1G1-v1-AODSIM_Skim2Ele/319d9d50ddc1c21c2a4623a85e06b6f6/output_86_1_QK8.root","open")); 
  for (int p = 0; p<nbFile; p++) {

    TString wholename((input[p])->GetName());
    cout << wholename << endl;
    int index1 = wholename.Index("FullTrees/")+10;
    int index2 = wholename.Index(".root");
    cout << index1 << " " << index2 << endl;
    cout << wholename(index1,index2-index1) << endl;
    TString name(wholename(index1,index2-index1));
    cout << name << endl;
    cout << name.Data() << endl;

    // Denominator for fake rate 
    denomEB.push_back(new TH1F(TString("denomEB").Append(name.Data()),"",nbin,0,ptmax)); 
    denomEE.push_back(new TH1F(TString("denomEE").Append(name.Data()),"",nbin,0,ptmax)); 
    // Numerator for fake rate
    numEB.push_back(new TH1F(TString("numEB").Append(name.Data()),"",nbin,0,ptmax)); 
    numEE.push_back(new TH1F(TString("numEE").Append(name.Data()),"",nbin,0,ptmax)); 

    cout << "file " << (input[p])->GetName() << endl;
    (input[p])->cd();
    TTree *thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");
    Init(thetree);
    Long64_t nentries = (*thetree).GetEntries();
    cout << nentries << " entries" << endl;
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      thetree->GetEntry(jentry);
      if ((jentry%50000) == 0 ) cout << "entry " << jentry << endl;

      // find out how many good gsf electrons with high enough et are in the events
      int numOfGoodGSFEle = 0;
      for (int m = 0; m < gsf_size; m++) {
        if (gsf_gsfet[m] > ptcut) numOfGoodGSFEle++;
      }
      // only events with a maximum number of gsf electron with high enough et should be looked at
      if (numOfGoodGSFEle > maxNumGSFEle) continue;

      // we have now only events with the required number of gsf electrons
      // loop over them
      for (int n = 0; n < gsf_size; n++) {
	if (gsf_gsfet[n] < ptcut) continue;

        if (gsf_isEB[n]) {
          if (gsf_sigmaIetaIeta[n] >=  maxSigmaIEtaIEtaEB) continue;
          if (gsf_hovere[n] >= maxHOverEEB ) continue;
          if (calomet >= maxPFMet ) continue;
          denomEB[p]->Fill(gsf_gsfet[n]);
          // the numerator has also to pass the HEEP selection
          if (IsHEEP(gsf_gsfet[n],
                     gsfsc_eta[n],
                     gsf_deltaeta[n],
                     gsf_deltaphi[n],
                     gsf_trackiso[n],
                     gsf_ecaliso[n],
                     gsf_hcaliso1[n],
                     gsf_hcaliso2[n],
                     gsf_sigmaIetaIeta[n],
                     gsf_hovere[n],
                     gsf_e2x5overe5x5[n],
                     gsf_e1x5overe5x5[n])) numEB[p]->Fill(gsf_gsfet[n]);
        }
        else if (gsf_isEE[n]) {
          if (gsf_sigmaIetaIeta[n] >=  maxSigmaIEtaIEtaEE) continue;
          if (gsf_hovere[n] >= maxHOverEEE ) continue;
          if (calomet >= maxPFMet ) continue;
          denomEE[p]->Fill(gsf_gsfet[n]);
          // the numerator has also to pass the HEEP selection
          if (IsHEEP(gsf_gsfet[n],
                     gsfsc_eta[n],
                     gsf_deltaeta[n],
                     gsf_deltaphi[n],
                     gsf_trackiso[n],
                     gsf_ecaliso[n],
                     gsf_hcaliso1[n],
                     gsf_hcaliso2[n],
                     gsf_sigmaIetaIeta[n],
                     gsf_hovere[n],
                     gsf_e2x5overe5x5[n],
                     gsf_e1x5overe5x5[n])) numEE[p]->Fill(gsf_gsfet[n]);
        }

      }  
  
    }
 
  }


  SingleEff(numEB, denomEB, numEE, denomEE, "c1", "fake rate", "", "");
  //SingleEff(hprobedEtaBNvtx,hprobeBNvtx,hprobedEtaENvtx,hprobeENvtx,"c2", "dEta cut eff","","");
  //SingleEff(hprobeTckIsoBNvtx,hprobeBNvtx,hprobeTckIsoENvtx,hprobeENvtx,"c3", "TckIso cut eff","","");
  //SingleEff(hprobeEcalHcal1IsoBNvtx,hprobeBNvtx,hprobeEcalHcal1IsoENvtx,hprobeENvtx,"c4", "EcalHcal1Iso cut eff","","");


}//end of method





void FakeRate::SingleEff(vector<TH1F *> histoB,vector<TH1F *> histobisB,vector<TH1F *> histoE,vector<TH1F *> histobisE,TString canvasname, TString Title, TString XTitle,TString YTitle )
{
  TGraphAsymmErrors *hgB[2]= {new TGraphAsymmErrors,new TGraphAsymmErrors};
  TGraphAsymmErrors *hgE[2]={new TGraphAsymmErrors,new TGraphAsymmErrors};
  TMultiGraph *mgB = new TMultiGraph(); 
  TMultiGraph *mgE = new TMultiGraph(); 

  hgB[0]->BayesDivide(histoB[0],histobisB[0]);
  hgB[0]->SetMarkerStyle(20);
  //hgB[1]->BayesDivide(histoB[1],histobisB[1]);
  //hgB[1]->SetMarkerStyle(0);
  hgE[0]->BayesDivide(histoE[0],histobisE[0]);
  hgE[0]->SetMarkerStyle(20);
  //hgE[1]->BayesDivide(histoE[1],histobisE[1]);
  //hgE[1]->SetMarkerStyle(0);
  mgB->Add(hgB[0]);
  //mgB->Add(hgB[1]);
  mgE->Add(hgE[0]);
  //mgE->Add(hgE[1]);
  //hgB[1]->SetLineColor(4);
  //hgE[1]->SetLineColor(4);


  TCanvas *c1 = new TCanvas(canvasname, Title, 100, 100, 1200, 725);
  TLegend *legend = new TLegend(0.4,0.2,0.7,0.5);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(hgB[0], "Data","lep");
  //legend->AddEntry(hgB[1], "DY MC","lp");


  c1->Divide(2); 
  c1->cd(1); 
  //mgB->SetMinimum(valmin); 
  //mgB->SetMaximum(1.03); 
  mgB->SetTitle(TString(Title).Append(" (Barrel)")); 
  mgB->Draw("AP");
  legend->Draw("sames");
  c1->cd(2);
  //mgE->SetMinimum(valmin); 
  //mgE->SetMaximum(1.03); 
  mgE->SetTitle(TString(Title).Append(" (Endcaps)")); 
  mgE->Draw("AP");
  legend->Draw("sames");
  if(saveplots){
    TImage *img = TImage::Create();
    img->FromPad(c1);
    TString pngcanvasfilename(TString(Title.Data()).Append(".png"));
    img->WriteImage(pngcanvasfilename.Data());
    
  }
}







bool FakeRate::IsVBTF(float pt,float eta, float dEta,float dPhi, float tckiso, float ecaliso, float hcaliso1,float hcaliso2,float sigmaietaieta, float hovere)
{ 

  bool isit = true; 
  //Barrel
  if(fabs(eta) <1.442 ){
    if((double)tckiso/(double)pt >0.09 )isit = false; 
    if((double)ecaliso/(double)pt>0.07) isit = false;
    if((double)hcaliso1/(double)pt +(double)hcaliso2/(double)pt >0.10 ) isit = false;
    if((double)sigmaietaieta>0.01) isit = false;
    if(fabs((double)dPhi) >0.06)isit = false;
    if(fabs((double)dEta) >0.004)isit = false;
    if((double)hovere>0.04)isit = false;
    if(dEta <0)cout << dEta << " dEta < 0 !!! " << endl;

  }
  //Endcaps
  if(fabs(eta) >1.442 ){
    if((double)tckiso/(double)pt >0.04 )isit = false; 
    if((double)ecaliso/(double)pt>0.05) isit = false;
    if((((double)hcaliso1 +(double)hcaliso2))/(double)pt >0.025 ) isit = false;
    if((double)sigmaietaieta>0.03) isit = false;
    if(fabs((double)dPhi) >0.03)isit = false;
    if(fabs((double)dEta) >0.007)isit = false;
    if((double)hovere>0.025)isit = false;
  }
  
  return isit;
}

bool FakeRate::IsHEEP(float pt,float eta, float dEta,float dPhi, float tckiso, float ecaliso, float hcaliso1,float hcaliso2,float sigmaietaieta, float hovere, float e2x5overe5x5,float e1x5overe5x5 )
{ 

  bool isit = true; 
  //Barrel
  if(fabs(eta) <1.442 ){
    if((double)tckiso > 7.5 )isit = false; 
    if(((double)ecaliso +(double)hcaliso1) >(2+0.03*(double)pt)) isit = false;
    if((double)e2x5overe5x5 < 0.94 && (double)e1x5overe5x5<0.83 ) isit = false;
    //if(sigmaietaieta>0.01) isit=false ; 
    if(fabs((double)dPhi) >0.09)isit = false;
       if(fabs((double)dEta) >0.005)isit = false;
    if((double)hovere>0.05)isit = false;
  }
  //Endcaps
  if(fabs(eta) >1.442 ){
    if((double)tckiso >15 )isit = false; 
    if((double)pt>50 && (((double)ecaliso +(double)hcaliso1 )>(1+0.03*(double)pt))) isit = false;
    if((double)pt<50 && (((double)ecaliso +(double)hcaliso1 )>2.5)) isit = false;
    if((double)hcaliso2 >0.5) isit = false;
    if((double)sigmaietaieta>0.03) isit = false;
    if(fabs((double)dPhi) >0.09)isit = false;
    if(fabs((double)dEta) >0.007)isit = false;
    if((double)hovere>0.05)isit = false;

  }
  
  return isit;
}
