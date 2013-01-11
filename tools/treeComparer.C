#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>

#include "TLorentzVector.h"
#include "TString.h"

//void TreeComparer::Loop()
void treeComparer(void)
{
  // parameters //////////////////////////////////////////////////////////////
  const int nBins = 100;
  const float massMin = 0.;
  const float massMax = 1500.;
  const float etMin = 0.;
  const float etMax = 500.;
  const float ptMin = 0.;
  const float ptMax = 500.;

  const int font = 42;
  ////////////////////////////////////////////////////////////////////////////

  TH1::SetDefaultSumw2(kTRUE);

  /////////////////////////////////////////////////////////////////////////
  // input files
  /////////////////////////////////////////////////////////////////////////
  vector<pair<TFile *, double> > input;
  vector<TString> inFileTag;
  vector<TString> legendName;
  //TFile *inRef = TFile::Open("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/mcsamples2012/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_35_3297045ev.root");
  //TFile *inRef = TFile::Open("file:////user/treis/mcsamples/ZprimePSIToEE_M-3000_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_25280ev.root");
  TFile *inRef = TFile::Open("file:////user/treis/mcsamples/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_28150723ev.root");
  input.push_back(make_pair(inRef, 1.));
  inFileTag.push_back("Ref");
  legendName.push_back("t#bar{t} powheg");

  //TFile *inTest = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1_AODSIM_gct1_33_3297045ev.root");
  //input.push_back(make_pair(inTest, 1.));
  //TFile *inTest = TFile::Open("file:////user/treis/mcsamples/DYToEE_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_gct1_45_99993ev.root");
  TFile *inTest = TFile::Open("file:////user/treis/mcsamples/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1+v2_AODSIM_gct1_46_16365457ev.root");
  input.push_back(make_pair(inTest, 1.));
  inFileTag.push_back("Test1");
  legendName.push_back("t#bar{t}+Jets MadGraph");
  ////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////////////
  // define test histograms
  /////////////////////////////////////////////////////////////////////////
  vector<TH1F *> hInvMass;
  vector<TH1F *> hGsfGsfEt;
  vector<TH1F *> hGsfPt;
  vector<TH1F *> hGsfEta;
  vector<TH1F *> hGsfPhi;
  vector<TH1I *> hGsfSize;
  vector<TH1I *> hGenMomPdgId;
  vector<TH1F *> hGenMomMass;
  vector<TH1F *> hGenMomPt;
  vector<TH1F *> hGenMomEta;
  vector<TH1F *> hGenMomPhi;
  /////////////////////////////////////////////////////////////////////////

  TCanvas *c0 = new TCanvas("invMass", "Invariant Mass", 100, 100, 800, 600);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFillColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(font);
  gStyle->SetLabelFont(font);
  gStyle->SetLegendFont(font);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.);
  gPad->SetTicks(1, 1);

  c0->Divide(1, 2);
  c0->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c0->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c0->cd(1);
  c0->SetBorderMode(0);
  c0->SetBorderSize(2);

 /////////////////////////////////////////////////////////////////////////
  // loop over files
  /////////////////////////////////////////////////////////////////////////
  for (unsigned int p = 0; p < input.size(); ++p) {
    string infile(input[p].first->GetName());
    cout << "file " << infile << endl;
    input[p].first->Cd("");
    TTree *thetree = (TTree*)input[p].first->Get("gsfcheckerjob/tree");

    // fill the histograms from the tree
    hInvMass.push_back(new TH1F("hInvMass" + inFileTag[p], "hInvMass" + inFileTag[p], nBins, massMin, massMax));
    cout << "Generating " << hInvMass.back()->GetName() << " Histogram.";
    thetree->Draw("heepHeepMass>>hInvMass" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGsfGsfEt.push_back(new TH1F("hGsfGsfEt" + inFileTag[p], "hGsfGsfEt" + inFileTag[p], nBins, etMin, etMax));
    cout << "Generating " << hGsfGsfEt.back()->GetName() << " Histogram.";
    thetree->Draw("gsf_gsfet>>hGsfGsfEt" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGsfPt.push_back(new TH1F("hGsfPt" + inFileTag[p], "hGsfPt" + inFileTag[p], nBins, ptMin, ptMax));
    cout << "Generating " << hGsfPt.back()->GetName() << " Histogram.";
    thetree->Draw("gsf_pt>>hGsfPt" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGsfEta.push_back(new TH1F("hGsfEta" + inFileTag[p], "hGsfEta" + inFileTag[p], nBins, -3, 3));
    cout << "Generating " << hGsfEta.back()->GetName() << " Histogram.";
    thetree->Draw("gsf_eta>>hGsfEta" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGsfPhi.push_back(new TH1F("hGsfPhi" + inFileTag[p], "hGsfPhi" + inFileTag[p], nBins, -4., 4.));
    cout << "Generating " << hGsfPhi.back()->GetName() << " Histogram.";
    thetree->Draw("gsf_phi>>hGsfPhi" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGsfSize.push_back(new TH1I("hGsfSize" + inFileTag[p], "hGsfSize" + inFileTag[p], 20, 0, 20));
    cout << "Generating " << hGsfSize.back()->GetName() << " Histogram.";
    thetree->Draw("gsf_size>>hGsfSize" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGenMomPdgId.push_back(new TH1I("hGenMomPdgId" + inFileTag[p], "hGenMomPdgId" + inFileTag[p], 100, -50, 50));
    cout << "Generating " << hGenMomPdgId.back()->GetName() << " Histogram.";
    thetree->Draw("genelemom_pdgid>>hGenMomPdgId" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGenMomMass.push_back(new TH1F("hGenMomMass" + inFileTag[p], "hGenMomMass" + inFileTag[p], nBins, massMin, massMax));
    cout << "Generating " << hGenMomMass.back()->GetName() << " Histogram.";
    thetree->Draw("genelemom_mass>>hGenMomMass" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGenMomPt.push_back(new TH1F("hGenMomPt" + inFileTag[p], "hGenMomPt" + inFileTag[p], nBins, ptMin, ptMax));
    cout << "Generating " << hGenMomPt.back()->GetName() << " Histogram.";
    thetree->Draw("genelemom_pt>>hGenMomPt" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGenMomEta.push_back(new TH1F("hGenMomEta" + inFileTag[p], "hGenMomEta" + inFileTag[p], nBins, -3., 3.));
    cout << "Generating " << hGenMomEta.back()->GetName() << " Histogram.";
    thetree->Draw("genelemom_eta>>hGenMomEta" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    hGenMomPhi.push_back(new TH1F("hGenMomPhi" + inFileTag[p], "hGenMomPhi" + inFileTag[p], nBins, -4., 4.));
    cout << "Generating " << hGenMomPhi.back()->GetName() << " Histogram.";
    thetree->Draw("genelemom_phi>>hGenMomPhi" + inFileTag[p], "gsf_gsfet[0]>25. && gsf_gsfet[1]>25.");
    cout << "       Done" << endl;

    // scale the test histogram to the ref
    if (p > 0) {
      hInvMass.back()->Scale(hInvMass[0]->Integral() / hInvMass.back()->Integral());
      hGsfGsfEt.back()->Scale(hGsfGsfEt[0]->Integral() / hGsfGsfEt.back()->Integral());
      hGsfPt.back()->Scale(hGsfPt[0]->Integral() / hGsfPt.back()->Integral());
      hGsfEta.back()->Scale(hGsfEta[0]->Integral() / hGsfEta.back()->Integral());
      hGsfPhi.back()->Scale(hGsfPhi[0]->Integral() / hGsfPhi.back()->Integral());
      hGsfSize.back()->Scale(hGsfSize[0]->Integral() / hGsfSize.back()->Integral());
      hGenMomPdgId.back()->Scale(hGenMomPdgId[0]->Integral() / hGenMomPdgId.back()->Integral());
      hGenMomMass.back()->Scale(hGenMomMass[0]->Integral() / hGenMomMass.back()->Integral());
      hGenMomPt.back()->Scale(hGenMomPt[0]->Integral() / hGenMomPt.back()->Integral());
      hGenMomEta.back()->Scale(hGenMomEta[0]->Integral() / hGenMomEta.back()->Integral());
      hGenMomPhi.back()->Scale(hGenMomPhi[0]->Integral() / hGenMomPhi.back()->Integral());
    }

  } // end of loop over input files

  TH1F * hInvMassDiff = (TH1F*)hInvMass[0]->Clone("hInvMassDiff");
  TH1F * hGsfGsfEtDiff = (TH1F*)hGsfGsfEt[0]->Clone("hGsfGsfEtDiff");
  TH1F * hGsfPtDiff = (TH1F*)hGsfPt[0]->Clone("hGsfPtDiff");
  TH1F * hGsfEtaDiff = (TH1F*)hGsfEta[0]->Clone("hGsfEtaDiff");
  TH1F * hGsfPhiDiff = (TH1F*)hGsfPhi[0]->Clone("hGsfPhiDiff");
  TH1F * hGsfSizeDiff = (TH1F*)hGsfSize[0]->Clone("hGsfSizeDiff");
  TH1F * hGenMomPdgIdDiff = (TH1F*)hGenMomPdgId[0]->Clone("hGenMomPdgIdDiff");
  TH1F * hGenMomMassDiff = (TH1F*)hGenMomMass[0]->Clone("hGenMomMassDiff");
  TH1F * hGenMomPtDiff = (TH1F*)hGenMomPt[0]->Clone("hGenMomPtDiff");
  TH1F * hGenMomEtaDiff = (TH1F*)hGenMomEta[0]->Clone("hGenMomEtaDiff");
  TH1F * hGenMomPhiDiff = (TH1F*)hGenMomPhi[0]->Clone("hGenMomPhiDiff");

  hInvMassDiff->Scale(-1.);
  hInvMassDiff->Add(hInvMass[1]);
  hInvMassDiff->Divide(hInvMass[0]);
  hGsfGsfEtDiff->Scale(-1.);
  hGsfGsfEtDiff->Add(hGsfGsfEt[1]);
  hGsfGsfEtDiff->Divide(hGsfGsfEt[0]);
  hGsfPtDiff->Scale(-1.);
  hGsfPtDiff->Add(hGsfPt[1]);
  hGsfPtDiff->Divide(hGsfPt[0]);
  hGsfEtaDiff->Scale(-1.);
  hGsfEtaDiff->Add(hGsfEta[1]);
  hGsfEtaDiff->Divide(hGsfEta[0]);
  hGsfPhiDiff->Scale(-1.);
  hGsfPhiDiff->Add(hGsfPhi[1]);
  hGsfPhiDiff->Divide(hGsfPhi[0]);
  hGsfSizeDiff->Scale(-1.);
  hGsfSizeDiff->Add(hGsfSize[1]);
  hGsfSizeDiff->Divide(hGsfSize[0]);
  hGenMomPdgIdDiff->Scale(-1.);
  hGenMomPdgIdDiff->Add(hGenMomPdgId[1]);
  hGenMomPdgIdDiff->Divide(hGenMomPdgId[0]);
  hGenMomMassDiff->Scale(-1.);
  hGenMomMassDiff->Add(hGenMomMass[1]);
  hGenMomMassDiff->Divide(hGenMomMass[0]);
  hGenMomPtDiff->Scale(-1.);
  hGenMomPtDiff->Add(hGenMomPt[1]);
  hGenMomPtDiff->Divide(hGenMomPt[0]);
  hGenMomEtaDiff->Scale(-1.);
  hGenMomEtaDiff->Add(hGenMomEta[1]);
  hGenMomEtaDiff->Divide(hGenMomEta[0]);
  hGenMomPhiDiff->Scale(-1.);
  hGenMomPhiDiff->Add(hGenMomPhi[1]);
  hGenMomPhiDiff->Divide(hGenMomPhi[0]);

  TCanvas *c1 = new TCanvas("gsfGsfEt", "gsfGsfEt", 100, 100, 800, 600);
  c1->Divide(1, 2);
  c1->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c1->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c1->cd(1);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  TCanvas *c2 = new TCanvas("gsfPt", "gsfPt", 100, 100, 800, 600);
  c2->Divide(1, 2);
  c2->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c2->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c2->cd(1);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  TCanvas *c3 = new TCanvas("gsfEta", "gsfEta", 100, 100, 800, 600);
  c3->Divide(1, 2);
  c3->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c3->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c3->cd(1);
  c3->SetBorderMode(0);
  c3->SetBorderSize(2);
  TCanvas *c4 = new TCanvas("gsfPhi", "gsfPhi", 100, 100, 800, 600);
  c4->Divide(1, 2);
  c4->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c4->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c4->cd(1);
  c4->SetBorderMode(0);
  c4->SetBorderSize(2);
  TCanvas *c5 = new TCanvas("gsfiSize", "gsfSize", 100, 100, 800, 600);
  c5->Divide(1, 2);
  c5->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c5->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c5->cd(1);
  c5->SetBorderMode(0);
  c5->SetBorderSize(2);
  TCanvas *c6 = new TCanvas("genMomPdgId", "genMomPdgId", 100, 100, 800, 600);
  c6->Divide(1, 2);
  c6->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c6->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c6->cd(1);
  c6->SetBorderMode(0);
  c6->SetBorderSize(2);
  TCanvas *c7 = new TCanvas("genMomMass", "genMomMass", 100, 100, 800, 600);
  c7->Divide(1, 2);
  c7->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c7->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c7->cd(1);
  c7->SetBorderMode(0);
  c7->SetBorderSize(2);
  TCanvas *c8 = new TCanvas("genMomPt", "genMomPt", 100, 100, 800, 600);
  c8->Divide(1, 2);
  c8->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c8->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c8->cd(1);
  c8->SetBorderMode(0);
  c8->SetBorderSize(2);
  TCanvas *c9 = new TCanvas("genMomEta", "genMomEta", 100, 100, 800, 600);
  c9->Divide(1, 2);
  c9->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c9->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c9->cd(1);
  c9->SetBorderMode(0);
  c9->SetBorderSize(2);
  TCanvas *c10 = new TCanvas("genMomPhi", "genMomPhi", 100, 100, 800, 600);
  c10->Divide(1, 2);
  c10->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
  c10->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
  c10->cd(1);
  c10->SetBorderMode(0);
  c10->SetBorderSize(2);

  ////////////////////////////////////////////////////////////////////////////
  c0->cd(1);
  gPad->SetLogy();
  hInvMass.front()->SetLineColor(kBlue);
  hInvMass.front()->SetMarkerColor(kBlue);
  hInvMass.front()->SetMarkerStyle(20);
  hInvMass.front()->SetMarkerSize(0.8);
  hInvMass.front()->GetXaxis()->SetTitle("m(ee) [GeV]");
  //hInvMass.front()->GetYaxis()->SetTitle("Events");
  hInvMass.front()->GetYaxis()->SetLabelFont(font);
  hInvMass.front()->Draw();

  TLegend *legend1 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend1->SetTextFont(font);
  legend1->SetTextSize(0.03);
  legend1->SetBorderSize(0);
  legend1->AddEntry(hInvMass.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hInvMass.at(p)->SetLineColor(kRed + p - 1);
    hInvMass.at(p)->SetMarkerColor(kRed + p - 1);
    hInvMass.at(p)->SetMarkerStyle(24);
    hInvMass.at(p)->SetMarkerSize(0.8);
    hInvMass.at(p)->Draw("sames");
    legend1->AddEntry(hInvMass.at(p), legendName[p].Data());
  }  
  legend1->Draw("sames");

  c0->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hInvMassDiff->Draw();
  TLine *zero = new TLine(massMin, 0., massMax, 0.);
  zero->Draw();
  hInvMassDiff->SetLineColor(kBlue);
  hInvMassDiff->SetMarkerColor(kBlue);
  hInvMassDiff->SetMarkerStyle(24);
  hInvMassDiff->SetMarkerSize(0.8);
  hInvMassDiff->GetYaxis()->SetLabelFont(font);
  hInvMassDiff->GetYaxis()->SetNdivisions(204);
  hInvMassDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hInvMassDiff->GetYaxis()->SetTitleFont(font);
  hInvMassDiff->GetYaxis()->SetTitleSize(0.15);
  hInvMassDiff->GetYaxis()->SetTitleOffset(0.2);
  hInvMassDiff->GetYaxis()->SetLabelSize(0.15);
  hInvMassDiff->GetXaxis()->SetTitle("m(ee) [GeV]");
  hInvMassDiff->GetXaxis()->SetTitleSize(0.15);
  hInvMassDiff->GetXaxis()->SetLabelSize(0.15);
  hInvMassDiff->Draw("sameaxis");
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c1->cd(1);
  gPad->SetLogy();
  hGsfGsfEt.front()->SetLineColor(kBlue);
  hGsfGsfEt.front()->SetMarkerColor(kBlue);
  hGsfGsfEt.front()->SetMarkerStyle(20);
  hGsfGsfEt.front()->SetMarkerSize(0.8);
  hGsfGsfEt.front()->GetXaxis()->SetTitle("Et [GeV]");
  hGsfGsfEt.front()->GetYaxis()->SetLabelFont(font);
  hGsfGsfEt.front()->Draw();

  TLegend *legend = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend->SetTextFont(font);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(hGsfGsfEt.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGsfGsfEt.at(p)->SetLineColor(kRed + p - 1);
    hGsfGsfEt.at(p)->SetMarkerColor(kRed + p - 1);
    hGsfGsfEt.at(p)->SetMarkerStyle(24);
    hGsfGsfEt.at(p)->SetMarkerSize(0.8);
    hGsfGsfEt.at(p)->Draw("sames");
    legend->AddEntry(hGsfGsfEt.at(p), legendName[p].Data());
  }  
  legend->Draw("sames");

  c1->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGsfGsfEtDiff->GetYaxis()->SetLabelFont(font);
  hGsfGsfEtDiff->Draw();
  zero->DrawLine(etMin, 0., etMax, 0.);
  hGsfGsfEtDiff->GetYaxis()->SetNdivisions(204);
  hGsfGsfEtDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGsfGsfEtDiff->GetYaxis()->SetTitleFont(font);
  hGsfGsfEtDiff->GetYaxis()->SetTitleSize(0.15);
  hGsfGsfEtDiff->GetYaxis()->SetTitleOffset(0.2);
  hGsfGsfEtDiff->GetYaxis()->SetLabelSize(0.15);
  hGsfGsfEtDiff->GetXaxis()->SetTitle("E_{T}^{RECO} [GeV]");
  hGsfGsfEtDiff->GetXaxis()->SetTitleSize(0.15);
  hGsfGsfEtDiff->GetXaxis()->SetLabelSize(0.15);
  hGsfGsfEtDiff->SetLineColor(kBlue);
  hGsfGsfEtDiff->SetMarkerColor(kBlue);
  hGsfGsfEtDiff->SetMarkerStyle(24);
  hGsfGsfEtDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c2->cd(1);
  gPad->SetLogy();
  hGsfPt.front()->SetLineColor(kBlue);
  hGsfPt.front()->SetMarkerColor(kBlue);
  hGsfPt.front()->SetMarkerStyle(20);
  hGsfPt.front()->SetMarkerSize(0.8);
  hGsfPt.front()->GetXaxis()->SetTitle("pt [GeV]");
  hGsfPt.front()->GetYaxis()->SetLabelFont(font);
  hGsfPt.front()->Draw();

  TLegend *legend2 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend2->SetTextFont(font);
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->AddEntry(hGsfPt.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGsfPt.at(p)->SetLineColor(kRed + p - 1);
    hGsfPt.at(p)->SetMarkerColor(kRed + p - 1);
    hGsfPt.at(p)->SetMarkerStyle(24);
    hGsfPt.at(p)->SetMarkerSize(0.8);
    hGsfPt.at(p)->Draw("sames");
    legend2->AddEntry(hGsfPt.at(p), legendName[p].Data());
  }  
  legend2->Draw("sames");

  c2->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGsfPtDiff->GetYaxis()->SetLabelFont(font);
  hGsfPtDiff->Draw();
  zero->DrawLine(ptMin, 0., ptMax, 0.);
  hGsfPtDiff->GetYaxis()->SetNdivisions(204);
  hGsfPtDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGsfPtDiff->GetYaxis()->SetTitleFont(font);
  hGsfPtDiff->GetYaxis()->SetTitleSize(0.15);
  hGsfPtDiff->GetYaxis()->SetTitleOffset(0.2);
  hGsfPtDiff->GetYaxis()->SetLabelSize(0.15);
  hGsfPtDiff->GetXaxis()->SetTitle("p_{T}^{RECO} [GeV]");
  hGsfPtDiff->GetXaxis()->SetTitleSize(0.15);
  hGsfPtDiff->GetXaxis()->SetLabelSize(0.15);
  hGsfPtDiff->SetLineColor(kBlue);
  hGsfPtDiff->SetMarkerColor(kBlue);
  hGsfPtDiff->SetMarkerStyle(24);
  hGsfPtDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c3->cd(1);
  gPad->SetLogy();
  hGsfEta.front()->SetLineColor(kBlue);
  hGsfEta.front()->SetMarkerColor(kBlue);
  hGsfEta.front()->SetMarkerStyle(20);
  hGsfEta.front()->SetMarkerSize(0.8);
  hGsfEta.front()->GetXaxis()->SetTitle("#eta_{RECO} [deg]");
  hGsfEta.front()->GetYaxis()->SetLabelFont(font);
  hGsfEta.front()->Draw();

  TLegend *legend3 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend3->SetTextFont(font);
  legend3->SetTextSize(0.03);
  legend3->SetBorderSize(0);
  legend3->AddEntry(hGsfEta.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGsfEta.at(p)->SetLineColor(kRed + p - 1);
    hGsfEta.at(p)->SetMarkerColor(kRed + p - 1);
    hGsfEta.at(p)->SetMarkerStyle(24);
    hGsfEta.at(p)->SetMarkerSize(0.8);
    hGsfEta.at(p)->Draw("sames");
    legend3->AddEntry(hGsfEta.at(p), legendName[p].Data());
  }  
  legend3->Draw("sames");

  c3->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGsfEtaDiff->GetYaxis()->SetLabelFont(font);
  hGsfEtaDiff->Draw();
  zero->DrawLine(-3., 0., 3., 0.);
  hGsfEtaDiff->GetYaxis()->SetNdivisions(204);
  hGsfEtaDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGsfEtaDiff->GetYaxis()->SetTitleFont(font);
  hGsfEtaDiff->GetYaxis()->SetTitleSize(0.15);
  hGsfEtaDiff->GetYaxis()->SetTitleOffset(0.2);
  hGsfEtaDiff->GetYaxis()->SetLabelSize(0.15);
  hGsfEtaDiff->GetXaxis()->SetTitle("#eta_{RECO} [deg]");
  hGsfEtaDiff->GetXaxis()->SetTitleSize(0.15);
  hGsfEtaDiff->GetXaxis()->SetLabelSize(0.15);
  hGsfEtaDiff->SetLineColor(kBlue);
  hGsfEtaDiff->SetMarkerColor(kBlue);
  hGsfEtaDiff->SetMarkerStyle(24);
  hGsfEtaDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c4->cd(1);
  gPad->SetLogy();
  hGsfPhi.front()->SetLineColor(kBlue);
  hGsfPhi.front()->SetMarkerColor(kBlue);
  hGsfPhi.front()->SetMarkerStyle(20);
  hGsfPhi.front()->SetMarkerSize(0.8);
  hGsfPhi.front()->GetXaxis()->SetTitle("#phi_{RECO} [deg]");
  hGsfPhi.front()->GetYaxis()->SetLabelFont(font);
  hGsfPhi.front()->Draw();

  TLegend *legend4 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend4->SetTextFont(font);
  legend4->SetTextSize(0.03);
  legend4->SetBorderSize(0);
  legend4->AddEntry(hGsfPhi.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGsfPhi.at(p)->SetLineColor(kRed + p - 1);
    hGsfPhi.at(p)->SetMarkerColor(kRed + p - 1);
    hGsfPhi.at(p)->SetMarkerStyle(24);
    hGsfPhi.at(p)->SetMarkerSize(0.8);
    hGsfPhi.at(p)->Draw("sames");
    legend4->AddEntry(hGsfPhi.at(p), legendName[p].Data());
  }  
  legend4->Draw("sames");

  c4->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGsfPhiDiff->GetYaxis()->SetLabelFont(font);
  hGsfPhiDiff->Draw();
  zero->DrawLine(-4., 0., 4., 0.);
  hGsfPhiDiff->GetYaxis()->SetNdivisions(204);
  hGsfPhiDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGsfPhiDiff->GetYaxis()->SetTitleFont(font);
  hGsfPhiDiff->GetYaxis()->SetTitleSize(0.15);
  hGsfPhiDiff->GetYaxis()->SetTitleOffset(0.2);
  hGsfPhiDiff->GetYaxis()->SetLabelSize(0.15);
  hGsfPhiDiff->GetXaxis()->SetTitle("#phi_{RECO} [deg]");
  hGsfPhiDiff->GetXaxis()->SetTitleSize(0.15);
  hGsfPhiDiff->GetXaxis()->SetLabelSize(0.15);
  hGsfPhiDiff->SetLineColor(kBlue);
  hGsfPhiDiff->SetMarkerColor(kBlue);
  hGsfPhiDiff->SetMarkerStyle(24);
  hGsfPhiDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c5->cd(1);
  gPad->SetLogy();
  hGsfSize.front()->SetLineColor(kBlue);
  hGsfSize.front()->SetMarkerColor(kBlue);
  hGsfSize.front()->SetMarkerStyle(20);
  hGsfSize.front()->SetMarkerSize(0.8);
  hGsfSize.front()->GetXaxis()->SetTitle("#");
  hGsfSize.front()->GetYaxis()->SetLabelFont(font);
  hGsfSize.front()->Draw();

  TLegend *legend5 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend5->SetTextFont(font);
  legend5->SetTextSize(0.03);
  legend5->SetBorderSize(0);
  legend5->AddEntry(hGsfSize.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGsfSize.at(p)->SetLineColor(kRed + p - 1);
    hGsfSize.at(p)->SetMarkerColor(kRed + p - 1);
    hGsfSize.at(p)->SetMarkerStyle(24);
    hGsfSize.at(p)->SetMarkerSize(0.8);
    hGsfSize.at(p)->Draw("sames");
    legend5->AddEntry(hGsfSize.at(p), legendName[p].Data());
  }  
  legend5->Draw("sames");

  c5->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGsfSizeDiff->GetYaxis()->SetLabelFont(font);
  hGsfSizeDiff->Draw();
  zero->DrawLine(0., 0., 20., 0.);
  hGsfSizeDiff->GetYaxis()->SetNdivisions(204);
  hGsfSizeDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGsfSizeDiff->GetYaxis()->SetTitleFont(font);
  hGsfSizeDiff->GetYaxis()->SetTitleSize(0.15);
  hGsfSizeDiff->GetYaxis()->SetTitleOffset(0.2);
  hGsfSizeDiff->GetYaxis()->SetLabelSize(0.15);
  hGsfSizeDiff->GetXaxis()->SetTitle("#");
  hGsfSizeDiff->GetXaxis()->SetTitleSize(0.15);
  hGsfSizeDiff->GetXaxis()->SetLabelSize(0.15);
  hGsfSizeDiff->SetLineColor(kBlue);
  hGsfSizeDiff->SetMarkerColor(kBlue);
  hGsfSizeDiff->SetMarkerStyle(24);
  hGsfSizeDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c6->cd(1);
  gPad->SetLogy();
  hGenMomPdgId.front()->SetLineColor(kBlue);
  hGenMomPdgId.front()->SetMarkerColor(kBlue);
  hGenMomPdgId.front()->SetMarkerStyle(20);
  hGenMomPdgId.front()->SetMarkerSize(0.8);
  hGenMomPdgId.front()->GetXaxis()->SetTitle("PDG id");
  hGenMomPdgId.front()->GetYaxis()->SetLabelFont(font);
  hGenMomPdgId.front()->Draw();

  TLegend *legend6 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend6->SetTextFont(font);
  legend6->SetTextSize(0.03);
  legend6->SetBorderSize(0);
  legend6->AddEntry(hGenMomPdgId.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGenMomPdgId.at(p)->SetLineColor(kRed + p - 1);
    hGenMomPdgId.at(p)->SetMarkerColor(kRed + p - 1);
    hGenMomPdgId.at(p)->SetMarkerStyle(24);
    hGenMomPdgId.at(p)->SetMarkerSize(0.8);
    hGenMomPdgId.at(p)->Draw("sames");
    legend6->AddEntry(hGenMomPdgId.at(p), legendName[p].Data());
  }  
  legend6->Draw("sames");

  c6->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGenMomPdgIdDiff->GetYaxis()->SetLabelFont(font);
  hGenMomPdgIdDiff->Draw();
  zero->DrawLine(-50., 0., 50., 0.);
  hGenMomPdgIdDiff->GetYaxis()->SetNdivisions(204);
  hGenMomPdgIdDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGenMomPdgIdDiff->GetYaxis()->SetTitleFont(font);
  hGenMomPdgIdDiff->GetYaxis()->SetTitleSize(0.15);
  hGenMomPdgIdDiff->GetYaxis()->SetTitleOffset(0.2);
  hGenMomPdgIdDiff->GetYaxis()->SetLabelSize(0.15);
  hGenMomPdgIdDiff->GetXaxis()->SetTitle("PDG id");
  hGenMomPdgIdDiff->GetXaxis()->SetTitleSize(0.15);
  hGenMomPdgIdDiff->GetXaxis()->SetLabelSize(0.15);
  hGenMomPdgIdDiff->SetLineColor(kBlue);
  hGenMomPdgIdDiff->SetMarkerColor(kBlue);
  hGenMomPdgIdDiff->SetMarkerStyle(24);
  hGenMomPdgIdDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c7->cd(1);
  gPad->SetLogy();
  hGenMomMass.front()->SetLineColor(kBlue);
  hGenMomMass.front()->SetMarkerColor(kBlue);
  hGenMomMass.front()->SetMarkerStyle(20);
  hGenMomMass.front()->SetMarkerSize(0.8);
  hGenMomMass.front()->GetXaxis()->SetTitle("m_{true} [GeV]");
  hGenMomMass.front()->GetYaxis()->SetLabelFont(font);
  hGenMomMass.front()->Draw();

  TLegend *legend7 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend7->SetTextFont(font);
  legend7->SetTextSize(0.03);
  legend7->SetBorderSize(0);
  legend7->AddEntry(hGenMomMass.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGenMomMass.at(p)->SetLineColor(kRed + p - 1);
    hGenMomMass.at(p)->SetMarkerColor(kRed + p - 1);
    hGenMomMass.at(p)->SetMarkerStyle(24);
    hGenMomMass.at(p)->SetMarkerSize(0.8);
    hGenMomMass.at(p)->Draw("sames");
    legend7->AddEntry(hGenMomMass.at(p), legendName[p].Data());
  }  
  legend7->Draw("sames");

  c7->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGenMomMassDiff->GetYaxis()->SetLabelFont(font);
  hGenMomMassDiff->Draw();
  zero->DrawLine(massMin, 0., massMax, 0.);
  hGenMomMassDiff->GetYaxis()->SetNdivisions(204);
  hGenMomMassDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGenMomMassDiff->GetYaxis()->SetTitleFont(font);
  hGenMomMassDiff->GetYaxis()->SetTitleSize(0.15);
  hGenMomMassDiff->GetYaxis()->SetTitleOffset(0.2);
  hGenMomMassDiff->GetYaxis()->SetLabelSize(0.15);
  hGenMomMassDiff->GetXaxis()->SetTitle("m_{true} [GeV]");
  hGenMomMassDiff->GetXaxis()->SetTitleSize(0.15);
  hGenMomMassDiff->GetXaxis()->SetLabelSize(0.15);
  hGenMomMassDiff->SetLineColor(kBlue);
  hGenMomMassDiff->SetMarkerColor(kBlue);
  hGenMomMassDiff->SetMarkerStyle(24);
  hGenMomMassDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c8->cd(1);
  gPad->SetLogy();
  hGenMomPt.front()->SetLineColor(kBlue);
  hGenMomPt.front()->SetMarkerColor(kBlue);
  hGenMomPt.front()->SetMarkerStyle(20);
  hGenMomPt.front()->SetMarkerSize(0.8);
  hGenMomPt.front()->GetXaxis()->SetTitle("p_{T}^{true} [GeV]");
  hGenMomPt.front()->GetYaxis()->SetLabelFont(font);
  hGenMomPt.front()->Draw();

  TLegend *legend8 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend8->SetTextFont(font);
  legend8->SetTextSize(0.03);
  legend8->SetBorderSize(0);
  legend8->AddEntry(hGenMomPt.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGenMomPt.at(p)->SetLineColor(kRed + p - 1);
    hGenMomPt.at(p)->SetMarkerColor(kRed + p - 1);
    hGenMomPt.at(p)->SetMarkerStyle(24);
    hGenMomPt.at(p)->SetMarkerSize(0.8);
    hGenMomPt.at(p)->Draw("sames");
    legend8->AddEntry(hGenMomPt.at(p), legendName[p].Data());
  }  
  legend8->Draw("sames");

  c8->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGenMomPtDiff->GetYaxis()->SetLabelFont(font);
  hGenMomPtDiff->Draw();
  zero->DrawLine(ptMin, 0., ptMax, 0.);
  hGenMomPtDiff->GetYaxis()->SetNdivisions(204);
  hGenMomPtDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGenMomPtDiff->GetYaxis()->SetTitleFont(font);
  hGenMomPtDiff->GetYaxis()->SetTitleSize(0.15);
  hGenMomPtDiff->GetYaxis()->SetTitleOffset(0.2);
  hGenMomPtDiff->GetYaxis()->SetLabelSize(0.15);
  hGenMomPtDiff->GetXaxis()->SetTitle("p_{T}^{true} [GeV]");
  hGenMomPtDiff->GetXaxis()->SetTitleSize(0.15);
  hGenMomPtDiff->GetXaxis()->SetLabelSize(0.15);
  hGenMomPtDiff->SetLineColor(kBlue);
  hGenMomPtDiff->SetMarkerColor(kBlue);
  hGenMomPtDiff->SetMarkerStyle(24);
  hGenMomPtDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c9->cd(1);
  gPad->SetLogy();
  hGenMomEta.front()->SetLineColor(kBlue);
  hGenMomEta.front()->SetMarkerColor(kBlue);
  hGenMomEta.front()->SetMarkerStyle(20);
  hGenMomEta.front()->SetMarkerSize(0.8);
  hGenMomEta.front()->GetXaxis()->SetTitle("#eta_{true} [deg]");
  hGenMomEta.front()->GetYaxis()->SetLabelFont(font);
  hGenMomEta.front()->Draw();

  TLegend *legend9 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend9->SetTextFont(font);
  legend9->SetTextSize(0.03);
  legend9->SetBorderSize(0);
  legend9->AddEntry(hGenMomEta.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGenMomEta.at(p)->SetLineColor(kRed + p - 1);
    hGenMomEta.at(p)->SetMarkerColor(kRed + p - 1);
    hGenMomEta.at(p)->SetMarkerStyle(24);
    hGenMomEta.at(p)->SetMarkerSize(0.8);
    hGenMomEta.at(p)->Draw("sames");
    legend9->AddEntry(hGenMomEta.at(p), legendName[p].Data());
  }  
  legend9->Draw("sames");

  c9->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGenMomEtaDiff->GetYaxis()->SetLabelFont(font);
  hGenMomEtaDiff->Draw();
  zero->DrawLine(-3., 0., 3., 0.);
  hGenMomEtaDiff->GetYaxis()->SetNdivisions(204);
  hGenMomEtaDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGenMomEtaDiff->GetYaxis()->SetTitleFont(font);
  hGenMomEtaDiff->GetYaxis()->SetTitleSize(0.15);
  hGenMomEtaDiff->GetYaxis()->SetTitleOffset(0.2);
  hGenMomEtaDiff->GetYaxis()->SetLabelSize(0.15);
  hGenMomEtaDiff->GetXaxis()->SetTitle("#eta_{true} [deg]");
  hGenMomEtaDiff->GetXaxis()->SetTitleSize(0.15);
  hGenMomEtaDiff->GetXaxis()->SetLabelSize(0.15);
  hGenMomEtaDiff->SetLineColor(kBlue);
  hGenMomEtaDiff->SetMarkerColor(kBlue);
  hGenMomEtaDiff->SetMarkerStyle(24);
  hGenMomEtaDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  c10->cd(1);
  gPad->SetLogy();
  hGenMomPhi.front()->SetLineColor(kBlue);
  hGenMomPhi.front()->SetMarkerColor(kBlue);
  hGenMomPhi.front()->SetMarkerStyle(20);
  hGenMomPhi.front()->SetMarkerSize(0.8);
  hGenMomPhi.front()->GetXaxis()->SetTitle("#phi_{true} [deg]");
  hGenMomPhi.front()->GetYaxis()->SetLabelFont(font);
  hGenMomPhi.front()->Draw();

  TLegend *legend10 = new TLegend(0.65, 0.68, 0.82, 0.82);
  legend10->SetTextFont(font);
  legend10->SetTextSize(0.03);
  legend10->SetBorderSize(0);
  legend10->AddEntry(hGenMomPhi.front(), legendName[0].Data());
  for (unsigned int p = 1; p < input.size(); ++p) {
    hGenMomPhi.at(p)->SetLineColor(kRed + p - 1);
    hGenMomPhi.at(p)->SetMarkerColor(kRed + p - 1);
    hGenMomPhi.at(p)->SetMarkerStyle(24);
    hGenMomPhi.at(p)->SetMarkerSize(0.8);
    hGenMomPhi.at(p)->Draw("sames");
    legend10->AddEntry(hGenMomPhi.at(p), legendName[p].Data());
  }  
  legend10->Draw("sames");

  c10->cd(2);
  gPad->SetTopMargin(0.);
  gPad->SetBottomMargin(0.3);
  gPad->SetBorderMode(0);
  hGenMomPhiDiff->GetYaxis()->SetLabelFont(font);
  hGenMomPhiDiff->Draw();
  zero->DrawLine(-4., 0., 4., 0.);
  hGenMomPhiDiff->GetYaxis()->SetNdivisions(204);
  hGenMomPhiDiff->GetYaxis()->SetTitle("(test- ref)/ref");
  hGenMomPhiDiff->GetYaxis()->SetTitleFont(font);
  hGenMomPhiDiff->GetYaxis()->SetTitleSize(0.15);
  hGenMomPhiDiff->GetYaxis()->SetTitleOffset(0.2);
  hGenMomPhiDiff->GetYaxis()->SetLabelSize(0.15);
  hGenMomPhiDiff->GetXaxis()->SetTitle("#phi_{true} [deg]");
  hGenMomPhiDiff->GetXaxis()->SetTitleSize(0.15);
  hGenMomPhiDiff->GetXaxis()->SetLabelSize(0.15);
  hGenMomPhiDiff->SetLineColor(kBlue);
  hGenMomPhiDiff->SetMarkerColor(kBlue);
  hGenMomPhiDiff->SetMarkerStyle(24);
  hGenMomPhiDiff->SetMarkerSize(0.8);
  ////////////////////////////////////////////////////////////////////////////


} //end of method



