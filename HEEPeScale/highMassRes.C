#include "highMassRes.h"

RooWorkspace * HighMassRes::cryBall(
                   TTree *tree,
                   const int regions,
                   double minMass,
                   double maxMass,
                   double trueMassMin,
                   double trueMassMax,
                   double cutoff_cb,
                   const char* plotOpt,
                   const int nBins,
                   const unsigned int fitModelType) // (0, 1) = (CB, double CB)
{
  setTDRStyle();

  // make a workspace to return
  RooWorkspace *w = new RooWorkspace("w","workspace");

  //RooRealVar trueMass("trueMass", "m_{true}", minMass, maxMass, "GeV");
  RooRealVar trueMass("trueMass", "m_{true}", trueMassMin, trueMassMax, "GeV");
  RooRealVar mass("mass", "m(ee)", 0, maxMass * 10, "GeV");
  //RooRealVar mass("massDrMatched", "m(ee)", 0, maxMass * 10, "GeV");
  RooRealVar evtRegion("evtRegion", "Detector region", 0., 2.);
  RooRealVar eleEta("eleEta", "eleEta", -3., 3.);

  // get invariant mass histo from file
  //TFile input(inFile, "read");
  //input.cd();
  //TTree *tree = (TTree *)gDirectory->Get(treeName);

  // Read dataset from tree
  RooDataSet *dataAll = new RooDataSet("dataAll", "complete dataset", RooArgSet(mass, trueMass, evtRegion, eleEta), Import(*tree));
  dataAll->Print();

  //input.Close(); 

  // add column with m_reco - m_true
  RooFormulaVar massDiff("massDiff", "m_{RECO} - m_{true} (GeV)", "@0-@1", RooArgSet(mass, trueMass));
  RooRealVar *mDiff = (RooRealVar *)dataAll->addColumn(massDiff);
  mDiff->setRange(-maxMass, maxMass);
  dataAll->Print("v");

  RooDataSet *data = new RooDataSet("data", "pruned dataset", RooArgSet(*mDiff));

  // prune dataset according to the desired regions. EE(4)/BE(2)/BB(1)
  if ((regions & 1) == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0"));
  if ((regions & 2) == 2) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 1"));
  if (data->numEntries() == 0) { // return if no falid region selected
    std::cout << "No event region selected. No data is imported." << std::endl;
    return w;
  }
  data->Print("v");

  float plotRangeMin = (minMass - maxMass) / 4.;
  float plotRangeMax = (maxMass - minMass) / 4.;
  if (plotRangeMax < 15.) {
    plotRangeMin = (minMass - maxMass) / 2.;
    plotRangeMax = (maxMass - minMass) / 2.;
  }
  RooPlot *plot = mDiff->frame(Name("plot"), Title("plot"), Range(plotRangeMin, plotRangeMax), Bins(nBins));
  data->plotOn(plot, Name("dataMDiff"));

  // Build p.d.f.
  ////////////////////////////////////////////////
  //             Parameters                     //
  ////////////////////////////////////////////////

  //  Signal p.d.f. parameters
  //  Parameters for a Crystal Ball and a double sided Crystal Ball lineshape
  RooRealVar cbBias ("#Deltam_{CB}", "CB Bias", -.5, -500, 500, "GeV");
  RooRealVar cbSigma("#sigma_{CB}", "CB Width", 40., 0.2, 500., "GeV");
  RooRealVar cbCut  ("a_{CB}", "CB Cut", 1., 0.1, 200.);
  RooRealVar cbPower("n_{CB}", "CB Power", 2., 0.5, 20.);
  cbCut.setVal(cutoff_cb);
  //cbCut.setConstant(kTRUE);  // fix the cbCut_off

  RooRealVar dCBBias ("#Deltam_{DCB}", "Double CB Bias", -.5, -500, 500, "GeV");
  RooRealVar dCBSigma ("#sigma_{DCB}", "Double CB Width", 4., 0.2, 500., "GeV");
  RooRealVar dCBCutL ("al_{DCB}", "Double CB Cut left", 1., 0.1, 50.);
  RooRealVar dCBCutR ("ar_{DCB}", "Double CB Cut right", 1., 0.1, 50.);
  RooRealVar dCBPowerL ("nl_{DCB}", "Double CB Power left", 2., 0.2, 50.);
  RooRealVar dCBPowerR ("nr_{DCB}", "Double CB Power right", 2., 0.2, 50.);

  RooRealVar gBias ("#Deltam_{Gauss}", "Gaussian Bias", -.5, -50, 50, "GeV");
  RooRealVar gSigma("#sigma_{Gauss}", "Gaussian Width", 4., 0.02, 50., "GeV");

  // fraction of signal
  RooRealVar  nsig("N_{S}", "#signal events", 524, 0.1, 10000000000.);
  nsig.setVal(data->numEntries());
  
  RooArgSet plotParam(nsig, "plotParam");
  if (fitModelType == 0) plotParam.add(RooArgSet(cbBias, cbSigma, cbCut, cbPower));
  else if (fitModelType == 1) plotParam.add(RooArgSet(dCBBias, dCBSigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR));
  else plotParam.add(RooArgSet(gBias, gSigma));
  ////////////////////////////////////////////////
  //               P.D.F.s                      //
  ////////////////////////////////////////////////
  RooCBShape cball("cball", "A Crystal Ball Lineshape", *mDiff, cbBias, cbSigma, cbCut, cbPower);

#ifdef __CINT__
  gROOT->ProcessLineSync(".x RooDCBShape.cxx+") ;
#endif  
  RooDCBShape dCBall("dCBall", "A double Crystal Ball lineshape", *mDiff, dCBBias, dCBSigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR);

  RooGaussian gauss("gauss", "A Gaussian", *mDiff, gBias, gSigma);

  // set fit range
  RooArgList fitModel;
  Float_t fitMin = (minMass - maxMass) / 4.;
  Float_t fitMax = (maxMass - minMass) / 4.;
  if (fitModelType == 0) {
    fitModel.add(cball);
  }
  else if (fitModelType == 1) {
    fitModel.add(dCBall);
    fitMin = (minMass - maxMass) / 2.;
    fitMax = (maxMass - minMass) / 2.;
  }
  else 
    fitModel.add(gauss);
  if (maxMass - minMass < 40.) {
    fitMin = (minMass - maxMass) / 2.;
    fitMax = (maxMass - minMass) / 2.;
  }
 
  // Di-Electron mass model p.d.f.
  RooAddPdf coarseModel("coarseModel", "signal", fitModel, RooArgList(nsig)); 
  RooAddPdf model("model", "signal", fitModel, RooArgList(nsig)); 
 
  TStopwatch sw;
  if (fitModelType == 0) {
    // do the coarse fit
    sw.Start();
    coarseModel.fitTo(*data, NumCPU(2), Timer(kTRUE), Range(fitMin, fitMax));
    sw.Stop();
    sw.Print();

    coarseModel.plotOn(plot, Name("coarseMod"), LineStyle(7));
    std::cout << "chi2 of the coarse fit: " << plot->chiSquare("coarseMod", "dataMDiff", 5+fitModelType*2) << std::endl;;

    // do the fine fit
    std::cout << "coarse fit found mean:  " << cbBias.getValV() << std::endl;
    std::cout << "coarse fit found sigma: " << cbSigma.getValV() << std::endl;
    std::cout << "coarse fit range was: " << fitMin << " to " << fitMax << std::endl;
    fitMin = cbBias.getValV() - 4.5 * cbSigma.getValV();
    fitMax = cbBias.getValV() + 3. * cbSigma.getValV();
    std::cout << "fine fit range: " << fitMin << " to " << fitMax << std::endl;
    sw.Reset();
  }
  sw.Start();
  model.fitTo(*data, NumCPU(2), Timer(kTRUE), Save(), Range(fitMin, fitMax));
  sw.Stop();
  sw.Print();

  // make a plot    
  model.plotOn(plot, Name("fineMod"));
  model.paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(plotParam), Layout(0.111, 0.455, 0.951));
  plot->getAttText("model_paramBox")->SetTextFont(font);
  plot->drawBefore("dataMDiff", "model_paramBox");

  // garnish it with some text
  TLatex *texLumi = new TLatex(0.625, 0.865, "CMS simulation preliminary");
  texLumi->SetName("texLumi");
  texLumi->SetTextColor(kBlack);
  texLumi->SetNDC();
  texLumi->SetTextFont(font);
  texLumi->SetLineWidth(2);
  texLumi->SetTextSize(0.035);
  TLatex *texReg = (TLatex *)texLumi->Clone("texReg");  
  texReg->SetTextSize(0.05);
  texReg->SetText(0.692, 0.678, "");
  TLatex *texRange = (TLatex *)texLumi->Clone("texRange");
  texRange->SetText(0.692, 0.765, Form("m_{true} [%.0f, %.0f] GeV", trueMassMin, trueMassMax));
  TLatex *texChi2Coarse = (TLatex *)texLumi->Clone("texChi2Coarse");
  if (fitModelType == 0) {
    texChi2Coarse->SetTextColor(kBlue);
    texChi2Coarse->SetText(0.13, 0.48, Form("#chi^{2} coarse fit = %.2f", plot->chiSquare("coarseMod", "dataMDiff", 5+fitModelType*2)));
  }
  TLatex *texChi2Fine = (TLatex *)texLumi->Clone("texChi2Fine");
  texChi2Fine->SetTextColor(kBlue);
  texChi2Fine->SetText(0.13, 0.435, Form("#chi^{2} fit = %.2f", plot->chiSquare("fineMod", "dataMDiff", 5+fitModelType*2)));

  plot->addObject(texLumi);
  plot->addObject(texReg);
  plot->addObject(texRange);
  if (fitModelType == 0) plot->addObject(texChi2Coarse);
  plot->addObject(texChi2Fine);

  //std::cout << "chi2 of the fine fit: " << plot->chiSquare("fineMod", "dataMDiff", 5+fitModelType*2) << std::endl;;

  //plot->Print("v");
  //plot->Draw();

  // Make a pull distribution
  RooHist *pullHist = plot->pullHist("dataMDiff", "fineMod");
  RooPlot *pullPlot = mDiff->frame(Name("pullPlot"), Title("pullPlot"), Range(plotRangeMin, plotRangeMax), Bins(nBins));
  pullPlot->addPlotable((RooPlotable *)pullHist, "p");
  //pullPlot->Draw();

  // define workspace
  //w->import(coarseModel);
  //w->import(model);
  //w->import(*data);  
  //w->import(*pullHist);
  w->import(*plot);
  w->import(*pullPlot);
  w->defineSet("plotParam", plotParam, kTRUE);
  //w->Print("v"); 
  return w;
}

void HighMassRes::RunCryBall()
{
  // make tree from all dy trees in file
  TFile input(inFile, "read");

  for (unsigned int reg = 0; reg < 3; ++reg) {
    if (!plotReg[reg]) continue;
    stringstream sStream;

    //========================================================================
    // fit for signal MC
    for (unsigned int iZp = 0; iZp < zPrimeGenMasses.size(); ++iZp) {
      Float_t minMass = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.5;
      Float_t maxMass = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.5;
      Float_t trueMassMin = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.01; // ~3sigma around the peak
      Float_t trueMassMax = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.01; // ~3sigma around the peak

      // get invariant mass histo from file
      input.cd();
      TTree *zpTree = (TTree *)gDirectory->Get(zPrimeGenMasses[iZp].second);

      // fit
      std::cout << "FIT " << zPrimeGenMasses[iZp].second << " mass: " << zPrimeGenMasses[iZp].first << std::endl;
      RooWorkspace *mcZpWkSpc = cryBall(zpTree, reg + 1, minMass, maxMass, trueMassMin, trueMassMax, cutoff_cb, plotOpt, nBins, fitModelType);
      //mcZpWkSpc->Print("v");

      // retrieve plot from workspace
      RooPlot *plot = (RooPlot *)mcZpWkSpc->genobj("plot");
      RooPlot *pullPlot = (RooPlot *)mcZpWkSpc->genobj("pullPlot");

      // plot the data
      TCanvas* c = new TCanvas("c" + regTxt[reg] += iZp, "Unbinned Invariant Mass Fit " + regTxt[reg] += iZp, 0, 0, 800, 600);
      c->cd(); 
      if (plotPull) {
        c->Divide(1, 2);
        c->GetPad(1)->SetPad(0.01, 0.2, 0.99, 0.99);
        c->GetPad(2)->SetPad(0.01, 0.01, 0.99, 0.2);
        c->SetBorderMode(0);
        c->SetBorderSize(2);
        c->cd(1);
        gPad->SetBottomMargin(0.039);
        plot->SetTitle("");
        ((TPaveText *)plot->findObject("model_paramBox"))->SetX2(0.4);
        ((TLatex *)plot->findObject("texLumi"))->SetX(0.674);
        ((TLatex *)plot->findObject("texLumi"))->SetY(0.867);
        ((TLatex *)plot->findObject("texRange"))->SetX(0.732);
        ((TLatex *)plot->findObject("texRange"))->SetY(0.763);
      }
      // plot the fit results
      plot->getAttLine("fineMod")->SetLineColor(fitColor);
      plot->getAttText("texChi2Fine")->SetTextColor(fitColor);
      ((TLatex *)plot->findObject("texReg"))->SetText(0.732, 0.678, regTxt[reg].Data());
      plot->Draw();

      // plot pull histogram
      if (plotPull) {
        c->cd(2);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.4);
        gPad->SetBorderMode(0);
        pullPlot->GetYaxis()->SetNdivisions(102);
        pullPlot->GetYaxis()->SetTitle("pull");
        pullPlot->GetYaxis()->SetTitleSize(0.17);
        pullPlot->GetYaxis()->SetTitleOffset(0.32);
        pullPlot->GetYaxis()->SetLabelSize(0.15);
        pullPlot->GetXaxis()->SetTitleSize(0.17);
        pullPlot->GetXaxis()->SetLabelSize(0.15);
        pullPlot->Draw();
      }

      // safe in various file formats
      sStream.str("");
      sStream << plotDir << "resSigMc" << regFileNameSuffix[reg];
      sStream << iZp << fileNameExtra << "_" << lumi << "pb-1";
      TString saveFileName = sStream.str();
      if (saveFitsAsPdf) c->Print(saveFileName + ".pdf", "pdf");
      if (saveFitsAsPng) c->Print(saveFileName + ".png", "png");
      if (saveFitsAsRoot) c->Print(saveFileName + ".root", "root");

      // get the fit results from workspace
      //RooRealVar *cbBiasMC = mcZpWkSpc->var(biasName.Data());
      RooRealVar *cbSigmaMC = mcZpWkSpc->var(sigmaName.Data());
      //RooRealVar *cbCutMC = mcZpWkSpc->var(cutLName.Data());
      RooRealVar *cbPowerMC = mcZpWkSpc->var(powerLName.Data());

      // fill the histograms with fit results from fit range
      Float_t denom = (Float_t)zPrimeGenMasses[iZp].first;
      sigmaHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * cbSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].first, 2)));
      sigmaHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * cbSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].second, 2)));
      //dmHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbBiasMC->getVal());
      //dmHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbBiasMC->getError());
      //acbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbCutMC->getVal());
      //acbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbCutMC->getError());
      //ncbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbPowerMC->getVal());
      //ncbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbPowerMC->getError());
    }
    //========================================================================

    //========================================================================
    // plot the the high energy resolution
    TCanvas* c1 = new TCanvas("c1" + regTxt[reg], "High mass resolution fit " + regTxt[reg], 0, 0, 800, 600);
    c1->cd();
    TFile* outFile = new TFile("resolutionFit.root", "recreate");
    // plot the data and fit the model
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", fitRangeMin, fitRangeMax);
    fitFunc->SetParameters(1., 0.001);
    fitFunc->SetParNames("N", "S");
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineColor(fitColorRes);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0.5);
    sigmaHistosZpPsi[reg]->Fit("fitFunc");
cout << "Chi^2 / NDF: " << fitFunc->GetChisquare() << " / " << fitFunc->GetNDF() << endl;
    fitFunc->Write();
    outFile->Close();
    if (zPrimeGenMasses.size() > 0) {
      sigmaHistosZpPsi[reg]->SetLineColor(fitColor);
      sigmaHistosZpPsi[reg]->SetMarkerColor(fitColor);
      sigmaHistosZpPsi[reg]->SetMarkerStyle(22);
      sigmaHistosZpPsi[reg]->Draw("e1");
    }
    sigmaHistosZpPsi[reg]->GetXaxis()->SetTitle("m(e#mu) [GeV]");
    sigmaHistosZpPsi[reg]->GetXaxis()->SetTitleSize(0.04);
    sigmaHistosZpPsi[reg]->GetXaxis()->SetLabelSize(0.035);
    sigmaHistosZpPsi[reg]->GetYaxis()->SetTitle("#sigma_{fit} (%)");
    sigmaHistosZpPsi[reg]->GetYaxis()->SetTitleSize(0.04);
    sigmaHistosZpPsi[reg]->GetYaxis()->SetLabelSize(0.035);

    //TLegend *legend = new TLegend(0.741, 0.758, 0.931, 0.916);
    //legend->SetTextSize(0.04);
    //legend->SetFillStyle(0);
    //legend->AddEntry(sigmaHistosZpPsi[reg], "Z'_{LFV} #rightarrow e#mu", "lep");
    //TLegendEntry *legEntry = legend->AddEntry("fitFunc", "fit", "l");
    //legEntry->SetLineColor(fitColorRes);
    //legEntry->SetLineWidth(2);
    //legend->Draw("sames");

    TLatex *tex = new TLatex(0.25, 0.85, "CMS simulation, 8 TeV");
    tex->SetNDC();
    tex->SetTextFont(font);
    tex->SetLineWidth(2);
    tex->SetTextSize(0.04);
    tex->Draw();
    tex->DrawLatex(0.25, 0.77, regTxt[reg].Data());
    tex->DrawLatex(0.65, 0.35, "P(m|p_{0},p_{1}) = p_{0} + p_{1}m");
    tex->DrawLatex(0.65, 0.22, Form("#splitline{p_{0} = %.3f #pm %.3f}{p_{1} = %.3e #pm %.1e}", fabs(fitFunc->GetParameter(0)), fitFunc->GetParError(0), fabs(fitFunc->GetParameter(1)), fitFunc->GetParError(1)));

    // safe in various file formats
    sStream.str("");
    sStream << plotDir << "highMassRes" << regFileNameSuffix[reg];
    sStream << fileNameExtra << "_" << lumi << "pb-1";
    TString saveFileName = sStream.str();
    if (saveResAsPdf) c1->Print(saveFileName + ".pdf", "pdf");
    if (saveResAsPng) c1->Print(saveFileName + ".png", "png");
    if (saveResAsRoot) c1->Print(saveFileName + ".root", "root");
    //========================================================================

    ////========================================================================
    //// plot the the high energy bias
    //TCanvas* c2 = new TCanvas("c2" + regTxt[reg], "High mass bias " + regTxt[reg], 0, 0, 800, 600);
    //c2->cd();
    //// plot the data and fit the model
    //gStyle->SetErrorX(0.5);
    //dmHistos[reg]->SetLineWidth(1);
    //dmHistos[reg]->SetLineColor(fitColor);
    //dmHistos[reg]->SetMarkerColor(fitColor);
    //dmHistos[reg]->SetMarkerStyle(20);
    //dmHistos[reg]->Draw("e1");
    //if (zPrimeGenMasses.size() > 0) {
    //  dmHistosZpPsi[reg]->SetLineColor(fitColor);
    //  dmHistosZpSsm[reg]->SetLineColor(fitColor);
    //  dmHistosZpPsi[reg]->SetMarkerColor(fitColor);
    //  dmHistosZpSsm[reg]->SetMarkerColor(fitColor);
    //  dmHistosZpPsi[reg]->SetMarkerStyle(22);
    //  dmHistosZpSsm[reg]->SetMarkerStyle(23);
    //  dmHistosZpPsi[reg]->Draw("e1sames");
    //  dmHistosZpSsm[reg]->Draw("e1sames");
    //}

    //dmHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    //dmHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    //dmHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    //dmHistos[reg]->GetYaxis()->SetTitle("#Deltam_{CB} (GeV)");
    //dmHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    //dmHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    //TLegend *legend2 = new TLegend(0.741, 0.758, 0.931, 0.916);
    //legend2->SetTextSize(0.04);
    //legend2->SetFillStyle(0);
    //legend2->AddEntry(dmHistos[reg], "DY #rightarrow ee", "lep");
    //if (zPrimeGenMasses.size() > 0) legend2->AddEntry(dmHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    //if (zPrimeGenMasses.size() > 7) legend2->AddEntry(dmHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    //legend2->Draw("sames");

    //TLatex *tex2 = new TLatex(0.25, 0.85, "#splitline{CMS preliminary}{#sqrt{s} = 8 TeV    #int L dt = 19.6 fb^{-1}}");
    //tex2->SetNDC();
    //tex2->SetTextFont(font);
    //tex2->SetLineWidth(2);
    //tex2->SetTextSize(0.04);
    //tex2->Draw();
    //tex2->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    ////========================================================================

    ////========================================================================
    //// plot the the high energy Crystall Ball cut off parameter
    //TCanvas* c3 = new TCanvas("c3" + regTxt[reg], "High mass CB cut off " + regTxt[reg], 0, 0, 800, 600);
    //c3->cd();
    //// plot the data and fit the model
    //gStyle->SetErrorX(0.5);
    //acbHistos[reg]->SetLineWidth(1);
    //acbHistos[reg]->SetLineColor(fitColor);
    //acbHistos[reg]->SetMarkerColor(fitColor);
    //acbHistos[reg]->SetMarkerStyle(20);
    //acbHistos[reg]->Draw("e1");
    //if (zPrimeGenMasses.size() > 0) {
    //  acbHistosZpPsi[reg]->SetLineColor(fitColor);
    //  acbHistosZpSsm[reg]->SetLineColor(fitColor);
    //  acbHistosZpPsi[reg]->SetMarkerColor(fitColor);
    //  acbHistosZpSsm[reg]->SetMarkerColor(fitColor);
    //  acbHistosZpPsi[reg]->SetMarkerStyle(22);
    //  acbHistosZpSsm[reg]->SetMarkerStyle(23);
    //  acbHistosZpPsi[reg]->Draw("e1sames");
    //  acbHistosZpSsm[reg]->Draw("e1sames");
    //}

    //acbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    //acbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    //acbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    //acbHistos[reg]->GetYaxis()->SetTitle("a_{CB}");
    //acbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    //acbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    //TLegend *legend3 = new TLegend(0.741, 0.758, 0.931, 0.916);
    //legend3->SetTextSize(0.04);
    //legend3->SetFillStyle(0);
    //legend3->AddEntry(acbHistos[reg], "DY #rightarrow ee", "lep");
    //if (zPrimeGenMasses.size() > 0) legend3->AddEntry(acbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    //if (zPrimeGenMasses.size() > 7) legend3->AddEntry(acbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    //legend3->Draw("sames");

    //TLatex *tex3 = new TLatex(0.25, 0.85, "#splitline{CMS preliminary}{#sqrt{s} = 8 TeV    #int L dt = 19.6 fb^{-1}}");
    //tex3->SetNDC();
    //tex3->SetTextFont(font);
    //tex3->SetLineWidth(2);
    //tex3->SetTextSize(0.04);
    //tex3->Draw();
    //tex3->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    ////========================================================================

    ////========================================================================
    //// plot the the high energy Crystall Ball power parameter
    //TCanvas* c4 = new TCanvas("c4" + regTxt[reg], "High mass CB power " + regTxt[reg], 0, 0, 800, 600);
    //c4->cd();
    //// plot the data and fit the model
    //gStyle->SetErrorX(0.5);
    //ncbHistos[reg]->SetLineWidth(1);
    //ncbHistos[reg]->SetLineColor(fitColor);
    //ncbHistos[reg]->SetMarkerColor(fitColor);
    //ncbHistos[reg]->SetMarkerStyle(20);
    //ncbHistos[reg]->Draw("e1");
    //if (zPrimeGenMasses.size() > 0) {
    //  ncbHistosZpPsi[reg]->SetLineColor(fitColor);
    //  ncbHistosZpSsm[reg]->SetLineColor(fitColor);
    //  ncbHistosZpPsi[reg]->SetMarkerColor(fitColor);
    //  ncbHistosZpSsm[reg]->SetMarkerColor(fitColor);
    //  ncbHistosZpPsi[reg]->SetMarkerStyle(22);
    //  ncbHistosZpSsm[reg]->SetMarkerStyle(23);
    //  ncbHistosZpPsi[reg]->Draw("e1sames");
    //  ncbHistosZpSsm[reg]->Draw("e1sames");
    //}

    //ncbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    //ncbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    //ncbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    //ncbHistos[reg]->GetYaxis()->SetTitle("n_{CB}");
    //ncbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    //ncbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    //TLegend *legend4 = new TLegend(0.741, 0.758, 0.931, 0.916);
    //legend4->SetTextSize(0.04);
    //legend4->SetFillStyle(0);
    //legend4->AddEntry(ncbHistos[reg], "DY #rightarrow ee", "lep");
    //if (zPrimeGenMasses.size() > 0) legend4->AddEntry(ncbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    //if (zPrimeGenMasses.size() > 7) legend4->AddEntry(ncbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    //legend4->Draw("sames");

    //TLatex *tex4 = new TLatex(0.25, 0.85, "#splitline{CMS preliminary}{#sqrt{s} = 8 TeV    #int L dt = 19.6 fb^{-1}}");
    //tex4->SetNDC();
    //tex4->SetTextFont(font);
    //tex4->SetLineWidth(2);
    //tex4->SetTextSize(0.04);
    //tex4->Draw();
    //tex4->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    ////========================================================================

  }
}

