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
  RooRealVar ele1Eta("ele1Eta", "ele1Eta", -3., 3.);
  RooRealVar ele2Eta("ele2Eta", "ele2Eta", -3., 3.);

  // get invariant mass histo from file
  //TFile input(inFile, "read");
  //input.cd();
  //TTree *tree = (TTree *)gDirectory->Get(treeName);

  // Read dataset from tree
  RooDataSet *dataAll = new RooDataSet("dataAll", "complete dataset", RooArgSet(mass, trueMass, evtRegion, ele1Eta, ele2Eta), Import(*tree));
  dataAll->Print();

  //input.Close(); 

  // add column with m_reco - m_true
  RooFormulaVar massDiff("massDiff", "m_{RECO} - m_{true} (GeV)", "@0-@1", RooArgSet(mass, trueMass));
  RooRealVar *mDiff = (RooRealVar *)dataAll->addColumn(massDiff);
  mDiff->setRange(-maxMass, maxMass);
  dataAll->Print("v");

  RooDataSet *data = new RooDataSet("data", "pruned dataset", RooArgSet(*mDiff));

  // prune dataset according to the desired regions. EE(4)/BE(2)/BB(1)
  if ((regions & 1) == 1) {
    if (eleEbReg == 1) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0 && abs(ele1Eta) < 0.7 && abs(ele2Eta) < 0.7")); // both electrons |eta| < 0.7
    else if (eleEbReg == 2) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0 && abs(ele1Eta) >= 0.7 && abs(ele2Eta) >= 0.7")); // both electrons |eta| >= 0.7
    else if (eleEbReg == 3) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0 && (abs(ele1Eta) < 0.7 || abs(ele2Eta) < 0.7)")); // one electron |eta| < 0.7
    else if (eleEbReg == 4) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0 && (abs(ele1Eta) >= 0.7 || abs(ele2Eta) >= 0.7)")); // one electron |eta| >= 0.7
    else if (eleEbReg == 5) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0 && ((abs(ele1Eta) >= 0.7 && abs(ele2Eta) < 0.7) || (abs(ele2Eta) >= 0.7 && abs(ele1Eta) < 0.7))")); // one electron |eta| < 0.7 and one |eta| >= 0.7
    else data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 0"));
  }
  if ((regions & 2) == 2) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 1"));
  if ((regions & 4) == 4) data->append(*(RooDataSet *)dataAll->reduce(RooArgSet(*mDiff), "evtRegion == 2"));
  if (data->numEntries() == 0) { // return if no falid region selected
    std::cout << "No event region selected. No data is imported." << std::endl;
    return w;
  }
  data->Print("v");

  float plotRangeMin = (minMass - maxMass) / 4.;
  float plotRangeMax = (maxMass - minMass) / 4.;
  if (regions == 4) {
    plotRangeMin = (minMass - maxMass) / 2.; // for EE-EE
    plotRangeMax = (maxMass - minMass) / 2.; // for EE-EE
  }
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
  RooRealVar cbBias ("#Deltam_{CB}", "CB Bias", -.5, -50, 50, "GeV");
  RooRealVar cbSigma("#sigma_{CB}", "CB Width", 4., 0.02, 50., "GeV");
  RooRealVar cbCut  ("#alpha_{CB}", "CB Cut", 1., 0.1, 20.);
  RooRealVar cbPower("n_{CB}", "CB Power", 2., 0.5, 20.);
  cbCut.setVal(cutoff_cb);
  //cbCut.setConstant(kTRUE);  // fix the cbCut_off

  RooRealVar dCBBias ("#Deltam_{DCB}", "Double CB Bias", -0.5, -50, 50, "GeV");
  RooRealVar dCBSigma ("#sigma_{DCB}", "Double CB Width", 4., 0.02, 50., "GeV");
  RooRealVar dCBCutL ("#alpha_{L}^{DCB}", "Double CB Cut left", 1., 0.1, 50.);
  RooRealVar dCBCutR ("#alpha_{R}^{DCB}", "Double CB Cut right", 1., 0.1, 50.);
  RooRealVar dCBPowerL ("n_{L}^{DCB}", "Double CB Power left", 2., 0.2, 50.);
  RooRealVar dCBPowerR ("n_{R}^{DCB}", "Double CB Power right", 2., 0.2, 50.);

  // avoid failing fit with different starting values
  if (eleEbReg == 1 && maxMass == 2500) {
    dCBBias.setVal(-15.);
    dCBSigma.setVal(15.);
  }

  // fraction of signal
  RooRealVar  nsig("N_{S}", "#signal events", 524, 0.1, 10000000000.);
  nsig.setVal(data->numEntries());
  
  RooArgSet plotParam(nsig, "plotParam");
  if (fitModelType == 0) plotParam.add(RooArgSet(cbBias, cbSigma, cbCut, cbPower));
  else plotParam.add(RooArgSet(dCBBias, dCBSigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR));
  ////////////////////////////////////////////////
  //               P.D.F.s                      //
  ////////////////////////////////////////////////
  RooCBShape cball("cball", "A Crystal Ball Lineshape", *mDiff, cbBias, cbSigma, cbCut, cbPower);

#ifdef __CINT__
  gROOT->ProcessLineSync(".x RooDCBShape.cxx+") ;
#endif  
  RooDCBShape dCBall("dCBall", "A double Crystal Ball lineshape", *mDiff, dCBBias, dCBSigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR);

  // set fit range
  RooArgList fitModel;
  Float_t fitMin = (minMass - maxMass) / 4.;
  Float_t fitMax = (maxMass - minMass) / 4.;
  if (regions == 4) {
    fitMin = (minMass - maxMass) / 2.; // for EE-EE
    fitMax = (maxMass - minMass) / 2.; // for EE-EE
  }
  if (fitModelType == 0) {
    fitModel.add(cball);
  }
  else {
    fitModel.add(dCBall);
    fitMin = (minMass - maxMass) / 2.;
    fitMax = (maxMass - minMass) / 2.;
    if (regions == 4) {
      fitMin = (minMass - maxMass); // for EE-EE
      fitMax = (maxMass - minMass); // for EE-EE
    }
  }
  if (maxMass - minMass < 40.) {
    fitMin = (minMass - maxMass) / 2.;
    fitMax = (maxMass - minMass) / 2.;
    if (regions == 4) {
      fitMin = (minMass - maxMass); // for EE-EE
      fitMax = (maxMass - minMass); // for EE-EE
    }
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
  TLatex *texLumi = new TLatex(0.625, 0.865, "CMS Simulation");
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
  texRange->SetText(0.692, 0.765, Form("m_{true} [%.0f, %.0f] GeV", minMass, maxMass));
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
  //input.cd();
  THashTable* dyTreeTable = new THashTable();
  std::vector<TTree*> inTrees(1, (TTree*)gDirectory->Get(dyTreeNames[0]));
  inTrees.reserve(dyTreeNames.size());
  for (unsigned int i = 1; i < dyTreeNames.size(); ++i) {
    inTrees.push_back((TTree*)gDirectory->Get(dyTreeNames[i]));
    dyTreeTable->Add(inTrees.back());
    cout << "Entries in " << dyTreeNames[i] << ": " << inTrees.back()->GetEntries() << endl;
  }
  TFile* f = new TFile("intermediate.root", "recreate");
  TTree* allDyTree = inTrees.front()->CloneTree();
  allDyTree->Merge(dyTreeTable);
  cout << "Entries in allDyTree: " << allDyTree->GetEntries() << ", Number of branches: " << allDyTree->GetNbranches() << endl;


  for (unsigned int reg = 0; reg < 7; ++reg) {
    if (!plotReg[reg]) continue;
    stringstream sStream;
    TString saveFileName = "";
    //========================================================================
    // fit for DY samples
    for (unsigned int iRange = 0; iRange < dyFitRanges.size() - 1; ++iRange) {
      Float_t minMass = dyFitRanges[iRange];
      Float_t maxMass = dyFitRanges[iRange + 1];
      //float meanMass = 0.5 * (dyFitRanges[iRange] + dyFitRanges[iRange + 1]);
      Float_t trueMassMin = dyFitRanges[iRange];
      Float_t trueMassMax = dyFitRanges[iRange + 1];
      //Float_t trueMassMin = meanMass - meanMass * 0.02;
      //Float_t trueMassMax = meanMass + meanMass * 0.02;

      // fit
      std::cout << "mass: " << dyFitRanges[iRange] << std::endl;
      RooWorkspace *mcWkSpc = cryBall(allDyTree, reg + 1, minMass, maxMass, trueMassMin, trueMassMax, cutoff_cb, plotOpt, nBins, fitModelType);
      //mcWkSpc->Print("v");

      // retrieve plot from workspace
      RooPlot *plot = (RooPlot *)mcWkSpc->genobj("plot");
      RooPlot *pullPlot = (RooPlot *)mcWkSpc->genobj("pullPlot");
      //plot->Print("v");
      //pullPlot->Print("v");

      // save workspace to file
      mcWkSpc->SetNameTitle("wDY"+regTxt[reg]+=iRange, "workspaceDY"+regTxt[reg]+=iRange);
      sStream.str("");
      sStream << plotDir << "workSpaces" << regFileNameSuffix[reg];
      sStream << fileNameExtra << "_" << lumi << "pb-1";
      saveFileName = sStream.str();
      if (iRange == 0) mcWkSpc->writeToFile(saveFileName + ".root", kTRUE);
      else mcWkSpc->writeToFile(saveFileName + ".root", kFALSE);

      // plot the data
      TCanvas* c = new TCanvas("c" + regTxt[reg] += iRange, "Unbinned Invariant Mass Fit " + regTxt[reg] += iRange, 0, 0, 800, 600);
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
        ((TLatex *)plot->findObject("texLumi"))->SetX(0.732);
        ((TLatex *)plot->findObject("texLumi"))->SetY(0.867);
        ((TLatex *)plot->findObject("texRange"))->SetX(0.732);
        ((TLatex *)plot->findObject("texRange"))->SetY(0.763);
      }
      // plot the fit results
      plot->getAttLine("fineMod")->SetLineColor(fitColorDy);
      plot->getAttText("texChi2Fine")->SetTextColor(fitColorDy);
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
      sStream << plotDir << "resDyMc" << regFileNameSuffix[reg];
      sStream << iRange << fileNameExtra << "_" << lumi << "pb-1";
      saveFileName = sStream.str();
      if (saveFitsAsPdf) c->Print(saveFileName + ".pdf", "pdf");
      if (saveFitsAsPng) c->Print(saveFileName + ".png", "png");
      if (saveFitsAsRoot) c->Print(saveFileName + ".root", "root");

      // get the fit results from workspace
      RooRealVar *cbBiasMC = mcWkSpc->var(biasName.Data());
      RooRealVar *cbSigmaMC = mcWkSpc->var(sigmaName.Data());
      RooRealVar *cbCutMC = mcWkSpc->var(cutLName.Data());
      RooRealVar *cbPowerMC = mcWkSpc->var(powerLName.Data());
      RooRealVar *dcbCutRMC = mcWkSpc->var(cutRName.Data());
      RooRealVar *dcbPowerRMC = mcWkSpc->var(powerRName.Data());

      // fill the histograms with fit results from fit range
      Float_t denom = (dyFitRanges[iRange] + dyFitRanges[iRange + 1]) / 2;
      sigmaHistos[reg]->SetBinContent(iRange + 1, sqrt(pow(100 * cbSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg].first, 2)));
      sigmaHistos[reg]->SetBinError(iRange + 1, sqrt(pow(100 * cbSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg].second, 2)));
cout << sqrt(pow(100 * cbSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg].first, 2)) << " +/- ";
cout << sqrt(pow(100 * cbSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg].second, 2)) << endl;
      dmHistos[reg]->SetBinContent(iRange + 1, cbBiasMC->getVal());
      dmHistos[reg]->SetBinError(iRange + 1, cbBiasMC->getError());
      acbHistos[reg]->SetBinContent(iRange + 1, cbCutMC->getVal());
      acbHistos[reg]->SetBinError(iRange + 1, cbCutMC->getError());
      ncbHistos[reg]->SetBinContent(iRange + 1, cbPowerMC->getVal());
      ncbHistos[reg]->SetBinError(iRange + 1, cbPowerMC->getError());
      if (fitModelType > 0) {
        ardcbHistos[reg]->SetBinContent(iRange + 1, dcbCutRMC->getVal());
        ardcbHistos[reg]->SetBinError(iRange + 1, dcbCutRMC->getError());
        nrdcbHistos[reg]->SetBinContent(iRange + 1, dcbPowerRMC->getVal());
        nrdcbHistos[reg]->SetBinError(iRange + 1, dcbPowerRMC->getError());
      }
    }
    //========================================================================

    //========================================================================
    // fit for signal MC
    for (unsigned int iZp = 0; iZp < zPrimeGenMasses.size(); ++iZp) {
      Float_t minMass = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.075;
      Float_t maxMass = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.075;
      Float_t trueMassMin = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.01; // ~3sigma around the peak
      Float_t trueMassMax = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.01; // ~3sigma around the peak
      if (iZp >= ssmStart) {
        trueMassMin = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.05; // ~3sigma around the peak
        trueMassMax = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.05; // ~3sigma around the peak
      }

      // get invariant mass histo from file
      //TFile input(inFile, "read");
      input.cd();
      TTree *zpTree = (TTree *)gDirectory->Get(zPrimeGenMasses[iZp].second);

      // fit
      std::cout << "FIT " << zPrimeGenMasses[iZp].second << " mass: " << zPrimeGenMasses[iZp].first << std::endl;
      RooWorkspace *mcZpWkSpc = cryBall(zpTree, reg + 1, minMass, maxMass, trueMassMin, trueMassMax, cutoff_cb, plotOpt, nBins, fitModelType);
      //mcZpWkSpc->Print("v");

      // retrieve plot from workspace
      RooPlot *plot = (RooPlot *)mcZpWkSpc->genobj("plot");
      RooPlot *pullPlot = (RooPlot *)mcZpWkSpc->genobj("pullPlot");

      // save workspace to file
      mcZpWkSpc->SetNameTitle("wZp"+regTxt[reg]+=iZp, "workspaceZp"+regTxt[reg]+=iZp);
      sStream.str("");
      sStream << plotDir << "workSpaces" << regFileNameSuffix[reg];
      sStream << fileNameExtra << "_" << lumi << "pb-1";
      saveFileName = sStream.str();
      mcZpWkSpc->writeToFile(saveFileName + ".root", kFALSE);

      // plot the data
      TCanvas* c = new TCanvas("c" + regTxt[reg] += (dyFitRanges.size() - 1 + iZp), "Unbinned Invariant Mass Fit " + regTxt[reg] += (dyFitRanges.size() - 1 + iZp), 0, 0, 800, 600);
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
        ((TLatex *)plot->findObject("texLumi"))->SetX(0.732);
        ((TLatex *)plot->findObject("texLumi"))->SetY(0.867);
        ((TLatex *)plot->findObject("texRange"))->SetX(0.732);
        ((TLatex *)plot->findObject("texRange"))->SetY(0.763);
      }
      // plot the fit results
      if (iZp < ssmStart) {
        plot->getAttLine("fineMod")->SetLineColor(fitColorZpPsi);
        plot->getAttText("texChi2Fine")->SetTextColor(fitColorZpPsi);
      } else {
        plot->getAttLine("fineMod")->SetLineColor(fitColorZpSsm);
        plot->getAttText("texChi2Fine")->SetTextColor(fitColorZpSsm);
      }
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
      saveFileName = sStream.str();
      if (saveFitsAsPdf) c->Print(saveFileName + ".pdf", "pdf");
      if (saveFitsAsPng) c->Print(saveFileName + ".png", "png");
      if (saveFitsAsRoot) c->Print(saveFileName + ".root", "root");

      // get the fit results from workspace
      RooRealVar *cbBiasSigMC = mcZpWkSpc->var(biasName.Data());
      RooRealVar *cbSigmaSigMC = mcZpWkSpc->var(sigmaName.Data());
      RooRealVar *cbCutSigMC = mcZpWkSpc->var(cutLName.Data());
      RooRealVar *cbPowerSigMC = mcZpWkSpc->var(powerLName.Data());
      RooRealVar *dcbCutRSigMC = mcZpWkSpc->var(cutRName.Data());
      RooRealVar *dcbPowerRSigMC = mcZpWkSpc->var(powerRName.Data());

      // fill the histograms with fit results from fit range
      Float_t denom = (Float_t)zPrimeGenMasses[iZp].first;
      if (iZp < ssmStart) {
        sigmaHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * cbSigmaSigMC->getVal() / denom, 2) + pow(sigmaExtras[reg].first, 2)));
        sigmaHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * cbSigmaSigMC->getError() / denom, 2) + pow(sigmaExtras[reg].second, 2)));
        dmHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbBiasSigMC->getVal());
        dmHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbBiasSigMC->getError());
        acbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbCutSigMC->getVal());
        acbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbCutSigMC->getError());
        ncbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbPowerSigMC->getVal());
        ncbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbPowerSigMC->getError());
        if (fitModelType > 0) {
          ardcbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, dcbCutRSigMC->getVal());
          ardcbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), dcbCutRSigMC->getError());
          nrdcbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, dcbPowerRSigMC->getVal());
          nrdcbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), dcbPowerRSigMC->getError());
        }
      } else {
        sigmaHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * cbSigmaSigMC->getVal() / denom, 2) + pow(sigmaExtras[reg].first, 2)));
        sigmaHistosZpSsm[reg]->SetBinError(sigmaHistosZpSsm[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * cbSigmaSigMC->getError() / denom, 2) + pow(sigmaExtras[reg].second, 2)));
        dmHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, cbBiasSigMC->getVal());
        dmHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbBiasSigMC->getError());
        acbHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, cbCutSigMC->getVal());
        acbHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbCutSigMC->getError());
        ncbHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, cbPowerSigMC->getVal());
        ncbHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbPowerSigMC->getError());
        if (fitModelType > 0) {
          ardcbHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, dcbCutRSigMC->getVal());
          ardcbHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), dcbCutRSigMC->getError());
          nrdcbHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, dcbPowerRSigMC->getVal());
          nrdcbHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), dcbPowerRSigMC->getError());
        }
      }
    }
    //========================================================================

    //========================================================================
    // plot the the high energy resolution
    TCanvas* c1 = new TCanvas("c1" + regTxt[reg], "High mass resolution fit " + regTxt[reg], 0, 0, 800, 600);
    c1->cd();
    // plot the data and fit the model
    TF1 *fitFunc = new TF1("fitFunc", "sqrt([0]^2/x^2 + [1]^2/x + [2]^2 + [3]^2*x)", dyFitRanges.front(), dyFitRanges.back());
    fitFunc->SetParameters(165., 10., 1.5, 0.01);
    if (!useRootTermForFit) fitFunc->FixParameter(3, 0.);
    fitFunc->SetParNames("N", "S", "C", "R");
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineColor(fitColorRes);
    gStyle->SetOptFit(0);
    sigmaHistos[reg]->Fit("fitFunc");
cout << "Chi^2 / NDF: " << fitFunc->GetChisquare() << " / " << fitFunc->GetNDF() << endl;
    gStyle->SetErrorX(0.5);
    sigmaHistos[reg]->SetLineWidth(1);
    sigmaHistos[reg]->SetLineColor(fitColorDy);
    sigmaHistos[reg]->SetMarkerColor(fitColorDy);
    sigmaHistos[reg]->SetMarkerStyle(20);
    sigmaHistos[reg]->Draw("e1sames");
    if (zPrimeGenMasses.size() > 0) {
      sigmaHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      sigmaHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      sigmaHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      sigmaHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      sigmaHistosZpPsi[reg]->SetMarkerStyle(22);
      sigmaHistosZpSsm[reg]->SetMarkerStyle(23);
      sigmaHistosZpPsi[reg]->Draw("e1sames");
      sigmaHistosZpSsm[reg]->Draw("e1sames");
    }

    sigmaHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    sigmaHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    sigmaHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    if (mcOnly) sigmaHistos[reg]->GetYaxis()->SetTitle(sigmaName+" (%)");
    else sigmaHistos[reg]->GetYaxis()->SetTitle("#sqrt{#sigma_{fit}^{2}+#sigma_{extra}^{2}} (%)");
    sigmaHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    sigmaHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->AddEntry(sigmaHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend->AddEntry(sigmaHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend->AddEntry(sigmaHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    TLegendEntry *legEntry = legend->AddEntry("fitFunc", "fit to DY", "l");
    legEntry->SetLineColor(fitColorRes);
    legEntry->SetLineWidth(2);
    legend->Draw("sames");

    TLatex *tex;
    if (mcOnly) tex = new TLatex(0.77, 0.965, "CMS Simulation");
    else tex = new TLatex(0.49, 0.965, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex->SetNDC();
    tex->SetTextFont(font);
    tex->SetLineWidth(2);
    tex->SetTextSize(0.04);
    tex->Draw();
    tex->DrawLatex(0.1, 0.965, regTxt[reg].Data());
    //tex->DrawLatex(0.3, 0.65, Form("#sigma_{extra} = %.2f%% #pm %.2f%%",sigmaExtras[reg].first,sigmaExtras[reg].second));
    if (eleEbReg == 1) tex->DrawLatex(0.3, 0.65, "two e with |#eta| < 0.7");
    else if (eleEbReg == 2) tex->DrawLatex(0.3, 0.65, "two e with |#eta| #geq 0.7");
    else if (eleEbReg == 3 || eleEbReg == 5) tex->DrawLatex(0.3, 0.65, "one e with |#eta| < 0.7");
    else if (eleEbReg == 4) tex->DrawLatex(0.3, 0.65, "one e with |#eta| #geq 0.7");
    if (eleEbReg == 5) tex->DrawLatex(0.3, 0.6, "one e with |#eta| #geq 0.7");
    if (useRootTermForFit) {
      tex->DrawLatex(0.6, 0.61, "#sqrt{#frac{N^{2}}{m^{2}} + #frac{S^{2}}{m} + C^{2} + R^{2}m}");
      tex->DrawLatex(0.6, 0.45, Form("#splitline{#splitline{N = %.3f #pm %.3f}{S = %.3f #pm %.3f}}{#splitline{C = %.3f #pm %.3f}{R = %.3f #pm %.3f}}", fabs(fitFunc->GetParameter(0)), fitFunc->GetParError(0), fabs(fitFunc->GetParameter(1)), fitFunc->GetParError(1), fabs(fitFunc->GetParameter(2)), fitFunc->GetParError(2), fabs(fitFunc->GetParameter(3)), fitFunc->GetParError(3)));
    }
    else {
      tex->DrawLatex(0.6, 0.61, "#sqrt{#frac{N^{2}}{m^{2}} + #frac{S^{2}}{m} + C^{2}}");
      tex->DrawLatex(0.6, 0.5, Form("#splitline{N = %.3f #pm %.3f}{#splitline{S = %.3f #pm %.3f}{C = %.3f #pm %.3f}}", fabs(fitFunc->GetParameter(0)), fitFunc->GetParError(0), fabs(fitFunc->GetParameter(1)), fitFunc->GetParError(1), fabs(fitFunc->GetParameter(2)), fitFunc->GetParError(2)));
    }

    // safe in various file formats
    sStream.str("");
    sStream << plotDir << "highMassRes" << regFileNameSuffix[reg];
    sStream << fileNameExtra << "_" << lumi << "pb-1";
    saveFileName = sStream.str();
    if (saveResAsPdf) c1->Print(saveFileName + ".pdf", "pdf");
    if (saveResAsPng) c1->Print(saveFileName + ".png", "png");
    if (saveResAsRoot) c1->Print(saveFileName + ".root", "root");
    //========================================================================

    //========================================================================
    // plot the the high energy bias
    TCanvas* c2 = new TCanvas("c2" + regTxt[reg], "High mass bias " + regTxt[reg], 0, 0, 800, 600);
    c2->cd();
    // plot the data and fit the model
    gStyle->SetErrorX(0.5);
    dmHistos[reg]->SetLineWidth(1);
    dmHistos[reg]->SetLineColor(fitColorDy);
    dmHistos[reg]->SetMarkerColor(fitColorDy);
    dmHistos[reg]->SetMarkerStyle(20);
    dmHistos[reg]->Draw("e1");
    if (zPrimeGenMasses.size() > 0) {
      dmHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      dmHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      dmHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      dmHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      dmHistosZpPsi[reg]->SetMarkerStyle(22);
      dmHistosZpSsm[reg]->SetMarkerStyle(23);
      dmHistosZpPsi[reg]->Draw("e1sames");
      dmHistosZpSsm[reg]->Draw("e1sames");
    }

    dmHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    dmHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    dmHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    dmHistos[reg]->GetYaxis()->SetTitle(biasName+" (GeV)");
    dmHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    dmHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend2 = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend2->SetTextSize(0.04);
    legend2->SetFillStyle(0);
    legend2->AddEntry(dmHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend2->AddEntry(dmHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend2->AddEntry(dmHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    legend2->Draw("sames");

    TLatex *tex2;
    if (mcOnly) tex2 = new TLatex(0.77, 0.965, "CMS Simulation");
    else tex2 = new TLatex(0.49, 0.965, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex2->SetNDC();
    tex2->SetTextFont(font);
    tex2->SetLineWidth(2);
    tex2->SetTextSize(0.04);
    tex2->Draw();
    tex2->DrawLatex(0.1, 0.965, regTxt[reg].Data());
    if (eleEbReg == 1) tex2->DrawLatex(0.72, 0.71, "two e with |#eta| < 0.7");
    else if (eleEbReg == 2) tex2->DrawLatex(0.72, 0.71, "two e with |#eta| #geq 0.7");
    else if (eleEbReg == 3 || eleEbReg == 5) tex2->DrawLatex(0.72, 0.71, "one e with |#eta| < 0.7");
    else if (eleEbReg == 4) tex2->DrawLatex(0.72, 0.71, "one e with |#eta| #geq 0.7");
    if (eleEbReg == 5) tex2->DrawLatex(0.72, 0.66, "one e with |#eta| #geq 0.7");

    // safe in various file formats
    sStream.str("");
    sStream << plotDir << "highMassBias" << regFileNameSuffix[reg];
    sStream << fileNameExtra << "_" << lumi << "pb-1";
    saveFileName = sStream.str();
    if (saveResAsPdf) c2->Print(saveFileName + ".pdf", "pdf");
    if (saveResAsPng) c2->Print(saveFileName + ".png", "png");
    if (saveResAsRoot) c2->Print(saveFileName + ".root", "root");
    //========================================================================

    //========================================================================
    // plot the the high energy Crystall Ball cut off parameter
    TCanvas* c3 = new TCanvas("c3" + regTxt[reg], "High mass CB cut off " + regTxt[reg], 0, 0, 800, 600);
    c3->cd();
    // plot the data and fit the model
    gStyle->SetErrorX(0.5);
    acbHistos[reg]->SetLineWidth(1);
    acbHistos[reg]->SetLineColor(fitColorDy);
    acbHistos[reg]->SetMarkerColor(fitColorDy);
    acbHistos[reg]->SetMarkerStyle(20);
    acbHistos[reg]->Draw("e1");
    if (zPrimeGenMasses.size() > 0) {
      acbHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      acbHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      acbHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      acbHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      acbHistosZpPsi[reg]->SetMarkerStyle(22);
      acbHistosZpSsm[reg]->SetMarkerStyle(23);
      acbHistosZpPsi[reg]->Draw("e1sames");
      acbHistosZpSsm[reg]->Draw("e1sames");
    }

    acbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    acbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    acbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    acbHistos[reg]->GetYaxis()->SetTitle(cutLName.Data());
    acbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    acbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend3 = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend3->SetTextSize(0.04);
    legend3->SetFillStyle(0);
    legend3->AddEntry(acbHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend3->AddEntry(acbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend3->AddEntry(acbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    legend3->Draw("sames");

    TLatex *tex3;
    if (mcOnly) tex3 = new TLatex(0.77, 0.965, "CMS Simulation");
    else tex3 = new TLatex(0.49, 0.965, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex3->SetNDC();
    tex3->SetTextFont(font);
    tex3->SetLineWidth(2);
    tex3->SetTextSize(0.04);
    tex3->Draw();
    tex3->DrawLatex(0.1, 0.965, regTxt[reg].Data());
    if (eleEbReg == 1) tex3->DrawLatex(0.72, 0.71, "two e with |#eta| < 0.7");
    else if (eleEbReg == 2) tex3->DrawLatex(0.72, 0.71, "two e with |#eta| #geq 0.7");
    else if (eleEbReg == 3 || eleEbReg == 5) tex3->DrawLatex(0.72, 0.71, "one e with |#eta| < 0.7");
    else if (eleEbReg == 4) tex3->DrawLatex(0.72, 0.71, "one e with |#eta| #geq 0.7");
    if (eleEbReg == 5) tex3->DrawLatex(0.72, 0.66, "one e with |#eta| #geq 0.7");

    // safe in various file formats
    sStream.str("");
    sStream << plotDir << "highMassCBCutoff" << regFileNameSuffix[reg];
    sStream << fileNameExtra << "_" << lumi << "pb-1";
    saveFileName = sStream.str();
    if (saveResAsPdf) c3->Print(saveFileName + ".pdf", "pdf");
    if (saveResAsPng) c3->Print(saveFileName + ".png", "png");
    if (saveResAsRoot) c3->Print(saveFileName + ".root", "root");
    //========================================================================

    //========================================================================
    // plot the the high energy Crystall Ball power parameter
    TCanvas* c4 = new TCanvas("c4" + regTxt[reg], "High mass CB power " + regTxt[reg], 0, 0, 800, 600);
    c4->cd();
    // plot the data and fit the model
    gStyle->SetErrorX(0.5);
    ncbHistos[reg]->SetLineWidth(1);
    ncbHistos[reg]->SetLineColor(fitColorDy);
    ncbHistos[reg]->SetMarkerColor(fitColorDy);
    ncbHistos[reg]->SetMarkerStyle(20);
    ncbHistos[reg]->Draw("e1");
    if (zPrimeGenMasses.size() > 0) {
      ncbHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      ncbHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      ncbHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      ncbHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      ncbHistosZpPsi[reg]->SetMarkerStyle(22);
      ncbHistosZpSsm[reg]->SetMarkerStyle(23);
      ncbHistosZpPsi[reg]->Draw("e1sames");
      ncbHistosZpSsm[reg]->Draw("e1sames");
    }

    ncbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    ncbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    ncbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    ncbHistos[reg]->GetYaxis()->SetTitle(powerLName.Data());
    ncbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    ncbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend4 = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend4->SetTextSize(0.04);
    legend4->SetFillStyle(0);
    legend4->AddEntry(ncbHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend4->AddEntry(ncbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend4->AddEntry(ncbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    legend4->Draw("sames");

    TLatex *tex4;
    if (mcOnly) tex4 = new TLatex(0.77, 0.965, "CMS Simulation");
    else tex4 = new TLatex(0.49, 0.965, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex4->SetNDC();
    tex4->SetTextFont(font);
    tex4->SetLineWidth(2);
    tex4->SetTextSize(0.04);
    tex4->Draw();
    tex4->DrawLatex(0.1, 0.965, regTxt[reg].Data());
    if (eleEbReg == 1) tex4->DrawLatex(0.72, 0.71, "two e with |#eta| < 0.7");
    else if (eleEbReg == 2) tex4->DrawLatex(0.72, 0.71, "two e with |#eta| #geq 0.7");
    else if (eleEbReg == 3 || eleEbReg == 5) tex4->DrawLatex(0.72, 0.71, "one e with |#eta| < 0.7");
    else if (eleEbReg == 4) tex4->DrawLatex(0.72, 0.71, "one e with |#eta| #geq 0.7");
    if (eleEbReg == 5) tex4->DrawLatex(0.72, 0.66, "one e with |#eta| #geq 0.7");

    // safe in various file formats
    sStream.str("");
    sStream << plotDir << "highMassCBPower" << regFileNameSuffix[reg];
    sStream << fileNameExtra << "_" << lumi << "pb-1";
    saveFileName = sStream.str();
    if (saveResAsPdf) c4->Print(saveFileName + ".pdf", "pdf");
    if (saveResAsPng) c4->Print(saveFileName + ".png", "png");
    if (saveResAsRoot) c4->Print(saveFileName + ".root", "root");
    //========================================================================

    if (fitModelType > 0) {
      //========================================================================
      // plot the the high energy Crystall Ball cut off parameter
      TCanvas* c5 = new TCanvas("c5" + regTxt[reg], "High mass DCB right cut off " + regTxt[reg], 0, 0, 800, 600);
      c5->cd();
      // plot the data and fit the model
      gStyle->SetErrorX(0.5);
      ardcbHistos[reg]->SetLineWidth(1);
      ardcbHistos[reg]->SetLineColor(fitColorDy);
      ardcbHistos[reg]->SetMarkerColor(fitColorDy);
      ardcbHistos[reg]->SetMarkerStyle(20);
      ardcbHistos[reg]->Draw("e1");
      if (zPrimeGenMasses.size() > 0) {
        ardcbHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
        ardcbHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
        ardcbHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
        ardcbHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
        ardcbHistosZpPsi[reg]->SetMarkerStyle(22);
        ardcbHistosZpSsm[reg]->SetMarkerStyle(23);
        ardcbHistosZpPsi[reg]->Draw("e1sames");
        ardcbHistosZpSsm[reg]->Draw("e1sames");
      }

      ardcbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
      ardcbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
      ardcbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
      ardcbHistos[reg]->GetYaxis()->SetTitle(cutRName.Data());
      ardcbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
      ardcbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

      TLegend *legend5 = new TLegend(0.741, 0.758, 0.931, 0.916);
      legend5->SetTextSize(0.04);
      legend5->SetFillStyle(0);
      legend5->AddEntry(ardcbHistos[reg], "DY #rightarrow ee", "lep");
      if (zPrimeGenMasses.size() > 0) legend5->AddEntry(ardcbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
      if (zPrimeGenMasses.size() > 7) legend5->AddEntry(acbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
      legend5->Draw("sames");

      TLatex *tex5;
      if (mcOnly) tex5 = new TLatex(0.77, 0.965, "CMS Simulation");
      else tex5 = new TLatex(0.49, 0.965, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
      tex5->SetNDC();
      tex5->SetTextFont(font);
      tex5->SetLineWidth(2);
      tex5->SetTextSize(0.04);
      tex5->Draw();
      tex5->DrawLatex(0.1, 0.965, regTxt[reg].Data());
      if (eleEbReg == 1) tex5->DrawLatex(0.72, 0.71, "two e with |#eta| < 0.7");
      else if (eleEbReg == 2) tex5->DrawLatex(0.72, 0.71, "two e with |#eta| #geq 0.7");
      else if (eleEbReg == 3 || eleEbReg == 5) tex5->DrawLatex(0.72, 0.71, "one e with |#eta| < 0.7");
      else if (eleEbReg == 4) tex5->DrawLatex(0.72, 0.71, "one e with |#eta| #geq 0.7");
      if (eleEbReg == 5) tex5->DrawLatex(0.72, 0.66, "one e with |#eta| #geq 0.7");

      // safe in various file formats
      sStream.str("");
      sStream << plotDir << "highMassDCBCutoffR" << regFileNameSuffix[reg];
      sStream << fileNameExtra << "_" << lumi << "pb-1";
      saveFileName = sStream.str();
      if (saveResAsPdf) c5->Print(saveFileName + ".pdf", "pdf");
      if (saveResAsPng) c5->Print(saveFileName + ".png", "png");
      if (saveResAsRoot) c5->Print(saveFileName + ".root", "root");
      //========================================================================

      //========================================================================
      // plot the the high energy Crystall Ball power parameter
      TCanvas* c6 = new TCanvas("c6" + regTxt[reg], "High mass DCB right power " + regTxt[reg], 0, 0, 800, 600);
      c6->cd();
      // plot the data and fit the model
      gStyle->SetErrorX(0.5);
      nrdcbHistos[reg]->SetLineWidth(1);
      nrdcbHistos[reg]->SetLineColor(fitColorDy);
      nrdcbHistos[reg]->SetMarkerColor(fitColorDy);
      nrdcbHistos[reg]->SetMarkerStyle(20);
      nrdcbHistos[reg]->Draw("e1");
      if (zPrimeGenMasses.size() > 0) {
        nrdcbHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
        nrdcbHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
        nrdcbHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
        nrdcbHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
        nrdcbHistosZpPsi[reg]->SetMarkerStyle(22);
        nrdcbHistosZpSsm[reg]->SetMarkerStyle(23);
        nrdcbHistosZpPsi[reg]->Draw("e1sames");
        nrdcbHistosZpSsm[reg]->Draw("e1sames");
      }

      nrdcbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
      nrdcbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
      nrdcbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
      nrdcbHistos[reg]->GetYaxis()->SetTitle(powerRName.Data());
      nrdcbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
      nrdcbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

      TLegend *legend6 = new TLegend(0.741, 0.758, 0.931, 0.916);
      legend6->SetTextSize(0.04);
      legend6->SetFillStyle(0);
      legend6->AddEntry(nrdcbHistos[reg], "DY #rightarrow ee", "lep");
      if (zPrimeGenMasses.size() > 0) legend6->AddEntry(nrdcbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
      if (zPrimeGenMasses.size() > 7) legend6->AddEntry(nrdcbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
      legend6->Draw("sames");

      TLatex *tex6;
      if (mcOnly) tex6 = new TLatex(0.77, 0.965, "CMS Simulation");
      else tex6 = new TLatex(0.49, 0.965, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
      tex6->SetNDC();
      tex6->SetTextFont(font);
      tex6->SetLineWidth(2);
      tex6->SetTextSize(0.04);
      tex6->Draw();
      tex6->DrawLatex(0.1, 0.965, regTxt[reg].Data());
      if (eleEbReg == 1) tex6->DrawLatex(0.72, 0.71, "two e with |#eta| < 0.7");
      else if (eleEbReg == 2) tex6->DrawLatex(0.72, 0.71, "two e with |#eta| #geq 0.7");
      else if (eleEbReg == 3 || eleEbReg == 5) tex6->DrawLatex(0.72, 0.71, "one e with |#eta| < 0.7");
      else if (eleEbReg == 4) tex6->DrawLatex(0.72, 0.71, "one e with |#eta| #geq 0.7");
      if (eleEbReg == 5) tex6->DrawLatex(0.72, 0.66, "one e with |#eta| #geq 0.7");

      // safe in various file formats
      sStream.str("");
      sStream << plotDir << "highMassDCBPowerR" << regFileNameSuffix[reg];
      sStream << fileNameExtra << "_" << lumi << "pb-1";
      saveFileName = sStream.str();
      if (saveResAsPdf) c6->Print(saveFileName + ".pdf", "pdf");
      if (saveResAsPng) c6->Print(saveFileName + ".png", "png");
      if (saveResAsRoot) c6->Print(saveFileName + ".root", "root");
      //========================================================================
    }
  }
}

void HighMassRes::CompareCryBall()
{
  // make tree from all dy trees in file
  TFile input(inFile, "read");
  input.cd();
  THashTable* dyTreeTable = new THashTable();
  for (unsigned int i = 0; i < dyTreeNames.size(); ++i) {
    TTree* inTree = (TTree*)gDirectory->Get(dyTreeNames[i]);
    dyTreeTable->Add(inTree);
  }
  TTree* allDyTree = new TTree("allDyTree", "allDyTree");
  allDyTree->Merge(dyTreeTable);

  for (unsigned int reg = 0; reg < 7; ++reg) {
    if (!plotReg[reg]) continue;
    stringstream sStream;
    //========================================================================
    // fit for DY samples
    for (unsigned int iRange = 0; iRange < dyFitRanges.size() - 1; ++iRange) {
      Float_t minMass = dyFitRanges[iRange];
      Float_t maxMass = dyFitRanges[iRange + 1];

      // fit
      std::cout << "mass: " << dyFitRanges[iRange] << std::endl;
      RooWorkspace *mcWkSpc = cryBall(allDyTree, reg + 1, minMass, maxMass, minMass, maxMass, cutoff_cb, plotOpt, nBins, fitModelType);
      RooWorkspace *mcWkSpc2 = cryBall(allDyTree, reg + 1, minMass, maxMass, minMass, maxMass, cutoff_cb, plotOpt, nBins, fitModelType2);
      //mcWkSpc->Print("v");

      // retrieve plot from workspace
      RooPlot *plot = (RooPlot *)mcWkSpc->genobj("plot");
      RooPlot *pullPlot = (RooPlot *)mcWkSpc->genobj("pullPlot");
      RooPlot *plot2 = (RooPlot *)mcWkSpc2->genobj("plot");
      RooPlot *pullPlot2 = (RooPlot *)mcWkSpc2->genobj("pullPlot");
      //plot->Print("v");
      //pullPlot->Print("v");

      // plot the data
      TCanvas* c = new TCanvas("c" + regTxt[reg] += iRange, "Unbinned Invariant Mass Fit " + regTxt[reg] += iRange, 0, 0, 800, 600);
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
        ((TPaveText *)plot2->findObject("model_paramBox"))->SetX2(0.4);
        ((TLatex *)plot->findObject("texLumi"))->SetX(0.732);
        ((TLatex *)plot->findObject("texLumi"))->SetY(0.867);
        ((TLatex *)plot->findObject("texRange"))->SetX(0.732);
        ((TLatex *)plot->findObject("texRange"))->SetY(0.763);
      }
      // plot the fit results
      plot->getAttLine("fineMod")->SetLineColor(fitColorDy);
      plot->getAttLine("fineMod")->SetLineStyle(7);
      plot->getAttText("texChi2Fine")->SetTextColor(fitColorDy);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetX1(0.69);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetY1(0.335);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetX2(0.978);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetY2(0.656);
      ((TLatex *)plot->findObject("texReg"))->SetText(0.732, 0.678, regTxt[reg].Data());
      plot->Draw();
      plot2->getAttLine("fineMod")->SetLineColor(fit2ColorDy);
      plot2->getAttLine("fineMod")->SetLineStyle(7);
      plot2->getAttText("texChi2Fine")->SetTextColor(fit2ColorDy);
      ((TLatex *)plot2->findObject("texChi2Fine"))->SetY(0.39);
      plot2->getAttText("texChi2Fine")->SetTextColor(fit2ColorDy);
      plot2->remove("dataMDiff");
      plot2->remove("texLumi");
      plot2->remove("texReg");
      plot2->remove("texRange");
      plot2->Draw("sames");

      TLegend *leg = new TLegend(0.760, 0.161, 0.951, 0.320);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->AddEntry((TObject *)((RooHist *)plot->findObject("dataMDiff")), "MC data", "p");
      leg->AddEntry(((RooCurve *)plot->findObject("coarseMod")), "1st CB fit", "l");
      leg->AddEntry(((RooCurve *)plot->findObject("fineMod")), "2nd CB fit", "l");
      leg->AddEntry(((RooCurve *)plot2->findObject("fineMod")), "double CB fit", "l");
      leg->Draw("sames");

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
        pullPlot->getAttMarker()->SetMarkerColor(fitColorDy);
        pullPlot->Draw();
        pullPlot2->GetYaxis()->SetNdivisions(102);
        pullPlot2->GetYaxis()->SetTitle("pull");
        pullPlot2->GetYaxis()->SetTitleSize(0.17);
        pullPlot2->GetYaxis()->SetTitleOffset(0.32);
        pullPlot2->GetYaxis()->SetLabelSize(0.15);
        pullPlot2->GetXaxis()->SetTitleSize(0.17);
        pullPlot2->GetXaxis()->SetLabelSize(0.15);
        pullPlot2->getAttMarker()->SetMarkerColor(fit2ColorDy);
        pullPlot2->getAttMarker()->SetMarkerStyle(22);
        pullPlot2->getAttLine()->SetLineColor(fit2ColorDy);
        pullPlot2->Draw("sames");
      }

      // safe in various file formats
      sStream.str("");
      sStream << plotDir << "resDyMc" << regFileNameSuffix[reg];
      sStream << iRange << fileNameExtra << "_" << lumi << "pb-1";
      TString saveFileName = sStream.str();
      if (saveFitsAsPdf) c->Print(saveFileName + ".pdf", "pdf");
      if (saveFitsAsPng) c->Print(saveFileName + ".png", "png");
      if (saveFitsAsRoot) c->Print(saveFileName + ".root", "root");

      // get the fit results from workspace
      RooRealVar *cbBiasMC = mcWkSpc->var(biasName.Data());
      RooRealVar *cbSigmaMC = mcWkSpc->var(sigmaName.Data());
      RooRealVar *cbCutMC = mcWkSpc->var(cutLName.Data());
      RooRealVar *cbPowerMC = mcWkSpc->var(powerLName.Data());

      RooRealVar *dCBSigmaMC = mcWkSpc2->var(sigmaName.Data());

      // fill the histograms with fit results from fit range
      Float_t denom = (dyFitRanges[iRange] + dyFitRanges[iRange + 1]) / 2;
      sigmaHistos[reg]->SetBinContent(iRange + 1, sqrt(pow(100 * cbSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].first, 2)));
      sigmaHistos[reg]->SetBinError(iRange + 1, sqrt(pow(100 * cbSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].second, 2)));
      dmHistos[reg]->SetBinContent(iRange + 1, cbBiasMC->getVal());
      dmHistos[reg]->SetBinError(iRange + 1, cbBiasMC->getError());
      acbHistos[reg]->SetBinContent(iRange + 1, cbCutMC->getVal());
      acbHistos[reg]->SetBinError(iRange + 1, cbCutMC->getError());
      ncbHistos[reg]->SetBinContent(iRange + 1, cbPowerMC->getVal());
      ncbHistos[reg]->SetBinError(iRange + 1, cbPowerMC->getError());

      sigmaDCBHistos[reg]->SetBinContent(iRange + 1, sqrt(pow(100 * dCBSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType2].first, 2)));
      sigmaDCBHistos[reg]->SetBinError(iRange + 1, sqrt(pow(100 * dCBSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType2].second, 2)));
    }
    //========================================================================

    //========================================================================
    // fit for signal MC
    for (unsigned int iZp = 0; iZp < zPrimeGenMasses.size(); ++iZp) {
      Float_t minMass = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.075;
      Float_t maxMass = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.075;
      Float_t trueMassMin = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.01; // ~3sigma around the peak
      Float_t trueMassMax = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.01; // ~3sigma around the peak
      if (iZp >= ssmStart) {
        trueMassMin = zPrimeGenMasses[iZp].first - zPrimeGenMasses[iZp].first * 0.05; // ~3sigma around the peak
        trueMassMax = zPrimeGenMasses[iZp].first + zPrimeGenMasses[iZp].first * 0.05; // ~3sigma around the peak
      }

      // get invariant mass histo from file
      //TFile input(inFile, "read");
      input.cd();
      TTree *zpTree = (TTree *)gDirectory->Get(zPrimeGenMasses[iZp].second);

      // fit
      std::cout << "mass: " << zPrimeGenMasses[iZp].first << std::endl;
      RooWorkspace *mcZpWkSpc = cryBall(zpTree, reg + 1, minMass, maxMass, trueMassMin, trueMassMax, cutoff_cb, plotOpt, nBins, fitModelType);
      RooWorkspace *mcZpWkSpc2 = cryBall(zpTree, reg + 1, minMass, maxMass, trueMassMin, trueMassMax, cutoff_cb, plotOpt, nBins, fitModelType2);
      //mcZpWkSpc->Print("v");

      // retrieve plot from workspace
      RooPlot *plot = (RooPlot *)mcZpWkSpc->genobj("plot");
      RooPlot *pullPlot = (RooPlot *)mcZpWkSpc->genobj("pullPlot");
      RooPlot *plot2 = (RooPlot *)mcZpWkSpc2->genobj("plot");
      RooPlot *pullPlot2 = (RooPlot *)mcZpWkSpc2->genobj("pullPlot");

      // plot the data
      TCanvas* c = new TCanvas("c" + regTxt[reg] += (dyFitRanges.size() - 1 + iZp), "Unbinned Invariant Mass Fit " + regTxt[reg] += (dyFitRanges.size() - 1 + iZp), 0, 0, 800, 600);
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
        ((TPaveText *)plot2->findObject("model_paramBox"))->SetX2(0.4);
        ((TLatex *)plot->findObject("texLumi"))->SetX(0.732);
        ((TLatex *)plot->findObject("texLumi"))->SetY(0.867);
        ((TLatex *)plot->findObject("texRange"))->SetX(0.732);
        ((TLatex *)plot->findObject("texRange"))->SetY(0.763);
      }
      // plot the fit results
      plot->getAttLine("fineMod")->SetLineStyle(7);
      if (iZp < ssmStart) {
        plot->getAttLine("fineMod")->SetLineColor(fitColorZpPsi);
        plot->getAttText("texChi2Fine")->SetTextColor(fitColorZpPsi);
      } else {
        plot->getAttLine("fineMod")->SetLineColor(fitColorZpSsm);
        plot->getAttText("texChi2Fine")->SetTextColor(fitColorZpSsm);
      }
      ((TPaveText *)plot->findObject("model_paramBox"))->SetX1(0.69);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetY1(0.335);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetX2(0.978);
      ((TPaveText *)plot->findObject("model_paramBox"))->SetY2(0.656);
      ((TLatex *)plot->findObject("texReg"))->SetText(0.732, 0.678, regTxt[reg].Data());
      plot->Draw();
      if (iZp < ssmStart) {
        plot2->getAttLine("fineMod")->SetLineColor(fit2ColorZpPsi);
        plot2->getAttText("texChi2Fine")->SetTextColor(fit2ColorZpPsi);
        plot2->getAttText("texChi2Fine")->SetTextColor(fit2ColorZpPsi);
      }
      else {
        plot2->getAttLine("fineMod")->SetLineColor(fit2ColorZpSsm);
        plot2->getAttText("texChi2Fine")->SetTextColor(fit2ColorZpSsm);
        plot2->getAttText("texChi2Fine")->SetTextColor(fit2ColorZpSsm);
      }
      plot2->getAttLine("fineMod")->SetLineStyle(7);
      ((TLatex *)plot2->findObject("texChi2Fine"))->SetY(0.39);
      plot2->remove("dataMDiff");
      plot2->remove("texLumi");
      plot2->remove("texReg");
      plot2->remove("texRange");
      plot2->Draw("sames");

      TLegend *leg = new TLegend(0.760, 0.161, 0.951, 0.320);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->AddEntry((TObject *)((RooHist *)plot->findObject("dataMDiff")), "MC data", "p");
      leg->AddEntry(((RooCurve *)plot->findObject("coarseMod")), "1st CB fit", "l");
      leg->AddEntry(((RooCurve *)plot->findObject("fineMod")), "2nd CB fit", "l");
      leg->AddEntry(((RooCurve *)plot2->findObject("fineMod")), "double CB fit", "l");
      leg->Draw("sames");


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
        if (iZp < ssmStart) pullPlot->getAttMarker()->SetMarkerColor(fitColorZpPsi);
        else pullPlot->getAttMarker()->SetMarkerColor(fitColorZpSsm);
        pullPlot->Draw();
        pullPlot2->GetYaxis()->SetNdivisions(102);
        pullPlot2->GetYaxis()->SetTitle("pull");
        pullPlot2->GetYaxis()->SetTitleSize(0.17);
        pullPlot2->GetYaxis()->SetTitleOffset(0.32);
        pullPlot2->GetYaxis()->SetLabelSize(0.15);
        pullPlot2->GetXaxis()->SetTitleSize(0.17);
        pullPlot2->GetXaxis()->SetLabelSize(0.15);
        pullPlot2->getAttMarker()->SetMarkerStyle(22);
        if (iZp < ssmStart) {
          pullPlot2->getAttMarker()->SetMarkerColor(fit2ColorZpPsi);
          pullPlot2->getAttLine()->SetLineColor(fit2ColorZpPsi);
        } 
        else {
          pullPlot2->getAttMarker()->SetMarkerColor(fit2ColorZpSsm);
          pullPlot2->getAttLine()->SetLineColor(fit2ColorZpSsm);
        }
        pullPlot2->Draw("sames");
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
      RooRealVar *cbBiasSigMC = mcZpWkSpc->var(biasName.Data());
      RooRealVar *cbSigmaSigMC = mcZpWkSpc->var(sigmaName.Data());
      RooRealVar *cbCutSigMC = mcZpWkSpc->var(cutLName.Data());
      RooRealVar *cbPowerSigMC = mcZpWkSpc->var(powerLName.Data());

      RooRealVar *dCBSigmaMC = mcZpWkSpc2->var(sigmaName.Data());

      // fill the histograms with fit results from fit range
      Float_t denom = (Float_t)zPrimeGenMasses[iZp].first;
      if (iZp < ssmStart) {
        sigmaHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * cbSigmaSigMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].first, 2)));
        sigmaHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * cbSigmaSigMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].second, 2)));
        dmHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbBiasSigMC->getVal());
        dmHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbBiasSigMC->getError());
        acbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbCutSigMC->getVal());
        acbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbCutSigMC->getError());
        ncbHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, cbPowerSigMC->getVal());
        ncbHistosZpPsi[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbPowerSigMC->getError());

        sigmaDCBHistosZpPsi[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * dCBSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType2].first, 2)));
        sigmaDCBHistosZpPsi[reg]->SetBinError(sigmaDCBHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * dCBSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType2].second, 2)));
      } else {
        sigmaHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * cbSigmaSigMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].first, 2)));
        sigmaHistosZpSsm[reg]->SetBinError(sigmaHistosZpSsm[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * cbSigmaSigMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType].second, 2)));
        dmHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, cbBiasSigMC->getVal());
        dmHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbBiasSigMC->getError());
        acbHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, cbCutSigMC->getVal());
        acbHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbCutSigMC->getError());
        ncbHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, cbPowerSigMC->getVal());
        ncbHistosZpSsm[reg]->SetBinError(sigmaHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), cbPowerSigMC->getError());

        sigmaDCBHistosZpSsm[reg]->Fill(zPrimeGenMasses[iZp].first, sqrt(pow(100 * dCBSigmaMC->getVal() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType2].first, 2)));
        sigmaDCBHistosZpSsm[reg]->SetBinError(sigmaDCBHistosZpPsi[reg]->FindBin(zPrimeGenMasses[iZp].first), sqrt(pow(100 * dCBSigmaMC->getError() / denom, 2) + pow(sigmaExtras[reg+7*fitModelType2].second, 2)));
      }
    }
    //========================================================================

    //========================================================================
    // plot the the high energy resolution
    TCanvas* c1 = new TCanvas("c1" + regTxt[reg], "High mass resolution fit " + regTxt[reg], 0, 0, 800, 600);
    c1->cd();
    // plot the data and fit the model
    TF1 *fitFunc = new TF1("fitFunc", "sqrt([0]^2 / x + [1]^2 / x^2 + [2]^2)", dyFitRanges.front(), dyFitRanges.back());
    fitFunc->SetParameters(10., 125., 1.);
    fitFunc->SetParNames("S", "N", "C");
    fitFunc->SetLineWidth(2);
    fitFunc->SetLineColor(fitColorRes);
    TF1 *fitFunc2 = new TF1("fitFunc2", "sqrt([0]^2 / x + [1]^2 / x^2 + [2]^2)", dyFitRanges.front(), dyFitRanges.back());
    fitFunc2->SetParameters(10., 125., 1.);
    fitFunc2->SetParNames("S_{DCB}", "N_{DCB}", "C_{DCB}");
    fitFunc2->SetLineWidth(2);
    fitFunc2->SetLineColor(fitColorRes2);
    //fitFunc2->SetLineStyle(7);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0.5);
    sigmaHistos[reg]->Fit("fitFunc");
    sigmaHistos[reg]->SetLineWidth(1);
    sigmaHistos[reg]->SetLineColor(fitColorDy);
    sigmaHistos[reg]->SetMarkerColor(fitColorDy);
    sigmaHistos[reg]->SetMarkerStyle(20);
    sigmaDCBHistos[reg]->Fit("fitFunc2");
    sigmaDCBHistos[reg]->SetLineWidth(1);
    sigmaDCBHistos[reg]->SetLineColor(fit2ColorDy);
    sigmaDCBHistos[reg]->SetMarkerColor(fit2ColorDy);
    sigmaDCBHistos[reg]->SetMarkerStyle(24);
    sigmaHistos[reg]->Draw();
    if (zPrimeGenMasses.size() > 0) {
      sigmaHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      sigmaHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      sigmaHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      sigmaHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      sigmaHistosZpPsi[reg]->SetMarkerStyle(22);
      sigmaHistosZpSsm[reg]->SetMarkerStyle(23);
      sigmaHistosZpPsi[reg]->Draw("e1sames");
      sigmaHistosZpSsm[reg]->Draw("e1sames");
      sigmaDCBHistosZpPsi[reg]->SetLineColor(fit2ColorZpPsi);
      sigmaDCBHistosZpSsm[reg]->SetLineColor(fit2ColorZpSsm);
      sigmaDCBHistosZpPsi[reg]->SetMarkerColor(fit2ColorZpPsi);
      sigmaDCBHistosZpSsm[reg]->SetMarkerColor(fit2ColorZpSsm);
      sigmaDCBHistosZpPsi[reg]->SetMarkerStyle(26);
      sigmaDCBHistosZpSsm[reg]->SetMarkerStyle(32);
      sigmaDCBHistosZpPsi[reg]->Draw("e1sames");
      sigmaDCBHistosZpSsm[reg]->Draw("e1sames");
    }
    sigmaDCBHistos[reg]->Draw("e1sames");
    sigmaHistos[reg]->Draw("samesaxis");

    sigmaHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    sigmaHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    sigmaHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    if (mcOnly) sigmaHistos[reg]->GetYaxis()->SetTitle("#sigma_{fit} (%)"); // for MC only
    else sigmaHistos[reg]->GetYaxis()->SetTitle("#sqrt{#sigma_{fit}^{2}+#sigma_{extra}^{2}} (%)");
    sigmaHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    sigmaHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend = new TLegend(0.577, 0.645, 0.951, 0.935);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    legend->AddEntry(sigmaHistos[reg], "DY #rightarrow ee (CB)", "lep");
    legend->AddEntry(sigmaDCBHistos[reg], "DY #rightarrow ee (double CB)", "lep");
    if (zPrimeGenMasses.size() > 0) {
      legend->AddEntry(sigmaHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee (CB)", "lep");
      legend->AddEntry(sigmaDCBHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee (double CB)", "lep");
    }
    if (zPrimeGenMasses.size() > 7) {
      legend->AddEntry(sigmaHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee  (CB)", "lep");
      legend->AddEntry(sigmaDCBHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee (double CB)", "lep");
    }
    TLegendEntry *legEntry = legend->AddEntry("fitFunc", "CB fit to DY", "l");
    legEntry->SetLineColor(fitColorRes);
    legEntry->SetLineWidth(2);
    TLegendEntry *legEntry2 = legend->AddEntry("fitFunc2", "double CB fit to DY", "l");
    legEntry2->SetLineColor(fitColorRes2);
    legEntry2->SetLineWidth(2);
    legend->Draw("sames");

    TLatex *tex;
    if (mcOnly) tex = new TLatex(0.25, 0.85, "CMS Simulation");
    else tex = new TLatex(0.25, 0.85, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex->SetNDC();
    tex->SetTextFont(font);
    tex->SetLineWidth(2);
    tex->SetTextSize(0.04);
    tex->Draw();
    tex->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    tex->DrawLatex(0.4, 0.5, "#sqrt{#frac{S^{2}}{m} + #frac{N^{2}}{m^{2}} + C^{2}}");
    tex->DrawLatex(0.6, 0.578, Form("#splitline{C_{CB} = %.3f #pm %.3f}{#splitline{S_{CB} = %.3f #pm %.3f}{N_{CB} = %.3f #pm %.3f}}", fitFunc->GetParameter(2), fitFunc->GetParError(2), fitFunc->GetParameter(0), fitFunc->GetParError(0), fitFunc->GetParameter(1), fitFunc->GetParError(1)));
    tex->DrawLatex(0.6, 0.424, Form("#splitline{C_{DCB} = %.3f #pm %.3f}{#splitline{S_{DCB} = %.3f #pm %.3f}{N_{DCB} = %.3f #pm %.3f}}", fitFunc2->GetParameter(2), fitFunc2->GetParError(2), fitFunc2->GetParameter(0), fitFunc2->GetParError(0), fitFunc2->GetParameter(1), fitFunc2->GetParError(1)));

    // safe in various file formats
    sStream.str("");
    sStream << plotDir << "highMassRes" << regFileNameSuffix[reg];
    sStream << fileNameExtra << "_" << lumi << "pb-1";
    TString saveFileName = sStream.str();
    if (saveResAsPdf) c1->Print(saveFileName + ".pdf", "pdf");
    if (saveResAsPng) c1->Print(saveFileName + ".png", "png");
    if (saveResAsRoot) c1->Print(saveFileName + ".root", "root");
    //========================================================================

    //========================================================================
    // plot the the high energy bias
    TCanvas* c2 = new TCanvas("c2" + regTxt[reg], "High mass bias " + regTxt[reg], 0, 0, 800, 600);
    c2->cd();
    // plot the data and fit the model
    gStyle->SetErrorX(0.5);
    dmHistos[reg]->SetLineWidth(1);
    dmHistos[reg]->SetLineColor(fitColorDy);
    dmHistos[reg]->SetMarkerColor(fitColorDy);
    dmHistos[reg]->SetMarkerStyle(20);
    dmHistos[reg]->Draw("e1");
    if (zPrimeGenMasses.size() > 0) {
      dmHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      dmHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      dmHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      dmHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      dmHistosZpPsi[reg]->SetMarkerStyle(22);
      dmHistosZpSsm[reg]->SetMarkerStyle(23);
      dmHistosZpPsi[reg]->Draw("e1sames");
      dmHistosZpSsm[reg]->Draw("e1sames");
    }

    dmHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    dmHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    dmHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    dmHistos[reg]->GetYaxis()->SetTitle("#Deltam_{CB} (GeV)");
    dmHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    dmHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend2 = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend2->SetTextSize(0.04);
    legend2->SetFillStyle(0);
    legend2->AddEntry(dmHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend2->AddEntry(dmHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend2->AddEntry(dmHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    legend2->Draw("sames");

    TLatex *tex2;
    if (mcOnly) tex2 = new TLatex(0.25, 0.85, "CMS Simulation");
    else tex2 = new TLatex(0.25, 0.85, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex2->SetNDC();
    tex2->SetTextFont(font);
    tex2->SetLineWidth(2);
    tex2->SetTextSize(0.04);
    tex2->Draw();
    tex2->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    //========================================================================

    //========================================================================
    // plot the the high energy Crystall Ball cut off parameter
    TCanvas* c3 = new TCanvas("c3" + regTxt[reg], "High mass CB cut off " + regTxt[reg], 0, 0, 800, 600);
    c3->cd();
    // plot the data and fit the model
    gStyle->SetErrorX(0.5);
    acbHistos[reg]->SetLineWidth(1);
    acbHistos[reg]->SetLineColor(fitColorDy);
    acbHistos[reg]->SetMarkerColor(fitColorDy);
    acbHistos[reg]->SetMarkerStyle(20);
    acbHistos[reg]->Draw("e1");
    if (zPrimeGenMasses.size() > 0) {
      acbHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      acbHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      acbHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      acbHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      acbHistosZpPsi[reg]->SetMarkerStyle(22);
      acbHistosZpSsm[reg]->SetMarkerStyle(23);
      acbHistosZpPsi[reg]->Draw("e1sames");
      acbHistosZpSsm[reg]->Draw("e1sames");
    }

    acbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    acbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    acbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    acbHistos[reg]->GetYaxis()->SetTitle("a_{CB}");
    acbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    acbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend3 = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend3->SetTextSize(0.04);
    legend3->SetFillStyle(0);
    legend3->AddEntry(acbHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend3->AddEntry(acbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend3->AddEntry(acbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    legend3->Draw("sames");

    TLatex *tex3;
    if (mcOnly) tex3 = new TLatex(0.25, 0.85, "CMS Simulation");
    else tex3 = new TLatex(0.25, 0.85, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex3->SetNDC();
    tex3->SetTextFont(font);
    tex3->SetLineWidth(2);
    tex3->SetTextSize(0.04);
    tex3->Draw();
    tex3->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    //========================================================================

    //========================================================================
    // plot the the high energy Crystall Ball power parameter
    TCanvas* c4 = new TCanvas("c4" + regTxt[reg], "High mass CB power " + regTxt[reg], 0, 0, 800, 600);
    c4->cd();
    // plot the data and fit the model
    gStyle->SetErrorX(0.5);
    ncbHistos[reg]->SetLineWidth(1);
    ncbHistos[reg]->SetLineColor(fitColorDy);
    ncbHistos[reg]->SetMarkerColor(fitColorDy);
    ncbHistos[reg]->SetMarkerStyle(20);
    ncbHistos[reg]->Draw("e1");
    if (zPrimeGenMasses.size() > 0) {
      ncbHistosZpPsi[reg]->SetLineColor(fitColorZpPsi);
      ncbHistosZpSsm[reg]->SetLineColor(fitColorZpSsm);
      ncbHistosZpPsi[reg]->SetMarkerColor(fitColorZpPsi);
      ncbHistosZpSsm[reg]->SetMarkerColor(fitColorZpSsm);
      ncbHistosZpPsi[reg]->SetMarkerStyle(22);
      ncbHistosZpSsm[reg]->SetMarkerStyle(23);
      ncbHistosZpPsi[reg]->Draw("e1sames");
      ncbHistosZpSsm[reg]->Draw("e1sames");
    }

    ncbHistos[reg]->GetXaxis()->SetTitle("m(ee) [GeV]");
    ncbHistos[reg]->GetXaxis()->SetTitleSize(0.04);
    ncbHistos[reg]->GetXaxis()->SetLabelSize(0.035);
    ncbHistos[reg]->GetYaxis()->SetTitle("n_{CB}");
    ncbHistos[reg]->GetYaxis()->SetTitleSize(0.04);
    ncbHistos[reg]->GetYaxis()->SetLabelSize(0.035);

    TLegend *legend4 = new TLegend(0.741, 0.758, 0.931, 0.916);
    legend4->SetTextSize(0.04);
    legend4->SetFillStyle(0);
    legend4->AddEntry(ncbHistos[reg], "DY #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 0) legend4->AddEntry(ncbHistosZpPsi[reg], "Z'_{#psi} #rightarrow ee", "lep");
    if (zPrimeGenMasses.size() > 7) legend4->AddEntry(ncbHistosZpSsm[reg], "Z'_{SSM} #rightarrow ee", "lep");
    legend4->Draw("sames");

    TLatex *tex4;
    if (mcOnly) tex4 = new TLatex(0.25, 0.85, "CMS Simulation");
    else tex4 = new TLatex(0.25, 0.85, "CMS Preliminary, 8 TeV, (19.7 #pm 0.5) fb^{-1}");
    tex4->SetNDC();
    tex4->SetTextFont(font);
    tex4->SetLineWidth(2);
    tex4->SetTextSize(0.04);
    tex4->Draw();
    tex4->DrawLatex(0.3, 0.7, regTxt[reg].Data());
    //========================================================================

  }
}

