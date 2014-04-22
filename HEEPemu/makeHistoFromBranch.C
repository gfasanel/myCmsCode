#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TString.h"
#include "TParameter.h"
#include "THashList.h"

// flags: [apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
TH1F *
MakeHistoFromBranch(TFile *input, const char * treeName, const char *brName, int signs, int region, vector<const char *> &cutVariables, vector<float> &cutLows, vector<float> &cutHighs, vector<float> &mcWeigthsForCutsRanges, vector<float> &binning, unsigned int flags, bool normToBinWidth = false, float userScale = 1.)
{
  TDirectory *dir = gDirectory->CurrentDirectory();
  // prepare the histogram
  TH1F *histo = new TH1F("dummy", "dummy", 3000, 0., 3000.);
  histo->Sumw2();
  TString histoName = treeName;
  histoName.Remove(0, 7);
  histoName.Prepend(brName);
  histo->SetName(histoName);
  histo->SetTitle(histoName);
  float *bins = &binning[0];
  histo->GetXaxis()->Set(binning.size() - 1, bins);
  map<const char *, float> cutVarMap;

  if (flags & 1<<7) userScale *= ((TParameter<float> *)input->Get("lumi"))->GetVal();
  if (flags & 1<<6 && cutVariables.size() == 0) {
    THashList *mcWeights = (THashList *)input->Get("mcWeights");
    unsigned int charOffset = 8;
    if (flags & 1<<9) charOffset += 2;
    userScale *= ((TParameter<float> *)mcWeights->FindObject(treeName + charOffset))->GetVal();
  }
  float lumiScaleFactorEB = ((TParameter<float> *)input->Get("lumiScaleFactorEB"))->GetVal();
  float lumiScaleFactorEE = ((TParameter<float> *)input->Get("lumiScaleFactorEE"))->GetVal();

  //cout << "Scalefactor = " << userScale << endl;

  // get the tree
  TTree *tree;
  tree = (TTree *)dir->Get(treeName);

  // get branches
  float var;
  bool passTrg;
  bool passHeep;
  float puWeight = 1.;
  int eCharge;
  int muCharge;
  int evtRegion;
  float fakeRate = 0.;
  float trgEff = 1.;
  float trgEffSf = 1.;
  float eleEffSf = 1.;
  float muEffSf = 1.;
  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus(brName,1);
  tree->SetBranchAddress(brName, &var);
  if (flags & 1<<8) {
    tree->SetBranchStatus("passTrg",1);
    tree->SetBranchAddress("passTrg", &passTrg);
  }
  if (signs != 0) {
    tree->SetBranchStatus("eCharge",1);
    tree->SetBranchStatus("muCharge",1);
    tree->SetBranchAddress("eCharge", &eCharge);
    tree->SetBranchAddress("muCharge", &muCharge);
  }
  if (region < 2) {
    tree->SetBranchStatus("evtRegion",1);
    tree->SetBranchAddress("evtRegion", &evtRegion);
  }
  if (flags & 1) {
    tree->SetBranchStatus("puWeight",1);
    tree->SetBranchAddress("puWeight", &puWeight);
  }
  for (unsigned int i = 0; i < cutVariables.size(); ++i) {
    if (cutVariables.at(i)[0] != '\0') {
      cutVarMap.insert(pair<const char *, float>(cutVariables[i], -1.));
      tree->SetBranchStatus(cutVariables[i],1);
      tree->SetBranchAddress(cutVariables[i], &cutVarMap.at(cutVariables[i]));
    }
  }
  if (flags & 1<<5) {
    tree->SetBranchStatus("trgEff",1);
    tree->SetBranchAddress("trgEff", &trgEff);
  }
  if (flags & 1<<4) {
    tree->SetBranchStatus("trgEffSf",1);
    tree->SetBranchAddress("trgEffSf", &trgEffSf);
  }
  if (flags & 1<<2) {
    tree->SetBranchStatus("eleEffSf",1);
    tree->SetBranchAddress("eleEffSf", &eleEffSf);
  }
  if (flags & 1<<1) {
    tree->SetBranchStatus("muEffSf",1);
    tree->SetBranchAddress("muEffSf", &muEffSf);
  }
  if (flags & 1<<9) {
    tree->SetBranchStatus("passHeep",1);
    tree->SetBranchStatus("fakeRate",1);
    tree->SetBranchAddress("passHeep", &passHeep);
    tree->SetBranchAddress("fakeRate", &fakeRate);
  }

  Long64_t nEntries = (*tree).GetEntries();
  for (unsigned int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // trigger fired?
    if ((flags & 1<<8) && passTrg == false) continue;

    // select electron region
    if (evtRegion == 0 && region == 1) continue;
    if (evtRegion == 1 && region == 0) continue;

    float scaleFactor = userScale;
    if (flags & 1<<5) scaleFactor *= trgEff;
    if (flags & 1<<4) scaleFactor *= trgEffSf;
    if (flags & 1<<2) scaleFactor *= eleEffSf;
    if (flags & 1<<1) scaleFactor *= muEffSf;
    // set lumi according to detector region
    if (evtRegion == 0 && flags & 1<<3) scaleFactor *= lumiScaleFactorEB;
    if (evtRegion == 1 && flags & 1<<3) scaleFactor *= lumiScaleFactorEE;

    // PU reweight
    if (flags & 1) scaleFactor *= puWeight;

    // get only the desired charge combination. Scheme emu -3 to +3: -+, +-, OS, ALL, SS, ++, --
    if (signs < 0 && (eCharge * muCharge) > 0) continue; // OS
    if (signs > 0 && (eCharge * muCharge) < 0) continue; // SS
    if (abs(signs) == 3 && eCharge > 0) continue; // e-mu+ or e-mu-
    if (abs(signs) == 2 && eCharge < 0) continue; // e+mu- or e+mu+

    // set correct event weight if there are cuts defined
    if (flags & 1<<6 && cutVariables.size() > 0) {
      float totMcWeight = 0.;
      for (unsigned int j = 0; j < cutVariables.size(); ++j) {
        if (cutVariables.at(j)[0] == '\0') {
          totMcWeight += 1./mcWeigthsForCutsRanges[j];
        } else {
          float cutVar = cutVarMap.at(cutVariables[j]);
          if (cutVar < cutLows[j] || cutVar >= cutHighs[j]) {
            continue;
          } else {
            totMcWeight += 1./mcWeigthsForCutsRanges[j];
          }
        }
      }
      // if totMcWeight is still 0 the event will be ignored
      if (totMcWeight == 0.) continue;

      totMcWeight = 1./totMcWeight;
      scaleFactor *= totMcWeight;
    }
    
    if (flags & 1<<9) {
      if (!passHeep) scaleFactor *= fakeRate / (1 - fakeRate);
      else continue;
    }

    if (normToBinWidth) scaleFactor /= histo->GetBinWidth(histo->FindBin(var));
    histo->Fill(var, scaleFactor);
  }

  //cout << "integral: " << histo->Integral() << "       overflow: " << histo->GetBinContent(histo->GetNbinsX() + 1) << endl;
  return histo;
}

