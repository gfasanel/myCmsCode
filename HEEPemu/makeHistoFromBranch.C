#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TString.h"
#include "TParameter.h"
#include "THashList.h"

// flags: [XOR multi samples | apply top reweighting | apply fake rate | pass Trigger | lumi | MC weight | trigger Eff | trigger MC to data SF | lumi SF | ele SF | mu SF | PU reweight]
TH1F *
MakeHistoFromBranch(TFile *input, const char *treeName, const char *shapeUncName, const char *brName, int signs, int region, vector<const char *> &cutVariables, vector<float> &cutLows, vector<float> &cutHighs, vector<float> &mcWeigthsForCutsRanges, vector<float> &binning, unsigned int flags, bool normToBinWidth = false, float userScale = 1.)
{
  TString histoName = treeName;
  histoName += shapeUncName;
  histoName.Remove(0, 7);
  histoName.Prepend(brName);

  // subtraction of e-mu+ - e+mu-
  if (signs == -4) {
    TH1F *histo1 = new TH1F("histo1", "histo1", 3000, 0., 3000.);
    histo1->Sumw2();
    histo1 = MakeHistoFromBranch(input, treeName, shapeUncName, brName, -3, region, cutVariables, cutLows, cutHighs, mcWeigthsForCutsRanges, binning, flags, normToBinWidth, userScale);
    TH1F *histo2 = new TH1F("histo2", "histo2", 3000, 0., 3000.);
    histo2->Sumw2();
    histo2 = MakeHistoFromBranch(input, treeName, shapeUncName, brName, -2, region, cutVariables, cutLows, cutHighs, mcWeigthsForCutsRanges, binning, flags, normToBinWidth, userScale);
    histo1->Add(histo2, -1.);
    //histo1->Divide(histo2);
    return histo1;
  }

  unsigned int fills = 0.;
  TDirectory *dir = gDirectory->CurrentDirectory();
  // get the tree
  TTree *tree;
  TString treeNameStr = treeName;
  treeNameStr += shapeUncName;
  tree = (TTree *)dir->Get((const char*)treeNameStr);
  // if the tree does not exist because the shape uncertainty tree was 
  // not produced try to fall back to the tree without shape uncertainty
  if (tree == NULL) {
    std::cout << "Could not find tree '" << treeNameStr << "' in root file directory '" << dir->GetPath() << "'. Try to fall back to base tree '" << treeName << "'." << std::endl;
    TDirectory *baseDir = dir->GetDirectory("/");
    treeNameStr = treeName;
    tree = (TTree *)baseDir->Get((const char*)treeName);
    if (tree == NULL) std::cerr << "Error: Fallback tree '" << treeName << "' not found." << std::endl;
  }
  dir->cd();

  // prepare the histogram
  TH1F *histo = new TH1F((const char*)histoName, (const char*)histoName, 3000, 0., 3000.);
  histo->Sumw2();
  float *bins = &binning[0];
  histo->GetXaxis()->Set(binning.size() - 1, bins);
  map<const char *, float> cutVarMap;

  if (flags & 1<<7) userScale *= ((TParameter<float> *)input->Get("lumi"))->GetVal();
  if (flags & 1<<6 && cutVariables.size() == 0) {
    THashList *mcWeights = (THashList *)input->Get("mcWeights");
    unsigned int charOffset = 8;
    if (flags & 1<<9) charOffset += 2;
    userScale *= ((TParameter<float> *)mcWeights->FindObject((const char*)treeNameStr + charOffset))->GetVal();
  }
  float lumiScaleFactorEB = ((TParameter<float> *)input->Get("lumiScaleFactorEB"))->GetVal();
  float lumiScaleFactorEE = ((TParameter<float> *)input->Get("lumiScaleFactorEE"))->GetVal();

  //cout << "Scalefactor = " << userScale << endl;

  float sf_sum = 0.;
  // get branches
  float var;
  bool passTrg;
  bool passHeep;
  float puWeight = 1.;
  int eCharge = 0;
  int muCharge = 0;
  int evtRegion = 0;
  float fakeRate = 0.;
  float trgEff = 1.;
  float trgEffSf = 1.;
  float eleEffSf = 1.;
  float muEffSf = 1.;
  float topRewSf = 1.;
  float topRewSf_sum = 0.;
  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus(brName,1);
  tree->SetBranchAddress(brName, &var);
  tree->SetBranchStatus("evtRegion",1);
  tree->SetBranchAddress("evtRegion", &evtRegion);
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
  if (flags & 1) {
    tree->SetBranchStatus("puWeight",1);
    tree->SetBranchAddress("puWeight", &puWeight);
  }
  for (unsigned int i = 0; i < cutVariables.size(); ++i) {
    if (cutVariables.at(i)[0] != '\0') {
      cutVarMap.insert(pair<const char *, float>(cutVariables[i], -1.));
      tree->SetBranchStatus(cutVariables[i],1);
      if (strcmp(cutVariables[i], brName) != 0)
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
  if (flags & 1<<10) {
    tree->SetBranchStatus("topRewSf",1);
    tree->SetBranchAddress("topRewSf", &topRewSf);
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
    sf_sum += scaleFactor;

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
          float cutVar;
          if (strcmp(cutVariables[j], brName) != 0) 
            cutVar = cutVarMap.at(cutVariables[j]);
          else
            cutVar = var;
          if (cutVar < cutLows[j] || cutVar >= cutHighs[j]) {
            // decide if samples should be ignored if one cut fails or not
            if (flags & 1<<11) {
              totMcWeight = 0.;
              break;
            } else { 
              continue;
            }
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

    if (flags & 1<<10) {
      topRewSf_sum += topRewSf;
      scaleFactor *= topRewSf;
    }

    if (normToBinWidth) scaleFactor /= histo->GetBinWidth(histo->FindBin(var));
    ++fills;
    histo->Fill(var, scaleFactor);
  }

  // fix the normalisation after top reweighting
  if (flags & 1<<10) {
    cout << "Top reweighting overall scale factor from " << fills << " events: " << fills/topRewSf_sum << endl;
    histo->Scale(fills/topRewSf_sum);
  }

  //cout << "integral: " << histo->Integral() << "       overflow: " << histo->GetBinContent(histo->GetNbinsX() + 1) << "       underflow: " << histo->GetBinContent(0) << endl;
  return histo;
}

