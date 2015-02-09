#!/usr/bin/env python

import ROOT
import math
#ROOT.gROOT.SetBatch(True)
from ROOT import TFile,TTree,TH1F,TH1D,TF1,TGraph,TGraphErrors
from ROOT import TH2F
from ROOT import TCanvas,TLegend,TLatex,TPaveStats
from ROOT import gROOT,gStyle,gPad
from array import array

savePlots = True
#savePlots = False
font = 42
markers = [ROOT.kFullCircle, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown]
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen]
gStyle.SetTitleFont(font)
gStyle.SetStatFont(font)
gStyle.SetTextFont(font)
gStyle.SetLabelFont(font)
gStyle.SetLegendFont(font)
gStyle.SetMarkerStyle(20)
gStyle.SetOptStat(0)
gStyle.SetOptFit(1)
gStyle.SetOptTitle(0)

#open files to comare
input = TFile("eScaleEvents19712pb-1_thesis_hovere.root")
treeNames = ['eleDY20Tree',
             'eleDY120Tree',
             'eleDY200Tree',
             'eleDY400Tree',
             'eleDY500Tree',
             'eleDY700Tree',
             'eleDY800Tree',
             'eleDY1000Tree',
             'eleDY1500Tree',
             'eleDY2000Tree']

trees = []

massRanges = [100.,
              120.,
              200.,
              300.,
              400.,
              500.,
              600.,
              700.,
              800.,
              1000.,
              1250.,
              1500.,
              1750.,
              2000.,
              2250.,
              2500.,
              3500.]

hEle1DrMatched_eta10to1p442_eta20to1p442_List = []
hEle1DrMatched_eta10to0p7_eta20to0p7_List = []
hEle1DrMatched_eta10p7to1p442_eta20p7to1p442_List = []
hEle1DrMatched_eta10to0p7_eta20p7to1p442_List = []
hEle2DrMatched_eta10to1p442_eta20to1p442_List = []
hEle2DrMatched_eta10to0p7_eta20to0p7_List = []
hEle2DrMatched_eta10p7to1p442_eta20p7to1p442_List = []
hEle2DrMatched_eta10to0p7_eta20p7to1p442_List = []

# define the histograms
for i in range(len(massRanges)-1):
    massLow = massRanges[i]
    massHigh = massRanges[i+1]
    namePostfix = '_m'+str(massRanges[i])+'tom'+str(massRanges[i+1])
    #print name
    hEle1DrMatched_eta10to1p442_eta20to1p442_List.append(TH1F('hEle1DrMatched_eta10to1p442_eta20to1p442'+namePostfix, 'hEle1DrMatched_eta10to1p442_eta20to1p442'+namePostfix, 100, 0., 0.05))
    hEle1DrMatched_eta10to0p7_eta20to0p7_List.append(TH1F('hEle1DrMatched_eta10to0p7_eta20to0p7'+namePostfix, 'hEle1DrMatched_eta10to0p7_eta20to0p7'+namePostfix, 100, 0., 0.05))
    hEle1DrMatched_eta10p7to1p442_eta20p7to1p442_List.append(TH1F('hEle1DrMatched_eta10p7to1p442_eta20p7to1p442'+namePostfix, 'hEle1DrMatched_eta10p7to1p442_eta20p7to1p442'+namePostfix, 100, 0., 0.05))
    hEle1DrMatched_eta10to0p7_eta20p7to1p442_List.append(TH1F('hEle1DrMatched_eta10to0p7_eta20p7to1p442'+namePostfix, 'hEle1DrMatched_eta10to0p7_eta20p7to1p442'+namePostfix, 100, 0., 0.05))
    hEle2DrMatched_eta10to1p442_eta20to1p442_List.append(TH1F('hEle2DrMatched_eta10to1p442_eta20to1p442'+namePostfix, 'hEle2DrMatched_eta10to1p442_eta20to1p442'+namePostfix, 100, 0., 0.05))
    hEle2DrMatched_eta10to0p7_eta20to0p7_List.append(TH1F('hEle2DrMatched_eta10to0p7_eta20to0p7'+namePostfix, 'hEle2DrMatched_eta10to0p7_eta20to0p7'+namePostfix, 100, 0., 0.05))
    hEle2DrMatched_eta10p7to1p442_eta20p7to1p442_List.append(TH1F('hEle2DrMatched_eta10p7to1p442_eta20p7to1p442'+namePostfix, 'hEle2DrMatched_eta10p7to1p442_eta20p7to1p442'+namePostfix, 100, 0., 0.05))
    hEle2DrMatched_eta10to0p7_eta20p7to1p442_List.append(TH1F('hEle2DrMatched_eta10to0p7_eta20p7to1p442'+namePostfix, 'hEle2DrMatched_eta10to0p7_eta20p7to1p442'+namePostfix, 100, 0., 0.05))

# define 2D histograms
h2Ele1DrMatched_eta10to0p7_eta20to0p7 = TH2F('h2Ele1DrMatched_eta10to0p7_eta20to0p7', 'h2Ele1DrMatched_eta10to0p7_eta20to0p7', len(massRanges)-1, array('d',massRanges), 30, 0., 0.01)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetXaxis().SetTitle('m(ee) (GeV)')
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetXaxis().SetTitleFont(font)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetXaxis().SetLabelFont(font)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetXaxis().SetRangeUser(massRanges[0], massRanges[len(massRanges)-1])
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetYaxis().SetTitle('H / E electron 1')
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetYaxis().SetTitleOffset(2.1)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetYaxis().SetTitleFont(font)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetYaxis().SetLabelFont(font)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetZaxis().SetTitleFont(font)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetZaxis().SetLabelFont(font)
h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetZaxis().SetRangeUser(0., 0.2)

h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442 = h2Ele1DrMatched_eta10to0p7_eta20to0p7.Clone('h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442')
h2Ele1DrMatched_eta10to0p7_eta20p7to1p442 = h2Ele1DrMatched_eta10to0p7_eta20to0p7.Clone('h2Ele1DrMatched_eta10to0p7_eta20p7to1p442')


#get the trees from the root files
input.cd()

for treeName in treeNames:
    trees.append(input.Get(treeName))

for tree in trees:
    print tree.GetName()
    nEv = tree.GetEntries()
    for ev in range(nEv):
        if ev%100000 == 0: 
            print ev
        tree.GetEntry(ev)
        massBin = 0
        #find the mass bin for this event
        for i in range(1,len(massRanges)):
            if tree.mass < massRanges[i]:
                massBin = i-1
                break
        if abs(tree.ele1EtaDrMatched) < 1.442 and abs(tree.ele2EtaDrMatched) < 1.442:
            hEle1DrMatched_eta10to1p442_eta20to1p442_List[massBin].Fill(tree.ele1HoEDrMatched)
            hEle2DrMatched_eta10to1p442_eta20to1p442_List[massBin].Fill(tree.ele2HoEDrMatched)
            if abs(tree.ele1EtaDrMatched) < 0.7 and abs(tree.ele2EtaDrMatched) < 0.7:
                hEle1DrMatched_eta10to0p7_eta20to0p7_List[massBin].Fill(tree.ele1HoEDrMatched)
                hEle2DrMatched_eta10to0p7_eta20to0p7_List[massBin].Fill(tree.ele2HoEDrMatched)
                h2Ele1DrMatched_eta10to0p7_eta20to0p7.Fill(tree.mass, tree.ele1HoEDrMatched)
            if abs(tree.ele1EtaDrMatched) >= 0.7 and abs(tree.ele2EtaDrMatched) >= 0.7:
                hEle1DrMatched_eta10p7to1p442_eta20p7to1p442_List[massBin].Fill(tree.ele1HoEDrMatched)
                hEle2DrMatched_eta10p7to1p442_eta20p7to1p442_List[massBin].Fill(tree.ele2HoEDrMatched)
                h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.Fill(tree.mass, tree.ele1HoEDrMatched)
            if abs(tree.ele1EtaDrMatched) < 0.7 and abs(tree.ele2EtaDrMatched) >= 0.7:
                hEle1DrMatched_eta10to0p7_eta20p7to1p442_List[massBin].Fill(tree.ele1HoEDrMatched)
                hEle2DrMatched_eta10to0p7_eta20p7to1p442_List[massBin].Fill(tree.ele2HoEDrMatched)
                h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.Fill(tree.mass, tree.ele1HoEDrMatched)

# make graphs
gEle1DrMatched_eta10to1p442_eta20to1p442_mean = TGraphErrors(len(massRanges)-1)
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetName('gEle1DrMatched_eta10to1p442_eta20to1p442_mean') 
gEle1DrMatched_eta10to0p7_eta20to0p7_mean = TGraphErrors(len(massRanges)-1)
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetName('gEle1DrMatched_eta10to0p7_eta20to0p7_mean') 
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean = TGraphErrors(len(massRanges)-1)
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetName('gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean') 
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean = TGraphErrors(len(massRanges)-1)
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetName('gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean') 

gEle2DrMatched_eta10to1p442_eta20to1p442_mean = TGraphErrors(len(massRanges)-1)
gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetName('gEle2DrMatched_eta10to1p442_eta20to1p442_mean') 
gEle2DrMatched_eta10to0p7_eta20to0p7_mean = TGraphErrors(len(massRanges)-1)
gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetName('gEle2DrMatched_eta10to0p7_eta20to0p7_mean') 
gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean = TGraphErrors(len(massRanges)-1)
gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetName('gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean') 
gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean = TGraphErrors(len(massRanges)-1)
gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetName('gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean') 
for i in range(len(massRanges)-1):
    gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle1DrMatched_eta10to1p442_eta20to1p442_List[i].GetMean())
    gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle1DrMatched_eta10to1p442_eta20to1p442_List[i].GetMean()/hEle1DrMatched_eta10to1p442_eta20to1p442_List[i].Integral())
    gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle1DrMatched_eta10to0p7_eta20to0p7_List[i].GetMean())
    gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle1DrMatched_eta10to0p7_eta20to0p7_List[i].GetMean()/hEle1DrMatched_eta10to0p7_eta20to0p7_List[i].Integral())
    gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle1DrMatched_eta10p7to1p442_eta20p7to1p442_List[i].GetMean())
    gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle1DrMatched_eta10p7to1p442_eta20p7to1p442_List[i].GetMean()/hEle1DrMatched_eta10p7to1p442_eta20p7to1p442_List[i].Integral())
    gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle1DrMatched_eta10to0p7_eta20p7to1p442_List[i].GetMean())
    gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle1DrMatched_eta10to0p7_eta20p7to1p442_List[i].GetMean()/hEle1DrMatched_eta10to0p7_eta20p7to1p442_List[i].Integral())

    gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle2DrMatched_eta10to1p442_eta20to1p442_List[i].GetMean())
    gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle2DrMatched_eta10to1p442_eta20to1p442_List[i].GetMean()/hEle2DrMatched_eta10to1p442_eta20to1p442_List[i].Integral())
    gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle2DrMatched_eta10to0p7_eta20to0p7_List[i].GetMean())
    gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle2DrMatched_eta10to0p7_eta20to0p7_List[i].GetMean()/hEle2DrMatched_eta10to0p7_eta20to0p7_List[i].Integral())
    gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle2DrMatched_eta10p7to1p442_eta20p7to1p442_List[i].GetMean())
    gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle2DrMatched_eta10p7to1p442_eta20p7to1p442_List[i].GetMean()/hEle2DrMatched_eta10p7to1p442_eta20p7to1p442_List[i].Integral())
    gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetPoint(i, (massRanges[i] + massRanges[i+1]) / 2., hEle2DrMatched_eta10to0p7_eta20p7to1p442_List[i].GetMean())
    gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetPointError(i, (massRanges[i+1] - massRanges[i]) / 2., hEle2DrMatched_eta10to0p7_eta20p7to1p442_List[i].GetMean()/hEle2DrMatched_eta10to0p7_eta20p7to1p442_List[i].Integral())

#normalise 2D histo
for xbin in range(1,h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetNbinsX()+1):
    int = h2Ele1DrMatched_eta10to0p7_eta20to0p7.Integral(xbin,xbin,1,h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetNbinsY()+1)
    for ybin in range(1, h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetNbinsY()+1):
        if int > 0.: 
            h2Ele1DrMatched_eta10to0p7_eta20to0p7.SetBinContent(xbin,ybin,h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetBinContent(xbin,ybin)/int)
            h2Ele1DrMatched_eta10to0p7_eta20to0p7.SetBinError(xbin,ybin,h2Ele1DrMatched_eta10to0p7_eta20to0p7.GetBinError(xbin,ybin)/int)
for xbin in range(1,h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.GetNbinsX()+1):
    int = h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.Integral(xbin,xbin,1,h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.GetNbinsY()+1)
    for ybin in range(1, h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.GetNbinsY()+1):
        if int > 0.: 
            h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.SetBinContent(xbin,ybin,h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.GetBinContent(xbin,ybin)/int)
            h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.SetBinError(xbin,ybin,h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.GetBinError(xbin,ybin)/int)
for xbin in range(1,h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.GetNbinsX()+1):
    int = h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.Integral(xbin,xbin,1,h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.GetNbinsY()+1)
    for ybin in range(1, h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.GetNbinsY()+1):
        if int > 0.: 
            h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.SetBinContent(xbin,ybin,h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.GetBinContent(xbin,ybin)/int)
            h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.SetBinError(xbin,ybin,h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.GetBinError(xbin,ybin)/int)

###############################################################################
## plot the histo
c0 = TCanvas('tests', 'tests', 100, 20, 600, 600)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

hEle1DrMatched_eta10to0p7_eta20to0p7_List[9].Draw('hist')

###############################################################################
c1 = TCanvas('eta10to1p442_eta20to1p442', 'eta10to1p442_eta20to1p442', 100, 20, 600, 600)
gPad.SetLeftMargin(0.18)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetMarkerStyle(markers[0])
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetMarkerColor(colors[0])
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetLineColor(colors[0])
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.SetLineWidth(2)
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.Draw('APZ')
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetXaxis().SetTitle('m(ee) (GeV)')
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetXaxis().SetTitleFont(font)
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetXaxis().SetLabelFont(font)
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetXaxis().SetRangeUser(massRanges[0], massRanges[len(massRanges)-1])
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetYaxis().SetTitle('#LTH / E#GT')
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetYaxis().SetTitleOffset(2.3)
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetYaxis().SetTitleFont(font)
gEle1DrMatched_eta10to1p442_eta20to1p442_mean.GetYaxis().SetLabelFont(font)

gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetMarkerStyle(markers[1])
gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetMarkerColor(colors[1])
gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetLineColor(colors[1])
gEle2DrMatched_eta10to1p442_eta20to1p442_mean.SetLineWidth(2)
gEle2DrMatched_eta10to1p442_eta20to1p442_mean.Draw('PZ')

legend1 = TLegend(0.23, 0.73, 0.44, 0.86)
legend1.SetTextFont(font)
legend1.SetTextSize(0.04)
legend1.SetBorderSize(0)
legend1.SetFillColor(19)
legend1.SetFillStyle(0)
legend1.AddEntry(gEle1DrMatched_eta10to1p442_eta20to1p442_mean, 'e_{1}', 'lp')
legend1.AddEntry(gEle2DrMatched_eta10to1p442_eta20to1p442_mean, 'e_{2}', 'lp')
legend1.Draw('same')

tex1 = TLatex()
tex1.SetNDC()
tex1.SetTextFont(font)
tex1.SetTextSize(0.04)
tex1.DrawLatex(0.18, 0.91, 'CMS Simulation, 8 TeV')
tex1.DrawLatex(0.5, 0.16, 'two e in barrel')

###############################################################################
c2 = TCanvas('eta10to0p7_eta20to0p7', 'eta10to0p7_eta20to0p7', 100, 20, 600, 600)
gPad.SetLeftMargin(0.18)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetMarkerStyle(markers[0])
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetMarkerColor(colors[0])
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetLineColor(colors[0])
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.SetLineWidth(2)
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.Draw('APZ')
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetXaxis().SetTitle('m(ee) (GeV)')
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetXaxis().SetTitleFont(font)
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetXaxis().SetLabelFont(font)
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetXaxis().SetRangeUser(massRanges[0], massRanges[len(massRanges)-1])
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetYaxis().SetTitle('#LTH / E#GT')
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetYaxis().SetTitleOffset(2.3)
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetYaxis().SetTitleFont(font)
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.GetYaxis().SetLabelFont(font)

gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetMarkerStyle(markers[1])
gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetMarkerColor(colors[1])
gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetLineColor(colors[1])
gEle2DrMatched_eta10to0p7_eta20to0p7_mean.SetLineWidth(2)
gEle2DrMatched_eta10to0p7_eta20to0p7_mean.Draw('PZ')

legend2 = TLegend(0.23, 0.73, 0.44, 0.86)
legend2.SetTextFont(font)
legend2.SetTextSize(0.04)
legend2.SetBorderSize(0)
legend2.SetFillColor(19)
legend2.SetFillStyle(0)
legend2.AddEntry(gEle1DrMatched_eta10to0p7_eta20to0p7_mean, 'e_{1}', 'lp')
legend2.AddEntry(gEle2DrMatched_eta10to0p7_eta20to0p7_mean, 'e_{2}', 'lp')
legend2.Draw('same')

tex2 = TLatex()
tex2.SetNDC()
tex2.SetTextFont(font)
tex2.SetTextSize(0.04)
tex2.DrawLatex(0.18, 0.91, 'CMS Simulation, 8 TeV')
tex2.DrawLatex(0.5, 0.16, 'two e with |#eta| < 0.7')

###############################################################################
c3 = TCanvas('eta10p7to1p442_eta20p7to1p442', 'eta10p7to1p442_eta20p7to1p442', 100, 20, 600, 600)
gPad.SetLeftMargin(0.18)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetMarkerStyle(markers[0])
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetMarkerColor(colors[0])
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetLineColor(colors[0])
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetLineWidth(2)
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.Draw('APZ')
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetXaxis().SetTitle('m(ee) (GeV)')
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetXaxis().SetTitleFont(font)
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetXaxis().SetLabelFont(font)
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetXaxis().SetRangeUser(massRanges[0], massRanges[len(massRanges)-1])
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetYaxis().SetTitle('#LTH / E#GT')
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetYaxis().SetTitleOffset(2.3)
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetYaxis().SetTitleFont(font)
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.GetYaxis().SetLabelFont(font)

gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetMarkerStyle(markers[1])
gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetMarkerColor(colors[1])
gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetLineColor(colors[1])
gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.SetLineWidth(2)
gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean.Draw('PZ')

legend3 = TLegend(0.23, 0.73, 0.44, 0.86)
legend3.SetTextFont(font)
legend3.SetTextSize(0.04)
legend3.SetBorderSize(0)
legend3.SetFillColor(19)
legend3.SetFillStyle(0)
legend3.AddEntry(gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean, 'e_{1}', 'lp')
legend3.AddEntry(gEle2DrMatched_eta10p7to1p442_eta20p7to1p442_mean, 'e_{2}', 'lp')
legend3.Draw('same')

tex3 = TLatex()
tex3.SetNDC()
tex3.SetTextFont(font)
tex3.SetTextSize(0.04)
tex3.DrawLatex(0.18, 0.91, 'CMS Simulation, 8 TeV')
tex3.DrawLatex(0.5, 0.16, 'two e with 0.7 #leq |#eta| < 1.442')

###############################################################################
c4 = TCanvas('eta10to0p7_eta20p7to1p442', 'eta10to0p7_eta20p7to1p442', 100, 20, 600, 600)
gPad.SetLeftMargin(0.18)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetMarkerStyle(markers[0])
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetMarkerColor(colors[0])
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetLineColor(colors[0])
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.SetLineWidth(2)
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.Draw('APZ')
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetXaxis().SetTitle('m(ee) (GeV)')
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetXaxis().SetTitleFont(font)
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetXaxis().SetLabelFont(font)
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetXaxis().SetRangeUser(massRanges[0], massRanges[len(massRanges)-1])
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetYaxis().SetTitle('#LTH / E#GT')
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetYaxis().SetTitleOffset(2.3)
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetYaxis().SetTitleFont(font)
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.GetYaxis().SetLabelFont(font)

gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetMarkerStyle(markers[1])
gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetMarkerColor(colors[1])
gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetLineColor(colors[1])
gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.SetLineWidth(2)
gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean.Draw('PZ')

legend4 = TLegend(0.23, 0.73, 0.44, 0.86)
legend4.SetTextFont(font)
legend4.SetTextSize(0.04)
legend4.SetBorderSize(0)
legend4.SetFillColor(19)
legend4.SetFillStyle(0)
legend4.AddEntry(gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean, 'e_{1}', 'lp')
legend4.AddEntry(gEle2DrMatched_eta10to0p7_eta20p7to1p442_mean, 'e_{2}', 'lp')
legend4.Draw('same')

tex4 = TLatex()
tex4.SetNDC()
tex4.SetTextFont(font)
tex4.SetTextSize(0.04)
tex4.DrawLatex(0.18, 0.91, 'CMS Simulation, 8 TeV')
tex4.DrawLatex(0.5, 0.2, 'one e with |#eta| < 0.7')
tex4.DrawLatex(0.5, 0.15, 'one e with 0.7 #leq |#eta| < 1.442')

###############################################################################
c21 = TCanvas('2d_eta10to0p7_eta20to0p7', '2d_eta10to0p7_eta20to0p7', 100, 20, 600, 600)
gStyle.SetPalette(53)
gStyle.SetNumberContours(99)
gPad.SetLeftMargin(0.16)
gPad.SetRightMargin(0.14)

h2Ele1DrMatched_eta10to0p7_eta20to0p7.Draw('COLZ')
gEle1DrMatched_eta10to0p7_eta20to0p7_mean.Draw('PZ')

legend21 = TLegend(0.19, 0.69, 0.60, 0.74)
legend21.SetTextFont(font)
legend21.SetTextSize(0.04)
legend21.SetTextColor(ROOT.kWhite)
legend21.SetBorderSize(0)
legend21.SetFillColor(19)
legend21.SetFillStyle(0)
legend21.AddEntry(gEle1DrMatched_eta10to0p7_eta20to0p7_mean, '#LTH / E#GT of m(ee) slice', 'lp')
legend21.Draw('same')

tex21 = TLatex()
tex21.SetNDC()
tex21.SetTextFont(font)
tex21.SetTextSize(0.04)
tex21.DrawLatex(0.16, 0.91, 'CMS Simulation, 8 TeV')
tex21.SetTextColor(ROOT.kWhite)
tex21.DrawLatex(0.19, 0.84, 'two e with |#eta| < 0.7')

###############################################################################
c22 = TCanvas('2d_eta10p7to1p442_eta20p7to1p442', '2d_eta10p7to1p442_eta20p7to1p442', 100, 20, 600, 600)
gPad.SetLeftMargin(0.16)
gPad.SetRightMargin(0.14)

h2Ele1DrMatched_eta10p7to1p442_eta20p7to1p442.Draw('COLZ')
gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean.Draw('PZ')

legend22 = TLegend(0.19, 0.69, 0.60, 0.74)
legend22.SetTextFont(font)
legend22.SetTextSize(0.04)
legend22.SetTextColor(ROOT.kWhite)
legend22.SetBorderSize(0)
legend22.SetFillColor(19)
legend22.SetFillStyle(0)
legend22.AddEntry(gEle1DrMatched_eta10p7to1p442_eta20p7to1p442_mean, '#LTH / E#GT of m(ee) slice', 'lp')
legend22.Draw('same')

tex22 = TLatex()
tex22.SetNDC()
tex22.SetTextFont(font)
tex22.SetTextSize(0.04)
tex22.DrawLatex(0.16, 0.91, 'CMS Simulation, 8 TeV')
tex22.SetTextColor(ROOT.kWhite)
tex22.DrawLatex(0.19, 0.84, 'two e with 0.7 #leq |#eta| < 1.442')

###############################################################################
c23 = TCanvas('2d_eta10to0p7_eta20p7to1p442', '2d_eta10to0p7_eta20p7to1p442', 100, 20, 600, 600)
gPad.SetLeftMargin(0.16)
gPad.SetRightMargin(0.14)

h2Ele1DrMatched_eta10to0p7_eta20p7to1p442.Draw('COLZ')
gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean.Draw('PZ')

legend23 = TLegend(0.19, 0.69, 0.60, 0.74)
legend23.SetTextFont(font)
legend23.SetTextSize(0.04)
legend23.SetTextColor(ROOT.kWhite)
legend23.SetBorderSize(0)
legend23.SetFillColor(19)
legend23.SetFillStyle(0)
legend23.AddEntry(gEle1DrMatched_eta10to0p7_eta20p7to1p442_mean, '#LTH / E#GT of m(ee) slice', 'lp')
legend23.Draw('same')

tex23 = TLatex()
tex23.SetNDC()
tex23.SetTextFont(font)
tex23.SetTextSize(0.04)
tex23.DrawLatex(0.16, 0.91, 'CMS Simulation, 8 TeV')
tex23.SetTextColor(ROOT.kWhite)
tex23.DrawLatex(0.19, 0.84, 'one e with |#eta| < 0.7')
tex23.DrawLatex(0.19, 0.79, 'one e with 0.7 #leq |#eta| < 1.442')

##############################################################################
# save canvases to root file
if savePlots:
    output = TFile('./thesisplots_HoverE/hovere_plots.root', 'recreate')
    plotDir = './thesisplots_HoverE/'
    output.cd()

    c1.Write(c1.GetName())
    c1.Print(plotDir+c1.GetName()+'.png', 'png')
    c1.Print(plotDir+c1.GetName()+'.pdf', 'pdf')
    c2.Write(c2.GetName())
    c2.Print(plotDir+c2.GetName()+'.png', 'png')
    c2.Print(plotDir+c2.GetName()+'.pdf', 'pdf')
    c3.Write(c3.GetName())
    c3.Print(plotDir+c3.GetName()+'.png', 'png')
    c3.Print(plotDir+c3.GetName()+'.pdf', 'pdf')
    c4.Write(c4.GetName())
    c4.Print(plotDir+c4.GetName()+'.png', 'png')
    c4.Print(plotDir+c4.GetName()+'.pdf', 'pdf')
    c21.Write(c21.GetName())
    c21.Print(plotDir+c21.GetName()+'.png', 'png')
    c21.Print(plotDir+c21.GetName()+'.pdf', 'pdf')
    c22.Write(c22.GetName())
    c22.Print(plotDir+c22.GetName()+'.png', 'png')
    c22.Print(plotDir+c22.GetName()+'.pdf', 'pdf')
    c23.Write(c23.GetName())
    c23.Print(plotDir+c23.GetName()+'.png', 'png')
    c23.Print(plotDir+c23.GetName()+'.pdf', 'pdf')

    output.Close()

# wait
if not ROOT.gROOT.IsBatch(): 
    raw_input("Press ENTER to quit.")

