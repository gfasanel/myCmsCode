#!/usr/bin/env python

import ROOT
import math
#ROOT.gROOT.SetBatch(True)
from ROOT import TFile,TTree,TH1F,TH1D,TF1,TGraphErrors
from ROOT import TCanvas,TLegend,TLatex,TPaveStats
from ROOT import gROOT,gStyle,gPad

savePlots = True
font = 42
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
input = TFile("emuSpec_19703pb-1.root")
treePrefix = 'emuTree_sigNoAccCuts'

#get the trees from the root files
input.cd()

# list with all canvases
cList = []

# check how many trees to look at are in the file
nTrees = 0
keyList = input.GetListOfKeys()
for key in keyList:
    if key.GetName()[:len(treePrefix)] == treePrefix:
       ++nTrees

res = TGraphErrors(nTrees)

# use all the signal trees in the file
for key in keyList:
    if key.GetName()[:len(treePrefix)] == treePrefix:
        print key.GetName()
        massPt = int(key.GetName()[len(treePrefix):])

        tree = input.Get(key.GetName())

        #brName = 'mass'
        brName = '(mass-emu_mass)/emu_mass'

        cList.append(TCanvas(key.GetName()+'_res', key.GetName()+'_res', 600, 600))
        i = len(cList)-1

        #get the histograms from the tree
        tree.Draw(brName+'>>{0}_hMassDiff(500)'.format(key.GetName()))
        hMassDiff = gROOT.FindObject('{0}_hMassDiff'.format(key.GetName()))

        #normalise
        hMassDiff.Sumw2()
        #hMassDiff.Scale(1./hMassDiff.Integral())

        #now draw
        gPad.SetLeftMargin(0.14)
        gPad.SetRightMargin(0.06)
        gPad.SetTickx(1)
        gPad.SetTicky(1)

        hMassDiff.SetMarkerStyle(20)
        hMassDiff.SetMarkerColor(ROOT.kBlue)
        #hMassDiff.SetLineWidth(2)
        hMassDiff.SetLineColor(ROOT.kBlue)
        hMassDiff.Draw('e0')
        hMassDiff.GetXaxis().SetTitle('(M_{e#mu}^{reco} - M_{e#mu}^{gen})/M_{e#mu}^{gen}')
        hMassDiff.GetXaxis().SetTitleOffset(1.05)
        hMassDiff.GetXaxis().SetTitleFont(font)
        hMassDiff.GetXaxis().SetLabelFont(font)
        hMassDiff.GetYaxis().SetTitle('Events')
        hMassDiff.GetYaxis().SetTitleOffset(1.5)
        hMassDiff.GetYaxis().SetTitleFont(font)
        hMassDiff.GetYaxis().SetLabelFont(font)

        gauss = TF1('gauss', 'gaus')
        gauss.SetLineColor(ROOT.kRed)
        peak = hMassDiff.GetBinCenter(hMassDiff.GetMaximumBin())
        gauss.SetParameter(1, peak)
        gauss.SetParameter(2, hMassDiff.GetRMS())
        factorFit1 = 1.
        factorFit2 = 1.8
        # raw and fine fit
        hMassDiff.Fit(gauss, '', '', peak-factorFit1*hMassDiff.GetRMS(), peak+factorFit1*hMassDiff.GetRMS())
        hMassDiff.Fit(gauss, '', '', gauss.GetParameter(1)-factorFit2*gauss.GetParameter(2), gauss.GetParameter(1)+factorFit2*gauss.GetParameter(2))

        #hMassDiff.GetXaxis().SetRangeUser(math.floor(gauss.GetParameter(1)-4.*gauss.GetParameter(2)), math.floor(gauss.GetParameter(1)+4.*gauss.GetParameter(2)))
        if massPt < 1000:
            hMassDiff.GetXaxis().SetRangeUser(-0.1, 0.1)
        elif massPt < 2000:
            hMassDiff.GetXaxis().SetRangeUser(-0.2, 0.2)
        else:
            hMassDiff.GetXaxis().SetRangeUser(-0.4, 0.4)
        cList[i].Modified()
        cList[i].Update()

        #res.SetPoint(i, massPt, gauss.GetParameter(2)/massPt)
        #res.SetPointError(i, 0., gauss.GetParError(2)/massPt)
        res.SetPoint(i, massPt, gauss.GetParameter(2))
        res.SetPointError(i, 0., gauss.GetParError(2))

        #legend = TLegend(0.16, 0.54, 0.38, 0.86)
        #legend.SetTextFont(font)
        #legend.SetTextSize(0.03)
        #legend.SetBorderSize(0)
        #legend.SetFillColor(19)
        #legend.SetFillStyle(0)
        #legend.AddEntry(hMassDiff, "#splitline{M_{Z'} = 1 TeV,}{M_{a'} = 1 TeV}")
        #legend.Draw('same')

        tex = TLatex()
        tex.SetNDC()
        tex.SetTextFont(font)
        tex.SetTextSize(0.04)
        tex.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')
        tex.SetTextSize(0.03)
        tex.DrawLatex(0.68, 0.72, "M_{Z'} = " + str(massPt) + ' GeV')

c1 = TCanvas('resolution', 'resolution', 600, 600)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)
res.SetMarkerStyle(20)
res.SetMarkerColor(ROOT.kBlue)
res.SetLineColor(ROOT.kBlue)
res.GetXaxis().SetTitle('M_{e#mu} (GeV)')
res.GetXaxis().SetTitleFont(font)
res.GetXaxis().SetLabelFont(font)
res.GetYaxis().SetTitle('#sigma(M_{e#mu})/M_{e#mu}')
res.GetYaxis().SetTitleOffset(1.55)
res.GetYaxis().SetTitleFont(font)
res.GetYaxis().SetLabelFont(font)
res.Draw('ap')
resFitFunc = TF1('resFitFunc', 'pol2', 250, 5000);
resFitFunc.SetLineColor(ROOT.kRed)
res.Fit(resFitFunc, '', '', 250, 5000)

# move the stats box
gPad.Update()
statsBox = res.GetListOfFunctions().FindObject('stats')
statsBox.SetX1NDC(0.18)
statsBox.SetY1NDC(0.66)
statsBox.SetX2NDC(0.54)
statsBox.SetY2NDC(0.86)

tex2 = TLatex()
tex2.SetNDC()
tex2.SetTextFont(font)
tex2.SetTextSize(0.04)
tex2.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')

c1.Modified()
c1.Update()

# save canvases to root file
if savePlots:
    output = TFile('./plots/resolutionPlots.root', 'recreate')
    output.cd()
    for canvas in cList:
        canvas.Write(canvas.GetName())
    c1.Write(c1.GetName())
    output.Close()

