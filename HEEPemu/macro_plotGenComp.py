#!/usr/bin/env python

import ROOT
import math
#ROOT.gROOT.SetBatch(True)
from ROOT import TMath,TFile,TDCacheFile,TTree,TH1F,TH1D
from ROOT import TCanvas,TLegend,TLatex
from ROOT import gROOT,gStyle,gPad

savePlots = True
brNames = ['genmu_eta', 'genele_eta', 'genmu_pt', 'genele_pt']
xAxisTitles = ['#eta_{#mu}^{gen}', '#eta_{e}^{gen}', 'p_{T #mu}^{gen}', 'p_{T e}^{gen}']
cuts = ['', 'genele_charge[0]<0', 'genele_charge[0]>0']
cutTitles = ['', '_e-#mu+', '_e+#mu-']
evtSelLabels = ['', 'e^{-}#mu^{+} events', 'e^{+}#mu^{-} events']

font = 42
gStyle.SetTitleFont(font)
gStyle.SetStatFont(font)
gStyle.SetTextFont(font)
gStyle.SetLabelFont(font)
gStyle.SetLegendFont(font)
gStyle.SetMarkerStyle(20)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

#open files to comare
inMzMa1000 = TDCacheFile("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/treeNursery/ZprimeToEMu_M1000_TuneZ2star_8TeV_madgraph_v2_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1_USER_99831ev.root");
inMz1000 = TDCacheFile("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/treeNursery/ZprimeToEMu_Mz-1000_Ma-20000_TuneZ2star_8TeV_madgraph_v1_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1-9fd29d4694dc4f5962d161df214899f8_USER_9999ev.root");
inMa1000 = TDCacheFile("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/treis/treeNursery/ZprimeToEMu_Ma-1000_Mz-20000_TuneZ2star_8TeV_madgraph_v1_treis-MCRECO_Private13_DR53X_PU_S10_START53_V19E-v1-9fd29d4694dc4f5962d161df214899f8_USER_9998ev.root");

#get the trees from the root files
inMzMa1000.cd()
treeMzMa1000 = inMzMa1000.Get('gsfcheckerjob/tree')
inMz1000.cd()
treeMz1000 = inMz1000.Get('gsfcheckerjob/tree')
inMa1000.cd()
treeMa1000 = inMa1000.Get('gsfcheckerjob/tree')

# list with all canvases
cList = []

for brName, xAxisTitle in zip(brNames, xAxisTitles):
    for cut, cutTitle, evtSelLabel in zip(cuts, cutTitles, evtSelLabels):
        cList.append(TCanvas(brName+cutTitle, brName+cutTitle, 100, 100, 600, 600))
        cList[len(cList)-1].cd()
        
        #get the histograms from the tree
        treeMzMa1000.Draw(brName+'[0]>>hMzMa1000_{0}{1}(30)'.format(brName, cutTitle), cut)
        hMzMa1000 = gROOT.FindObject('hMzMa1000_{0}{1}'.format(brName, cutTitle))
        treeMz1000.Draw(brName+'[0]>>hMz1000_{0}{1}({2},{3},{4})'.format(brName, cutTitle, hMzMa1000.GetNbinsX(), hMzMa1000.GetXaxis().GetXmin(), hMzMa1000.GetXaxis().GetXmax()), cut)
        hMz1000 = gROOT.FindObject('hMz1000_{0}{1}'.format(brName, cutTitle))
        treeMa1000.Draw(brName+'[0]>>hMa1000_{0}{1}({2},{3},{4})'.format(brName, cutTitle, hMzMa1000.GetNbinsX(), hMzMa1000.GetXaxis().GetXmin(), hMzMa1000.GetXaxis().GetXmax()), cut)
        hMa1000 = gROOT.FindObject('hMa1000_{0}{1}'.format(brName, cutTitle))
        
        #normalise
        hMzMa1000.Sumw2()
        hMz1000.Sumw2()
        hMa1000.Sumw2()
        print 'N_hMzMa1000: {0}, N_hMz1000: {1}, N_hMa1000: {2}'.format(hMzMa1000.Integral(), hMz1000.Integral(), hMa1000.Integral())
        hMzMa1000.Scale(1./hMzMa1000.Integral())
        hMz1000.Scale(1./hMz1000.Integral())
        hMa1000.Scale(1./hMa1000.Integral())
        #print 'N_hMzMa1000: {0}, N_hMz1000: {1}, N_hMa1000: {2}'.format(hMzMa1000.Integral(), hMz1000.Integral(), hMa1000.Integral())
        
        #get the proper values for the y axis
        minimum = hMzMa1000.GetMinimum()
        if hMz1000.GetMinimum() < minimum: minimum = hMz1000.GetMinimum()
        if hMa1000.GetMinimum() < minimum: minimum = hMa1000.GetMinimum()
        maximum = hMzMa1000.GetMaximum()
        if hMz1000.GetMaximum() > maximum: maximum = hMz1000.GetMaximum()
        if hMa1000.GetMaximum() > maximum: maximum = hMa1000.GetMaximum()
        
        #now draw
        gPad.SetLeftMargin(0.14)
        gPad.SetRightMargin(0.06)
        gPad.SetTickx(1)
        gPad.SetTicky(1)
        
        hMzMa1000.SetMarkerStyle(21)
        hMzMa1000.SetMarkerColor(1)
        hMzMa1000.SetLineWidth(2)
        hMzMa1000.SetLineColor(ROOT.kBlack)
        hMzMa1000.Draw('e0')
        hMzMa1000.GetYaxis().SetRangeUser(0.95*minimum, 1.05*maximum)
        hMzMa1000.GetXaxis().SetTitle(xAxisTitle)
        hMzMa1000.GetXaxis().SetTitleFont(font)
        hMzMa1000.GetXaxis().SetLabelFont(font)
        hMzMa1000.GetYaxis().SetTitle('Normalised distribution')
        hMzMa1000.GetYaxis().SetTitleOffset(1.5)
        hMzMa1000.GetYaxis().SetTitleFont(font)
        hMzMa1000.GetYaxis().SetLabelFont(font)
        hMz1000.SetMarkerStyle(20)
        hMz1000.SetMarkerColor(2)
        hMz1000.SetLineWidth(2)
        hMz1000.SetLineColor(ROOT.kRed)
        hMz1000.Draw('e0same')
        hMa1000.SetMarkerStyle(24)
        hMa1000.SetMarkerColor(4)
        hMa1000.SetLineWidth(2)
        hMa1000.SetLineColor(ROOT.kBlue)
        hMa1000.Draw('e0same')
        
        legend = TLegend(0.16, 0.54, 0.38, 0.86)
        legend.SetTextFont(font)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillColor(19)
        legend.SetFillStyle(0)
        legend.AddEntry(hMzMa1000, "#splitline{M_{Z'} = 1 TeV,}{M_{a'} = 1 TeV}")
        legend.AddEntry(hMz1000, "#splitline{M_{Z'} = 1 TeV,}{M_{a'} = 20 TeV}")
        legend.AddEntry(hMa1000, "#splitline{M_{Z'} = 20 TeV,}{M_{a'} = 1 TeV}")
        legend.Draw('same')
        
        tex = TLatex()
        tex.SetNDC()
        tex.SetTextFont(font)
        tex.SetTextSize(0.04)
        tex.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')
        tex.DrawLatex(0.767, 0.91, evtSelLabel)

        legend.Draw('same')

        cList[len(cList)-1].Modified()
        cList[len(cList)-1].Update()

# save canvases to root file
if savePlots:
    output = TFile('./plots/signalComparisonPlots.root', 'recreate')
    output.cd()
    for canvas in cList:
        canvas.Write(canvas.GetName())
    output.Close()

