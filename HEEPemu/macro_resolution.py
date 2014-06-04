#!/usr/bin/env python

import ROOT
import math
ROOT.gROOT.SetBatch(True)
from ROOT import TFile,TTree,TH1F,TH1D,TF1,TGraph,TGraphErrors
from ROOT import TCanvas,TLegend,TLatex,TPaveStats
from ROOT import gROOT,gStyle,gPad

savePlots = True
font = 42
markers = [ROOT.kFullTriangleUp, ROOT.kFullTriangleDown]
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen]
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
#input = TFile("emuSpec_19703pb-1.root")
input = TFile("emuSpec_singleMuTrg_19706pb-1.root")
treePrefices = ['emuTree_sigNoAccCuts', 'emuTree_sigTagV7C2_']
#treePrefices = ['emuTree_sigNoAccCuts']

#get the trees from the root files
input.cd()

# list with all canvases
cList = []

maxMass = 4000

# check how many trees to look at are in the file
nTrees = 0
keyList = input.GetListOfKeys()
for key in keyList:
    keyName = key.GetName()
    if keyName[:len(treePrefices[0])] == treePrefices[0]:
       if int(keyName[len(treePrefices[0]):]) <= maxMass:
           nTrees += 1

resList = []

for i, treePrefix in enumerate(treePrefices):
    resList.append(TGraphErrors(nTrees))
    pt = 0
    # use all the signal trees in the file
    for key in keyList:
        keyName = key.GetName()
        if keyName[:len(treePrefix)] == treePrefix and keyName[len(treePrefix)].isdigit() and keyName[len(keyName)-1].isdigit():
            massPt = int(keyName[len(treePrefix):])
            # go only up to masses of 4 TeV
            if massPt > 4000:
                continue
    
            print keyName
            tree = input.Get(keyName)
    
            brName = '(mass-emu_mass)/emu_mass'
            #brName = '(mass-genEmuMass)/genEmuMass'
    
            cList.append(TCanvas(keyName+'_res', keyName+'_res', 100, 20, 600, 600))
            j = len(cList)-1
    
            #get the histograms from the tree
            tree.Draw(brName+'>>{0}_hMassDiff(500)'.format(keyName))
            hMassDiff = gROOT.FindObject('{0}_hMassDiff'.format(keyName))
    
            #normalise
            hMassDiff.Sumw2()
            #hMassDiff.Scale(1./hMassDiff.Integral())
    
            #now draw
            gPad.SetLeftMargin(0.14)
            gPad.SetRightMargin(0.06)
            gPad.SetTickx(1)
            gPad.SetTicky(1)
    
            hMassDiff.SetMarkerStyle(markers[i])
            hMassDiff.SetMarkerColor(colors[0])
            #hMassDiff.SetLineWidth(2)
            hMassDiff.SetLineColor(colors[0])
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
            gauss.SetLineColor(colors[1])
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
            cList[j].Modified()
            cList[j].Update()
    
            #resList[i].SetPoint(pt, massPt, gauss.GetParameter(2)/massPt)
            #resList[i].SetPointError(pt, 0., gauss.GetParError(2)/massPt)
            resList[i].SetPoint(pt, massPt, gauss.GetParameter(2))
            resList[i].SetPointError(pt, 0., gauss.GetParError(2))
    
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

            pt += 1

##############################################################################
# plot the resolution
c1 = TCanvas('resolution', 'resolution', 100, 20, 600, 600)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

resFitFuncList = []
statsBoxList = []
statsBoxCoord = [[0.18, 0.66, 0.54, 0.86], [0.54, 0.14, 0.90, 0.34], [0.18, 0.76, 0.54, 0.86]]

for i, treePrefix in enumerate(treePrefices):
    resList[i].SetMarkerStyle(markers[i])
    resList[i].SetMarkerColor(colors[i])
    resList[i].SetLineColor(colors[i])
    if i == 0:
        resList[i].GetXaxis().SetTitle('M_{e#mu} (GeV)')
        resList[i].GetXaxis().SetTitleFont(font)
        resList[i].GetXaxis().SetLabelFont(font)
        resList[i].GetYaxis().SetTitle('#sigma(M_{e#mu})/M_{e#mu}')
        resList[i].GetYaxis().SetTitleOffset(1.55)
        resList[i].GetYaxis().SetTitleFont(font)
        resList[i].GetYaxis().SetLabelFont(font)
        resList[i].Draw('ap')
    else:
        resList[i].Draw('psame')

    resFitFuncList.append(TF1('resFitFunc{0}'.format(i), 'pol2', 150, 5000))
    resFitFuncList[i].SetLineColor(colors[i])
    resList[i].Fit(resFitFuncList[i], '', '', 150, 5000)

    # move the stats box
    gPad.Update()
    statsBoxList.append(resList[i].GetListOfFunctions().FindObject('stats'))
    statsBoxList[i].SetX1NDC(statsBoxCoord[i][0])
    statsBoxList[i].SetY1NDC(statsBoxCoord[i][1])
    statsBoxList[i].SetX2NDC(statsBoxCoord[i][2])
    statsBoxList[i].SetY2NDC(statsBoxCoord[i][3])
    statsBoxList[i].SetLineColor(colors[i])
    statsBoxList[i].SetTextColor(colors[i])

tex2 = TLatex()
tex2.SetNDC()
tex2.SetTextFont(font)
tex2.SetTextSize(0.04)
tex2.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')

## old AN-13-422 resolution function
#anFitFunc = TF1('anFitFunc', 'pol2', 150, 5000)
#anFitFunc.SetParameters(0.013, 1.8e-5, -8.9e-10)
#anFitFunc.Draw('psame')

legend = TLegend(0.18, 0.52, 0.47, 0.65)
legend.SetTextFont(font)
legend.SetTextSize(0.03)
legend.SetBorderSize(0)
legend.SetFillColor(19)
legend.SetFillStyle(0)
legend.AddEntry(resList[0], 'START53_V7C1', 'lp')
legend.AddEntry(resList[1], 'START53_V7C2', 'lp')
#legend.AddEntry(anFitFunc, 'AN-13-422', 'l')
legend.Draw('same')
 
c1.Modified()
c1.Update()

##############################################################################
# if there are two functions extract the error from the difference
if len(treePrefices) > 1:
    c2 = TCanvas('resError', 'resError', 100, 20, 600, 600)
    gPad.SetLeftMargin(0.14)
    gPad.SetRightMargin(0.06)
    gPad.SetTickx(1)
    gPad.SetTicky(1)
    
    resErr = TGraphErrors(resList[0].GetN())
    g1x = resList[0].GetX()
    g1y = resList[0].GetY()
    g2y = resList[1].GetY()
    for i in range(resList[0].GetN()):
        num = abs(g1y[i]-g2y[i])
        resErr.SetPoint(i, g1x[i], num/g1y[i])
        numErr2 = resList[0].GetErrorY(i)**2 + resList[1].GetErrorY(i)**2
        resErr.SetPointError(i, 0., num/g1y[i] * math.sqrt((resList[0].GetErrorY(i)/g1y[i])**2 + numErr2/num**2)) # error propagation for |A-B|/A
    resErr.GetXaxis().SetTitle('M_{e#mu} (GeV)')
    resErr.GetXaxis().SetTitleFont(font)
    resErr.GetXaxis().SetLabelFont(font)
    resErr.GetYaxis().SetTitle('|#sigma(M_{e#mu})_{V7C1}-#sigma(M_{e#mu})_{V7C2}|/#sigma(M_{e#mu})_{V7C1}')
    resErr.GetYaxis().SetTitleOffset(1.55)
    resErr.GetYaxis().SetTitleFont(font)
    resErr.GetYaxis().SetLabelFont(font)
    resErr.Draw('ap')

    resErrFitFuncLow = TF1('resErrFitFuncLow', 'pol0', 150, 700)
    resErrFitFuncLow.SetLineColor(colors[0])
    resErr.Fit(resErrFitFuncLow, '', '', 150, 700)

    ## relative difference of the two reolution fit functions drawn on plot
    #resErrFunc = TF1('resErrFunc', 'pol2 / pol2(3)', 150, 5000)
    #resErrFunc.SetParameters(resFitFuncList[0].GetParameter(0)-resFitFuncList[1].GetParameter(0), resFitFuncList[0].GetParameter(1)-resFitFuncList[1].GetParameter(1), resFitFuncList[0].GetParameter(2)-resFitFuncList[1].GetParameter(2), resFitFuncList[0].GetParameter(0), resFitFuncList[0].GetParameter(1), resFitFuncList[0].GetParameter(2))
    #resErrFunc.SetLineColor(colors[1])
    #resErrFunc.Draw('lsame')
    #resErrFunc.Print()

    tex2.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')
    # move the stats box
    gPad.Update()
    statsBoxList.append(resErr.GetListOfFunctions().FindObject('stats').Clone('statsPol0'))
    statsBoxList[len(statsBoxList)-1].SetX1NDC(statsBoxCoord[2][0])
    statsBoxList[len(statsBoxList)-1].SetY1NDC(statsBoxCoord[2][1])
    statsBoxList[len(statsBoxList)-1].SetX2NDC(statsBoxCoord[2][2])
    statsBoxList[len(statsBoxList)-1].SetY2NDC(statsBoxCoord[2][3])
    statsBoxList[len(statsBoxList)-1].SetLineColor(colors[0])
    statsBoxList[len(statsBoxList)-1].SetTextColor(colors[0])

    resErrFitFuncHigh = TF1('resErrFitFuncHigh', 'pol2', 600, 5000)
    resErrFitFuncHigh.SetLineColor(colors[0])
    resErr.Fit(resErrFitFuncHigh, '', '', 600, 5000)

    gPad.Update()
    statsBoxList.append(resErr.GetListOfFunctions().FindObject('stats'))
    statsBoxList[len(statsBoxList)-1].SetX1NDC(statsBoxCoord[1][0])
    statsBoxList[len(statsBoxList)-1].SetY1NDC(statsBoxCoord[1][1])
    statsBoxList[len(statsBoxList)-1].SetX2NDC(statsBoxCoord[1][2])
    statsBoxList[len(statsBoxList)-1].SetY2NDC(statsBoxCoord[1][3])
    statsBoxList[len(statsBoxList)-1].SetLineColor(colors[0])
    statsBoxList[len(statsBoxList)-1].SetTextColor(colors[0])

    resErrFitFuncLow.Draw('lsame')
    statsBoxList[len(statsBoxList)-2].Draw('same')

    c2.Modified()
    c2.Update()

##############################################################################
# save canvases to root file
if savePlots:
    #output = TFile('./plots/resolutionPlots.root', 'recreate')
    output = TFile('./plots/resolutionPlots_singleMu.root', 'recreate')
    plotDir = './plots/plots/'
    output.cd()
    for canvas in cList:
        canvas.Write(canvas.GetName())
        canvas.Print(plotDir+canvas.GetName()+'.png', 'png')
        canvas.Print(plotDir+canvas.GetName()+'.pdf', 'pdf')
    c1.Write(c1.GetName())
    c1.Print(plotDir+c1.GetName()+'.png', 'png')
    c1.Print(plotDir+c1.GetName()+'.pdf', 'pdf')
    for func in resFitFuncList:
        func.Write()
    #anFitFunc.Write()
    if len(treePrefices) > 1:
        c2.Write(c2.GetName())
        c2.Print(plotDir+c2.GetName()+'.png', 'png')
        c2.Print(plotDir+c2.GetName()+'.pdf', 'pdf')
        resErrFitFuncLow.Write()
        resErrFitFuncHigh.Write()
        #resErrFunc.Write()
    output.Close()

