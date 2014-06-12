#!/usr/bin/env python

import ROOT
import math
ROOT.gROOT.SetBatch(True)
from ROOT import TFile,TTree,TH1F,TH1D,TF1,TGraph,TGraphErrors
from ROOT import TCanvas,TLegend,TLatex,TPaveStats, TList
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
#input = TFile("emuSpec_singleMuTrg_topxsect245p8_19706pb-1.root")
treePrefices = ['emuTree_sigNoAccCuts', 'emuTree_sigTagV7C2_']
#treePrefices = ['emuTree_sigTagV7C2_']
recoTags = ['START53_V7C1', 'START53_V7C2']

muPt_ranges = [75., 150., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1125., 1250., 1375., 1500., 1625., 1750., 2000., 2400.]

#get the trees from the root files
input.cd()

# list with all canvases
cList = []

keyList = input.GetListOfKeys()

resList = []

hPtDiffs = []
helperCanvas = TCanvas('helper', 'helper', 100, 20, 600, 600)
for i, treePrefix in enumerate(treePrefices):
    resList.append(TGraphErrors(len(muPt_ranges)))
    pt = 0
    j = 0
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

            lowerPt = 0.
            for k, upperPt in enumerate(muPt_ranges):
                #brName = '(muPt-genMuPt)/genMuPt'
                brName = '(genMuPt-muPt)/muPt'
                #brName = '(genMuPt*muCharge-muPt*genMuCharge)/muPt'
                #brName = 'genMuPt'

                axisRange = 0.2
                if lowerPt > 500.:
                    axisRange = 0.5
                if lowerPt > 1500.:
                    axisRange = 0.8
                #get the histograms from the tree
                helperCanvas.cd()
                tree.Draw(brName+'>>{0}_hPtDiff_{1}_{2}(100, -{3}, {3})'.format(treePrefices[i], keyName, k, axisRange), 'genMuPt>={0} && genMuPt<{1}'.format(lowerPt, upperPt))
                hPtDiff_part = gROOT.FindObject('{0}_hPtDiff_{1}_{2}'.format(treePrefices[i], keyName, k))

                if j == 0:
                    hPtDiffs.append(hPtDiff_part)
                else:
                    l = k + i*len(muPt_ranges)
                    hPtDiffs[l].Add(hPtDiff_part)

                lowerPt = upperPt
            j += 1

    lowerPt = 35.
    for k, upperPt in enumerate(muPt_ranges):
        l = k + i*len(muPt_ranges)
        cList.append(TCanvas(treePrefices[i]+'_res_genMuPt_%.0f_to_%.0fGeV'%(lowerPt, upperPt), treePrefices[i]+'_res_genMuPt_%.0f_to_%.0fGeV'%(lowerPt, upperPt), 100, 20, 600, 600))

        print 'Total events pT range [{0}, {1}]: {2}'.format(lowerPt, upperPt, hPtDiffs[l].GetEntries()) 
        print l
        #normalise
        hPtDiffs[l].Sumw2()
        #hPtDiffs[l].Scale(1./hPtDiffs[l].Integral())

        #now draw
        gPad.SetLeftMargin(0.14)
        gPad.SetRightMargin(0.06)
        gPad.SetTickx(1)
        gPad.SetTicky(1)

        hPtDiffs[l].SetMarkerStyle(markers[i])
        hPtDiffs[l].SetMarkerColor(colors[0])
        #hPtDiffs[l].SetLineWidth(2)
        hPtDiffs[l].SetLineColor(colors[0])
        hPtDiffs[l].Draw('e0')
        #hPtDiffs[l].GetXaxis().SetTitle('(p^{#mu reco}_{T} - p^{#mu gen}_{T})/p^{#mu gen}_{T}')
        hPtDiffs[l].GetXaxis().SetTitle('(p^{#mu gen}_{T} - p^{#mu reco}_{T})/p^{#mu reco}_{T}')
        hPtDiffs[l].GetXaxis().SetTitleOffset(1.05)
        hPtDiffs[l].GetXaxis().SetTitleFont(font)
        hPtDiffs[l].GetXaxis().SetLabelFont(font)
        hPtDiffs[l].GetYaxis().SetTitle('Events')
        hPtDiffs[l].GetYaxis().SetTitleOffset(1.5)
        hPtDiffs[l].GetYaxis().SetTitleFont(font)
        hPtDiffs[l].GetYaxis().SetLabelFont(font)

        gauss = TF1('gauss', 'gaus')
        gauss.SetLineColor(colors[1])
        peak = hPtDiffs[l].GetBinCenter(hPtDiffs[l].GetMaximumBin())
        gauss.SetParameter(1, peak)
        gauss.SetParameter(2, hPtDiffs[l].GetRMS())
        factorFit1 = 1.
        factorFit2 = 2.
        if hPtDiffs[l].GetEntries() < 400:
            factorFit2 = 3. 
        # raw and fine fit
        hPtDiffs[l].Fit(gauss, '', '', peak-factorFit1*hPtDiffs[l].GetRMS(), peak+factorFit1*hPtDiffs[l].GetRMS())
        hPtDiffs[l].Fit(gauss, '', '', gauss.GetParameter(1)-factorFit2*gauss.GetParameter(2), gauss.GetParameter(1)+factorFit2*gauss.GetParameter(2))
    
        #hPtDiffs[l].GetXaxis().SetRangeUser(math.floor(gauss.GetParameter(1)-4.*gauss.GetParameter(2)), math.floor(gauss.GetParameter(1)+4.*gauss.GetParameter(2)))
        #if massPt < 1000:
        #    hPtDiffs[l].GetXaxis().SetRangeUser(-0.1, 0.1)
        #elif massPt < 2000:
        #    hPtDiffs[l].GetXaxis().SetRangeUser(-0.2, 0.2)
        #else:
        #    hPtDiffs[l].GetXaxis().SetRangeUser(-0.4, 0.4)

        resList[i].SetPoint(k, (upperPt+lowerPt)/2., gauss.GetParameter(2))
        resList[i].SetPointError(k, (upperPt-lowerPt)/2., gauss.GetParError(2))

        tex = TLatex()
        tex.SetNDC()
        tex.SetTextFont(font)
        tex.SetTextSize(0.04)
        tex.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')
        tex.SetTextSize(0.03)
        tex.DrawLatex(0.18, 0.85, "p_{T}^{#mu gen} = ["+str(lowerPt)+', '+str(upperPt)+'] GeV')
        tex.DrawLatex(0.18, 0.80, recoTags[i])
        cList[l].Modified()
        cList[l].Update()

        lowerPt = upperPt


##############################################################################
# plot the resolution
c1 = TCanvas('mu_pt_resolution', 'mu_pt_resolution', 100, 20, 600, 600)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.06)
gPad.SetTickx(1)
gPad.SetTicky(1)

resFitFuncList = []
statsBoxList = []
statsBoxCoord = [[0.18, 0.66, 0.54, 0.86], [0.54, 0.14, 0.90, 0.34]]

for i, treePrefix in enumerate(treePrefices):
    resList[i].SetMarkerStyle(markers[i])
    resList[i].SetMarkerColor(colors[i])
    resList[i].SetLineColor(colors[i])
    if i == 0:
        resList[i].GetXaxis().SetTitle('p^{#mu gen}_{T} (GeV)')
        resList[i].GetXaxis().SetTitleFont(font)
        resList[i].GetXaxis().SetLabelFont(font)
        resList[i].GetYaxis().SetTitle('#sigma(p^{#mu}_{T})/p^{#mu}_{T}')
        resList[i].GetYaxis().SetTitleOffset(1.55)
        resList[i].GetYaxis().SetTitleFont(font)
        resList[i].GetYaxis().SetLabelFont(font)
        resList[i].Draw('ap')
    else:
        resList[i].Draw('psame')

    #resFitFuncList.append(TF1('resFitFunc{0}'.format(i), 'pol2', 0, 3000))
    resFitFuncList.append(TF1('resFitFunc{0}'.format(i), '[0]+[1]*(1-exp(-x/[2]))', 0, 3000))
    resFitFuncList[i].SetParameters(0.01, 0.22, 2300.)
    resFitFuncList[i].SetLineColor(colors[i])
    resList[i].Fit(resFitFuncList[i], '', '', 0, 3000)

    # move the stats box
    gPad.Update()
    statsBoxList.append(resList[i].GetListOfFunctions().FindObject('stats'))
    statsBoxList[i].SetX1NDC(statsBoxCoord[i][0])
    statsBoxList[i].SetY1NDC(statsBoxCoord[i][1])
    statsBoxList[i].SetX2NDC(statsBoxCoord[i][2])
    statsBoxList[i].SetY2NDC(statsBoxCoord[i][3])
    statsBoxList[i].SetLineColor(colors[i])
    statsBoxList[i].SetTextColor(colors[i])

#function = TF1('function', '[0]+[1]*(1-exp(-x/[2]))', 0, 3000)
#function.SetParameters(1.35863e-2, 5.54127e-1, 6.23052e+3)
#function.SetLineStyle(ROOT.kDashed)
#function.SetLineColor(colors[2])
#function.Draw('lsame')

tex2 = TLatex()
tex2.SetNDC()
tex2.SetTextFont(font)
tex2.SetTextSize(0.04)
tex2.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')
tex2.SetTextSize(0.03)
tex2.DrawLatex(0.20, 0.16, 'f(x|p_{0},p_{1},p_{2}) = p_{0}+p_{1}#left(1-e^{-x/p_{2}}#right)')

legend = TLegend(0.18, 0.56, 0.47, 0.65)
legend.SetTextFont(font)
legend.SetTextSize(0.03)
legend.SetBorderSize(0)
legend.SetFillColor(19)
legend.SetFillStyle(0)
for i, res in enumerate(resList):
    legend.AddEntry(res, recoTags[i], 'lp')
legend.Draw('same')
 
c1.Modified()
c1.Update()

##############################################################################
# if there are two functions extract the error from the difference
if len(treePrefices) > 1:
    c2 = TCanvas('mu_pt_resError', 'mu_pt_resError', 100, 20, 600, 600)
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
        resErr.SetPointError(i, resList[0].GetErrorX(i), num/g1y[i] * math.sqrt((resList[0].GetErrorY(i)/g1y[i])**2 + numErr2/num**2)) # error propagation for |A-B|/A
    resErr.GetXaxis().SetTitle('p^{#mu gen}_{T} (GeV)')
    resErr.GetXaxis().SetTitleFont(font)
    resErr.GetXaxis().SetLabelFont(font)
    resErr.GetYaxis().SetTitle('|#sigma(p^{#mu V7C1}_{T})-#sigma(p^{#mu V7C2}_{T})|/#sigma(p^{#mu V7C1}_{T})')
    resErr.GetYaxis().SetTitleOffset(1.55)
    resErr.GetYaxis().SetTitleFont(font)
    resErr.GetYaxis().SetLabelFont(font)
    resErr.Draw('ap')

    resErrFitFunc = TF1('resErrFitFunc', 'pol2', 35, 2800)
    resErrFitFunc.SetParameters(-0.086, 5.e-4, -1.2e-7)
    resErrFitFunc.SetLineColor(colors[0])
    resErr.Fit(resErrFitFunc, '', '', 35, 2800)

    #resErrFunc = TF1('resErrFunc', 'pol2 / pol2(3)', 250, 5000)
    #resErrFunc.SetParameters(resFitFuncList[0].GetParameter(0)-resFitFuncList[1].GetParameter(0), resFitFuncList[0].GetParameter(1)-resFitFuncList[1].GetParameter(1), resFitFuncList[0].GetParameter(2)-resFitFuncList[1].GetParameter(2), resFitFuncList[0].GetParameter(0), resFitFuncList[0].GetParameter(1), resFitFuncList[0].GetParameter(2))
    #resErrFunc.SetLineColor(colors[1])
    #resErrFunc.Draw('lsame')
    #resErrFunc.Print()

    tex2.DrawLatex(0.14, 0.91, 'CMS Simulation, 8 TeV')
    # move the stats box
    gPad.Update()
    statsBoxList.append(resErr.GetListOfFunctions().FindObject('stats'))
    statsBoxList[len(statsBoxList)-1].SetX1NDC(statsBoxCoord[1][0])
    statsBoxList[len(statsBoxList)-1].SetY1NDC(statsBoxCoord[1][1])
    statsBoxList[len(statsBoxList)-1].SetX2NDC(statsBoxCoord[1][2])
    statsBoxList[len(statsBoxList)-1].SetY2NDC(statsBoxCoord[1][3])
    statsBoxList[len(statsBoxList)-1].SetLineColor(colors[0])
    statsBoxList[len(statsBoxList)-1].SetTextColor(colors[0])

    c2.Modified()
    c2.Update()

##############################################################################
# save canvases to root file
if savePlots:
    #output = TFile('./plots/muPtResolutionPlots.root', 'recreate')
    output = TFile('./plots/muPtResolutionPlots_singleMu.root', 'recreate')
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
    if len(treePrefices) > 1:
        c2.Write(c2.GetName())
        c2.Print(plotDir+c2.GetName()+'.png', 'png')
        c2.Print(plotDir+c2.GetName()+'.pdf', 'pdf')
        resErrFitFunc.Write()
        #resErrFunc.Write()
    output.Close()

