#!/usr/bin/env python

import os
import shutil
import math
import array
import numpy
import ROOT
#ROOT.gROOT.SetBatch(True)
from ROOT import TMath,TFile,TH1F,TH1D,TGraph,TGraphAsymmErrors,TF1
from ROOT import TCanvas,TLegend,TLatex
from ROOT import gROOT,gStyle,gPad
from optparse import OptionParser

# set the propper axis style ans labels
def setAxisLabels(graph, min, max, inFb, fontStyle):
    graph.GetXaxis().SetTitle("M_{e#mu} (GeV)")
    #graph.GetXaxis().SetTitleOffset(1.)
    graph.GetXaxis().SetTitleFont(fontStyle)
    graph.GetXaxis().SetLabelFont(fontStyle)
    graph.GetXaxis().SetRangeUser(min, max)
    if inFb:
        graph.GetYaxis().SetTitle('#sigma x BR (fb)')
    else:
        graph.GetYaxis().SetTitle('#sigma x BR (pb)')
    #graph.GetYaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetTitleFont(fontStyle)
    graph.GetYaxis().SetLabelFont(fontStyle)

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
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# command line options
parser = OptionParser(usage="usage: %prog [options]", description="Plot limits.\n")
parser.add_option("-w" ,"--workDir1", dest="workDir1", default=os.getcwd()+'/', help="Path to working directory one.")
parser.add_option("-W" ,"--workDir2", dest="workDir2", default=os.getcwd()+'/', help="Path to working directory two.")
parser.add_option("-n" ,"--name1", dest="name1", default='', help="Display name for legend one.")
parser.add_option("-N" ,"--name2", dest="name2", default='', help="Display name for legend two.")
parser.add_option("-o" ,"--obs1", dest="drawObs1", default=False, action="store_true", help="Draw observed limit one.")
parser.add_option("-O" ,"--obs2", dest="drawObs2", default=False, action="store_true", help="Draw observed limit two.")
parser.add_option("-e" ,"--exp1", dest="drawExp1", default=False, action="store_true", help="Draw expected limit one.")
parser.add_option("-E" ,"--exp2", dest="drawExp2", default=False, action="store_true", help="Draw expected limit two.")
parser.add_option("-1" ,"--1sigma1", dest="draw1sigma1", default=False, action="store_true", help="Draw one sigma band one.")
parser.add_option("-2" ,"--2sigma1", dest="draw2sigma1", default=False, action="store_true", help="Draw two sigma band one.")
parser.add_option("-3" ,"--3sigma1", dest="draw3sigma1", default=False, action="store_true", help="Draw three sigma band one.")
parser.add_option("-4" ,"--1sigma2", dest="draw1sigma2", default=False, action="store_true", help="Draw one sigma band two.")
parser.add_option("-5" ,"--2sigma2", dest="draw2sigma2", default=False, action="store_true", help="Draw two sigma band two.")
parser.add_option("-6" ,"--3sigma2", dest="draw3sigma2", default=False, action="store_true", help="Draw three sigma band two.")
parser.add_option("-t" ,"--theory1", dest="drawTheory1", default=False, action="store_true", help="Draw theory curve one.")
parser.add_option("-T" ,"--theory2", dest="drawTheory2", default=False, action="store_true", help="Draw theory curve one.")
parser.add_option("-s" ,"--save", dest="savePlots", default=False, action="store_true", help="Save plots to root file in results directory.")
parser.add_option("-d"  ,"--dots1", dest="showMarker1", default=False, action="store_true", help="Plot a marker at the calculated mass point for limit one.")
parser.add_option("-D"  ,"--dots2", dest="showMarker2", default=False, action="store_true", help="Plot a marker at the calculated mass point for limit two.")
(options, args) = parser.parse_args()

if not options.workDir1[-1:] == '/':
    options.workDir1 = options.workDir1+'/'
if not options.workDir2[-1:] == '/':
    options.workDir2 = options.workDir2+'/'

# graphs and curves to plot
if not os.path.isfile(options.workDir1+'results/limitPlot.root'):
    raise Exception('File {} does not exist.'.format(options.workDir1+'results/limitPlot.root'))
file1 = TFile(options.workDir1+'results/limitPlot.root')
obs_lim1 = file1.Get('obs_lim')
median_3sigma1 = file1.Get('median_3sigma')
median_2sigma1 = file1.Get('median_2sigma')
median_1sigma1 = file1.Get('median_1sigma')
theory_curve1 = file1.Get('theory_curve')
obs_lim1.SetName('obs_lim1')
median_3sigma1.SetName('median_3sigma1')
median_2sigma1.SetName('median_2sigma1')
median_1sigma1.SetName('median_1sigma1')
theory_curve1.SetName('theory_curve1')

if not os.path.isfile(options.workDir2+'results/limitPlot.root'):
    raise Exception('File {} does not exist.'.format(options.workDir2+'results/limitPlot.root'))
file2 = TFile(options.workDir2+'results/limitPlot.root')
obs_lim2 = file2.Get('obs_lim')
median_3sigma2 = file2.Get('median_3sigma')
median_2sigma2 = file2.Get('median_2sigma')
median_1sigma2 = file2.Get('median_1sigma')
theory_curve2 = file2.Get('theory_curve')
obs_lim2.SetName('obs_lim2')
median_3sigma2.SetName('median_3sigma2')
median_2sigma2.SetName('median_2sigma2')
median_1sigma2.SetName('median_1sigma2')
theory_curve2.SetName('theory_curve2')

# clone the median graph for plotting the expected limit
median1 = median_1sigma1.Clone('median1')
median1.SetLineWidth(2)
median1.SetLineColor(ROOT.kBlue)
median1.SetLineStyle(ROOT.kDashed)
median2 = median_1sigma2.Clone('median2')
median2.SetLineWidth(2)
median2.SetLineColor(ROOT.kRed)
median2.SetLineStyle(ROOT.kDotted)

##############################################################################
#now draw the limits
c = TCanvas('limitComparison')
c.SetLogy(True)

optionAxis = 'A'
#draw expected limits
rd1 = gROOT.GetColor(ROOT.kRed-3)
rd1.SetAlpha(0.5)
median_3sigma1.SetFillColor(ROOT.kRed-3)
median_3sigma1.SetLineColor(ROOT.kRed-3)
axisHisto = median_3sigma1.GetHistogram()
if options.draw3sigma1:
    median_3sigma1.Draw(optionAxis+'3')
    optionAxis = ''
rd2 = gROOT.GetColor(ROOT.kRed-6)
rd2.SetAlpha(0.5)
median_3sigma2.SetFillColor(ROOT.kRed-6)
median_3sigma2.SetLineColor(ROOT.kRed-6)
if options.draw3sigma2:
    median_3sigma2.Draw(optionAxis+'3')
    optionAxis = ''

ye1 = gROOT.GetColor(ROOT.kYellow+1)
ye1.SetAlpha(0.5)
median_2sigma1.SetFillColor(ROOT.kYellow+1)
median_2sigma1.SetLineColor(ROOT.kYellow+1)
if options.draw2sigma1:
    median_2sigma1.Draw(optionAxis+'3')
    optionAxis = ''
ye2 = gROOT.GetColor(ROOT.kYellow)
ye2.SetAlpha(0.5)
median_2sigma2.SetFillColor(ROOT.kYellow)
median_2sigma2.SetLineColor(ROOT.kYellow)
if options.draw2sigma2:
    median_2sigma2.Draw(optionAxis+'3')
    optionAxis = ''

gr1 = gROOT.GetColor(ROOT.kGreen+1)
#gr1.SetAlpha(0.5)
median_1sigma1.SetFillColor(ROOT.kGreen+1)
median_1sigma1.SetLineColor(ROOT.kGreen+1)
if options.draw1sigma1:
    median_1sigma1.Draw(optionAxis+'3')
    optionAxis = ''
gr2 = gROOT.GetColor(ROOT.kGreen)
gr2.SetAlpha(0.5)
median_1sigma2.SetFillColor(ROOT.kGreen)
median_1sigma2.SetLineColor(ROOT.kGreen)
if options.draw1sigma2:
    median_1sigma2.Draw(optionAxis+'3')
    optionAxis = ''

if options.drawExp1:
    median1.Draw(optionAxis+'LX')
    optionAxis = ''
if options.drawExp2:
    median2.Draw(optionAxis+'LX')
    optionAxis = ''

#draw observed limits
obs_lim1.SetLineWidth(2)
obs_lim1.SetMarkerStyle(20)
drawStyle1 = 'l'
if options.drawObs1:
    if options.showMarker1:
        drawStyle1 += 'p'
    obs_lim1.Draw(optionAxis+drawStyle1)
    optionAxis = ''
obs_lim2.SetLineWidth(2)
obs_lim2.SetLineColor(ROOT.kRed)
obs_lim2.SetMarkerStyle(ROOT.kCircle)
obs_lim2.SetMarkerColor(ROOT.kRed)
drawStyle2 = 'l'
if options.drawObs2:
    if options.showMarker2:
        drawStyle2 += 'p'
    obs_lim2.Draw(optionAxis+drawStyle2)
    optionAxis = ''

#draw theory curve
theory_curve1.SetLineWidth(2)
theory_curve1.SetLineColor(ROOT.kCyan+1)
if options.drawTheory1:
    if optionAxis == '':
        theory_curve1.Draw('Lsame')
    else:
        theory_curve1.Draw('L')
theory_curve2.SetLineWidth(2)
theory_curve2.SetLineColor(ROOT.kCyan-1)
theory_curve2.SetLineStyle(ROOT.kDashed)
if options.drawTheory2:
    if optionAxis == '':
        theory_curve2.Draw('Lsame')
    else:
        theory_curve2.Draw('L')

#legend
legend = TLegend(0.32, 0.50, 0.89, 0.88)
legend.SetTextFont(font)
legend.SetTextSize(0.04)
legend.SetBorderSize(0)
legend.SetFillColor(19)
legend.SetFillStyle(0)
legend.SetNColumns(2)
leftCol = True
if options.name1 != '':
    if not leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(options.name1, options.name1, '')
    leftCol = False
if options.name2 != '':
    if leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(options.name2, options.name2, '')
    leftCol = True
if options.drawObs1:
    if not leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(obs_lim1, '95% CL limit', drawStyle1)
    leftCol = False
if options.drawObs2:
    if leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(obs_lim2, '95% CL limit', drawStyle2)
    leftCol = True
if options.drawExp1:
    if not leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median1, 'median exp. limit', 'l')
    leftCol = False
if options.drawExp2:
    if leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median2, 'median exp. limit', 'l')
    leftCol = True
if options.draw1sigma1:
    if not leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median_1sigma1, '68% expected', 'f')
    leftCol = False
if options.draw1sigma2:
    if leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median_1sigma2, '68% expected', 'f')
    leftCol = True
if options.draw2sigma1:
    if not leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median_2sigma1, '95% expected', 'f')
    leftCol = False
if options.draw2sigma2:
    if leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median_2sigma2, '95% expected', 'f')
    leftCol = True
if options.draw3sigma1:
    if not leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median_3sigma1, '99% expected', 'f')
    leftCol = False
if options.draw3sigma2:
    if leftCol:
        legend.AddEntry('empty', '', '')
    legend.AddEntry(median_3sigma2, '99% expected', 'f')
    leftCol = True
if options.drawTheory1:
    if not leftCol:
        legend.AddEntry('empty', '', '')
    sigText = "Z' signal (LO)"
    #sigText = "#splitline{{Z' signal (LO)}}{{#kappa = {:.1g} #times M_{{Z'}} / 100 TeV}}".format(options.kappa)
    legend.AddEntry(theory_curve1, sigText, 'l')
    leftCol = False
if options.drawTheory2:
    if leftCol:
        legend.AddEntry('empty', '', '')
    sigText = "Z' signal (LO)"
    #sigText = "#splitline{{Z' signal (LO)}}{{#kappa = {:.1g} #times M_{{Z'}} / 100 TeV}}".format(options.kappa)
    legend.AddEntry(theory_curve2, sigText, 'l')
    leftCol = True
legend.Draw('same')

tex = TLatex()
tex.SetNDC()
tex.SetTextFont(font)
tex.SetTextSize(0.04)
tex.DrawLatex(0.1, 0.91, 'CMS Preliminary')
tex.DrawLatex(0.717, 0.91, '19.7 fb^{-1} (8 TeV)')

axisHisto.Draw("sameaxis")
##############################################################################
# save canvases to root file
if options.savePlots:
    output = TFile(options.workDir1+'results/limitComparisonPlot.root', 'recreate')
    c.Write(c.GetName())
    output.Close()
    print 'Limit comparion plot written to: '+options.workDir1+'results/limitComparisonPlot.root'

# wait
raw_input("Press ENTER to quit.")

