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
parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
parser.add_option("-m"  ,"--minMass", dest="minMass", default="100", type="float", help="Minimal resonance mass.")
parser.add_option("-M"  ,"--maxMass", dest="maxMass", default="100", type="float", help="Maximal resonance mass.")
parser.add_option("-S"  ,"--massStep", dest="massStep", default="50", type="float", help="Mass step size.")
parser.add_option("-k"  ,"--kappa", dest="kappa", default="1.", type="float", help="Kappa factor for coupling.")
parser.add_option("-o"  ,"--obs", dest="drawObs", default=False, action="store_true", help="Draw observed limit.")
parser.add_option("-e"  ,"--exp", dest="drawExp", default=False, action="store_true", help="Draw expected limit.")
parser.add_option("-1"  ,"--1sigma", dest="draw1sigma", default=False, action="store_true", help="Draw one sigma band.")
parser.add_option("-2"  ,"--2sigma", dest="draw2sigma", default=False, action="store_true", help="Draw two sigma band.")
parser.add_option("-3"  ,"--3sigma", dest="draw3sigma", default=False, action="store_true", help="Draw three sigma band.")
parser.add_option("-t"  ,"--theory", dest="drawTheory", default=False, action="store_true", help="Draw theory curve.")
parser.add_option("-f"  ,"--fb", dest="inFemtoBarn", default=False, action="store_true", help="Output cross section times BR in fb instead of pb.")
parser.add_option("-s"  ,"--save", dest="savePlots", default=False, action="store_true", help="Save plots to root file in results directory.")
parser.add_option("-d"  ,"--dots", dest="showMarker", default=False, action="store_true", help="Plot a marker at the calculated mass point.")
(options, args) = parser.parse_args()

if not options.workDir[-1:] == '/':
    options.workDir = options.workDir+'/'

# mass points to draw
massPoints = numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep)

# theory cross section curve
xsec_x_br_sig = TF1('xsec_x_br_sig', '[0] * exp(-x/[1] + [2]/x - x**2/[3])', 0.5*options.minMass, 1.2*options.maxMass)
xsec_x_br_sig.SetParameters(3.98360e-4, 397.414, 323.610, 1.05017e7)

# factor in case the y axis should be in fb
femtoFact = 1.
if options.inFemtoBarn:
    femtoFact = 1000.

# graphs and curves to plot
obs_lim = TGraph(len(massPoints))
obs_lim.SetName('obs_lim')

median_3sigma = TGraphAsymmErrors(len(massPoints))
median_3sigma.SetName('median_3sigma')
median_2sigma = median_3sigma.Clone('median_2sigma')
median_1sigma = median_3sigma.Clone('median_1sigma')
theory_curve = xsec_x_br_sig.Clone('theory_curve')
theory_curve.SetParameter(0, options.kappa**2 * xsec_x_br_sig.GetParameter(0)*femtoFact)

print '----------------------------------------------------------------------------------------------'
print '  Mass      |   Observed limit  |         Expected limit (+/- 1 sigma)        |  Theory'
print '----------------------------------------------------------------------------------------------'
# loop over the mass points
oi = -1
ei = -1
for mass in massPoints:
    scale = options.kappa**2 * xsec_x_br_sig.Eval(mass) * femtoFact

    #get the observed limit from the root file
    resultsFile = TFile(options.workDir+'results/higgsCombineObserved.MarkovChainMC.mH120_m{:.0f}.root'.format(mass))
    tree = resultsFile.Get('limit')
    obs_lim_value = 0.
    if tree == None:
        pass
        #print "No data found in file '{0}'. Will set observed limit to 0 for mass={1} GeV.".format(resultsFile.GetName(), mass)
    else:
        if tree.GetEntry(0) > 0:
            obs_lim_value = tree.limit * scale
        else:
            pass
            #print "No limit found in file '{0}'. Will set observed limit to 0 for mass={1} GeV.".format(resultsFile.GetName(), mass)
    if obs_lim_value > 0.:
        oi += 1
        obs_lim.SetPoint(oi, mass, obs_lim_value)

    #get the 1-, 2-, 3-sigma quantiles and the median from the toys for the expected limits
    expectedFile = TFile(options.workDir+'results/higgsCombineExpected.MarkovChainMC.mH120_m{:.0f}.root'.format(mass))
    expTree = expectedFile.Get('limit')
    if expTree == None: 
        #print "No data found in file '{0}'. Will set expected limit to 0 for mass={1} GeV.".format(resultsFile.GetName(), mass)
        quantiles = array.array('d', [0., 0., 0., 0., 0., 0., 0.])
    else:
        expLimit = array.array('d')
        for row in expTree:
            expLimit.append(row.limit)

        #calculate quantiles
        prob = array.array('d', [0.0013, 0.0228, 0.1587, 0.5, 0.8413, 0.9772, 0.9987])
        quantiles = array.array('d', [0., 0., 0., 0., 0., 0., 0.])
        TMath.Quantiles(expTree.GetEntries(), 7, expLimit, quantiles, prob, False)

    # fill the graphs
    if quantiles[3] > 0:
        ei += 1
        median_3sigma.SetPoint(ei, mass, quantiles[3] * scale)
        median_3sigma.SetPointError(ei, 0., 0., (quantiles[3]-quantiles[0]) * scale, (quantiles[6]-quantiles[3]) * scale)
        median_2sigma.SetPoint(ei, mass, quantiles[3] * scale)
        median_2sigma.SetPointError(ei, 0., 0., (quantiles[3]-quantiles[1]) * scale, (quantiles[5]-quantiles[3]) * scale)
        median_1sigma.SetPoint(ei, mass, quantiles[3] * scale)
        median_1sigma.SetPointError(ei, 0., 0., (quantiles[3]-quantiles[2]) * scale, (quantiles[4]-quantiles[3]) * scale)

    # text output of the limits
    if options.inFemtoBarn:
        print "{:>6.1f} GeV  |  {:>12.5f} fb  |  {:>12.5f} ({:>+11.5g} {:>+11.5g}) fb  |  {:g} fb".format(mass, obs_lim_value, quantiles[3] * scale, (quantiles[4]-quantiles[3]) * scale, (quantiles[2]-quantiles[3]) * scale, theory_curve.Eval(mass))
    else:
        print "{:>6.1f} GeV  |  {:>12.5f} pb  |  {:>12.5f} ({:>+11.5g} {:>+11.5g}) pb  |  {:g} pb".format(mass, obs_lim_value, quantiles[3] * scale, (quantiles[4]-quantiles[3]) * scale, (quantiles[2]-quantiles[3]) * scale, theory_curve.Eval(mass))

# remove needless points
while median_3sigma.GetN() > oi+1:
    obs_lim.RemovePoint(obs_lim.GetN()-1)
while median_3sigma.GetN() > ei+1:
    median_3sigma.RemovePoint(median_3sigma.GetN()-1)
    median_2sigma.RemovePoint(median_2sigma.GetN()-1)
    median_1sigma.RemovePoint(median_1sigma.GetN()-1)

# clone the median graph for plotting the expected limit
median = median_1sigma.Clone('median')
median.SetLineWidth(2)
median.SetLineColor(ROOT.kBlue)
median.SetLineStyle(ROOT.kDashed)

##############################################################################
#now draw the limits
c = TCanvas('limit')
c.SetLogy(True)

optionAxis = 'A'
#draw expected limits
median_3sigma.SetFillColor(ROOT.kRed-6)
median_3sigma.SetLineColor(ROOT.kRed-6)
setAxisLabels(median_3sigma, options.minMass, options.maxMass, options.inFemtoBarn, font)
if options.draw3sigma:
    median_3sigma.Draw(optionAxis+'3')
    optionAxis = ''

median_2sigma.SetFillColor(ROOT.kYellow)
median_2sigma.SetLineColor(ROOT.kYellow)
setAxisLabels(median_2sigma, options.minMass, options.maxMass, options.inFemtoBarn, font)
if options.draw2sigma:
    median_2sigma.Draw(optionAxis+'3')
    optionAxis = ''

median_1sigma.SetFillColor(ROOT.kGreen)
median_1sigma.SetLineColor(ROOT.kGreen)
setAxisLabels(median_1sigma, options.minMass, options.maxMass, options.inFemtoBarn, font)
if options.draw1sigma:
    median_1sigma.Draw(optionAxis+'3')
    optionAxis = ''

setAxisLabels(median, options.minMass, options.maxMass, options.inFemtoBarn, font)
if options.drawExp:
    median.Draw(optionAxis+'LX')
    optionAxis = ''

#draw observed limits
obs_lim.SetLineWidth(2)
obs_lim.SetMarkerStyle(20)
setAxisLabels(obs_lim, options.minMass, options.maxMass, options.inFemtoBarn, font)
drawStyle = 'l'
if options.drawObs:
    if options.showMarker:
        drawStyle += 'p'
    obs_lim.Draw(optionAxis+drawStyle)
    optionAxis = ''

#draw theory curve
theory_curve.SetLineWidth(2)
theory_curve.SetLineColor(ROOT.kCyan+1)
setAxisLabels(theory_curve, options.minMass, options.maxMass, options.inFemtoBarn, font)
if options.drawTheory:
    if optionAxis == '':
        theory_curve.Draw('Lsame')
    else:
        theory_curve.Draw('L')
theory_curve.SetRange(options.minMass, options.maxMass)

#legend
legend = TLegend(0.54, 0.50, 0.89, 0.88)
legend.SetTextFont(font)
legend.SetTextSize(0.04)
legend.SetBorderSize(0)
legend.SetFillColor(19)
legend.SetFillStyle(0)
if options.drawObs:
    legend.AddEntry(obs_lim, '95% CL limit', drawStyle)
if options.drawExp:
    legend.AddEntry(median, 'median expected limit', 'l')
if options.draw1sigma:
    legend.AddEntry(median_1sigma, '68% expected', 'f')
if options.draw2sigma:
    legend.AddEntry(median_2sigma, '95% expected', 'f')
if options.draw3sigma:
    legend.AddEntry(median_3sigma, '99% expected', 'f')
if options.drawTheory:
    if options.kappa == 1.:
        sigText = "Z' signal (LO)"
    else:
        sigText = "#splitline{{Z' signal (LO)}}{{#kappa = {:.1g} #times M_{{Z'}} / 100 TeV}}".format(options.kappa)
    legend.AddEntry(theory_curve, sigText, 'l')
legend.Draw('same')

tex = TLatex()
tex.SetNDC()
tex.SetTextFont(font)
tex.SetTextSize(0.04)
tex.DrawLatex(0.1, 0.91, 'CMS Preliminary, 8 TeV, 19.7 fb^{-1}')

##############################################################################
# save canvases to root file
if options.savePlots:
    output = TFile(options.workDir+'results/limitPlot.root', 'recreate')
    c.Write(c.GetName())
    obs_lim.Write()
    median_1sigma.Write()
    median_2sigma.Write()
    median_3sigma.Write()
    theory_curve.Write()
    output.Close()
    print 'Limit plot written to: '+options.workDir+'results/limitPlot.root'

# wait
raw_input("Press ENTER to quit.")

