#!/usr/bin/env python

import os
import shutil
import math
import numpy
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TFile,TF1,TH1F,TH1D,TGraph,TF1
from optparse import OptionParser

def main():
    parser = OptionParser(usage="usage: %prog [options]", description="Generate files with shape histograms for limit setting.\n")
    parser.add_option("-d"  ,"--dataSet", dest="dataset", default='MuEG', help="Dataset used to generate input histograms. Can be MuEG or SingleMu.")
    parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
    parser.add_option("-i"  ,"--input", dest="input", default=os.getcwd()+'/data/histograms_1gevbins.root', help="Path to input root file.")
    parser.add_option("-m"  ,"--minMass", dest="minMass", default="100", type="float", help="Minimal resonance mass.")
    parser.add_option("-M"  ,"--maxMass", dest="maxMass", default="100", type="float", help="Maximal resonance mass.")
    parser.add_option("-S"  ,"--massStep", dest="massStep", default="50", type="float", help="Mass step size.")
    parser.add_option("-k"  ,"--kappa", dest="kappa", default="1.", type="float", help="Kappa factor for coupling.")
    (options, args) = parser.parse_args()

    if options.dataset != 'MuEG' and options.dataset != 'SingleMu':
        raise Exception('Unknown option '+options.dataset+" for dataset. Use 'MuEG' or 'SingleMu'.")

    if not options.workDir[-1:] == '/':
        options.workDir = options.workDir+'/'
    if not os.path.exists(options.workDir):
        os.makedirs(options.workDir)
    if not os.path.exists(options.workDir+'data'):
        os.makedirs(options.workDir+'data')
        
    shm = ShapeHMaker(options.dataset, options.workDir, options.input, options.kappa)
    massPoints = numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep)
    obs_lim = TGraph(len(massPoints))
    for i, mass in enumerate(massPoints):
        #prepare shape histograms for this mass point
        shm.run(mass)

class ShapeHMaker:
    def __init__(self, dataset, work_dir, histo_file, kappa):
        print 'Setup for {0} dataset.'.format(dataset)
        self.workDir = work_dir
        self.histoFile = histo_file
        self.kappa = kappa
        self.dataHisto_name = "mass_data_obs"
        self.outfile_name = "shapeHistos_m"
        self.outfile_dir = self.workDir+"data/"
        self.relMassRes = TF1('relMassRes', '[0] + [1]*x + [2]*x**2', 0., 1.e6)
        self.accTimesEff = TF1('accTimesEff', '[0] + [1]/(x+[2]) + [3]*x', 0., 1.e6)
        self.xsec_x_br_sig = TF1('xsec_x_br_sig', '[0] * exp(-x/[1] + [2]/x - x**2/[3])', 1.e-6, 1.e6)
        self.xsec_x_br_sig.SetParameters(3.98360e-4, 397.414, 323.610, 1.05017e7)
        # define resolution error function
        self.resErrFunc = TF1('resErrFunc', 'abs(([0]+[1]*x+[2]*x**2)/([3]+[4]*x+[5]*x**2))') # not used
        self.resErrFuncLow = TF1('resErrFuncLow', 'pol0')
        self.resErrFuncHigh = TF1('resErrFuncHigh', 'pol2')
        self.resErrIntersection = 0.
        if dataset == 'SingleMu':  # singleMu dataset
            self.lumi = 19706.
            self.accTimesEff.SetParameters(0.739589, -141.339, 165.642, -2.6972e-5)
            self.relMassRes.SetParameters(0.01323, 1.434e-5, 3.288e-10)
            self.resErrFunc.SetParameters(-0.00069275, 1.55093e-6, 1.78443e-9, 0.013229, 1.43445e-5, 3.28827e-10)
            self.resErrFuncLow.SetParameter(0, 0.01849)
            self.resErrFuncHigh.SetParameters(-0.1173, 2.097e-4, -1.786e-8)
            self.resErrIntersection = 687.8
        else:  # MuEG dataset
            self.lumi = 19703.
            self.accTimesEff.SetParameters(0.764791, -129.597, 151.976, -2.93389e-5)
            self.relMassRes.SetParameters(0.0138464, 1.39666e-5, 3.99985e-10)
            self.resErrFunc.SetParameters(-0.000570474, 1.48494e-6, 1.80443e-9, 0.0138464, 1.39666e-5, 3.99985e-10)
            self.resErrFuncLow.SetParameter(0, 0.02567)
            self.resErrFuncHigh.SetParameters(-0.1103, 2.0363e-4, -1.66997e-8)
            self.resErrIntersection = 709.
        self.shapes = ['sigRes', 
                       'eleScale', 
                       'muScale', 
                       #'muonRes', 
                       'muonResSmear', 
                       'topPtReweight']
        self.massWinFactorInSigmas = 3.

    def run(self, mass):
        self.outfile = TFile(self.outfile_dir + self.outfile_name + "{:.0f}.root".format(mass), 'recreate')
        self.infile = TFile(self.histoFile)
        self.baseNames = []
        self.allNames = []
        self.blackList = []
        self.factorMin = 1.
        self.factorMax = 1.
        # find a proper range for this mass point
        # this sets self.factorMin properly
        self.findRange(mass)
        # make the signal shapes for sigRes
        self.makeSigShapeHisto(mass, 500000)
        self.makeSigShapeHisto(mass, 500000, -1.)
        self.makeSigShapeHisto(mass, 500000, 1.)
        # make the bkg shapes from the input file
        keyList = self.infile.GetListOfKeys()
        for name in self.allNames:
            self.makeShapeHisto(name, mass)
        print "Created shapes file: {0}".format(self.outfile.GetName())
        self.infile.Close()
        self.outfile.Close()

    def getSigExpEvts(self, mass):
        nExp = self.lumi * self.kappa**2 * self.xsec_x_br_sig.Eval(mass) * self.accTimesEff.Eval(mass)
        return nExp

    def getResErr(self, mass):
        #return self.resErrFunc.Eval(mass)
        if mass < self.resErrIntersection:
            return self.resErrFuncLow.Eval(mass)
        else:
            return self.resErrFuncHigh.Eval(mass)

    def makeSigShapeHisto(self, mass, nEvts, resScale=0):
        nExp = self.getSigExpEvts(mass)
        minMass = mass - self.massWinFactorInSigmas*self.factorMin*mass*self.relMassRes.Eval(mass)
        # no upper limit above a mass of 800 GeV
        if mass <= 800.:
            maxMass = mass + self.massWinFactorInSigmas*self.factorMax*mass*self.relMassRes.Eval(mass)
        else:
            self.infile.cd()
            dataHist = self.infile.Get('mass_data_obs')
            maxMass = dataHist.GetXaxis().GetXmax()
        self.outfile.cd()
        shapeHisto = TH1F('mass_sig', 'mass_sig', int(math.ceil(maxMass)-math.floor(minMass)), math.floor(minMass), math.ceil(maxMass))
        gaussShape = TF1('gaussShape', 'gaus', math.floor(minMass), math.ceil(maxMass))
        gaussShape.SetParameters(1., mass, (mass-minMass)/self.massWinFactorInSigmas)
        # upscale or downscale the width of the Gaussian if requested
        if resScale < 0.:
            shapeHisto.SetName(shapeHisto.GetName()+'_sigResDown')
            gaussShape.SetParameter(2, gaussShape.GetParameter(2) + resScale*gaussShape.GetParameter(2)*self.getResErr(mass))
        elif resScale > 0.:
            shapeHisto.SetName(shapeHisto.GetName()+'_sigResUp')
            gaussShape.SetParameter(2, gaussShape.GetParameter(2) + resScale*gaussShape.GetParameter(2)*self.getResErr(mass))
        # fill the histogram
        shapeHisto.FillRandom('gaussShape', nEvts)
        # normalise
        shapeHisto.Scale(nExp/shapeHisto.Integral())
        shapeHisto.Write()

    def makeShapeHisto(self, name, mass):
        # check if the histogram is blacklisted. If so use the base histogram
        for badName in self.blackList:
            blackListed = False
            if name == badName:
                blacklisted = True
                for baseName in self.baseNames:
                    if name[:len(baseName)] == baseName:
                        print 'Use {:s} instead of the blacklisted {:s}.'.format(baseName, name)
                        name = baseName
                        break
            if blackListed:
                break
 
        nExp = 0.
        self.infile.cd()
        bkgHist = self.infile.Get(name)
        if bkgHist:
            minMass = mass - self.massWinFactorInSigmas*self.factorMin*mass*self.relMassRes.Eval(mass)
            # no upper limit above a mass of 800 GeV
            if mass <= 800.: 
                maxMass = mass + self.massWinFactorInSigmas*self.factorMax*mass*self.relMassRes.Eval(mass)
            else:
                maxMass = bkgHist.GetXaxis().GetXmax()
            self.outfile.cd()
            shapeHisto = TH1F(name, name, int(math.ceil(maxMass)-math.floor(minMass)), math.floor(minMass), math.ceil(maxMass))
            for m in range(int(math.floor(minMass)), int(math.ceil(maxMass))):
                shapeHisto.SetBinContent(shapeHisto.FindBin(m), bkgHist.GetBinContent(bkgHist.FindBin(m)))
                shapeHisto.SetBinError(shapeHisto.FindBin(m), bkgHist.GetBinError(bkgHist.FindBin(m)))
            shapeHisto.Rebuild()
            shapeHisto.Write()
        else:
            raise Exception("I did not find histogram '{0}' in file {1}.".format(self.dataHisto_name, self.histoFile))

    def findRange(self, mass):
        self.infile.cd()
        keyList = self.infile.GetListOfKeys()
        # list of base histograms
        shapeNamess = []
        baseHistos = []
        shapeHistoss = []
        for keyB in keyList:
            if keyB.GetName()[:5] == 'mass_' and keyB.GetName()[-4:] != 'Down' and keyB.GetName()[-2:] != 'Up':
                self.baseNames.append(keyB.GetName())
                baseHistos.append(self.infile.Get(keyB.GetName()))
        # complete the list of histograms corresponding to baseHisto
        for i, baseHisto in enumerate(baseHistos):
            baseName = self.baseNames[i]
            shapeNames = []
            shapeHistos = []
            for keyS in keyList:
                keySName = keyS.GetName()
                if keySName != baseName and keySName[:len(baseName)] == baseName:
                    # check if the shape is in the list of good shapes
                    for goodShName in self.shapes:
                        index = keySName.rfind(goodShName)
                        if index > 0 and (keySName[index+len(goodShName):] == 'Up' or keySName[index+len(goodShName):] == 'Down'):
                            shapeNames.append(keySName)
                            shapeHistos.append(self.infile.Get(keySName))
                            break
            shapeNamess.append(shapeNames)
            shapeHistoss.append(shapeHistos)

        # make one global list of all histograms
        for i, bn in enumerate(self.baseNames):
            self.allNames.append(bn)
            for sn in shapeNamess[i]:
                self.allNames.append(sn)

        self.factorMin = 1.
        self.factorMax = 1.
        rangeChanged = True
        conflictDet = False
        counter = 0
        firstBadHistoName = ''
        while rangeChanged or conflictDet:
            rangeChanged = False
            conflictDet = False
            # catch infinite loops
            counter += 1
            if counter > 100:
                print 'It seems I am stuck in an infinite loop. Break.'
                break

            # define new range minimum
            minMass = mass - self.massWinFactorInSigmas*self.factorMin*mass*self.relMassRes.Eval(mass)
            # no upper limit above a mass of 800 GeV, so define range maximum from base histogram an case
            if mass <= 800.: 
                maxMass = mass + self.massWinFactorInSigmas*self.factorMax*mass*self.relMassRes.Eval(mass)
            else:
                maxMass = baseHistos[0].GetXaxis().GetXmax()
            minMassBin = baseHistos[0].FindBin(minMass)
            maxMassBin = baseHistos[0].FindBin(maxMass)
            if minMassBin > maxMassBin:
                print 'Error: Lower edge of window above upper edge.'
            minMassFromBin = baseHistos[0].GetBinLowEdge(minMassBin)
            maxMassFromBin = baseHistos[0].GetBinLowEdge(maxMassBin) + baseHistos[0].GetBinWidth(maxMassBin)
            print 'Trying range [{:.1f}, {:.1f}] GeV for mass point {:.1f} GeV.'.format(minMassFromBin, maxMassFromBin, mass)

            # loop over all base histograms
            for i, baseHisto in enumerate(baseHistos):
                baseName = self.baseNames[i]
                #print 'baseHisto: {:<12s}'.format(baseName)
                baseH_int = baseHisto.Integral(minMassBin, maxMassBin)
                # save the previous factors in case this histogram has to be blacklisted
                prevFactorMin = self.factorMin
                prevFactorMax = self.factorMax
                prevRangeChanged = rangeChanged
                # loop over shape histograms corresponding to base histogram
                for shapeHisto in shapeHistoss[i]:
                    shapeName = shapeHisto.GetName()
                    #check if the histogram is on the blacklist
                    blackListed = False
                    for badName in self.blackList:
                        if shapeName == badName:
                            blackListed = True
                    if blackListed:
                        continue

                    histoRangeChanged = True
                    while histoRangeChanged:
                        histoRangeChanged = False
                        # define new range minimum
                        minMass = mass - self.massWinFactorInSigmas*self.factorMin*mass*self.relMassRes.Eval(mass)
                        minMassBin = baseHistos[0].FindBin(minMass)
                        # define new range maximum
                        if mass <= 800.: 
                            maxMass = mass + self.massWinFactorInSigmas*self.factorMax*mass*self.relMassRes.Eval(mass)
                        else:
                             maxMass = baseHistos[0].GetXaxis().GetXmax()
                        maxMassBin = baseHistos[0].FindBin(maxMass)
                        shapeH_int = shapeHisto.Integral(minMassBin, maxMassBin)
                        #print 'baseHisto: {:<12s}, Integral={:>8g}; shapeHisto: {:<25s}, Integral={:>8g}'.format(baseHisto.GetName(), baseH_int, shapeName, shapeH_int)
                        nbins = maxMassBin-minMassBin+1
                        # windows have to be enlarged
                        if baseH_int > 0. and shapeH_int == 0.:
                            # remember the name of the first histogram that needed to change the range
                            if firstBadHistoName == '':
                                firstBadHistoName = shapeName
                            # if the range was enlarged by a previous histogram we have a loop condition if this histogram is trying to narrow it
                            # raise flag so that the initially enlarging histogram is blacklisted on the next iteration
                            if prevFactorMin < 1.:
                                conflictDet = True
                            # enlarge the range
                            if mass <= 800.:
                                # decide if widening from the upper or lower edge by looking which half has more events (advantage for lower side if odd number of bins)
                                if baseHisto.Integral(minMassBin, minMassBin+int(math.ceil(nbins/2))-1) > baseHisto.Integral(minMassBin+int(math.ceil(nbins/2)), maxMassBin):
                                    self.factorMin += 0.01
                                else:
                                    self.factorMax += 0.01
                            else:
                                self.factorMin += 0.01
                            histoRangeChanged = True
                            rangeChanged = True
                        # windows have to be narrowed
                        elif baseH_int == 0. and shapeH_int > 0.:
                            # remember the name of the first histogram that needed to change the range
                            if firstBadHistoName == '':
                                firstBadHistoName = shapeName
                            # if the range was enlarged by a previous histogram we have a loop condition if this histogram is trying to narrow it
                            # raise flag so that the initially narrowing histogram is blacklisted on the next iteration
                            if prevFactorMin > 1.:
                                conflictDet = True
                            # narrow the range
                            if mass <= 800.:
                                # decide if narrowing from the upper or lower edge by looking which half has more events (advantage for lower side if odd number of bins)
                                if shapeHisto.Integral(minMassBin, minMassBin+int(math.ceil(nbins/2))-1) > shapeHisto.Integral(minMassBin+int(math.ceil(nbins/2)), maxMassBin):
                                    self.factorMin -= 0.01
                                else:
                                    self.factorMax -= 0.01
                            else:
                                self.factorMin -= 0.01
                            histoRangeChanged = True
                            rangeChanged = True

                    if conflictDet == True:
                        print 'Conflict detected between {:s} and {:s}. Blacklist first histogram to change the range: {:s}.'.format(firstBadHistoName, shapeName, firstBadHistoName)
                        self.blackList.append(firstBadHistoName)
                        firstBadHistoName = ''
                        self.factorMin = 1
                        self.factorMax = 1
                        break

                    # in case the factor is too large blacklist the histogram and use the base histogram instead
                    if abs(self.factorMin-1.) > 0.2 or abs(self.factorMax-1.) > 0.2:
                        print 'Range modification factors exceed 20% of the nominal value. Blacklist {:s}.'.format(shapeName)
                        self.blackList.append(shapeName)
                        self.factorMin = prevFactorMin
                        self.factorMax = prevFactorMax
                        rangeChanged = prevRangeChanged
                        continue

                    # if one histogram triggered a range change stop the check for this range and start with the new range
                    if rangeChanged:
                        if self.factorMin > 1.:
                            print 'Lower edge shifted by {:.0f}% to enlarge the window of {:s} histogram from nominal width.'.format(100*(self.factorMin-1.), shapeName)
                        elif self.factorMin < 1.:
                            print 'Lower edge shifted by {:.0f}% to narrow the window of {:s} histogram from nominal width.'.format(100*(self.factorMin-1.), shapeName)
                        if self.factorMax > 1.:
                            print 'Upper edge shifted by {:.0f}% to enlarge the window of {:s} histogram from nominal width.'.format(100*(self.factorMax-1.), shapeName)
                        elif self.factorMax < 1.:
                            print 'Upper edge shifted by {:.0f}% to narrow the window of {:s} histogram from nominal width.'.format(100*(self.factorMax-1.), shapeName)
                        break 

                if rangeChanged or conflictDet:
                    break


if __name__=='__main__': main()
