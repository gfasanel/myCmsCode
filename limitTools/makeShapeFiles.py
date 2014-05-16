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
    parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
    parser.add_option("-i"  ,"--input", dest="input", default=os.getcwd()+'/data/histograms_1gevbins.root', help="Path to input root file.")
    parser.add_option("-m"  ,"--minMass", dest="minMass", default="100", type="float", help="Minimal resonance mass.")
    parser.add_option("-M"  ,"--maxMass", dest="maxMass", default="100", type="float", help="Maximal resonance mass.")
    parser.add_option("-S"  ,"--massStep", dest="massStep", default="50", type="float", help="Mass step size.")
    parser.add_option("-k"  ,"--kappa", dest="kappa", default="1.", type="float", help="Kappa factor for coupling.")
    (options, args) = parser.parse_args()

    if not options.workDir[-1:] == '/':
        options.workDir = options.workDir+'/'
    if not os.path.exists(options.workDir):
        os.makedirs(options.workDir)
    if not os.path.exists(options.workDir+'data'):
        os.makedirs(options.workDir+'data')
        
    shm = ShapeHMaker(options.workDir, options.input, options.kappa)
    massPoints = numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep)
    obs_lim = TGraph(len(massPoints))
    for i, mass in enumerate(massPoints):
        #prepare shape histograms for this mass point
        shm.run(mass)

class ShapeHMaker:
    def __init__(self, work_dir, histo_file, kappa):
        self.workDir = work_dir
        self.histoFile = histo_file
        self.kappa = kappa
        self.lumi = 19703.
        self.relMassRes = TF1('relMassRes', '[0] + [1]*x + [2]*x**2', 0., 1.e6)
        self.relMassRes.SetParameters(0.0132055, 1.53824e-5, -9.33481e-11)
        self.accTimesEff = TF1('accTimesEff', '[0] + [1]/(x+[2]) + [3]*x', 0., 1.e6)
        self.accTimesEff.SetParameters(0.74, -131.8, 151.9, -2.81e-5)
        self.xsec_x_br_sig = TF1('xsec_x_br_sig', '[0] * exp(-x/[1] + [2]/x - x**2/[3])', 1.e-6, 1.e6)
        self.xsec_x_br_sig.SetParameters(3.98360e-4, 397.414, 323.610, 1.05017e7)
        # define resolution error function
        self.resErrFunc = TF1('resErrFunc', 'abs(([0]+[1]*x+[2]*x**2)/([3]+[4]*x+[5]*x**2))')
        self.resErrFunc.SetParameters(-0.00162388, 3.74716e-6, 1.04324e-9, 0.0132055, 1.53824e-5, -9.33481e-11)
        self.dataHisto_name = "mass_data_obs"
        self.outfile_name = "shapeHistos_m"
        self.outfile_dir = self.workDir+"data/"

    def run(self, mass):
        self.outfile = TFile(self.outfile_dir + self.outfile_name + "{:.0f}.root".format(mass), 'recreate')
        self.infile = TFile(self.histoFile)
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
        for key in keyList:
            if key.GetName()[:5] == 'mass_': 
                self.makeBkgShapeHisto(key.GetName(), mass)
        print "Created shapes file: {0}".format(self.outfile.GetName())
        self.infile.Close()
        self.outfile.Close()

    def getSigExpEvts(self, mass):
        nExp = self.lumi * self.kappa**2 * self.xsec_x_br_sig.Eval(mass) * self.accTimesEff.Eval(mass)
        return nExp

    def makeSigShapeHisto(self, mass, nEvts, resScale=0):
        nExp = self.getSigExpEvts(mass)
        minMass = mass - 3*self.factorMin*mass*self.relMassRes.Eval(mass)
        # no upper limit above a mass of 800 GeV
        if mass <= 800.:
            maxMass = mass + 3*self.factorMax*mass*self.relMassRes.Eval(mass)
        else:
            self.infile.cd()
            dataHist = self.infile.Get('mass_data_obs')
            maxMass = dataHist.GetXaxis().GetXmax()
        self.outfile.cd()
        shapeHisto = TH1F('mass_sig', 'mass_sig', int(math.ceil(maxMass)-math.floor(minMass)), math.floor(minMass), math.ceil(maxMass))
        gaussShape = TF1('gaussShape', 'gaus', math.floor(minMass), math.ceil(maxMass))
        gaussShape.SetParameters(1., mass, (mass-minMass)/3.)
        # upscale or downscale the width of the Gaussian if requested
        if resScale < 0.:
            shapeHisto.SetName(shapeHisto.GetName()+'_sigResDown')
            gaussShape.SetParameter(2, gaussShape.GetParameter(2) + resScale*gaussShape.GetParameter(2)*self.resErrFunc.Eval(mass))
        elif resScale > 0.:
            shapeHisto.SetName(shapeHisto.GetName()+'_sigResUp')
            gaussShape.SetParameter(2, gaussShape.GetParameter(2) + resScale*gaussShape.GetParameter(2)*self.resErrFunc.Eval(mass))
        # fill the histogram
        shapeHisto.FillRandom('gaussShape', nEvts)
        # normalise
        shapeHisto.Scale(nExp/shapeHisto.Integral())
        shapeHisto.Write()

    def makeBkgShapeHisto(self, bkg_name, mass):
        nExp = 0.
        self.infile.cd()
        bkgHist = self.infile.Get(bkg_name)
        if bkgHist:
            minMass = mass - 3*self.factorMin*mass*self.relMassRes.Eval(mass)
            # no upper limit above a mass of 800 GeV
            if mass <= 800.: 
                maxMass = mass + 3*self.factorMax*mass*self.relMassRes.Eval(mass)
            else:
                maxMass = bkgHist.GetXaxis().GetXmax()
            self.outfile.cd()
            shapeHisto = TH1F(bkg_name, bkg_name, int(math.ceil(maxMass)-math.floor(minMass)), math.floor(minMass), math.ceil(maxMass))
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
        baseNames = []
        shapeNamess = []
        baseHistos = []
        shapeHistoss = []
        for keyB in keyList:
            if keyB.GetName()[:5] == 'mass_' and keyB.GetName()[-4:] != 'Down' and keyB.GetName()[-2:] != 'Up':
                baseNames.append(keyB.GetName())
                baseHistos.append(self.infile.Get(keyB.GetName()))
        # complete the list of histograms corresponding to baseHisto
        for i, baseHisto in enumerate(baseHistos):
            baseName = baseNames[i]
            shapeNames = []
            shapeHistos = []
            for keyS in keyList:
                if keyS.GetName() != baseName and keyS.GetName()[:len(baseName)] == baseName:
                    shapeNames.append(keyS.GetName())
                    shapeHistos.append(self.infile.Get(keyS.GetName()))
            shapeNamess.append(shapeNames)
            shapeHistoss.append(shapeHistos)

        self.factorMin = 1.
        rangeChanged = True
        while rangeChanged:
            rangeChanged = False
            # define new range minimum
            minMass = mass - 3*self.factorMin*mass*self.relMassRes.Eval(mass)
            # no upper limit above a mass of 800 GeV, so define range maximum from base histogram an case
            if mass <= 800.: 
                maxMass = mass + 3*self.factorMax*mass*self.relMassRes.Eval(mass)
            else:
                maxMass = baseHistos[0].GetXaxis().GetXmax()
            minMassBin = baseHistos[0].FindBin(minMass)
            maxMassBin = baseHistos[0].FindBin(maxMass)
            minMassFromBin = baseHistos[0].GetBinLowEdge(minMassBin)
            maxMassFromBin = baseHistos[0].GetBinLowEdge(maxMassBin) + baseHistos[0].GetBinWidth(maxMassBin)
            print 'Trying range [{:.1f}, {:.1f}] GeV for mass point {:.1f} GeV.'.format(minMassFromBin, maxMassFromBin, mass)

            # loop over all base histograms
            for i, baseHisto in enumerate(baseHistos):
                baseName = baseNames[i]
                #print 'baseHisto: {:<12s}'.format(baseName)
                baseH_int = baseHisto.Integral(minMassBin, maxMassBin)
                 # loop over shape histograms corresponding to base histogram
                for shapeHisto in shapeHistoss[i]:
                    histoRangeChanged = True
                    while histoRangeChanged:
                        histoRangeChanged = False
                        # define new range minimum
                        minMass = mass - 3*self.factorMin*mass*self.relMassRes.Eval(mass)
                        minMassBin = baseHistos[0].FindBin(minMass)
                        shapeH_int = shapeHisto.Integral(minMassBin, maxMassBin)
                        #print 'baseHisto: {:<12s}, Integral={:>8g}; shapeHisto: {:<25s}, Integral={:>8g}'.format(baseHisto.GetName(), baseH_int, shapeHisto.GetName(), shapeH_int)
                        # windows have to be enlarged
                        if baseH_int > 0. and shapeH_int == 0.:
                            self.factorMin += 0.01
                            histoRangeChanged = True
                            rangeChanged = True
                        # windows have to be narrowed
                        elif baseH_int == 0. and shapeH_int > 0.:
                            self.factorMin -= 0.01
                            histoRangeChanged = True
                            rangeChanged = True

                # if one histogram triggered a range change stop the check for this range and start with the new range
                if rangeChanged:
                    if self.factorMin > 1.:
                        print 'Lower edge shifted by {:.0f}% to enlarge the window of {:s}* histograms from nominal width.'.format(100*(self.factorMin-1.), baseName)
                    elif self.factorMin < 1.:
                        print 'Lower edge shifted by {:.0f}% to narrow  the window of {:s}* histograms from nominal width.'.format(100*(self.factorMin-1.), baseName)
                    break 

if __name__=='__main__': main()
