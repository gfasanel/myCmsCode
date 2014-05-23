#!/usr/bin/env python

import os
import shutil
import math
import numpy
import ROOT
#ROOT.gROOT.SetBatch(True)
from ROOT import TFile,TH1F,TH1D,TGraph,TCanvas,TF1
from optparse import OptionParser

def main():
    parser = OptionParser(usage="usage: %prog [options]", description="Generate cards for limit setting.\n")
    parser.add_option("-d"  ,"--dataSet", dest="dataset", default='MuEG', help="Dataset used to generate input histograms. Can be MuEG or SingleMu.")
    parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
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
    if not os.path.exists(options.workDir+'cards'):
        os.makedirs(options.workDir+'cards')

    dcm = DcMaker(options.dataset, options.workDir, options.kappa)
    massPoints = numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep)
    obs_lim = TGraph(len(massPoints))
    for i, mass in enumerate(massPoints):
        #prepare card for this mass point
        dcm.run(mass)

class DcMaker:
    divider_str = "------------------------------------------------------------------------------- \n"
    
    def __init__(self, dataset, work_dir, kappa):
        self.workDir = work_dir
        self.kappa = kappa
        self.outfile_name = "emuLimitCard_m"
        self.outfile_dir = self.workDir+"cards/"
        self.accTimesEff = TF1('accTimesEff', '[0] + [1]/(x+[2]) + [3]*x', 0., 1.e6)
        self.xsec_x_br_sig = TF1('xsec_x_br_sig', '[0] * exp(-x/[1] + [2]/x - x**2/[3])', 1.e-6, 1.e6)
        self.xsec_x_br_sig.SetParameters(3.98360e-4, 397.414, 323.610, 1.05017e7)
        self.pdf_unc = TF1('pdf_unc', '[0] + [1]*x + [2]*x**2', 0., 1.e6)
        self.pdf_unc.SetParameters(1.045, 3.4e-5, 1.5e-8)
        if dataset == 'SingleMu':  # singleMu dataset
            self.lumi = 19706.
            self.accTimesEff.SetParameters(0.748776, -161.321, 206.203, -3.11382e-5)
        else:  # MuEG dataset
            self.lumi = 19703.
            self.accTimesEff.SetParameters(0.755094, -129.296, 153.844, -2.77505e-5)
 
    def run(self, mass):
        self.mass = mass
        self.outfile = open(self.outfile_dir + self.outfile_name + "{:.0f}.txt".format(self.mass), 'w')
        self.writeSetup()
        self.writeBins()
        self.writeExpected()
        self.writeUncert()
        print "Created card file: {0}".format(self.outfile.name)
        self.outfile.close()

    def getNEvts(self, name):
        nEvt = 0.
        file = TFile(self.workDir+'data/shapeHistos_m{:.0f}.root'.format(self.mass))
        hist = file.Get("mass_"+name)
        if hist:
            nEvt = hist.Integral()
            return nEvt
        else:
            raise Exception("I did not find histogram 'mass_{:s}' in file {:s}./shapeHistos_m{:.0f}.root.".format(name, self.workDir, self.mass))

    def getSigExpEvts(self):
        nExp = self.lumi * self.kappa**2 * self.xsec_x_br_sig.Eval(self.mass) * self.accTimesEff.Eval(self.mass)
        return nExp

    def getStatErr(self, name):
        statErr = 1.
        file = TFile(self.workDir+'data/shapeHistos_m{:.0f}.root'.format(self.mass))
        hist = file.Get("mass_"+name)
        if hist:
            nEvt = hist.Integral()
            binErr = 0.
            for bin in range(1, hist.GetNbinsX()+1):
                binErr += hist.GetBinError(bin)**2
            if nEvt > 0.:
                statErr += math.sqrt(binErr)/nEvt
                return statErr
            else:
                return '-'
        else:
            raise Exception("I did not find histogram 'mass_{:s}' in file {:s}./shapeHistos_m{:.0f}.root.".format(name, self.workDir, self.mass))


    def writeSetup(self):
        imax_str = "imax    {0}     number of channels \n".format(1)
        jmax_str = "jmax    {0}     number of backgrounds \n".format(9)
        kmax_str = "kmax    *     number of nuisance parameters \n"
        shape_str = "shapes * * ./data/shapeHistos_m{:.0f}.root mass_$PROCESS mass_$PROCESS_$SYSTEMATIC \n".format(self.mass)
        self.outfile.write(imax_str+jmax_str+kmax_str+shape_str+self.divider_str)

    def writeBins(self):
        bin_str = "bin              zp_lfv \n"
        observation_str = "observation {:11g} \n".format(self.getNEvts('data_obs'))
        self.outfile.write(bin_str+observation_str+self.divider_str)

    def writeExpected(self):
        bin_str = "bin   {:>17} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('zp_lfv','zp_lfv','zp_lfv','zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv')
        procName_str = "process   {:>13} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('sig', 'ttbar', 'ww', 'wz', 'zz', 'tw', 'ztautau', 'zmumu', 'zee', 'qcd')
        procNum_str = "process   {:>13} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        expYield_str = "rate   {:>16.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} \n".format(self.getSigExpEvts(), self.getNEvts('ttbar'), self.getNEvts('ww'), self.getNEvts('wz'), self.getNEvts('zz'), self.getNEvts('tw'), self.getNEvts('ztautau'), self.getNEvts('zmumu'), self.getNEvts('zee'), self.getNEvts('qcd'))
        self.outfile.write(bin_str+procName_str+procNum_str+expYield_str+self.divider_str)

    def writeUncert(self):
        self.outfile.write("lumi     lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026, '-'))
        self.outfile.write("bkgXsec  lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1.036, 1.035, 1.038, 1.025, 1.069, 1.054, 1.054, 1.054, '-'))
        self.outfile.write("fakeRate lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', '-', '-', '-', '-', '-', '-', '-', '-', 1.3))
        pdf_unc = self.pdf_unc.Eval(self.mass)
        self.outfile.write("pdf      lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, '-'))
        self.outfile.write("accXeff  lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, '-'))
        self.outfile.write("statMC   lnN   {:>8} {:>9.6} {:>9.6} {:>9.6} {:>9.6} {:>9.6} {:>9.6} {:>9.6} {:>9.6} {:>9} \n".format('-', self.getStatErr('ttbar'), self.getStatErr('ww'), self.getStatErr('wz'), self.getStatErr('zz'), self.getStatErr('tw'), self.getStatErr('ztautau'), self.getStatErr('zmumu'), self.getStatErr('zee'), '-'))
        self.outfile.write("sigRes   shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(1, '-', '-', '-', '-', '-', '-', '-', '-', '-'))
        self.outfile.write("eleScale shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1, 1, 1, 1, 1, 1, 1, 1, '-'))
        self.outfile.write("muScale  shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1, 1, 1, 1, 1, 1, 1, 1, '-'))
        self.outfile.write("muonRes  shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1, 1, 1, 1, 1, 1, 1, 1, '-'))

if __name__=='__main__': main()
