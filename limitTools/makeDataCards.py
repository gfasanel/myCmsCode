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
    if not os.path.exists(options.workDir+'cards'):
        os.makedirs(options.workDir+'cards')

    dcm = DcMaker(options.workDir, options.input, options.kappa)
    massPoints = numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep)
    obs_lim = TGraph(len(massPoints))
    for i, mass in enumerate(massPoints):
        #prepare card for this mass point
        dcm.run(mass)

class DcMaker:
    divider_str = "------------------------------------------------------------------------------- \n"
    
    def __init__(self, work_dir, histo_file, kappa):
        self.workDir = work_dir
        self.histoFile = histo_file
        self.kappa = kappa
        self.lumi = 19703.
        self.accTimesEff = TF1('accTimesEff', '[0] + [1]/(x+[2]) + [3]*x', 0., 1.e6)
        self.accTimesEff.SetParameters(0.74, -131.8, 151.9, -2.81e-5)
        self.xsec_x_br_sig = TF1('xsec_x_br_sig', '[0] * exp(-x/[1] + [2]/x - x**2/[3])', 1.e-6, 1.e6)
        self.xsec_x_br_sig.SetParameters(3.98360e-4, 397.414, 323.610, 1.05017e7)
        self.outfile_name = "emuLimitCard_m"
        self.outfile_dir = self.workDir+"cards/"

    def run(self, mass):
        outfile = open(self.outfile_dir + self.outfile_name + "{:.0f}.txt".format(mass), 'w')
        print "Created card file: {0}".format(outfile.name)
        self.writeSetup(outfile, mass)
        self.writeBins(outfile, mass)
        self.writeExpected(outfile, mass)
        self.writeUncert(outfile, mass)
        outfile.close()

    def getNEvts(self, name, mass):
        nEvt = 0
        file = TFile(self.workDir+'data/shapeHistos_m{:.0f}.root'.format(mass))
        hist = file.Get("mass_"+name)
        if hist:
            nEvt = hist.Integral()
            return nEvt
        else:
            raise Exception("I did not find histogram 'mass_{0}' in file {1}.".format(name, self.histoFile))

    def getSigExpEvts(self, mass):
        nExp = self.lumi * self.kappa**2 * self.xsec_x_br_sig.Eval(mass) * self.accTimesEff.Eval(mass)
        return nExp

    def writeSetup(self, outfile, mass):
        imax_str = "imax    {0}     number of channels \n".format(1)
        jmax_str = "jmax    {0}     number of backgrounds \n".format(9)
        kmax_str = "kmax    *     number of nuisance parameters \n"
        shape_str = "shapes * * ./data/shapeHistos_m{:.0f}.root mass_$PROCESS mass_$PROCESS_$SYSTEMATIC \n".format(mass)
        outfile.write(imax_str+jmax_str+kmax_str+shape_str+self.divider_str)

    def writeBins(self, outfile, mass):
        bin_str = "bin              zp_lfv \n"
        observation_str = "observation {:11g} \n".format(self.getNEvts('data_obs', mass))
        outfile.write(bin_str+observation_str+self.divider_str)

    def writeExpected(self, outfile, mass):
        bin_str = "bin   {:>17} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('zp_lfv','zp_lfv','zp_lfv','zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv', 'zp_lfv')
        procName_str = "process   {:>13} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('sig', 'ttbar', 'ww', 'wz', 'zz', 'tw', 'ztautau', 'zmumu', 'zee', 'qcd')
        procNum_str = "process   {:>13} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        expYield_str = "rate   {:>16.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} {:>9.3f} \n".format(self.getSigExpEvts(mass), self.getNEvts('ttbar', mass), self.getNEvts('ww', mass), self.getNEvts('wz', mass), self.getNEvts('zz', mass), self.getNEvts('tw', mass), self.getNEvts('ztautau', mass), self.getNEvts('zmumu', mass), self.getNEvts('zee', mass), self.getNEvts('qcd', mass))
        outfile.write(bin_str+procName_str+procNum_str+expYield_str+self.divider_str)

    def writeUncert(self, outfile, mass):
        unc_str = "lumi     lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026, 1.026)
        outfile.write(unc_str)
        unc_str = "bkgXsec  lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1.036, 1.035, 1.038, 1.025, 1.069, 1.054, 1.054, 1.054, '-')
        outfile.write(unc_str)
        unc_str = "fakeRate lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', '-', '-', '-', '-', '-', '-', '-', '-', 1.3)
        outfile.write(unc_str)
        pdf_unc = 1.045 + 3.4e-5*mass + 1.5e-8*mass**2
        unc_str = "pdf      lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, pdf_unc, '-')
        outfile.write(unc_str)
        unc_str = "accXeff  lnN   {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05)
        outfile.write(unc_str)
        unc_str = "sigRes   shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format(1, '-', '-', '-', '-', '-', '-', '-', '-', '-')
        outfile.write(unc_str)
        unc_str = "eleScale shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1, 1, 1, 1, 1, 1, 1, 1, '-')
        outfile.write(unc_str)
        unc_str = "muScale  shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1, 1, 1, 1, 1, 1, 1, 1, '-')
        outfile.write(unc_str)
        unc_str = "muonRes  shape {:>8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} \n".format('-', 1, 1, 1, 1, 1, 1, 1, 1, '-')
        outfile.write(unc_str)

if __name__=='__main__': main()
