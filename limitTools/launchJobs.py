#!/usr/bin/env python

import os
import numpy
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options]", description="Launch jobs to calculate the expected limits.\n")
parser.add_option("-j"  ,"--jobs", dest="nJobs", default="1", type="int", help="Number of jobs per mass point to be launched.")
parser.add_option("-q"  ,"--queue", dest="queue", default="8nm", help="Queue to submit the jobs to.")
parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
parser.add_option("-m"  ,"--minMass", dest="minMass", default="100", type="float", help="Minimal resonance mass.")
parser.add_option("-M"  ,"--maxMass", dest="maxMass", default="100", type="float", help="Maximal resonance mass.")
parser.add_option("-S"  ,"--massStep", dest="massStep", default="100", type="float", help="Mass step size.")
parser.add_option("-t"  ,"--toys", dest="nToys", default="-1", type="int", help="Number of toy experiments. Calculate observed limts if < 1 or not specified.")
parser.add_option("-T"  ,"--tries", dest="nTries", default="10" , type="int",help="Number of tries for MCMC.")
parser.add_option("-s"  ,"--seed", dest="seed", default="-1", type="int", help="Seed for pseudo random number generator. (-1 for random seed)")
parser.add_option("-v"  ,"--verbosity", dest="verbosity", default="0", type="int", help="Verbosity level. -1 (very quite) to 2 (debug)")
(options, args) = parser.parse_args()

if not options.workDir[:1] == '/':
    if not options.workDir[:2] == './':
        options.workDir = os.getcwd()+'/'+options.workDir
if not options.workDir[-1:] == '/':
    options.workDir = options.workDir+'/'

cwd = os.getcwd()

if options.nToys > 0:
    for mass in numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep):
        for job in range(options.nJobs):
            os.system('bsub -R "pool>100" -cwd %s -q %s -J j_exp_m%d_%d %s/calcExpLimit.py -w %s -m %f -t %d -T %d -s %d -v %d -b'%(options.workDir, options.queue, mass, job, cwd, options.workDir, mass, options.nToys, options.nTries, options.seed, options.verbosity))
else:
    for mass in numpy.arange(options.minMass, options.maxMass+options.massStep, options.massStep):
        os.system('bsub -R "pool>100" -cwd %s -q %s -J j_obs_m%d %s/calcObsLimit.py -w %s -m %f -T %d -s %d -v %d -b'%(options.workDir, options.queue, mass, cwd, options.workDir, mass, options.nTries, options.seed, options.verbosity))

os.system('bjobs')

