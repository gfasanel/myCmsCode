#!/usr/bin/env python

import os
import shutil
import random
import glob

from optparse import OptionParser
from string import ascii_letters, digits

parser = OptionParser(usage="usage: %prog [options]", description="Calculate the observed limits.\n")
parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
parser.add_option("-m"  ,"--mass", dest="mass", default="100", type="float", help="Resonance mass.")
parser.add_option("-T"  ,"--tries", dest="nTries", default="10" , type="int",help="Number of tries for MCMC.")
parser.add_option("-s"  ,"--seed", dest="seed", default="-1", type="int", help="Seed for pseudo random number generator. (-1 for random seed)")
parser.add_option("-v"  ,"--verbosity", dest="verbosity", default="0", type="int", help="Verbosity level. -1 (very quite) to 2 (debug)")
parser.add_option("-b"  ,"--batch", dest="batch", default=False, action="store_true", help="Run on cluster.")
(options, args) = parser.parse_args()

if not options.workDir[:1] == '/':
    if not options.workDir[:2] == './':
        options.workDir = os.getcwd()+'/'+options.workDir
if not options.workDir[-1:] == '/':
    options.workDir = options.workDir+'/'
if not os.path.exists(options.workDir+'results'):
    os.makedirs(options.workDir+'results')

#make the limit tool available by setting up cmssw
shCommands = ['true']
if options.batch:
    shCommands.append('source $VO_CMS_SW_DIR/cmsset_default.sh')
    shCommands.append('cd ' + options.workDir)
    shCommands.append('eval `scram runtime -sh`')
uid = ''.join([random.choice(ascii_letters + digits) for n in xrange(6)]) 
shCommands.append("combine -M MarkovChainMC -H ProfileLikelihood --tries %d  -n Observed_%s -s %d -v %d %scards/emuLimitCard_m%d.txt" % (options.nTries, uid, options.seed, options.verbosity, options.workDir, options.mass))
os.system('; '.join(shCommands))
#copy output root file results
rootFile_paths = glob.glob(options.workDir + 'higgsCombineObserved_%s.MarkovChainMC.mH120.*.root' % (uid))
if len(rootFile_paths) > 0:
    rootFile_name = rootFile_paths[0][len(options.workDir):]
    resultsDir = options.workDir + 'results/'
    #shutil.move(rootFile_paths[0], resultsDir + rootFile_name.replace('MC.mH120', 'MC.mH120_m{0}'.format(int(options.mass))))
    shutil.move(rootFile_paths[0], resultsDir + 'higgsCombineObserved.MarkovChainMC.mH120_m{0}.root'.format(int(options.mass)))

