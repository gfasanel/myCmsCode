#!/usr/bin/env python

import os
import re
import glob
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options]", description="Merge root files with expected limits from different jobs.\n")
parser.add_option("-w"  ,"--workDir", dest="workDir", default=os.getcwd()+'/', help="Path to working directory.")
parser.add_option("-d"  ,"--destDir", dest="destDir", default='',help="Destination for merged files.")
parser.add_option("-m"  ,"--mass", dest="massStr", default="*", help="Mass point to be merged. (* for all)")
(options, args) = parser.parse_args()

if not options.workDir[-1:] == '/':
    options.workDir = options.workDir+'/'

# in case no destination is specified save in the working directory
if options.destDir == '':
    options.destDir = options.workDir+'results/'

massDirs = glob.glob(options.workDir + 'results/toys_m%s'%(options.massStr))

for massDir in massDirs:
    mass = int(re.search(r'\d+', massDir[massDir.rfind('/')+1:]).group())
    os.system('hadd %shiggsCombineExpected.MarkovChainMC.mH120_m%d.root %s/higgsCombineExpected_??????.MarkovChainMC.mH120_m*.root'%(options.destDir, mass, massDir))

