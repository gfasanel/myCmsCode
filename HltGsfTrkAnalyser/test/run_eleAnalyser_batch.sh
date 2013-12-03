#!/bin/bash

##############################################################################
# DoubleEle33 + Ele80 paths
cmsRun hltgsfeleanalyser_cfg.py \
globTag=FT_53_V21_AN6::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim.root \
outFile=histoele_dataZskim.root

cmsRun hltgsfeleanalyser_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_8tev50ns.root \
outFile=histoele_dy_8tev50ns.root

cmsRun hltgsfeleanalyser_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputQCD_8tev50ns.root \
outFile=histoele_qcd_8tev50ns.root

##############################################################################
# DoubleEle33
cmsRun hltgsfeleanalyser_doubleele33only_cfg.py \
globTag=FT_53_V21_AN6::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim.root \
outFile=histoele_dataZskim_doubleele33only.root

cmsRun hltgsfeleanalyser_doubleele33only_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_8tev50ns.root \
outFile=histoele_dy_8tev50ns_doubleele33only.root

cmsRun hltgsfeleanalyser_doubleele33only_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputQCD_8tev50ns.root \
outFile=histoele_qcd_8tev50ns_doubleele33only.root

##############################################################################
# Ele80 paths
cmsRun hltgsfeleanalyser_ele80only_cfg.py \
globTag=FT_53_V21_AN6::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim.root \
outFile=histoele_dataZskim_ele80only.root

cmsRun hltgsfeleanalyser_ele80only_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_8tev50ns.root \
outFile=histoele_dy_8tev50ns_ele80only.root

cmsRun hltgsfeleanalyser_ele80only_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputQCD_8tev50ns.root \
outFile=histoele_qcd_8tev50ns_ele80only.root

##############################################################################
# Ele25 paths
cmsRun hltgsfeleanalyser_lowEt_cfg.py \
globTag=FT_53_V21_AN6::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim_lowEt.root \
outFile=histoele_dataZskim_lowEtMenu.root

cmsRun hltgsfeleanalyser_lowEt_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_lowEt_8tev50ns.root \
outFile=histoele_dy_8tev50ns_lowEtMenu.root

cmsRun hltgsfeleanalyser_lowEt_cfg.py \
globTag=START53_V19D::All \
dataFile=file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputQCD_lowEt_8tev50ns.root \
outFile=histoele_qcd_8tev50ns_lowEtMenu.root

