#!/bin/bash
cmsRun hltgsftrkanalyser_cfg.py testNcSuffix=NC2 outFile=recoVsNc2_pt10_dr0p5_passHltDEle33.root
cmsRun hltgsftrkanalyser_cfg.py testNcSuffix=NC3 outFile=recoVsNc3_pt10_dr0p5_passHltDEle33.root
cmsRun hltgsftrkanalyser_cfg.py testNcSuffix=NC4 outFile=recoVsNc4_pt10_dr0p5_passHltDEle33.root
cmsRun hltgsftrkanalyser_cfg.py testNcSuffix=NC5 outFile=recoVsNc5_pt10_dr0p5_passHltDEle33.root
cmsRun hltgsftrkanalyser_cfg.py outFile=recoVsNc6_pt10_dr0p5_passHltDEle33.root

