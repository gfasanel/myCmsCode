#!/bin/bash
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root recoVsNc2_pt10_dr0p5_passHltDEle33.root

sed -i -e 's/NC2/NC3/g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root recoVsNc3_pt10_dr0p5_passHltDEle33.root

sed -i -e 's/NC3/NC4/g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root recoVsNc4_pt10_dr0p5_passHltDEle33.root

sed -i -e 's/NC4/NC5/g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root recoVsNc5_pt10_dr0p5_passHltDEle33.root

sed -i -e 's/NC5//g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root recoVsNc6_pt10_dr0p5_passHltDEle33.root

sed -i -e "s/testNcSuffix = ''/testNcSuffix = 'NC2'/g" ./hltgsftrkanalyser_cfg.py

