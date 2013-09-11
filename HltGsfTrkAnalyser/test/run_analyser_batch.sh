#!/bin/bash
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root nc6VsNc2_pt10_sharedHits3_passHltDEle33.root

sed -i -e 's/NC2/NC3/g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root nc6VsNc3_pt10_sharedHits3_passHltDEle33.root

sed -i -e 's/NC3/NC4/g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root nc6VsNc4_pt10_sharedHits3_passHltDEle33.root

sed -i -e 's/NC4/NC5/g' ./hltgsftrkanalyser_cfg.py
cmsRun hltgsftrkanalyser_cfg.py
mv histo.root nc6VsNc5_pt10_sharedHits3_passHltDEle33.root

sed -i -e 's/NC5/NC2/g' ./hltgsftrkanalyser_cfg.py

