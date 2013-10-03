from hltgsftrkanalyser_cfg import *

process.source.fileNames = cms.untracked.vstring(
   'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim_lowEt.root'
)

# track collections
process.HltGsfTrkAna.referenceTrackCollTag = cms.InputTag('hltL1SeededElectronGsfTracks' + options.refNcSuffix, '', 'GSFTEST'),
process.HltGsfTrkAna.testTrackCollTag = cms.InputTag('hltCtfL1SeededWithMaterialTracks', '', 'GSFTEST')

# Deta and Dphi
process.HltGsfTrkAna.referenceTrackVarsTag = cms.InputTag('hltElectronL1SeededDetaDphi', '', 'GSFTEST')
process.HltGsfTrkAna.testTrackVarsTag = cms.InputTag('hltElectronL1SeededDetaDphi', '', 'GSFTEST')
process.HltGsfTrkAna.referenceEcalCandTagForDEta = cms.InputTag('hltEle25CaloIdVTOneOEMinusOneOPFilter', '', 'GSFTEST')
process.HltGsfTrkAna.testEcalCandTagForDEta = cms.InputTag('hltEle25CaloIdVTOneOEMinusOneOPFilter', '', 'GSFTEST')
process.HltGsfTrkAna.referenceEcalCandTagForDPhi = cms.InputTag('hltEle25CaloIdVTTrkIdTDetaFilter', '', 'GSFTEST')
process.HltGsfTrkAna.testEcalCandTagForDPhi = cms.InputTag('hltEle25CaloIdVTTrkIdTDetaFilter', '', 'GSFTEST')
process.HltGsfTrkAna.referenceDEtaCut = cms.untracked.double(0.008)
process.HltGsfTrkAna.referenceDPhiCut = cms.untracked.double(0.07)
process.HltGsfTrkAna.testDEtaCut = cms.untracked.double(0.008)
process.HltGsfTrkAna.testDPhiCut = cms.untracked.double(0.07)

# trigger bits
process.HltGsfTrkAna.referenceTrigger.triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT')

