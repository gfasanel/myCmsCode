from hltgsftrkanalyser_lowEt_kfTrk_cfg import *

process.HltGsfTrkAna.referenceTrackCollTag = cms.InputTag('electronGsfTracks', '', 'RECO')

process.HltGsfTrkAna.matchMethod = cms.untracked.uint32(0)      # for matching by dR
process.HltGsfTrkAna.discrCutForMatch = cms.untracked.double(0.5)
process.HltGsfTrkAna.largerThan = cms.untracked.bool(False)

