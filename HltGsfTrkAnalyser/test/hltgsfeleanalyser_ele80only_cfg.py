from hltgsfeleanalyser_cfg import *

process.HltGsfEleAna = cms.EDAnalyzer('HltGsfEleAnalyser',
    # track collections
    electronCollTag = cms.InputTag('gsfElectrons', '', 'RECO'),
    trgResultsTag = cms.InputTag('TriggerResults', '', 'GSFTEST'),
    trgEventTag = cms.InputTag('hltTriggerSummaryAOD', '', 'GSFTEST'),

    # triggers 
    triggers = cms.untracked.VPSet(
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC2'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC3'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC4'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC5'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_v2'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_TrkIdT'),
            refTriggerName = cms.untracked.string('HLT_Ele80_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(80.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
    ),
)


