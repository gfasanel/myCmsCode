from hltgsfeleanalyser_cfg import *

process.HltGsfEleAna = cms.EDAnalyzer('HltGsfEleAnalyser',
    # track collections
    electronCollTag = cms.InputTag('gsfElectrons', '', 'RECO'),
    trgResultsTag = cms.InputTag('TriggerResults', '', 'GSFTEST'),
    trgEventTag = cms.InputTag('hltTriggerSummaryAOD', '', 'GSFTEST'),

    # triggers 
    triggers = cms.untracked.VPSet(
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC2'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC3'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC4'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC5'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdT_v10'),
            refTriggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(2),
            minEts = cms.untracked.vdouble(35., 35.),
            tpTriggerName = cms.untracked.string('HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v7'),
            tagFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),
            probeFilterName = cms.untracked.string('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
        ),
    ),
)


