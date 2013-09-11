import FWCore.ParameterSet.Config as cms

process = cms.Process("GSFANA")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_P_V42::All'
process.GlobalTag.globaltag = 'START53_V19D::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputAForPP.root'
    )
)

refNcSuffix = ''
testNcSuffix = 'NC2'
process.HltGsfTrkAna = cms.EDAnalyzer('HltGsfTrkAnalyser',
    # track collections
    #referenceTrackCollTag = cms.InputTag('electronGsfTracks', '', 'RECO'),
    referenceTrackCollTag = cms.InputTag('hltActivityElectronGsfTracks' + refNcSuffix, '', 'GSFTEST'),
    #referenceTrackCollTag = cms.InputTag('hltL1SeededElectronGsfTracks' + refNcSuffix, '', 'GSFTEST'),
    testTrackCollTag = cms.InputTag('hltActivityElectronGsfTracks' + testNcSuffix, '', 'GSFTEST'),
    #testTrackCollTag = cms.InputTag('hltL1SeededElectronGsfTracks' + testNcSuffix, '', 'GSFTEST'),

    # Deta and Dphi
    referenceTrackVarsTag = cms.InputTag('hltActivityGsfTrackVars'+ refNcSuffix, '', 'GSFTEST'),
    testTrackVarsTag = cms.InputTag('hltActivityGsfTrackVars' + testNcSuffix, '', 'GSFTEST'),
    referenceEcalCandTagForDEta = cms.InputTag('hltDiEle33CaloIdLPixelMatchDoubleFilter', '', 'GSFTEST'),
    testEcalCandTagForDEta = cms.InputTag('hltDiEle33CaloIdLPixelMatchDoubleFilter', '', 'GSFTEST'),
    referenceEcalCandTagForDPhi = cms.InputTag('hltDiEle33CaloIdLGsfTrkIdVLDEtaDoubleFilter' + refNcSuffix, '', 'GSFTEST'),
    testEcalCandTagForDPhi = cms.InputTag('hltDiEle33CaloIdLGsfTrkIdVLDEtaDoubleFilter' + testNcSuffix, '', 'GSFTEST'),
    referenceDEtaCut = cms.untracked.double(0.02),
    referenceDPhiCut = cms.untracked.double(0.15),
    testDEtaCut = cms.untracked.double(0.02),
    testDPhiCut = cms.untracked.double(0.15),

    # track matching
    #matchMethod = cms.untracked.uint32(0),      # for matching by dR
    #discrCutForMatch = cms.untracked.double(0.5),
    #largerThan = cms.untracked.bool(False),
    #matchMethod = cms.untracked.uint32(1),      # for matching by pT
    #discrCutForMatch = cms.untracked.double(1000.),
    #largerThan = cms.untracked.bool(False),
    #matchMethod = cms.untracked.uint32(2),      # for matching by dz
    #discrCutForMatch = cms.untracked.double(500.),
    #largerThan = cms.untracked.bool(False),
    matchMethod = cms.untracked.uint32(3),      # for matching by shared hits
    discrCutForMatch = cms.untracked.double(3.),
    largerThan = cms.untracked.bool(True),

    # pt cuts for _ptCut_ histograms
    referenceMinPt = cms.untracked.double(10.),
    testMinPt = cms.untracked.double(0.),

    # trigger bits
    referenceTrigger = cms.untracked.PSet(
        triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_v14'),
        #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC2'),
        #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC3'),
        #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC4'),
        #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC5'),
        #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7'),
        #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_v2'),
        #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC5'),
        #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC4'),
        #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC3'),
        #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC2'),
        #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_TrkIdT'),
        invertBit = cms.untracked.bool(False)
    ),
    #testTrigger = cms.untracked.PSet(
    #    triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC2'),
    #    #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC3'),
    #    #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC4'),
    #    #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_nC5'),
    #    #triggerName = cms.untracked.string('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7'),
    #    #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_v2'),
    #    #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC5'),
    #    #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC4'),
    #    #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC3'),
    #    #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_GsfTrkIdT_nC2'),
    #    #triggerName = cms.untracked.string('HLT_Ele80_CaloIdVT_TrkIdT'),
    #    invertBit = cms.untracked.bool(False)
    #),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)

process.p = cms.Path(
    process.HltGsfTrkAna
)

#process.writeOutput = cms.OutputModule( "PoolOutputModule",
#    fileName = cms.untracked.string( "outputTest.root" ),
#    outputCommands = cms.untracked.vstring( 'drop *',
#      'keep *_*_*_GSFANA',
#    )
#)
#
#process.pEnd = cms.EndPath(process.writeOutput)

