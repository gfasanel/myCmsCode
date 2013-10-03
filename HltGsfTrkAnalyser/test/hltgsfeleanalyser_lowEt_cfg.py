import FWCore.ParameterSet.Config as cms

process = cms.Process("GSFANA")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'FT_53_V21_AN6::All'
process.GlobalTag.globaltag = 'START53_V19D::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputQCD_lowEt_8tev50ns.root'
        #'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_lowEt_8tev50ns.root'
        #'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim_lowEt.root'
    ),
    secondaryFileNames = cms.untracked.vstring(),
)

process.HltGsfEleAna = cms.EDAnalyzer('HltGsfEleAnalyser',
    # track collections
    electronCollTag = cms.InputTag('gsfElectrons', '', 'RECO'),
    trgResultsTag = cms.InputTag('TriggerResults', '', 'GSFTEST'),
    trgEventTag = cms.InputTag('hltTriggerSummaryAOD', '', 'GSFTEST'),

    # triggers 
    triggers = cms.untracked.VPSet(
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT_GsfTrkIdT_nC2'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT_GsfTrkIdT_nC3'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT_GsfTrkIdT_nC4'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT_GsfTrkIdT_nC5'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT_GsfTrkIdT_v2'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
        cms.untracked.PSet(
            triggerName = cms.untracked.string('HLT_Ele25_CaloIdVT_TrkIdT'),
            refTriggerName = cms.untracked.string('HLT_Ele25_CaloIdVT'),
            invertBit = cms.untracked.bool(False),
            minEle = cms.untracked.uint32(1),
            minEts = cms.untracked.vdouble(25.),
            tpTriggerName = cms.untracked.string('HLT_Ele27_WP80_v13'),
            tagFilterName = cms.untracked.string('hltEle27WP80TrackIsoFilter'),
        ),
    ),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histoele.root')
)

# The next three lines are for rho computation (energy density, highly correlated to PU), see here :
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRecipesFor2011#FastJet_based_pile_up_isolation
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

process.p = cms.Path(
    process.kt6PFJets * 
    process.HltGsfEleAna
)

#process.writeOutput = cms.OutputModule( "PoolOutputModule",
#    fileName = cms.untracked.string( "outputTest.root" ),
#    outputCommands = cms.untracked.vstring( 'drop *',
#      'keep *_*_*_GSFANA',
#    )
#)
#
#process.pEnd = cms.EndPath(process.writeOutput)

