import FWCore.ParameterSet.Config as cms

# initalize the command line parser
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('globTag', 'START53_V19D::All', # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global tag")
options.register('dataFile', 'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_8tev50ns.root', # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Data file")
options.register('outFile', 'histoele.root', # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Output file name")
options.parseArguments()

# the cms process
process = cms.Process("GSFANA")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'FT_53_V21_AN6::All'
process.GlobalTag.globaltag = options.globTag

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.dataFile, 
        #'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputQCD_8tev50ns.root'
        #'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputDY_8tev50ns.root'
        #'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputZskim.root'
        #'file:/afs/cern.ch/work/t/treis/hlt/gsftracking/CMSSW_6_2_0_patch1/src/HLTrigger/Configuration/test/outputAForPP.root'
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

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFile)
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

