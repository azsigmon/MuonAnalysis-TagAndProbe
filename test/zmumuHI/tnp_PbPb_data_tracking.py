import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:Dimuon_Skim_Mass40_450.root'),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")

process.GlobalTag.globaltag = cms.string('75X_dataRun2_v12')

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")

## ==== Filters ====
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("hiSelectedVertex"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),
)

process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')

process.fastFilter = cms.Sequence(process.goodVertexFilter * process.hfCoincFilter3)

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    mergeCaloMuons = cms.bool(False), # AOD
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("hiGeneralTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("pt > 3 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 3"),
    tracksCut    = cms.string("pt > 3"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")

process.patMuonsWithoutTrigger.pvSrc = "hiSelectedVertex"
process.patTriggerFull.l1GtReadoutRecordInputTag = cms.InputTag("gtDigis","","RECO")                 
process.patTrigger.collections.remove("hltL3MuonCandidates")
process.patTrigger.collections.append("hltHIL3MuonCandidates")
process.muonMatchHLTL3.matchedCuts = cms.string('coll("hltHIL3MuonCandidates")')

## ==== Tag and probe definitions
from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")
process.muonImpactParameter.vertices = cms.InputTag("hiSelectedVertex")
process.muonDxyPVdzmin.vertices = cms.InputTag("hiSelectedVertex")
process.muonImpactParameterTags.vertices = cms.InputTag("hiSelectedVertex")
process.muonDxyPVdzminTags.vertices = cms.InputTag("hiSelectedVertex")
from MuonAnalysis.TagAndProbe.heavyIon_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.heavyIon_modules_cff")

InAcceptance = '(abs(eta)<2.4 && pt>=15)'
TightId =  "isGlobalMuon && globalTrack.normalizedChi2 < 10 && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1 && track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0 && abs(dB) < 0.2"

TrackHISel1  = "track.isNonnull && track.quality('highPurity') && track.ptError/track.pt < 0.3"
TrackHISel2  = "track.isNonnull && track.quality('highPurity') && track.ptError/track.pt < 0.1 && track.numberOfValidHits > 10 && track.chi2/track.ndof/track.hitPattern.trackerLayersWithMeasurement < 0.15"

HighPtTriggerFlags = cms.PSet(
   HIL2Mu7_NHitQ10 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu7_NHitQ10_v*',1,0).empty()"),
   HIL3Mu7_NHitQ15 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu7_NHitQ15_v*',1,0).empty()"),
   HIL2Mu15 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*',1,0).empty()"),
   HIL3Mu15 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu15_v*',1,0).empty()"),
   HIL2Mu20 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu20_v*',1,0).empty()"),
   HIL3Mu20 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu20_v*',1,0).empty()"),
)

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(InAcceptance + ' && ' + TightId),
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        KinematicVariables,
	IsolationVariables,
	MuonIDVariables,
	MuonCaloVariables,
	TrackQualityVariables,
	GlobalTrackQualityVariables,
	StaOnlyVariables,
	dxyBS = cms.InputTag("muonDxyPVdzmin","dxyBS"),
	dxyPVdzmin = cms.InputTag("muonDxyPVdzmin","dxyPVdzmin"),
	dxyPVerr = cms.InputTag("muonDxyPVdzmin","dxyPVerr"),
	dzPV = cms.InputTag("muonDxyPVdzmin","dzPV"),
	dzPVerr = cms.InputTag("muonDxyPVdzmin","dzPVerr"),
	tkDxy = cms.InputTag("muonImpactParameter","Dxy"),
	tkDxyError = cms.InputTag("muonImpactParameter","DxyError"),
	tkDxyOverDxyError = cms.InputTag("muonImpactParameter","DxyOverDxyError"),
	tkDz = cms.InputTag("muonImpactParameter","Dz"),
	tkDzError = cms.InputTag("muonImpactParameter","DzError"),
	tkDzOverDzError = cms.InputTag("muonImpactParameter","DzOverDzError"),
	tkChi2OverNdofOverNlayer = cms.string("? track.isNull ? -1 : track.chi2/track.ndof/track.hitPattern.trackerLayersWithMeasurement"),
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       TightHI = cms.string(TightId),
       Track_HISel_v1  = cms.string(TrackHISel1),
       Track_HISel_v2  = cms.string(TrackHISel2),
    ),
    tagVariables = cms.PSet(
        KinematicVariables,
	IsolationVariables,
	MuonIDVariables,
	MuonCaloVariables,
	TrackQualityVariables,
	GlobalTrackQualityVariables,
	StaOnlyVariables,
        CentralityVariables,
	dxyBS = cms.InputTag("muonDxyPVdzminTags","dxyBS"),
	dxyPVdzmin = cms.InputTag("muonDxyPVdzminTags","dxyPVdzmin"),
	dxyPVerr = cms.InputTag("muonDxyPVdzminTags","dxyPVerr"),
	dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
	dzPVerr = cms.InputTag("muonDxyPVdzminTags","dzPVerr"),
	tkDxy = cms.InputTag("muonImpactParameterTags","Dxy"),
	tkDxyError = cms.InputTag("muonImpactParameterTags","DxyError"),
	tkDxyOverDxyError = cms.InputTag("muonImpactParameterTags","DxyOverDxyError"),
	tkDz = cms.InputTag("muonImpactParameterTags","Dz"),
	tkDzError = cms.InputTag("muonImpactParameterTags","DzError"),
	tkDzOverDzError = cms.InputTag("muonImpactParameterTags","DzOverDzError"),
    ),
    tagFlags = cms.PSet(
        HighPtTriggerFlags,
    ),
    pairVariables = cms.PSet(
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        probeMultiplicity = cms.InputTag("probeMultiplicity"),
        probeMultiplicity_TMGM = cms.InputTag("probeMultiplicityTMGM"),
        probeMultiplicity_Pt10_M60140 = cms.InputTag("probeMultiplicityPt10M60140"),
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.extraProbeVariablesSeq = cms.Sequence(
    process.splitTrackTagger +
    process.muonDxyPVdzmin + 
    process.muonImpactParameter
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.extraProbeVariablesSeq +
    process.probeMultiplicities + 
    process.bestPairByZMass + 
    process.muonDxyPVdzminTags +
    process.muonImpactParameterTags +
    process.centralityInfo +
    process.centralityBinInfo +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
      process.fastFilter
    * process.centralityBin
    * process.mergedMuons
    * process.patMuonsWithTriggerSequence
    * process.tnpSimpleSequence
)

##    _____               _    _             
##   |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##     | || '__/ _` |/ __| |/ / | '_ \ / _` |
##     | || | | (_| | (__|   <| | | | | (_| |
##     |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                                     |___/ 

## Then make another collection for standalone muons, using standalone track to define the 4-momentum
process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)
## Match to trigger, to measure the efficiency of HLT tracking
from PhysicsTools.PatAlgos.tools.helpers import *
process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(process, process.patMuonsWithTriggerSequence, "Sta")
process.muonMatchHLTL2Sta.maxDeltaR = 0.5
process.muonMatchHLTL3Sta.maxDeltaR = 0.5
massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta, "mergedMuons", "muonsSta")

## Define probes and T&P pairs
process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerSta"),
    cut = cms.string("outerTrack.isNonnull"), # no real cut now
)

process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = '40 < mass < 150')

process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))

process.pCutTracks.src = cms.InputTag("hiGeneralTracks")
process.staToTkMatch.maxDeltaR     = 0.3
process.staToTkMatch.maxDeltaPtRel = 2.
process.staToTkMatchNoZ.maxDeltaR     = 0.3
process.staToTkMatchNoZ.maxDeltaPtRel = 2.

process.load("MuonAnalysis.TagAndProbe.tracking_reco_info_cff")

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    arbitration   = "OneProbe",
    variables = cms.PSet(
        KinematicVariables, 
        StaOnlyVariables,
        ## track matching variables
        tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ","deltaR"),
        tk_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ","deltaEta"),
    ),
    flags = cms.PSet(
        outerValidHits = cms.string("outerTrack.numberOfValidHits > 0"),
        TM  = cms.string("isTrackerMuon"),
        Glb = cms.string("isGlobalMuon"),
        Tk  = cms.string("track.isNonnull"),
        StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
    ),
    tagVariables = cms.PSet(
        CentralityVariables,
        pt = cms.string("pt"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
        chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
        neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
        photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
        combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
    ),
    tagFlags = cms.PSet(
        HighPtTriggerFlags
    ),
    pairVariables = cms.PSet(
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
    ),
    pairFlags = cms.PSet(),
)

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta   +
    process.tpPairsSta      +
    process.onePairSta      +
    process.staToTkMatchSequenceZ +
    process.centralityInfo +
    process.centralityBinInfo +
    process.tpTreeSta
)

## Add extra RECO-level info
if True: 
    process.tracksNoMuonSeeded = cms.EDFilter("TrackSelector",
      src = cms.InputTag("hiGeneralTracks"),
      cut = cms.string(" || ".join("isAlgoInMask('%s')" % a for a in [
                    'initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep',
                    'mixedTripletStep', 'pixelLessStep', 'tobTecStep', 'jetCoreRegionalStep' ] ) )
    )
    process.pCutTracks0 = process.pCutTracks.clone(src = 'tracksNoMuonSeeded')
    process.tkTracks0 = process.tkTracks.clone(src = 'pCutTracks0')
    process.tkTracksNoZ0 = process.tkTracksNoZ.clone(src = 'tkTracks0')
    process.preTkMatchSequenceZ.replace(
            process.tkTracksNoZ, process.tkTracksNoZ +
            process.tracksNoMuonSeeded + process.pCutTracks0 + process.tkTracks0 + process.tkTracksNoZ0)
    process.staToTkMatch0 = process.staToTkMatch.clone(matched = 'tkTracks0')
    process.staToTkMatchNoZ0 = process.staToTkMatchNoZ.clone(matched = 'tkTracksNoZ0')
    process.staToTkMatchSequenceZ.replace( process.staToTkMatch, process.staToTkMatch + process.staToTkMatch0 )
    process.staToTkMatchSequenceZ.replace( process.staToTkMatchNoZ, process.staToTkMatchNoZ + process.staToTkMatchNoZ0 )
    process.tpTreeSta.variables.tk0_deltaR     = cms.InputTag("staToTkMatch0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta   = cms.InputTag("staToTkMatch0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ0","deltaEta")

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.centralityBin +
    process.muonsSta +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

##    _____     _          ____       _            
##   |  ___|_ _| | _____  |  _ \ __ _| |_ ___  ___ 
##   | |_ / _` | |/ / _ \ | |_) / _` | __/ _ \/ __|
##   |  _| (_| |   <  __/ |  _ < (_| | ||  __/\__ \
##   |_|  \__,_|_|\_\___| |_| \_\__,_|\__\___||___/
##                                                 
##   
process.load("MuonAnalysis.TagAndProbe.fakerate_all_cff")

process.fakeRateZPlusProbeTree = process.tpTree.clone(
    tagProbePairs = 'zPlusProbe',
    arbitration   = 'None', 
    tagVariables  = process.ZPlusProbeTagVariables,
    tagFlags      = cms.PSet(),
    pairVariables = cms.PSet(), 
    pairFlags     = cms.PSet(), 
)
process.fakeRateZPlusProbe = cms.Path(
    process.fastFilter +
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.zPlusProbeSequence +
    process.fakeRateZPlusProbeTree
)

process.schedule = cms.Schedule(
   process.tagAndProbe, 
   process.tagAndProbeSta, 
   process.fakeRateZPlusProbe
)

process.RandomNumberGeneratorService.tkTracksNoZ = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoZ0 = cms.PSet( initialSeed = cms.untracked.uint32(81) )

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_Data_PbPb_AOD.root"))
