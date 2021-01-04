import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo",eras.Run2_2018_highBetaStar)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# the following lines are difined in EnergyLossProducer_cff.py ...be careful!
##process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff')
##process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
##process.load('RecoTracker.TrackProducer.TrackRefitter_cfi')
##TrackRefitter.src = 'generalTracks'

# V0 ...Luiz
##process.load('RecoVertex/V0Producer/generalV0Candidates_cff')
process.load('RecoVertex/V0Producer/generalV0Candidates_cfi')
process.generalV0Candidates.tkPtCut = cms.double(0.0)
####...the following needs TreckRefitters!
####process.generalV0Candidates.trackRecoAlgorithm = 'TrackRefitter'
##process.generalV0Candidates = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

####process.load("UserCode.EnergyLossPID.EnergyLossProducer_cff")
###############################################################################

from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag.globaltag = "101X_dataRun2_Express_v8"
process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v11"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(XXX)
#,
#lumisToProcess = cms.untracked.VLuminosityBlockRange('319104:1-319104:10',
#'319104:15-319104:185','319124:91-319124:277','319125:1-319125:208','319159:125-319159:618',
#'319174:1-319174:77','319175:1-319175:139','319176:1-319176:1803','319177:1-319177:232',
#'319190:1-319190:317','319222:108-319222:294','319223:1-319223:131','319254:115-319254:263',
#'319255:1-319255:164','319256:1-319256:726','319262:10-319262:10','319262:15-319262:16',
#'319262:20-319262:23','319262:29-319262:34','319262:39-319262:40','319262:46-319262:58',
#'319262:61-319262:78','319262:82-319262:123','319262:129-319262:362','319263:1-319263:367',
#'319264:1-319264:57','319265:1-319265:396','319266:1-319266:26','319267:1-319267:204',
#'319268:1-319268:467','319270:1-319270:206','319300:1-319300:1132','319311:1-319311:1733'
#)
)

process.options = cms.untracked.PSet(
)

process.demo = cms.EDAnalyzer('PromptAnalyzer'
# Ferenc
    ,tracks = cms.InputTag('generalTracks')
####    ,tracks = cms.InputTag('TrackRefitter') ??
####    ,tracks = cms.InputTag('refitterForEnergyLoss')
#
    ,dedxs = cms.InputTag('dedxHarmonic2')
    ,dedxPIXs = cms.InputTag('dedxPixelHarmonic2')
# Ferenc
####    ,dedxPIXs = cms.InputTag('energyLossProducer','energyLossAllHits')
#
    ,dedxpixels = cms.InputTag('dedxHitInfo')
    ,RPtracks = cms.InputTag('ctppsLocalTrackLiteProducer')
    ,vertices = cms.InputTag('offlinePrimaryVertices')
    ,beamspot = cms.InputTag('offlineBeamSpot')
    ,triggers = cms.InputTag('TriggerResults','','HLT')
###4    ,pflows = cms.InputTag('particleFlow')
###5    ,muons = cms.InputTag('muons')
#     ,clusters = cms.InputTag('siPixelClusters')
# Luiz
    ,kshorts = cms.InputTag('generalV0Candidates','Kshort')
    ,lambdas = cms.InputTag('generalV0Candidates','Lambda')                      
)

#
##process.TFileService = cms.Service("TFileService",
##    fileName = cms.string('data_YYY.root')
##)
# ...Luiz
process.TFileService = cms.Service("TFileService",
            fileName = cms.string("output.root"),
            closeFileFast = cms.untracked.bool(False)
)

# ...Luiz
# Configure the object that writes an output file
##process.out = cms.OutputModule("PoolOutputModule",
##    fileName = cms.untracked.string("output.root")
##)

# produce new dEdx by Ferenc
#
####process.energyLossProducer.tag = cms.string('totem')

####process.reco = cms.Path(process.MeasurementTrackerEvent
####                      * process.siPixelClusterShapeCache)

# turned off : creates huge .err file
####process.Tracer = cms.Service("Tracer")

process.p = cms.Path( process.generalV0Candidates * process.demo )
####process.p = cms.Path( process.TrackRefitter * process.generalV0Candidates * process.demo )
####process.p = cms.Path( process.demo )

# ...Luiz
##process.ep = cms.EndPath(process.out)


####process.schedule = cms.Schedule(process.p)
##process.schedule = cms.Schedule(process.p,process.ep)
####process.schedule = cms.Schedule(process.reco,process.produceEnergyLoss,process.p)
#####process.schedule = cms.Schedule(process.generalV0Candidates,process.p)

##process.schedule = cms.Schedule(process.reco,process.produceEnergyLoss,process.p,process.ep)

#do not add changes to your config after this point (unless you know what you are doing)
##from FWCore.ParameterSet.Utilities import convertToUnscheduled
##process=convertToUnscheduled(process)

# copied from FErenc
# Add early deletion of temporary data products to reduce peak memory need
##from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
##process = customiseEarlyDelete(process)
