import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
from Configuration.AlCa.GlobalTag import GlobalTag
import sys

process = cms.Process("PVAnalysis")

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag =  GlobalTag(process.GlobalTag, "106X_mc2017_realistic_v9")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) ) 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

# myinfile = "file:/user/kakang/Analysis/CMSSW_10_6_30_patch1/src/tuples/PAT/output_53.root"
# myoutfile = "file:test.root"
# myinfile = sys.argv[2]
# myoutfile = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
    # fileNames = cms.untracked.vstring(myinfile)
)

process.RecoAnalyzer = cms.EDAnalyzer('RecoAnalyzer',
    pfCands = cms.InputTag('packedPFCandidates'),
    muons = cms.InputTag('slimmedMuons'),
    primvtx = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamspot = cms.InputTag("offlineBeamSpot"),

    TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        minPt = cms.double(0.5),
        maxEta = cms.double(5.0), # 2.4
        maxD0Significance = cms.double(4.0),
        maxNormalizedChi2 = cms.double(10.0),
        minSiliconLayersWithHits = cms.int32(5),
        minPixelLayersWithHits = cms.int32(2),
        trackQuality = cms.string("any")
    ),

    TkClusParameters = cms.PSet(
        algorithm = cms.string("DA_vect"),
        TkDAClusParameters = cms.PSet(coolingFactor = cms.double(0.6),  #  moderate annealing speed
            Tmin = cms.double(2.0),            #  end of annealing
            Tpurge = cms.double(2.0),         # cleaning
            Tstop = cms.double(0.5),          # end of annealing
            uniquetrkweight = cms.double(0.8), # require at least two tracks with this weight at T=Tpurge
            zmerge = cms.double(1e-2),        # merge intermediat clusters separated by less than zmerge
            vertexSize = cms.double(0.006),    #  ~ resolution / sqrt(Tmin)
            d0CutOff = cms.double(3.),        # downweight high IP tracks 
            dzCutOff = cms.double(3.)         # outlier rejection after freeze-out (T<Tmin)
        )
    ),

    VxFitterParameters = cms.PSet(algorithm=cms.string('AVF'),
        minNdof = cms.double(0.0),
        maxDistanceToBeam = cms.double(1.0)
    ),

    VxFitterBSParameters = cms.PSet(algorithm=cms.string('AVFBS'),
        minNdof = cms.double(2.0),
        maxDistanceToBeam = cms.double(1.0)
    ),
)

process.TFileService = cms.Service(
    # "TFileService", fileName=cms.string(myoutfile), closeFileFast = cms.untracked.bool(True)
    "TFileService", fileName=cms.string("output.root"), closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.RecoAnalyzer)
