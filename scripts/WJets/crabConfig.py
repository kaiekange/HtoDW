from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MC2017UL_RecoBKG_WJets_20250630'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.psetName = '/user/kakang/Analysis/CMSSW_10_6_30_patch1/src/EDAnalyzers/RecoAnalyzer/test/RecoBestAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
# config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting='FileBased'
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 10 
config.Data.publication = False
config.Data.inputDataset = '/WJetsToLNu_012JetsNLO_34JetsLO_EWNLOcorr_13TeV-sherpa/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v4/MINIAODSIM'

config.Data.outputDatasetTag = 'MC2017UL_RecoBKG_WJets'
config.Data.outLFNDirBase = '/store/user/kakang/Analysis/Reconstruction/'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'

config.section_("User")
config.User.voGroup = 'becms'

config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
