from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'H2DW_BKG_WJets'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/test/FullReco.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
##config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 3000 # min

config.section_('Data')
config.Data.splitting='FileBased'
config.Data.totalUnits = -1
#config.Data.unitsPerJob = 1
#config.Data.unitsPerJob = 100
config.Data.unitsPerJob = 20

#config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.inputDataset = '/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'Run3Summer22_MINIAODSIM'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = '/store/user/kakang/H2DW/Datasets/20260119/Background'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'
