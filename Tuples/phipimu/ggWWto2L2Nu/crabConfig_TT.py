from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'ggWWto2L2Nu_TT_H_phipimu'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/BKG_H_phipimu.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
##config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2750 #3000 # min

config.section_('Data')
config.Data.splitting='FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 10

#config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.inputDataset = '/ggWWto2L2Nu_OS_PolarizationTT_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'ggWWto2L2Nu_TT_H_phipimu'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = '/store/user/kakang/H2DW/Analysis/MC/2022'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
config.Site.storageSite = 'T2_BE_IIHE'
