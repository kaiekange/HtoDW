from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'Muon_Run2022C_H_phipimu_SS'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/DATA_H_phipimu_SS.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
##config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2750 #3000 # min

config.section_('Data')
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 150
config.Data.lumiMask = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/Data/JSON/Cert_Collisions2022_355100_362760_Golden.json'

#config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.inputDataset = '/Muon/Run2022C-22Sep2023-v1/MINIAOD'
config.Data.outputDatasetTag = 'Muon_Run2022C_H_phipimu_SS'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = '/store/user/kakang/H2DW/Analysis/Data/2022'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
config.Site.storageSite = 'T2_BE_IIHE'
