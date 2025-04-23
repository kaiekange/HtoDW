from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'MC_H2DW_2017UL_20250417'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.psetName = './LHEGEN_cfg.py'
config.JobType.pluginName = 'PrivateMC'

config.section_("Data")
config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting='EventBased'
config.Data.unitsPerJob = 200
config.Data.totalUnits = 100000
config.Data.publication = False
config.Data.outputDatasetTag = 'H2DW_2017UL_20250417'
config.Data.outLFNDirBase = '/store/user/kakang/Analysis/Simulation/'

config.section_("User")
config.User.voGroup = 'becms'

config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
