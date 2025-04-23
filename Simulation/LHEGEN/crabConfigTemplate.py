from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'REQUESTNAME'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.psetName = './LHEGEN_cfg.py'
config.JobType.pluginName = 'PrivateMC'

config.section_("Data")
config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting='EventBased'
config.Data.unitsPerJob = UNITSPERJOB
config.Data.totalUnits = TOTALUNITS
config.Data.publication = False
config.Data.outputDatasetTag = 'OUTPUTDATASETTAG'
config.Data.outLFNDirBase = 'OUTLFNDIRBASE'

config.section_("User")
config.User.voGroup = 'becms'

config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
