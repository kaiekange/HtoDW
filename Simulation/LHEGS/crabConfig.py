from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = "MC_ggHtoDsW_KKpimu_2022_20260408"
config.General.transferLogs = True

config.section_("JobType")
config.JobType.psetName = "LHEGS_ggHtoDsW_KKpimu_cfg.py"
config.JobType.pluginName = "PrivateMC"

config.section_("Data")
config.Data.outputPrimaryDataset = "LHEGS"
config.Data.splitting="EventBased"
config.Data.unitsPerJob = 500
config.Data.totalUnits = 500000
config.Data.publication = False
config.Data.outputDatasetTag = "2022"
config.Data.outLFNDirBase = "/store/user/kakang/Analysis/ggHtoDsW_KKpimu/20260408/Simulation"

config.section_("User")
config.User.voGroup = "becms"

config.section_("Site")
config.Site.storageSite = "T2_BE_IIHE"