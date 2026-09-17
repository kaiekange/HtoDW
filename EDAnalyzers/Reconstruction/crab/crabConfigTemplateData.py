# Template -- placeholders are substituted by submitData.sh into a scratch crabConfig.py.
# Do not edit the generated crabConfig.py; it is overwritten on every submission.
#
# Placeholders are fenced in double underscores (matching the old Tuples/Data/ template):
# the pset's option is itself named GLOBALTAG, so a bare 'GLOBALTAG=GLOBALTAG' would sed
# into '130X_dataRun3_v2=130X_dataRun3_v2'. submitData.sh rejects any generated config
# that still contains a fenced token.
#
# Splitting/lumi-mask/storage-site settings are copied verbatim from
# Tuples/Data/crabConfigTemplateData.py (the old H_phipimu-tree production) -- only the
# pset, output path and dataset table changed for the NanoAOD workflow.
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = '__REQUESTNAME__'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '__PSETNAME__'
config.JobType.pluginName = 'Analysis'
# Per-dataset GlobalTag, from Tuples/Data/DatasetList/data_2022.txt (shared with the old
# workflow -- it's a pure DAS/GT table, not release-specific). SS selects the OS/SS
# tuplizer variant. The pset has no default for GLOBALTAG and raises at parse time if
# it's missing.
# YEAR (2022Pre/2022Post, the Year column) selects the era ParticleNet WPs; required by the pset.
config.JobType.pyCfgParams = ['GLOBALTAG=__GLOBALTAG_VALUE__', 'SS=__SS_VALUE__', 'YEAR=__YEAR_VALUE__']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2750  # min
# 2026-09-15: the C/D tasks ran with the CRAB default 2000 MB and jobs were killed at
# 2250-2500 MB (exit 50660), about half of Run2022C at T2_FR_GRIF. Completed D jobs peaked
# at 2215 MB.
# 3000 is the CRAB maximum for a single-core job (3500 is rejected at submission).
config.JobType.maxMemoryMB = 3000

config.section_('Data')
config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
# 150 lumis/job (inherited from the old tuplizer) gave 3 h average and 6.4 h longest jobs on
# Run2022D; E/F/G have more events per lumi. 50 keeps single jobs short.
config.Data.unitsPerJob = 50
config.Data.lumiMask = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/crab/JSON/Cert_Collisions2022_355100_362760_Golden.json'

config.Data.publication = False
config.Data.inputDataset = '__INPUTDATASET__'
config.Data.outputDatasetTag = '__OUTPUTDATASETTAG__'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
# Lands on EOS as
#   /eos/home-k/kakang/H2DW/Analysis/NanoMVA/Data/<era>/Muon/<outputDatasetTag>/<ts>/0000/*.root
# -- a Data/ sibling of NanoMVA/MC (see condor_submit.py's EOS_BASE), same 2022/2022EE era
# names. The era is in the path, not just a comment: every Lumi_reweight downstream is
# normalised per-era, so preEE and postEE data must never land in the same directory.
# submitData.sh fills __ERADIR__ from the Year column of the dataset list.
config.Data.outLFNDirBase = '/store/user/kakang/H2DW/Analysis/NanoMVA/Data/__ERADIR__'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
config.Site.storageSite = 'T3_CH_CERNBOX'
