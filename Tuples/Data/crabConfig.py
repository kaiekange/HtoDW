# Template -- placeholders are substituted by submitData.sh into a scratch crabConfig.py.
# Do not edit the generated crabConfig.py; it is overwritten on every submission.
#
# Placeholders are fenced in double underscores rather than written as bare uppercase
# words: the pset's option is itself named GLOBALTAG, so a bare 'GLOBALTAG=GLOBALTAG'
# would sed into '130X_dataRun3_v2=130X_dataRun3_v2'. The fencing also rules out one
# placeholder matching inside another. submitData.sh rejects any generated config that
# still contains a fenced token, so keep the fence out of prose in this file.
#
# Every setting other than the placeholders is copied verbatim from the
# crabConfig_Run2022{C,D}{,_SS}.py files that produced the 2026-08-03 data, so the new
# eras run under the same splitting, lumi mask and storage as the old ones.
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'Muon_Run2022D_H_phipimu_SS'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/DATA_H_phipimu_SS.py'
config.JobType.pluginName = 'Analysis'
# Per-dataset era, straight from DatasetList/data_2022.txt. The pset has no default for
# either and raises at parse time if one is missing or unknown.
config.JobType.pyCfgParams = ['GLOBALTAG=130X_dataRun3_v2', 'YEAR=2022Pre']
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
config.Data.inputDataset = '/Muon/Run2022D-22Sep2023-v1/MINIAOD'
config.Data.outputDatasetTag = 'Muon_Run2022D_H_phipimu_SS'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
# Lands on EOS as
#   /eos/home-k/kakang/H2DW/Analysis/Data/<era>/Muon/<outputDatasetTag>/<ts>/0000/*.root
# i.e. a Data/ sibling of the existing Analysis/{MC,MVA}, using the same 2022 / 2022EE era
# names as the MC tree.
#
# The era is in the path, not just a comment: all background MC is per-era and every
# Lumi_reweight is normalised to that era's luminosity, so preEE and postEE data must never
# be merged. submitData.sh fills 2022 from the Year column of the dataset list
# (2022Pre -> 2022, 2022Post -> 2022EE), so re-enabling the Run2022C/D rows there is all it
# takes to produce preEE into its own directory.
config.Data.outLFNDirBase = '/store/user/kakang/H2DW/Analysis/Data/2022'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
# CERNBOX rather than IIHE, so the data tuples sit next to the MC on EOS.
# T3_CH_CERNBOX maps outLFNDirBase's /store/user/kakang/... onto /eos/user/k/kakang/...,
# which is the same inode as /eos/home-k/kakang/... where the MC already lives.
config.Site.storageSite = 'T3_CH_CERNBOX'
