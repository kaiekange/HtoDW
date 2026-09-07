# Inclusive WtoLNu-2Jets, ERA=2022EE, opposite-sign.
#
# One of four near-identical configs (2022|2022EE x OS|SS) that differ only in
# inputDataset, pyCfgParams and the output names. Everything else is deliberately
# identical so the four tasks are comparable.
#
# The pset is the CRAB sibling of WtoLNu2Jets_H_phipimu.py: its source is empty and CRAB
# fills fileNames from Data.inputDataset. Era drives BOTH the GlobalTag and the JEC/veto-map
# year inside the pset -- see its ERAS dict; there is no default, so a missing ERA fails at
# parse time rather than silently producing preEE numbers.
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'WtoLNu2Jets_2022EE_H_phipimu'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/EDAnalyzers/Reconstruction/test/WtoLNu2Jets_crab_H_phipimu.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['ERA=2022EE', 'SS=0']
config.JobType.allowUndistributedCMSSW = True
##config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2750 #3000 # min

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
# 10 files/job, inherited from the never-submitted crabConfig.py in this directory.
# Files average ~3.8 GB, so a job reads ~38 GB. UNVERIFIED at this size -- if jobs hit
# maxJobRuntimeMin, lower this and resubmit.
config.Data.unitsPerJob = 10

#config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.inputDataset = '/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'
config.Data.outputDatasetTag = 'WtoLNu2Jets_2022EE_H_phipimu'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
# Lands on EOS as
#   /eos/home-k/kakang/H2DW/Analysis/MC/2022EE/<primaryDS>/WtoLNu2Jets_2022EE_H_phipimu/<ts>/0000/*.root
# Note the extra <primaryDS> level compared with the HTCondor jet-binned layout
# (.../MC/<era>/WtoLNu2Jets_<bin>_H_phipimu/), which CRAB always inserts.
#
# Era is in the path because every Lumi_reweight is normalised to that era's luminosity
# (7990 /pb preEE vs 26680 /pb postEE); preEE and postEE must never be merged.
config.Data.outLFNDirBase = '/store/user/kakang/H2DW/Analysis/MC/2022EE'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
# T3_CH_CERNBOX maps /store/user/kakang/... onto /eos/user/k/kakang/..., the same inode as
# /eos/home-k/kakang/... where the rest of the analysis lives -- so the output is readable
# with POSIX paths and needs no xrootd for hadd/dedup.
config.Site.storageSite = 'T3_CH_CERNBOX'
