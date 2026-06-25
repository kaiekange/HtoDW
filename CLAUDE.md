# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working Convention

At the end of each session, update the "Last Session" section below.

## Environment Setup

This is a CMSSW 15.0.5 package at CERN. Before any build or run command, the CMSSW environment must be active:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/
cmsenv
```

## Build

```bash
scram b -j8
```

To clean and rebuild a single package:
```bash
scram b -j8 EDAnalyzers/Reconstruction
```

## Running Analyzers Locally

```bash
cmsRun EDAnalyzers/Reconstruction/test/RecoAll_DsW_KKpimu.py
cmsRun EDAnalyzers/GenAnalyzer/test/GenMatch_DsW_KKpimu.py
```

Config files accept `JOBID`, `GROUPSIZE`, and `OUTFILE` arguments:
```bash
cmsRun EDAnalyzers/Reconstruction/test/RecoAll_DsW_KKpimu.py JOBID=0 GROUPSIZE=1 OUTFILE=file:output.root
```

## HTCondor Submission

Each subdirectory under `Tuples/` corresponds to one analyzer job type. To submit:
```bash
cd Tuples/RecoAll_DsW_KKpimu/
python3 condor_submit.py
```

The `run.sh` in each `Tuples/` subdirectory is the HTCondor executable. It sources `cmsenv`, runs `cmsRun`, and writes output to `/eos/home-k/kakang/H2DW/Analysis/MC/2022/<JobType>/`.

## Post-processing ROOT Macros

Analysis macros in `Tuples/<JobType>/` (e.g., `bestcand.cc`, `compare.cc`, `PrepareSB.cc`) and `analysis/macros/` are standalone ROOT scripts run with `root -l -q macro.cc`.

## Architecture

### Physics Goal
Search for H→DsW→(KKπ)(μν). Also studies the Ds* variant. The decay chain is:
- H(125) → Ds±(→K+K-π±) + W∓(→μ∓ν)

### EDAnalyzers Package Structure

**`EDAnalyzers/GenAnalyzer/`** — gen-level and gen-reco matching analyzers:
- `GenMatch_DsW_KKpimu` / `GenMatch_DsstW_KKpimu`: gen-level particle matching with PF candidates
- `GenRecoDsW` / `GenRecoDsstW`: gen+reco combined ntuples
- `TrigAnalyzer`, `Count_DsW_KKpimu`: utility analyzers

**`EDAnalyzers/JMEcorr/`** — jet/MET energy correction producers (must run before the main analyzer):
- `JecAppliedJetProducer`: applies JEC to slimmedJetsPuppi → `CorrectedAK4`
- `JecAppliedMetProducer`: applies JEC-propagated MET → `CorrectedMet`
- `JvmAppliedEventFilter`: jet veto map filter
- `tightJetId`: produces `ValueMap<Int_t>` tight jet ID

The cmsRun path is always: `jecAppliedJetProducer * jecAppliedMetProducer * jvmAppliedEventFilter * tightJetId * <MainAnalyzer>`

**`EDAnalyzers/Reconstruction/`** — full reco analyzers:
- `RecoAll_DsW_KKpimu` / `RecoAll_DsstW_KKpimu`: full event reconstruction including kinematic vertex fit

### Analyzer Data Flow (RecoAll_DsW_KKpimu)

1. **Gen filter**: require exactly 1H + 1Ds + 1W + 1μ + 0e (hard process), then find gen K+, K-, π daughters. Returns early if not found.
2. **Gen matching**: match PF candidates to gen particles via `packedPFCandidateToGenAssociation` + ΔR<0.05 + Δpt/pt<0.3. Match muon via `muon.genParticle()`.
3. **Trigger**: require `HLT_IsoMu24_v` or `HLT_Mu50_v`. Returns early if not fired.
4. **Primary vertex**: take `offlineSlimmedPrimaryVertices.front()`.
5. **Best muon**: highest-pT tight muon with pT>24 GeV.
6. **Ds candidate loop**: triple nested loop over PF candidates (charged hadrons with `hasTrackDetails`, `trackHighPurity`, pT>1 GeV). Pre-cuts on Ds pT>5 and |M-MDs|<0.3. Kinematic vertex fit via `KinematicParticleVertexFitter`. Post-fit cuts: |M-MDs|<0.15, pT>7, dira>0.99, (Ds+μ) mass<125.
7. **Output**: per-event scalars (gen info, PV, MET, triggers, muon) + per-Ds-candidate vectors. `nCandidates` = number of candidates passing all cuts. Tree filled once per event.

### Tree Structure Convention
- `interface/<Name>_Tree.h` declares all branches as public members
- `src/<Name>_Tree.cc` implements `CreateBranches()` and `Init()` (reset to defaults)
- `null` sentinel = `-777` (defined as `#define null -777` in each tree header)
- Per-candidate quantities are `std::vector<T>`, indexed by candidate number
- Per-event quantities are scalars

### Key Constants
- GlobalTag (2022 MC): `130X_mcRun3_2022_realistic_v5`
- Ds mass: 1.96835 GeV, φ mass: 1.019460 GeV, K: 0.493677 GeV, π: 0.13957039 GeV, μ: 0.1056583755 GeV
- Input MiniAOD on XRootD: `root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsW_KKpimu/`
- Output ntuples: `/eos/home-k/kakang/H2DW/Analysis/MC/2022/`

## Last Session (2026-05-26)

### DNN_Dsmu (v1 重新训练完成)
- 路径：`analysis/DNN_Dsmu/`；模型 128→64→32→16→1
- 特征：61→46（移除 9 个 `_p` 冗余变量、6 个 `_ipxysig`/`FDxysig`）；新增 dira sentinel 过滤（`> -1`）
- 产物：`/eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/v1/Prediction/` 下的 scored/best 文件
- 性能：多候选 Ds 选对率 DsstW 85.3%，DsW 83.7%（vs baseline Ds_pt 60.1%/58.8%，提升 +25%）

### DNN_gamma (v1 训练完成)
- 路径：`analysis/DNN_gamma/`；模型 **32→16→8→1**（缩小）
- 特征：15→13（移除 dM 防 mass sculpting，移除 Dsstmu_p 冗余）
- 输入文件从 `MVA/Dsmu/v1/Prediction/` 读取（硬编码 v1）
- 产物：`/eos/home-k/kakang/H2DW/Analysis/MVA/gamma/v1/Prediction/` 下的 best_scored/best_gamma 文件
- 性能：有真光子时选中率 95.9%（vs baseline 最高pT 94.0%）；仅 13.4% DsstW 事例有真光子入选
- run.sh 包含全流程（训练+predict+best_cand，3个样本串行）
