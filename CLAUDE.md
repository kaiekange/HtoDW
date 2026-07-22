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
cmsRun EDAnalyzers/Reconstruction/test/SIG_H_phipimu.py        # DsW signal
cmsRun EDAnalyzers/Reconstruction/test/SIG_DsstW_H_phipimu.py  # DsstW signal
cmsRun EDAnalyzers/Reconstruction/test/BKG_H_phipimu.py        # background (doTruth=False)
```

Config files accept `JOBID`, `GROUPSIZE`, and `OUTFILE` arguments:
```bash
cmsRun EDAnalyzers/Reconstruction/test/SIG_H_phipimu.py JOBID=0 GROUPSIZE=1 OUTFILE=file:output.root
```

## HTCondor Submission

Each subdirectory under `Tuples/` corresponds to one sample. To submit:
```bash
cd Tuples/Signal/
python3 condor_submit.py
```

The `run.sh` in each `Tuples/` subdirectory is the HTCondor executable. It sources `cmsenv`, runs `cmsRun`, and writes output to `/eos/home-k/kakang/H2DW/Analysis/MC/2022/<JobType>/`.

## Post-processing ROOT Macros

Analysis macros in `Tuples/<JobType>/` (e.g., `bestcand.cc`, `compare.cc`, `PrepareSB.cc`) and `analysis/macros/` are standalone ROOT scripts run with `root -l -q macro.cc`.

## Architecture

### Physics Goal
Search for H→DsW, with **Ds± → φπ±, φ → K+K−** only. The non-resonant Ds→KKπ mode
was dropped on 2026-07-16 and its code deleted; recover it from git history if ever
needed. Also studies the Ds* variant. The decay chains are:
- H(125) → Ds±(→φπ± → K+K−π±) + W∓(→ℓ∓ν, ℓ=e/μ/τ)
- H(125) → Ds*±(→Ds±γ or Ds±π0(→γγ)) + W∓(→ℓ∓ν, ℓ=e/μ/τ)

Note the final state is the same K+K−π± either way; what makes this the φπ mode is
the φ requirement — the gen filter demands the φ ancestry and the candidate loop
applies a ±30 MeV window on M(KK).

**Reco/selection is μ-only regardless of the gen-level W flavor**: the analyzer always
requires a reco tight muon with pT>24 GeV before filling the tree (see Data Flow below).
So on disk, only prompt W→μ and W→τ→μ events ever produce rows; W→e and W→τ→(e/hadronic)
pass the gen filter (and increment `nEventsGen_`/`nEventsLepton_`) but never reach `Fill()`.
This was a deliberate 2026-07-22 decision — gen truth generalized, reco untouched.

### EDAnalyzers Package Structure

**`EDAnalyzers/JMEcorr/`** — jet/MET energy correction producers (must run before the main analyzer):
- `JecAppliedJetProducer`: applies JEC to slimmedJetsPuppi → `CorrectedAK4`
- `JecAppliedMetProducer`: applies JEC-propagated MET → `CorrectedMet`
- `JvmAppliedEventFilter`: jet veto map filter
- `tightJetId`: produces `ValueMap<Int_t>` tight jet ID

The cmsRun path is always: `jecAppliedJetProducer * jecAppliedMetProducer * jvmAppliedEventFilter * tightJetId * <MainAnalyzer>`

**`EDAnalyzers/Reconstruction/`** — the analyzer:
- `H_phipimu`: the only analyzer. Handles both channels via `channel = "DsW" | "DsstW"`,
  and both signal and background via `doTruth`. Full reconstruction including the
  kinematic vertex fit.

### Analyzer Data Flow (H_phipimu)

1. **Gen filter** (`doTruth` only): require exactly 1H + 1Ds(*) + 1W + 1 hard-process W-daughter lepton (e/μ/τ, any flavor — `Gen_W_lep_pdgId`), then find the status-1 gen K+, K−, π with the full φ ancestry chain (`TwoBodyDecay(gp, 333, 321, -321)` etc.). Returns early if not found. DsstW additionally requires the Ds* to go via γ **or** π0(→γγ) (`Gen_Dsst_mode`: 0=γ, 1=π0), both ancestry-checked the same way. If the W-lepton is a τ, a second pass over `isDirectHardProcessTauDecayProductFinalState()` particles classifies `Gen_tau_decay_mode` (0=hadronic, 1=e, 2=μ) and — only for τ→μ — supplies the muon used below. **The tree stores no gen-level kinematics** (no Gen_H/Ds/Dsst/W/Kp/Km/pi/mu/gamma/pi0/tau 4-vectors, no PV/SV, no Gen MET, no gen dR family) — only the categorical tags (`Gen_W_lep_pdgId`, `Gen_tau_decay_mode`, `Gen_Dsst_mode`) and the match booleans below. This was a deliberate 2026-07-22 simplification; the local `reco::GenParticle` objects needed for matching still exist in `analyze()`, just aren't written out.
2. **Gen matching**: match PF candidates via `packedPFCandidateToGenAssociation` + ΔR<0.05 + Δpt/pt<0.3 (ancestry checks accept either Ds* submode). The γ (DsstW, γ-mode) is matched geometrically instead — ΔR<0.02, Δpt/pt<0.5 — because it has no association; in π0-mode the two π0 daughter photons are matched the same way, producing `Gen_pi0gamma1_match`/`Gen_pi0gamma2_match`. Muon via `muon.genParticle()`; for a τ→μ event the ancestry check uses `isDirectHardProcessTauDecayProductFinalState()` instead of the direct W→μν chain.
3. **Trigger**: require `HLT_IsoMu24_v` or `HLT_Mu50_v`. Returns early if not fired.
4. **Primary vertex**: take `offlineSlimmedPrimaryVertices.front()`.
5. **Best muon**: highest-pT tight muon with pT>24 GeV.
6. **Ds candidate loop**: triple nested loop over the charged-hadron pool (`hasTrackDetails`, `trackHighPurity`, pT>1, |pdgId|=211). Per-track: K± pT>2, K+ and π opposite in charge to the muon, K− same as the muon. Pre-fit: (KKπ) pT>5, |M−MDs|<0.3. Kinematic vertex fit via `KinematicParticleVertexFitter`. Post-fit: vertex prob>0.001, |M−MDs|<0.15, pT>7, dira>0.99, (Ds+μ) mass<125, **φ: pT(KK)>3 and |M(KK)−Mφ|<0.03**.
7. **Ds* (DsstW only)**: photon pool is |pdgId|=22, pT>0.5; (Ds*+μ) mass<125.
8. **Output**: per-event scalars (gen info, PV, MET, jets, triggers, muon) + per-candidate vectors. `nCandidates` = number passing all cuts. Tree filled once per event.

### Tree Structure Convention
- `interface/<Name>_Tree.h` declares all branches as public members
- `src/<Name>_Tree.cc` implements `CreateBranches()` and `Init()` (reset to defaults)
- `null` sentinel = `-777` (defined as `#define null -777` in each tree header)
- Per-candidate quantities are `std::vector<T>`, indexed by candidate number
- Per-event quantities are scalars

### Key Constants
- GlobalTag (2022 MC): `130X_mcRun3_2022_realistic_v5`
- Ds mass: 1.96835 GeV, φ mass: 1.019460 GeV, K: 0.493677 GeV, π: 0.13957039 GeV, μ: 0.1056583755 GeV
- Input MiniAOD on XRootD: `root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsW_phipi/` (DsstW: `ggHtoDsstW_phipi/`; `output_188` there is bad and is skipped in the config)
- Output ntuples: `/eos/home-k/kakang/H2DW/Analysis/MC/2022/`

## Last Session (2026-07-22)

### H_phipimu gen truth generalized to e/μ/τ and Ds*→γ/π0 (gen-only; reco untouched)

新信号样本里 W 已经是 e/μ/τ inclusive（`24:onIfMatch` 三个都开着，フラグメント本来就这样，不是新改的），DsstW 的 Ds* 现在也按真实 BR 分两支：Ds*→Dsγ（94.2%）和 Ds*→Dsπ0(→γγ)（5.8%）。原来的代码只认 W→μ（`num_Gen_mu==1 && num_Gen_e==0`，τ 隐性被当场 return 掉）和 Ds*→Dsγ 单一模式。这次把 gen filter 泛化到能识别全部分支，**但 reco/选择完全没动**——用户明确要求只改 gen 部分。

**改了什么**（`H_phipimu.cc` 两个 channel 块 + `H_phipimu_Tree.h`/`.cc`）：
- `Gen_W_lep_pdgId`（±11/13/15，W 的直接 hard-process 带电轻子）；若是 τ，用 `isDirectHardProcessTauDecayProductFinalState()`（CMSSW 内建 flag，`GenParticle.h:85`）找它的直接衰变产物，得到 `Gen_tau_decay_mode`（0=hadronic/1=lep-e/2=lep-mu）
- PF/muon truth-match 的祖先检查按 τ 分支：τ→μ 时用上面那个 flag，不再走"muon 的 mother 必须是 W"那条链（原来这条链对 τ→μ 会静默漏判，是这次一起修的 bug，不是新引入的）
- DsstW：`Gen_Dsst_mode`（0=γ/1=π0）；K+/K−/π 的祖先检查两种模式都认（原来只认 γ 模式，会导致 π0 事件在 KKπ finder 那步就被吞掉）；π0 模式新增 `Gen_pi0gamma1_match`/`Gen_pi0gamma2_match`（PF dR<0.02 匹配两个光子）
- **Gen 树上不再存任何动力学量**（H/Ds/Ds*/W/K+/K−/π/μ/γ/π0/τ 的四矢量、PV/SV、Gen MET、dR 族、Dsmu/Dsstmu 全删）——这是中途用户明确要求的简化："只要知道这个事例是 mu 还是 e 还是 tau，然后能 truth-match，其他动力学数据可以删掉"。gen-finding 用到的本地 `reco::GenParticle` 对象（`Kp_gp`/`Km_gp`/`pi_gp`/`Ds_gp`/`mu_gp`/`gamma_gp`/`pi0gamma1_gp`/`pi0gamma2_gp` 等）全部保留，只删了写到 tree 上的那些行；truth-match 的产出仍然是原有的那批 Bool 分支（`Kp_match`…`event_match`）+ 新的两个 pi0gamma match。
- Ds→KKπ 的判定逻辑没变（φ→KK 本来就是 EvtGen 强制的 Br=1）

**关键：reco 端的 muon gate（`if (best_mu_idx == -1) return;`）在 Fill() 之前，一行没动**。所以 W→e、W→τ→e、W→τ→hadronic 的事件会在 gen filter 通过、`nEventsGen_`/`nEventsLepton_` 计数之后，在 reco 选择阶段被挡掉，不会写进 tree。只有 prompt W→μ 和 W→τ→μ 真正落盘。

**没验证的地方**（新 MiniAOD 还没产生，只能靠 CMSSW 文档保证语义）：
1. τ 作为 hard-process 粒子时 `num_Gen_Wlep==1` 会不会像 μ 一样干净地计一次——如果 Pythia 给 τ 弄出多个 hard-process copy，τ 支路会在计数 return 那一步被静默清零。样本出来后第一件事：画 `num_Gen_Wlep` 分布确认 τ 事件也稳定是 1。
2. `nEventsGen_`/`nEventsLepton_` 现在把 e/μ/τ 全部计进去了（以前这两个计数器只算 W→μ，e/τ 都被旧的 return 挡在外面）。如果下游拿它当 μ-channel 分母做归一化会被摊高，类似 [[project_gencount_normalization]] 那次的坑。**这个问题问了用户两次，用户都还没直接回答**——下次话题回到这里时记得追问，别让它再滑过去。

编译过（`-Werror=unused-but-set-variable` 逼着把改动过程中变成 dead code 的 `H_gp`/`Ds_gp`（DsstW 那份）/`Gen_H_idx`/`Gen_pi0_idx` 一起清掉了），没有跑 cmsRun（没有新 MiniAOD）。

**同一 session 后续又删了 GenJet/GenMET**（token、`consumes<>`、include、两处 fetch 全删）：`GenJetToken_` 在原版代码里本来就从没被 `getByToken` 过，是死代码，删它不改变任何行为；`GenMETToken_` 原来在两个 channel 里各有一次 `if (GenMETHandle->empty()) return;`，是个事件级 gate，这次跟着 fills 一起删掉了——**这是唯一一处真的改变了选择逻辑的删除**（虽然标准 MiniAOD 里 GenMET collection 几乎不可能是空的，实际影响预期为零，但如果哪天诡异地空过，行为跟以前不一样了）。

> 之前几次 session 的详细记录（geometric matching 研究、lostTracks recovery 研究、MVA v1/v2/v3 训练配比）已经归档进 memory（见 `[[project_geometric_matching]]`、`[[project_lostrack_recovery]]`、`[[project_phipimu_mva]]`），不再重复贴在这里。
