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
cd Tuples/HtoDsW/
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
requires a reco tight muon with pT>26 GeV before filling the tree (see Data Flow below).
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
5. **Best muon**: highest-pT tight muon with pT>26 GeV (raised from 24 on 2026-07-22 — above the `HLT_IsoMu24` turn-on so the offline selection sits on the trigger efficiency plateau).
6. **Ds candidate loop**: triple nested loop over the charged-hadron pool (`hasTrackDetails`, `trackHighPurity`, pT>1, |pdgId|=211). Per-track: K± pT>2, K+ and π opposite in charge to the muon, K− same as the muon. Pre-fit: (KKπ) pT>5, |M−MDs|<0.5; **φ pre-fit** M(KK)<Mφ+0.08 (loose, pre-vertex-fit speed cut on raw KK mass). Kinematic vertex fit via `KinematicParticleVertexFitter`. Post-fit: vertex prob>0.001, |M−MDs|<0.3, pT>7, dira>0.99, (Ds+μ) mass<125, **φ: pT(KK)>3 and M(KK)<Mφ+0.05** (upper-cut only — the KK threshold, Mφ−32 MeV, physically bounds the lower side, so no lower cut is applied). Ds/φ windows widened 2026-07-22 (post-fit Ds ±0.15→±0.3, pre-fit KKπ ±0.3→±0.5) to give the Ds-mass sideband more lever arm for an SR/CR split; the pre-fit window is intentionally wider than post-fit so it never truncates the sideband.
7. **Ds* (DsstW only)**: photon pool is |pdgId|=22, pT>0.5; (Ds*+μ) mass<125. Mass window unchanged (±0.15) — not yet revisited pending a DsstW rerun to look at its resolution.
8. **Output**: the tree is row-per-candidate (flat ntuple), not one row per event — `tree->Fill()` is called once per Ds×gamma combination (or once per Ds if no gamma passes), so an event with N candidates writes N rows. Event-level fields (gen info, PV, MET, `jet_ht`/`jet_leadpt`/`jet_maxbscore`, triggers, muon, electron/muon os/ss lepton counts) are set once before the candidate loop and repeat identically across that event's rows. `Ds_nearestjet_dR`/`bscore` and `n_b_jets_{loose,medium,tight}{,_pt20}_dR{03,04}` (b-tagged central jets excluding jets within dR03/dR04 of the Ds, added 2026-07-23) are genuinely per-candidate: recomputed inside the loop from each row's own `Ds_mtm`. `nCandidates` = number of rows passing all cuts.

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

## Last Session (2026-07-23)

### Same-sign (SS) control-region tuplizer

新增 `H_phipimu_SS`（`plugins/H_phipimu_SS.cc`），跟 `H_phipimu` 逐字一样，唯一 delta 是 K+K−π 候选循环里 pi 的电荷 cut 从 `!= mu.charge()`（要求反号，即真信号的电荷关联）翻成 `== mu.charge()` 排除、也就是**只要求 pi 与 mu 同号**——因为 Ds 和 W 电荷相反，真信号 pi 必然跟 mu 反号，所以 SS 样本在结构上排除了真信号，是纯 combinatoric/wrong-charge 本底的干净 control region。Kp/Km 的 cut 没动（它们只是"哪颗 K 叫 Kp/Km"的记账，跟 SS/OS 选择正交，不影响 phase space）。诊断用的 `num_os_mu_*`/`num_os_e_*` lepton 计数公式也跟着翻了一下，保持它们"相对 pi 电荷"的语义在 OS/SS 两个文件里一致（不翻的话同名分支在两个文件里会代表相反的量）。

在旧信号 MC 上验证过：`Ds_match`/`event_match` ≈0%（SS 确实排除了真信号，唯一一个 match 是电荷误重建的 per-mille 量级噪声，不是 bug）；`mu_match` 跟 OS 基本一样（~99.7-99.9% vs 99.98-99.99%，印证了 tight muon 选择不受 pi 电荷翻转影响，fake 只在强子侧）；`phi_match`~41.5%，说明假候选里接近一半是"真 φ+错电荷 π"。**这次验证是信号 MC 的自污染检查，不等于 OS/SS 本底对称性验证**——后者需要拿背景 MC 分别跑 OS/SS 两个 tuplizer，比较产额比例和 M(KKπ) 等拟合变量的 shape，这一步还没做。

配套建了 `test/{BKG,SIG,SIG_DsstW}_H_phipimu_SS.py`（跟各自 OS 版本一样，只换 `EDAnalyzer` 类型）、`Tuples/HtoDsW/{condor_submit_SS.py,run_SS.sh}`、`Tuples/HtoDsstW/{condor_submit_SS.py,run_SS.sh}`、`Tuples/WtoLNu2Jets/crabConfig_SS.py`——都是复制现有 OS 版本再改一两行，不是参数化开关（跟这个仓库一直以来"cp 原版再加 delta"的规矩一致）。

**踩坑记录**：新建的 plugin 编译成功后本地 cmsRun 报 `PluginNotFound`，原因是 `$CMSSW_BASE/lib/<arch>/.edmplugincache` 在 `scram b` 之后没有自动收录新插件；修法是 `edmPluginRefresh -p $CMSSW_BASE/lib/el9_amd64_gcc12`。

### Tuples 目录改名

`Tuples/Signal`→`Tuples/HtoDsW`，`Tuples/Signal_DsstW`→`Tuples/HtoDsstW`，`Tuples/WJets`→`Tuples/WtoLNu2Jets`——按 dataset 名字去掉连字符的规矩（跟 `TbarQtoLNu`/`TTto2L2Nu` 等现有目录的命名风格对齐，不加下划线）。`CLAUDE.md` 里的示例路径同步改了。

### 新信号样本 (20260717) 上线，处理还没跑完的问题

发现新的 DsW 信号样本在 `root://maite.iihe.ac.be//store/user/kakang/Analysis/ggHtoDsW_phipi/20260717/`——这才是真正的 inclusive W→e/μ/τ + Ds*→γ/π0 fragment（`[[project_h_phipimu_multiflavor_gen]]` 记录的旧样本 `20260624` 从未测试过 τ/π0 代码路径，这批新样本才是能测的）。生产还在跑（condor batch `MC_ggHtoDsW_phipi_AODSIM_2022_20260717`），预计 500 个文件，目前只有部分落盘且中间有零散缺口。

`test/SIG_H_phipimu.py`/`SIG_H_phipimu_SS.py` 的 `INFILE` 改成运行时现查 `xrdfs maite.iihe.ac.be ls <PAT目录>`，只把真实存在的文件塞进 `PoolSource`，不存在的自动跳过——每次重跑都会现查，以后生产陆续补齐也会自动捡到新文件，不用手动改排除列表（20260624 那行保留作注释参考）。

**发现并修了一个资源浪费的 bug**：`Tuples/HtoDsW/condor_submit.py` 的 `njobs = 1000/GroupSize` 是照抄旧样本（~1000 个文件）的规模，但新样本上限只有 500，导致 JOBID 25-49 那一半 job（每批 50 个里的 25 个，OS+SS 两批共 50 个空 job）分到的文件号根本不存在，产出全是 16967 字节的空 ROOT 文件（有 tree header 但 0 entries，不是崩溃，`cmsRun` exit 0）。修法是 `njobs` 也现查 `xrdfs` 拿到的最大文件编号动态算。顺便把 `GroupSize` 从 20 改成 10，缩短单个 job 的跑时。

### H_phipimu / H_phipimu_SS：n_b_jets 改成按候选算、排除 Ds 自己所在的 jet

`n_b_jets_{loose,medium,tight}{,_pt20}` 原来是候选循环之前算一次的事例级标量，跟具体哪个 Ds 候选无关。现在搬到 K+K−π 循环内部、`Ds_mtm` 算出来之后重新算，新增按 ΔR(jet,Ds) 门槛区分的两套（`_dR03`/`_dR04`，各自排除跟这一行 Ds 距离在 0.3/0.4 以内的 jet），一共 12 个新分支替换掉原来的 6 个——避免把"Ds 自己衰变径迹聚出来的那个 jet"错当成额外的 b-jet 计数进去。两个 plugin（OS/SS）改法完全同步。这次改动只编译过，没有跑 cmsRun 验证。

### 顺带发现（未处理）：数据侧的 dataset

如果之后要跑 data，用 `dasgoclient` 查过，2022 Run3 单/双 muon 触发路径统一在 `/Muon/Run2022<era>-<reprocessing>/MINIAOD`；现有信号 MC 是 `Run3Summer22MiniAODv4`（preEE），对应数据应选 `Run2022C`+`Run2022D`（`22Sep2023` 版 ReReco），不要混进 postEE 的 E/F/G。Data GlobalTag 和 2022 Golden JSON 仓库里完全没有，是全新要加的东西。

> 之前几次 session 的详细记录（geometric matching 研究、lostTracks recovery 研究、MVA v1/v2/v3 训练配比、2026-07-22 的 gen 多味化+jet/b-tag批次）已经归档进 memory（见 `[[project_geometric_matching]]`、`[[project_lostrack_recovery]]`、`[[project_phipimu_mva]]`、`[[project_h_phipimu_multiflavor_gen]]`），不再重复贴在这里。
