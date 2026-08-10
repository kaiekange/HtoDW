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

## Last Session (2026-08-10)

### MVA 背景从 W+jets MC 换成 SS data，并把特征集对齐后重训

**起因是一个无效的比较。** 上一轮把训练背景换成 SS data 训出新模型，跟旧的 W+jets 模型对比得出"相同信号效率下背景低 36–54%"。这个数字**当时不成立**：旧模型 40 特征、新模型 43 特征（多了 `jet_ht`、`n_b_jets_loose_dR04`、`n_b_jets_medium_dR04`），两者差的不只是背景。

根因是 `FEATURE_NAMES` 在 `DNN.py`/`predict.py`/`BDT/*.py`/macros 里各有一份拷贝，只靠注释 "must stay identical" 同步。

**修法：`analysis/DNN/common.py` 作为单一事实来源**，装 `FEATURE_NAMES`(43)、`LOG_FEATURES`(16)、`SENTINEL_FIX`、`SELECTION`/`B_VETO`、`load()`、`preprocess()`、`phipimuNet`、`TRAININGS` 表。`DNN.py`/`predict.py` 全部 import，自己不再声明。**`BDT/` 故意不接**（钉在更老的 42 特征集、已声明不可比，接进来会导致下次加特征时静默改变 BDT）。

训练改成参数化，两个配置只差背景：
```
python3 condor_submit.py ssdata    # signal(OS) vs Data_2022_SS，偶数 event_id
python3 condor_submit.py wjets     # signal(OS) vs WtoLNu2Jets 0J+1J+2J
```
产物在 `dsw_dsstw_{ssdata,wjets}_f43`（`_f43` 是特征集代际标签）。**wjets 刻意用 `parity=None`**——parity 切分只为保护 SS 的非循环 sculpting 检验，W+jets 没这需求，加了会让重训相对旧版同时改两件事。

**防漂移守卫**（特征顺序跟 `scaler_*.npy` 和 `nn.Linear` 是位置绑定的，改顺序不重训不会抛异常、只会静默出垃圾分数）：训练时写 `out_dir/feature_names.txt`，`predict.py` 打分前调 `check_model_dir()` 校验 stamp/scaler尺寸/输入层宽度。三条分支都实测过：旧40特征目录（无stamp）、顺序调换、长度不符，都给出可读拒绝。

**回归验证**：重构后重跑 ssdata，跟重构前**逐项一致**——scaler 逐位相同、`X_sig`/`X_bkg` 数组相同、early stopping 同在 epoch 98、test AUC 同为 0.9977。

### 有效的对比结果（OS sideband data，两个模型都没见过）

| ε_sig | wjets 模型 | ssdata 模型 | ratio |
|---|---|---|---|
| 90% | 7,789 | 4,990 | 0.641 ± 0.007 |
| 70% | 1,917 | 1,104 | 0.576 ± 0.015 |
| 50% | 751 | 388 | 0.517 ± 0.023 |

**结论在特征对齐后依然成立**：SS 背景让背景低 36–48%（原来声称 36–54%，方向对、幅度略小）。误差用"两模型选出的事例集高度重叠"这一点正确传播（拆成 A-only/B-only/共同 C 三块独立 Poisson），不是 sqrt(1/N1+1/N2)。**事例身份必须用 `(run, lumi, event)`——`event_id` 只在单个 run 内唯一，Data_2022 有 646 个跨 run 碰撞。**

**sculpting 也更轻**（SR/sideband 真-Ds 富集，越低越好）：ε_sig=50% 时 SS 样本上 0.455→0.318、W+jets 样本上 0.497→0.405。即 ssdata 模型**同时**在拒绝率和 sculpting 两条轴上更好。

**重训 wjets 的 AUC 是 0.9951，跟旧的 40 特征版一模一样**——新增那 3 个特征对 W+jets 完全无效，复现了当初删掉它们的理由。注意 AUC 是各自对自己背景类测的，**不是跨模型指标**，判据只能是上面 OS sideband 的拒绝率。

### 顺带修掉的 bug

- v3 的图全盖 "Simulation" 但背景是真实 data → CMS 标注和图例跟 `cfg["is_data"]`/`cfg["label"]` 走
- 光度权重守卫 `[check] weighted <num_jets_pt30>` 在 SS 改版时被删（唯一能抓 `WJETS_WEIGHTS` 静默失效的检查）→ 恢复，实测 1.080 vs raw-concat 1.55
- 旧的 SS 对比图没做 parity 切分，"trained on" 那条曲线含训练用的那一半（自证自洁）→ 训练后所有 SS 图强制 `event_id % 2 == 1`

### 脚本与图表清理

删了 14 个脚本：2 个指向已不存在的路径（`DNN/plot.py`→`MVA/Dsmu/`、`DNN/check_score.py`→`MVA/phipimu/v1/`，且各带一份 `FEATURE_NAMES` 拷贝）、3 个读已不再产出的 inclusive `WtoLNu2Jets.root`/去重前信号路径、4 个各画一张 Ds_mass 叠图、5 个各自重复实现同一张 score-cut-vs-质量谱图。

`analysis/macros/` 现在只剩：`compare_models.py`（取代那5个）、`draw_features.py`（43特征 signal/W+jets/OS data/SS data 四线对比）、`check_os_vs_ss_features.py`（OS-vs-SS 电荷位形检验）、`draw_nbjets_veto.py` 和 `screen_new_jet_vars.py`（b-veto 工作点和 43 特征选型的决策记录），加两个生产工具。图 19 张→5 张。

EOS 上删了 `dsw_wjets_test`(2.9G) 和 `dsw_dsstw_wjets/Prediction`(2.8G)；`dsw_dsstw_wjets` 的模型本体保留作 40 特征存档。

### 下一步没做的事

- **真 charm 富集是个待处理的系统性**：两个模型、两个背景样本上都随 cut 变紧长出真 Ds 峰（SR/sideband 从 0.14 涨到 0.32–0.50），sideband 外推会系统性低估峰下本底。要么在本底模型里显式加一个真 charm 成分，要么削弱分类器对真 charm 的偏好。
- `BDT/` 仍在旧的 42 特征集上，跟 DNN 不可比。

> 之前几次 session 的详细记录见 memory（`[[project_mva_test_build_v1]]`、`[[project_h_phipimu_ss_tuplizer]]`、`[[project_mva_sculpting_finding]]` 等）。2026-07-23 那次的 SS tuplizer / Tuples 改名 / 新信号样本上线记录见 `[[project_h_phipimu_ss_tuplizer]]` 和 git history。
