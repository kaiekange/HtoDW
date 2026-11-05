# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working Convention

At the end of each session, update the "Last Session" section below. Keep only the two most
recent sessions here; move older ones verbatim to the top of `SESSION_LOG.md` (not auto-loaded).

## Environment Setup / Repository Layout

As of 2026-11-03 this checkout holds **only Python analysis code** (`analysis/`); there is no
C++ package to build here any more (the old private tuplizer `EDAnalyzers/`, `Tuples/`,
`GenCount/` were deleted, recover from git history if needed).

- **NanoAOD production** (C++ producers, cmsRun psets, condor/CRAB submission) lives in
  `/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/EDAnalyzers/Reconstruction/`
  (el8 release): `plugins/HDsPhiPiMu{,GenTag}Producer.cc` + `BestPhiPiMuCandFilter.cc`,
  `python/nano_phipimu_cff.py`, `test/NANO_{,DATA_}H_phipimu.py`, `condor/condor_submit.py`
  (MC), `crab/submitData.sh` (data).
- **Corrections** (friend trees, b-tag efficiency maps, normalization, signal gencount) live in
  `/afs/cern.ch/work/k/kakang/H2DW/Corrections/`. It uses this CMSSW_15_0_5 release only as the
  el9 Python stack (correctionlib/uproot), inside `cmssw-el9`:

```bash
cmssw-el9   # on an el8 node
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/ && cmsenv
```

- **`analysis/`** here: `macros/make_samples_nano.py` (merged skims -> `NanoMVA/Samples/`),
  `DNN/` + `BDT/` (training/scoring, single source of truth `DNN/common.py`; run under conda
  `ds_dnn`, submit with `condor_submit.py`), and the plotting/fit macros in `analysis/macros/`
  (matplotlib + mplhep, figures via `common.fig_path`).

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

### Candidate Producer Data Flow (13_0_13 `HDsPhiPiMuProducer`, copied from the old `H_phipimu`)

Differences from the old tuplizer described below: jets come from the standard-NANO
`updatedJetsPuppiWithUserData` (no JMEcorr chain, no pT threshold in the producer loop); the
muon pre-cut is 24 GeV on the uncorrected pT (the analysis cut, corrected pT>26, is applied
offline); the post-fit φ cut is `|M(KK)−Mφ|<20 MeV` (since 2026-10-26); `BestPhiPiMuCandFilter`
keeps only the highest-Ds_pt candidate, so the nano has one row per event. GenTag variant:
`HDsPhiPiMuGenTagProducer` (background truth tags `cand_gen_*`).


1. **Gen filter** (`doTruth` only): require exactly 1H + 1Ds(*) + 1W + 1 hard-process W-daughter lepton (e/μ/τ, any flavor — `Gen_W_lep_pdgId`), then find the status-1 gen K+, K−, π with the full φ ancestry chain (`TwoBodyDecay(gp, 333, 321, -321)` etc.). Returns early if not found. DsstW additionally requires the Ds* to go via γ **or** π0(→γγ) (`Gen_Dsst_mode`: 0=γ, 1=π0), both ancestry-checked the same way. If the W-lepton is a τ, a second pass over `isDirectHardProcessTauDecayProductFinalState()` particles classifies `Gen_tau_decay_mode` (0=hadronic, 1=e, 2=μ) and — only for τ→μ — supplies the muon used below. **The tree stores no gen-level kinematics** (no Gen_H/Ds/Dsst/W/Kp/Km/pi/mu/gamma/pi0/tau 4-vectors, no PV/SV, no Gen MET, no gen dR family) — only the categorical tags (`Gen_W_lep_pdgId`, `Gen_tau_decay_mode`, `Gen_Dsst_mode`) and the match booleans below. This was a deliberate 2026-07-22 simplification; the local `reco::GenParticle` objects needed for matching still exist in `analyze()`, just aren't written out.
2. **Gen matching**: match PF candidates via `packedPFCandidateToGenAssociation` + ΔR<0.05 + Δpt/pt<0.3 (ancestry checks accept either Ds* submode). The γ (DsstW, γ-mode) is matched geometrically instead — ΔR<0.02, Δpt/pt<0.5 — because it has no association; in π0-mode the two π0 daughter photons are matched the same way, producing `Gen_pi0gamma1_match`/`Gen_pi0gamma2_match`. Muon via `muon.genParticle()`; for a τ→μ event the ancestry check uses `isDirectHardProcessTauDecayProductFinalState()` instead of the direct W→μν chain.
3. **Trigger**: require `HLT_IsoMu24_v` or `HLT_Mu50_v`. Returns early if not fired.
4. **Primary vertex**: take `offlineSlimmedPrimaryVertices.front()`.
5. **Best muon**: highest-pT tight muon with pT>26 GeV (raised from 24 on 2026-07-22 — above the `HLT_IsoMu24` turn-on so the offline selection sits on the trigger efficiency plateau).
6. **Ds candidate loop**: triple nested loop over the charged-hadron pool (`hasTrackDetails`, `trackHighPurity`, pT>1, |pdgId|=211). Per-track: K± pT>2, K+ and π opposite in charge to the muon, K− same as the muon. Pre-fit: (KKπ) pT>5, |M−MDs|<0.5; **φ pre-fit** M(KK)<Mφ+0.08 (loose, pre-vertex-fit speed cut on raw KK mass). Kinematic vertex fit via `KinematicParticleVertexFitter`. Post-fit: vertex prob>0.001, |M−MDs|<0.3, pT>7, dira>0.99, (Ds+μ) mass<125, **φ: pT(KK)>3 and M(KK)<Mφ+0.05** (upper-cut only — the KK threshold, Mφ−32 MeV, physically bounds the lower side, so no lower cut is applied). Ds/φ windows widened 2026-07-22 (post-fit Ds ±0.15→±0.3, pre-fit KKπ ±0.3→±0.5) to give the Ds-mass sideband more lever arm for an SR/CR split; the pre-fit window is intentionally wider than post-fit so it never truncates the sideband.
7. **Ds* (DsstW only)**: photon pool is |pdgId|=22, pT>0.5; (Ds*+μ) mass<125. Mass window unchanged (±0.15) — not yet revisited pending a DsstW rerun to look at its resolution.
8. **Output**: the tree is row-per-candidate (flat ntuple), not one row per event — `tree->Fill()` is called once per Ds×gamma combination (or once per Ds if no gamma passes), so an event with N candidates writes N rows. Event-level fields (gen info, PV, MET, `jet_ht`/`jet_leadpt`/`jet_maxbscore`, triggers, muon, electron/muon os/ss lepton counts) are set once before the candidate loop and repeat identically across that event's rows. `Ds_nearestjet_dR`/`bscore` and `n_b_jets_{loose,medium,tight}{,_pt20}_dR{03,04}` (b-tagged central jets excluding jets within dR03/dR04 of the Ds, added 2026-07-23) are genuinely per-candidate: recomputed inside the loop from each row's own `Ds_mtm`. `nCandidates` = number of rows passing all cuts.

### Key Constants
- GlobalTag (2022 MC): `130X_mcRun3_2022_realistic_v5`
- Ds mass: 1.96835 GeV, φ mass: 1.019460 GeV, K: 0.493677 GeV, π: 0.13957039 GeV, μ: 0.1056583755 GeV
- Input MiniAOD on XRootD: `root://maite.iihe.ac.be//store/user//kakang/Analysis/ggHtoDsW_phipi/` (DsstW: `ggHtoDsstW_phipi/`; `output_188` there is bad and is skipped in the config)
- Nano outputs: `/eos/home-k/kakang/H2DW/Analysis/NanoMVA/{MC,Data}/`; merged skims: `NanoMVA/Samples/`

## Last Session (2026-11-04 晚 ~ 11-05)

### JME 完整链进 friend 并成为分析默认；全样本 friend / skim / merge 重做

**方案（用户定）**：JES（data residual 按 run）+ MC JER smear，Type-1 PuppiMET 传 JER + MET xy；μ 修正不传 MET；nano 原值保留。
系统误差（MC）26 个：Reduced JES 11 源（`_year` 源跨 era 不相关）×2、JER×2、unclustered×2（只动 MET）；JES shift 在 JER smear 前（tutorial 顺序）。
- Corrections：`jme_direct.Jme.variations()`、`jme.full_columns()`（`jme_*` 列，闭合到 10-22 验证过的 `jet_pt`/`met`+`met_xy`）、`btv` 接受 jet pt、
  效率图 `--jme`（JME-pt，`btag_efficiency_jme/` 36 组，旧图保留）。friend 222→261 列，约 2 倍大。
- `make_samples_nano.py`：7 个 jet/MET `*_corr` = JME，nano 版改名 `*_nano`；cut 用 `jme_jvm_veto`/`jme_n_bjet`，权重 `jme_w_btag_nom`；
  只在某个 JME 变化下过 b-veto 的 MC 行进 `EventsSyst`（`weight`=NaN，`weight_nobtag` × `jme_w_btag_<v>` = 变化权重）；ttcr 拒绝运行；merge 拒收旧 chunk。
  jet veto map 只有 nominal。
- 结果：data+MC friend、SR+zmumu skim、18 个 merge 全部重做并验收；事例数比旧版多 0.1–0.9%（JME jet 的 jet veto map veto 更少，data 3.68%→3.31%）。

**Z→μμ（JME 样本）**：DY jet-bin SF 0.807/0.763/1.168（≈ nano 版；`common.DY_JETBIN_SF` 未改，仍 0.802/0.765/1.155）。
**MET data/MC 变差**：nano 1.07→0.80，JME 1.23→0.43（MC 均值 +2.4 GeV）；N_jet≥4、m_T 高端也变差。= 10-22 已知的「MET 传 JER 变硬」。
**待定**：MET 传不传 JER（选项：保持 / `propagate_jer=False` 全部重跑 / 先本地拆原因——我推荐）；DY SF 是否更新；JES/JER 对 N_jet 的影响未拆。
旧 DNN/BDT 模型与 Prediction 不匹配新样本（MVA 暂停）；`draw_zmumu_jme.py`/`jme_zmumu_study.py` 的 nano-vs-JME 对比因改名失效。

**其他**：删孤儿 friend 105 个（MC 104 + data Run2022G `nano_52`）。apply_corrections 的 per-job log 目录满（AFS 单目录 ~2.2 万文件，"File too large"）
→ 2022EE W+jets 1J_SS 583 次失败、用户让 `condor_rm` 全部；`run_apply_corrections.sh` 改为只留失败 job 的 log，旧 log 全删。hadd/merge 改回由我做（用户撤销旧规定）。

### 我犯的错
1. 定 JME 方案时没翻出 10-22 的「MET 传 JER 变硬」提醒用户，全样本跑完才在 Z→μμ 看到。
2. 设计变化列时没追到 skim 的 b-veto（并集问题），advisor 指出、用户问「为什么之前没想过」。
3. 没查 log 目录文件数就交几千 job；还把未扫描的 per-job log 一起删了，只好另交 condor 补查。
4. 例行用 `eos ls` 核对、回复用了英文，用户不满。

## 上一 Session (2026-11-04)

### POG 修正逐类对照官方例子 + Corrections 按 POG 拆模块（MUO 用 muonscarekit .cc、BTV 改 C++）

**对照**（用户把 BTV/EGM/JME/LUM/MUO examples + muonscarekit 放在 ~/work/tmp，之后删了；地址和 commit 记在 Corrections `external/POG_EXAMPLES.txt`）：
- LUM 和例子一致。
- MUO：随机数用 JSON `RandomSmearing`（muonscarekit 新写法，旧 MUO examples 的 SeedSequence 只用 event 号）；muonscarekit 的 .py 和 .cc 不一致的两处按 .cc（用户定）。
- BTV：公式和例子一致；heavy/light 一起变化，没有 correlated/uncorrelated 拆分，**未改**。
- **EGM 完全没做**（ID/Reco SF、scale/smearing 都没有，只影响 eμ/ttcr）。
- JME：不做 JEC/JER（10-21 定）；jet ID 用 nano `Jet_jetId` 位（NanoV12 在 |η|>2.7 偏松，POG 用 `jetid.json`），**未改**；MET xy 没做。
- muon η 截到 ±2.399 没有依据（JSON 是 flow=clamp），**已删**（3 处）。

**拆分**（Corrections 仓库；列名和内容不变，下游不用改）：`lum.py` / `muo.py` / `btv.py` / `jme.py`，`apply_corrections.py` 只负责读 nano、调用、闭合、写 friend。
- `muo.py` 用 ROOT 编译 `external/muonscarekit/MuonScaRe.cc`（@1c7426b，原样）。`muo_scalesmearing.py` 已 `git rm`（连同 11-03 未提交的 RNG 改动，被 .cc 取代）。
- BTV 只有 Python 版权重（`lib/btag_sf.py`，要 Python ≥3.10，CMSSW_15_0_5 是 3.9），用户定改 C++：`corrections/btag_weight.cc` 逐行移植 `compute_btag_event_weight_multi_wp`。
  `checks/validate_btag_weight_cc.py` 和官方逐事例比（单/多 WP、各种 variation），差 ≤1.3e-15。`btv.py` 发现 jet 落在 ε=0（tagged）或 ε=1（untagged）的 bin 就 exit（效率图过时）。
- `make_samples_nano.py` 的额外 μ 改调 `muo.py`。
- 验收：signal/tt̄/data（2022）各 1 个文件，拆分前后 friend 222 列 ≤1e-12；`n_xmu≥1 ≡ ~pass_muveto` 成立。**2022EE 和 condor worker 上没跑过**（worker 要 JIT 两个 .cc），正式提交前先交 1 个 job。
- 代价：每个 job 多约 30 s（import ROOT）；MC 中 event>2³¹ 的事例 smearing 随机数和旧移植不同（.cc 用 32 位 int，同 RDF）。

**待定**：EGM、jet ID、BTV 系统误差拆分、MET xy 要不要做；2022EE nano 补齐后先重做效率图；φ20 新 nano skim 要 `--correct`，Z→μμ 的 DY SF / r12 重做。

### 我犯的错
1. 用户问「是按这个来的吗」我反问「这个是什么」，没自己去找参考对比（memory `feedback_compare_to_reference_code`）。
2. 只要代码对比时我去跑数值验证，被打断。
3. `btag_weight.cc` 文件头的 BTV commit 号先凭记忆写错，核对后改正。
