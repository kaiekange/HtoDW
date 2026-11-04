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

## Last Session (2026-11-04)

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

## 上一 Session (2026-11-03)

### 仓库清理：15_0_5 只剩 `analysis/`（旧 tuplizer 全删），13_0_13 不再依赖 15_0_5

**13_0_13**（`7a29ec3`、`58b6a31`）：删 `test/nano_base_preEE.py`（cmsDriver 骨架，无人引用）；
W+jets 文件列表 `filelist_{0,1,2}J_{2022,2022EE}.txt` 移到 `condor/das/`（`condor_submit.py` 的 `WJ_FILELIST_DIR`），
`data_2022.txt` 移到 `crab/DatasetList/`（`submitData.sh`），逐字节相同、仍被 `*.txt` gitignore；2023/2024 golden JSON 拷进 `crab/JSON/`。
`test/SIG_H_phipimu.py`、`WtoLNu2Jets_H_phipimu.py`、`plugins/H_phipimu*.cc` 是死代码（import 不存在的 JMEcorr），**还没删**。
未跟踪的 `plugins/_fillcand_body.txt`（09-10 草稿）、根目录 `BKG_datasets.txt` 已删。删 `H_phipimu*.cc` 要 el8 `scram b` 重编译，等 condor job 排空再做；`H_phipimu*_Tree` 类不能删（producer 拿它当 scratch struct）。

**15_0_5**：删 `Tuples/`（含未进 git 的 condor_logs）、`GenCount/`（被 Corrections/gencount 取代）、`BKGsets/`（未进 git）、
`EDAnalyzers/`（旧 `H_phipimu` + JMEcorr；删后 el9 容器里 cmsenv + correctionlib 2.6.4 正常）。
`analysis/`：删 corrections pilot 旧副本 10 个（都在 Corrections 仓库，`make_samples_nano.py` 本来就 import 那边的）、
只能配合已删模型（osdata/oswjets、wjetsres f24/nostage1、ssdata f55）或旧列名（`Lumi_reweight`/`event_id`）的脚本、
已结束的一次性研究（09-25 系列、tt̄ CR、OS−SS、BR/DY τ 分数、FD/顶点图、MET raw vs corr 等）、`functions/`（ROOT C++ 风格，无人用）、
`DNN/score_only.sh`、DNN/BDT 全部 scheduler log、`macros/condor_logs/` 里旧 apply_corrections/btag log。
**`common.TRAININGS` 删了 wjets/oswjets/osdata/osdata_phi20/wjetsres_nostage1/wjetsres_s1seed***；
`ssdata`、`ssdata_phi20`、`wjetsres` 模型已删但**保留**，因为 `ssdata_f45`、`resdy`、`wjetsres_v2` 由它们派生。
**所有默认 key 改成 `ssdata_f45`**（原 `ssdata` / `ssdata_phi20`，`phi_win` 由 f45 继承、行为不变）：
`draw_stage1_mass_by_score_cut`、`draw_cr_zmumu`、`draw_zmumu_feature`、`ds_mass_fit_io`、`run_ds_mass_fit.sh`、`fit_dy_jetbin_zmumu`、
`fom_scan_stage1`、`fit_sideband_scales_stage1`、`draw_stage1_dy_mass_by_cut`、`draw_tau3pi_masses_by_cut`、`draw_stage1_dsmass_by_era`。
验收：全部 `.py` py_compile、`.sh` bash -n、ds_dnn 下 `import common` + 三个现存模型目录存在 + `stage1_threshold(ssdata_f45, 0.90)=0.9596`；**脚本没有逐个实跑**。
留着的进行中工作：`check_{data,mc}_phi20nano_vs_old`（MC 版未提交）、`jme_zmumu_study`/`draw_zmumu_jme`、`fit_zmumu_r12_*`、`draw_zmumu_feature`、
`check_dsmunu_mass`+`compare_dsmunu_signal_gentag`（m(Dsμν) 进不进 MVA 未定）、`compare_stage2_dy_vs_nody`。用户未提交的 `draw_wp_mass.py --ds-only` 没带进 commit。

**下午追加（用户批准「φ20 之前用到的全删」）**：EOS `Samples/*_merged_prephi20*`、`skim_prephi20/`、`MC/2022/btag_efficiency_prephi20/`（~34 GB）删；
EOS postEE 11 个 `*_pre_phi20nano` 本底目录（63 GB，`eos rm -r`，在 EOS 回收站）删 → **新旧 nano 对比再也做不了**，
postEE WJ0/WJ1/DY2J/DY1J SS 新 nano 补齐前没有完整 postEE nano（现有 `*_merged.root` 仍含旧数据）。
13_0_13 `condor_logs/` 的 40 个 `*_pre_phi20nano.log`、`filelists_before_phi20nano/` 删；15_0_5 `macros/condor_logs/` 只剩 `make_samples_nano/lists/`（merge 要读）。
删 `check_data_phi20nano_vs_old.py`、`run_check_data_phi20nano.sh`、`fit_zmumu_r12_quick.py`（留 `fit_zmumu_r12_cutscan.py`）。
`check_mc_phi20nano_vs_old.py` + `run_check_mc_phi20nano.sh` 被我一并删了（未跟踪），用户要提交，从 AFS 备份卷 `work.kakang.backup`（11-02 18:23 快照）恢复后提交；
它们读的 `_pre_phi20nano` 已删，**现在跑不了**，只作记录。`draw_wp_mass.py --ds-only`（只画 Ds_mass panel）一并提交。

### 我犯的错
1. 每删一组就 commit + push（一天 7 个小 commit），用户恼火 → memory `feedback_batch_commits`。
2. 删了用户未提交的两个脚本（虽在批准的选项里），用户随后要提交 → 从 AFS 备份卷恢复。
3. 又用 shell 变量拼 `rm` 路径被安全检查拦；删除一律写字面路径。
