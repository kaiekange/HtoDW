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

- **`analysis/`** here: `macros/make_samples_nano.py` (merged skims -> `NanoMVA/Samples/<year>/`, `--year=2022|2023`),
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
- Nano outputs: `/eos/home-k/kakang/H2DW/Analysis/NanoMVA/{MC,Data}/`; merged skims: `NanoMVA/Samples/<year>/` (since 2026-11-12: `2022/` = 2022+2022EE, all groups; `2023/` = 2023+2023BPix, data only; group/file names identical, `common.SAMPLES_DIR` = `Samples/2022`)

## Last Session (2026-11-13)

### 2024 data CRAB 查状态 + 一次 resubmit（无代码改动）

- 32 个任务（`CMSSW_15_0_2/.../crab/crab_Muon{0,1}_Run2024{C..I}[_v2]_NanoH_phipimu[_SS]`，用户 1112 交）。我的 shell 里 `crab status` 能用了
  （cmsenv 15_0_2 + `crab.sh prod` + `X509_USER_PROXY`）。
- 上午：71 个 failed（8028 读输入为主；60328 ASO→EOS HTTP 500；50660 内存 3500>3000 MB ×5；−1 ×1）→ 用户让重交，22 个任务 `crab resubmit`
  （Muon0 D OS、Muon0 H SS、Muon1 G SS 加 `--maxmemory 4000`），复查全部回到 idle/running；Muon0 D SS 新出 2 个 8901 也重交。
- 15:44：12 个任务 100%；又有 26 个 failed（10 个任务；Muon0 D SS、Muon0 F OS 是 50660，重交要加 `--maxmemory 4000`），**用户：先不交**。
  2533 个 idle 集中在 8 个任务；每个 dataset disk 副本只有 2-3 个站点（Muon0 I = Vanderbilt+Bari，完成的全在 Bari，OS/SS 几小时 0 running）；
  CRAB 警告 3000 MB 只保证 2500。原因未定，未处理。proxy 当时只剩约 44 min，需用户续。
- 2024 生产代码（`CMSSW_15_0_2`）仍不在 git 里。详见 memory `project_data_datasets_2022_2026`。

## 上一 Session (2026-11-12)

### recoil calibration 进整套修正（CDF 映射，在 skim 里应用）；friend 删 nano 列；friend + SR skim 全量重跑

**用户定**：nominal = CDF 映射（JME POG 原文："map the simulated and measured cumulative distributions, which also
accounts for non-Gaussian shapes"）；按 era 分开推导；q_T 不卡；参考 boson：DY = Z、W+jets = W、**信号 = H**
（Ds 当衰变产物扣掉）；tt̄、data 不修；Z→μμ 的 map 直接用到其他 MC。recoil 在 skim（make_samples_nano worker）里应用。
**不确定度和各 JME 变化的 map_v 推后**（POG 要求在 JEC/uncl 变化下重做推导，不要重复计算）。W+jets sideband 外推检查撤回。
DY jet-bin SF 推导时照乘现值。nano jet/MET 列删掉（friend 也重跑）。

- `recoil_zmumu.py`：data/TT 用 reco μμ，DY 用 gen V（V = hard-process 轻子和，信号用 H；V_vis = V − hard-process ν）；
  q_T[0,10,20,30,50,80,∞) × N_jet{≤1,2,≥3}，x_∥ = u_∥ + q_T 和 u_⊥ 各自 CDF 映射到 data − TT；map `Samples/CR/zmumu_recoil/recoil_<era>.npz`。
  闭合（偶推奇验）：MET 60–150 GeV data/MC 0.51/0.47 → 0.86/0.91（2022/2022EE），m_T 1.00/0.99；分格图 `cr_zmumu_recoil_cells_*`。
  缩放法（中位数+宽度）尾巴只到 0.73/0.87，弃用。
- `make_samples_nano.py`：`gen_boson_cols`、`recoil_correct`；DY/W+jets/信号（OS、SS）的 `met_*_corr`、`W_tmass_corr` = recoil 后；
  新列 `*_prerecoil`、`gen_V*`；26 个变化 MET 改名 `jme_met_*_<v>_prerecoil`（**未做 recoil，和 nominal 不可比**）；`*_nano` 删；meta 记 map sha1。
- Corrections：`apply_corrections.py`/`jme.py` 删 13 个 nano 列，μ 变化 `W_tmass_*` 改用 JME MET。`check_dsmunu_mass.py` 改读 jme_ 列。
- 重跑：SR + friend 506 job 全 0；SR merge 9 个完成（W+jets ⟨MET⟩ 46.72→44.43）。
- 晚上：zmumu skim 506 job 全 0 → merge 9 个完成 → `*_prerecoil` 与旧 zmumu MET/m_T 逐行一致（差 ≤1.6e-9，浮点舍入），
  DY gen_V* 与 genv 逐位一致 → `recoil_zmumu.py` 复推两个 era 的 map **与 skim 所用的逐位相同**，DY met_pt_corr 与 map 差 ≤5e-4 GeV。
  已删：`zmumu_jme_steps/`、nojer/`_jer` npz、`genv/`、`old_zmumu_met/`、旧 `MC/2022/btag_efficiency/`、`CR/zmumu_jme/`、
  `jetpt_corr_ratio_eta.py`、`draw_signal_met_jme.py`；JME 拆解两个 .py 最终版 commit `5cf2caa` 后 `git rm`（连同 run_*.sh，`e8b7075`）。
  所有已验收的 condor log 删掉（15_0_5、Corrections、13_0_13；只留 `make_samples_nano/lists/` 和 13_0_13 的坏 JOBID 记录）。
  详见 memory `project_recoil_into_chain_20261112`。

### 晚上：2023 data 下游（不需要 MC）→ `Samples/2023/` SR + zmumu
- 2023 data 生产：24 任务 23 个齐；Muon1 C v4 OS job 7（FileReadError×6）放弃，SS 保留（我误删后从 EOS recycle 恢复）；golden = 2026-02 版。
- **用户定**：按年分目录 `Samples/<年>/`、组名不改；region SR + zmumu；lumi 2023C 17.981、2023D 9.516 fb⁻¹（显示 2 位），eras.py 加 3/4。
- Corrections（`b4e497f`）：2023 b-veto 用 RobustParT（`bem.ERAS[*]["tagger"]`、`btv.TAGGERS`，WP_T 0.7969/0.7994）；
  JME 2023 payload（2026-06-05，`Summer23Prompt23`/`Summer23BPixPrompt23`），MET xy 文件名按 era；**2023BPix L2Relative 带 JetPhi**
  （`jme_direct.jes` 加 phi）。2023 MC 的 JER（按 run 分 tag）、b-tag 权重/效率图未接，会报错。回归：同机新旧代码 2022 data/MC friend 逐位相同。
- skim（`e8b7075`）：`ERA_YEAR`、`GROUPS_BY_YEAR`、`DATA_TASKS`/`data_pd`（Muon0/Muon1、`_vN`）、`--year`；fom_scan_stage1/fom_map_stage2_run3
  的总光度改为只取 era 1+2。EOS `Samples/*` 全部 mv 进 `Samples/2022/`，recoil_zmumu 复推 map 仍逐位相同。
- condor：SR `--correct` 144 job + zmumu 144 job 全 0；merged：SR OS 1232998 / SS 1275562，zmumu OS 98220 / SS 97839；
  era 只有 3/4，OS SR 窗口 0 行，列与 2022 相同。**每 fb⁻¹ 产额比 2022 高约 30%（SR 和 Z 窗都是）**，原因未查（b-veto PNet→ParT、pileup/组合候选？）；已排除重复事例（merged 里 (run,lumi,event) 全唯一）。
- 1112 晚三个分支都已 push（combined_channel `e8b7075`、corrections `b4e497f`、CMSSW_13_0_13 `c531f25` = 2023 data 生产代码）。
  未跟踪的 `check_mu_tkiso_eff.py`/`run_check_mu_tkiso_eff.sh` 用户说先留着。
- 下一步：2024 data 交 CRAB（用户；证书约 12-01 到期）；2023 MC 生产；之后 MVA。详见 memory `project_2023_data_downstream_port`。

### 我犯的错
1. memory 转述 POG 时漏了 "accounts for non-Gaussian shapes"，先做了缩放法；被问"原文怎么写"才去会话记录翻原文。
2. 让用户选方法时只给数字不解释（"我为什么要确认这个方案，你什么都没解释"）。
3. 提出 W+jets sideband 外推检查，被问才想到 sideband MC 缺 QCD、得不出结论，撤回。
4. 154 个 genv job 没先单 job 上 condor 测就整批交（advisor 指出）；`pkill -f` 的模式匹配到自己的 shell。
5. 用户说"失败的不管，删了就行"，我套用 MC 的 OS/SS 两边删规则删了 2023 data SS 的好文件（已从 recycle 恢复）。
