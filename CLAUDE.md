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

## Last Session (2026-08-14)

### seed 噪声底 + 工作点定在 ε_sig=70% + BDT 全对齐接入 common.py

**先量了一件之前从没量过的东西：同一配置换 seed 重训，结果会晃多少。** seed 42（现有模型）
+ 101/102/103 三个新作业，配置逐字相同。OS sideband N_kept 的 RMS：**90% 处 1.09%、70% 处 1.48%、
50% 处 1.85%**。四个模型打分的是同一批 686,497 个事例，所以这里面**没有独立 Poisson 成分，
100% 是模型自身涨落**（别拿 √N 的 1.4/3.0/5.1% 去解释它）。n=4 时 RMS 自身有 41% 不确定度，
当阈值用要往上取整：**90% 处 ≳2%、70% 处 ≳3%、50% 处 ≳4% 才算超出噪声。**
副产物：四次训练 early stop 在 98/105/118/151，val AUC 全在 0.9975-0.9977 —— **AUC 平台上停在哪
基本是随机的**。

**「本底太多要不要先手工筛一遍」→ 不用，BCE 已经在做。** 背景事例在 BCE 下的梯度**恰好等于它
自己的 score**，所以 score 分布就是学习信号本身。实测：score<0.01 的 83.3% 只贡献 **3.2%** 梯度，
最难的 3% 贡献 68.5%。想筛掉的那些一共只值 3%，而手切的 cut 打分时还得照施。**真正的瓶颈是难区
事例数不够**（训练集当量 90% 处 ~808、70% 处 ~94、50% 处 ~14），而预筛只会减少行数。

**工作点：ε_sig = 70%。** 决定性理由是 **sideband 是拟合本底的地方**——70% 剩 1,108 个事例，
50% 只剩 386，而灵敏度从 ~75% 往下在误差内是平的（FOM 归一化：90% 0.605±0.023 / 70% 0.945±0.090
/ 50% 1.000±0.153，**别读成"50% 最高"**，差 0.055 而合并误差 ~0.18）。稳定性上 **50%→40% 之间
有清楚断裂**（产额 RMS 1.85%→5.09%→8.56%），跟训练统计量对得上。脚本 `scratchpad/wp_scan2.py`。
注：不同 seed 选出的**事例集** Jaccard 只有 0.75(90%)/0.65(50%)，但那量的是边界有多"软"，
**不是不稳定性**——拟合只关心产额和质量形状，这两个到 50% 都稳。

**BDT 全面对齐并接入 `common.py`。** 之前 BDT 跟 DNN 三处同时不同（42 特征的 fork、W+jets 背景、
单信号），"BDT vs DNN"什么都测不出。现在共享 43 特征/`load()`/`preprocess()`/`SIGNALS`/`TRAININGS`
/类权重/70-10-20 split/seed，**只剩算法不同**。XGBoost 超参和 `eval_metric="auc"` 一个没动
（依据：上一版 1085 轮里 490 轮落在最佳 val AUC 的 1e-4 内，停在哪轮几乎不改变模型）。
`common.py` 当初拒绝接 BDT 的顾虑用 `check_feature_stamp` 守卫**中和**了，不是无视；四条分支
提交前全部实测通过。`predict.py` 里第二份手抄的 42 特征表和那句已失效的 "Must stay identical to
DNN/DNN.py" 一并删掉。基础设施同时加了 `out_dir_for(key, tool, seed)`（**`cfg["out_dir"]` 原样
保留**，7 个 macro 靠它）、`DNN.py`/`BDT.py` 的可选 seed、`predict.py --samples`、
`compare_models.py` 的 `[TOOL:]bkg[:seedN]` 泛化（回归验证：旧命令逐项复现旧结果）。

**结果：BDT 在未见的 OS sideband data 上少 11%/18%/34% 背景**（0.894±0.009 / 0.822±0.021 /
0.657±0.033），远超噪声底。**已排除 overfit 假象**：用从未参与训练的 test-only 信号定阈值，
比值 0.892/0.822/0.658，几乎不变。

### 判据被纠正：SR/sideband 富集不是问题，形状对不上才是

用户指出：data 和 W+jets 里确实有真 D/Ds，score cut 富集它们是物理上该发生的；**问题在于 cut
之后 data 和 W+jets 的形状对不上**。**以后别再拿 SR/sideband 富集当模型好坏的判据。**

改用低/高双比值（归一无关），配对 bootstrap：no cut 1.016±0.005 → 90% DNN 1.291±0.058 vs BDT
1.280±0.060（0.3σ）→ 70% 1.777±0.165 vs 1.644±0.177（1.0σ）→ 50% 2.568±0.448 vs 1.924±0.393
（1.6σ）。**测不出区别。** DNN 那列逐项复现 0811 记录。

**这是一条新的独立证据：tilt 与算法无关。** 0811 是换训练背景证明与分类器无关，现在 BDT 是
**换算法**（GBDT vs NN）却给出几乎相同的 tilt。**换 MVA 修不了它。**

**净结果：BDT 在拒绝率上赢 18%（ε_sig=70%），在真正要紧的形状闭合上打平。**

### 加权假说被否掉（问得好，答案干净）

问"全范围的 data/MC 差距会不会是 jet-bin 权重错了"。**不是，因为权重不随 score cut 改变。**
b-veto 区 OS sideband 的 data/W+jets 比值，0-jet 行 1.108 → 0.973/0.862/0.816，
≥4 jet 行 3.280 → 0.976/1.073/1.079。0J 权重若错 11%，0-jet 那行每一格都该是 1.108。
**能被分类器切掉的是物理成分，不是归一化因子。** cut 后残留 ~±10% 且 N_jet 上平 → 只剩"三个 bin
公共的整体归一化"这一种可能。低/高双比值更是比值的比值，权重精确抵消。
（≥4 我得 3.280 而 0811 记 4.761：我限制在 M_Ds±0.3 sideband，`check_jetbin_closure.py` 不加质量
窗口。窗口不同，不是矛盾。）

### 两个自己犯的错，记下来

1. **拿 b-tagged 区的 171/128 当论据是错的。** 分析已经 veto 掉那个区，171 不能量化 veto 之后的
   残留，而残留才是唯一有意义的量；且 b-veto 区超出只有 1.44，权重误差原则上能解释那个量级。
   **b-tagged 只能说明"存在一个 MC 没有的过程"，不能当定量论据。**
   **veto 后的 tt̄ 残留现在量不出来**——要 tt̄ MC（不在当前 `Samples/`），或一个本身也需要 tt̄ MC
   才能定的 transfer factor。
2. **"50% 时 W+jets 在 m≈2.05 有个峰" —— 没有峰。** 判据是看 70% 那档（统计量 3 倍）：完全光滑。
   50% 时的假象来自 2.018 bin 偏低 2σ 造出假谷，加上 2.058 里 3 个 0J 事例各带权重 2.33——它裸
   事例只有 38（比邻居的 54 还少），画出来高度却相当。**图上是加权高度，眼睛读不出背后几个事例。**

### 下一步

- **BDT 自己的 seed 噪声底没量**（1.09/1.48/1.85% 是拿 DNN 量的）。18% 要被 seed 涨落吃掉需要
  BDT 的 floor 比 DNN 大十倍，不太可能，但没实测。`BDT.py ssdata 101/102/103` 即可。
- **`check_data_vs_wjets.py` 等 6 个 macro 仍只认 DNN**（读 `cfg["out_dir"]`）。本 session 的
  data/W+jets 对比是另写脚本算的，应并回去做 `out_dir_for` + `TOOL:key` 泛化。
- **SR closure 仍未解决**，且现在多一条证据说明换分类器解决不了。真正的前提问题没变：
  最终本底是拟合 data sideband 还是用 MC 给形状。

## 上一 Session (2026-08-11)

### 追查"DNN cut 越紧 data 跟 W+jets 差别越大"，结论跟原来的假设相反

只看 ssdata 模型（W+jets MC 在这里不是训练背景，只是被比较的绝对预言）。三个新/改脚本：
`analysis/macros/check_jetbin_closure.py`（新）、`check_feature_tail.py`（新）、
`check_data_vs_wjets.py`（加了 SR/CR 双比值、N_eff、四质量窗口、D± excess）。

**先纠正两处读图。** `compare_models_mass.png` **不是 data/MC 闭合图**：绿色是 OS W+jets 全
质量范围、蓝点是 SS data 的 odd parity 一半、黑点是 OS data 且 blind 窗内为空，三种归一互不
可比（docstring 里"按 sideband 归一"是过期的假话，已改）。另外 **1.44 是 post-b-veto 的数**：
裸的 CR 比值 OS 2.19 / SS 2.29，加 1 tight muon 变 1.92/1.98，再加 tight b-veto 才是
1.44/1.46 —— b-veto 砍掉 25.3% 的 data 但只砍掉 0.7% 的 W+jets MC。

**关键发现：1.44 这个缺口在质量上平，在 jet 多重数上完全不平。**

| N_jet(pT>30) | OS b-veto | SS b-veto | OS b-tagged |
|---|---|---|---|
| 0 | 1.108±0.005 | 1.118±0.006 | 2.19±0.17 |
| 1 | 1.138±0.004 | 1.192±0.005 | 20.2±0.7 |
| 2 | 1.787±0.007 | 1.810±0.009 | 72.2±1.8 |
| 3 | 2.762±0.015 | 2.713±0.019 | 104.7±3.2 |
| ≥4 | **4.761±0.038** | 4.649±0.045 | **171.5±6.3** |

χ²/ndf = 26941/4。OS 与 SS 一致 → 电荷无关。b-tagged 区的 171 是 tt̄ 的直接证据（W+jets MC
在那里几乎什么都预言不出来）。21 万事例的超出有 80% 坐在 ≥2 jet。所以过去那句"比值平 ⇒ 形状对、
只缺归一"是错的推论 —— 它只在被看过的那一个变量上平。

**DNN 把这些超出几乎全部去掉了。** ε_sig=90%、OS b-veto CR：ε(W+jets MC)=1.105%，
ε(未建模超出 := data−MC)= **−0.132±0.046%**；70% 时 0.233% vs −0.003±0.021%。
缺失本底比 W+jets 本身好拒约 10 倍 —— 这正是"用含 tt/DY/QCD 的 SS data 当训练背景"应该买到的
东西，现在是量出来的而不是假设的。

**而且不是尾部放大效应。** ε_data/ε_MC 对 DNN score **在 ε_MC 从 1e-1 到 5e-4 三个数量级上
恒为 0.65-0.68**。完全去掉超出应给 1/1.44=0.694，实测 0.65-0.68 → **超出去除解释了这个因子的
约 95%，剩下约 5% 是真实的特征失配残差**（ε(未建模)=−0.132±0.046% 比 0 低 2.9σ，即分类器拒掉的
data 比"只拒超出"应有的还多一点）。所以之前"CR 回到 1.0 是两个大小相当的错误相消"那个读法**量级上
是错的**——主项是超出去除，不是相消；但 **CR 闭合到 ~5%，不是精确闭合**，引用时别把这 5% 丢掉。

**没有单个特征能复现它**，而且有两个特征本身失配、方向相反：ε_MC=1e-2 时 DNN 是 0.66，
`Ds_pt` 单独切是 2.52（到 5e-4 时 4.57），`Ds_ipsig` 单独切是 0.54（到 0.35）。
**坑**：`Ds_IsoR04_sumChargedHadronPt` / `Ds_IsoR03_RelIso` 有巨大的单值尖峰，阈值切不动，
ε_MC 卡在 9.7e-2（`jet_ht` 卡在 5.7e-3），它们那几行"平"是假象。**必须打印实际达到的 ε，不是目标 ε。**

### 真正没修好的两件事

**1. SR 特有超出在涨，不是常数。** SS SR/CR 双比值（归一无关）：1.01±0.01 → 1.15±0.04(99%)
→ 1.28±0.12(90%) → 1.57±0.36(70%) → 1.82±0.48(65%)。单调。**不给 p 值**：21 个扫描点是
嵌套 cut、共享事例，跨点做 χ² trend test 无效，要显著性得改用不相交的 score bin。
MC 的 N_eff 在 ε_sig≈50% 掉破 30，之后是噪声。

**2. 一个 OS 独有的低质量结构。** 四个质量窗口（都不碰 blind 区）：

| ε_sig | OS <1.830 | OS D± 1.83-1.91 | OS >2.008 | SS <1.830 |
|---|---|---|---|---|
| no cut | 1.450 | 1.439 | 1.427 | 1.468 |
| 90% | 1.055±0.033 | 1.006±0.043 | 0.817±0.025 | 0.928±0.059 |
| 70% | 1.355±0.091 | 0.997±0.083 | 0.763±0.051 | 0.750±0.139 |
| 50% | 1.787±0.217 | 1.039±0.142 | 0.696±0.078 | 0.782±0.291 |

**OS 和 SS 都在变差，而且方向相反。** 光度加权 W+jets 是绝对预言（σ×L，无自由归一），
所以每个窗口的比值本身就有意义，不需要挑基准：

- OS：高质量端 MC **多**预言且越切越多（0.817→0.763→0.696→0.391@20%）；低质量端 MC **少**
  预言且越切越少（1.055→1.355→1.787→2.100）
- SS：反过来 —— 低质量端 0.928→0.750→0.782，高质量端 1.060→1.233→2.006

显著性分清楚：**OS 那边是硬的**（低/高双比值 90% 时 1.272±0.050，5.4σ），**SS 那边单点只有
1-1.8σ**，是趋势不是已确立效应；真正硬的是两者之差（2.6-3.7σ）。也**不是精确镜像**：90% 时
OS 1.272 而 SS 的倒数是 1.05，OS 动得更大。

**不是 `Ds_pt` 失配投影到质量上的假象**（做过对照，因为 Ds_pt 既是重要性第一、又是单特征失配最
严重的那个，而 M(KKπ) 跟它对 combinatorial 候选并不独立）：固定 Ds_pt band 重算，pt>18 段
90% 时 OS 1.226±0.050 vs SS 0.978±0.084，70% 时 1.464±0.125 vs 0.784±0.163，50% 时
1.863±0.265 vs 0.650±0.240。pre-cut 时三段在 OS/SS 都是 ~1.0。注意 70% 之后只有 pt>18 段
还有统计量 —— cut 本身就在选高 Ds_pt。

**机制未知，不要硬凑。** 曾经写过"超出在 D± 峰以下 ⇒ 部分重建 charm ⇒ 电荷关联真 charm"这条
链，**它撑不住**：(a) D± 窗口 79-94% 是连续本底（峰只占 5.9%@no-cut / 21%@90%），所以"D± 窗口
~1.0 ⇒ MC 全重建 charm 没问题"是过度解读，已撤回；(b) W+jets MC **本来就有** Ds→φρ⁺/φℓν/
D*→Dπ，所以问题只能是"量不对"，而这一点没有独立证据；(c) 它完全解释不了 SS 为什么反号，也解释
不了 OS 高质量端为什么 MC 会多预言到 2.6 倍。**当前状态：现象测清楚了，机制说不出来。**

**D± 峰的 excess 抽取失败、不可用**：两种 below-peak 线性外推方案只在不切时一致
（2.42±0.63 vs 1.89±0.48），一切就分家（50% 时 0.67±0.36 vs 0.22±0.37），外推背景从峰窗密度的
89% 掉到 ~40%。跨 ~100 MeV 的线性外推在 cut 重塑低质量谱之后不可靠 —— 而那正是本研究发现它会做的事。
用上面那个 1.83-1.91 窗口比值代替（不需要外推）。对称局部 sideband 方案在这里另有一个禁忌：
上侧带会跟 blind 区重叠。

**无法做的事**：`Ds_match` 在 W+jets tuple 里全是 False（背景 tuplizer 是 `doTruth=False`），
所以不能逐事例把 MC 拆成真 Ds / combinatorial，只能靠质量 sideband 相减。

### OS − SS：把问题定位到电荷关联成分（`analysis/macros/check_os_minus_ss.py`）

减掉 SS 就消掉了一切不把 π 电荷跟 μ 电荷关联起来的本底，剩下电荷关联那部分。

- **DNN 把电荷关联成分从 13% 富集到 ~77%，MC 跟得上**：f_cc=(OS−SS)/OS，data
  0.126→0.615(90%)→0.747(70%)→0.774(50%)，MC 0.139→0.645→0.759→0.850。
- **它的总量 cut 后是对的**：D(data)/D(MC) = 1.306±0.029 → 0.904±0.035 → 0.978±0.062
  → 0.954±0.093。
- **它的质量形状不对且单调恶化**：M<1.830 是 1.14±0.07 / 1.64±0.17 / 2.15±0.37，
  M>2.008 是 0.70±0.05 / 0.65±0.07 / 0.55±0.10（高端偏离 1 达 6σ/5σ/4.5σ，每个点都硬）。
  D± 窗口一直 ~1.0。
- **减掉 SS 之后 tilt 变大**（纯 OS 低/高 1.29@90% → OS−SS 1.63）→ **tilt 就住在电荷关联
  成分里**，这是本次唯一真正的定位进展。
- **限定**：OS−SS 剩下的是"一切偏好 OS 的东西"，既含真 W+c charm 也含 combinatorial 自身的
  电荷不对称（jet 带净电荷）；不切时 f_cc 才 12.6%，未必都是真 charm。ε_sig=50% 那栏 MC
  直方图目视很尖（统计不足），读形状以 90% 为准。
- **仍然没有机制**：逐 bin 图上是横跨全质量范围的平滑斜率，不是局部鼓包，所以"部分重建 charm
  形成低质量肩膀"这个具体图像也对不上。

### OS/SS score response：上升的原因（`analysis/macros/draw_score_os_ss.py`）

归一化后的 score 谱（只看形状 —— OS 产额本来就比 SS 多 14%，不除掉会把效应埋掉），
OS vs SS，data 用 CR、W+jets 用全范围。SS data 必须用 odd 半边（模型训练用了 even 半边）；
归一化让这个一半/全量因子自动约掉。

**OS/SS 形状比值随 score 单调上升**，data 积分是 1.92(>0.9) / 2.58(>0.98)。三步定位：

1. **D± 峰毫无贡献；Ds 峰贡献很小且只有 MC 能测。** 按窗口拆开（脚本现在会打印这张表）：

   | | >0.9 | >0.98 |
   |---|---|---|
   | W+jets 全范围 | 2.324±0.041 | 3.041±0.085 |
   | W+jets 挖 Ds | 2.167±0.041 | 2.702±0.082 |
   | W+jets 挖 D± | 2.325±0.044 | 3.080±0.093 |
   | data 挖 Ds（默认） | 1.917±0.050 | 2.591±0.116 |
   | data 两个都挖 | 1.920±0.056 | 2.584±0.130 |

   挖 D± 两边都没变（data 1.917→1.920，MC 2.324→2.325）。挖 Ds 让 MC 掉约 7%，即相对 1 的
   超出掉 12% —— **即便在 MC 里也有约 88% 的上升是非峰成分**。**data 的 Ds 窗根本没法测**：
   OS 在那里 blind，任何保留 Ds 窗的行都是拿被挖空的 OS 分子比完整的 SS 分母（所以 data
   "全范围"那行读出 1.808 是假数，不可引用）。

   **踩过一次的坑**：默认 data 选择本来就已挖掉 Ds 窗（blind，且对 SS 同样施加以保持对称），
   而默认 W+jets 用全质量范围。所以拿默认 vs `--veto-peaks` 相减，等于**只从 data 挖了 D±、
   却从 MC 挖了两个峰**，得出的"data 平、MC 掉 10%"看起来像 data/MC 差异，其实是两件不同的
   比较并排放。要用逐窗口那张表，不要用两个模式相减。
2. **不是训练假象。** 隐患是结构性的：ssdata 模型训的是 signal(OS) vs SS data，所以**任何能区分
   OS/SS 位形的运动学量在训练时都是完美判别量**，与是不是信号无关，OS 本底会仅因为是 OS 而拿高分。
   对照：用 **wjets 模型**（signal(OS) vs W+jets **OS**，从未见过 OS/SS 区别）给同一批 data 打分，
   得到**同样甚至更强**的上升：2.031±0.060(>0.9) / 2.973±0.180(>0.98)，对比 ssdata 的
   1.920±0.056 / 2.584±0.130。**分类器无关 ⇒ 假象假说排除。**
   （wjets 只当量具，不是提议换回去；且只有它的 **data** 数可用 —— OS W+jets 是它的训练集、
   SS W+jets 不是，所以它的 MC OS/SS 比值被系统性压低，不能引用。）
3. 43 特征的 OS/SS 一维边缘分布基本重合，唯一可见差别是 OS 在 dR(track,μ)≈0 处的尖峰更高，
   而那是**背景样**方向，会把 OS 分数往低推。旁证，但一维边缘重合不代表尾部联合分布重合。

**结论：OS 本底里真的有一类 SS 结构上没有的 signal-like 成分，它电荷关联且不峰化。**
最自然的候选是 **W+c**（charm hadron 电荷与 W 的关联方式跟信号一样，且是挨着孤立 μ 的真 displaced
charm，在用到的每个特征上都像信号）——**但这一步仍是推断**。已测到的是三条性质：电荷关联、不峰化、
分类器无关。

**实际后果：换训练背景修不了这个**（OS-vs-OS 训的模型上升更大），只能在本底模型里处理。

### 下一步

- SR 闭合仍未解决。现在知道伴随它的质量依赖结构住在电荷关联成分里，但机制不明。
- **质量 sideband 这条路到头了**。truth 能回答，但**走不了现有 tuplizer**：它的 gen matching
  是沿 H→Ds(*)W 的 ancestry 往上走的，W+jets 里的 Ds 不来自 H，所以把 matching 从 gen filter
  解绑没有用，得另写"匹配到任意 gen Ds"的逻辑（之前那次成分研究是用 FWLite 直接读 MiniAOD 做的）。
- **先决问题**：最终本底是拟合 data sideband 还是用 MC 给形状？前者的话 MC 形状错不直接伤拟合
  （伤拟合的是 data 自己 sideband 里的 peaking 成分），后者的话上面这个斜率就是直接系统。
  这决定要不要继续追。
- `BDT/` 仍在旧的 42 特征集上，跟 DNN 不可比。

## 上一 Session (2026-08-10)

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

### 追加：W+jets 到底能不能代表 data（`analysis/macros/check_data_vs_wjets.py`）

上面"真 charm 富集"那条系统性追下去，做了个正式的 closure test——不重归一化，直接比较光度加权 W+jets 的绝对预言跟 data 的绝对事例数，分三个区间：OS CR（唯一能测的 OS 区间，SR blind）、SS CR、SS SR（SS 没 blind，唯一能测 SR 的地方）。

- **没切之前，三个区间比值一致**：1.44 / 1.46 / 1.47——缺口是常数（缺整个 process：tt/DY/QCD），能被一个 transfer factor 吸收。
- **wjets 模型切下去，比值单调爬升、三个区间开始分家**：OS CR 从 1.44 一路到 90%时2.09、70%时2.87、50%时3.83、最紧的20%时4.51±0.69——分类器专挑 MC 建模得了的那部分背景拒绝，留下的越来越不是 MC 能描述的。
- **ssdata 模型把 OS CR/SS CR 都摁回 1 附近，而且在 99%→~55% 这一大段都稳**（21点细扫：OS CR 0.93–1.05）。**但 SS SR 从来没回到过 1**，同样随 cut 单调恶化：1.19(99%)→1.33(95%)→1.64(70%)→~2.0(55–50%)。**换背景解决了 sideband 闭合，没解决 signal region 闭合**——这才是更要命的那个数字，之前的"背景低36-48%"不能读成"背景问题解决了"。
- ε_sig<45% 之后 SS SR 统计量掉到 <15 个事例，误差 >30-50%，那之后的点不能读结构。
- **逐 bin 图上一个反直觉的坑**：远离 Ds/D± 峰的地方比值贴着1，个别高质量端的散点"翘起来"看着像效应，实测把那几个 bin 加总起来算真实统计量反而是 0.80（低于1不是高于1）——单 bin 统计量太小时，眼睛会被涨落骗，加总回去再看。

结论一句话：**换 SS 背景修好了 sideband，没修好 signal region**，这个是下一步要处理的真正问题，不是"already better"。

### 下一步没做的事

- **SR closure 没解决**：见上一节，是当前最大的开放问题。方向大致两条：本底模型里显式加真 charm 成分，或者削弱分类器对真 charm 的偏好（也可能两条都要）。
- `BDT/` 仍在旧的 42 特征集上，跟 DNN 不可比。

> 之前几次 session 的详细记录见 memory（`[[project_mva_test_build_v1]]`、`[[project_h_phipimu_ss_tuplizer]]`、`[[project_mva_sculpting_finding]]` 等）。2026-07-23 那次的 SS tuplizer / Tuples 改名 / 新信号样本上线记录见 `[[project_h_phipimu_ss_tuplizer]]` 和 git history。
