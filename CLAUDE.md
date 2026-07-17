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

Each subdirectory under `Tuples/phipimu/` corresponds to one sample. To submit:
```bash
cd Tuples/phipimu/Signal/
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
- H(125) → Ds±(→φπ± → K+K−π±) + W∓(→μ∓ν)
- H(125) → Ds*±(→Ds±γ) + W∓(→μ∓ν)

Note the final state is the same K+K−π± either way; what makes this the φπ mode is
the φ requirement — the gen filter demands the φ ancestry and the candidate loop
applies a ±30 MeV window on M(KK).

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

1. **Gen filter** (`doTruth` only): require exactly 1H + 1Ds + 1W + 1μ + 0e (hard process), then find the status-1 gen K+, K−, π with the full φ ancestry chain (`TwoBodyDecay(gp, 333, 321, -321)` etc.). Returns early if not found.
2. **Gen matching**: match PF candidates via `packedPFCandidateToGenAssociation` + ΔR<0.05 + Δpt/pt<0.3. The γ (DsstW) is matched geometrically instead — ΔR<0.02, Δpt/pt<0.5 — because it has no association. Muon via `muon.genParticle()`.
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

## Last Session (2026-07-16)

### association vs geometric truth-matching（phipimu）

> 代码已移出本 branch，全部在 **`study/truthmatch-lostrack`**（已 push）。
> `git checkout study/truthmatch-lostrack` 即可取回。

问题：现在 truth matching 用 `packedPFCandidateToGenAssociation`，换成 geometric matching to genParticle 会更好还是更差？

**结论：两者在同一 working point 下等价。建议不换，理由只是「等价、没必要动 association」**（它是 hit-based truth，更有原则，且是现成的）——**不是** geometric 有害。这是个不强的理由。

三元组（K+/K−/π 全 match）层级。两种方法跑在**同一批 event** 上（相关 ⇒ 只有 discordant 带信息，用 McNemar）：

| 通道 | assoc-only | geo-only | geo − assoc |
|---|---|---|---|
| DsW | 58 | 57 | −0.03 ± 0.31 pp（z=−0.09） |
| DsstW | 65 | 68 | +0.08 ± 0.33 pp（z=+0.26） |

#### ⚠️ geometric matching 的正确写法（第一版写错了）

**正确**（= 原版 `H_phipimu` gamma matching 的写法，`H_phipimu.cc:548-556`）：
```
对每个 gen 粒子：
  收集所有同电荷 charged-hadron PF 中 dR<0.05 且 |dpt|<0.3 的
  在通过的候选里取 min dR
```
两个 cut 都在 loop **里面**。

**错的写法**（第一版）：dR<0.3 里取 min dR（loop 里没有 dpt cut），拿到唯一候选，离线再加 cut。
⇒ 最近的那条是软 track 时，会选中它然后离线被砍 → 记为 no match；正确写法会跳过它、继续找次近的合格 track → match 上。**人为削弱了 geometric**。

**教训**：cut 决定「选哪个候选」时**不能**推迟到离线。ungated 记录只对 association 成立（它靠 map 认粒子，cut 纯属事后质量检查）；套到 geometric 上逻辑不成立。这个错误也无法离线补救——每个粒子只存了一个候选。

#### 修正后的数字（第一版全部作废）

| | 错的写法 | 正确写法 |
|---|---|---|
| purity cost | 0.9–1.8% | **0.00–0.03%**（一致率 99.97–100%） |
| K+/π 撞同一条 track | 7–8% | **0.06–0.17%** |
| assoc-null bucket dpt | −0.65~−0.78（据此误判为软垃圾） | **−0.01~−0.07, RMS 0.11 = 和真 match 无异** |
| bucket 大小 | 210–333 | 22–32 |

**原先「证伪」的假说其实是对的**：assoc-null bucket 里那 ~25 条 **是真的 track**，是 association 被 `quickTrackAssociatorByHits` 的 `Cut_RecoToSim=0.75`（须与 TrackingParticle 共享 >75% hits，否则完全没有 association）丢掉的。只是量少（~0.7%），而 geometric 自己也漏掉差不多同样多，净值为零。

#### 源码核实的事实（CMSSW 15.0.5，非推测）
- **association 是 hit-based MC truth**：`packedPFCandidateToGenAssociation` → `prunedTrackMCMatch`(`MCTrackMatcher`) → `quickPrunedTrackAssociatorByHits`（cluster→TrackingParticle→genParticle）。geometric 没有 sim hit，一般不可能更正确。
- **Association 建在 `prunedGenParticles`**（`PackedCandidateGenAssociationProducer.cc:78,85,134`）＝ tuplizer `GenPartToken_` 的同一集合 ⇒ `genRef.key()` 可直接比 `Gen_<粒子>_idx`。**实测 ~6000 条 0 mismatch**。
- **TP-pruning 假说错了**：`prunedTrackingParticles` = `keep++ (400<|pdgId|<600)` + `drop status!=1` + `drop charge==0` ⇒ Ds/Ds* 后代（K±/π±）**按设计保留**。
- **assoc 的 `dR<0.05` 完全不起作用**：hit-based 的 dR 最大只有 ~0.012。只有 `|dpt|<0.3` 去掉 ~2.5%。
- **geometric 必须 gen 驱动**：`prunedGenParticles` 有 `drop status==1` + 选择性 re-keep，不是完整末态清单。若写成 PF 驱动（PF→最近 gen）必须改用 `packedGenParticles`，否则会把 PF 错判给我们的粒子、系统性高估 geo 效率。

**验证**：key 对拍 0 mismatch；association 离线重加 cut 后与原版 `LostTrack_phipimu` **逐项相等**（101/102/101），证明语义没被动过；geometric 与 hit-based truth 一致率 99.97–100%（新代码没写错的凭据）。

**没做**：没提 condor（本地 25 file × 双通道 ≈ 7000 event，差值已 consistent with zero）；M(KKπ) 闭合图没做（dpt 残差已够定论，除非要写 note/talk）。

### lostTracks recovery 研究（phipimu）

> 代码已移出本 branch，全部在 **`study/truthmatch-lostrack`**（已 push）。

问题：MiniAOD 里有 `lostTracks` 和 `lostTracksToGenAssociation`，把它们加进 Ds 的 track 池能回收多少信号、组合本底涨多少？

**结论：不做。S/√B 净降约 1%**（信号 +1.7% 且质量差 2 倍，真本底 +5.2%）。

两个新分析器，都由 `H_phipimu` **完整复制**而来再删改（不手写重实现，见 memory `feedback_copy_original_not_rewrite`）：

- **`LostTrack_phipimu`**（gen 上限，砍掉 reco 起；不加载 JMEcorr）：每个粒子的 PF match rate ~84%，PF 没 match 上的粒子只有 **7.8%** 能在 lostTracks 找回；池子层面的上限 **+2.0%**。匹配质量：lostTracks 的 Δpt/pt σ68 = **6.5%** vs PF **1.5%**（差 4.3 倍）。
- **`LostTrackReco_phipimu`**（完整 reco，保留 JMEcorr）：信号 **+1.73%**（DsW）/ **+1.40%**（DsstW）；组合本底 **+3.72%**；pT 排序纯度 99.925%→99.917%（**没被破坏**，现有选候选策略安全）；含 LT 候选的 M(KKπ) 分辨差 **2.11 倍**。

**lostTracks 关键事实**：约 30/事例，`hasTrackDetails()` 恒 100%，pdgId 只有 ±211 和极少 ±13——**没有光子**（它装的是 generalTracks 里没被 PF 用掉的径迹，结构上不可能有中性粒子）。所以对 DsstW 的 gamma 完全帮不上忙。pt>1 的只有 ~66% 过 `trackHighPurity`。

**双来源池子用边界法**：`idx_pf` 类型/名字不变、仍存原生下标；PF 先 push 并记 `const Int_t n_pf_pool = idx_pf.size()`，再 push lostTracks；`Kp_isLost = (iKp >= n_pf_pool)`。**池子绝不可重排**（如按 pT 排序）否则身份判定静默崩坏。必改三处：
1. `idx_pf.size()<3` 挪到**合并后**判（否则 2PF+1LT 的 recovery 事例被提前扔掉）
2. 孤立度排除加 `!isLost` 保护（原版 `i==idx_Kp` 在双来源下会误伤 PF 里同号候选）；**求和仍只用 PF**，鉴别变量不能中途改定义
3. `Ds_match` 每个粒子按各自集合比对

**不变式**：`n_lost==0` 的候选 = 只用 PF 池子的候选（已对同文件的原版 H_phipimu 逐项验证：80 行/79 事例/74 gen-matched 全等）。所以**跑一遍同时得到 PF-only 和 PF+LT**。

**最终裁决**：WJets 真本底（CRAB `260716_135152`，99/99 成功，产出在 IIHE `/store/user/kakang/.../WJets_LostTrackReco_phipimu/output_{1..99}.root`）：去重事例 **39995 → 42067，+5.18% ± 0.11%**——比信号事例内的假候选 +3.7% 还高，是信号增幅的 3 倍。

| 道 | S | B | S/√B |
|---|---|---|---|
| DsW | ×1.0173 | ×1.0518 | **0.992（−0.8%）** |
| DsstW | ×1.0140 | ×1.0518 | **0.989（−1.1%）** |

两个道都变差约 1%，且未计入多出信号质量差 2 倍的贬值。WJets 占加权本底 77%，可代表总本底。本底 yield = **事例数**（每事例只留一个候选），不是候选数。

**另一条线（与 lostTracks 无关，仍开着）**：DsstW 的真正瓶颈是**光子**——gamma PF 匹配率仅 **28%**，vs 带电腿 84%。要提效率得往光子重建找（放宽 gamma 匹配的 dR<0.02 / Δpt/pt<0.5，或看 slimmedPhotons / ECAL）。

图：`/eos/.../Signal_LostTrack_phipimu/match_quality_dpt.png`、`/eos/.../Signal_LostTrackReco_phipimu/{pt,mass}_resolution.png`

### phipimu MVA 训练配比研究：v1 / v2 / v3

三个 version 共用同一套 **37 特征**（旧 46 特征会 sculpt 拟合变量，已废弃）和同一份代码，只有 `CONFIG[model_name]` 的背景组成 + 采样模式不同。旧 v1/v2/v3 产物已删除，原 v4(flat) 改名为 v3。

| version | 训练背景 | 采样 |
|---------|----------|------|
| v1 | 只 WtoLNu_2Jets | — |
| **v2** | 4 个 | **proportional**（38715:3418:3570:2024） |
| v3 | 4 个 | flat（等量） |

4 个训练背景 = WtoLNu_2Jets, TTtoLNu2Q, TTto2L2Nu, DYto2L_MLL50。
产物：`/eos/home-k/kakang/H2DW/Analysis/MVA/{BDT,phipimu}/{v1,v2,v3}/`，各 36 个 Prediction。

#### 裁决：**v2 (proportional) + BDT 最优**

指标 = 全 34 背景加权残留（基线 B0=14,444,602）。BDT 残留 @ 50%/30%/10% 信号效率：

| version | 50% | 30% | 10% |
|---------|-----|-----|-----|
| v1 | 20,823 | 7,498 | 2,077 |
| **v2** | **11,417** | **3,105** | **788** |
| v3 | 14,580 | 3,970 | 940 |

排序稳定 **v2 < v3 < v1**；BDT 全面优于 DNN。

**反直觉的机制**：per-background 看 v3(flat) 几乎压每个背景都最狠，唯独 W+jets 是 v1 最狠（0.041% vs v3 0.114%）。但 W+jets 占加权背景 **77%**（基线 11,147,414），v3 在它上面多放走的量超过在 TT/DY 上省下的量 → 总残留输给 v2。要优化的是加权总背景，不是 per-sample 存活率。

#### 对比图
- `/eos/home-k/kakang/H2DW/Analysis/MVA/version_comparison.png`（双面板 BDT|DNN）
- `/eos/home-k/kakang/H2DW/Analysis/MVA/version_comparison_overlay.png`（6 条同图，实线 BDT / 虚线 DNN）
- x = signal efficiency %，y = background retention %（log）

#### 未完成
- v2-BDT 的 mass sculpting 复查未做
- 本轮代码改动（CONFIG 三模式 + proportional 采样 + condor 参数化）**未 commit**
