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

## Last Session (2026-07-16)

### association vs geometric truth-matching（phipimu）

问题：现在 truth matching 用 `packedPFCandidateToGenAssociation`，换成 geometric matching to genParticle 会更好还是更差？

**结论：不要换。工作点上效率与 association 统计上无法区分，而 geometric 白白多出纯度代价。**

三元组（三条腿全中）层级，两条臂跑在同一批事例上（相关 ⇒ 只有 discordant 带信息，用 McNemar）：

| 通道 | assoc-only | geo-only | geo − assoc |
|---|---|---|---|
| DsW | 58 | 53 | **-0.14 ± 0.30 pp**（z=-0.47） |
| DsstW | 68 | 67 | **-0.03 ± 0.33 pp**（z=-0.09） |

⚠️ 别把 -0.1pp 说成"geometric 输"——统计上就是**打平**。真正的理由是打平之余还多出代价。

**源码核实的关键事实**（非推测）：
- **association 是 hit-based MC truth**：`packedPFCandidateToGenAssociation` → `prunedTrackMCMatch`(`MCTrackMatcher`) → `quickPrunedTrackAssociatorByHits`（cluster→TrackingParticle→genParticle）。geometric 没有 sim hit，一般不可能更正确。
- **Association 建在 `prunedGenParticles`**（`PackedCandidateGenAssociationProducer.cc:78,85,134`）＝ tuplizer `GenPartToken_` 的同一集合 ⇒ `genRef.key()` 可直接比 `Gen_<leg>_idx`。**实测 ~6000 条腿 0 mismatch**。
- **TP-pruning 假说错了**：`prunedTrackingParticles` = `keep++ (400<|pdgId|<600)` + `drop status!=1` + `drop charge==0` ⇒ Ds/Ds* 后代（K±/π±）**按设计保留**。
- **assoc 的 `dR<0.05` 门完全不起作用**：hit-based 的 dR 最大只有 ~0.012。只有 `|dpt|<0.3` 切掉 ~2.5%。

**geometric 为什么赢不了**：裸效率 geo 94.5% vs assoc 85%（+9pp），但加质量门后收益全蒸发（ungated +9pp → dR<0.05 +1.4pp → 再加 |dpt|<0.3 ≈ **0**）。
决定性证据：多出的 9pp 里 ~70% 是 `assoc_key == null`，**该桶 dpt 均值 -0.65~-0.78**（真匹配是 -0.01，RMS 0.11）——geo 捡到的 PF 只带 gen 腿 ~25% 的 pT，根本不是同一个粒子。真腿没被重建成 PF 时，geometric 不返回"无匹配"，而是抓最近的软 PU/UE 径迹，几乎全被 `|dpt|<0.3` 杀掉。
⟹ **association 那 16% 缺失是真实重建低效，不是记账假象**（与 lostTracks 研究一致：PF 漏掉的腿只有 7.8% 能在 lostTracks 找回）。

**换过去的额外代价**：纯度——assoc 成功的腿上 geo 有 **0.9-1.8%** 选中不同 PF；**同号撞车 7-8%**（`Kp_geo_idx==pi_geo_idx`，K+/π 同号，association 靠真身份不会犯）。

**geo 臂没写错**：两臂都匹配的腿上，geo 与 hit-based truth 选同一个 PF 达 **98.2-99.1%**——有 bug 不可能符合到 99%。

产物：`plugins/GeoMatch_phipimu.cc` + `interface/GeoMatch_phipimu_Tree.h` + `src/GeoMatch_phipimu_Tree.cc`（由 `LostTrack_phipimu` **完整复制**而来）、`test/SIG_{,DsstW_}GeoMatch_phipimu.py`、`Tuples/phipimu/Signal_GeoMatch/geomatch.C`（一个宏出全部数字）。

**设计要点（可复用）**：两条臂同树 + **ungated 记录**（门离线加）⟹ 一次跑扫遍工作点、零重跑。assoc 臂去门后离线重加门，计数与原版 `LostTrack_phipimu` **逐项相等**（101/102/101），证明语义没被动过。

**没做**：没提 condor（本地 25 文件 × 双通道 ≈ 7000 事例，且差值已 consistent with zero，加统计动不了结论）；计划里的 M(KKπ) 闭合图没做——dpt 残差已经是不依赖 association 的判据，够定论了，除非要写 note/talk。

### lostTracks recovery 研究（phipimu）

问题：MiniAOD 里有 `lostTracks` 和 `lostTracksToGenAssociation`，把它们加进 Ds 的径迹池能回收多少信号、组合本底涨多少？

**结论：不做。S/√B 净降约 1%**（信号 +1.7% 且质量差 2 倍，真本底 +5.2%）。

两个新分析器，都由 `H_phipimu` **完整复制**而来再删改（不手写重实现，见 memory `feedback_copy_original_not_rewrite`）：

- **`LostTrack_phipimu`**（gen 上限，砍掉 reco 起；不加载 JMEcorr）：每条腿 PF 匹配率 ~84%，PF 缺失的腿只有 **7.8%** 能在 lostTracks 找回；池子级上限 **+2.0%**。匹配质量：lostTracks 的 Δpt/pt σ68 = **6.5%** vs PF **1.5%**（差 4.3 倍）。
- **`LostTrackReco_phipimu`**（完整 reco，保留 JMEcorr）：信号 **+1.73%**（DsW）/ **+1.40%**（DsstW）；组合本底 **+3.72%**；pT 排序纯度 99.925%→99.917%（**没被破坏**，现有选候选策略安全）；含 LT 候选的 M(KKπ) 分辨差 **2.11 倍**。

**lostTracks 关键事实**：约 30/事例，`hasTrackDetails()` 恒 100%，pdgId 只有 ±211 和极少 ±13——**没有光子**（它装的是 generalTracks 里没被 PF 用掉的径迹，结构上不可能有中性粒子）。所以对 DsstW 的 gamma 完全帮不上忙。pt>1 的只有 ~66% 过 `trackHighPurity`。

**双来源池子用边界法**：`idx_pf` 类型/名字不变、仍存原生下标；PF 先 push 并记 `const Int_t n_pf_pool = idx_pf.size()`，再 push lostTracks；`Kp_isLost = (iKp >= n_pf_pool)`。**池子绝不可重排**（如按 pT 排序）否则身份判定静默崩坏。必改三处：
1. `idx_pf.size()<3` 挪到**合并后**判（否则 2PF+1LT 的 recovery 事例被提前扔掉）
2. 孤立度排除加 `!isLost` 保护（原版 `i==idx_Kp` 在双来源下会误伤 PF 里同号候选）；**求和仍只用 PF**，鉴别变量不能中途改定义
3. `Ds_match` 每条腿按各自集合比对

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
