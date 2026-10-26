"""Does adding DY mu-from-tau to the stage-2 training help? resdy vs wjetsres_v2 (2026-10-26).

The two BDTs differ ONLY in the background: resdy = W+jets modes 1-5 + DY (Muon_genPartFlav
15), wjetsres_v2 = W+jets modes 1-5. Same features (prodtopo_v2), same stage-1 pre-filter
(ssdata_f45 eps90), same weights, seed and split recipe.

Compared at MATCHED eps_tot (stage-1 eps90 AND stage-2, counted directly on truth-matched
signal, weighted): the targets are resdy's WP2 90/70/50% (thr2 = unweighted quantile, as in
draw_stage1_mass_by_score_cut.py); the control's thr2 is chosen to reach the same eps_tot.

Rows per column (MC x br_reweight, DY also x dy_jetbin_sf):
  DY mu<-tau, test split : DY flav-15 rows resdy did NOT train on (its 20% test split,
                           reproduced and checked bit-for-bit against its X_*.npy). All
                           masses (the training population), and the Ds SR.
  DY mu<-tau, all rows   : in-sample for resdy -> biased in its favour, reference only.
  TT 2L2Nu mode 0        : tau 3-prong in ttbar, neither model trained on it.
  W+jets modes 1-5 (SR)  : the resonant rejection stage-2 exists for (in-sample for both).
  Data OS low sideband   : M < 1.830 GeV; Data OS sideband = all non-blind rows.
  Punzi                  : eps_tot / (1.5 + sqrt(B)), B = W+jets + DY + TT in the Ds SR.
Ratio = resdy / control (< 1: resdy removes more). Errors split the two kept sets into
both / resdy-only / control-only pieces (independent), sum w^2 per piece.

    python3 compare_stage2_dy_vs_nody.py
"""
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "DNN"))
import common

K_A, K_B = "resdy", "wjetsres_v2"
M_DS, SR_HALF, LOW_HI = 1.96835, 0.040, 1.830
PHI_WIN = 0.020
SEED = 42
cfgA = common.resolve(K_A)
assert {k: v for k, v in cfgA.items() if k not in ("samples", "label", "out_name", "out_dir")} == \
       {k: v for k, v in common.resolve(K_B).items() if k not in ("samples", "label", "out_name", "out_dir")}
DIR_A, DIR_B = (common.out_dir_for(k, tool="BDT") for k in (K_A, K_B))
FLAV = np.load(f"{common.MVA_ROOT}/aux/DY_merged_mu_genPartFlav.npz")
_o = np.argsort(FLAV["rowkey"])
FK, FV = FLAV["rowkey"][_o], FLAV["flav"][_o]


def rowkey(a):
    return (a["luminosityBlock"].astype(np.int64) << 32) | a["event"].astype(np.int64)


def flav_of(a):
    rk = rowkey(a)
    pos = np.minimum(np.searchsorted(FK, rk), len(FK) - 1)
    return np.where(FK[pos] == rk, FV[pos], -1)


# ---------- 1. reproduce resdy's test split; DY test-row keys ----------
thr1, _ = common.stage1_threshold(cfgA["stage1_cut"])
S1 = dict(cfgA["stage1_cut"], thr=thr1)
feats = cfgA["features"]
xs, _ = common.load("Signal_merged", feats=feats, truth_filter=True, stage1=S1, extra_branches=("weight",))
xw, _ = common.load("WtoLNu2Jets_merged", feats=feats, gen_mode_filter=(1, 2, 3, 4, 5), stage1=S1,
                    extra_branches=("weight",))
xd, exd = common.load("DY_merged", feats=feats, stage1=S1,
                      extra_branches=("run", "luminosityBlock", "event"))
keep = flav_of(exd) == 15
xd, dykey = xd[keep], rowkey(exd)[keep]
X = np.vstack([xs, xw, xd])
y = np.concatenate([np.ones(len(xs)), np.zeros(len(xw) + len(xd))])
src = np.concatenate([np.zeros(len(xs)), np.ones(len(xw)), np.full(len(xd), 2)])
dyk_all = np.concatenate([np.zeros(len(xs) + len(xw), np.int64), dykey])
X, m = common.preprocess(X, feats)
X, y, src, dyk_all = X[m], y[m], src[m], dyk_all[m]
# log columns differ by ~1e-16 between the condor node and here (libm), so not bit-for-bit
for _c, _f in ((1, "X_sig.npy"), (0, "X_bkg.npy")):
    _ref = np.load(f"{DIR_A}/{_f}")
    assert _ref.shape == X[y == _c].shape and np.allclose(X[y == _c], _ref, rtol=0, atol=1e-12), \
        f"matrix != resdy {_f}"
idx = np.random.default_rng(SEED).permutation(len(y))
test_i = idx[:int(0.2 * len(y))]
assert np.array_equal(y[test_i], np.load(f"{DIR_A}/y_test.npy")), "test split != resdy y_test.npy"
DY_TEST = set(dyk_all[test_i][src[test_i] == 2].tolist())
print(f"resdy split reproduced: {len(xd):,} DY flav-15 rows, {len(DY_TEST):,} in its test split")


# ---------- 2. read both Prediction files, aligned ----------
SEL = list(common.SELECTION)
BASE = ["Ds_mass", "phi_mass", "score_dsmu", "score_dsmu_res", "weight", "run", "luminosityBlock",
        "event", "cand_gen_mode", "cand_gen_common_hf_pdgId"] + SEL


def read(name, extra=()):
    br = list(dict.fromkeys(BASE + list(extra)))
    is_mc_bkg = not (name.startswith("Data") or name.startswith("Signal"))
    if not is_mc_bkg:   # no truth tag; signal is one BR by construction (Punzi is BR-free)
        br = [b for b in br if not b.startswith("cand_gen")]
    a = uproot.open(f"{DIR_A}/Prediction/{name}.root")["Events"].arrays(br, library="np")
    b = uproot.open(f"{DIR_B}/Prediction/{name}.root")["Events"].arrays(
        ["score_dsmu", "score_dsmu_res", "run", "luminosityBlock", "event"], library="np")
    for k in ("run", "luminosityBlock", "event", "score_dsmu"):
        assert np.array_equal(a[k], b[k]), f"{name}: {k} differs between the two Prediction files"
    a["s2A"], a["s2B"] = a.pop("score_dsmu_res"), b["score_dsmu_res"]
    m = (common.selection_mask(a) & (a["score_dsmu"] > -1) & (a["s2A"] > -1) & (a["s2B"] > -1)
         & (np.abs(a["phi_mass"] - common.M_PHI) < PHI_WIN))
    a = {k: v[m] for k, v in a.items()}
    w = a["weight"].astype(np.float64)
    if is_mc_bkg:
        w = w * common.br_reweight(a)
    if name.startswith("DY"):
        w = w * common.dy_jetbin_sf(a)
    a["w"] = w
    a["sr"] = np.abs(a["Ds_mass"] - M_DS) < SR_HALF
    print(f"  {name}: {m.sum():,} rows")
    return a


sig = read("Signal_merged", ("Ds_match", "mu_match"))
sig = {k: v[sig["Ds_match"].astype(bool) & sig["mu_match"].astype(bool)] for k, v in sig.items()}
wj = read("WtoLNu2Jets_merged")
dy = read("DY_merged", ("jet_bin",))
tt = read("TT_merged", ("tt_decay",))
do = read("Data_OS_merged")
dy["flav"] = flav_of(dy)
dy["test"] = np.isin(rowkey(dy), np.fromiter(DY_TEST, np.int64))
n_found = (dy["test"] & (dy["flav"] == 15) & (dy["score_dsmu"] > thr1)).sum()
assert n_found == len(DY_TEST), f"DY test rows found {n_found} != {len(DY_TEST)}"

# ---------- 3. matched eps_tot thresholds ----------
s1p = sig["score_dsmu"] > thr1
W = sig["w"].sum()


def thr_for_eps(s2, eps):
    """Stage-2 cut giving weighted eps_tot = eps on signal passing stage-1."""
    v, ww = s2[s1p], sig["w"][s1p]
    o = np.argsort(-v)
    c = np.cumsum(ww[o]) / W
    return float(v[o][np.searchsorted(c, eps)])


COLS = [("stage-1 eps90 only", -np.inf, -np.inf)]
for e in (0.90, 0.70, 0.50):
    tA = float(np.quantile(sig["s2A"][sig["score_dsmu"] > -1], 1 - e))
    eps = sig["w"][s1p & (sig["s2A"] > tA)].sum() / W
    tB = thr_for_eps(sig["s2B"], eps)
    epsB = sig["w"][s1p & (sig["s2B"] > tB)].sum() / W
    COLS.append((f"WP2 {e:.0%} (eps_tot {eps:.3f}/{epsB:.3f})", tA, tB))


def passes(a, tA, tB):
    s1 = a["score_dsmu"] > thr1
    return s1 & (a["s2A"] > tA), s1 & (a["s2B"] > tB)


def ratio(a, sel, tA, tB):
    pA, pB = passes(a, tA, tB)
    pA, pB = pA & sel, pB & sel
    w = a["w"]
    C, Ao, Bo = (w[pA & pB], w[pA & ~pB], w[~pA & pB])
    nA, nB = C.sum() + Ao.sum(), C.sum() + Bo.sum()
    vC, vA, vB = (x @ x for x in (C, Ao, Bo))
    R = nA / nB
    err = np.sqrt(vA / nB**2 + vB * (nA / nB**2)**2 + vC * ((Bo.sum() - Ao.sum()) / nB**2)**2)
    return nA, nB, R, err


ROWS = [
    ("DY mu<-tau test split, all M", dy, (dy["flav"] == 15) & dy["test"]),
    ("DY mu<-tau test split, SR", dy, (dy["flav"] == 15) & dy["test"] & dy["sr"]),
    ("DY mu<-tau all rows, SR (in-sample A)", dy, (dy["flav"] == 15) & dy["sr"]),
    ("DY prompt mu (flav 1), SR", dy, (dy["flav"] == 1) & dy["sr"]),
    ("TT 2L2Nu mode 0, all M", tt, (tt["tt_decay"] == 0) & (tt["cand_gen_mode"] == 0)),
    ("TT 2L2Nu mode 0, SR", tt, (tt["tt_decay"] == 0) & (tt["cand_gen_mode"] == 0) & tt["sr"]),
    ("W+jets modes 1-5, SR", wj, np.isin(wj["cand_gen_mode"], [1, 2, 3, 4, 5]) & wj["sr"]),
    ("W+jets mode 0, SR", wj, (wj["cand_gen_mode"] == 0) & wj["sr"]),
    ("Data OS low sideband", do, do["Ds_mass"] < LOW_HI),
    ("Data OS sideband (all)", do, np.ones(len(do["w"]), bool)),
]

print(f"\nresdy (A) vs wjetsres_v2 (B); stage-1 thr {thr1:.6f}; ratio = A/B\n")
for lab, tA, tB in COLS:
    print(f"== {lab}  thr2 A {tA:.4f}  B {tB:.4f}")
    for rl, a, sel in ROWS:
        nA, nB, R, err = ratio(a, sel, tA, tB)
        print(f"   {rl:40s} A {nA:10.2f}  B {nB:10.2f}  A/B {R:.3f} +- {err:.3f}")
    eA = sig["w"][passes(sig, tA, tB)[0]].sum() / W
    eB = sig["w"][passes(sig, tA, tB)[1]].sum() / W
    BA = sum(a["w"][passes(a, tA, tB)[0] & a["sr"]].sum() for a in (wj, dy, tt))
    BB = sum(a["w"][passes(a, tA, tB)[1] & a["sr"]].sum() for a in (wj, dy, tt))
    pA, pB = eA / (1.5 + np.sqrt(max(BA, 0))), eB / (1.5 + np.sqrt(max(BB, 0)))
    print(f"   {'Punzi (B = all MC in SR)':40s} A {pA:.5f} (B={BA:.1f})  B {pB:.5f} (B={BB:.1f})  "
          f"A/B {pA / pB:.3f}")
