"""Single source of truth for the phipimu DNN: feature set, preprocessing chain,
network definition, and the background-configuration table.

Why this module exists: FEATURE_NAMES used to be copy-pasted into DNN.py and
predict.py, kept in sync only by a comment saying "must stay identical". They
drifted anyway -- the 40-feature v2 (W+jets background) and the 43-feature v3 (SS
data background) were compared head-to-head on 2026-08-10 as if the only
difference were the training background, which invalidated the comparison. Import
from here instead; do not re-declare any of it downstream.

BDT/ is deliberately NOT wired into this module. It is pinned to an older
42-feature set and is already declared non-comparable to the DNN; importing from
here would silently change the BDT the next time a feature is added or reordered.

FEATURE ORDER IS LOAD-BEARING. It is positionally bound to scaler_mean.npy,
scaler_scale.npy, and nn.Linear(n_in, 128). Reordering the list without retraining
produces garbage scores and raises no exception -- check_model_dir() below is what
turns that into a readable failure.
"""

import os
import sys

import numpy as np
import torch
import torch.nn as nn
import uproot

SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples"
MVA_DIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/DNN"
TREE = "Events"

# ========= 特征定义 =========
# Feature set excludes every decay-kinematics variable: those reconstruct
# phi_mass / Ds_mass to R^2 = 0.98 / 0.94, which would sculpt the resonance peaks
# used for the final fit. Dsmu_mass is a feature here, not a fit variable -- it has
# no resonant peak, and keeping it out of the fit is what lets Ds_pt stay in the MVA.
FEATURE_NAMES = [
    # lifetime / vertex -- prompt-from-H vs displaced-from-B (the real-Ds background)
    "Ds_ipsig",
    "Ds_FDsig",
    "dira",
    "Ds_prob",
    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "mu_ipsig",
    # 2D lifetime/vertex: transverse resolution beats longitudinal.
    # FDxysig separates DY (fake Ds from prompt tracks -> no flight distance);
    # ipxysig separates TT (real Ds from a displaced B -> momentum misses the PV).
    "Ds_FDxysig",
    "Ds_ipxysig",
    # Per-track 2D IP significance. Added 2026-08-04: the same "transverse beats
    # longitudinal" argument that motivated Ds_ipxysig applies to the tracks, and it
    # had simply never been carried down to them. Measured on this sample they beat
    # their 3D counterparts (AUC 0.725/0.723/0.780 vs 0.710/0.709/0.774) and are only
    # weakly correlated with them (r = 0.14/0.19/0.27) -- the 3D versions are diluted
    # by the much worse z resolution, so the two carry largely complementary
    # information rather than being two views of the same quantity.
    # mu_ipxysig deliberately left out: AUC 0.504, i.e. nothing -- the muon comes from
    # a prompt W and has no displacement to measure (its 3D version, mu_ipsig, is
    # equally useless at 0.503 but is kept for continuity with the earlier feature set).
    "Kp_ipxysig",
    "Km_ipxysig",
    "pi_ipxysig",
    # isolation -- the components, not just RelIso. sumChargedHadronPt is the single
    # strongest variable against every background (AUC 0.91 vs TT): a Ds buried in a
    # b jet is surrounded by charged hadrons. RelIso divides by Ds_pt, which dilutes it.
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
    "Ds_IsoR04_RelIso",
    "Ds_IsoR03_RelIso",
    "mu_IsoR04_RelIso",
    # MET / event. num_jets_pt20..pt60 are near-duplicates (AUC all ~0.80), so only
    # one threshold is kept; num_jets (no threshold) adds the soft-jet activity.
    "met_pt",
    "met_sig",
    "met_sumEt",
    "W_tmass",
    "num_jets",
    "num_jets_pt30",
    # Ds + mu kinematics -- hard from a 125 GeV H vs soft inside a jet
    "Ds_pt",
    "Dsmu_pt",
    "Dsmu_eta",
    "mu_pt",
    "mu_eta",
    "dR_Ds_mu",
    "dR_Kp_mu",
    "dR_Km_mu",
    "dR_pi_mu",
    "dR_phi_mu",
    # Ds+mu mass: H-mass constrained for signal, continuum for background
    "Dsmu_mass",
    # --- jet/b-tag environment, added 2026-08-03, re-screened 2026-08-10 for the
    # SS-background switch (screen_new_jet_vars.py, b-veto applied, SS as the
    # reference instead of W+jets MC -- see new_jet_vars_screening.png).
    "Ds_nearestjet_bscore",
    "Ds_nearestjet_dR",
    "jet_leadpt",
    # jet_ht: against W+jets MC this scored ~0 importance and was dropped 2026-08-04
    # -- W+jets has no tt-like hadronic activity to separate from. Against SS data
    # it keeps a real (if reduced post-veto) shape difference from signal, so it's
    # back in. jet_ht_pt20 is the same variable at a second threshold, dropped as a
    # near-duplicate (same precedent as num_jets_pt20..pt60 above).
    "jet_ht",
    # n_b_jets_{loose,medium}_dR04: dropped 2026-08-04 after ranking dead last
    # (41st/42nd of 42, importance 0.0016/0.0020) against W+jets MC -- which has no
    # extra soft b-jet activity to distinguish. Re-added: against SS (the actual
    # candidate background now) they show a real gap in the 0-jet vs 1-jet fraction
    # that never showed up against W+jets. jet_maxbscore was also re-checked and
    # dropped again -- its pre-veto separation turned out to come almost entirely
    # from events the tight b-veto already removes, so little is left post-veto.
    "n_b_jets_loose_dR04",
    "n_b_jets_medium_dR04",
]

LOG_FEATURES = [
    "Kp_ipsig",
    "Km_ipsig",
    "pi_ipsig",
    "Ds_ipsig",
    "mu_ipsig",
    "Ds_FDsig",
    "Ds_FDxysig",
    "Ds_ipxysig",
    "Kp_ipxysig",
    "Km_ipxysig",
    "pi_ipxysig",
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
    "jet_ht",
]

# Sentinel handling -- needed here but NOT in the BDT: trees just split around a
# -777, but StandardScaler would let it wreck the mean/variance of the whole
# feature. Each sentinel is mapped to the physically sensible edge of its range
# ("no jet nearby" / "tagger gave nothing"), keeping the input bounded.
#   Ds_nearestjet_bscore  -777 (no jet) and ~-10 (tagger unavailable, ~14% of rows)
#   Ds_nearestjet_dR      -777 (no jet)  -> far away
#   jet_leadpt            -777 (no jet)  -> zero pt
SENTINEL_FIX = {
    "Ds_nearestjet_bscore": -0.1,   # just below the [0,1] tagger range
    "Ds_nearestjet_dR": 5.0,        # beyond any real dR in the sample
    "jet_leadpt": 0.0,
}

# Exactly one tight muon in the event (the best muon itself, no extra of either
# charge), plus a b-jet veto. The veto is a CUT, not a feature. With W+jets MC as
# background its job was removing a process the background model didn't contain at
# all (tt); now that SS data does contain tt, its stated purpose is what the
# 2026-08-10 re-screening showed directly: the loose/medium b-jet-count separation
# against SS survives the tight veto, so the veto and those features are not
# redundant. Tight WP keeps 99.3% of signal while removing 25% of data; dR04
# matches the AK4 jet cone radius. See BDT/BDT.py for the full reasoning.
B_VETO = "n_b_jets_tight_dR04"
SELECTION = ["num_os_mu_tight", "num_ss_mu_tight", B_VETO]
LOAD_BRANCHES = list(dict.fromkeys(FEATURE_NAMES + SELECTION))

DIRA_IDX = FEATURE_NAMES.index("dira")
LOG_IDX = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]

# ========= 训练配置 =========
# Both signal channels, each event at weight 1. The two productions were sized in
# proportion to their cross sections (500 vs ~992 PAT files, N_processed 1.54M vs
# 2.86M = 1:1.86), so N_gen is already proportional to sigma x BR and the
# equivalent-luminosity weight per event is the same for both -- concatenating them
# reproduces the physical DsW:DsstW mixture, and their selection efficiencies are
# near identical anyway (9.32% vs 9.49%), so it survives the selection too.
# If the two channels ever need genuinely unequal weights, put them here.
SIGNALS = {
    "HtoDsW":   1.0,
    "HtoDsstW": 1.0,
}

# WtoLNu2Jets is split into three jet-multiplicity bins of the SAME physical
# process, weight-merged using the full statistics of all three files. Weights are
# the equivalent-luminosity reweighting factors (external -- these files have no
# embedded Lumi_reweight branch).
WJETS_WEIGHTS = {
    "WtoLNu2Jets_0J": 2.3293047885,
    "WtoLNu2Jets_1J": 0.5382864853,
    "WtoLNu2Jets_2J": 0.2748710631,
}

# The two background configurations. `_f43` is a feature-set generation tag: a
# future feature change gets its own suffix rather than clobbering these.
#
# parity: SS data holds out the odd-event_id half so a later sculpting check on
#   Data_2022_SS is not circular (train-set self-validation). W+jets gets None on
#   purpose -- it has no such check to protect, and the original 40-feature v2 used
#   no parity split either. Adding one would change two things at once relative to
#   v2 (feature set AND training statistics), which is the exact defect this
#   refactor exists to remove. The head-to-head metric is OS sideband data, which
#   neither model ever sees.
# is_data: drives hep.cms.label -- real data gets "Preliminary"/data=True,
#   simulation gets "Simulation"/data=False. The pre-refactor v3 run stamped every
#   panel "Simulation" while training against real SS data; that is the bug this
#   field exists to prevent.
TRAININGS = {
    "ssdata": dict(
        samples={"Data_2022_SS": 1.0},
        ss=True,
        parity=0,
        weighted=False,          # real data: no per-process lumi reweighting
        label="Data SS",
        is_data=True,
        out_dir=f"{MVA_DIR}/dsw_dsstw_ssdata_f43",
    ),
    "wjets": dict(
        samples=WJETS_WEIGHTS,
        ss=False,
        parity=None,
        weighted=True,
        label="W+jets MC",
        is_data=False,
        out_dir=f"{MVA_DIR}/dsw_dsstw_wjets_f43",
    ),
}


def resolve(key):
    """Look up a background config, failing loudly on a typo."""
    if key not in TRAININGS:
        sys.exit(f"unknown background '{key}'; expected one of: {', '.join(TRAININGS)}")
    return TRAININGS[key]


# ========= 数据加载 =========
def load(name, *, truth_filter=False, ss=False, parity=None, extra_branches=()):
    """Load one Samples/ file, apply the analysis selection, return the feature matrix.

    ss:      select the same-sign muon configuration instead of opposite-sign.
    parity:  keep only rows with event_id % 2 == parity (None = keep all).
    Returns (X, extras) where extras is a dict of any extra_branches requested,
    already masked in step with X.
    """
    branches = list(LOAD_BRANCHES)
    if truth_filter:
        branches += ["Ds_match", "mu_match"]
    if parity is not None:
        branches += ["event_id"]
    branches += list(extra_branches)
    branches = list(dict.fromkeys(branches))

    with uproot.open(f"{SAMPLES_DIR}/{name}.root")[TREE] as t:
        a = t.arrays(branches, library="np")

    n_in = len(a[B_VETO])
    if ss:
        m = (a["num_ss_mu_tight"] == 1) & (a["num_os_mu_tight"] == 0)
    else:
        m = (a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
    m &= a[B_VETO] == 0
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if parity is not None:
        m &= (a["event_id"] % 2) == parity

    X = np.column_stack([a[f][m] for f in FEATURE_NAMES])
    print(f"  {name:18s} {n_in:>9,} -> {m.sum():>9,} rows after selection")
    return X, {b: a[b][m] for b in extra_branches}


# ========= 预处理 =========
def preprocess(X):
    """The one implementation of the chain: sentinel fix -> finite mask -> log.

    Returns (X, mask). X is modified in place for the sentinel step and the log
    transform is applied only to rows inside the mask (so log(1.01 - (-777)) is
    never evaluated). Callers differ in what they do with ~mask: DNN.py drops those
    rows, predict.py keeps them and writes score_dsmu = -777.

    Order matters. Sentinels must be regularized BEFORE the scaler ever sees them,
    and before the log, or the transform and the fitted mean/variance are garbage.
    """
    for feat, fill in SENTINEL_FIX.items():
        i = FEATURE_NAMES.index(feat)
        bad = X[:, i] < -1
        if bad.any():
            print(f"  sentinel: {feat:24s} {bad.sum():>9,} rows ({100*bad.mean():5.2f}%) -> {fill}")
            X[bad, i] = fill

    mask = np.isfinite(X).all(axis=1) & (X[:, DIRA_IDX] > -1)

    Xm = X[mask]
    Xm[:, LOG_IDX] = np.log1p(Xm[:, LOG_IDX])
    Xm[:, DIRA_IDX] = np.log(1.01 - Xm[:, DIRA_IDX])
    X[mask] = Xm
    return X, mask


# ========= 模型 =========
class phipimuNet(nn.Module):
    def __init__(self, n_in, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_in, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 64),  nn.BatchNorm1d(64),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(64, 32),   nn.BatchNorm1d(32),  nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, 16),   nn.BatchNorm1d(16),  nn.ReLU(),
            nn.Linear(16, 1),
        )

    def forward(self, x):
        return self.net(x)


# ========= 防漂移守卫 =========
FEATURE_FILE = "feature_names.txt"


def save_feature_names(out_dir):
    """Stamp the model dir with the exact feature list (and order) it was trained on."""
    with open(f"{out_dir}/{FEATURE_FILE}", "w") as f:
        f.write("\n".join(FEATURE_NAMES) + "\n")


def check_model_dir(model_dir):
    """Refuse to score with a model whose feature set does not match this module.

    Catches the two silent failure modes: a model dir from before the feature list
    was stamped (the frozen 40-feature v2), and a list that was reordered without
    retraining -- the latter raises no exception on its own because feature order is
    only implicitly bound to the scaler and the first Linear layer.
    """
    path = f"{model_dir}/{FEATURE_FILE}"
    if not os.path.exists(path):
        sys.exit(
            f"{path} not found.\n"
            f"This model dir predates feature_names.txt -- it is a frozen older build "
            f"(e.g. the 40-feature v2 at dsw_dsstw_wjets). Do NOT re-score it with the "
            f"current {len(FEATURE_NAMES)}-feature code: the scores would be silently wrong. "
            f"Retrain into a new _f43 dir instead."
        )

    with open(path) as f:
        stored = [ln.strip() for ln in f if ln.strip()]
    if stored != FEATURE_NAMES:
        n = min(len(stored), len(FEATURE_NAMES))
        first = next((i for i in range(n) if stored[i] != FEATURE_NAMES[i]), n)
        sys.exit(
            f"feature mismatch vs {path}\n"
            f"  stored: {len(stored)} features, common.py: {len(FEATURE_NAMES)}\n"
            f"  first difference at index {first}: "
            f"stored={stored[first] if first < len(stored) else '<end>'} "
            f"vs common.py={FEATURE_NAMES[first] if first < len(FEATURE_NAMES) else '<end>'}"
        )

    mean = np.load(f"{model_dir}/scaler_mean.npy")
    scale = np.load(f"{model_dir}/scaler_scale.npy")
    if mean.size != len(FEATURE_NAMES) or scale.size != len(FEATURE_NAMES):
        sys.exit(f"scaler size {mean.size}/{scale.size} != {len(FEATURE_NAMES)} features "
                 f"in {model_dir} -- the model was trained on a different feature set.")

    state = torch.load(f"{model_dir}/phipimu_model.pt", map_location="cpu")
    in_features = state["net.0.weight"].shape[1]
    if in_features != len(FEATURE_NAMES):
        sys.exit(f"model input layer expects {in_features} features but common.py "
                 f"declares {len(FEATURE_NAMES)} -- {model_dir} is a different feature set.")
    return state
