"""Single source of truth for the phipimu DNN: feature set, preprocessing chain,
network definition, and the background-configuration table.

Why this module exists: FEATURE_NAMES used to be copy-pasted into DNN.py and
predict.py, kept in sync only by a comment saying "must stay identical". They
drifted anyway -- the 40-feature v2 (W+jets background) and the 43-feature v3 (SS
data background) were compared head-to-head on 2026-08-10 as if the only
difference were the training background, which invalidated the comparison. Import
from here instead; do not re-declare any of it downstream.

BDT/ imports from here as of 2026-08-14. It used to be deliberately left out, on
the grounds that "importing would silently change the BDT the next time a feature is
added or reordered" -- a real hazard, but the fix for it is the stamp guard below,
not staying disconnected. check_feature_stamp() makes a feature edit a loud failure
at scoring time for the BDT exactly as check_model_dir() already did for the DNN, so
the objection is neutralized rather than overridden. Before this, BDT was pinned to a
different 42-feature set (its own copy, plus a second hand-typed copy in
BDT/predict.py), which is precisely the drift this module exists to prevent, and it
made "BDT vs DNN" unmeasurable: the two differed in feature set, background and
signal all at once.

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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eras import (B_VETO_BRANCH, ERAS, LEPTON_BRANCHES,  # noqa: F401  (re-exported)
                  SIGNAL_ERA_WEIGHT, TAG_TO_ERA)

MVA_ROOT = "/eos/home-k/kakang/H2DW/Analysis/MVA"
# Moved out of MC/2022/ on 2026-08-24: this directory holds real data as well as MC,
# so living under MC/ was a misnomer, and it is an MVA input rather than a tuplizer
# output. analysis/Samples/{run.sh,condor_submit.py} write here -- keep all three in step.
SAMPLES_DIR = f"{MVA_ROOT}/Samples"
MVA_DIR = f"{MVA_ROOT}/DNN"
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
# Defined in analysis/eras.py, which imports nothing heavy, so the pure-uproot
# producers (highest_pt_dedup.py, make_samples.py) can share these names instead of
# hand-copying them -- this module pulls in torch, which they must not.
B_VETO = B_VETO_BRANCH
SELECTION = [*LEPTON_BRANCHES, B_VETO]
LOAD_BRANCHES = list(dict.fromkeys(FEATURE_NAMES + SELECTION))

DIRA_IDX = FEATURE_NAMES.index("dira")
LOG_IDX = [FEATURE_NAMES.index(f) for f in LOG_FEATURES]

# ========= 训练配置 =========
# The one seed the canonical models were trained with. Lives here rather than in
# DNN.py because out_dir_for() needs it to decide when to add a _seed suffix.
SEED_DEFAULT = 42

# Both signal channels, each event at weight 1. The two productions were sized in
# proportion to their cross sections (500 vs ~992 PAT files, N_processed 1.60M vs
# 2.98M = 1:1.86 after the 2026-08-17 reproduction), so N_gen is already
# proportional to sigma x BR and the
# equivalent-luminosity weight per event is the same for both -- concatenating them
# reproduces the physical DsW:DsstW mixture, and their selection efficiencies are
# near identical anyway (9.32% vs 9.49%), so it survives the selection too.
# If the two channels ever need genuinely unequal weights, put them here.
# The _2022 suffix is the era tag, added 2026-08-24 to match the Data_2022 /
# Data_2022EE convention.
#
# postEE signal MC arrived 2026-09-03 and is ADDED here, not substituted. The two eras
# were generated with the same event count (DsW 1.6M, DsstW 3M each), so N_gen carries
# no era information and the weight is purely the luminosity ratio:
#
#     w_postEE / w_preEE = 26.68 / 7.99 = 3.3392
#
# preEE stays at 1.0 as the anchor -- signal is not on an absolute scale anyway (the
# backgrounds carry sigma*L/sum(w) while signal is per-event 1.0), so the anchor choice
# is physically arbitrary and this one leaves every preEE-based result undisturbed.
# The channel-relative factor stays 1.0: N_gen/BR is 1.000e11 for both DsW (1.6M /
# 1.6e-5) and DsstW (3M / 3.0e-5), so N_gen is proportional to sigma x BR by construction
# and the two channels are already correctly mixed at equal weight.
#
# This closes the era bias documented on TRAININGS below: weighted, the signal is now
# 23.0% preEE / 77.0% postEE, matching the background's 77% postEE.
#
# NOTE the weight lives ONLY here -- the four Samples/ files are indistinguishable on
# disk (no Lumi_reweight branch, by the --no-lumi-weight convention), so any macro that
# reads them without going through SIGNALS sees both eras at equal weight, with no error.
# 2026-09-04: era is a BRANCH now, not a filename suffix, and the era weight is a
# Lumi_reweight COLUMN in the file rather than a value in this dict. So SIGNALS is a
# plain list of the two channels -- one file each, both eras inside.
#
# This is what closes the hazard the previous comment ended on: the weight used to
# live only here, so a macro reading the files directly saw both eras at equal weight
# with no error. It cannot now -- the weight travels with the rows.
#
# Read the weight with common.load(..., extra_branches=("Lumi_reweight",)). For signal
# it is RELATIVE (era 1 = 1.0, era 2 = L2/L1), not an absolute yield; see eras.py.
SIGNALS = ["HtoDsW", "HtoDsstW"]

# Everything in Samples/, i.e. everything predict.py scores. ONE list, imported by both
# DNN/predict.py and BDT/predict.py -- they used to keep separate literal copies that
# BDT's comment claimed were "identical to the list in DNN/predict.py" while being 13
# samples short, which is why the BDT could not reproduce the DNN's signal definition.
# make_samples.py asserts its table matches this, so the producer and the consumer
# cannot drift either.
#
# The four axes that are still FILENAMES rather than branches, and must stay that way:
# OS vs SS (different candidate reconstructions, 12.9% event overlap), incl vs the jet
# bins (alternative descriptions of one process, never components), GenTag (the same
# events as the bins, different branch set), and the process itself.
ALL_SAMPLES = [
    "HtoDsW", "HtoDsstW",
    "Data", "Data_SS",
    "WtoLNu2Jets", "WtoLNu2Jets_SS",
    "WtoLNu2Jets_incl", "WtoLNu2Jets_incl_SS",
    "WtoLNu2Jets_GenTag",
]

# WtoLNu2Jets is split into three jet-multiplicity bins of the SAME physical process.
# The bins are a `jet_bin` BRANCH since the 2026-09-04 merge, so this is the list of bin
# ids -- not a weight table.
#
# The weight table that used to live here (WJETS_WEIGHTS: three sigma*L/N constants, one
# per preEE bin, applied uniformly) was deleted with the per-bin files on 2026-09-04. It
# was wrong twice over: (a) the samples carry a per-event gen_weight (amcatnloFXFX, two
# valued +-|w|) and the negative fraction varies 3.4x across bins, so no single constant
# per bin can reproduce the post-selection yield -- old and new agree only on the TOTAL,
# diverging by (1-2f)/(1-2f_sel) once anything is cut; and (b) it only ever had preEE
# keys, so every macro that iterated it was silently preEE-only even after the postEE
# jet bins landed in Samples/ on 2026-08-25.
#
# The weight now travels with the rows: read Lumi_reweight from the file.
WJETS_JET_BINS = (0, 1, 2)

# The two background configurations. `_f43` is a feature-set generation tag: a
# future feature change gets its own suffix rather than clobbering these.
#
# out_name vs out_dir: `out_name` is the bare directory name, which out_dir_for()
#   composes under any tool (DNN/, BDT/) and optional seed suffix. `out_dir` stays
#   as the literal DNN path it has always been -- seven macros read cfg["out_dir"]
#   directly (check_data_vs_wjets, check_feature_tail, check_jetbin_closure,
#   check_os_minus_ss, compare_models, draw_score_os_ss, draw_lepton_bjet_vs_cut),
#   so it is kept, not renamed. The two must agree; the assertion below enforces it.
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
# Feature groups used by the pileup control models below. Named here rather than
# inline in TRAININGS so the era-shift macro and the training config cannot disagree
# about what "the isolation features" means.
#
# ISO_ACTIVITY: everything that measures how busy the event is around the candidate.
# This is where pileup enters, and it holds the three largest data-vs-MC era residuals
# (Ds_IsoR04_sumPUPt 0.276, sumNeutralHadronEt 0.180, n_b_jets_loose_dR04 0.141;
# see macros/draw_era_feature_shift.py).
#
# VTX_IP: displacement observables. They also carry a measurable era residual
# (Ds_prob 0.082, the *_ipsig 0.056-0.078) because pileup degrades the primary vertex,
# but they are the analysis's core discriminants, so they are neutralised only in the
# second, more aggressive control.
#
# MET is deliberately in NEITHER: its residual is 0.009-0.029, and it carries real
# W -> mu nu information, so removing it would weaken the classifier for no pileup reason.
ISO_ACTIVITY = ("Ds_IsoR04_sumChargedHadronPt", "Ds_IsoR04_sumNeutralHadronEt",
                "Ds_IsoR04_sumPhotonEt", "Ds_IsoR04_sumPUPt", "Ds_IsoR04_RelIso",
                "Ds_IsoR03_RelIso", "mu_IsoR04_RelIso", "num_jets", "num_jets_pt30",
                "jet_ht", "jet_leadpt", "n_b_jets_loose_dR04", "n_b_jets_medium_dR04",
                "Ds_nearestjet_dR", "Ds_nearestjet_bscore")
VTX_IP = ("Ds_prob", "Ds_ipsig", "Ds_ipxysig", "Kp_ipsig", "Km_ipsig", "pi_ipsig",
          "mu_ipsig", "Kp_ipxysig", "Km_ipxysig", "pi_ipxysig")


TRAININGS = {
    "ssdata": dict(
        # Both eras of SS data, at weight 1.0 each: data carries no luminosity weight, and
        # more SS statistics is exactly what the 70%-working-point study identified as the
        # bottleneck (2026-08-14: ~94 hard events of training equivalent at eps_sig=70%).
        #
        # The era-imbalance bias noted here on 2026-08-21 is GONE, and this comment used to
        # claim otherwise long after it stopped being true -- it said the signal was 100%
        # preEE, which was correct only until the postEE signal MC was produced (2026-08-25
        # .. 09-04) and merged into Samples/. Measured 2026-09-08, the era composition of the
        # two training classes is identical to the decimal:
        #     signal (weighted by Lumi_reweight, which is what training uses)  23.2% / 76.8%
        #     background Data_SS at parity=0                                   23.2% / 76.8%
        # so era carries no class information and there is nothing for the network to gain
        # by learning it. DNN.py reads the signal's Lumi_reweight unconditionally; the
        # `weighted` flag below governs only the background branch.
        #
        # THE BIAS THAT IS REAL is a different one, and it does not need any era imbalance:
        # signal is MC and background is real data, so every data/MC mismodelling is a free
        # discriminant. Pileup is the measured case -- <Ds_IsoR04_sumPUPt> is 9.65 GeV in
        # signal MC versus 11.99 (preEE) / 14.03 (postEE) GeV in data, and no pileup
        # reweighting exists anywhere in this analysis. Because the MC is era-flat while the
        # data is not, the gap is 4% in preEE and 22% in postEE, which is why the DNN cut
        # removes ~25-38% more postEE data than preEE data. See the memory note
        # project_era_pileup_mismodelling and analysis/macros/check_era_split.py.
        #
        # One file, both eras inside (era branch). Data carries no luminosity weight,
        # so the two eras really are at weight 1.0 each -- nothing is lost by merging.
        samples={"Data_SS": 1.0},
        ss=True,
        parity=0,
        weighted=False,          # real data: no per-process lumi reweighting
        label="Data SS",
        is_data=True,
        out_name="dsw_dsstw_ssdata_f43",
        out_dir=f"{MVA_DIR}/dsw_dsstw_ssdata_f43",
    ),
    "wjets": dict(
        # One file, both eras and all three jet bins inside. `weighted=True` means the
        # weight is read per event from the file's Lumi_reweight rather than taken from
        # this dict -- which is what replaces the deleted WJETS_WEIGHTS constants, and
        # is the fix the 2026-08-21 note asked for before these models are retrained.
        # The value here is unused when weighted=True; it stays 1.0 for shape.
        samples={"WtoLNu2Jets": 1.0},
        ss=False,
        parity=None,
        weighted=True,
        label="W+jets MC",
        is_data=False,
        out_name="dsw_dsstw_wjets_f43",
        out_dir=f"{MVA_DIR}/dsw_dsstw_wjets_f43",
    ),
}

# Pileup control models. Same samples and same 43 inputs as "ssdata"; the listed
# features are blanked to a constant so the network cannot use them. They exist to
# answer one question: if the classifier cannot see the mismodelled pileup information,
# does the postEE/preEE data asymmetry go away?
#
# Two of them, because a null result from the first alone would be ambiguous -- the
# vertex/IP features carry a smaller but real era residual and could route the same
# information through. Reading agreed in advance:
#   noiso alone fixes it        -> the mechanism is the isolation/jet group
#   only noiso_novtx fixes it   -> it travels mainly through the vertex/IP observables
#   neither fixes it            -> pileup is not the main cause and the diagnosis is wrong
for _name, _drop in (("ssdata_noiso", ISO_ACTIVITY),
                     ("ssdata_noiso_novtx", ISO_ACTIVITY + VTX_IP)):
    TRAININGS[_name] = dict(TRAININGS["ssdata"], neutralize=_drop,
                            label=f"Data SS [{_name.replace('ssdata_', '')}]",
                            out_name=f"dsw_dsstw_{_name}_f43",
                            out_dir=f"{MVA_DIR}/dsw_dsstw_{_name}_f43")

for _k, _c in TRAININGS.items():
    assert _c["out_dir"] == f"{MVA_DIR}/{_c['out_name']}", f"{_k}: out_dir/out_name disagree"


def resolve(key):
    """Look up a background config, failing loudly on a typo."""
    if key not in TRAININGS:
        sys.exit(f"unknown background '{key}'; expected one of: {', '.join(TRAININGS)}")
    return TRAININGS[key]


TOOLS = ("DNN", "BDT")


def out_dir_for(key, tool="DNN", seed=None):
    """Model directory for one (background, tool, seed) combination.

    seed=None or the default 42 gives the canonical path -- i.e. the existing
    dsw_dsstw_*_f43 dirs, which the seed-spread study must not touch. Any other
    seed gets a _seed{N} suffix so several runs of the same config coexist.
    """
    if tool not in TOOLS:
        sys.exit(f"unknown tool '{tool}'; expected one of: {', '.join(TOOLS)}")
    name = resolve(key)["out_name"]
    if seed is not None and seed != SEED_DEFAULT:
        name = f"{name}_seed{seed}"
    return f"{MVA_ROOT}/{tool}/{name}"


# Figures live on EOS next to the model that produced them, never in the source
# tree (they used to all land in one flat Tuples/HtoDsW/figures dir on AFS).
# Categories are the subdirectory names; keep this list short and add to it only
# when a genuinely new study appears, not per script.
FIG_CATEGORIES = (
    "training",       # DNN.py / BDT.py own diagnostics
    "closure",        # data vs MC closure and yield scans
    "os_ss",          # OS-vs-SS charge-configuration studies
    "misid",          # K->pi reflection / misid veto
    "era",            # preEE vs postEE splits
    "gentag",         # GenTag truth-composition plots
    "features",       # input-feature distributions and screening
    "preselection",   # pre-MVA selection studies
    "comparisons",    # anything overlaying several models
)

# Model-free figures (pre-MVA, feature screening, multi-model overlays) have no
# single model to belong to, so they go here rather than into an arbitrary one.
SHARED_FIG_DIR = f"{MVA_ROOT}/figures"


def fig_path(category, name, *, key=None, tool="DNN", seed=None):
    """Absolute path for one figure, creating its directory.

    key=None      -> SHARED_FIG_DIR/<category>/, for figures that do not belong
                     to a single model (pre-MVA studies, multi-model overlays).
    key set       -> out_dir_for(key, tool, seed)/figures/<category>/, so a new
                     training run's figures follow it into its own directory
                     instead of piling up under whichever model was hardcoded.

    Pass tool/seed through from whatever the caller actually scored with; they
    default to the canonical DNN/seed-42 model the same way out_dir_for does.
    """
    if category not in FIG_CATEGORIES:
        sys.exit(f"unknown figure category '{category}'; "
                 f"expected one of: {', '.join(FIG_CATEGORIES)}")
    base = SHARED_FIG_DIR if key is None else f"{out_dir_for(key, tool, seed)}/figures"
    d = f"{base}/{category}"
    os.makedirs(d, exist_ok=True)
    return f"{d}/{name}"


# ========= 数据加载 =========
def load(name, *, era=None, jet_bin=None, truth_filter=False, ss=False, parity=None,
         extra_branches=()):
    """Load one Samples/ file, apply the analysis selection, return the feature matrix.

    era:     keep only this era id (see eras.ERAS). None = every era in the file.
    jet_bin: keep only this jet bin. None = every bin. Errors if the file has no
             jet_bin branch, rather than silently returning everything.
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
    if era is not None:
        branches += ["era"]
    if jet_bin is not None:
        branches += ["jet_bin"]
    branches += list(extra_branches)
    branches = list(dict.fromkeys(branches))

    with uproot.open(f"{SAMPLES_DIR}/{name}.root")[TREE] as t:
        if jet_bin is not None and "jet_bin" not in t.keys():
            sys.exit(f"load({name!r}, jet_bin={jet_bin}): that file has no jet_bin branch. "
                     f"Only the jet-binned W+jets samples do -- incl and GenTag-merged "
                     f"files are separate productions, not a fourth bin.")
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
    if era is not None:
        m &= a["era"] == era
    if jet_bin is not None:
        m &= a["jet_bin"] == jet_bin

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
def neutralize(X, feats):
    """Blank the named feature columns in place, and return the names actually blanked.

    The columns are set to a constant (their own mean), which carries exactly zero
    information, rather than being DELETED. That is the whole point: the feature count,
    the scaler size, the input-layer width and the feature_names.txt stamp all stay at
    43, so every drift guard and every downstream consumer keeps working untouched.
    Deleting columns instead would mean making FEATURE_NAMES per-training, which has
    20+ call sites across DNN.py, predict.py, BDT and six macros -- far too much blast
    radius for a control experiment.

    The constant is a literal 0.0, NOT the column mean. That distinction is the whole
    correctness argument: predict.py must feed the network exactly what training fed it,
    and a column mean is a property of whichever array is in hand -- the training set
    for DNN.py, but each scored sample separately for predict.py. Using the mean would
    hand the network a different constant at scoring time, and a different one per
    sample, injecting a spurious sample-dependent offset through the first-layer weight
    that training never removed (a constant input is absorbed into the bias, so that
    weight keeps its initialisation). Zero is the same number everywhere.

    A constant is used rather than a random permutation for the same reason: predict.py
    has to reproduce it exactly, and a permutation is not reproducible.

    StandardScaler maps an all-zero column to mean_ = 0, scale_ = 1 (it substitutes 1
    for a zero variance rather than dividing by zero), so the scaled input is exactly
    0.0 in both training and scoring.
    """
    if not feats:
        return ()
    idx = [FEATURE_NAMES.index(f) for f in feats]      # loud on a typo
    X[:, idx] = 0.0
    return tuple(feats)


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


def check_feature_stamp(model_dir):
    """Refuse to score with a model whose feature stamp does not match this module.

    Catches the two silent failure modes shared by every algorithm: a model dir from
    before the feature list was stamped (the frozen 40-feature v2, or the frozen
    42-feature BDT at dsw_wjets_test), and a list that was reordered without
    retraining -- the latter raises no exception on its own because feature order is
    only implicitly bound to the scaler, the first Linear layer, and the tree splits.

    Split out of check_model_dir so BDT/predict.py can reuse it without importing
    torch machinery it has no use for.
    """
    path = f"{model_dir}/{FEATURE_FILE}"
    if not os.path.exists(path):
        sys.exit(
            f"{path} not found.\n"
            f"This model dir predates feature_names.txt -- it is a frozen older build "
            f"(e.g. the 40-feature v2 at dsw_dsstw_wjets, or the 42-feature BDT at "
            f"dsw_wjets_test). Do NOT re-score it with the current "
            f"{len(FEATURE_NAMES)}-feature code: the scores would be silently wrong. "
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


def check_n_features(model_dir, n_in, what):
    """The algorithm-specific half: the trained object's own input width."""
    if n_in != len(FEATURE_NAMES):
        sys.exit(f"{what} expects {n_in} features but common.py declares "
                 f"{len(FEATURE_NAMES)} -- {model_dir} is a different feature set.")


def check_model_dir(model_dir):
    """check_feature_stamp plus the DNN-specific scaler and input-layer checks.

    Returns the loaded state dict, so the file that was checked is the file loaded.
    """
    check_feature_stamp(model_dir)

    mean = np.load(f"{model_dir}/scaler_mean.npy")
    scale = np.load(f"{model_dir}/scaler_scale.npy")
    if mean.size != len(FEATURE_NAMES) or scale.size != len(FEATURE_NAMES):
        sys.exit(f"scaler size {mean.size}/{scale.size} != {len(FEATURE_NAMES)} features "
                 f"in {model_dir} -- the model was trained on a different feature set.")

    state = torch.load(f"{model_dir}/phipimu_model.pt", map_location="cpu")
    check_n_features(model_dir, state["net.0.weight"].shape[1], "model input layer")
    return state
