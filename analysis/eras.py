"""Era table and the selection branch names -- the one place both halves of the
analysis agree on what "2022EE" means.

This module deliberately imports NOTHING but the standard library. It exists because
the two natural homes for these constants cannot both reach each other:

  * analysis/DNN/common.py imports torch at module level, so a pure-uproot macro
    cannot import it just to learn a branch name. That is why
    analysis/macros/highest_pt_dedup.py used to carry a hand-copied B_VETO_BRANCH
    with a "if it changes there, change it here too" comment -- exactly the kind of
    duplicated constant that goes stale silently.
  * The macros each grew their own copy of the era table. Before this module,
    "7.99"/"26.68" existed only as display strings in three macros, and the
    matching 7990/26680 only in a common.py comment and an argparse help string.
    Nothing tied the two representations together.

So: this is the upstream of both. common.py re-exports what it needs; the macros
import from here directly.

ERA IDs ARE A WIRE FORMAT. They are written into the `era` branch of every file in
Samples/, so a value here can never be reused for a different era. Append only.
"""

# id -> era. `tag` is the string that appears in the old per-era filenames and in
# the production paths (/eos/.../MC/<tag>/...); `label` is for plot titles;
# `lumi_pb` is the integrated luminosity in /pb, which is the unit --lumi wants.
ERAS = {
    1: dict(tag="2022",   label="2022 (preEE)",    lumi_pb=7990.0,  lumi_fb="7.99"),
    2: dict(tag="2022EE", label="2022EE (postEE)", lumi_pb=26680.0, lumi_fb="26.68"),
    # 3: 2023, 4: 2023BPix, 5: 2024 -- append here, never renumber.
}

TAG_TO_ERA = {v["tag"]: k for k, v in ERAS.items()}

# Measured data pileup profile per era, from pileupCalc.py over that era's certified
# lumi sections. Histogram name inside each file is "pileup".
#
# This is the NUMERATOR of the pileup reweighting; the denominator is the nTrueInt_all
# histogram the tuplizer accumulates over every generated event. Both therefore describe
# THE BEAM, not any selection -- see make_samples.py's pu_weight_lookup().
#
# The 2022 profiles are 99 bins over [0, 99] while the tuplizer's are 100 over [0, 100].
# The edges coincide on the overlap and the MC is empty above 99, but that is asserted
# rather than assumed, because a silent one-bin offset would misweight every event.
PU_PROFILE_DIR = "/eos/home-k/kakang/H2DW/Analysis/MVA/Samples/pileup"
PU_PROFILE_HIST = "pileup"
PU_PROFILE = {
    1: f"{PU_PROFILE_DIR}/PU_2022.root",
    2: f"{PU_PROFILE_DIR}/PU_2022EE.root",
    # 3: PU_2023.root, 4: PU_2023BPix.root, 5: PU_2024.root -- present on disk already.
}

# The relative per-era weight for the signal samples. The two signal productions were
# generated with the SAME event count in every era (DsW 1.6M, DsstW 3M each), so N_gen
# carries no era information and the weight is purely the luminosity ratio, anchored on
# era 1. Derived rather than written out, so it cannot drift from the lumi values above:
# it used to be the literal 3.3392 in common.py's SIGNALS, which is this number rounded.
#
# RELATIVE, NOT ABSOLUTE. Signal is not on an absolute yield scale in this analysis --
# the backgrounds carry sigma*L/sum(w) while signal is anchored at 1.0 -- so summing
# these weights does NOT give an expected event count. See make_samples.py's Meta tree,
# which records weight_mode per file so this distinction is machine-readable.
SIGNAL_ERA_WEIGHT = {e: ERAS[e]["lumi_pb"] / ERAS[1]["lumi_pb"] for e in ERAS}

# Branch names the row selection is written in terms of. common.py derives its
# B_VETO / SELECTION from these; highest_pt_dedup.py and make_samples.py use them
# directly. One definition, two consumers.
B_VETO_BRANCH = "n_b_jets_tight_dR04"
LEPTON_BRANCHES = ("num_os_mu_tight", "num_ss_mu_tight")

# The Ds mass window that --blind removes from OS data.
M_DS = 1.96835
SR_HALF_WIDTH = 0.04
