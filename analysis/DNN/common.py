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
from eras import B_VETO_BRANCH, ERAS, LEPTON_BRANCHES, TAG_TO_ERA  # noqa: F401  (re-exported)

# MVA_ROOT still holds model output (DNN/, BDT/) and figures -- only the sample source
# moved. NanoMVA/Samples/ is the custom-NanoAOD + POG-corrections production
# (make_samples_nano.py), which replaced the private-tuplizer MVA/Samples/ production
# as of the 2026-09-17 migration (project_nano_merged_samples memory). The old
# MVA/Samples/*.root files are left on disk (not deleted) for reference, but nothing
# in this module reads them any more.
MVA_ROOT = "/eos/home-k/kakang/H2DW/Analysis/MVA"
SAMPLES_DIR = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Samples"
MVA_DIR = f"{MVA_ROOT}/DNN"
TREE = "Events"

# ========= 特征定义 =========
# Feature set excludes every decay-kinematics variable: those reconstruct
# phi_mass / Ds_mass to R^2 = 0.98 / 0.94, which would sculpt the resonance peaks
# used for the final fit. Dsmu_mass is a feature here, not a fit variable -- it has
# no resonant peak, and keeping it out of the fit is what lets Ds_pt stay in the MVA.
#
# NanoMVA migration (2026-09-17): every feature that has a POG-corrected counterpart
# in NanoMVA/Samples/ now uses it (met_pt_corr, W_tmass_corr, num_jets_pt30_corr,
# Dsmu_pt_corr, Dsmu_eta_corr, Dsmu_mass_corr, mu_pt_corr, jet_leadpt_corr), per the
# 0916 "default to corrected quantities" convention (feedback_corrected_by_default
# memory). jet_ht and num_jets (the two no-pt-threshold sums) are DROPPED outright,
# not given a same-threshold corrected replacement -- there isn't one: the custom
# NanoAOD Jet table only stores jets with pt > 15 GeV (jme_corrections.py docstring),
# so 9.6% of the old jet_ht (the sub-15-GeV tail) is permanently unrecoverable from
# any stored branch, corrected or not, and num_jets is the same problem at threshold
# zero. Swapping in a different-threshold corrected variable (jet_ht_pt20_corr,
# num_jets_pt30_corr is already a separate feature below) would silently change what
# is being tested, so it was rejected rather than substituted. 43 features -> 41.
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
    # mu_ipxysig was left out on 2026-08-04 as "AUC 0.504, i.e. nothing". ADDED BACK
    # 2026-09-21: that number does not survive the NanoMVA migration. Re-measured on
    # Signal_merged vs Data_SS_merged it is 0.566 (mu_ipsig 0.575, Kp_ipxysig 0.692,
    # pi_ipxysig 0.744) -- weak, but not nothing. Caveat on the old number: it is dated
    # 2026-08-04 while the background switch to SS data happened 08-10, so 0.504 was
    # very likely measured against W+jets MC and is not the same comparison.
    #
    # The reason to add it does NOT rest on that AUC, which is confounded (signal is MC,
    # background is real data, so part of any separation is instrumental IP-reconstruction
    # difference -- the same confound that invalidated the data/MC TVD study, see
    # project_mva_feature_set_revision memory). It rests on the redundancy measurement,
    # which is taken WITHIN one sample and carries no such confound: of the four
    # 2D/3D IP pairs, mu is the LEAST redundant -- Spearman r 0.669, vs 0.829/0.829/0.857
    # for Kp/Km/pi (same ordering under Pearson raw and Pearson log1p). So the original
    # "2D carries information complementary to 3D" argument, which is what put the three
    # track 2D features in, applies to the muon MORE strongly than to the tracks, not less.
    # (The r = 0.14/0.19/0.27 quoted above for the tracks could NOT be reproduced under
    # any of those three conventions -- measured 0.78-0.87 -- so do not cite it.)
    #
    # No sentinel and no log-domain problem: min = 0, zero negative, zero non-finite in
    # both Signal_merged and Data_SS_merged, so it needs a LOG_FEATURES entry (below)
    # like every other IP significance and nothing in SENTINEL_FIX.
    "Kp_ipxysig",
    "Km_ipxysig",
    "pi_ipxysig",
    "mu_ipxysig",
    # isolation -- the components, not just RelIso. sumChargedHadronPt is the single
    # strongest variable against every background (AUC 0.91 vs TT): a Ds buried in a
    # b jet is surrounded by charged hadrons. RelIso divides by Ds_pt, which dilutes it.
    # No corrected counterpart exists for any of these (isolation sums are not JES/JER
    # quantities), so they are unchanged from the pre-migration feature set.
    #
    # 2026-09-22: ALL of Ds_IsoR03_*, Ds_IsoR04_* and mu_IsoR04_* -- 18 branches, 11 of
    # them new. What was here before was the four R04 Ds components plus three RelIso, an
    # asymmetry that was inherited rather than argued: the paragraph above says "the
    # components, not just RelIso", yet R03 kept ONLY RelIso, the very form it calls
    # diluted. mu_IsoR03_* is deliberately NOT included (not requested).
    #
    # R03 is not redundant with R04. R03 is the INNER cone, so R04 = R03 + the 0.3-0.4
    # annulus, and measured on signal vs SS data: Ds_IsoR03_sumChargedHadronPt alone is
    # AUC 0.825 (only R04_RelIso 0.844 and R04_sumChargedHadronPt 0.837 beat it, and it
    # was not in the set); the annulus is AUC 0.606 and correlates only r = 0.164 with
    # R03 inside signal, so core and annulus are near-independent and handing the model
    # only their sum threw the decomposition away; the existing Ds_IsoR03_RelIso is NOT a
    # proxy for it (r = 0.332). r(R03, R04) also differs by class -- 0.782 in signal vs
    # 0.942 in SS data -- so "do the two cones track each other" discriminates in its own
    # right, and is only constructible when both cones are present.
    #
    # KNOWN REDUNDANCY, KEPT DELIBERATELY (user decision 2026-09-22). Verified to float32
    # precision on all three families:
    #     PFIso  == sumCH + max(0, sumNH + sumPh - 0.5*sumPU)     (max rel dev 0, bitwise)
    #     RelIso == PFIso / pt
    # so the three PFIso and the two Ds RelIso are determined by the twelve sums plus
    # Ds_pt. They stay because a NET CANNOT EASILY LEARN A RATIO -- handing it the
    # combination costs one input and saves it from reconstructing a division out of
    # log1p'd sums and a linear pt. (For a BDT the same columns would be dead weight.)
    # mu_IsoR04_RelIso additionally carries something no other feature does: the
    # UNCORRECTED mu_pt. That is not a convention slip -- the POG iso WP and its SF are
    # defined on the uncorrected pt, which is why apply_corrections.py leaves RelIso as
    # stored. mu_IsoR04_PFIso is the one near-miss: rebuilding it from its four sums is
    # off by 7.4e-03, so muon isolation uses the official pfIsolationR04() rather than a
    # plain dBeta sum.
    #
    # Numerics: all 11 new columns have min = 0, zero negative and zero non-finite across
    # all four merged samples (9,452,569 rows after selection), so each one goes into
    # LOG_FEATURES and none needs a SENTINEL_FIX entry.
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
    "Ds_IsoR04_RelIso",
    "Ds_IsoR04_PFIso",
    "Ds_IsoR03_sumChargedHadronPt",
    "Ds_IsoR03_sumNeutralHadronEt",
    "Ds_IsoR03_sumPhotonEt",
    "Ds_IsoR03_sumPUPt",
    "Ds_IsoR03_RelIso",
    "Ds_IsoR03_PFIso",
    "mu_IsoR04_sumChargedHadronPt",
    "mu_IsoR04_sumNeutralHadronEt",
    "mu_IsoR04_sumPhotonEt",
    "mu_IsoR04_sumPUPt",
    "mu_IsoR04_RelIso",
    "mu_IsoR04_PFIso",
    # MET / event. num_jets_pt20..pt60 are near-duplicates (AUC all ~0.80), so only
    # one threshold is kept -- and that threshold's corrected counterpart is what is
    # used (num_jets_pt30_corr). num_jets (no threshold) is dropped, see the module
    # comment above.
    #
    # met_sig DROPPED 2026-09-21, same category as jet_ht: a JES-dependent quantity for
    # which no corrected version can be built. metSignificance() is MET^T V^-1 MET with V
    # assembled by PAT's METSignificance from jet pt x resolution, so the stored value is
    # computed from pre-JEC inputs; and the Puppi V is not in the custom NanoAOD (only the
    # PF MET covariance, MET_covXX/XY/YY, which is a ~10x different quantity), so it cannot
    # be recomputed. Concrete cost of keeping it: apply_corrections.py does not recompute it
    # for _jesUp/_jesDown/_jerUp/_jerDown, so it stays frozen under JME variations and the
    # JES systematic on the score would be UNDERESTIMATED.
    #
    # met_sumEt KEPT, and the distinction is the whole point: it is the scalar sum of
    # candidate ET. Type-1 propagation moves the MET VECTOR, not the scalar sum, and JEC
    # calibrates jet-level objects while this sums candidates -- so it is JES-independent
    # by construction and freezing it under a JME variation is CORRECT, not a deficiency.
    # "Has no correction" means "correction does not apply" here, exactly as it does for
    # the isolation and IP features, and is not grounds for dropping it. (There is a
    # separate argument for dropping it -- a whole-event activity sum is the most PU/UE-
    # tune-coupled input in a signal=MC vs background=SS-data training -- but that one has
    # never been measured and no clean data handle exists for it without a tt MC sample.
    # User decision 2026-09-21: keep.)
    "met_pt_corr",
    "met_sumEt",
    "W_tmass_corr",
    "num_jets_pt30_corr",
    # Ds + mu kinematics -- hard from a 125 GeV H vs soft inside a jet. Ds_pt/mu_eta
    # have no corrected counterpart (the muon correction only rescales pt, not
    # direction; Ds_pt is a decay-vertex-fit quantity JME does not touch).
    "Ds_pt",
    "Dsmu_pt_corr",
    "Dsmu_eta_corr",
    "mu_pt_corr",
    "mu_eta",
    "dR_Ds_mu",
    "dR_Kp_mu",
    "dR_Km_mu",
    "dR_pi_mu",
    "dR_phi_mu",
    # Ds+mu mass: H-mass constrained for signal, continuum for background
    "Dsmu_mass_corr",
    # --- jet/b-tag environment, added 2026-08-03, re-screened 2026-08-10 for the
    # SS-background switch (screen_new_jet_vars.py, b-veto applied, SS as the
    # reference instead of W+jets MC -- see new_jet_vars_screening.png).
    "Ds_nearestjet_bscore",
    "Ds_nearestjet_dR",
    "jet_leadpt_corr",
    # jet_ht dropped outright (see module comment above) rather than moved to
    # jet_ht_pt20_corr, a different threshold -- unlike num_jets_pt30 above, jet_ht had
    # no same-threshold corrected sibling already occupying a second feature slot.
    # n_b_jets_{loose,medium}_dR04: dropped 2026-08-04 after ranking dead last
    # (41st/42nd of 42, importance 0.0016/0.0020) against W+jets MC -- which has no
    # extra soft b-jet activity to distinguish. Re-added: against SS (the actual
    # candidate background now) they show a real gap in the 0-jet vs 1-jet fraction
    # that never showed up against W+jets. jet_maxbscore was also re-checked and
    # dropped again -- its pre-veto separation turned out to come almost entirely
    # from events the tight b-veto already removes, so little is left post-veto.
    # No corrected counterpart in NanoMVA/Samples/ (jet-quality counts, not JES/JER
    # quantities), so unchanged.
    "n_b_jets_loose_dR04",
    "n_b_jets_medium_dR04",
    # --- DERIVED features, added 2026-09-21. These are NOT branches: they are computed
    # in build_matrix() from DERIVED below. Appended at the END on purpose so the first
    # 41 positions keep the indices they had -- `names` is positionally bound to
    # scaler_mean/scale.npy and the first nn.Linear, so a stable prefix makes any future
    # diff of feature_names.txt readable.
    #
    # dphi(MET, X) folded to [0, pi]. The motivation is that the signal topology is a
    # 125 GeV H recoiling into Ds + W(->mu nu): the neutrino (hence MET) is produced with
    # the muon on the W side, back-to-back against the Ds. Combinatorial background has no
    # such angular relation, and for W+jets the MET is a real W neutrino but the Ds is a
    # fragmentation product with no fixed orientation to it. met_pt/W_tmass already carry
    # the MET magnitude; none of the existing 41 features carry its DIRECTION relative to
    # anything -- met_phi was never a feature (a bare azimuth is meaningless on its own).
    # Untested as of the change: no AUC or correlation measurement was made before adding.
    "dphi_met_mu",
    "dphi_met_Ds",
    "dphi_met_Dsmu",
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
    "mu_ipxysig",
    # Isolation: the heavy-tailed ABSOLUTE quantities -- the twelve cone sums and the
    # three PFIso. The three RelIso stay OUT: they are ratios, already O(1), which is the
    # convention this list has always followed. Added 2026-09-22 along with the eleven new
    # isolation features; every one has min = 0 / zero negative / zero non-finite across
    # the four merged samples, so log1p is safe and no sentinel is involved.
    "Ds_IsoR04_sumChargedHadronPt",
    "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt",
    "Ds_IsoR04_sumPUPt",
    "Ds_IsoR04_PFIso",
    "Ds_IsoR03_sumChargedHadronPt",
    "Ds_IsoR03_sumNeutralHadronEt",
    "Ds_IsoR03_sumPhotonEt",
    "Ds_IsoR03_sumPUPt",
    "Ds_IsoR03_PFIso",
    "mu_IsoR04_sumChargedHadronPt",
    "mu_IsoR04_sumNeutralHadronEt",
    "mu_IsoR04_sumPhotonEt",
    "mu_IsoR04_sumPUPt",
    "mu_IsoR04_PFIso",
]

# Sentinel handling -- needed here but NOT in the BDT: trees just split around a
# -777, but StandardScaler would let it wreck the mean/variance of the whole
# feature. Each sentinel is mapped to the physically sensible edge of its range
# ("no jet nearby" / "tagger gave nothing"), keeping the input bounded.
#   Ds_nearestjet_bscore  -777 (no jet) and ~-10 (tagger unavailable, ~14% of rows)
#   Ds_nearestjet_dR      -777 (no jet)  -> far away
#   jet_leadpt_corr       -777 (no jet)  -> zero pt (same sentinel semantics as the
#                          pre-migration jet_leadpt; apply_corrections.py writes NULL
#                          = -777.0 for an event with no tight jet, same convention)
SENTINEL_FIX = {
    "Ds_nearestjet_bscore": -0.1,   # just below the [0,1] tagger range
    "Ds_nearestjet_dR": 5.0,        # beyond any real dR in the sample
    "jet_leadpt_corr": 0.0,
}


def _dphi(phi1, phi2):
    """|phi1 - phi2| folded into [0, pi]."""
    d = np.abs(np.asarray(phi1, dtype=np.float64) - np.asarray(phi2, dtype=np.float64))
    d = np.mod(d, 2.0 * np.pi)
    return np.minimum(d, 2.0 * np.pi - d)


# Features that are COMPUTED rather than read: {name: (fn(arrays_dict), source_branches)}.
# build_matrix() is the single place they are evaluated, so common.load() (training) and
# DNN/BDT predict.py (scoring) cannot drift -- before this existed all three independently
# did `column_stack([data[b] for b in names])`, which works only while every feature name
# is also a branch name.
#
# The source branches are what goes into FeatureSpec.load_branches; the derived NAMES are
# deliberately kept out of it, or uproot would raise on a key the file does not have.
# All four sources exist in every one of ALL_SAMPLES and in a stage-1 Prediction/ file
# (checked 2026-09-21), which is the other input path load() can be pointed at.
#
# Dsmu_phi_corr (not Dsstmu_phi_corr) for the DsstW channel too, matching how every other
# Ds+mu feature here (Dsmu_mass_corr / Dsmu_pt_corr / Dsmu_eta_corr) is defined.
DERIVED = {
    "dphi_met_mu":   (lambda a: _dphi(a["met_phi_corr"], a["mu_phi"]),
                      ("met_phi_corr", "mu_phi")),
    "dphi_met_Ds":   (lambda a: _dphi(a["met_phi_corr"], a["Ds_phi"]),
                      ("met_phi_corr", "Ds_phi")),
    "dphi_met_Dsmu": (lambda a: _dphi(a["met_phi_corr"], a["Dsmu_phi_corr"]),
                      ("met_phi_corr", "Dsmu_phi_corr")),
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
# producers (make_samples_nano.py) can share these names instead of hand-copying
# them -- this module pulls in torch, which they must not. B_VETO_BRANCH is
# n_bjet_corr as of the NanoMVA migration -- see eras.py for why.
B_VETO = B_VETO_BRANCH

# Muon selection in CORRECTED quantities, 2026-09-21. This REPLACES the old
# `num_os_mu_tight == 1 & num_ss_mu_tight == 0` condition, which was never a muon
# selection in the sense its name suggests: the producer counts ALL tight muons with
# no pt, eta or isolation threshold at all, so requiring "exactly one" was in effect an
# uncalibrated heavy-flavour tag on soft muons inside jets (measured composition of the
# extra muons: 86% below 10 GeV, 81% non-isolated, 58% inside a stored jet). It rejected
# 7.57% of OS data against 1.98% of signal, with no POG scale factor available for any
# of it -- MUO ID/ISO SFs start at 15 GeV and the tag lives below that. The corrections
# work deliberately replaced it with pass_muveto (prompt second-muon veto: tight ID,
# |eta| < 2.4, corrected pt >= 15, Loose PF iso), which make_samples_nano.py already
# applied when it skimmed NanoMVA/Samples/.
#
# The three values are NOT free parameters -- each is fixed by a correctionlib key's own
# binning, see feedback_sf_key_defines_the_cut memory:
#   pt  > 26    MUO trigger SF pt edges start at 26 with flow="error" (it RAISES below,
#               it does not clamp), and the SF must be evaluated on the corrected pt.
#   |eta| < 2.4 the SF eta binning ends there.
#   relIso04 < 0.15 = PFIsoTight, forced by the trigger SF denominator: every trigger key
#               in muon_Z.json has an isolation requirement in its DEN, so tight ID
#               leaves DEN_CutBasedIdTight_and_PFIsoTight as the only option. Measured
#               equivalent to nano Muon_pfIsoId >= 4 on 3334/3334 muons.
# All four cuts are already baked into every row of NanoMVA/Samples/ (verified: 100.0000%
# pass, zero failures, in all four merged files), so stating them here is a no-op
# re-affirmation exactly like B_VETO above -- not a new cut. They are stated rather than
# assumed so that a sample which was NOT pre-skimmed cannot silently skip them.
MU_PT_MIN, MU_ETA_MAX, MU_ISO_MAX = 26.0, 2.4, 0.15
MU_SEL_BRANCHES = ["mu_pt_corr", "mu_eta", "mu_isTight", "mu_IsoR04_RelIso"]

# mu_charge/pi_charge: the ss= consistency check in load(). NOT a selection -- see there.
CHARGE_BRANCHES = ["mu_charge", "pi_charge"]

# LEPTON_BRANCHES is still loaded, but NOTHING in this module selects on it any more.
# It stays only because eight macros do `[*common.SELECTION]` and then index
# num_os_mu_tight/num_ss_mu_tight themselves; dropping it here would turn those into a
# KeyError. THOSE MACROS ARE NOW STALE -- they still apply the retired veto by hand:
#   check_feature_datamc_gap, check_jetbin_closure_nanomva, draw_features,
#   draw_stage1_mass_by_score_cut, draw_stage1_dsmass_by_era,
#   draw_ds_mass_signal_data_wjets, check_os_vs_ss_features, screen_new_jet_vars
# (check_os_vs_ss_features is the one legitimate user: an OS-vs-SS charge study.)
SELECTION = list(dict.fromkeys(
    [*MU_SEL_BRANCHES, *CHARGE_BRANCHES, B_VETO, *LEPTON_BRANCHES]))


def selection_mask(a, *, bveto=True):
    """The analysis row selection: corrected muon cuts + b-jet veto.

    One function so load() and stage1_threshold() cannot drift -- they used to carry
    two hand-written copies of the same condition, and only one of them would have been
    noticed if it had been changed alone.

    bveto=False drops ONLY the b-jet veto and keeps the muon cuts. Two macros study the
    b-veto itself -- draw_nbjets_veto.py picks its working point and draw_lepton_bjet_vs_cut.py
    plots b-jet multiplicity with the veto lifted -- so applying it would collapse their own
    observable to a spike at zero. They still need the MUON selection to come from here
    rather than a hand-written copy, which is the entire point of this function.
    """
    m = ((a["mu_pt_corr"] > MU_PT_MIN)
         & (np.abs(a["mu_eta"]) < MU_ETA_MAX)
         & a["mu_isTight"].astype(bool)
         & (a["mu_IsoR04_RelIso"] < MU_ISO_MAX))
    return (m & (a[B_VETO] == 0)) if bveto else m


def check_charge_config(a, ss, where):
    """Assert the file's charge configuration matches what the caller asked for.

    This is an ASSERTION, not a selection. Which charge configuration a row has is a
    property of the FILE, fixed by the producer's sameSignPi flag: HDsPhiPiMuProducer.cc
    sets `pi_charge_ref = sameSignPi ? mu.charge() : -mu.charge()`, so the pion's charge
    relative to the muon is decided before any candidate is built. Measured:
    mu_charge == pi_charge is 100.0000% in Data_SS_merged and 0.0000% in all three OS
    files -- there is nothing to select, only a wrong-file mistake to catch.

    Worth having because that mistake is currently SILENT: the old condition applied to
    a mismatched file returned ~1% of its rows, a plausible-looking subset, with no error.
    """
    same_sign = a["mu_charge"] == a["pi_charge"]
    if not (same_sign.all() if ss else (~same_sign).all()):
        frac = 100.0 * same_sign.mean()
        sys.exit(f"{where}: ss={ss} does not match this file's charge configuration "
                 f"(mu_charge == pi_charge in {frac:.4f}% of rows; expected "
                 f"{'100%' if ss else '0%'}). Wrong sample for this config?")


class FeatureSpec:
    """A feature set plus everything positionally bound to it.

    Until a second classifier appeared, FEATURE_NAMES was a module global with
    LOG_FEATURES / SENTINEL_FIX / DIRA_IDX / LOG_IDX / LOAD_BRANCHES all derived
    from it -- implicitly "the one feature set". The stage-2 resonant classifier
    (feature_set "prodtopo_v1") makes the feature set a per-training property, so
    this bundles the six and every function that used a global now takes `feats=`.

    `names` IS LOAD-BEARING: positionally bound to scaler_mean/scale.npy, the first
    nn.Linear(n_in, 128) and the feature_names.txt stamp. Reordering it without
    retraining produces garbage scores and raises nothing on its own --
    check_model_dir() / check_feature_stamp() below turn that into a readable
    failure. `names` is NOT copied, so FEATURE_SETS["f55"].names is FEATURE_NAMES
    itself and the back-compat shim below cannot drift from it.

    `dira` is a feature name or None. None (stage-2 has no `dira`) drops both the
    `X[:, dira_idx] > -1` finite-mask term and the log(1.01 - dira) transform in
    preprocess() -- so the non-finite fraction a stage differs legitimately.
    """

    def __init__(self, names, *, log=(), sentinel=None, dira=None, derived=None):
        self.names = names
        self.log = list(log)
        self.sentinel = dict(sentinel or {})
        self.dira = dira
        self.derived = dict(derived or {})
        self.log_idx = [names.index(f) for f in self.log]
        self.dira_idx = names.index(dira) if dira is not None else None
        # Derived names are NOT branches: read their sources instead, or uproot raises
        # on a key the file does not have. Sources may repeat across derived features
        # (met_phi_corr feeds all three dphi), hence the dict.fromkeys dedup.
        _branch = [f for f in names if f not in self.derived]
        _src = [s for _, srcs in self.derived.values() for s in srcs]
        self.load_branches = list(dict.fromkeys(_branch + _src + SELECTION))
        for f in (list(self.sentinel) + self.log + list(self.derived)
                  + ([dira] if dira is not None else [])):
            if f not in names:
                raise ValueError(f"FeatureSpec: {f!r} is not in names")


# Stage-2 (resonant W+jets) production-topology feature set. Branch-backed only.
# Excludes ALL decay kinematics (would sculpt the fit peaks, same reason as f41)
# AND all lifetime/IP: W+jets true-Ds is ~89% prompt like signal, while the missing
# tt background is displaced b->Ds, so an IP-heavy stage-2 would rank displaced
# background HIGHER -- wrong direction. Revisit only with a tt MC sample. Also drops
# Ds_IsoR04_sumPUPt / sumNeutralHadronEt / sumPhotonEt and n_b_jets_* (pileup
# sensitive, and the W+jets PU reweight is still provisional). See the stage-2 plan.
#
# NanoMVA migration: Dsmu_mass/Dsmu_pt/num_jets_pt30 switched to their _corr
# counterparts, same as f41 (default-to-corrected convention). jet_ht dropped
# outright, same reason as f41 -- no same-threshold corrected replacement exists,
# and it is 11 features -> 10.
PRODTOPO_V1 = [
    # 2026-09-22, USER FRAMEWORK: four classes, and one rule deciding membership --
    # no daughter-level quantities, only the Ds and things above it. That is what this
    # set always claimed ("production topology and not decay kinematics") but did not
    # follow: Kp_ipsig/Km_ipsig were added 2026-09-21 in direct contradiction of it.
    #
    # The rule has a measured justification, not just symmetry: of the whole displacement
    # family ONLY Kp_ipsig sculpts the fit variable. Within mode 3-5 the SR/sideband
    # density ratio moves 1.119 -> 1.305 -> 1.482 (+-0.031/0.045/0.071, 6.8 sigma) as the
    # cut tightens, while Ds_FDsig is flat (0.974 +- 0.053, 0.5 sigma) and Ds_ipsig and
    # isolation are flat to under 1 sigma. Mechanism: 93% of mode 3-5 is mode 5 (K->pi
    # misassignment), so the track ASSIGNED as K+ has its parameters coupled to the
    # misreconstructed mass; Ds_FDsig is a whole-Ds quantity with no such coupling.
    # HOW to measure this correctly is in the feedback_sculpting_measure memory: an
    # SR/sideband ratio taken across MIXED components measures composition, not shape.

    # --- class 1: isolation (8). The four sums are a BASIS -- PFIso == sumCH + max(0,
    # sumNH + sumPh - 0.5*sumPU) and RelIso == PFIso / Ds_pt, both verified bitwise, so
    # once the sums and Ds_pt are present the two PFIso and two RelIso are determined.
    # Dropped: for a BDT a redundant column is dead weight. (stage-1 keeps its own for
    # the OPPOSITE reason -- a net cannot easily learn a ratio, so handing it the
    # combination is worth one input.)
    "Ds_IsoR04_sumChargedHadronPt", "Ds_IsoR04_sumNeutralHadronEt",
    "Ds_IsoR04_sumPhotonEt", "Ds_IsoR04_sumPUPt",
    "Ds_IsoR03_sumChargedHadronPt", "Ds_IsoR03_sumNeutralHadronEt",
    "Ds_IsoR03_sumPhotonEt", "Ds_IsoR03_sumPUPt",

    # --- class 2: prompt Ds (3). Kp_ipsig/Km_ipsig REMOVED: daughter-level, sculpting
    # (above), and r = 0.754 with Ds_FDsig, which measures the same flight directly and
    # beats them inside every Ds_pt bin. Ds_ipsig ADDED: it asks whether the Ds MOMENTUM
    # points back at the PV (the b->Ds question) rather than how far it flew, and is
    # nearly independent of the rest (r = 0.180 with Ds_FDsig, 0.146 with Kp_ipsig,
    # -0.380 with dira). AUC 0.519 is low BY CONSTRUCTION -- the ceiling for a variable
    # tagging f ~ 11% of the background is 0.5 + f/2 = 0.555 -- and it is the only
    # displacement variable that does NOT collapse inside the SR (0.5187 -> 0.5172, vs
    # Kp_ipsig 0.6701 -> 0.5943 and Ds_FDsig 0.6781 -> 0.6118). The collapse happens
    # because the SR is 82.9% real Ds while the sideband is 92% D/reflection, so the
    # ctau-based Ds-vs-D separation these variables provide is work the mass window has
    # already done.
    "Ds_FDsig", "Ds_ipsig", "dira",

    # --- class 3: jet / event activity (5).
    "met_sumEt", "num_jets_pt30_corr", "Ds_nearestjet_dR",
    "n_b_jets_loose_dR04", "n_b_jets_medium_dR04",

    # --- class 4: Ds-mu back-to-back (8). Ds_pt ADDED: post-eps90 AUC 0.7189, the highest
    # of anything measured on this problem, and a whole-Ds quantity (Kp_pt + Km_pt +
    # dR_Kp_Km together rebuild m(KK); Ds_pt alone does not, so it cannot sculpt the fit
    # peaks). CAVEAT: CLAUDE.md 2026-08-11 records Ds_pt as the worst data/MC-mismatched
    # feature in the set (single-cut eps_data/eps_MC 2.52, rising to 4.57) -- harmless
    # while training MC-vs-MC, but it propagates when the model is applied to data.
    # mu_pt_corr ADDED in place of Dsmu_Ht_corr, which was exactly Ds_pt + mu_pt_corr to
    # float32: with both originals present that sum is a fixed 1:1 combination the BDT can
    # build itself and weight as it likes. Dsmu_pt_corr does NOT substitute for it -- that
    # is the VECTOR sum, so given Ds_pt it does not determine mu_pt (a dphi is missing).
    "dR_Ds_mu", "dphi_met_mu", "dphi_met_Ds", "dphi_met_Dsmu",
    "Dsmu_mass_corr", "Dsmu_pt_corr", "Ds_pt", "mu_pt_corr",
]

FEATURE_SETS = {
    "f55": FeatureSpec(FEATURE_NAMES, log=LOG_FEATURES, sentinel=SENTINEL_FIX,
                       dira="dira", derived=DERIVED),
    "prodtopo_v1": FeatureSpec(
        PRODTOPO_V1,
        # log1p the heavy-tailed isolation SUMS; the RelIso / mass / dR / met / dphi /
        # jet-count features are already O(1) or bounded. The three sums re-added in
        # 2026-09-21 get the same treatment they have in LOG_FEATURES for stage-1, so a
        # variable is not scaled one way in one stage and another way in the other.
        # (Cosmetic while stage-2 is a BDT -- trees are invariant under a monotone
        # transform -- but it is not BDT-only by construction, and preprocess() is shared.)
        log=["Ds_IsoR04_sumChargedHadronPt", "Ds_IsoR04_sumPUPt",
             "Ds_IsoR04_sumNeutralHadronEt", "Ds_IsoR04_sumPhotonEt",
             "Ds_IsoR03_sumChargedHadronPt", "Ds_IsoR03_sumPUPt",
             "Ds_IsoR03_sumNeutralHadronEt", "Ds_IsoR03_sumPhotonEt",
             "Ds_FDsig", "Ds_ipsig"],
        # Ds_nearestjet_dR IS in this set (Ds_nearestjet_bscore and jet_leadpt_corr
        # are not) and needs the same "-777 (no jet) -> far away" sentinel as f55.
        sentinel={"Ds_nearestjet_dR": 5.0},
        # dira="dira" (was None): dira is in this set as of 2026-09-21, so it needs the
        # same log(1.01 - dira) transform stage-1 gives it. Without this it stays a raw
        # cosine piled up against 1.0 and is numerically almost a constant.
        dira="dira",
        # The SAME DERIVED dict object as f55, so there is one definition of the dphi
        # formulas for both stages instead of a second copy that could drift. Sharing it
        # works only because all three of its keys are in this set's names too --
        # FeatureSpec REQUIRES every derived key to appear in `names` and raises
        # otherwise, so a future set that wants only some of them needs a subset dict,
        # not this one.
        derived=DERIVED),
}

# Back-compat shims: the stage-1 macros import these names directly and must keep
# working untouched. They are the f41 spec's attributes, by identity for `names`.
STAGE1 = FEATURE_SETS["f55"]
LOAD_BRANCHES = STAGE1.load_branches
DIRA_IDX = STAGE1.dira_idx
LOG_IDX = STAGE1.log_idx

# XGBoost hyperparameters per feature set. f41 keeps the values tuned for the
# (pre-migration 43-, now 41-) feature stage-1 problem on a much larger sample
# (unchanged since 2026-08-14 -- two fewer features does not warrant a re-tune).
# prodtopo_v1 (stage-2) trains on N_eff ~ O(10^3) background with 28 features, so it
# is shallower and more regularised: depth 3, larger min_child_weight. n_estimators
# raised 800 -> 1500 on 2026-09-09: the first seed-42 run hit the 800 cap without
# early-stopping (val AUC still rising ~1e-4/round at round 799), i.e. underfit.
BDT_PARAMS = {
    "f55": dict(n_estimators=2000, max_depth=5, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8, min_child_weight=10),
    # n_estimators history: 800 -> 1500 (2026-09-09) -> 3000 (2026-09-21). Each raise
    # was forced by the SAME symptom, not by tuning: the run ends with best_iteration at
    # the cap and early_stopping_rounds=30 never fires, i.e. val AUC is still climbing
    # when the trees run out. At 1500 with the new 28-feature set both stage-2 configs
    # stopped at 1497/1499 with a tail slope of ~3e-6/round (extrapolating to about
    # +0.003 AUC over another 1000 rounds), so the 10-feature-era cap was underfitting a
    # set 2.8x wider. Read a best_iteration equal to n_estimators-1 as "raise this
    # again", never as convergence.
    "prodtopo_v1": dict(n_estimators=3000, max_depth=3, learning_rate=0.03,
                        subsample=0.7, colsample_bytree=0.9, min_child_weight=40,
                        reg_lambda=5.0),
}

# ========= 训练配置 =========
# The one seed the canonical models were trained with. Lives here rather than in
# DNN.py because out_dir_for() needs it to decide when to add a _seed suffix.
SEED_DEFAULT = 42

# Signal is ONE file (NanoMVA/Samples/Signal_merged.root) since the 2026-09-17 NanoMVA
# migration, not two per-channel files -- DsW/DsstW are now a `channel` BRANCH
# (0/1) inside it, same as era. SIGNALS stays a list (not a bare string) because every
# consumer (stage1_threshold(), DNN.py, BDT.py) already loops `for name in SIGNALS`
# to concatenate signal pieces; a one-element list keeps that loop unchanged.
#
# The weight is now ABSOLUTE (sigma*L/sum(w) per (channel, era) piece, computed by
# make_samples_nano.py), not the old RELATIVE anchored-at-era-1 SIGNAL_ERA_WEIGHT --
# see eras.py. It already encodes the correct DsW:DsstW and preEE:postEE mixture, so
# there is no separate per-channel or per-era factor to apply here any more; read it
# with common.load(..., extra_branches=("weight",)).
SIGNALS = ["Signal_merged"]

# Everything in NanoMVA/Samples/, i.e. everything predict.py scores. ONE list,
# imported by both DNN/predict.py and BDT/predict.py so they cannot drift apart (see
# the historical note in project_samples_merge_era_branch memory -- BDT used to carry
# a 13-sample-short hand copy).
#
# OS vs SS stay separate FILES (different candidate reconstructions, ~13% event
# overlap in the pre-migration samples) -- everything else that used to be a
# filename axis (channel, era, jet_bin, GenTag's cand_gen_mode) is now a BRANCH
# inside one of these four files.
ALL_SAMPLES = ["Signal_merged", "WtoLNu2Jets_merged", "Data_OS_merged", "Data_SS_merged"]

# WtoLNu2Jets_merged holds all three jet-multiplicity bins of the SAME physical
# process as a `jet_bin` branch, and its per-event `weight` column already carries
# the amcatnloFXFX sign and the correct sigma*L/sum(w) normalization -- no per-bin
# weight table exists or is needed (see the pre-migration WJETS_WEIGHTS post-mortem
# in project_genweight_missing memory for why a per-bin constant cannot work at all).
WJETS_JET_BINS = (0, 1, 2)

# The two background configurations. `_f55` is a feature-set generation tag: a
# future feature change gets its own suffix rather than clobbering these.
#
# out_name vs out_dir: `out_name` is the bare directory name, which out_dir_for()
#   composes under any tool (DNN/, BDT/) and optional seed suffix. `out_dir` stays
#   as the literal DNN path it has always been -- seven macros read cfg["out_dir"]
#   directly (check_data_vs_wjets, check_feature_tail, check_jetbin_closure,
#   check_os_minus_ss, compare_models, draw_score_os_ss, draw_lepton_bjet_vs_cut),
#   so it is kept, not renamed. The two must agree; the assertion below enforces it.
#
# parity: SS data holds out the odd-event half so a later sculpting check on
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
        # NanoMVA migration (2026-09-17): file renamed Data_SS -> Data_SS_merged, still
        # one file with both eras inside (era branch). Data carries no luminosity
        # weight (weight is 1.0 for every row by the production convention), so the
        # `weighted=False` / flat-1.0 behavior below is unaffected by the migration.
        #
        # Everything below this line describes the pre-migration ssdata/Data_SS.root
        # measurements (era balance, the pileup-mismodelling bias). Re-measure after
        # the first post-migration retrain -- the corrections pipeline now includes a
        # proper w_pu_nom weight, which did not exist when this was written, and the
        # underlying jet/MET quantities changed (see project_nanoaod_migration and
        # project_corrections_pilot_implementation memory).
        #
        # More SS statistics is exactly what the 70%-working-point study identified as
        # the bottleneck (2026-08-14: ~94 hard events of training equivalent at
        # eps_sig=70%). Measured 2026-09-08 (pre-migration), the era composition of the
        # two training classes was identical to the decimal (23.2%/76.8% both), so era
        # carried no class information. THE BIAS THAT WAS REAL: signal is MC and
        # background is real data, so every data/MC mismodelling is a free
        # discriminant -- pileup was the measured case (<Ds_IsoR04_sumPUPt> 9.65 GeV MC
        # vs 11.99/14.03 GeV preEE/postEE data, no pileup reweighting existed anywhere).
        # See project_era_pileup_mismodelling and analysis/macros/check_era_split.py.
        samples={"Data_SS_merged": 1.0},
        ss=True,
        parity=0,
        weighted=False,          # real data: no per-process lumi reweighting
        label="Data SS",
        is_data=True,
        out_name="dsw_dsstw_ssdata_f55",
        out_dir=f"{MVA_DIR}/dsw_dsstw_ssdata_f55",
    ),
    "wjets": dict(
        # One file, both eras and all three jet bins inside (WtoLNu2Jets_merged, was
        # WtoLNu2Jets pre-migration). `weighted=True` reads the per-event `weight`
        # column (was `Lumi_reweight`) -- sigma*L/sum(w), carrying the amcatnloFXFX
        # sign. The value here is unused when weighted=True; it stays 1.0 for shape.
        samples={"WtoLNu2Jets_merged": 1.0},
        ss=False,
        parity=None,
        weighted=True,
        label="W+jets MC",
        is_data=False,
        out_name="dsw_dsstw_wjets_f55",
        out_dir=f"{MVA_DIR}/dsw_dsstw_wjets_f55",
    ),
    "wjetsres": dict(
        # Stage 2 of the two-stage MVA: signal MC vs W+jets RESONANT background (real
        # D/Ds), which the stage-1 f41 classifier cannot remove -- in the f41 feature
        # set a true W+jets Ds->phipi and a true signal Ds->phipi are the same object.
        # Background = WtoLNu2Jets_merged rows with cand_gen_mode in {1,2,3,4,5}
        # (1-2 true Ds, 3-4 true D, 5 D->Kpipi misID reflection), MC-vs-MC. The
        # cand_gen_mode branch lives directly in WtoLNu2Jets_merged as of the NanoMVA
        # migration -- no separate GenTag file/production any more.
        #
        # stage1_cut: training rows must first pass stage-1 at eff_sig -- the 80% of
        # modes 1-5 that stage-1 removes is a DIFFERENT population in the production-
        # topology features stage-2 uses, so training on the unfiltered sample trains
        # on the wrong problem. eff_sig = 0.90 (pre-migration value, 2026-09-09
        # fom_scan_stage1.py) -- MUST be re-derived once stage-1 is retrained on the
        # NanoMVA sample (Step 3 of the migration plan), the training population and
        # statistics are not the same any more. common.load() reads the training rows
        # from stage-1's Prediction/ dir and applies score_dsmu > threshold; predict.py
        # does the same input redirect but applies NO threshold (cut goes downstream).
        #
        # feature_set prodtopo_v1: branch-backed production-topology features only
        # (10 as of the migration, jet_ht dropped -- see FEATURE_SETS above). Excludes
        # ALL lifetime/IP -- W+jets true-Ds is ~89% prompt like signal, the missing tt
        # background is displaced b->Ds, so an IP-heavy stage-2 would rank displaced
        # background HIGHER. Revisit only with a tt MC sample.
        #
        # score_branch score_dsmu_res so the stage-1 score_dsmu and this coexist in one
        # Prediction file. Algorithm is BDT (analysis/BDT/) -- N_eff ~ O(10^3), tabular.
        samples={"WtoLNu2Jets_merged": 1.0},   # unused when weighted=True
        ss=False,
        parity=None,
        weighted=True,
        gen_mode_filter=(1, 2, 3, 4, 5),
        stage1_cut=dict(key="ssdata", tool="DNN", seed=42, eff_sig=0.90),
        feature_set="prodtopo_v1",
        score_branch="score_dsmu_res",
        label="W+jets resonant (GenTag 1-5, post stage-1 90%)",
        is_data=False,                          # MC-vs-MC -> "Simulation"
        # DNN-convention path (the out_dir == MVA_DIR/out_name assertion below pins it).
        # A BDT-only entry still needs the field; every stage-2 macro must resolve the
        # real dir with common.out_dir_for("wjetsres", tool="BDT"), never cfg["out_dir"].
        out_name="dsw_dsstw_wjetsres_f24_wp90",
        out_dir=f"{MVA_DIR}/dsw_dsstw_wjetsres_f24_wp90",
    ),
    "wjetsres_nostage1": dict(
        # Comparison arm for Step 4 of the NanoMVA migration plan: same background
        # population (WtoLNu2Jets_merged, cand_gen_mode in {1,2,3,4,5}) and the same
        # prodtopo_v1 feature set as "wjetsres" above, but stage1_cut=None -- trains
        # on ALL of modes 1-5, not just the 20% that survives stage-1 at eff_sig=0.90.
        # Tests whether pre-filtering on stage-1 helps stage-2 learn the production-
        # topology problem, or whether it just throws away training statistics stage-2
        # could have used.
        #
        # stage1_cut=None means common.load()/predict.py read WtoLNu2Jets_merged
        # straight from SAMPLES_DIR, with NO dependency on any stage-1 model existing
        # or having been scored -- this config can train and be scored independently
        # of "ssdata"/stage-1, in parallel. The cost: its own Prediction/ files carry
        # score_dsmu_res but NOT score_dsmu (stage-1's score never enters this path),
        # so the joint stage-1 x stage-2 rejection has to be recovered later by
        # joining this config's Prediction/WtoLNu2Jets_merged.root against stage-1's
        # own Prediction/WtoLNu2Jets_merged.root on (run, luminosityBlock, event) --
        # not by reading through stage-1's directory the way "wjetsres" does.
        samples={"WtoLNu2Jets_merged": 1.0},   # unused when weighted=True
        ss=False,
        parity=None,
        weighted=True,
        gen_mode_filter=(1, 2, 3, 4, 5),
        stage1_cut=None,
        feature_set="prodtopo_v1",
        score_branch="score_dsmu_res",
        label="W+jets resonant (GenTag 1-5, no stage-1 pre-filter)",
        is_data=False,
        out_name="dsw_dsstw_wjetsres_nostage1_f24",
        out_dir=f"{MVA_DIR}/dsw_dsstw_wjetsres_nostage1_f24",
    ),
}

# Stage-2 seed-coupling check (2026-09-16): "wjetsres" above always conditions on
# stage-1 SEED 42's WP1=90% threshold, regardless of what seed stage-2 itself trains
# with -- so a stage-2 seed-spread study run behind it only measures stage-2's own
# training noise, not the fact that a different stage-1 seed passes a ~25%-different
# set of events at the same nominal WP1=90% (Jaccard 0.75, vs 1.09% on the aggregate
# yield -- membership, not yield, is what defines stage-2's training population).
# These three entries are identical to "wjetsres" except stage1_cut points at a
# different stage-1 seed's own Prediction/ dir and threshold; train stage-2 behind
# each (itself at seed=42, so only the stage-1 side varies) and compare the resulting
# rejection spread to the fixed-stage-1 spread already measured from wjetsres seeds
# 101/102/103. Requires dsw_dsstw_ssdata_f55_seed{101,102,103} to exist first -- the
# pre-migration f43 seed models this study originally measured were deleted with the
# rest of the pre-NanoMVA model directories (2026-09-17), so this needs re-running,
# not just re-pointing.
for _s1seed in (101, 102, 103):
    _name = f"wjetsres_s1seed{_s1seed}"
    TRAININGS[_name] = dict(
        TRAININGS["wjetsres"],
        stage1_cut=dict(key="ssdata", tool="DNN", seed=_s1seed, eff_sig=0.90),
        label=f"W+jets resonant (GenTag 1-5, post stage-1 seed{_s1seed} 90%)",
        out_name=f"dsw_dsstw_{_name}_f24_wp90",
        out_dir=f"{MVA_DIR}/dsw_dsstw_{_name}_f24_wp90",
    )

for _k, _c in TRAININGS.items():
    assert _c["out_dir"] == f"{MVA_DIR}/{_c['out_name']}", f"{_k}: out_dir/out_name disagree"

# Every config resolves to a FeatureSpec and an output score-branch name. Existing
# entries omit both and get the stage-1 defaults; a stage-2 entry sets feature_set
# and score_branch explicitly. Not `.get()` at the call site: a future entry that
# forgets score_branch should fail loudly rather than silently write into stage-1's
# "score_dsmu".
for _c in TRAININGS.values():
    _c.setdefault("feature_set", "f55")
    _c["features"] = FEATURE_SETS[_c["feature_set"]]
    _c.setdefault("score_branch", "score_dsmu")


def resolve(key):
    """Look up a background config, failing loudly on a typo."""
    if key not in TRAININGS:
        sys.exit(f"unknown background '{key}'; expected one of: {', '.join(TRAININGS)}")
    return TRAININGS[key]


TOOLS = ("DNN", "BDT")


def out_dir_for(key, tool="DNN", seed=None):
    """Model directory for one (background, tool, seed) combination.

    seed=None or the default 42 gives the canonical path -- i.e. the existing
    dsw_dsstw_*_f55 dirs, which the seed-spread study must not touch. Any other
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
def load(name, *, feats=STAGE1, era=None, jet_bin=None, truth_filter=False, ss=False,
         parity=None, gen_mode_filter=None, stage1=None, extra_branches=()):
    """Load one Samples/ file, apply the analysis selection, return the feature matrix.

    feats:   which FeatureSpec to build the matrix for (default = stage-1 f41).
    era:     keep only this era id (see eras.ERAS). None = every era in the file.
    jet_bin: keep only this jet bin. None = every bin. Errors if the file has no
             jet_bin branch, rather than silently returning everything.
    ss:      select the same-sign muon configuration instead of opposite-sign.
    parity:  keep only rows with event % 2 == parity (None = keep all).
    gen_mode_filter: keep only rows whose cand_gen_mode is in this iterable (the
             GenTag truth tag). Errors if the file has no cand_gen_mode branch.
    stage1:  a stage1_cut dict {key, tool, seed, thr, ...}. When set, the file is
             read from that stage-1 model's Prediction/ dir (a superset of Samples/,
             carrying score_dsmu) instead of Samples/, and rows are additionally
             required to have score_dsmu > stage1["thr"]. TRAINING-TIME ONLY --
             predict.py does its own input redirect and applies no threshold.
    Returns (X, extras) where extras is a dict of any extra_branches requested,
    already masked in step with X.
    """
    branches = list(feats.load_branches)
    if truth_filter:
        branches += ["Ds_match", "mu_match"]
    if parity is not None:
        branches += ["event"]
    if era is not None:
        branches += ["era"]
    if jet_bin is not None:
        branches += ["jet_bin"]
    if gen_mode_filter is not None:
        branches += ["cand_gen_mode"]
    if stage1 is not None:
        branches += ["score_dsmu"]
    branches += list(extra_branches)
    branches = list(dict.fromkeys(branches))

    if stage1 is None:
        src = f"{SAMPLES_DIR}/{name}.root"
    else:
        src = (f"{out_dir_for(stage1['key'], stage1.get('tool', 'DNN'), stage1.get('seed'))}"
               f"/Prediction/{name}.root")

    with uproot.open(src)[TREE] as t:
        keys = t.keys()
        if jet_bin is not None and "jet_bin" not in keys:
            sys.exit(f"load({name!r}, jet_bin={jet_bin}): that file has no jet_bin branch. "
                     f"Only the jet-binned W+jets samples do -- incl and GenTag-merged "
                     f"files are separate productions, not a fourth bin.")
        if gen_mode_filter is not None and "cand_gen_mode" not in keys:
            sys.exit(f"load({name!r}, gen_mode_filter=...): {src} has no cand_gen_mode "
                     f"branch -- only the GenTag production carries the truth tag.")
        a = t.arrays(branches, library="np")

    n_in = len(a[B_VETO])
    check_charge_config(a, ss, f"load({name!r})")
    m = selection_mask(a)
    if truth_filter:
        m &= a["Ds_match"].astype(bool) & a["mu_match"].astype(bool)
    if parity is not None:
        m &= (a["event"] % 2) == parity
    if era is not None:
        m &= a["era"] == era
    if jet_bin is not None:
        m &= a["jet_bin"] == jet_bin
    if gen_mode_filter is not None:
        m &= np.isin(a["cand_gen_mode"], list(gen_mode_filter))
    if stage1 is not None:
        m &= a["score_dsmu"] > stage1["thr"]

    X = build_matrix(a, feats, mask=m)
    print(f"  {name:18s} {n_in:>9,} -> {m.sum():>9,} rows after selection")
    return X, {b: a[b][m] for b in extra_branches}


def build_matrix(a, feats=STAGE1, mask=None):
    """Assemble the feature matrix from a dict of arrays, evaluating derived features.

    THE single place a feature matrix is built. load() (training) and both predict.py
    (scoring) call it, so the column order, the derived-feature formulas and the dtype
    are identical on every path by construction. Before this existed each of the three
    did its own `column_stack([data[b] for b in feats.names])`, which silently assumed
    every feature name is also a branch name -- true until the first derived feature.

    Computed on the full column and masked afterwards. dphi is elementwise so the two
    orders agree; fixing one order here keeps it that way for any future derived
    feature that is not.

    float64 explicitly: the two predict.py already cast, load() relied on column_stack's
    dtype promotion, so training and scoring could in principle have differed.
    """
    cols = []
    for f in feats.names:
        v = feats.derived[f][0](a) if f in feats.derived else a[f]
        v = np.asarray(v, dtype=np.float64)
        cols.append(v if mask is None else v[mask])
    return np.column_stack(cols)


def stage1_threshold(stage1_cut):
    """The score_dsmu cut value for a stage-1 signal efficiency.

    stage1_cut: {key, tool, seed, eff_sig}. Reads the two signals from that stage-1
    model's Prediction/ dir, applies OS + b-veto + truth match + score_dsmu > -1,
    and returns (thr, n_sig_rows) where thr is the UNWEIGHTED (1 - eff_sig) quantile
    of the pooled signal score -- unweighted to match draw_gentag_mass_modes.py and
    every working-point macro. The weighted quantile is printed alongside as a check.

    This defines the stage-2 TRAINING POPULATION; DNN.py / BDT.py stamp the returned
    value into stage1_filter.txt, and it moves if stage-1 is ever re-scored.
    """
    pred = (f"{out_dir_for(stage1_cut['key'], stage1_cut.get('tool', 'DNN'), stage1_cut.get('seed'))}"
            f"/Prediction")
    # Same selection as load(), via selection_mask() -- this used to be a second
    # hand-written copy of the condition, so the WP1 threshold could silently be derived
    # on a different signal subset than the one the model was trained on.
    want = ["score_dsmu", *SELECTION, "Ds_match", "mu_match", "weight"]
    want = list(dict.fromkeys(want))
    sc, w = [], []
    for name in SIGNALS:
        with uproot.open(f"{pred}/{name}.root")[TREE] as t:
            a = t.arrays(want, library="np")
        check_charge_config(a, False, f"stage1_threshold({name!r})")
        m = (selection_mask(a)
             & a["Ds_match"].astype(bool) & a["mu_match"].astype(bool) & (a["score_dsmu"] > -1))
        sc.append(a["score_dsmu"][m])
        w.append(a["weight"][m])
    sc = np.concatenate(sc)
    w = np.concatenate(w)
    q = 1.0 - stage1_cut["eff_sig"]
    thr = float(np.quantile(sc, q))
    order = np.argsort(sc)
    cw = np.cumsum(w[order])
    thr_w = float(sc[order][np.searchsorted(cw, q * cw[-1])])
    print(f"  stage1_threshold({stage1_cut['key']}, eff_sig={stage1_cut['eff_sig']}): "
          f"thr={thr:.6f} (weighted {thr_w:.6f}), from {len(sc):,} truth-matched signal rows")
    return thr, len(sc)


# ========= 预处理 =========
def preprocess(X, feats=STAGE1):
    """The one implementation of the chain: sentinel fix -> finite mask -> log.

    Returns (X, mask). X is modified in place for the sentinel step and the log
    transform is applied only to rows inside the mask (so log(1.01 - (-777)) is
    never evaluated). Callers differ in what they do with ~mask: DNN.py drops those
    rows, predict.py keeps them and writes the score = -777.

    Order matters. Sentinels must be regularized BEFORE the scaler ever sees them,
    and before the log, or the transform and the fitted mean/variance are garbage.

    feats with dira=None (stage-2) skips the dira finite-mask term and the
    log(1.01 - dira) transform -- the feature is simply not in the set.
    """
    for feat, fill in feats.sentinel.items():
        i = feats.names.index(feat)
        bad = X[:, i] < -1
        if bad.any():
            print(f"  sentinel: {feat:24s} {bad.sum():>9,} rows ({100*bad.mean():5.2f}%) -> {fill}")
            X[bad, i] = fill

    mask = np.isfinite(X).all(axis=1)
    if feats.dira_idx is not None:
        mask &= (X[:, feats.dira_idx] > -1)

    Xm = X[mask]
    if feats.log_idx:
        Xm[:, feats.log_idx] = np.log1p(Xm[:, feats.log_idx])
    if feats.dira_idx is not None:
        Xm[:, feats.dira_idx] = np.log(1.01 - Xm[:, feats.dira_idx])
    X[mask] = Xm
    return X, mask


# ========= 模型 =========
def neutralize(X, cols, feats=STAGE1):
    """Blank the named feature columns in place, and return the names actually blanked.

    The columns are set to a constant (their own mean), which carries exactly zero
    information, rather than being DELETED. That is the whole point: the feature count,
    the scaler size, the input-layer width and the feature_names.txt stamp all stay at
    43, so every drift guard and every downstream consumer keeps working untouched.
    Deleting columns instead would mean a per-training feature list threaded through
    every call site -- which is now exactly what FeatureSpec does. The nine ssdata_*
    pileup-control models that used this were removed (config and EOS) on 2026-09-16;
    no current TRAININGS entry sets `neutralize`.

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
    if not cols:
        return ()
    idx = [feats.names.index(f) for f in cols]      # loud on a typo
    X[:, idx] = 0.0
    return tuple(cols)


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
# Every guard takes feats=STAGE1 as a DEFAULT. A forgotten feats= on a stage-2 call
# is safe here BY CONSTRUCTION: every mismatch direction (stage-2 stamp vs f41, or
# f41 model vs stage-2 spec) is a length or content difference, never a silent pass.
FEATURE_FILE = "feature_names.txt"


def save_feature_names(out_dir, feats=STAGE1):
    """Stamp the model dir with the exact feature list (and order) it was trained on."""
    with open(f"{out_dir}/{FEATURE_FILE}", "w") as f:
        f.write("\n".join(feats.names) + "\n")


def check_feature_stamp(model_dir, feats=STAGE1):
    """Refuse to score with a model whose feature stamp does not match feats.names.

    Catches the two silent failure modes shared by every algorithm: a model dir from
    before the feature list was stamped (the frozen 40-feature v2, or the frozen
    42-feature BDT at dsw_wjets_test), and a list that was reordered without
    retraining -- the latter raises no exception on its own because feature order is
    only implicitly bound to the scaler, the first Linear layer, and the tree splits.

    Split out of check_model_dir so BDT/predict.py can reuse it without importing
    torch machinery it has no use for.
    """
    names = feats.names
    path = f"{model_dir}/{FEATURE_FILE}"
    if not os.path.exists(path):
        sys.exit(
            f"{path} not found.\n"
            f"This model dir predates feature_names.txt -- it is a frozen older build "
            f"(e.g. the 40-feature v2 at dsw_dsstw_wjets, or the 42-feature BDT at "
            f"dsw_wjets_test). Do NOT re-score it with the current "
            f"{len(names)}-feature code: the scores would be silently wrong. "
            f"Retrain into a new dir instead."
        )

    with open(path) as f:
        stored = [ln.strip() for ln in f if ln.strip()]
    if stored != names:
        n = min(len(stored), len(names))
        first = next((i for i in range(n) if stored[i] != names[i]), n)
        sys.exit(
            f"feature mismatch vs {path}\n"
            f"  stored: {len(stored)} features, common.py: {len(names)}\n"
            f"  first difference at index {first}: "
            f"stored={stored[first] if first < len(stored) else '<end>'} "
            f"vs common.py={names[first] if first < len(names) else '<end>'}"
        )


def check_n_features(model_dir, n_in, what, feats=STAGE1):
    """The algorithm-specific half: the trained object's own input width."""
    if n_in != len(feats.names):
        sys.exit(f"{what} expects {n_in} features but common.py declares "
                 f"{len(feats.names)} -- {model_dir} is a different feature set.")


def check_model_dir(model_dir, feats=STAGE1):
    """check_feature_stamp plus the DNN-specific scaler and input-layer checks.

    Returns the loaded state dict, so the file that was checked is the file loaded.
    """
    check_feature_stamp(model_dir, feats)
    n = len(feats.names)

    mean = np.load(f"{model_dir}/scaler_mean.npy")
    scale = np.load(f"{model_dir}/scaler_scale.npy")
    if mean.size != n or scale.size != n:
        sys.exit(f"scaler size {mean.size}/{scale.size} != {n} features "
                 f"in {model_dir} -- the model was trained on a different feature set.")

    state = torch.load(f"{model_dir}/phipimu_model.pt", map_location="cpu")
    check_n_features(model_dir, state["net.0.weight"].shape[1], "model input layer", feats)
    return state
