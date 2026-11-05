"""BTV: b-jet counting and the fixed-WP b-tag event weight.

The weight is btag_weight.cc, a C++ port of compute_btag_event_weight_multi_wp from
BTV/examples lib/btag_sf.py (BTV has no C++ version; the Python one needs Python >= 3.10),
compiled through ROOT. One WP ("T"), PNet B. For one WP it is
Pi_tagged SF * Pi_untagged (1-SF*eps)/(1-eps). Heavy and light jets are varied together
(variation_heavy = variation_light), as before. up/down are evaluated with
ignore_warnings=True, as in the official scripts/validate_btag_sf.py.

Jets: tight ID, |eta| < 2.5, dR(Ds, jet) > 0.4, nano Jet_pt >= 20 (the b-veto population;
the efficiency map is built on the same jets, build_btag_efficiency_map.py).
  w_btag_*      the b-veto weight: NaN where n_bjet > 0
  w_btag_any_*  valid at any n_bjet (ttcr); equals w_btag_* where n_bjet == 0
"""
import os
import sys

import awkward as ak
import correctionlib
import numpy as np

BTAG_WEIGHT_CC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "btag_weight.cc")

POG_DIR = "/cvmfs/cms-griddata.cern.ch/cat/metadata"
BTV_VERSION = "2025-08-20"
ALGO_KEY = "particleNet"          # correctionlib key of btagPNetB
EFF_KEY = "particleNet_efficiency"
WP = "T"
JET_PT_MIN = 20.0

BRANCHES = ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_btagPNetB", "Ds_eta", "Ds_phi"]
MC_BRANCHES = ["Jet_hadronFlavour"]

_ROOT = None


def _load_root():
    global _ROOT
    if _ROOT is None:
        import ROOT
        correctionlib.register_pyroot_binding()
        ROOT.gROOT.ProcessLine(f'#include "{BTAG_WEIGHT_CC}"')
        if not hasattr(ROOT, "BtagWeight"):
            sys.exit(f"{BTAG_WEIGHT_CC} did not compile")
        _ROOT = ROOT
    return _ROOT


def sf_path(campaign):
    return f"{POG_DIR}/BTV/{campaign}/{BTV_VERSION}/btagging.json.gz"


def sf_set(campaign, wp_tight):
    """BTV CorrectionSet; asserts its tight WP equals the one used to count b jets."""
    btv = correctionlib.CorrectionSet.from_file(sf_path(campaign))
    wp_json = btv[f"{ALGO_KEY}_wp_values"].evaluate(WP)
    if abs(wp_json - wp_tight) > 1e-6:
        sys.exit(f"tight WP mismatch: JSON {wp_json} vs {wp_tight}")
    return btv


def scored_jets(a, pt=None):
    """Mask of the jets the b-veto examines (jet pt = nano Jet_pt, or `pt` if given)."""
    pt = a["Jet_pt"] if pt is None else pt
    dphi = abs(a["Jet_phi"] - a["Ds_phi"])
    dphi = ak.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    dR = np.sqrt((a["Jet_eta"] - a["Ds_eta"]) ** 2 + dphi ** 2)
    return ((abs(a["Jet_eta"]) < 2.5) & (a["Jet_jetId"] >= 2) & (dR > 0.4)
            & (pt >= JET_PT_MIN))


def n_bjet(a, wp_tight, pt=None):
    m = scored_jets(a, pt)
    return ak.to_numpy(ak.sum(m & (a["Jet_btagPNetB"] > wp_tight), axis=1))


def weights(a, campaign, btv, eff_json, nb, pt=None, sf_vars=("nom", "up", "down"), stale="exit"):
    """{"w_btag_nom", ..., "w_btag_any_down"} for MC. pt: jet pt (default nano Jet_pt).
    sf_vars: subset of nom/up/down to compute. stale="count": do not exit on jets in eps=0/1
    bins (JME variations: jets move across pt bins of the nominal map; the official function
    then sets that event's weight to 1); the count is returned as cols["n_stale"]."""
    jpt = a["Jet_pt"] if pt is None else pt
    m = scored_jets(a, jpt)
    pt = np.ascontiguousarray(ak.to_numpy(ak.flatten(jpt[m])), dtype=np.float64)
    eta = np.ascontiguousarray(ak.to_numpy(ak.flatten(a["Jet_eta"][m])), dtype=np.float64)
    flav = np.ascontiguousarray(ak.to_numpy(ak.flatten(a["Jet_hadronFlavour"][m])), dtype=np.int32)
    disc = np.ascontiguousarray(ak.to_numpy(ak.flatten(a["Jet_btagPNetB"][m])), dtype=np.float64)
    counts = np.ascontiguousarray(ak.to_numpy(ak.num(jpt[m], axis=1)), dtype=np.int32)

    # A map built on this sample's own jets has no untagged jet in an eps=1 bin and no tagged
    # jet in an eps=0 bin. Either one means the map is stale (nano changed after the map was
    # built); the BTV function would then return 1 for the whole event or fail its check.
    eps = correctionlib.CorrectionSet.from_file(eff_json)[EFF_KEY].evaluate(WP, flav, np.abs(eta), pt)
    tag = disc > btv[f"{ALGO_KEY}_wp_values"].evaluate(WP)
    n_bad = int(((eps >= 1.0) & ~tag).sum() + ((eps <= 0.0) & tag).sum())
    if n_bad and stale == "exit":
        sys.exit(f"{n_bad} jets in eps=0/1 bins of {eff_json}: efficiency map stale, rebuild it")

    R = _load_root()
    wps = R.std.vector["std::string"]()
    wps.push_back(WP)
    bw = R.BtagWeight(sf_path(campaign), eff_json, ALGO_KEY, EFF_KEY, wps)
    rv = R.VecOps.AsRVec
    cols = {} if stale == "exit" else {"n_stale": n_bad}
    for suffix, s in (("nom", "central"), ("up", "up"), ("down", "down")):
        if suffix not in sf_vars:
            continue
        w = np.array(bw.weights(rv(pt), rv(eta), rv(flav), rv(disc), rv(counts), s, s, s != "central"),
                     dtype=float)
        if len(w) != len(nb):
            sys.exit(f"btag weight: {len(w)} events, expected {len(nb)}")
        cols[f"w_btag_any_{suffix}"] = w
        cols[f"w_btag_{suffix}"] = np.where(nb == 0, w, np.nan)
    return cols
