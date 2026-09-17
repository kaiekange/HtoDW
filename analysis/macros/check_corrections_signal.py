#!/usr/bin/env python3
"""Pilot: read the 6 POG JSONs (BTV/EGM/JME/LUM/MUO) with correctionlib and evaluate
nominal+up/down weights on the full preEE statistics of ONE signal sample.
TAU is skipped (see below).

Usage: check_corrections_signal.py <DsW_preEE|DsstW_preEE>

One sample per run, one condor job per sample (see run_corrections.sh /
condor_submit_corrections.py) -- matches this repo's convention (Samples/, DNN/, the
tuplizer condor scripts all submit one job per sample, never a combined job). The b-tag
efficiency map is therefore also per-sample, which is preferable anyway: DsW and DsstW
have different jet environments (DsstW adds a photon), and the official jetBtag.py
docstring itself recommends computing the efficiency map "at least separately for
signal and background processes." No cross-sample merge step exists or is needed.

Not production code: does not touch common.py / make_samples.py / eras.py, does not
write an analysis/corrections.py, does not propagate JES/JER/muon-scale variations into
the derived DNN features. See ~/.claude/plans/quiet-hatching-lerdorf.md for the full plan.

Must run inside the cmssw-el9 container (correctionlib 2.6.4). Measured, JSON by JSON,
against the el8-native correctionlib 2.1.0 in CMSSW_13_0_13: JME (jet_jerc, jetvetomaps,
jetid, met_xy), LUM puWeights and BTV btagging DO parse there, but all five MUO JSONs
fail ("Invalid edges array type" / "Unrecognized Content object nodetype"), as do EGM
electron.json and photon.json ("binning edges are not monotone increasing") and TAU
tau.json. Since this script needs MUO and EGM, el9 is required for the whole of it --
an earlier version of this docstring claimed none of these JSONs parse on el8, which is
too broad.

  /cvmfs/cms.cern.ch/common/cmssw-el9 -- bash -c '
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src && eval `scramv1 runtime -sh`
    python3 analysis/macros/check_corrections_signal.py DsW_preEE
  '

postEE / W+jets / data are deliberately NOT included yet (different era JSON tag, WPs,
GT) -- preEE DsW/DsstW only, per plan.
"""
import glob
import math
import os
import sys

import awkward as ak
import correctionlib
import numpy as np
import ROOT
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jme_corrections as jme  # noqa: E402  (JEC/JER via the official tutorial modules)

FILE_GROUPS = {
    "DsW_preEE": sorted(glob.glob(
        "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/Signal_DsW_H_phipimu/output_*.root")),
    "DsstW_preEE": sorted(glob.glob(
        "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/Signal_DsstW_H_phipimu/output_*.root")),
    "WtoLNu2Jets_0J_preEE": sorted(glob.glob(
        "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/WtoLNu2Jets_0J_H_phipimu_GenTag/output_*.root")),
    "WtoLNu2Jets_1J_preEE": sorted(glob.glob(
        "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/WtoLNu2Jets_1J_H_phipimu_GenTag/output_*.root")),
    "WtoLNu2Jets_2J_preEE": sorted(glob.glob(
        "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/WtoLNu2Jets_2J_H_phipimu_GenTag/output_*.root")),
}
TREE = "Events"

# Jet-bin cross sections (pb) and preEE luminosity (pb^-1): same values already
# established in analysis/macros/make_samples.py (XSEC_JETBIN / ERAS[1]["lumi_pb"]) --
# reused here, not re-derived. amcatnloFXFX W+jets carries genWeight = +/-|w| (a sizable
# negative fraction, ~10-35% depending on jet bin per memory); Lumi_reweight_i =
# genWeight_i * (xsec * L / sum(genWeight)) folds in both the per-event sign and the
# per-bin normalization, and IS the weight used below for W+jets samples (not raw
# genWeight) -- per-bin correction means are unaffected (Lumi_reweight is genWeight times
# one constant per bin, so a weighted mean/std within one bin is invariant to it), but it
# is the physically correct weight to report, and it is what makes combining the three
# jet bins into one W+jets picture afterward valid.
XSEC_JETBIN_PB = {0: 51846.0, 1: 8860.0, 2: 3284.0}
LUMI_PREEE_PB = 7990.0

POG_DIR = "/cvmfs/cms-griddata.cern.ch/cat/metadata"
JSON_LUM = f"{POG_DIR}/LUM/Run3-22CDSep23-Summer22-NanoAODv12/2024-01-31/puWeights.json.gz"
JSON_MUO_Z = f"{POG_DIR}/MUO/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-18/muon_Z.json.gz"
JSON_MUO_SS = f"{POG_DIR}/MUO/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-18/muon_scalesmearing.json.gz"
# JME (JEC, JER, jet veto map, Type-1 MET) JSONs and tags are resolved by the official
# tutorial config through jme_corrections.py, not declared here.
JSON_BTV = f"{POG_DIR}/BTV/Run3-22CDSep23-Summer22-NanoAODv12/2025-08-20/btagging.json.gz"
JSON_EGM = f"{POG_DIR}/EGM/Run3-22CDSep23-Summer22-NanoAODv12/2025-12-15/electron.json.gz"
# Per-sample efficiency maps from build_btag_efficiency_map.py -- each sample is
# corrected with ITS OWN tagging efficiency (see that script's docstring). An earlier
# version applied one W+jets map to every sample, signal included.
_BTAG_EFF_DIR = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/btag_efficiency"
BTAG_EFF_MAP_JSON = {
    "DsW_preEE": f"{_BTAG_EFF_DIR}/DsW_preEE.json",
    "DsstW_preEE": f"{_BTAG_EFF_DIR}/DsstW_preEE.json",
    "WtoLNu2Jets_0J_preEE": f"{_BTAG_EFF_DIR}/WtoLNu2Jets_preEE.json",
    "WtoLNu2Jets_1J_preEE": f"{_BTAG_EFF_DIR}/WtoLNu2Jets_preEE.json",
    "WtoLNu2Jets_2J_preEE": f"{_BTAG_EFF_DIR}/WtoLNu2Jets_preEE.json",
}
# BTV SFs are measured for 20 <= pt < 1000 only; outside that range particleNet_comb /
# particleNet_light silently return a flow constant instead of raising (measured at
# |eta|=1, T WP: comb 0.837 at pt 19.9 vs 1.011 at pt 20). The b-veto is therefore
# defined on jets with corrected pt >= 20 -- the SF key fixes the cut, as for the muon.
BTAG_JET_PT_MIN = 20.0

PU_KEY = "Collisions2022_355100_357900_eraBCD_GoldenJson"
MUO_ID_KEY = "NUM_TightID_DEN_TrackerMuons"
MUO_ISO_KEY = "NUM_TightPFIso_DEN_TightID"
MUO_TRIG_KEY = "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight"
EGM_YEAR = "2022Re-recoBCD"

BTAG_WP_TIGHT = 0.6734  # producer's hardcoded value; cross-checked below against the JSON
# PFIsoTight = relIso04(dBeta) < 0.15. Not a free choice: with tight ID the only trigger-SF
# DEN available is CutBasedIdTight_and_PFIsoTight, which fixes the iso WP (see section 3b).
ISO_WP_TIGHT = 0.15
NULL = -777

BRANCHES = [
    "mu_pt", "mu_eta", "mu_phi", "mu_charge", "mu_IsoR04_RelIso",
    "n_b_jets_tight_dR04", "n_b_jets_tight_pt20_dR04", "Ds_eta", "Ds_phi",
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_tightId", "Muon_nTrackerLayers",
    "Muon_pfIsoId", "num_os_mu_tight", "num_ss_mu_tight",
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_btagPNetB", "Jet_hadronFlavour",
    "Electron_pt", "Electron_eta", "Electron_deltaEtaSC",
    "Pileup_nTrueInt", "Rho_fixedGridRhoFastjetAll", "genWeight",
    "event", "luminosityBlock",
    # stored producer jet quantities, for the jet-variable closure self-check
    "num_jets", "num_jets_pt20", "num_jets_pt30", "jet_ht", "jet_ht_pt20", "jet_leadpt",
] + jme.JET_BRANCHES + jme.JVM_BRANCHES + jme.MET_BRANCHES + ["PuppiMET_pt"]

MUO_LOW_PT_THRESHOLD = 26.0  # matches cms-cat/nanoAOD-tools-modules MuonScaRe default


def wmean(x, w):
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    return float(np.average(x, weights=w)), float(np.sqrt(np.average((x - np.average(x, weights=w)) ** 2, weights=w)))


def delta_phi(phi1, phi2):
    d = phi1 - phi2
    return np.where(d > np.pi, d - 2 * np.pi, np.where(d < -np.pi, d + 2 * np.pi, d))


# ---- MUO momentum scale/smearing: faithful port of the OFFICIAL recipe ----
# Source: cms-cat/nanoAOD-tools-modules src/MuonScaRe.cc + interface/MuonScaRe.h +
# python/modules/muonScaleRes.py, itself imported from (per that file's own header
# comment) gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/773fd624/examples/
# MuonScaRe.cc -- i.e. the same repo that publishes the muon_scalesmearing.json this
# script reads. This is a line-by-line port, not a reinterpretation: pt_scale here
# reproduces `return 1./(m/pt + charge*a);` and pt_resol reproduces
# `ptc = pt*(1 + k*std*rndm)` with k/std/rndm built exactly as in get_k/get_std/
# get_rndm there (CrystalBall inverse-CDF sampled with a TRandom3 seeded via the same
# custom SeedSequence mixing of (event, luminosityBlock, phi)).
import scipy.special


class _CrystalBall:
    def __init__(self, mean, sigma, alpha, n):
        self.m, self.s, self.a, self.n = mean, sigma, alpha, n
        pi = 3.14159
        self.sqrtPiOver2 = math.sqrt(pi / 2.0)
        self.sqrt2 = math.sqrt(2.0)
        fa = abs(alpha)
        ex = math.exp(-fa * fa / 2)
        A = (n / fa) ** n * ex
        C1 = n / fa / (n - 1) * ex
        D1 = 2 * self.sqrtPiOver2 * math.erf(fa / self.sqrt2)
        self.B = n / fa - fa
        self.C = (D1 + 2 * C1) / C1
        self.D = (D1 + 2 * C1) / 2
        self.N = 1.0 / sigma / (D1 + 2 * C1)
        self.k = 1.0 / (n - 1)
        self.NA = self.N * A
        self.Ns = self.N * sigma
        self.NC = self.Ns * C1
        self.F = 1 - fa * fa / n
        self.G = sigma * n / fa
        self.cdfMa = self._cdf(mean - alpha * sigma)
        self.cdfPa = self._cdf(mean + alpha * sigma)

    def _cdf(self, x):
        d = (x - self.m) / self.s
        if d < -self.a:
            return self.NC / (self.F - self.s * d / self.G) ** (self.n - 1)
        if d > self.a:
            return self.NC * (self.C - (self.F + self.s * d / self.G) ** (1 - self.n))
        return self.Ns * (self.D - self.sqrtPiOver2 * math.erf(-d / self.sqrt2))

    def invcdf(self, u):
        if u < self.cdfMa:
            return self.m + self.G * (self.F - (self.NC / u) ** self.k)
        if u > self.cdfPa:
            return self.m - self.G * (self.F - (self.C - u / self.NC) ** (-self.k))
        return self.m - self.sqrt2 * self.s * scipy.special.erfinv((self.D - u / self.Ns) / self.sqrtPiOver2)


def _seed_sequence_one(seeds):
    """Port of the custom SeedSequence::generate() for n=1 output word."""
    mult = 0x9E3779B9
    mix_const = 0x85EBCA6B
    mask = 0xFFFFFFFF
    buf = 0x8B8B8B8B
    # n=1, s=len(seeds)=3 -> the "for i in range(min(n,s))" loop runs once (i=0)
    buf = (buf ^ ((seeds[0] + mult * 0) & mask)) & mask
    # mixing loop, single element: z = buf ^ (buf >> 27); buf = (z*mix_const) ^ (buf<<13)
    z = (buf ^ (buf >> 27)) & mask
    buf = ((z * mix_const) & mask) ^ ((buf << 13) & mask)
    return buf & mask


def muo_get_rndm(cset, eta, phi, nL, event_number, lumi_number):
    mean = cset["cb_params"].evaluate(abs(eta), nL, 0)
    sigma = cset["cb_params"].evaluate(abs(eta), nL, 1)
    n = cset["cb_params"].evaluate(abs(eta), nL, 2)
    alpha = cset["cb_params"].evaluate(abs(eta), nL, 3)
    cb = _CrystalBall(mean, sigma, alpha, n)
    phi_seed = int((phi / math.pi) * ((1 << 31) - 1)) & 0xFFF
    seed = _seed_sequence_one([
        int(event_number) & 0xFFFFFFFF,
        int(lumi_number) & 0xFFFFFFFF,
        phi_seed & 0xFFFFFFFF,
    ])
    rnd = ROOT.TRandom3()
    rnd.SetSeed(seed)
    return cb.invcdf(rnd.Rndm())


def muo_get_std(cset, pt, nL, abseta):
    p0 = cset["poly_params"].evaluate(abseta, nL, 0)
    p1 = cset["poly_params"].evaluate(abseta, nL, 1)
    p2 = cset["poly_params"].evaluate(abseta, nL, 2)
    sigma = p0 + p1 * pt + p2 * pt * pt
    return max(sigma, 0.0)


def muo_get_k(cset, abseta, var="nom"):
    k_data = cset["k_data"].evaluate(abseta, var)
    k_mc = cset["k_mc"].evaluate(abseta, var)
    return math.sqrt(k_data * k_data - k_mc * k_mc) if k_mc < k_data else 0.0


def muo_pt_scale(cset, is_data, pt, eta, phi, charge):
    a = cset["a_data" if is_data else "a_mc"].evaluate(eta, phi, "nom")
    m = cset["m_data" if is_data else "m_mc"].evaluate(eta, phi, "nom")
    if pt < MUO_LOW_PT_THRESHOLD:
        return pt
    return 1.0 / (m / pt + charge * a)


def muo_pt_resol(cset, pt, eta, phi, nL, event_number, lumi_number):
    # The official C++ (MuonScaRe.cc) has no guard against degenerate cb_params (e.g.
    # alpha=0 for some (abseta, nTrackerLayers) bin) -- it doesn't need one, because C++
    # float division by zero silently yields inf/nan, which the `isnan(ptc)` line below
    # already catches. Python's plain `/` raises ZeroDivisionError instead of returning
    # nan, so without this try/except a muon landing in such a bin crashes the whole job
    # instead of falling back to the unsmeared pt like the reference implementation does.
    # Seen in practice: WtoLNu2Jets_2J_preEE hit this on a real muon (not a signal muon,
    # signal never triggered it in 473358 rows).
    try:
        rndm = muo_get_rndm(cset, eta, phi, nL, event_number, lumi_number)
        std = muo_get_std(cset, pt, nL, abs(eta))
        k = muo_get_k(cset, abs(eta))
        ptc = pt * (1 + k * std * rndm)
    except (ZeroDivisionError, ValueError, OverflowError):
        ptc = float("nan")
    if math.isnan(ptc):
        ptc = pt
    if ptc / pt > 2 or ptc / pt < 0.1 or ptc < 0 or pt < MUO_LOW_PT_THRESHOLD or pt > 200:
        ptc = pt
    return ptc


def muo_pt_resol_var(cset, pt_woresol, pt_wresol, eta, updn):
    """Smearing systematic. Ported from MuonScaRe.cc::pt_resol_var / the official
    MuonScaRe.py::pt_resol_var (cross-checked, both agree on this formula)."""
    k = muo_get_k(cset, abs(eta))
    if k == 0:
        return pt_wresol
    k_unc = cset["k_mc"].evaluate(abs(eta), "stat")
    std_x_rndm = (pt_wresol / pt_woresol - 1) / k
    if updn == "up":
        pt_var = pt_woresol * (1 + (k + k_unc) * std_x_rndm)
    elif updn == "dn":
        pt_var = pt_woresol * (1 + (k - k_unc) * std_x_rndm)
    else:
        raise ValueError("updn must be 'up' or 'dn'")
    if pt_var / pt_woresol > 2 or pt_var / pt_woresol < 0.1 or pt_var < 0:
        pt_var = pt_woresol
    return pt_var


def muo_pt_scale_var(cset, pt, eta, phi, charge, updn):
    """Scale systematic. Ported from MuonScaRe.cc::pt_scale_var / the official
    MuonScaRe.py::pt_scale_var (cross-checked, both agree on this formula). Per the
    official calling convention (muonScaleRes.py::getPtVarScale), `pt` here is the
    FINAL corrected pt (after scale+smear), not the raw pt."""
    stat_a = cset["a_mc"].evaluate(eta, phi, "stat")
    stat_m = cset["m_mc"].evaluate(eta, phi, "stat")
    stat_rho = cset["m_mc"].evaluate(eta, phi, "rho_stat")
    unc = pt * pt * math.sqrt(stat_m * stat_m / (pt * pt) + stat_a * stat_a
                               + 2 * charge * stat_rho * stat_m / pt * stat_a)
    if updn == "up":
        return pt + unc
    elif updn == "dn":
        return pt - unc
    raise ValueError("updn must be 'up' or 'dn'")


def check_muon_reid(a):
    """Locate the selected muon inside standard Muon_* by (pt>26 & tightId, argmax pt).
    Returns (idx array, match bool array)."""
    n = len(a)
    idx = np.full(n, -1, dtype=np.int64)
    match = np.zeros(n, dtype=bool)
    for i in range(n):
        mpt = np.asarray(a["Muon_pt"][i])
        tid = np.asarray(a["Muon_tightId"][i])
        meta = np.asarray(a["Muon_eta"][i])
        mch = np.asarray(a["Muon_charge"][i])
        sel = (mpt > 26) & tid
        if not sel.any():
            continue
        j = int(np.argmax(np.where(sel, mpt, -1)))
        idx[i] = j
        match[i] = (
            abs(mpt[j] - a["mu_pt"][i]) < 1e-3
            and abs(meta[j] - a["mu_eta"][i]) < 1e-3
            and mch[j] == a["mu_charge"][i]
        )
    return idx, match


def reco_examined_jet_mask(a, i):
    """Boolean mask, within event i's Jet_* array, of the jets the b-veto EXAMINES
    (central, tight ID, dR(Ds,jet)>0.4) -- no tag requirement. This is the population
    the fixed-WP b-tag SF product must run over."""
    jeta = np.asarray(a["Jet_eta"][i])
    jphi = np.asarray(a["Jet_phi"][i])
    jid = np.asarray(a["Jet_jetId"][i])
    dphi = delta_phi(jphi, a["Ds_phi"][i])
    dR = np.sqrt((jeta - a["Ds_eta"][i]) ** 2 + dphi ** 2)
    return (np.abs(jeta) < 2.5) & (jid >= 2) & (dR > 0.4)


def reco_tagged_jet_mask(a, i):
    """Same as reco_examined_jet_mask, plus the tight-WP tag requirement -- this is
    exactly what n_b_jets_tight_dR04 counts, used only for the closure self-check."""
    jb = np.asarray(a["Jet_btagPNetB"][i])
    return reco_examined_jet_mask(a, i) & (jb > BTAG_WP_TIGHT)


# Standard nano stores several Jet_* branches with lossy, reduced-precision ROOT floats:
# Jet_btagPNetB precision=10 (jetsAK4_Puppi_cff.py), Jet_eta/Jet_phi precision=12 (the
# shared CandVars template in common_cff.py). n_b_jets_tight_dR04 was computed by the
# producer from the FULL-precision PAT-level values, before any of that nano compression,
# so a jet whose true value sits just inside a cut can get rounded to just outside it (or
# vice versa) in the standard branches -- confirmed on two DIFFERENT variables so far:
#   - WtoLNu2Jets_0J row 20753/65169: btag stored as 0.67333984 vs WP 0.6734 (~6e-5 gap,
#     inside the ~5e-4 step for precision=10 at this magnitude)
#   - WtoLNu2Jets_1J row 5023: eta stored as EXACTLY 2.500000 (the |eta|<2.5 cut boundary
#     itself), true value was presumably just under 2.5 before precision=12 rounding
# This is not a bug in the producer or this script -- it is real, understood information
# loss in the standard branches. It never happened on signal (0/473358 rows) because that
# sample's jets rarely sit this close to any of the three cut boundaries. Tolerances below
# are a >=2x margin on each variable's quantization step at the relevant magnitude.
BOUNDARY_TOL = dict(eta=0.001, btag=0.001, dR=0.001)


def check_btag_closure(a):
    """Returns (exact_match, explained_mismatch) bool arrays. A mismatch counts as
    "explained" only if flipping the pass/fail status of some near-boundary jet(s) --
    within BOUNDARY_TOL of the |eta|<2.5, btag>WP, or dR>0.4 cut -- could account for the
    count difference. An unexplained mismatch is a real bug, not precision noise, and
    must still fail loudly."""
    n = len(a)
    exact = np.zeros(n, dtype=bool)
    explained = np.zeros(n, dtype=bool)
    for i in range(n):
        n_reco = int(reco_tagged_jet_mask(a, i).sum())
        n_stored = int(a["n_b_jets_tight_dR04"][i])
        if n_reco == n_stored:
            exact[i] = True
            continue
        jeta = np.asarray(a["Jet_eta"][i])
        jphi = np.asarray(a["Jet_phi"][i])
        jid = np.asarray(a["Jet_jetId"][i])
        jb = np.asarray(a["Jet_btagPNetB"][i])
        dphi = delta_phi(jphi, a["Ds_phi"][i])
        dR = np.sqrt((jeta - a["Ds_eta"][i]) ** 2 + dphi ** 2)
        # jetId>=2 is required outright (an int branch, not precision-truncated). For a
        # near-boundary jet's flip to actually change the tagged COUNT, every OTHER cut
        # must already pass -- e.g. a jet near the eta edge only matters if it also
        # clears the btag and dR cuts (row 5023: eta=2.500000 exactly, btag=0.91, dR=2.49
        # -- only eta was ambiguous, the other two were nowhere near their boundaries).
        base = jid >= 2
        passes_eta = np.abs(jeta) < 2.5
        passes_btag = jb > BTAG_WP_TIGHT
        passes_dR = dR > 0.4
        near_eta = base & passes_btag & passes_dR & (np.abs(np.abs(jeta) - 2.5) < BOUNDARY_TOL["eta"])
        near_btag = base & passes_eta & passes_dR & (np.abs(jb - BTAG_WP_TIGHT) < BOUNDARY_TOL["btag"])
        near_dR = base & passes_eta & passes_btag & (np.abs(dR - 0.4) < BOUNDARY_TOL["dR"])
        n_near = int((near_eta | near_btag | near_dR).sum())
        if n_near >= abs(n_reco - n_stored):
            explained[i] = True
    return exact, explained


def check_jet_variable_closure(a):
    """Which stored jet quantities can be rebuilt from nano `Jet_*`, and which cannot.

    Every jet loop in HDsPhiPiMuProducer.cc cuts only on userInt("tightId")!=0 (== nano
    `Jet_jetId` bit 2) and, in the candidate loop, |eta|<2.5 -- there is NO pt floor
    anywhere. nano's jetTable stores finalJetsPuppi = pt>15 GeV. So the producer counted
    jets nano does not contain, and the pt-thresholdless quantities cannot close.

    Reports rather than asserts: the failures here are a known, quantified property of
    the NanoAOD, not a bug to halt on. What matters downstream is only that the
    JME-dependent quantities we actually recompute (the pt>20 family, jet_leadpt, MET)
    do close -- those are asserted by the caller.
    """
    tight = (a["Jet_jetId"] & 2) == 2
    pt = a["Jet_pt"]
    n = len(a)
    nano_njet = ak.to_numpy(ak.sum(tight, axis=1)).astype(int)
    deficit = ak.to_numpy(a["num_jets"]).astype(int) - nano_njet

    def frac_bad(mine, stored, rtol=None):
        m = ak.to_numpy(mine).astype(float)
        s = ak.to_numpy(stored).astype(float)
        bad = (np.abs(m - s) > rtol * np.maximum(np.abs(s), 1e-9)) if rtol else (m != s)
        return float(bad.mean())

    out = {
        "num_jets_pt20": frac_bad(ak.sum(tight & (pt > 20), axis=1), a["num_jets_pt20"]),
        "num_jets_pt30": frac_bad(ak.sum(tight & (pt > 30), axis=1), a["num_jets_pt30"]),
        "jet_ht_pt20": frac_bad(ak.sum(ak.where(tight & (pt > 20), pt, 0.0), axis=1),
                                a["jet_ht_pt20"], rtol=2e-3),
        "num_jets": frac_bad(nano_njet, a["num_jets"]),
        "jet_ht": frac_bad(ak.sum(ak.where(tight, pt, 0.0), axis=1), a["jet_ht"], rtol=2e-3),
    }
    lead = ak.where(ak.sum(tight, axis=1) > 0,
                    ak.fill_none(ak.max(ak.where(tight, pt, -1.0), axis=1), -1.0), float(NULL))
    out["jet_leadpt"] = frac_bad(lead, a["jet_leadpt"], rtol=2e-3)

    ht_stored = ak.to_numpy(a["jet_ht"]).astype(float)
    ht_nano = ak.to_numpy(ak.sum(ak.where(tight, pt, 0.0), axis=1)).astype(float)
    ok = ht_stored > 0
    missing = dict(
        n_evt_frac=float((deficit > 0).mean()),
        mean_njet=float(deficit.mean()),
        mean_ht=float((ht_stored - ht_nano).mean()),
        frac_of_ht=float(np.mean((ht_stored - ht_nano)[ok] / ht_stored[ok])),
    )
    return out, missing, n


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in FILE_GROUPS:
        sys.exit(f"usage: {sys.argv[0]} <{'|'.join(FILE_GROUPS)}>")
    sample = sys.argv[1]
    files = FILE_GROUPS[sample]

    print(f"Reading {len(files)} files for {sample} : {TREE}")
    a = uproot.concatenate([f"{f}:{TREE}" for f in files], filter_name=BRANCHES)
    n = len(a)
    gen_weight = np.asarray(a["genWeight"])
    print(f"{sample}: n rows = {n}")

    if sample.startswith("WtoLNu2Jets_"):
        jet_bin = int(sample.split("_")[1][0])  # "0J"/"1J"/"2J" -> 0/1/2
        sum_gen_weight = float(gen_weight.sum())
        neg_frac = float((gen_weight < 0).mean())
        norm = XSEC_JETBIN_PB[jet_bin] * LUMI_PREEE_PB / sum_gen_weight
        genw = gen_weight * norm
        print(f"[W+jets weighting] jet_bin={jet_bin} xsec={XSEC_JETBIN_PB[jet_bin]} pb,"
              f" L={LUMI_PREEE_PB} pb^-1, sum(genWeight)={sum_gen_weight:.6g},"
              f" norm={norm:.6g}, neg-weight frac={neg_frac:.4f}")
        print(f"  sum(Lumi_reweight) for this bin = {genw.sum():.6g}  (needed to combine"
              f" the three jet bins afterward: weight each bin's sums by this, not by row count)")
    else:
        genw = gen_weight

    # ---- self-checks (must pass before trusting anything below) ----
    idx_mu, muon_match = check_muon_reid(a)
    if not muon_match.all():
        bad = np.where(~muon_match)[0]
        raise SystemExit(f"muon re-id closure FAILED on {len(bad)}/{n} rows, e.g. {bad[:5]}")
    print(f"[self-check] muon re-id: {muon_match.sum()}/{n} OK")

    btag_exact, btag_explained = check_btag_closure(a)
    btag_match = btag_exact | btag_explained
    if not btag_match.all():
        bad = np.where(~btag_match)[0]
        raise SystemExit(f"b-veto jet-set closure FAILED (unexplained) on {len(bad)}/{n} rows,"
                          f" e.g. {bad[:5]} -- this is NOT the known nano precision-quantization"
                          f" effect (see BOUNDARY_TOL comment), investigate before proceeding")
    print(f"[self-check] b-veto jet-set closure: {btag_exact.sum()}/{n} exact,"
          f" {btag_explained.sum()}/{n} explained by nano branch precision quantization"
          f" (Jet_eta/Jet_phi precision=12 or Jet_btagPNetB precision=10) near a cut"
          f" boundary (see comment on BOUNDARY_TOL) -- 0 unexplained")

    jetvar, missing, _ = check_jet_variable_closure(a)
    print("[self-check] jet-variable closure vs nano Jet_* (mismatch fraction):")
    for k in ("num_jets_pt20", "num_jets_pt30", "jet_ht_pt20", "jet_leadpt"):
        print(f"    {k:16s} {jetvar[k]*100:7.3f}%   <- JME-dependent, must be rebuildable")
    for k in ("num_jets", "jet_ht"):
        print(f"    {k:16s} {jetvar[k]*100:7.3f}%   <- no pt threshold in the producer")
    print(f"    events where the producer saw a tightId jet nano does not store:"
          f" {missing['n_evt_frac']*100:.1f}%  (mean {missing['mean_njet']:.2f} jets,"
          f" {missing['mean_ht']:.1f} GeV = {missing['frac_of_ht']*100:.1f}% of stored jet_ht)")
    for k in ("num_jets_pt20", "num_jets_pt30", "jet_ht_pt20", "jet_leadpt"):
        if jetvar[k] > 0.02:
            raise SystemExit(f"jet-variable closure FAILED: {k} mismatches {jetvar[k]*100:.2f}%"
                             f" of rows (>2%) -- this quantity has a pt threshold above nano's"
                             f" 15 GeV storage cut and should rebuild almost exactly; investigate"
                             f" before trusting any JME-corrected quantity below")
    print("    -> jet_ht and num_jets are NOT recoverable (soft jets nano dropped); they are")
    print("       the only JME-dependent DNN feature (jet_ht) / a JME-invariant count (num_jets).")

    iso_fail_frac = float((np.asarray(a["mu_IsoR04_RelIso"]) >= ISO_WP_TIGHT).mean())
    bveto = np.asarray(a["n_b_jets_tight_dR04"]) == 0
    iso_fail_frac_bveto = float((np.asarray(a["mu_IsoR04_RelIso"])[bveto] >= ISO_WP_TIGHT).mean())
    print(f"[note] frac(mu_IsoR04_RelIso >= {ISO_WP_TIGHT}): all rows={iso_fail_frac:.4f}"
          f"  b-veto rows={iso_fail_frac_bveto:.4f}  (PFIsoTight, fixed by the trigger-SF DEN;"
          f" applied offline in section 3b together with |eta|<2.4 and corrected pt>26)")

    # ---- 1. LUM: pileup ----
    cset = correctionlib.CorrectionSet.from_file(JSON_LUM)
    c = cset[PU_KEY]
    nti_raw = np.asarray(a["Pileup_nTrueInt"])
    nti = np.clip(nti_raw, 0, 98)  # guard against out-of-range
    w_pu = {v: np.array([c.evaluate(float(x), v) for x in nti]) for v in ("nominal", "up", "down")}
    print("\n[LUM pileup]")
    print(f"  nTrueInt in this file: min={nti_raw.min():.2f} max={nti_raw.max():.2f}"
          f" mean={nti_raw.mean():.2f}")
    for v in ("nominal", "up", "down"):
        m, s = wmean(w_pu[v], genw)
        print(f"  w_pu[{v}]: mean={m:.4f} std={s:.4f}")
    top = np.argsort(w_pu["nominal"])[-3:]
    print(f"  largest nominal weights come from low-nTrueInt outliers (few MC events there):"
          f" {[(round(float(nti_raw[i]), 1), round(float(w_pu['nominal'][i]), 2)) for i in top]}"
          f"  -- the std is dominated by this tail, not a flat per-event spread")

    # Does our private MC's PU profile match the MC denominator puWeights.json assumes?
    # <w_pu> ~ 1 is only weak evidence: the weights are BUILT so their mean sits near 1
    # over the assumed profile, so a mean near 1 is close to unfalsifiable. The
    # discriminating quantity is the SHAPE -- where our events sit relative to where the
    # weight deviates from 1, and how much of the total weight comes from each region.
    # If our profile matched the assumed one, the weight would be distributed around 1
    # with our population concentrated where w ~ 1; a profile offset shows up as our bulk
    # sitting systematically on one side.
    print("  [profile check] our nTrueInt distribution vs where puWeights deviates from 1:")
    edges = np.array([0, 10, 20, 30, 40, 50, 60, 70, 99], dtype=float)
    tot_w = float(np.abs(genw).sum())
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (nti_raw >= lo) & (nti_raw < hi)
        if not m.any():
            continue
        frac = float(np.abs(genw[m]).sum() / tot_w)
        wm, _ = wmean(w_pu["nominal"][m], genw[m])
        print(f"    nTrueInt [{lo:4.0f},{hi:4.0f}): {frac*100:5.2f}% of events, "
              f"<w_pu>={wm:6.3f}")
    bulk = (nti_raw >= 20) & (nti_raw < 70)
    wb, _ = wmean(w_pu["nominal"][bulk], genw[bulk])
    print(f"    bulk (20<=nTrueInt<70, {float(np.abs(genw[bulk]).sum()/tot_w)*100:.1f}% of"
          f" events): <w_pu>={wb:.4f}")
    # An earlier version printed an interpretation here ("bulk <w_pu> near 1 = healthy")
    # that does not hold: w_pu IS the data/MC profile ratio, so per-bin structure far from
    # 1 is exactly what it should show, and a mean near 1 only says the normalisation is
    # roughly right. It cannot tell whether our private MC was generated with the same
    # profile puWeights.json divides by. On top of that, these rows are AFTER the event
    # selection, and selection is known to shift pileup-correlated quantities (memory:
    # feedback_pu_reweight_denominator), so no test on these rows is unbiased. The
    # table above is descriptive only.
    print("    Descriptive only: w_pu is a data/MC ratio, so per-bin structure is expected;")
    print("    rows are post-selection, so this cannot validate the generation profile.")

    # ---- 2. MUO: inputs. The SFs are NOT evaluated here -- see section 3b. ----
    # The muon SFs are functions of the muon pt, and the pt we are entitled to hand them
    # is the CORRECTED one, which does not exist until section 3 has run. Evaluating them
    # here on the stored pt (as an earlier version did) meant every SF was looked up in
    # the wrong pt bin for the rows near a bin edge, and -- worse -- it hid that the
    # trigger SF is only defined from pt=26 upwards.
    mu_pt_raw = np.asarray(a["mu_pt"])
    mu_eta_signed = np.asarray(a["mu_eta"])  # SIGNED: see section 3b, this matters

    # ---- 3. MUO: momentum scale + smearing ----
    # Combination formula source: cms-cat/nanoAOD-tools-modules src/MuonScaRe.cc (imported
    # from gitlab.cern.ch/cms-nanoAOD/jsonpog-integration itself) -- see the muo_pt_scale /
    # muo_pt_resol functions above for the faithful port. is_data=False (signal MC): scale
    # uses m_mc/a_mc, smearing is applied on top using k_mc/k_data/cb_params/poly_params +
    # a deterministic (event, luminosityBlock, phi)-seeded random draw.
    cset = correctionlib.CorrectionSet.from_file(JSON_MUO_SS)
    mu_phi = np.asarray(a["mu_phi"])
    mu_eta_clip = np.clip(np.asarray(a["mu_eta"]), -2.399, 2.399)
    mu_charge_arr = np.asarray(a["mu_charge"])
    event_arr = np.asarray(a["event"])
    lumi_arr = np.asarray(a["luminosityBlock"])

    pt_after_scale = np.empty(n)
    pt_after_smear = np.empty(n)
    pt_scale_up = np.empty(n)
    pt_scale_dn = np.empty(n)
    pt_smear_up = np.empty(n)
    pt_smear_dn = np.empty(n)
    for i in range(n):
        nL = float(np.asarray(a["Muon_nTrackerLayers"][i])[idx_mu[i]])
        eta_i, phi_i, charge_i = float(mu_eta_clip[i]), float(mu_phi[i]), int(mu_charge_arr[i])
        pt_after_scale[i] = muo_pt_scale(cset, False, mu_pt_raw[i], eta_i, phi_i, charge_i)
        pt_after_smear[i] = muo_pt_resol(cset, pt_after_scale[i], eta_i, phi_i, nL,
                                          event_arr[i], lumi_arr[i])
        # official calling convention (muonScaleRes.py): scale variations are built
        # around the FINAL (scale+smear) pt; smear variations compare scale-only vs
        # final pt. Same try/except rationale as muo_pt_resol -- Python raises where
        # C++ would silently give nan/inf on a degenerate parameter.
        try:
            pt_scale_up[i] = muo_pt_scale_var(cset, pt_after_smear[i], eta_i, phi_i, charge_i, "up")
            pt_scale_dn[i] = muo_pt_scale_var(cset, pt_after_smear[i], eta_i, phi_i, charge_i, "dn")
        except (ZeroDivisionError, ValueError, OverflowError):
            pt_scale_up[i] = pt_scale_dn[i] = pt_after_smear[i]
        try:
            pt_smear_up[i] = muo_pt_resol_var(cset, pt_after_scale[i], pt_after_smear[i], eta_i, "up")
            pt_smear_dn[i] = muo_pt_resol_var(cset, pt_after_scale[i], pt_after_smear[i], eta_i, "dn")
        except (ZeroDivisionError, ValueError, OverflowError):
            pt_smear_up[i] = pt_smear_dn[i] = pt_after_smear[i]
    scale_shift = pt_after_scale - mu_pt_raw
    total_shift = pt_after_smear - mu_pt_raw
    print("\n[MUO momentum scale/smearing]  (formula source: cms-cat/nanoAOD-tools-modules"
          " MuonScaRe.cc/.py, see comment above)")
    print(f"  scale only:  mean shift={scale_shift.mean():.5f} GeV std={scale_shift.std():.5f} GeV"
          f"  ({100 * scale_shift.mean() / mu_pt_raw.mean():.4f}% of <mu_pt>)")
    print(f"  scale+smear: mean shift={total_shift.mean():.5f} GeV std={total_shift.std():.5f} GeV"
          f"  ({100 * total_shift.mean() / mu_pt_raw.mean():.4f}% of <mu_pt>)")
    n_reverted = int(((pt_after_smear == pt_after_scale) & (mu_pt_raw >= MUO_LOW_PT_THRESHOLD)).sum())
    print(f"  rows where smearing hit a guard rail and reverted to the scale-only value:"
          f" {n_reverted}/{n}")
    for label, up, dn in (("scaleUnc", pt_scale_up, pt_scale_dn), ("smearUnc", pt_smear_up, pt_smear_dn)):
        du = (up - pt_after_smear)
        dd = (dn - pt_after_smear)
        print(f"  {label}: up mean={du.mean():.5f} GeV std={du.std():.5f} GeV,"
              f" down mean={dd.mean():.5f} GeV std={dd.std():.5f} GeV")

    # ---- 3b. MUO: the analysis muon selection, then ID x ISO x trigger SF ----
    # The selection below is not a free choice -- every piece of it is pinned by what the
    # SF keys will accept, measured against the JSON itself rather than assumed:
    #
    #   pt > 26   The trigger SF's pt binning is [26,30,40,50,60,120,inf] with
    #             flow="error". Measured: evaluate(eta, 25.9, ...) raises
    #             "Index below bounds", it does NOT clamp. 26 is the SF's own floor.
    #             Applied to the CORRECTED pt, because that is what we pass to evaluate().
    #   |eta|<2.4 The eta binning is [-2.4,...,2.4], also flow="error".
    #   PFIsoTight  Every trigger key's DEN carries an iso requirement; there is no
    #             "tight ID, no iso" denominator in this JSON at all. With tight ID the
    #             only consistent DEN is CutBasedIdTight_and_PFIsoTight, which fixes the
    #             iso WP to Tight = relIso04(dBeta) < 0.15. Measured: nano's official
    #             Muon_pfIsoId>=4 and mu_IsoR04_RelIso<0.15 agree on 3334/3334 rows, and
    #             mu_IsoR04_RelIso equals nano Muon_pfRelIso04_all to 1.2e-7, so the
    #             producer's hand-rolled iso IS the POG's quantity.
    #
    # SIGNED eta, and only for the trigger key it actually matters: measured, ID and ISO
    # are exactly symmetric under eta -> -eta (0.000% over a scan of |eta| in 0.5-2.3 and
    # pt in 27-100), but the TRIGGER SF is not -- 17.8% between +2.3 and -2.3 at pt=27,
    # and still 1.3-1.5% at |eta| ~ 1-2 near threshold. An earlier version passed
    # np.abs(mu_eta) to all three keys, which left the trigger SF wrong by up to 17.8% in
    # the forward region. Pass signed eta.
    mu_pt_corr = pt_after_smear
    sel_eta = np.abs(mu_eta_signed) < 2.4
    sel_pt = mu_pt_corr > 26.0
    sel_iso = np.asarray(a["mu_IsoR04_RelIso"]) < ISO_WP_TIGHT
    sel_mu = sel_eta & sel_pt & sel_iso
    print("\n[MUO analysis muon selection]  (each cut set by the SF keys, see comment)")
    print(f"  |eta|<2.4          : {sel_eta.sum():6d}/{n} ({100*sel_eta.mean():6.2f}%)")
    print(f"  corrected pt>26    : {sel_pt.sum():6d}/{n} ({100*sel_pt.mean():6.2f}%)"
          f"   <- the trigger SF's own floor, applied to the corrected pt")
    print(f"  PFIsoTight (<{ISO_WP_TIGHT})  : {sel_iso.sum():6d}/{n} ({100*sel_iso.mean():6.2f}%)")
    print(f"  all three          : {sel_mu.sum():6d}/{n} ({100*sel_mu.mean():6.2f}%)")

    print("\n[MUO ID/ISO/trigger SF]  at the corrected pt, signed eta, on selected rows")
    cset_z = correctionlib.CorrectionSet.from_file(JSON_MUO_Z)
    w_mu = {}
    for label, key in (("id", MUO_ID_KEY), ("iso", MUO_ISO_KEY), ("trig", MUO_TRIG_KEY)):
        c = cset_z[key]
        for v in ("nominal", "systup", "systdown"):
            w = np.ones(n)
            # No clamping: these JSONs are flow="error", so a clamp would only convert a
            # loud out-of-range into a silently wrong bin. Rows outside the selection
            # keep weight 1 and are excluded from the mean below instead.
            for i in np.where(sel_mu)[0]:
                w[i] = c.evaluate(float(mu_eta_signed[i]), float(mu_pt_corr[i]), v)
            w_mu[(label, v)] = w
            m, s = wmean(w[sel_mu], genw[sel_mu])
            print(f"  w_mu_{label}[{v}]: mean={m:.4f} std={s:.4f}")
    w_mu_total = w_mu[("id", "nominal")] * w_mu[("iso", "nominal")] * w_mu[("trig", "nominal")]
    m, s = wmean(w_mu_total[sel_mu], genw[sel_mu])
    print(f"  w_mu_total[nominal] = id x iso x trig: mean={m:.4f} std={s:.4f}")
    print("  (the three chain exactly: DEN(iso)=NUM(id)=TightID, DEN(trig)="
          "CutBasedIdTight_and_PFIsoTight=NUM(id)xNUM(iso))")
    print("  CAVEAT: our HLT is IsoMu24||Mu50, but the JSON has no two-path OR key --"
          " the one used adds CascadeMu100/HighPtTkMu100, which only contribute above"
          " ~100 GeV, so the difference is confined to the high-pt tail.")

    # Select-before-correct: the producer required pt>26 on the UNCORRECTED pt. Migration
    # across 26 can only go DOWNWARD: the official MuonScaRe (MuonScaRe.cc pt_scale /
    # pt_resol, low_pt_threshold = 26 by default; MuonScaRe.py filter_boundaries likewise)
    # returns the input pt unchanged for pt < 26, so a muon below 26 before correction is
    # still below 26 after it and could never have passed the corrected-pt cut. The events
    # the producer dropped therefore contain nothing this cut would keep. The downward half
    # -- rows above 26 that the correction pushes below -- is removed by the cut above.
    #
    # An earlier version of this comment claimed an unrecoverable "mirror half" and argued
    # for a 20 GeV producer pre-cut from the maximum upward shift (1.0856 over 65k muons).
    # That shift was measured on muons already above 26 and is irrelevant to crossing 26
    # from below; the argument overlooked the recipe's own low-pt threshold.
    n_drop = int((mu_pt_corr < MUO_LOW_PT_THRESHOLD).sum())
    print(f"  [select-before-correct] rows whose corrected pt falls below"
          f" {MUO_LOW_PT_THRESHOLD}: {n_drop}/{n} ({100*n_drop/n:.3f}%) -- removed by the"
          f" cut above. No upward migration exists: MuonScaRe leaves pt<26 uncorrected.")

    # ---- 3c. MUO: prompt second-muon veto ----
    # Replaces the producer's (num_os_mu_tight==1 & num_ss_mu_tight==0), which counts every
    # tight muon at any pt with no isolation. Measured on the existing nano: the extra muons
    # it vetoes are 77-86% below 10 GeV, 81-89% non-isolated and mostly inside jets --
    # heavy-flavour decay products, not a second prompt lepton (which is what a DY / Z->mumu
    # / ttbar-dilepton veto is for). The soft, in-jet component is a separate, possible
    # ttbar handle to be studied once ttbar or data nano exist; it has no POG SF.
    #
    # Prompt veto = no OTHER muon with tight ID, |eta|<2.4, corrected pt>=15, Loose PF iso:
    #   pt>=15    muon_Z.json ID SF floor (the SF key sets the cut, as for the selected muon)
    #   |eta|<2.4 SF eta edges
    #   Loose iso NUM_LoosePFIso_DEN_TightID exists; for a veto the looser WP is the safer
    #             choice, and it differs from Tight by 0.02% of 2J events.
    # Charge-blind: the selected muon is itself counted in num_os (OS mode) or num_ss (SS
    # mode), so "no extra muon" is the equivalent condition in both modes.
    # No veto SF is applied: the veto removes 0-0.1% of signal / W+jets events, so the SF
    # effect is that fraction times (SF-1). Revisit when DY / ttbar MC enter.
    veto_pass = np.ones(n, dtype=bool)
    for i in range(n):
        mpt = np.asarray(a["Muon_pt"][i])
        cand = np.where((np.asarray(a["Muon_tightId"][i]) == 1)
                        & (np.abs(np.asarray(a["Muon_eta"][i])) < 2.4)
                        & (np.asarray(a["Muon_pfIsoId"][i]) >= 2)
                        & (mpt >= 15.0))[0]
        for j in cand:
            if j == idx_mu[i]:
                continue
            eta_j = float(np.clip(a["Muon_eta"][i][j], -2.399, 2.399))
            phi_j = float(a["Muon_phi"][i][j])
            pt_j = muo_pt_scale(cset, False, float(mpt[j]), eta_j, phi_j, int(a["Muon_charge"][i][j]))
            pt_j = muo_pt_resol(cset, pt_j, eta_j, phi_j, float(a["Muon_nTrackerLayers"][i][j]),
                                event_arr[i], lumi_arr[i])
            if pt_j >= 15.0:
                veto_pass[i] = False
                break
    old_veto = ((np.asarray(a["num_os_mu_tight"]) == 1) & (np.asarray(a["num_ss_mu_tight"]) == 0))
    print(f"\n[MUO prompt second-muon veto]  (tight ID, |eta|<2.4, corrected pt>=15, Loose PF iso)")
    print(f"  pass: {veto_pass.sum()}/{n} ({100*veto_pass.mean():.3f}%)   producer's"
          f" num_os==1&num_ss==0 would pass: {old_veto.sum()}/{n} ({100*old_veto.mean():.3f}%)")
    print(f"  selected muon & veto: {(sel_mu & veto_pass).sum()}/{n}"
          f" ({100*(sel_mu & veto_pass).mean():.2f}%)")

    # ---- 4. JME: jet veto map -- moved into section 5 ----
    # It needs the corrected jet pt (see jme_corrections.jet_veto_flag), which only
    # exists after the JEC/JER chain has run.

    # ---- 5. JME: the actual JEC/JER chain, applied ----
    # Not a lookup any more. nano `Jet_pt` carries the GlobalTag's JEC payload and NO JER
    # smearing at all (there is no smearing module anywhere in jetsAK4_Puppi_cff.py), so
    # the analysis jet pt is: undo rawFactor -> reapply Summer22_22Sep2023_V4_MC ->
    # smear with Summer22_22Sep2023_JRV2_MC. See jme_corrections.py for the full
    # derivation, the ordering's primary source, and which quantities this can and
    # cannot rebuild.
    jme_refs = jme.JmeRefs(year=jme.YEAR_PREEE, is_data=False)
    jet_pt_var = {v: jme.corrected_jet_pt(a, jme_refs, v)
                  for v in ("nominal", "jes_up", "jes_down", "jer_up", "jer_down")}
    jet_pt_corr = jet_pt_var["nominal"]

    nano_pt_flat = ak.to_numpy(ak.flatten(a["Jet_pt"]))
    corr_pt_flat = ak.to_numpy(ak.flatten(jet_pt_corr))
    gen_matched = ak.to_numpy(ak.flatten(a["Jet_genJetIdx"])) >= 0
    print("\n[JME JEC+JER applied]  (tag Summer22_22Sep2023_V4_MC / JRV2_MC)")
    print(f"  jets: {len(nano_pt_flat)}, gen-matched {gen_matched.mean()*100:.1f}%"
          f" (matched jets take the deterministic smear branch, the rest the stochastic one)")
    print(f"  corrected/nano Jet_pt: mean={(corr_pt_flat/nano_pt_flat).mean():.5f}"
          f" median={np.median(corr_pt_flat/nano_pt_flat):.5f}")

    # Jet veto map, official JvmApplication.apply_jvm_flag, evaluated per variation on
    # that variation's corrected pt (as applyJecAndJvm.C does, once per pass). An earlier
    # version looked up EVERY jet with no selection; on Signal_DsW output_0 that vetoed
    # 161 events instead of the official 106 -- 52 of the extra in-map jets fail
    # jetId>=6 and 7 fall below pt 15 after correction.
    print("  jet veto map (official apply_jvm_flag; jetId>=6, pt>=15 corrected,"
          " chEmEF+neEmEF<=0.9, |eta|<=5.191):")
    jvm_veto = {}
    for v in ("nominal", "jes_up", "jes_down", "jer_up", "jer_down"):
        jvm_veto[v] = ak.to_numpy(jme.jet_veto_flag(a, jet_pt_var[v]))
        flips = int((jvm_veto[v] != jvm_veto["nominal"]).sum())
        print(f"    {v:9s} vetoed {jvm_veto[v].sum():7d}/{n} ({100*jvm_veto[v].mean():5.2f}%)"
              f"  flips vs nominal: {flips}")

    tight_jet = (a["Jet_jetId"] & 2) == 2
    print("  event-level quantities recomputed from the corrected pt:")
    print(f"    {'variation':10s} {'<jet_ht_pt20>':>14s} {'<num_jets_pt30>':>16s} {'<jet_leadpt>':>13s}")
    for v in ("nominal", "jes_up", "jes_down", "jer_up", "jer_down"):
        pt = jet_pt_var[v]
        sel = tight_jet & (pt > 20)
        ht = ak.to_numpy(ak.sum(ak.where(sel, pt, 0.0), axis=1))
        nj = ak.to_numpy(ak.sum(tight_jet & (pt > 30), axis=1))
        ld = ak.to_numpy(ak.fill_none(ak.max(ak.where(tight_jet, pt, 0.0), axis=1), 0.0))
        print(f"    {v:10s} {wmean(ht, genw)[0]:14.4f} {wmean(nj, genw)[0]:16.4f}"
              f" {wmean(ld, genw)[0]:13.4f}")
    sel0 = tight_jet & (a["Jet_pt"] > 20)
    ht0 = ak.to_numpy(ak.sum(ak.where(sel0, a["Jet_pt"], 0.0), axis=1))
    nj0 = ak.to_numpy(ak.sum(tight_jet & (a["Jet_pt"] > 30), axis=1))
    print(f"    {'(nano)':10s} {wmean(ht0, genw)[0]:14.4f} {wmean(nj0, genw)[0]:16.4f}"
          f"   <- uncorrected, for reference")
    print("    jet_ht (the DNN feature) is deliberately absent: it sums pt with no")
    print("    threshold, so ~10% of it sits in jets nano never stored -- see the")
    print("    jet-variable closure self-check above. jet_ht_pt20 is its usable form.")

    # Type-1 MET, JES/JER propagated. get_corrected_met rebuilds its own CorrT1+AK4 jet
    # chain internally, so it takes the raw nano columns rather than jet_pt_corr.
    print("  Type-1 MET (official get_corrected_met, CorrT1+AK4, JES/JER propagated):")
    met_nano = np.asarray(a["PuppiMET_pt"]) if "PuppiMET_pt" in a.fields else None
    for v in ("nominal", "jes_up", "jes_down", "jer_up", "jer_down"):
        m = ak.to_numpy(jme.corrected_met_pt(a, jme_refs, v))
        print(f"    met_pt[{v:9s}] mean={wmean(m, genw)[0]:8.4f}")
    if met_nano is not None:
        print(f"    met_pt[{'(nano)':9s}] mean={wmean(met_nano, genw)[0]:8.4f}"
              f"   <- uncorrected PuppiMET_pt, for reference")

    # ---- 6. BTV: b-tag SF on the exact veto jet set (self-check #2's mask) ----
    cset = correctionlib.CorrectionSet.from_file(JSON_BTV)
    wp_json = cset["particleNet_wp_values"].evaluate("T")
    print("\n[BTV WP cross-check]")
    print(f"  particleNet_wp_values('T') = {wp_json}  vs producer hardcoded {BTAG_WP_TIGHT}"
          f"  -> {'MATCH' if abs(wp_json - BTAG_WP_TIGHT) < 1e-6 else 'MISMATCH'}")

    c_comb, c_light = cset["particleNet_comb"], cset["particleNet_light"]
    # This sample's own efficiency map (build_btag_efficiency_map.py), consumed as the
    # official compute_btag_event_weight_multi_wp() does: .evaluate(wp, flav, abseta, pt).
    eff_map = correctionlib.CorrectionSet.from_file(
        BTAG_EFF_MAP_JSON[sample])["particleNet_efficiency"]

    # The jets the veto examines: tight ID, |eta|<2.5, dR(Ds,jet)>0.4 -- plus, for the
    # analysis veto, CORRECTED pt >= BTAG_JET_PT_MIN (the BTV SF validity floor).
    dphi_j = abs(a["Jet_phi"] - a["Ds_phi"])
    dphi_j = ak.where(dphi_j > np.pi, 2 * np.pi - dphi_j, dphi_j)
    dR_j = np.sqrt((a["Jet_eta"] - a["Ds_eta"]) ** 2 + dphi_j ** 2)
    examined = (abs(a["Jet_eta"]) < 2.5) & (a["Jet_jetId"] >= 2) & (dR_j > 0.4)
    tagged = examined & (a["Jet_btagPNetB"] > BTAG_WP_TIGHT)

    # Closure: the producer's own pt>20 count, on nano (uncorrected) pt. Confirms the
    # definition below differs from the stored pt20 branch only by using corrected pt.
    nb_pt20_nano = ak.to_numpy(ak.sum(tagged & (a["Jet_pt"] > 20.0), axis=1))
    nb_pt20_stored = np.asarray(a["n_b_jets_tight_pt20_dR04"])
    print(f"\n[BTV b-veto definition]")
    print(f"  closure n_b_jets_tight_pt20_dR04 (nano pt>20) vs stored:"
          f" {(nb_pt20_nano != nb_pt20_stored).sum()}/{n} mismatches")

    # Analysis veto: no tagged examined jet with corrected pt >= 20.
    nb_corr = ak.to_numpy(ak.sum(tagged & (jet_pt_corr >= BTAG_JET_PT_MIN), axis=1))
    veto_pass = nb_corr == 0
    old_pass = np.asarray(a["n_b_jets_tight_dR04"]) == 0
    print(f"  veto on corrected pt>={BTAG_JET_PT_MIN:g}: pass {veto_pass.sum()}/{n}"
          f" ({100*veto_pass.mean():.2f}%);  stored n_b_jets_tight_dR04==0 (no pt cut):"
          f" pass {old_pass.sum()}/{n};  events whose veto decision changes:"
          f" {(veto_pass != old_pass).sum()}")

    # w_btag = Pi_j (1-SF_j*eps_j)/(1-eps_j) over the examined jets of veto-passing
    # events (all untagged by construction, so only the untagged branch of the fixed-WP
    # formula is needed). btagging.json has zero compound_corrections; this is the BTV
    # recipe as implemented in the official lib/btag_sf.py. SF and eps are evaluated at
    # the corrected pt. No try/except: |eta|<2.5 and pt>=20 are guaranteed by the mask,
    # and above pt 1000 the JSON returns a flow constant rather than raising, so those
    # jets are counted and reported instead of being silently swallowed.
    scored = examined & (jet_pt_corr >= BTAG_JET_PT_MIN)
    scored = scored & ak.broadcast_arrays(ak.Array(veto_pass), scored)[0]
    counts = ak.to_numpy(ak.num(scored[scored], axis=1))
    pt_s = ak.to_numpy(ak.flatten(jet_pt_corr[scored])).astype(float)
    abseta_s = np.abs(ak.to_numpy(ak.flatten(a["Jet_eta"][scored]))).astype(float)
    flav_s = ak.to_numpy(ak.flatten(a["Jet_hadronFlavour"][scored])).astype(int)
    n_above = int((pt_s >= 1000.0).sum())
    eps = eff_map.evaluate("T", flav_s, abseta_s, pt_s)
    heavy = (flav_s == 4) | (flav_s == 5)
    w_btag = {}
    for v, key in (("nom", "central"), ("up", "up"), ("down", "down")):
        sf = np.empty_like(pt_s)
        sf[heavy] = c_comb.evaluate(key, "T", flav_s[heavy], abseta_s[heavy], pt_s[heavy])
        sf[~heavy] = c_light.evaluate(key, "T", flav_s[~heavy], abseta_s[~heavy], pt_s[~heavy])
        per_jet = np.where(eps < 1.0, (1.0 - sf * eps) / (1.0 - eps), 1.0)
        w = ak.to_numpy(ak.prod(ak.unflatten(per_jet, counts), axis=1))
        w_btag[v] = np.where(veto_pass, w, np.nan)
    print(f"\n[BTV b-tag SF]  (this sample's efficiency map: {BTAG_EFF_MAP_JSON[sample]})")
    print("  SF and efficiency evaluated at the JME-corrected jet pt; jets pt>=20 only.")
    print(f"  veto-passing events: {veto_pass.sum()}/{n}, jets scored: {len(pt_s)},"
          f" of which pt>=1000 (JSON flow constant, no measured SF): {n_above}")
    for v in ("nom", "up", "down"):
        m, s = wmean(w_btag[v][veto_pass], genw[veto_pass])
        print(f"  w_btag[{v}]: mean={m:.4f} std={s:.4f}  (over veto-passing events only)")

    # ---- 7. EGM: NOT APPLICABLE to this analysis ----
    # Earlier this was recorded as "blocked pending the electron-veto cut definition".
    # That framing was wrong: there is no pending decision, because the analysis applies
    # no electron selection at all. The selection is, in full (analysis/DNN/common.py
    # load(), with the branch names coming from analysis/eras.py):
    #     num_os_mu_tight == 1 & num_ss_mu_tight == 0 & n_b_jets_tight_dR04 == 0
    # No electron cut, no electron veto, and FEATURE_NAMES contains no electron variable
    # either. An ID scale factor corrects the efficiency of a cut that was applied; with
    # no cut applied there is nothing for it to correct, exactly as for TAU.
    #
    # Checked for an official combination rule anyway, in case one existed to adopt:
    # electron.json.gz has zero compound_corrections, and cms-analysis-corrections/EGM's
    # examples (electronExample.py, electronHltExample.py, photonExample.py) only ever do
    # single-object lookups -- the "veto" in that repo's README is the photon CSEV/PixVeto
    # discriminator, an unrelated object-level quantity. So no event-level electron-veto
    # weight is defined by EGM for us to apply even if we wanted one.
    #
    # The raw per-electron lookups are still printed below as a live check that the JSON
    # loads and the supercluster-eta convention (Electron_eta + Electron_deltaEtaSC) is
    # right -- so that if an electron veto is ever ADDED to the selection, this is ready.
    # No event weight is assembled, and none should be until such a cut exists.
    cset = correctionlib.CorrectionSet.from_file(JSON_EGM)
    c = cset["Electron-ID-SF"]
    sf_vals = []
    n_with_ele = 0
    for i in range(n):
        pts = np.asarray(a["Electron_pt"][i])
        etas = np.asarray(a["Electron_eta"][i]) + np.asarray(a["Electron_deltaEtaSC"][i])
        if len(pts) == 0:
            continue
        n_with_ele += 1
        for p, e in zip(pts, etas):
            sf_vals.append(c.evaluate(EGM_YEAR, "sf", "Veto", float(e), float(max(p, 10.0))))
    sf_vals = np.array(sf_vals)
    print("\n[EGM] NOT APPLICABLE -- the analysis selection applies no electron cut at all")
    print("  (num_os_mu_tight==1 & num_ss_mu_tight==0 & n_b_jets_tight_dR04==0; no electron")
    print("  term, and no electron feature in FEATURE_NAMES). Lookups below are a liveness")
    print("  check on the JSON and the supercluster-eta convention only -- no weight is built.")
    print(f"  events with >=1 electron: {n_with_ele}/{n}, electrons scored: {len(sf_vals)}")
    if len(sf_vals):
        print(f"  sf: mean={sf_vals.mean():.4f} range=[{sf_vals.min():.4f},{sf_vals.max():.4f}]")

    print("\n[TAU] NOT APPLICABLE -- no slimmedTaus/patTau/DeepTau reference anywhere in the"
          " producer; the analysis never reads a reco hadronic-tau object.")


if __name__ == "__main__":
    main()
