"""Muon momentum scale and resolution (MuonScaRe) plus the MUO scale-factor keys.

Moved verbatim out of the corrections pilot (CMSSW_15_0_5/src/analysis/macros/
check_corrections_signal.py, lines 105-107, 129, 143-305) so that apply_corrections.py no
longer imports the pilot. No code change; only the surrounding module is new.
"""
import math

import ROOT

MUO_ID_KEY = "NUM_TightID_DEN_TrackerMuons"
MUO_ISO_KEY = "NUM_TightPFIso_DEN_TightID"
MUO_TRIG_KEY = "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight"

MUO_LOW_PT_THRESHOLD = 26.0  # matches cms-cat/nanoAOD-tools-modules MuonScaRe default


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
