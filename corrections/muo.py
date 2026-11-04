"""MUO: muon momentum scale/smearing and muon scale factors.

Scale/smearing: the official muonscarekit MuonScaRe.cc (external/muonscarekit, unmodified),
compiled through ROOT exactly as in muonscarekit scripts/apply_corrections_rdf.py: a global
`cset` holding muon_scalesmearing.json, then `#include "MuonScaRe.cc"`. The functions below
only call it; the one adaptation is that the event number is passed as a 32-bit int, as an
RDataFrame ULong64_t `event` column passed to the .cc's `int evtNumber` would be.
On top of the .cc (kept from the previous port): a non-finite scale/resolution VARIATION
falls back to the nominal corrected pt (.cc pt_scale_var gives NaN when its sqrt argument
is negative).

Scale factors: muon_Z.json, evaluated as in MUO/examples muonExample.py
(signed eta, pt, "nominal"/"systup"/"systdown").
"""
import os
import sys

import awkward as ak
import correctionlib
import numpy as np

POG_DIR = "/cvmfs/cms-griddata.cern.ch/cat/metadata"
MUO_VERSION = "2026-06-18"
SCARE_CC = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                                         "external", "muonscarekit", "MuonScaRe.cc"))

MUO_ID_KEY = "NUM_TightID_DEN_TrackerMuons"
MUO_ISO_KEY = "NUM_TightPFIso_DEN_TightID"
MUO_LOOSE_ISO_KEY = "NUM_LoosePFIso_DEN_TightID"
MUO_TRIG_KEY = "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight"

MU_PRECUT = 24.0       # HDsPhiPiMuProducer: best muon = highest-pt tight muon, pt >= 24
MU_PT_MIN = 26.0       # trigger-SF floor, on corrected pt
MU_ETA_MAX = 2.4
VETO_MU_PT_MIN = 15.0  # ID-SF floor

MU_VARS = ("nom", "muScaleUp", "muScaleDown", "muResUp", "muResDown")


def json_dir(campaign):
    return f"{POG_DIR}/MUO/{campaign}/{MUO_VERSION}"


def sf_set(campaign):
    return correctionlib.CorrectionSet.from_file(f"{json_dir(campaign)}/muon_Z.json.gz")


# ---------------------------------------------------------------- MuonScaRe.cc
_ROOT = None
_LOADED = None


def load_scalesmearing(campaign):
    """Compile MuonScaRe.cc against this campaign's muon_scalesmearing.json (once per process:
    the .cc reads one global `cset`)."""
    global _ROOT, _LOADED
    path = f"{json_dir(campaign)}/muon_scalesmearing.json.gz"
    if _LOADED == path:
        return
    if _LOADED is not None:
        sys.exit(f"MuonScaRe.cc already loaded with {_LOADED}; cannot switch to {path}")
    import ROOT
    correctionlib.register_pyroot_binding()
    ROOT.gROOT.ProcessLine(f'auto cset = correction::CorrectionSet::from_file("{path}");')
    ROOT.gROOT.ProcessLine(f'#include "{SCARE_CC}"')
    for f in ("pt_scale", "pt_resol", "pt_scale_var", "pt_resol_var"):
        if not hasattr(ROOT, f):
            sys.exit(f"MuonScaRe.cc did not compile ({f} missing)")
    _ROOT, _LOADED = ROOT, path


def _i32(x):
    """uint64 -> C++ int (two's complement wrap), as RDataFrame passes `event` to the .cc."""
    return int(((int(x) + 2 ** 31) % 2 ** 32) - 2 ** 31)


def pt_scale(is_data, pt, eta, phi, charge):
    return _ROOT.pt_scale(bool(is_data), float(pt), float(eta), float(phi), int(charge))


def pt_resol(pt, eta, phi, nL, event, lumi):
    return _ROOT.pt_resol(float(pt), float(eta), float(phi), float(nL), _i32(event), _i32(lumi))


def pt_corrected(is_data, pt, eta, phi, charge, nL, event, lumi):
    """Nominal corrected pt: scale, plus resolution smearing on MC."""
    p = pt_scale(is_data, pt, eta, phi, charge)
    return p if is_data else pt_resol(p, eta, phi, nL, event, lumi)


def _var(x, fallback):
    return x if np.isfinite(x) else fallback


# ---------------------------------------------------------------- selected muon
def locate_selected_muon(a):
    """Index in Muon_* of the producer's best muon; asserts it reproduces mu_pt/eta/charge."""
    sel = (a["Muon_pt"] >= MU_PRECUT) & a["Muon_tightId"]
    idx = ak.argmax(ak.where(sel, a["Muon_pt"], -1.0), axis=1, keepdims=True)
    pt, eta, ch = (ak.to_numpy(ak.flatten(a[b][idx])) for b in ("Muon_pt", "Muon_eta", "Muon_charge"))
    ok = ((np.abs(pt - ak.to_numpy(a["mu_pt"])) < 1e-3) & (np.abs(eta - ak.to_numpy(a["mu_eta"])) < 1e-3)
          & (ch == ak.to_numpy(a["mu_charge"]).astype(int)))
    if not ok.all():
        sys.exit(f"selected-muon re-id FAILED on {int((~ok).sum())}/{len(ok)} rows")
    return ak.to_numpy(ak.flatten(idx))


def muon_momenta(a, idx, is_data=False):
    """Corrected pt of the selected muon for every MU_VARS entry (data: "nom" only, scale only)."""
    n = len(a)
    out = {v: np.empty(n) for v in (("nom",) if is_data else MU_VARS)}
    pt_raw = ak.to_numpy(a["mu_pt"]).astype(float)
    eta = ak.to_numpy(a["mu_eta"]).astype(float)
    phi = ak.to_numpy(a["mu_phi"]).astype(float)
    ch = ak.to_numpy(a["mu_charge"]).astype(int)
    nL = ak.to_numpy(ak.flatten(a["Muon_nTrackerLayers"][idx[:, None]])).astype(float)
    evt, lumi = ak.to_numpy(a["event"]), ak.to_numpy(a["luminosityBlock"])
    R = _ROOT
    for i in range(n):
        p_scale = pt_scale(is_data, pt_raw[i], eta[i], phi[i], ch[i])
        if is_data:
            out["nom"][i] = p_scale          # no resolution smearing on data
            continue
        p = pt_resol(p_scale, eta[i], phi[i], nL[i], evt[i], lumi[i])
        out["nom"][i] = p
        out["muScaleUp"][i] = _var(R.pt_scale_var(p, eta[i], phi[i], int(ch[i]), "up"), p)
        out["muScaleDown"][i] = _var(R.pt_scale_var(p, eta[i], phi[i], int(ch[i]), "dn"), p)
        out["muResUp"][i] = _var(R.pt_resol_var(p_scale, p, eta[i], "up"), p)
        out["muResDown"][i] = _var(R.pt_resol_var(p_scale, p, eta[i], "dn"), p)
    return out


def prompt_muon_veto(a, idx, is_data=False):
    """True = passes (no other prompt-like muon: tight ID, |eta|<2.4, loose PF iso, nominal
    corrected pt >= 15)."""
    n = len(a)
    passed = np.ones(n, dtype=bool)
    evt, lumi = ak.to_numpy(a["event"]), ak.to_numpy(a["luminosityBlock"])
    cand = ((a["Muon_tightId"] == 1) & (abs(a["Muon_eta"]) < MU_ETA_MAX)
            & (a["Muon_pfIsoId"] >= 2) & (a["Muon_pt"] >= VETO_MU_PT_MIN)
            & (ak.local_index(a["Muon_pt"]) != idx))
    for i in np.where(ak.to_numpy(ak.any(cand, axis=1)))[0]:
        for j in np.where(ak.to_numpy(cand[i]))[0]:
            pt_j = pt_corrected(is_data, a["Muon_pt"][i][j], a["Muon_eta"][i][j], a["Muon_phi"][i][j],
                                a["Muon_charge"][i][j], a["Muon_nTrackerLayers"][i][j], evt[i], lumi[i])
            if pt_j >= VETO_MU_PT_MIN:
                passed[i] = False
                break
    return passed


# ---------------------------------------------------------------- scale factors
def sf(cset_z, key, eta, pt, var, sel):
    """SF on the rows in `sel`, NaN elsewhere (the keys have flow="error": no clamping)."""
    w = np.full(len(eta), np.nan)
    w[sel] = cset_z[key].evaluate(eta[sel], pt[sel], var)
    return w
