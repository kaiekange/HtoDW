"""phi helicity angle: why the signal M(pi+pi-) of (pi, K-) has two peaks.

    python3 draw_phi_helicity.py

Ds (spin 0) -> phi (spin 1) pi (spin 0) leaves the phi with helicity 0, so in the phi rest frame the
K- direction relative to the phi flight direction follows cos^2(theta_h). M^2(pi K-) is linear in
cos(theta_h), so the two populated ends cos ~ +-1 become the two M(pipi) peaks.
theta_h: angle between K- and the phi flight direction (= -pi direction) in the phi rest frame, kaon mass
for K+/K-. Panels: (a) truth-matched signal cos(theta_h), no cut and ssdata eps90, with (3/2)cos^2;
(b) signal M(pi+pi-) of (pi, K-) vs cos(theta_h); (c) DY muon-from-tau / prompt muon, ssdata eps90.
DY muon origin from the draw_dy_mu_origin.py cache. Output: os_ss/phi_helicity_signal_dy.png (shared dir).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "DNN"))
import common  # noqa: E402

KEY = "ssdata"
M_PI, M_K = 0.13957039, 0.493677
CACHE = f"{common.MVA_ROOT}/aux/DY_merged_mu_genPartFlav.npz"
TRK = [f"{p}_{v}" for p in ("Kp", "Km", "pi") for v in ("pt", "eta", "phi")]
WANT = list(dict.fromkeys(["Ds_mass", "weight", "score_dsmu", "luminosityBlock", "event", "mu_charge",
                           "pi_charge", *TRK, *common.SELECTION]))


def p4(a, p, m):
    pt, eta, phi = a[f"{p}_pt"], a[f"{p}_eta"], a[f"{p}_phi"]
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return np.stack([np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m * m), px, py, pz])


def boost(v, b):
    """Boost 4-vectors v (4,N) by velocity b (3,N)."""
    b2 = (b ** 2).sum(0)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = (b * v[1:]).sum(0)
    g2 = np.where(b2 > 0, (g - 1.0) / b2, 0.0)
    E = g * (v[0] - bp)
    p = v[1:] + (g2 * bp - g * v[0]) * b
    return np.vstack([E, p])


def helicity(a):
    kp, km, pi = p4(a, "Kp", M_K), p4(a, "Km", M_K), p4(a, "pi", M_PI)
    phi = kp + km
    b = phi[1:] / phi[0]
    km_r, pi_r = boost(km, b), boost(pi, b)
    u1 = km_r[1:] / np.linalg.norm(km_r[1:], axis=0)
    u2 = -pi_r[1:] / np.linalg.norm(pi_r[1:], axis=0)   # phi flight direction in the phi frame
    return (u1 * u2).sum(0)


def m_pipi(a):
    v = p4(a, "pi", M_PI) + p4(a, "Km", M_PI)
    return np.sqrt(np.maximum(v[0] ** 2 - (v[1:] ** 2).sum(0), 0.0))


def load(fname, extra=()):
    with uproot.open(f"{common.out_dir_for(KEY, 'DNN')}/Prediction/{fname}.root")[common.TREE] as t:
        a = t.arrays(list(dict.fromkeys(WANT + list(extra))), library="np")
    common.check_charge_config(a, False, fname)
    return a


thr = common.stage1_threshold({"key": KEY, "tool": "DNN", "eff_sig": 0.90})[0]
sg = load("Signal_merged", ["Ds_match", "mu_match"])
sg_sel = common.selection_mask(sg) & sg["Ds_match"].astype(bool) & sg["mu_match"].astype(bool)
sg_cos, sg_m, sg_w = helicity(sg)[sg_sel], m_pipi(sg)[sg_sel], sg["weight"][sg_sel].astype(float)
sg_cut = sg["score_dsmu"][sg_sel] > thr
dy = load("DY_merged")
c = np.load(CACHE)
assert np.array_equal(c["rowkey"], (dy["luminosityBlock"].astype(np.int64) << 32) | dy["event"].astype(np.int64))
dy_sel = common.selection_mask(dy) & (dy["score_dsmu"] > thr)
dy_cos, dy_w, dy_fl = helicity(dy)[dy_sel], dy["weight"][dy_sel].astype(float), c["flav"][dy_sel]

bins = np.linspace(-1, 1, 41)
bc = 0.5 * (bins[1:] + bins[:-1])
plt.style.use(hep.style.CMS)
fig, axes = plt.subplots(1, 3, figsize=(27, 8.5))


def shape(x, w, ax, **kw):
    h, _ = np.histogram(x, bins=bins, weights=w)
    e = np.sqrt(np.histogram(x, bins=bins, weights=w ** 2)[0])
    norm = w.sum() * np.diff(bins)
    hep.histplot(h / norm, bins, ax=ax, yerr=e / norm, **kw)


ax = axes[0]
shape(sg_cos, sg_w, ax, histtype="step", color="#c0392b", lw=2.2, label="signal, truth-matched, no cut")
shape(sg_cos[sg_cut], sg_w[sg_cut], ax, histtype="errorbar", color="black", marker="o", markersize=4,
      label=rf"signal, {KEY} $\varepsilon_{{sig}}$=90%")
ax.plot(bc, 1.5 * bc ** 2, color="green", ls="--", lw=2, label=r"ideal $\frac{3}{2}\cos^2\theta_h$")
ax.set_xlabel(r"$\cos\theta_h$ (K$^-$ vs $\phi$ flight, $\phi$ frame)")
ax.set_ylabel("normalized to unit area")
ax.set_ylim(0, 2.6)
ax.legend(fontsize=13, loc="upper center")

ax = axes[1]
hh, xe, ye = np.histogram2d(sg_cos, sg_m, bins=[np.linspace(-1, 1, 41), np.linspace(0.6, 1.5, 46)], weights=sg_w)
im = ax.pcolormesh(xe, ye, hh.T, cmap="viridis")
fig.colorbar(im, ax=ax, label="signal weight")
ax.set_xlabel(r"$\cos\theta_h$")
ax.set_ylabel(r"$M(\pi^+\pi^-)$ of $\pi K^-$ [GeV]")
ax.set_title("signal, truth-matched, no cut", fontsize=15)

ax = axes[2]
shape(sg_cos[sg_cut], sg_w[sg_cut], ax, histtype="step", color="#c0392b", lw=2.2, ls="--",
      label=rf"signal, {KEY} $\varepsilon_{{sig}}$=90%")
for lab, fv, col in (("DY, muon from tau", 15, "steelblue"), ("DY, prompt muon", 1, "grey")):
    m = dy_fl == fv
    shape(dy_cos[m], dy_w[m], ax, histtype="errorbar", color=col, marker="o", markersize=4,
          label=rf"{lab}, {KEY} $\varepsilon_{{sig}}$=90%")
ax.set_xlabel(r"$\cos\theta_h$")
ax.set_ylabel("normalized to unit area")
ax.set_ylim(0, 2.6)
ax.legend(fontsize=13, loc="upper center")

hep.cms.label(ax=axes[0], rlabel="2022 (13.6 TeV)", data=False, text="Preliminary", fontsize=12)
plt.tight_layout()
out = common.fig_path("os_ss", "phi_helicity_signal_dy.png")
plt.savefig(out, dpi=90)
plt.close(fig)
w2 = lambda x, w, m: w[m].sum() / w.sum()  # noqa: E731
print(f"signal truth-matched: {sg_sel.sum():,} rows; |cos|>0.7 share {w2(sg_cos, sg_w, np.abs(sg_cos) > 0.7):.3f} "
      f"(ideal cos^2: {1 - 0.7 ** 3:.3f}, flat: 0.300)")
for lab, fv in (("DY mu<-tau", 15), ("DY prompt", 1)):
    m = dy_fl == fv
    print(f"{lab} eps90: |cos|>0.7 share {w2(dy_cos[m], dy_w[m], np.abs(dy_cos[m]) > 0.7):.3f}")
print(f"Saved {out}")
