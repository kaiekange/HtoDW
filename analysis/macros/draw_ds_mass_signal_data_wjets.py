"""Ds_mass (M(K+K-pi+-)) and phi_mass (M(K+K-)) overlay, side by side: truth-matched
signal, OS data, W+jets MC.

No MVA cut -- base analysis selection only (OS, b-veto), same as the "no cut" column
of compare_models.py's mass panel. Signal is truth-matched (Ds_match & mu_match), same
definition as draw_wp_mass.py / draw_signal_dsst_mass_truth.py. Signal and W+jets use
their absolute per-row `weight` (sigma*L/sum(w), see project_nano_merged_samples
memory); OS data is unweighted counts.

Blinding is a property of the EVENT, not of which mass variable is on the x-axis (same
convention as draw_wp_mass.py): the Ds_mass SR mask is computed once from OS data's
Ds_mass and applied to both the Ds_mass and phi_mass panels for OS data. W+jets and
signal are drawn over the full range in both panels, unblinded
(feedback_blind_sideband_care memory: masked with NaN, never forced to 0).

Signal's absolute yield is ~1e4-1e5x below data/W+jets in this mass window (a few
hundred pb-scale signal events across all of Run 2022 vs a bulk QCD-adjacent
background), so it is invisible on a shared linear axis at its true scale. The script
prints two different scale factors, computed once from the Ds_mass panel and reused
for phi_mass too (one uniform multiplier for both panels, same convention as
compare_dsmunu_signal_gentag.py):
  (a) integral ratio  -- signal_yield / wjets_yield (full range) and
                          signal_yield / data_sideband_yield (SR-excluded range)
                          for the record; NOT what is used to draw the plot.
  (b) peak-bin ratio  -- W+jets height in the M_Ds bin (unblinded) / signal height in
                          the same bin. This is the multiplier applied to the drawn
                          signal curve, because it is what makes the shapes visually
                          comparable on the page; it is a display choice, not a yield
                          statement, and the legend says so explicitly.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_ds_mass_signal_data_wjets.py
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, M_PHI, SR_HALF = 1.96835, 1.019460, 0.040
LUMI = common.ERAS[1]["lumi_pb"] + common.ERAS[2]["lumi_pb"]
SEL = [*common.SELECTION]  # num_os_mu_tight, num_ss_mu_tight, n_bjet_corr

VARS = [   # (branch, xlabel, bin edges, blind half-width or None, ref line or None)
    ("phi_mass", r"$M(K^+K^-)$ [GeV]",
     np.arange(0.984, 1.074 + 1e-9, 0.002), None, M_PHI),
    ("Ds_mass", r"$M(K^+K^-\pi^\pm)$ [GeV]",
     np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010), SR_HALF, None),
]


def ossel(a):
    return ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
             & (a[common.B_VETO] == 0))


def read(name, extra):
    path = f"{common.SAMPLES_DIR}/{name}.root"
    return uproot.open(path)["Events"].arrays(
        ["Ds_mass", "phi_mass"] + SEL + extra, library="np")


# ---- truth-matched signal (OS + b-veto + Ds_match & mu_match), weighted ----
s = read("Signal_merged", ["weight", "Ds_match", "mu_match"])
sm = ossel(s) & s["Ds_match"].astype(bool) & s["mu_match"].astype(bool)
sig_ds, sig_phi, sig_w = s["Ds_mass"][sm], s["phi_mass"][sm], s["weight"][sm]

# ---- W+jets MC (OS + b-veto), weighted, full mass range ----
w = read("WtoLNu2Jets_merged", ["weight"])
wm = ossel(w)
wj_ds, wj_phi, wj_w = w["Ds_mass"][wm], w["phi_mass"][wm], w["weight"][wm]

# ---- OS data (OS + b-veto), unweighted; SR blind is an EVENT-level mask, reused for
# both panels (see module docstring) ----
d = read("Data_OS_merged", [])
dm = ossel(d)
d_ds, d_phi = d["Ds_mass"][dm], d["phi_mass"][dm]
d_blind_event = np.abs(d_ds - M_DS) < SR_HALF   # event-level, not bin-level

print(f"selected rows: signal(truth-matched) {sm.sum():,}  W+jets {wm.sum():,}  "
      f"OS data {dm.sum():,}  (OS data in Ds SR, excluded from both panels: "
      f"{d_blind_event.sum():,})")

# ---- scale factor, computed once from Ds_mass (regardless of panel order below)
# and reused for phi_mass -- one uniform multiplier for both panels ----
ds_bins = np.arange(M_DS - 0.3, M_DS + 0.3 + 1e-9, 0.010)
ds_bc = 0.5 * (ds_bins[1:] + ds_bins[:-1])
h_sig_ds, _ = np.histogram(sig_ds, bins=ds_bins, weights=sig_w)
h_wj_ds, _ = np.histogram(wj_ds, bins=ds_bins, weights=wj_w)
h_d_ds, _ = np.histogram(d_ds, bins=ds_bins)
peak_bin = np.argmin(np.abs(ds_bc - M_DS))
factor_peak = h_wj_ds[peak_bin] / h_sig_ds[peak_bin]
factor_int_wj = h_wj_ds.sum() / h_sig_ds.sum()
factor_int_data = h_d_ds.sum() / h_sig_ds.sum()
factor_round = float(f"{factor_peak:.1g}")   # 1 sig fig, for a clean legend
print(f"signal peak-bin height (M_Ds bin, unscaled): {h_sig_ds[peak_bin]:.4g}  "
      f"W+jets peak-bin height: {h_wj_ds[peak_bin]:.4g}")
print(f"scale factor (a) integral, signal->W+jets (full range): {factor_int_wj:,.0f}")
print(f"scale factor (a) integral, signal->OS data (Ds SR excluded): "
      f"{factor_int_data:,.0f}")
print(f"scale factor (b) peak-bin, signal->W+jets (used to draw both panels): "
      f"{factor_peak:,.0f}  (rounded {factor_round:,.0f})")

fig = plt.figure(figsize=(17.0, 8.5))
outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.24)

for col, (branch, xlabel, bins, blind_half, refline) in enumerate(VARS):
    inner = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[col], height_ratios=[3, 1], hspace=0.05)
    ax, rax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
    sig_vals = sig_ds if branch == "Ds_mass" else sig_phi
    wj_vals = wj_ds if branch == "Ds_mass" else wj_phi
    d_vals = d_ds if branch == "Ds_mass" else d_phi

    bc = 0.5 * (bins[1:] + bins[:-1])
    width = bins[1] - bins[0]
    vis = ~d_blind_event

    h_sig, _ = np.histogram(sig_vals, bins=bins, weights=sig_w)
    e_sig = np.sqrt(np.histogram(sig_vals, bins=bins, weights=sig_w ** 2)[0])

    h_wj, _ = np.histogram(wj_vals, bins=bins, weights=wj_w)
    e_wj = np.sqrt(np.histogram(wj_vals, bins=bins, weights=wj_w ** 2)[0])

    h_d, _ = np.histogram(d_vals[vis], bins=bins)
    assert h_d.sum() == vis.sum(), (
        f"{branch}: histogram dropped {vis.sum() - h_d.sum()} OS data events as "
        f"under/overflow -- bin range [{bins[0]:.4f},{bins[-1]:.4f}] does not cover "
        f"the data.")
    e_d = np.sqrt(h_d)

    ax.bar(bc, h_wj, width=width, color="seagreen", alpha=0.55, edgecolor="none",
           label=f"W+jets MC  N={wm.sum():,}  yield={h_wj.sum():,.1f}")
    ax.bar(bc, 2 * e_wj, width=width, bottom=h_wj - e_wj, color="none",
           hatch="////", edgecolor="0.35", linewidth=0.0, zorder=3,
           label="W+jets MC stat")

    # h_d was built only from the visible (non-blind-event) rows -- rebuild onto the
    # full bin grid so blind bins on the Ds_mass panel stay NaN, not silently 0.
    h_d_full, _ = np.histogram(d_vals, bins=bins)
    e_d_full = np.sqrt(h_d_full)
    blind_bins = (np.abs(bc - M_DS) < blind_half) if blind_half else np.zeros_like(bc, dtype=bool)
    d_plot = np.where(blind_bins, np.nan, h_d_full)
    ax.errorbar(bc, d_plot, yerr=np.where(blind_bins, 0, e_d_full), fmt="o", ms=4,
                color="black", elinewidth=1.2,
                label=f"OS data{' (SR blind)' if blind_half else ''}  N={h_d.sum():,}")
    if blind_half:
        ax.axvspan(M_DS - blind_half, M_DS + blind_half, color="grey", alpha=0.15,
                   zorder=0)
    if refline:
        ax.axvline(refline, color="green", ls=":", lw=1.4)

    ax.errorbar(bc, h_sig * factor_round, yerr=e_sig * factor_round, fmt="^", ms=5,
                color="firebrick", elinewidth=1.2,
                label=(f"truth-matched signal $\\times${factor_round:,.0f}  "
                       f"(Ds_match & mu_match)  N={sm.sum():,}"))

    top = np.nanmax([h_wj.max(), np.nanmax(d_plot), (h_sig * factor_round).max()])
    ax.set_ylim(0, top * 1.55)
    ax.set_xlim(bins[0], bins[-1])
    ax.tick_params(labelbottom=False)
    if col == 0:
        ax.set_ylabel(f"Events / {width*1000:.0f} MeV")
    ax.legend(fontsize=8.5, loc="upper right",
              title="no MVA cut -- OS + b-veto only", title_fontsize=8.5)
    if col == 0:
        hep.cms.label(ax=ax, data=True, text="Preliminary",
                      rlabel=rf"{LUMI/1000:.1f} fb$^{{-1}}$ (13.6 TeV)", fontsize=10)

    # ---- data / W+jets ratio subpanel ----
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(blind_bins | (h_wj <= 0), np.nan, d_plot / h_wj)
        re = np.where(blind_bins | (h_wj <= 0), np.nan, np.sqrt(h_d_full) / h_wj)
        mc_band = np.where(h_wj > 0, e_wj / h_wj, np.nan)
    rax.bar(bc, 2 * mc_band, width=width, bottom=1 - mc_band, color="none",
            hatch="////", edgecolor="0.35", linewidth=0.0)
    rax.errorbar(bc, r, yerr=re, fmt="o", ms=3.5, color="black", elinewidth=1.0)
    rax.axhline(1.0, color="red", lw=1.0)
    if blind_half:
        rax.axvspan(M_DS - blind_half, M_DS + blind_half, color="grey", alpha=0.15)
    if refline:
        rax.axvline(refline, color="green", ls=":", lw=1.4)
    # y-range covers every CR (non-blinded) ratio point plus its error bar, not a
    # fixed window -- the two panels have very different CR spreads.
    lo, hi = np.nanmin(r - re), np.nanmax(r + re)
    pad = 0.05 * (hi - lo)
    rax.set_ylim(lo - pad, hi + pad)
    rax.set_xlim(bins[0], bins[-1])
    rax.set_xlabel(xlabel)
    if col == 0:
        rax.set_ylabel("data / MC")

out = common.fig_path("preselection", "ds_mass_signal_data_wjets.png", key=None)
plt.savefig(out, dpi=120, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")
