"""How the postEE data/MC gap was traced to pileup, in one figure.

The question this answers: postEE data loses far more events to the DNN cut than preEE
data does (era double ratio 0.75 / 0.71 / 0.69 at eps_sig = 90/70/50%), while the
W+jets MC shows almost none of that. Which of the 43 inputs is responsible?

THE QUANTITY THAT MATTERS IS THE RESIDUAL, not either shift on its own. For every
feature, measure the preEE -> postEE shift separately in data and in MC, in units of
the preEE DATA sigma so both are on the same ruler:

    d(data) - d(MC)

A feature both agree on is a real era difference the MC reproduces -- fine, no
mismodelling. A feature that moves in data and not in MC is a condition the MC does
not know about, and it is exactly what a classifier trained on MC-vs-data will exploit.

The features are coloured by a grouping declared FROM PHYSICS before looking at the
answer, so that the separation between groups is a result rather than a post-hoc
label. The kinematic/angular group is the negative control: those observables have no
reason to care about pileup, and if they had moved, the measurement itself would be
suspect (e.g. the two era samples not being comparable in the first place).

No DNN cut is applied anywhere here -- this is the population the cut acts on, so the
diagnosis must not itself depend on the model.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 draw_era_feature_shift.py
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

plt.style.use(hep.style.CMS)

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis")
sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common
from eras import PU_PROFILE, PU_PROFILE_HIST

# Read Samples/, not Prediction/: this measurement uses no DNN score, so it must not
# depend on which model happens to be current. The 43 features are identical in both.
SAMPLES = "/eos/home-k/kakang/H2DW/Analysis/MVA/Samples"
M_DS, SR_HALF = 1.96835, 0.04
F = common.FEATURE_NAMES

# Declared from physics BEFORE the measurement, so that "the groups separate" is a
# result and not a label applied afterwards.
ISO_PU = {"Ds_IsoR04_sumChargedHadronPt", "Ds_IsoR04_sumNeutralHadronEt",
          "Ds_IsoR04_sumPhotonEt", "Ds_IsoR04_sumPUPt", "Ds_IsoR04_RelIso",
          "Ds_IsoR03_RelIso", "mu_IsoR04_RelIso", "num_jets", "num_jets_pt30",
          "jet_ht", "jet_leadpt", "n_b_jets_loose_dR04", "n_b_jets_medium_dR04",
          "Ds_nearestjet_dR", "Ds_nearestjet_bscore", "met_pt", "met_sig", "met_sumEt"}
VTX_IP = {"Ds_prob", "dira", "Ds_ipsig", "Ds_ipxysig", "Ds_FDsig", "Ds_FDxysig",
          "Kp_ipsig", "Km_ipsig", "pi_ipsig", "mu_ipsig",
          "Kp_ipxysig", "Km_ipxysig", "pi_ipxysig"}
GROUP = {"isolation / jets / MET": ("#C0392B", ISO_PU),
         "vertex & impact parameter": ("#E67E22", VTX_IP),
         "kinematics & angles (control)": ("#2C7BB6", set(F) - ISO_PU - VTX_IP)}


def group_of(f):
    for name, (_c, members) in GROUP.items():
        if f in members:
            return name
    return "kinematics & angles (control)"


def grab(name, is_data):
    br = list(dict.fromkeys(F + ["Ds_mass", "era", "num_os_mu_tight",
                                 "num_ss_mu_tight", common.B_VETO]
                            + ([] if is_data else ["Lumi_reweight"])))
    a = uproot.open(f"{SAMPLES}/{name}.root")["Events"].arrays(br, library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a[common.B_VETO] == 0)
         & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALF))       # sideband on BOTH sides
    X = np.stack([a[f][m].astype(np.float64) for f in F], axis=1)
    X, ok = common.preprocess(X)
    w = np.ones(int(ok.sum())) if is_data else a["Lumi_reweight"][m][ok].astype(float)
    return X[ok], a["era"][m][ok], w


print("loading data ...");      Xd, ed, wd = grab("Data", True)
print("loading W+jets MC ..."); Xm, em, wm = grab("WtoLNu2Jets", False)
# NOTE: WtoLNu2Jets in Samples/ is NOT pileup-reweighted -- it has no nTrueInt yet, so
# make_samples.py refuses to reweight it. That is exactly the state this figure
# documents: the diagnosis, before the correction.


def wmean(x, w):
    return float(np.sum(w * x) / np.sum(w))


rows = []
for i, f in enumerate(F):
    d1, d2 = Xd[ed == 1, i], Xd[ed == 2, i]
    s0 = np.std(d1)
    if s0 <= 0:
        continue
    dd = (np.mean(d2) - np.mean(d1)) / s0
    dm = (wmean(Xm[em == 2, i], wm[em == 2]) - wmean(Xm[em == 1, i], wm[em == 1])) / s0
    rows.append((f, dd, dm, dd - dm))
rows.sort(key=lambda r: abs(r[3]))            # smallest residual at the bottom

print(f"\n{'feature':<32}{'d(data)':>9}{'d(MC)':>9}{'resid':>9}")
for f, dd, dm, r in reversed(rows):
    print(f"{f:<32}{dd:>9.3f}{dm:>9.3f}{r:>9.3f}")

# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(21, 13))
gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1.0], height_ratios=[1.0, 1.0],
                      wspace=0.28, hspace=0.32)
axL = fig.add_subplot(gs[:, 0])
axR1 = fig.add_subplot(gs[0, 1])
axR2 = fig.add_subplot(gs[1, 1])

# ---- left: the ranked shifts ----
y = np.arange(len(rows))
for j, (f, dd, dm, r) in enumerate(rows):
    c = GROUP[group_of(f)][0]
    axL.plot([dm, dd], [j, j], color=c, lw=2.4, alpha=0.55, zorder=1,
             solid_capstyle="round")
    axL.plot(dm, j, marker="s", ms=7, mfc="white", mec=c, mew=1.8, zorder=2)
    axL.plot(dd, j, marker="o", ms=8, color=c, zorder=3)
axL.axvline(0, color="black", lw=1.0, alpha=0.6)
axL.set_yticks(y)
axL.set_yticklabels([r[0] for r in rows], fontsize=10)
for t in axL.get_yticklabels():
    t.set_color(GROUP[group_of(t.get_text())][0])
axL.set_ylim(-1, len(rows))
axL.set_xlabel(r"preEE $\to$ postEE shift  [preEE data $\sigma$]", fontsize=14)
axL.set_title("open square = MC     filled circle = data\n"
              "the GAP between them is the mismodelling", fontsize=14, pad=14)
h = ([plt.Line2D([], [], color=c, lw=3, label=n) for n, (c, _m) in GROUP.items()]
     + [plt.Line2D([], [], marker="o", color="0.35", ls="none", ms=8, label="data"),
        plt.Line2D([], [], marker="s", color="0.35", ls="none", ms=7, mfc="white",
                   mew=1.8, label="W+jets MC")])
axL.legend(handles=h, fontsize=11, loc="lower right", framealpha=0.95)
axL.grid(axis="x", alpha=0.25)

# ---- top right: nTrueInt profiles, the quantity the top feature is a proxy for ----
mc = uproot.open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_H_phipimu/"
                 "HtoDsW_2022.root")["Hphipimu/nTrueInt_all"]
mv, me = mc.values(), mc.axis().edges()
prof = {}
for era, lab, col in ((1, "2022 preEE data", "#2C7BB6"), (2, "2022EE postEE data", "#C0392B")):
    hh = uproot.open(PU_PROFILE[era])[PU_PROFILE_HIST]
    v, e = hh.values(), hh.axis().edges()
    prof[era] = (v / v.sum(), e, lab, col)
n = min(len(mv), len(prof[1][0]))
mcn = mv[:n] / mv[:n].sum()
cen = 0.5 * (me[:n] + me[1:n + 1])
for era in (1, 2):
    v, e, lab, col = prof[era]
    hep.histplot(v[:n], me[:n + 1], ax=axR1, histtype="step", lw=2.4, color=col,
                 label=f"{lab}   $\\langle$PU$\\rangle$={np.sum(v[:n]*cen)/v[:n].sum():.1f}")
hep.histplot(mcn, me[:n + 1], ax=axR1, histtype="fill", alpha=0.35, color="seagreen",
             label=f"MC (Run3Summer22)   $\\langle$PU$\\rangle$={np.sum(mcn*cen):.1f}")
axR1.set_xlim(0, 80)
axR1.set_xlabel("nTrueInt", fontsize=13)
axR1.set_ylabel("normalised", fontsize=13)
axR1.legend(fontsize=10.5, loc="upper right")
axR1.set_title("what the top feature is a proxy for:\n"
               "one MC profile, two different data eras", fontsize=13, pad=10)

# ---- bottom right: the two extremes, before/after, as a bar ----
top = sorted(rows, key=lambda r: -abs(r[3]))[:6]
ctrl = [r for r in rows if group_of(r[0]).startswith("kinematics")][:4]
show = top + ctrl
yy = np.arange(len(show))[::-1]
axR2.barh(yy, [r[1] for r in show], height=0.36, color="0.25", label="data")
axR2.barh(yy - 0.38, [r[2] for r in show], height=0.36, color="seagreen",
          label="W+jets MC")
axR2.set_yticks(yy - 0.19)
axR2.set_yticklabels([r[0] for r in show], fontsize=9.5)
for t, r in zip(axR2.get_yticklabels(), show):
    t.set_color(GROUP[group_of(r[0])][0])
axR2.axvline(0, color="black", lw=1.0)
axR2.set_xlabel(r"preEE $\to$ postEE shift  [preEE data $\sigma$]", fontsize=13)
axR2.legend(fontsize=11, loc="lower right")
axR2.set_title("the six largest residuals, and four controls", fontsize=13, pad=10)

hep.cms.label(ax=axL, rlabel="2022 + 2022EE (13.6 TeV)", data=True,
              text="Preliminary", fontsize=12)
out = common.fig_path("features", "era_feature_shift.png")
plt.savefig(out, dpi=110, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved {out}")

# ---------------------------------------------------------------------------
# Second figure: the distributions themselves, for the features with the largest
# data-vs-MC era residual. Four curves per panel, each normalised to unit area --
# the question is shape, and the two eras differ in luminosity by 3.3x anyway:
#
#     data preEE / data postEE   (solid)   -> these separate
#     MC   preEE / MC   postEE   (dashed)  -> these lie on top of each other
#
# Raw values, not the preprocessed/log ones used for the ranking, so the axes are the
# physical quantity. Re-read here for only the plotted branches rather than carrying a
# second full 43-column copy of both samples in memory.
# ---------------------------------------------------------------------------
plot_feats = [r[0] for r in sorted(rows, key=lambda r: -abs(r[3]))[:8]]


def grab_raw(name, feats):
    br = list(dict.fromkeys(list(feats) + ["Ds_mass", "era", "num_os_mu_tight",
                                           "num_ss_mu_tight", common.B_VETO,
                                           "Lumi_reweight"]))
    a = uproot.open(f"{SAMPLES}/{name}.root")["Events"].arrays(br, library="np")
    m = ((a["num_os_mu_tight"] == 1) & (a["num_ss_mu_tight"] == 0)
         & (a[common.B_VETO] == 0)
         & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALF))
    return ({f: a[f][m] for f in feats}, a["era"][m], a["Lumi_reweight"][m])


print("\nre-reading raw values for the distribution figure ...")
Rd, Rde, _ = grab_raw("Data", plot_feats)
Rm, Rme, Rmw = grab_raw("WtoLNu2Jets", plot_feats)

NCOL, NROW = 4, 2
fig, axes = plt.subplots(NROW, NCOL, figsize=(6.0 * NCOL, 5.0 * NROW), squeeze=False)
# Colour = era, line = data vs MC. The MC pair is drawn last, in a lighter shade and
# with a long dash, because the whole point is that the two MC curves coincide -- and
# where they also coincide with data preEE (they nearly do, by construction) three
# curves would otherwise hide under one line and the panel would look like it holds two.
STYLE = [("data preEE",   1, "data", "#0B4F8A", "-",              3.2, 1.00, 2),
         ("data postEE",  2, "data", "#A01F1F", "-",              3.2, 1.00, 2),
         ("MC preEE",     1, "mc",   "#5FA8DC", (0, (7, 2.5)),    2.2, 0.95, 3),
         ("MC postEE",    2, "mc",   "#E8836F", (0, (2.5, 2.5)),  2.2, 0.95, 3)]

for j, f in enumerate(plot_feats):
    ax = axes[j // NCOL][j % NCOL]
    x = Rd[f]
    fin = np.isfinite(x) & (x > -100)
    uniq = np.unique(x[fin])
    if len(uniq) <= 12 and np.allclose(uniq, np.round(uniq)):
        bins = np.arange(uniq.min() - 0.5, uniq.max() + 1.5, 1.0)
    else:
        lo, hi = np.percentile(x[fin], [0.5, 99.0])
        bins = np.linspace(lo, hi, 41)

    drawn = []
    for lab, era, kind, col, ls, lw, al, zo in STYLE:
        if kind == "data":
            h, _ = np.histogram(Rd[f][Rde == era], bins=bins)
        else:
            km = Rme == era
            h, _ = np.histogram(Rm[f][km], bins=bins, weights=Rmw[km])
        h = h / (h.sum() if h.sum() else 1)
        drawn.append(h)
        hep.histplot(h, bins, ax=ax, histtype="step", lw=lw, ls=ls, color=col,
                     alpha=al, zorder=zo, label=lab)

    # Log y only where one spike would otherwise flatten everything else (the isolation
    # sums and the IP significances all pile up at zero).
    pos = np.concatenate([h[h > 0] for h in drawn])
    if pos.size and pos.max() / np.median(pos) > 20:
        ax.set_yscale("log")
        ax.set_ylim(max(pos.min() * 0.5, 1e-5), pos.max() * 4)
    ax.set_xlabel(f, fontsize=15)
    ax.set_ylabel("normalised", fontsize=13)
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.2)

handles, labels = axes[0][0].get_legend_handles_labels()
fig.legend(handles, labels, ncol=4, fontsize=17, loc="upper center",
           bbox_to_anchor=(0.5, 1.045), frameon=False)
plt.tight_layout()
out2 = common.fig_path("features", "era_feature_distributions.png")
plt.savefig(out2, dpi=110, bbox_inches="tight")
plt.close(fig)
print(f"Saved {out2}")
