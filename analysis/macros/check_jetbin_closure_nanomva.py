"""Re-measure the 2026-08-11 jet-bin closure finding (project_wjets_background_incomplete
memory) on the post-migration NanoMVA/Samples/ production, to check whether it still
holds after the POG corrections + NanoAOD switch.

Original (pre-migration, private tuplizer) OS b-veto CR data/W+jets ratio by
N_jet(pT>30): 0->1.108, 1->1.138, 2->1.787, 3->2.762, >=4->4.761 -- monotonically
rising, read as data containing a multi-jet (largely tt-like, per the b-tagged-region
cross-check) component W+jets MC does not model.

Same definition here: OS + b-veto (n_bjet_corr==0, NOT the jet count itself) + Ds mass
sideband (CR; SR is blind in OS data by construction), jet count from
num_jets_pt30_corr (the corrected-pt >=30 GeV jet count -- the direct NanoMVA
successor of the old N_jet(pT>30) branch). W+jets uses its absolute `weight`
(sigma*L/sum(w)); OS data is unweighted counts.

    source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh
    conda activate ds_dnn
    python3 check_jetbin_closure_nanomva.py
"""
import sys

import numpy as np
import uproot

sys.path.insert(0, "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/analysis/DNN")
import common

M_DS, SR_HALF = 1.96835, 0.040
SEL = [*common.SELECTION]


def ossel(a):
    # common.selection_mask (the 2026-09-21 selection), NOT a hand-written copy: the old
    # num_os/ss_mu_tight condition was an uncalibrated soft-in-jet-muon tag, not a muon
    # selection, and NanoMVA/Samples/ already has the corrected cuts baked in. See common.py.
    return common.selection_mask(a)


def load(name, extra):
    a = uproot.open(f"{common.SAMPLES_DIR}/{name}.root")["Events"].arrays(
        SEL + ["Ds_mass", "num_jets_pt30_corr"] + extra, library="np")
    m = ossel(a) & (np.abs(a["Ds_mass"] - M_DS) >= SR_HALF)
    return {k: v[m] for k, v in a.items()}


d = load("Data_OS_merged", [])
w = load("WtoLNu2Jets_merged", ["weight"])
w_w = w["weight"]

print(f"OS data (CR) {len(d['Ds_mass']):,} rows | W+jets (CR) {len(w['Ds_mass']):,} "
      f"rows, {w_w.sum():,.0f} weighted\n")

njet_d = np.clip(d["num_jets_pt30_corr"], 0, 4)
njet_w = np.clip(w["num_jets_pt30_corr"], 0, 4)

print(f"{'N_jet(pT>30)':>14s} {'OS data':>12s} {'W+jets (wtd)':>14s} {'ratio':>10s}")
for j in range(5):
    label = f"{j}" if j < 4 else ">=4"
    nd = int((njet_d == j).sum())
    mw = w_w[njet_w == j]
    nw, nw2 = mw.sum(), (mw ** 2).sum()
    r = nd / nw if nw > 0 else np.nan
    err = r * np.sqrt(1 / max(nd, 1) + nw2 / nw ** 2) if nw > 0 else np.nan
    print(f"{label:>14s} {nd:12,} {nw:14,.1f} {r:7.3f}+-{err:.3f}")

print("\noriginal (pre-migration): 0->1.108  1->1.138  2->1.787  3->2.762  >=4->4.761")
