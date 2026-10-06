"""Does the stage-1 DNN prefer the charge configuration it was NOT trained against?

    python3 compare_os_ss_retention.py [key ...]      (default: ssdata oswjets osdata)

ssdata (background = SS data) keeps OS-configured W+jets combinatorial (GenTag mode 0)
~2.3x more than SS-configured at the same cut (2026-10-06). oswjets / osdata are the
same network trained against an OS background; if the preference is "reject whatever
charge configuration the training background had", their eps_OS/eps_SS flips below 1.

Per model: score_dsmu thresholds at eps_sig = 90/70/50% (that model's own
stage1_threshold), then eps_OS/eps_SS for
  W+jets OS vs SS, per GenTag mode group      (MC weights; error from sum w^2)
  data OS vs SS, sideband windows only        (Poisson; OS data is blinded in the Ds window)
in the low (<1.830) and high (>2.008) Ds_mass windows. ODD events only everywhere: every
model here trained on the even half of its background, so odd is out-of-sample for all.
"""
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "DNN"))
import common  # noqa: E402

KEYS = sys.argv[1:] or ["ssdata", "oswjets", "osdata"]
EFFS = (0.90, 0.70, 0.50)
WINDOWS = {"low <1.830": (0.0, 1.830), "high >2.008": (2.008, 9.0)}
MODE_GROUPS = {"mode 0 (comb.)": [0], "mode 8": [8], "mode 1-5 (true charm)": [1, 2, 3, 4, 5],
               "all modes": list(range(11))}
WANT = ["Ds_mass", "weight", "score_dsmu", "event", "mu_charge", "pi_charge", *common.SELECTION]


def load(key, name, ss, mc):
    path = f"{common.out_dir_for(key, 'DNN')}/Prediction/{name}.root"
    cols = list(dict.fromkeys(WANT + (["cand_gen_mode"] if mc else [])))
    with uproot.open(path)[common.TREE] as t:
        a = t.arrays(cols, library="np")
    common.check_charge_config(a, ss, name)
    m = common.selection_mask(a) & ((a["event"] % 2) == 1)
    return {k: v[m] for k, v in a.items()}


def eff(a, base, thr, mc):
    w = a["weight"] if mc else np.ones(len(a["weight"]))
    k = base & (a["score_dsmu"] > thr)
    n0, n = w[base].sum(), w[k].sum()
    if n0 <= 0 or n <= 0:
        return np.nan, np.nan
    e = n / n0
    err = (np.sqrt((w[k] ** 2).sum()) if mc else np.sqrt(k.sum())) / n0
    return e, err


def ratio(eo, so):
    (a, ea), (b, eb) = eo, so
    r = a / b
    return r, abs(r) * np.hypot(ea / a, eb / b)


for key in KEYS:
    cfg = common.resolve(key)
    thr = {e: common.stage1_threshold({"key": key, "tool": "DNN", "eff_sig": e})[0] for e in EFFS}
    print(f"\n######## {key}  (trained against: {cfg['label']})   thresholds "
          + "  ".join(f"eps{int(100 * e)}={t:.4f}" for e, t in thr.items()))
    S = {("wj", "OS"): load(key, "WtoLNu2Jets_merged", False, True),
         ("wj", "SS"): load(key, "WtoLNu2Jets_SS_merged", True, True),
         ("data", "OS"): load(key, "Data_OS_merged", False, False),
         ("data", "SS"): load(key, "Data_SS_merged", True, False)}
    for wn, (lo, hi) in WINDOWS.items():
        print(f"  --- {wn}: eps_OS / eps_SS at " + " / ".join(f"eps{int(100 * e)}" for e in EFFS))
        rows = [(f"W+jets {g}", "wj", modes) for g, modes in MODE_GROUPS.items()] + [("data", "data", None)]
        for label, src, modes in rows:
            out = []
            for e in EFFS:
                r = []
                for side in ("OS", "SS"):
                    a = S[(src, side)]
                    base = (a["Ds_mass"] >= lo) & (a["Ds_mass"] < hi)
                    if modes is not None:
                        base &= np.isin(a["cand_gen_mode"], modes)
                    r.append(eff(a, base, thr[e], src == "wj"))
                v, ev = ratio(*r)
                out.append(f"{v:6.2f}±{ev:<5.2f}")
            print(f"    {label:28s} " + "  ".join(out))
