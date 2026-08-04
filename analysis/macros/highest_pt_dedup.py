"""Dedup an H_phipimu (OS or sameSignPi) flat ntuple to one row per event: keep the
candidate with the highest Ds_pt, breaking ties on gamma_pt (rows are per-(Ds,gamma)
candidate: Ds_pt is identical across every gamma row of the same Ds, so when the
winning Ds has several gamma rows they all tie on Ds_pt -- gamma_pt picks the best one).

Works on any H_phipimu(_SS) output, signal or background, OS or SS -- just point
it at the right TDirectory/tree path (the module label differs by cfg, e.g.
"HphipimuDsW/Events" for the signal cfgs, "Hphipimu/Events" for WtoLNu2Jets).

Usage:
    python3 highest_pt_dedup.py <input.root> <tree_path> <output.root>

Example:
    python3 highest_pt_dedup.py \\
        /eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_H_phipimu/HtoDsW.root \\
        HphipimuDsW/Events best_HtoDsW.root
"""

import sys
import numpy as np
import uproot

input_file = sys.argv[1]
tree_path  = sys.argv[2]
out_file   = sys.argv[3]
out_tree   = tree_path.split("/")[-1]

print(f"Reading {input_file} [{tree_path}]...")
with uproot.open(input_file)[tree_path] as t:
    all_data = t.arrays(library="np")

n_total = len(all_data["event_id"])
print(f"Total rows: {n_total}")

# Group by (run_id, lumi_id, event_id) -- event_id alone is not safe across runs/lumis.
run_id, lumi_id, event_id = all_data["run_id"], all_data["lumi_id"], all_data["event_id"]
order = np.lexsort((event_id, lumi_id, run_id))
r, l, e = run_id[order], lumi_id[order], event_id[order]
group_start = np.ones(len(order), dtype=bool)
group_start[1:] = (r[1:] != r[:-1]) | (l[1:] != l[:-1]) | (e[1:] != e[:-1])
starts = np.flatnonzero(group_start)
ends = np.append(starts[1:], len(order))

ds_pt = all_data["Ds_pt"]
gamma_pt = all_data["gamma_pt"]
best = np.empty(len(starts), dtype=np.int64)
for i, (s, en) in enumerate(zip(starts, ends)):
    idx = order[s:en]
    pt = ds_pt[idx]
    top = idx[pt == pt.max()]
    best[i] = top[0] if top.size == 1 else top[np.argmax(gamma_pt[top])]
best.sort()

n_kept = best.size
print(f"Unique events kept: {n_kept}  (dropped {n_total - n_kept} lower-Ds_pt rows)")

out_branches = {k: v[best] for k, v in all_data.items()}
branch_types = {k: v.dtype for k, v in out_branches.items()}

print(f"Writing {out_file} [{out_tree}]...")
with uproot.recreate(out_file) as f_out:
    f_out.mktree(out_tree, branch_types)
    f_out[out_tree].extend(out_branches)

print("Done.")
