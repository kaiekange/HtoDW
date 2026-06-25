"""Select the best gamma per event (highest score_gamma).
Input: scored file from predict.py (one Ds per event, multiple gamma rows).
Output: one row per event — the best gamma for the best Ds.

Usage:
    python3 best_cand.py <input_scored.root> <tree_name> <output.root>
"""

import sys
import numpy as np
import uproot

input_file = sys.argv[1]
tree_name  = sys.argv[2]
out_file   = sys.argv[3]

print(f"Reading {input_file} [{tree_name}]...")
with uproot.open(input_file)[tree_name] as t:
    all_data = t.arrays(library="np")

event_id = all_data["event_id"]
score    = all_data["score_gamma"]
n_total  = len(event_id)
print(f"Total rows: {n_total}")

# For each event pick the row with highest score_gamma
ev_best_idx = {}   # event_id -> row index of best gamma
ev_best_sc  = {}
for i in range(n_total):
    ev = int(event_id[i])
    sc = float(score[i])
    if ev not in ev_best_sc or sc > ev_best_sc[ev]:
        ev_best_sc[ev]  = sc
        ev_best_idx[ev] = i

keep = np.array(list(ev_best_idx.values()), dtype=int)
keep.sort()
print(f"Events (one row each): {len(keep)}")

out_branches = {k: v[keep] for k, v in all_data.items()}
branch_types = {k: v.dtype for k, v in out_branches.items()}

print(f"Writing {out_file}...")
with uproot.recreate(out_file) as f_out:
    f_out.mktree(tree_name, branch_types)
    f_out[tree_name].extend(out_branches)

print("Done.")
