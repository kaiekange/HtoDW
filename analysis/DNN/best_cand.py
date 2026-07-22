"""Select the best Ds candidate per event (highest score_dsmu) and keep all
gamma rows for that Ds.

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
cand_id  = all_data["cand_id"]
score    = all_data["score_dsmu"]
n_total  = len(event_id)
print(f"Total rows: {n_total}")

# For each event find the cand_id with the highest score_dsmu.
# Same (event_id, cand_id) always has the same score, so any row suffices.
ev_best_cand = {}   # event_id -> best cand_id
ev_best_score = {}  # event_id -> best score (to break ties)
for i in range(n_total):
    ev = int(event_id[i])
    sc = float(score[i])
    if ev not in ev_best_score or sc > ev_best_score[ev]:
        ev_best_score[ev] = sc
        ev_best_cand[ev]  = int(cand_id[i])

# Keep all rows where cand_id matches the best one for that event
mask = np.array([
    int(cand_id[i]) == ev_best_cand[int(event_id[i])]
    for i in range(n_total)
], dtype=bool)

n_kept = mask.sum()
print(f"Kept rows (all gammas of best Ds): {n_kept}")
print(f"Unique events: {len(ev_best_cand)}")

# Write filtered branches
out_branches = {k: v[mask] for k, v in all_data.items()}
branch_types = {k: v.dtype for k, v in out_branches.items()}

print(f"Writing {out_file}...")
with uproot.recreate(out_file) as f_out:
    f_out.mktree(tree_name, branch_types)
    f_out[tree_name].extend(out_branches)

print(f"Done.")
