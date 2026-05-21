"""Pick the highest-score candidate per event and write a new ROOT file.

Usage:
    python3 best_cand.py <input_scored.root> [tree_name] [output_root]

Output:
    one row per event — the candidate with the highest DNN score.
    Events whose every candidate has score == -777 are dropped.
"""

import sys
import numpy as np
import uproot
from collections import defaultdict

input_file = sys.argv[1]
tree_name  = sys.argv[2] if len(sys.argv) > 2 else "Events"
out_root   = input_file.replace("_scored.root", "_best.root")

# ========= load =========
print(f"Reading {input_file} [{tree_name}]...")
with uproot.open(input_file) as f:
    data = f[tree_name].arrays(library="np")

n_in = len(data["event_id"])
print(f"Loaded {n_in} candidates")

# ========= best-candidate selection =========
event_ids = data["event_id"]
scores    = data["score"]

evt_to_rows = defaultdict(list)
for i, eid in enumerate(event_ids):
    evt_to_rows[eid].append(i)

selected = []
for eid, rows in evt_to_rows.items():
    rows_arr = np.array(rows)
    s = scores[rows_arr]
    valid = s > -777
    if not valid.any():
        continue
    best = rows_arr[np.argmax(np.where(valid, s, -np.inf))]
    selected.append(best)

selected = np.array(sorted(selected))
print(f"Selected {len(selected)} events from {len(evt_to_rows)} unique event IDs")

# ========= write output =========
out_branches = {k: v[selected] for k, v in data.items()}
branch_types = {k: v.dtype for k, v in out_branches.items()}

print(f"Writing {out_root}...")
with uproot.recreate(out_root) as f_out:
    f_out.mktree(tree_name, branch_types)
    f_out[tree_name].extend(out_branches)

print(f"Done. {len(selected)} entries saved to {out_root}")
