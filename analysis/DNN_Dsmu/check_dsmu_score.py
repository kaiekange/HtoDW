"""
For each event, deduplicate by (event_id, cand_id) — same Ds with different
gamma rows has the same score_dsmu — then check whether the highest-score Ds
is the matched one (Ds_match && mu_match).
"""

import numpy as np
import uproot

FILES = [
    ("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_DsstW_H_KKpimu/Signal_DsstW_H_KKpimu_scored.root",
     "HKKpimuDsstW/Events", "DsstW"),
    ("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_H_KKpimu/Signal_H_KKpimu_scored.root",
     "HKKpimuDsW/Events", "DsW"),
]

BRANCHES = ["event_id", "cand_id", "Ds_match", "mu_match", "score_dsmu"]

for fpath, tree, label in FILES:
    print(f"\n{'='*50}")
    print(f"{label}: {fpath}")

    with uproot.open(fpath)[tree] as t:
        d = t.arrays(BRANCHES, library="np")

    event_id  = d["event_id"]
    cand_id   = d["cand_id"]
    ds_match  = d["Ds_match"].astype(bool)
    mu_match  = d["mu_match"].astype(bool)
    score     = d["score_dsmu"]
    matched   = ds_match & mu_match

    print(f"Total rows: {len(event_id)}")

    # Deduplicate: keep first occurrence of each (event_id, cand_id) pair.
    # Score is the same for all rows of the same Ds candidate regardless of gamma.
    _, first_idx = np.unique(np.column_stack([event_id, cand_id]), axis=0, return_index=True)
    event_id_u = event_id[first_idx]
    cand_id_u  = cand_id[first_idx]
    score_u    = score[first_idx]
    matched_u  = matched[first_idx]

    print(f"Unique (event, Ds) candidates: {len(event_id_u)}")

    # For each event: find cand with highest score, check if it's matched
    n_events_with_match = 0
    n_best_is_match     = 0
    n_multi_cand        = 0
    n_multi_best_match  = 0

    unique_events = np.unique(event_id_u)
    for ev in unique_events:
        ev_mask  = event_id_u == ev
        ev_score = score_u[ev_mask]
        ev_match = matched_u[ev_mask]

        if not ev_match.any():
            continue  # no matched Ds in this event, skip
        n_events_with_match += 1

        best_idx    = np.argmax(ev_score)
        best_match  = ev_match[best_idx]
        n_cands     = ev_mask.sum()

        if n_cands > 1:
            n_multi_cand += 1
            if best_match:
                n_multi_best_match += 1

        if best_match:
            n_best_is_match += 1

    print(f"Events with >= 1 matched Ds:         {n_events_with_match}")
    print(f"  Best score Ds is matched:           {n_best_is_match} / {n_events_with_match}"
          f"  ({100*n_best_is_match/n_events_with_match:.1f}%)")
    print(f"  Multi-Ds events (with match):       {n_multi_cand}")
    print(f"  Multi-Ds: best score is matched:    {n_multi_best_match} / {n_multi_cand}"
          f"  ({100*n_multi_best_match/n_multi_cand:.1f}%)" if n_multi_cand else "  (none)")
