"""
For each event, deduplicate by (event_id, cand_id) — same Ds with different
gamma rows has the same score_dsmu — then check whether the highest-score Ds
is the matched one (Ds_match && mu_match).
"""

import numpy as np
import uproot

model_name = "v1"
PRED = f"/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/{model_name}/Prediction"

FILES = [
    (f"{PRED}/Signal_DsstW_H_phipimu_scored.root", "HphipimuDsstW/Events", "DsstW"),
    (f"{PRED}/Signal_H_phipimu_scored.root",        "HphipimuDsW/Events",   "DsW"),
]

BRANCHES = ["event_id", "cand_id", "Ds_match", "mu_match", "score_dsmu", "Ds_pt"]

def eval_selector(event_id_u, score_u, matched_u, label):
    n_events_with_match = 0
    n_best_is_match     = 0
    n_multi_cand        = 0
    n_multi_best_match  = 0
    for ev in np.unique(event_id_u):
        ev_mask  = event_id_u == ev
        ev_score = score_u[ev_mask]
        ev_match = matched_u[ev_mask]
        if not ev_match.any():
            continue
        n_events_with_match += 1
        best_match = ev_match[np.argmax(ev_score)]
        n_cands = ev_mask.sum()
        if n_cands > 1:
            n_multi_cand += 1
            if best_match:
                n_multi_best_match += 1
        if best_match:
            n_best_is_match += 1
    pct     = 100 * n_best_is_match / n_events_with_match
    pct_m   = (100 * n_multi_best_match / n_multi_cand) if n_multi_cand else float("nan")
    print(f"  [{label}] overall: {n_best_is_match}/{n_events_with_match} ({pct:.1f}%)"
          f"  |  multi-Ds: {n_multi_best_match}/{n_multi_cand} ({pct_m:.1f}%)")

for fpath, tree, label in FILES:
    print(f"\n{'='*50}")
    print(f"{label}: {fpath}")

    with uproot.open(fpath)[tree] as t:
        d = t.arrays(BRANCHES, library="np")

    event_id = d["event_id"]
    cand_id  = d["cand_id"]
    matched  = d["Ds_match"].astype(bool) & d["mu_match"].astype(bool)
    score    = d["score_dsmu"]
    ds_pt    = d["Ds_pt"]

    print(f"Total rows: {len(event_id)}")

    # Deduplicate by (event_id, cand_id)
    _, first_idx = np.unique(np.column_stack([event_id, cand_id]), axis=0, return_index=True)
    event_id_u = event_id[first_idx]
    matched_u  = matched[first_idx]
    score_u    = score[first_idx]
    pt_u       = ds_pt[first_idx]

    print(f"Unique (event, Ds) candidates: {len(event_id_u)}")
    eval_selector(event_id_u, score_u, matched_u, "DNN score")
    eval_selector(event_id_u, pt_u,    matched_u, "Ds_pt    ")
