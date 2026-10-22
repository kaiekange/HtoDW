"""JME study on the Z->mumu CR (2026-10-22): JME-corrected jet/MET columns as a row-aligned
friend of the existing zmumu skim. The event set is NOT changed -- selection (b-veto, jet
veto map, ...) stays on nano jets exactly as in the skim; only the plotted jet/MET quantities
are recomputed.

    python3 jme_zmumu_study.py submit <group> [--dry-run]        (lxplus, one job per zmumu chunk)
    python3 jme_zmumu_study.py worker <chunk list.txt> <out.root>  (inside cmssw-el9 / CMSSW_15_0_5)
    python3 jme_zmumu_study.py merge <group>                       (inside cmssw-el9 / CMSSW_15_0_5)
    <group>: Data_OS_merged | TT_merged | WtoLNu2Jets_merged | DY_merged

Corrections: Corrections/corrections/jme_direct.py (JME POG JSONs 2026-06-05, JER-Smearing
2026-09-18; data residual at the event's own run). Variants (column suffix):
  _jer       JES + JER (MC) jets; Type-1 PuppiMET with JER propagated (MC)
  _nojer     same jets; Type-1 PuppiMET with JES only (data: identical to _jer)
  _jer_xy, _nojer_xy   the above + MET xy correction (PuppiMET, nom)
  jets: jme_{num_jets_pt30, jet_leadpt, jet_ht_pt20} with the nano tight-ID bit (as
  apply_corrections.py: tight ID, no eta cut), and *_idfix with the NanoAODv12 jet-ID bugfix
  recipe (JetID13p6TeV). MET does not depend on jet ID.
W_tmass = corrected best muon (skim mu_pt_corr, nano mu_phi) + that MET.
Fail-loud checks per chunk: every skim row matched exactly once by (run, lumi, event), and
num_jets_pt30 rebuilt from nano Jet_pt/Jet_jetId equals the skim's num_jets_pt30_corr row by row.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import make_samples_nano as msn

STUDY_DIR = f"{msn.OUT_DIR}/CR/zmumu_jme"
GROUPS = ("Data_OS_merged", "TT_merged", "WtoLNu2Jets_merged", "DY_merged")
NULL = -777.0
VARIANTS = ("jer", "nojer", "jer_xy", "nojer_xy")
LOGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "condor_logs", "jme_zmumu_study")


def out_path(era, sample, k):
    return f"{STUDY_DIR}/skim/{era}/{sample}/chunk_{k}.root"


def chunk_jobs(group):
    """(list, out) for every zmumu chunk of the group whose skim root exists."""
    jobs = []
    for era, sample in msn.GROUPS[group]:
        lists = msn.piece_chunks(era, sample)
        if not lists:
            sys.exit(f"no chunk lists for {era}/{sample}")
        for lp in lists:
            k = int(lp.rsplit("_chunk", 1)[1][:-len(".txt")])
            sp = msn.skim_path(era, sample, k, "zmumu")
            with open(sp[:-len(".root")] + ".json") as f:
                m = json.load(f)
            if m["root_written"]:
                jobs.append((lp, out_path(era, sample, k)))
    return jobs


def jet_cols(pt, tight, prefix):
    import awkward as ak
    import numpy as np
    return {f"{prefix}num_jets_pt30": ak.to_numpy(ak.sum(tight & (pt > 30), axis=1)).astype(np.int32),
            f"{prefix}jet_leadpt": ak.to_numpy(ak.fill_none(ak.max(pt[tight], axis=1), NULL)),
            f"{prefix}jet_ht_pt20": ak.to_numpy(ak.sum(ak.where(tight & (pt > 20), pt, 0.0), axis=1))}


def worker(list_path, out):
    import awkward as ak
    import numpy as np
    import uproot
    sys.path.insert(0, os.path.dirname(msn.APPLY_CORRECTIONS))
    import jme_direct as jd

    with open(list_path) as f:
        nano_files = [x.strip() for x in f if x.strip()]
    era, sample, _, is_data = msn.parse_nano_path(nano_files[0])
    k = int(list_path.rsplit("_chunk", 1)[1][:-len(".txt")])
    sp = msn.skim_path(era, sample, k, "zmumu")
    with uproot.open(sp) as f:
        sk = f["Events"].arrays(["run", "luminosityBlock", "event", "mu_pt_corr", "mu_phi",
                                 "num_jets_pt30_corr"], library="np")
    n_sk = len(sk["event"])
    keys = list(zip(sk["run"].tolist(), sk["luminosityBlock"].tolist(), sk["event"].tolist()))
    row_of = {key: i for i, key in enumerate(keys)}
    if len(row_of) != n_sk:
        sys.exit(f"{sp}: duplicate (run, lumi, event) keys")

    J = jd.Jme(era, is_data)
    cols = {}
    filled = np.zeros(n_sk, dtype=np.int32)
    br = jd.BRANCHES + ["luminosityBlock"] + ([] if is_data else jd.MC_BRANCHES)
    for nano_path in nano_files:
        with uproot.open(nano_path) as fn:
            t = fn["Events"]
            ids = t.arrays(["run", "luminosityBlock", "event"], library="np")
            hit = np.array([key in row_of for key in zip(ids["run"].tolist(),
                            ids["luminosityBlock"].tolist(), ids["event"].tolist())])
            if not hit.any():
                continue
            a = t.arrays(br)[hit]
        rows = np.array([row_of[(int(r), int(l), int(e))] for r, l, e in
                         zip(ak.to_numpy(a.run), ak.to_numpy(a.luminosityBlock), ak.to_numpy(a.event))])
        filled[rows] += 1
        c = {}
        tight_nano = (a.Jet_jetId & 2) > 0
        c["chk_num_jets_pt30_nano"] = jet_cols(a.Jet_pt, tight_nano, "chk_")["chk_num_jets_pt30"]
        pt = J.jet_pt(a)
        c.update(jet_cols(pt, tight_nano, "jme_"))
        tight_fix, _ = jd.jet_id_v12_fix(a)
        c.update({f"{key}_idfix": v for key, v in jet_cols(pt, tight_fix, "jme_").items()})
        npv = ak.to_numpy(a.PV_npvsGood)
        for prop, name in ((True, "jer"), (False, "nojer")):
            m, p = J.met(a, propagate_jer=prop)
            c[f"jme_met_pt_{name}"], c[f"jme_met_phi_{name}"] = m, p
            mx, px = J.met_xy(m, p, npv)
            c[f"jme_met_pt_{name}_xy"], c[f"jme_met_phi_{name}_xy"] = mx, px
        if not cols:
            for key, v in c.items():
                cols[key] = np.full(n_sk, NULL, dtype=v.dtype)
        for key, v in c.items():
            cols[key][rows] = v

    if not np.all(filled == 1):
        sys.exit(f"{sp}: {int((filled == 0).sum())} skim rows unmatched, "
                 f"{int((filled > 1).sum())} matched more than once")
    bad = cols["chk_num_jets_pt30_nano"] != sk["num_jets_pt30_corr"]
    if bad.any():
        sys.exit(f"{sp}: nano num_jets_pt30 rebuilt != skim on {int(bad.sum())}/{n_sk} rows")
    del cols["chk_num_jets_pt30_nano"]
    for v in VARIANTS:
        cols[f"jme_W_tmass_{v}"] = np.sqrt(2.0 * sk["mu_pt_corr"] * cols[f"jme_met_pt_{v}"]
                                           * (1.0 - np.cos(sk["mu_phi"] - cols[f"jme_met_phi_{v}"])))
    for key in ("run", "luminosityBlock", "event"):
        cols[key] = sk[key]
    for key, v in cols.items():
        if not np.all(np.isfinite(v)):
            sys.exit(f"{sp}: non-finite values in {key}")

    os.makedirs(os.path.dirname(out), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        local = os.path.join(tmp, "chunk.root")
        with uproot.recreate(local) as f:
            f["Events"] = cols
        shutil.copyfile(local, out)
    with uproot.open(out) as f:
        back = f["Events"].arrays(["event"], library="np")["event"]
    if not np.array_equal(back, sk["event"]):
        sys.exit(f"read-back of {out} does not match")
    meta = dict(era=era, sample=sample, list=os.path.basename(list_path), skim=sp, n_rows=n_sk,
                jme_json_dir=jd.JSON_DIR[era], jer_smear_json=jd.JER_SMEAR_JSON)
    with open(out[:-len(".root")] + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"{era}/{sample} chunk {k}: {n_sk} rows -> {out}")


def submit(group, dry_run):
    jobs = chunk_jobs(group)
    print(f"{group}: {len(jobs)} jobs")
    if dry_run:
        return
    for _, o in jobs:
        for old in (o, o[:-len(".root")] + ".json"):
            if os.path.exists(old):
                os.remove(old)
    os.makedirs(LOGDIR, exist_ok=True)
    job_table = f"{LOGDIR}/{group}_jobs.txt"
    with open(job_table, "w") as f:
        f.write("\n".join(f"{lp} {o}" for lp, o in jobs) + "\n")
    sched = f"condor_logs/HTCondor_jme_zmumu_study_{group}.log"
    if os.path.exists(sched):
        os.remove(sched)
    desc = f"""
executable = run_jme_zmumu_study.sh
arguments = $(LIST) $(OUT)
log = {sched}
JobBatchName = jme_zmumu_{group}
request_cpus = 1
request_memory = 4G
request_disk = 2G
+JobFlavour = "longlunch"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue LIST,OUT from {job_table}
"""
    subprocess.run(["condor_submit"], input=desc.encode(), check=True)


def merge(group):
    import numpy as np
    import uproot
    jobs = chunk_jobs(group)
    missing = [o for _, o in jobs if not os.path.exists(o[:-len(".root")] + ".json")]
    if missing:
        sys.exit(f"merge refused: {len(missing)} chunks missing, first {missing[0]}")
    parts = []
    for _, o in jobs:
        with uproot.open(o) as f:
            parts.append(f["Events"].arrays(library="np"))
    cols = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    with uproot.open(f"{msn.region_dir('zmumu')}/{group}.root") as f:
        ref = f["Events"].arrays(["run", "luminosityBlock", "event"], library="np")
    for key in ref:
        if not np.array_equal(ref[key], cols[key]):
            sys.exit(f"merge refused: {key} not row-aligned with the zmumu {group}")
    target = f"{STUDY_DIR}/{group}.root"
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        local = os.path.join(tmp, "m.root")
        with uproot.recreate(local) as f:
            f["Events"] = cols
        shutil.copyfile(local, target)
    print(f"wrote {target} ({len(cols['event'])} rows, row-aligned with zmumu {group})")


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == "worker":
        worker(a[1], a[2])
    elif len(a) in (2, 3) and a[0] == "submit" and a[1] in GROUPS and a[2:] in ([], ["--dry-run"]):
        submit(a[1], dry_run=len(a) == 3)
    elif len(a) == 2 and a[0] == "merge" and a[1] in GROUPS:
        merge(a[1])
    else:
        sys.exit(__doc__)
