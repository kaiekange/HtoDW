"""JME step-by-step study on the Z->mumu CR (2026-11-06): which JME step makes the MC MET harder.
Output is a row-aligned friend of the current zmumu skim. The event set is the skim's (selection
on jme_jvm_veto / jme_n_bjet) and is the same for every step; only jet/MET quantities change.

    python3 jme_zmumu_study.py submit <group> [--dry-run]        (lxplus, one job per zmumu chunk)
    python3 jme_zmumu_study.py worker <chunk list.txt> <out.root>  (inside cmssw-el9 / CMSSW_15_0_5)
    python3 jme_zmumu_study.py merge <group>                       (inside cmssw-el9 / CMSSW_15_0_5)
    <group>: Data_OS_merged | DY_merged | TT_merged   (W+jets < 1% of the CR in every MET bin)

Steps, built from Corrections/corrections/jme_direct.py pieces (JES L1L2L3, data L2L3Residual at
the event's run; JER = Jme.jer, the official v2 smearing):
  jets  j0 nano Jet_pt | j1 re-JES | j2 j1 + MC JER smear (= the analysis default *_corr)
        columns j<s>_{num_jets_pt30, jet_leadpt, jet_ht_pt20}, nano tight-ID bit.
  MET   (Type-1 PuppiMET from RawPuppiMET; MET jets = CorrT1METJet + AK4 Jet, muon-subtracted,
        kept if pt>15, |eta|<5.2, EM fraction<0.9; columns m<s>_{pt,phi})
        m0      nano PuppiMET                 m0xy  m0 + MET xy
        m1      Type-1 with re-JES (also the MET of step j2: JER on jets only)
        m3      m1 + JER on the AK4 Jet entries only
        m4      m1 + JER on all MET jets (AK4 + CorrT1)   m5  m4 + MET xy (= analysis default)
        diagnostics of m4: m4match (JER only on gen-matched jets, stochastic ones at JES pt),
        m4nocross (15 GeV threshold decided on the JES pt), m4jf (AK4 entries smeared with the
        jet-collection factor pt_j2/pt_j1 instead of re-evaluating on the muon-subtracted pt);
        counts ncross_{up,dn}_{ct1,ak4}: MET jets that cross 15 GeV only after / before smearing.
  npv = PV_npvsGood (input of the MET xy correction).
  data: no JER, so j2 = j1 and m3 = m4 = m4* = m1, m5 = m1 + xy (data parameters).
Fail-loud checks per chunk: every skim row matched once by (run, lumi, event); j0 == skim *_nano,
j2 == skim *_corr, m0 == skim met_*_nano, m5 == skim met_*_corr; for gen-unmatched jets the JER
factor does not depend on the gen jet (checks the replicated match condition).
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import make_samples_nano as msn

STUDY_DIR = f"{msn.OUT_DIR}/CR/zmumu_jme_steps"
GROUPS = ("Data_OS_merged", "DY_merged", "TT_merged")
NULL = -777.0
MET_STEPS = ("m0", "m0xy", "m1", "m3", "m4", "m4match", "m4nocross", "m4jf", "m5")
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


def gen_match(jd, J, pt, eta, phi, rho, gen_idx, gen_pt, gen_eta, gen_phi):
    """The gen-match condition of Jme.jer, replicated (checked against Jme.jer in step_columns)."""
    import awkward as ak
    import numpy as np
    rho_j = ak.broadcast_arrays(rho, eta)[0]
    reso = jd._flat_eval(J.reso, eta, pt, rho_j)
    valid = (gen_idx >= 0) & (gen_idx < ak.num(gen_pt))
    safe = ak.where(valid, gen_idx, 0)

    def pad(x, v):
        return ak.fill_none(ak.pad_none(x, 1, axis=1), v)
    gpt, geta, gphi = pad(gen_pt, 0.0)[safe], pad(gen_eta, 999.0)[safe], pad(gen_phi, 0.0)[safe]
    dphi = np.arctan2(np.sin(phi - gphi), np.cos(phi - gphi))
    dr = np.sqrt((eta - geta) ** 2 + dphi ** 2)
    return valid & (dr < jd.JER_MATCH_DR) & (abs(pt - gpt) < 3.0 * reso * pt)


def step_columns(jd, J, a, is_data):
    """{column: numpy array} for every jet/MET step (see the module docstring)."""
    import awkward as ak
    import numpy as np
    rho, run, ev = a["Rho_fixedGridRhoFastjetAll"], a["run"], a["event"]
    tight = (a["Jet_jetId"] & 2) > 0
    c = {}
    # ---- jets ----
    _, j1 = J.jes(a["Jet_pt"] * (1.0 - a["Jet_rawFactor"]), a["Jet_area"], a["Jet_eta"], rho, run)
    j2 = j1 if is_data else J.jer(j1, a["Jet_eta"], a["Jet_phi"], rho, ev, a["Jet_genJetIdx"],
                                   a["GenJet_pt"], a["GenJet_eta"], a["GenJet_phi"])
    for s, pt in (("j0", a["Jet_pt"]), ("j1", j1), ("j2", j2)):
        c.update(jet_cols(pt, tight, f"{s}_"))
    # ---- MET jets, per collection: (raw, musub, eta, phi, area, em) ----
    coll = {"ct1": (a["CorrT1METJet_rawPt"], a["CorrT1METJet_muonSubtrFactor"], a["CorrT1METJet_eta"],
                    a["CorrT1METJet_phi"], a["CorrT1METJet_area"], ak.zeros_like(a["CorrT1METJet_eta"])),
            "ak4": (a["Jet_pt"] * (1.0 - a["Jet_rawFactor"]), a["Jet_muonSubtrFactor"], a["Jet_eta"],
                    a["Jet_phi"], a["Jet_area"], a["Jet_chEmEF"] + a["Jet_neEmEF"])}
    if not is_data:
        gen = (a["GenJet_pt"], a["GenJet_eta"], a["GenJet_phi"])
        ge, gp = a["GenJet_eta"], a["GenJet_phi"]   # CorrT1: closest GenJet within dR<0.2 (as jme_direct)
        ce, cp = a["CorrT1METJet_eta"], a["CorrT1METJet_phi"]
        dph = np.arctan2(np.sin(cp[:, :, None] - gp[:, None, :]), np.cos(cp[:, :, None] - gp[:, None, :]))
        dr2 = (ce[:, :, None] - ge[:, None, :]) ** 2 + dph ** 2
        amin = ak.fill_none(ak.argmin(dr2, axis=2), -1)
        mind = ak.fill_none(ak.min(dr2, axis=2), 1e9)
        gidx = {"ct1": ak.values_astype(ak.where(mind < jd.JER_MATCH_DR ** 2, amin, -1), np.int64),
                "ak4": ak.values_astype(a["Jet_genJetIdx"], np.int64)}
    M = {}
    for name, (raw, musub, eta, phi, area, em) in coll.items():
        l1, jes = J.jes(raw * (1.0 - musub), area, eta, rho, run)
        m = dict(l1=l1, jes=jes, eta=eta, phi=phi, ok=(abs(eta) < 5.2) & (em < 0.9))
        if is_data:
            m.update(jer=jes, match=ak.zeros_like(jes, dtype=bool), jf=jes)
        else:
            m["jer"] = J.jer(jes, eta, phi, rho, ev, gidx[name], *gen)
            m["match"] = gen_match(jd, J, jes, eta, phi, rho, gidx[name], *gen)
            nogen = J.jer(jes, eta, phi, rho, ev, ak.full_like(gidx[name], -1), *gen)
            if not ak.all(ak.flatten(ak.where(m["match"], True, m["jer"] == nogen))):
                sys.exit(f"replicated gen match disagrees with Jme.jer ({name})")
            if name == "ak4":
                if not ak.all(ak.flatten(ak.where(musub == 0, (jes == j1) & (m["jer"] == j2), True))):
                    sys.exit("MET AK4 entries without muon != jet collection pt")
                m["jf"] = jes * (j2 / j1)
            else:
                m["jf"] = m["jer"]
        M[name] = m
        c[f"ncross_up_{name}"] = ak.to_numpy(ak.sum(m["ok"] & (jes <= 15.0) & (m["jer"] > 15.0), axis=1)).astype(np.int16)
        c[f"ncross_dn_{name}"] = ak.to_numpy(ak.sum(m["ok"] & (jes > 15.0) & (m["jer"] <= 15.0), axis=1)).astype(np.int16)

    def cat(key):
        return ak.concatenate([M["ct1"][key], M["ak4"][key]], axis=1)
    l1, jes, jer, ok, phi = cat("l1"), cat("jes"), cat("jer"), cat("ok"), cat("phi")
    raw_px = a["RawPuppiMET_pt"] * np.cos(a["RawPuppiMET_phi"])
    raw_py = a["RawPuppiMET_pt"] * np.sin(a["RawPuppiMET_phi"])

    def type1(pt_after, pt_thr):
        dpt = ak.where(ok & (pt_thr > 15.0), pt_after - l1, 0.0)
        px = ak.to_numpy(raw_px - ak.sum(dpt * np.cos(phi), axis=1))
        py = ak.to_numpy(raw_py - ak.sum(dpt * np.sin(phi), axis=1))
        return np.hypot(px, py), np.arctan2(py, px)
    npv = ak.to_numpy(a["PV_npvsGood"])
    c["npv"] = npv.astype(np.int32)
    met = {"m0": (ak.to_numpy(a["PuppiMET_pt"]).astype(np.float64), ak.to_numpy(a["PuppiMET_phi"]).astype(np.float64))}
    met["m0xy"] = J.met_xy(*met["m0"], npv)
    met["m1"] = type1(jes, jes)
    met["m3"] = type1(ak.concatenate([M["ct1"]["jes"], M["ak4"]["jer"]], axis=1),
                      ak.concatenate([M["ct1"]["jes"], M["ak4"]["jer"]], axis=1))
    met["m4"] = type1(jer, jer)
    sm = ak.where(cat("match"), jer, jes)
    met["m4match"] = type1(sm, sm)
    met["m4nocross"] = type1(jer, jes)
    jf = cat("jf")
    met["m4jf"] = type1(jf, jf)
    met["m5"] = J.met_xy(*met["m4"], npv)
    for s, (pt, ph) in met.items():
        c[f"{s}_pt"], c[f"{s}_phi"] = np.asarray(pt, dtype=np.float64), np.asarray(ph, dtype=np.float64)
    return c


def close(name, x, ref, rtol, phi=False):
    """Fail loud unless x == ref within rtol (phi: wrapped absolute difference < rtol)."""
    import numpy as np
    d = np.abs(np.arctan2(np.sin(x - ref), np.cos(x - ref))) if phi else \
        np.abs(x - ref) / np.maximum(np.abs(ref), 1.0)
    bad = d > rtol
    if bad.any():
        sys.exit(f"closure FAILED {name}: {int(bad.sum())}/{len(x)} rows, max diff {d.max():.3g}")


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
    ref_cols = {"j0_num_jets_pt30": "num_jets_pt30_nano", "j0_jet_leadpt": "jet_leadpt_nano",
                "j0_jet_ht_pt20": "jet_ht_pt20_nano", "j2_num_jets_pt30": "num_jets_pt30_corr",
                "j2_jet_leadpt": "jet_leadpt_corr", "j2_jet_ht_pt20": "jet_ht_pt20_corr",
                "m0_pt": "met_pt_nano", "m0_phi": "met_phi_nano", "m5_pt": "met_pt_corr", "m5_phi": "met_phi_corr"}
    with uproot.open(sp) as f:
        sk = f["Events"].arrays(["run", "luminosityBlock", "event"] + list(ref_cols.values()), library="np")
    n_sk = len(sk["event"])
    keys = list(zip(sk["run"].tolist(), sk["luminosityBlock"].tolist(), sk["event"].tolist()))
    row_of = {key: i for i, key in enumerate(keys)}
    if len(row_of) != n_sk:
        sys.exit(f"{sp}: duplicate (run, lumi, event) keys")

    J = jd.Jme(era, is_data)
    cols = {}
    filled = np.zeros(n_sk, dtype=np.int32)
    br = jd.BRANCHES + ["luminosityBlock", "PuppiMET_pt", "PuppiMET_phi"] + ([] if is_data else jd.MC_BRANCHES)
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
        c = step_columns(jd, J, a, is_data)
        if not cols:
            for key, v in c.items():
                cols[key] = np.full(n_sk, NULL, dtype=v.dtype)
        for key, v in c.items():
            cols[key][rows] = v

    if not np.all(filled == 1):
        sys.exit(f"{sp}: {int((filled == 0).sum())} skim rows unmatched, "
                 f"{int((filled > 1).sum())} matched more than once")
    for mine, theirs in ref_cols.items():
        if mine.endswith("num_jets_pt30"):
            if not np.array_equal(cols[mine], sk[theirs]):
                sys.exit(f"{sp}: closure FAILED {mine} != skim {theirs}")
        else:   # m0 vs the producer's MET: 1e-3 as jme.full_columns; the rest is the same code
            close(f"{sp}: {mine} vs skim {theirs}", cols[mine], sk[theirs].astype(np.float64),
                  1e-3 if mine.startswith("m0") else 1e-6, phi=mine.endswith("_phi"))
    for key in list(cols):
        if cols[key].dtype == np.float64:
            cols[key] = cols[key].astype(np.float32)
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
