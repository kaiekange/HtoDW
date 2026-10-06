#!/usr/bin/env python3
"""Merged custom-NanoAOD samples: signal (DsW+DsstW), W+jets GenTag (0J/1J/2J; OS and same-sign-pi SS),
ttbar GenTag (2L2Nu + LNu2Q), DY GenTag (0J/1J/2J), and Muon data OS / SS (Run2022 C-G), each merged over 2022+2022EE into one file.
NanoAOD-pipeline analogue of the old tuplizer-based make_samples.py.

Two modes, driven by condor_submit_make_samples_nano.py:

  worker <list.txt> <out.root> [--correct]
      One condor job = one chunk (<=50 nano files of ONE piece, listed in list.txt).
      Per file: read nano + its apply_corrections.py friend (row i <-> row i), assert
      the join on run/luminosityBlock/event, select, weight, tag; the whole chunk is
      written ONCE (local temp, then copy), so the merged file ends up with one basket
      per branch per chunk instead of one per nano file. Then <out>.json is written
      last: its existence is the success marker the merge checks.
      --correct: first run the Corrections repo's apply_corrections.py on every file of
      the chunk (subprocess, unchanged script, writes the usual friend), then skim. One
      submission instead of two; the friends are still written because MC systematics
      need unselected rows and the skim gets re-run far more often than the corrections.

  merge <group>   (Signal_merged | WtoLNu2Jets_merged | WtoLNu2Jets_SS_merged | TT_merged | DY_merged | Data_OS_merged | Data_SS_merged)
      Reads the chunk lists the submit script wrote (never re-globs), requires every
      chunk's .root+.json, checks that the chunks cover exactly the files on disk and
      share one branch list, hadds them in list order, writes <group>_meta.json (per
      piece sums) and asserts merged entries == sum of selected rows.

SELECTION (friend columns, same combination as every corrections validator):
    pass_mu_nom & pass_muveto & ~jvm_veto_nom & (n_bjet_nom == 0)

WEIGHT (background, i.e. W+jets, ttbar and DY):
    weight = genWeight * norm * w_pu_nom * w_mu_nom * w_btag_nom
`norm` is the piece's normalization.json value (xsec*L / (sum genWeight * acceptance
correction), Corrections repo normalization.py).

WEIGHT (signal) -- changed 2026-09-17, genWeight dropped:
    weight = norm * w_pu_nom * w_mu_nom * w_btag_nom
Found 7/869,385 post-selection signal rows (49/1.6M pre-selection) with genWeight =
-23.1645 instead of +23.1645 -- a technical two-valued weight with no cross-section
meaning for this LO-ish EvtGen+Pythia generation (unlike the real amcatnloFXFX NLO
subtraction term W+jets has). It crashed DNN.py's WeightedRandomSampler (torch.multinomial
rejects any negative probability) the first time stage-1 was retrained on this sample.
`norm` for signal now divides by genEventCount (plain N_gen), not sum(genWeight), to stay
consistent with genWeight being absent from the numerator too -- see normalization.py.

SIGNAL IS ON THE SAME ABSOLUTE SCALE AS BACKGROUND here -- unlike the old Samples/
(weight_mode="signal_era", relative, anchored at 1.0; see eras.py SIGNAL_ERA_WEIGHT)
summing `weight` gives an expected yield.

DATA: no genWeight/normalization/SFs, so `weight` = 1.0 (present only so MC and data read
the same way), no `norm`. OS data is BLINDED here: rows with |Ds_mass - M_DS| <
SR_HALF_WIDTH (eras.py, 1.96835 +- 0.04) are dropped after the selection, the same
definition as highest_pt_dedup.apply_row_cuts(blind=True). SS data is not blinded.

PROVENANCE BRANCHES: era (eras.py wire format, 1=2022, 2=2022EE, data too); signal file: channel (0=DsW, 1=DsstW);
W+jets files (OS and SS): jet_bin (0/1/2);
ttbar file: tt_decay (0=TTto2L2Nu, 1=TTtoLNu2Q); DY file: jet_bin (0/1/2). Data has no extra tag -- `run` identifies the run period.

CARRIED: every custom analyzer branch (everything not in a stock NanoAOD collection this
analysis never reads, see _STD_PREFIXES) + run/luminosityBlock/event/genWeight/nTrueInt;
from the friend: w_pu_nom/w_mu_nom/w_btag_nom, `norm`, and the corrected nominal
kinematics renamed *_corr (FRIEND_CORR_KIN). The raw nano versions are kept too.
NOT CARRIED: the up/down systematic columns -- not asked for; they stay in the friends.

Runs inside cmssw-el9 + CMSSW_15_0_5 (uproot, awkward, hadd).
"""
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eras import M_DS, SR_HALF_WIDTH, TAG_TO_ERA  # noqa: E402

MC_BASE = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"
DATA_BASE = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Data"
OUT_DIR = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Samples"
SKIM_DIR = f"{OUT_DIR}/skim"
LIST_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "condor_logs", "make_samples_nano", "lists")
FILES_PER_JOB = 50
APPLY_CORRECTIONS = "/afs/cern.ch/work/k/kakang/H2DW/Corrections/corrections/apply_corrections.py"

SIGNAL_CHANNEL = {"Signal_DsW_H_phipimu": 0, "Signal_DsstW_H_phipimu": 1}
WJETS_BIN = {
    "WtoLNu2Jets_0J_H_phipimu_GenTag": 0,
    "WtoLNu2Jets_1J_H_phipimu_GenTag": 1,
    "WtoLNu2Jets_2J_H_phipimu_GenTag": 2,
}
# same-sign-pi W+jets (control region for the OS-minus-SS comparison), same jet_bin tag
WJETS_SS_BIN = {f"{s}_SS": b for s, b in WJETS_BIN.items()}
# ttbar decay-filtered datasets -> tt_decay tag (normalized separately, each with its own xsec)
TT_DECAY = {
    "TTto2L2Nu_H_phipimu_GenTag": 0,
    "TTtoLNu2Q_H_phipimu_GenTag": 1,
}
DY_BIN = {
    "DYto2L2Jets_0J_H_phipimu_GenTag": 0,
    "DYto2L2Jets_1J_H_phipimu_GenTag": 1,
    "DYto2L2Jets_2J_H_phipimu_GenTag": 2,
}
DATA_RUNS = {"2022": ("C", "D"), "2022EE": ("E", "F", "G")}


def data_task(run, same_sign):
    return f"Muon_Run2022{run}_NanoH_phipimu" + ("_SS" if same_sign else "")


# data task -> blind? (OS only)
DATA_TASK_BLIND = {data_task(r, ss): not ss
                   for runs in DATA_RUNS.values() for r in runs for ss in (False, True)}

# group -> ordered pieces (era tag, sample); this order is the row order of the merged file
GROUPS = {
    "Signal_merged": [(e, s) for e in ("2022", "2022EE") for s in SIGNAL_CHANNEL],
    "WtoLNu2Jets_merged": [(e, s) for e in ("2022", "2022EE") for s in WJETS_BIN],
    "WtoLNu2Jets_SS_merged": [(e, s) for e in ("2022", "2022EE") for s in WJETS_SS_BIN],
    "TT_merged": [(e, s) for e in ("2022", "2022EE") for s in TT_DECAY],
    "DY_merged": [(e, s) for e in ("2022", "2022EE") for s in DY_BIN],
    "Data_OS_merged": [(e, data_task(r, False)) for e, runs in DATA_RUNS.items() for r in runs],
    "Data_SS_merged": [(e, data_task(r, True)) for e, runs in DATA_RUNS.items() for r in runs],
}


_STD_PREFIXES = (
    "HLT_", "L1_", "L1Reco_step", "L1simulation_step", "Muon_", "Electron_", "Tau_",
    "FatJet_", "Photon_", "Jet_", "LowPtElectron_", "Flag_", "boostedTau_", "SV_",
    "GenJet", "GenPart", "GenVisTau", "GenDressedLepton", "GenIsolatedPhoton",
    "GenProton", "GenVtx", "SubGenJet", "SubJet_", "SoftActivityJet", "CorrT1METJet",
    "IsoTrack", "OtherPV", "PSWeight", "LHE", "PV_", "Pileup_", "RawPuppiMET", "RawMET",
    "PuppiMET", "MET_", "ChsMET", "DeepMETResolutionTune", "DeepMETResponseTune",
    "TrigObj", "btagWeight", "fixedGridRho", "Rho_", "TkMET", "CaloMET", "FsrPhoton_",
    "HTXS_", "Generator_", "HLTriggerFinalPath", "HLTriggerFirstPath", "BeamSpot_",
    "GenMET_", "bunchCrossing", "genTtbarId", "cand_id", "nCorrT1METJet", "nFatJet",
    "nFsrPhoton", "nGenDressedLepton", "nGenIsolatedPhoton", "nGenJet", "nGenJetAK8",
    "nGenPart", "nGenProton", "nGenVisTau", "nIsoTrack", "nLHEPart", "nLHEPdfWeight",
    "nLHEReweightingWeight", "nLHEScaleWeight", "nLowPtElectron", "nOtherPV", "nPSWeight",
    "nPUInt", "nPhoton", "nSV", "nSoftActivityJet", "nSubGenJetAK8", "nSubJet", "nTau",
    "nTrigObj", "nboostedTau",
    "PPSLocalTrack_", "Proton_", "nPPSLocalTrack", "nProton_",  # data only
)
_ALWAYS = ("run", "luminosityBlock", "event", "genWeight", "nTrueInt")

FRIEND_JOIN_COLS = ["run", "luminosityBlock", "event"]
FRIEND_WEIGHT_COLS = ["w_pu_nom", "w_mu_nom", "w_btag_nom"]
FRIEND_CUT_COLS = ["pass_mu_nom", "pass_muveto", "jvm_veto_nom", "n_bjet_nom"]
FRIEND_CORR_KIN = {
    "mu_pt_nom": "mu_pt_corr",
    "jet_ht_pt20_nom": "jet_ht_pt20_corr",
    "jet_leadpt_nom": "jet_leadpt_corr",
    "met_pt_nom": "met_pt_corr",
    "met_phi_nom": "met_phi_corr",
    "num_jets_pt30_nom": "num_jets_pt30_corr",
    "W_tmass_nom": "W_tmass_corr",
    "n_bjet_nom": "n_bjet_corr",
    "Dsmu_pt_nom": "Dsmu_pt_corr", "Dsmu_eta_nom": "Dsmu_eta_corr",
    "Dsmu_phi_nom": "Dsmu_phi_corr", "Dsmu_mass_nom": "Dsmu_mass_corr",
    "Dsmu_Ht_nom": "Dsmu_Ht_corr",
    "Dsstmu_pt_nom": "Dsstmu_pt_corr", "Dsstmu_eta_nom": "Dsstmu_eta_corr",
    "Dsstmu_phi_nom": "Dsstmu_phi_corr", "Dsstmu_mass_nom": "Dsstmu_mass_corr",
    "Dsstmu_Ht_nom": "Dsstmu_Ht_corr",
}

_IDX_RE = re.compile(r"(?:output|nano)_(\d+)\.root$")


def nano_branch_list(keys):
    return sorted({k for k in keys
                   if k in _ALWAYS or not any(k.startswith(p) for p in _STD_PREFIXES)})


def piece_nano_files(era, sample):
    if sample in DATA_TASK_BLIND:
        # CRAB layout <task>/<timestamp>/<0000|0001|...>/nano_N.root; the friend dir is flat
        files = glob.glob(f"{DATA_BASE}/{era}/Muon/{sample}/*/*/nano_*.root")
        names = [os.path.basename(f) for f in files]
        if len(set(names)) != len(names):
            sys.exit(f"{era}/{sample}: duplicate nano_N basenames across CRAB subdirectories")
    else:
        files = glob.glob(f"{MC_BASE}/{era}/{sample}/output_*.root")
    return sorted(files, key=lambda p: int(_IDX_RE.search(p).group(1)))


def chunk_name(era, sample, k):
    return f"{era}_{sample}_chunk{k}"


def piece_chunks(era, sample):
    """Chunk list files of one piece as written by the submit script, in chunk order."""
    lists = glob.glob(f"{LIST_DIR}/{era}_{sample}_chunk*.txt")
    return sorted(lists, key=lambda p: int(re.search(r"_chunk(\d+)\.txt$", p).group(1)))


def skim_path(era, sample, k):
    return f"{SKIM_DIR}/{era}/{sample}/chunk_{k}.root"


def parse_nano_path(path):
    """-> (era, sample, friend path, is_data); same friend layout as apply_corrections.py."""
    parts = os.path.abspath(path).split(os.sep)
    if "Data" in parts:
        i = len(parts) - 1 - parts[::-1].index("Data")
        era, task = parts[i + 1], parts[i + 3]
        if (era not in TAG_TO_ERA or parts[i + 2] != "Muon" or task not in DATA_TASK_BLIND
                or len(parts) != i + 7 or not _IDX_RE.search(parts[-1])):
            sys.exit(f"not a nano data path: {path}")
        return era, task, os.sep.join(parts[:i + 2] + ["Corrections", task, parts[-1]]), True
    era, sample = parts[-3], parts[-2]
    if parts[-4] != "MC" or era not in TAG_TO_ERA or not _IDX_RE.search(parts[-1]):
        sys.exit(f"not a nano MC path: {path}")
    if not any(sample in t for t in (SIGNAL_CHANNEL, WJETS_BIN, WJETS_SS_BIN, TT_DECAY, DY_BIN)):
        sys.exit(f"sample {sample!r} belongs to no group")
    friend = os.sep.join(parts[:-2] + ["Corrections", sample, parts[-1]])
    return era, sample, friend, False


# ---------------------------------------------------------------- worker
def worker(list_path, out_path, correct=False):
    # numpy/uproot only here and in merge(): the submit script imports this module from
    # lxplus's native python, which has neither.
    import numpy as np
    import uproot

    with open(list_path) as f:
        nano_files = [line.strip() for line in f if line.strip()]
    if not nano_files:
        sys.exit(f"empty list {list_path}")
    pieces = {parse_nano_path(p)[:2] for p in nano_files}
    if len(pieces) != 1:
        sys.exit(f"{list_path} mixes pieces: {sorted(pieces)}")
    era, sample, _, is_data = parse_nano_path(nano_files[0])

    tags = {"era": TAG_TO_ERA[era]}
    if is_data:
        norm = xsec_pb = lumi_pb = None
        blind = DATA_TASK_BLIND[sample]
        friend_cols = FRIEND_JOIN_COLS + FRIEND_CUT_COLS + list(FRIEND_CORR_KIN)
    else:
        with open(f"{MC_BASE}/{era}/Corrections/{sample}/normalization.json") as f:
            nd = json.load(f)
        norm, xsec_pb, lumi_pb = nd["norm"], nd["xsec_pb"], nd["lumi_pb"]
        blind = False
        friend_cols = (FRIEND_JOIN_COLS + FRIEND_WEIGHT_COLS + FRIEND_CUT_COLS
                       + list(FRIEND_CORR_KIN))
        if sample in SIGNAL_CHANNEL:
            tags["channel"] = SIGNAL_CHANNEL[sample]
        elif sample in TT_DECAY:
            tags["tt_decay"] = TT_DECAY[sample]
        elif sample in DY_BIN:
            tags["jet_bin"] = DY_BIN[sample]
        elif sample in WJETS_SS_BIN:
            tags["jet_bin"] = WJETS_SS_BIN[sample]
        else:
            tags["jet_bin"] = WJETS_BIN[sample]

    branches = None
    parts = []
    n_total = n_sel_total = n_blinded = 0
    wsum = 0.0
    for nano_path in nano_files:
        friend_path = parse_nano_path(nano_path)[2]
        if correct:
            subprocess.run([sys.executable, APPLY_CORRECTIONS, nano_path], check=True)
        with uproot.open(nano_path) as fn:
            tn = fn["Events"]
            if branches is None:
                branches = nano_branch_list(tn.keys())
            nano = tn.arrays(branches, library="np")
        with uproot.open(friend_path) as ff:
            fr = ff["Events"].arrays(friend_cols, library="np")
        n = len(nano["event"])
        n_total += n
        if len(fr["event"]) != n:
            sys.exit(f"row count mismatch {nano_path} ({n}) vs {friend_path} ({len(fr['event'])})")
        for c in FRIEND_JOIN_COLS:
            if not np.array_equal(nano[c], fr[c]):
                sys.exit(f"join mismatch on {c}: {nano_path} vs {friend_path}")

        sel = (fr["pass_mu_nom"] & fr["pass_muveto"] & ~fr["jvm_veto_nom"]
               & (fr["n_bjet_nom"] == 0))
        if blind:
            in_sr = sel & (np.abs(nano["Ds_mass"] - M_DS) < SR_HALF_WIDTH)
            n_blinded += int(in_sr.sum())
            sel &= ~in_sr
        n_sel = int(sel.sum())
        if n_sel == 0:
            continue
        if is_data:
            w = np.ones(n_sel, dtype=np.float64)
        else:
            # Signal drops genWeight (see module docstring, WEIGHT (signal)): the
            # ~23.1645-magnitude technical weight isn't a real cross-section factor for
            # this generation, and a rare sign flip (7/869,385 post-selection) crashes
            # torch's WeightedRandomSampler downstream. normalization.py's `norm` for
            # signal is computed to match (genEventCount denominator, not sum(genWeight)).
            gw = (np.ones(n_sel, dtype=np.float64) if sample in SIGNAL_CHANNEL
                  else nano["genWeight"][sel].astype(np.float64))
            w = gw * norm * fr["w_pu_nom"][sel] * fr["w_mu_nom"][sel] * fr["w_btag_nom"][sel]
            if not np.all(np.isfinite(w)):
                sys.exit(f"non-finite weight in {nano_path} after selection")
        out = {b: nano[b][sel] for b in branches}
        out["weight"] = w
        if not is_data:
            out["norm"] = np.full(n_sel, norm, dtype=np.float64)
        for fcol, ocol in FRIEND_CORR_KIN.items():
            out[ocol] = fr[fcol][sel]
        for k, v in tags.items():
            out[k] = np.full(n_sel, v, dtype=np.int32)
        parts.append(out)
        n_sel_total += n_sel
        wsum += float(w.sum())

    meta = dict(era=era, sample=sample, list=os.path.basename(list_path),
                files=[os.path.basename(p) for p in nano_files], n_files=len(nano_files),
                n_rows_total=n_total, n_rows_selected=n_sel_total, n_rows_blinded=n_blinded,
                weight_sum=wsum,
                norm=norm, xsec_pb=xsec_pb, lumi_pb=lumi_pb,
                branches_sha1=None, root_written=False)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if n_sel_total:
        cols = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        meta["branches_sha1"] = hashlib.sha1("\n".join(sorted(cols)).encode()).hexdigest()
        with tempfile.TemporaryDirectory() as tmp:
            local = os.path.join(tmp, "chunk.root")
            with uproot.recreate(local) as f:
                f["Events"] = cols
            shutil.copyfile(local, out_path)
        with uproot.open(out_path) as f:
            t = f["Events"]
            back = t.arrays(["event", "weight"], library="np")
            if (t.num_entries != n_sel_total or not np.array_equal(back["event"], cols["event"])
                    or not np.array_equal(back["weight"], cols["weight"])):
                sys.exit(f"read-back of {out_path} does not match")
        meta["root_written"] = True
    elif os.path.exists(out_path):
        os.remove(out_path)
    with open(out_path[:-len(".root")] + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"{era}/{sample} {os.path.basename(list_path)}: files={len(nano_files)} "
          f"rows={n_total} sel={n_sel_total} blinded={n_blinded} weight_sum={wsum:.6g}"
          f" -> {out_path}")


# ---------------------------------------------------------------- merge
def merge(group):
    import uproot

    if group not in GROUPS:
        sys.exit(f"unknown group {group!r}, choose from {sorted(GROUPS)}")
    roots, piece_meta, shas, problems = [], [], set(), []
    for era, sample in GROUPS[group]:
        lists = piece_chunks(era, sample)
        if not lists:
            problems.append(f"{era}/{sample}: no chunk lists in {LIST_DIR}")
            continue
        listed = []
        agg = dict(era=era, sample=sample, n_chunks=len(lists), n_files=0, n_rows_total=0,
                   n_rows_selected=0, n_rows_blinded=0, weight_sum=0.0,
                   norm=None, xsec_pb=None, lumi_pb=None)
        ks = [int(re.search(r"_chunk(\d+)\.txt$", lp).group(1)) for lp in lists]
        if ks != list(range(len(lists))):
            problems.append(f"{era}/{sample}: chunk numbers {ks} are not 0..{len(lists) - 1}")
            continue
        for k, lp in zip(ks, lists):
            with open(lp) as f:
                names = [os.path.basename(x.strip()) for x in f if x.strip()]
            listed += names
            rp = skim_path(era, sample, k)
            jp = rp[:-len(".root")] + ".json"
            if not os.path.exists(jp):
                problems.append(f"missing {jp}")
                continue
            with open(jp) as f:
                m = json.load(f)
            if m["files"] != names:
                problems.append(f"{jp}: file list differs from {lp}")
                continue
            if m["root_written"]:
                if not os.path.exists(rp):
                    problems.append(f"missing {rp}")
                    continue
                roots.append(rp)
                shas.add(m["branches_sha1"])
            for key in ("n_files", "n_rows_total", "n_rows_selected", "n_rows_blinded",
                        "weight_sum"):
                agg[key] += m[key]
            for key in ("norm", "xsec_pb", "lumi_pb"):
                agg[key] = m[key]
        on_disk = [os.path.basename(p) for p in piece_nano_files(era, sample)]
        if sorted(listed) != sorted(on_disk) or len(set(listed)) != len(listed):
            problems.append(f"{era}/{sample}: chunk lists ({len(listed)}) do not cover the "
                            f"{len(on_disk)} nano files on disk exactly")
        piece_meta.append(agg)
    if len(shas) > 1:
        problems.append(f"chunks disagree on the branch list ({len(shas)} variants)")
    if problems:
        sys.exit("merge refused:\n  " + "\n  ".join(problems))

    target = f"{OUT_DIR}/{group}.root"
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        local = os.path.join(tmp, f"{group}.root")
        subprocess.run(["hadd", "-f", local] + roots, check=True)
        shutil.copyfile(local, target)
    n_expect = sum(p["n_rows_selected"] for p in piece_meta)
    with uproot.open(target) as f:
        n_got = f["Events"].num_entries
    if n_got != n_expect:
        sys.exit(f"{target}: {n_got} entries, chunk JSONs sum to {n_expect}")
    with open(f"{OUT_DIR}/{group}_meta.json", "w") as f:
        json.dump(piece_meta, f, indent=1)
    for p in piece_meta:
        xl = "" if p["xsec_pb"] is None else f" xsec*L={p['xsec_pb'] * p['lumi_pb']:.6g}"
        print(f"  {p['era']}/{p['sample']}: chunks={p['n_chunks']} files={p['n_files']} "
              f"rows={p['n_rows_total']} sel={p['n_rows_selected']} "
              f"blinded={p['n_rows_blinded']} weight_sum={p['weight_sum']:.6g}{xl}")
    print(f"wrote {target} ({n_got} entries from {len(roots)} chunks)")


def main():
    if len(sys.argv) in (4, 5) and sys.argv[1] == "worker" and sys.argv[4:] in ([], ["--correct"]):
        worker(sys.argv[2], sys.argv[3], correct=len(sys.argv) == 5)
    elif len(sys.argv) == 3 and sys.argv[1] == "merge":
        merge(sys.argv[2])
    else:
        sys.exit("usage: make_samples_nano.py worker <list.txt> <out.root> [--correct]\n"
                 "       make_samples_nano.py merge <" + "|".join(GROUPS) + ">")


if __name__ == "__main__":
    main()
