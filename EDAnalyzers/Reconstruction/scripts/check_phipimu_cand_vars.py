#!/usr/bin/env python3
"""Verify nano_phipimu_cff.py's variable maps against the tree classes they mirror.

The _phipimu_vars / _phipimu_gentag_vars dicts in python/nano_phipimu_cff.py must contain
exactly one same-named entry per branch that src/H_phipimu*_Tree.cc::CreateBranches() writes,
minus the five NanoAOD-native drops, plus nCandidates. This re-derives the expected sets
from the .cc/.h and diffs them against the cff. Run from the package root:

    python3 scripts/check_phipimu_cand_vars.py
"""
import re
import sys
import os

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP = {"event_id", "run_id", "lumi_id", "gen_weight", "jet_vetomap_pass"}


def tree_branches(cc_path, hdr_path):
    hdr = open(hdr_path).read()
    tmap = {}
    for m in re.finditer(r"(Float_t|Double_t|Int_t|UInt_t|ULong64_t|Bool_t)\s+([A-Za-z_]\w*)\s*;", hdr):
        tmap[m.group(2)] = m.group(1)
    body = open(cc_path).read().split("CreateBranches()")[1]
    order, seen = [], set()
    for b, mem in re.findall(r'tree->Branch\("([^"]+)",\s*&([A-Za-z_]\w*)\)', body):
        if b in seen:
            continue
        seen.add(b)
        if b in SKIP:
            continue
        order.append((b, "userFloat" if tmap[b] in ("Float_t", "Double_t") else "userInt"))
    return order


def cff_vars(dict_name):
    src = open(os.path.join(PKG, "python", "nano_phipimu_cff.py")).read()
    out = {}
    for name, acc in re.findall(dict_name + r"\[[\"']([\w]+)[\"']\]\s*=\s*Var\(\"(userFloat|userInt)\('", src):
        out[name] = acc
    return out


def check(label, cc, hdr, dict_name):
    expect = dict(tree_branches(os.path.join(PKG, cc), os.path.join(PKG, hdr)))
    expect["nCandidates"] = "userInt"
    got = cff_vars(dict_name)
    ok = True
    miss = set(expect) - set(got)
    extra = set(got) - set(expect)
    badacc = {k for k in set(expect) & set(got) if expect[k] != got[k]}
    if miss:
        print(f"[{label}] MISSING from cff: {sorted(miss)}"); ok = False
    if extra:
        print(f"[{label}] EXTRA in cff: {sorted(extra)}"); ok = False
    if badacc:
        print(f"[{label}] accessor mismatch: {sorted(badacc)}"); ok = False
    if ok:
        print(f"[{label}] OK ({len(expect)} variables, 1:1 with CreateBranches)")
    return ok


a = check("nominal", "src/H_phipimu_Tree.cc", "interface/H_phipimu_Tree.h", "_phipimu_vars")
b = check("gentag", "src/H_phipimu_GenTag_Tree.cc", "interface/H_phipimu_GenTag_Tree.h", "_phipimu_gentag_vars")
sys.exit(0 if (a and b) else 1)
