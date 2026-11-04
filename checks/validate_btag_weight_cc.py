#!/usr/bin/env python3
"""Check corrections/btag_weight.cc against the official BTV compute_btag_event_weight_multi_wp
(external/btv_examples/lib/btag_sf.py), event by event, on one MC nano file.

Usage (cmssw-el9 + CMSSW_15_0_5 cmsenv):
    python3 checks/validate_btag_weight_cc.py [nano MC file] [campaign]

An L/M/T efficiency map is built on the file's own b-veto jets with the official
compute_btagging_efficiency, then both implementations are evaluated for several WP lists and
variations. Exits 1 if any event differs by more than 1e-12.
btag_sf.py annotates `str | list` (Python >= 3.10); it is compiled here with postponed
annotations so that Python 3.9 can run it unchanged.
"""
import __future__
import os
import sys
import tempfile
import types

import awkward as ak
import correctionlib
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "corrections"))
import btv  # noqa: E402

REF = os.path.join(HERE, os.pardir, "external", "btv_examples", "lib", "btag_sf.py")
NANO = sys.argv[1] if len(sys.argv) > 1 else \
    "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/TTto2L2Nu_H_phipimu_GenTag/output_42.root"
CAMPAIGN = sys.argv[2] if len(sys.argv) > 2 else "Run3-22CDSep23-Summer22-NanoAODv12"

mod = types.ModuleType("btag_sf")
with open(REF) as f:
    exec(compile(f.read(), REF, "exec", flags=__future__.annotations.compiler_flag, dont_inherit=True),
         mod.__dict__)

a = uproot.open(NANO)["Events"].arrays(btv.BRANCHES + btv.MC_BRANCHES)
m = btv.scored_jets(a)
jets = ak.zip({"pt": ak.values_astype(a["Jet_pt"][m], np.float64),
               "eta": ak.values_astype(a["Jet_eta"][m], np.float64),
               "hadronFlavour": ak.values_astype(a["Jet_hadronFlavour"][m], np.int64),
               "btagPNetB": ak.values_astype(a["Jet_btagPNetB"][m], np.float64)})
sfs = correctionlib.CorrectionSet.from_file(btv.sf_path(CAMPAIGN))
wpv = {w: sfs["particleNet_wp_values"].evaluate(w) for w in "LMT"}
eff_cs = mod.compute_btagging_efficiency(jets, "btagPNetB", wpv)

with tempfile.TemporaryDirectory() as tmp:
    eff_json = os.path.join(tmp, "eff_LMT.json")
    with open(eff_json, "w") as f:
        f.write(eff_cs._data)
    eff = correctionlib.CorrectionSet.from_file(eff_json)["particleNet_efficiency"]

    R = btv._load_root()
    rv = R.VecOps.AsRVec

    def flat(x, t):
        return np.ascontiguousarray(ak.to_numpy(ak.flatten(x)), dtype=t)

    pt, eta = flat(jets.pt, np.float64), flat(jets.eta, np.float64)
    fl, di = flat(jets.hadronFlavour, np.int32), flat(jets.btagPNetB, np.float64)
    cnt = np.ascontiguousarray(ak.to_numpy(ak.num(jets.pt)), dtype=np.int32)
    worst = 0.0
    for wps in (["T"], ["L", "M", "T"], ["M", "L"]):
        v = R.std.vector["std::string"]()
        for w in wps:
            v.push_back(w)
        bw = R.BtagWeight(btv.sf_path(CAMPAIGN), eff_json, "particleNet", "particleNet_efficiency", v)
        for vh, vl in (("central", "central"), ("up", "up"), ("down", "down"), ("up", "down"),
                       ("up_correlated", "central")):
            ign = (vh, vl) != ("central", "central")
            py = ak.to_numpy(mod.compute_btag_event_weight_multi_wp(
                jets, "btagPNetB", sfs, eff, list(wps), vh, vl, ignore_warnings=ign)).astype(float)
            cc = np.array(bw.weights(rv(pt), rv(eta), rv(fl), rv(di), rv(cnt), vh, vl, ign))
            d = float(np.max(np.abs(py - cc)))
            worst = max(worst, d)
            print(f"{','.join(wps):6s} heavy={vh:14s} light={vl:8s} events {len(py)}  max|diff| {d:.2e}")
print(f"worst {worst:.2e}")
sys.exit(1 if worst > 1e-12 else 0)
