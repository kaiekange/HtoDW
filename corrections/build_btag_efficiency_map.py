#!/usr/bin/env python3
"""Build the MC b-tagging efficiency map for ONE sample group, as a correctionlib
CorrectionSet on disk.

Usage: build_btag_efficiency_map.py <DsW|DsstW|WtoLNu2Jets>_<2022|2022EE>
Output: /eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/btag_efficiency/<group>.json
        (one directory for both eras; the era is in the file name)

WHY PER SAMPLE. The fixed-WP weight Pi (1-SF*eps)/(1-eps) needs the MC tagging
efficiency of the sample being corrected -- eps is a property of the simulation, the SF
carries the data/MC difference, and with SF=1 the weight is 1 for any eps. Using a
sample's own map on itself is therefore what the method requires, not a circularity.
The bias to avoid is the opposite one: an earlier version built a single W+jets map and
the pilot applied it to signal too. BTV's own guidance is to compute the map at least
separately for signal and background. DsW and DsstW get their own; the three W+jets
jet-bin samples are pooled into one map because they are one process sliced by jet
multiplicity. Signal has few genuine b jets, so its b bins are sparse; an empty bin
gets eps=0, i.e. a factor (1-SF*0)/(1-0) = 1 (no correction), as in the official code.

WHICH JETS. Exactly the population the b-veto examines, so eps is measured on the
same jets the weight is applied to:
  tight jet ID, |eta| < 2.5, dR(Ds, jet) > 0.4, and CORRECTED pt >= 20.
- corrected pt: JEC+JER nominal via jme_corrections.py (nano Jet_pt carries the GT's
  stale JEC and no JER smearing). An earlier version binned on nano Jet_pt.
- pt >= 20: BTV SFs are measured only for 20 <= pt < 1000, and outside that range
  particleNet_comb/light do NOT raise -- they silently return a flow constant (measured
  at |eta|=1, T WP: comb 0.837 at pt 19.9 vs 1.011 at pt 20). The veto is therefore
  defined on pt >= 20 jets, and so is this map.
- dR > 0.4: the veto excludes the jet the Ds sits in; an earlier version of this map
  included it.
Events: every event in the nano files, i.e. after the producer's event selection and
BEFORE the b-veto. Measuring after the veto would force eps -> 0 by construction.

compute-and-write logic is adapted from cms-analysis-corrections/BTV/examples
lib/btag_sf.py (same _ABSETA_EDGES/_PT_EDGES/_FLAVORS, same schemav2 structure, same
empty-bin convention). The one structural change is that counts are accumulated file
by file rather than from one in-memory jet array: W+jets is 1049 files and the JEC/JER
chain needs ~15 jet/genjet branches per file. Counts are unweighted, as in the official
code (a sign(genWeight)-weighted variant was tried earlier and gave negative bin counts).
"""
import glob
import os
import sys

import awkward as ak
import correctionlib
import correctionlib.schemav2 as cs
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jme_corrections as jme  # noqa: E402

NANO = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"
# Era -> JME year key of the jerc-application-tutorial config, and the ParticleNet tight WP
# (BTV btagging.json particleNet_wp_values, 2025-08-20; must match nano_phipimu_cff.PNET_WP).
ERAS = {
    "2022":   dict(jme_year="2022Pre",  wp_tight=0.6734),
    "2022EE": dict(jme_year="2022Post", wp_tight=0.6915),
}
# Group -> (era, input directories). One map per sample per era.
GROUPS = {}
for _era in ERAS:
    GROUPS[f"DsW_{_era}"] = (_era, [f"{NANO}/{_era}/Signal_DsW_H_phipimu"])
    GROUPS[f"DsstW_{_era}"] = (_era, [f"{NANO}/{_era}/Signal_DsstW_H_phipimu"])
    GROUPS[f"WtoLNu2Jets_{_era}"] = (_era, [f"{NANO}/{_era}/WtoLNu2Jets_{jb}J_H_phipimu_GenTag" for jb in (0, 1, 2)])
OUT_DIR = f"{NANO}/2022/btag_efficiency"

JET_PT_MIN = 20.0  # BTV SF validity floor, see docstring
_ABSETA_EDGES = [0, 1.305, 2.5]
_PT_EDGES = [20., 30., 50., 70., 100., 140., 200., 300., 600., 1000.]
_FLAVORS = [0, 4, 5]

BRANCHES = sorted(set(jme.JET_BRANCHES + [
    "Jet_jetId", "Jet_btagPNetB", "Jet_hadronFlavour", "Ds_eta", "Ds_phi"]))


def veto_jet_mask(a, pt_corr):
    """The jets the b-veto examines: tight ID, |eta|<2.5, dR(Ds,jet)>0.4, corrected pt>=20."""
    dphi = abs(a["Jet_phi"] - a["Ds_phi"])
    dphi = ak.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    dR = np.sqrt((a["Jet_eta"] - a["Ds_eta"]) ** 2 + dphi ** 2)
    return ((abs(a["Jet_eta"]) < 2.5) & (a["Jet_jetId"] >= 2) & (dR > 0.4)
            & (pt_corr >= JET_PT_MIN))


def accumulate(files, refs, wp_tight):
    """Per-flavour (total, pass) counts on the (abseta, pt) grid, summed over files."""
    shape = (len(_ABSETA_EDGES) - 1, len(_PT_EDGES) - 1)
    total = {f: np.zeros(shape) for f in _FLAVORS}
    passed = {f: np.zeros(shape) for f in _FLAVORS}
    n_events = 0
    for i, fn in enumerate(files):
        a = uproot.open(fn)["Events"].arrays(BRANCHES)
        if len(a) == 0:
            continue
        n_events += len(a)
        pt_corr = jme.corrected_jet_pt(a, refs, "nominal")
        m = veto_jet_mask(a, pt_corr)
        pt = ak.to_numpy(ak.flatten(pt_corr[m]))
        abseta = np.abs(ak.to_numpy(ak.flatten(a["Jet_eta"][m])))
        flav = ak.to_numpy(ak.flatten(a["Jet_hadronFlavour"][m]))
        tagged = ak.to_numpy(ak.flatten(a["Jet_btagPNetB"][m])) > wp_tight
        for f in _FLAVORS:
            fm = flav == f
            total[f] += np.histogram2d(abseta[fm], pt[fm], bins=[_ABSETA_EDGES, _PT_EDGES])[0]
            passed[f] += np.histogram2d(abseta[fm & tagged], pt[fm & tagged],
                                        bins=[_ABSETA_EDGES, _PT_EDGES])[0]
        if (i + 1) % 50 == 0 or i + 1 == len(files):
            print(f"  {i + 1}/{len(files)} files, {n_events} events", flush=True)
    return total, passed, n_events


def build_correction_set(total, passed, group):
    """Adapted from lib/btag_sf.py compute_btagging_efficiency (counts already filled)."""
    flav_content = []
    for flav in _FLAVORS:
        with np.errstate(divide="ignore", invalid="ignore"):
            # empty bins -> 0.0, exactly as the official code: factor (1-SF*0)/(1-0) = 1
            eff_vals = np.where(total[flav] > 0, passed[flav] / total[flav], 0.0)
        abseta_content = []
        for ieta in range(len(_ABSETA_EDGES) - 1):
            abseta_content.append(cs.Binning.parse_obj({
                "nodetype": "binning", "input": "pt", "edges": list(_PT_EDGES),
                "content": [float(eff_vals[ieta, ipt]) for ipt in range(len(_PT_EDGES) - 1)],
                "flow": "clamp",
            }))
        flav_content.append({"key": int(flav), "value": cs.Binning.parse_obj({
            "nodetype": "binning", "input": "abseta", "edges": list(_ABSETA_EDGES),
            "content": abseta_content, "flow": "clamp",
        })})
        n_pass, n_tot = int(passed[flav].sum()), int(total[flav].sum())
        n_empty = int((total[flav] == 0).sum())
        print(f"  wp=T flav={flav}: {n_pass}/{n_tot} jets pass"
              f" ({100 * n_pass / n_tot if n_tot else 0:.3f}%),"
              f" empty bins {n_empty}/{total[flav].size}")

    correction = cs.Correction.parse_obj({
        "version": 1,
        "name": "particleNet_efficiency",
        "description": (f"MC b-tagging efficiency for particleNet, {group}. Jets: tight ID,"
                        f" |eta|<2.5, dR(Ds,jet)>0.4, JEC+JER-corrected pt>=20."),
        "inputs": [
            {"name": "wp", "type": "string", "description": "Working point: T"},
            {"name": "flav", "type": "int", "description": "Jet hadron flavor: 0=udsg, 4=c, 5=b"},
            {"name": "abseta", "type": "real", "description": "Jet |eta|"},
            {"name": "pt", "type": "real", "description": "JEC+JER-corrected jet pt [GeV]"},
        ],
        "output": {"name": "efficiency", "type": "real", "description": "MC b-tagging efficiency"},
        "data": cs.Category.parse_obj({"nodetype": "category", "input": "wp", "content": [
            {"key": "T", "value": cs.Category.parse_obj(
                {"nodetype": "category", "input": "flav", "content": flav_content})},
        ]}),
    })
    return cs.CorrectionSet.parse_obj({
        "schema_version": 2,
        "description": f"MC b-tagging efficiency map, {group}",
        "corrections": [correction],
    })


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in GROUPS:
        sys.exit(f"usage: {sys.argv[0]} <{'|'.join(GROUPS)}>")
    group = sys.argv[1]
    era, dirs = GROUPS[group]
    files = sum([sorted(glob.glob(f"{d}/output_*.root")) for d in dirs], [])
    print(f"{group}: {len(files)} files")

    refs = jme.JmeRefs(year=ERAS[era]["jme_year"], is_data=False)
    total, passed, n_events = accumulate(files, refs, ERAS[era]["wp_tight"])
    print(f"events: {n_events}, veto-population jets:"
          f" {int(sum(t.sum() for t in total.values()))}")

    cset_model = build_correction_set(total, passed, group)
    os.makedirs(OUT_DIR, exist_ok=True)
    out_json = f"{OUT_DIR}/{group}.json"
    with open(out_json, "w") as f:
        f.write(cset_model.model_dump_json(exclude_unset=True))
    print(f"Wrote {out_json}")

    check = correctionlib.CorrectionSet.from_file(out_json)["particleNet_efficiency"]
    for flav, pt, eta in [(5, 40.0, 1.0), (4, 40.0, 1.0), (0, 40.0, 1.0), (5, 150.0, 2.0)]:
        print(f"  sanity eval: wp=T flav={flav} pt={pt} abseta={eta} -> "
              f"eff={check.evaluate('T', flav, eta, pt):.4f}")


if __name__ == "__main__":
    main()
