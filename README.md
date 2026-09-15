# H→Ds(*)W corrections

POG corrections (JME, MUO, BTV, LUM) for the custom H→Ds(*)W NanoAOD, written as one
friend tree per nano output file. No CMSSW code is used: the scripts only need Python with
correctionlib ≥ 2.6, uproot, awkward, numpy, scipy and ROOT (for `TRandom3` in MuonScaRe).

The nano itself is produced separately (branch `CMSSW_13_0_13`, el8). This branch holds
only the correction step.

## Layout

| path | content |
|---|---|
| `corrections/apply_corrections.py` | per-file friend tree: muon scale/resolution + ID×ISO×trigger SF, prompt muon veto, PU weights, JEC+JER (Total up/down), Type-1 MET, jet veto map, b-veto + fixed-WP b-tag SF |
| `corrections/build_btag_efficiency_map.py` | MC b-tag efficiency map per sample and era (input to the b-tag SF) |
| `corrections/jme_corrections.py` | wrapper around the official JERC application tutorial |
| `corrections/muo_scalesmearing.py` | port of the official MuonScaRe recipe |
| `external/jerc-application-tutorial/` | unmodified subset of the official tutorial; upstream commit in `UPSTREAM` |
| `condor/` | HTCondor submission (one job per nano file / per efficiency-map group) |
| `checks/validate_friends.py` | per-file acceptance of a finished friend-tree batch |

## Running

Order: nano finished → efficiency maps → corrections.

```bash
cd condor
python3 condor_submit_btag_eff.py DsW_2022 DsstW_2022
python3 condor_submit_apply_corrections.py 2022/Signal_DsW_H_phipimu 2022/Signal_DsstW_H_phipimu
```

A single file interactively (writes to the default output unless a second path is given):

```bash
python3 corrections/apply_corrections.py /eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/2022/Signal_DsW_H_phipimu/output_0.root [out.root]
```

Input  `/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/<era>/<sample>/output_N.root`
Output `/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/<era>/Corrections/<sample>/output_N.root`
(tree `Events`, row-aligned with the nano file, carries run/luminosityBlock/event).

Environment: the condor scripts currently run inside `cmssw-el9` with the CMSSW_15_0_5
Python stack. The planned switch is a plain LCG view on el9 (`LCG_109`, correctionlib
2.6.4) with `PYTHONNOUSERSITE=1`, after a column-by-column comparison.

## Contract with the nano producer

These constants are duplicated here and must match `CMSSW_13_0_13`
(`HDsPhiPiMuProducer.cc`, `nano_phipimu_cff.py`):

- muon pre-cut `pt >= 24` (`MU_PRECUT`) — checked: the selected-muon re-identification
  must close on every row;
- ParticleNet tight WP per era (`build_btag_efficiency_map.ERAS`) — checked against the
  BTV JSON, not against the producer config;
- branch names read in `BRANCHES`.

## Reproducibility note

Running the same code on an lxplus node and on a condor worker gives bit-identical output
except for the last floating-point bit of eta/phi columns (~4e-16) and `W_tmass`
(~5e-11), from CPU-dependent `arctan2`/`arcsinh`.
