# H→Ds(*)W corrections

POG corrections (MUO, BTV, LUM, JME) for the custom H→Ds(*)W NanoAOD, written as one
friend tree per nano output file. No CMSSW code is used: the scripts only need Python with
correctionlib ≥ 2.6, uproot, awkward, numpy, scipy and ROOT (compiles muonscarekit `MuonScaRe.cc`).

The nano itself is produced separately (branch `CMSSW_13_0_13`, el8). This branch holds
only the correction step.

## Layout

| path | content |
|---|---|
| `corrections/apply_corrections.py` | per-file friend tree: reads the nano, calls the four POG modules below, runs the closures, writes the friend. Two jet/MET sets: the nano ones (`*_nom`, nano jets, kept unchanged) and, since 2026-11-04, the JME ones (`jme_*`: JES + MC JER, Type-1 PuppiMET with JER propagated + MET xy, plus 26 MC variations: Reduced JES 11×2, JER ×2, unclustered ×2). The analysis skim uses the `jme_*` set |
| `corrections/lum.py` | pileup weights |
| `corrections/muo.py` | muon scale/resolution (muonscarekit `MuonScaRe.cc` via ROOT) + ID×ISO×trigger SF, selected muon, prompt muon veto |
| `corrections/btv.py` | b-veto count + fixed-WP b-tag event weight via `btag_weight.cc` (C++ port of BTV examples `lib/btag_sf.py::compute_btag_event_weight_multi_wp`, commit in its header; compiled through ROOT) |
| `corrections/jme.py` | jet veto map + nano jet/MET columns (`columns`); JME columns and variations (`full_columns`, via `jme_direct.Jme.variations`, b-tag weight on the JME-pt efficiency map) |
| `corrections/build_btag_efficiency_map.py` | MC b-tag efficiency map per sample and era (input to the b-tag SF); `--jme`: binned on JME jet pt, written to `btag_efficiency_jme/` |
| `corrections/jme_corrections.py` | wrapper around the official JERC application tutorial; only its jet veto map is used by `apply_corrections.py` now |
| `corrections/jme_direct.py` | JEC/JER/Type-1 MET/MET xy straight from the JME correctionlib JSONs (no tutorial); `variations()` gives nominal + systematic jet pt and MET |
| `corrections/normalization.py` | per-sample normalisation (σ·L / Σw; signal: gen-level acceptance from `gencount/`) |
| `gencount/` | signal generator-level counts on MiniAOD (FWLite), input to `normalization.py` |
| `external/jerc-application-tutorial/` | unmodified subset of the official tutorial; upstream commit in `UPSTREAM` |
| `external/muonscarekit/` | unmodified `MuonScaRe.cc` from cms-muonPOG/muonscarekit; commit in `UPSTREAM` |
| `external/btv_examples/` | unmodified `lib/btag_sf.py` from BTV/examples, reference only (not imported) for `btag_weight.cc`; commit in `UPSTREAM` |
| `external/POG_EXAMPLES.txt` | URLs and commits of the six POG example repos compared against on 2026-11-04 |
| `checks/validate_btag_weight_cc.py` | `btag_weight.cc` vs the official Python function, event by event |
| `condor/` | HTCondor submission (one job per nano file / per efficiency-map group). `run_apply_corrections.sh` keeps a per-job log only if the job fails (~2e4 logs filled the AFS directory on 2026-11-04) |
| `checks/compare_{wjets,ss_os}_btag_maps.py` | is pooling the W+jets jet bins / sharing the OS map with SS justified (cell-by-cell efficiency comparison) |

## Running

Order: nano finished → efficiency maps → corrections.

```bash
cd condor
python3 condor_submit_btag_eff.py DsW_2022 DsstW_2022
python3 condor_submit_btag_eff.py --jme DsW_2022 DsstW_2022   # JME-pt maps, needed by the jme_* b-tag weight
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
