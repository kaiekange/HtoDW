"""2025 smoke-test checks on one NANO_DATA_H_phipimu.py output (CMSSW_15_0_10).

  python3 check_smoke_2025.py nano_2025C_os.root

1. jet ID: producer num_jets_pt20 vs nano sum(Jet_pt>20 && jetid.json AK4PUPPI_Tight), per event.
   (num_jets itself is not comparable: nano Jet has pt>15, the producer collection has no cut.)
2. b-score: producer jet_maxbscore (tight-ID, |eta|<2.5, any pt) vs max nano Jet_btagUParTAK4B
   over the same jets with pt>15 (nano maps <=0 to -1). Differences can only come from pt<15 jets.
3. b-score range / sentinel: distribution of Ds_nearestjet_bscore and jet_maxbscore.
"""
import sys
import numpy as np
import uproot
import awkward as ak
import correctionlib

JETID = "/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_10/src/EDAnalyzers/Reconstruction/data/jetid_2025_Summer24.json"

ev = uproot.open(sys.argv[1])["Events"]
n = ev.num_entries
print(f"Events: {n}")
a = ev.arrays(["Jet_pt", "Jet_eta", "Jet_chHEF", "Jet_neHEF", "Jet_chEmEF", "Jet_neEmEF", "Jet_muEF",
               "Jet_chMultiplicity", "Jet_neMultiplicity", "Jet_btagUParTAK4B",
               "num_jets_pt20", "jet_maxbscore", "Ds_nearestjet_bscore", "Ds_nearestjet_dR"])

cj = correctionlib.CorrectionSet.from_file(JETID)["AK4PUPPI_Tight"]
counts = ak.num(a.Jet_pt)
flat = lambda x: ak.to_numpy(ak.flatten(x))
chm = flat(a.Jet_chMultiplicity).astype(np.int64)
nem = flat(a.Jet_neMultiplicity).astype(np.int64)
jid = cj.evaluate(flat(a.Jet_eta), flat(a.Jet_chHEF), flat(a.Jet_neHEF), flat(a.Jet_chEmEF),
                  flat(a.Jet_neEmEF), flat(a.Jet_muEF), chm, nem, chm + nem) > 0.5
jid = ak.unflatten(jid, counts)
print(f"nano jets: {len(chm)}, pass tight ID: {int(ak.sum(jid))} ({ak.sum(jid)/len(chm):.3f})")

# 1. jet ID
nano_n20 = ak.to_numpy(ak.sum(jid & (a.Jet_pt > 20), axis=1))
prod_n20 = ak.to_numpy(a.num_jets_pt20)
d = prod_n20 - nano_n20
print(f"[jetID] num_jets_pt20 producer - nano: identical in {np.sum(d == 0)}/{n} events; "
      f"diff values {dict(zip(*np.unique(d, return_counts=True)))}")

# 2. b-score
sel = jid & (np.abs(a.Jet_eta) < 2.5)
nano_max = ak.to_numpy(ak.fill_none(ak.max(ak.where(sel, a.Jet_btagUParTAK4B, -999.0), axis=1), -999.0))
prod_max = ak.to_numpy(a.jet_maxbscore)
both = (nano_max > 0) & (prod_max > 0)
same = np.isclose(prod_max[both], nano_max[both], rtol=0, atol=2e-3)  # nano precision=12 bits
print(f"[bscore] jet_maxbscore vs nano max (both >0): equal in {same.sum()}/{both.sum()}; "
      f"producer > nano in {(prod_max[both][~same] > nano_max[both][~same]).sum()} (pt<15 jet allowed)")
bad = both.copy(); bad[both] = ~same
if (prod_max[bad] < nano_max[bad]).any():
    print("  !! producer max below nano max -> inputs differ:", list(zip(prod_max[bad], nano_max[bad]))[:10])

# 3. range / sentinel
for name in ("jet_maxbscore", "Ds_nearestjet_bscore"):
    v = ak.to_numpy(a[name])
    print(f"[range] {name}: min {v.min():.4g} max {v.max():.4g}; "
          f"<0 values {dict(zip(*np.unique(np.round(v[v < 0], 4), return_counts=True)))}; "
          f"in [0,1]: {np.sum((v >= 0) & (v <= 1))}/{n}")
raw = flat(a.Jet_btagUParTAK4B)
print(f"[range] nano Jet_btagUParTAK4B: min {raw.min():.4g} max {raw.max():.4g}; "
      f"==-1 for |eta|>=2.5: {np.sum(raw[np.abs(flat(a.Jet_eta)) >= 2.5] == -1)}/{np.sum(np.abs(flat(a.Jet_eta)) >= 2.5)}")
