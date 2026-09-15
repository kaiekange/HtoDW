import glob, os, re, uproot, numpy as np
N="/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"; L="/afs/cern.ch/work/k/kakang/H2DW/Corrections/condor/condor_logs/apply_corrections"
tot_bad=0
for era in ("2022","2022EE"):
  for s in ("Signal_DsW_H_phipimu","Signal_DsstW_H_phipimu"):
    nanos=sorted(glob.glob(f"{N}/{era}/{s}/output_*.root")); bad=[]; rows=0; selw=[]
    for f in nanos:
        stem=os.path.basename(f)[:-5]; fr=f"{N}/{era}/Corrections/{s}/{stem}.root"
        log=f"{L}/{era}_{s}_{stem}.log"
        lt=open(log).read() if os.path.exists(log) else ""
        if "Wrote " not in lt or "FAILED" in lt or "Traceback" in lt: bad.append((stem,"log")); continue
        try:
            a=uproot.open(f)["Events"].arrays(["run","luminosityBlock","event"],library="np")
            b=uproot.open(fr)["Events"].arrays(["run","luminosityBlock","event","pass_mu_nom","pass_muveto","jvm_veto_nom","n_bjet_nom","w_pu_nom","w_mu_nom","w_btag_nom","Dsmu_phi_nom","Dsstmu_Ht_nom"],library="np")
        except Exception as e: bad.append((stem,f"open {e}")); continue
        if not all(np.array_equal(a[k],b[k]) for k in a): bad.append((stem,"keys")); continue
        sel=b["pass_mu_nom"]&b["pass_muveto"]&~b["jvm_veto_nom"]&(b["n_bjet_nom"]==0)
        w=b["w_pu_nom"]*b["w_mu_nom"]*b["w_btag_nom"]
        if np.isnan(w[sel]).any(): bad.append((stem,"NaN in sel"))
        rows+=len(a["event"]); selw.append(w[sel])
    sw=np.concatenate(selw) if selw else np.array([])
    print(f"{era} {s}: nano {len(nanos)} ok {len(nanos)-len(bad)} rows {rows} sel {len(sw)} <w> {sw.mean():.4f} bad {bad[:5]}")
    tot_bad+=len(bad)
print("TOTAL BAD", tot_bad)
