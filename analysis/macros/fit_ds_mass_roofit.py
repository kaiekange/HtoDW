"""Ds_mass sideband fit of MC scale factors with RooFit/HistFactory (run under cmsenv).

    cmsenv; python3 fit_ds_mass_roofit.py <fit dir>

<fit dir> holds inputs.npz + meta.json written by `ds_mass_fit_io.py inputs` (see its docstring
for the model). This script only needs ROOT + numpy, so it runs under cmsenv, which has no torch.

Likelihood: binned Poisson over the Ds-mass SIDEBAND bins only (mapped to a compact bin-index
observable, so no split fit ranges), one HistFactory channel per fitted region. The regions are
meta.json's fit_regions: the disjoint score slices (joint=slices) or the nested cuts (joint=nested;
these share events, so every event above the tightest cut enters the likelihood once per region
containing it -- errors are underestimated and the tight end is over-weighted). Each
floated factor is one NormFactor shared by name across channels, bounded to [0, 10]; fixed
components enter at their MC yield. MC statistics: Barlow-Beeston-lite StatError (one gamma per
bin per channel, Poisson constraint, threshold 0 -> every bin), so a few heavy DY / negative-weight
rows cannot make the fit overconfident. Errors: HESSE, plus MINOS on every factor.

Only the shared-factor joint fit is done (user, 2026-10-16: no per-region fits). Printed alongside:
the old linear chi2 (sigma^2 = data + sum MC sumw2) on the same regions, as a cross-check.
Writes <fit dir>/result.json, read by `ds_mass_fit_io.py draw`.
"""
import json
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
HF = ROOT.RooStats.HistFactory

FITDIR = sys.argv[1]
NO_STAT = "--no-stat" in sys.argv   # diagnostic: drop the Barlow-Beeston-lite gammas
A = np.load(f"{FITDIR}/inputs.npz")
with open(f"{FITDIR}/meta.json") as f:
    META = json.load(f)
VIS = A["vis"].astype(bool)
NB = int(VIS.sum())
COMPS = [(g, p) for g, p in META["comps"]]
FLOAT = META["float"]
REGIONS = META["fit_regions"]
LO, HI = 0.0, 10.0


def th1(name, w, w2=None):
    h = ROOT.TH1D(name, name, NB, 0, NB)
    h.SetDirectory(0)
    for i, (v, e2) in enumerate(zip(w[VIS], (w2 if w2 is not None else w)[VIS])):
        h.SetBinContent(i + 1, v)
        h.SetBinError(i + 1, np.sqrt(max(e2, 0.0)))
    return h


def build(regions):
    """HistFactory workspace, one channel per region; every floated factor is one parameter shared by all."""
    meas = HF.Measurement("meas", "meas")
    meas.SetExportOnly(True)
    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.01)
    meas.AddConstantParam("Lumi")
    hold = []
    names = set()
    for r in regions:
        ch = HF.Channel(r)
        hd = th1(f"{r}_data", A[f"{r}_data"]); hold.append(hd)
        ch.SetData(hd)
        ch.SetStatErrorConfig(0.0, "Poisson")
        for g, p in COMPS:
            h = th1(f"{r}_{g}", A[f"{r}_{g}_w"], A[f"{r}_{g}_w2"]); hold.append(h)
            s = HF.Sample(f"{r}_{g}")
            s.SetHisto(h)
            s.SetNormalizeByTheory(False)
            if not NO_STAT:
                s.ActivateStatError()
            if p is not None:
                s.AddNormFactor(p, 1.0, LO, HI)
                names.add(p)
            ch.AddSample(s)
        meas.AddChannel(ch)
    meas.SetPOI(sorted(names)[0])
    ws = HF.HistoToWorkspaceFactoryFast.MakeCombinedModel(meas)
    return ws, sorted(names), hold


def fit(regions, minos=False):
    ws, names, _hold = build(regions)
    mc = ws.obj("ModelConfig")
    pdf, data = mc.GetPdf(), ws.data("obsData")
    nll = pdf.createNLL(data, ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),
                        ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()), ROOT.RooFit.Offset(False))
    m = ROOT.RooMinimizer(nll)
    m.setPrintLevel(-1)
    m.setPrintEvalErrors(-1)   # MIGRAD probing regions with negative bin yields: transient, not a failure
    m.setStrategy(2)
    st = m.minimize("Minuit2", "Migrad")
    m.hesse()
    pars = [ws.var(n) for n in names]
    if minos:
        m.minos(ROOT.RooArgSet(*pars))
    res = m.save()
    out = dict(status=int(st), covQual=int(res.covQual()), nll=float(nll.getVal()), params={})
    for v in pars:
        q = dict(val=v.getVal(), err=v.getError(), at_bound=bool(v.getVal() < LO + 1e-3 or v.getVal() > HI - 1e-3))
        if minos and v.hasAsymError():
            q.update(minos_lo=v.getAsymErrorLo(), minos_hi=v.getAsymErrorHi())
        out["params"][v.GetName()] = q
    out["corr"] = {f"{a}|{b}": res.correlation(a, b) for i, a in enumerate(names) for b in names[i + 1:]}
    gam = [(x.GetName(), x.getVal()) for x in ws.allVars() if x.GetName().startswith("gamma_stat")]
    out["gamma_min"] = min(gam, key=lambda t: t[1]) if gam else None
    out["gamma_max"] = max(gam, key=lambda t: t[1]) if gam else None
    # goodness of fit, same definition as fit_sideband_scales_stage1.py and the draw step:
    # chi2 = sum (n - mu)^2 / (n + sum_c f_c^2 sumw2_c), mu = factor-scaled templates (gammas at 1)
    c2, nb = 0.0, 0
    for r in regions:
        f = {g: 1.0 if p is None else out["params"][p]["val"] for g, p in COMPS}
        n = A[f"{r}_data"][VIS]
        mu = sum(f[g] * A[f"{r}_{g}_w"][VIS] for g, _p in COMPS)
        var = n + sum(f[g] ** 2 * A[f"{r}_{g}_w2"][VIS] for g, _p in COMPS)
        c2 += float((((n - mu) ** 2) / np.maximum(var, 1.0)).sum())
        nb += NB
    out["chi2"], out["ndf"] = c2, nb - len(names)
    return out


def linear_chi2(regions):
    """fit_sideband_scales_stage1.py's linear chi2, same inputs, shared factors."""
    ys, sigs, cols = [], [], {p: [] for p in FLOAT}
    for r in regions:
        dh = A[f"{r}_data"]
        sigs.append(np.sqrt(np.maximum(dh + sum(A[f"{r}_{g}_w2"] for g, _p in COMPS), 1.0))[VIS])
        ys.append((dh - sum(A[f"{r}_{g}_w"] for g, p in COMPS if p is None))[VIS])
        for p in FLOAT:
            cols[p].append(sum(A[f"{r}_{g}_w"] for g, pp in COMPS if pp == p)[VIS])
    y, sig = np.concatenate(ys), np.concatenate(sigs)
    X = np.column_stack([np.concatenate(cols[p]) for p in FLOAT])
    Am = X / sig[:, None]
    sol = np.linalg.lstsq(Am, y / sig, rcond=None)[0]
    err = np.sqrt(np.diag(np.linalg.inv(Am.T @ Am)))
    return dict(zip(FLOAT, zip(sol, err))), float((((y - X @ sol) / sig) ** 2).sum()), len(y) - len(FLOAT)


def show(tag, F):
    s = "  ".join(f"{n}={q['val']:.3f}" + (f"(+{q['minos_hi']:.3f}/{q['minos_lo']:.3f})" if "minos_lo" in q
                                         else f"±{q['err']:.3f}") + ("[BOUND]" if q["at_bound"] else "")
                  for n, q in F["params"].items())
    print(f"{tag:34s} {s}  | chi2 {F['chi2']:.1f}/{F['ndf']}  status {F['status']} covQual {F['covQual']}")


print(f"fit dir {FITDIR}\nkey {META['key']}  phi_win {META['phi_win']}  joint {META['joint']} {META['effs']}"
      f"  float {FLOAT}")
R = {"meta": META}
R["shared"] = fit(REGIONS, minos=True)
show(f"SHARED ({META['joint']} {'/'.join(f'{100 * e:.0f}' for e in META['effs'])})", R["shared"])
print("  correlations: " + "  ".join(f"{k} {v:+.2f}" for k, v in R["shared"]["corr"].items()))
print(f"  gamma range: {R['shared']['gamma_min']} .. {R['shared']['gamma_max']}")
lin, c2, ndf = linear_chi2(REGIONS)
print("  linear chi2 cross-check: " + "  ".join(f"{p}={v:.3f}±{e:.3f}" for p, (v, e) in lin.items())
      + f"  chi2/ndf {c2:.1f}/{ndf}")
R["linear_chi2"] = dict(params={p: dict(val=v, err=e) for p, (v, e) in lin.items()}, chi2=c2, ndf=ndf)
with open(f"{FITDIR}/result.json", "w") as f:
    json.dump(R, f, indent=1)
print(f"Wrote {FITDIR}/result.json")
