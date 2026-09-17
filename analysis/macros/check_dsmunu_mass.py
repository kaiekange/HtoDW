#!/usr/bin/env python3
"""Signal test of the full H mass m(Ds mu nu), neutrino pz from the W-mass constraint.

Usage: check_dsmunu_mass.py [era ...]     (default: 2022 2022EE)

Inputs per event, all at the NOMINAL corrections (friend trees from apply_corrections.py):
  Ds     = stored (Ds+mu) 4-vector minus the stored (uncorrected) muon -- exactly the Ds
           the producer added to the muon. For DsstW also Ds* = stored (Ds*+mu) - mu.
  mu     = stored muon direction, corrected pt (mu_pt_nom)
  nu     = px, py from the Type-1 corrected MET (met_pt_nom, met_phi_nom);
           pz from (mu + nu)^2 = mW^2:
             a    = mW^2/2 + px_mu px_nu + py_mu py_nu
             pz   = (a pz_mu +- E_mu sqrt(a^2 - pT_mu^2 pT_nu^2)) / pT_mu^2
           complex (negative discriminant) -> real part; two real -> smaller |pz|.
Selection: pass_mu_nom & pass_muveto & !jvm_veto_nom & n_bjet_nom == 0.

Gen cross-checks (W->mu events only, |Gen_W_lep_pdgId| == 13): the status-1 nu_mu
fromHardProcess gives the true pz, so we report
  - how often the smaller-|pz| choice is the solution closer to the truth,
  - m with the TRUE nu pz (reco MET px/py): isolates the pz ambiguity,
  - m with the full gen nu (reco Ds, reco mu): isolates the MET resolution.
Not used for any choice -- picking the solution closest to 125 would sculpt background.
"""
import glob
import os
import sys

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "DNN"))
import common  # noqa: E402  (fig_path)

NANO = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"
M_W = 80.377   # PDG 2022
M_MU = 0.1056583755
NULL = -777.0
FROM_HARD_PROCESS = 1 << 8

NANO_BR = ["run", "luminosityBlock", "event", "mu_px", "mu_py", "mu_pz", "mu_pt",
           "Dsmu_px", "Dsmu_py", "Dsmu_pz", "Dsmu_mass",
           "Dsstmu_px", "Dsstmu_py", "Dsstmu_pz", "Dsstmu_mass",
           "Ds_match", "Dsst_match", "Gen_W_lep_pdgId",
           "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags",
           "GenPart_pt", "GenPart_eta", "GenPart_phi"]
FRIEND_BR = ["run", "luminosityBlock", "event", "mu_pt_nom", "met_pt_nom", "met_phi_nom",
             "pass_mu_nom", "pass_muveto", "jvm_veto_nom", "n_bjet_nom"]


def four(px, py, pz, m):
    return np.stack([px, py, pz, np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2)])


def inv_mass(p):
    return np.sqrt(np.maximum(p[3] ** 2 - p[0] ** 2 - p[1] ** 2 - p[2] ** 2, 0.0))


def nu_pz(mu, nux, nuy):
    """Both W-constraint solutions (real part when complex) and the complex flag."""
    pt2 = mu[0] ** 2 + mu[1] ** 2
    a = 0.5 * M_W ** 2 + mu[0] * nux + mu[1] * nuy
    disc = a ** 2 - pt2 * (nux ** 2 + nuy ** 2)
    root = np.sqrt(np.maximum(disc, 0.0))
    return (a * mu[2] + mu[3] * root) / pt2, (a * mu[2] - mu[3] * root) / pt2, disc < 0


def load(era, sample):
    cols = {}
    for f in sorted(glob.glob(f"{NANO}/{era}/{sample}/output_*.root")):
        fr = f"{NANO}/{era}/Corrections/{sample}/{os.path.basename(f)}"
        a = uproot.open(f)["Events"].arrays(NANO_BR)
        b = uproot.open(fr)["Events"].arrays(FRIEND_BR, library="np")
        for k in ("run", "luminosityBlock", "event"):
            if not np.array_equal(ak.to_numpy(a[k]), b[k]):
                sys.exit(f"friend key mismatch in {fr}")
        # true nu_mu pz: status-1 |pdgId|==14 fromHardProcess (exactly one expected in W->mu)
        g = ((abs(a["GenPart_pdgId"]) == 14) & (a["GenPart_status"] == 1)
             & ((a["GenPart_statusFlags"] & FROM_HARD_PROCESS) != 0))
        ng = ak.to_numpy(ak.sum(g, axis=1))
        gpt = ak.to_numpy(ak.fill_none(ak.firsts(a["GenPart_pt"][g]), np.nan))
        geta = ak.to_numpy(ak.fill_none(ak.firsts(a["GenPart_eta"][g]), np.nan))
        gphi = ak.to_numpy(ak.fill_none(ak.firsts(a["GenPart_phi"][g]), np.nan))
        row = {k: ak.to_numpy(a[k]) for k in NANO_BR if not k.startswith("GenPart_")}
        row.update({k: b[k] for k in FRIEND_BR})
        row.update(gen_nnu=ng, gen_pt=gpt, gen_eta=geta, gen_phi=gphi)
        for k, v in row.items():
            cols.setdefault(k, []).append(v)
    return {k: np.concatenate(v) for k, v in cols.items()}


def q(x):
    x = x[np.isfinite(x)]
    p16, p50, p84 = np.percentile(x, [16, 50, 84])
    return (f"median {p50:6.1f}  16-84% [{p16:6.1f},{p84:6.1f}] (half-width {0.5*(p84-p16):5.1f})"
            f"  in [110,140]: {100*np.mean((x > 110) & (x < 140)):4.1f}%  n={len(x)}")


def main():
    eras = sys.argv[1:] or ["2022", "2022EE"]
    hists = {}
    for era in eras:
        for sample in ("Signal_DsW_H_phipimu", "Signal_DsstW_H_phipimu"):
            c = load(era, sample)
            sel = c["pass_mu_nom"] & c["pass_muveto"] & ~c["jvm_veto_nom"] & (c["n_bjet_nom"] == 0)
            f64 = {k: v.astype(float) for k, v in c.items() if v.dtype.kind == "f"}
            mu_raw = four(f64["mu_px"], f64["mu_py"], f64["mu_pz"], M_MU)
            k = f64["mu_pt_nom"] / f64["mu_pt"]
            mu = four(k * f64["mu_px"], k * f64["mu_py"], k * f64["mu_pz"], M_MU)
            dsmu = four(f64["Dsmu_px"], f64["Dsmu_py"], f64["Dsmu_pz"], f64["Dsmu_mass"])
            ds = dsmu - mu_raw
            nux = f64["met_pt_nom"] * np.cos(f64["met_phi_nom"])
            nuy = f64["met_pt_nom"] * np.sin(f64["met_phi_nom"])
            pz1, pz2, cplx = nu_pz(mu, nux, nuy)
            pz = np.where(np.abs(pz1) < np.abs(pz2), pz1, pz2)
            nu = four(nux, nuy, pz, 0.0)
            m_dsmunu = inv_mass(ds + mu + nu)

            wmu = np.abs(c["Gen_W_lep_pdgId"]) == 13
            dsm = c["Ds_match"].astype(bool)
            print(f"\n===== {era} {sample}: selected {sel.sum()}/{len(sel)}")
            print(f"  complex nu pz (real part used): {100*cplx[sel].mean():.1f}%")
            print(f"  m(Ds mu nu), all selected           : {q(m_dsmunu[sel])}")
            print(f"  m(Ds mu nu), W->mu & Ds matched     : {q(m_dsmunu[sel & wmu & dsm])}")
            print(f"  m(Ds mu nu), W->tau->mu & Ds matched: {q(m_dsmunu[sel & ~wmu & dsm])}")
            print(f"  Dsmu_mass (current feature), W->mu & Ds matched: {q(f64['Dsmu_mass'][sel & wmu & dsm])}")

            # gen cross-checks, W->mu, exactly one hard-process nu_mu
            gm = sel & wmu & dsm & (c["gen_nnu"] == 1)
            gpz = f64["gen_pt"] * np.sinh(f64["gen_eta"])
            closer_small = (np.abs(pz - gpz) <= np.abs(np.where(pz == pz1, pz2, pz1) - gpz))
            real = gm & ~cplx
            print(f"  [gen] W->mu with one hard nu_mu: {gm.sum()} (of {(sel & wmu & dsm).sum()});"
                  f" smaller-|pz| is the closer solution in {100*closer_small[real].mean():.1f}% of"
                  f" real-solution events")
            m_truepz = inv_mass(ds + mu + four(nux, nuy, gpz, 0.0))
            gnu = four(f64["gen_pt"] * np.cos(f64["gen_phi"]), f64["gen_pt"] * np.sin(f64["gen_phi"]), gpz, 0.0)
            m_gennu = inv_mass(ds + mu + gnu)
            print(f"  [gen] m with reco MET px,py + TRUE pz : {q(m_truepz[gm])}")
            print(f"  [gen] m with full gen nu               : {q(m_gennu[gm])}")

            hists[(era, sample, "all")] = m_dsmunu[sel]
            hists[(era, sample, "wmu")] = m_dsmunu[sel & wmu & dsm]
            hists[(era, sample, "gennu")] = m_gennu[gm]
            if sample == "Signal_DsstW_H_phipimu":
                ok = f64["Dsstmu_mass"] != NULL
                dsst = four(f64["Dsstmu_px"], f64["Dsstmu_py"], f64["Dsstmu_pz"], f64["Dsstmu_mass"]) - mu_raw
                m_dsst = inv_mass(dsst + mu + nu)
                dsstm = c["Dsst_match"].astype(bool)
                print(f"  m(Ds* mu nu), Ds* candidate present   : {q(m_dsst[sel & ok])}")
                print(f"  m(Ds* mu nu), W->mu & Ds* matched     : {q(m_dsst[sel & ok & wmu & dsstm])}")
                hists[(era, sample, "dsst_wmu")] = m_dsst[sel & ok & wmu & dsstm]

    bins = np.linspace(0, 300, 61)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, sample, title in ((axes[0], "Signal_DsW_H_phipimu", "DsW"),
                              (axes[1], "Signal_DsstW_H_phipimu", "DsstW")):
        keys = [("all", "selected, m(Ds mu nu)"), ("wmu", "W->mu & matched, m(Ds mu nu)"),
                ("gennu", "W->mu & matched, gen nu")]
        if title == "DsstW":
            keys.append(("dsst_wmu", "W->mu & Ds* matched, m(Ds* mu nu)"))
        for key, label in keys:
            x = np.concatenate([hists[(e, sample, key)] for e in eras])
            ax.hist(x, bins=bins, histtype="step", density=True, label=label)
        ax.axvline(125, color="gray", ls="--", lw=1)
        ax.set_xlabel("m [GeV]")
        ax.set_ylabel("normalised")
        ax.set_title(f"{title}, {'+'.join(eras)}")
        ax.legend(fontsize=8)
    out = common.fig_path("features", "dsmunu_mass_signal.png")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
