#!/usr/bin/env python3
"""Per-sample MC normalisation: sum of generator weights before any selection.

Usage: normalization.py <era>/<sample dir> [<era>/<sample dir> ...]
       e.g. normalization.py 2022/WtoLNu2Jets_0J_H_phipimu_GenTag
Output: /eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC/<era>/Corrections/<sample dir>/normalization.json

    Lumi_reweight_i = genWeight_i * norm,   norm = sigma * L / sum(genWeight) = L / L_eff

This is the BACKGROUND (amcatnloFXFX) convention. SIGNAL is different as of 2026-09-17 --
see the SIGNAL section below: signal's weight formula drops genWeight entirely, so its
`norm` denominator is genEventCount, not genEventSumw. `norm_denominator` in the output
JSON records which one was actually used.

Run once per sample after its nano production has finished (the sum is over every file
present). It does not depend on the per-file friend trees, and they do not depend on it;
the two are combined only when the analysis samples are built.

DENOMINATOR (background). `Runs.genEventSumw`, filled by genWeightsTable in the nano job's
sumw_step path on every generated event, i.e. before the producer's selection. Each nano
job writes its own file (no hadd), so the sum is over the Runs entries of all files. A
missing or failed job removes its events from numerator and denominator alike: the
normalisation stays unbiased and only statistics is lost.

CHECK. The NTrueIntCounterTableProducer counter is filled on the same events with the
same weight, so for every file sum(nTrueIntSumw) + overflow + noInfo must equal
genEventSumw. A mismatch means the denominator does not cover the full generated sample.

CONSTANTS. Copied from CMSSW_15_0_5/src/analysis (this repository does not import it):
luminosity from eras.py ERAS[*]["lumi_pb"], W+jets cross sections from make_samples.py
XSEC_JETBIN (the 0.929803403 rescaling is already inside them -- do not apply it again).

SIGNAL. The cross sections (given by the analyst) are for the final state
    Ds -> phi pi, phi -> K+K-,  W -> mu nu  or  W -> tau nu, tau -> mu nu nu
  DsW  : 2 * 59.9 * 1.6e-5 * 0.0225 * (0.1063 + 0.1138*0.1737)                       = 5.437e-6 pb
  DsstW: 2 * 59.9 * 3.0e-5 * (0.936 + 0.0577*0.98823) * 0.0225 * (0.1063 + 0.1138*0.1737) = 10.13e-6 pb
  (2: charge conjugate; 59.9 pb: gg->H; then BR(H->Ds(*)W), BR(Ds*->Ds gamma) +
  BR(Ds*->Ds pi0) BR(pi0->gg), BR(Ds->phi pi) BR(phi->KK), BR(W->mu) + BR(W->tau) BR(tau->mu).)
The samples are generated with W -> e/mu/tau all open and the Ds(*) chain forced, so only
a fraction of the generated events lies in that final state:
    norm = sigma * L / (genEventCount * f_W_gen * f_Ds)
  genEventCount, NOT sum(genWeight): signal's per-event weight (make_samples_nano.py) does
           not carry a genWeight factor -- the ~49/1.6M events with genWeight=-23.1645
           instead of +23.1645 (found 2026-09-17) are a technical two-valued weight with
           no cross-section meaning for this LO-ish EvtGen+Pythia generation, unlike the
           amcatnloFXFX case above, so it is dropped rather than propagated. Using
           genEventCount here instead of genEventSumw is what keeps norm consistent with
           that: dropping genWeight from the numerator but not this denominator would
           silently rescale every signal yield by <genWeight> ~= 23.16.
  f_W_gen  fraction of generated W decays giving mu or tau->mu, from the branching ratios
           the GENERATOR used (Pythia 8.306 ParticleData.xml), not the PDG values in sigma:
           the numerator is the physical yield, the denominator counts generated events.
  f_Ds     fraction of generated events whose signal Ds(*) really decayed through the forced
           phi pi chain (EvtGen does not always apply the forcing), counted on the full
           MiniAOD by gencount/gencount.py. Its job i covers the same MiniAOD files as nano
           output_i, which is checked event count by event count.
"""
import glob
import json
import os
import sys

import numpy as np
import uproot

NANO = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"
LUMI_PB = {"2022": 7990.0, "2022EE": 26680.0}
XSEC_PB = {
    "WtoLNu2Jets_0J_H_phipimu_GenTag": 51846.0,
    "WtoLNu2Jets_1J_H_phipimu_GenTag": 8860.0,
    "WtoLNu2Jets_2J_H_phipimu_GenTag": 3284.0,
    "Signal_DsW_H_phipimu": 5.437e-6,
    "Signal_DsstW_H_phipimu": 10.13e-6,
}
SIGNALS = ("Signal_DsW_H_phipimu", "Signal_DsstW_H_phipimu")
# Pythia 8.306 (CMSSW_12_4_14_patch3, the generation release) ParticleData.xml
BR_GEN = dict(W_e=0.108166, W_mu=0.108166, W_tau=0.108087, tau_mu=0.1731072)
F_W_GEN = ((BR_GEN["W_mu"] + BR_GEN["W_tau"] * BR_GEN["tau_mu"])
           / (BR_GEN["W_e"] + BR_GEN["W_mu"] + BR_GEN["W_tau"]))
RUNS_BRANCHES = ["genEventCount", "genEventSumw", "genEventSumw2",
                 "nTrueIntSumw", "nTrueIntSumwOverflow", "nTrueIntSumwNoInfo"]


def sum_file(path):
    r = uproot.open(path)["Runs"].arrays(RUNS_BRANCHES, library="np")
    if len(r["genEventSumw"]) == 0:
        sys.exit(f"{path}: empty Runs tree")
    counter = (sum(float(np.sum(x)) for x in r["nTrueIntSumw"])
               + float(np.sum(r["nTrueIntSumwOverflow"])) + float(np.sum(r["nTrueIntSumwNoInfo"])))
    sumw = float(np.sum(r["genEventSumw"]))
    # A job killed mid-event can leave the two counters one event apart (seen on the
    # truncated 2022EE W+jets files: the difference was exactly one |genWeight|). That is
    # harmless for the denominator, so warn instead of stopping; a large mismatch still is
    # a real inconsistency and stops the run.
    rel = abs(counter - sumw) / max(abs(sumw), 1.0)
    if rel > 1e-3:
        sys.exit(f"{path}: nTrueInt counter {counter!r} != genEventSumw {sumw!r} (rel {rel:.3g})")
    if rel > 1e-6:
        print(f"  [warn] {os.path.basename(path)}: counter vs genEventSumw differ by"
              f" {counter - sumw:+.6g} (rel {rel:.3g}) -- truncated file, using genEventSumw")
    return dict(count=float(np.sum(r["genEventCount"])), sumw=sumw,
                sumw2=float(np.sum(r["genEventSumw2"])))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for target in sys.argv[1:]:
        era, sample = target.strip("/").split("/")
        if era not in LUMI_PB or sample not in XSEC_PB:
            sys.exit(f"unknown era/sample {target!r}")
        files = sorted(glob.glob(f"{NANO}/{era}/{sample}/output_*.root"))
        if not files:
            sys.exit(f"no output_*.root under {NANO}/{era}/{sample}")
        out_dir = f"{NANO}/{era}/Corrections/{sample}"
        tot = dict(count=0.0, sumw=0.0, sumw2=0.0)
        gen = {}
        for f in files:
            per = sum_file(f)
            for k, v in per.items():
                tot[k] += v
            if sample in SIGNALS:
                idx = os.path.basename(f)[len("output_"):-len(".root")]
                cpath = f"{out_dir}/gencount/count_{idx}.json"
                if not os.path.exists(cpath):
                    sys.exit(f"{target}: missing gencount {cpath}")
                cnt = json.load(open(cpath))
                if cnt["total"] != per["count"]:
                    sys.exit(f"{target}: {cpath} counted {cnt['total']} events but"
                             f" {os.path.basename(f)} has genEventCount {per['count']:.0f}")
                for k, v in cnt.items():
                    if isinstance(v, int):
                        gen[k] = gen.get(k, 0) + v
        if tot["sumw"] <= 0:
            sys.exit(f"{target}: sum(genWeight) = {tot['sumw']} -- refusing to normalise")

        # DENOMINATOR, signal vs background. Background (W+jets) is amcatnloFXFX: its
        # +-|genWeight| encodes a real NLO subtraction term, and Lumi_reweight_i =
        # genWeight_i * norm needs the SIGNED sum(genWeight) in the denominator for the
        # numerator/denominator cancellation to be self-consistent (the standard CMS
        # convention, unchanged here).
        #
        # Signal is EvtGen+Pythia (not NLO-matched): 2026-09-17, 7/869,385 post-selection
        # rows (49/1.6M pre-selection) were found to carry genWeight = -23.1645 instead of
        # +23.1645, a technical two-valued weight with no corresponding physical
        # cross-section meaning at LO here (unlike the amcatnloFXFX case) -- confirmed by
        # the user, who reads genWeight as not meaningful for this signal generation and
        # asked for it out of the signal weight entirely. So for SIGNALS, genWeight is
        # dropped from BOTH the per-event weight (make_samples_nano.py) and this
        # denominator: using genEventCount (plain N_gen) instead of genEventSumw keeps
        # numerator (weight = norm * w_pu * w_mu * w_btag, no genWeight factor) and
        # denominator consistent. Background is untouched.
        denom = tot["count"] if sample in SIGNALS else tot["sumw"]

        xsec, lumi = XSEC_PB[sample], LUMI_PB[era]
        out = dict(
            era=era, sample=sample, n_files=len(files),
            files=[os.path.basename(f) for f in files],
            genEventCount=tot["count"], genEventSumw=tot["sumw"], genEventSumw2=tot["sumw2"],
            n_eff=tot["sumw"] ** 2 / tot["sumw2"],
            lumi_pb=lumi, xsec_pb=xsec,
        )
        acc = 1.0
        if sample in SIGNALS:
            n = gen["total"]
            f_ds = gen["chain"] / n
            f_w_meas = (gen["W_mu"] + gen["W_tau_mu"]) / n
            acc = F_W_GEN * f_ds
            out.update(
                gencount=gen, br_gen=BR_GEN, f_W_gen=F_W_GEN,
                f_Ds=f_ds, f_Ds_err=(f_ds * (1 - f_ds) / n) ** 0.5,
                f_W_measured=f_w_meas, f_W_measured_err=(f_w_meas * (1 - f_w_meas) / n) ** 0.5,
                f_joint_measured=gen["chain_and_mu_final"] / n,
            )
        out.update(lumi_eff_pb=denom * acc / xsec, norm=xsec * lumi / (denom * acc),
                   norm_denominator=("genEventCount" if sample in SIGNALS else "genEventSumw"))
        os.makedirs(out_dir, exist_ok=True)
        with open(f"{out_dir}/normalization.json", "w") as fh:
            json.dump(out, fh, indent=1)
        extra = ""
        if sample in SIGNALS:
            extra = (f" f_Ds={out['f_Ds']:.5f}+-{out['f_Ds_err']:.5f}"
                     f" f_W_gen={F_W_GEN:.5f} f_W_meas={out['f_W_measured']:.5f}+-{out['f_W_measured_err']:.5f}")
        print(f"{target}: files={len(files)} genEventCount={tot['count']:.0f}"
              f" sumw={tot['sumw']:.6g} n_eff={out['n_eff']:.6g}{extra} norm={out['norm']:.6g}"
              f" -> {out_dir}/normalization.json")


if __name__ == "__main__":
    main()
