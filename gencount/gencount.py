#!/usr/bin/env python3
"""Generator-level composition of the H->Ds(*)W signal samples, for the normalisation.

Usage (CMSSW env with FWLite; the only part of this repository that needs CMSSW):
    gencount.py <DsW|DsstW> <filelist.txt> <JOBID> <GROUPSIZE> <out.json>

The signal cross sections are given for the final state
    Ds -> phi pi, phi -> K+K-,  W -> mu nu  or  W -> tau nu, tau -> mu nu nu
while the samples are generated with W -> e/mu/tau all open (tau decayed by Pythia at
its natural branching ratios) and the Ds(*) chain FORCED in EvtGen:
    DsW  : Ds -> phi pi (1.0), phi -> K+K- (1.0)
    DsstW: Ds* -> Ds gamma (0.942) / Ds pi0 (0.058), then the same Ds chain
The forcing is not always effective (evtgen_if picks one forced particle per event at
random; the old samples measured 94.4% / 97.2%), so the fraction of generated events
inside the cross-section's final state has to be counted on the full generated sample.
The nano cannot be used for this: it stores only events passing the reco selection.

Per event, from prunedGenParticles:
  chain   the hard-process Ds(*) (last copy) decays exactly as forced; extra photons from
          PHOTOS final-state radiation are allowed in every decay
  W       daughter lepton of the hard-process W: e / mu / tau; for tau, a final-state
          e or mu flagged isDirectHardProcessTauDecayProduct -> tau_e / tau_mu, else tau_had
Output: JSON with the counts, including the joint count chain & (W->mu or tau->mu).
"""
import json
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Events, Handle  # noqa: E402

RESID = {"DsW": 431, "DsstW": 433}


def last_copy(gp):
    """Follow same-pdgId daughters down to the last copy."""
    while True:
        nxt = None
        for j in range(gp.numberOfDaughters()):
            d = gp.daughter(j)
            if d.pdgId() == gp.pdgId():
                nxt = d
                break
        if nxt is None:
            return gp
        gp = nxt


def daughters_no_fsr(gp):
    """Daughter pdgIds (sorted), with FSR photons removed and the daughter objects by pdgId."""
    ids, objs = [], {}
    for j in range(gp.numberOfDaughters()):
        d = gp.daughter(j)
        if d.pdgId() == 22:
            continue
        ids.append(d.pdgId())
        objs.setdefault(d.pdgId(), d)
    return sorted(ids), objs


def ds_chain_ok(ds):
    """Ds -> phi pi(same charge as Ds), phi -> K+ K-."""
    ds = last_copy(ds)
    q = 1 if ds.pdgId() > 0 else -1
    ids, objs = daughters_no_fsr(ds)
    if ids != sorted([333, q * 211]):
        return False
    phi = last_copy(objs[333])
    ids, _ = daughters_no_fsr(phi)
    return ids == [-321, 321]


def dsst_chain_ok(dsst):
    """Ds* -> Ds gamma or Ds pi0, then the Ds chain. Returns (ok, mode) mode in gamma/pi0."""
    dsst = last_copy(dsst)
    q = 1 if dsst.pdgId() > 0 else -1
    ids = []
    ds = None
    has_pi0 = False
    n_gamma = 0
    for j in range(dsst.numberOfDaughters()):
        d = dsst.daughter(j)
        pid = d.pdgId()
        if pid == q * 431:
            ds = d
        elif pid == 111:
            has_pi0 = True
        elif pid == 22:
            n_gamma += 1
        else:
            ids.append(pid)
    if ds is None or ids:
        return False, None
    if has_pi0:
        mode = "pi0"
    elif n_gamma >= 1:
        mode = "gamma"
    else:
        return False, None
    return ds_chain_ok(ds), mode


def main():
    channel, filelist, jobid, groupsize, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    resid = RESID[channel]
    urls = [l.strip() for l in open(filelist) if l.strip()]
    files = urls[jobid * groupsize:(jobid + 1) * groupsize]
    keys = ["total", "no_signal_resonance", "no_hard_W", "chain",
            "W_e", "W_mu", "W_tau_e", "W_tau_mu", "W_tau_had", "W_other",
            "chain_and_mu_final", "mode_gamma", "mode_pi0", "chain_gamma", "chain_pi0"]
    c = dict.fromkeys(keys, 0)
    h = Handle("std::vector<reco::GenParticle>")
    for event in Events(files):
        c["total"] += 1
        event.getByLabel("prunedGenParticles", h)
        gps = h.product()
        res = w = None
        tau_e = tau_mu = False
        for gp in gps:
            pid = abs(gp.pdgId())
            flags = gp.statusFlags()
            if pid == resid and flags.isHardProcess() and res is None:
                res = gp
            elif pid == 24 and flags.isHardProcess() and w is None:
                w = gp
            if gp.status() == 1 and flags.isDirectHardProcessTauDecayProduct():
                if pid == 11:
                    tau_e = True
                elif pid == 13:
                    tau_mu = True
        if res is None:
            c["no_signal_resonance"] += 1
            continue
        if w is None:
            c["no_hard_W"] += 1
            continue

        if channel == "DsW":
            ok = ds_chain_ok(res)
        else:
            ok, mode = dsst_chain_ok(res)
            if mode:
                c[f"mode_{mode}"] += 1
                if ok:
                    c[f"chain_{mode}"] += 1
        c["chain"] += int(ok)

        wl = last_copy(w)
        lep = None
        for j in range(wl.numberOfDaughters()):
            a = abs(wl.daughter(j).pdgId())
            if a in (11, 13, 15):
                lep = a
        if lep == 11:
            cat = "W_e"
        elif lep == 13:
            cat = "W_mu"
        elif lep == 15:
            cat = "W_tau_mu" if tau_mu else ("W_tau_e" if tau_e else "W_tau_had")
        else:
            cat = "W_other"
        c[cat] += 1
        if ok and cat in ("W_mu", "W_tau_mu"):
            c["chain_and_mu_final"] += 1

    c.update(channel=channel, files=files)
    with open(out, "w") as fh:
        json.dump(c, fh, indent=1)
    print(json.dumps({k: v for k, v in c.items() if k != "files"}))


if __name__ == "__main__":
    main()
