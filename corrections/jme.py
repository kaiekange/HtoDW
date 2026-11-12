"""JME: jet veto map and the JME-corrected jet/MET columns. The jet veto map is evaluated with the
jerc-application-tutorial's JvmApplication.apply_jvm_flag (via jme_corrections.jet_veto_flag) on the
JME jet pt. Jet ID is the nano Jet_jetId bit (bit 2 = tight).
2026-11-12 (user): the nano-jet/MET columns (*_nom: jvm_veto, jet_ht_pt20, num_jets_pt30, jet_leadpt,
met_pt, met_phi; written since 2026-10-21) are removed; only the JME ones below remain."""
import sys

import awkward as ak

import jme_corrections

# ---------------------------------------------------------------------------------------
# JME-corrected jets/MET (2026-11-04, user decision: these are the analysis default).
# jme_direct.Jme.variations(): JES (L1L2L3, data
# residual at the event's run) + MC JER smear; Type-1 PuppiMET rebuilt from RawPuppiMET with
# JER propagated, then MET xy. Systematics (MC): Reduced JES set (11 sources; "_year" sources
# and JER are uncorrelated across eras -> split by era downstream), JER (1 per era),
# unclustered energy (MET only). Jet ID: nano Jet_jetId bit.
#   jme_jvm_veto, jme_num_jets_pt30, jme_jet_ht_pt20, jme_jet_leadpt, jme_n_bjet,
#   jme_met_pt, jme_met_phi, jme_W_tmass (nominal muon), MC: jme_w_btag[_any]_{nom,up,down}
#   per variation v (MC): jme_{num_jets_pt30,jet_leadpt,n_bjet,w_btag,met_pt,met_phi}_<v>
#   (w_btag_<v>: central SF on the nominal-JME map at the varied pt, NaN where n_bjet_<v> > 0);
#   uncl{Up,Down}: met_pt/met_phi only. Variation columns are float32/int8.
# ---------------------------------------------------------------------------------------
import numpy as np  # noqa: E402

import btv  # noqa: E402
import jme_direct  # noqa: E402

FULL_BRANCHES = sorted(set(jme_direct.BRANCHES + jme_direct.UNCL_BRANCHES + jme_corrections.JVM_BRANCHES))
FULL_MC_BRANCHES = jme_direct.MC_BRANCHES


def _jet_cols(pt, tight, null):
    return (ak.to_numpy(ak.sum(tight & (pt > 30), axis=1)),
            ak.to_numpy(ak.sum(ak.where(tight & (pt > 20), pt, 0.0), axis=1)),
            ak.to_numpy(ak.fill_none(ak.max(pt[tight], axis=1), null)))


def full_columns(a, era, is_data, jme_year, tagger, wp_tight, campaign, btv_set, eff_json_jme,
                 mu_pt, mu_phi, null):
    J = jme_direct.Jme(era, is_data)
    var = J.variations(a)
    tight = (a["Jet_jetId"] & 2) == 2
    pt, met, mphi = var["nom"]

    # closures: nominal = the 10-22 validated chain (jet_pt; met(propagate_jer=True) + met_xy)
    if not ak.all(ak.flatten(pt) == ak.flatten(J.jet_pt(a))):
        sys.exit("JME closure FAILED: variations() nominal jet pt != Jme.jet_pt")
    m0, p0 = J.met_xy(*J.met(a, propagate_jer=True), ak.to_numpy(a["PV_npvsGood"]))
    if not (np.array_equal(m0, met) and np.array_equal(p0, mphi)):
        sys.exit("JME closure FAILED: variations() nominal MET != met(propagate_jer=True) + met_xy")
    d = np.abs(ak.to_numpy(a["PuppiMET_pt"]) - ak.to_numpy(a["met_pt"]))
    if (d > 1e-3 * np.maximum(ak.to_numpy(a["met_pt"]), 1.0)).any():
        sys.exit(f"nano met_pt != PuppiMET_pt on {int((d > 1e-3).sum())} rows")

    c = {}
    nj, ht, lead = _jet_cols(pt, tight, null)
    nb = btv.n_bjet(a, wp_tight, tagger, pt)
    c.update({"jme_jvm_veto": ak.to_numpy(jme_corrections.jet_veto_flag(a, pt, year=jme_year)),
              "jme_num_jets_pt30": nj, "jme_jet_ht_pt20": ht, "jme_jet_leadpt": lead,
              "jme_n_bjet": nb, "jme_met_pt": met, "jme_met_phi": mphi,
              "jme_W_tmass": np.sqrt(2.0 * mu_pt * met * (1.0 - np.cos(mu_phi - mphi)))})
    if is_data:
        return c
    c.update({f"jme_{k}": x for k, x in btv.weights(a, campaign, btv_set, eff_json_jme, nb, tagger, pt=pt).items()})
    n_stale = 0
    for v, (vpt, vmet, vphi) in var.items():
        if v == "nom":
            continue
        if vpt is not None:
            nj, _, lead = _jet_cols(vpt, tight, null)
            nbv = btv.n_bjet(a, wp_tight, tagger, vpt)
            w = btv.weights(a, campaign, btv_set, eff_json_jme, nbv, tagger, pt=vpt, sf_vars=("nom",), stale="count")
            n_stale += w["n_stale"]
            c[f"jme_num_jets_pt30_{v}"] = nj.astype(np.int8)
            c[f"jme_jet_leadpt_{v}"] = lead.astype(np.float32)
            c[f"jme_n_bjet_{v}"] = nbv.astype(np.int8)
            c[f"jme_w_btag_{v}"] = w["w_btag_nom"].astype(np.float32)
        c[f"jme_met_pt_{v}"] = np.asarray(vmet, dtype=np.float32)
        c[f"jme_met_phi_{v}"] = np.asarray(vphi, dtype=np.float32)
    print(f" JME: {len(var) - 1} variations; jets in eps=0/1 bins of the JME map summed over "
          f"variations: {n_stale} (their events get w_btag = 1 in that variation)")
    return c
