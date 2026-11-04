"""JME: jet veto map and the jet/MET columns. No jet energy correction (2026-10-21, user
decision): jets and MET are the nano values (Jet_pt with the production GT's JEC, no JER in
MC; slimmedMETsPuppi). The jet veto map is evaluated with the jerc-application-tutorial's
JvmApplication.apply_jvm_flag (via jme_corrections.jet_veto_flag) on nano Jet_pt.
Jet ID is the nano Jet_jetId bit (bit 2 = tight)."""
import awkward as ak

import jme_corrections

BRANCHES = sorted(set(["Jet_pt", "Jet_jetId", "met_pt", "met_phi"] + jme_corrections.JVM_BRANCHES))


def columns(a, jme_year, null):
    """{"jvm_veto_nom", "jet_ht_pt20_nom", "num_jets_pt30_nom", "jet_leadpt_nom",
    "met_pt_nom", "met_phi_nom"}."""
    pt = a["Jet_pt"]
    tight = (a["Jet_jetId"] & 2) == 2
    return {
        "jvm_veto_nom": ak.to_numpy(jme_corrections.jet_veto_flag(a, pt, year=jme_year)),
        "jet_ht_pt20_nom": ak.to_numpy(ak.sum(ak.where(tight & (pt > 20), pt, 0.0), axis=1)),
        "num_jets_pt30_nom": ak.to_numpy(ak.sum(tight & (pt > 30), axis=1)),
        "jet_leadpt_nom": ak.to_numpy(ak.fill_none(ak.max(pt[tight], axis=1), null)),
        "met_pt_nom": ak.to_numpy(a["met_pt"]).astype(float),
        "met_phi_nom": ak.to_numpy(a["met_phi"]).astype(float),
    }
