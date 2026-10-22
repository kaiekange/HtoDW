"""JME corrections straight from the JME POG correctionlib JSONs (no tutorial wrapper).

Inputs (pinned, not `latest`):
  jet_jerc / met_xyCorrections : JME/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05   (2022)
                                 JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05 (2022EE)
  JER smearing formula         : JME/JER-Smearing/2026-09-18 (v2: smeared pt floored at 0.01 GeV)

Chain (AK4 Puppi):
  jet : pt_raw = Jet_pt*(1-Jet_rawFactor) -> L1FastJet -> L2Relative -> L3Absolute
        -> (data) L2L3Residual(run) -> (MC) JER smear
  MET : RawPuppiMET - sum over (CorrT1METJet + Jet) of (pt_corr - pt_L1) (cos, sin),
        jets muon-subtracted, kept if pt_corr>15, |eta|<5.2, chEmEF+neEmEF<0.9;
        MC JER propagated or not (propagate_jer); optional MET xy correction on top.

Differences from jme_corrections.py (the tutorial wrapper) that are intended:
  * data residual is evaluated at the event's own run. The tutorial uses one representative
    run per year (2022Post: 359022 = Run E), but the 2022EE residual has three run bins
    (E/F/G) that differ by up to ~10%.
  * JER smearing uses the official v2 formula (floor at 0.01 GeV); the tutorial copy is v1.
"""

import awkward as ak
import correctionlib
import numpy as np

JME_BASE = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME"
JSON_DIR = {
    "2022": f"{JME_BASE}/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05",
    "2022EE": f"{JME_BASE}/Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05",
}
JER_SMEAR_JSON = f"{JME_BASE}/JER-Smearing/2026-09-18/jer_smear.json.gz"
TAG = {"2022": "Summer22_22Sep2023", "2022EE": "Summer22EE_22Sep2023"}
# tutorial's representative runs, only for closure against jme_corrections.py
TUTORIAL_RUN = {"2022": 355065.0, "2022EE": 359022.0}

JER_MATCH_DR = 0.2  # AK4: R/2


def _flat_eval(corr, *args):
    """Evaluate a correction on jagged arrays of identical structure; returns jagged."""
    counts = ak.num(args[0])
    flat = [ak.to_numpy(ak.flatten(x)) for x in args]
    return ak.unflatten(corr.evaluate(*flat), counts)


class Jme:
    def __init__(self, era, is_data, smear_json=JER_SMEAR_JSON):
        self.era, self.is_data = era, is_data
        cs = correctionlib.CorrectionSet.from_file(f"{JSON_DIR[era]}/jet_jerc.json.gz")
        lvl = "DATA" if is_data else "MC"
        p = f"{TAG[era]}_V4_{lvl}"
        self.l1 = cs[f"{p}_L1FastJet_AK4PFPuppi"]
        self.l2 = cs[f"{p}_L2Relative_AK4PFPuppi"]
        self.l3 = cs[f"{p}_L3Absolute_AK4PFPuppi"]
        self.res = cs[f"{p}_L2L3Residual_AK4PFPuppi"] if is_data else None
        if not is_data:
            self.reso = cs[f"{TAG[era]}_JRV2_MC_PtResolution_AK4PFPuppi"]
            self.sf = cs[f"{TAG[era]}_JRV2_MC_ScaleFactor_AK4PFPuppi"]
            self.smear = correctionlib.CorrectionSet.from_file(smear_json)["JERSmear"]
        self.metxy = correctionlib.CorrectionSet.from_file(
            f"{JSON_DIR[era]}/met_xyCorrections_2022_{era}.json.gz")["met_xy_corrections"]

    # ---------------- JES ----------------
    def jes(self, pt_raw, area, eta, rho, run, residual_run="event"):
        """Return (pt_after_L1, pt_after_full_JES). rho, run are per event."""
        rho_j = ak.broadcast_arrays(rho, eta)[0]
        pt_l1 = pt_raw * _flat_eval(self.l1, area, eta, pt_raw, rho_j)
        pt = pt_l1 * _flat_eval(self.l2, eta, pt_l1)
        pt = pt * _flat_eval(self.l3, eta, pt)
        if self.is_data:
            if residual_run == "event":
                run_j = ak.values_astype(ak.broadcast_arrays(run, eta)[0], np.float64)
            elif residual_run == "tutorial":
                run_j = ak.ones_like(eta) * TUTORIAL_RUN[self.era]
            else:
                raise ValueError(residual_run)
            pt = pt * _flat_eval(self.res, run_j, eta, pt)
        return pt_l1, pt

    # ---------------- JER ----------------
    def jer(self, pt, eta, phi, rho, event, gen_idx, gen_pt, gen_eta, gen_phi):
        """MC only: hybrid JER smearing (gen-matched scaling, else stochastic)."""
        rho_j = ak.broadcast_arrays(rho, eta)[0]
        reso = _flat_eval(self.reso, eta, pt, rho_j)
        sf = _flat_eval(self.sf, eta, pt)
        valid = (gen_idx >= 0) & (gen_idx < ak.num(gen_pt))
        safe = ak.where(valid, gen_idx, 0)
        # events without GenJets: pad one dummy so indexing is defined (never used: valid=False)
        def pad(x, v):
            return ak.fill_none(ak.pad_none(x, 1, axis=1), v)
        gpt_p, geta_p, gphi_p = pad(gen_pt, 0.0), pad(gen_eta, 999.0), pad(gen_phi, 0.0)
        gpt, geta, gphi = gpt_p[safe], geta_p[safe], gphi_p[safe]
        dphi = np.arctan2(np.sin(phi - gphi), np.cos(phi - gphi))
        dr = np.sqrt((eta - geta) ** 2 + dphi ** 2)
        match = valid & (dr < JER_MATCH_DR) & (abs(pt - gpt) < 3.0 * reso * pt)
        genpt = ak.where(match, gpt, -1.0)
        eid = ak.Array((ak.to_numpy(ak.values_astype(event, np.uint64)) & np.uint64(0x7FFFFFFF))
                       .astype(np.int64))
        eid_j = ak.broadcast_arrays(eid, eta)[0]
        s = _flat_eval(self.smear, pt, eta, genpt, rho_j, eid_j, reso, sf)
        return pt * s

    # ---------------- jets ----------------
    def jet_pt(self, a, residual_run="event"):
        pt_raw = a["Jet_pt"] * (1.0 - a["Jet_rawFactor"])
        _, pt = self.jes(pt_raw, a["Jet_area"], a["Jet_eta"], a["Rho_fixedGridRhoFastjetAll"],
                         a["run"], residual_run)
        if self.is_data:
            return pt
        return self.jer(pt, a["Jet_eta"], a["Jet_phi"], a["Rho_fixedGridRhoFastjetAll"], a["event"],
                        a["Jet_genJetIdx"], a["GenJet_pt"], a["GenJet_eta"], a["GenJet_phi"])

    # ---------------- Type-1 MET ----------------
    def met(self, a, propagate_jer=True, residual_run="event"):
        """Type-1 PuppiMET (pt, phi) rebuilt from RawPuppiMET."""
        ct1_raw = a["CorrT1METJet_rawPt"]
        ak4_raw = a["Jet_pt"] * (1.0 - a["Jet_rawFactor"])
        raw = ak.concatenate([ct1_raw, ak4_raw], axis=1)
        musub = ak.concatenate([a["CorrT1METJet_muonSubtrFactor"], a["Jet_muonSubtrFactor"]], axis=1)
        eta = ak.concatenate([a["CorrT1METJet_eta"], a["Jet_eta"]], axis=1)
        phi = ak.concatenate([a["CorrT1METJet_phi"], a["Jet_phi"]], axis=1)
        area = ak.concatenate([a["CorrT1METJet_area"], a["Jet_area"]], axis=1)
        em = ak.concatenate([ak.zeros_like(a["CorrT1METJet_eta"]), a["Jet_chEmEF"] + a["Jet_neEmEF"]],
                            axis=1)
        pt_mu = raw * (1.0 - musub)
        rho = a["Rho_fixedGridRhoFastjetAll"]
        pt_l1, pt = self.jes(pt_mu, area, eta, rho, a["run"], residual_run)
        if not self.is_data and propagate_jer:
            # CorrT1 jets: closest GenJet within dR<0.2; AK4 jets: nano genJetIdx
            ge, gp = a["GenJet_eta"], a["GenJet_phi"]
            ce, cp = a["CorrT1METJet_eta"], a["CorrT1METJet_phi"]
            dph = np.arctan2(np.sin(cp[:, :, None] - gp[:, None, :]), np.cos(cp[:, :, None] - gp[:, None, :]))
            dr2 = (ce[:, :, None] - ge[:, None, :]) ** 2 + dph ** 2
            amin = ak.fill_none(ak.argmin(dr2, axis=2), -1)
            mind = ak.fill_none(ak.min(dr2, axis=2), 1e9)
            ct1_idx = ak.values_astype(ak.where(mind < JER_MATCH_DR ** 2, amin, -1), np.int64)
            gidx = ak.concatenate([ct1_idx, ak.values_astype(a["Jet_genJetIdx"], np.int64)], axis=1)
            pt = self.jer(pt, eta, phi, rho, a["event"], gidx, a["GenJet_pt"], a["GenJet_eta"],
                          a["GenJet_phi"])
        keep = (pt > 15.0) & (abs(eta) < 5.2) & (em < 0.9)
        dpt = ak.where(keep, pt - pt_l1, 0.0)
        px = a["RawPuppiMET_pt"] * np.cos(a["RawPuppiMET_phi"]) - ak.sum(dpt * np.cos(phi), axis=1)
        py = a["RawPuppiMET_pt"] * np.sin(a["RawPuppiMET_phi"]) - ak.sum(dpt * np.sin(phi), axis=1)
        return ak.to_numpy(np.hypot(px, py)), ak.to_numpy(np.arctan2(py, px))

    def met_xy(self, met_pt, met_phi, npv_good):
        """MET xy correction (PuppiMET, nominal). Returns (pt, phi)."""
        dtmc = "DATA" if self.is_data else "MC"
        npv = np.asarray(npv_good, dtype=np.float64)
        pt = self.metxy.evaluate("pt", "PuppiMET", self.era, dtmc, "nom", met_pt, met_phi, npv)
        phi = self.metxy.evaluate("phi", "PuppiMET", self.era, dtmc, "nom", met_pt, met_phi, npv)
        return pt, phi


def jet_id_v12_fix(a):
    """(tight, tightLepVeto) per jet: JME JetID13p6TeV bugfix recipe for NanoAODv12."""
    aeta = abs(a["Jet_eta"])
    tight_bit = (a["Jet_jetId"] & 2) > 0
    tight = ak.where(aeta <= 2.7, tight_bit,
                     ak.where(aeta <= 3.0, tight_bit & (a["Jet_neHEF"] < 0.99),
                              tight_bit & (a["Jet_neEmEF"] < 0.4)))
    lepveto = ak.where(aeta <= 2.7, tight & (a["Jet_muEF"] < 0.8) & (a["Jet_chEmEF"] < 0.8), tight)
    return tight, lepveto


BRANCHES = [
    "run", "event", "Rho_fixedGridRhoFastjetAll", "PV_npvsGood",
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Jet_muonSubtrFactor",
    "Jet_chEmEF", "Jet_neEmEF", "Jet_neHEF", "Jet_muEF", "Jet_jetId",
    "CorrT1METJet_rawPt", "CorrT1METJet_eta", "CorrT1METJet_phi", "CorrT1METJet_area",
    "CorrT1METJet_muonSubtrFactor", "RawPuppiMET_pt", "RawPuppiMET_phi",
]
MC_BRANCHES = ["Jet_genJetIdx", "GenJet_pt", "GenJet_eta", "GenJet_phi"]
