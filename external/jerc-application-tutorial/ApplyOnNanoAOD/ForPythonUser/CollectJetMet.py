#!/usr/bin/env python3
# CollectJetMet.py
#
# Central helper for AK4/AK8 jet selection logic used by applyJecAndJvm.py.

from dataclasses import dataclass
import numpy as np
import awkward as ak
from typing import Optional
import PrintHelper

@dataclass
class JetSelections:
    """
    Convenience container for jet-path selections.

    - ak8_sel          : AK8 jets used in the AK8 jet path (mode-dependent).
    - basic_ak4        : Baseline AK4 selection used in AK4 jet path.
    - ak4_nonoverlap   : AK4 jets passing baseline and non-overlapping
                         with selected AK8 jets (used only in MODE_BOTH).
                         None in other modes.
    - ak4_mask_for_debug : Mask used to locate a representative debug event.
    """
    ak8_sel: ak.Array
    basic_ak4: ak.Array
    ak4_nonoverlap: ak.Array
    ak4_mask_for_debug: ak.Array


def build_jet_selections(arrs, mode_ak4: bool, mode_ak8: bool, mode_both: bool) -> JetSelections:
    """
    Construct the AK4/AK8 selection masks used in the jet paths.

    Parameters
    ----------
    arrs       : dict-like NanoAOD column container.
    mode_ak4   : True if AK4-only mode.
    mode_ak8   : True if AK8-only mode.
    mode_both  : True if combined (AK4+AK8) mode.

    Returns
    -------
    JetSelections
        Container with all masks needed by applyJecAndJvm.py.
    """
    # AK8 inputs
    ak8_pt  = arrs["FatJet_pt"]
    ak8_eta = arrs["FatJet_eta"]
    ak8_phi = arrs["FatJet_phi"]

    # AK4 inputs
    ak4_pt  = arrs["Jet_pt"]
    ak4_eta = arrs["Jet_eta"]
    ak4_phi = arrs["Jet_phi"]

    # ---------------------------
    # AK8 preselection
    # ---------------------------
    if mode_ak8 or mode_both:
        ak8_sel = (ak8_pt >= 100.0) & (abs(ak8_eta) <= 5.2)
    else:
        # Same structure, but no AK8 jet path
        ak8_sel = ak.zeros_like(ak8_pt, dtype=np.bool_)

    # ---------------------------
    # AK4 baseline preselection
    # ---------------------------
    basic_ak4 = (ak4_pt >= 15.0) & (abs(ak4_eta) <= 5.2)

    # ---------------------------
    # AK4 non-overlap w.r.t. selected AK8 (MODE_BOTH)
    # ---------------------------
    if mode_both:
        # Only compare to AK8 jets that passed the AK8 jet-path selection
        ak8_eta_sel = ak8_eta[ak8_sel]
        ak8_phi_sel = ak8_phi[ak8_sel]

        # Broadcast: (nEv, nAK4, 1) vs (nEv, 1, nAK8_sel)
        deta = ak4_eta[:, :, None] - ak8_eta_sel[:, None, :]
        dphi_ = ak4_phi[:, :, None] - ak8_phi_sel[:, None, :]
        dphi = np.arctan2(np.sin(dphi_), np.cos(dphi_)) 

        dr2 = deta * deta + dphi * dphi

        # Use a large initial value so events with 0 AK8 jets do not blow up
        min_dr = np.sqrt(ak.min(dr2, axis=2, initial=1.0e12))
        overlap = ak.fill_none((min_dr < 0.6), False)

        ak4_nonoverlap = basic_ak4 & (~overlap)
    else:
        ak4_nonoverlap = None

    # ---------------------------
    # Debug event mask
    # ---------------------------
    if mode_ak4:
        ak4_mask_for_debug = basic_ak4
    elif mode_both:
        ak4_mask_for_debug = ak4_nonoverlap
    else:
        ak4_mask_for_debug = ak.zeros_like(ak4_pt, dtype=np.bool_)

    return JetSelections(
        ak8_sel=ak8_sel,
        basic_ak4=basic_ak4,
        ak4_nonoverlap=ak4_nonoverlap,
        ak4_mask_for_debug=ak4_mask_for_debug,
    )


# ---------------------------
# Collect inputs for MET
# ---------------------------
@dataclass
class JetForMetColumns:
    """
    Columnar equivalent of the C++ JetForMet POD, but with jagged arrays:

      struct JetForMet {
        double phi, eta, area, rawPt;
        double muonSubtrFactor;
        double chEmEf, neEmEf;
      };
    """
    phi: ak.Array
    eta: ak.Array
    area: ak.Array
    rawPt: ak.Array
    muonSubtrFactor: ak.Array
    chEmEf: ak.Array
    neEmEf: ak.Array

    @property
    def rawPt_mu_sub(self):
        """Raw pT after muon subtraction (per jet)."""
        return self.rawPt * (1.0 - self.muonSubtrFactor)


@dataclass
class MetInputs:
    """
    struct MetInputs { double metPt; double metPhi; };
    """
    pt: ak.Array
    phi: ak.Array

    @property
    def px(self):
        return self.pt * np.cos(self.phi)

    @property
    def py(self):
        return self.pt * np.sin(self.phi)

def uses_puppi_met(year: str) -> bool:
    return year in {"2022Pre", "2022Post", "2023Pre", "2023Post", "2024", "2025"}


def build_met_inputs(arrs, year: str, print_evt: Optional[int]) -> MetInputs:
    """Extract raw MET (Puppi or PF) into a small POD-like object."""
    met_raw_pt  = arrs["RawPuppiMET_pt"]  if uses_puppi_met(year) else arrs["RawMET_pt"]
    met_raw_phi = arrs["RawPuppiMET_phi"] if uses_puppi_met(year) else arrs["RawMET_phi"]

    PrintHelper.debug_values("[MET] Nano MET_pt :", arrs["MET_pt"],
                             indent=4, event_index=print_evt)
    PrintHelper.debug_values("[MET] Raw MET pt :", met_raw_pt,
                             indent=4, event_index=print_evt)

    return MetInputs(pt=met_raw_pt, phi=met_raw_phi)


def build_corrt1_jets_for_met(arrs, print_evt: Optional[int]) -> JetForMetColumns:
    """Build JetForMetColumns for CorrT1 jets (counted only)."""
    ct1_pt_raw = arrs["CorrT1METJet_rawPt"]
    ct1_eta    = arrs["CorrT1METJet_eta"]
    ct1_phi    = arrs["CorrT1METJet_phi"]
    ct1_area   = arrs["CorrT1METJet_area"]
    ct1_muSub  = arrs["CorrT1METJet_muonSubtrFactor"]

    # No EM fractions for CorrT1 in NanoAOD → zeros
    ct1_chEm = ak.zeros_like(ct1_eta)
    ct1_neEm = ak.zeros_like(ct1_eta)

    jets = JetForMetColumns(
        phi=ct1_phi,
        eta=ct1_eta,
        area=ct1_area,
        rawPt=ct1_pt_raw,
        muonSubtrFactor=ct1_muSub,
        chEmEf=ct1_chEm,
        neEmEf=ct1_neEm,
    )

    PrintHelper.debug_values("[MET] CorrT1 raw pt:", jets.rawPt,
                             indent=6, event_index=print_evt)
    return jets


def build_ak4_jets_for_met(arrs, print_evt: Optional[int]) -> JetForMetColumns:
    """Build JetForMetColumns for AK4 jets."""
    ak4_pt   = arrs["Jet_pt"]
    ak4_eta  = arrs["Jet_eta"]
    ak4_phi  = arrs["Jet_phi"]
    ak4_area = arrs["Jet_area"]
    ak4_rawF = arrs["Jet_rawFactor"]
    ak4_muSub= arrs["Jet_muonSubtrFactor"]
    ak4_chEm = arrs["Jet_chEmEF"]
    ak4_neEm = arrs["Jet_neEmEF"]

    # Reconstruct raw pT from pT and rawFactor
    ak4_pt_raw_all = ak4_pt * (1.0 - ak4_rawF)

    jets = JetForMetColumns(
        phi=ak4_phi,
        eta=ak4_eta,
        area=ak4_area,
        rawPt=ak4_pt_raw_all,
        muonSubtrFactor=ak4_muSub,
        chEmEf=ak4_chEm,
        neEmEf=ak4_neEm,
    )

    PrintHelper.debug_values("[MET] AK4   raw pt:", jets.rawPt,
                             indent=6, event_index=print_evt)
    return jets


def combine_met_jets(ct1: JetForMetColumns,
                      ak4: JetForMetColumns) -> JetForMetColumns:
    """Concatenate CorrT1 and AK4 jet collections into a single columnar POD."""
    return JetForMetColumns(
        phi=ak.concatenate([ct1.phi,   ak4.phi],   axis=1),
        eta=ak.concatenate([ct1.eta,   ak4.eta],   axis=1),
        area=ak.concatenate([ct1.area,  ak4.area], axis=1),
        rawPt=ak.concatenate([ct1.rawPt, ak4.rawPt], axis=1),
        muonSubtrFactor=ak.concatenate([ct1.muonSubtrFactor, ak4.muonSubtrFactor], axis=1),
        chEmEf=ak.concatenate([ct1.chEmEf, ak4.chEmEf], axis=1),
        neEmEf=ak.concatenate([ct1.neEmEf, ak4.neEmEf], axis=1),
    )


