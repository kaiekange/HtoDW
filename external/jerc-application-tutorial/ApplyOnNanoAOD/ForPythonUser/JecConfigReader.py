import json
from functools import lru_cache
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import correctionlib._core as correction


# ============================================================
# Struct-like containers 
# ============================================================

@dataclass(frozen=True)
class ConfigPaths:
    ak4: str = "../../JecConfigAK4.json"
    ak8: str = "../../JecConfigAK8.json"


CONFIG_PATHS = ConfigPaths()


@dataclass
class Tags:
    """
    Container for the nominal tag names and the base JERC JSON path.
    """
    jercJsonPath: str
    tagNameL1FastJet: str
    tagNameL2Relative: str
    # Optional because:
    #   - Data has L2L3Residual
    #   - MC does not use residual corrections
    tagNameL2L3Residual: Optional[str]
    # Optional because JER is only relevant for MC
    tagNamePtResolution: Optional[str]
    tagNameJerScaleFactor: Optional[str]
    tagNameJerSFUncertainty: Optional[str]


@dataclass
class JesDataRefs:
    cs: correction.CorrectionSet
    l1FastJet: correction.Correction
    l2Relative: correction.Correction
    l2L3Residual: correction.Correction

@dataclass
class JesMcRefs:
    cs: correction.CorrectionSet
    l1FastJet: correction.Correction
    l2Relative: correction.Correction

@dataclass
class JerMcRefs:
    cs: correction.CorrectionSet
    ptResolution: correction.Correction
    scaleFactor: correction.Correction
    sfUncertainty: Optional[correction.Correction]

@dataclass
class JerBin:
    label: str
    etaMin: float
    etaMax: float
    ptMin: float
    ptMax: float

@dataclass
class JesUncSetRefs:
    """
    Keys are CMS names (e.g. 'CMS_scale_j_AbsoluteMPFBias').
    """
    full: Dict[str, correction.Correction]
    reduced: Dict[str, correction.Correction]
    total: Dict[str, correction.Correction]

@dataclass
class JerUncSets:
    """
    Keys are the CMS labels for the bins
    (e.g. 'CMS_res_j_eta0to1p93_pt0toInf_2016').
    """
    full: Dict[str, JerBin]
    total: Dict[str, JerBin]


# ============================================================
# JER smearing set (external JSON, year-agnostic node)
# ============================================================

@lru_cache(maxsize=1)
def get_jer_smear_cset() -> correction.CorrectionSet:
    return correction.CorrectionSet.from_file("../../jer_smear.json.gz")

def get_jer_smear() -> correction.Correction:
    cs = get_jer_smear_cset()
    try:
        return cs["JERSmear"]
    except Exception as e:
        raise RuntimeError(
            "Could not get 'JERSmear' from '../../jer_smear.json.gz': "
            f"{e}"
        )

# ============================================================
# Small JSON loader with caching
# ============================================================

@lru_cache(maxsize=None)
def load_json_cached(path: str) -> dict:
    """Read and cache small JSON configuration files by absolute path."""
    with open(path, "r") as f:
        return json.load(f)


# ============================================================
# Global cache: config JSONs and CorrectionRefs
# ============================================================

class CorrCache:
    """
    Caches:
      - Parsed JEC config JSONs per 'kind' (e.g. 'AK4', 'AK8')
      - CorrectionRefs objects per (kind, year, isData, era)
    """

    _cfg: Dict[str, dict] = {}  # kind -> parsed config JSON
    _refs: Dict[Tuple[str, str, bool, str], "CorrectionRefs"] = {}

    @classmethod
    def get_cfg(cls, kind: str) -> dict:
        """
        Return the parsed config JSON for a given 'kind'
        (e.g. 'AK4', 'AK8', ...).
        """
        if kind not in cls._cfg:
            kind_upper = kind.upper()
            if kind_upper == "AK4":
                path = CONFIG_PATHS.ak4
            elif kind_upper == "AK8":
                path = CONFIG_PATHS.ak8
            else:
                # fallback / extension point for future kinds
                path = f"../../JecConfig{kind}.json"

            cls._cfg[kind] = load_json_cached(path)
        return cls._cfg[kind]

    @classmethod
    def get_refs(
        cls,
        kind: str,
        year: str,
        isData: bool,
        era: Optional[str],
    ) -> "CorrectionRefs":
        """
        Return a CorrectionRefs object for a given (kind, year, isData, era).

        * kind: 'AK4', 'AK8', etc.  (maps to JecConfigAK4.json, JecConfigAK8.json)
        * year: e.g. '2016Pre'
        * isData: True for data, False for MC
        * era: e.g. 'Era2016PreBCD', only used for data
        """
        key = (kind, year, bool(isData), era or "")
        if key not in cls._refs:
            cfg = cls.get_cfg(kind)

            # Fallback for AK8 tags (reuse AK4 tags if AK8 not present for this year)
            #if kind == "AK8" and year not in cfg:
            #    cfg = cls.get_cfg("AK4")

            tags = get_tag_names(cfg, year, isData, era)
            cls._refs[key] = CorrectionRefs(tags, isData=isData)
        return cls._refs[key]


# ============================================================
# Tag container for nominal corrections
# ============================================================

def get_tag_names(base: dict, year: str, isData: bool, era: Optional[str]) -> Tags:
    """
    Expected layout (per year):
      {
        "2016Pre": {
          "jercJsonPath": "...jet_jerc.json.gz",
          "ApplyOnData": {
            "JesNominal": {
              "Era2016PreBCD": {
                "tagNameL1FastJet": "...",
                "tagNameL2Relative": "...",
                "tagNameL2L3Residual": "..."
              },
              "Era2016PreEF": { ... }
            }
          },
          "ApplyOnMC": {
            "JesNominal": {
              "tagNameL1FastJet": "...",
              "tagNameL2Relative": "..."
            },
            "JerNominal": {
              "tagNamePtResolution": "...",
              "tagNameJerScaleFactor": "...",
              "tagNameJerSFUncertainty": "..."
            }
          }
        }
      }
    """
    if year not in base:
        raise RuntimeError(f"Year key not found in JSON: {year}")

    y = base[year]
    jercJsonPath = y["jercJsonPath"]

    if isData:
        apply_data = y.get("ApplyOnData")
        if apply_data is None:
            raise RuntimeError(
                f"Requested data but no 'ApplyOnData' section for year: {year}"
            )

        jes_nom = apply_data.get("JesNominal")
        if jes_nom is None:
            raise RuntimeError(
                f"No 'JesNominal' under 'ApplyOnData' for year: {year}"
            )

        if era is None:
            raise RuntimeError(
                f"Data requested but no era provided for year: {year}"
            )

        if era not in jes_nom:
            raise RuntimeError(
                f"Era '{era}' not found under ApplyOnData.JesNominal "
                f"for year: {year}"
            )

        e_cfg = jes_nom[era]
        tagNameL1FastJet = e_cfg["tagNameL1FastJet"]
        tagNameL2Relative = e_cfg["tagNameL2Relative"]
        tagNameL2L3Residual = e_cfg.get("tagNameL2L3Residual")

        # Data: we normally do not use JER resolution / scale factor
        tagNamePtResolution = None
        tagNameJerScaleFactor = None
        tagNameJerSFUncertainty = None

    else:
        apply_mc = y.get("ApplyOnMC")
        if apply_mc is None:
            raise RuntimeError(
                f"Requested MC but no 'ApplyOnMC' section for year: {year}"
            )

        jes_nom = apply_mc.get("JesNominal")
        if jes_nom is None:
            raise RuntimeError(
                f"No 'JesNominal' under 'ApplyOnMC' for year: {year}"
            )

        jer_nom = apply_mc.get("JerNominal")  # may be None if JER disabled

        tagNameL1FastJet = jes_nom["tagNameL1FastJet"]
        tagNameL2Relative = jes_nom["tagNameL2Relative"]
        # MC: no residual corrections in the new layout
        tagNameL2L3Residual = None

        tagNamePtResolution = (
            jer_nom["tagNamePtResolution"] if jer_nom is not None else None
        )
        tagNameJerScaleFactor = (
            jer_nom["tagNameJerScaleFactor"] if jer_nom is not None else None
        )
        tagNameJerSFUncertainty = (
            jer_nom.get("tagNameJerSFUncertainty") if jer_nom is not None else None
        )

    return Tags(
        jercJsonPath=jercJsonPath,
        tagNameL1FastJet=tagNameL1FastJet,
        tagNameL2Relative=tagNameL2Relative,
        tagNameL2L3Residual=tagNameL2L3Residual,
        tagNamePtResolution=tagNamePtResolution,
        tagNameJerScaleFactor=tagNameJerScaleFactor,
        tagNameJerSFUncertainty=tagNameJerSFUncertainty,
    )


# ============================================================
# CorrectionRefs: wrap correctionlib objects
# ============================================================


class CorrectionRefs:
    """
    Lightweight handle on the correctionlib nodes needed for JEC/JER.

    Backwards-compatible scalar attributes:

      * cs        : CorrectionSet
      * cL1       : L1FastJet
      * cL2       : L2Relative
      * cL2L3Res  : residual corrections (data only; may be None for MC)
      * cReso     : JER pT resolution (MC; None for data)
      * cJerSF    : JER scale factors (MC; None for data)
      * cJerSFUnc : JER SF uncertainties (MC; None for data / older payloads)
      * cJerSmear : year-agnostic JER smearing node (jer_smear.json.gz)

    Struct-like views 

      * jes : JesDataRefs  (for data)  or JesMcRefs (for MC)
      * jer : JerMcRefs    (for MC)    or None      (for data)
    """

    def __init__(self, tags: Tags, isData: bool):
        # --- Load the base CorrectionSet ---
        self.cs: correction.CorrectionSet = correction.CorrectionSet.from_file(
            tags.jercJsonPath
        )

        # --- JEC pieces ---
        self.cL1: correction.Correction = self.cs[tags.tagNameL1FastJet]
        self.cL2: correction.Correction = self.cs[tags.tagNameL2Relative]

        # Residual corrections (data only) may be absent for MC
        if tags.tagNameL2L3Residual:
            self.cL2L3Res: Optional[correction.Correction] = self.cs[
                tags.tagNameL2L3Residual
            ]
        else:
            self.cL2L3Res = None

        # --- JER pieces (MC only) ---
        if tags.tagNamePtResolution:
            self.cReso: Optional[correction.Correction] = self.cs[
                tags.tagNamePtResolution
            ]
        else:
            self.cReso = None

        if tags.tagNameJerScaleFactor:
            self.cJerSF: Optional[correction.Correction] = self.cs[
                tags.tagNameJerScaleFactor
            ]
        else:
            self.cJerSF = None

        if tags.tagNameJerSFUncertainty:
            try:
                self.cJerSFUnc: Optional[correction.Correction] = self.cs[
                    tags.tagNameJerSFUncertainty
                ]
            except Exception as e:
                print(
                    "[JecConfig] Warning: JER SF uncertainty correction "
                    f"'{tags.tagNameJerSFUncertainty}' is not available in "
                    f"{tags.jercJsonPath}. JER up/down variations will fall "
                    f"back to the scale-factor payload only: {e}"
                )
                self.cJerSFUnc = None
        else:
            self.cJerSFUnc = None

        # Year-agnostic JER smearing node (separate JSON)
        self.cJerSmear: correction.Correction = get_jer_smear()

        # ----------------------------------------------------
        # Struct-like views (ergonomic access)
        # ----------------------------------------------------
        if isData:
            # Data: has L1, L2, L2L3Residual; no JER
            if self.cL2L3Res is None:
                raise RuntimeError(
                    "Data configuration without L2L3Residual is unexpected."
                )
            self.jes: JesDataRefs = JesDataRefs(
                cs=self.cs,
                l1FastJet=self.cL1,
                l2Relative=self.cL2,
                l2L3Residual=self.cL2L3Res,
            )
        else:
            # MC: L1, L2; residuals unused in new layout
            self.jes: JesMcRefs = JesMcRefs(
                cs=self.cs,
                l1FastJet=self.cL1,
                l2Relative=self.cL2,
            )

        if self.cReso is not None and self.cJerSF is not None:
            self.jer: Optional[JerMcRefs] = JerMcRefs(
                cs=self.cs,
                ptResolution=self.cReso,
                scaleFactor=self.cJerSF,
                sfUncertainty=self.cJerSFUnc,
            )
        else:
            self.jer = None


# ============================================================
# Systematic tag helpers (new JSON layout)
# ============================================================


def get_syst_tag_names(
    year_json: dict,
    year: str,
) -> Dict[str, List[Tuple[str, str]]]:
    """
    Read JES uncertainty tag names from the new JSON layout.

    Input: the full config JSON (e.g. CorrCache.get_cfg("AK4")), and a year key.

    New structure:

      year_json[year]["ApplyOnMC"]["JesUncertaintySet"] = {
        "JesUncertaintySetFull": {
          "CMS_scale_j_AbsoluteMPFBias": "...fullTagName...",
          ...
        },
        "JesUncertaintySetReduced": {
          "CMS_scale_j_FlavorQCD": "...",
          ...
        },
        "JesUncertaintySetTotal": {
          "CMS_scale_j_Total": "..."
        }
      }

    Returns:
      dict:
        setName -> list of (fullTagName, cmsBaseName)
        e.g. "JesUncertaintySetFull" -> [
                 ("Summer19UL16APV_V7_MC_AbsoluteMPFBias_AK4PFchs",
                  "CMS_scale_j_AbsoluteMPFBias"),
                 ...
              ]
    """
    out: Dict[str, List[Tuple[str, str]]] = {}

    if year not in year_json:
        return out

    y = year_json[year]
    apply_mc = y.get("ApplyOnMC", {})
    jes_sets = apply_mc.get("JesUncertaintySet", {})

    for setName, mapping in jes_sets.items():
        # mapping: cmsName -> fullTagName
        pairs: List[Tuple[str, str]] = []
        for cms_name, full_tag in mapping.items():
            # We want (fullTagName, cmsName), matching earlier interface
            pairs.append((full_tag, cms_name))
        out[setName] = pairs

    return out


def get_jer_sets(year_json: dict, year: str) -> Dict[str, List[JerBin]]:
    """
    Read JER uncertainty bin definitions from the new JSON layout.

    New structure:

      year_json[year]["ApplyOnMC"]["JerUncertaintySet"] = {
        "JerUncertaintySetFull": {
          "CMS_res_j_eta0to1p93_pt0toInf_2016": [0.0, 1.93, 0, 10000],
          ...
        },
        "JerUncertaintySetTotal": {
          "CMS_res_j_2016": [0.0, 6.0, 0, 10000]
        }
      }

    Returns:
      dict:
        setName -> list[JerBin]
    """
    out: Dict[str, List[JerBin]] = {}

    if year not in year_json:
        return out

    y = year_json[year]
    apply_mc = y.get("ApplyOnMC", {})
    jer_sets = apply_mc.get("JerUncertaintySet", {})

    for setName, obj in jer_sets.items():
        bins: List[JerBin] = []
        for label, arr in obj.items():
            if isinstance(arr, list) and len(arr) == 4:
                etaMin, etaMax, ptMin, ptMax = arr
                bins.append(JerBin(label, etaMin, etaMax, ptMin, ptMax))
        out[setName] = bins

    return out


# ============================================================
# Helpers to build C++-style "uncertainty" structs
# ============================================================


def build_jes_unc_set_refs(
    cs: correction.CorrectionSet,
    year_json: dict,
    year: str,
) -> JesUncSetRefs:
    """
    Build a JesUncSetRefs POD from the config for a given year.

    'full', 'reduced', 'total' are keyed by the CMS uncertainty name,
    with values being correctionlib Correction objects.
    """
    tag_sets = get_syst_tag_names(year_json, year)

    full: Dict[str, correction.Correction] = {}
    reduced: Dict[str, correction.Correction] = {}
    total: Dict[str, correction.Correction] = {}

    for setName, pairs in tag_sets.items():
        # Decide which bucket this set belongs to
        if "Full" in setName:
            bucket = full
        elif "Reduced" in setName:
            bucket = reduced
        elif "Total" in setName:
            bucket = total
        else:
            # Unknown naming; you can extend this logic if needed
            continue

        for full_tag, cms_name in pairs:
            bucket[cms_name] = cs[full_tag]

    return JesUncSetRefs(full=full, reduced=reduced, total=total)


def build_jer_unc_sets(
    year_json: dict,
    year: str,
) -> JerUncSets:
    """
    Build a JerUncSets POD from the config for a given year.

    'full' and 'total' are maps:
      label (CMS name) -> JerBin
    """
    raw_sets = get_jer_sets(year_json, year)

    full: Dict[str, JerBin] = {}
    total: Dict[str, JerBin] = {}

    for setName, bins in raw_sets.items():
        if "Full" in setName:
            bucket = full
        elif "Total" in setName:
            bucket = total
        else:
            continue

        for b in bins:
            bucket[b.label] = b

    return JerUncSets(full=full, total=total)
