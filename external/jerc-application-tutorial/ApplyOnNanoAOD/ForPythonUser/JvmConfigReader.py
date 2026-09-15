"""
JvmConfigReader
----------------
Light-weight helper module to load Jet Veto Map (JVM) corrections.

Design goals:
  * Provide a *single* point where JVM config files are read.
  * Cache all loaded JSONs and CorrectionSets to avoid repeated disk I/O.
  * Expose a minimal interface (get_jvm) returning (jvm_ref, jvm_key)
    for a given year.
  * Fail gracefully if a year is not configured.

"""

import json
from functools import lru_cache
import correctionlib._core as correction


# ----------------------------------------------------------------------
# JSON caching
# ----------------------------------------------------------------------
@lru_cache(maxsize=None)
def load_json_cached(path: str) -> dict:
    """
    Load a JSON file once and cache the result permanently.

    Parameters
    ----------
    path : str
        Path to the JSON configuration file.

    Returns
    -------
    dict
        Parsed JSON dictionary.
    """
    with open(path, "r") as f:
        return json.load(f)


# ----------------------------------------------------------------------
# JVM cache manager
# ----------------------------------------------------------------------
class CorrCache:
    """
    CorrCache
    ---------
    Cache for:
      * JVM configuration JSON (via load_json_cached)
      * correctionlib CorrectionSets for JVM
      * per-year (jvm_ref, jvm_key) tuples

    Usage:
        jvm_ref, jvm_key = CorrCache.get_jvm("2022")

    If a year is not configured, (None, None) is returned.
    """

    # year → (Correction object, key string)
    _jvm = {}

    @classmethod
    def get_jvm(cls, year: str):
        """
        Retrieve JVM correction for a given year.

        Parameters
        ----------
        year : str
            Year label, e.g. "2022", "2023", "2024".

        Returns
        -------
        (correction.Correction, str)  or  (None, None)
            * jvm_ref : callable correction from correctionlib
            * jvm_key : category key used during evaluation
        """

        # Already cached? Return fast.
        if year in cls._jvm:
            return cls._jvm[year]

        # Load the JSON config file defining JVM per year
        cfgJvm = load_json_cached("../../JvmConfig.json")

        if year in cfgJvm:
            entry = cfgJvm[year]

            # Load the CorrectionSet from the JSON file on disk
            cs = correction.CorrectionSet.from_file(entry["jvmFilePath"])

            # Extract the correction object and the category name
            jvm_ref = cs[entry["jvmTagName"]]
            jvm_key = entry["jvmKeyName"]

            cls._jvm[year] = (jvm_ref, jvm_key)
        else:
            # No JVM entry for this year → return empty markers
            cls._jvm[year] = (None, None)

        return cls._jvm[year]

