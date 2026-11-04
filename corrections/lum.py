"""LUM: pileup weights, as in LUM/examples lumExample.py (Pileup_nTrueInt, "nominal"/"up"/"down")."""
import correctionlib

POG_DIR = "/cvmfs/cms-griddata.cern.ch/cat/metadata"
LUM_VERSION = "2024-01-31"
PU_KEY = {
    "2022": "Collisions2022_355100_357900_eraBCD_GoldenJson",
    "2022EE": "Collisions2022_359022_362760_eraEFG_GoldenJson",
}


def pu_weights(campaign, era, n_true_int):
    """{"nom", "up", "down"} -> per-event weight (JSON flow=clamp, nTrueInt passed as is)."""
    pu = correctionlib.CorrectionSet.from_file(
        f"{POG_DIR}/LUM/{campaign}/{LUM_VERSION}/puWeights.json.gz")[PU_KEY[era]]
    return {s: pu.evaluate(n_true_int, v) for s, v in (("nom", "nominal"), ("up", "up"), ("down", "down"))}
