
import ColumnarHelper

def apply_jvm_flag(
    ak4_eta,
    ak4_phi,
    ak4_pt,
    ak4_id,
    ak4_chEm,
    ak4_neEm,
    jvm_ref,
    jvm_key,
):
    """
    Apply Jet Veto Map. 
    Always uses ALL AK4 jets, independent of AK4/AK8 mode.

    Parameters
    ----------
    ak4_eta, ak4_phi, ak4_pt, ak4_id, ak4_chEm, ak4_neEm : awkward Arrays
        AK4 jet collections.
    jvm_ref : correctionlib object or None
        Jet Veto Map reference (correction) object.
    jvm_key : str or None
        Category key for correctionlib evaluation.

    Returns
    -------
    vetoed : awkward.Array (boolean per event)
        True for events that contain at least one vetoed jet.
        If JVM is disabled, returns None.
    """

    import numpy as np
    import awkward as ak

    # JVM enabled?
    use_jvm = jvm_ref is not None and jvm_key is not None
    if not use_jvm:
        return None

    # -------------------------------------------------------
    # Jet selection for veto map
    # -------------------------------------------------------
    msel = (
        (abs(ak4_eta) < 5.191) &
        (abs(ak4_phi) <= np.pi) &
        (ak4_id >= 6) &
        (ak4_pt >= 15.0) &
        ((ak4_chEm + ak4_neEm) <= 0.9)
    )

    # -------------------------------------------------------
    # Evaluate JVM values
    # -------------------------------------------------------
    jvm_vals = ColumnarHelper.eval_with_string_first(
        jvm_ref,
        ak4_eta[msel],
        ak4_phi[msel],
        cat=jvm_key,
    )

    # -------------------------------------------------------
    # flagged jets → per-event boolean
    # -------------------------------------------------------
    flagged = ak.unflatten(
        ak.to_numpy(ak.flatten(jvm_vals > 0.0, axis=None)),
        ak.num(ak4_eta[msel])
    )

    vetoed = ak.any(flagged, axis=1)

    # -------------------------------------------------------
    # Logging
    # -------------------------------------------------------
    print(f"   Number of events vetoed due to JetVetoMap: {int(ak.sum(vetoed))}")

    return vetoed

