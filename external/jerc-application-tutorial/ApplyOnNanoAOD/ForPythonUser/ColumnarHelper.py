from typing import Optional, Tuple
import math

import awkward as ak
import numpy as np

"""
ColumnarHelper
--------------
Helper utilities for:

  * Selecting a good "debug event" in the columnar world 
  * Safely evaluating correctionlib corrections on awkward arrays, even when
    the correctionlib build does NOT support NumPy vectorization.
  * Broadcasting mixed inputs (scalars, lists, awkward arrays) to a common
    jagged shape.

These helpers are deliberately standalone so they can be reused in other
columnar examples or tutorials.
"""


# ======================================================================
# Debug-event helper
# ======================================================================

def find_debug_event_index(
    ak4_sel,
    ak8_sel,
    apply_only_ak4: bool,
    apply_only_ak8: bool,
    apply_both: bool,
) -> Optional[int]:
    """
    Return the first event index suitable for debug printing.

    The logic: 

      * MODE_BOTH  : event must have >= 1 AK8 AND >= 1 AK4
                     (caller should pass the non-overlap AK4 selection here)
      * MODE_AK4   : event must have >= 1 AK4
      * MODE_AK8   : event must have >= 1 AK8

    Parameters
    ----------
    ak4_sel :
        Awkward array of per-jet booleans, jagged by event, or None.
        Typically the jet selection mask used for AK4.
    ak8_sel :
        Awkward array of per-jet booleans, jagged by event, or None.
        Typically the jet selection mask used for AK8.
    apply_only_ak4 : bool
        True if running only the AK4 jet path.
    apply_only_ak8 : bool
        True if running only the AK8 jet path.
    apply_both : bool
        True if running both AK4 and AK8 jet paths.

    Returns
    -------
    int or None
        The first event index matching the requested mode, or None
        if no suitable event is found.
    """

    # Determine nominal event count from whichever selector is available
    n4 = len(ak4_sel) if ak4_sel is not None else None
    n8 = len(ak8_sel) if ak8_sel is not None else None
    n_events = n4 if n4 is not None else n8

    if n_events is None or n_events == 0:
        return None

    # If both provided but lengths differ, align conservatively
    if (n4 is not None and n8 is not None) and (n4 != n8):
        n_events = min(n4, n8)
        if ak4_sel is not None:
            ak4_sel = ak4_sel[:n_events]
        if ak8_sel is not None:
            ak8_sel = ak8_sel[:n_events]

    def event_mask(sel, n: int) -> np.ndarray:
        """
        Turn a jagged boolean per-jet selection (or None) into a flat
        per-event boolean mask.
        """
        if sel is None:
            return np.zeros(n, dtype=bool)

        try:
            # ak.any over jets per event; fill None→False; then convert to numpy
            m = ak.fill_none(ak.any(sel, axis=1), False)
            return ak.to_numpy(m)
        except Exception:
            # Defensive fallback: if anything goes wrong, pretend nothing passes
            return np.zeros(n, dtype=bool)

    ak4_mask = event_mask(ak4_sel, n_events)
    ak8_mask = event_mask(ak8_sel, n_events)

    # Mode-dependent requirement
    if apply_both:
        candidate_mask = ak4_mask & ak8_mask
    elif apply_only_ak4:
        candidate_mask = ak4_mask
    elif apply_only_ak8:
        candidate_mask = ak8_mask
    else:
        # Defensive fallback: any event with at least one selected jet
        candidate_mask = ak4_mask | ak8_mask

    idx = np.flatnonzero(candidate_mask)
    return int(idx[0]) if idx.size else None


# ======================================================================
# correctionlib evaluation helpers
# ======================================================================

def eval_with_string_last(corr, *num_cols, cat, blocksize: int = 200_000):
    """
    Evaluate a correctionlib node where the string categorical argument
    comes LAST in the signature.

    Example node signature:
        (eta, pt, 'up'/'down')

    Parameters
    ----------
    corr :
        correctionlib Correction or CorrectionRef.
    *num_cols :
        Numeric inputs (arrays, awkward arrays, or scalars).
    cat :
        String category, either a single string or an array-like of strings.
    blocksize : int, optional
        Chunk size for processing entries. Smaller chunks reduce peak memory,
        larger chunks are faster.

    Returns
    -------
    awkward.Array
        Awkward array with the same jagged structure as the broadcast inputs.
    """
    bcols, counts = broadcast_to_common_shape(*num_cols)
    flats = [
        ak.to_numpy(ak.flatten(x, axis=None)).astype(np.float64, copy=False)
        for x in bcols
    ]

    n = flats[0].size if flats else 0
    if n == 0:
        # Return an empty array with the right shape
        return ak.zeros_like(bcols[0] if bcols else ak.Array([[]]), dtype=np.float64)

    # Handle categorical argument (single string or per-entry array)
    if isinstance(cat, str):
        cats = np.full(n, cat, dtype=object)
    else:
        bcat, _ = broadcast_to_common_shape(cat, bcols[0])
        cats = ak.to_numpy(ak.flatten(bcat[0], axis=None)).astype(object)

    out = np.empty(n, dtype=np.float64)

    # Chunked scalar evaluation (guaranteed to work in all correctionlib builds)
    for start in range(0, n, blocksize):
        end = min(start + blocksize, n)
        for i in range(start, end):
            args = [a[i] for a in flats] + [cats[i]]
            out[i] = float(corr.evaluate(*args))

    return ak.unflatten(out, counts)


def eval_with_string_first(corr, *num_cols, cat, blocksize: int = 200_000):
    """
    Evaluate a correctionlib node where the string categorical argument
    comes FIRST in the signature.

    Example node signature:
        ('map_key', eta, phi)

    Parameters
    ----------
    corr :
        correctionlib Correction or CorrectionRef.
    *num_cols :
        Numeric inputs (arrays, awkward arrays, or scalars).
    cat :
        String category, either a single string or an array-like of strings.
    blocksize : int, optional
        Chunk size for processing entries.

    Returns
    -------
    awkward.Array
        Awkward array with the same jagged structure as the broadcast inputs.
    """
    bcols, counts = broadcast_to_common_shape(*num_cols)
    flats = [
        ak.to_numpy(ak.flatten(x, axis=None)).astype(np.float64, copy=False)
        for x in bcols
    ]

    n = flats[0].size if flats else 0
    if n == 0:
        # Return an empty array with the right shape
        return ak.zeros_like(bcols[0] if bcols else ak.Array([[]]), dtype=np.float64)

    # Handle categorical argument
    if isinstance(cat, str):
        cats = np.full(n, cat, dtype=object)
    else:
        bcat, _ = broadcast_to_common_shape(cat, bcols[0])
        cats = ak.to_numpy(ak.flatten(bcat[0], axis=None)).astype(object)

    out = np.empty(n, dtype=np.float64)

    # Chunked scalar evaluation
    for start in range(0, n, blocksize):
        end = min(start + blocksize, n)
        for i in range(start, end):
            args = [cats[i]] + [a[i] for a in flats]
            out[i] = float(corr.evaluate(*args))

    return ak.unflatten(out, counts)


def _promote_to_ak(x, template: Optional[ak.Array]) -> ak.Array:
    """
    Convert `x` to an awkward Array.

    If `x` is a scalar (Python or 0-D numpy scalar), broadcast it to the
    jagged shape of `template`. If no template is available, fall back to
    a length-1 jagged array [[x]].

    Parameters
    ----------
    x :
        Input object (scalar, list, numpy array, awkward array, ...).
    template : awkward.Array or None
        Reference jagged array defining the desired shape for broadcasting.

    Returns
    -------
    awkward.Array
        Converted / broadcasted awkward array.
    """
    if isinstance(x, ak.Array):
        return x

    # Convert 0-D numpy array to Python scalar
    if isinstance(x, np.ndarray) and x.ndim == 0:
        x = x.item()

    # Python / numpy scalar → broadcast to template
    if np.isscalar(x):
        if template is None:
            # Rare path: no template available; create a trivial jagged array
            return ak.Array([[x]])
        return ak.full_like(template, x)

    # list / ndarray / other sequence → regular awkward.Array
    return ak.Array(x)


def eval_numeric(corr, *cols, blocksize: int = 200_000):
    """
    Evaluate a numeric-only correctionlib node on awkward inputs.

    Strategy:
      1. Broadcast all inputs to a common jagged shape.
      2. Flatten to NumPy arrays of float64.
      3. Process in chunks:
         * Try vectorized `corr.evaluate(*chunk)` if available.
         * On failure, fall back to scalar loop for that chunk.
      4. Unflatten back to the original jagged structure.

    Parameters
    ----------
    corr :
        correctionlib Correction or CorrectionRef.
    *cols :
        Numeric inputs (arrays, awkward arrays, or scalars).
    blocksize : int, optional
        Chunk size for processing entries.

    Returns
    -------
    awkward.Array
        Awkward array with the same jagged structure as the broadcast inputs.
    """
    bcols, counts = broadcast_to_common_shape(*cols)

    flats = [
        ak.to_numpy(ak.flatten(x, axis=None)).astype(np.float64, copy=False)
        for x in bcols
    ]
    n = flats[0].size if flats else 0
    if n == 0:
        # Return empty array with correct type/shape
        return ak.zeros_like(bcols[0] if bcols else ak.Array([[]]), dtype=np.float64)

    out = np.empty(n, dtype=np.float64)

    # Process in chunks to keep memory/timing sane
    for start in range(0, n, blocksize):
        end = min(start + blocksize, n)
        sl = slice(start, end)
        chunk_args = [a[sl] for a in flats]

        try:
            # Preferred: vectorized evaluate if the build supports it
            out[sl] = np.asarray(corr.evaluate(*chunk_args), dtype=np.float64)
        except Exception:
            # Fallback: scalar loop (pybind always supports scalar evaluate)
            for i in range(start, end):
                out[i] = float(corr.evaluate(*(a[i] for a in flats)))

    return ak.unflatten(out, counts)


# ======================================================================
# Broadcasting helper
# ======================================================================

def broadcast_to_common_shape(*cols) -> Tuple[Tuple[ak.Array, ...], ak.Array]:
    """
    Broadcast a set of inputs to a common jagged shape.

    This helper:
      * Accepts a mix of scalars, lists, numpy arrays, and awkward arrays.
      * Promotes any scalar to match the jagged shape of the first
        array-like column (used as a template).
      * Uses awkward.broadcast_arrays to align everything.
      * Returns both the broadcasted columns and per-event counts.

    Parameters
    ----------
    *cols :
        Inputs to broadcast.

    Returns
    -------
    (tuple[awkward.Array, ...], awkward.Array)
        * bcols  : broadcast awkward arrays.
        * counts : per-event counts (ak.num along axis=1) from the first column.
    """
    # Find a template jagged array:
    #   - prefer an existing awkward.Array,
    #   - otherwise try to convert the first array-like to awkward.
    template = None
    for c in cols:
        if isinstance(c, ak.Array):
            template = c
            break
        if isinstance(c, (list, tuple, np.ndarray)):
            try:
                candidate = ak.Array(c)
                template = candidate
                break
            except Exception:
                # If conversion fails, keep looking
                pass

    # Convert everything with scalar promotion
    ak_cols = [_promote_to_ak(c, template) for c in cols]

    # Awkward broadcasting to a common jagged structure
    bcols = ak.broadcast_arrays(*ak_cols)

    # Per-event counts (works for any jagged array)
    counts = ak.num(bcols[0], axis=1)

    return bcols, counts

