
from typing import Dict, List, Optional, Tuple
import awkward as ak
import numpy as np

from tqdm.auto import tqdm  # type: ignore
try:
    from tqdm.auto import tqdm  # type: ignore
    TQDM_AVAILABLE = True
    def _progress_print(message: str) -> None:
        """Print messages without breaking tqdm progress bars."""
        tqdm.write(message)

except Exception:  # pragma: no cover - fallback when tqdm is unavailable
    tqdm = None  # type: ignore
    TQDM_AVAILABLE = False

    def _progress_print(message: str) -> None:
        print(message)


def progress_print(message: str) -> None:
    """Print messages without breaking tqdm progress bars."""
    tqdm.write(message)



def _extract_event_slice(arr, event_index: int, mask=None):
    """Return the slice of ``arr`` corresponding to ``event_index`` and optional ``mask``."""

    if event_index is None:
        return arr

    data = arr
    if isinstance(arr, ak.Array):
        if len(arr) <= event_index:
            return arr
        data = arr[event_index]
        if mask is not None:
            mask_arr = mask[event_index]
            try:
                data = data[mask_arr]
            except Exception:
                pass
    else:
        np_arr = np.asarray(arr)
        if np_arr.ndim == 0:
            data = np_arr
        elif event_index < np_arr.shape[0]:
            data = np_arr[event_index]
            if mask is not None:
                mask_arr = np.asarray(mask)
                try:
                    mask_slice = mask_arr[event_index]
                    data = data[mask_slice]
                except Exception:
                    pass
        else:
            data = np_arr
    return data


def preview_values(arr, limit: int = 5):
    """Return a small list of representative values from an array-like object."""
    if limit <= 0:
        return []

    if isinstance(arr, ak.Array):
        sample = arr[: max(1, limit)]
        try:
            flat = ak.flatten(sample, axis=None)
            np_arr = ak.to_numpy(flat)
        except Exception:
            try:
                as_list = ak.to_list(sample)
            except Exception:
                return []

            collected = []

            def collect(obj):
                if len(collected) >= limit:
                    return
                if isinstance(obj, list):
                    for item in obj:
                        if len(collected) >= limit:
                            break
                        collect(item)
                else:
                    collected.append(obj)

            collect(as_list)
            return collected
    else:
        np_arr = np.asarray(arr)

    if np_arr.size == 0:
        return []

    flat_np = np_arr.reshape(-1)
    if flat_np.size == 0:
        return []
    return [round(v, 6) for v in flat_np[:limit].tolist()]


isDebug = True

def printDebug(*args, indent: int = 0, enable: Optional[bool] = None) -> None:
    """Conditionally print debug information with optional indentation."""
    if enable is None:
        enable = isDebug
    if not enable:
        return

    prefix = " " * max(indent, 0)
    message = " ".join(str(arg) for arg in args)
    if prefix:
        print(prefix + message)
    else:
        print(message)


def debug_values(
    label: str,
    arr,
    indent: int = 0,
    limit: int = 5,
    event_index: Optional[int] = None,
    mask=None,
) -> None:
    """Helper to print a label and a preview of array values when debugging."""
    values = preview_values(_extract_event_slice(arr, event_index, mask=mask), limit=limit)
    printDebug(f"{label}", values, indent=indent)


def take_by_index(arr, idx, fill=0.0):
    valid = (idx >= 0) & (idx < ak.num(arr, axis=1))
    idx_safe = ak.where(valid, idx, 0)
    picked = arr[idx_safe]
    return ak.where(valid, picked, fill), valid

# ---------------------------
# Helper: apply a value only on selected positions, keep original elsewhere
# ---------------------------
def only_nonoverlap(sel, updated, original):
    # sel, updated, original are jagged-with-identical layout
    return ak.where(sel, updated, original)



