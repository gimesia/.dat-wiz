"""Siemens ``.dat`` loading and ICE-dimension auto-detection.

This module wraps ``mapvbvd`` for raw TWIX parsing (``src/utils`` does not
provide a ``.dat`` loader) and exposes the data as a single labelled array
together with helpers for the UI to slice out one k-space plane.

Canonical per-plane shape handed to the recon layer: ``(Cha, Lin, Col)``
i.e. ``(n_coils, n_phase_encode_lines, n_readout_samples)``, complex64.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

try:
    import mapvbvd
except ImportError as exc:  # pragma: no cover - environment guard
    raise ImportError(
        "mapvbvd is required for .dat loading. Install with `pip install pymapvbvd` "
        "(import name is `mapvbvd`)."
    ) from exc


# K-space axes (never offered as selectors) and the two dedicated sliders.
_READOUT_AXIS = "Col"   # readout samples  -> kx
_PHASE_AXIS = "Lin"     # phase-encode      -> ky
_COIL_AXIS = "Cha"      # receive channels  -> dedicated coil slider
_SLICE_AXIS = "Sli"     # slices            -> dedicated slice slider

# Axes that are part of the k-space grid / dedicated sliders, not ICE selectors.
_NON_SELECTOR = {_READOUT_AXIS, _PHASE_AXIS, _COIL_AXIS, _SLICE_AXIS}


@dataclass
class LoadedData:
    """A labelled raw-data array plus metadata for slicing.

    Attributes
    ----------
    data:
        Complex64 array. ``dims`` gives the name of each axis in order.
    dims:
        Axis names, e.g. ``["Cha", "Lin", "Col"]`` or with extra ICE dims.
    sizes:
        Mapping ``dim_name -> size``.
    selector_dims:
        Non-singleton ICE dimensions excluding coil/slice/readout/phase, in
        display order (these become the left-column selectors).
    calib:
        Navigator / phase-correction lines ``(Cha, Lin, Col)`` complex, or
        ``None`` if the file has no ``phasecor`` data.
    """

    data: np.ndarray
    dims: List[str]
    sizes: Dict[str, int] = field(default_factory=dict)
    selector_dims: List[str] = field(default_factory=list)
    calib: np.ndarray | None = None

    # -- convenience accessors ------------------------------------------------
    @property
    def n_coils(self) -> int:
        return self.sizes[_COIL_AXIS]

    @property
    def n_lines(self) -> int:
        return self.sizes[_PHASE_AXIS]

    @property
    def n_samples(self) -> int:
        return self.sizes[_READOUT_AXIS]

    @property
    def n_slices(self) -> int:
        return self.sizes[_SLICE_AXIS]

    def default_state(self) -> Dict[str, int]:
        """A selector state with every index at 0 (incl. ``Sli``)."""
        state = {d: 0 for d in self.selector_dims}
        state[_SLICE_AXIS] = 0
        return state

    def get_plane(self, state: Dict[str, int]) -> np.ndarray:
        """Return the ``(Cha, Lin, Col)`` plane for the given selector state.

        ``state`` must contain an index for every entry in ``selector_dims``
        and for ``"Sli"``. Coil is *not* selected here — the full coil axis is
        always returned so the recon layer can do per-coil or SOS.
        """
        index = []
        kept: List[str] = []
        for name in self.dims:
            if name in (_COIL_AXIS, _PHASE_AXIS, _READOUT_AXIS):
                index.append(slice(None))
                kept.append(name)
            elif name == _SLICE_AXIS:
                index.append(int(state.get(_SLICE_AXIS, 0)))
            else:
                index.append(int(state.get(name, 0)))

        plane = self.data[tuple(index)]  # axes follow `kept` order
        order = [kept.index(n) for n in (_COIL_AXIS, _PHASE_AXIS, _READOUT_AXIS)]
        plane = np.transpose(plane, order)
        return np.ascontiguousarray(plane.astype(np.complex64))


def _ensure_axis(data: np.ndarray, dims: List[str], name: str) -> tuple[np.ndarray, List[str]]:
    """Insert a length-1 axis named ``name`` at the front if missing."""
    if name not in dims:
        data = data[np.newaxis, ...]
        dims = [name] + dims
    return data, dims


def load_dat(path: str) -> LoadedData:
    """Load a Siemens ``.dat`` file and auto-detect its ICE dimensions.

    Readout oversampling is kept (``flagRemoveOS = False``) because the
    ramp-sampled EPI trajectory needs the full readout.

    Singleton ICE dimensions are squeezed out and hidden from the UI, except
    ``Sli`` and ``Cha`` which are always materialised as dedicated axes.
    """
    twix = mapvbvd.mapVBVD(path, quiet=True)
    if isinstance(twix, (list, tuple)):
        twix = twix[-1]

    img = twix.image
    img.flagRemoveOS = False  # keep readout oversampling for ramp NUFFT

    dims = list(img.sqzDims)            # squeezed axis names, in array order
    data = np.asarray(img[""].squeeze())

    # Guarantee the two dedicated-slider axes exist even when singleton.
    data, dims = _ensure_axis(data, dims, _SLICE_AXIS)
    data, dims = _ensure_axis(data, dims, _COIL_AXIS)

    sizes = {name: int(data.shape[ax]) for ax, name in enumerate(dims)}

    selector_dims = [
        d for d in dims if d not in _NON_SELECTOR and sizes[d] > 1
    ]

    return LoadedData(
        data=data.astype(np.complex64),
        dims=dims,
        sizes=sizes,
        selector_dims=selector_dims,
        calib=_load_calib(twix),
    )


def _load_calib(twix) -> np.ndarray | None:
    """Load navigator/phase-correction lines as ``(Cha, Lin, Col)`` or None.

    The ``phasecor`` object often squeezes to a different axis order than
    ``image`` (and may carry extra Rep/Eco/Sli/Ave dims); we reduce any extra
    axes to index 0 and reorder to ``(Cha, Lin, Col)``.
    """
    if not hasattr(twix, "phasecor"):
        return None
    try:
        pc = twix.phasecor
        pc.flagRemoveOS = False
        for flag in ("flagDoRawDataCorrect", "flagDoAverage"):
            if hasattr(pc, flag):
                try:
                    setattr(pc, flag, False)
                except Exception:
                    pass

        dims = list(pc.sqzDims)
        arr = np.asarray(pc[""].squeeze())
        keep = {_COIL_AXIS, _PHASE_AXIS, _READOUT_AXIS}
        if arr.ndim == 0 or not keep.issubset(dims):
            return None

        # Collapse any non-(Cha/Lin/Col) axis to its first index.
        index = [slice(None) if name in keep else 0 for name in dims]
        arr = arr[tuple(index)]
        kept = [name for name in dims if name in keep]
        order = [kept.index(n) for n in (_COIL_AXIS, _PHASE_AXIS, _READOUT_AXIS)]
        arr = np.transpose(arr, order)
        return np.ascontiguousarray(arr.astype(np.complex64))
    except Exception:
        return None
