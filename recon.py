"""NUFFT reconstruction + SOS coil combination + in-memory caching.

``src/utils`` does not provide a NUFFT recon, so the NUFFT itself is done here
with ``mrinufft`` (mirroring the established hospital-session recon scripts).
The SOS coil combination and the EPI odd-line flip *are* reused from
``src/utils``.

Caching (in-memory, two-level)
------------------------------
* ``_coil_cache[plane_key]`` -> complex coil images ``(C, Ny, Nx)``. The NUFFT
  adjoint is run once per distinct k-space plane (slice + ICE selectors); SOS
  and single-coil views are both derived from this without recomputing.
* ``_image_cache[full_key]`` -> displayed magnitude image ``(Ny, Nx)``.
  ``full_key = plane_key + (coil_index,)`` when SOS is off, or
  ``plane_key + ("SOS",)`` when SOS is on, so toggling SOS never collides with
  per-coil entries (per spec section 7).
"""

from __future__ import annotations

import os
from typing import Dict, Tuple

import numpy as np
import torch
import mrinufft

from ._kspace_utils import combine_sos  # image-domain SOS, expects (Ny, C, Nx)
from ._kspace_utils import flip_alternating_lines  # EPI odd-line flip
from .trajectory import Trajectory


PlaneKey = Tuple


class ReconEngine:
    """Reconstructs a selected k-space plane via NUFFT with caching."""

    def __init__(
        self,
        trajectory: Trajectory,
        recon_shape: Tuple[int, int] | None = None,
        backend: str | None = None,
    ) -> None:
        self.trajectory = trajectory
        # Nominal (square) recon matrix. Readout oversampling lives in the
        # trajectory, not the image grid, so NX == NY == n_lines by default.
        self.recon_shape = recon_shape or (trajectory.n_lines, trajectory.n_lines)
        self.backend = backend or os.environ.get("NUFFT_BACKEND", "finufft")

        # Operators keyed by (n_coils, acquired-line signature) since the
        # trajectory differs between full and partial-Fourier planes.
        self._operators: Dict[Tuple, object] = {}
        self._coil_cache: Dict[PlaneKey, np.ndarray] = {}   # plane -> (C,Ny,Nx)
        self._image_cache: Dict[PlaneKey, np.ndarray] = {}  # full -> (Ny,Nx)

    # -- operator management --------------------------------------------------
    def _get_operator(self, n_coils: int, acquired: Tuple[int, ...]):
        op = self._operators.get((n_coils, acquired))
        if op is None:
            samples = torch.from_numpy(self.trajectory.samples_for_lines(acquired))
            op = mrinufft.get_operator(
                self.backend,
                samples=samples,
                shape=self.recon_shape,
                n_coils=n_coils,
                density=True,
            )
            self._operators[(n_coils, acquired)] = op
        return op

    @staticmethod
    def _acquired_lines(plane: np.ndarray) -> Tuple[int, ...]:
        """Indices of PE lines that carry data (non-zero).

        Un-acquired partial-Fourier lines are exactly zero-filled by mapvbvd,
        so any line with non-zero magnitude was actually acquired.
        """
        energy = np.abs(plane).sum(axis=(0, 2))  # (Lin,)
        lines = np.nonzero(energy > 0)[0]
        if lines.size == 0:  # degenerate guard
            lines = np.arange(plane.shape[1])
        return tuple(int(i) for i in lines)

    @staticmethod
    def _coil_key(plane_key: PlaneKey, toggles: dict | None) -> Tuple:
        """Cache key for the coil images: plane + any recon-affecting toggles.

        Toggles (e.g. odd-line reversal) change the coil images themselves, so
        each toggle combination caches separately. Extend by adding entries to
        the ``toggles`` dict — they fold into the key automatically.
        """
        items = tuple(sorted((k, bool(v)) for k, v in (toggles or {}).items()))
        return plane_key + items

    # -- core recon -----------------------------------------------------------
    def _coil_images(
        self, plane: np.ndarray, coil_key: Tuple, toggles: dict
    ) -> np.ndarray:
        """NUFFT-adjoint the plane to complex coil images ``(C, Ny, Nx)``."""
        cached = self._coil_cache.get(coil_key)
        if cached is not None:
            return cached

        # EPI odd-line reversal (toggle) so all lines share the same kx
        # direction. On -> flipped data matches the zigzag trajectory (correct
        # image); off -> raw data with the zigzag trajectory (N/2-ghost "raw").
        if toggles.get("reverse_odd_lines", True):
            data_src = flip_alternating_lines(plane)  # (C, Lin, Col)
        else:
            data_src = plane
        n_coils = data_src.shape[0]

        # Use only the acquired lines so the trajectory matches the sampled
        # k-space (partial Fourier: drop zero-filled lines from both).
        acquired = self._acquired_lines(plane)
        data = data_src[:, acquired, :]  # (C, n_acq, Col)
        kdata = torch.from_numpy(data.reshape(n_coils, -1).astype(np.complex64))

        op = self._get_operator(n_coils, acquired)
        coil_imgs = op.adj_op(kdata)
        if isinstance(coil_imgs, torch.Tensor):
            coil_imgs = coil_imgs.cpu().numpy()
        coil_imgs = np.asarray(coil_imgs).reshape(n_coils, *self.recon_shape)

        self._coil_cache[coil_key] = coil_imgs
        return coil_imgs

    def reconstruct(
        self,
        plane: np.ndarray,
        plane_key: PlaneKey,
        sos: bool,
        coil_index: int,
        toggles: dict | None = None,
    ) -> Tuple[np.ndarray, bool]:
        """Reconstruct one plane.

        Returns ``(magnitude_image, was_cached)``.
        """
        toggles = toggles or {}
        coil_key = self._coil_key(plane_key, toggles)
        full_key = coil_key + (("SOS",) if sos else (int(coil_index),))

        cached_img = self._image_cache.get(full_key)
        if cached_img is not None:
            return cached_img, True

        coil_imgs = self._coil_images(plane, coil_key, toggles)  # (C, Ny, Nx)

        if sos:
            # combine_sos expects (Ny, C, Nx); coil_imgs is (C, Ny, Nx).
            img = combine_sos(np.transpose(coil_imgs, (1, 0, 2)), coil_axis=1)
        else:
            c = int(np.clip(coil_index, 0, coil_imgs.shape[0] - 1))
            img = np.abs(coil_imgs[c])

        img = np.asarray(img, dtype=np.float64)
        self._image_cache[full_key] = img
        return img, False

    def has_coils(self, plane_key: PlaneKey, toggles: dict | None = None) -> bool:
        """True if the all-coil recon for this plane+toggles is already computed.

        Lets the UI re-derive SOS / single-coil views instantly (cheap combine)
        without triggering an expensive NUFFT when the user only toggles a
        display option.
        """
        return self._coil_key(plane_key, toggles) in self._coil_cache

    def clear_cache(self) -> None:
        self._coil_cache.clear()
        self._image_cache.clear()
