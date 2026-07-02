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
        calib: np.ndarray | None = None,
    ) -> None:
        self.trajectory = trajectory
        # Nominal (square) recon matrix. Readout oversampling lives in the
        # trajectory, not the image grid, so NX == NY == n_lines by default.
        self.recon_shape = recon_shape or (trajectory.n_lines, trajectory.n_lines)
        self.backend = backend or os.environ.get("NUFFT_BACKEND", "finufft")
        self.calib = calib  # navigator lines (Cha, Lin, Col) for ghost correction

        # Operators keyed by (n_coils, acquired-line signature) since the
        # trajectory differs between full and partial-Fourier planes.
        self._operators: Dict[Tuple, object] = {}
        self._coil_cache: Dict[PlaneKey, np.ndarray] = {}   # plane -> (C,Ny,Nx)
        self._image_cache: Dict[PlaneKey, np.ndarray] = {}  # full -> (Ny,Nx)
        # Ghost-correction: 1-D readout NUFFT + navigator phase map per coil count.
        self._gc_ops: Dict[int, object] = {}
        self._phase_corr: Dict[int, object] = {}

    @property
    def can_ghost_correct(self) -> bool:
        """True if navigator data is available for ghost correction."""
        return self.calib is not None and self.calib.shape[1] >= 2

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

    # -- ghost correction (navigator / Nunes) --------------------------------
    def _gc_operator(self, n_coils: int):
        """1-D readout NUFFT operator (n_samples -> NX) for ghost correction."""
        op = self._gc_ops.get(n_coils)
        if op is None:
            kx = self.trajectory.kx_norm.astype(np.float32)[:, None]  # +kx ramp
            op = mrinufft.get_operator(
                self.backend,
                samples=torch.from_numpy(kx),
                shape=(self.recon_shape[1],),
                n_coils=n_coils,
                density=True,
            )
            self._gc_ops[n_coils] = op
        return op

    def _get_phase_corr(self, n_coils: int):
        """Navigator phase-correction map ``exp(i(phi_neg - phi_pos))`` (NX,).

        Estimated once (per coil count) from the navigators the way the recon
        pipelines do: grid each nav to hybrid (x) space with the +kx operator
        (mapvbvd already reverses the negative-kx readouts), then take the
        coil-weighted phase difference of the odd (neg) nav vs the averaged even
        (pos) navs.
        """
        pc = self._phase_corr.get(n_coils)
        if pc is not None:
            return pc
        if not self.can_ghost_correct:
            return None
        op = self._gc_operator(n_coils)
        nav = torch.from_numpy(self.calib.astype(np.complex64))  # (C, L, Col)
        n_lines = nav.shape[1]
        hybrids = []
        for n in range(n_lines):
            h = op.adj_op(nav[:, n, :])
            hybrids.append(h if isinstance(h, torch.Tensor) else torch.as_tensor(h))
        hybrid_pos = torch.stack(hybrids[0::2]).mean(0)  # avg even navs (C, NX)
        hybrid_neg = hybrids[1]                          # odd nav        (C, NX)
        phase_diff = hybrid_neg * hybrid_pos.conj()
        w = hybrid_pos.abs()
        pc = (phase_diff * w).sum(0) / (w.sum(0) + 1e-8)  # (NX,)
        pc = pc / (pc.abs() + 1e-8)                       # unit phasor
        self._phase_corr[n_coils] = pc
        return pc

    def _apply_ghost_correction(self, plane: np.ndarray) -> np.ndarray:
        """Phase-correct the odd (originally negative-kx) lines in hybrid space.

        For each odd PE line: NUFFT-adjoint to image (x) space, multiply by
        ``phase_corr.conj()``, NUFFT-forward back to k-space. Zero-filled lines
        stay zero. Operates on the raw (pre-flip) plane so the correction is
        baked in before the usual reversal + zigzag NUFFT.
        """
        pc = self._get_phase_corr(plane.shape[0])
        if pc is None:
            return plane
        op = self._gc_operator(plane.shape[0])
        corr = pc.conj()
        result = torch.from_numpy(plane.astype(np.complex64)).clone()  # (C,Lin,Col)
        for line in range(1, result.shape[1], 2):
            hybrid = op.adj_op(result[:, line, :])       # (C, NX)
            if not isinstance(hybrid, torch.Tensor):
                hybrid = torch.as_tensor(hybrid)
            back = op.op(hybrid * corr)                  # (C, Col)
            result[:, line, :] = back if isinstance(back, torch.Tensor) else torch.as_tensor(back)
        return result.numpy()

    # -- core recon -----------------------------------------------------------
    def _coil_images(
        self, plane: np.ndarray, coil_key: Tuple, toggles: dict
    ) -> np.ndarray:
        """NUFFT-adjoint the plane to complex coil images ``(C, Ny, Nx)``."""
        cached = self._coil_cache.get(coil_key)
        if cached is not None:
            return cached

        # Acquired (non-zero) lines from the ORIGINAL plane, before any GC edits.
        acquired = self._acquired_lines(plane)

        # Navigator ghost correction (toggle) — applied to raw data before flip.
        if toggles.get("ghost_correction", False) and self.can_ghost_correct:
            plane = self._apply_ghost_correction(plane)

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
