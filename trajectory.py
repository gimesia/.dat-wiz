"""k-space trajectory extraction from a Pulseq ``.seq`` file.

Design note (flagged assumption)
--------------------------------
For a standard EPI readout the *only* genuinely non-uniform sampling is along
the ramp-sampled readout (kx). The phase-encode positions (ky) sit on a uniform
Cartesian grid set by constant PE blips. We therefore build the 2-D trajectory
from:

* **kx** — the ramp-sampled readout trajectory of a single ADC line, obtained
  via ``utils_kspace.get_kx_and_ramp_info`` (reused from ``src/utils``);
* **ky** — a uniform grid of ``n_lines`` positions in ``[-0.5, 0.5)``.

This is robust across echoes / repetitions / diffusion directions, because the
EPI readout gradient (hence kx) does not change between them — only the
diffusion-weighting gradients (played before readout) differ, and those do not
shift the readout k-space. The trajectory is extracted **once per session**.

If a future dataset genuinely varies its readout trajectory per selector
(e.g. blip-up/blip-down), specialise :meth:`Trajectory.for_selector`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
import pypulseq as pp

from ._kspace_utils import get_kx_and_ramp_info


@dataclass
class Trajectory:
    """A normalised 2-D EPI trajectory plus the metadata used to build it.

    Attributes
    ----------
    kx_norm:
        Readout k positions, shape ``(n_samples,)``, normalised to ``[-0.5, 0.5]``.
    n_lines:
        Number of phase-encode lines (ky grid size).
    n_samples:
        Number of readout samples per line.
    ramp_sampled:
        ``True`` if a ramp-sampled trapezoid readout was detected in the .seq;
        ``False`` if we fell back to a uniform (Cartesian) kx.
    """

    kx_norm: np.ndarray
    n_lines: int
    n_samples: int
    ramp_sampled: bool

    def _ky_norm(self) -> np.ndarray:
        idx = np.arange(self.n_lines)
        return (idx - self.n_lines // 2) / self.n_lines  # [-0.5, 0.5)

    def samples_for_lines(self, line_indices) -> np.ndarray:
        """Trajectory for a subset of PE lines, shape ``(len*n_samples, 2)``.

        ``line_indices`` are indices into the **full** ``n_lines`` grid, so each
        selected line keeps its true ky position and its zigzag parity. This is
        what makes partial-Fourier work: only the acquired (non-zero) lines are
        passed, and the trajectory matches exactly those ky positions — the
        un-acquired lines are excluded from the NUFFT rather than fed as zeros.

        The readout is a **zigzag**: even PE lines traverse +kx, odd lines −kx
        (EPI acquisition order). The recon layer flips the odd *data* lines to
        match — the combination verified to reconstruct cleanly against the
        reference pipeline. A monotonic kx with flipped data (or zigzag with
        un-flipped data) produces a half-FOV readout artefact.

        Sample ordering is line-major to match the data plane flattened with
        ``reshape(-1)``. Columns are ``[kx, ky]``, dtype float32.
        """
        ky = self._ky_norm()
        rows = []
        for line in line_indices:
            kxl = self.kx_norm if (int(line) % 2 == 0) else self.kx_norm[::-1]
            rows.append(np.stack([kxl, np.full(self.n_samples, ky[int(line)])], axis=-1))
        return np.concatenate(rows, axis=0).astype(np.float32)

    def samples(self) -> np.ndarray:
        """Full-grid trajectory, shape ``(n_lines*n_samples, 2)``."""
        return self.samples_for_lines(range(self.n_lines))

    def for_selector(self, state: Dict[str, int]) -> "Trajectory":
        """Hook for per-selector trajectories. Currently returns ``self``.

        The readout trajectory is shared across all ICE selectors for standard
        EPI (see module docstring). Override here if a dataset re-derives the
        trajectory per repetition / diffusion direction.
        """
        return self


def _resample_to(arr: np.ndarray, length: int) -> np.ndarray:
    """Linearly resample a 1-D array to ``length`` points."""
    if len(arr) == length:
        return arr
    src = np.linspace(0.0, 1.0, len(arr))
    dst = np.linspace(0.0, 1.0, length)
    return np.interp(dst, src, arr)


def extract_trajectory(seq_path: str, n_lines: int, n_samples: int) -> Trajectory:
    """Read a ``.seq`` file and build the shared EPI trajectory.

    Parameters
    ----------
    seq_path:
        Path to the Pulseq ``.seq`` file.
    n_lines, n_samples:
        Phase-encode and readout sizes from the loaded ``.dat`` (so the
        trajectory is matched to the data grid; kx is resampled if the .seq
        ADC sample count differs).
    """
    seq = pp.Sequence()
    seq.read(seq_path)

    try:
        kx_adc, *_ = get_kx_and_ramp_info(seq)
        kx = np.asarray(kx_adc, dtype=np.float64)
        ramp_sampled = True
    except Exception:
        # Fall back to a uniform readout (Cartesian kx) so the tool still works.
        kx = np.linspace(-1.0, 1.0, n_samples)
        ramp_sampled = False

    kx = _resample_to(kx, n_samples)

    # get_kx_and_ramp_info integrates only the readout gradient, so kx runs
    # 0 -> kmax (monotonic, all positive). The actual EPI readout is prewound
    # to start at -kmax, so re-centre on the midpoint before normalising.
    kx = kx - 0.5 * (kx.max() + kx.min())

    max_abs = np.max(np.abs(kx))
    if max_abs > 0:
        kx_norm = kx / max_abs * 0.5  # -> [-0.5, 0.5]
    else:
        kx_norm = kx

    return Trajectory(
        kx_norm=kx_norm.astype(np.float32),
        n_lines=int(n_lines),
        n_samples=int(n_samples),
        ramp_sampled=ramp_sampled,
    )
