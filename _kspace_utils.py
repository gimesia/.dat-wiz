"""Vendored numeric helpers — self-contained copies so this tool has no
dependency on the parent repo's ``src/utils``.

These three functions were copied verbatim from the parent project:

* ``combine_sos``            — ``src/utils/utils_coils.py``
* ``flip_alternating_lines`` — ``src/utils/utils_kspace.py``
* ``get_kx_and_ramp_info``   — ``src/utils/utils_kspace.py``

They depend only on ``numpy`` (and, for ``get_kx_and_ramp_info``, the attributes
of the passed pypulseq ``Sequence`` object). Copying them here keeps
``kspace_viewer`` importable without pulling in MRzeroCore / matplotlib / scipy.
"""

from __future__ import annotations

import numpy as np


def combine_sos(img_coils, coil_axis=1):
    """
    Sum-of-Squares coil combination.

    Parameters
    ----------
    img_coils : ndarray
        Complex coil images with shape (Ny, Ncoil, Nx).
    coil_axis : int
        Axis corresponding to coil dimension.

    Returns
    -------
    ndarray
        Magnitude image (Ny, Nx).
    """
    return np.sqrt(np.sum(np.abs(img_coils) ** 2, axis=coil_axis))


def flip_alternating_lines(kspace: np.ndarray) -> np.ndarray:
    """Flip every odd-indexed ky line along the readout (kx) axis.

    In EPI acquisitions the readout gradient alternates direction each line:
    even lines are acquired in the +kx direction, odd lines in the −kx
    direction.  Reversing the odd lines corrects this so that all lines share
    the same kx orientation before further processing.

    Parameters
    ----------
    kspace:
        Shape ``(C, Ny, Nx)``, complex.

    Returns
    -------
    np.ndarray
        Same shape and dtype; lines 1, 3, 5, … are reversed along axis 2.
    """
    if kspace.ndim != 3:
        raise ValueError(f"Expected (C, Ny, Nx) array, got shape {kspace.shape}.")
    result = kspace.copy()
    result[:, 1::2, :] = result[:, 1::2, ::-1]
    return result


def get_kx_and_ramp_info(seq, gamma=42.577e6):
    """
    Returns:
        kx_adc        : (N_adc,) full kx trajectory including ramp samples
        n_ramp_left   : number of ADC samples on the rising ramp
        n_ramp_right  : number of ADC samples on the falling ramp
        n_flat        : number of ADC samples on the flat top
    """
    grad_raster = seq.system.grad_raster_time

    for block_idx in seq.block_events:
        block = seq.get_block(block_idx)
        if block.adc is None:
            continue
        if not hasattr(block, "gx") or block.gx is None:
            continue

        gx = block.gx
        adc = block.adc

        # Build waveform
        if gx.type == "trap":
            n_rise = round(gx.rise_time / grad_raster)
            n_flat = round(gx.flat_time / grad_raster)
            n_fall = round(gx.fall_time / grad_raster)
            waveform = np.concatenate(
                [
                    np.linspace(0, gx.amplitude, n_rise),
                    np.full(n_flat, gx.amplitude),
                    np.linspace(gx.amplitude, 0, n_fall),
                ]
            )
            # Time axis of the waveform (in seconds)
            t_waveform = np.arange(len(waveform)) * grad_raster

            # Time axis of the ADC samples (relative to block start)
            t_adc = adc.delay + np.arange(adc.num_samples) * adc.dwell

            # kx at each ADC sample via interpolation of cumulative integral
            kx_full = np.cumsum(waveform) * grad_raster * gamma
            kx_adc = np.interp(t_adc, t_waveform, kx_full)

            # ── classify each ADC sample as ramp or flat ──────────────────
            t_rise_end = gx.rise_time
            t_flat_end = gx.rise_time + gx.flat_time

            on_left_ramp = t_adc < t_rise_end
            on_right_ramp = t_adc >= t_flat_end
            on_flat = ~on_left_ramp & ~on_right_ramp

            n_ramp_left = int(np.sum(on_left_ramp))
            n_ramp_right = int(np.sum(on_right_ramp))
            n_flat = int(np.sum(on_flat))

        return kx_adc, n_ramp_left, n_ramp_right, n_flat

    raise RuntimeError("No readout block found.")
