"""Calibration / navigator diagnostics window.

Shows the EPI navigator (phasecor) lines the way the recon pipelines do
(cf. Section 5b of the hospital-session recon scripts):

* raw navigator k-space magnitude per readout sample,
* navigator magnitude in image (x) space,
* the phase difference between the negative- and positive-kx navigators,
  which should look like a roughly linear ramp (the source of the EPI
  N/2-ghost phase error).

The image-space transform here is a plain centred 1-D IFFT along the readout.
The navigators are ramp-sampled, so this is an *approximate* hybrid transform
(the full pipeline grids with a 1-D NUFFT); it is faithful enough for a visual
diagnostic and needs no trajectory.
"""

from __future__ import annotations

import numpy as np
from PyQt5 import QtWidgets

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


def _ift_readout(x: np.ndarray) -> np.ndarray:
    """Centred 1-D inverse FFT along the last (readout) axis."""
    return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(x, axes=-1), axis=-1), axes=-1)


def _mask_low_signal(phase: np.ndarray, mag: np.ndarray, frac: float = 0.08) -> np.ndarray:
    """Blank (NaN) phase where magnitude is below ``frac`` of its peak.

    Phase is meaningless in near-zero-signal regions; masking keeps the phase
    panels readable instead of full of edge noise.
    """
    out = np.asarray(phase, dtype=float).copy()
    out[mag < frac * (mag.max() + 1e-12)] = np.nan
    return out


class CalibrationWindow(QtWidgets.QMainWindow):
    """A standalone window plotting the navigator/calibration diagnostics."""

    def __init__(self, calib: np.ndarray, coil_index: int = 0, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Calibration / navigator diagnostics")
        self.calib = calib  # (Cha, Lin, Col)
        self.resize(820, 720)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        v = QtWidgets.QVBoxLayout(central)

        # -- controls -----------------------------------------------------
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Coil"))
        self.coil_spin = QtWidgets.QSpinBox()
        self.coil_spin.setRange(0, calib.shape[0] - 1)
        self.coil_spin.setValue(int(np.clip(coil_index, 0, calib.shape[0] - 1)))
        self.coil_spin.valueChanged.connect(self._replot)
        row.addWidget(self.coil_spin)
        self.sos_check = QtWidgets.QCheckBox("Coil-combined (SOS)")
        self.sos_check.stateChanged.connect(self._replot)
        row.addWidget(self.sos_check)
        row.addStretch(1)
        v.addLayout(row)

        # -- figure -------------------------------------------------------
        self.fig = Figure(figsize=(8, 7), tight_layout=True)
        self.canvas = FigureCanvas(self.fig)
        v.addWidget(NavigationToolbar(self.canvas, self))
        v.addWidget(self.canvas, stretch=1)

        self._replot()

    # --------------------------------------------------------------------
    def _replot(self, *_) -> None:
        calib = self.calib
        n_coils, n_lines, n_col = calib.shape
        sos = self.sos_check.isChecked()          # coil-combine the magnitudes
        c = self.coil_spin.value()                # coil for the phase panels
        hyb = _ift_readout(calib)                 # (C, L, Col), image space

        self.fig.clear()
        # Left column = magnitude analysis, right column = phase analysis.
        ax_kmag = self.fig.add_subplot(2, 2, 1)   # top-left
        ax_imag = self.fig.add_subplot(2, 2, 3)   # bottom-left
        ax_iph = self.fig.add_subplot(2, 2, 2)    # top-right
        ax_diff = self.fig.add_subplot(2, 2, 4)   # bottom-right

        # --- LEFT: magnitudes -------------------------------------------
        for n in range(n_lines):
            if sos:
                km = np.sqrt((np.abs(calib[:, n, :]) ** 2).sum(0))
                im = np.sqrt((np.abs(hyb[:, n, :]) ** 2).sum(0))
            else:
                km = np.abs(calib[c, n, :])
                im = np.abs(hyb[c, n, :])
            ax_kmag.plot(km, label=f"nav {n}")
            ax_imag.plot(im, label=f"nav {n}")

        ax_kmag.set_title("k-space magnitude")
        ax_kmag.set_xlabel("readout sample")
        ax_kmag.grid(True, alpha=0.3)
        ax_kmag.legend(fontsize=8)

        ax_imag.set_title("image-space magnitude  (|IFFT| vs x)")
        ax_imag.set_xlabel("x pixel")
        ax_imag.grid(True, alpha=0.3)
        ax_imag.legend(fontsize=8)

        # --- RIGHT: phases (selected coil) ------------------------------
        # Per-navigator image-space phase, masked to the signal region.
        for n in range(n_lines):
            hc = hyb[c, n, :]
            ax_iph.plot(_mask_low_signal(np.angle(hc), np.abs(hc)),
                        ".-", ms=3, label=f"nav {n}")
        ax_iph.set_title(f"image-space phase (coil {c})")
        ax_iph.set_xlabel("x pixel")
        ax_iph.set_ylabel("phase [rad]")
        ax_iph.set_ylim(-np.pi, np.pi)
        ax_iph.grid(True, alpha=0.3)
        ax_iph.legend(fontsize=8)

        # Phase difference: negative-kx nav (1) vs positive-kx navs (evens).
        if n_lines >= 2:
            if sos:
                pos = hyb[:, 0::2, :].mean(1)                # (C, Col)
                neg = hyb[:, 1, :]
                w = np.abs(pos)
                pc = (neg * np.conj(pos) * w).sum(0) / (w.sum(0) + 1e-8)
                pos_mag = np.sqrt((np.abs(pos) ** 2).sum(0))
                lbl = "coil-combined"
            else:
                pos = hyb[c, 0::2, :].mean(0)                # (Col,)
                neg = hyb[c, 1, :]
                pc = neg * np.conj(pos)
                pos_mag = np.abs(pos)
                lbl = f"coil {c}"
            pc = pc / (np.abs(pc) + 1e-8)                    # unit phasor
            ax_diff.plot(_mask_low_signal(np.angle(pc), pos_mag), "b.-", ms=3)
            ax_diff.set_ylim(-np.pi, np.pi)
            ax_diff.set_title(f"phase difference: neg vs pos ({lbl})")
            ax_diff.set_xlabel("x pixel")
            ax_diff.set_ylabel("phase [rad]")
            ax_diff.grid(True, alpha=0.3)
        else:
            ax_diff.text(0.5, 0.5, "phase difference needs ≥ 2 navigator lines",
                         ha="center", va="center", transform=ax_diff.transAxes)
            ax_diff.axis("off")

        energy = np.sqrt((np.abs(calib) ** 2).sum(axis=(0, 2)))
        energy = energy / (energy.max() + 1e-12)
        mag_where = "SOS" if sos else f"coil {c}"
        self.fig.suptitle(
            f"Magnitude (left, {mag_where})  ·  Phase (right, coil {c})   —   "
            f"{n_lines} navs, {n_coils} coils, {n_col} samples   "
            f"energy: {np.round(energy, 2).tolist()}",
            fontsize=10,
        )
        self.canvas.draw_idle()
