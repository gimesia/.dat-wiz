"""PyQt5 + pyqtgraph UI for the k-space viewer / NUFFT recon tool."""

from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore, QtWidgets

from .loader import LoadedData
from .recon import ReconEngine
from .trajectory import Trajectory

# Display arrays as arr[row, col] (matplotlib-like), not transposed.
pg.setConfigOptions(imageAxisOrder="row-major")

# Recon-affecting toggles: (key, label, default). The key must match what
# ReconEngine understands (see recon._coil_images). Add future toggles here and
# they appear in the panel and fold into the recon cache key automatically.
RECON_TOGGLES = [
    ("reverse_odd_lines", "Reverse odd lines", True),
]


def _gray_colormap() -> pg.ColorMap:
    return pg.ColorMap([0.0, 1.0], [(0, 0, 0, 255), (255, 255, 255, 255)])


def _make_imageview() -> pg.ImageView:
    iv = pg.ImageView()
    iv.ui.roiBtn.hide()
    iv.ui.menuBtn.hide()
    iv.setColorMap(_gray_colormap())
    return iv


class MainWindow(QtWidgets.QMainWindow):
    def __init__(
        self,
        data: LoadedData,
        trajectory: Trajectory,
        engine: ReconEngine,
        dat_name: str = "",
        seq_name: str = "",
    ) -> None:
        super().__init__()
        self.data = data
        self.trajectory = trajectory
        self.engine = engine

        title = "K-space Viewer / NUFFT Recon"
        if dat_name:
            title += f"  —  {dat_name}"
        self.setWindowTitle(title)

        self._selectors: Dict[str, QtWidgets.QSpinBox] = {}
        self._cancel_requested = False
        self._kspace_fitted = False
        self._recon_fitted = False
        self._build_ui()
        self.update_kspace_display()
        self._update_recon(force=False)

    # ------------------------------------------------------------------ build
    def _build_ui(self) -> None:
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        outer = QtWidgets.QVBoxLayout(central)

        # -- TOP: shared progress bar (hidden when idle) ----------------------
        prog_row = QtWidgets.QHBoxLayout()
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setTextVisible(True)
        prog_row.addWidget(self.progress_bar, stretch=1)
        self.progress_cancel = QtWidgets.QPushButton("Cancel")
        self.progress_cancel.setVisible(False)
        self.progress_cancel.clicked.connect(self._on_cancel_progress)
        prog_row.addWidget(self.progress_cancel)
        outer.addLayout(prog_row)

        # -- MAIN: three columns ----------------------------------------------
        root = QtWidgets.QHBoxLayout()
        outer.addLayout(root, stretch=1)

        # -- LEFT: ICE-dimension selectors ------------------------------------
        left = QtWidgets.QVBoxLayout()

        self.open_button = QtWidgets.QPushButton("Open .dat…")
        self.open_button.setToolTip("Load a different .dat file for a new analysis.")
        self.open_button.clicked.connect(self.on_open_dataset)
        left.addWidget(self.open_button)

        left.addWidget(QtWidgets.QLabel("<b>Dimensions</b>"))
        for dim in self.data.selector_dims:
            row = QtWidgets.QHBoxLayout()
            size = self.data.sizes[dim]
            row.addWidget(QtWidgets.QLabel(f"{dim} (0–{size - 1})"))
            spin = QtWidgets.QSpinBox()
            spin.setRange(0, size - 1)
            spin.valueChanged.connect(self.on_selector_change)
            row.addWidget(spin)
            self._selectors[dim] = spin
            left.addLayout(row)
        if not self.data.selector_dims:
            left.addWidget(QtWidgets.QLabel("<i>(no extra ICE dims)</i>"))
        if not self.trajectory.ramp_sampled:
            warn = QtWidgets.QLabel(
                "<span style='color:#c80'>⚠ No ramp readout found in .seq;\n"
                "using uniform kx (Cartesian).</span>"
            )
            warn.setWordWrap(True)
            left.addWidget(warn)
        left.addStretch(1)

        # Batch: reconstruct every selector/slice combination up front.
        self.recon_all_button = QtWidgets.QPushButton("Reconstruct All")
        self.recon_all_button.setToolTip(
            "Run the all-coil NUFFT for every dimension/slice combination "
            "and cache the results so browsing is instant."
        )
        self.recon_all_button.clicked.connect(self.on_reconstruct_all)
        left.addWidget(self.recon_all_button)

        # -- IMAGES: equal-size grid (col 0 = k-space, col 1 = recon) ---------
        # Both ImageViews share grid row 2 (the only stretch row) and two
        # equal-stretch columns, so the panels are always the same size.
        grid = QtWidgets.QGridLayout()
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        grid.setRowStretch(2, 1)

        # row 0 — titles (+ SOS toggle on the recon side)
        grid.addWidget(QtWidgets.QLabel("<b>k-space (magnitude)</b>"), 0, 0)
        recon_title = QtWidgets.QHBoxLayout()
        recon_title.addWidget(QtWidgets.QLabel("<b>Reconstruction</b>"))
        recon_title.addStretch(1)
        self.sos_checkbox = QtWidgets.QCheckBox("SOS (all coils)")
        self.sos_checkbox.stateChanged.connect(self.on_sos_toggle)
        recon_title.addWidget(self.sos_checkbox)
        grid.addLayout(recon_title, 0, 1)

        # row 1 — coil slider above the k-space panel
        coil_row = QtWidgets.QHBoxLayout()
        coil_row.addWidget(QtWidgets.QLabel("Coil"))
        self.coil_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.coil_slider.setRange(0, self.data.n_coils - 1)
        self.coil_slider.valueChanged.connect(self.on_coil_change)
        coil_row.addWidget(self.coil_slider)
        self.coil_label = QtWidgets.QLabel(f"0 / {self.data.n_coils - 1}")
        coil_row.addWidget(self.coil_label)
        grid.addLayout(coil_row, 1, 0)

        # row 2 — the two equal image panels
        self.kspace_view = _make_imageview()
        grid.addWidget(self.kspace_view, 2, 0)
        self.recon_view = _make_imageview()
        grid.addWidget(self.recon_view, 2, 1)

        # column 2 — recon toggles panel (extensible), beside the recon image
        grid.setColumnStretch(2, 0)  # keep it narrow
        toggles_group = QtWidgets.QGroupBox("Toggles")
        tg_layout = QtWidgets.QVBoxLayout(toggles_group)
        self._toggle_boxes: Dict[str, QtWidgets.QCheckBox] = {}
        for key, label, default in RECON_TOGGLES:
            cb = QtWidgets.QCheckBox(label)
            cb.setChecked(default)                       # set before connect
            cb.stateChanged.connect(self.on_toggle_change)
            tg_layout.addWidget(cb)
            self._toggle_boxes[key] = cb
        tg_layout.addStretch(1)
        grid.addWidget(toggles_group, 2, 2)

        # row 3 — slice slider (k-space side) | cached indicator (recon side)
        slice_row = QtWidgets.QHBoxLayout()
        slice_row.addWidget(QtWidgets.QLabel("Slice"))
        self.slice_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slice_slider.setRange(0, self.data.n_slices - 1)
        self.slice_slider.valueChanged.connect(self.on_slice_change)
        slice_row.addWidget(self.slice_slider)
        self.slice_label = QtWidgets.QLabel(f"0 / {self.data.n_slices - 1}")
        slice_row.addWidget(self.slice_label)
        grid.addLayout(slice_row, 3, 0)

        self.cache_label = QtWidgets.QLabel("")
        self.cache_label.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.cache_label, 3, 1)

        # row 4 — log toggle + Run button (k-space side)
        ctrl_row = QtWidgets.QHBoxLayout()
        self.log_checkbox = QtWidgets.QCheckBox("Log scale")
        self.log_checkbox.setChecked(True)
        self.log_checkbox.stateChanged.connect(self.update_kspace_display)
        ctrl_row.addWidget(self.log_checkbox)
        ctrl_row.addStretch(1)
        self.recon_button = QtWidgets.QPushButton("Run NUFFT Recon")
        self.recon_button.clicked.connect(self.on_run_recon)
        ctrl_row.addWidget(self.recon_button)
        grid.addLayout(ctrl_row, 4, 0)

        root.addLayout(left, stretch=0)
        root.addLayout(grid, stretch=1)

    # ------------------------------------------------------------------ state
    def _state(self) -> Dict[str, int]:
        state = {dim: spin.value() for dim, spin in self._selectors.items()}
        state["Sli"] = self.slice_slider.value()
        return state

    def _toggles(self) -> Dict[str, bool]:
        """Current recon toggle states, e.g. ``{"reverse_odd_lines": True}``."""
        return {k: cb.isChecked() for k, cb in self._toggle_boxes.items()}

    def _plane_key_for(self, state: Dict[str, int]) -> Tuple:
        return tuple((d, state[d]) for d in self.data.selector_dims) + (
            ("Sli", state["Sli"]),
        )

    def _plane_key(self) -> Tuple:
        return self._plane_key_for(self._state())

    def _all_states(self) -> List[Dict[str, int]]:
        """Every (selectors × slices) combination — one per distinct k-space."""
        ranges = [range(self.data.sizes[d]) for d in self.data.selector_dims]
        states: List[Dict[str, int]] = []
        for combo in itertools.product(*ranges, range(self.data.n_slices)):
            *sel_vals, sli = combo
            state = dict(zip(self.data.selector_dims, sel_vals))
            state["Sli"] = sli
            states.append(state)
        return states

    # --------------------------------------------------------------- handlers
    def on_selector_change(self, *_) -> None:
        self.update_kspace_display()
        self._update_recon(force=False)

    def on_coil_change(self, value: int) -> None:
        self.coil_label.setText(f"{value} / {self.data.n_coils - 1}")
        self.update_kspace_display()
        if not self.sos_checkbox.isChecked():
            self._update_recon(force=False)

    def on_slice_change(self, value: int) -> None:
        self.slice_label.setText(f"{value} / {self.data.n_slices - 1}")
        self.update_kspace_display()
        self._update_recon(force=False)

    def on_sos_toggle(self, *_) -> None:
        # Coil index is irrelevant to the SOS result -> disable the slider.
        self.coil_slider.setEnabled(not self.sos_checkbox.isChecked())
        # Re-derive the display from the already-computed coil images (instant).
        self._update_recon(force=False)

    def on_toggle_change(self, *_) -> None:
        # A recon-affecting toggle changed -> reconstruct the current plane so
        # the effect is shown immediately (cheap if that state is cached).
        self._update_recon(force=True)

    def on_run_recon(self) -> None:
        # Full NUFFT of all coils for the current slice / echo / rep / …
        self._update_recon(force=True)

    def on_reconstruct_all(self) -> None:
        """Reconstruct every k-space (all selector/slice combinations)."""
        states = self._all_states()
        sos = self.sos_checkbox.isChecked()
        coil = self.coil_slider.value()
        toggles = self._toggles()
        n = len(states)

        self.recon_all_button.setEnabled(False)
        self.recon_button.setEnabled(False)
        self._progress_begin(n, text=f"Reconstructing all — 0/{n}", cancelable=True)
        try:
            for i, state in enumerate(states):
                if self._cancel_requested:
                    break
                key = self._plane_key_for(state)
                if not self.engine.has_coils(key, toggles):
                    plane = self.data.get_plane(state)
                    self.engine.reconstruct(plane, key, sos=sos, coil_index=coil, toggles=toggles)
                self._progress_set(i + 1, f"Reconstructing all — {i + 1}/{n}")
        finally:
            self._progress_end()
            self.recon_all_button.setEnabled(True)
            self.recon_button.setEnabled(True)

        # Show the current selection (now cached, unless cancelled early).
        self._update_recon(force=False)

    # --------------------------------------------------------------- displays
    def update_kspace_display(self, *_) -> None:
        plane = self.data.get_plane(self._state())  # (C, Lin, Col)
        coil = min(self.coil_slider.value(), plane.shape[0] - 1)
        mag = np.abs(plane[coil])
        if self.log_checkbox.isChecked():
            mag = np.log(mag + 1e-10)
        # Fit the image to the panel on first show; keep the user's zoom after.
        self.kspace_view.setImage(mag, autoLevels=True, autoRange=not self._kspace_fitted)
        self._kspace_fitted = True

    def _update_recon(self, force: bool) -> None:
        """Update the recon panel.

        ``force=True`` (Run button) always runs the full all-coil NUFFT for the
        current plane. ``force=False`` (SOS toggle, coil/slice/selector change)
        only refreshes when the plane's coil images are already computed, so a
        display toggle is instant and never triggers an unexpected NUFFT.
        """
        plane_key = self._plane_key()
        toggles = self._toggles()
        if not force and not self.engine.has_coils(plane_key, toggles):
            return
        heavy = not self.engine.has_coils(plane_key, toggles)  # a real NUFFT will run
        plane = self.data.get_plane(self._state())
        if heavy:
            self._progress_begin(0, text="Reconstructing…")  # busy indicator
        try:
            img, cached = self.engine.reconstruct(
                plane,
                plane_key,
                sos=self.sos_checkbox.isChecked(),
                coil_index=self.coil_slider.value(),
                toggles=toggles,
            )
        finally:
            if heavy:
                self._progress_end()
        self._show_recon(img, cached)

    # ------------------------------------------------------------- progress
    def _on_cancel_progress(self) -> None:
        self._cancel_requested = True

    def _progress_begin(self, maximum: int, text: str = "", cancelable: bool = False) -> None:
        """Show the top progress bar. ``maximum=0`` gives a busy/indeterminate bar."""
        self._cancel_requested = False
        self.progress_bar.setRange(0, maximum)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat(text or "%p%")
        self.progress_bar.setVisible(True)
        self.progress_cancel.setVisible(cancelable)
        QtWidgets.QApplication.processEvents()

    def _progress_set(self, value: int, text: str | None = None) -> None:
        self.progress_bar.setValue(value)
        if text is not None:
            self.progress_bar.setFormat(text)
        QtWidgets.QApplication.processEvents()

    def _progress_end(self) -> None:
        self.progress_bar.setVisible(False)
        self.progress_cancel.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setFormat("%p%")

    def _show_recon(self, img: np.ndarray, cached: bool) -> None:
        self.recon_view.setImage(img, autoLevels=True, autoRange=not self._recon_fitted)
        self._recon_fitted = True
        self.cache_label.setText(
            "<span style='color:#0a0'>cached</span>" if cached else ""
        )

    # ------------------------------------------------------------- dataset
    def on_open_dataset(self) -> None:
        """Pick a new .dat (with auto-located .seq) and open it in a fresh window."""
        from . import main as _main  # lazy import to avoid a cycle

        win = _main.open_new_dataset(parent=self)
        if win is not None:
            win.show()
            self.close()
