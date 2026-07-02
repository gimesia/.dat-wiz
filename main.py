"""Entry point: CLI args, file dialogs, and app wiring.

Usage
-----
Simple run (browse for the .dat; the .seq is auto-located):

    python src/kspace_viewer/main.py

Or with explicit paths:

    python -m src.kspace_viewer.main --dat <path.dat> [--seq <path.seq>]

If ``--dat`` is omitted a browse dialog opens. The matching ``.seq`` is found
automatically from the ``.dat`` name/location; if it can't be found a browse
dialog opens for it too.
"""

from __future__ import annotations

import argparse
import difflib
import glob
import os
import sys

# Support running as a script (python main.py) in addition to -m.
if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from kspace_viewer.loader import load_dat
    from kspace_viewer.trajectory import extract_trajectory
    from kspace_viewer.recon import ReconEngine
    from kspace_viewer.ui import MainWindow
else:
    from .loader import load_dat
    from .trajectory import extract_trajectory
    from .recon import ReconEngine
    from .ui import MainWindow

from PyQt5 import QtWidgets

# Keep references to open windows so Python doesn't garbage-collect them.
_OPEN_WINDOWS: list = []


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="K-space viewer / NUFFT recon tool")
    p.add_argument("--dat", help="Siemens raw data file (.dat)")
    p.add_argument("--seq", help="Companion Pulseq sequence file (.seq); "
                                 "auto-located from --dat if omitted")
    return p.parse_args()


def _pick_file(caption: str, file_filter: str, start_dir: str = "") -> str:
    path, _ = QtWidgets.QFileDialog.getOpenFileName(None, caption, start_dir, file_filter)
    if not path:
        raise SystemExit(f"No file selected for: {caption}")
    return path


def _find_seq_for_dat(dat_path: str) -> str | None:
    """Best-effort auto-location of the .seq that matches a .dat.

    Searches the .dat's own directory and a sibling ``in/`` directory (the
    hospital-session layout puts raw ``out/*.dat`` next to ``in/*.seq``), then
    fuzzy-matches on the filename tokens. Returns ``None`` if nothing plausible
    is found.
    """
    dat_dir = os.path.dirname(os.path.abspath(dat_path))
    stem = os.path.splitext(os.path.basename(dat_path))[0]

    search_dirs = [dat_dir, os.path.join(os.path.dirname(dat_dir), "in")]
    candidates: list[str] = []
    for d in search_dirs:
        if os.path.isdir(d):
            candidates.extend(glob.glob(os.path.join(d, "*.seq")))
    if not candidates:
        return None

    # Score by shared tokens; ".dat" names use "_" separators, ".seq" often "-".
    dat_tokens = set(stem.replace("-", "_").lower().split("_"))

    def score(seq_path: str) -> float:
        seq_stem = os.path.splitext(os.path.basename(seq_path))[0]
        seq_tokens = set(seq_stem.replace("-", "_").lower().split("_"))
        overlap = len(dat_tokens & seq_tokens)
        ratio = difflib.SequenceMatcher(None, stem, seq_stem).ratio()
        return overlap + ratio

    best = max(candidates, key=score)
    # Require at least some token overlap to avoid a spurious match.
    best_stem = os.path.splitext(os.path.basename(best))[0]
    if dat_tokens & set(best_stem.replace("-", "_").lower().split("_")):
        return best
    return None


def _build_window(dat_path: str, seq_path: str) -> MainWindow:
    """Load the data + trajectory and construct a MainWindow (not shown)."""
    print(f"Loading .dat: {dat_path}")
    data = load_dat(dat_path)
    print(
        f"  dims={data.dims}  coils={data.n_coils}  lines={data.n_lines}  "
        f"samples={data.n_samples}  slices={data.n_slices}"
    )
    print(f"  ICE selectors: {data.selector_dims or '(none)'}")

    print(f"Extracting trajectory from .seq: {seq_path}")
    trajectory = extract_trajectory(seq_path, data.n_lines, data.n_samples)
    print(
        f"  ramp_sampled={trajectory.ramp_sampled}  "
        f"recon grid will be {trajectory.n_lines}x{trajectory.n_lines}"
    )

    engine = ReconEngine(trajectory, calib=data.calib)
    window = MainWindow(
        data,
        trajectory,
        engine,
        dat_name=os.path.basename(dat_path),
        seq_name=os.path.basename(seq_path),
    )
    window.resize(1300, 700)
    _OPEN_WINDOWS.append(window)
    return window


def open_new_dataset(parent=None) -> "MainWindow | None":
    """Prompt for a new .dat (auto-locate its .seq) and build a window.

    Returns ``None`` if the user cancels — non-fatal, unlike startup, so the
    current session keeps running.
    """
    dat_path, _ = QtWidgets.QFileDialog.getOpenFileName(
        parent, "Select Siemens .dat file", "", "TWIX raw (*.dat)"
    )
    if not dat_path:
        return None

    seq_path = _find_seq_for_dat(dat_path)
    if seq_path:
        print(f"Auto-located .seq: {seq_path}")
    else:
        seq_path, _ = QtWidgets.QFileDialog.getOpenFileName(
            parent, "Select Pulseq .seq file (auto-locate failed)",
            os.path.dirname(os.path.abspath(dat_path)), "Pulseq (*.seq)",
        )
        if not seq_path:
            return None

    return _build_window(dat_path, seq_path)


def main() -> int:
    args = _parse_args()

    # QApplication must exist before any dialog / widget.
    app = QtWidgets.QApplication(sys.argv)

    dat_path = args.dat or _pick_file("Select Siemens .dat file", "TWIX raw (*.dat)")

    seq_path = args.seq
    if not seq_path:
        seq_path = _find_seq_for_dat(dat_path)
        if seq_path:
            print(f"Auto-located .seq: {seq_path}")
        else:
            seq_path = _pick_file(
                "Select Pulseq .seq file (auto-locate failed)",
                "Pulseq (*.seq)",
                start_dir=os.path.dirname(os.path.abspath(dat_path)),
            )

    window = _build_window(dat_path, seq_path)
    window.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
