# K-space Viewer & NUFFT Reconstruction Tool

A standalone PyQt5 desktop tool to inspect raw Siemens `.dat` k-space, browse
the stored ICE dimensions (echo, repetition, average, …), and reconstruct the
selected slice with NUFFT using the trajectory from a companion Pulseq `.seq`.

## Launch

* **Double-click** `Run_KSpace_Viewer.bat` — runs in the `pulses-gpu` env,
  browse for the `.dat`, the `.seq` is auto-located.
* **Run the file:** `python src/kspace_viewer/main.py`
* **Explicit paths:**

  ```bash
  conda activate pulses-gpu
  python -m src.kspace_viewer.main --dat <path.dat> [--seq <path.seq>]
  ```

`--dat` opens a browse dialog if omitted. The matching `.seq` is auto-located
from the `.dat` (searches its folder and a sibling `in/` folder, fuzzy-matched
on filename); a browse dialog opens only if that fails. See
[MANUAL.md](MANUAL.md) for full usage.

Dev/test dataset:

```bash
python -m src.kspace_viewer.main \
  --dat "src/hospital sessions/hdl/0416/out/meas_MID01073_FID422509_1_19_te80_nudir_lbl.dat" \
  --seq "src/hospital sessions/hdl/0416/in/1-19-te80-nudir-lbl.seq"
```

## Environment

Inside this repo it runs in the **`pulses-gpu`** conda env, which has the full
recon stack (`mapvbvd`, `mrinufft`, `finufft`, `pypulseq`, `torch`, `pyqtgraph`,
`PyQt5`). The NUFFT backend defaults to `finufft` (CPU); override with
`NUFFT_BACKEND=<backend>`.

## Sharing this tool / standalone install

This folder is **self-contained** — it does not import anything from the parent
repo (the few numeric helpers it needs are vendored in `_kspace_utils.py`). To
give it to someone else, copy the `kspace_viewer/` folder out on its own; that
folder is a complete project (zip it, or `git init` inside it).

**For a first-time setup, follow the step-by-step [INSTALL.md](INSTALL.md).**
The short version:

```bash
python -m venv venv
venv\Scripts\activate            # Windows  (source venv/bin/activate on macOS/Linux)
pip install -r requirements.txt
python -m kspace_viewer.main     # run from the folder that contains kspace_viewer/
```

On Windows you can instead double-click `Run_KSpace_Viewer.bat` (set
`KSPACE_PYTHON` to a specific `python.exe` if `python` on PATH isn't the right
env). Dependencies are listed in [requirements.txt](requirements.txt); note
`pymapvbvd` installs under that name but imports as `mapvbvd`, and
`finufft`/`mrinufft` need binary wheels or a compiler.

## UI

```
┌──────────┬──────────────────────────────┬──────────────────────────┐
│ Dims     │            Coil  [────●────]  │  Reconstruction   □ SOS  │
│  Eco ▢   │   ┌────────────────────────┐  │  ┌────────────────────┐  │
│  Rep ▢   │   │   k-space (magnitude)  │  │  │   recon image      │  │
│   …      │   │      grayscale         │  │  │   grayscale        │  │
│          │   └────────────────────────┘  │  └────────────────────┘  │
│          │           Slice [──●──────]   │        (cached)          │
│          │   ☑ Log scale  [Run NUFFT…]   │                          │
└──────────┴──────────────────────────────┴──────────────────────────┘
```

* **Left** — one spinbox per non-singleton ICE dim (excluding coil/slice).
  Singleton dims are auto-hidden.
* **Coil slider** (above k-space) and **slice slider** (below) are dedicated.
* **k-space** display is log-magnitude by default (toggle to linear).
* **Run NUFFT Recon** reconstructs the current selection. With **SOS** off it
  reconstructs the active coil; with SOS on it combines all coils
  (sum-of-squares, image domain). The coil slider is disabled under SOS.
* Both image panes support zoom/pan/levels (pyqtgraph).
* A **cached** indicator appears when a recon comes from the in-memory cache.

## Architecture

| File | Responsibility |
|------|----------------|
| `main.py` | CLI args, file dialogs, app wiring |
| `loader.py` | `.dat` loading (`mapvbvd`) + ICE-dim auto-detection |
| `trajectory.py` | `.seq` trajectory extraction (built from the ramp-sampled readout) |
| `recon.py` | NUFFT recon (`mrinufft`) + SOS + two-level cache |
| `ui.py` | PyQt5 + pyqtgraph widgets and layout |
| `_kspace_utils.py` | Vendored numeric helpers (see below) — keeps the package self-contained |

### Numeric helpers (`_kspace_utils.py`)

Three small helpers were originally reused from the parent repo's `src/utils`
and are now **vendored** here so the package has no dependency on the parent
project (and doesn't drag in MRzeroCore / matplotlib):

* `combine_sos` — image-domain sum-of-squares coil combination
  (from `utils_coils.py`);
* `flip_alternating_lines` — EPI odd-line reversal (from `utils_kspace.py`);
* `get_kx_and_ramp_info` — ramp-sampled readout kx (from `utils_kspace.py`).

Everything else — `.dat` loading via `mapvbvd`, the NUFFT via `mrinufft`,
ICE-dim auto-detection, caching, and the UI — is original to this module and
mirrors the established hospital-session recon scripts.

## Caching

In-memory only, cleared on close. Two levels:

* per **plane** (slice + ICE selectors) → complex coil images, so the NUFFT
  adjoint runs once and SOS / single-coil views are derived without recompute;
* per **display** keyed by the full tuple including the coil index, or `"SOS"`
  when SOS is on — so toggling SOS never collides with per-coil entries.

Disk caching is a possible future addition (not built).

## Flagged assumptions

* **Trajectory model.** For standard EPI the only non-uniform sampling is the
  ramp-sampled readout (kx); phase-encode (ky) is a uniform Cartesian grid. The
  trajectory is built once per session from `get_kx_and_ramp_info` as a
  **zigzag** (even lines +kx, odd lines −kx) and the recon flips the odd data
  lines to match — the combination verified to reconstruct cleanly against the
  reference pipeline (`…/0416/reconstructed/…_recon.py`). A monotonic kx with
  flipped data (or zigzag with un-flipped data) yields a half-FOV readout
  artefact.
* **Shared trajectory across selectors.** The readout trajectory does not change
  with echo/rep/diffusion-direction (only the pre-readout diffusion gradients
  differ, which don't shift readout k-space), so it is extracted once.
  `Trajectory.for_selector(state)` is the hook to specialise per selector if a
  future dataset (e.g. blip-up/blip-down) needs it.
* **Recon grid.** Nominal square matrix `(n_lines, n_lines)`; readout
  oversampling lives in the trajectory, not the image grid. Override via
  `ReconEngine(recon_shape=…)`.
* **SOS in image domain** after per-coil NUFFT (matches `combine_sos`).
* **Partial-Fourier echoes** (e.g. TE1). mapvbvd zero-fills the un-acquired PE
  lines in the `Lin` axis. The recon detects the acquired (non-zero) lines per
  plane and passes **only those** to the NUFFT, with a trajectory built from
  exactly their ky positions — so the trajectory matches the sampled k-space
  rather than treating the zero lines as if they were acquired. (No Hermitian
  fill is applied; the asymmetric ky coverage is reconstructed as-is.)
  For the dev dataset TE1 resolves to 72 acquired lines (indices 24–95),
  matching the reference pipeline.

## Validation

Verified against the dev dataset (`1_19_te80_nudir`): single-coil, SOS, and
partial-Fourier (TE1) reconstructions all match the reference pipeline's clean
circular-phantom output. The k-space display correctly shows the partial-Fourier
zero band on TE1 and the full centred DC on TE2/TE3.
