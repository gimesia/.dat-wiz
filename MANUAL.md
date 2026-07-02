# K-space Viewer — User Manual

A step-by-step guide to using the tool. For architecture and design notes see
[README.md](README.md).

---

## 1. Before you start

You need:

* a Siemens raw data file (`.dat`),
* the companion Pulseq sequence file (`.seq`) that was used to acquire it.

The `.seq` is only used to read the k-space trajectory — no reconstruction
parameters are taken from it beyond the readout geometry.

Activate the environment that has all dependencies:

```bash
conda activate pulses-gpu
```

---

## 2. Starting the tool

### Simplest — double-click

Double-click **`Run_KSpace_Viewer.bat`** (in this folder). It launches in the
`pulses-gpu` env, opens a **browse dialog for the `.dat`**, and auto-locates the
matching `.seq`. A console window stays open showing the startup log.

### Simple — run the file

```bash
python src/kspace_viewer/main.py
```

Browse for the `.dat`; the `.seq` is found automatically (see below).

### The `.seq` is auto-located

After you pick a `.dat`, the tool looks for the matching `.seq` in the `.dat`'s
own folder and in a sibling `in/` folder, matching on the filename. If it can't
find one, a browse dialog opens for the `.seq`. You can always override with
`--seq`.

### With explicit paths

```bash
python -m src.kspace_viewer.main --dat "path/to/scan.dat" [--seq "path/to/scan.seq"]
```

### Example (dev dataset)

```bash
python -m src.kspace_viewer.main \
  --dat "src/hospital sessions/hdl/0416/out/meas_MID01073_FID422509_1_19_te80_nudir_lbl.dat" \
  --seq "src/hospital sessions/hdl/0416/in/1-19-te80-nudir-lbl.seq"
```

On startup the terminal prints what was detected, e.g.:

```
Loading .dat: ...meas_MID01073...nudir_lbl.dat
  dims=['Sli','Col','Cha','Lin','Eco','Rep']  coils=58  lines=96  samples=160  slices=1
  ICE selectors: ['Eco', 'Rep']
Extracting trajectory from .seq: ...1-19-te80-nudir-lbl.seq
  ramp_sampled=True  recon grid will be 96x96
```

If you see `ramp_sampled=False`, the tool could not find a ramp readout in the
`.seq` and fell back to a uniform (Cartesian) kx — a warning also shows in the
window. Reconstructions will still run but may not match a ramp-sampled scan.

---

## 3. The window at a glance

```
┌──────── progress bar (only while working) ────────────────────────────┐
├──────────┬──────────────────────────────┬──────────────────────────┤
│[Open .dat│  k-space (magnitude)         │  Reconstruction   ☑ SOS  │
│ Dims     │  Coil  [────●────]  0/17     │                          │
│  Eco ▢   │   ┌────────────────────────┐ │  ┌────────────────────┐  │
│  Rep ▢   │   │   k-space  (equal      │ │  │   recon  (equal    │  │
│          │   │    panel size)         │ │  │    panel size)     │  │
│          │   └────────────────────────┘ │  └────────────────────┘  │
│[Recon All│  Slice [──●──────]  0/2      │        cached            │
│          │  ☑ Log scale  [Run NUFFT…]   │                          │
└──────────┴──────────────────────────────┴──────────────────────────┘
```

Both image panels are always the same size. The progress bar appears across the
top only while a reconstruction runs.

| Area | What it does |
|------|--------------|
| **Open .dat…** (top-left) | Load a different `.dat` for a new analysis — opens in a fresh window, its `.seq` auto-located. |
| **Dimensions** (left) | One selector per non-singleton ICE dimension (echo, repetition, …). Singleton and coil/slice dims are not shown here. |
| **Coil slider** (top) | Picks the receive channel shown in the k-space panel and used for single-coil recon. |
| **k-space panel** (center-left) | Magnitude of the selected k-space plane, grayscale. |
| **Slice slider** (bottom) | Picks the slice. |
| **Log scale** | Toggles log vs linear magnitude for the k-space panel. |
| **Run NUFFT Recon** | Reconstructs the current selection into the right panel. |
| **SOS** (top-right) | Combine all coils (sum-of-squares) instead of a single coil. |
| **Reconstruction panel** (right) | The reconstructed image, grayscale. |
| **Toggles** (far right) | Recon options that change the reconstruction. **Reverse odd lines** (on by default) flips odd EPI lines to remove the N/2 ghost; turn it off to see the raw (ghosted) image. Toggling reconstructs immediately. |
| **cached** | Lights up when the shown recon came from the in-memory cache. |

---

## 4. Browsing k-space

The k-space panel updates **live** — no recon needed — whenever you change:

* any left-column **dimension** selector,
* the **coil** slider,
* the **slice** slider,
* the **Log scale** toggle.

Tips:

* **Log scale** (default on) makes low-signal structure visible. Turn it off to
  judge relative magnitude.
* On a **partial-Fourier** echo you will see a black band of un-acquired lines
  at one edge of k-space — that is expected.
* **Zoom / pan:** scroll to zoom, drag to pan. Right-click a panel for
  autoscale and export options (pyqtgraph). Double-clicking the little
  histogram resets levels.

---

## 5. Reconstructing an image

1. Set the **dimensions**, **slice**, and (for single-coil) the **coil** to what
   you want.
2. Choose the mode with the **SOS** checkbox:
   * **SOS off** → reconstructs the single coil shown by the coil slider.
   * **SOS on** → reconstructs every coil and combines them (sum-of-squares).
     The coil slider is greyed out because it no longer affects the result.
3. Click **Run NUFFT Recon**. The image appears in the right panel.

Partial-Fourier echoes are handled automatically: only the acquired lines are
used, matched to the trajectory (no manual setting needed).

### Reconstruct everything at once

The **Reconstruct All** button at the bottom of the left column runs the
all-coil NUFFT for **every** dimension/slice combination and caches the results.
A progress bar across the **top of the window** shows how many k-spaces remain,
with a **Cancel** button to stop early (whatever finished stays cached). Once it
finishes, switching between any echo / rep / slice / coil — and toggling SOS —
is instant, because every plane is already reconstructed. Use it when you want
to browse the whole dataset without waiting on each selection.

> The same top progress bar also appears (as a moving *busy* indicator) whenever
> a single **Run NUFFT Recon** actually computes — so you always get feedback
> that work is happening. Instant cache hits show no bar.

---

## 6. The cache

Every reconstruction is kept in memory, keyed by the full selection:

* single-coil results are keyed by coil index,
* SOS results are keyed separately, so toggling SOS on/off never overwrites your
  per-coil images.

When you return to a selection you already reconstructed, the image is shown
**instantly** and the **cached** indicator appears — no recompute. Changing a
selector to something new leaves the previous recon in place until you press
**Run NUFFT Recon** again.

The cache is in-memory only and is cleared when you close the tool.

---

## 7. Typical workflows

**Compare echoes for one slice/coil**

1. Coil slider → your channel, SOS off.
2. Set `Eco = 0`, click **Run NUFFT Recon**.
3. Set `Eco = 1`, **Run** again; `Eco = 2`, **Run** again.
4. Flip back and forth between echoes — each shows instantly from cache.

**Single-coil vs combined**

1. **Run** with SOS off (single coil).
2. Tick **SOS**, **Run** again to see the combined image.
3. Untick SOS — your single-coil image returns from cache instantly.

**Step through diffusion directions / repetitions**

1. SOS on for a clean combined image.
2. Increment the `Rep` selector and **Run** for each direction.

---

## 8. Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `ModuleNotFoundError` (mapvbvd / mrinufft / pyqtgraph) | Wrong env — `conda activate pulses-gpu`. |
| Startup warning "No ramp readout found" | The `.seq` has no trapezoid readout block; kx falls back to uniform. Check you passed the matching `.seq`. |
| Recon looks wrapped/mirrored horizontally | Usually a `.dat`/`.seq` mismatch — confirm both come from the same acquisition. |
| Coil slider does nothing | You are in **SOS** mode; the coil index is irrelevant there. |
| Nothing happens on selector change (recon panel) | By design — the recon panel only updates on **Run NUFFT Recon** or on a cache hit. |
| Slow first recon per selection | Building the NUFFT operator; subsequent identical selections are cached. Set `NUFFT_BACKEND=finufft` (default) for CPU. |

---

## 9. Quick reference

| Action | How |
|--------|-----|
| Change echo/rep/etc. | Left-column spinboxes |
| Change coil | Coil slider (top) |
| Change slice | Slice slider (bottom) |
| Log ↔ linear k-space | **Log scale** checkbox |
| Single-coil ↔ combined | **SOS** checkbox |
| Reconstruct current | **Run NUFFT Recon** |
| Reconstruct all planes | **Reconstruct All** (bottom-left) |
| Open a different `.dat` | **Open .dat…** (top-left) |
| Zoom / pan | Scroll / drag on a panel |
| Reset levels | Right-click → auto, or the histogram |
