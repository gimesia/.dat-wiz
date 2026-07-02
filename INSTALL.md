# Getting Started — Install & Run

Two ways to set up. Once it's running, see [MANUAL.md](MANUAL.md) for how to use
the interface.

- **Option A — Conda, one-click (Windows):** easiest if you have Miniconda/Anaconda.
- **Option B — pip + venv:** works everywhere.

You need the `kspace_viewer/` folder and either conda or **Python 3.10–3.13**.

---

# Option A — Conda, one-click (Windows)

Requires Miniconda or Anaconda installed
(https://docs.conda.io/en/latest/miniconda.html).

1. **Double-click `kspace_viewer/Setup_Conda.bat`** — once. It creates a conda
   env named `kspace-viewer` (or reuses it if it already exists), installs all
   dependencies, and remembers that env's Python for the launcher. (Takes a few
   minutes; leave the window open.) Safe to re-run — it also updates the deps.
2. **Double-click `kspace_viewer/Run_KSpace_Viewer.bat`** — every time you want
   to use the tool. It opens a browse dialog for your `.dat`.

That's it — no typing, no environment variables to set.

### Conda on macOS/Linux (no .bat there)

```bash
conda create -n kspace-viewer python=3.11 -y
conda activate kspace-viewer
pip install -r kspace_viewer/requirements.txt
python -m kspace_viewer.main          # run from the folder containing kspace_viewer/
```

---

# Option B — pip + venv

You need **Python 3.10–3.13** and the `kspace_viewer/` folder.

## 1. Open a terminal in the project folder

Go to the folder that **contains** `kspace_viewer/` (not inside it).

```
your-project/
  kspace_viewer/        <-- the tool
    main.py  requirements.txt  ...
```

- **Windows:** open the folder in File Explorer, type `cmd` in the address bar, press Enter.
- **macOS/Linux:** open Terminal and `cd` into the folder.

Check Python is available (should print 3.10 or newer):

```
python --version
```

If that fails on macOS/Linux, use `python3` everywhere below.

---

## 2. Create and activate a virtual environment

```bash
python -m venv venv
```

Activate it:

- **Windows (cmd):** `venv\Scripts\activate`
- **Windows (PowerShell):** `venv\Scripts\Activate.ps1`
- **macOS/Linux:** `source venv/bin/activate`

Your prompt now starts with `(venv)`. (Run the activate step again each new terminal.)

---

## 3. Install the dependencies

```bash
pip install -r kspace_viewer/requirements.txt
```

This pulls numpy, torch, pypulseq, pymapvbvd, mrinufft, finufft, PyQt5 and
pyqtgraph. First install can take a few minutes (torch is large).

---

## 4. Run it

```bash
python -m kspace_viewer.main
```

A file dialog opens — pick your Siemens **`.dat`** file. The matching **`.seq`**
is located automatically (it looks in the `.dat`'s folder and a sibling `in/`
folder); if it can't find one, a second dialog asks for the `.seq`.

The main window then opens. → Continue with [MANUAL.md](MANUAL.md).

### Other ways to launch

- **Explicit files (skip the dialogs):**
  ```bash
  python -m kspace_viewer.main --dat "path/to/scan.dat" --seq "path/to/scan.seq"
  ```
- **Windows double-click:** `kspace_viewer/Run_KSpace_Viewer.bat`. It uses, in
  order: a `KSPACE_PYTHON` you set, then `env_python.txt` (written by
  `Setup_Conda.bat`), then `python` on PATH. With a venv you can point it at the
  venv once:
  ```
  set KSPACE_PYTHON=C:\full\path\to\your-project\venv\Scripts\python.exe
  ```

---

## 5. Troubleshooting

| Problem | Fix |
|---------|-----|
| `python: command not found` | Install Python from python.org (tick "Add to PATH"); on macOS/Linux use `python3`. |
| `No module named kspace_viewer` | Run from the folder that **contains** `kspace_viewer/`, and make sure the venv is activated. |
| `No module named mapvbvd` | The pip package is `pymapvbvd` but imports as `mapvbvd` — `pip install pymapvbvd`. |
| pip fails building `finufft` / `mrinufft` | It's trying to compile from source. Upgrade pip (`python -m pip install -U pip`) so it grabs a prebuilt wheel; on an unusual platform install a C/C++ compiler. |
| Qt / "xcb" / display errors on Linux | Install system Qt libs, e.g. `sudo apt install libxcb-xinerama0 libgl1`. |
| Window opens but recon is slow the first time | Normal — the NUFFT operator is built on first use, then cached. Use **Reconstruct All** to pre-compute everything. |

The tool runs on CPU by default (`finufft` backend); no GPU is required.
