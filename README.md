## HV-SWD-DFA: H/V Spectral Ratio and Dispersion Modeling (Diffuse Wavefield Assumption)

Fortran implementation of H/V spectral ratio and surface-wave dispersion, with a thin Python wrapper. The
Python API mirrors the original CLI behavior (no normalization) and offers simple functions for HV, 
dispersion and components.

- **Original authors**: HV‑INV project team (see headers in `HV.f90`)
- **Modifications and API**: Shihao Yuan (`syuan@mines.edu`)

### Repository layout
- `HV.f90`: Original CLI program (reference)
- `hv_api.f90`: Public Fortran API module used by Python
- Core Fortran sources: `modules.f90`, `aux_procedures.f90`, `dispersion_solver.f90`, `root_solver.f90`, `dispersion_equation.f90`, `greens_rayleigh.f90`, `greens_love.f90`, `rayleigh_mode.f90`, `love_mode.f90`, `bodywave_integrals.f90`
- Python build: `setup.py` (builds compiled module `HVSWDpy`), `hvdfa.pyf` (explicit f2py interface)
- Python wrapper module: `hvswdpy.py` (imports `HVSWDpy.hv_dfa_api` and provides a friendly API)
- Examples:
  - `examples/hv_quickstart.py`: basic H/V usage
  - `examples/compare_model_both.py`: compare HV and dispersion (CLI vs API) using `examples/model.txt`
- Outputs: `results/` (plots, comparisons)

### Requirements
- gfortran (OpenMP-capable recommended)
- Python 3.9+ with NumPy
- macOS/Linux (tested on macOS ARM with Conda)

### Build / Install
- Original CLI only:
  ```bash
  make hv_orig
  ```
- Python extension and wrapper (recommended):
  ```bash
  # Option A: Makefile convenience
  make python

  # Option B: pip install from source
  python -m pip install .
  ```

Build options (env vars):
- `USE_OPENMP=0|1` (default 1) to disable/enable OpenMP in Fortran
- `DEBUG_F2PY=1` to add debug flags (`-O0 -g -fcheck=all -fbacktrace`)

### Model format (API)
- API arrays:
  - `vp` (m/s), `vs` (m/s), `rho` (kg/m³), `thickness` (m) for all layers except halfspace
  - At least 2 layers (including halfspace) so `thickness` has length `nlayers-1`
- CLI `model.txt` (used by examples):
  - First line: `N_LAYERS`
  - Next lines: `THICKNESS VP VS RHO`, with `THICKNESS=0` for the halfspace

### Original CLI usage
```bash
make hv_orig
./hv_orig -f examples/model.txt -fmin 0.1 -fmax 100 -nf 100 -logsam -nmr 3 -nml 3 -prec 1.0 -nks 0 -ph -hv > HV.dat
# Outputs: Rph.dat (Rayleigh slowness), Lph.dat (Love slowness), HV.dat (freq, hv)
```

### Python API (no normalization; matches CLI default)
  ```python
  import numpy as np
  import hvswdpy as hv
  vp = np.array([300., 1500.])
  vs = np.array([150., 800.])
  rho = np.array([1800., 2200.])
  thickness = np.array([20.])  # nlayers-1 values (halfspace excluded)
  f = np.logspace(-1, 2, 100)

  # H/V spectral ratio
  hv_curve, status = hv.hv(
      frequencies_hz=f,
      vp=vp, vs=vs, rho=rho, thickness=thickness,
      n_rayleigh_modes=1, n_love_modes=0, precision_percent=1.0,
  )
  ```
- Components
  ```python
  comps = hv.hv_components(f, vp, vs, rho, thickness)
  ```
- Dispersion (convert slowness to velocity via 1/slowness)
  ```python
  disp = hv.dispersion(
      frequencies_hz=f,
      vp=vp, vs=vs, rho=rho, thickness=thickness,
      n_rayleigh_modes=1, n_love_modes=0, precision_percent=1.0,
  )
  mask = (disp.rayleigh_valid[:, 0] != 0)
  rayleigh_vel_mode1 = 1.0 / disp.rayleigh_slowness[mask, 0]

Troubleshooting imports in notebooks:
- If running from a subfolder (e.g., `examples/`), make sure the project root is on `sys.path` or `PYTHONPATH` so `import hvswdpy` works:
  ```python
  import os, sys
  ROOT = os.path.abspath(os.path.join(os.getcwd(), ".."))
  if ROOT not in sys.path:
      sys.path.insert(0, ROOT)
  import hvswdpy
  ```
  ```

### Examples
- Quickstart (HV):
  ```bash
  python examples/hv_quickstart.py
  ```
- Compare HV and dispersion (CLI vs API) — generates plots into `results/`:
  ```bash
  python examples/compare_model_both.py
  ```
  - `results/compare_hv.png`
  - `results/compare_rayleigh_dispersion.png` (phase velocity vs frequency)
  - `results/compare_love_dispersion.png` (if Love modes requested)

- Bayesian inversion (BayesBay) example notebook:
  - `examples/BayesBay/BayesBay.ipynb` (and `BayesBay_fixed.ipynb`): demonstrates joint inversion of Rayleigh
    dispersion and HVSR using a Voronoi1D Vs model.
  - Notes:
    - Use `from bayesbay.likelihood import Target, LogLikelihood` (current BayesBay API)
    - Forward functions must pass `thickness[:-1]` (exclude halfspace)
    - To keep output quiet during long runs: set `verbose=False` in `inv.run(...)`
    - When running from the notebook folder, add the project root to `sys.path` (see snippet above)

### Notes
- f2py auto-generated files are build artifacts; don’t edit or commit.
- `make clean` removes executables, build artifacts, and Python extension outputs.

### License
See `LICENSE`.
