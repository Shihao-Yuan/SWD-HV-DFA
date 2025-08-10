## HV-SWD-DFA: H/V Spectral Ratio and Dispersion Modeling Based on Diffuse Wavefield Assumption

Fortran implementation of H/V spectral ratio and surface-wave dispersion, with a Python wrapper. The Python API matches the original CLI results, and includes scripts to compare and visualize outputs.

- **Original authors**: HV‑INV project team (see headers in `HV.f90`)
- **Modifications and API**: Shihao Yuan (syuan@mines.edu)

### Repository layout
- `HV.f90`: Original CLI program (reference)
- `hv_api.f90`: Public API module used by Python
- Core Fortran sources: `modules.f90`, `aux_procedures.f90`, `dispersion_solver.f90`, `root_solver.f90`, `dispersion_equation.f90`, `greens_rayleigh.f90`, `greens_love.f90`, `rayleigh_mode.f90`, `love_mode.f90`, `bodywave_integrals.f90`
- Python build: `setup.py` (builds `HVSWDpy`), `hvdfa.pyf` (explicit f2py interface)
- Examples:
  - `examples/hv_quickstart.py`: basic H/V usage
  - `examples/compare_model_both.py`: compare HV and dispersion (CLI vs API) using `examples/model.txt`
- Outputs: `results/` (plots, comparisons)

### Requirements
- gfortran (OpenMP-capable recommended)
- Python 3.9+ with NumPy
- macOS/Linux (tested on macOS ARM with Conda)

### Build
- Original CLI:
  ```bash
  make hv_orig.exe
  ```
- Python wrapper (recommended):
  ```bash
  make python
  ```
- Both:
  ```bash
  make all
  ```

### Model format
- API arrays:
  - `vp` (m/s), `vs` (m/s), `rho` (kg/m³), `thickness` (m) for all layers except halfspace
  - At least 2 layers (including halfspace) so `thickness` has length `nlayers-1`
- CLI `model.txt` (used by examples):
  - First line: `N_LAYERS`
  - Next lines: `THICKNESS VP VS RHO`, with `THICKNESS=0` for the halfspace

### Original CLI usage
```bash
make hv_orig.exe
./hv_orig.exe -f examples/model.txt -fmin 0.1 -fmax 100 -nf 100 -logsam -nmr 3 -nml 3 -prec 1.0 -nks 0 -ph -hv > HV.dat
# Outputs: Rph.dat (Rayleigh slowness), Lph.dat (Love slowness), HV.dat (freq, hv)
```

### Python API
- H/V (no normalization; matches CLI default)
  ```python
  import numpy as np
  import hvswdpy as hv
  vp = np.array([300., 1500.]); vs = np.array([150., 800.])
  rho = np.array([1800., 2200.]); thickness = np.array([20.])
  f = np.logspace(-1, 2, 100)  # log-spaced to match CLI -logsam
  hv_curve, status = hv.hv(f, vp, vs, rho, thickness)
  ```
- Components
  ```python
  comps = hv.hv_components(f, vp, vs, rho, thickness)
  ```
- Dispersion (convert slowness to velocity via 1/slowness)
  ```python
  disp = hv.dispersion(f, vp, vs, rho, thickness)
  rayleigh_vel_mode1 = 1.0 / disp.rayleigh_slowness[disp.rayleigh_valid[:, 0] != 0, 0]
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

### Notes
- f2py auto-generated files (`hvdfa-f2pywrappers.f`, `hvdfa-f2pywrappers2.f90`) are build artifacts; don’t edit or commit.
- `make clean` removes executables, build artifacts, and Python extension outputs.

### License
See `LICENSE`.
