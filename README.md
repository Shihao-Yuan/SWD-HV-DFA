## HV-SWD-DFA: H/V Spectral Ratio and Surface Wave Dispersion Modeling (Diffuse Wavefield Assumption)

Fortran implementation of H/V spectral ratio and surface-wave dispersion, with a thin Python wrapper. The Python API mirrors the original CLI and offers simple functions.

- **Original authors**: HV‑INV project team (see headers in `HV.f90`)
- **Modifications and API**: Shihao Yuan (`syuan@mines.edu`)

### DISCLAIMER:
This is a development build. The code may contain errors or unstable functionality. Contributions and feedback are welcome.

### Repository layout
- `src/hvswdpy/`: Python package (wrapper)
- `src/hvswdpy/_fortran/`: Fortran sources and f2py interface (`hvdfa.pyf`)
- `src/cli/`: Fortran CLI sources (`HV.f90`, `cli_args*.f90`, optional drivers)
- `src/Makefile`: build targets (`hv_orig`, `python`, `python-dev`)
- `bin/`: built CLI executable (`bin/hv_orig`)
- `examples/`: scripts and notebooks
- `examples/results/`: output plots

### Requirements
- gfortran 
- Python 3.9+ with NumPy
- macOS/Linux (tested on macOS ARM with Conda)

### Build / Install
- From `src/` (recommended):
  - Original CLI:
    ```bash
    cd src
    make hv_orig        
    ```
  - Python extension and wrapper:
    ```bash
    cd src
    make python          
    ```
- From `src/` as an alternative:
  ```bash
  cd src
  python -m pip install -e .  
  ```

### Model format (API)
- API arrays:
  - `vp` (m/s), `vs` (m/s), `rho` (kg/m³), `thickness` (m) for all layers except halfspace
  - At least 2 layers (including halfspace) so `thickness` has length `nlayers-1`
- CLI `model.txt` (used by examples):
  - First line: `N_LAYERS`
  - Next lines: `THICKNESS VP VS RHO`, with `THICKNESS=0` for the halfspace

### Original CLI usage
```bash
# Build from src/
cd src && make hv_orig  

# Run from repo root
bin/hv_orig -f examples/model.txt -fmin 0.1 -fmax 100 -nf 100 -logsam -nmr 3 -nml 3 -prec 1.0 -nks 0 -ph -hv > examples/HV.dat
# Outputs in examples/: Rph.dat (Rayleigh slowness), Lph.dat (Love slowness), HV.dat (freq, hv)
```

### Python API 
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
  ```

Troubleshooting imports in notebooks:
- If running from a subfolder (e.g., `examples/`), make sure the project root is on `sys.path` or `PYTHONPATH` so `import hvswdpy` works:
  ```python
  import os, sys
  ROOT = os.path.abspath(os.path.join(os.getcwd(), ".."))
  if ROOT not in sys.path:
      sys.path.insert(0, ROOT)
  import hvswdpy
  ```

### Examples
- Quickstart (HV):
  ```bash
  python examples/hv_quickstart.py
  ```
- Compare HV and dispersion (CLI vs API) — generates plots into `examples/results/`:
  ```bash
  python examples/compare_API_CLI.py
  ```
  - `examples/results/compare_hv.png`
  - `examples/results/compare_rayleigh_dispersion.png` (phase velocity vs frequency)
  - `examples/results/compare_love_dispersion.png` (if Love modes requested)


### License
This project is licensed under the MIT License. See `LICENSE` for details.
