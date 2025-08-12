from __future__ import annotations

import os
import sys
from pathlib import Path

from numpy.distutils.core import Extension, setup


def _bool_env(name: str, default: bool) -> bool:
    val = os.environ.get(name)
    if val is None:
        return default
    return val.lower() in {"1", "true", "yes", "on"}


ROOT = Path(__file__).parent.resolve()

# Fortran sources (refactored API + core algorithms)
sources = [
    "modules.f90",
    "aux_procedures.f90",
    "dispersion_solver.f90",
    "root_solver.f90",
    "dispersion_equation.f90",
    "greens_rayleigh.f90",
    "greens_love.f90",
    "rayleigh_mode.f90",
    "love_mode.f90",
    "bodywave_integrals.f90",
    "hv_api.f90",
]

use_openmp = _bool_env("USE_OPENMP", True)

extra_compile_args = []
extra_link_args = []

if use_openmp:
    # Default to GNU OpenMP flags. You can override via OMP_FLAG/OMP_LIB env vars.
    omp_flag = os.environ.get("OMP_FLAG", "-fopenmp")
    extra_compile_args.append(omp_flag)
    extra_link_args.append(omp_flag)
    # Do not force a specific OpenMP runtime by default; gfortran typically
    # links the correct runtime automatically. Allow override via OMP_LIB.
    omp_lib = os.environ.get("OMP_LIB", "")
    if omp_lib:
        extra_link_args.append(omp_lib)

# Ensure free-form parsing for any file misdetected as fixed-form
extra_compile_args.append("-ffree-form")

# Optional debug flags for runtime checks (enable with DEBUG_F2PY=1)
if _bool_env("DEBUG_F2PY", False):
    extra_compile_args.extend(["-O0", "-g", "-fcheck=all", "-fbacktrace"]) 

FORTRAN_DIR = ROOT / "hvswdpy" / "_fortran"

ext = Extension(
    name="HVSWDpy",
    sources=[str(FORTRAN_DIR / "hvdfa.pyf")] + [str(FORTRAN_DIR / s) for s in sources],
    extra_f77_compile_args=extra_compile_args,
    extra_f90_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
)


setup(
    name="hvswdpy",
    version="0.1.0",
    description="HV-SWD-DFA Fortran API and Python wrapper for dispersion and H/V",
    author="Shihao Yuan",
    author_email="syuan@mines.edu",
    url="",
    package_dir={"": "."},
    packages=["hvswdpy"],
    ext_modules=[ext],
    include_package_data=True,
)


