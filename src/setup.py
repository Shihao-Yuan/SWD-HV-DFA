from __future__ import annotations

import os
from pathlib import Path

try:
    from numpy.distutils.core import Extension, setup
    HAS_NUMPY_DISTUTILS = True
except ImportError:

    from setuptools import setup
    from setuptools.extension import Extension
    HAS_NUMPY_DISTUTILS = False


def _bool_env(name: str, default: bool) -> bool:
    val = os.environ.get(name)
    if val is None:
        return default
    return val.lower() in {"1", "true", "yes", "on"}


ROOT = Path(__file__).parent.resolve()

# Fortran sources
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
    
    omp_flag = os.environ.get("OMP_FLAG", "-fopenmp")
    extra_compile_args.append(omp_flag)
    extra_link_args.append(omp_flag)
    omp_lib = os.environ.get("OMP_LIB", "")
    if omp_lib:
        extra_link_args.append(omp_lib)

extra_compile_args.append("-ffree-form")

if _bool_env("DEBUG_F2PY", False):
    extra_compile_args.extend(["-O0", "-g", "-fcheck=all", "-fbacktrace"]) 

FORTRAN_DIR = ROOT / "hvswdpy" / "_fortran"

if HAS_NUMPY_DISTUTILS:
    ext = Extension(
        name="hvswdpy.HVSWDpy", 
        sources=[str(FORTRAN_DIR / "hvdfa.pyf")] + [str(FORTRAN_DIR / s) for s in sources],
        extra_f77_compile_args=extra_compile_args,
        extra_f90_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
    ext_modules = [ext]
else:
    ext_modules = []

setup(
    name="hvswdpy",
    version="0.1.0",
    description="HV-SWD-DFA Fortran API and Python wrapper for dispersion and H/V",
    author="Shihao Yuan",
    author_email="syuan@mines.edu",
    package_dir={"": "."},
    packages=["hvswdpy"],
    ext_modules=ext_modules if ext_modules else None,
    include_package_data=True,
)
