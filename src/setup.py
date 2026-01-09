from __future__ import annotations

import os
import sys
import shutil
import subprocess
from pathlib import Path
from setuptools import setup
from setuptools.command.build_ext import build_ext


def _bool_env(name: str, default: bool) -> bool:
    val = os.environ.get(name)
    return val.lower() in {"1", "true", "yes", "on"} if val else default


ROOT = Path(__file__).parent.resolve()
FORTRAN_DIR = ROOT / "hvswdpy" / "_fortran"

SOURCES = [
    "modules.f90", "aux_procedures.f90", "dispersion_solver.f90",
    "root_solver.f90", "dispersion_equation.f90", "greens_rayleigh.f90",
    "greens_love.f90", "rayleigh_mode.f90", "love_mode.f90",
    "bodywave_integrals.f90", "hv_api.f90",
]

extra_compile_args = ["-ffree-form"]
extra_link_args = []

if _bool_env("USE_OPENMP", True):
    omp_flag = os.environ.get("OMP_FLAG", "-fopenmp")
    extra_compile_args.append(omp_flag)
    extra_link_args.append(omp_flag)
    if omp_lib := os.environ.get("OMP_LIB"):
        extra_link_args.append(omp_lib)

if _bool_env("DEBUG_F2PY", False):
    extra_compile_args.extend(["-O0", "-g", "-fcheck=all", "-fbacktrace"])


class F2PyBuildExt(build_ext):
    def run(self):
        # Override run() to ensure build_extensions() is called even when ext_modules=[]
        # The parent run() skips build_extensions() when there are no extensions
        self.build_extensions()
    
    def build_extensions(self):
        import numpy.f2py
        pyf_file = FORTRAN_DIR / "hvdfa.pyf"
        # Ensure output dirs are absolute (our f2py build runs with cwd=FORTRAN_DIR)
        output_dirs = [(ROOT / self.build_lib / "hvswdpy").resolve()]
        if self.inplace:
            output_dirs.append((ROOT / "hvswdpy").resolve())
        for output_dir in output_dirs:
            output_dir.mkdir(parents=True, exist_ok=True)
            self._build_with_f2py(pyf_file, output_dir)
    
    def _build_with_f2py(self, pyf_file, output_dir):
        # Use f2py compile mode (-c) so it actually builds the extension.
        # In the dascore conda toolchain we must force free-form + no line-length truncation,
        # otherwise .f90 sources may be compiled as fixed-form (72-col) and fail.
        f90flags = " ".join(extra_compile_args + ["-ffree-line-length-none"])
        f77flags = f90flags
        f2py_args = [
            sys.executable, "-m", "numpy.f2py",
            "-c",
            f"--f90flags={f90flags}",
            f"--f77flags={f77flags}",
            str(pyf_file),
            *[str(FORTRAN_DIR / s) for s in SOURCES],
            "-m", "HVSWDpy",
            "--quiet",
        ]
        
        env = os.environ.copy()
        if extra_link_args:
            existing_ldflags = env.get("LDFLAGS", "")
            env["LDFLAGS"] = f"{existing_ldflags} {' '.join(extra_link_args)}".strip()
        result = subprocess.run(
            f2py_args,
            check=False,
            capture_output=True,
            text=True,
            env=env,
            cwd=str(FORTRAN_DIR),
        )
        if result.returncode != 0:
            print(f"f2py failed with return code {result.returncode}")
            if result.stderr:
                print(result.stderr)
            if result.stdout:
                print(result.stdout)
            raise RuntimeError(f"f2py failed with return code {result.returncode}")

        so_files = list(FORTRAN_DIR.glob("HVSWDpy*.so"))
        if not so_files:
            so_files = [f for f in FORTRAN_DIR.glob("*.so") if "HVSWDpy" in f.name.lower()]
        for so_file in so_files:
            output_dir.mkdir(parents=True, exist_ok=True)
            dest = output_dir / so_file.name
            if dest.exists():
                dest.unlink()
            shutil.move(so_file, dest)
            print(f"Built extension: {dest}")

setup(
    name="hvswdpy",
    version="0.1.0",
    description="HV-SWD-DFA Fortran API and Python wrapper for surface wave dispersion and H/V forward modeling",
    author="Shihao Yuan",
    author_email="syuan@mines.edu",
    package_dir={"": "."},
    packages=["hvswdpy"],
    package_data={"hvswdpy": ["HVSWDpy*.so"]},
    ext_modules=[], 
    cmdclass={"build_ext": F2PyBuildExt},
    include_package_data=True,
)
