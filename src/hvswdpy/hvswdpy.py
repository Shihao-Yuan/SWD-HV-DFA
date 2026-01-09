from __future__ import annotations

"""
Python wrapper for the HV-SWD-DFA Fortran extension.

Author (modifications and API): Shihao Yuan <syuan@mines.edu>
Original authorship for the Fortran algorithms is retained in the
corresponding source headers (see HV.f90 and related modules).
"""

from typing import NamedTuple, Optional, Tuple

import numpy as np

# Lazy import: extension may not be built yet after pip install -e
_api = None

def _get_api():
    """Get the Fortran API module, importing it lazily."""
    global _api
    if _api is None:
        try:
            from .HVSWDpy import hv_dfa_api as _api
        except ImportError as e:
            raise ImportError(
                "hvswdpy.HVSWDpy extension module not found. "
                "The extension needs to be built. Run:\n"
                "  python setup.py build_ext --inplace\n"
                f"Original error: {e}"
            ) from e
    return _api


class HVComponents(NamedTuple):
    img11_pihalf_total: np.ndarray
    img33_total: np.ndarray
    img11_pihalf_rayleigh: np.ndarray
    img11_pihalf_love: np.ndarray
    imvv: np.ndarray
    imhpsv: np.ndarray
    imhsh: np.ndarray
    status: Optional[int]


class DispersionCurves(NamedTuple):
    rayleigh_slowness: np.ndarray
    rayleigh_valid: np.ndarray
    love_slowness: np.ndarray
    love_valid: np.ndarray
    status: Optional[int]


def hv(
    frequencies_hz: np.ndarray,
    vp: np.ndarray,
    vs: np.ndarray,
    rho: np.ndarray,
    thickness: np.ndarray,
    *,
    n_rayleigh_modes: int = 3,
    n_love_modes: int = 3,
    precision_percent: float = 1.0,
    nks: int = 0,
    sh_damp: Optional[float] = None,
    psv_damp: Optional[float] = None,
) -> Tuple[np.ndarray, Optional[int]]:
    """Compute H/V spectral ratio for a layered model.

    Sizes are inferred from arrays; parameters mirror the Fortran API with
    simpler names and defaults.
    """
    api = _get_api()
    nf = int(np.asarray(frequencies_hz).size)
    nl = int(np.asarray(vp).size)
    return api.hv_compute_f2py(
        nf=nf,
        nl=nl,
        frequencies_hz=frequencies_hz,
        vp_in=vp,
        vs_in=vs,
        rho_in=rho,
        thickness_in=thickness,
        n_rayleigh_modes=n_rayleigh_modes,
        n_love_modes=n_love_modes,
        precision_percent=precision_percent,
        nks=nks,
        sh_damp=sh_damp,
        psv_damp=psv_damp,
    )


def hv_components(
    frequencies_hz: np.ndarray,
    vp: np.ndarray,
    vs: np.ndarray,
    rho: np.ndarray,
    thickness: np.ndarray,
    *,
    n_rayleigh_modes: int = 3,
    n_love_modes: int = 3,
    precision_percent: float = 1.0,
    nks: int = 0,
) -> HVComponents:
    """Compute individual H/V contributions (Rayleigh/Love, components)."""
    api = _get_api()
    nf = int(np.asarray(frequencies_hz).size)
    nl = int(np.asarray(vp).size)
    return HVComponents(
        *api.hv_compute_components_f2py(
            nf=nf,
            nl=nl,
            frequencies_hz=frequencies_hz,
            vp_in=vp,
            vs_in=vs,
            rho_in=rho,
            thickness_in=thickness,
            n_rayleigh_modes=n_rayleigh_modes,
            n_love_modes=n_love_modes,
            precision_percent=precision_percent,
            nks=nks,
        )
    )


def dispersion(
    frequencies_hz: np.ndarray,
    vp: np.ndarray,
    vs: np.ndarray,
    rho: np.ndarray,
    thickness: np.ndarray,
    *,
    n_rayleigh_modes: int = 3,
    n_love_modes: int = 3,
    precision_percent: float = 1.0,
) -> DispersionCurves:
    """Compute Rayleigh/Love slowness curves and validity masks."""
    api = _get_api()
    nf = int(np.asarray(frequencies_hz).size)
    nl = int(np.asarray(vp).size)
    rayleigh_slowness, rayleigh_valid_int, love_slowness, love_valid_int, status = api.hv_compute_dispersion_f2py(
        nf=nf,
        nl=nl,
        frequencies_hz=frequencies_hz,
        vp_in=vp,
        vs_in=vs,
        rho_in=rho,
        thickness_in=thickness,
        n_rayleigh_modes=n_rayleigh_modes,
        n_love_modes=n_love_modes,
        precision_percent=precision_percent,
    )
    # Convert integer arrays (0/1) to boolean arrays
    rayleigh_valid = rayleigh_valid_int.astype(bool)
    love_valid = love_valid_int.astype(bool)
    return DispersionCurves(
        rayleigh_slowness=rayleigh_slowness,
        rayleigh_valid=rayleigh_valid,
        love_slowness=love_slowness,
        love_valid=love_valid,
        status=status,
    )


