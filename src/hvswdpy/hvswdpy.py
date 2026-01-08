from __future__ import annotations

"""
Python wrapper for the HV-SWD-DFA Fortran extension.

Author (modifications and API): Shihao Yuan <syuan@mines.edu>
Original authorship for the Fortran algorithms is retained in the
corresponding source headers (see HV.f90 and related modules).
"""

from typing import NamedTuple, Optional, Tuple

import numpy as np

from .HVSWDpy import hv_dfa_api as _api


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
    nf = int(np.asarray(frequencies_hz).size)
    nl = int(np.asarray(vp).size)
    return _api.hv_compute_f2py(
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
    nf = int(np.asarray(frequencies_hz).size)
    nl = int(np.asarray(vp).size)
    return HVComponents(
        *_api.hv_compute_components_f2py(
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
    nf = int(np.asarray(frequencies_hz).size)
    nl = int(np.asarray(vp).size)
    return DispersionCurves(
        *_api.hv_compute_dispersion_f2py(
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
    )


