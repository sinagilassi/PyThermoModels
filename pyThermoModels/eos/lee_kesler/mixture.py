from dataclasses import dataclass
from typing import Any, Sequence, cast

import numpy as np

from .model import LeeKeslerResult, lee_kesler_pure, Branch


@dataclass(frozen=True)
class LeeKeslerMixtureResult:
    """
    Lee-Kesler mixture result container.

    Attributes
    ----------
    result : LeeKeslerResult
        The result of the Lee-Kesler pure fluid evaluation at the pseudo-critical state.
    pseudo_critical_temperature : float
        The pseudo-critical temperature of the mixture.
    pseudo_critical_pressure : float
        The pseudo-critical pressure of the mixture.
    pseudo_critical_volume : float
        The pseudo-critical volume of the mixture.
    mixture_acentric_factor : float
        The acentric factor of the mixture.
    composition : tuple[float, ...]
        The mole fraction composition of the mixture.
    metadata : dict[str, Any]
        Additional metadata related to the mixture calculation.
    """
    result: LeeKeslerResult
    pseudo_critical_temperature: float
    pseudo_critical_pressure: float
    pseudo_critical_volume: float
    mixture_acentric_factor: float
    composition: tuple[float, ...]
    metadata: dict[str, Any]


def _normalize_composition(x) -> np.ndarray:
    # SECTION: Composition validation
    arr = np.asarray(x, dtype=float)
    # ! Ploecker mixing rules require an ordered mole-fraction vector.
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("x must be a non-empty one-dimensional composition.")
    if np.any(arr < 0):
        raise ValueError("x must not contain negative values.")
    total = float(np.sum(arr))
    if total <= 0:
        raise ValueError("x must have a positive sum.")
    # NOTE: Normalize so amount-like feed vectors can be accepted safely.
    return arr / total


def _critical_volumes(Tc: np.ndarray, Pc: np.ndarray, Zc: np.ndarray | None) -> np.ndarray:
    # SECTION: Critical-volume helper
    from ...configs import R_CONST

    # ? If Zc is unavailable, use the documented LKP-style default and report it in metadata.
    zc = np.full_like(Tc, 0.2905) if Zc is None else np.asarray(
        Zc, dtype=float)
    # ! Zc order must match the component order.
    if zc.shape != Tc.shape:
        raise ValueError("Zc must match Tc/Pc shape.")
    return zc * R_CONST * Tc / Pc


def lee_kesler_ploecker_mixture(
    P: float,
    T: float,
    components: Sequence[str],
    x,
    Tc,
    Pc,
    acentric_factors,
    Zc=None,
    k_ij=None,
    phase: str = "vapor",
) -> LeeKeslerMixtureResult:
    # SECTION: Ordered component inputs
    comp = tuple(components)
    mole = _normalize_composition(x)
    n = len(comp)
    if mole.size != n:
        raise ValueError("components and x must have the same length.")

    # SECTION: Pure-component property arrays
    Tc_arr = np.asarray(Tc, dtype=float)
    Pc_arr = np.asarray(Pc, dtype=float)
    omega_arr = np.asarray(acentric_factors, dtype=float)
    # ! All property arrays must be aligned with the component order.
    if Tc_arr.shape != (n,) or Pc_arr.shape != (n,) or omega_arr.shape != (n,):
        raise ValueError(
            "Tc, Pc, and acentric_factors must match component count.")

    # SECTION: Binary interaction matrix
    # ? Missing k_ij is allowed but made visible in result metadata.
    if k_ij is None:
        kij = np.zeros((n, n), dtype=float)
        default_kij = True
    else:
        kij = np.asarray(k_ij, dtype=float)
        if kij.shape != (n, n):
            raise ValueError(f"k_ij must have shape {(n, n)}.")
        default_kij = False

    # SECTION: Ploecker pseudo-critical reducing rules
    Vc = _critical_volumes(
        Tc_arr, Pc_arr, None if Zc is None else np.asarray(Zc, dtype=float))
    vc_mix = 0.0
    tc_vc_mix = 0.0
    for i in range(n):
        for j in range(n):
            vc_ij = 0.125 * (Vc[i] ** (1.0 / 3.0) + Vc[j] ** (1.0 / 3.0)) ** 3
            tc_ij = (Tc_arr[i] * Tc_arr[j]) ** 0.5 * (1.0 - kij[i, j])
            weight = mole[i] * mole[j]
            vc_mix += weight * vc_ij
            tc_vc_mix += weight * tc_ij * vc_ij

    # ! Reducing volume must stay physically positive.
    if vc_mix <= 0:
        raise ValueError("Pseudo-critical volume must be positive.")

    # SECTION: Pseudo-critical state
    Tc_mix = tc_vc_mix / vc_mix
    omega_mix = float(np.dot(mole, omega_arr))
    zc_mix = float(np.dot(mole, np.full(n, 0.2905)
                   if Zc is None else np.asarray(Zc, dtype=float)))
    from ...configs import R_CONST

    Pc_mix = zc_mix * R_CONST * Tc_mix / vc_mix

    # SECTION: Corresponding pure-fluid evaluation
    # NOTE: The mixture reduces to a pure Lee-Kesler call at the pseudo-critical state.
    pure_like = lee_kesler_pure(
        P=P,
        T=T,
        Tc=float(Tc_mix),
        Pc=float(Pc_mix),
        acentric_factor=omega_mix,
        phase=cast(Branch, phase),
    )
    # SECTION: Result packaging
    return LeeKeslerMixtureResult(
        result=pure_like,
        pseudo_critical_temperature=float(Tc_mix),
        pseudo_critical_pressure=float(Pc_mix),
        pseudo_critical_volume=float(vc_mix),
        mixture_acentric_factor=omega_mix,
        composition=tuple(float(v) for v in mole),
        metadata={
            "model": "Lee-Kesler-Ploecker",
            "components": comp,
            "composition_normalized": not np.isclose(np.sum(np.asarray(x, dtype=float)), 1.0),
            "k_ij_defaulted_to_zero": default_kij,
            "Zc_defaulted_to_0.2905": Zc is None,
        },
    )
