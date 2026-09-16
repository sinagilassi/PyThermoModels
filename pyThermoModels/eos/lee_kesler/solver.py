from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import optimize

from .constants import LeeKeslerCoefficients


Branch = Literal["vapor", "liquid", "stable"]


@dataclass(frozen=True)
class LeeKeslerRoot:
    reduced_volume: float
    compressibility_factor: float
    roots: tuple[float, ...]
    branch: str
    converged: bool
    iterations: int


def _temperature_terms(Tr: float, coeffs: LeeKeslerCoefficients) -> tuple[float, float, float]:
    # SECTION: Temperature-dependent BWR coefficients
    B = coeffs.b1 - coeffs.b2 / Tr - coeffs.b3 / Tr**2 - coeffs.b4 / Tr**3
    C = coeffs.c1 - coeffs.c2 / Tr + coeffs.c3 / Tr**3
    D = coeffs.d1 + coeffs.d2 / Tr
    return B, C, D


def z_from_reduced_volume(Tr: float, vr: float, coeffs: LeeKeslerCoefficients) -> float:
    # SECTION: Modified BWR compressibility
    # ! Reduced temperature and volume must be positive.
    if Tr <= 0 or vr <= 0:
        raise ValueError("Tr and reduced volume must be positive.")
    B, C, D = _temperature_terms(Tr, coeffs)
    inv_v = 1.0 / vr
    exponential = np.exp(-coeffs.gamma * inv_v**2)
    return float(
        1.0
        + B * inv_v
        + C * inv_v**2
        + D * inv_v**5
        + coeffs.c4 / Tr**3 * inv_v**2 * (coeffs.beta + coeffs.gamma * inv_v**2) * exponential
    )


def _pressure_residual(vr: float, Tr: float, Pr: float, coeffs: LeeKeslerCoefficients) -> float:
    # SECTION: Implicit reduced-pressure residual
    return Pr * vr / Tr - z_from_reduced_volume(Tr, vr, coeffs)


def _unique_positive(values: list[float], rel_tol: float = 1e-8) -> tuple[float, ...]:
    # SECTION: Root cleanup
    roots: list[float] = []
    for value in sorted(values):
        # ! Only positive finite roots are physical reduced volumes.
        if value <= 0 or not np.isfinite(value):
            continue
        # NOTE: Collapse repeated roots discovered from adjacent brackets.
        if not roots or abs(value - roots[-1]) > rel_tol * max(1.0, abs(value)):
            roots.append(float(value))
    return tuple(roots)


def solve_reduced_volume(
    Tr: float,
    Pr: float,
    coeffs: LeeKeslerCoefficients,
    branch: Branch = "vapor",
    tolerance: float = 1e-10,
    max_iterations: int = 100,
) -> LeeKeslerRoot:
    # SECTION: Solver input validation
    if Tr <= 0 or Pr <= 0:
        raise ValueError("Tr and Pr must be positive.")
    # ? Branch policy selects the smallest or largest physical reduced-volume root.
    if branch not in ("vapor", "liquid", "stable"):
        raise ValueError("branch must be 'vapor', 'liquid', or 'stable'.")

    # SECTION: Bracket scan
    # NOTE: A logarithmic grid is robust across vapor-like and liquid-like roots.
    grid = np.logspace(-4, 5, 450)
    residuals = np.asarray([_pressure_residual(v, Tr, Pr, coeffs) for v in grid])
    roots: list[float] = []
    iterations = 0

    # SECTION: Root solve
    for left, right, f_left, f_right in zip(grid[:-1], grid[1:], residuals[:-1], residuals[1:]):
        if not np.isfinite(f_left) or not np.isfinite(f_right):
            continue
        if f_left == 0.0:
            roots.append(float(left))
            continue
        if f_left * f_right > 0.0:
            continue
        result = optimize.root_scalar(
            _pressure_residual,
            args=(Tr, Pr, coeffs),
            bracket=(float(left), float(right)),
            xtol=tolerance,
            rtol=tolerance,
            maxiter=max_iterations,
        )
        iterations += result.iterations
        if result.converged:
            roots.append(float(result.root))

    # SECTION: Branch selection
    unique_roots = _unique_positive(roots)
    if not unique_roots:
        raise RuntimeError("Lee-Kesler reduced-volume solver did not find a physical root.")

    if branch == "liquid":
        vr = unique_roots[0]
    else:
        vr = unique_roots[-1]
    Z = Pr * vr / Tr

    # SECTION: Result packaging
    return LeeKeslerRoot(
        reduced_volume=float(vr),
        compressibility_factor=float(Z),
        roots=unique_roots,
        branch=branch,
        converged=True,
        iterations=iterations,
    )
