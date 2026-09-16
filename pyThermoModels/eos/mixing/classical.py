from typing import Optional
from math import sqrt

import numpy as np

from .base import MixingResult


def _as_composition(xi) -> np.ndarray:
    # SECTION: Composition validation
    x = np.asarray(xi, dtype=float)
    # ! Composition must be a flat vector for all mixing rules.
    if x.ndim != 1 or x.size == 0:
        raise ValueError("xi must be a non-empty one-dimensional composition.")
    # ! Negative mole fractions are physically invalid.
    if np.any(x < 0):
        raise ValueError("xi must not contain negative mole fractions.")
    total = float(np.sum(x))
    if total <= 0:
        raise ValueError("xi must have a positive sum.")
    # NOTE: Normalize defensively so callers can pass feed-like amounts.
    return x / total


def _as_kij(k_ij: Optional[np.ndarray | list], size: int) -> tuple[np.ndarray, bool]:
    # SECTION: Binary interaction parameter validation
    # ? Missing k_ij means the documented zero-interaction approximation.
    if k_ij is None:
        return np.zeros((size, size), dtype=float), True
    kij = np.asarray(k_ij, dtype=float)
    # ! k_ij must align exactly with component order.
    if kij.shape != (size, size):
        raise ValueError(f"k_ij must have shape {(size, size)}.")
    return kij, False


def classical_quadratic_mixing(
    xi,
    params_list: list[dict],
    k_ij: Optional[np.ndarray | list] = None,
) -> MixingResult:
    """
    Classical quadratic attraction and linear covolume mixing rule.

    This preserves the legacy EOSModels.eos_mixing_rule behavior while returning
    a named result object for newer mixing-rule implementations.
    """
    # SECTION: Inputs
    x = _as_composition(xi)
    n = len(params_list)
    if x.size != n:
        raise ValueError("xi length must match params_list length.")

    # SECTION: Pure-component arrays
    kij, used_default_kij = _as_kij(k_ij, n)
    ai = np.asarray([params["a"] for params in params_list], dtype=float)
    bi = np.asarray([params["b"] for params in params_list], dtype=float)
    Ai = np.asarray([params["A"] for params in params_list], dtype=float)
    Bi = np.asarray([params["B"] for params in params_list], dtype=float)

    # SECTION: Quadratic attraction matrices
    a_ij = np.zeros((n, n), dtype=float)
    A_ij = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            # NOTE: Preserve the legacy geometric-mean attraction rule.
            a_ij[i, j] = (1.0 - kij[i, j]) * sqrt(ai[i] * ai[j])
            A_ij[i, j] = (1.0 - kij[i, j]) * sqrt(Ai[i] * Ai[j])

    # SECTION: Mixture parameters
    a_mix = float(np.sum(x[:, None] * x[None, :] * a_ij))
    A_mix = float(np.sum(x[:, None] * x[None, :] * A_ij))
    b_mix = float(np.dot(x, bi))
    B_mix = float(np.dot(x, Bi))

    # SECTION: Result packaging
    return MixingResult(
        a_mix=a_mix,
        b_mix=b_mix,
        a_ij=a_ij,
        A_mix=A_mix,
        B_mix=B_mix,
        metadata={
            "rule": "classical_quadratic",
            "composition_normalized": not np.isclose(np.sum(np.asarray(xi, dtype=float)), 1.0),
            "k_ij_defaulted_to_zero": used_default_kij,
        },
    )
