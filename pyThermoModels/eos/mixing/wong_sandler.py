from typing import Callable, Optional
from math import sqrt

import numpy as np

from ...configs import R_CONST
from .base import MixingResult
from .classical import _as_composition, _as_kij


PR_WONG_SANDLER_OMEGA = -0.62323


def wong_sandler_pr_nrtl_mixing(
    T: float,
    xi,
    params_list: list[dict],
    excess_gibbs_over_rt: Callable[[float, np.ndarray], float],
    k_ij: Optional[np.ndarray | list] = None,
) -> MixingResult:
    """
    Wong-Sandler mixing rule for Peng-Robinson coupled to a G^E/RT callback.

    First supported scope is PR + NRTL. The callback must return dimensionless
    G^E/(RT) at the supplied temperature and composition.
    """
    # SECTION: Input validation
    # ! Current implementation is intentionally limited to PR + NRTL-compatible G^E/RT.
    if T <= 0:
        raise ValueError("T must be positive.")

    x = _as_composition(xi)
    n = len(params_list)
    if x.size != n:
        raise ValueError("xi length must match params_list length.")
    for params in params_list:
        if params.get("eos-model") != "PR":
            raise ValueError("Wong-Sandler implementation currently supports PR only.")

    # SECTION: Pure-component arrays
    kij, used_default_kij = _as_kij(k_ij, n)
    ai = np.asarray([params["a"] for params in params_list], dtype=float)
    bi = np.asarray([params["b"] for params in params_list], dtype=float)
    Ai = np.asarray([params["A"] for params in params_list], dtype=float)
    Bi = np.asarray([params["B"] for params in params_list], dtype=float)
    # ! PR covolume parameters must be positive for the WS formulas.
    if np.any(bi <= 0):
        raise ValueError("All pure-component b parameters must be positive.")

    # SECTION: Second-virial cross term
    b_minus_a_rt_ij = np.zeros((n, n), dtype=float)
    a_ij = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            # NOTE: This form follows the PR Wong-Sandler first-scope implementation.
            a_ij[i, j] = (1.0 - kij[i, j]) * sqrt(ai[i] * ai[j])
            b_minus_a_rt_ij[i, j] = 0.5 * (bi[i] + bi[j]) - a_ij[i, j] / (R_CONST * T)

    # SECTION: Wong-Sandler mixture calculation
    q_value = float(np.sum(x[:, None] * x[None, :] * b_minus_a_rt_ij))
    # ? Activity model is injected so NRTL/UNIQUAC/etc. can be validated separately.
    ge_rt = float(excess_gibbs_over_rt(T, x))
    denominator = 1.0 - np.sum(x * ai / (bi * R_CONST * T)) - ge_rt / PR_WONG_SANDLER_OMEGA
    # ! Avoid silently returning singular or explosive parameters.
    if abs(denominator) < 1e-14:
        raise ZeroDivisionError("Wong-Sandler denominator is too close to zero.")

    b_mix = q_value / denominator
    a_mix = b_mix * (float(np.sum(x * ai / bi)) + (R_CONST * T * ge_rt) / PR_WONG_SANDLER_OMEGA)
    A_mix = float(a_mix / (R_CONST * T) ** 2)
    B_mix = float(np.dot(x, Bi))

    # SECTION: Result packaging
    return MixingResult(
        a_mix=float(a_mix),
        b_mix=float(b_mix),
        a_ij=a_ij,
        A_mix=A_mix,
        B_mix=B_mix,
        metadata={
            "rule": "wong_sandler",
            "eos_model": "PR",
            "activity_model": "NRTL-compatible",
            "ge_rt": ge_rt,
            "omega_constant": PR_WONG_SANDLER_OMEGA,
            "q": q_value,
            "k_ij_defaulted_to_zero": used_default_kij,
        },
    )
