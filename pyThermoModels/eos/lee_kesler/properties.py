from dataclasses import dataclass

from ...configs import R_CONST
from .constants import LeeKeslerCoefficients
from .solver import _temperature_terms


@dataclass(frozen=True)
class LeeKeslerDeparture:
    helmholtz_over_rt: float
    enthalpy_departure: float
    entropy_departure: float
    gibbs_departure: float


def residual_helmholtz_over_rt(Tr: float, vr: float, coeffs: LeeKeslerCoefficients) -> float:
    # SECTION: Residual Helmholtz expression
    # ! Reduced temperature and volume must be positive.
    if Tr <= 0 or vr <= 0:
        raise ValueError("Tr and reduced volume must be positive.")
    B, C, D = _temperature_terms(Tr, coeffs)
    y = 1.0 / vr
    exp_term = __import__("math").exp(-coeffs.gamma * y**2)
    return float(
        B * y
        + 0.5 * C * y**2
        + 0.2 * D * y**5
        - coeffs.c4 / (2.0 * coeffs.gamma * Tr**3)
        * (coeffs.beta + 1.0 + coeffs.gamma * y**2)
        * exp_term
        + coeffs.c4 / (2.0 * coeffs.gamma * Tr**3) * (coeffs.beta + 1.0)
    )


def departure_properties(
    T: float,
    Tr: float,
    vr: float,
    Z: float,
    coeffs: LeeKeslerCoefficients,
) -> LeeKeslerDeparture:
    # SECTION: Base residual Helmholtz value
    alpha = residual_helmholtz_over_rt(Tr, vr, coeffs)
    tau = 1.0 / Tr
    # NOTE: Finite difference is scoped to temperature derivative only.
    step = max(1e-5, 1e-5 * tau)

    def alpha_tau(tau_value: float) -> float:
        # ? Hold reduced volume fixed while differentiating with respect to tau.
        return residual_helmholtz_over_rt(1.0 / tau_value, vr, coeffs)

    # SECTION: Departure identities
    d_alpha_d_tau = (alpha_tau(tau + step) - alpha_tau(tau - step)) / (2.0 * step)
    delta_alpha_delta = Z - 1.0
    h_over_rt = tau * d_alpha_d_tau + delta_alpha_delta
    s_over_r = tau * d_alpha_d_tau - alpha

    # SECTION: SI unit conversion
    enthalpy = R_CONST * T * h_over_rt
    entropy = R_CONST * s_over_r
    gibbs = enthalpy - T * entropy

    # SECTION: Result packaging
    return LeeKeslerDeparture(
        helmholtz_over_rt=float(alpha),
        enthalpy_departure=float(enthalpy),
        entropy_departure=float(entropy),
        gibbs_departure=float(gibbs),
    )
