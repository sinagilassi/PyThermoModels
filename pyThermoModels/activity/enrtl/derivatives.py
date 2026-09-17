"""Analytic temperature derivatives for explicit ENRTL/NRTL tau correlations."""

from typing import Any, Optional

import numpy as np

from ...utils.utility import TauCorrelation


def tau_temperature_derivative(
    temperature: float,
    tau_correlation: TauCorrelation,
    *,
    dg_ij: Optional[Any] = None,
    a_ij: Optional[Any] = None,
    b_ij: Optional[Any] = None,
    c_ij: Optional[Any] = None,
    d_ij: Optional[Any] = None,
    d_tau_dT: Optional[Any] = None,
    gas_constant: float = 8.31446261815324,
) -> np.ndarray:
    """Return analytic ``d tau_ij / dT`` for a supported explicit correlation."""
    T = float(temperature)
    if not np.isfinite(T) or T <= 0.0:
        raise ValueError("temperature must be finite and greater than 0 K")
    if gas_constant <= 0.0 or not np.isfinite(gas_constant):
        raise ValueError("gas_constant must be finite and positive")

    # SECTION: Direct tau requires caller-supplied derivative information.
    if tau_correlation == "direct_tau":
        if d_tau_dT is None:
            raise NotImplementedError(
                "direct_tau has no temperature derivative without d_tau_dT"
            )
        return _finite_array(d_tau_dT, "d_tau_dT")
    if tau_correlation == "gibbs_energy":
        return -_required_array(dg_ij, "dg_ij") / (gas_constant * T**2)

    b = _required_array(b_ij, "b_ij")
    derivative = -b / T**2
    if tau_correlation == "inverse_temperature":
        return derivative
    if tau_correlation == "inverse_temperature_squared":
        return derivative - 2.0 * _required_array(c_ij, "c_ij") / T**3
    if tau_correlation == "inverse_log_temperature":
        return derivative + _required_array(c_ij, "c_ij") / T
    if tau_correlation == "extended_temperature":
        return derivative + _required_array(c_ij, "c_ij") / T + _required_array(d_ij, "d_ij")
    raise ValueError(f"Unsupported tau_correlation: {tau_correlation}")


def _required_array(value: Optional[Any], name: str) -> np.ndarray:
    if value is None:
        raise ValueError(f"{name} is required for this tau correlation")
    return _finite_array(value, name)


def _finite_array(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array
