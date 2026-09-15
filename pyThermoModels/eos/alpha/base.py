"""
Alpha functions for the cubic equation-of-state layer.

This module centralizes the temperature-dependent alpha term used to build the
pure-component attraction parameter ``a(T)`` for the cubic EOS models currently
supported by PyThermoModels. It deliberately preserves the legacy Phase-1
behavior of ``EOSModels.eos_parameters()``:

* ``vdW`` uses ``alpha = 1``.
* ``RK`` uses ``alpha = Tr**-0.5``.
* ``SRK`` and ``PR`` use the existing Soave-style squared expression.

The public helpers accept absolute temperature ``T`` and critical temperature
``Tc`` rather than reduced temperature ``Tr`` so callers have one consistent
interface for values and derivatives. The derivative helpers return analytic
``dalpha/dT`` and ``d2alpha/dT2`` values with respect to absolute temperature
in kelvin. ``attraction_parameter*`` helpers simply scale those alpha results
by the critical attraction parameter ``a_c`` for equations of the form
``a(T) = a_c * alpha(T)``.

Important Phase-1 compatibility note: for PR/SRK, the caller still passes the
EOS parameter-table value named ``omega`` as ``acentric_factor``. That matches
the pre-existing numerical behavior and intentionally does not switch to
component databank acentric factors in this step.
"""

# SECTION: Imports
from dataclasses import dataclass
from math import pow


# SECTION: Result Containers
@dataclass(frozen=True)
class AlphaResult:
    """Alpha value and its first two temperature derivatives."""

    alpha: float
    dalpha_dT: float
    d2alpha_dT2: float


# SECTION: Internal Validation And Model Selection
def _validate_temperatures(temperature: float, critical_temperature: float) -> None:
    # ! Alpha correlations are undefined for non-positive absolute temperatures.
    if temperature <= 0:
        raise ValueError("temperature must be greater than zero")
    if critical_temperature <= 0:
        raise ValueError("critical_temperature must be greater than zero")


def _normalize_model(model: str) -> str:
    # NOTE: Keep public spelling compatible with the existing eos model names.
    model_key = model.strip()
    if model_key.upper() == "VDW":
        return "vdW"
    if model_key.upper() in {"RK", "SRK", "PR"}:
        return model_key.upper()
    raise ValueError(f"Unsupported alpha model: {model}")


def _m_coefficient(model: str, acentric_factor: float | None) -> float:
    # NOTE: Phase 1 preserves the legacy PR/SRK alpha input convention.
    omega = 0.0 if acentric_factor is None else acentric_factor
    if model == "SRK":
        return 0.480 + 1.574 * omega - 0.176 * pow(omega, 2)
    if model == "PR":
        return 0.37464 + 1.54226 * omega - 0.26992 * pow(omega, 2)
    raise ValueError(f"Model {model} does not use a Soave alpha coefficient")


# SECTION: Alpha Values And Temperature Derivatives
def alpha_result(
    temperature: float,
    critical_temperature: float,
    acentric_factor: float | None = None,
    *,
    model: str,
) -> AlphaResult:
    """
    Return alpha and temperature derivatives for one cubic EOS alpha model.

    Parameters are in SI-compatible scalar form: ``temperature`` and
    ``critical_temperature`` are kelvin, and ``acentric_factor`` is the
    dimensionless alpha correlation input. For this Phase-1 implementation,
    PR/SRK callers pass the legacy EOS parameter-table ``omega`` value here.

    The returned derivatives are analytic derivatives with respect to absolute
    temperature, not reduced temperature.
    """
    _validate_temperatures(temperature, critical_temperature)
    model_key = _normalize_model(model)

    # ? vdW has no temperature-dependent attraction correction.
    if model_key == "vdW":
        return AlphaResult(alpha=1.0, dalpha_dT=0.0, d2alpha_dT2=0.0)

    # ? RK uses the current reduced-temperature expression exactly.
    if model_key == "RK":
        alpha = pow(temperature / critical_temperature, -0.5)
        dalpha = -0.5 * alpha / temperature
        d2alpha = 0.75 * alpha / pow(temperature, 2)
        return AlphaResult(alpha=alpha, dalpha_dT=dalpha, d2alpha_dT2=d2alpha)

    # NOTE: PR/SRK share the same squared Soave-style alpha shape.
    m_value = _m_coefficient(model_key, acentric_factor)
    sqrt_tr = pow(temperature / critical_temperature, 0.5)
    g_value = 1 + m_value * (1 - sqrt_tr)
    alpha = pow(g_value, 2)

    # NOTE: Derivatives are with respect to absolute temperature [K].
    dg_dT = -m_value / (2 * critical_temperature * sqrt_tr)
    d2g_dT2 = m_value / (4 * pow(critical_temperature, 2) * pow(sqrt_tr, 3))
    dalpha = 2 * g_value * dg_dT
    d2alpha = 2 * pow(dg_dT, 2) + 2 * g_value * d2g_dT2
    return AlphaResult(alpha=alpha, dalpha_dT=dalpha, d2alpha_dT2=d2alpha)


# SECTION: Convenience Wrappers
def alpha_value(
    temperature: float,
    critical_temperature: float,
    acentric_factor: float | None = None,
    *,
    model: str,
) -> float:
    """Return only the dimensionless alpha value for ``model``."""
    return alpha_result(
        temperature,
        critical_temperature,
        acentric_factor,
        model=model,
    ).alpha


def dalpha_dT(
    temperature: float,
    critical_temperature: float,
    acentric_factor: float | None = None,
    *,
    model: str,
) -> float:
    """Return ``dalpha/dT`` for ``model`` in 1/K."""
    return alpha_result(
        temperature,
        critical_temperature,
        acentric_factor,
        model=model,
    ).dalpha_dT


def d2alpha_dT2(
    temperature: float,
    critical_temperature: float,
    acentric_factor: float | None = None,
    *,
    model: str,
) -> float:
    """Return ``d2alpha/dT2`` for ``model`` in 1/K^2."""
    return alpha_result(
        temperature,
        critical_temperature,
        acentric_factor,
        model=model,
    ).d2alpha_dT2


# SECTION: Attraction-Parameter Scaling
def attraction_parameter(
    critical_attraction_parameter: float,
    temperature: float,
    critical_temperature: float,
    acentric_factor: float | None = None,
    *,
    model: str,
) -> float:
    """Return ``a(T) = a_c * alpha(T)`` for a pure cubic EOS component."""
    return critical_attraction_parameter * alpha_value(
        temperature,
        critical_temperature,
        acentric_factor,
        model=model,
    )


def attraction_parameter_derivatives(
    critical_attraction_parameter: float,
    temperature: float,
    critical_temperature: float,
    acentric_factor: float | None = None,
    *,
    model: str,
) -> tuple[float, float]:
    """Return ``da/dT`` and ``d2a/dT2`` for ``a(T) = a_c * alpha(T)``."""
    alpha = alpha_result(
        temperature,
        critical_temperature,
        acentric_factor,
        model=model,
    )
    return (
        critical_attraction_parameter * alpha.dalpha_dT,
        critical_attraction_parameter * alpha.d2alpha_dT2,
    )
