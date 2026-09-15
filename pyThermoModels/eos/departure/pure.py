# SECTION: Imports
from dataclasses import dataclass
from math import log

from ..alpha import attraction_parameter, attraction_parameter_derivatives
from ...configs import R_CONST


# SECTION: Result Containers
@dataclass(frozen=True)
class PureCubicResidualProperties:
    """
    Container for the residual thermodynamic properties of a pure cubic EOS fluid.
    All properties are evaluated for a single homogeneous phase at the specified
    temperature and pressure.

    Attributes
    ----------
    residual_enthalpy : float
        Residual enthalpy of the fluid.
    residual_entropy : float
        Residual entropy of the fluid.
    residual_gibbs : float
        Residual Gibbs free energy of the fluid.
    residual_internal_energy : float
        Residual internal energy of the fluid.
    residual_isobaric_heat_capacity : float
        Residual isobaric heat capacity of the fluid.
    residual_isochoric_heat_capacity : float
        Residual isochoric heat capacity of the fluid.
    compressibility_factor : float
        Compressibility factor of the fluid.
    molar_volume : float
        Molar volume of the fluid.
    phase : str
        Phase of the fluid (e.g., "liquid" or "vapor").
    eos_model : str
        EOS model used (e.g., "PR" or "SRK").
    convention : str, optional
        Description of the residual property convention, by default "real-fluid minus ideal-gas at same T, P, and composition".
    """
    residual_enthalpy: float
    residual_entropy: float
    residual_gibbs: float
    residual_internal_energy: float
    residual_isobaric_heat_capacity: float
    residual_isochoric_heat_capacity: float
    compressibility_factor: float
    molar_volume: float
    phase: str
    eos_model: str
    convention: str = "real-fluid minus ideal-gas at same T, P, and composition"


# SECTION: Internal Validation And EOS Helpers
def _validate_model(model: str) -> str:
    model_key = model.strip().upper()
    if model_key not in {"PR", "SRK"}:
        raise ValueError(
            "Pure residual properties are implemented for PR and SRK only")
    return model_key


def _log_term(Z: float, B: float, sigma: float, epsilon: float) -> float:
    # ! Invalid logarithm arguments usually mean the chosen root is not valid for this state.
    numerator = Z + sigma * B
    denominator = Z + epsilon * B
    if numerator <= 0 or denominator <= 0:
        raise ValueError(
            "Invalid logarithm argument for cubic residual property")
    return log(numerator / denominator)


# SECTION: Residual H/S/G/U/Cv Core
def _base_residual_properties(
    *,
    P: float,
    T: float,
    Z: float,
    a: float,
    b: float,
    da_dT: float,
    d2a_dT2: float,
    B: float,
    sigma: float,
    epsilon: float,
) -> tuple[float, float, float, float, float]:
    # ! These low-level equations assume one valid homogeneous EOS root.
    if P <= 0:
        raise ValueError("pressure must be greater than zero")
    if T <= 0:
        raise ValueError("temperature must be greater than zero")
    if b <= 0:
        raise ValueError("covolume parameter b must be greater than zero")
    if Z <= B:
        raise ValueError("compressibility factor must be greater than B")

    delta = sigma - epsilon
    if delta == 0:
        raise ValueError("sigma and epsilon must differ")

    # NOTE: Generic cubic attraction logarithm for the selected root.
    log_attraction = _log_term(Z, B, sigma, epsilon)
    attraction_factor = b * delta

    # NOTE: Residual convention is real-fluid minus ideal-gas at same T, P, composition.
    residual_enthalpy = (
        R_CONST * T * (Z - 1)
        + (T * da_dT - a) * log_attraction / attraction_factor
    )
    residual_entropy = (
        R_CONST * log(Z - B)
        + da_dT * log_attraction / attraction_factor
    )

    # NOTE: Lower-risk properties are evaluated from thermodynamic identities.
    residual_gibbs = residual_enthalpy - T * residual_entropy
    residual_internal_energy = residual_enthalpy - R_CONST * T * (Z - 1)
    residual_isochoric_heat_capacity = (
        T * d2a_dT2 * log_attraction / attraction_factor
    )

    return (
        residual_enthalpy,
        residual_entropy,
        residual_gibbs,
        residual_internal_energy,
        residual_isochoric_heat_capacity,
    )


# SECTION: Residual Cp Identity
def _residual_cp_from_cv(
    *,
    T: float,
    molar_volume: float,
    a: float,
    b: float,
    da_dT: float,
    residual_cv: float,
    sigma: float,
    epsilon: float,
) -> float:
    # ? Cp - Cv is evaluated from EOS pressure derivatives for the same root.
    d_value = (molar_volume + epsilon * b) * (molar_volume + sigma * b)
    d_dv = 2 * molar_volume + (sigma + epsilon) * b
    dP_dT_v = R_CONST / (molar_volume - b) - da_dT / d_value
    dP_dV_t = (
        -R_CONST * T / ((molar_volume - b) ** 2)
        + a * d_dv / (d_value ** 2)
    )

    # NOTE: Subtract ideal-gas R to convert real Cp-Cv into residual Cp from residual Cv.
    cp_minus_cv_real = -T * (dP_dT_v ** 2) / dP_dV_t
    return residual_cv + cp_minus_cv_real - R_CONST


# SECTION: Public Pure-Fluid Residual API
def residual_properties_cubic_pure(
    *,
    P: float,
    T: float,
    Z: float,
    Tc: float,
    Pc: float,
    alpha_acentric_input: float | None,
    eos_model: str,
    sigma: float,
    epsilon: float,
    psi: float,
    omega: float,
    phase: str,
) -> PureCubicResidualProperties:
    """
    Return pure-fluid residual/departure properties for PR/SRK.

    This is a root-explicit low-level helper: the caller must choose the EOS
    root before calling and pass that compressibility factor as ``Z``. The
    helper does not solve for roots, does not infer vapor versus liquid, and
    does not check phase equilibrium. The ``phase`` argument is carried only as
    metadata in the returned result so the selected root can be identified by
    the calling layer.

    H, S, G, U, Cp, and Cv are residual values with the convention real-fluid
    property minus ideal-gas property at the same T, P, and composition.
    Heat-capacity values are single homogeneous-root values.
    """
    model_key = _validate_model(eos_model)

    # SECTION: Pure-Component EOS Parameters
    ac_value = psi * (R_CONST ** 2) * (Tc ** 2) / Pc
    a_value = attraction_parameter(
        ac_value,
        T,
        Tc,
        alpha_acentric_input,
        model=model_key,
    )
    da_dT, d2a_dT2 = attraction_parameter_derivatives(
        ac_value,
        T,
        Tc,
        alpha_acentric_input,
        model=model_key,
    )
    b_value = omega * R_CONST * Tc / Pc
    B = b_value * P / (R_CONST * T)
    molar_volume = Z * R_CONST * T / P

    # SECTION: Residual Property Evaluation
    (
        residual_enthalpy,
        residual_entropy,
        residual_gibbs,
        residual_internal_energy,
        residual_isochoric_heat_capacity,
    ) = _base_residual_properties(
        P=P,
        T=T,
        Z=Z,
        a=a_value,
        b=b_value,
        da_dT=da_dT,
        d2a_dT2=d2a_dT2,
        B=B,
        sigma=sigma,
        epsilon=epsilon,
    )
    residual_isobaric_heat_capacity = _residual_cp_from_cv(
        T=T,
        molar_volume=molar_volume,
        a=a_value,
        b=b_value,
        da_dT=da_dT,
        residual_cv=residual_isochoric_heat_capacity,
        sigma=sigma,
        epsilon=epsilon,
    )

    # SECTION: Result Packaging
    return PureCubicResidualProperties(
        residual_enthalpy=residual_enthalpy,
        residual_entropy=residual_entropy,
        residual_gibbs=residual_gibbs,
        residual_internal_energy=residual_internal_energy,
        residual_isobaric_heat_capacity=residual_isobaric_heat_capacity,
        residual_isochoric_heat_capacity=residual_isochoric_heat_capacity,
        compressibility_factor=Z,
        molar_volume=molar_volume,
        phase=phase,
        eos_model=model_key,
    )
