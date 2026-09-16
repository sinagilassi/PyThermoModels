from dataclasses import dataclass
from typing import Any

from ...configs import R_CONST
from .constants import REFERENCE_ACENTRIC_FACTOR, REFERENCE_FLUID, SIMPLE_FLUID
from .properties import departure_properties
from .solver import Branch, solve_reduced_volume


@dataclass(frozen=True)
class LeeKeslerResult:
    """
    Lee-Kesler pure fluid result container.

    Attributes
    ----------
    compressibility_factor : float
        The compressibility factor of the fluid.
    molar_volume : float
        The molar volume of the fluid.
    enthalpy_departure : float
        The enthalpy departure of the fluid.
    entropy_departure : float
        The entropy departure of the fluid.
    gibbs_departure : float
        The Gibbs free energy departure of the fluid.
    reduced_temperature : float
        The reduced temperature of the fluid.
    reduced_pressure : float
        The reduced pressure of the fluid.
    roots : tuple[float, ...]
        The roots of the reduced volume equation.
    metadata : dict[str, Any]
        Additional metadata related to the pure fluid calculation.
    """
    compressibility_factor: float
    molar_volume: float
    enthalpy_departure: float
    entropy_departure: float
    gibbs_departure: float
    reduced_temperature: float
    reduced_pressure: float
    roots: tuple[float, ...]
    metadata: dict[str, Any]


def _interpolate(simple: float, reference: float, acentric_factor: float) -> float:
    # SECTION: Three-parameter corresponding-states interpolation
    # NOTE: The simple fluid has omega=0 and the reference fluid uses n-octane omega.
    return simple + (acentric_factor / REFERENCE_ACENTRIC_FACTOR) * (reference - simple)


def lee_kesler_pure(
    P: float,
    T: float,
    Tc: float,
    Pc: float,
    acentric_factor: float,
    phase: Branch = "vapor",
) -> LeeKeslerResult:
    # SECTION: Input validation
    # ! Lee-Kesler reduced variables require positive dimensional inputs.
    if min(P, T, Tc, Pc) <= 0:
        raise ValueError("P, T, Tc, and Pc must be positive.")

    # SECTION: Reduced state
    Tr = T / Tc
    Pr = P / Pc

    # SECTION: Simple and reference fluid solves
    simple_root = solve_reduced_volume(Tr, Pr, SIMPLE_FLUID, branch=phase)
    reference_root = solve_reduced_volume(
        Tr, Pr, REFERENCE_FLUID, branch=phase)

    # SECTION: Compressibility and volume interpolation
    Z = _interpolate(
        simple_root.compressibility_factor,
        reference_root.compressibility_factor,
        acentric_factor,
    )
    molar_volume = Z * R_CONST * T / P

    # SECTION: Departure-property evaluations
    simple_dep = departure_properties(
        T=T,
        Tr=Tr,
        vr=simple_root.reduced_volume,
        Z=simple_root.compressibility_factor,
        coeffs=SIMPLE_FLUID,
    )
    reference_dep = departure_properties(
        T=T,
        Tr=Tr,
        vr=reference_root.reduced_volume,
        Z=reference_root.compressibility_factor,
        coeffs=REFERENCE_FLUID,
    )

    # SECTION: Result packaging
    # NOTE: Fugacity support is deliberately not claimed for Lee-Kesler in this step.
    return LeeKeslerResult(
        compressibility_factor=float(Z),
        molar_volume=float(molar_volume),
        enthalpy_departure=float(_interpolate(
            simple_dep.enthalpy_departure,
            reference_dep.enthalpy_departure,
            acentric_factor,
        )),
        entropy_departure=float(_interpolate(
            simple_dep.entropy_departure,
            reference_dep.entropy_departure,
            acentric_factor,
        )),
        gibbs_departure=float(_interpolate(
            simple_dep.gibbs_departure,
            reference_dep.gibbs_departure,
            acentric_factor,
        )),
        reduced_temperature=float(Tr),
        reduced_pressure=float(Pr),
        roots=simple_root.roots,
        metadata={
            "model": "Lee-Kesler",
            "source": "Lee and Kesler 1975, AIChE J. 21(3), 510-527, DOI 10.1002/aic.690210313",
            "reference_acentric_factor": REFERENCE_ACENTRIC_FACTOR,
            "phase_branch": phase,
            "simple_solver": simple_root,
            "reference_solver": reference_root,
            "fugacity_coefficient_supported": False,
        },
    )
