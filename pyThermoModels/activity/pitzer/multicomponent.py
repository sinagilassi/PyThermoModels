"""Multicomponent aqueous-ion Pitzer v2 thermodynamic equations."""

from collections.abc import Mapping
from itertools import combinations
from math import exp, sqrt

from pythermodb_settings.models import Component

from .components import normalize_molalities, normalize_pitzer_components
from .functions import calc_debye_huckel_f, calc_g, calc_g_prime
from .interactions import (
    make_binary_interaction_key, make_psi_c_aa_key, make_psi_cc_a_key,
    make_theta_key, require_parameter,
)
from .mixing import calc_Phi, calc_Phi_phi, calc_Phi_prime
from .parameters import (
    PitzerMulticomponentBinaryParameters, PitzerParameterSet,
    PitzerPsiParameter, PitzerThetaParameter,
)


# SECTION: Composition functions
def calc_multicomponent_ionic_strength(components: list[Component], molalities: Mapping[str, float]) -> float:
    """Return ``I = 0.5 sum(m_i z_i^2)`` from Component metadata."""
    normalized = normalize_molalities(components, molalities)
    return 0.5 * sum(normalized[c.get_key("Formula-State")] * c.get_net_charge() ** 2 for c in components)


def calc_Z(components: list[Component], molalities: Mapping[str, float]) -> float:
    """Return Pitzer's ``Z = sum(m_i |z_i|)`` quantity."""
    normalized = normalize_molalities(components, molalities)
    return sum(normalized[c.get_key("Formula-State")] * abs(c.get_net_charge()) for c in components)


def calc_multicomponent_net_charge(components: list[Component], molalities: Mapping[str, float]) -> float:
    """Return molality-weighted net charge for a true-species state."""
    normalized = normalize_molalities(components, molalities)
    return sum(normalized[c.get_key("Formula-State")] * c.get_net_charge() for c in components)


def validate_multicomponent_electroneutrality(components: list[Component], molalities: Mapping[str, float], atol: float = 1.0e-12) -> bool:
    """Require a Pitzer v2 species state to be electrically neutral."""
    charge = calc_multicomponent_net_charge(components, molalities)
    if abs(charge) > atol:
        raise ValueError(f"Pitzer ionic state is not electroneutral: net molal charge = {charge:.6e}")
    return True


# SECTION: Binary interaction functions
def calc_B(ionic_strength: float, params: PitzerMulticomponentBinaryParameters) -> float:
    """Return the standard Pitzer binary activity interaction ``B``."""
    root_i = sqrt(ionic_strength)
    result = params.beta0 + params.beta1 * calc_g(params.alpha1 * root_i)
    if params.beta2 != 0.0:
        # ? alpha2 is validated by the parameter container when beta2 is active.
        result += params.beta2 * calc_g(float(params.alpha2) * root_i)
    return result


def calc_B_prime(ionic_strength: float, params: PitzerMulticomponentBinaryParameters) -> float:
    """Return derivative of B with respect to ionic strength."""
    if ionic_strength == 0.0:
        return 0.0
    root_i = sqrt(ionic_strength)
    result = params.beta1 * calc_g_prime(params.alpha1 * root_i) * params.alpha1 / (2.0 * root_i)
    if params.beta2 != 0.0:
        result += params.beta2 * calc_g_prime(float(params.alpha2) * root_i) * float(params.alpha2) / (2.0 * root_i)
    return result


def calc_B_phi(ionic_strength: float, params: PitzerMulticomponentBinaryParameters) -> float:
    """Return binary interaction used in the osmotic coefficient."""
    root_i = sqrt(ionic_strength)
    result = params.beta0 + params.beta1 * exp(-params.alpha1 * root_i)
    if params.beta2 != 0.0:
        result += params.beta2 * exp(-float(params.alpha2) * root_i)
    return result


def calc_C(c_phi: float, z_cation: int, z_anion: int) -> float:
    """Convert ``C_phi`` to the charge-scaled Pitzer ``C`` parameter."""
    return float(c_phi) / (2.0 * sqrt(abs(z_cation * z_anion)))


def _groups(components: list[Component]) -> tuple[list[Component], list[Component]]:
    return [c for c in components if c.is_cation()], [c for c in components if c.is_anion()]


def _m(component: Component, values: Mapping[str, float]) -> float:
    return values[component.get_key("Formula-State")]


def _binary(params: PitzerParameterSet, cation: Component, anion: Component, strict: bool) -> PitzerMulticomponentBinaryParameters | None:
    parameter = require_parameter(params.binary, make_binary_interaction_key(cation, anion), "binary", strict)
    if parameter is not None and not isinstance(parameter, PitzerMulticomponentBinaryParameters):
        raise TypeError("binary parameter mappings must contain PitzerMulticomponentBinaryParameters")
    return parameter


def _theta(params: PitzerParameterSet, left: Component, right: Component, strict: bool) -> PitzerThetaParameter | None:
    parameter = require_parameter(params.theta, make_theta_key(left, right), "theta", strict)
    if parameter is not None and not isinstance(parameter, PitzerThetaParameter):
        raise TypeError("theta parameter mappings must contain PitzerThetaParameter")
    return parameter


def _psi_cc_a(params: PitzerParameterSet, c1: Component, c2: Component, anion: Component, strict: bool) -> PitzerPsiParameter | None:
    parameter = require_parameter(params.psi_cc_a, make_psi_cc_a_key(c1, c2, anion), "psi_cc_a", strict)
    if parameter is not None and not isinstance(parameter, PitzerPsiParameter):
        raise TypeError("psi_cc_a parameter mappings must contain PitzerPsiParameter")
    return parameter


def _psi_c_aa(params: PitzerParameterSet, cation: Component, a1: Component, a2: Component, strict: bool) -> PitzerPsiParameter | None:
    parameter = require_parameter(params.psi_c_aa, make_psi_c_aa_key(cation, a1, a2), "psi_c_aa", strict)
    if parameter is not None and not isinstance(parameter, PitzerPsiParameter):
        raise TypeError("psi_c_aa parameter mappings must contain PitzerPsiParameter")
    return parameter


def calc_F(components: list[Component], values: Mapping[str, float], params: PitzerParameterSet, ionic_strength: float, strict_parameters: bool = True) -> float:
    """Return common Pitzer ``F`` term in individual-ion coefficients."""
    cations, anions = _groups(components)
    result = calc_debye_huckel_f(ionic_strength, params.A_phi, params.b)
    for cation in cations:
        for anion in anions:
            binary = _binary(params, cation, anion, strict_parameters)
            if binary is not None:
                result += _m(cation, values) * _m(anion, values) * calc_B_prime(ionic_strength, binary)
    for group in (cations, anions):
        for left, right in combinations(group, 2):
            theta = _theta(params, left, right, strict_parameters)
            if theta is not None:
                result += _m(left, values) * _m(right, values) * calc_Phi_prime(ionic_strength, left.get_net_charge(), right.get_net_charge(), params.A_phi)
    return result


# SECTION: Thermodynamic properties
def calc_ln_activity_coefficients(components: list[Component], molalities: Mapping[str, float], params: PitzerParameterSet, strict_parameters: bool = True) -> dict[str, float]:
    """Calculate conventional individual-ion Pitzer ``ln(gamma_i)`` values."""
    normalize_pitzer_components(components)
    values = normalize_molalities(components, molalities)
    validate_multicomponent_electroneutrality(components, values)
    ionic_strength = calc_multicomponent_ionic_strength(components, values)
    if ionic_strength == 0.0:
        return {c.get_key("Formula-State"): 0.0 for c in components}
    cations, anions = _groups(components)
    F = calc_F(components, values, params, ionic_strength, strict_parameters)
    result: dict[str, float] = {}
    for ion, opposites, same, is_cation in ((c, anions, cations, True) for c in cations):
        total = ion.get_net_charge() ** 2 * F
        for opposite in opposites:
            binary = _binary(params, ion, opposite, strict_parameters)
            if binary is not None:
                # NOTE: 3/2 converts C_phi to the individual-ion convention.
                total += _m(opposite, values) * (2.0 * calc_B(ionic_strength, binary) + 1.5 * calc_Z(components, values) * calc_C(binary.c_phi, ion.get_net_charge(), opposite.get_net_charge()))
        for partner in same:
            if partner is ion:
                continue
            theta = _theta(params, ion, partner, strict_parameters)
            if theta is not None:
                total += 2.0 * _m(partner, values) * calc_Phi(theta.theta, ionic_strength, ion.get_net_charge(), partner.get_net_charge(), params.A_phi)
            for opposite in opposites:
                psi = _psi_cc_a(params, ion, partner, opposite, strict_parameters)
                if psi is not None:
                    total += _m(partner, values) * _m(opposite, values) * psi.psi
        if is_cation:
            for a1, a2 in combinations(anions, 2):
                psi = _psi_c_aa(params, ion, a1, a2, strict_parameters)
                if psi is not None:
                    total += _m(a1, values) * _m(a2, values) * psi.psi
        result[ion.get_key("Formula-State")] = total
    for ion in anions:
        total = ion.get_net_charge() ** 2 * F
        for opposite in cations:
            binary = _binary(params, opposite, ion, strict_parameters)
            if binary is not None:
                total += _m(opposite, values) * (2.0 * calc_B(ionic_strength, binary) + 1.5 * calc_Z(components, values) * calc_C(binary.c_phi, opposite.get_net_charge(), ion.get_net_charge()))
        for partner in anions:
            if partner is ion:
                continue
            theta = _theta(params, ion, partner, strict_parameters)
            if theta is not None:
                total += 2.0 * _m(partner, values) * calc_Phi(theta.theta, ionic_strength, ion.get_net_charge(), partner.get_net_charge(), params.A_phi)
            for opposite in cations:
                psi = _psi_c_aa(params, opposite, ion, partner, strict_parameters)
                if psi is not None:
                    total += _m(partner, values) * _m(opposite, values) * psi.psi
        for c1, c2 in combinations(cations, 2):
            psi = _psi_cc_a(params, c1, c2, ion, strict_parameters)
            if psi is not None:
                total += _m(c1, values) * _m(c2, values) * psi.psi
        result[ion.get_key("Formula-State")] = total
    return result


def calc_activity_coefficients(*args, **kwargs) -> dict[str, float]:
    """Exponentiate :func:`calc_ln_activity_coefficients` componentwise."""
    return {key: exp(value) for key, value in calc_ln_activity_coefficients(*args, **kwargs).items()}


def calc_osmotic_coefficient_multicomponent(components: list[Component], molalities: Mapping[str, float], params: PitzerParameterSet, strict_parameters: bool = True) -> float:
    """Calculate the ionic multicomponent Pitzer osmotic coefficient."""
    normalize_pitzer_components(components)
    values = normalize_molalities(components, molalities)
    validate_multicomponent_electroneutrality(components, values)
    total_m = sum(values.values())
    if total_m == 0.0:
        return 1.0
    ionic_strength, Z = calc_multicomponent_ionic_strength(components, values), calc_Z(components, values)
    cations, anions = _groups(components)
    term = -params.A_phi * ionic_strength ** 1.5 / (1.0 + params.b * sqrt(ionic_strength))
    for cation in cations:
        for anion in anions:
            binary = _binary(params, cation, anion, strict_parameters)
            if binary is not None:
                term += _m(cation, values) * _m(anion, values) * (calc_B_phi(ionic_strength, binary) + Z * calc_C(binary.c_phi, cation.get_net_charge(), anion.get_net_charge()))
    for group, is_cation in ((cations, True), (anions, False)):
        for left, right in combinations(group, 2):
            theta = _theta(params, left, right, strict_parameters)
            if theta is not None:
                term += _m(left, values) * _m(right, values) * calc_Phi_phi(theta.theta, ionic_strength, left.get_net_charge(), right.get_net_charge(), params.A_phi)
            opposites = anions if is_cation else cations
            for opposite in opposites:
                psi = _psi_cc_a(params, left, right, opposite, strict_parameters) if is_cation else _psi_c_aa(params, opposite, left, right, strict_parameters)
                if psi is not None:
                    term += _m(left, values) * _m(right, values) * _m(opposite, values) * psi.psi
    return 1.0 + 2.0 * term / total_m


def calc_water_activity_multicomponent(components: list[Component], molalities: Mapping[str, float], params: PitzerParameterSet, strict_parameters: bool = True, water_molar_mass: float = 0.01801528) -> float:
    """Return water activity from multicomponent osmotic coefficient."""
    values = normalize_molalities(components, molalities)
    phi = calc_osmotic_coefficient_multicomponent(components, values, params, strict_parameters)
    return exp(-float(water_molar_mass) * sum(values.values()) * phi)


def calc_mean_ionic_activity_coefficient(ion_activity_coefficients: Mapping[str, float], stoichiometry: Mapping[str, int]) -> float:
    """Return a stoichiometrically weighted mean ionic activity coefficient."""
    if not stoichiometry:
        raise ValueError("stoichiometry must not be empty")
    total_nu, logarithm = 0, 0.0
    from math import log
    for component_id, nu in stoichiometry.items():
        if component_id not in ion_activity_coefficients or int(nu) <= 0:
            raise ValueError("stoichiometry must reference coefficients with positive integers")
        gamma = float(ion_activity_coefficients[component_id])
        if gamma <= 0.0:
            raise ValueError("activity coefficients must be positive")
        total_nu += int(nu)
        logarithm += int(nu) * log(gamma)
    return exp(logarithm / total_nu)
