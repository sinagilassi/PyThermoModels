"""Independent Pitzer electrolyte activity-model package."""

from .binary import (
    calc_B_gamma,
    calc_B_phi,
    calc_C_gamma,
    calc_f_gamma,
    calc_f_phi,
    calc_gamma_mean_binary,
    calc_ln_gamma_mean_binary,
    calc_osmotic_coefficient_binary,
    calc_water_activity,
)
from .core import (
    build_binary_ion_molalities,
    calc_ionic_strength,
    calc_net_charge,
    validate_binary_electrolyte_v1,
    validate_electroneutrality,
)
from .model import Pitzer
from .parameters import (
    PitzerBinaryParameters, PitzerMulticomponentBinaryParameters,
    PitzerParameterSet, PitzerPsiParameter, PitzerThetaParameter,
)
from .components import PitzerSpeciesState, build_pitzer_species_state, normalize_molalities, normalize_pitzer_components
from .interactions import make_binary_interaction_key, make_psi_c_aa_key, make_psi_cc_a_key, make_theta_key
from .validation import validate_component_molality_mapping, validate_components, validate_molalities
from .multicomponent import (
    calc_B, calc_B_phi as calc_B_phi_multicomponent, calc_B_prime, calc_C,
    calc_F, calc_Z, calc_activity_coefficients, calc_ln_activity_coefficients,
    calc_mean_ionic_activity_coefficient, calc_multicomponent_ionic_strength,
    calc_multicomponent_net_charge, calc_osmotic_coefficient_multicomponent,
    calc_water_activity_multicomponent, validate_multicomponent_electroneutrality,
)

__all__ = [
    "Pitzer",
    "PitzerBinaryParameters",
    "PitzerMulticomponentBinaryParameters",
    "PitzerThetaParameter",
    "PitzerPsiParameter",
    "PitzerParameterSet",
    "PitzerSpeciesState",
    "normalize_pitzer_components",
    "normalize_molalities",
    "build_pitzer_species_state",
    "make_binary_interaction_key",
    "make_theta_key",
    "make_psi_cc_a_key",
    "make_psi_c_aa_key",
    "validate_components",
    "validate_molalities",
    "validate_component_molality_mapping",
    "build_binary_ion_molalities",
    "calc_ionic_strength",
    "calc_net_charge",
    "validate_binary_electrolyte_v1",
    "validate_electroneutrality",
    "calc_f_phi",
    "calc_B_phi",
    "calc_f_gamma",
    "calc_B_gamma",
    "calc_C_gamma",
    "calc_osmotic_coefficient_binary",
    "calc_ln_gamma_mean_binary",
    "calc_gamma_mean_binary",
    "calc_water_activity",
    "calc_multicomponent_ionic_strength",
    "calc_multicomponent_net_charge",
    "validate_multicomponent_electroneutrality",
    "calc_Z",
    "calc_B",
    "calc_B_phi_multicomponent",
    "calc_B_prime",
    "calc_C",
    "calc_F",
    "calc_ln_activity_coefficients",
    "calc_activity_coefficients",
    "calc_osmotic_coefficient_multicomponent",
    "calc_water_activity_multicomponent",
    "calc_mean_ionic_activity_coefficient",
]
