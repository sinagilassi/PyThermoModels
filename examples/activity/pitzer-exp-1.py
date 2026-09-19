"""Run the component-centric multicomponent Pitzer v2 aqueous-ion example.

The example uses a true-species molality state, not apparent salts.  Parameter
values below are illustrative only; production work must use a validated,
temperature-specific Pitzer parameter source.
"""

from itertools import combinations
from pprint import pprint

from pythermodb_settings.models import Component, Temperature

from pyThermoModels.activity.pitzer import (
    Pitzer,
    PitzerMulticomponentBinaryParameters,
    PitzerParameterSet,
    PitzerPsiParameter,
    PitzerThetaParameter,
    make_psi_c_aa_key,
    make_psi_cc_a_key,
    make_theta_key,
)


# SECTION: True aqueous species and their molality state
components = [
    Component(name="Sodium ion", formula="Na{+}", state="aq"),
    Component(name="Potassium ion", formula="K{+}", state="aq"),
    Component(name="Magnesium ion", formula="Mg{2+}", state="aq"),
    Component(name="Chloride", formula="Cl{-}", state="aq"),
    Component(name="Sulfate", formula="SO4{2-}", state="aq"),
]
molalities = {
    "Na{+}-aq": 0.40,
    "K{+}-aq": 0.10,
    "Mg{2+}-aq": 0.05,
    "Cl{-}-aq": 0.40,
    "SO4{2-}-aq": 0.10,
}
cations, anions = components[:3], components[3:]


# SECTION: Explicit Formula-keyed parameters
# NOTE: Values are examples, not a literature parameter set.
binary_parameters = {
    (cation.get_key("Formula"), anion.get_key("Formula")):
    PitzerMulticomponentBinaryParameters(
        beta0=0.05,
        beta1=0.15,
        c_phi=0.001,
    )
    for cation in cations
    for anion in anions
}

# ! Zero-valued interactions are supplied explicitly because v2 is strict.
theta_parameters = {
    make_theta_key(left, right): PitzerThetaParameter(theta=0.0)
    for group in (cations, anions)
    for left, right in combinations(group, 2)
}
psi_cc_a_parameters = {
    make_psi_cc_a_key(cation_1, cation_2, anion): PitzerPsiParameter(psi=0.0)
    for cation_1, cation_2 in combinations(cations, 2)
    for anion in anions
}
psi_c_aa_parameters = {
    make_psi_c_aa_key(cation, anions[0], anions[1]): PitzerPsiParameter(psi=0.0)
    for cation in cations
}
parameters = PitzerParameterSet(
    A_phi=0.3915,
    b=1.2,
    binary=binary_parameters,
    theta=theta_parameters,
    psi_cc_a=psi_cc_a_parameters,
    psi_c_aa=psi_c_aa_parameters,
)


# SECTION: Calculate ionic Pitzer properties
# ? Temperature may use any PyCUC-supported unit; the model converts it to K.
model = Pitzer(components, formulation="multicomponent_pitzer_v2")
result, details = model.cal(
    {
        "temperature": Temperature(value=25.0, unit="C"),
        "molalities": molalities,
        "parameters": parameters,
    }
)

# Normalized molalities, I, Z, gamma_i, phi, and water activity are in details.
pprint(result)
pprint(details)
