"""Component-centric multicomponent Pitzer v2 regression coverage."""

import pytest
from pythermodb_settings.models import Component, Temperature

from pyThermoModels.activity.pitzer import (
    Pitzer, PitzerMulticomponentBinaryParameters, PitzerParameterSet,
    PitzerPsiParameter, PitzerThetaParameter, calc_Z,
    calc_multicomponent_ionic_strength, calc_multicomponent_net_charge,
    make_psi_c_aa_key, make_psi_cc_a_key, make_theta_key,
)


def _components():
    return [
        Component(name="Sodium", formula="Na{+}", state="aq"),
        Component(name="Potassium", formula="K{+}", state="aq"),
        Component(name="Magnesium", formula="Mg{2+}", state="aq"),
        Component(name="Chloride", formula="Cl{-}", state="aq"),
        Component(name="Sulfate", formula="SO4{2-}", state="aq"),
    ]


def _parameters(components):
    cations, anions = components[:3], components[3:]
    binary = {(c.get_key("Formula"), a.get_key("Formula")): PitzerMulticomponentBinaryParameters(.05, .15, c_phi=.001) for c in cations for a in anions}
    theta = {make_theta_key(a, b): PitzerThetaParameter(0.0) for group in (cations, anions) for a, b in __import__('itertools').combinations(group, 2)}
    psi_cc_a = {make_psi_cc_a_key(c1, c2, a): PitzerPsiParameter(0.0) for c1, c2 in __import__('itertools').combinations(cations, 2) for a in anions}
    psi_c_aa = {make_psi_c_aa_key(c, anions[0], anions[1]): PitzerPsiParameter(0.0) for c in cations}
    return PitzerParameterSet(A_phi=.3915, b=1.2, binary=binary, theta=theta, psi_cc_a=psi_cc_a, psi_c_aa=psi_c_aa)


def test_component_composition_functions_and_temperature_normalization():
    components = _components()
    molalities = {"Na{+}-aq": .4, "K{+}-aq": .1, "Mg{2+}-aq": .05, "Cl{-}-aq": .4, "SO4{2-}-aq": .1}
    assert calc_multicomponent_ionic_strength(components, molalities) == pytest.approx(.75)
    assert calc_Z(components, molalities) == pytest.approx(1.2)
    assert calc_multicomponent_net_charge(components, molalities) == pytest.approx(0.0)
    model = Pitzer(components, formulation="multicomponent_pitzer_v2")
    parameters = _parameters(components)
    _, celsius = model.cal({"temperature": Temperature(value=25, unit="C"), "molalities": molalities, "parameters": parameters})
    _, kelvin = model.cal({"temperature": Temperature(value=298.15, unit="K"), "molalities": molalities, "parameters": parameters})
    assert celsius["temperature_K"] == pytest.approx(298.15)
    assert celsius["activity_coefficients"] == pytest.approx(kelvin["activity_coefficients"])


def test_v2_binary_reduces_to_existing_nacl_regression():
    components = [Component(name="Sodium", formula="Na{+}", state="aq"), Component(name="Chloride", formula="Cl{-}", state="aq")]
    parameters = PitzerParameterSet(A_phi=.3915, b=1.2, binary={("Na{+}", "Cl{-}"): PitzerMulticomponentBinaryParameters(.0765, .2664, c_phi=.00127)}, theta={}, psi_cc_a={}, psi_c_aa={})
    result, other = Pitzer(components, formulation="multicomponent_pitzer_v2").cal({"temperature": Temperature(value=298.15, unit="K"), "molalities": {"Na{+}-aq": 1, "Cl{-}-aq": 1}, "parameters": parameters})
    assert result["value"]["Na{+}-aq"] == pytest.approx(.6555080908595792, rel=1e-12)
    assert other["osmotic_coefficient"] == pytest.approx(.9358687739996882, rel=1e-12)


def test_v2_requires_explicit_parameters_and_rejects_neutral_component():
    na, cl = Component(name="Sodium", formula="Na{+}", state="aq"), Component(name="Chloride", formula="Cl{-}", state="aq")
    parameters = PitzerParameterSet(A_phi=.3915, b=1.2, binary={}, theta={}, psi_cc_a={}, psi_c_aa={})
    with pytest.raises(ValueError, match="Missing required binary"):
        Pitzer([na, cl], formulation="multicomponent_pitzer_v2").cal({"temperature": Temperature(value=298.15, unit="K"), "molalities": {"Na{+}-aq": 1, "Cl{-}-aq": 1}, "parameters": parameters})
    water = Component(name="Water", formula="H2O", state="aq")
    with pytest.raises(NotImplementedError, match="lambda/zeta"):
        Pitzer([na, cl, water], formulation="multicomponent_pitzer_v2")
