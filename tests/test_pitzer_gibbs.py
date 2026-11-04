"""Thermodynamic consistency checks for the canonical Gibbs kernel."""
import numpy as np
import pytest
from autograd import grad, hessian
from pythermodb_settings.models import Component

from pyThermoModels.activity.pitzer import PitzerMultiComponentBinaryParameters, PitzerParameterSet
from pyThermoModels.activity.pitzer.gibbs import _j, evaluate_gibbs, excess_gibbs_molal, prepare_terms


def test_model_validates_mixed_formulation_inputs():
    from pythermodb_settings.models import Temperature
    from pyThermoModels.activity.pitzer import Pitzer

    components = [Component(name="Sodium", formula="Na{+}", state="aq"),
                  Component(name="Chloride", formula="Cl{-}", state="aq"),
                  Component(name="Carbon dioxide", formula="CO2", state="aq")]
    params = PitzerParameterSet(.3915, 1.2, {("Na{+}", "Cl{-}"): PitzerMultiComponentBinaryParameters(.0765, .2664, c_phi=.00127)}, {}, {}, {})
    model = Pitzer(components, formulation="multicomponent_pitzer_v2")
    state = {"temperature": Temperature(value=298.15, unit="K"),
             "molalities": {"Na{+}-aq": 1., "Cl{-}-aq": 1., "CO2-aq": 0.},
             "parameters": params, "strict_parameters": False}
    with pytest.raises(TypeError, match="strict_parameters must be a boolean"):
        model.cal({**state, "strict_parameters": "false"})
    for bad in (float("nan"), float("inf"), "1", True):
        with pytest.raises(ValueError, match="positive integers"):
            model.cal({**state, "stoichiometry": {"Na{+}-aq": bad}})
    with pytest.raises(ValueError, match="ions only"):
        model.cal({**state, "stoichiometry": {"Na{+}-aq": 1, "Cl{-}-aq": 1, "CO2-aq": 1}})
    _, details = model.cal({**state, "stoichiometry": {"Na{+}-aq": 1, "Cl{-}-aq": 1}})
    assert details["gamma_mean"] == pytest.approx(np.sqrt(details["activity_coefficients"]["Na{+}-aq"] * details["activity_coefficients"]["Cl{-}-aq"]))


def test_gibbs_derivative_and_reciprocity():
    components = [Component(name="Sodium", formula="Na{+}", state="aq"), Component(name="Chloride", formula="Cl{-}", state="aq")]
    params = PitzerParameterSet(.3915, 1.2, {("Na{+}", "Cl{-}"): PitzerMultiComponentBinaryParameters(.0765, .2664, c_phi=.00127)}, {}, {}, {})
    terms = prepare_terms(components, {}, params)
    charges = np.array([1., -1.])
    m = np.array([.8, .8])
    f = lambda x: excess_gibbs_molal(x, charges, params, terms)
    gradient = grad(f)(m)
    step = 1e-5
    for i in range(2):
        delta = np.eye(2)[i]*step
        assert gradient[i] == pytest.approx((f(m+delta)-f(m-delta))/(2*step), abs=1e-8)
    matrix = hessian(f)(m)
    assert matrix == pytest.approx(matrix.T, abs=1e-10)
    G, derivatives, phi, ln_aw = evaluate_gibbs(components, {"Na{+}-aq": .8, "Cl{-}-aq": .8}, params)
    assert phi == pytest.approx(1 + (m @ derivatives-G)/sum(m))
    assert ln_aw == pytest.approx(-.01801528*sum(m)*phi)


def test_j_pitzer_1975_reference():
    # Wolery et al. 2010, Appendix A, table A-3, Pitzer (1975) approximation.
    assert _j(2.) == pytest.approx(.2973357, abs=1e-7)


def test_neutral_self_pair_and_pure_water_limit():
    co2 = Component(name="Carbon dioxide", formula="CO2", state="aq")
    params = PitzerParameterSet(.3915, 1.2, {}, {}, {}, {}, lambda_nn={("CO2", "CO2"): .1})
    G, derivative, phi, ln_aw = evaluate_gibbs([co2], {"CO2-aq": 2.}, params)
    assert G == pytest.approx(.4)
    assert derivative[0] == pytest.approx(.4)
    assert phi == pytest.approx(1.2)
    assert ln_aw == pytest.approx(-.01801528*2*1.2)
    assert evaluate_gibbs([co2], {"CO2-aq": 0.}, params)[2:] == (1., 0.)
