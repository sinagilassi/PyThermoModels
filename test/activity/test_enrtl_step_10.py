import numpy as np
import pytest

from pyThermoModels.activity import ENRTL
from pyThermoModels.activity.enrtl.caloric import excess_enthalpy_not_available
from pyThermoModels.activity.enrtl.derivatives import tau_temperature_derivative


class FakeComponent:
    def __init__(self, key, charge):
        self.key = key
        self.charge = charge
        self.species_type = []

    def get_formula_state(self):
        return self.key

    def get_net_charge(self):
        return self.charge

    def is_radical(self):
        return False


def _neutral_model_input():
    return {
        "composition_representation": "true_species",
        "temperature": [298.15, "K"],
        "mole_fraction": {"A-l": 0.4, "B-l": 0.6},
        "tau_ij": [[0.0, 0.2], [0.1, 0.0]],
        "alpha_ij": [[0.0, 0.3], [0.3, 0.0]],
        "long_range": {
            "model": "pitzer_debye_huckel",
            "basis": "mole_fraction",
            "A_phi": 0.392,
        },
    }


def _ionic_model_input():
    return {
        "composition_representation": "true_species",
        "temperature": [298.15, "K"],
        "mole_fraction": {"H2O-l": 0.98, "Na{+}-aq": 0.01, "Cl{-}-aq": 0.01},
        "molality": {"H2O-l": 0.0, "Na{+}-aq": 0.1, "Cl{-}-aq": 0.1},
        "tau_ij": [[0.0, 0.10, 0.20], [0.15, 0.0, 0.30], [0.25, 0.40, 0.0]],
        "alpha_ij": [[0.0, 0.2, 0.2], [0.2, 0.0, 0.2], [0.2, 0.2, 0.0]],
        "long_range": {
            "model": "pitzer_debye_huckel",
            "basis": "molality",
            "A_phi": 0.392,
        },
    }


def test_neutral_nrtl_limit_excess_gibbs_matches_activity_identity():
    model = ENRTL([FakeComponent("A-l", 0), FakeComponent("B-l", 0)])
    _, details = model.cal(_neutral_model_input())

    expected = np.dot([0.4, 0.6], details["ln_gamma_local_composition"])
    assert details["excess_gibbs_RT"] == pytest.approx(expected)
    assert details["excess_gibbs_contributions_RT"]["total_status"] == "complete"
    assert details["reference_state"] == "true_species_activity_model"
    assert details["activity_convention"] == "natural_log_activity_coefficients"
    assert details["species_basis"] == "true_species"


def test_charged_local_excess_gibbs_matches_chen_evans_kernel_diagnostics():
    model = ENRTL([
        FakeComponent("H2O-l", 0),
        FakeComponent("Na{+}-aq", 1),
        FakeComponent("Cl{-}-aq", -1),
    ])
    _, details = model.cal(_ionic_model_input())

    contributions = details["excess_gibbs_contributions_RT"]
    assert contributions["local_composition"] == pytest.approx(
        details["local_composition_diagnostics"]["gE_local_composition_RT"]
    )
    assert details["excess_gibbs_RT"] is None
    assert contributions["total_status"] == "unavailable"
    assert contributions["long_range"] is None


def test_excess_gibbs_rejects_apparent_species_and_implicit_ionic_identity():
    model = ENRTL([
        FakeComponent("Na{+}-aq", 1),
        FakeComponent("Cl{-}-aq", -1),
    ])

    with pytest.raises(NotImplementedError, match="composition_representation"):
        model.excess_gibbs_free_energy(
            [0.5, 0.5], [0.1, 0.1], composition_representation="apparent_species"
        )
    with pytest.raises(NotImplementedError, match="not a validated ionic"):
        model.excess_gibbs_free_energy([0.5, 0.5], [0.1, 0.1])

    result = model.excess_gibbs_free_energy(
        [0.5, 0.5], [0.1, 0.1], allow_ionic_identity=True
    )
    assert result["value"] == pytest.approx(0.1)


@pytest.mark.parametrize(
    ("correlation", "kwargs"),
    [
        ("gibbs_energy", {"dg_ij": [[1200.0]]}),
        ("inverse_temperature", {"a_ij": [[0.2]], "b_ij": [[350.0]]}),
        ("inverse_temperature_squared", {"a_ij": [[0.2]], "b_ij": [[350.0]], "c_ij": [[4000.0]]}),
        ("inverse_log_temperature", {"a_ij": [[0.2]], "b_ij": [[350.0]], "c_ij": [[0.05]]}),
        ("extended_temperature", {"a_ij": [[0.2]], "b_ij": [[350.0]], "c_ij": [[0.05]], "d_ij": [[0.001]]}),
    ],
)
def test_tau_temperature_derivatives_match_finite_difference(correlation, kwargs):
    temperature = 325.0
    step = 1e-3
    analytic = tau_temperature_derivative(temperature, correlation, **kwargs)
    numerical = (_tau_value(temperature + step, correlation, **kwargs) - _tau_value(temperature - step, correlation, **kwargs)) / (2.0 * step)

    assert analytic == pytest.approx(numerical, rel=1e-8, abs=1e-10)


def test_direct_tau_derivative_requires_explicit_data():
    with pytest.raises(NotImplementedError, match="d_tau_dT"):
        tau_temperature_derivative(298.15, "direct_tau")
    assert tau_temperature_derivative(298.15, "direct_tau", d_tau_dT=[[0.01]]) == pytest.approx([[0.01]])


def test_excess_enthalpy_remains_gated():
    model = ENRTL([FakeComponent("A-l", 0), FakeComponent("B-l", 0)])
    with pytest.raises(NotImplementedError, match="not available"):
        model.excess_enthalpy()
    with pytest.raises(NotImplementedError, match="not available"):
        excess_enthalpy_not_available()


def test_enrtl_rejects_invalid_composition_before_gibbs_diagnostics():
    model = ENRTL([FakeComponent("A-l", 0), FakeComponent("B-l", 0)])
    model_input = _neutral_model_input()
    model_input["mole_fraction"] = {"A-l": 0.4, "B-l": 0.5}

    with pytest.raises(ValueError, match="sum to 1.0"):
        model.cal(model_input)


def test_enrtl_rejects_non_electroneutral_state_before_gibbs_diagnostics():
    model = ENRTL([
        FakeComponent("H2O-l", 0),
        FakeComponent("Na{+}-aq", 1),
        FakeComponent("Cl{-}-aq", -1),
    ])
    model_input = _ionic_model_input()
    model_input["mole_fraction"] = {
        "H2O-l": 0.97,
        "Na{+}-aq": 0.02,
        "Cl{-}-aq": 0.01,
    }

    with pytest.raises(ValueError, match="not electrically neutral"):
        model.cal(model_input)

def _tau_value(temperature, correlation, **kwargs):
    if correlation == "gibbs_energy":
        return np.asarray(kwargs["dg_ij"], dtype=float) / (8.31446261815324 * temperature)
    value = np.asarray(kwargs["a_ij"], dtype=float) + np.asarray(kwargs["b_ij"], dtype=float) / temperature
    if correlation == "inverse_temperature":
        return value
    if correlation == "inverse_temperature_squared":
        return value + np.asarray(kwargs["c_ij"], dtype=float) / temperature**2
    if correlation == "inverse_log_temperature":
        return value + np.asarray(kwargs["c_ij"], dtype=float) * np.log(temperature)
    return value + np.asarray(kwargs["c_ij"], dtype=float) * np.log(temperature) + np.asarray(kwargs["d_ij"], dtype=float) * temperature
