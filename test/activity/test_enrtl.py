import numpy as np
import pytest

from pyThermoModels.activity import ActivityCore, ENRTL


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


def test_activity_core_exposes_enrtl_for_legacy_string_components():
    activity = ActivityCore(
        datasource={},
        equationsource={},
        components=["A", "B"],
    )

    assert isinstance(activity.enrtl, ENRTL)
    assert isinstance(activity.select("ENRTL"), ENRTL)


def test_enrtl_neutral_calculation_uses_log_contribution_sum():
    model = ENRTL([
        FakeComponent("A-l", 0),
        FakeComponent("B-l", 0),
    ])

    res, details = model.cal({
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
    })

    assert res["model"] == "ENRTL"
    assert res["formulation"] == "chen_evans_1986"
    assert np.allclose(
        details["ln_gamma_total"],
        details["ln_gamma_long_range"] + details["ln_gamma_local_composition"],
    )
    assert details["net_charge"] == pytest.approx(0.0)


def test_enrtl_mean_ionic_activity_coefficient():
    model = ENRTL([
        FakeComponent("Na{+}-aq", 1),
        FakeComponent("Cl{-}-aq", -1),
    ])

    value = model.mean_ionic_activity_coefficient(
        gamma={"Na{+}-aq": 0.8, "Cl{-}-aq": 0.9},
        cation="Na{+}-aq",
        anion="Cl{-}-aq",
    )

    assert value == pytest.approx((0.8 * 0.9) ** 0.5)


def test_enrtl_charged_calculation_combines_chen_evans_contributions():
    model = ENRTL([
        FakeComponent("H2O-l", 0),
        FakeComponent("Na{+}-aq", 1),
        FakeComponent("Cl{-}-aq", -1),
    ])

    res, details = model.cal({
        "composition_representation": "true_species",
        "temperature": [298.15, "K"],
        "mole_fraction": {
            "H2O-l": 0.98,
            "Na{+}-aq": 0.01,
            "Cl{-}-aq": 0.01,
        },
        "molality": {
            "H2O-l": 0.0,
            "Na{+}-aq": 0.1,
            "Cl{-}-aq": 0.1,
        },
        "tau_ij": [
            [0.0, 0.10, 0.20],
            [0.15, 0.00, 0.30],
            [0.25, 0.40, 0.00],
        ],
        "alpha_ij": [
            [0.0, 0.2, 0.2],
            [0.2, 0.0, 0.2],
            [0.2, 0.2, 0.0],
        ],
        "long_range": {
            "model": "pitzer_debye_huckel",
            "basis": "molality",
            "A_phi": 0.392,
        },
    })

    assert res["model"] == "ENRTL"
    assert np.all(np.isfinite(res["value"]))
    assert np.allclose(
        details["ln_gamma_total"],
        details["ln_gamma_long_range"] + details["ln_gamma_local_composition"],
    )
    assert details["local_composition_diagnostics"]["formulation"] == "chen_evans_1986"
