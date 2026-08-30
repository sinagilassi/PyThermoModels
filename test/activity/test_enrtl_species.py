import pytest

from pyThermoModels.activity.enrtl.component_adapter import ENRTLComponentAdapter


class FakeComponent:
    def __init__(self, key, charge, species_type=None):
        self.key = key
        self.charge = charge
        self.species_type = species_type or []

    def get_formula_state(self):
        return self.key

    def get_net_charge(self):
        return self.charge

    def is_radical(self):
        return "radical" in self.species_type


def test_enrtl_component_adapter_consumes_component_charge_metadata():
    components = [
        FakeComponent("H2O-l", 0),
        FakeComponent("Na{+}-aq", 1),
        FakeComponent("Cl{-}-aq", -1),
    ]

    adapter = ENRTLComponentAdapter(components)

    assert adapter.component_keys == ["H2O-l", "Na{+}-aq", "Cl{-}-aq"]
    assert adapter.charges.tolist() == [0, 1, -1]
    assert adapter.charges_dict == {
        "H2O-l": 0,
        "Na{+}-aq": 1,
        "Cl{-}-aq": -1,
    }


def test_enrtl_component_adapter_rejects_conflicting_charge_override():
    components = [FakeComponent("Na{+}-aq", 1)]

    with pytest.raises(ValueError, match="disagrees with component metadata"):
        ENRTLComponentAdapter(
            components,
            charge_overrides={"Na{+}-aq": -1},
        )


def test_enrtl_component_adapter_rejects_zwitterions_for_v1():
    adapter = ENRTLComponentAdapter(
        [FakeComponent("glycine-aq", 0, species_type=["zwitterion"])]
    )

    with pytest.raises(ValueError, match="zwitterion"):
        adapter.validate_species_support()
