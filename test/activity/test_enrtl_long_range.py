import pytest

from pyThermoModels.activity.enrtl_long_range import ENRTLLongRange


def test_pitzer_debye_huckel_infinite_dilution_limit():
    model = ENRTLLongRange()

    result = model.cal_ln_gamma_long_range(
        ionic_strength=0.0,
        charges={"Na{+}-aq": 1, "Cl{-}-aq": -1, "H2O-l": 0},
        model="pitzer_debye_huckel",
        A_phi=0.392,
    )

    assert result == {
        "Na{+}-aq": 0.0,
        "Cl{-}-aq": 0.0,
        "H2O-l": 0.0,
    }


def test_debye_huckel_zero_charge_species_has_zero_contribution():
    model = ENRTLLongRange()

    result = model.cal_ln_gamma_long_range(
        ionic_strength=0.1,
        charges={"H2O-l": 0},
        model="debye_huckel",
        A=0.5,
        B=0.3,
    )

    assert result["H2O-l"] == 0.0


def test_debye_huckel_requires_ion_size_for_charged_species():
    model = ENRTLLongRange()

    with pytest.raises(ValueError, match="Missing ion_size"):
        model.cal_ln_gamma_long_range(
            ionic_strength=0.1,
            charges={"Na{+}-aq": 1},
            model="debye_huckel",
            A=0.5,
            B=0.3,
        )
