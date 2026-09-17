import numpy as np
import pytest

from pyThermoModels.activity.pitzer import (
    Pitzer,
    PitzerBinaryParameters,
    build_binary_ion_molalities,
    calc_ionic_strength,
    calc_ln_gamma_mean_binary,
    calc_osmotic_coefficient_binary,
    validate_electroneutrality,
)


NACL = PitzerBinaryParameters(
    beta0=0.0765,
    beta1=0.2664,
    c_phi=0.00127,
    alpha=2.0,
    A_phi=0.3915,
    b=1.2,
)


def _nacl_input(molality=1.0):
    return {
        "temperature": [298.15, "K"],
        "salt_molality": molality,
        "charges": [1, -1],
        "stoichiometry": [1, 1],
        "beta0": 0.0765,
        "beta1": 0.2664,
        "c_phi": 0.00127,
        "alpha": 2.0,
        "A_phi": 0.3915,
        "b": 1.2,
    }


def test_nacl_ionic_strength_equals_salt_molality():
    ion_molalities = build_binary_ion_molalities(1.0, 1, 1)

    assert calc_ionic_strength(ion_molalities, np.asarray([1, -1])) == pytest.approx(1.0)


def test_binary_state_requires_electroneutrality():
    validate_electroneutrality(np.asarray([1.0, 1.0]), np.asarray([1, -1]))

    with pytest.raises(ValueError, match="not electrically neutral"):
        validate_electroneutrality(np.asarray([1.0, 0.8]), np.asarray([1, -1]))


def test_nacl_1m_binary_equations_regression():
    phi = calc_osmotic_coefficient_binary(1.0, 1.0, 1, -1, 1, 1, NACL)
    ln_gamma = calc_ln_gamma_mean_binary(1.0, 1.0, 1, -1, 1, 1, NACL)

    assert phi == pytest.approx(0.9358687739996882, rel=1e-12, abs=1e-12)
    assert ln_gamma == pytest.approx(-0.4223446328193487, rel=1e-12, abs=1e-12)


def test_zero_molality_returns_exact_ideal_limit():
    model = Pitzer(["Na+", "Cl-"])
    result, details = model.cal(_nacl_input(0.0))

    assert details["ionic_strength"] == 0.0
    assert details["osmotic_coefficient"] == 1.0
    assert details["ln_gamma_mean"] == 0.0
    assert result["value"] == 1.0
    assert details["water_activity"] == 1.0


def test_nacl_1m_model_regression():
    model = Pitzer(["Na+", "Cl-"])
    result, details = model.cal(_nacl_input())

    assert result["model"] == "PITZER"
    assert result["formulation"] == "binary_single_alpha_v1"
    assert result["value"] == pytest.approx(0.6555080908595792, rel=1e-12)
    assert details["ionic_strength"] == pytest.approx(1.0)
    assert details["osmotic_coefficient"] == pytest.approx(0.9358687739996882, rel=1e-12)
    assert details["water_activity"] == pytest.approx(0.9668423024271156, rel=1e-12)


def test_rejects_unsupported_binary_scope():
    model = Pitzer(["Mg+2", "SO4-2"])
    model_input = _nacl_input()
    model_input["charges"] = [2, -2]

    with pytest.raises(NotImplementedError, match="2:2"):
        model.cal(model_input)

    model_input = _nacl_input()
    model_input["stoichiometry"] = [1, 2]
    with pytest.raises(ValueError, match="not electrically neutral"):
        Pitzer(["Ca+2", "Cl-"]).cal(model_input)
