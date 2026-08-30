import math

import numpy as np

from pyThermoModels.activity.uniquac import UNIQUAC
from pyThermoModels.activity.uniquac.local_composition import UNIQUACLocalComposition
from pyThermoModels.activity.nrtl.local_composition import NRTLLocalComposition


def test_uniquac_local_composition_tau_methods_are_exponential():
    model = UNIQUACLocalComposition(
        components=["A", "B"],
        component_idx={"A": 0, "B": 1},
    )
    temperature = 300.0
    a_ij = np.array([[0.0, 0.2], [0.4, 0.0]])
    b_ij = np.array([[0.0, 30.0], [60.0, 0.0]])
    c_ij = np.array([[0.0, 0.01], [0.02, 0.0]])
    d_ij = np.array([[0.0, 0.001], [0.002, 0.0]])
    dU_ij = np.array([[0.0, 1000.0], [2000.0, 0.0]])

    tau_m1, _ = model.cal_tau_ij_M1(temperature, dU_ij)
    tau_m2, _ = model.cal_tau_ij_M2(
        temperature, a_ij, b_ij, c_ij, d_ij
    )
    tau_m3, _ = model.cal_tau_ij_M3(temperature, a_ij)
    tau_m4, _ = model.cal_tau_ij_M4(temperature, a_ij, b_ij)
    tau_m5, _ = model.cal_tau_ij_M5(temperature, a_ij, b_ij, c_ij)
    tau_m6, _ = model.cal_tau_ij_M6(temperature, a_ij, b_ij, c_ij)

    assert np.allclose(np.diag(tau_m1), [1.0, 1.0])
    assert np.allclose(np.diag(tau_m2), [1.0, 1.0])
    assert np.allclose(np.diag(tau_m3), [1.0, 1.0])
    assert np.allclose(np.diag(tau_m4), [1.0, 1.0])
    assert np.allclose(np.diag(tau_m5), [1.0, 1.0])
    assert np.allclose(np.diag(tau_m6), [1.0, 1.0])

    assert np.isclose(tau_m1[0, 1], math.exp(-1000.0 / (8.314 * temperature)))
    assert np.isclose(tau_m2[0, 1], math.exp(
        0.2 + 30.0 / temperature + 0.01 * math.log(temperature)
        + 0.001 * temperature
    ))
    assert np.isclose(tau_m3[0, 1], math.exp(-0.2 / temperature))
    assert np.isclose(tau_m4[0, 1], math.exp(0.2 + 30.0 / temperature))
    assert np.isclose(tau_m5[0, 1], math.exp(
        0.2 + 30.0 / temperature + 0.01 / temperature**2
    ))
    assert np.isclose(tau_m6[0, 1], math.exp(
        0.2 + 30.0 / temperature + 0.01 * math.log(temperature)
    ))


def test_uniquac_parameter_builder_exposes_tau_methods():
    uniquac = UNIQUAC(components=["A", "B"])
    a_ij = {"A | A": 0.0, "A | B": 0.2, "B | A": 0.4, "B | B": 0.0}
    b_ij = {"A | A": 0.0, "A | B": 30.0, "B | A": 60.0, "B | B": 0.0}
    c_ij = {"A | A": 0.0, "A | B": 0.01, "B | A": 0.02, "B | B": 0.0}

    tau_ij, tau_ij_comp = uniquac.uniquac_parameter_builder.cal_tau_ij_M6(
        temperature=300.0,
        a_ij=a_ij,
        b_ij=b_ij,
        c_ij=c_ij,
    )

    assert np.allclose(np.diag(tau_ij), [1.0, 1.0])
    assert tau_ij_comp["A | A"] == 1.0
    assert tau_ij_comp["B | B"] == 1.0
    assert np.isclose(tau_ij[0, 1], math.exp(
        0.2 + 30.0 / 300.0 + 0.01 * math.log(300.0)
    ))


def test_nrtl_tau_diagonal_remains_zero():
    model = NRTLLocalComposition(
        components=["A", "B"],
        component_idx={"A": 0, "B": 1},
    )
    dg_ij = np.array([[0.0, 1000.0], [2000.0, 0.0]])

    tau_ij, tau_ij_comp = model.cal_tau_ij_M1(
        temperature=300.0,
        dg_ij=dg_ij,
    )

    assert np.allclose(np.diag(tau_ij), [0.0, 0.0])
    assert tau_ij_comp["A | A"] == 0
    assert tau_ij_comp["B | B"] == 0
