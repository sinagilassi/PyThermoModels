from pyThermoModels.core.activity_parameter_calculators import (
    _activity_tau_configuration,
)
from pyThermoModels.core.activity_methods import (
    _activity_tau_correlation_default,
)
from pyThermoModels.activity.nrtl.parameter_builder import NRTLParameterBuilder
from pyThermoModels.activity.uniquac.parameter_builder import UNIQUACParameterBuilder
from pyThermoModels.utils.utility import (
    map_tau_correlation_to_method,
    map_uniquac_tau_correlation_to_method,
)

import numpy as np


def test_nrtl_tau_correlation_maps_to_local_composition_method_ids():
    assert map_tau_correlation_to_method("direct_tau") == "M0"
    assert map_tau_correlation_to_method("gibbs_energy") == "M1"
    assert map_tau_correlation_to_method("extended_temperature") == "M2"
    assert map_tau_correlation_to_method("inverse_temperature") == "M3"
    assert map_tau_correlation_to_method("inverse_temperature_squared") == "M4"
    assert map_tau_correlation_to_method("inverse_log_temperature") == "M5"


def test_uniquac_tau_correlation_maps_to_local_composition_method_ids():
    assert map_uniquac_tau_correlation_to_method("direct_tau") == "M0"
    assert map_uniquac_tau_correlation_to_method("gibbs_energy") == "M1"
    assert map_uniquac_tau_correlation_to_method("extended_temperature") == "M2"
    assert map_uniquac_tau_correlation_to_method("inverse_temperature") == "M4"
    assert map_uniquac_tau_correlation_to_method("inverse_temperature_squared") == "M5"
    assert map_uniquac_tau_correlation_to_method("inverse_log_temperature") == "M6"


def test_activity_tau_configuration_passes_descriptive_names():
    nrtl_config = _activity_tau_configuration(
        model_name="NRTL",
        tau_correlation="inverse_temperature",
    )
    uniquac_config = _activity_tau_configuration(
        model_name="UNIQUAC",
        tau_correlation="inverse_temperature",
    )

    assert nrtl_config["tau_correlation"] == "inverse_temperature"
    assert uniquac_config["tau_correlation"] == "inverse_temperature"


def test_activity_tau_correlation_default_uses_direct_tau_for_nrtl_tau_source():
    tau_correlation = _activity_tau_correlation_default(
        model_name="NRTL",
        datasource={
            "A|B": {
                "tau": [[0.0, 0.2], [0.4, 0.0]],
                "dg": [[0.0, 100.0], [200.0, 0.0]],
                "alpha": [[0.0, 0.3], [0.3, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
    )

    assert tau_correlation == "direct_tau"


def test_activity_tau_correlation_default_honors_nrtl_custom_source_priority():
    tau_correlation = _activity_tau_correlation_default(
        model_name="NRTL",
        datasource={
            "A|B": {
                "tau": [[0.0, 0.2], [0.4, 0.0]],
                "dg": [[0.0, 100.0], [200.0, 0.0]],
                "a": [[0.0, 0.1], [0.2, 0.0]],
                "b": [[0.0, 10.0], [20.0, 0.0]],
                "c": [[0.0, 0.01], [0.02, 0.0]],
                "d": [[0.0, 0.001], [0.002, 0.0]],
                "alpha": [[0.0, 0.3], [0.3, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
        tau_source_priority=("coefficients", "dg", "tau"),
    )

    assert tau_correlation == "extended_temperature"


def test_activity_tau_correlation_default_uses_gibbs_energy_for_nrtl_dg_source():
    tau_correlation = _activity_tau_correlation_default(
        model_name="NRTL",
        datasource={
            "A|B": {
                "dg": [[0.0, 100.0], [200.0, 0.0]],
                "alpha": [[0.0, 0.3], [0.3, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
    )

    assert tau_correlation == "gibbs_energy"


def test_activity_tau_correlation_default_uses_extended_temperature_for_nrtl_coefficients():
    tau_correlation = _activity_tau_correlation_default(
        model_name="NRTL",
        datasource={
            "A|B": {
                "a": [[0.0, 0.1], [0.2, 0.0]],
                "b": [[0.0, 10.0], [20.0, 0.0]],
                "c": [[0.0, 0.01], [0.02, 0.0]],
                "d": [[0.0, 0.001], [0.002, 0.0]],
                "alpha": [[0.0, 0.3], [0.3, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
    )

    assert tau_correlation == "extended_temperature"


def test_activity_tau_correlation_default_uses_direct_tau_for_uniquac_tau_source():
    tau_correlation = _activity_tau_correlation_default(
        model_name="UNIQUAC",
        datasource={
            "A|B": {
                "tau": [[1.0, 0.2], [0.4, 1.0]],
                "dU": [[0.0, 100.0], [200.0, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
    )

    assert tau_correlation == "direct_tau"


def test_activity_tau_correlation_default_honors_uniquac_custom_source_priority():
    tau_correlation = _activity_tau_correlation_default(
        model_name="UNIQUAC",
        datasource={
            "A|B": {
                "tau": [[1.0, 0.2], [0.4, 1.0]],
                "dU": [[0.0, 100.0], [200.0, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
        tau_source_priority=("dU", "tau", "coefficients"),
    )

    assert tau_correlation == "gibbs_energy"


def test_activity_tau_correlation_default_uses_gibbs_energy_for_uniquac_du_source():
    tau_correlation = _activity_tau_correlation_default(
        model_name="UNIQUAC",
        datasource={
            "A|B": {
                "dU": [[0.0, 100.0], [200.0, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
    )

    assert tau_correlation == "gibbs_energy"


def test_activity_tau_correlation_default_keeps_extended_temperature_for_coefficients():
    tau_correlation = _activity_tau_correlation_default(
        model_name="UNIQUAC",
        datasource={
            "A|B": {
                "a": [[0.0, 0.1], [0.2, 0.0]],
                "b": [[0.0, 10.0], [20.0, 0.0]],
                "c": [[0.0, 0.01], [0.02, 0.0]],
                "d": [[0.0, 0.001], [0.002, 0.0]],
            },
        },
        mixture_ids={"Name": "A|B"},
    )

    assert tau_correlation == "extended_temperature"


def test_nrtl_parameter_builder_maps_descriptive_name_internally():
    builder = NRTLParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        temperature=[300.0, "K"],
        tau_correlation="inverse_temperature",
        include_alpha=False,
        model_input={
            "a_ij": [[0.0, 0.2], [0.4, 0.0]],
            "b_ij": [[0.0, 30.0], [60.0, 0.0]],
        },
    )

    assert inputs["tau_ij"][0, 0] == 0
    assert np.isclose(inputs["tau_ij"][0, 1], 0.2 + 30.0 / 300.0)


def test_nrtl_parameter_builder_accepts_direct_tau_name():
    builder = NRTLParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        tau_correlation="direct_tau",
        include_alpha=False,
        model_input={
            "tau_ij": [[0.0, 0.2], [0.4, 0.0]],
        },
    )

    assert inputs["tau_ij"][0, 0] == 0.0
    assert inputs["tau_ij"][0, 1] == 0.2


def test_nrtl_parameter_builder_can_prefer_coefficients_over_direct_tau():
    builder = NRTLParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        temperature=[300.0, "K"],
        tau_correlation="extended_temperature",
        include_alpha=False,
        model_input={
            "tau_ij": [[0.0, 99.0], [99.0, 0.0]],
            "a_ij": [[0.0, 0.2], [0.4, 0.0]],
            "b_ij": [[0.0, 30.0], [60.0, 0.0]],
            "c_ij": [[0.0, 0.01], [0.02, 0.0]],
            "d_ij": [[0.0, 0.001], [0.002, 0.0]],
        },
    )

    expected = 0.2 + 30.0 / 300.0 + 0.01 * np.log(300.0) + 0.001 * 300.0
    assert np.isclose(inputs["tau_ij"][0, 1], expected)


def test_uniquac_parameter_builder_uses_default_extended_temperature_method():
    builder = UNIQUACParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        temperature=[300.0, "K"],
        include_pure_component_parameters=False,
        model_input={
            "a_ij": [[0.0, 0.2], [0.4, 0.0]],
            "b_ij": [[0.0, 30.0], [60.0, 0.0]],
            "c_ij": [[0.0, 0.01], [0.02, 0.0]],
            "d_ij": [[0.0, 0.001], [0.002, 0.0]],
        },
    )

    expected = np.exp(
        0.2 + 30.0 / 300.0 + 0.01 * np.log(300.0) + 0.001 * 300.0
    )
    assert inputs["tau_ij"][0, 0] == 1.0
    assert np.isclose(inputs["tau_ij"][0, 1], expected)


def test_uniquac_parameter_builder_maps_descriptive_name_internally():
    builder = UNIQUACParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        temperature=[300.0, "K"],
        tau_correlation="inverse_temperature",
        include_pure_component_parameters=False,
        model_input={
            "a_ij": [[0.0, 0.2], [0.4, 0.0]],
            "b_ij": [[0.0, 30.0], [60.0, 0.0]],
        },
    )

    assert inputs["tau_ij"][0, 0] == 1.0
    assert np.isclose(inputs["tau_ij"][0, 1], np.exp(0.2 + 30.0 / 300.0))


def test_uniquac_parameter_builder_accepts_direct_tau_name():
    builder = UNIQUACParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        temperature=[300.0, "K"],
        tau_correlation="direct_tau",
        include_pure_component_parameters=False,
        model_input={
            "tau_ij": [[1.0, 0.2], [0.4, 1.0]],
        },
    )

    assert inputs["tau_ij"][0, 0] == 1.0
    assert inputs["tau_ij"][0, 1] == 0.2


def test_uniquac_parameter_builder_can_prefer_energy_over_direct_tau():
    builder = UNIQUACParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    inputs = builder.inputs_generator(
        temperature=[300.0, "K"],
        tau_correlation="gibbs_energy",
        include_pure_component_parameters=False,
        model_input={
            "tau_ij": [[1.0, 99.0], [99.0, 1.0]],
            "dU_ij": [[0.0, 100.0], [200.0, 0.0]],
        },
    )

    expected = np.exp(-100.0 / (8.314 * 300.0))
    assert np.isclose(inputs["tau_ij"][0, 1], expected)


def test_uniquac_parameter_builder_reports_mismatched_gibbs_energy_source():
    builder = UNIQUACParameterBuilder(
        components=["A", "B"],
        comp_idx={"A": 0, "B": 1},
        datasource={},
        equationsource={},
    )

    try:
        builder.inputs_generator(
            temperature=[300.0, "K"],
            tau_correlation="gibbs_energy",
            include_pure_component_parameters=False,
            model_input={
                "a_ij": [[0.0, 0.2], [0.4, 0.0]],
                "b_ij": [[0.0, 30.0], [60.0, 0.0]],
                "c_ij": [[0.0, 0.01], [0.02, 0.0]],
                "d_ij": [[0.0, 0.001], [0.002, 0.0]],
            },
        )
    except Exception as exc:
        assert "requires dU_ij for UNIQUAC" in str(exc)
    else:
        raise AssertionError("Expected gibbs_energy with coefficient matrices to fail")


def test_uniquac_tau_correlation_rejects_raw_method_id():
    try:
        map_uniquac_tau_correlation_to_method("M4")
    except ValueError as exc:
        assert "Unsupported UNIQUAC tau correlation" in str(exc)
    else:
        raise AssertionError("Expected raw UNIQUAC method id to be rejected")
