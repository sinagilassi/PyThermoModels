from pyThermoModels.core.activity_parameter_calculators import (
    _activity_tau_configuration,
)
from pyThermoModels.activity.nrtl_parameter_builder import NRTLParameterBuilder
from pyThermoModels.utils.utility import (
    map_tau_correlation_to_method,
    map_uniquac_tau_correlation_to_method,
)

import numpy as np


def test_nrtl_tau_correlation_maps_to_local_composition_method_ids():
    assert map_tau_correlation_to_method("gibbs_energy") == "M1"
    assert map_tau_correlation_to_method("extended_temperature") == "M2"
    assert map_tau_correlation_to_method("inverse_temperature") == "M3"
    assert map_tau_correlation_to_method("inverse_temperature_squared") == "M4"
    assert map_tau_correlation_to_method("inverse_log_temperature") == "M5"


def test_uniquac_tau_correlation_maps_to_local_composition_method_ids():
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


def test_uniquac_tau_correlation_rejects_raw_method_id():
    try:
        map_uniquac_tau_correlation_to_method("M4")
    except ValueError as exc:
        assert "Unsupported UNIQUAC tau correlation" in str(exc)
    else:
        raise AssertionError("Expected raw UNIQUAC method id to be rejected")
