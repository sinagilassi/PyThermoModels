# import packages/modules
from examples.source.model_source_1 import model_source, model_source_dict, C3H8
import sys
from pathlib import Path
from typing import Dict
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
from pyThermoModels.core import (
    calc_gas_fugacity,
    calc_residual_properties,
    check_component_eos_roots,
)

# SECTION: local example imports
# NOTE: make this file runnable directly from the examples/eos-models folder.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# ! model source & components

# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# NOTE: examples
# eos model
eos_model = "SRK"

# phase
phase = "VAPOR"

# temperature
temperature = Temperature(value=300.1, unit='K')
# pressure
pressure = Pressure(value=9.99, unit='bar')

# ------------------------------------------------
# NOTE: eos root analysis
# ------------------------------------------------
# ! new method
res = check_component_eos_roots(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name=eos_model,
    component_key='Name-State',
    phase=phase,
    mode='log',
)
print(res)

# ------------------------------------------------
# NOTE: calculation
# ------------------------------------------------
# ! new method
res = calc_gas_fugacity(
    component=C3H8,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name=eos_model,
    component_key='Name-State',
    phase=phase,
    mode='log',
)
print(res)
# ------------------------------------------------
# NOTE: residual/departure properties
# ------------------------------------------------
# ! Phase 1 pure-fluid PR/SRK residual properties
# NOTE: This API is explicit about the homogeneous EOS root through phase.
res = calc_residual_properties(
    component=C3H8,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name=eos_model,
    component_key='Name-State',
    phase=phase,
    mode='log',
)
print(res)
