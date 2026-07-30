# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
from pyThermoModels.core import calc_gas_fugacity, check_component_eos_roots
# ! model source & components
from examples.source.model_source_1 import model_source, C3H8

# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# ! LOAD THERMODB
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, '..', 'thermodb')
print(thermodb_dir)

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# NOTE: examples
# phase
phase = "VAPOR-LIQUID"

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
    component_key='Name-State',
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
    component_key='Name-State'
)
print(res)
