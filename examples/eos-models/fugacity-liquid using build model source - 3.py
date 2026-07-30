# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
from pyThermoModels.core import calc_liquid_fugacity, check_component_eos_roots
# locals
from examples.source.model_source_2 import model_source, C3H8

# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================

# NOTE: Example 3.13s
# temperature
temperature = Temperature(value=340, unit='K')
# pressure
pressure = Pressure(value=30, unit='bar')

# ------------------------------------------------
# NOTE: eos root analysis
# ------------------------------------------------
# ! new method
res = check_component_eos_roots(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
)
print(res)

# ------------------------------------------------
# NOTE: calculation
# ------------------------------------------------
# ! new method
res = calc_liquid_fugacity(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    liquid_fugacity_mode='Poynting',
)
print(res)
