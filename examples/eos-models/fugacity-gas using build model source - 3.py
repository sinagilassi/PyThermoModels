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
from examples.source.model_source_1 import model_source, model_source_dict, C3H8

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


# ========================================
# NOTE: INITIALIZE OBJECT
# ========================================
tm = ptm.init()
# log
print("tm: ", tm)

# NOTE: check reference
# all reference
res_ = tm.references
print(res_)

# eos model
eos_model = 'PR'

res_ = tm.check_fugacity_reference(eos_model)
print(res_)

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# NOTE: examples
# Example 3.9 (page 100), Introduction to Chemical Engineering Thermodynamics
# by J.M. Smith, H.C. Van Ness, M.M. Abbott

# Example 10.9, Introduction to Chemical Engineering Thermodynamics

# Example 3.11 (page 68), The Thermodynamics of Phase and Reaction Equilibria

# NOTE: init eos
eos = ptm.eos()
print(ptm.eos.metadata)

# NOTE: check details
print(eos.cal_fugacity.metadata)

# model input
# eos model
eos_model = 'PR'

# component phase
phase = "VAPOR"

# NOTE: Example 3.9
# component
component = "n-butane"
# temperature [K]
T = 350
# pressure bar
P = 9.653800

# NOTE: Example 10.9
component = "1-butene"
# temperature [C]
T = 200
# pressure bar
P = 70

# NOTE: Example 3.11
component = "propane-g"
# phase
phase = "VAPOR-LIQUID"
# temperature [K]
T = 300.1
# pressure [bar]
P = 9.99

# NOTE: Example 3.13
# # component
# component = "CH4"
# # temperature [K]
# T = 340
# # pressure [bar]
# P = 30


# SECTION: model input
# temperature
temperature = Temperature(value=T, unit='K')
# pressure
pressure = Pressure(value=P, unit='bar')

model_input = {
    # "phase": phase,
    "component": component,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# ------------------------------------------------
# NOTE: eos root analysis
# ------------------------------------------------
# ! old method
res = eos.check_eos_roots_single_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source_dict
)
print(res)

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
# NOTE: gas fugacity calculation method
# ! old method
res = eos.cal_fugacity(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source_dict
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
