# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
from pyThermoModels.core import calc_mixture_fugacity, check_multi_component_eos_roots
# ! v
from examples.source.model_source_2 import model_source, CO2, C4H10

# version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# ! CALCULATE FUGACITY FOR MULTI COMPONENT
# =======================================
# NOTE: Reference
# NOTE: Example 9.5 (page 257), The Thermodynamics of Phase and Reaction Equilibria
# eos model
eos_model = 'RK'
# feed spec
N0s = {
    'CO2-g': 0.15,
    'n-butane-g': 0.85,
}
# temperature [K]
temperature = Temperature(value=444, unit='K')

# pressure [bar]
pressure = Pressure(value=10, unit='bar')

# binary interaction parameter
k_ij = [[0, 0.18],
        [0.18, 0]]


# SECTION: components
# ! CO2
CO2 = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g',
    mole_fraction=N0s.get('CO2-g', 0)
)
# ! n-butane
C4H10 = Component(
    name='n-butane',
    formula='C4H10',
    state='g',
    mole_fraction=N0s.get('n-butane-g', 0)
)
# ! components
components = [CO2, C4H10]

# =======================================
# EOS ROOT ANALYSIS
# =======================================
# ! new method using build model source
res_ = check_multi_component_eos_roots(
    components=components,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    component_key='Name-State'  # ! component key (optional)
)
print(f"new method using build model source:")
print(res_)

# =======================================
# FUGACITY CALCULATION
# =======================================
# ! new method using build model source and component key
res = calc_mixture_fugacity(
    components=components,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    component_key='Name-State'
)
print(res)
