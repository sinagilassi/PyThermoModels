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

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(thermodb_dir, 'carbon dioxide-g.pkl')

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(thermodb_dir, 'acetylene-g.pkl')

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-g.pkl')

# ! ethanol
# thermodb file name
ethanol_thermodb_file = os.path.join(thermodb_dir, 'ethanol-l.pkl')

# ! methanol
# thermodb file name
methanol_thermodb_file = os.path.join(thermodb_dir, 'methanol-g.pkl')

# ! 1-butene
# thermodb file name
butene_thermodb_file = os.path.join(thermodb_dir, '1-butene-g.pkl')

# ! propane
# thermodb file name
propane_thermodb_file = os.path.join(thermodb_dir, 'propane-g.pkl')

# ! methane
# thermodb file name
methane_thermodb_file = os.path.join(thermodb_dir, 'methane-g.pkl')

# =======================================
# SECTION: COMPONENTS THERMODB SOURCE
# =======================================
# NOTE: carbon dioxide
CO2_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='carbon dioxide',
        formula='CO2',
        state='g'
    ),
    source=CO2_thermodb_file
)

# NOTE: acetylene
acetylene_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='acetylene',
        formula='C2H2',
        state='g'
    ),
    source=acetylene_thermodb_file
)

# NOTE: n-butane
n_butane_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='n-butane',
        formula='C4H10',
        state='g'
    ),
    source=n_butane_thermodb_file
)

# NOTE: ethanol
ethanol_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='ethanol',
        formula='C2H5OH',
        state='l'
    ),
    source=ethanol_thermodb_file
)

# NOTE: methanol
methanol_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='methanol',
        formula='CH3OH',
        state='g'
    ),
    source=methanol_thermodb_file
)

# NOTE: 1-butene
butene_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='1-butene',
        formula='C4H8',
        state='g'
    ),
    source=butene_thermodb_file
)

# NOTE: propane
C3H8_Comp = Component(
    name='propane',
    formula='C3H8',
    state='g'
)
propane_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=C3H8_Comp,
    source=propane_thermodb_file
)

# NOTE: methane
CH4_Comp = Component(
    name='methane',
    formula='CH4',
    state='g'
)
methane_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CH4_Comp,
    source=methane_thermodb_file
)

# NOTE: component thermodb source
_component_thermodb = [
    CO2_component_thermodb,
    acetylene_component_thermodb,
    n_butane_component_thermodb,
    ethanol_component_thermodb,
    methanol_component_thermodb,
    butene_component_thermodb,
    propane_component_thermodb,
    methane_component_thermodb,
]
# =======================================
# SECTION: BUILD THERMODB MODEL SOURCE
# =======================================
# update thermodb rule
thermodb_rules: Dict[str, Dict[str, ComponentRule]] = {
    'ALL': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa'
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    },
    'CH4-g': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa'
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    },
    'Methane-g': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa'
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    }
}

model_source_: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=_component_thermodb,
    rules=thermodb_rules,
)
print(model_source_)

# get data source and equation source
datasource = model_source_.data_source
equationsource = model_source_.equation_source

# ------------------------------------------------
# ! THERMODYNAMIC PROPERTIES
# ------------------------------------------------
# vapor pressure
VaPr = equationsource['propane-g']['VaPr'].cal(T=300.1)
print(VaPr)

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
phase = "LIQUID"
# temperature [K]
T = 300.1
# pressure [bar]
P = 10

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

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}


# ------------------------------------------------
# NOTE: eos root analysis
# ------------------------------------------------
# ! old method
res = eos.check_eos_roots_single_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source
)
print(res)

# ! new method
res = check_component_eos_roots(
    component=C3H8_Comp,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source_,
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
    model_source=model_source,
    liquid_fugacity_mode='Poynting'
)
print(res)

# ! new method
res = eos.calc_fugacity(
    component=C3H8_Comp,
    model_input=model_input,
    model_source=model_source_,
    liquid_fugacity_mode='Poynting'
)
print(res)

# ------------------------------------------------
# NOTE: calculation
# ------------------------------------------------
# ! new method
res = calc_liquid_fugacity(
    component=C3H8_Comp,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source_,
    liquid_fugacity_mode='Poynting'
)
print(res)
