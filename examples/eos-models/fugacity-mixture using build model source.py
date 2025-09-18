# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure

# version
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
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)

# check
print(CO2_thermodb.check())

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(thermodb_dir, 'acetylene-g.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-g.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)

# ! ethanol
# thermodb file name
ethanol_thermodb_file = os.path.join(thermodb_dir, 'ethanol-l.pkl')
# load
ethanol_thermodb = ptdb.load_thermodb(ethanol_thermodb_file)

# ! methanol
# thermodb file name
methanol_thermodb_file = os.path.join(thermodb_dir, 'methanol-g.pkl')
# load
methanol_thermodb = ptdb.load_thermodb(methanol_thermodb_file)

# ! 1-butene
# thermodb file name
butene_thermodb_file = os.path.join(thermodb_dir, '1-butene-g.pkl')
# load
butene_thermodb = ptdb.load_thermodb(butene_thermodb_file)

# ! propane
# thermodb file name
propane_thermodb_file = os.path.join(thermodb_dir, 'propane-g.pkl')
# load
propane_thermodb = ptdb.load_thermodb(propane_thermodb_file)

# ! methane
# thermodb file name
methane_thermodb_file = os.path.join(thermodb_dir, 'methane-g.pkl')
# load
CH4_thermodb = ptdb.load_thermodb(methane_thermodb_file)

# ! N2
# thermodb file name
N2_thermodb_file = os.path.join(thermodb_dir, 'nitrogen-g.pkl')
# load
N2_thermodb = ptdb.load_thermodb(N2_thermodb_file)

# ! C2H6
# thermodb file name
C2H6_thermodb_file = os.path.join(thermodb_dir, 'ethane-g.pkl')
# load
C2H6_thermodb = ptdb.load_thermodb(C2H6_thermodb_file)

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

# ========================================
# NOTE: INITIALIZE PYTHERMOMODELS
# ========================================
tm = ptm.init()
# log
print("tm: ", tm)

# NOTE: CHECK REFERENCES
eos_model = 'PR'  # example: 'PR', 'SRK', 'RK'
# check reference
res_ = tm.check_fugacity_reference(eos_model)
print(res_)

# =======================================
# ! CALCULATE FUGACITY FOR MULTI COMPONENT
# =======================================
# NOTE: Reference
# Example 10.7 (page 378) in Introduction to Chemical Engineering Thermodynamics
# by J.M. Smith, H.C. Van Ness, M.M. Abbott
# Example 15.2 (page 587) in Introductory Chemical Engineering Thermodynamics

# init
eos = ptm.eos()
print("eos: ", eos)

# model input
# eos model
eos_model = 'PR'

# component phase
phase = "VAPOR"

# NOTE: example 1
# feed spec
N0s = {
    'N2-g': 0.40,
    'CH4-g': 0.60,
}
# temperature [K]
T = 200
# pressure [bar]
P = 30

# NOTE: Example 9.2. (page 252), The Thermodynamics of Phase and Reaction Equilibria
# feed spec
N0s = {
    'CH4-g': 0.35,
    'C2H6-g': 0.65,
}
# # temperature [K]
T = 373.15
# # pressure [bar]
P = 30

# NOTE: Example 9.5 (page 257), The Thermodynamics of Phase and Reaction Equilibria
# eos model
eos_model = 'RK'
# feed spec
N0s = {
    'CO2-g': 0.15,
    'n-butane-g': 0.85,
}
# temperature [K]
T = 444
# pressure [bar]
P = 10
# binary interaction parameter
k_ij = [[0, 0.18],
        [0.18, 0]]

# SECTION: model input
model_input = {
    "feed-specification": N0s,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# SECTION: components
# ! CO2
CO2_Comp = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g',
    mole_fraction=N0s.get('CO2-g', 0)
)
# ! n-butane
C4H10_Comp = Component(
    name='n-butane',
    formula='C4H10',
    state='g',
    mole_fraction=N0s.get('n-butane-g', 0)
)
# ! components
components = [
    CO2_Comp, C4H10_Comp
]

# =======================================
# EOS ROOT ANALYSIS
# =======================================
# eos root analysis
res_ = eos.check_eos_roots_multi_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source
)
print(res_)

# =======================================
# FUGACITY CALCULATION
# =======================================
# calculate fugacity
# ! old method
res = eos.cal_fugacity_mixture(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source,
)
print(res)

# ! new method using build model source
res = eos.calc_fugacity_mixture(
    components=components,
    model_input=model_input,
    model_source=model_source_
)
print(res)
