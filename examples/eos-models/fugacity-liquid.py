# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# check version
print(ptm.__version__)
# check version
print(ptdb.__version__)
# check version
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
print(type(CO2_thermodb))
print(CO2_thermodb.check())

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(thermodb_dir, 'acetylene-g.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)
print(type(acetylene_thermodb))
print(acetylene_thermodb.check())

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-g.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)
print(type(n_butane_thermodb))
print(n_butane_thermodb.check())

# ! ethanol
# thermodb file name
ethanol_thermodb_file = os.path.join(thermodb_dir, 'ethanol-l.pkl')
# load
ethanol_thermodb = ptdb.load_thermodb(ethanol_thermodb_file)
print(type(ethanol_thermodb))
print(ethanol_thermodb.check())

# ! methanol
# thermodb file name
methanol_thermodb_file = os.path.join(thermodb_dir, 'methanol-g.pkl')
# load
methanol_thermodb = ptdb.load_thermodb(methanol_thermodb_file)
print(type(methanol_thermodb))
print(methanol_thermodb.check())

# ! propane
# thermodb file name
propane_thermodb_file = os.path.join(thermodb_dir, 'propane-g.pkl')
# load
propane_thermodb_file = ptdb.load_thermodb(propane_thermodb_file)
print(type(propane_thermodb_file))
print(propane_thermodb_file.check())

# ========================================
# ! INITIALIZE OBJECT
# ========================================
tm = ptm.init()
# log
print("tm: ", tm)

# NOTE: check reference
eos_model = 'PR'
res_ = tm.check_fugacity_reference(eos_model)
print(res_)

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
# thub1.add_thermodb('CO2', CO2_thermodb)
# thub1.add_thermodb('EtOH', ethanol_thermodb)
# thub1.add_thermodb('MeOH', methanol_thermodb)
# acetylene
thub1.add_thermodb('acetylene', acetylene_thermodb)
# propane
thub1.add_thermodb('propane', propane_thermodb_file)

# * add thermodb rule
thermodb_config_file = os.path.join(
    parent_dir,
    'thermodb_config_link.yml'
)

# one component
# thub1.config_thermodb_rule(thermodb_config_file, name='EtOH')
# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# init
eos = ptm.eos()
print("eos: ", eos)

# model input
# eos model
eos_model = 'PR'

# component phase
phase = "LIQUID"

# NOTE: example
# component
component = 'acetylene'
# temperature [K]
T = 250
# pressure [bar]
P = 20

# NOTE: example
# component
component = 'propane'
# temperature [K]
T = 300.1
# pressure [bar]
P = 10

# SECTION: model input content
model_input_content = """
phase: LIQUID
component: propane
pressure: [10, 'bar']
temperature: [300.1, 'K']
"""

# parse model input content
model_inputs_parsed = eos.parse_model_inputs(model_input_content)

# SECTION: model input
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
res = eos.check_eos_roots_single_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source)
print(res)

# ------------------------------------------------
# NOTE: eos
# ------------------------------------------------
# NOTE: liquid fugacity calculation method
res = eos.cal_fugacity(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source,
    liquid_fugacity_mode='Poynting')
print(res)
