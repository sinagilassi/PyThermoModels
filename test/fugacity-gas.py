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
# NOTE: thermodb directory
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(thermodb_dir, 'carbon dioxide-1.pkl')
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(thermodb_dir, 'acetylene-1.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-1.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)

# ! ethanol
# thermodb file name
ethanol_thermodb_file = os.path.join(thermodb_dir, 'ethanol-1.pkl')
# load
ethanol_thermodb = ptdb.load_thermodb(ethanol_thermodb_file)

# ! methanol
# thermodb file name
methanol_thermodb_file = os.path.join(thermodb_dir, 'methanol-1.pkl')
# load
methanol_thermodb = ptdb.load_thermodb(methanol_thermodb_file)

# ! 1-butene
# thermodb file name
butene_thermodb_file = os.path.join(thermodb_dir, '1-butene-1.pkl')
# load
butene_thermodb = ptdb.load_thermodb(butene_thermodb_file)

# ! propane
# thermodb file name
propane_thermodb_file = os.path.join(thermodb_dir, 'propane-1.pkl')
# load
propane_thermodb = ptdb.load_thermodb(propane_thermodb_file)

# ! methane
# thermodb file name
methane_thermodb_file = os.path.join(thermodb_dir, 'methane-1.pkl')
# load
methane_thermodb = ptdb.load_thermodb(methane_thermodb_file)

# ========================================
# ! INITIALIZE OBJECT
# ========================================
tm = ptm.init()
# log
print("tm: ", tm)

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('CO2', CO2_thermodb)
thub1.add_thermodb('CH4', methane_thermodb)
thub1.add_thermodb('EtOH', ethanol_thermodb)
thub1.add_thermodb('MeOH', methanol_thermodb)
thub1.add_thermodb('acetylene', acetylene_thermodb)
thub1.add_thermodb('1-butene', butene_thermodb)
thub1.add_thermodb('n-butane', n_butane_thermodb)
thub1.add_thermodb('propane', propane_thermodb)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')
# one component
# thub1.config_thermodb_rule(thermodb_config_file, name='EtOH')
# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# ------------------------------------------------
# ! THERMODYNAMIC PROPERTIES
# ------------------------------------------------
# vapor pressure
VaPr = equationsource['propane']['VaPr'].cal(T=300.1)
print(VaPr)

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# NOTE:
# Example 3.9 (page 100), Introduction to Chemical Engineering Thermodynamics
# by J.M. Smith, H.C. Van Ness, M.M. Abbott

# Example 10.9, Introduction to Chemical Engineering Thermodynamics

# Example 3.11 (page 68), The Thermodynamics of Phase and Reaction Equilibria


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
component = "propane"
# phase
phase = "VAPOR-LIQUID"
# temperature [K]
T = 300.1
# pressure [bar]
P = 9.99

# NOTE: Example 3.13
# component
component = "CH4"
# temperature [K]
T = 340
# pressure [bar]
P = 30


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
# NOTE: check reference
# ------------------------------------------------
# res_ = tm.check_fugacity_reference(eos_model)
# print(res_)

# ------------------------------------------------
# NOTE: eos root analysis
# ------------------------------------------------
res = tm.check_eos_roots_single_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source)
print(res)

# ------------------------------------------------
# NOTE: calculation
# ------------------------------------------------
# NOTE: gas fugacity calculation method
res = tm.cal_fugacity(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source)
print(res)
