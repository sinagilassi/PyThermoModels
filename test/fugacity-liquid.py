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
print(type(CO2_thermodb))
print(CO2_thermodb.check())

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(thermodb_dir, 'acetylene-1.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)
print(type(acetylene_thermodb))
print(acetylene_thermodb.check())

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-1.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)
print(type(n_butane_thermodb))
print(n_butane_thermodb.check())

# ! ethanol
# thermodb file name
ethanol_thermodb_file = os.path.join(thermodb_dir, 'ethanol-1.pkl')
# load
ethanol_thermodb = ptdb.load_thermodb(ethanol_thermodb_file)
print(type(ethanol_thermodb))
print(ethanol_thermodb.check())

# ! methanol
# thermodb file name
methanol_thermodb_file = os.path.join(thermodb_dir, 'methanol-1.pkl')
# load
methanol_thermodb = ptdb.load_thermodb(methanol_thermodb_file)
print(type(methanol_thermodb))
print(methanol_thermodb.check())

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
# thub1.add_thermodb('CO2', CO2_thermodb)
# thub1.add_thermodb('EtOH', ethanol_thermodb)
# thub1.add_thermodb('MeOH', methanol_thermodb)
# acetylene
thub1.add_thermodb('acetylene', acetylene_thermodb)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')
# one component
# thub1.config_thermodb_rule(thermodb_config_file, name='EtOH')
# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# model input
# eos model
eos_model = 'SRK'

# component phase
phase = "LIQUID"

# component
N0s = {
    "acetylene": 1.0
}

# temperature [K]
T = 250

# pressure [MPa]
P = 2

# model input
model_input = {
    "phase": phase,
    "feed-specification": N0s,
    "pressure": [P, 'MPa'],
    "temperature": [T, 'K'],
}

# model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}
# ------------------------------------------------
# NOTE: check reference
# ------------------------------------------------
res_ = tm.check_fugacity_reference(eos_model)
print(res_)

# ------------------------------------------------
# NOTE: eos root analysis
# ------------------------------------------------
res = tm.check_eos_roots(
    model_name=eos_model, model_input=model_input, model_source=model_source)
print(res)

# ------------------------------------------------
# NOTE: eos
# ------------------------------------------------
# NOTE: liquid fugacity calculation method
res = tm.cal_fugacity(model_name=eos_model, model_input=model_input,
                      model_source=model_source)

# NOTE: gas fugacity calculation method
# res = tm.cal_fugacity(
#     model_name=eos_model, model_input=model_input)

Z, Phi, eos_parms, phi_parms = res
# res
print(f"Z: {Z}")
print('-'*50)
print(f"Phi: {Phi}")
print('-'*50)
# print(eos_parms)
# print('-'*50)
# print(phi_parms)
# print('-'*50)
