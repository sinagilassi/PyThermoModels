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

# ! N2
# thermodb file name
N2_thermodb_file = os.path.join(thermodb_dir, 'nitrogen-1.pkl')
# load
N2_thermodb = ptdb.load_thermodb(N2_thermodb_file)

# ! methane
# thermodb file name
CH4_thermodb_file = os.path.join(thermodb_dir, 'methane-1.pkl')
# load
CH4_thermodb = ptdb.load_thermodb(CH4_thermodb_file)

# ========================================
# ! INITIALIZE FUGACITY OBJECT
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
thub1.add_thermodb('N2', N2_thermodb)
thub1.add_thermodb('CH4', CH4_thermodb)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! CALCULATE FUGACITY FOR MULTI COMPONENT
# =======================================
# NOTE: Reference
# Example 10.7 (page 378) in Introduction to Chemical Engineering Thermodynamics
# by J.M. Smith, H.C. Van Ness, M.M. Abbott
# Example 15.2 (page 587) in Introductory Chemical Engineering Thermodynamics

# model input
# eos model
eos_model = 'PR'

# component phase
phase = "VAPOR"

# feed spec
N0s = {
    'N2': 0.50,
    'CH4': 0.50
}

# temperature [K]
T = 200
# pressure [bar]
P = 30

# # temperature [K]
# T = 100
# # pressure [MPa]
# P = 0.4119

# model input
model_input = {
    "feed-specification": N0s,
    "pressure": [P, 'MPa'],
    "temperature": [T, 'K'],
}

# model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# =======================================
# ! CHECK REFERENCES
# =======================================
# check reference
# res_ = tm.check_fugacity_reference(eos_model)
# print(res_)

# =======================================
# ! EOS ROOT ANALYSIS
# =======================================
# eos root analysis
res_ = tm.check_eos_roots_multi_component(model_name=eos_model,
                                          model_input=model_input,
                                          model_source=model_source)
print(res_)

# =======================================
# ! FUGACITY CALCULATION
# =======================================
# calculate fugacity
res = tm.cal_fugacity_mixture(model_name=eos_model,
                              model_input=model_input,
                              model_source=model_source)
print(res)
