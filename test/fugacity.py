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

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(os.getcwd(), 'test', 'Carbon Dioxide.pkl')
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)
print(type(CO2_thermodb))
print(CO2_thermodb.check())

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(os.getcwd(), 'test', 'acetylene.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)
print(type(acetylene_thermodb))
print(acetylene_thermodb.check())

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(os.getcwd(), 'test', 'n-butane.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)
print(type(n_butane_thermodb))
print(n_butane_thermodb.check())

# ========================================
# ! INITIALIZE OBJECT
# ========================================
tm = ptm.init()
# log
print("tm: ", tm)

# =======================================
# SECTION: THERMODB CONFIGURATION
# =======================================
# # ! CO2
# # add CO2 thermodb
# tm.add_thermodb('CO2', CO2_thermodb)

# # ! acetylene
# # add acetylene thermodb
# tm.add_thermodb('acetylene', acetylene_thermodb)

# # ! n-butane
# # add n-butane thermodb
# tm.add_thermodb('n-butane', n_butane_thermodb)


# # * add thermodb rule
# thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config.yml')
# tm.config_thermodb('acetylene', thermodb_config_file)
# tm.config_thermodb('CO2', thermodb_config_file)

# # check thermodb
# print(tm.check_thermodb())

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
# thub1.add_thermodb('EtOH', EtOH_thermodb)
# thub1.add_thermodb('MeOH', MeOH_thermodb)
thub1.add_thermodb('CO2', CO2_thermodb)

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
phase = "VAPOR"

# component list
# comp_list = ["CO2"]

# mole fraction
# MoFri = []

# component
N0s = {
    "CO2": 1.0
}

# temperature [K]
T = 350

# pressure [Pa]
P = 9.4573*1e5

# model input
model_input = {
    "eos-model": eos_model,
    "phase": phase,
    "feed-spec": N0s,
    "operating-conditions": {
        "pressure": [P, 'Pa'],
        "temperature": [T, 'K'],
    },
    "datasource": datasource,
    "equationsource": equationsource
}

# ------------------------------------------------
# check reference
# ------------------------------------------------
# print(tm.fugacity_check_reference(eos_model))

# ------------------------------------------------
# eos
# ------------------------------------------------
# method 1
# res1, res2 = tm.fugacity_cal_init(model_input)
# print(res1)
# print('-'*50)
# print(res2)


# method 2
# liquid fugacity calculation method
# Z, Phi, eos_parms = tm.cal_fugacity_coefficient(
#     model_input, root_analysis_set=1, liquid_fugacity_calculation_method='EOS')

# gas fugacity calculation method
Z, Phi, eos_parms = tm.cal_fugacity_coefficient(model_input)

# res
print(Z)
print('-'*50)
print(Phi)
print('-'*50)
print(eos_parms)
print('-'*50)

# ------------------------------------------------
# eos root analysis
# ------------------------------------------------
# res = tm.check_eos_roots(model_input)
# print(res)

# ------------------------------------------------
# thermo lib
# ------------------------------------------------
# t_lib = ptm.thermo_lib()

# calculate molar volume
# Vm = t_lib.cal_molar_volume(P, T, Z[0])*1e6
# print(Vm)
