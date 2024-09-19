# import packages/modules
from pprint import pprint as pp
import pyThermoModels as ptm
import pyThermoDB as ptdb
import os

# check version
print(ptm.__version__)
# check version
print(ptdb.__version__)

# =======================================
# LOAD THERMODB
# =======================================

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(os.getcwd(), 'test', 'Carbon Dioxide.pkl')
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)
print(type(CO2_thermodb))

# CO2_thermodb
CO2_thermodb

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(os.getcwd(), 'test', 'acetylene.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)
print(type(acetylene_thermodb))

acetylene_thermodb

# ========================================
# INITIALIZE FUGACITY OBJECT
# ========================================
fugacity_obj = ptm.fugacity()
# log
print("fugacity_obj: ", fugacity_obj)

# =======================================
# THERMODB CONFIGURATION
# =======================================

# ! CO2
# add CO2 thermodb
fugacity_obj.add_thermodb('CO2', CO2_thermodb)

# ! acetylene
# add acetylene thermodb
fugacity_obj.add_thermodb('acetylene', acetylene_thermodb)


# * add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config.yml')
fugacity_obj.config_thermodb('acetylene', thermodb_config_file)
fugacity_obj.config_thermodb('CO2', thermodb_config_file)


# check thermodb
print(fugacity_obj.check_thermodb())


# =======================================
# CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# model input
# eos model
eos_model = 'SRK'

# component phase
phase = "GAS"

# component list
comp_list = ["CO2"]
# required component input
# MW,Tc,Pc,w,Zc,Vc

# mole fraction
MoFri = []

# temperature [K]
T = 310

# pressure [Pa]
P = 8*1e5

# model input
model_input = {
    "eos-model": eos_model,
    "phase": phase,
    "components": comp_list,
    "mole-fraction": MoFri,
    "operating_conditions": {
        "pressure": P,
        "temperature": T,
    },
}

# check reference
# pp(fugacity_obj.fugacity_check_reference(eos_model))
# eos
# res = fugacity_obj.fugacity_cal_init(model_input)
# pp(res)

# method 2
res = fugacity_obj.fugacity_cal(model_input)
pp(res)
