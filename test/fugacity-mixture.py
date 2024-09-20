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

# ! ethanol
# thermodb file name
EtOH_thermodb_file = os.path.join(os.getcwd(), 'test', 'ethanol.pkl')
# load
EtOH_thermodb = ptdb.load_thermodb(EtOH_thermodb_file)
print(type(EtOH_thermodb))

# CO2_thermodb
EtOH_thermodb

# ! methanol
# thermodb file name
MeOH_thermodb_file = os.path.join(os.getcwd(), 'test', 'methanol.pkl')
# load
MeOH_thermodb = ptdb.load_thermodb(MeOH_thermodb_file)
print(type(MeOH_thermodb))

MeOH_thermodb

# ========================================
# INITIALIZE FUGACITY OBJECT
# ========================================
fugacity_obj = ptm.fugacity()
# log
print("fugacity_obj: ", fugacity_obj)

# =======================================
# THERMODB CONFIGURATION
# =======================================

# ! EtOH
# add EtOH thermodb
fugacity_obj.add_thermodb('EtOH', EtOH_thermodb)

# ! MeOH
# add acetylene thermodb
fugacity_obj.add_thermodb('MeOH', MeOH_thermodb)


# * add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config.yml')
fugacity_obj.config_thermodb('EtOH', thermodb_config_file)
fugacity_obj.config_thermodb('MeOH', thermodb_config_file)


# check thermodb
print(fugacity_obj.check_thermodb())


# =======================================
# CALCULATE FUGACITY FOR MULTI COMPONENT
# =======================================
# model input
# eos model
eos_model = 'SRK'

# component phase
phase = "VAPOR"

# component list
comp_list = ["EtOH", "MeOH"]

# mole fraction
MoFri = [0.75, 0.25]

# temperature [K]
T = 300

# pressure [Pa]
P = 1.2*1e5

# model input
model_input = {
    "eos-model": eos_model,
    "phase": phase,
    "components": comp_list,
    "mole-fraction": MoFri,
    "operating_conditions": {
        "pressure": [P, 'Pa'],
        "temperature": [T, 'K'],
    },
}

# =======================================
# EOS ROOT ANALYSIS
# =======================================
# eos root analysis
# res = fugacity_obj.check_eos_roots(model_input)
# pp(res)

# =======================================
# CHECK REFERENCES
# =======================================
# check reference
# pp(fugacity_obj.fugacity_check_reference(eos_model))

# =======================================
# FUGACITY CALCULATION
# =======================================
# method 2
Phi, Z, _ = fugacity_obj.fugacity_cal(model_input)
pp(Z)
pp(Phi)
