# import packages/modules
import PyThermoModels as ptm
import pyThermoDB as ptdb
import os

# check version
print(ptm.__version__)
# check version
print(ptdb.__version__)

# =======================================
# LOAD THERMODB
# =======================================
# thermodb file name
CO2_thermodb_file = os.path.join(os.getcwd(), 'test', 'Carbon Dioxide.pkl')
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)
print(type(CO2_thermodb))


# ========================================
# INITIALIZE FUGACITY OBJECT
# ========================================
fugacity_obj = ptm.fugacity()
# log
print("fugacity_obj: ", fugacity_obj)

# =======================================
# THERMODB CONFIGURATION
# =======================================
# add CO2 thermodb
fugacity_obj.add_thermodb('CO2', CO2_thermodb)

# add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config.yml')
fugacity_obj.config_thermodb('CO2', thermodb_config_file)

# check thermodb
print(fugacity_obj.check_thermodb())

# =======================================
# CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# model input
# eos model
eos_model = 'PR'

# component phase
phase = "gas"

# component list
comp_list = ["CO2"]
# required component input
# MW,Tc,Pc,w,Zc,Vc

# mole fraction
MoFri = []

# temperature [K]
T = 200 + 273.15

# pressure [Pa]
P = 15.55*1e5

# model input
model_input = {
    "eos-model": eos_model,
    "phase": phase,
    "components": comp_list,
    "mole-fraction": MoFri,
    "operating_conditions": {
        "P": P,
        "T": T,
    },
}

# eos
fugacity_obj.fugacity_cal_init(model_input)
