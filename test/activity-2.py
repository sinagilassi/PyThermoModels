# import packages/modules
from pprint import pprint as pp
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import os

# check version
print(ptm.__version__)
# check version
print(ptdb.__version__)

# =======================================
# ! LOAD THERMODB
# =======================================

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(os.getcwd(), 'test', 'Carbon Dioxide.pkl')
# load
CO2_thermodb = ptdb.load_thermodb(CO2_thermodb_file)
print(type(CO2_thermodb))

# CO2_thermodb
pp(CO2_thermodb.check())

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(os.getcwd(), 'test', 'acetylene.pkl')
# load
acetylene_thermodb = ptdb.load_thermodb(acetylene_thermodb_file)
print(type(acetylene_thermodb))

acetylene_thermodb

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(os.getcwd(), 'test', 'n-butane.pkl')
# load
n_butane_thermodb = ptdb.load_thermodb(n_butane_thermodb_file)
print(type(n_butane_thermodb))

n_butane_thermodb

# ========================================
# ! INITIALIZE ACTIVITY OBJECT
# ========================================
activity_obj = ptm.activity_lib()
# log
print("activity_obj: ", activity_obj)

# ========================================
# ! CHECK ACTIVITY REFERENCE
# ========================================


# =======================================
# ! THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.thermodb_hub()
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
thub1.config_thermodb_rule(thermodb_config_file, name='ALL')

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! CALCULATE FUGACITY FOR PURE COMPONENT
# =======================================
# model input
# activity model
activity_model = 'Margules'

# component phase
phase = "LIQUID"

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
    "eos-model": activity_model,
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
# pp(fugacity_obj.fugacity_check_reference(eos_model))

# ------------------------------------------------
# eos
# ------------------------------------------------
# method 1
# res1, res2 = fugacity_obj.fugacity_cal_init(model_input)
# pp(res1)
# print('-'*50)
# pp(res2)


# method 2
# liquid fugacity calculation method
# Z, Phi, eos_parms = fugacity_obj.cal_fugacity_coefficient(
#     model_input, root_analysis_set=1, liquid_fugacity_calculation_method='EOS')

# gas fugacity calculation method
res = fugacity_obj.cal_fugacity_coefficient(
    model_input)
pp(res)
# res
# pp(Z)
# print('-'*50)
# pp(Phi)
# print('-'*50)
# pp(eos_parms)
# print('-'*50)

# ------------------------------------------------
# eos root analysis
# ------------------------------------------------
# res = fugacity_obj.check_eos_roots(model_input)
# pp(res)

# ------------------------------------------------
# thermo lib
# ------------------------------------------------
# t_lib = ptm.thermo_lib()

# calculate molar volume
# Vm = t_lib.cal_molar_volume(P, T, Z[0])*1e6
# pp(Vm)
