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


# ! nrtl ethanol-butyl-methyl-ether
nrtl_path = os.path.join(
    thermodb_dir, 'thermodb_nrtl_ethanol_butyl-methyl-ether_1.pkl')
# load
thermodb_nrtl_1 = ptdb.load_thermodb(nrtl_path)
# check
print(thermodb_nrtl_1.check())

# ========================================
# ! INITIALIZE PYTHERMOMODELS
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

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# ! CALCULATE ACTIVITY
# =======================================
# NOTE: Reference

# model input
# eos model
eos_model = 'NRTL'

# NOTE: Example 9.2. (page 252), The Thermodynamics of Phase and Reaction Equilibria
# feed spec
N0s = {
    'ethanol': 0.4,
    'butyl-methyl-ether': 0.6
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
    'CO2': 0.15,
    'n-butane': 0.85,
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

# =======================================
# ! CHECK REFERENCES
# =======================================
# check reference
# res_ = tm.check_fugacity_reference(eos_model)
# print(res_)

# =======================================
# ! FUGACITY CALCULATION
# =======================================
# calculate fugacity
res = tm.cal_activity(model_name=eos_model,
                      model_input=model_input,
                      model_source=model_source)
print(res)
