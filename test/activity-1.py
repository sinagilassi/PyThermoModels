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
thub1.add_thermodb('nrtl', thermodb_nrtl_1)

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
# SECTION: configure activity model
# model input
activity_model = 'NRTL'

# NOTE: Example: Ethanol-Butyl-Methyl-Ether
# feed spec
N0s = {
    'ethanol': 0.4,
    'butyl-methyl-ether': 0.6
}
# # temperature [K]
T = 373.15
# # pressure [bar]
P = 30

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
# NOTE CHECK REFERENCES
# =======================================
# check reference
# res_ = tm.check_activity_reference(activity_model)
# print(res_)

# =======================================
# NOTE ACTIVITY CALCULATION
# =======================================
# calculate fugacity
res = tm.cal_activity(model_name=activity_model,
                      model_input=model_input,
                      model_source=model_source)
print(res)
