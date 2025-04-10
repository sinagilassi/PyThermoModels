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

# ! ethanol
# thermodb file name
EtOH_thermodb_file = os.path.join(thermodb_dir, 'ethanol-1.pkl')
# load
EtOH_thermodb = ptdb.load_thermodb(EtOH_thermodb_file)
print(type(EtOH_thermodb))
print(EtOH_thermodb.check())

# ! methanol
# thermodb file name
MeOH_thermodb_file = os.path.join(thermodb_dir, 'methanol-1.pkl')
# load
MeOH_thermodb = ptdb.load_thermodb(MeOH_thermodb_file)
print(type(MeOH_thermodb))
print(MeOH_thermodb.check())

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
thub1.add_thermodb('EtOH', EtOH_thermodb)
thub1.add_thermodb('MeOH', MeOH_thermodb)

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
# model input
# eos model
eos_model = 'SRK'

# component phase
phase = "VAPOR"

# feed spec
N0s = {
    'EtOH': 0.75,
    'MeOH': 0.25
}

# temperature [K]
T = 300

# pressure [Pa]
P = 1.2*1e5

# model input
model_input = {
    "phase": phase,
    "feed-spec": N0s,
    "operating-conditions": {
        "pressure": [P, 'Pa'],
        "temperature": [T, 'K'],
    },
    "datasource": datasource,
    "equationsource": equationsource,
}

# =======================================
# EOS ROOT ANALYSIS
# =======================================
# eos root analysis
res_ = tm.check_eos_roots(model_name=eos_model, model_input=model_input)
print(type(res_))
print(res_)

# =======================================
# CHECK REFERENCES
# =======================================
# check reference
res_ = tm.check_fugacity_reference(eos_model)
print(type(res_))
print(res_)

# =======================================
# FUGACITY CALCULATION
# =======================================
# method 2
res = tm.cal_fugacity(model_name=eos_model, model_input=model_input,
                      root_analysis_set=1, liquid_fugacity_mode='EOS')
print(type(res))
print(res)
