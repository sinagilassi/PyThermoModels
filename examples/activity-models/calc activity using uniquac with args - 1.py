# import packages/modules
import os
import numpy as np
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
from pyThermoModels.activity import UNIQUAC
import pyThermoLinkDB as ptdblink
import numpy as np
from pyThermoModels.core import calc_activity_coefficient_using_uniquac_model
from pythermodb_settings.models import Component, Pressure, Temperature

# version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# ! LOAD THERMODB
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, '..', 'thermodb')
print(thermodb_dir)

# ! nrtl methanol-ethanol
nrtl_path = os.path.join(
    thermodb_dir,
    'mixture methanol-ethanol.pkl'
)
# load
thermodb_nrtl_1 = ptdb.load_thermodb(nrtl_path)
# check
print(thermodb_nrtl_1.check())

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('uniquac', thermodb_nrtl_1)

# * add thermodb rule
thermodb_config_file = os.path.join(
    parent_dir,
    'thermodb_config_link.yml'
)

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()

# =======================================
# NOTE: INITIALIZE ACTIVITY
# =======================================
# SECTION: configure activity model
# components
components = ['methanol', 'ethanol']

# model input
activity_model = 'UNIQUAC'

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# activity model
activity_uniquac = ptm.activities(
    components=components,
    model_name=activity_model
)
# check uniquac model
if not isinstance(activity_uniquac, UNIQUAC):
    raise ValueError("Activity model is not UNIQUAC")
print(activity_uniquac)

# =======================================
# NOTE CHECK REFERENCES
# =======================================
# check reference
# res_ = tm.check_activity_reference(activity_model)
# print(res_)

# ========================================
# NOTE ACTIVITY CALCULATION
# ========================================
# NOTE: Example: Methanol-Ethanol

# NOTE: components
ethanol = Component(
    name='ethanol',
    formula='C2H5OH',
    state='l',
    mole_fraction=0.6
)

# methanol
methanol = Component(
    name='methanol',
    formula='CH3OH',
    state='l',
    mole_fraction=0.4
)


# mole fraction
mole_fraction = {
    'methanol': 0.4,
    'ethanol': 0.6
}

# molecular id
# MeOH-EtOH
r_i = [1.4311, 2.1055]
q_i = [1.4320, 1.8920]

# as dict
r_i_comp = {
    'methanol': 1.4311,
    'ethanol': 2.1055
}

q_i_comp = {
    'methanol': 1.4320,
    'ethanol': 1.8920
}

# NOTE: non-randomness parameters
# binary energy of interaction parameters
tau_ij = [
    [1, 1.031995],
    [1.309036, 1]
]

# tau_ij as dict
tau_ij_comp = {
    'methanol | methanol': 1,
    'methanol | ethanol': 1.031995,
    'ethanol | methanol': 1.309036,
    'ethanol | ethanol': 1
}

# ! > to numpy array>
r_i = np.array(r_i)
q_i = np.array(q_i)
tau_ij = np.array(tau_ij)

# NOTE: operating conditions
# temperature [K]
T = 350.1546257
# pressure [bar]
P = 30

# NOTE: operating conditions
# temperature
temperature = Temperature(value=323.15, unit='K')
# pressure
pressure = Pressure(value=30, unit='bar')


# SECTION: model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "r_i": r_i,
    "q_i": q_i
}

# ! in case tau_ij, and r_i, q_i are not provided but existing in the datasource
# model_input = {
#     "mole_fraction": mole_fraction,
#     "temperature": [T, 'K'],
# }

# NOTE: calculate activity
res_, others_ = activity_uniquac.cal(model_input=model_input)
# print(res_)

# print the results
print("res_:")
print(res_)
print("-" * 50)
print("others_:")
print(others_)
print("-" * 50)

# ! new method
res_2, others_2, Gx_2 = calc_activity_coefficient_using_uniquac_model(
    components=[ethanol, methanol],
    pressure=pressure,
    temperature=temperature,
    tau_ij=tau_ij_comp,
    r_i=r_i_comp,
    q_i=q_i_comp
)
print("res_2:")
print(res_2)
print("-" * 50)
print("others_2:")
print(others_2)
print("-" * 50)
print("Gx_2:")
print(Gx_2)
