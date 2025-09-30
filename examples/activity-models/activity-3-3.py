# import packages/modules
import os
import numpy as np
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
from pyThermoModels.activity import UNIQUAC
import pyThermoLinkDB as ptdblink
import numpy as np

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

# ! nrtl ethanol-butyl-methyl-ether
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

# mole fraction
mole_fraction = {
    'methanol': 0.4,
    'ethanol': 0.6
}

# molecular id
# MeOH-EtOH
r_i = [1.4311, 2.1055]
q_i = [1.4320, 1.8920]

# NOTE: non-randomness parameters
# binary energy of interaction parameters
tau_ij = [
    [1, 1.031995],
    [1.309036, 1]
]

# ! >> to numpy array
r_i = np.array(r_i)
q_i = np.array(q_i)
tau_ij = np.array(tau_ij)

# NOTE: operating conditions
# temperature [K]
T = 350.1546257
# pressure [bar]
P = 30

# SECTION: model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "r_i": r_i,
    "q_i": q_i
}

# ! in case tau_ij, and r_i, q_i are not provided but existing in the datasource
model_input = {
    "mole_fraction": mole_fraction,
    "temperature": [T, 'K']
}

# NOTE: calculate activity
res_, others_ = activity_uniquac.cal(model_input=model_input)
# print(res_)

# print the results
print(f"res_: {res_}")
print("-" * 50)
print(f"others_: {others_}")
print("-" * 50)


# NOTE: excess gibbs free energy
gibbs_energy = activity_uniquac.excess_gibbs_free_energy(
    mole_fraction=mole_fraction,
    tau_ij=tau_ij,
    r_i=r_i,
    q_i=q_i
)
print(f"excess gibbs free energy 1: {gibbs_energy}")
print("-" * 50)
gibbs_energy = activity_uniquac.excess_gibbs_free_energy()
print(f"excess gibbs free energy 2: {gibbs_energy}")
print("-" * 50)
