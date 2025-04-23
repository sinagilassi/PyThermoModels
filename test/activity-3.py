# import packages/modules
import os
import numpy as np
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
    thermodb_dir, 'thermodb_nrtl_methanol_ethanol_1.pkl')
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
thub1.add_thermodb('nrtl', thermodb_nrtl_1)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')

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
activity = ptm.activity(components=components, model_name=activity_model)
print(activity)

# select uniquac model
activity_uniquac = activity.uniquac
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
r_i = np.array([1.4311, 2.1055])
q_i = np.array([1.4320, 1.8920])

# NOTE: non-randomness parameters
# binary energy of interaction parameters
tau_ij = np.array([
    [1, 1.031995],
    [1.309036, 1]
])

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
    mole_fraction=mole_fraction, tau_ij=tau_ij, r_i=r_i, q_i=q_i)
print(f"excess gibbs free energy: {gibbs_energy}")
print("-" * 50)
