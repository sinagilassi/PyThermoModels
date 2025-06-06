# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL, UNIQUAC
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

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# init thermodb hub
thub1 = ptdblink.init()
print(type(thub1))

# add component thermodb
thub1.add_thermodb('NRTL', thermodb_nrtl_1)

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
components = ['ethanol', 'butyl-methyl-ether']

# model input
activity_model = 'NRTL'

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# activity model
activity_nrtl = ptm.activities(
    components=components,
    model_name=activity_model,
    model_source=model_source)
# check activity model
if not isinstance(activity_nrtl, NRTL):
    raise TypeError(
        f"activity_nrtl is not NRTL class, but {type(activity_nrtl)}")

print(activity_nrtl)


# =======================================
# NOTE CHECK REFERENCES
# =======================================
# check reference
# res_ = tm.check_activity_reference(activity_model)
# print(res_)

# ========================================
# NOTE ACTIVITY CALCULATION
# ========================================
# NOTE: Example: Ethanol-Butyl-Methyl-Ether

# feed spec
mole_fraction = {
    'ethanol': 0.4,
    'butyl-methyl-ether': 0.6
}

# NOTE: non-randomness parameters
non_randomness_parameters = thermodb_nrtl_1.select('non-randomness-parameters')
print(type(non_randomness_parameters))
# symbols
print(non_randomness_parameters.matrix_symbol)

# dg_ij
dg_ij = non_randomness_parameters.ijs(
    f"dg | {components[0]} | {components[1]}")
print(type(dg_ij))
print(dg_ij)

# or
dg_ij = non_randomness_parameters.mat('dg', components)
print(dg_ij)


# alpha_ij
alpha_ij = non_randomness_parameters.ijs(
    f"alpha | {components[0]} | {components[1]}")
print(type(alpha_ij))
print(alpha_ij)

# or
alpha_ij = non_randomness_parameters.mat('alpha', components)
print(alpha_ij)


# NOTE: operating conditions
# temperature [K]
T = 323.15
# pressure [bar]
P = 30

# NOTE: calculate the interaction parameter matrix (tau_ij)
tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(temperature=T, dg_ij=dg_ij)
print(f"tau_ij: {tau_ij}")
print(f"tau_ij_comp: {tau_ij_comp}")

# SECTION: model input
model_input = {
    "mole_fraction": mole_fraction,
    "temperature": [T, 'K'],
}

# NOTE: calculate activity
res_, others_ = activity_nrtl.cal(model_input=model_input)
# print(res_)

# print the results
print(f"res_: {res_}")
print("-" * 50)
G_ij = others_['G_ij']
print(f"G_ij: {G_ij}")
print("-" * 50)

# NOTE: excess gibbs free energy
gibbs_energy = activity_nrtl.excess_gibbs_free_energy(
    mole_fraction=mole_fraction, G_ij=G_ij, tau_ij=tau_ij)
print(f"excess gibbs free energy method 1: {gibbs_energy}")
print("-" * 50)
gibbs_energy = activity_nrtl.excess_gibbs_free_energy()
print(f"excess gibbs free energy method 1: {gibbs_energy}")
print("-" * 50)
