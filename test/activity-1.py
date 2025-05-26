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
components = ['ethanol', 'butyl-methyl-ether']

# model input
activity_model = 'NRTL'

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# activity model
activity = ptm.activity(components=components, model_name=activity_model)
print(activity)

# select nrtl
activity_nrtl = activity.nrtl
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

# dg_ij
dg_ij = non_randomness_parameters.ijs(
    f"dg | {components[0]} | {components[1]}")
print(type(dg_ij))
print(dg_ij)

# alpha_ij
alpha_ij = non_randomness_parameters.ijs(
    f"alpha | {components[0]} | {components[1]}")
print(type(alpha_ij))
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
    "tau_ij": tau_ij,
    "alpha_ij": alpha_ij
}

# NOTE: calculate activity
res_, others_ = activity_nrtl.cal(model_input=model_input)
# print(res_)

# print the results
print(f"res_: {res_}")
print("-" * 50)
# activity coefficients
activity_coefficients = others_['AcCo_i_comp']
print(f"activity coefficients: {activity_coefficients}")
print("-" * 50)
G_ij = others_['G_ij']
print(f"G_ij: {G_ij}")
print("-" * 50)

# SECTION: calculate excess gibbs free energy
# NOTE: excess gibbs free energy
gibbs_energy = activity_nrtl.excess_gibbs_free_energy(
    mole_fraction=mole_fraction, G_ij=G_ij, tau_ij=tau_ij)
print(f"excess gibbs free energy method 1: {gibbs_energy}")
print("-" * 50)
gibbs_energy = activity_nrtl.excess_gibbs_free_energy()
print(f"excess gibbs free energy method 1: {gibbs_energy}")
print("-" * 50)

# NOTE: general excess gibbs free energy
gibbs_energy = activity.general_excess_molar_gibbs_free_energy(
    mole_fraction=mole_fraction, activity_coefficients=activity_coefficients)
print(f"general excess gibbs free energy: {gibbs_energy}")
print("-" * 50)
