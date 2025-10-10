# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL, UNIQUAC
import pyThermoDB as ptdb
from pyThermoDB.core import TableMatrixData
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
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, '..', 'thermodb')
print(thermodb_dir)

# ! nrtl ethanol-methanol
nrtl_path = os.path.join(
    thermodb_dir, 'mixture methanol-ethanol.pkl')
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
activity_model = 'NRTL'

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# activity model
activity_nrtl = ptm.activities(
    components=components,
    model_name=activity_model
)

# check
if not isinstance(activity_nrtl, NRTL):
    raise TypeError(
        f"activity_nrtl is not an instance of NRTL: {type(activity_nrtl)}"
    )

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
# NOTE: Example: Methanol-Ethanol

# mole fraction
mole_fraction = {
    'methanol': 0.0316,
    'ethanol': 0.9684
}

# NOTE: non-randomness parameters
# non-randomness parameters
non_randomness_parameters = thermodb_nrtl_1.select(
    'CUSTOM-REF-1::NRTL Non-randomness parameters-2'
)
# >> check
if not isinstance(non_randomness_parameters, TableMatrixData):
    raise ValueError("non_randomness_parameters is not TableMatrixData")

print(type(non_randomness_parameters))

# a_ij
a_ij = non_randomness_parameters.ijs(f"a | {components[0]} | {components[1]}")
print(type(a_ij))
print(a_ij)

# b_ij
b_ij = non_randomness_parameters.ijs(f"b | {components[0]} | {components[1]}")
print(type(b_ij))
print(b_ij)

# c_ij
c_ij = non_randomness_parameters.ijs(f"c | {components[0]} | {components[1]}")
print(type(c_ij))
print(c_ij)

# alpha_ij
alpha_ij = non_randomness_parameters.ijs(
    f"alpha | {components[0]} | {components[1]}")
print(type(alpha_ij))
print(alpha_ij)

# NOTE: operating conditions
# temperature [K]
T = 350.1546257
# pressure [bar]
P = 30

# NOTE: calculate the interaction parameter matrix (dg_ij)
dg_ij, dg_ij_comp = activity_nrtl.cal_dg_ij_M1(
    temperature=T,
    a_ij=a_ij,
    b_ij=b_ij,
    c_ij=c_ij
)
print(f"dg_ij: {dg_ij}")
print(f"dg_ij_comp: {dg_ij_comp}")
print("-" * 50)

# NOTE: calculate the interaction parameter matrix (tau_ij)
tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(
    temperature=T,
    dg_ij=dg_ij
)
print(f"tau_ij: {tau_ij}")
print(f"tau_ij_comp: {tau_ij_comp}")
print("-" * 50)

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
G_ij = others_['G_ij']
print(f"G_ij: {G_ij}")
print("-" * 50)

# NOTE: excess gibbs free energy
gibbs_energy = activity_nrtl.excess_gibbs_free_energy(
    mole_fraction=mole_fraction,
    G_ij=G_ij,
    tau_ij=tau_ij
)
print(f"excess gibbs free energy method 1: {gibbs_energy}")
print("-" * 50)
