# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL, UNIQUAC
import pyThermoDB as ptdb
from pyThermoDB.core import TableMatrixData
import pyThermoLinkDB as ptdblink
from pyThermoModels.core import calc_activity_coefficient_using_nrtl_model
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

# ! nrtl ethanol-butyl-methyl-ether
nrtl_path = os.path.join(
    thermodb_dir,
    'thermodb_nrtl_ethanol_butyl-methyl-ether.pkl'
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
components = ['ethanol', 'butyl-methyl-ether']

# model input
activity_model = 'NRTL'

# SECTION: model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# activity model
activity = ptm.activity(
    components=components,
    model_name=activity_model
)
print(activity)
print(ptm.activity.metadata)

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

# NOTE: components
ethanol = Component(
    name='ethanol',
    formula='C2H6O',
    state='l',
    mole_fraction=0.4
)

butyl_methyl_ether = Component(
    name='butyl-methyl-ether',
    formula='C5H12O',
    state='l',
    mole_fraction=0.6
)

# NOTE: operating conditions
# temperature
temperature = Temperature(value=323.15, unit='K')
# pressure
pressure = Pressure(value=30, unit='bar')

# feed spec
mole_fraction = {
    'ethanol': 0.4,
    'butyl-methyl-ether': 0.6
}

# NOTE: non-randomness parameters
non_randomness_parameters = thermodb_nrtl_1.select('non-randomness-parameters')
# >> check
if not isinstance(non_randomness_parameters, TableMatrixData):
    raise ValueError("non_randomness_parameters is not TableMatrixData")

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

# alpha_ij_comp
alpha_ij_comp = non_randomness_parameters.ijs(
    f"alpha | {components[0]} | {components[1]}",
    res_format='alphabetic')
print(type(alpha_ij_comp))
print(alpha_ij_comp)
if not isinstance(alpha_ij_comp, dict):
    raise ValueError("alpha_ij_comp is not dict")

# NOTE: operating conditions
# temperature [K]
T = 323.15
# pressure [bar]
P = 30

# NOTE: calculate the interaction parameter matrix (tau_ij)
tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(
    temperature=T,
    dg_ij=dg_ij
)
print(f"tau_ij: {tau_ij}")
print(f"tau_ij_comp: {tau_ij_comp}")


# SECTION: model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "alpha_ij": alpha_ij
}

# NOTE: calculate activity
res_, others_ = activity_nrtl.cal(
    model_input=model_input
)
print("res_:")
print(res_)

# ! new method
res_2, others_2, Gx_2 = calc_activity_coefficient_using_nrtl_model(
    components=[ethanol, butyl_methyl_ether],
    pressure=pressure,
    temperature=temperature,
    tau_ij=tau_ij_comp,
    alpha_ij=alpha_ij_comp
)
print("res_2:")
print(res_2)
print("Gx_2:")
print(Gx_2)
