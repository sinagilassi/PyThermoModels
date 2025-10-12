# import libs
import os
from pythermodb_settings.models import Component, Temperature, Pressure
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL, UNIQUAC
import pyThermoDB as ptdb
from pyThermoDB.core import TableMatrixData
import pyThermoLinkDB as ptdblink
from pyThermoModels.activity import (
    calc_dU_ij_using_uniquac_model,
    calc_tau_ij_with_dU_ij_using_uniquac_model
)

# check version
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
# NOTE: components
methanol = Component(name='methanol', formula='CH3OH', state='l')
ethanol = Component(name='ethanol', formula='C2H5OH', state='l')

# components
components = [methanol.name, ethanol.name]

# ========================================
# NOTE ACTIVITY CALCULATION
# ========================================

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
a_ij = non_randomness_parameters.ijs(
    property=f"a | {components[0]} | {components[1]}",
    res_format='alphabetic'
)
if not isinstance(a_ij, dict):
    raise ValueError("a_ij is not a dictionary")
print(type(a_ij))
print(a_ij)

# b_ij
b_ij = non_randomness_parameters.ijs(
    property=f"b | {components[0]} | {components[1]}",
    res_format='alphabetic'
)
if not isinstance(b_ij, dict):
    raise ValueError("b_ij is not a dictionary")
print(type(b_ij))
print(b_ij)

# c_ij
c_ij = non_randomness_parameters.ijs(
    property=f"c | {components[0]} | {components[1]}",
    res_format='alphabetic'
)
if not isinstance(c_ij, dict):
    raise ValueError("c_ij is not a dictionary")
print(type(c_ij))
print(c_ij)


# NOTE: operating conditions
# temperature [K]
temperature = Temperature(value=350.1546257, unit='K')
# pressure [bar]
P = Pressure(value=30, unit='bar')

# NOTE: calculate the interaction parameter matrix (dU_ij)
dU_ij, dU_ij_comp, dU_ij_comp_upd = calc_dU_ij_using_uniquac_model(
    components=[methanol, ethanol],
    temperature=temperature,
    a_ij=a_ij,
    b_ij=b_ij,
    c_ij=c_ij
)
print(dU_ij)
print(dU_ij_comp)
print(dU_ij_comp_upd)

# # NOTE: calculate the interaction parameter matrix (tau_ij)
tau_ij, tau_ij_comp, tau_ij_comp_upd = calc_tau_ij_with_dU_ij_using_uniquac_model(
    components=[methanol, ethanol],
    temperature=temperature,
    dU_ij=dU_ij_comp,
)
print(tau_ij)
print(tau_ij_comp)
print(tau_ij_comp_upd)
