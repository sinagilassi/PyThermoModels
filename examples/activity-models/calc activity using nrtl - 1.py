# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL
import pyThermoDB as ptdb
from pyThermoDB.core import TableMatrixData
import pyThermoLinkDB as ptdblink
from pyThermoModels.core import calc_activity_coefficient
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoLinkDB.models import ModelSource

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

# ! nrtl ethanol-butyl-methyl-ether
nrtl_path = os.path.join(
    thermodb_dir, 'thermodb_nrtl_ethanol_butyl-methyl-ether.pkl')
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
    parent_dir,
    'thermodb_config_link.yml'
)

# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()


# =======================================
# SECTION: INITIALIZE INPUTS
# =======================================
# NOTE" create model source
model_source = ModelSource(
    data_source=datasource,
    equation_source=equationsource
)

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

# =======================================
# SECTION: CALCULATION
# =======================================
# NOTE: calculate activity
res_, others_, _ = calc_activity_coefficient(
    components=[ethanol, butyl_methyl_ether],
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name='NRTL',
    verbose=True
)
# print(res_)

# print the results
print(f"res_: {res_}")
print("-" * 50)
G_ij = others_['G_ij']
print(f"G_ij: {G_ij}")
print("-" * 50)

# # NOTE: excess gibbs free energy
# gibbs_energy = activity_nrtl.excess_gibbs_free_energy(
#     mole_fraction=mole_fraction, G_ij=G_ij, tau_ij=tau_ij)
# print(f"excess gibbs free energy method 1: {gibbs_energy}")
# print("-" * 50)
# gibbs_energy = activity_nrtl.excess_gibbs_free_energy()
# print(f"excess gibbs free energy method 1: {gibbs_energy}")
# print("-" * 50)
