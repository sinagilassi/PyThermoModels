# import packages/modules
import os
from typing import Any, Dict, List, Optional, Tuple
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL
import pyThermoDB as ptdb
from pyThermoDB.core import TableMatrixData
import pyThermoLinkDB as ptdblink
from pyThermoModels.core import calc_activity_coefficient
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoLinkDB.models import ModelSource
from pyThermoDB import (
    build_mixture_thermodb_from_reference,
    MixtureThermoDB
)

# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# ====================================
# ☑️ CUSTOM REFERENCES
# ====================================
REFERENCE_CONTENT = """
REFERENCES:
    NRTL:
      DATABOOK-ID: 1
      TABLES:
        Non-randomness parameters of the NRTL equation:
          TABLE-ID: 1
          DESCRIPTION:
            This table provides the NRTL non-randomness parameters for the NRTL equation.
          MATRIX-SYMBOL:
            - alpha constant: alpha
            - binary interaction parameter: dg
          STRUCTURE:
            COLUMNS: [No.,Mixture,Name,Formula,State,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
            SYMBOL: [None,None,None,None,None,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
            UNIT: [None,None,None,None,None,1,1,1,1]
          VALUES:
            - [1,ethanol|butyl-methyl-ether,ethanol,C2H5OH,l,0,0.680715,0,3268.884433]
            - [2,ethanol|butyl-methyl-ether,butyl-methyl-ether,C5H12O,l,0.680715,0,1768.662389,0]
"""
# =======================================
# ☑️ DIRECTORY SETUP
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, '..', 'thermodb')
print(thermodb_dir)

# ====================================
# ☑️ COMPONENTS
# ====================================
# NOTE: components
ethanol = Component(
    name='ethanol',
    formula='C2H5OH',
    state='l',
    mole_fraction=0.4
)

butyl_methyl_ether = Component(
    name='butyl-methyl-ether',
    formula='C5H12O',
    state='l',
    mole_fraction=0.6
)

# ! binary mixture
binary_mixture = [ethanol, butyl_methyl_ether]

# ====================================
# SECTION: BUILD THERMODB
# ====================================
# NOTE: normal build
# ! >> binary mixture
thermodb_components_: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=binary_mixture,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir
)
print(f"thermodb_components_: {type(thermodb_components_)}")

# >> thermodb
if thermodb_components_ is None:
    raise ValueError("thermodb_components_ is None!")

# ! nrtl ethanol-butyl-methyl-ether
thermodb_nrtl_1 = thermodb_components_.thermodb
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

# NOTE: operating conditions
# temperature
temperature = Temperature(value=323.15, unit='K')
# pressure
pressure = Pressure(value=30, unit='bar')

# =======================================
# SECTION: CALCULATION
# =======================================
# NOTE: calculate activity
res_, others_ = calc_activity_coefficient(
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
