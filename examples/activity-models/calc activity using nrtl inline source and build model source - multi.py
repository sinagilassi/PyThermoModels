# import packages/modules
import os

from rich import print

import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import pyThermoModels as ptm
from pyThermoDB import MixtureThermoDB, build_mixture_thermodb_from_reference
from pyThermoLinkDB import build_mixture_model_source, build_model_source
from pyThermoLinkDB.models import MixtureModelSource, ModelSource
from pyThermoModels.core import calc_activity_coefficient
from pythermodb_settings.models import Component, Pressure, Temperature


# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# ====================================
# CUSTOM REFERENCES
# ====================================
REFERENCE_CONTENT = """
REFERENCES:
    NRTL:
      DATABOOK-ID: 2
      TABLES:
        NRTL parameters for methanol-ethanol-butyl-methyl-ether:
          TABLE-ID: 1
          DESCRIPTION:
            This table provides ternary NRTL non-randomness and interaction
            energy data encoded as binary pair rows.
          MATRIX-SYMBOL:
            - alpha constant: alpha
            - binary interaction parameter: dg
          STRUCTURE:
            COLUMNS: [No.,Mixture,Name,Formula,State,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
            SYMBOL: [None,None,None,None,None,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
            UNIT: [None,None,None,None,None,1,1,J/mol,J/mol]
          VALUES:
            - [1,methanol|ethanol,methanol,CH3OH,l,0,0.30,0,900]
            - [2,methanol|ethanol,ethanol,C2H5OH,l,0.30,0,850,0]
            - [1,methanol|butyl-methyl-ether,methanol,CH3OH,l,0,0.31,0,750]
            - [2,methanol|butyl-methyl-ether,butyl-methyl-ether,C5H12O,l,0.31,0,650,0]
            - [1,ethanol|butyl-methyl-ether,ethanol,C2H5OH,l,0,0.28,0,700]
            - [2,ethanol|butyl-methyl-ether,butyl-methyl-ether,C5H12O,l,0.28,0,500,0]
"""

# =======================================
# DIRECTORY SETUP
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, '..', 'thermodb')
print(thermodb_dir)

# ====================================
# COMPONENTS
# ====================================
methanol = Component(
    name='methanol',
    formula='CH3OH',
    state='l',
    mole_fraction=0.30,
)

ethanol = Component(
    name='ethanol',
    formula='C2H5OH',
    state='l',
    mole_fraction=0.45,
)

butyl_methyl_ether = Component(
    name='butyl-methyl-ether',
    formula='C5H12O',
    state='l',
    mole_fraction=0.25,
)

# ! multi-component mixture
multi_component_mixture = [methanol, ethanol, butyl_methyl_ether]

# ====================================
# SECTION: BUILD THERMODB
# ====================================
# NOTE: build mixture thermodb from inline NRTL source
thermodb_components_: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=multi_component_mixture,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
print(f"thermodb_components_: {type(thermodb_components_)}")

# >> thermodb
if thermodb_components_ is None:
    raise ValueError("thermodb_components_ is None!")

thermodb_nrtl_1 = thermodb_components_.thermodb

# check
print(thermodb_nrtl_1.check())

# ====================================
# BUILD MIXTURE MODEL SOURCE
# ====================================
# NOTE: build mixture model source
mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=thermodb_components_,
)
print(f"mixture_model_source: {mixture_model_source}")

# SECTION: build model source
model_source: ModelSource = build_model_source(
    source=[mixture_model_source],
)
print(f"model_source: {model_source}")

# =======================================
# SECTION: INITIALIZE INPUTS
# =======================================
# NOTE: operating conditions
temperature = Temperature(value=323.15, unit='K')
pressure = Pressure(value=30, unit='bar')

# =======================================
# SECTION: CALCULATION
# =======================================
# NOTE: calculate activity
res_, others_, G_ex = calc_activity_coefficient(
    components=multi_component_mixture,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name='NRTL',
    verbose=True,
)

# print the results
print("res_:")
print(res_)
print("-" * 50)
print("others_:")
print(others_)
print("-" * 50)
print("G_ex:")
print(G_ex)
