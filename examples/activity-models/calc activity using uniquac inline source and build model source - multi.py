# import packages/modules
import os
from typing import List

from rich import print

import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import pyThermoModels as ptm
from pyThermoDB import (
    build_component_thermodb_from_reference,
    build_mixture_thermodb_from_reference,
)
from pyThermoLinkDB import (
    build_components_model_source,
    build_mixture_model_source,
    build_model_source,
)
from pyThermoLinkDB.models import (
    ComponentModelSource,
    MixtureModelSource,
    ModelSource,
)
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
    CUSTOM-REFERENCE-1:
      DATABOOK-ID: 1
      TABLES:
        UNIQUAC parameters for methanol-ethanol-butyl-methyl-ether:
          TABLE-ID: 1
          DESCRIPTION:
            This table provides ternary UNIQUAC interaction data encoded as
            binary pair rows.
          MATRIX-SYMBOL:
            - a constant: a
            - b constant: b
            - c constant: c
            - d constant: d
          STRUCTURE:
            COLUMNS: [No.,Mixture,Name,Formula,State,a_i_1,a_i_2,b_i_1,b_i_2,c_i_1,c_i_2,d_i_1,d_i_2]
            SYMBOL: [None,None,None,None,None,a_i_1,a_i_2,b_i_1,b_i_2,c_i_1,c_i_2,d_i_1,d_i_2]
            UNIT: [None,None,None,None,None,1,1,1,1,1,1,1,1]
          VALUES:
            - [1,methanol|ethanol,methanol,CH3OH,l,0,0.300492719,0,1.564200272,0,35.05450323,0,0]
            - [2,methanol|ethanol,ethanol,C2H5OH,l,0.380229054,0,-20.63243601,0,0.059982839,0,0,0]
            - [1,methanol|butyl-methyl-ether,methanol,CH3OH,l,0,0.1201,0,2.25,0,18.4,0,0]
            - [2,methanol|butyl-methyl-ether,butyl-methyl-ether,C5H12O,l,0.2152,0,-8.75,0,0.041,0,0,0]
            - [1,ethanol|butyl-methyl-ether,ethanol,C2H5OH,l,0,0.1803,0,3.268,0,22.6,0,0]
            - [2,ethanol|butyl-methyl-ether,butyl-methyl-ether,C5H12O,l,0.2457,0,-12.48,0,0.052,0,0,0]
        general-data:
          TABLE-ID: 2
          DESCRIPTION:
            This table provides pure component UNIQUAC volume and surface-area
            parameters.
          DATA: []
          STRUCTURE:
            COLUMNS: [No.,Name,Formula,State,Molecular-Weight,Critical-Temperature,Critical-Pressure,Critical-Molar-Volume,Critical-Compressibility-Factor,Acentric-Factor,Enthalpy-of-Formation,Gibbs-Energy-of-Formation,Volume-Parameter,Surface-Area-Parameter]
            SYMBOL: [None,None,None,None,MW,Tc,Pc,Vc,Zc,AcFa,EnFo,GiEnFo,r,q]
            UNIT: [None,None,None,None,g/mol,K,MPa,m3/kmol,None,None,kJ/mol,kJ/mol,None,None]
            CONVERSION: [None,None,None,None,1,1,1,1,1,1,1,1,1,1]
          VALUES:
            - [1,'methanol','CH3OH','l',32.04,512.5,8.084,0.117,0.222,0.5658,-200.7,-162,1.4311,1.4320]
            - [2,'ethanol','C2H5OH','l',46.068,514,6.137,0.168,0.241,0.6436,-277.70,-174.80,2.1055,1.8920]
            - [3,'butyl-methyl-ether','C5H12O','l',88.15,497.15,3.43,0.329,0.273,0.2662,-283.5,-117.5,4.0672,3.4920]
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
# NOTE: build mixture thermodb from inline UNIQUAC source
thermodb_components_ = build_mixture_thermodb_from_reference(
    components=multi_component_mixture,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
print(f"thermodb_components_: {type(thermodb_components_)}")

# >> thermodb
if thermodb_components_ is None:
    raise ValueError("thermodb_components_ is None!")

thermodb_uniquac_1 = thermodb_components_.thermodb

# check
print(thermodb_uniquac_1.check())

# NOTE: component thermodbs provide r and q for UNIQUAC
methanol_thermodb = build_component_thermodb_from_reference(
    component_name=methanol.name,
    component_formula=methanol.formula,
    component_state=methanol.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
if methanol_thermodb is None:
    raise ValueError("methanol_thermodb is None!")
print(methanol_thermodb.thermodb.check())

ethanol_thermodb = build_component_thermodb_from_reference(
    component_name=ethanol.name,
    component_formula=ethanol.formula,
    component_state=ethanol.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
if ethanol_thermodb is None:
    raise ValueError("ethanol_thermodb is None!")
print(ethanol_thermodb.thermodb.check())

butyl_methyl_ether_thermodb = build_component_thermodb_from_reference(
    component_name=butyl_methyl_ether.name,
    component_formula=butyl_methyl_ether.formula,
    component_state=butyl_methyl_ether.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
if butyl_methyl_ether_thermodb is None:
    raise ValueError("butyl_methyl_ether_thermodb is None!")
print(butyl_methyl_ether_thermodb.thermodb.check())

# ====================================
# BUILD MIXTURE MODEL SOURCE
# ====================================
# NOTE: build mixture model source
mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=thermodb_components_,
)
print(f"mixture_model_source: {mixture_model_source}")

# NOTE: build component model source
components_model_source: List[ComponentModelSource] = build_components_model_source(
    components_thermodb=[
        methanol_thermodb,
        ethanol_thermodb,
        butyl_methyl_ether_thermodb,
    ],
)
print(f"components_model_source: {components_model_source}")

# SECTION: build model source
source: list = []
source.extend(components_model_source)
source.append(mixture_model_source)

model_source: ModelSource = build_model_source(
    source=source,
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
    model_name='UNIQUAC',
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
