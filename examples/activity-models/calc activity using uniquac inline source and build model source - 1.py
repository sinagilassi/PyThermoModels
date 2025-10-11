# import packages/modules
import os
from typing import Any, Dict, List, Optional, Tuple
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL
import pyThermoDB as ptdb
from pyThermoDB.core import TableMatrixData
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB import (
    build_model_source,
    build_mixture_model_source,
    build_component_model_source,
    build_components_model_source
)
from pyThermoModels.core import calc_activity_coefficient
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoLinkDB.models import ModelSource, MixtureModelSource, ComponentModelSource
from pyThermoDB import (
    build_mixture_thermodb_from_reference,
    build_component_thermodb_from_reference,
    MixtureThermoDB,
    ComponentThermoDB
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
    CUSTOM-REFERENCE-1:
      DATABOOK-ID: 1
      TABLES:
        Non-randomness parameters of the UNIQUAC equation:
          TABLE-ID: 1
          DESCRIPTION:
            This table provides the UNIQUAC non-randomness parameters for the UNIQUAC equation.
          MATRIX-SYMBOL:
            - a constant: a
            - b constant: b
            - c constant: c
            - d constant: d
            - alpha constant: alpha
          STRUCTURE:
            COLUMNS: [No.,Mixture,Name,Formula,State,a_i_1,a_i_2,b_i_1,b_i_2,c_i_1,c_i_2,d_i_1,d_i_2,alpha_i_1,alpha_i_2]
            SYMBOL: [None,None,None,None,None,a_i_1,a_i_2,b_i_1,b_i_2,c_i_1,c_i_2,d_i_1,d_i_2,alpha_i_1,alpha_i_2]
            UNIT: [None,None,None,None,None,1,1,1,1,1,1,1,1,1,1]
          VALUES:
            - [1,methanol|ethanol,methanol,CH3OH,l,0,0.300492719,0,1.564200272,0,35.05450323,0,0,0,4.481683583]
            - [2,methanol|ethanol,ethanol,C2H5OH,l,0.380229054,0,-20.63243601,0,0.059982839,0,0,0,4.481683583,0]
        general-data:
          TABLE-ID: 2
          DESCRIPTION:
            This table provides the general data of different chemical species participating in the CO2 hydrogenation reaction and includes molecular weight (MW) in g/mol, critical temperature (Tc) in K, critical pressure (Pc) in MPa, and critical molar volume (Vc) in m3/kmol. The table also includes the critical compressibility factor (Zc), acentric factor (AcFa), enthalpy of formation (EnFo) in kJ/mol, and Gibbs energy of formation (GiEnFo) in kJ/mol. Moreover,
            The chemical state of the species is also provided in the table and hence the enthalpy of formation and Gibbs energy of formation are provided for the ideal gas and liquid state are designated as EnFo_IG, GiEnFo_IG, EnFo_LIQ, and GiEnFo_LIQ, respectively.
          DATA: []
          STRUCTURE:
            COLUMNS: [No.,Name,Formula,State,Molecular-Weight,Critical-Temperature,Critical-Pressure,Critical-Molar-Volume,Critical-Compressibility-Factor,Acentric-Factor,Enthalpy-of-Formation,Gibbs-Energy-of-Formation,Volume-Parameter,Surface-Area-Parameter]
            SYMBOL: [None,None,None,None,MW,Tc,Pc,Vc,Zc,AcFa,EnFo,GiEnFo,r,q]
            UNIT: [None,None,None,None,g/mol,K,MPa,m3/kmol,None,None,kJ/mol,kJ/mol,None,None]
            CONVERSION: [None,None,None,None,1,1,1,1,1,1,1,1,1,1]
          VALUES:
            - [1,'methanol','CH3OH','l',32.04,512.5,8.084,0.117,0.222,0.5658,-200.7,-162,1.4311,1.4320]
            - [2,'ethanol','C2H5OH','l',46.068,514,6.137,0.168,0.241,0.6436,-277.70,-174.80,2.1055,1.8920]
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
    mole_fraction=0.6
)

# methanol
methanol = Component(
    name='methanol',
    formula='CH3OH',
    state='l',
    mole_fraction=0.4
)

# ! binary mixture
binary_mixture = [ethanol, methanol]

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

# NOTE: components in thermodb
# ! methanol
methanol_thermodb: ComponentThermoDB | None = build_component_thermodb_from_reference(
    component_name=methanol.name,
    component_formula=methanol.formula,
    component_state=methanol.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir
)
# check
if methanol_thermodb is None:
    raise ValueError("methanol_thermodb is None!")
# thermodb
print(methanol_thermodb.thermodb.check())

# ! ethanol
ethanol_thermodb: ComponentThermoDB | None = build_component_thermodb_from_reference(
    component_name=ethanol.name,
    component_formula=ethanol.formula,
    component_state=ethanol.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir
)
# check
if ethanol_thermodb is None:
    raise ValueError("ethanol_thermodb is None!")
# thermodb
print(ethanol_thermodb.thermodb.check())

# ====================================
# ☑️ BUILD MIXTURE MODEL SOURCE
# ====================================
# NOTE: build mixture model source
mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=thermodb_components_,
)
print(f"mixture_model_source: {mixture_model_source}")

# NOTE: build component model source
components_model_source: List[ComponentModelSource] = build_components_model_source(
    components_thermodb=[methanol_thermodb, ethanol_thermodb]
)
print(f"components_model_source: {components_model_source}")

# SECTION: build model source
# all model source
source: list = []
source.extend(components_model_source)
source.append(mixture_model_source)

# build model source
model_source: ModelSource = build_model_source(
    source=source,
)
print(f"model_source: {model_source}")

# =======================================
# SECTION: INITIALIZE INPUTS
# =======================================
# NOTE: operating conditions
# temperature
temperature = Temperature(value=323.15, unit='K')
# pressure
pressure = Pressure(value=30, unit='bar')

# =======================================
# SECTION: CALCULATION
# =======================================
# NOTE: calculate activity
res_, others_, G_ex = calc_activity_coefficient(
    components=[ethanol, methanol],
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name='UNIQUAC',
    verbose=True
)
# print the results
print("res_:")
print(res_)
print("-" * 50)
print("G_ex:")
print(G_ex)
