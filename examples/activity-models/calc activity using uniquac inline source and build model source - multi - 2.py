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
# ! reference
from examples.references.uniquac_reference_1 import REFERENCE_CONTENT


# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)


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
# ! methanol
methanol = Component(
    name='methanol',
    formula='CH3OH',
    state='l',
    mole_fraction=0.30,
)

# ! ethanol
ethanol = Component(
    name='ethanol',
    formula='C2H5OH',
    state='l',
    mole_fraction=0.45,
)

# dimethyl-carbonate
dimethyl_carbonate = Component(
    name='dimethyl-carbonate',
    formula='C3H6O3',
    state='l',
    mole_fraction=0.25,
)

# ! multi-component mixture
multi_component_mixture = [methanol, ethanol, dimethyl_carbonate]

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

dimethyl_carbonate_thermodb = build_component_thermodb_from_reference(
    component_name=dimethyl_carbonate.name,
    component_formula=dimethyl_carbonate.formula,
    component_state=dimethyl_carbonate.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
if dimethyl_carbonate_thermodb is None:
    raise ValueError("dimethyl_carbonate_thermodb is None!")
print(dimethyl_carbonate_thermodb.thermodb.check())

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
        dimethyl_carbonate_thermodb,
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
