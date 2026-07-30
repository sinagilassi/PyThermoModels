# import packages/modules
import os
import time
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from typing import Callable, Dict, Optional, Union, List, Any
from pyThermoLinkDB import (
    build_component_model_source,
    build_components_model_source,
    build_model_source
)
from pyThermoLinkDB.models import ComponentModelSource, ModelSource
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoDB import ComponentThermoDB
from pyThermoDB import build_component_thermodb_from_reference
from pyThermoModels.core import calc_liquid_fugacity
# locals
from examples.references.reference_3 import REFERENCE_CONTENT

# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# ====================================
# BUILD COMPONENT THERMODB
# ====================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, 'thermodb')
print(thermodb_dir)

# NOTE: create component
# ! propane
# component info
# =======================================
# SECTION: COMPONENTS THERMODB SOURCE
# =======================================
# NOTE: carbon dioxide
CO2 = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g'
)

# NOTE: acetylene
C2H2 = Component(
    name='acetylene',
    formula='C2H2',
    state='g'
)

# NOTE: n-butane
C4H10 = Component(
    name='n-butane',
    formula='C4H10',
    state='g'
)

# NOTE: ethanol
C2H5OH = Component(
    name='ethanol',
    formula='C2H5OH',
    state='l'
)

# NOTE: methanol
CH3OH = Component(
    name='methanol',
    formula='CH3OH',
    state='g'
)

# NOTE: 1-butene
C4H8 = Component(
    name='1-butene',
    formula='C4H8',
    state='g'
)

# NOTE: propane
C3H8 = Component(
    name='propane',
    formula='C3H8',
    state='g'
)

# NOTE: methane
CH4 = Component(
    name='methane',
    formula='CH4',
    state='g'
)

components = [CO2, C3H8, CH3OH, C4H10]


# NOTE: ignore state properties
ignore_state_props = ['MW', 'VaPr']

# ====================================================
# SECTION: build components thermodb
# ====================================================
thermodb_components: List[ComponentThermoDB] = []

_thermodb_t0 = time.perf_counter()
for comp in components:
    thermodb_component = build_component_thermodb_from_reference(
        component_name=comp.name,
        component_formula=comp.formula,
        component_state=comp.state,
        reference_content=REFERENCE_CONTENT,
        ignore_state_props=ignore_state_props,
        thermodb_save=True,
        thermodb_save_path=thermodb_dir,
    )
    if thermodb_component is None:
        raise ValueError(f"thermodb_component for {comp.name} is None")
    thermodb_components.append(thermodb_component)
_thermodb_t1 = time.perf_counter()
print(
    f"[bold cyan]Timing[/bold cyan] build_component_thermodb_from_reference (total): "
    f"{(_thermodb_t1 - _thermodb_t0) * 1000.0:.2f} ms"
)

# ====================================================
# SECTION: build model source
# ====================================================
# NOTE: with partially matched rules
_build_t0 = time.perf_counter()
component_model_source: List[ComponentModelSource] = build_components_model_source(
    components_thermodb=thermodb_components,
    rules=None,
)

# model source
model_source: ModelSource = build_model_source(
    source=component_model_source,
)
_build_t1 = time.perf_counter()
print(
    f"[bold cyan]Timing[/bold cyan] build_components_model_source + build_model_source: "
    f"{(_build_t1 - _build_t0) * 1000.0:.2f} ms"
)
# ====================================================
# SECTION: THERMODB LINK CONFIGURATION
# ====================================================

# build datasource & equationsource
datasource = model_source.data_source
equationsource = model_source.equation_source

# dict
model_source_dict = {
    "datasource": datasource,
    "equationsource": equationsource
}
print(model_source_dict)
