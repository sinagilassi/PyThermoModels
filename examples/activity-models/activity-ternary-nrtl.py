# import packages/modules
import os
import runpy
import numpy as np

from rich import print

import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoModels.core import calc_activity_coefficient
from pyThermoLinkDB.models import MixtureModelSource, ModelSource
from pythermodb_settings.models import Component, Pressure, Temperature
from pythermodb_settings.utils import create_mixture_id


# check version
print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)

# ====================================
# LOAD TABLE DATA FROM REFERENCE FILE
# ====================================
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

reference_file = os.path.join(parent_dir, '..', 'references', 'reference_2.py')
reference_data = runpy.run_path(reference_file)

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

ternary_mixture = [methanol, ethanol, butyl_methyl_ether]

# matrix data stored in the reference example file
alpha_ij = reference_data['NRTL_TERNARY_ALPHA_DICT']
dg_ij = reference_data['NRTL_TERNARY_DG_DICT']

# build the exact mixture ids that the core resolves at runtime
mixture_name = create_mixture_id(
    components=ternary_mixture,
    mixture_key='Name',
    delimiter='|',
)
mixture_formula = create_mixture_id(
    components=ternary_mixture,
    mixture_key='Formula',
    delimiter='|',
)

mixture_data_source = {
    mixture_name: {
        'alpha': alpha_ij,
        'dg': dg_ij,
    },
    mixture_formula: {
        'alpha': alpha_ij,
        'dg': dg_ij,
    },
}

mixture_equation_source = {
    mixture_name: {},
    mixture_formula: {},
}

mixture_model_source = MixtureModelSource(
    components=ternary_mixture,
    data_source=mixture_data_source,
    equation_source=mixture_equation_source,
    check_labels=True,
    label_link=False,
)
print(f"mixture_model_source: {mixture_model_source}")

model_source = ModelSource(
    data_source=mixture_model_source.data_source,
    equation_source=mixture_model_source.equation_source,
)
print(f"model_source: {model_source}")

# ====================================
# INPUTS
# ====================================
temperature = Temperature(value=323.15, unit='K')
pressure = Pressure(value=30, unit='bar')

# ====================================
# CALCULATION
# ====================================
res_, others_, G_ex = calc_activity_coefficient(
    components=ternary_mixture,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name='NRTL',
    verbose=True,
)

print('res_:')
print(res_)
print('-' * 50)
print('others_:')
print(others_)
print('-' * 50)
print('G_ex:')
print(G_ex)
