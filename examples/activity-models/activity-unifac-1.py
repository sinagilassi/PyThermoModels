# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import UNIFAC
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
# local
from private.unifac_data_loader import load_unifac_data

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

# =======================================
# SECTION: THERMODB LINK CONFIGURATION
# =======================================
# ! no need
# build datasource & equationsource
datasource, equationsource = {}, {}

# =======================================
# NOTE: INITIALIZE ACTIVITY
# =======================================
# SECTION: configure activity model
# components
components = ['acetone', 'n_heptane']

# model input
activity_model = 'UNIFAC'

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

# select unifac
activity_unifac = activity.unifac
print(activity_unifac)

# =======================================
# SECTION: Load UNIFAC parameters
# =======================================
# load unifac data
group_parameters_dict, group_interaction_dict = load_unifac_data()
# log
print("UNIFAC group parameters and interaction parameters loaded.")
print(len(group_parameters_dict), "group parameters loaded.")
print(len(group_interaction_dict), "group interaction parameters loaded.")

# NOTE: set unifac data
activity_unifac.load_data(
    group_data=group_parameters_dict,
    interaction_data=group_interaction_dict
)

print("UNIFAC data set successfully.")
# =======================================
# NOTE COMPONENT GROUP ASSIGNMENT
# =======================================
# Component Definitions: List of {Subgroup ID: Count}
acetone_structure = {"1": 1.0, "18": 1.0}

# Component 2:
n_heptane_structure = {"1": 2.0, "2": 3.0}

components_group_data = {
    "acetone": acetone_structure,
    "n_heptane": n_heptane_structure
}

# set component group data
activity_unifac.set_component_groups(
    component_groups=components_group_data
)
# ========================================
# NOTE ACTIVITY CALCULATION
# ========================================
# NOTE: Example: Calculate activity coefficients using UNIFAC model for acetone and n-heptane mixture at 307 K

# feed spec
mole_fraction = {
    'acetone': 0.047,
    'n_heptane': 0.953
}


# SECTION: model input
model_inputs = {
    "mole_fraction": mole_fraction,
    "temperature": [307, "K"]
}

# NOTE: calculate activity
res_, other_values = activity_unifac.cal(
    model_inputs=model_inputs
)
print(res_)
print(other_values)


# print the results
activity_coefficients = res_['value']

# NOTE: general excess gibbs free energy
gibbs_energy = activity.general_excess_molar_gibbs_free_energy(
    mole_fraction=mole_fraction,
    activity_coefficients=activity_coefficients
)
print(f"general excess gibbs free energy: {gibbs_energy}")
print("-" * 50)
