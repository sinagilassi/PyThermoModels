import os

import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
import pyThermoModels as ptm
from pyThermoDB import MixtureThermoDB, build_mixture_thermodb_from_reference
from pyThermoLinkDB import build_mixture_model_source, build_model_source
from pyThermoLinkDB.models import ModelSource, MixtureModelSource
from pythermodb_settings.models import Component
from rich import print


print(ptm.__version__)
print(ptdb.__version__)
print(ptdblink.__version__)


REFERENCE_CONTENT = """
REFERENCES:
    ENRTL:
      DATABOOK-ID: 1
      TABLES:
        Chen-Evans ENRTL local-composition parameters:
          TABLE-ID: 1
          DESCRIPTION:
            Example ENRTL local-composition parameters for aqueous NaCl.
            These values are illustrative for exercising the API, not a
            validated literature parameter set.
          MATRIX-SYMBOL:
            - non-randomness parameter: alpha
            - local interaction parameter: tau
          STRUCTURE:
            COLUMNS: [No.,Mixture,Name,Formula,State,alpha_i_1,alpha_i_2,tau_i_1,tau_i_2]
            SYMBOL: [None,None,None,None,None,alpha_i_1,alpha_i_2,tau_i_1,tau_i_2]
            UNIT: [None,None,None,None,None,1,1,1,1]
          VALUES:
            - [1,sodium|water,sodium,"Na{+}",aq,0.0,0.2,0.00,0.15]
            - [2,sodium|water,water,H2O,l,0.2,0.0,0.10,0.00]
            - [3,chloride|water,chloride,"Cl{-}",aq,0.0,0.2,0.00,0.25]
            - [4,chloride|water,water,H2O,l,0.2,0.0,0.20,0.00]
            - [5,chloride|sodium,chloride,"Cl{-}",aq,0.0,0.2,0.00,0.40]
            - [6,chloride|sodium,sodium,"Na{+}",aq,0.2,0.0,0.30,0.00]
"""


parent_dir = os.path.dirname(os.path.abspath(__file__))
thermodb_dir = os.path.join(parent_dir, "..", "thermodb")


water = Component(
    name="water",
    formula="H2O",
    state="l",
    mole_fraction=0.98,
)
sodium = Component(
    name="sodium",
    formula="Na{+}",
    state="aq",
    mole_fraction=0.01,
)
chloride = Component(
    name="chloride",
    formula="Cl{-}",
    state="aq",
    mole_fraction=0.01,
)

true_species = [water, sodium, chloride]


thermodb_components: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=true_species,
    reference_content=REFERENCE_CONTENT,
    check_labels=False,
    ignore_state_all_props=True,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)
if thermodb_components is None:
    raise ValueError("thermodb_components is None")

mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=thermodb_components,
)
model_source: ModelSource = build_model_source(
    source=[mixture_model_source],
)
model_source_dict = {
    "datasource": model_source.data_source,
    "equationsource": model_source.equation_source,
}


base_model_input = {
    "composition_representation": "true_species",
    "temperature": [298.15, "K"],
    "mole_fraction": {
        "H2O-l": 0.98,
        "Na{+}-aq": 0.01,
        "Cl{-}-aq": 0.01,
    },
    "molality": {
        "H2O-l": 0.0,
        "Na{+}-aq": 0.1,
        "Cl{-}-aq": 0.1,
    },
    "long_range": {
        "model": "pitzer_debye_huckel",
        "basis": "molality",
        "A_phi": 0.392,
    },
}


activity_from_source = ptm.activity(
    components=true_species,
    model_name="ENRTL",
    model_source=model_source_dict,
)
enrtl_from_source = activity_from_source.enrtl

res_source, details_source = enrtl_from_source.cal(
    model_input=base_model_input,
    tau_correlation="direct_tau",
)

print("\nENRTL from inline YAML model source")
print(res_source)
print("tau_ij:")
print(details_source["tau_ij"])
print("alpha_ij:")
print(details_source["alpha_ij"])
print("local diagnostics:")
print(details_source["local_composition_diagnostics"])


activity_from_inputs = ptm.activity(
    components=true_species,
    model_name="ENRTL",
)
enrtl_from_inputs = activity_from_inputs.enrtl

model_input_with_parameters = {
    **base_model_input,
    "tau_ij": [
        [0.00, 0.10, 0.20],
        [0.15, 0.00, 0.30],
        [0.25, 0.40, 0.00],
    ],
    "alpha_ij": [
        [0.0, 0.2, 0.2],
        [0.2, 0.0, 0.2],
        [0.2, 0.2, 0.0],
    ],
}

res_inputs, details_inputs = enrtl_from_inputs.cal(
    model_input=model_input_with_parameters,
)
g_ex_inputs = enrtl_from_inputs.excess_gibbs_free_energy(
    mole_fraction=model_input_with_parameters["mole_fraction"],
    ln_gamma=details_inputs["ln_gamma_total"],
)

print("\nENRTL from direct model_input parameters")
print(res_inputs)
print("mean ionic activity coefficient:")
print(
    enrtl_from_inputs.mean_ionic_activity_coefficient(
        gamma=details_inputs["activity_coefficients_comp"],
        cation="Na{+}-aq",
        anion="Cl{-}-aq",
    )
)
print("G_ex / RT:")
print(g_ex_inputs)
