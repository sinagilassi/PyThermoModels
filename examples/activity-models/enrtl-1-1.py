import numpy as np
from pythermodb_settings.models import Component

import pyThermoModels as ptm


components = [
    Component(name="water", formula="H2O", state="l"),
    Component(name="sodium", formula="Na{+}", state="aq"),
    Component(name="chloride", formula="Cl{-}", state="aq"),
]

activity = ptm.activity(
    components=components,
    model_name="ENRTL",
)
enrtl = activity.enrtl

model_input = {
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
    "tau_ij": np.zeros((3, 3)),
    "alpha_ij": np.zeros((3, 3)),
    "long_range": {
        "model": "pitzer_debye_huckel",
        "basis": "molality",
        "A_phi": 0.392,
    },
}

try:
    result, details = enrtl.cal(model_input=model_input)
    print(result)
    print(details)
except NotImplementedError as exc:
    print("Chen-Evans 1986 local-composition implementation is incomplete.")
    print(exc)
