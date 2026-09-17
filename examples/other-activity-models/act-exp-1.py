"""Calculate Wilson activity coefficients for a ternary liquid mixture."""

from rich import print
import numpy as np

import pyThermoModels as ptm


# SECTION: Mixture definition
components = ["acetone", "methanol", "water"]
mole_fraction = {
    "acetone": 0.25,
    "methanol": 0.45,
    "water": 0.30,
}

# SECTION: Wilson interaction parameters
# NOTE: Lambda_ij is dimensionless; the diagonal is fixed to one by the model.
lambda_ij = np.array([
    [1.0, 0.92, 1.18],
    [1.11, 1.0, 0.76],
    [0.84, 1.29, 1.0],
])

# ! Every supplied Wilson Lambda value must be positive.
# ? Change the composition or Lambda values to explore another liquid mixture.
model_input = {
    "mole_fraction": mole_fraction,
    "lambda_ij": lambda_ij,
}

# SECTION: Activity-coefficient calculation
activity = ptm.activity(components=components, model_name="WILSON")
res, other_values = activity.wilson.cal(model_input=model_input)

# Activity coefficients are returned in the standard PyThermoModels result format.
print(res)
print({"G^E/RT": other_values["excess_gibbs_RT"]})

# SECTION: Excess Gibbs free energy
gibbs_energy = activity.wilson.excess_gibbs_free_energy(
    mole_fraction=mole_fraction,
    lambda_ij=lambda_ij,
)
print(gibbs_energy)
