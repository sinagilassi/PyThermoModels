"""Calculate Redlich-Kister activity coefficients for a binary mixture."""

from rich import print

import pyThermoModels as ptm


# SECTION: Mixture and Redlich-Kister coefficients
components = ["acetone", "chloroform"]
mole_fraction = {
    "acetone": 0.35,
    "chloroform": 0.65,
}

# NOTE: a[k] is the coefficient of (x1 - x2)^k, beginning with a[0].
redlich_kister_input = {
    "mole_fraction": mole_fraction,
    "a": [-0.82, 0.31, -0.08],
}

# ! Redlich-Kister is implemented for binary systems only.
# ? Add coefficients to the list when a higher-order expansion is required.

# SECTION: Activity-coefficient calculation
activity = ptm.activity(components=components, model_name="REDLICH_KISTER")
res, other_values = activity.redlich_kister.cal(model_input=redlich_kister_input)

# The calculated values remain aligned with the component order above.
print(res)
print({"G^E/RT": other_values["excess_gibbs_RT"]})

# SECTION: Excess Gibbs free energy
gibbs_energy = activity.redlich_kister.excess_gibbs_free_energy(
    mole_fraction=mole_fraction,
    a=redlich_kister_input["a"],
)
print(gibbs_energy)
