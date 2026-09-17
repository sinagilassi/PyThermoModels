"""Pitzer v1 example: NaCl at 298.15 K and 1 mol/kg."""
from rich import print
from pyThermoModels.activity.pitzer import Pitzer


# SECTION: Direct parameters for the documented 298.15 K NaCl regression.
model = Pitzer(["Na+", "Cl-"])
result, details = model.cal(
    {
        "temperature": [298.15, "K"],
        "salt_molality": 1.0,
        "charges": [1, -1],
        "stoichiometry": [1, 1],
        "beta0": 0.0765,
        "beta1": 0.2664,
        "c_phi": 0.00127,
        "alpha": 2.0,
        "A_phi": 0.3915,
        "b": 1.2,
    }
)

# NOTE: Pitzer v1 reports only gamma_pm, never individual-ion coefficients.
print(result)
print(details)
