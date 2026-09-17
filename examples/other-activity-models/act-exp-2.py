"""Compare Margules and van Laar activity models for a binary mixture."""

from rich import print

import pyThermoModels as ptm


# SECTION: Mixture definition
components = ["ethanol", "water"]
mole_fraction = {
    "ethanol": 0.40,
    "water": 0.60,
}

# SECTION: Three-suffix Margules model
# NOTE: A12 multiplies x1 and A21 multiplies x2 in G^E/RT.
margules_input = {
    "mole_fraction": mole_fraction,
    "A12": 1.35,
    "A21": 0.85,
}

# ! Margules is implemented for binary systems only.
# ? Use calculation_mode="two_suffix" with a single A value for a symmetric fit.
margules = ptm.activity(components=components, model_name="MARGULES").margules
margules_res, margules_other = margules.cal(
    model_input=margules_input,
    calculation_mode="three_suffix",
)

print("Margules")
print(margules_res)
print({"G^E/RT": margules_other["excess_gibbs_RT"]})

# The direct method produces the standard excess-Gibbs result dictionary.
margules_gibbs = margules.excess_gibbs_free_energy(
    mole_fraction=mole_fraction,
    A12=margules_input["A12"],
    A21=margules_input["A21"],
    calculation_mode="three_suffix",
)
print(margules_gibbs)

# SECTION: van Laar model
van_laar_input = {
    "mole_fraction": mole_fraction,
    "a1": 1.40,
    "a2": 0.90,
}

van_laar = ptm.activity(components=components, model_name="VAN_LAAR").van_laar
van_laar_res, van_laar_other = van_laar.cal(model_input=van_laar_input)

print("van Laar")
print(van_laar_res)
print({"G^E/RT": van_laar_other["excess_gibbs_RT"]})

van_laar_gibbs = van_laar.excess_gibbs_free_energy(
    mole_fraction=mole_fraction,
    a1=van_laar_input["a1"],
    a2=van_laar_input["a2"],
)
print(van_laar_gibbs)
