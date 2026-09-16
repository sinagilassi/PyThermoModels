# import packages/modules
from pyThermoModels.eos.lee_kesler import lee_kesler_ploecker_mixture
from pathlib import Path
import sys
from rich import print

# SECTION: Local package bootstrap
# NOTE: Allows running this file directly from the examples folder.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# SECTION: Lee-Kesler-Ploecker mixture example
# NOTE: Methane/ethane critical properties are supplied explicitly in SI units.
components = ["methane", "ethane"]
x = [0.25, 0.75]
P = 1.0e5  # Pa
T = 330.0  # K
Tc = [190.564, 305.32]  # K
Pc = [4.5992e6, 4.872e6]  # Pa
acentric_factors = [0.011, 0.099]
Zc = [0.286, 0.279]

# SECTION: Binary interaction parameters
# ? This example leaves k_ij unspecified, so the model records the zero-default assumption.
k_ij = None

# SECTION: Model calculation
result = lee_kesler_ploecker_mixture(
    P=P,
    T=T,
    components=components,
    x=x,
    Tc=Tc,
    Pc=Pc,
    acentric_factors=acentric_factors,
    Zc=Zc,
    k_ij=k_ij,
    phase="vapor",
)

# SECTION: Result summary
lk_result = result.result
print(f"components: {components}")
print(f"composition: {result.composition}")
print(
    f"pseudo-critical temperature [K]: {result.pseudo_critical_temperature:.8f}")
print(f"pseudo-critical pressure [Pa]: {result.pseudo_critical_pressure:.8f}")
print(f"pseudo-critical volume [m3/mol]: {result.pseudo_critical_volume:.8e}")
print(f"mixture acentric factor: {result.mixture_acentric_factor:.8f}")
print(f"Z: {lk_result.compressibility_factor:.8f}")
print(f"molar volume [m3/mol]: {lk_result.molar_volume:.8e}")
print(f"enthalpy departure [J/mol]: {lk_result.enthalpy_departure:.8f}")
print(f"entropy departure [J/(mol.K)]: {lk_result.entropy_departure:.8f}")
print(f"gibbs departure [J/mol]: {lk_result.gibbs_departure:.8f}")

# ! Metadata makes default assumptions visible.
print(f"k_ij defaulted to zero: {result.metadata['k_ij_defaulted_to_zero']}")
