# import packages/modules
from pyThermoModels.eos.lee_kesler import lee_kesler_pure
from pathlib import Path
import sys
from rich import print

# SECTION: Local package bootstrap
# NOTE: Allows running this file directly from the examples folder.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# SECTION: Pure-fluid Lee-Kesler example
# NOTE: Methane critical properties are provided explicitly in SI units.
component = "methane"
P = 1.0e5  # Pa
T = 300.0  # K
Tc = 190.564  # K
Pc = 4.5992e6  # Pa
acentric_factor = 0.011

# SECTION: Model calculation
# ? phase="vapor" selects the largest physical reduced-volume root.
result = lee_kesler_pure(
    P=P,
    T=T,
    Tc=Tc,
    Pc=Pc,
    acentric_factor=acentric_factor,
    phase="vapor",
)

# SECTION: Result summary
# NOTE: Lee-Kesler is used here for volumetric and departure properties.
print(f"component: {component}")
print(f"Z: {result.compressibility_factor:.8f}")
print(f"molar volume [m3/mol]: {result.molar_volume:.8e}")
print(f"enthalpy departure [J/mol]: {result.enthalpy_departure:.8f}")
print(f"entropy departure [J/(mol.K)]: {result.entropy_departure:.8f}")
print(f"gibbs departure [J/mol]: {result.gibbs_departure:.8f}")
print(f"reduced temperature: {result.reduced_temperature:.8f}")
print(f"reduced pressure: {result.reduced_pressure:.8f}")

# ! Fugacity coefficient support is not claimed by this Lee-Kesler step.
print(
    f"fugacity coefficient supported: {result.metadata['fugacity_coefficient_supported']}")
