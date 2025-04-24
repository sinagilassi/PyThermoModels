# 🧪 Examples

This documentation provides examples of how to use the PyThermoModels library for various thermodynamic calculations including fugacity coefficients and activity coefficients.

## 📚 Table of Contents

1. [Fugacity Calculations](#fugacity-calculations)
   - [Gas Phase Fugacity](#gas-phase-fugacity)
   - [Liquid Phase Fugacity](#liquid-phase-fugacity)
   - [Mixture Fugacity](#mixture-fugacity)
2. [Activity Coefficient Calculations](#activity-coefficient-calculations)
   - [NRTL Model - Binary Mixture (Ethanol/Butyl-Methyl-Ether)](#nrtl-model---binary-mixture-ethanolbutyl-methyl-ether)
   - [NRTL Model - Binary Mixture (Methanol/Ethanol)](#nrtl-model---binary-mixture-methanolethanol)
   - [UNIQUAC Model](#uniquac-model)

## 💨 Fugacity Calculations

### 🌬️ Gas Phase Fugacity

This example demonstrates how to calculate the fugacity coefficient for a pure component in the gas phase using various equation of state (EOS) models.

```python
# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# Load thermodynamic database files
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# Load component data
CO2_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'carbon dioxide-1.pkl'))
acetylene_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'acetylene-1.pkl'))
n_butane_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'n-butane-1.pkl'))
ethanol_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'ethanol-1.pkl'))
methanol_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'methanol-1.pkl'))
butene_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, '1-butene-1.pkl'))
propane_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'propane-1.pkl'))
methane_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'methane-1.pkl'))

# Initialize thermodb hub
thub1 = ptdblink.init()

# Add component thermodb
thub1.add_thermodb('CO2', CO2_thermodb)
thub1.add_thermodb('CH4', methane_thermodb)
thub1.add_thermodb('EtOH', ethanol_thermodb)
thub1.add_thermodb('MeOH', methanol_thermodb)
thub1.add_thermodb('acetylene', acetylene_thermodb)
thub1.add_thermodb('1-butene', butene_thermodb)
thub1.add_thermodb('n-butane', n_butane_thermodb)
thub1.add_thermodb('propane', propane_thermodb)

# Add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build datasource & equationsource
datasource, equationsource = thub1.build()

# Initialize pyThermoModels
tm = ptm.init()
eos = ptm.eos()

# Example: Propane at 300.1 K and 9.99 bar (based on literature example)
component = "propane"
T = 300.1  # Temperature [K]
P = 9.99   # Pressure [bar]

# Model input
model_input = {
    "component": component,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# Model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# Check reference
res_ = tm.check_fugacity_reference('PR')  # Check Peng-Robinson EOS references

# Check EOS roots
res = eos.check_eos_roots_single_component(
    model_name='PR',
    model_input=model_input,
    model_source=model_source)

# Calculate fugacity
result = eos.cal_fugacity(
    model_name='PR',
    model_input=model_input,
    model_source=model_source)
```

The example above calculates the fugacity coefficient for propane at 300.1 K and 9.99 bar using the Peng-Robinson equation of state. The code demonstrates:

1. Loading thermodynamic data from pickle files
2. Setting up the thermodynamic database configuration
3. Setting operating conditions (temperature and pressure)
4. Checking the EOS roots to determine the phase state
5. Calculating the fugacity coefficient

### 💧 Liquid Phase Fugacity

This example shows how to calculate the fugacity coefficient for a pure component in the liquid phase.

```python
# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# Load thermodynamic database files
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# Load component data
CO2_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'carbon dioxide-1.pkl'))
acetylene_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'acetylene-1.pkl'))
n_butane_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'n-butane-1.pkl'))
ethanol_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'ethanol-1.pkl'))
methanol_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'methanol-1.pkl'))
propane_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'propane-1.pkl'))

# Initialize thermodb hub
thub1 = ptdblink.init()

# Add component thermodb (example with two components)
thub1.add_thermodb('acetylene', acetylene_thermodb)
thub1.add_thermodb('propane', propane_thermodb)

# Add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build datasource & equationsource
datasource, equationsource = thub1.build()

# Initialize pyThermoModels
tm = ptm.init()
eos = ptm.eos()

# Example: Propane at 300.1 K and 10 bar
component = 'propane'
T = 300.1  # Temperature [K]
P = 10     # Pressure [bar]

# Model input
model_input = {
    "component": component,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# Model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# Check reference
res_ = tm.check_fugacity_reference('PR')

# Check EOS roots
res = eos.check_eos_roots_single_component(
    model_name='PR',
    model_input=model_input,
    model_source=model_source)

# Calculate liquid phase fugacity using Poynting correction
result = eos.cal_fugacity(
    model_name='PR',
    model_input=model_input,
    model_source=model_source,
    liquid_fugacity_mode='Poynting')
```

The liquid phase fugacity calculation uses the Poynting correction to account for the pressure effect on liquid fugacity. This example uses propane at 300.1 K and 10 bar.

### 🔄 Mixture Fugacity

This example demonstrates how to calculate fugacity coefficients for multi-component mixtures.

```python
# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# Load thermodynamic database files
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# Load component data
N2_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'nitrogen-1.pkl'))
CH4_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'methane-1.pkl'))
C2H6_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'ethane-1.pkl'))
n_butane_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'n-butane-1.pkl'))
CO2_thermodb = ptdb.load_thermodb(os.path.join(thermodb_dir, 'carbon dioxide-1.pkl'))

# Initialize pyThermoModels
tm = ptm.init()
eos = ptm.eos()

# Initialize thermodb hub
thub1 = ptdblink.init()

# Add component thermodb
thub1.add_thermodb('N2', N2_thermodb)
thub1.add_thermodb('CH4', CH4_thermodb)
thub1.add_thermodb('C2H6', C2H6_thermodb)
thub1.add_thermodb('n-butane', n_butane_thermodb)
thub1.add_thermodb('CO2', CO2_thermodb)

# Add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build datasource & equationsource
datasource, equationsource = thub1.build()

# Example: CO2/n-butane mixture (Example 9.5, The Thermodynamics of Phase and Reaction Equilibria)
eos_model = 'RK'  # Redlich-Kwong equation of state

# Mixture composition
N0s = {
    'CO2': 0.15,
    'n-butane': 0.85,
}

# Operating conditions
T = 444    # Temperature [K]
P = 10     # Pressure [bar]

# Binary interaction parameters
k_ij = [[0, 0.18],
        [0.18, 0]]

# Model input
model_input = {
    "feed-specification": N0s,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# Model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# Check reference
res_ = tm.check_fugacity_reference(eos_model)

# Check EOS roots for mixture
res_ = eos.check_eos_roots_multi_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source)

# Calculate mixture fugacity
result = eos.cal_fugacity_mixture(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source,
    k_ij=k_ij)
```

This example calculates fugacity coefficients for a CO2/n-butane mixture at 444 K and 10 bar using the Redlich-Kwong equation of state. The calculation includes binary interaction parameters (k_ij) to account for the non-ideal behavior of the mixture.

## 🧪 Activity Coefficient Calculations

### 🍸 NRTL Model - Binary Mixture (Ethanol/Butyl-Methyl-Ether)

This example demonstrates how to use the NRTL (Non-Random Two-Liquid) model to calculate activity coefficients for a binary mixture of ethanol and butyl-methyl-ether.

```python
# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL, UNIQUAC
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# Load thermodynamic database files
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# Load NRTL parameters for ethanol/butyl-methyl-ether
nrtl_path = os.path.join(thermodb_dir, 'thermodb_nrtl_ethanol_butyl-methyl-ether_1.pkl')
thermodb_nrtl_1 = ptdb.load_thermodb(nrtl_path)

# Initialize thermodb hub
thub1 = ptdblink.init()

# Add component thermodb
thub1.add_thermodb('nrtl', thermodb_nrtl_1)

# Add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build datasource & equationsource
datasource, equationsource = thub1.build()

# Configure activity model
components = ['ethanol', 'butyl-methyl-ether']
activity_model = 'NRTL'

# Initialize activity model
activity = ptm.activity(components=components, model_name=activity_model)
activity_nrtl = activity.nrtl

# Mixture composition
mole_fraction = {
    'ethanol': 0.4,
    'butyl-methyl-ether': 0.6
}

# Get NRTL parameters from the database
non_randomness_parameters = thermodb_nrtl_1.select('non-randomness-parameters')

# Extract required parameters
dg_ij = non_randomness_parameters.ijs(f"dg | {components[0]} | {components[1]}")
alpha_ij = non_randomness_parameters.ijs(f"alpha | {components[0]} | {components[1]}")

# Operating conditions
T = 323.15  # Temperature [K]
P = 30      # Pressure [bar]

# Calculate interaction parameters at specified temperature
tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(temperature=T, dg_ij=dg_ij)

# Model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "alpha_ij": alpha_ij
}

# Calculate activity coefficients
res_, others_ = activity_nrtl.cal(model_input=model_input)

# Calculate excess Gibbs free energy
G_ij = others_['G_ij']
gibbs_energy = activity_nrtl.excess_gibbs_free_energy(
    mole_fraction=mole_fraction, G_ij=G_ij, tau_ij=tau_ij)
```

This example calculates activity coefficients and excess Gibbs free energy for an ethanol/butyl-methyl-ether mixture using the NRTL model.

### 🧪 NRTL Model - Binary Mixture (Methanol/Ethanol)

This example demonstrates the NRTL model for a methanol/ethanol mixture.

```python
# import packages/modules
import os
from rich import print
import pyThermoModels as ptm
from pyThermoModels import NRTL, UNIQUAC
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# Load thermodynamic database files
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# Load NRTL parameters for methanol/ethanol
nrtl_path = os.path.join(thermodb_dir, 'thermodb_nrtl_methanol_ethanol_1.pkl')
thermodb_nrtl_1 = ptdb.load_thermodb(nrtl_path)

# Initialize thermodb hub
thub1 = ptdblink.init()

# Add component thermodb
thub1.add_thermodb('nrtl', thermodb_nrtl_1)

# Add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build datasource & equationsource
datasource, equationsource = thub1.build()

# Configure activity model
components = ['methanol', 'ethanol']
activity_model = 'NRTL'

# Initialize activity model
activity = ptm.activity(components=components, model_name=activity_model)
activity_nrtl = activity.nrtl

# Mixture composition
mole_fraction = {
    'methanol': 0.0316,
    'ethanol': 0.9684
}

# Get NRTL parameters from the database
non_randomness_parameters = thermodb_nrtl_1.select('non-randomness-parameters')

# Extract required parameters
a_ij = non_randomness_parameters.ijs(f"a | {components[0]} | {components[1]}")
b_ij = non_randomness_parameters.ijs(f"b | {components[0]} | {components[1]}")
c_ij = non_randomness_parameters.ijs(f"c | {components[0]} | {components[1]}")
alpha_ij = non_randomness_parameters.ijs(f"alpha | {components[0]} | {components[1]}")

# Operating conditions
T = 350.1546257  # Temperature [K]
P = 30           # Pressure [bar]

# Calculate interaction parameters at specified temperature
dg_ij, dg_ij_comp = activity_nrtl.cal_dg_ij_M1(
    temperature=T, a_ij=a_ij, b_ij=b_ij, c_ij=c_ij)
tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(temperature=T, dg_ij=dg_ij)

# Model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "alpha_ij": alpha_ij
}

# Calculate activity coefficients
res_, others_ = activity_nrtl.cal(model_input=model_input)
G_ij = others_['G_ij']

# Calculate excess Gibbs free energy
gibbs_energy = activity_nrtl.excess_gibbs_free_energy(
    mole_fraction=mole_fraction, G_ij=G_ij, tau_ij=tau_ij)
```

This example demonstrates the more general approach to NRTL calculations where the interaction parameters are calculated from temperature-dependent expressions.

### 🧮 UNIQUAC Model

This example demonstrates how to use the UNIQUAC (Universal Quasi-Chemical) model to calculate activity coefficients.

```python
# import packages/modules
import os
import numpy as np
from rich import print
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink

# Load thermodynamic database files
thermodb_dir = os.path.join(os.getcwd(), 'test', 'thermodb')

# Load NRTL parameters for methanol/ethanol
nrtl_path = os.path.join(thermodb_dir, 'thermodb_nrtl_methanol_ethanol_1.pkl')
thermodb_nrtl_1 = ptdb.load_thermodb(nrtl_path)

# Initialize thermodb hub
thub1 = ptdblink.init()

# Add component thermodb
thub1.add_thermodb('nrtl', thermodb_nrtl_1)

# Add thermodb rule
thermodb_config_file = os.path.join(os.getcwd(), 'test', 'thermodb_config_link.yml')
thub1.config_thermodb_rule(thermodb_config_file)

# Build datasource & equationsource
datasource, equationsource = thub1.build()

# Configure activity model
components = ['methanol', 'ethanol']
activity_model = 'UNIQUAC'

# Model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# Initialize activity model
activity = ptm.activity(components=components, model_name=activity_model)
activity_uniquac = activity.uniquac

# Mixture composition
mole_fraction = {
    'methanol': 0.4,
    'ethanol': 0.6
}

# UNIQUAC parameters
# Volume and surface area parameters
r_i = [1.4311, 2.1055]  # MeOH, EtOH
q_i = [1.4320, 1.8920]  # MeOH, EtOH

# Binary interaction parameters
tau_ij = [
    [1, 1.031995],
    [1.309036, 1]
]

# Operating conditions
T = 350.1546257  # Temperature [K]
P = 30           # Pressure [bar]

# Model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "r_i": r_i,
    "q_i": q_i
}

# Calculate activity coefficients
res_, others_ = activity_uniquac.cal(model_input=model_input)

# Calculate excess Gibbs free energy
gibbs_energy = activity_uniquac.excess_gibbs_free_energy(
    mole_fraction=mole_fraction, tau_ij=tau_ij, r_i=r_i, q_i=q_i)
```

This example demonstrates activity coefficient calculations using the UNIQUAC model for a methanol/ethanol mixture. The UNIQUAC model requires volume and surface area parameters (r_i and q_i) for each component in addition to binary interaction parameters.