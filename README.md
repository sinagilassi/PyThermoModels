# Python Thermodynamic Models

![PyThermoModels](./static/header1.png)

![Downloads](https://img.shields.io/pypi/dm/pyThermoModels) ![PyPI](https://img.shields.io/pypi/v/pyThermoModels) ![Python Version](https://img.shields.io/pypi/pyversions/pyThermoModels.svg) ![License](https://img.shields.io/pypi/l/pyThermoModels) ![Read the Docs](https://img.shields.io/readthedocs/pyThermoModels)


`PyThermoModels` is an open-source Python package designed to facilitate thermodynamic modeling and calculations. This package provides a comprehensive and user-friendly interface to popular thermodynamic models, enabling quick and accurate estimation of key properties.

**Important**: PyThermoModels is now suitable for production and critical applications. However, users are encouraged to thoroughly test the package in their specific use cases.

---

## 📦 Dependency Notice

PyThermoModels is dependent on the `PyThermoDB` Python package. All input data required for calculations should be provided through this package. Ensure that `PyThermoDB` is installed and properly configured before using PyThermoModels.

### 🔗 PyThermoLinkDB

`pyThermoLinkDB` is a Python package that acts as a bridge between `PyThermoDB` and `PyThermoModels`. It is a necessary dependency for enabling seamless integration and data exchange between the two packages. Ensure that `pyThermoLinkDB` is installed and properly configured to use the full functionality of `PyThermoModels`.

---

## ✨ Features

PyThermoModels offers a wide range of features to support thermodynamic modeling and calculations. Below are some of the key features:

### 1️⃣ Equations of State (EOS)

PyThermoModels provides support for popular equations of state (EOS) models, including:

- **Peng-Robinson (PR)**: Suitable for non-polar and mildly polar compounds.
- **Soave-Redlich-Kwong (SRK)**: Commonly used for hydrocarbon systems.
- **Redlich-Kwong (RK)**: A simpler EOS for gases and liquids.
- **van der Waals (vdW)**: A classical EOS for understanding basic thermodynamic behavior.

These EOS models are used to calculate key thermodynamic properties such as:

- Fugacity for pure components and mixtures.
- Phase equilibrium properties.

### 2️⃣ Activity Coefficient Models

For liquid-phase systems, PyThermoModels supports activity coefficient models such as:

- **Non-Random Two-Liquid (NRTL)**: Ideal for highly non-ideal liquid mixtures.
- **UNIQUAC (Universal Quasi-Chemical)**: Suitable for a wide range of liquid mixtures.

These models are essential for:

- Predicting phase equilibria in liquid mixtures.
- Calculating activity coefficients, which are crucial for understanding non-ideal behavior in solutions.

### 🌟 Role of Thermodynamic Properties in Modeling

Thermodynamic properties such as fugacity, activity coefficients, and phase equilibria play a critical role in:

- Designing chemical processes.
- Simulating and optimizing industrial operations.
- Understanding the behavior of chemical systems under various conditions.

By providing robust implementations of these models, PyThermoModels enables researchers and engineers to perform accurate and efficient thermodynamic calculations.

---

## 📥 Installation

Install PyThermoModels with pip:

```python
import pyThermoModels as ptm
# check version
print(ptm.__version__)
```

Import dependent packages:

```python
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
```

---

## 🌡️ Usage Examples

### 1️⃣ Calculating Thermodynamic Properties Using Equations of State (EOS)

PyThermoModels supports the calculation of thermodynamic properties such as fugacity for pure components and mixtures using popular equations of state (EOS) like Peng-Robinson (PR) and Redlich-Kwong (RK).

#### Example: Fugacity Calculation for a Pure Component

```python
import pyThermoModels as ptm
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
```

### 🛠️ Creating a Thermodynamic Database with `PyThermoDB`

Use `PyThermoDB` to create a thermodynamic database for a component:

```python
# propane: thermodb file name
propane_thermodb_file = os.path.join(thermodb_dir, 'propane-1.pkl')
propane_thermodb = ptdb.load_thermodb(propane_thermodb_file)
```

### 🔗 Linking Thermodynamic Data with `PyThermoDBLink`

Use `PyThermoDBLink` to link the thermodynamic data to `PyThermoModels`:

```python
# init thermodb hub
thub1 = ptdblink.init()

# register
thub1.add_thermodb('propane', propane_thermodb)

# * add thermodb rule
thermodb_config_file = os.path.join(
    os.getcwd(), 'test', 'thermodb_config_link.yml')
# all components
thub1.config_thermodb_rule(thermodb_config_file)

# build datasource & equationsource
datasource, equationsource = thub1.build()
```

### 🔑 DataSource and EquationSource in PyThermoModels

`datasource` and `equationsource` are key components in the `PyThermoModels` library. They are used to extract the necessary data and equations required for property calculations.

- **📂 DataSource**: This is responsible for providing the raw data needed for thermodynamic property calculations. It could be a database, a file, or an API that contains experimental or reference data such as temperature, pressure, and material properties.

- **📐 EquationSource**: This provides the mathematical models or equations that describe the relationships between thermodynamic properties. These equations are used to compute properties like enthalpy, entropy, or specific heat based on the data retrieved from the `datasource`.

Together, `datasource` and `equationsource` enable `PyThermoModels` to perform accurate and flexible property calculations by combining empirical data with theoretical models.

#### 🛠️ Example Workflow

1. The `📂 datasource` fetches the required input data (e.g., temperature and pressure values).
2. The `📐 equationsource` supplies the equations or models that define how the properties are calculated.
3. `PyThermoModels` uses both to compute the desired thermodynamic properties.

This modular approach ensures that the library can be easily extended to support new data sources or equations as needed.

### 📊 Calculating Fugacity Using an EOS Model

Choose an equation of state (EOS) to calculate fugacity:

```python
# init
eos = ptm.eos()

# eos model
eos_model = 'PR'
# component
component = "propane"
# phase (optional)
phase = "VAPOR-LIQUID"
# temperature [K]
T = 300.1
# pressure [bar]
P = 9.99

# model input
model_input = {
    # "phase": phase,
    "component": component,
    "pressure": [P, 'bar'],
    "temperature": [T, 'K'],
}

# model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}
```

### 🌡️ Phase Formation Summary

PyThermoModels allows you to summarize the phase formed at a desired temperature and pressure. This feature is particularly useful for understanding phase behavior under specific conditions.

```python
# eos root analysis
res = eos.check_eos_roots_single_component(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source)
print(res)

    # {
    #     'component_name': 'propane',
    #     'pressure': 999000.0,
    #     'pressure_unit': 'Pa',
    #     'temperature': 300.1,
    #     'temperature_unit': 'K',
    #     'root': 3,
    #     'root-no': '1 real root (vapor)',
    #     'phase': 'VAPOR',
    #     'vapor_pressure': 999723.1044,
    #     'vapor_pressure_unit': 'Pa',
    #     'critical_temperature': 369.83,
    #     'critical_temperature_unit': 'K',
    #     'critical_pressure': 4248000.0,
    #     'critical_pressure_unit': 'Pa',
    #     'tolerance': 0.1,
    #     'vapor_pressure_check': 723.1043999999529,
    #     'temperature_equality_value': 69.72999999999996,
    #     'pressure_equality_check': False,
    #     'temperature_equality_check': False
    # }
```

The analysis of `propane` under the given conditions reveals that the system is in the vapor phase at a temperature of 300.1 K and a pressure of 999,000 Pa. The root analysis identifies three roots, with one real root corresponding to the vapor phase. The calculated vapor pressure is 999,723.10 Pa, which is very close to the system pressure, confirming the stability of the vapor phase. These results provide a clear understanding of propane's phase behavior under the specified conditions.

### 🧪 Calculating Fugacity Using EOS Model

To calculate fugacity using an EOS model, follow the steps below:

```python
# Initialize EOS and calculate fugacity
res = eos.cal_fugacity(
    model_name=eos_model,
    model_input=model_input,
    model_source=model_source)
print(res)

    # {
    #     'phase': ['vapor'],
    #     'component': ['propane'],
    #     'vapor': {
    #         'mole_fraction': 1.0,
    #         'temperature': {'value': 300.1, 'unit': 'K', 'symbol': 'T'},
    #         'pressure': {'value': 999000.0, 'unit': 'Pa', 'symbol': 'P'},
    #         'molar_volume': {'value': 0.0020494767707729513, 'unit': 'm3/mol', 'symbol': 'MoVo'},
    #         'compressibility_coefficient': {'value': 0.8205552301471566, 'unit': 'dimensionless', 'symbol': 'Z'},
    #         'fugacity_coefficient': {'value': 0.8493394106824814, 'unit': 'dimensionless', 'symbol': 'phi'},
    #         'fugacity': {'value': 848490.0712717989, 'unit': 'Pa', 'symbol': 'Fug_PURE'},
    #         'mode': 'SINGLE',
    #         'phase': 'VAPOR',
    #         'eos_model': 'PR'
    #     }
    # }
```

This method leverages the selected equation of state (e.g., Peng-Robinson or Redlich-Kwong) to compute the fugacity of a component under specified conditions. The result provides insights into the thermodynamic behavior of the system.

### 2️⃣ Calculating Activity Coefficients Using NRTL and UNIQUAC

PyThermoModels provides support for activity coefficient models such as Non-Random Two-Liquid (NRTL) and UNIQUAC. These models are useful for phase equilibrium calculations in liquid mixtures.

#### Example: Activity Coefficient Calculation Using NRTL

### 🍷 Calculating Activity Coefficients Using NRTL

Initialize the NRTL model as:

```python
# components
components = ['ethanol', 'butyl-methyl-ether']

# feed spec
mole_fraction = {
    'ethanol': 0.4,
    'butyl-methyl-ether': 0.6
}

# model input
activity_model = 'NRTL'

# 📂 model source
model_source = {
    "datasource": datasource,
    "equationsource": equationsource
}

# activity model
activity = ptm.activity(components=components, model_name=activity_model)
print(activity)

# select nrtl
activity_nrtl = activity.nrtl
print(activity_nrtl)
```

### 🧪 Start the Calculation of Activity Coefficients

```python
# 🌡️ operating conditions
# temperature [K]
T = 323.15
# pressure [bar]
P = 30

# 🧮 calculate the interaction parameter matrix (tau_ij)
tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(temperature=T, dg_ij=dg_ij)
print(f"tau_ij: {tau_ij}")
print(f"tau_ij_comp: {tau_ij_comp}")

# 📋 model input
model_input = {
    "mole_fraction": mole_fraction,
    "tau_ij": tau_ij,
    "alpha_ij": alpha_ij
}
```

```python
# 🔍 calculate activity
res_, _ = activity_nrtl.cal(model_input=model_input)
print(res_)

# {'property_name': 'activity coefficients',
# 'components': ['ethanol', 'butyl-methyl-ether'],
# 'mole_fraction': [0.4, 0.6],
# 'value': [1.6102904900236217, 1.187152671740749],
# 'unit': 1,
# 'symbol': 'AcCo_i',
# 'message': 'Calculate activity coefficients for ethanol, butyl-methyl-ether using NRTL model'}

# 📈 excess molar gibbs free energy
gibbs_energy = activity_nrtl.excess_gibbs_free_energy()
print(f"excess gibbs free energy: {gibbs_energy}")

#  {'property_name': 'Excess Molar Gibbs Free Energy',
# 'components': ['ethanol', 'butyl-methyl-ether'],
# 'mole_fraction': [0.4, 0.6],
# 'value': 0.2935004728361861,
# 'unit': 1,
# 'symbol': 'ExMoGiFrEn',
# 'message': 'Excess Gibbs Free Energy for ethanol, butyl-methyl-ether'}
```

---

## 📚 Documentation

For comprehensive documentation, tutorials, and more detailed examples, please visit:

[PyThermoModels Documentation](https://pythermomodels.readthedocs.io/en/latest/)

The documentation includes:

- Detailed API references
- Step-by-step tutorials
- Advanced usage examples
- Integration guides with PyThermoDB and PyThermoLinkDB

We encourage all users to explore the documentation to fully leverage the capabilities of PyThermoModels in their projects.

---

## 📝 License

This project is licensed under the MIT License. You are free to use, modify, and distribute this software in your own applications or projects. However, if you choose to use this app in another app or software, please ensure that my name, Sina Gilassi, remains credited as the original author. This includes retaining any references to the original repository or documentation where applicable. By doing so, you help acknowledge the effort and time invested in creating this project.

---

## ❓ FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/)