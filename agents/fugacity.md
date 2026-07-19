# Fugacity Calculation Guide for Agents

This guide summarizes how to calculate fugacity in `PyThermoModels` using the
build-model-source workflow shown in:

- `examples/eos-models/fugacity-gas using build model source - 2.py`
- `examples/eos-models/fugacity-gas using build model source - 3.py`
- `examples/eos-models/fugacity-liquid using build model source - 2.py`
- `examples/eos-models/fugacity-mixture using build model source.py`

Prefer the newer core functions when writing agent-generated code:

- `calc_gas_fugacity`
- `calc_liquid_fugacity`
- `calc_mixture_fugacity`
- `check_component_eos_roots`
- `check_multi_component_eos_roots`

The older `ptm.eos()` methods are still useful for compatibility examples, but
new examples should use `Component`, `Temperature`, `Pressure`, and
`ModelSource` objects directly.

## Required Imports

```python
import os
from typing import Dict

import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import (
    Component,
    ComponentRule,
    ComponentThermoDBSource,
    Pressure,
    Temperature,
)
from pyThermoModels.core import (
    calc_gas_fugacity,
    calc_liquid_fugacity,
    calc_mixture_fugacity,
    check_component_eos_roots,
    check_multi_component_eos_roots,
)
```

## Shared Setup: Build `ModelSource`

Every fugacity calculation needs a `ModelSource` containing component data and
equations. The examples now show two supported ways to obtain it:

- Build the source inline from existing ThermoDB pickle files.
- Import a reusable source module that builds ThermoDB data from reference
  content and exports `model_source`, `model_source_dict`, and components.

### Required Thermodynamic Inputs

Before preparing `ModelSource`, make sure the source ThermoDB data can provide
the thermodynamic inputs required by the EOS and phase logic:

| Input | Model-source symbol | Required for |
| --- | --- | --- |
| Critical pressure | `Pc` | EOS parameter construction and root analysis |
| Critical temperature | `Tc` | EOS parameter construction and phase/root analysis |
| Acentric factor | `AcFa` | EOS model data completeness and model-source compatibility |
| Vapor-pressure equation | `VaPr` | automatic phase detection, EOS root analysis, and Poynting liquid fugacity |
| Ideal-gas heat-capacity equation | `Cp_IG` | optional for fugacity itself, but commonly included in EOS examples and shared model sources |

The minimum practical rule set for fugacity examples is:

```python
thermodb_rules = {
    "ALL": {
        "DATA": {
            "critical-pressure": "Pc",
            "critical-temperature": "Tc",
            "acentric-factor": "AcFa",
        },
        "EQUATIONS": {
            "CUSTOM-REF-1::vapor-pressure": "VaPr",
            "CUSTOM-REF-1::ideal-gas-heat-capacity": "Cp_IG",
        },
    }
}
```

For pure gas fugacity with a manually supplied vapor phase, the EOS path mainly
needs `Pc` and `Tc`. In normal agent workflows, still include `VaPr` because
automatic phase detection and root checks depend on it. For liquid fugacity with
`liquid_fugacity_mode="Poynting"`, `VaPr` is mandatory.

### Option 1: Build Inline From ThermoDB Files

Build the source from ThermoDB pickle files with
`ptdblink.load_and_build_model_source`.

```python
parent_dir = os.path.dirname(os.path.abspath(__file__))
thermodb_dir = os.path.join(parent_dir, "..", "thermodb")

CO2 = Component(name="carbon dioxide", formula="CO2", state="g")
C3H8 = Component(name="propane", formula="C3H8", state="g")
C4H10 = Component(name="n-butane", formula="C4H10", state="g")

component_thermodb = [
    ComponentThermoDBSource(
        component=CO2,
        source=os.path.join(thermodb_dir, "carbon dioxide-g.pkl"),
    ),
    ComponentThermoDBSource(
        component=C3H8,
        source=os.path.join(thermodb_dir, "propane-g.pkl"),
    ),
    ComponentThermoDBSource(
        component=C4H10,
        source=os.path.join(thermodb_dir, "n-butane-g.pkl"),
    ),
]

thermodb_rules: Dict[str, Dict[str, ComponentRule]] = {
    "ALL": {
        "DATA": {
            "critical-pressure": "Pc",
            "critical-temperature": "Tc",
            "acentric-factor": "AcFa",
        },
        "EQUATIONS": {
            "CUSTOM-REF-1::vapor-pressure": "VaPr",
            "CUSTOM-REF-1::ideal-gas-heat-capacity": "Cp_IG",
        },
    }
}

model_source: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=component_thermodb,
    rules=thermodb_rules,
)
```

### Option 2: Import a Reusable Model Source Module

The `fugacity-gas using build model source - 3.py` example uses
`examples/source/model_source_1.py`. That module builds component ThermoDB data
from `examples.references.reference_2.REFERENCE_CONTENT`, then builds a
`ModelSource` with:

- `build_components_model_source(...)`
- `build_model_source(...)`

Use this pattern when an example or agent workflow should centralize component
source construction in one file.

```python
from examples.source.model_source_1 import (
    C3H8,
    model_source,
    model_source_dict,
)
```

Then use:

- `model_source` with the newer core helpers such as `calc_gas_fugacity`.
- `model_source_dict` with older `ptm.eos()` methods that expect a dictionary
  containing `datasource` and `equationsource`.

```python
temperature = Temperature(value=300.1, unit="K")
pressure = Pressure(value=9.99, unit="bar")

root_result = check_component_eos_roots(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    component_key="Name-State",
)

fugacity_result = calc_gas_fugacity(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    component_key="Name-State",
)
```

Important: importing `examples.source.model_source_1` executes its build steps.
It can create or update ThermoDB files in `examples/thermodb` because it calls
`build_component_thermodb_from_reference(..., thermodb_save=True, ...)`.

### Component Keys

`component_key` controls how the component is matched against `model_source`.

- Use `"Name-State"` for keys like `propane-g`, `carbon dioxide-g`,
  `n-butane-g`.
- Use `"Formula-State"` for keys like `C3H8-g`, `CO2-g`, `C4H10-g`.

The `Component` object must contain the matching `name`, `formula`, and `state`.
For mixtures, each `Component` must also include `mole_fraction`.

## Gas Fugacity: Pure Component

Use `calc_gas_fugacity` for a pure-component vapor/gas calculation. Optionally
run `check_component_eos_roots` first to confirm the phase selected by the EOS
root analysis.

```python
temperature = Temperature(value=300.1, unit="K")
pressure = Pressure(value=9.99, unit="bar")

root_result = check_component_eos_roots(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name="PR",
    component_key="Name-State",
)

fugacity_result = calc_gas_fugacity(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name="PR",
    component_key="Name-State",
)
```

Expected output is a `ComponentGasFugacityResult`-style object containing the
phase, component, fugacity coefficient, fugacity, compressibility coefficient,
molar volume, temperature, pressure, and EOS model.

## Liquid Fugacity: Pure Component

Use `calc_liquid_fugacity` for pure-liquid fugacity. The liquid helper defaults
the phase to `"LIQUID"` internally, but passing `phase="LIQUID"` is acceptable
when clarity is useful.

```python
temperature = Temperature(value=340, unit="K")
pressure = Pressure(value=30, unit="bar")

root_result = check_component_eos_roots(
    component=CO2,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name="PR",
    component_key="Name-State",
)

fugacity_result = calc_liquid_fugacity(
    component=CO2,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name="PR",
    liquid_fugacity_mode="Poynting",
    component_key="Name-State",
)
```

`liquid_fugacity_mode` can be:

- `"EOS"`: calculate liquid fugacity from the liquid EOS root.
- `"Poynting"`: use the Poynting correction path shown in the liquid example.

## Mixture Fugacity

Use `calc_mixture_fugacity` for multi-component fugacity. Each component must
carry its mole fraction. The mole fractions should sum to 1.0.

```python
CO2_mix = Component(
    name="carbon dioxide",
    formula="CO2",
    state="g",
    mole_fraction=0.15,
)
C4H10_mix = Component(
    name="n-butane",
    formula="C4H10",
    state="g",
    mole_fraction=0.85,
)
components = [CO2_mix, C4H10_mix]

temperature = Temperature(value=444, unit="K")
pressure = Pressure(value=10, unit="bar")

# Optional binary interaction parameters. Omit this block when unavailable.
k_ij = [
    [0, 0.18],
    [0.18, 0],
]

root_result = check_multi_component_eos_roots(
    components=components,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name="RK",
    component_key="Name-State",
)

fugacity_result = calc_mixture_fugacity(
    components=components,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    model_name="RK",
    component_key="Name-State",
    k_ij=k_ij,
)
```

Expected output is a `MixtureFugacityResult`-style object with per-component
fugacity coefficients and fugacities for the detected phase or phases.

`k_ij` is optional for mixture fugacity. Pass it only when binary interaction
parameters are known and should be included in the EOS mixing rule. If no
`k_ij` is supplied, the calculation uses the model's default interaction
parameter handling.

## Old API to New API Mapping

The examples include both old and new call styles. When preparing new code,
prefer the new API.

| Task | Old API | New API |
| --- | --- | --- |
| Pure EOS root analysis | `eos.check_eos_roots_single_component(...)` | `check_component_eos_roots(...)` |
| Pure gas fugacity | `eos.cal_fugacity(...)` | `calc_gas_fugacity(...)` |
| Pure liquid fugacity | `eos.cal_fugacity(..., liquid_fugacity_mode=...)` | `calc_liquid_fugacity(...)` |
| Mixture EOS root analysis | `eos.check_eos_roots_multi_component(...)` | `check_multi_component_eos_roots(...)` |
| Mixture fugacity | `eos.cal_fugacity_mixture(...)` | `calc_mixture_fugacity(...)` |

When using the reusable source module:

- Old API calls should receive `model_source=model_source_dict`.
- New API calls should receive `model_source=model_source`.

## Common Agent Checklist

1. Define `Component` objects with correct `name`, `formula`, and `state`.
2. Wrap ThermoDB pickle paths in `ComponentThermoDBSource`.
3. Build `ModelSource` inline with rules mapping `Pc`, `Tc`, `AcFa`, `VaPr`,
   and optionally `Cp_IG`, or import a prepared `model_source` from
   `examples.source.model_source_1`.
4. Use `Temperature(value=..., unit="K")` and `Pressure(value=..., unit="bar")`
   or another supported unit.
5. Pick `model_name`, commonly `"PR"` for Peng-Robinson or `"RK"` for
   Redlich-Kwong examples.
6. Run the relevant root-analysis helper before fugacity when phase behavior is
   unclear.
7. Use the same `component_key` in root analysis and fugacity calculation.
8. For mixtures, set `mole_fraction` on each `Component`; optionally pass
   `k_ij` through `**kwargs` when binary interaction parameters are available.

## Common Failure Points

- Component key mismatch: `component_key="Name-State"` requires model-source
  keys based on names, while `"Formula-State"` requires formula-based keys.
- Missing vapor-pressure equation: liquid and phase-root checks depend on
  `VaPr`; include the vapor-pressure rule in `thermodb_rules`.
- Missing critical properties: EOS fugacity requires `Pc`, `Tc`, and `AcFa`.
- Incorrect state suffix: examples use gas-state ThermoDB files like
  `propane-g.pkl` and component states like `state="g"`.
- Incomplete mixture data: mixture calculations need `mole_fraction` values on
  the `Component` objects, not just a separate feed dictionary.
- Unit ambiguity: always construct `Temperature` and `Pressure` with explicit
  units.
- Import side effects: reusable source modules may build and save ThermoDB files
  at import time; account for that when using them in tests or scripts.

## Minimal Decision Guide

- Pure gas or vapor phase: `calc_gas_fugacity`.
- Pure liquid phase: `calc_liquid_fugacity`.
- Mixture with mole fractions: `calc_mixture_fugacity`.
- Unsure phase/root behavior: call the matching `check_*_eos_roots` helper
  first.
