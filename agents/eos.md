# Fugacity Calculation Guide for Agents

This guide summarizes the current fugacity workflow in `PyThermoModels` using
the build-model-source API shown in:

- `examples/eos-models/fugacity-gas using build model source - 4.py`
- `examples/eos-models/fugacity-liquid using build model source - 3.py`
- `examples/eos-models/fugacity-mixture using build model source - 1.py`
- `examples/source/model_source_1.py`
- `examples/source/model_source_2.py`

Use the `pyThermoModels.core` helpers with a `ModelSource` object. Do not use
the older `ptm.eos()` object API for new agent-generated fugacity examples.

## Current API

Import the core helpers directly:

```python
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoModels.core import (
    calc_gas_fugacity,
    calc_liquid_fugacity,
    calc_mixture_fugacity,
    check_component_eos_roots,
    check_multi_component_eos_roots,
)
```

The current helper signatures use structured objects:

- `component` or `components`: `Component` object or list of `Component`
  objects.
- `temperature`: `Temperature(value=..., unit=...)`.
- `pressure`: `Pressure(value=..., unit=...)`.
- `model_source`: a `pyThermoLinkDB.models.ModelSource` instance.
- `model_name`: optional EOS name, default is `"SRK"`.
- `component_key`: optional component lookup mode, default is `"Name-State"`.

## Function Argument Reference

When writing agent guidance or generated examples, describe every argument used
by the selected function, including optional arguments that rely on defaults.
Defaults are part of the API contract and should be chosen intentionally.

### `check_component_eos_roots(...)`

Use for pure-component EOS root and phase analysis before pure gas or liquid
fugacity calculations.

- `component`: required `Component`. Must include `name`, `formula`, and
  `state` matching the selected `component_key`.
- `pressure`: required `Pressure`. Use explicit units, for example
  `Pressure(value=9.99, unit="bar")`.
- `temperature`: required `Temperature`. Use explicit units, for example
  `Temperature(value=300.1, unit="K")`.
- `model_source`: required `ModelSource`. Pass the object exported by the
  source module, not `model_source_dict`.
- `model_name`: optional EOS model. Default is `"SRK"`; accepted values are
  `"SRK"`, `"PR"`, `"vdW"`, and `"RK"`.
- `component_key`: optional lookup mode. Default is `"Name-State"`; use
  `"Formula-State"` only when the source keys are formula-state keys.
- `**kwargs`: optional calculation controls. Common keys include `phase` for an
  explicit phase hint and `tolerance` for numerical/root comparison tolerance.

### `check_multi_component_eos_roots(...)`

Use for mixture EOS root and phase analysis before mixture fugacity.

- `components`: required list of `Component` objects. Each component must carry
  `mole_fraction` for mixture calculations.
- `pressure`: required `Pressure` with explicit units.
- `temperature`: required `Temperature` with explicit units.
- `model_source`: required `ModelSource` containing all mixture components.
- `model_name`: optional EOS model. Default is `"SRK"`; accepted values are
  `"SRK"`, `"PR"`, `"vdW"`, and `"RK"`.
- `bubble_point_pressure_mode`: optional bubble-point method. Default is
  `"Raoult"`.
- `dew_point_pressure_mode`: optional dew-point method. Default is `"Raoult"`.
- `component_key`: optional lookup mode. Default is `"Name-State"`.
- `**kwargs`: optional calculation controls such as `tolerance`.

### `calc_gas_fugacity(...)`

Use for pure-component vapor/gas fugacity.

- `component`: required `Component` for the pure component.
- `pressure`: required `Pressure` with explicit units.
- `temperature`: required `Temperature` with explicit units.
- `model_source`: required `ModelSource`; do not pass the old dictionary form.
- `model_name`: optional EOS model. Default is `"SRK"`; accepted values are
  `"SRK"`, `"PR"`, `"RK"`, and `"vdW"`.
- `solver_method`: optional numerical solver. Default is `"ls"`; accepted
  values are `"ls"`, `"newton"`, `"fsolve"`, and `"root"`.
- `component_key`: optional lookup mode. Default is `"Name-State"`.
- `phase_names`: optional list of phases considered by the calculation. Default
  is `["VAPOR", "LIQUID", "SUPERCRITICAL", "VAPOR-LIQUID"]`.
- `**kwargs`: optional calculation controls. Common keys include `phase` and
  `tolerance`.

### `calc_liquid_fugacity(...)`

Use for pure-component liquid fugacity.

- `component`: required `Component` for the pure component.
- `pressure`: required `Pressure` with explicit units.
- `temperature`: required `Temperature` with explicit units.
- `model_source`: required `ModelSource`; do not pass the old dictionary form.
- `model_name`: optional EOS model. Default is `"SRK"`; accepted values are
  `"SRK"`, `"PR"`, `"RK"`, and `"vdW"`.
- `solver_method`: optional numerical solver. Default is `"ls"`; accepted
  values are `"ls"`, `"newton"`, `"fsolve"`, and `"root"`.
- `liquid_fugacity_mode`: optional liquid method. Default is `"EOS"`; use
  `"Poynting"` for the Poynting correction path shown in the current liquid
  example.
- `component_key`: optional lookup mode. Default is `"Name-State"`.
- `phase_names`: optional list of phases considered by the calculation. Default
  is `["VAPOR", "LIQUID", "SUPERCRITICAL", "VAPOR-LIQUID"]`.
- `**kwargs`: optional calculation controls. Common keys include `phase` and
  `tolerance`.

### `calc_mixture_fugacity(...)`

Use for multi-component fugacity.

- `components`: required list of `Component` objects. Each component must carry
  `mole_fraction`; the fractions should sum to 1.0.
- `pressure`: required `Pressure` with explicit units.
- `temperature`: required `Temperature` with explicit units.
- `model_source`: required `ModelSource` containing every component in the
  mixture.
- `model_name`: optional EOS model. Default is `"SRK"`; accepted values are
  `"SRK"`, `"PR"`, `"RK"`, and `"vdW"`.
- `solver_method`: optional numerical solver. Default is `"ls"`; accepted
  values are `"ls"`, `"newton"`, `"fsolve"`, and `"root"`.
- `liquid_fugacity_mode`: optional liquid treatment for liquid-phase mixture
  paths. Default is `"EOS"`; `"Poynting"` is available where the target path
  supports it.
- `component_key`: optional lookup mode. Default is `"Name-State"`.
- `phase_names`: optional list of phases considered by the calculation. Default
  is `["VAPOR", "LIQUID", "SUPERCRITICAL", "VAPOR-LIQUID"]`.
- `**kwargs`: optional calculation controls such as `tolerance`; pass binary
  interaction data only when the called path supports and needs it.

## Building Model Sources

Reusable model-source modules should build ThermoDB components from reference
content, then convert those ThermoDB objects into a `ModelSource`.

Use this API:

```python
from pyThermoDB import build_component_thermodb_from_reference
from pyThermoLinkDB import build_components_model_source, build_model_source
from pyThermoLinkDB.models import ComponentModelSource, ModelSource
from pythermodb_settings.models import Component
```

Pattern:

```python
thermodb_components = []

for comp in components:
    thermodb_component = build_component_thermodb_from_reference(
        component_name=comp.name,
        component_formula=comp.formula,
        component_state=comp.state,
        reference_content=REFERENCE_CONTENT,
        ignore_state_props=["MW", "VaPr"],
        thermodb_save=True,
        thermodb_save_path=thermodb_dir,
    )
    if thermodb_component is None:
        raise ValueError(f"thermodb_component for {comp.name} is None")
    thermodb_components.append(thermodb_component)

component_model_source: list[ComponentModelSource] = build_components_model_source(
    components_thermodb=thermodb_components,
    rules=None,
)

model_source: ModelSource = build_model_source(source=component_model_source)
```

### Source-Builder Arguments

The reusable source modules use these builder functions before any EOS helper is
called.

`build_component_thermodb_from_reference(...)` builds one `ComponentThermoDB`
object from reference content.

- `component_name`: required component name, normally `comp.name` from the
  `Component` object.
- `component_formula`: required chemical formula, normally `comp.formula`.
- `component_state`: required state suffix, normally `comp.state` such as `"g"`
  or `"l"`.
- `reference_content`: required reference-data object, such as
  `REFERENCE_CONTENT` imported from an `examples.references` module.
- `ignore_state_props`: optional list of state properties to ignore while
  building the component ThermoDB. Current source modules pass `["MW", "VaPr"]`.
- `thermodb_save`: optional save flag. Current source modules pass `True`, which
  means imports may write ThermoDB files.
- `thermodb_save_path`: required when `thermodb_save=True`; directory where the
  generated ThermoDB files are written.

`build_components_model_source(...)` converts one or more component ThermoDB
objects into component model-source objects.

- `components_thermodb`: required list of `ComponentThermoDB` objects produced
  by `build_component_thermodb_from_reference(...)`.
- `rules`: optional mapping that controls property/equation selection. Current
  source modules pass `None` so the builder uses its default matching behavior.

`build_model_source(...)` combines component model-source objects into the
`ModelSource` required by the new EOS helpers.

- `source`: required list of `ComponentModelSource` objects returned by
  `build_components_model_source(...)`.

Important points:

- Prefer `rules=None` with `build_components_model_source(...)` when following
  `model_source_1.py` and `model_source_2.py`.
- Export `model_source` from reusable source modules and pass it directly to
  the new core helpers.
- `model_source_dict = {"datasource": ..., "equationsource": ...}` is only for
  older compatibility code. New fugacity examples should not need it.
- Importing these source modules can create or update files under the selected
  `thermodb_save_path` because `thermodb_save=True`.

## Reusable Source Modules

Use `examples.source.model_source_1` for the pure propane gas example:

```python
from examples.source.model_source_1 import C3H8, model_source
```

`model_source_1.py` currently builds a source for:

- `C3H8` as propane gas.

Use `examples.source.model_source_2` for liquid and mixture examples:

```python
from examples.source.model_source_2 import C3H8, model_source
from examples.source.model_source_2 import CO2, C4H10, model_source
```

`model_source_2.py` currently builds a source for:

- `CO2` as carbon dioxide gas.
- `C3H8` as propane gas.
- `CH3OH` as methanol gas.
- `C4H10` as n-butane gas.

When a reusable module is imported, use its exported `Component` objects when
possible. For mixtures, recreate the components with `mole_fraction` values,
as shown below.

## Component Keys

`component_key` controls how components are looked up in `model_source`.

- `"Name-State"` uses keys such as `propane-g`, `carbon dioxide-g`,
  `n-butane-g`.
- `"Formula-State"` uses keys such as `C3H8-g`, `CO2-g`, `C4H10-g`.

The current fugacity examples use `"Name-State"` explicitly for gas and mixture
and rely on the same default for liquid.

## Gas Fugacity

Use `check_component_eos_roots` first when root or phase behavior matters, then
use `calc_gas_fugacity`.

```python
from pythermodb_settings.models import Pressure, Temperature
from pyThermoModels.core import calc_gas_fugacity, check_component_eos_roots
from examples.source.model_source_1 import C3H8, model_source

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

`model_name` is optional. If omitted, the helper uses `"SRK"`.

## Liquid Fugacity

Use `calc_liquid_fugacity` for pure-component liquid fugacity. The current
liquid example uses propane from `model_source_2.py` and the Poynting mode.

```python
from pythermodb_settings.models import Pressure, Temperature
from pyThermoModels.core import calc_liquid_fugacity, check_component_eos_roots
from examples.source.model_source_2 import C3H8, model_source

temperature = Temperature(value=340, unit="K")
pressure = Pressure(value=30, unit="bar")

root_result = check_component_eos_roots(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
)

fugacity_result = calc_liquid_fugacity(
    component=C3H8,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    liquid_fugacity_mode="Poynting",
)
```

`liquid_fugacity_mode` can be:

- `"EOS"`: calculate liquid fugacity from the liquid EOS root.
- `"Poynting"`: use the Poynting correction path.

## Mixture Fugacity

Use `calc_mixture_fugacity` for multi-component fugacity. Each component in the
calculation must include `mole_fraction`.

```python
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoModels.core import (
    calc_mixture_fugacity,
    check_multi_component_eos_roots,
)
from examples.source.model_source_2 import model_source

temperature = Temperature(value=444, unit="K")
pressure = Pressure(value=10, unit="bar")

N0s = {
    "CO2-g": 0.15,
    "n-butane-g": 0.85,
}

CO2 = Component(
    name="carbon dioxide",
    formula="CO2",
    state="g",
    mole_fraction=N0s.get("CO2-g", 0),
)
C4H10 = Component(
    name="n-butane",
    formula="C4H10",
    state="g",
    mole_fraction=N0s.get("n-butane-g", 0),
)
components = [CO2, C4H10]

root_result = check_multi_component_eos_roots(
    components=components,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    component_key="Name-State",
)

fugacity_result = calc_mixture_fugacity(
    components=components,
    temperature=temperature,
    pressure=pressure,
    model_source=model_source,
    component_key="Name-State",
)
```

The mixture helper also accepts optional settings such as `model_name`,
`solver_method`, `liquid_fugacity_mode`, `phase_names`, and keyword arguments.
Only pass optional binary interaction data when the target helper path supports
and needs it.

## Old API to New API Mapping

| Task | Old API | New API |
| --- | --- | --- |
| Pure EOS root analysis | `eos.check_eos_roots_single_component(...)` | `check_component_eos_roots(...)` |
| Pure gas fugacity | `eos.cal_fugacity(...)` | `calc_gas_fugacity(...)` |
| Pure liquid fugacity | `eos.cal_fugacity(..., liquid_fugacity_mode=...)` | `calc_liquid_fugacity(...)` |
| Mixture EOS root analysis | `eos.check_eos_roots_multi_component(...)` | `check_multi_component_eos_roots(...)` |
| Mixture fugacity | `eos.cal_fugacity_mixture(...)` | `calc_mixture_fugacity(...)` |

## Agent Checklist

1. Import a reusable `model_source` from `examples.source.model_source_1` or
   `examples.source.model_source_2`, or build one with
   `build_components_model_source(..., rules=None)` and `build_model_source`.
2. Use `Temperature` and `Pressure` objects with explicit units.
3. Use `Component` objects with matching `name`, `formula`, and `state`.
4. For mixtures, set `mole_fraction` on each `Component`; do not rely only on a
   separate feed dictionary.
5. Pass the `ModelSource` object directly to the new core helper.
6. Keep `component_key` consistent between root checks and fugacity calls.
7. Add `model_name` only when the example needs a specific EOS; otherwise use
   the helper default.

## Common Failure Points

- Passing `model_source_dict` to a new core helper. New helpers require a
  `ModelSource` object.
- Component-key mismatch between the `Component` object and `model_source`.
- Missing `mole_fraction` values in mixture components.
- Ambiguous units. Always use `Temperature(value=..., unit=...)` and
  `Pressure(value=..., unit=...)`.
- Import side effects from source modules that save ThermoDB files.
- Copying old examples that instantiate `ptm.eos()` instead of importing from
  `pyThermoModels.core`.
