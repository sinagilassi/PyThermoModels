# Activity Coefficient Calculation Guide for Agents

This guide summarizes how to calculate liquid-mixture activity coefficients in
`PyThermoModels` using the inline-reference and build-model-source workflow
shown in:

- `examples/activity-models/calc activity using nrtl inline source and build model source - 1.py`
- `examples/activity-models/calc activity using uniquac inline source and build model source - 1.py`

Prefer the newer core helper when writing agent-generated code:

- `calc_activity_coefficient`

The examples use `PyThermoDB` to build ThermoDB objects directly from inline
reference content, `PyThermoLinkDB` to build a `ModelSource`, and
`PyThermoModels` to calculate activity coefficients and excess Gibbs energy.

## Required Imports

```python
import os

from pyThermoDB import (
    build_component_thermodb_from_reference,
    build_mixture_thermodb_from_reference,
    ComponentThermoDB,
    MixtureThermoDB,
)
from pyThermoLinkDB import (
    build_component_model_source,
    build_components_model_source,
    build_mixture_model_source,
    build_model_source,
)
from pyThermoLinkDB.models import (
    ComponentModelSource,
    MixtureModelSource,
    ModelSource,
)
from pythermodb_settings.models import Component, Pressure, Temperature
from pyThermoModels.core import calc_activity_coefficient
```

Use the component-source imports only when the selected activity model needs
pure-component data in addition to mixture interaction data. In the shown
examples, NRTL uses only a mixture source, while UNIQUAC uses both mixture and
component sources.

## Shared Setup

Activity models require at least two liquid components. Each component should
include `name`, `formula`, `state`, and `mole_fraction`.

```python
ethanol = Component(
    name="ethanol",
    formula="C2H5OH",
    state="l",
    mole_fraction=0.4,
)

butyl_methyl_ether = Component(
    name="butyl-methyl-ether",
    formula="C5H12O",
    state="l",
    mole_fraction=0.6,
)

components = [ethanol, butyl_methyl_ether]

temperature = Temperature(value=323.15, unit="K")
pressure = Pressure(value=30, unit="bar")
```

The mole fractions should sum to 1.0. The core helper internally prepares:

- component ids such as `ethanol-l` or `C2H5OH-l`
- mixture ids such as `ethanol|butyl-methyl-ether`
- mole-fraction feed data keyed by component name
- pressure and temperature as `[value, unit]`

## Required Thermodynamic Inputs

Before preparing `ModelSource`, make sure the source ThermoDB data can provide
the model-specific inputs.

| Model | Required mixture data | Required component data |
| --- | --- | --- |
| `NRTL` | `alpha` and either `dg` or direct `tau` | none in the shown build-source example |
| `UNIQUAC` | interaction parameters that can generate `tau`, such as `a`, `b`, `c`, `d`, or direct `tau` | `r` and `q` for each component |

For NRTL:

- `alpha`: non-randomness parameter matrix
- `dg`: interaction energy parameter matrix used to calculate `tau`
- optional `tau`: direct binary interaction parameter matrix

For UNIQUAC:

- `a`, `b`, `c`, `d`: interaction constants used to calculate `tau`
- optional `tau`: direct binary interaction parameter matrix
- `r`: UNIQUAC volume parameter for each pure component
- `q`: UNIQUAC surface-area parameter for each pure component

The inline references must use model-source symbols that the activity models can
resolve, such as `alpha`, `dg`, `a`, `b`, `c`, `d`, `tau`, `r`, and `q`.

## NRTL: Build Model Source

The NRTL example builds a mixture ThermoDB from inline reference content and then
converts it to a `ModelSource`.

The reference must contain a mixture table with matrix symbols for:

- `alpha`: NRTL non-randomness matrix
- `dg`: NRTL interaction energy matrix

Minimal build pattern:

```python
thermodb_dir = os.path.join(parent_dir, "..", "thermodb")

thermodb_components: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=components,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)

if thermodb_components is None:
    raise ValueError("thermodb_components is None")

mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=thermodb_components,
)

model_source: ModelSource = build_model_source(
    source=[mixture_model_source],
)
```

Then calculate:

```python
res, others, G_ex = calc_activity_coefficient(
    components=components,
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name="NRTL",
    verbose=True,
)
```

Expected output:

- `res`: activity coefficient result with components, mole fractions, values,
  unit, symbol, and message
- `others`: intermediate NRTL values, including `tau_ij`, `alpha_ij`, `G_ij`,
  and component-keyed dictionaries
- `G_ex`: excess Gibbs energy result

## UNIQUAC: Build Model Source

The UNIQUAC example needs two source types:

- a mixture source for UNIQUAC interaction parameters
- component sources for pure-component `r` and `q`

The inline reference must include:

- a mixture interaction table containing symbols such as `a`, `b`, `c`, `d`
- a component/general-data table containing `r` and `q`

Minimal build pattern:

```python
thermodb_dir = os.path.join(parent_dir, "..", "thermodb")

mixture_thermodb: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=components,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)

if mixture_thermodb is None:
    raise ValueError("mixture_thermodb is None")

methanol_thermodb: ComponentThermoDB | None = build_component_thermodb_from_reference(
    component_name=methanol.name,
    component_formula=methanol.formula,
    component_state=methanol.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)

ethanol_thermodb: ComponentThermoDB | None = build_component_thermodb_from_reference(
    component_name=ethanol.name,
    component_formula=ethanol.formula,
    component_state=ethanol.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)

if methanol_thermodb is None or ethanol_thermodb is None:
    raise ValueError("component ThermoDB build failed")

mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=mixture_thermodb,
)

components_model_source: list[ComponentModelSource] = build_components_model_source(
    components_thermodb=[methanol_thermodb, ethanol_thermodb],
)

source = []
source.extend(components_model_source)
source.append(mixture_model_source)

model_source: ModelSource = build_model_source(source=source)
```

Then calculate:

```python
res, others, G_ex = calc_activity_coefficient(
    components=[ethanol, methanol],
    pressure=pressure,
    temperature=temperature,
    model_source=model_source,
    model_name="UNIQUAC",
    verbose=True,
)
```

Expected output:

- `res`: activity coefficient result with component values
- `others`: intermediate UNIQUAC values, including `tau_ij`, `r_i`, `q_i`, and
  component-keyed dictionaries
- `G_ex`: excess Gibbs energy result

## Component and Mixture Keys

`calc_activity_coefficient` creates several IDs internally:

- `Name-State`: `ethanol-l`
- `Formula-State`: `C2H5OH-l`
- `Name`: `ethanol`
- `Formula`: `C2H5OH`
- mixture by name: `ethanol|butyl-methyl-ether`
- mixture by formula: `C2H5OH|C5H12O`

Default options:

```python
component_key = "Name-State"
mixture_key = "Name"
separator_symbol = "-"
delimiter = "|"
```

Use the defaults unless the source data was built with formula-based mixture
keys. The activity helper uses `mixture_key` to select the mixture id from the
prepared source and uses component ids to select pure-component data when the
model requires it.

## Calculation Flow

The core helper performs this flow:

```text
calc_activity_coefficient(...)
  -> validate components, pressure, temperature, model_source
  -> build component ids and mixture ids
  -> build mole_fraction from Component.mole_fraction
  -> convert ModelSource to {"datasource", "equationsource"}
  -> ThermoModelCore().select_activities(...)
  -> NRTL.cal(...) or UNIQUAC.cal(...)
  -> activity_models.excess_gibbs_free_energy()
  -> return res, others, G_ex
```

## Alternative: Direct Parameter Helpers

The package also provides direct helpers that do not require `ModelSource`:

- `calc_activity_coefficient_using_nrtl_model`
- `calc_activity_coefficient_using_uniquac_model`

Use these when the agent already has matrices or dictionaries for the model
inputs.

NRTL direct helper requires:

- `tau_ij`
- `alpha_ij`

UNIQUAC direct helper requires:

- `tau_ij`
- `r_i`
- `q_i`

For the build-model-source examples requested here, prefer
`calc_activity_coefficient` because it reads the parameters from `ModelSource`.

## Common Agent Checklist

1. Define at least two liquid `Component` objects.
2. Set `mole_fraction` on every component and verify the sum is 1.0.
3. Define `Temperature` and `Pressure` with explicit units.
4. For NRTL, include mixture interaction data with `alpha` and either `dg` or
   `tau`.
5. For UNIQUAC, include mixture interaction data plus component data for `r`
   and `q`.
6. Build a `MixtureThermoDB` with `build_mixture_thermodb_from_reference`.
7. For UNIQUAC, also build each component ThermoDB with
   `build_component_thermodb_from_reference`.
8. Convert ThermoDB objects to model-source objects with
   `build_mixture_model_source` and `build_components_model_source`.
9. Combine all source objects with `build_model_source`.
10. Call `calc_activity_coefficient(..., model_name="NRTL")` or
    `calc_activity_coefficient(..., model_name="UNIQUAC")`.

## Common Failure Points

- Missing mole fractions: activity calculations require mixture composition on
  the `Component` objects.
- Single component input: activity models require at least two components.
- Missing NRTL `alpha`: NRTL needs a non-randomness parameter matrix.
- Missing NRTL interaction data: provide `dg`, `tau`, or the supported
  constants used to generate `tau`.
- Missing UNIQUAC `r` and `q`: component model sources are required when those
  values are not passed directly.
- Mixture key mismatch: the generated mixture id must match the source key, for
  example `ethanol|methanol` versus `methanol|ethanol`.
- Component order mismatch: matrices are interpreted in the component order used
  by the model.
- Import/build side effects: the examples use `thermodb_save=True`, so scripts
  may create or update files in `examples/thermodb`.

## Minimal Decision Guide

- Use `NRTL` when the source contains NRTL `alpha` and `dg` or `tau`.
- Use `UNIQUAC` when the source contains UNIQUAC interaction data and
  pure-component `r` and `q`.
- Use `calc_activity_coefficient` for build-model-source workflows.
- Use direct model-specific helpers only when the matrices are already prepared
  and no `ModelSource` is needed.
