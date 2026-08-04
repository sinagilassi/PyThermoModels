# Activity Coefficient Calculation Guide for Agents

This guide summarizes how to calculate liquid-mixture activity coefficients and
binary interaction parameters (`tau_ij`) in `PyThermoModels` using the
inline-reference and build-model-source workflow
shown in:

- `examples/activity-models/calc activity using nrtl inline source and build model source - 1.py`
- `examples/activity-models/calc activity using uniquac inline source and build model source - 1.py`
- `examples/activity-models/calc activity using nrtl inline source and build model source - multi.py`
- `examples/activity-models/calc activity using uniquac inline source and build model source - multi.py`
- `examples/activity-models/calc tau_ij using nrtl inline source and build model source - 1.py`
- `examples/activity-models/calc tau_ij using uniquac inline source and build model source - 1.py`
- `examples/activity-models/calc tau_ij using nrtl inline source and build model source - multi.py`
- `examples/activity-models/calc activity using uniquac inline source and build model source - multi - 2.py`

Prefer these core helpers when writing agent-generated code:

- `calc_activity_coefficient`
- `calc_tau_ij_using_nrtl_model`
- `calc_tau_ij_using_uniquac_model`

The examples use `PyThermoDB` to build ThermoDB objects directly from inline
reference content, `PyThermoLinkDB` to build a `ModelSource`, and
`PyThermoModels` to calculate activity coefficients, `tau_ij`, and excess
Gibbs energy.

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
from pyThermoModels.core import (
    calc_activity_coefficient,
    calc_tau_ij_using_nrtl_model,
    calc_tau_ij_using_uniquac_model,
)
```

Use the component-source imports only when the selected activity model needs
pure-component data in addition to mixture interaction data. In the shown
examples, NRTL uses only a mixture source, UNIQUAC activity calculations use
both mixture and component sources, and UNIQUAC tau-only calculations use only
the mixture source.

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

For the multi-component NRTL example and the original multi-component UNIQUAC
example, the component list is:

```python
methanol = Component(
    name="methanol",
    formula="CH3OH",
    state="l",
    mole_fraction=0.30,
)

ethanol = Component(
    name="ethanol",
    formula="C2H5OH",
    state="l",
    mole_fraction=0.45,
)

butyl_methyl_ether = Component(
    name="butyl-methyl-ether",
    formula="C5H12O",
    state="l",
    mole_fraction=0.25,
)

multi_component_mixture = [methanol, ethanol, butyl_methyl_ether]
```

The requested UNIQUAC `multi - 2.py` example uses `dimethyl-carbonate` instead
of `butyl-methyl-ether`:

```python
dimethyl_carbonate = Component(
    name="dimethyl-carbonate",
    formula="C3H6O3",
    state="l",
    mole_fraction=0.25,
)

multi_component_mixture = [methanol, ethanol, dimethyl_carbonate]
```

## Required Thermodynamic Inputs

Before preparing `ModelSource`, make sure the source ThermoDB data can provide
the model-specific inputs.

| Model/workflow | Required mixture data | Required component data |
| --- | --- | --- |
| `NRTL` activity | `alpha` plus one tau route: direct `tau`, `dg`, or coefficient matrices required by `tau_correlation` | none |
| `NRTL` tau-only | one tau route: direct `tau`, `dg`, or coefficient matrices required by `tau_correlation` | none |
| `UNIQUAC` activity | one tau route: direct `tau`, `dU`, or coefficient matrices required by `tau_correlation` | `r` and `q` for each component |
| `UNIQUAC` tau-only | one tau route: direct `tau`, `dU`, or coefficient matrices required by `tau_correlation` | none |

For NRTL:

- `alpha`: non-randomness parameter matrix
- `dg`: interaction energy parameter matrix used to calculate `tau`
- `a`, `b`, `c`, `d`: coefficient matrices used to calculate `tau` when
  `dg` is not provided
- `tau`: direct binary interaction parameter matrix
- Full NRTL activity needs `alpha`; NRTL tau-only does not need `alpha`.

For UNIQUAC:

- `dU`: interaction energy parameter matrix used to calculate `tau`
- `a`, `b`, `c`, `d`: interaction constants used to calculate `tau`
- optional `tau`: direct binary interaction parameter matrix
- `r`: UNIQUAC volume parameter for each pure component
- `q`: UNIQUAC surface-area parameter for each pure component
- Full UNIQUAC activity needs `r` and `q`; UNIQUAC tau-only does not need
  `r` or `q`.

The inline references must use model-source symbols that the activity models can
resolve, such as `alpha`, `dg`, `dU`, `a`, `b`, `c`, `d`, `tau`, `r`, and `q`.

For NRTL build-model-source workflows, the source should include `alpha` plus
one interaction route: direct `tau`, `dg`, or all four coefficient matrices
`a`, `b`, `c`, `d`.

For NRTL tau-only build-model-source workflows, the source does not need
`alpha`; it only needs a valid tau route. The shown NRTL tau-only examples use
`dg` with `tau_correlation="gibbs_energy"`.

## Tau Route Equations and Required Symbols

Use this table to choose the source symbols for `tau_ij` generation. `T` is the
temperature in Kelvin and `R = 8.314`.

| Model | `tau_correlation` | Required mixture symbols | Equation |
| --- | --- | --- | --- |
| `NRTL` | `direct_tau` | `tau` or `tau_ij` | `tau_ij = tau_ij` |
| `NRTL` | `gibbs_energy` | `dg` or `dg_ij` | `tau_ij = dg_ij / (R*T)` |
| `NRTL` | `extended_temperature` | `a`, `b`, `c`, `d` | `tau_ij = a_ij + b_ij/T + c_ij*ln(T) + d_ij*T` |
| `NRTL` | `inverse_temperature` | `a`, `b` | `tau_ij = a_ij + b_ij/T` |
| `NRTL` | `inverse_temperature_squared` | `a`, `b`, `c` | `tau_ij = a_ij + b_ij/T + c_ij/T^2` |
| `NRTL` | `inverse_log_temperature` | `a`, `b`, `c` | `tau_ij = a_ij + b_ij/T + c_ij*ln(T)` |
| `UNIQUAC` | `direct_tau` | `tau` or `tau_ij` | `tau_ij = tau_ij` |
| `UNIQUAC` | `gibbs_energy` | `dU` or `dU_ij` | `tau_ij = exp(-dU_ij/(R*T))` |
| `UNIQUAC` | `extended_temperature` | `a`, `b`, `c`, `d` | `ln(tau_ij) = a_ij + b_ij/T + c_ij*ln(T) + d_ij*T`; `tau_ij = exp(ln(tau_ij))` |
| `UNIQUAC` | `inverse_temperature` | `a`, `b` | `ln(tau_ij) = a_ij + b_ij/T`; `tau_ij = exp(ln(tau_ij))` |
| `UNIQUAC` | `inverse_temperature_squared` | `a`, `b`, `c` | `ln(tau_ij) = a_ij + b_ij/T + c_ij/T^2`; `tau_ij = exp(ln(tau_ij))` |
| `UNIQUAC` | `inverse_log_temperature` | `a`, `b`, `c` | `ln(tau_ij) = a_ij + b_ij/T + c_ij*ln(T)`; `tau_ij = exp(ln(tau_ij))` |

Self-interaction terms differ by model:

- NRTL generated `tau_ii` terms are `0`.
- UNIQUAC generated `tau_ii` terms are `1`.

The parameter builders also expose energy-correlation helpers:

- NRTL `dg_ij = a_ij + b_ij*T + c_ij*T^2`.
- UNIQUAC `dU_ij = a_ij + b_ij*T + c_ij*T^2`.

Those helpers are separate from the public `tau_correlation` routes above. In
the build-model-source examples, the source usually provides `dg` or `dU`
directly instead of calculating them from `a`, `b`, and `c`.

## Two-Column Inline Matrix Pattern

In the inline ThermoDB references, each binary-pair row stores one component
row and two matrix columns for each interaction symbol. The suffixes `_i_1` and
`_i_2` mean "interaction from the current row component `i` to pair component
1 or 2".

For a binary pair `ethanol|butyl-methyl-ether`, an NRTL tau-only source using
`dg` can use only these interaction columns:

```yaml
MATRIX-SYMBOL:
  - binary interaction parameter: dg
STRUCTURE:
  COLUMNS: [No.,Mixture,Name,Formula,State,dg_i_1,dg_i_2]
  SYMBOL: [None,None,None,None,None,dg_i_1,dg_i_2]
```

For full NRTL activity, add the two `alpha` columns:

```yaml
MATRIX-SYMBOL:
  - alpha constant: alpha
  - binary interaction parameter: dg
STRUCTURE:
  COLUMNS: [No.,Mixture,Name,Formula,State,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
  SYMBOL: [None,None,None,None,None,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
```

For any coefficient route, include the required two-column pair for every
required symbol:

- `extended_temperature`: `a_i_1`, `a_i_2`, `b_i_1`, `b_i_2`, `c_i_1`,
  `c_i_2`, `d_i_1`, `d_i_2`
- `inverse_temperature`: `a_i_1`, `a_i_2`, `b_i_1`, `b_i_2`
- `inverse_temperature_squared`: `a_i_1`, `a_i_2`, `b_i_1`, `b_i_2`,
  `c_i_1`, `c_i_2`
- `inverse_log_temperature`: `a_i_1`, `a_i_2`, `b_i_1`, `b_i_2`, `c_i_1`,
  `c_i_2`

The same two-column matrix pattern applies to UNIQUAC `dU`, `tau`, and
coefficient symbols. UNIQUAC activity still needs separate component/general
data rows for `r` and `q`; UNIQUAC tau-only does not.

## Tau Correlation Names

Use public descriptive `tau_correlation` values, not raw internal method ids.
The supported names are:

- `direct_tau`: use a supplied `tau`/`tau_ij` table directly.
- `gibbs_energy`: calculate from NRTL `dg` or UNIQUAC `dU`.
- `extended_temperature`: calculate from `a`, `b`, `c`, and `d`.
- `inverse_temperature`: calculate from `a` and `b`.
- `inverse_temperature_squared`: calculate from `a`, `b`, and `c`.
- `inverse_log_temperature`: calculate from `a`, `b`, and `c`.

For NRTL, these map internally to `M0` through `M5`. For UNIQUAC, the same
public names are used, but coefficient correlations map to UNIQUAC-specific
internal method ids. Do not pass raw `M0`, `M1`, etc.

Defaults differ by helper:

- `calc_activity_coefficient`: if `tau_correlation` is omitted, the helper
  inspects the source and prefers direct `tau`, then energy data, then
  coefficient data.
- `calc_tau_ij_using_nrtl_model`: default `tau_correlation="gibbs_energy"`.
- `calc_tau_ij_using_uniquac_model`: default
  `tau_correlation="extended_temperature"`.

When the source has UNIQUAC `dU` energy data, pass
`tau_correlation="gibbs_energy"` explicitly, as shown in
`calc tau_ij using uniquac inline source and build model source - 1.py`.

## Multi-Component Inline References

For multi-component solutions, keep the `Component` list in the same order used
for the calculation and encode every binary pair needed by that mixture in the
inline reference.

For the ternary examples:

```python
multi_component_mixture = [methanol, ethanol, butyl_methyl_ether]
```

the inline mixture table stores three binary-pair groups:

- `methanol|ethanol`
- `methanol|butyl-methyl-ether`
- `ethanol|butyl-methyl-ether`

Each binary pair has two rows, one row per component in that pair. Therefore a
ternary source encoded as binary pair rows has six `VALUES` rows. A
four-component source would need all six binary pairs, each with two rows.

The `Mixture` value must match the default `mixture_key="Name"` behavior unless
the calculation is explicitly configured to use formula-based keys. With the
default delimiter, use compact pipe-separated names such as
`methanol|ethanol`; do not add spaces inside the stored mixture id unless the
source builder is configured to expect them.

For NRTL multi-component sources, each binary-pair row supplies `alpha` and one
interaction route. The multi example uses `dg`:

```yaml
MATRIX-SYMBOL:
  - alpha constant: alpha
  - binary interaction parameter: dg
STRUCTURE:
  COLUMNS: [No.,Mixture,Name,Formula,State,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
  SYMBOL: [None,None,None,None,None,alpha_i_1,alpha_i_2,dg_i_1,dg_i_2]
```

For UNIQUAC multi-component sources, the mixture table supplies one interaction
route and component tables supply pure-component `r` and `q`. The multi example
uses the coefficient route:

```yaml
MATRIX-SYMBOL:
  - a constant: a
  - b constant: b
  - c constant: c
  - d constant: d
STRUCTURE:
  COLUMNS: [No.,Mixture,Name,Formula,State,a_i_1,a_i_2,b_i_1,b_i_2,c_i_1,c_i_2,d_i_1,d_i_2]
  SYMBOL: [None,None,None,None,None,a_i_1,a_i_2,b_i_1,b_i_2,c_i_1,c_i_2,d_i_1,d_i_2]
```

The UNIQUAC `general-data` table must include one row for every component in the
multi-component mixture. At minimum, the row must resolve symbols `r` and `q`;
the example uses full general-data columns and maps
`Volume-Parameter -> r` and `Surface-Area-Parameter -> q`.

## Function Arguments in the Multi Examples

Use these argument contracts when generating or modifying the build-source
multi-component scripts.

`Component(...)`

- `name`: component name used for `Name` ids and default mixture ids.
- `formula`: formula used for `Formula` ids.
- `state`: phase/state string; activity examples use liquid state `"l"`.
- `mole_fraction`: component mole fraction in the liquid mixture. Required for
  `calc_activity_coefficient`; all component mole fractions should sum to 1.0.

`Temperature(...)`

- `value`: numeric temperature value.
- `unit`: unit string, such as `"K"`. Use explicit units.

`Pressure(...)`

- `value`: numeric pressure value.
- `unit`: unit string, such as `"bar"`. Use explicit units.

`build_mixture_thermodb_from_reference(...)`

- `components`: list of `Component` objects for the full mixture. For
  multi-component solutions, pass the full list, not only one binary pair.
- `reference_content`: inline reference string containing the mixture
  interaction table. For multi-component references encoded as binary rows,
  include every binary pair for the component set.
- `thermodb_save`: optional boolean controlling whether the generated ThermoDB
  is saved to disk. The examples set `True`.
- `thermodb_save_path`: optional output directory used when `thermodb_save=True`.
  The examples use `examples/thermodb`.
- Return value: a mixture ThermoDB wrapper or `None`. Always check for `None`
  before using `.thermodb`.

`build_component_thermodb_from_reference(...)`

- `component_name`: component name to extract from the inline component/general
  data table.
- `component_formula`: component formula to match the requested component.
- `component_state`: component state to match the requested component.
- `reference_content`: inline reference string containing component/general
  data.
- `thermodb_save`: optional boolean controlling whether the generated component
  ThermoDB is saved. The UNIQUAC multi example sets `True`.
- `thermodb_save_path`: optional output directory used when `thermodb_save=True`.
- Return value: a component ThermoDB wrapper or `None`. Always check each
  component result before building component model sources.
- Use this function for UNIQUAC multi-component workflows because each component
  must provide `r` and `q`. It is not needed in the shown NRTL multi example.

`build_mixture_model_source(...)`

- `mixture_thermodb`: mixture ThermoDB wrapper returned by
  `build_mixture_thermodb_from_reference`.
- Return value: `MixtureModelSource`, carrying mixture-level `data_source` and
  `equation_source` data.

`build_components_model_source(...)`

- `components_thermodb`: list of component ThermoDB wrappers, one for each
  component that must provide pure-component data. For the original UNIQUAC
  ternary example, pass
  `[methanol_thermodb, ethanol_thermodb, butyl_methyl_ether_thermodb]`. For
  the `multi - 2.py` example, pass
  `[methanol_thermodb, ethanol_thermodb, dimethyl_carbonate_thermodb]`.
- `rules`: optional mapping that can control source/property selection when the
  builder supports it. The activity multi examples omit it, so the builder uses
  its default matching behavior.
- Return value: list of `ComponentModelSource` objects.

`build_model_source(...)`

- `source`: list containing all model-source objects. For NRTL and tau-only
  examples, pass `[mixture_model_source]`. For UNIQUAC activity examples, pass
  all component model sources plus the mixture model source.
- Return value: `ModelSource`, which is the object passed to
  `calc_activity_coefficient`.

`calc_activity_coefficient(...)`

- Argument order:
  `components`, `pressure`, `temperature`, `model_source`, `model_name`,
  `tau_correlation`, `tau_source_priority`, `component_key`, `mixture_key`,
  `separator_symbol`, `delimiter`, `message`, `verbose`, `**kwargs`.
- `components`: list of `Component` objects in the target mixture. Keep the
  order consistent with how matrices should be interpreted.
- `pressure`: `Pressure` object. Activity models may not always use pressure
  directly, but the core helper requires it.
- `temperature`: `Temperature` object. Required when `tau_ij` must be calculated
  from `dg`, `dU`, or coefficient matrices.
- `model_source`: `ModelSource` built from the thermodb/model-source workflow.
- `model_name`: activity model name. Use `"NRTL"` for the NRTL examples and
  `"UNIQUAC"` for the UNIQUAC examples.
- `tau_correlation`: optional public tau-correlation name. Use `"direct_tau"`
  for source `tau`, `"gibbs_energy"` for NRTL `dg` or UNIQUAC `dU`, or one of
  the coefficient correlations for coefficient sources. If omitted, the helper
  selects a default from available source data.
- `tau_source_priority`: optional ordered source preference used only when
  `tau_correlation` is omitted. Default is
  `("tau", "energy", "coefficients")`; `dg` and `dU` are accepted aliases for
  `energy`.
- `component_key`: optional component id mode. Default is `"Name-State"`, which
  creates ids such as `methanol-l`. Use `"Formula-State"` only when the source
  should be resolved by formula-state ids such as `CH3OH-l`.
- `mixture_key`: optional mixture id mode. Default is `"Name"`, which creates
  mixture ids from component names. Use `"Formula"` only when the source stores
  formula-based mixture ids.
- `separator_symbol`: optional separator between component id and state.
  Default is `"-"`, producing ids such as `methanol-l`.
- `delimiter`: optional delimiter between components in mixture ids. Default is
  `"|"`, producing ids such as `methanol|ethanol`.
- `message`: optional custom result/log message. Default is `None`, so the core
  helper generates the message.
- `verbose`: optional boolean for diagnostic output. Default is `False`; the
  examples set `True`.
- `**kwargs`: optional extra keyword arguments forwarded through the core helper
  to the selected activity calculation path. Use only when a lower-level model
  option is needed.
- Return value: `(res, others, G_ex)`, where `res` contains activity
  coefficients, `others` contains model intermediates, and `G_ex` contains
  excess Gibbs energy.

`calc_tau_ij_using_nrtl_model(...)`

- Argument order:
  `components`, `temperature`, `model_source`, `tau_correlation`,
  `component_key`, `mixture_key`, `separator_symbol`, `delimiter`, `message`,
  `verbose`, `output_format`, `**kwargs`.
- `components`: list of `Component` objects in the target mixture.
- `temperature`: `Temperature` object used to evaluate calculated tau routes.
- `model_source`: `ModelSource` built from the mixture ThermoDB source.
- `tau_correlation`: default is `"gibbs_energy"`. Use `"gibbs_energy"` when
  the NRTL source contains `dg`, `"direct_tau"` when it contains `tau`, and a
  coefficient correlation when it contains compatible coefficient matrices.
- `component_key`, `mixture_key`, `separator_symbol`, and `delimiter`: same
  lookup controls as `calc_activity_coefficient`.
- `message`: optional custom log message.
- `verbose`: optional boolean for diagnostic output.
- `output_format`: component-key style for the returned dictionary. The
  examples use `"Name-State"`; the helper default is `"Name"`.
- `**kwargs`: optional lower-level options. The examples pass `mode="log"`.
- Return value: `(tau_ij_array, tau_ij_dict)`.

`calc_tau_ij_using_uniquac_model(...)`

- Argument order:
  `components`, `temperature`, `model_source`, `tau_correlation`,
  `component_key`, `mixture_key`, `separator_symbol`, `delimiter`, `message`,
  `verbose`, `output_format`, `**kwargs`.
- `components`: list of `Component` objects in the target mixture.
- `temperature`: `Temperature` object used to evaluate calculated tau routes.
- `model_source`: `ModelSource` built from the mixture ThermoDB source. The
  tau-only helper does not require UNIQUAC component sources for `r` and `q`.
- `tau_correlation`: default is `"extended_temperature"`. Pass
  `"gibbs_energy"` when the UNIQUAC source contains `dU`, `"direct_tau"` when
  it contains `tau`, and a coefficient correlation when it contains compatible
  coefficient matrices.
- `component_key`, `mixture_key`, `separator_symbol`, and `delimiter`: same
  lookup controls as `calc_activity_coefficient`.
- `message`: optional custom log message.
- `verbose`: optional boolean for diagnostic output.
- `output_format`: component-key style for the returned dictionary. The
  examples use `"Name-State"`; the helper default is `"Name"`.
- `**kwargs`: optional lower-level options. The examples pass `mode="log"`.
- Return value: `(tau_ij_array, tau_ij_dict)`.

## NRTL: Build Model Source

The NRTL example builds a mixture ThermoDB from inline reference content and then
converts it to a `ModelSource`.

The reference must contain a mixture table with matrix symbols for:

- `alpha`: NRTL non-randomness matrix
- `dg`: NRTL interaction energy matrix

Alternatively, the NRTL source can provide the temperature-dependent coefficient
set instead of `dg`:

- `a`: coefficient matrix
- `b`: coefficient matrix
- `c`: coefficient matrix
- `d`: coefficient matrix

`nrtl.py` gets or calculates `tau_ij` by one of these source routes:

- direct `tau`
- from `dg`: `tau_ij = dg_ij / (R * T)`
- from `a`, `b`, `c`, `d`:
  `tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T`

Accepted source key aliases are `alpha_ij`/`alpha`, `dg_ij`/`dg`,
`tau_ij`/`tau`, and `a_ij`/`a`, `b_ij`/`b`, `c_ij`/`c`, `d_ij`/`d`.
If both source data and `model_input` values are provided, `model_input`
overrides matching source keys.

`inputs_generator(...)` treats the literal string `"None"` the same as a
missing `tau`/`tau_ij` source, so it can still calculate `tau_ij` from `dg` or
from the full `a`/`b`/`c`/`d` coefficient set. A positive temperature is
required only for calculated `tau_ij`; direct `tau` input does not need
temperature.

Matrix sources accepted by NRTL input generation are:

- `TableMatrixData`
- list-of-lists
- `numpy.ndarray`
- component-pair dictionaries keyed with the selected delimiter, such as
  `ethanol | methanol`

This applies to direct `tau`, `alpha`, `dg`, and all four coefficient matrices
`a`, `b`, `c`, and `d`.

Minimal build pattern:

```python
thermodb_dir = os.path.join(parent_dir, "..", "thermodb")

thermodb_components: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=multi_component_mixture,
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
    components=multi_component_mixture,
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

To calculate only NRTL `tau_ij` from the same build-source pattern:

```python
tau_ij, tau_ij_comp = calc_tau_ij_using_nrtl_model(
    components=multi_component_mixture,
    temperature=temperature,
    model_source=model_source,
    tau_correlation="gibbs_energy",
    verbose=False,
    output_format="Name-State",
    mode="log",
)
```

For the NRTL tau examples, `tau_correlation="gibbs_energy"` matches the inline
`dg` source. The tau-only helper returns only the tau matrix and a
component-pair dictionary; it does not return `alpha_ij`, activity
coefficients, or excess Gibbs energy.

## UNIQUAC: Build Model Source

The UNIQUAC example needs two source types:

- a mixture source for UNIQUAC interaction parameters
- component sources for pure-component `r` and `q`

The inline reference must include:

- a mixture interaction table containing one interaction route: direct `tau`,
  `dU`, or all of `a`, `b`, `c`, `d`
- a component/general-data table containing `r` and `q`

`uniquac.py` gets or calculates `tau_ij` by one of these source routes:

- direct `tau`
- from `dU`: `tau_ij = exp(-dU_ij / (R * T))`
- from `a`, `b`, `c`, `d`:
  `tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T`

Accepted source key aliases are `r_i`/`r`, `q_i`/`q`, `dU_ij`/`dU`,
`tau_ij`/`tau`, and `a_ij`/`a`, `b_ij`/`b`, `c_ij`/`c`, `d_ij`/`d`.
If both source data and `model_input` values are provided, `model_input`
overrides matching source keys.

`inputs_generator(...)` treats the literal string `"None"` the same as a
missing `tau`/`tau_ij` source, so it can still calculate `tau_ij` from `dU` or
from the full `a`/`b`/`c`/`d` coefficient set. A positive temperature is
required only for calculated `tau_ij`; direct `tau` input does not need
temperature.

Matrix sources accepted by UNIQUAC input generation are:

- `TableMatrixData`
- list-of-lists
- `numpy.ndarray`
- component-pair dictionaries keyed with the selected delimiter, such as
  `ethanol | methanol`

This applies to direct `tau`, `dU`, and all four coefficient matrices `a`, `b`,
`c`, and `d`. Pure-component `r` and `q` can be lists, `numpy.ndarray` values,
or component dictionaries keyed by component name.

Minimal build pattern:

```python
thermodb_dir = os.path.join(parent_dir, "..", "thermodb")

mixture_thermodb: MixtureThermoDB | None = build_mixture_thermodb_from_reference(
    components=multi_component_mixture,
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

butyl_methyl_ether_thermodb: ComponentThermoDB | None = build_component_thermodb_from_reference(
    component_name=butyl_methyl_ether.name,
    component_formula=butyl_methyl_ether.formula,
    component_state=butyl_methyl_ether.state,
    reference_content=REFERENCE_CONTENT,
    thermodb_save=True,
    thermodb_save_path=thermodb_dir,
)

if (
    methanol_thermodb is None or
    ethanol_thermodb is None or
    butyl_methyl_ether_thermodb is None
):
    raise ValueError("component ThermoDB build failed")

mixture_model_source: MixtureModelSource = build_mixture_model_source(
    mixture_thermodb=mixture_thermodb,
)

components_model_source: list[ComponentModelSource] = build_components_model_source(
    components_thermodb=[
        methanol_thermodb,
        ethanol_thermodb,
        butyl_methyl_ether_thermodb,
    ],
)

source = []
source.extend(components_model_source)
source.append(mixture_model_source)

model_source: ModelSource = build_model_source(source=source)
```

Then calculate:

```python
res, others, G_ex = calc_activity_coefficient(
    components=multi_component_mixture,
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

To calculate only UNIQUAC `tau_ij` from a mixture source containing `dU`:

```python
tau_ij, tau_ij_comp = calc_tau_ij_using_uniquac_model(
    components=binary_mixture,
    temperature=temperature,
    model_source=model_source,
    tau_correlation="gibbs_energy",
    verbose=False,
    output_format="Name-State",
    mode="log",
)
```

For UNIQUAC tau-only workflows, component ThermoDB/model-source objects are not
needed unless the next step is an activity coefficient calculation. If the
source uses the coefficient route rather than `dU`, either omit
`tau_correlation` to use the UNIQUAC tau-only default
`"extended_temperature"` or pass the intended coefficient correlation
explicitly.

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

Direct `NRTL.cal(...)` also requires `mole_fraction`. If both `tau_ij` and
`alpha_ij` are already present in `model_input`, they are consumed directly as
the final matrices. Accepted shapes for these final direct inputs are numpy
arrays, list-of-lists, `TableMatrixData`, or component-pair dictionaries keyed
with the selected delimiter, such as `ethanol | methanol`.

If `tau_ij` or `alpha_ij` is missing from `model_input`, `NRTL.cal(...)` calls
`inputs_generator(...)`. `temperature` is required only when `tau_ij` must be
calculated from `dg` or coefficients. The same source rules apply: provide
`alpha`/`alpha_ij` plus `tau`/`tau_ij`, `dg`/`dg_ij`, or the complete
`a`/`b`/`c`/`d` coefficient set. Use Python `None`, omit the `tau` key, or use
the literal string `"None"` when a source contains a placeholder and tau should
be calculated instead.

UNIQUAC direct helper requires:

- `tau_ij`
- `r_i`
- `q_i`

Direct `UNIQUAC.cal(...)` also requires `mole_fraction`. If `tau_ij`, `r_i`,
or `q_i` is missing from `model_input`, `UNIQUAC.cal(...)` calls
`inputs_generator(...)`. `temperature` is required only when `tau_ij` must be
calculated from `dU` or coefficients. Use Python `None`, omit the `tau` key, or
use the literal string `"None"` when a source contains a placeholder and tau
should be calculated instead.

For the build-model-source examples requested here, prefer
`calc_activity_coefficient` because it reads the parameters from `ModelSource`.
For tau-only build-model-source examples, prefer
`calc_tau_ij_using_nrtl_model` or `calc_tau_ij_using_uniquac_model`.

## Common Agent Checklist

1. Define at least two liquid `Component` objects.
2. Set `mole_fraction` on every component and verify the sum is 1.0.
3. Define `Temperature` and `Pressure` with explicit units.
4. For NRTL build-source workflows, include mixture interaction data with
   `alpha` and one of `tau`, `dg`, or all of `a`, `b`, `c`, `d`.
5. For UNIQUAC activity coefficients, include mixture interaction data with one
   of `tau`, `dU`, or all of `a`, `b`, `c`, `d`, plus component data for `r`
   and `q`. For UNIQUAC tau-only workflows, the mixture interaction data is
   enough.
6. Build a `MixtureThermoDB` with `build_mixture_thermodb_from_reference`.
7. For UNIQUAC activity coefficients, also build each component ThermoDB with
   `build_component_thermodb_from_reference`.
8. Convert ThermoDB objects to model-source objects with
   `build_mixture_model_source` and `build_components_model_source`.
9. Combine all source objects with `build_model_source`.
10. For full activity coefficients, call
    `calc_activity_coefficient(..., model_name="NRTL")` or
    `calc_activity_coefficient(..., model_name="UNIQUAC")`.
11. For tau-only workflows, call `calc_tau_ij_using_nrtl_model(...)` or
    `calc_tau_ij_using_uniquac_model(...)` with a `tau_correlation` that
    matches the source route.

## Common Failure Points

- Missing mole fractions: activity calculations require mixture composition on
  the `Component` objects.
- Single component input: activity models require at least two components.
- Missing NRTL `alpha`: NRTL needs a non-randomness parameter matrix.
- Missing NRTL interaction data in a build-source workflow: provide direct
  `tau`, `dg`, or the full `a`, `b`, `c`, `d` coefficient set.
- Placeholder NRTL `tau`: use `None`, omit the key, or set `"None"` only when
  `dg` or all four coefficients are available and temperature is provided.
- Missing UNIQUAC `r` and `q`: component model sources are required when those
  values are not passed directly.
- Missing UNIQUAC interaction data in a build-source workflow: provide direct
  `tau`, `dU`, or the full `a`, `b`, `c`, `d` coefficient set.
- Placeholder UNIQUAC `tau`: use `None`, omit the key, or set `"None"` only
  when `dU` or all four coefficients are available and temperature is provided.
- Wrong `tau_correlation`: use `"gibbs_energy"` for NRTL `dg` and UNIQUAC
  `dU`; use `"direct_tau"` only for direct `tau`; use coefficient correlations
  only when the required coefficient matrices are present.
- Raw method ids: pass public names such as `"gibbs_energy"`, not internal
  values such as `"M1"`.
- Mixture key mismatch: the generated mixture id must match the source key, for
  example `ethanol|methanol` versus `methanol|ethanol`.
- Component order mismatch: matrices are interpreted in the component order used
  by the model.
- Import/build side effects: the examples use `thermodb_save=True`, so scripts
  may create or update files in `examples/thermodb`.

## Minimal Decision Guide

- Use `NRTL` when the source contains NRTL `alpha` and one of `tau`, `dg`, or
  all of `a`, `b`, `c`, `d`.
- Use `UNIQUAC` when the source contains UNIQUAC `r` and `q` plus one of
  `tau`, `dU`, or all of `a`, `b`, `c`, `d`.
- Use `calc_activity_coefficient` for build-model-source workflows.
- Use `calc_tau_ij_using_nrtl_model` or `calc_tau_ij_using_uniquac_model` when
  only `tau_ij` is needed from a build-model-source workflow.
- Use direct model-specific helpers only when the matrices are already prepared
  and no `ModelSource` is needed.
