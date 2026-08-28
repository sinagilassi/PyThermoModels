# ENRTL Activity Calculation Guide for Agents

This guide documents the current `PyThermoModels` ENRTL input contract. Treat
the implementation as the source of truth and keep the thermodynamic guidance
separate from user requests.

`PyThermoModels.ENRTL` is the public Electrolyte NRTL activity-coefficient
model. Version 1 is scoped to the Chen and Evans (1986) formulation name:

```python
formulation = "chen_evans_1986"
```

The model accepts a supplied true-species state. It does not perform
apparent-to-true species conversion, salt dissociation, acid/base speciation, or
reaction-equilibrium solving inside `ENRTL.cal()`.

## Current Implementation Status

Implemented infrastructure:

- component normalization through `ENRTLComponentAdapter`
- charge consumption from component metadata or checked legacy overrides
- radical and zwitterion guardrails
- `composition_representation="true_species"` validation
- mole-fraction validation and electroneutrality check
- ionic-strength calculation on an explicit basis
- long-range electrostatic contribution
- neutral molecular NRTL limiting path
- log-space contribution summation
- mean ionic activity coefficient helper
- excess Gibbs energy from supplied `ln_gamma`

Important limitation:

The charged-species Chen-Evans local-composition equations are not yet
implemented. If any component has nonzero charge and the local-composition mode
is `"chen_evans_1986"`, `ENRTL.cal()` raises `NotImplementedError`. This is
intentional; do not document the current model as production-complete for ionic
ENRTL activity coefficients.

For neutral-only mixtures, ENRTL falls back to the NRTL local-composition limit.
For ionic mixtures, the current code can validate inputs and defines the
long-range machinery, but the full calculation stops at the strict
Chen-Evans guardrail.

## Required Imports

Use structured components when possible:

```python
import numpy as np
from pythermodb_settings.models import Component

import pyThermoModels as ptm
```

Legacy string components can construct `ActivityCore.enrtl`, but a real ENRTL
calculation with `ENRTL.cal()` requires charge metadata or a checked
`model_input["charges"]` override.

## Component Inputs

Prefer `pythermodb_settings.models.Component` objects:

```python
components = [
    Component(name="water", formula="H2O", state="l"),
    Component(name="sodium", formula="Na{+}", state="aq"),
    Component(name="chloride", formula="Cl{-}", state="aq"),
]
```

ENRTL component keys are produced by `component.get_formula_state()` when that
method exists. Typical keys are:

- `H2O-l`
- `Na{+}-aq`
- `Cl{-}-aq`

Charges are read from `component.get_net_charge()` when available, then from a
`.charge` attribute if present. ENRTL does not parse formula strings to infer
charges.

If using legacy string-like components, provide a charge dictionary in
`model_input["charges"]` keyed by ENRTL component key:

```python
"charges": {
    "H2O-l": 0,
    "Na{+}-aq": 1,
    "Cl{-}-aq": -1,
}
```

If both component metadata and charge overrides are present, the values must
agree. A mismatch raises `ValueError`.

ENRTL v1 rejects radicals and zwitterions when component metadata exposes
`is_radical()` or `species_type`. Do not use `Component.is_ionic()` as the
ionic-strength criterion; ionic strength uses net charge `z_i`.

## Creating the Model

Direct model creation:

```python
enrtl = ptm.activity(
    components=components,
    model_name="ENRTL",
).enrtl
```

Equivalent lower-level construction:

```python
from pyThermoModels.activity import ENRTL

enrtl = ENRTL(components)
```

Constructor arguments:

- `components`: required list of true thermodynamic species, preferably
  `Component` objects.
- `datasource`: optional dictionary used by the ENRTL/NRTL parameter builders.
- `equationsource`: optional dictionary used by the parameter builders.
- `formulation`: optional formulation selector. The only accepted value is
  `"chen_evans_1986"`.
- `charges`: optional keyword-only legacy charge override dictionary. During
  construction missing charges default to zero only for initialization; during
  `cal()` charges are required.
- `mixture_id`: optional metadata string.
- Other `**kwargs`: forwarded to the parameter builder, matching existing
  NRTL-style infrastructure.

## Minimal `model_input`

Use this structure for the current direct-input workflow:

```python
model_input = {
    "composition_representation": "true_species",
    "temperature": [298.15, "K"],
    "mole_fraction": {
        "H2O-l": 0.98,
        "Na{+}-aq": 0.01,
        "Cl{-}-aq": 0.01,
    },
    "molality": {
        "H2O-l": 0.0,
        "Na{+}-aq": 0.1,
        "Cl{-}-aq": 0.1,
    },
    "tau_ij": np.zeros((3, 3)),
    "alpha_ij": np.zeros((3, 3)),
    "long_range": {
        "model": "pitzer_debye_huckel",
        "basis": "molality",
        "A_phi": 0.392,
    },
}
```

Required keys for `ENRTL.cal()`:

- `composition_representation`: optional in syntax, but defaults to
  `"true_species"`. Any other value raises `NotImplementedError`.
- `temperature`: required as `[value, unit]`; it is converted to kelvin with
  `pycuc`.
- `mole_fraction`: required dict, list, or numpy array. Values must be finite,
  non-negative, have one value per component, and sum to `1.0`.
- `tau_ij`: required unless it can be generated from `datasource`/
  `equationsource`.
- `alpha_ij`: required unless it can be generated from `datasource`/
  `equationsource`.
- `long_range`: required in practice because the default long-range model needs
  model-specific constants.
- `molality` or `molarity`: required when `long_range["basis"]` is `"molality"`
  or `"molarity"`.

Optional keys:

- `charges`: checked charge overrides keyed by component key.
- `activity_basis`: fallback basis used when `long_range["basis"]` is omitted.
  If both are omitted, the basis defaults to `"molality"`.
- `local_composition`: optional dictionary. The current recognized `mode`
  values are `"chen_evans_1986"` and `"neutral_nrtl_limit"`.

## Composition and Charge Rules

`mole_fraction` is always required because ENRTL validates the true-species
liquid composition and uses it for the local-composition term.

For dict inputs, keys must exactly match `enrtl.components`:

```python
["H2O-l", "Na{+}-aq", "Cl{-}-aq"]
```

For list or numpy-array inputs, values are interpreted in the same order as
`components`.

The mole-fraction composition must be electrically neutral:

```text
sum(x_i * z_i) = 0
```

If the net charge magnitude is greater than `1e-12`, `cal()` raises
`ValueError`.

Ionic strength is calculated as:

```text
I = 0.5 * sum(composition_i * z_i**2)
```

The `composition_i` vector is selected from the requested basis:

- `basis="mole_fraction"` uses `model_input["mole_fraction"]`
- `basis="molality"` uses `model_input["molality"]`
- `basis="molarity"` uses `model_input["molarity"]`

Do not call a mole-fraction-based ionic-strength value molality. Preserve the
basis in diagnostics and examples.

## Long-Range Arguments

`model_input["long_range"]` is resolved by
`ENRTLParameterBuilder.resolve_long_range()`.

Accepted keys:

- `model`: `"pitzer_debye_huckel"` or `"debye_huckel"`. Default is
  `"pitzer_debye_huckel"`.
- `basis`: `"molality"`, `"molarity"`, or `"mole_fraction"`. Default is
  `model_input["activity_basis"]`, then `"molality"`.
- `A_phi`: required for `"pitzer_debye_huckel"`.
- `b`: optional positive Pitzer-Debye-Huckel parameter. Default is `1.2`.
- `A`: required for `"debye_huckel"`.
- `B`: required for `"debye_huckel"`.
- `ion_size`: required for charged components when `model="debye_huckel"`;
  keys must be component keys.

The long-range methods return natural-log activity-coefficient contributions:
`ln(gamma_i^LR)`.

For `model="debye_huckel"`, charged components need `ion_size`:

```python
"long_range": {
    "model": "debye_huckel",
    "basis": "molality",
    "A": 0.509,
    "B": 0.328,
    "ion_size": {
        "Na{+}-aq": 4.0,
        "Cl{-}-aq": 3.5,
    },
}
```

Neutral species have zero ionic long-range contribution in the current
long-range equations.

## Local-Composition Arguments

`tau_ij` and `alpha_ij` may be supplied directly as:

- `numpy.ndarray`
- list-of-lists
- `TableMatrixData`
- component-pair dictionaries accepted by `ComponentParameterMixin.to_matrix_ij`

Matrices must have shape `(component_count, component_count)` and finite
values.

Dictionary pair keys use the selected symbol delimiter:

- with `symbol_delimiter="|"`, keys look like `A-l | B-l`
- with `symbol_delimiter="_"`, keys look like `A-l_B-l`

If either `tau_ij` or `alpha_ij` is missing, `ENRTL.cal()` calls
`inputs_generator(...)`, which currently delegates to the NRTL parameter
builder. Supported `tau_correlation` values match NRTL:

- `"direct_tau"`
- `"gibbs_energy"`
- `"extended_temperature"`
- `"inverse_temperature"`
- `"inverse_temperature_squared"`
- `"inverse_log_temperature"`

For source-generated values, provide NRTL-style source symbols such as `tau`,
`dg`, `a`, `b`, `c`, `d`, and `alpha`.

Local-composition modes:

- `"chen_evans_1986"` is the default. It raises for ionic mixtures until the
  full electrolyte-specific local-composition equations are implemented.
- `"neutral_nrtl_limit"` is valid only when every component charge is zero.

Do not force ionic ENRTL parameters into an ordinary unrestricted NRTL matrix
unless the code has been extended with the selected electrolyte formulation and
validated against that formulation.

## `cal()` Arguments

Signature:

```python
res, details = enrtl.cal(
    model_input,
    tau_correlation="gibbs_energy",
    symbol_delimiter="|",
    message=None,
    **kwargs,
)
```

Arguments:

- `model_input`: required dictionary containing composition, temperature,
  charges if needed, local-composition parameters, and long-range settings.
- `tau_correlation`: controls how missing `tau_ij` values are generated from
  sources. Default is `"gibbs_energy"`.
- `symbol_delimiter`: controls pair-dictionary parsing. Use `"|"` for keys like
  `A | B` or `"_"` for keys like `A_B`.
- `message`: optional custom result message.
- `**kwargs`: forwarded to the parameter builder when `tau_ij` or `alpha_ij`
  must be generated.

Return values:

- `res`: public result dictionary containing `property_name`, `model`,
  `formulation`, `components`, `mole_fraction`, `value`, `unit`, `symbol`, and
  `message`.
- `details`: diagnostic dictionary containing component-keyed values,
  contribution arrays, ionic strength, net charge, charge map, parameter
  matrices, `G_ij`, selected long-range model, selected local-composition mode,
  and composition representation.

Common `details` keys:

- `activity_coefficients_comp`
- `ln_gamma_total`
- `ln_gamma_total_comp`
- `ln_gamma_long_range`
- `ln_gamma_long_range_comp`
- `ln_gamma_local_composition`
- `ln_gamma_local_composition_comp`
- `ionic_strength`
- `ionic_strength_basis`
- `net_charge`
- `charges`
- `charges_comp`
- `tau_ij`
- `tau_ij_comp`
- `alpha_ij`
- `alpha_ij_comp`
- `G_ij`
- `G_ij_comp`
- `long_range_model`
- `local_composition_mode`
- `composition_representation`

## Neutral-Only Example

This path currently completes because all charges are zero:

```python
from pyThermoModels.activity import ENRTL


class FakeComponent:
    def __init__(self, key, charge):
        self.key = key
        self.charge = charge
        self.species_type = []

    def get_formula_state(self):
        return self.key

    def get_net_charge(self):
        return self.charge

    def is_radical(self):
        return False


model = ENRTL([
    FakeComponent("A-l", 0),
    FakeComponent("B-l", 0),
])

res, details = model.cal({
    "composition_representation": "true_species",
    "temperature": [298.15, "K"],
    "mole_fraction": {"A-l": 0.4, "B-l": 0.6},
    "tau_ij": [[0.0, 0.2], [0.1, 0.0]],
    "alpha_ij": [[0.0, 0.3], [0.3, 0.0]],
    "long_range": {
        "model": "pitzer_debye_huckel",
        "basis": "mole_fraction",
        "A_phi": 0.392,
    },
})
```

## Ionic True-Species Example

This input is structurally correct, but currently reaches the
Chen-Evans local-composition guardrail:

```python
components = [
    Component(name="water", formula="H2O", state="l"),
    Component(name="sodium", formula="Na{+}", state="aq"),
    Component(name="chloride", formula="Cl{-}", state="aq"),
]

activity = ptm.activity(components=components, model_name="ENRTL")
enrtl = activity.enrtl

model_input = {
    "composition_representation": "true_species",
    "temperature": [298.15, "K"],
    "mole_fraction": {
        "H2O-l": 0.98,
        "Na{+}-aq": 0.01,
        "Cl{-}-aq": 0.01,
    },
    "molality": {
        "H2O-l": 0.0,
        "Na{+}-aq": 0.1,
        "Cl{-}-aq": 0.1,
    },
    "tau_ij": np.zeros((3, 3)),
    "alpha_ij": np.zeros((3, 3)),
    "long_range": {
        "model": "pitzer_debye_huckel",
        "basis": "molality",
        "A_phi": 0.392,
    },
}

try:
    res, details = enrtl.cal(model_input=model_input)
except NotImplementedError:
    print("Ionic Chen-Evans local-composition equations are incomplete.")
```

## Mean Ionic Activity Coefficient

Use this helper after activity coefficients are available:

```python
gamma_pm = enrtl.mean_ionic_activity_coefficient(
    gamma={"Na{+}-aq": 0.8, "Cl{-}-aq": 0.9},
    cation="Na{+}-aq",
    anion="Cl{-}-aq",
    nu_cation=1,
    nu_anion=1,
)
```

Arguments:

- `gamma`: component-keyed activity coefficient dictionary.
- `cation`: component key for the cation.
- `anion`: component key for the anion.
- `nu_cation`: positive stoichiometric coefficient for the cation.
- `nu_anion`: positive stoichiometric coefficient for the anion.

The helper returns:

```text
exp((nu_cation*ln(gamma_cation) + nu_anion*ln(gamma_anion))
    / (nu_cation + nu_anion))
```

## Excess Gibbs Energy

`ENRTL.excess_gibbs_free_energy(...)` accepts:

- `mole_fraction`: dict, list, or numpy array in component order.
- `ln_gamma`: list or numpy array with one value per component.

It returns `G^E / RT`:

```text
sum(x_i * ln(gamma_i))
```

Use `details["ln_gamma_total"]` from `cal()` when available.

## Common Failure Points

- Passing `composition_representation="apparent_species"`: not implemented.
- Supplying apparent salts such as `NaCl` directly and expecting automatic
  dissociation into `Na{+}` and `Cl{-}`.
- Missing charge metadata for components.
- Charge override values that disagree with component metadata.
- Mole fractions that do not sum to `1.0`.
- Mole-fraction composition that is not electrically neutral.
- Missing `molality` or `molarity` when the long-range basis requires it.
- Missing `A_phi` for `pitzer_debye_huckel`.
- Missing `A`, `B`, or charged-component `ion_size` for `debye_huckel`.
- Expecting an ionic ENRTL calculation to complete before the
  Chen-Evans local-composition equations are implemented.
- Using `neutral_nrtl_limit` with charged species.
- Mismatching component order between `components`, list inputs, and matrices.

## Agent Checklist

1. Use true species, not apparent electrolyte names, for direct `ENRTL.cal()`.
2. Prefer `Component` objects so charges come from structured metadata.
3. Confirm every component key used in dictionaries matches
   `enrtl.components`.
4. Provide mole fractions that are finite, non-negative, sum to `1.0`, and are
   electrically neutral.
5. Choose an ionic-strength basis and provide the matching composition vector.
6. Provide long-range constants for the chosen long-range model.
7. Provide `tau_ij` and `alpha_ij` directly or ensure the source can generate
   them with the selected `tau_correlation`.
8. State the current ionic local-composition limitation whenever showing ionic
   ENRTL examples or tests.
