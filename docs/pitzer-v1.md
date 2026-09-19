# Pitzer Activity Model

`pyThermoModels.activity.pitzer` contains two explicit formulations.

- `binary_single_alpha_v1` is the default, compatibility-preserving binary
  salt API. It accepts the original salt molality, charge, and stoichiometry
  inputs and returns the measurable mean ionic coefficient.
- `multicomponent_pitzer_v2` uses shared `Component` metadata and true aqueous
  species molalities. It returns individual-ion Pitzer coefficients, osmotic
  coefficient, and water activity.

## Component-centric multicomponent v2

V2 requires aqueous cation/anion `Component` objects, a unit-bearing shared
`Temperature`, and an explicit `PitzerParameterSet`. Charges and composition
identifiers are derived from `Component`; callers must not provide a separate
charge map. Temperatures are converted to kelvin internally.

```python
from pythermodb_settings.models import Component, Temperature
from pyThermoModels.activity.pitzer import (
    Pitzer, PitzerMulticomponentBinaryParameters, PitzerParameterSet,
)

na = Component(name="Sodium", formula="Na{+}", state="aq")
cl = Component(name="Chloride", formula="Cl{-}", state="aq")
parameters = PitzerParameterSet(
    A_phi=0.3915,
    b=1.2,
    binary={
        ("Na{+}", "Cl{-}"): PitzerMulticomponentBinaryParameters(
            beta0=0.0765, beta1=0.2664, c_phi=0.00127,
        ),
    },
    theta={}, psi_cc_a={}, psi_c_aa={},
)
model = Pitzer([na, cl], formulation="multicomponent_pitzer_v2")
result, details = model.cal({
    "temperature": Temperature(value=25.0, unit="C"),
    "molalities": {"Na{+}-aq": 1.0, "Cl{-}-aq": 1.0},
    "parameters": parameters,
})
```

Interaction keys use Formula identifiers. Same-sign theta keys and the
same-sign pair inside psi keys are canonicalized, so their input order does not
change parameter lookup. Parameter coverage is strict by default: active
interactions must be represented explicitly, including interactions assigned a
zero value.

V2 is an ionic formulation. Neutral `lambda`/`zeta` interactions, radicals,
and zwitterions are intentionally rejected. Chemical dissociation and
speciation remain the responsibility of the reaction/equilibrium layer; Pitzer
consumes the resulting true-species molalities.
