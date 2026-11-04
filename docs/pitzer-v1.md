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
    Pitzer, PitzerMultiComponentBinaryParameters, PitzerParameterSet,
)

na = Component(name="Sodium", formula="Na{+}", state="aq")
cl = Component(name="Chloride", formula="Cl{-}", state="aq")
parameters = PitzerParameterSet(
    A_phi=0.3915,
    b=1.2,
    binary={
        ("Na{+}", "Cl{-}"): PitzerMultiComponentBinaryParameters(
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

Both formulations calculate properties from a shared molality-basis excess
Gibbs function. The multicomponent formulation also accepts neutral aqueous
solutes (water remains the solvent), with `lambda_ni`, `lambda_nn`, and
`zeta_g` mappings on `PitzerParameterSet`. Keys are formula tuples: neutral
then ion, sorted neutral pair (including self), and neutral/cation/anion,
respectively. Values are canonical coefficients in kg/mol for lambda and
(kg/mol)^2 for zeta; `zeta_g` is the coefficient of the Gibbs monomial, not
an unverified external database value. Strict coverage requires explicit zero
entries for omitted interactions.

Unequal-charge same-sign pairs require explicit `j_strategy="pitzer_1975"`
for the Pitzer (1975) equation 47 approximation, or
`"reduced_no_electrostatic_mixing"` for a deliberately reduced model. The
former is compared with Wolery et al. (2010), Appendix A, table A-3:
https://wipp.energy.gov/information_repository/cra/CRA-2014/References/Others/Wolery_Xiong_Long_2010_ERMS550239.pdf.
The choice must match the parameter library. Chemical speciation remains the
responsibility of a separate equilibrium layer.

Parameters supplied directly are already resolved at the requested temperature
and pressure. The current implementation does not yet resolve ModelSource
tables, validate parameter ranges or pressure, convert external conventions,
or provide validated MU/ETA adapters. Do not interpret direct input as a
verified external parameter library. Scientific mixed-electrolyte and neutral
reference comparisons remain outstanding; the code is not a validated complete
Pitzer implementation.
