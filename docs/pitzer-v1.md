# Pitzer v1 Activity Model

`pyThermoModels.activity.pitzer` implements a narrow, independent Pitzer
framework for one fully dissociated binary aqueous electrolyte. It uses a
molality basis and returns the thermodynamically measurable mean molal ionic
activity coefficient, osmotic coefficient, and water activity.

## Supported v1 scope

- one cation and one anion, supplied in that order;
- electroneutral dissociated-salt stoichiometry;
- at least one monovalent ion;
- single-alpha `beta0`, `beta1`, and `C_phi` interaction form;
- caller-supplied parameters at the specified temperature.

## Excluded

Mixed electrolytes, 2:2 electrolytes, beta2, theta, psi, neutral-species
terms, individual-ion conventions, speciation, parameter databases, and
caloric derivatives are intentionally unavailable in v1.

## NaCl at 298.15 K

```python
from pyThermoModels.activity.pitzer import Pitzer

model = Pitzer(["Na+", "Cl-"])
result, details = model.cal({
    "temperature": [298.15, "K"],
    "salt_molality": 1.0,
    "charges": [1, -1],
    "stoichiometry": [1, 1],
    "beta0": 0.0765,
    "beta1": 0.2664,
    "c_phi": 0.00127,
    "alpha": 2.0,
    "A_phi": 0.3915,
    "b": 1.2,
})
```

This returns `ionic_strength = 1.0 mol/kg`, `gamma_pm = 0.6555080909`,
`osmotic_coefficient = 0.9358687740`, and `water_activity = 0.9668423024`.
