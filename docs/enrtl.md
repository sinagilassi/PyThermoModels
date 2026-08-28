# Electrolyte NRTL

`PyThermoModels.ENRTL` is reserved for the Chen and Evans (1986)
Electrolyte NRTL formulation:

Chau-Chyun Chen and L. B. Evans, "A local composition model for the
excess Gibbs energy of aqueous electrolyte systems", AIChE Journal,
32(3), 444-454, 1986. DOI: `10.1002/aic.690320311`.

The v1 implementation accepts true-species states and keeps apparent-to-true
mapping, dissociation, and reaction/speciation solving outside `ENRTL.cal()`.

Supported infrastructure:

- component identity from `pythermodb_settings.models.Component`
- charge metadata consumption without formula reparsing
- electroneutrality checks
- ionic strength with explicit basis
- isolated long-range electrostatic contribution
- log-space contribution summation
- mean ionic activity coefficient helper
- neutral molecular NRTL limiting path

Strict formulation rule:

The ionic Chen-Evans local-composition equations must be implemented in
`enrtl_local_composition.py` before ionic ENRTL calculations are treated as
production-complete. The current code fails explicitly for ionic
local-composition calculations instead of substituting ordinary NRTL.
