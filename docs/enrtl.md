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
- Chen-Evans local-composition interaction classification
- charged true-species local-composition calculations

Formulation rule:

The ionic Chen-Evans local-composition contribution is implemented in
`pyThermoModels/activity/enrtl/local_composition.py`. It uses explicit
interaction classes, excludes like-ion local-composition terms, keeps ordinary
NRTL only as the all-neutral limiting path, and exposes local-composition
diagnostics in the `other_values["local_composition_diagnostics"]` result.

Validation boundary:

The charged Chen-Evans calculation path is now available for true-species
systems with supplied parameters. Quantitative production use still requires
case-specific validation against published electrolyte activity-coefficient
data and parameter sets.
