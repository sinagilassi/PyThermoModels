# individual methods
from .eos_methods import (
    check_component_eos_roots,
    check_multi_component_eos_roots,
    calc_gas_fugacity,
    calc_liquid_fugacity,
    calc_mixture_fugacity
)


__all__ = [
    'check_component_eos_roots',
    'check_multi_component_eos_roots',
    'calc_gas_fugacity',
    'calc_liquid_fugacity',
    'calc_mixture_fugacity'
]
