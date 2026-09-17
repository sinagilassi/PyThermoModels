"""Validated parameter container for Pitzer binary-electrolyte v1."""

from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class PitzerBinaryParameters:
    """Temperature-specific parameters for one binary Pitzer electrolyte."""

    beta0: float
    beta1: float
    c_phi: float
    alpha: float = 2.0
    A_phi: float = 0.3915
    b: float = 1.2

    def __post_init__(self) -> None:
        # SECTION: Validate caller-supplied, temperature-specific parameters.
        for name in ("beta0", "beta1", "c_phi", "alpha", "A_phi", "b"):
            value = float(getattr(self, name))
            if not isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)

        # ! The electrostatic and single-alpha parameters must be positive.
        for name in ("alpha", "A_phi", "b"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"{name} must be positive")
