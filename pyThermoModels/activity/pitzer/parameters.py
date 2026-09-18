"""Validated parameter container for Pitzer binary-electrolyte v1."""
from collections.abc import Mapping
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


# ! extract Pitzer binary parameters
def _extract_pitzer_binary_parameters(
        model_input: Mapping
) -> PitzerBinaryParameters:
    # SECTION: Validation
    beta0 = model_input.get("beta0")
    # >> check
    if beta0 is None:
        raise ValueError("beta0 must be provided in model_input")
    beta1 = model_input.get("beta1")
    # >> check
    if beta1 is None:
        raise ValueError("beta1 must be provided in model_input")

    c_phi = model_input.get("c_phi")
    # >> check
    if c_phi is None:
        raise ValueError("c_phi must be provided in model_input")

    # NOTE: default values for optional parameters
    alpha = model_input.get("alpha", 2.0)
    A_phi = model_input.get("A_phi", 0.3915)
    b = model_input.get("b", 1.2)

    return PitzerBinaryParameters(
        beta0=beta0,
        beta1=beta1,
        c_phi=c_phi,
        alpha=model_input.get("alpha", 2.0),
        A_phi=model_input.get("A_phi", 0.3915),
        b=model_input.get("b", 1.2),
    )
