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


# SECTION: Multicomponent Pitzer v2 parameter containers
@dataclass(frozen=True)
class PitzerMulticomponentBinaryParameters:
    """Temperature-specific cation--anion parameters for Pitzer v2.

    ``alpha2`` is required only when ``beta2`` is non-zero.  The parameters
    are deliberately supplied by the caller; database/source resolution is a
    separate concern from the numerical Pitzer formulation.
    """

    beta0: float
    beta1: float = 0.0
    beta2: float = 0.0
    c_phi: float = 0.0
    alpha1: float = 2.0
    alpha2: float | None = None

    def __post_init__(self) -> None:
        # ! Parameters must be finite before they enter the equation kernels.
        for name in ("beta0", "beta1", "beta2", "c_phi", "alpha1"):
            value = float(getattr(self, name))
            if not isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        if self.alpha1 <= 0.0:
            raise ValueError("alpha1 must be positive")
        if self.alpha2 is not None:
            alpha2 = float(self.alpha2)
            if not isfinite(alpha2) or alpha2 <= 0.0:
                raise ValueError("alpha2 must be finite and positive")
            object.__setattr__(self, "alpha2", alpha2)
        if self.beta2 != 0.0 and self.alpha2 is None:
            raise ValueError("alpha2 must be supplied when beta2 is non-zero")


@dataclass(frozen=True)
class PitzerThetaParameter:
    """Same-sign unlike-ion Pitzer theta interaction parameter."""

    theta: float

    def __post_init__(self) -> None:
        if not isfinite(float(self.theta)):
            raise ValueError("theta must be finite")
        object.__setattr__(self, "theta", float(self.theta))


@dataclass(frozen=True)
class PitzerPsiParameter:
    """Three-ion Pitzer psi interaction parameter."""

    psi: float

    def __post_init__(self) -> None:
        if not isfinite(float(self.psi)):
            raise ValueError("psi must be finite")
        object.__setattr__(self, "psi", float(self.psi))


@dataclass(frozen=True)
class PitzerParameterSet:
    """Explicit parameter source for one multicomponent Pitzer calculation."""

    A_phi: float
    b: float
    binary: Mapping
    theta: Mapping
    psi_cc_a: Mapping
    psi_c_aa: Mapping

    def __post_init__(self) -> None:
        # NOTE: Mappings use Formula-based canonical interaction identities.
        for name in ("A_phi", "b"):
            value = float(getattr(self, name))
            if not isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)
        for name in ("binary", "theta", "psi_cc_a", "psi_c_aa"):
            if not isinstance(getattr(self, name), Mapping):
                raise TypeError(f"{name} must be a mapping")
