from dataclasses import dataclass


# SECTION: Lee-Kesler reference-fluid identity
# NOTE: The reference fluid acentric factor is for the Lee-Kesler reference fluid.
REFERENCE_ACENTRIC_FACTOR = 0.3978


@dataclass(frozen=True)
class LeeKeslerCoefficients:
    # SECTION: Modified BWR coefficient container
    b1: float
    b2: float
    b3: float
    b4: float
    c1: float
    c2: float
    c3: float
    c4: float
    d1: float
    d2: float
    beta: float
    gamma: float


# SECTION: Simple-fluid constants
# NOTE: Lee and Kesler, AIChE J. 21(3), 510-527 (1975), DOI 10.1002/aic.690210313.
# ? Values are also reproduced in NIST teqp LKP documentation and open implementations.
SIMPLE_FLUID = LeeKeslerCoefficients(
    b1=0.1181193,
    b2=0.265728,
    b3=0.154790,
    b4=0.030323,
    c1=0.0236744,
    c2=0.0186984,
    c3=0.0,
    c4=0.042724,
    d1=0.155488e-4,
    d2=0.623689e-4,
    beta=0.65392,
    gamma=0.060167,
)

# SECTION: Reference-fluid constants
# NOTE: Same coefficient order as SIMPLE_FLUID for direct interpolation.
REFERENCE_FLUID = LeeKeslerCoefficients(
    b1=0.2026579,
    b2=0.331511,
    b3=0.027655,
    b4=0.203488,
    c1=0.0313385,
    c2=0.0503618,
    c3=0.016901,
    c4=0.041577,
    d1=0.48736e-4,
    d2=0.0740336e-4,
    beta=1.226,
    gamma=0.03754,
)
