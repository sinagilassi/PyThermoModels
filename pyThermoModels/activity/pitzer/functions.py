"""Pure mathematical functions used by multicomponent Pitzer equations."""

from math import exp, log, sqrt


def calc_g(x: float) -> float:
    """Return the standard Pitzer exponential function ``g(x)``."""
    x = float(x)
    if abs(x) < 1.0e-8:
        return 1.0 - (2.0 / 3.0) * x
    return 2.0 * (1.0 - (1.0 + x) * exp(-x)) / (x * x)


def calc_g_prime(x: float) -> float:
    """Return derivative of ``g`` with respect to its argument."""
    x = float(x)
    if abs(x) < 1.0e-6:
        return -2.0 / 3.0
    return 2.0 * ((x * x + 2.0 * x + 2.0) * exp(-x) - 2.0) / (x ** 3)


def calc_debye_huckel_f(ionic_strength: float, A_phi: float, b: float) -> float:
    """Return the Pitzer Debye-Huckel activity-coefficient term ``f(I)``."""
    if ionic_strength == 0.0:
        return 0.0
    root_i = sqrt(ionic_strength)
    return -A_phi * (root_i / (1.0 + b * root_i) + 2.0 * log(1.0 + b * root_i) / b)
