"""Same-sign Pitzer mixing functions for multicomponent ionic systems."""

from math import exp, sqrt


def calc_J(x: float) -> float:
    """Return the accepted Pitzer approximation to the electrostatic J function."""
    x = abs(float(x))
    if x == 0.0:
        return 0.0
    denominator = 4.0 + 4.581 * x ** -0.7237 * exp(-0.0120 * x ** 0.528)
    return x / denominator


def calc_J_prime(x: float) -> float:
    """Return a stable central-difference derivative of the J approximation."""
    x = abs(float(x))
    step = max(1.0e-6, x * 1.0e-5)
    return (calc_J(x + step) - calc_J(max(0.0, x - step))) / (2.0 * step)


def calc_etheta(ionic_strength: float, z_i: int, z_j: int, A_phi: float) -> float:
    """Return electrostatic unlike same-sign contribution ``Etheta``."""
    if ionic_strength == 0.0 or z_i == z_j:
        return 0.0
    root_i = sqrt(ionic_strength)
    x_ij = 6.0 * abs(z_i * z_j) * A_phi * root_i
    x_ii = 6.0 * z_i * z_i * A_phi * root_i
    x_jj = 6.0 * z_j * z_j * A_phi * root_i
    return z_i * z_j * (calc_J(x_ij) - 0.5 * calc_J(x_ii) - 0.5 * calc_J(x_jj)) / (4.0 * ionic_strength)


def calc_etheta_prime(ionic_strength: float, z_i: int, z_j: int, A_phi: float) -> float:
    """Return numerical derivative of ``Etheta`` with respect to ionic strength."""
    if ionic_strength == 0.0 or z_i == z_j:
        return 0.0
    step = max(1.0e-7, ionic_strength * 1.0e-5)
    return (calc_etheta(ionic_strength + step, z_i, z_j, A_phi) - calc_etheta(max(0.0, ionic_strength - step), z_i, z_j, A_phi)) / (2.0 * step)


def calc_Phi(theta: float, ionic_strength: float, z_i: int, z_j: int, A_phi: float) -> float:
    """Return same-sign interaction ``Phi = theta + Etheta``."""
    return theta + calc_etheta(ionic_strength, z_i, z_j, A_phi)


def calc_Phi_prime(ionic_strength: float, z_i: int, z_j: int, A_phi: float) -> float:
    """Return the ionic-strength derivative of Phi."""
    return calc_etheta_prime(ionic_strength, z_i, z_j, A_phi)


def calc_Phi_phi(theta: float, ionic_strength: float, z_i: int, z_j: int, A_phi: float) -> float:
    """Return same-sign contribution used in the osmotic coefficient."""
    return theta + calc_etheta(ionic_strength, z_i, z_j, A_phi) + ionic_strength * calc_etheta_prime(ionic_strength, z_i, z_j, A_phi)
