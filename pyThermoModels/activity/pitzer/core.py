"""Molality and charge bookkeeping for Pitzer binary-electrolyte v1."""

from numpy.typing import ArrayLike

import numpy as np


CHARGE_BALANCE_ATOL = 1e-12


# SECTION: Binary electrolyte state construction
def build_binary_ion_molalities(
    salt_molality: float,
    nu_cation: int,
    nu_anion: int,
) -> np.ndarray:
    """Return cation/anion molalities from an analytical salt molality."""
    m = float(salt_molality)
    if not np.isfinite(m) or m < 0.0:
        raise ValueError("salt_molality must be finite and non-negative")
    if not isinstance(nu_cation, int) or not isinstance(nu_anion, int):
        raise TypeError("stoichiometric coefficients must be integers")
    if nu_cation <= 0 or nu_anion <= 0:
        raise ValueError("stoichiometric coefficients must be positive")
    return np.asarray([nu_cation * m, nu_anion * m], dtype=float)


def calc_net_charge(molalities: ArrayLike, charges: ArrayLike) -> float:
    """Return the molality-weighted ionic charge for a binary ionic state."""
    m, z = _validate_vector_pair(molalities, charges)
    return float(np.dot(m, z))


def validate_electroneutrality(
    molalities: ArrayLike,
    charges: ArrayLike,
    atol: float = CHARGE_BALANCE_ATOL,
) -> bool:
    """Require the molal ionic state to be electrically neutral."""
    charge = calc_net_charge(molalities, charges)
    if abs(charge) > atol:
        raise ValueError(
            "Pitzer binary electrolyte is not electrically neutral: "
            f"net charge = {charge:.6e}"
        )
    return True


def calc_ionic_strength(molalities: ArrayLike, charges: ArrayLike) -> float:
    """Return molal ionic strength ``I = 0.5 sum(m_i z_i^2)``."""
    m, z = _validate_vector_pair(molalities, charges)
    return float(0.5 * np.sum(m * np.square(z)))


def validate_binary_electrolyte_v1(
    charges: ArrayLike,
    stoichiometry: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate the intentionally narrow Pitzer v1 binary-ion scope."""
    raw_charges = np.asarray(charges, dtype=float)
    raw_stoichiometry = np.asarray(stoichiometry, dtype=float)
    if raw_charges.shape != (2,) or raw_stoichiometry.shape != (2,):
        raise ValueError("Pitzer v1 requires exactly one cation and one anion")
    if not np.all(np.isfinite(raw_charges)) or not np.all(np.isfinite(raw_stoichiometry)):
        raise ValueError("Pitzer v1 charges and stoichiometry must be finite")
    if not np.all(raw_charges == np.floor(raw_charges)):
        raise ValueError("Pitzer v1 charges must be integers")
    if not np.all(raw_stoichiometry == np.floor(raw_stoichiometry)):
        raise ValueError("Pitzer v1 stoichiometry must contain integers")
    z = raw_charges.astype(int)
    nu = raw_stoichiometry.astype(int)
    if z[0] <= 0 or z[1] >= 0:
        raise ValueError("Pitzer v1 requires charges ordered as [cation, anion]")
    if np.any(nu <= 0):
        raise ValueError("Pitzer v1 stoichiometry must contain positive integers")
    if abs(int(z[0])) != z[0] or abs(int(z[1])) != -z[1]:
        raise ValueError("Pitzer v1 charges must be nonzero integers")
    if abs(int(z[0])) == 2 and abs(int(z[1])) == 2:
        raise NotImplementedError("Pitzer v1 does not support 2:2 electrolytes")
    if abs(int(z[0])) != 1 and abs(int(z[1])) != 1:
        raise NotImplementedError("Pitzer v1 requires at least one monovalent ion")

    # NOTE: stoichiometry and charge must define a neutral dissociated salt.
    ion_molalities = nu.astype(float)
    validate_electroneutrality(ion_molalities, z)
    return z, nu


def _validate_vector_pair(
    molalities: ArrayLike,
    charges: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    m = np.asarray(molalities, dtype=float)
    z = np.asarray(charges, dtype=float)
    if m.shape != (2,) or z.shape != (2,):
        raise ValueError("Pitzer v1 requires exactly two ionic species")
    if not np.all(np.isfinite(m)) or np.any(m < 0.0):
        raise ValueError("molalities must be finite and non-negative")
    if not np.all(np.isfinite(z)):
        raise ValueError("charges must be finite")
    return m, z
