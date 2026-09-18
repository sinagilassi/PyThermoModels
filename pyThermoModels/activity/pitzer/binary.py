"""Single-alpha binary Pitzer equations on the molality basis."""

from math import exp, isfinite, log, sqrt

from .parameters import PitzerBinaryParameters


# SECTION: Debye-Huckel and binary interaction terms
def calc_f_phi(ionic_strength: float, params: PitzerBinaryParameters) -> float:
    """Return the Pitzer Debye-Huckel osmotic term ``f_phi``."""
    I = _validate_ionic_strength(ionic_strength)
    if I == 0.0:
        return 0.0
    root_I = sqrt(I)
    return -params.A_phi * root_I / (1.0 + params.b * root_I)


def calc_B_phi(ionic_strength: float, params: PitzerBinaryParameters) -> float:
    """Return the single-alpha osmotic binary interaction ``B_phi``."""
    I = _validate_ionic_strength(ionic_strength)
    return params.beta0 + params.beta1 * exp(-params.alpha * sqrt(I))


def calc_f_gamma(ionic_strength: float, params: PitzerBinaryParameters) -> float:
    """Return the Pitzer Debye-Huckel mean-activity term ``f_gamma``."""
    I = _validate_ionic_strength(ionic_strength)
    if I == 0.0:
        return 0.0
    root_I = sqrt(I)
    bracket = root_I / (1.0 + params.b * root_I) + 2.0 * \
        log(1.0 + params.b * root_I) / params.b
    return -params.A_phi * bracket


def calc_B_gamma(ionic_strength: float, params: PitzerBinaryParameters) -> float:
    """Return the single-alpha mean-activity binary interaction ``B_gamma``."""
    I = _validate_ionic_strength(ionic_strength)
    if I == 0.0:
        # ? Retain this exact limit to avoid the direct expression's division by I.
        # NOTE: limiting form of the provided single-alpha expression.
        return 2.0 * params.beta0 + params.beta1
    root_I = sqrt(I)
    alpha_root_I = params.alpha * root_I
    bracket = 1.0 - (1.0 + alpha_root_I - 0.5 *
                     alpha_root_I**2) * exp(-alpha_root_I)
    return 2.0 * params.beta0 + 2.0 * params.beta1 * bracket / (params.alpha**2 * I)


def calc_C_gamma(params: PitzerBinaryParameters) -> float:
    """Return ``C_gamma = 3/2 C_phi`` for the v1 binary formulation."""
    return 1.5 * params.c_phi


# SECTION: Binary-electrolyte properties
def calc_osmotic_coefficient_binary(
    salt_molality: float,
    ionic_strength: float,
    z_cation: int,
    z_anion: int,
    nu_cation: int,
    nu_anion: int,
    params: PitzerBinaryParameters,
) -> float:
    """Return the binary Pitzer osmotic coefficient on the molality basis."""
    m = _validate_molality(salt_molality)
    if m == 0.0:
        return 1.0
    I = _validate_ionic_strength(ionic_strength)
    charge_product = abs(int(z_cation) * int(z_anion))
    binary_factor, c_factor = _stoichiometric_factors(nu_cation, nu_anion)
    return 1.0 + charge_product * calc_f_phi(I, params) + m * binary_factor * calc_B_phi(I, params) + m**2 * c_factor * params.c_phi


def calc_ln_gamma_mean_binary(
    salt_molality: float,
    ionic_strength: float,
    z_cation: int,
    z_anion: int,
    nu_cation: int,
    nu_anion: int,
    params: PitzerBinaryParameters,
) -> float:
    """Return ``ln(gamma_pm)``; no individual-ion convention is applied."""
    m = _validate_molality(salt_molality)
    if m == 0.0:
        return 0.0
    I = _validate_ionic_strength(ionic_strength)
    charge_product = abs(int(z_cation) * int(z_anion))
    binary_factor, c_factor = _stoichiometric_factors(nu_cation, nu_anion)
    return charge_product * calc_f_gamma(I, params) + m * binary_factor * calc_B_gamma(I, params) + m**2 * c_factor * calc_C_gamma(params)


def calc_gamma_mean_binary(
    salt_molality: float,
    ionic_strength: float,
    z_cation: int,
    z_anion: int,
    nu_cation: int,
    nu_anion: int,
    params: PitzerBinaryParameters,
) -> float:
    """Return the thermodynamically measurable binary mean ionic coefficient."""
    return exp(
        calc_ln_gamma_mean_binary(
            salt_molality=salt_molality,
            ionic_strength=ionic_strength,
            z_cation=z_cation,
            z_anion=z_anion,
            nu_cation=nu_cation,
            nu_anion=nu_anion,
            params=params,
        )
    )


def calc_water_activity(
    salt_molality: float,
    osmotic_coefficient: float,
    nu_cation: int,
    nu_anion: int,
    water_molar_mass: float = 0.01801528,
) -> float:
    """Return solvent water activity from the Pitzer osmotic coefficient."""
    m = _validate_molality(salt_molality)
    phi = float(osmotic_coefficient)
    M_w = float(water_molar_mass)
    if not isfinite(phi) or not isfinite(M_w) or M_w <= 0.0:
        raise ValueError(
            "osmotic_coefficient and water_molar_mass must be finite; water_molar_mass must be positive")
    if nu_cation <= 0 or nu_anion <= 0:
        raise ValueError("stoichiometric coefficients must be positive")
    return exp(-M_w * (nu_cation + nu_anion) * m * phi)


def _stoichiometric_factors(nu_cation: int, nu_anion: int) -> tuple[float, float]:
    if nu_cation <= 0 or nu_anion <= 0:
        raise ValueError("stoichiometric coefficients must be positive")
    nu = nu_cation + nu_anion
    product = nu_cation * nu_anion
    return 2.0 * product / nu, 2.0 * product**1.5 / nu


def _validate_molality(value: float) -> float:
    m = float(value)
    if not isfinite(m) or m < 0.0:
        raise ValueError("salt_molality must be finite and non-negative")
    return m


def _validate_ionic_strength(value: float) -> float:
    I = float(value)
    if not isfinite(I) or I < 0.0:
        raise ValueError("ionic_strength must be finite and non-negative")
    return I
