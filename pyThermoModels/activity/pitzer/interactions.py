"""Canonical Formula-based Pitzer interaction identities and lookup."""

from pythermodb_settings.models import Component


def make_binary_interaction_key(cation: Component, anion: Component) -> tuple[str, str]:
    """Return canonical ``(cation Formula, anion Formula)`` parameter key."""
    if not cation.is_cation() or not anion.is_anion():
        raise ValueError("Binary Pitzer interaction requires a cation then an anion")
    return cation.get_key("Formula"), anion.get_key("Formula")


def make_theta_key(component_1: Component, component_2: Component) -> tuple[str, str]:
    """Return an order-independent Formula key for unlike same-sign ions."""
    if not ((component_1.is_cation() and component_2.is_cation()) or (component_1.is_anion() and component_2.is_anion())):
        raise ValueError("Theta requires two cations or two anions")
    return tuple(sorted((component_1.get_key("Formula"), component_2.get_key("Formula"))))


def make_psi_cc_a_key(cation_1: Component, cation_2: Component, anion: Component) -> tuple[tuple[str, str], str]:
    """Return canonical key for a cation-cation-anion psi interaction."""
    if not (cation_1.is_cation() and cation_2.is_cation() and anion.is_anion()):
        raise ValueError("psi_cc_a requires two cations and one anion")
    return tuple(sorted((cation_1.get_key("Formula"), cation_2.get_key("Formula")))), anion.get_key("Formula")


def make_psi_c_aa_key(cation: Component, anion_1: Component, anion_2: Component) -> tuple[str, tuple[str, str]]:
    """Return canonical key for a cation-anion-anion psi interaction."""
    if not (cation.is_cation() and anion_1.is_anion() and anion_2.is_anion()):
        raise ValueError("psi_c_aa requires one cation and two anions")
    return cation.get_key("Formula"), tuple(sorted((anion_1.get_key("Formula"), anion_2.get_key("Formula"))))


def require_parameter(mapping: dict, key: object, family: str, strict: bool) -> object | None:
    """Resolve a parameter, never silently treating a missing value as zero."""
    if key in mapping:
        return mapping[key]
    if strict:
        raise ValueError(f"Missing required {family} Pitzer parameter for {key!r}")
    return None
