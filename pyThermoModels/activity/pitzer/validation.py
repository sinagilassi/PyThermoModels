"""Public validation helpers for component-centric Pitzer v2 inputs."""

from collections.abc import Mapping

from pythermodb_settings.models import Component

from .components import normalize_molalities, normalize_pitzer_components
from .multicomponent import validate_multicomponent_electroneutrality


# SECTION: Pitzer v2 input validation facade
def validate_components(components: list[Component]) -> bool:
    """Validate that Components are unique, aqueous, and supported ionic species."""
    normalize_pitzer_components(components)
    return True


def validate_molalities(components: list[Component], molalities: Mapping[str, float]) -> dict[str, float]:
    """Validate and normalize user molalities to Formula-State identifiers."""
    # NOTE: This returns the canonical mapping used by all numerical kernels.
    return normalize_molalities(components, molalities)


def validate_component_molality_mapping(components: list[Component], molalities: Mapping[str, float]) -> bool:
    """Require complete normalized molalities and an electroneutral ionic state."""
    normalized = normalize_molalities(components, molalities)
    validate_multicomponent_electroneutrality(components, normalized)
    return True
