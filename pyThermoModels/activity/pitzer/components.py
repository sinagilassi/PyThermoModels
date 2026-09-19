"""Component-derived state preparation for multicomponent Pitzer v2."""

from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np
from pythermodb_settings.models import Component


# SECTION: Validated numerical state
@dataclass(frozen=True)
class PitzerSpeciesState:
    """Private numerical view derived from validated shared Components."""

    component_ids: tuple[str, ...]
    formula_ids: tuple[str, ...]
    molalities: np.ndarray
    charges: np.ndarray
    cation_indices: np.ndarray
    anion_indices: np.ndarray


def normalize_pitzer_components(components: list[Component]) -> None:
    """Validate the aqueous ionic species currently supported by Pitzer v2."""
    if not isinstance(components, list) or not components:
        raise ValueError("Pitzer v2 requires a non-empty list of Components")
    ids: set[str] = set()
    for component in components:
        if not isinstance(component, Component):
            raise TypeError("Pitzer v2 components must be pythermodb_settings Component instances")
        component_id = component.get_key("Formula-State")
        if component_id in ids:
            raise ValueError(f"Duplicate Pitzer component identifier '{component_id}'")
        ids.add(component_id)
        # ! This formulation has no neutral lambda/zeta interaction terms.
        if not component.is_aqueous():
            raise ValueError("Pitzer v2 currently requires aqueous components only")
        if component.is_zwitterion():
            raise NotImplementedError("Pitzer v2 does not support zwitterion interactions")
        if component.is_radical() or component.is_radical_ion():
            raise NotImplementedError("Pitzer v2 does not support radical species")
        if component.is_neutral():
            raise NotImplementedError("Pitzer v2 neutral interactions require lambda/zeta parameters")
        if not (component.is_cation() or component.is_anion()):
            raise ValueError(f"Unsupported Pitzer ionic species '{component_id}'")


def normalize_molalities(
    components: list[Component], values: Mapping[str, float], component_key: str = "Formula-State"
) -> dict[str, float]:
    """Normalize supported Component identifiers to Formula-State molalities."""
    if not isinstance(values, Mapping):
        raise TypeError("molalities must be a mapping")
    aliases: dict[str, str] = {}
    for component in components:
        canonical = component.get_key("Formula-State")
        for key in ("Formula-State", "Name-State", "Formula", "Name", "Name-Formula", "Name-Formula-State", "Formula-Name-State"):
            aliases[component.get_key(key)] = canonical
        aliases[canonical] = canonical
    normalized: dict[str, float] = {}
    for raw_key, raw_value in values.items():
        key = str(raw_key)
        if key not in aliases:
            raise ValueError(f"Unknown Pitzer molality component '{key}'")
        canonical = aliases[key]
        if canonical in normalized:
            raise ValueError(f"Molality supplied more than once for '{canonical}'")
        value = float(raw_value)
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"Molality for '{canonical}' must be finite and non-negative")
        normalized[canonical] = value
    expected = {component.get_key("Formula-State") for component in components}
    missing = sorted(expected.difference(normalized))
    if missing:
        raise ValueError(f"Missing molalities for: {', '.join(missing)}")
    return normalized


def build_pitzer_species_state(components: list[Component], molalities: Mapping[str, float]) -> PitzerSpeciesState:
    """Build an array state after Component and molality validation."""
    normalize_pitzer_components(components)
    normalized = normalize_molalities(components, molalities)
    ids = tuple(component.get_key("Formula-State") for component in components)
    charges = np.asarray([component.get_net_charge() for component in components], dtype=int)
    values = np.asarray([normalized[key] for key in ids], dtype=float)
    return PitzerSpeciesState(
        component_ids=ids,
        formula_ids=tuple(component.get_key("Formula") for component in components),
        molalities=values,
        charges=charges,
        cation_indices=np.flatnonzero(charges > 0),
        anion_indices=np.flatnonzero(charges < 0),
    )
