# import libs
from typing import Any, Dict, List, Optional

import numpy as np


class ENRTLComponentAdapter:
    """Normalize ENRTL component identity without reparsing chemistry."""

    def __init__(
        self,
        components: List[Any],
        charge_overrides: Optional[Dict[str, int]] = None,
        require_charges: bool = True,
    ) -> None:
        """
        Adapt a list of upstream `Component` objects (or plain strings) into
        canonical ENRTL component keys and a validated charge vector.

        Parameters
        ----------
        components : List[Any]
            Component objects (e.g. `pythermodb_settings.models.Component`) or
            plain string identifiers.
        charge_overrides : Optional[Dict[str, int]]
            Legacy/manual charge values keyed by component key. Validated
            against component metadata when both are available.
        require_charges : bool
            If `True`, raise when a component has no resolvable charge;
            if `False`, default unresolved charges to `0`.

        Raises
        ------
        TypeError
            If `components` is not a list.
        """
        if not isinstance(components, list):
            raise TypeError("components must be a list")

        self.original_components = components
        self.components = components
        self.component_keys = [self._component_key(
            component) for component in components]
        self.charges = self._build_charges(
            charge_overrides=charge_overrides,
            require_charges=require_charges,
        )

    def _component_key(self, component: Any) -> str:
        """
        Resolve the canonical, charge-preserving key for a single component.

        Prefers `component.get_formula_state()`; falls back to `str(component)`
        for plain string identifiers.

        Raises
        ------
        ValueError
            If the resolved key is empty.
        """
        if hasattr(component, "get_formula_state"):
            key = component.get_formula_state()
        else:
            key = str(component).strip()

        if not key:
            raise ValueError("ENRTL component keys cannot be empty")

        return key

    def _component_charge(self, component: Any, key: str) -> Optional[int]:
        """
        Resolve a component's net charge from structured metadata.

        Prefers `component.get_net_charge()`, falls back to a plain `charge`
        attribute, and returns `None` if neither is available.
        """
        if hasattr(component, "get_net_charge"):
            return int(component.get_net_charge())

        charge = getattr(component, "charge", None)
        if charge is not None:
            return int(charge)

        return None

    def _build_charges(
        self,
        charge_overrides: Optional[Dict[str, int]] = None,
        require_charges: bool = True,
    ) -> np.ndarray:
        """
        Build the resolved charge vector, cross-checking metadata against any
        legacy override.

        Parameters
        ----------
        charge_overrides : Optional[Dict[str, int]]
            Legacy/manual charge values keyed by component key.
        require_charges : bool
            If `True`, raise for components with no resolvable charge;
            if `False`, default them to `0`.

        Returns
        -------
        np.ndarray
            Net charge per component, in `component_keys` order, dtype `int`.

        Raises
        ------
        TypeError
            If `charge_overrides` is not a dictionary.
        ValueError
            If a charge cannot be resolved (and `require_charges` is `True`),
            or if an override disagrees with component metadata.
        """
        charge_overrides = {} if charge_overrides is None else charge_overrides
        if not isinstance(charge_overrides, dict):
            raise TypeError("charge_overrides must be a dictionary")

        charges = []
        for component, key in zip(self.original_components, self.component_keys):
            metadata_charge = self._component_charge(component, key)
            override_charge = charge_overrides.get(key, None)

            # SECTION: no charge available from either source
            if metadata_charge is None and override_charge is None:
                if not require_charges:
                    charges.append(0)
                    continue
                raise ValueError(
                    f"Missing charge metadata for ENRTL component '{key}'. "
                    "Use pythermodb_settings.models.Component or provide a checked "
                    "legacy charge override."
                )

            # NOTE: cross-validate a legacy override against structured metadata
            if metadata_charge is not None and override_charge is not None:
                override_charge = int(override_charge)
                if override_charge != metadata_charge:
                    raise ValueError(
                        f"Charge override for ENRTL component '{key}' "
                        f"({override_charge}) disagrees with component metadata "
                        f"({metadata_charge})."
                    )

            resolved_charge = metadata_charge if metadata_charge is not None else override_charge
            if resolved_charge is None:
                raise ValueError(
                    f"Missing charge metadata for ENRTL component '{key}'.")

            charges.append(int(resolved_charge))

        return np.asarray(charges, dtype=int)

    @property
    def charges_dict(self) -> Dict[str, int]:
        """Return the resolved charges as a `component_key -> charge` mapping."""
        return {
            self.component_keys[i]: int(self.charges[i])
            for i in range(len(self.component_keys))
        }

    def validate_species_support(self) -> None:
        """
        Enforce the ENRTL v1 species support policy.

        Raises
        ------
        ValueError
            If a component is a radical (unsupported in v1) or a zwitterion
            (unsupported for long-range electrostatics without validation).
        """
        for component, key in zip(self.original_components, self.component_keys):
            # 1. radicals are rejected by default in ENRTL v1
            if hasattr(component, "is_radical") and component.is_radical():
                raise ValueError(
                    f"ENRTL v1 does not support radical species '{key}' without "
                    "a validated Chen-Evans parameterization."
                )

            species_type = getattr(component, "species_type", None)
            if species_type is None:
                continue

            species_values = (
                species_type
                if isinstance(species_type, list)
                else [species_type]
            )
            species_values = [str(value).lower() for value in species_values]

            # 2. zwitterions are excluded from ionic-strength/long-range treatment
            if "zwitterion" in species_values:
                raise ValueError(
                    f"ENRTL v1 does not support zwitterion long-range "
                    f"electrostatics for '{key}' without explicit validation."
                )
