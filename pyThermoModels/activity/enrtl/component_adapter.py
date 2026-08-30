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
        if not isinstance(components, list):
            raise TypeError("components must be a list")

        self.original_components = components
        self.components = components
        self.component_keys = [self._component_key(component) for component in components]
        self.charges = self._build_charges(
            charge_overrides=charge_overrides,
            require_charges=require_charges,
        )

    def _component_key(self, component: Any) -> str:
        if hasattr(component, "get_formula_state"):
            key = component.get_formula_state()
        else:
            key = str(component).strip()

        if not key:
            raise ValueError("ENRTL component keys cannot be empty")

        return key

    def _component_charge(self, component: Any, key: str) -> Optional[int]:
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
        charge_overrides = {} if charge_overrides is None else charge_overrides
        if not isinstance(charge_overrides, dict):
            raise TypeError("charge_overrides must be a dictionary")

        charges = []
        for component, key in zip(self.original_components, self.component_keys):
            metadata_charge = self._component_charge(component, key)
            override_charge = charge_overrides.get(key, None)

            if metadata_charge is None and override_charge is None:
                if not require_charges:
                    charges.append(0)
                    continue
                raise ValueError(
                    f"Missing charge metadata for ENRTL component '{key}'. "
                    "Use pythermodb_settings.models.Component or provide a checked "
                    "legacy charge override."
                )

            if metadata_charge is not None and override_charge is not None:
                override_charge = int(override_charge)
                if override_charge != metadata_charge:
                    raise ValueError(
                        f"Charge override for ENRTL component '{key}' "
                        f"({override_charge}) disagrees with component metadata "
                        f"({metadata_charge})."
                    )

            charges.append(metadata_charge if metadata_charge is not None else int(override_charge))

        return np.asarray(charges, dtype=int)

    @property
    def charges_dict(self) -> Dict[str, int]:
        return {
            self.component_keys[i]: int(self.charges[i])
            for i in range(len(self.component_keys))
        }

    def validate_species_support(self) -> None:
        for component, key in zip(self.original_components, self.component_keys):
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

            if "zwitterion" in species_values:
                raise ValueError(
                    f"ENRTL v1 does not support zwitterion long-range "
                    f"electrostatics for '{key}' without explicit validation."
                )
