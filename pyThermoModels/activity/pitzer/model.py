"""Pitzer v1 model for one fully dissociated binary aqueous electrolyte."""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
# locals
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from .binary import (
    calc_gamma_mean_binary,
    calc_ln_gamma_mean_binary,
    calc_osmotic_coefficient_binary,
    calc_water_activity,
)
from .core import (
    build_binary_ion_molalities,
    calc_ionic_strength,
    calc_net_charge,
    validate_binary_electrolyte_v1,
    validate_electroneutrality,
)
from .parameters import PitzerBinaryParameters, _extract_pitzer_binary_parameters


class Pitzer:
    """Binary, single-alpha Pitzer electrolyte model on a molality basis.

    Version 1 reports only the measurable mean molal ionic activity coefficient.
    It does not define individual-ion activity coefficients or speciation.
    """

    formulation = "binary_single_alpha_v1"

    def __init__(
        self,
        components: List[Any],
        datasource: Optional[Dict[str, Any]] = None,
        equationsource: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> None:
        if not isinstance(components, list) or len(components) != 2:
            raise ValueError(
                "Pitzer v1 requires exactly [cation, anion] components")
        self.components = [self._component_key(
            component) for component in components]
        self.datasource = {} if datasource is None else datasource
        self.equationsource = {} if equationsource is None else equationsource

    @add_attributes(metadata=ACTIVITY_MODELS["PITZER"])
    def cal(
        self,
        model_input: Dict[str, Any]
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Calculate binary Pitzer mean activity, osmotic coefficient, and water activity."""
        if not isinstance(model_input, dict):
            raise TypeError("model_input must be a dictionary")

        # SECTION: Build the validated binary ionic state.
        salt_molality = model_input.get("salt_molality")
        if salt_molality is None:
            raise ValueError("salt_molality must be provided in model_input")
        # >>> set
        salt_molality = float(salt_molality)

        # ? charges
        charges = model_input.get("charges")
        if charges is None or not isinstance(charges, (list, tuple)) or len(charges) != 2:
            raise ValueError(
                "charges must be provided as a list or tuple of length 2 in model_input")
        charges = [int(charge) for charge in charges]

        # ? stoichiometry
        stoichiometry = model_input.get("stoichiometry")
        if stoichiometry is None or not isinstance(stoichiometry, (list, tuple)) or len(stoichiometry) != 2:
            raise ValueError(
                "stoichiometry must be provided as a list or tuple of length 2 in model_input")
        stoichiometry = [int(nu) for nu in stoichiometry]

        charges, stoichiometry = validate_binary_electrolyte_v1(
            charges=charges,
            stoichiometry=stoichiometry,
        )
        ion_molalities = build_binary_ion_molalities(
            salt_molality=salt_molality,
            nu_cation=int(stoichiometry[0]),
            nu_anion=int(stoichiometry[1]),
        )
        validate_electroneutrality(ion_molalities, charges)
        ionic_strength = calc_ionic_strength(ion_molalities, charges)

        # SECTION: Parameters are direct, caller-supplied values at this temperature.
        # NOTE: extract parameters
        params: PitzerBinaryParameters = _extract_pitzer_binary_parameters(
            model_input=model_input
        )

        phi = calc_osmotic_coefficient_binary(
            salt_molality=salt_molality,
            ionic_strength=ionic_strength,
            z_cation=int(charges[0]),
            z_anion=int(charges[1]),
            nu_cation=int(stoichiometry[0]),
            nu_anion=int(stoichiometry[1]),
            params=params,
        )
        ln_gamma_mean = calc_ln_gamma_mean_binary(
            salt_molality=salt_molality,
            ionic_strength=ionic_strength,
            z_cation=int(charges[0]),
            z_anion=int(charges[1]),
            nu_cation=int(stoichiometry[0]),
            nu_anion=int(stoichiometry[1]),
            params=params,
        )
        gamma_mean = calc_gamma_mean_binary(
            salt_molality=salt_molality,
            ionic_strength=ionic_strength,
            z_cation=int(charges[0]),
            z_anion=int(charges[1]),
            nu_cation=int(stoichiometry[0]),
            nu_anion=int(stoichiometry[1]),
            params=params,
        )
        water_activity = calc_water_activity(
            salt_molality=salt_molality,
            osmotic_coefficient=phi,
            nu_cation=int(stoichiometry[0]),
            nu_anion=int(stoichiometry[1]),
            water_molar_mass=float(model_input.get(
                "water_molar_mass", 0.01801528)),
        )

        # ! gamma_mean is not an individual single-ion activity coefficient.
        res = {
            "property_name": "mean ionic activity coefficient",
            "model": "PITZER",
            "formulation": self.formulation,
            "components": self.components,
            "value": gamma_mean,
            "unit": 1,
            "symbol": "gamma_pm",
        }
        other_values = {
            "salt_molality": salt_molality,
            "molality_unit": "mol/kg",
            "ion_molalities": ion_molalities.tolist(),
            "charges": charges.tolist(),
            "stoichiometry": stoichiometry.tolist(),
            "ionic_strength": ionic_strength,
            "ionic_strength_unit": "mol/kg",
            "net_charge": calc_net_charge(ion_molalities, charges),
            "ln_gamma_mean": ln_gamma_mean,
            "gamma_mean": gamma_mean,
            "osmotic_coefficient": phi,
            "water_activity": water_activity,
            "water_molar_mass": float(model_input.get("water_molar_mass", 0.01801528)),
            "parameter_temperature": model_input.get("temperature"),
            "parameters": {
                "beta0": params.beta0,
                "beta1": params.beta1,
                "c_phi": params.c_phi,
                "alpha": params.alpha,
                "A_phi": params.A_phi,
                "b": params.b,
            },
        }
        return res, other_values

    @staticmethod
    def _component_key(component: Any) -> str:
        if hasattr(component, "get_formula_state"):
            key = str(component.get_formula_state())
        else:
            key = str(component).strip()
        if not key:
            raise ValueError("Pitzer component keys cannot be empty")
        return key
