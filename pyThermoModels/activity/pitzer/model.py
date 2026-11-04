"""Pitzer v1 model for one fully dissociated binary aqueous electrolyte."""

from typing import Any, Dict, List, Optional, Tuple, cast
from numbers import Real
from math import isfinite

import numpy as np
import pycuc
from pythermodb_settings.models import Component, Temperature
# locals
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from .binary import calc_water_activity
from .core import (
    build_binary_ion_molalities,
    calc_ionic_strength,
    calc_net_charge,
    validate_binary_electrolyte_v1,
    validate_electroneutrality,
)
from .parameters import PitzerBinaryParameters, _extract_pitzer_binary_parameters
from .parameters import PitzerParameterSet
from .parameters import PitzerMultiComponentBinaryParameters
from .gibbs import evaluate_gibbs, excess_gibbs_molal, gibbs_gradient
from .components import normalize_molalities, normalize_pitzer_components
from .multicomponent import (
    calc_Z,
    calc_multicomponent_ionic_strength,
    calc_multicomponent_net_charge,
    validate_multicomponent_electroneutrality,
)


def _is_finite_integer(value: Any) -> bool:
    """Return whether a numeric input represents a finite, non-boolean integer."""
    return bool(isinstance(value, Real) and not isinstance(value, (bool, np.bool_))
                and isfinite(value) and value % 1 == 0)


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
        formulation: str = "binary_single_alpha_v1",
        **kwargs: Any,
    ) -> None:
        # SECTION: Formulation-specific component validation.
        if formulation not in ("binary_single_alpha_v1", "multicomponent_pitzer_v2"):
            raise ValueError(f"Unsupported Pitzer formulation '{formulation}'")
        self._species_components: list[Component] = []
        self.components: list[str] | list[Component]
        if formulation == "binary_single_alpha_v1":
            if not isinstance(components, list) or len(components) != 2:
                raise ValueError("Pitzer v1 requires exactly [cation, anion] components")
            self.components = [self._component_key(component) for component in components]
        else:
            normalize_pitzer_components(components)
            self._species_components = cast(list[Component], components)
            self.components = self._species_components
        self.formulation = formulation
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

        if self.formulation == "multicomponent_pitzer_v2":
            return self._cal_multicomponent_v2(model_input)

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
        if any(not _is_finite_integer(charge) for charge in charges):
            raise ValueError("charges must be finite integers")
        charges = [int(charge) for charge in charges]

        # ? stoichiometry
        stoichiometry = model_input.get("stoichiometry")
        if stoichiometry is None or not isinstance(stoichiometry, (list, tuple)) or len(stoichiometry) != 2:
            raise ValueError(
                "stoichiometry must be provided as a list or tuple of length 2 in model_input")
        if any(not _is_finite_integer(nu) for nu in stoichiometry):
            raise ValueError("stoichiometry must contain finite integers")
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

        resolved = PitzerParameterSet(params.A_phi, params.b, {}, {}, {}, {})
        binary = PitzerMultiComponentBinaryParameters(params.beta0, params.beta1, c_phi=params.c_phi, alpha1=params.alpha)
        terms = ([(0, 1, binary)], [], [], [], [])
        if salt_molality == 0:
            phi, ln_gamma_mean, gibbs = 1.0, 0.0, 0.0
        else:
            kernel = lambda m: excess_gibbs_molal(m, charges, resolved, terms)
            gibbs = float(kernel(ion_molalities))
            derivatives = gibbs_gradient(kernel, ion_molalities)
            phi = 1 + (float(np.dot(ion_molalities, derivatives)) - gibbs) / sum(ion_molalities)
            ln_gamma_mean = float(np.dot(stoichiometry, derivatives) / sum(stoichiometry))
        gamma_mean = float(np.exp(ln_gamma_mean))
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
            "gibbs_excess_per_kg_rt": gibbs,
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

    # SECTION: Component-centric multicomponent Pitzer v2 orchestration
    def _cal_multicomponent_v2(
        self, model_input: Dict[str, Any]
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Calculate ionic multicomponent Pitzer properties from true species."""
        components = self._species_components
        temperature = model_input.get("temperature")
        if not isinstance(temperature, Temperature):
            raise TypeError("Pitzer v2 temperature must be a pythermodb_settings Temperature")
        molalities = model_input.get("molalities")
        if molalities is None:
            raise ValueError("Pitzer v2 requires molalities in model_input")
        parameters = model_input.get("parameters")
        if model_input.get("model_source") is not None:
            raise NotImplementedError("Pitzer ModelSource parameter resolution is not implemented; supply a resolved PitzerParameterSet")
        if model_input.get("pressure") is not None:
            raise NotImplementedError("Pitzer pressure-dependent parameter resolution is not implemented")
        if not isinstance(parameters, PitzerParameterSet):
            raise TypeError("Pitzer v2 parameters must be a PitzerParameterSet")
        strict_parameters = model_input.get("strict_parameters", True)
        if not isinstance(strict_parameters, bool):
            raise TypeError("strict_parameters must be a boolean")
        water_molar_mass = float(model_input.get("water_molar_mass", 0.01801528))
        if not np.isfinite(water_molar_mass) or water_molar_mass <= 0.0:
            raise ValueError("water_molar_mass must be finite and positive")

        # NOTE: Temperature remains metadata while kernels use kelvin.
        temperature_K = pycuc.convert_from_to(temperature.value, temperature.unit, "K")
        if not np.isfinite(temperature_K) or temperature_K <= 0:
            raise ValueError("temperature must be finite and positive")
        normalized = normalize_molalities(components, molalities)
        validate_multicomponent_electroneutrality(components, normalized)
        ionic_strength = calc_multicomponent_ionic_strength(components, normalized)
        Z = calc_Z(components, normalized)
        gibbs, derivatives, phi, ln_aw = evaluate_gibbs(components, normalized, parameters, strict_parameters, water_molar_mass)
        ln_gamma = {c.get_key("Formula-State"): float(v) for c, v in zip(components, derivatives)}
        gamma = {key: float(np.exp(value)) for key, value in ln_gamma.items()}
        water_activity = float(np.exp(ln_aw))
        component_ids = [component.get_key("Formula-State") for component in components]

        res = {
            "property_name": "activity coefficients",
            "model": "PITZER",
            "formulation": self.formulation,
            "components": component_ids,
            "value": gamma,
            "unit": 1,
            "symbol": "gamma_i",
        }
        other_values = {
            "temperature": temperature.model_dump(),
            "temperature_K": temperature_K,
            "molalities": normalized,
            "ionic_strength": ionic_strength,
            "ionic_strength_unit": "mol/kg",
            "Z": Z,
            "net_charge": calc_multicomponent_net_charge(components, normalized),
            "ln_activity_coefficients": ln_gamma,
            "activity_coefficients": gamma,
            "osmotic_coefficient": phi,
            "water_activity": water_activity,
            "ln_water_activity": ln_aw,
            "gibbs_excess_per_kg_rt": gibbs,
            "parameter_provenance": parameters.provenance,
            "warnings": (["Missing interactions may have been treated as zero because strict_parameters=False"] if not strict_parameters else []),
            "water_molar_mass": water_molar_mass,
            "strict_parameters": strict_parameters,
            "parameters": parameters,
            "component_metadata": {
                component.get_key("Formula-State"): {
                    "formula": component.formula,
                    "name": component.name,
                    "state": component.state,
                    "charge": component.get_net_charge(),
                    "species_type": component.species_type,
                } for component in components
            },
        }
        if "stoichiometry" in model_input:
            stoichiometry = model_input["stoichiometry"]
            if not isinstance(stoichiometry, dict) or not stoichiometry:
                raise ValueError("stoichiometry must be a nonempty component-to-integer mapping")
            component_by_id = {c.get_key("Formula-State"): c for c in components}
            if any(key not in component_by_id or not _is_finite_integer(nu) or nu <= 0 for key, nu in stoichiometry.items()):
                raise ValueError("stoichiometry must contain known components and positive integers")
            if any(component_by_id[key].get_net_charge() == 0 for key in stoichiometry):
                raise ValueError("mean ionic stoichiometry must contain ions only")
            if sum(component_by_id[k].get_net_charge()*v for k, v in stoichiometry.items()) != 0:
                raise ValueError("stoichiometry must be electroneutral")
            ln_mean = sum(ln_gamma[k]*v for k, v in stoichiometry.items())/sum(stoichiometry.values())
            other_values["ln_gamma_mean"] = ln_mean
            other_values["gamma_mean"] = float(np.exp(ln_mean))
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
