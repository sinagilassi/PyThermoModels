from math import exp, log
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
import pycuc
from pyThermoDB import TableMatrixData
# ! locals
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from ...utils.utility import TauCorrelation
from ..component_parameter_mixin import ComponentParameterMixin
from .component_adapter import ENRTLComponentAdapter
from .core import ActivityBasis, CompositionRepresentation, ENRTLCore
from .local_composition import ENRTLFormulation, ENRTLLocalComposition
from .long_range import ENRTLLongRange
from .parameter_builder import ENRTLParameterBuilder


class ENRTL(ENRTLCore):
    """Electrolyte NRTL activity-coefficient model.

    Version 1 is formulated for Chen and Evans (1986) and accepts a supplied
    true-species state. It intentionally keeps speciation and apparent-to-true
    mapping outside this model.
    """

    def __init__(
        self,
        components: List[Any],
        datasource: Optional[Dict] = None,
        equationsource: Optional[Dict] = None,
        formulation: ENRTLFormulation = "chen_evans_1986",
        **kwargs: Any,
    ) -> None:
        if formulation != "chen_evans_1986":
            raise ValueError(
                "ENRTL v1 supports only formulation 'chen_evans_1986'")

        datasource = {} if datasource is None else datasource
        equationsource = {} if equationsource is None else equationsource
        if not isinstance(datasource, dict):
            raise TypeError("datasource must be a dict")
        if not isinstance(equationsource, dict):
            raise TypeError("equationsource must be a dict")

        adapter = ENRTLComponentAdapter(
            components=components,
            charge_overrides=kwargs.get("charges"),
            require_charges=False,
        )

        self.original_components = components
        self.components = adapter.component_keys
        self.comp_num = len(self.components)
        self.comp_idx = {
            self.components[i]: i
            for i in range(self.comp_num)
        }
        self.datasource = datasource
        self.equationsource = equationsource
        self.formulation = formulation
        self._mixture_id = kwargs.get("mixture_id", "")

        self.enrtl_parameter_builder = ENRTLParameterBuilder(
            components=self.components,
            comp_idx=self.comp_idx,
            datasource=self.datasource,
            equationsource=self.equationsource,
            **kwargs,
        )
        self.inputs_generator = self.enrtl_parameter_builder.inputs_generator

        self.long_range = ENRTLLongRange()
        self.local_composition = ENRTLLocalComposition(
            components=self.components,
            comp_idx=self.comp_idx,
            formulation=formulation,
        )
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx,
        )
        self.to_dict_ij = self.component_parameter_mixin.to_dict_ij
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij
        self.to_dict_i = self.component_parameter_mixin.to_dict_i

    @add_attributes(metadata=ACTIVITY_MODELS["ENRTL"])
    def cal(
        self,
        model_input: Dict[str, Any],
        tau_correlation: TauCorrelation = "gibbs_energy",
        symbol_delimiter: Literal["|", "_"] = "|",
        message: Optional[str] = None,
        **kwargs: Any,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        inputs = self.check_and_build_inputs(
            model_input=model_input,
            tau_correlation=tau_correlation,
            symbol_delimiter=symbol_delimiter,
            **kwargs,
        )

        x = inputs["mole_fraction_array"]
        charges = inputs["charges"]
        charge_map = inputs["charges_comp"]
        long_range_inputs = inputs["long_range"]
        basis = long_range_inputs["basis"]

        ionic_composition = inputs["ionic_strength_composition"]
        ionic_strength = self.ionic_strength(
            composition=ionic_composition,
            charges=charges,
            basis=basis,
        )

        ln_gamma_lr_comp = self.long_range.cal_ln_gamma_long_range(
            ionic_strength=ionic_strength,
            charges=charge_map,
            model=long_range_inputs["model"],
            ion_size=long_range_inputs["ion_size"],
            A=long_range_inputs["A"],
            B=long_range_inputs["B"],
            A_phi=long_range_inputs["A_phi"],
            b=long_range_inputs["b"],
        )
        ln_gamma_lr = np.asarray(
            [ln_gamma_lr_comp[component] for component in self.components],
            dtype=float,
        )

        local_mode = inputs["local_composition"].get("mode", "chen_evans_1986")
        ln_gamma_lc = self.local_composition.cal_ln_gamma_lc(
            mole_fraction=x,
            charges=charges,
            tau_ij=inputs["tau_ij"],
            alpha_ij=inputs["alpha_ij"],
            mode=local_mode,
        )

        ln_gamma_total = self.combine_ln_gamma(
            ln_gamma_long_range=ln_gamma_lr,
            ln_gamma_local_composition=ln_gamma_lc,
        )
        gamma = np.exp(ln_gamma_total)
        if not np.all(np.isfinite(gamma)) or np.any(gamma <= 0):
            raise ValueError(
                "ENRTL activity coefficients must be finite and positive")

        gamma_list = [float(value) for value in gamma]
        gamma_comp = {
            self.components[i]: gamma_list[i]
            for i in range(self.comp_num)
        }
        ln_gamma_lc_comp = {
            self.components[i]: float(ln_gamma_lc[i])
            for i in range(self.comp_num)
        }
        ln_gamma_total_comp = {
            self.components[i]: float(ln_gamma_total[i])
            for i in range(self.comp_num)
        }
        G_ij, G_ij_comp = self.local_composition.build_local_composition_diagnostics(
            tau_ij=inputs["tau_ij"],
            alpha_ij=inputs["alpha_ij"],
            symbol_delimiter=symbol_delimiter,
        )

        message = (
            "Calculate activity coefficients using ENRTL "
            f"({self.formulation})"
            if message is None
            else message
        )

        res = {
            "property_name": "activity coefficients",
            "model": "ENRTL",
            "formulation": self.formulation,
            "components": self.components,
            "mole_fraction": x.tolist(),
            "value": gamma_list,
            "unit": 1,
            "symbol": "AcCo_i",
            "message": message,
        }
        other_values = {
            "activity_coefficients_comp": gamma_comp,
            "ln_gamma_total": ln_gamma_total,
            "ln_gamma_total_comp": ln_gamma_total_comp,
            "ln_gamma_long_range": ln_gamma_lr,
            "ln_gamma_long_range_comp": ln_gamma_lr_comp,
            "ln_gamma_local_composition": ln_gamma_lc,
            "ln_gamma_local_composition_comp": ln_gamma_lc_comp,
            "ionic_strength": ionic_strength,
            "ionic_strength_basis": basis,
            "net_charge": self.net_charge(composition=x, charges=charges),
            "charges": charges,
            "charges_comp": charge_map,
            "tau_ij": inputs["tau_ij"],
            "tau_ij_comp": self.to_dict_ij(inputs["tau_ij"], symbol_delimiter=symbol_delimiter),
            "alpha_ij": inputs["alpha_ij"],
            "alpha_ij_comp": self.to_dict_ij(inputs["alpha_ij"], symbol_delimiter=symbol_delimiter),
            "G_ij": G_ij,
            "G_ij_comp": G_ij_comp,
            "long_range_model": long_range_inputs["model"],
            "local_composition_mode": local_mode,
            "local_composition_diagnostics": self.local_composition.last_diagnostics,
            "composition_representation": inputs["composition_representation"],
        }
        return res, other_values

    def check_and_build_inputs(
        self,
        model_input: Dict[str, Any],
        tau_correlation: TauCorrelation = "gibbs_energy",
        symbol_delimiter: Literal["|", "_"] = "|",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        if not isinstance(model_input, dict):
            raise TypeError("model_input must be dict")

        representation: CompositionRepresentation = model_input.get(
            "composition_representation",
            "true_species",
        )
        if representation != "true_species":
            raise NotImplementedError(
                "ENRTL.cal() v1 accepts only composition_representation='true_species'. "
                "Apparent-to-true mapping belongs in a separate speciation layer."
            )

        adapter = ENRTLComponentAdapter(
            components=self.original_components,
            charge_overrides=model_input.get("charges"),
            require_charges=True,
        )
        adapter.validate_species_support()
        charges = adapter.charges
        charge_map = {
            self.components[i]: int(charges[i])
            for i in range(self.comp_num)
        }

        mole_fraction = model_input.get("mole_fraction")
        # >> check
        if mole_fraction is None:
            raise ValueError("model_input must contain 'mole_fraction'")

        x = self._composition_to_array(mole_fraction, "mole_fraction")
        if abs(float(np.sum(x)) - 1.0) > self.COMPOSITION_ATOL:
            raise ValueError("ENRTL mole fractions must sum to 1.0")
        self.check_electroneutrality(composition=x, charges=charges)

        temperature = self._validate_temperature(
            model_input.get("temperature"))

        tau_ij, alpha_ij = self._build_tau_alpha(
            model_input=model_input,
            tau_correlation=tau_correlation,
            symbol_delimiter=symbol_delimiter,
            **kwargs,
        )

        long_range = self.enrtl_parameter_builder.resolve_long_range(
            model_input)
        basis: ActivityBasis = long_range["basis"]
        ionic_strength_composition = self._ionic_strength_composition(
            model_input=model_input,
            basis=basis,
        )

        local_composition = model_input.get("local_composition", {})
        if local_composition is None:
            local_composition = {}
        if not isinstance(local_composition, dict):
            raise TypeError(
                "model_input['local_composition'] must be a dictionary")

        return {
            "mole_fraction": mole_fraction,
            "mole_fraction_array": x,
            "temperature": temperature,
            "charges": charges,
            "charges_comp": charge_map,
            "tau_ij": tau_ij,
            "alpha_ij": alpha_ij,
            "long_range": long_range,
            "local_composition": local_composition,
            "ionic_strength_composition": ionic_strength_composition,
            "composition_representation": representation,
        }

    def mean_ionic_activity_coefficient(
        self,
        gamma: Dict[str, float],
        cation: str,
        anion: str,
        nu_cation: int = 1,
        nu_anion: int = 1,
    ) -> float:
        if nu_cation <= 0 or nu_anion <= 0:
            raise ValueError("stoichiometric coefficients must be positive")
        if cation not in gamma or anion not in gamma:
            raise KeyError("cation and anion must both exist in gamma")

        gamma_cation = float(gamma[cation])
        gamma_anion = float(gamma[anion])
        if gamma_cation <= 0 or gamma_anion <= 0:
            raise ValueError("activity coefficients must be positive")

        ln_gamma_mean = (
            nu_cation * log(gamma_cation)
            + nu_anion * log(gamma_anion)
        ) / (nu_cation + nu_anion)
        return float(exp(ln_gamma_mean))

    def excess_gibbs_free_energy(
        self,
        mole_fraction: Dict[str, float] | List[float] | np.ndarray,
        ln_gamma: List[float] | np.ndarray,
    ) -> Dict[str, Any]:
        x = self._composition_to_array(mole_fraction, "mole_fraction")
        ln_gamma_array = np.asarray(ln_gamma, dtype=float)
        if ln_gamma_array.shape != (self.comp_num,):
            raise ValueError(f"ln_gamma must have shape ({self.comp_num},)")

        gE_RT = float(np.sum(x * ln_gamma_array))
        return {
            "property_name": "Excess Molar Gibbs Free Energy (G^E/RT)",
            "model": "ENRTL",
            "formulation": self.formulation,
            "components": self.components,
            "mole_fraction": x.tolist(),
            "value": gE_RT,
            "unit": 1,
            "symbol": "ExMoGiFrEn",
        }

    def _build_tau_alpha(
        self,
        model_input: Dict[str, Any],
        tau_correlation: TauCorrelation,
        symbol_delimiter: Literal["|", "_"],
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if "tau_ij" in model_input:
            tau_ij = self._matrix_to_array(
                model_input["tau_ij"], "tau_ij", symbol_delimiter)
        else:
            generated = self.inputs_generator(
                temperature=model_input.get("temperature"),
                tau_correlation=tau_correlation,
                symbol_delimiter=symbol_delimiter,
                model_input=model_input,
                **kwargs,
            )
            tau_ij = np.asarray(generated["tau_ij"], dtype=float)

        if "alpha_ij" in model_input:
            alpha_ij = self._matrix_to_array(
                model_input["alpha_ij"],
                "alpha_ij",
                symbol_delimiter,
            )
        else:
            generated = self.inputs_generator(
                temperature=model_input.get("temperature"),
                tau_correlation=tau_correlation,
                symbol_delimiter=symbol_delimiter,
                model_input=model_input,
                **kwargs,
            )
            alpha_ij = np.asarray(generated["alpha_ij"], dtype=float)

        return tau_ij, alpha_ij

    def _matrix_to_array(
        self,
        data: TableMatrixData | np.ndarray | Dict[str, float] | List[List[float]],
        property_name: str,
        symbol_delimiter: Literal["|", "_"],
    ) -> np.ndarray:
        if isinstance(data, np.ndarray):
            matrix = data.astype(float)
        elif isinstance(data, list):
            matrix = np.asarray(data, dtype=float)
        else:
            matrix = self.to_matrix_ij(
                data=data,
                property_name=property_name.replace("_ij", ""),
                symbol_delimiter=symbol_delimiter,
            )

        if matrix.shape != (self.comp_num, self.comp_num):
            raise ValueError(
                f"{property_name} must have shape ({self.comp_num}, {self.comp_num})")
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{property_name} contains non-finite values")
        return matrix

    def _composition_to_array(
        self,
        composition: Dict[str, float] | List[float] | np.ndarray,
        name: str,
    ) -> np.ndarray:
        if isinstance(composition, dict):
            missing = [
                component for component in self.components if component not in composition]
            if missing:
                raise KeyError(
                    f"Missing {name} value for component(s): {', '.join(missing)}")
            values = [composition[component] for component in self.components]
            array = np.asarray(values, dtype=float)
        elif isinstance(composition, list):
            array = np.asarray(composition, dtype=float)
        elif isinstance(composition, np.ndarray):
            array = composition.astype(float)
        else:
            raise TypeError(f"{name} must be a dict, list, or numpy array")

        if array.shape != (self.comp_num,):
            raise ValueError(f"{name} must have shape ({self.comp_num},)")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} contains non-finite values")
        if np.any(array < 0):
            raise ValueError(f"{name} values must be non-negative")
        return array

    def _ionic_strength_composition(
        self,
        model_input: Dict[str, Any],
        basis: ActivityBasis,
    ) -> np.ndarray:
        if basis == "mole_fraction":
            return self._composition_to_array(
                model_input["mole_fraction"],
                "mole_fraction",
            )

        if basis not in ("molality", "molarity"):
            raise ValueError(
                "long_range basis must be 'molality', 'molarity', or 'mole_fraction'"
            )

        if basis not in model_input:
            raise KeyError(
                f"model_input['{basis}'] is required for ionic strength")

        return self._composition_to_array(model_input[basis], basis)

    def _validate_temperature(self, temperature: List[float | str] | None) -> float:
        if temperature is None:
            raise KeyError("temperature is required in model_input")
        if not isinstance(temperature, list) or len(temperature) != 2:
            raise TypeError(
                "temperature must be a list formatted as [value, unit]")

        value = float(temperature[0])
        unit = str(temperature[1])
        temperature_k = pycuc.convert_from_to(value, unit, "K")
        if temperature_k is None or temperature_k <= 0:
            raise ValueError("temperature must be greater than 0 K")
        return float(temperature_k)
