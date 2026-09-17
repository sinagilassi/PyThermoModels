# import libs
import json
from math import exp
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
from pyThermoDB import TableMatrixData

# local
from .component_parameter_mixin import ComponentParameterMixin


class ClassicalActivityBase:
    """
    Shared helpers for small classical activity-coefficient models.
    """

    def __init__(
        self,
        components: List[str],
        datasource: Optional[Dict] = None,
        equationsource: Optional[Dict] = None,
        **kwargs
    ):
        # SECTION: validate configuration
        datasource = {} if datasource is None else datasource
        equationsource = {} if equationsource is None else equationsource

        if not isinstance(datasource, dict):
            raise TypeError("datasource must be a dict")
        if not isinstance(equationsource, dict):
            raise TypeError("equationsource must be a dict")
        if not isinstance(components, list):
            raise TypeError("components must be a list")

        # NOTE: direct model_input parameters are the Step 09 public path
        self.datasource = datasource
        self.equationsource = equationsource
        self.components = [component.strip() for component in components]
        self.comp_num = len(self.components)
        self.comp_idx = {
            self.components[i]: i for i in range(self.comp_num)
        }

        self._mixture_id: Optional[str] = kwargs.get("mixture_id", None)
        self.__mole_fraction = None
        self.__xi = None

        # SECTION: component parameter mixin
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx,
        )
        # ! set component access methods
        self.to_ij = self.component_parameter_mixin.to_ij
        self.to_i = self.component_parameter_mixin.to_i
        self.to_dict_ij = self.component_parameter_mixin.to_dict_ij
        self.to_dict_i = self.component_parameter_mixin.to_dict_i
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij
        self.to_dict_ij_ext = self.component_parameter_mixin.to_dict_ij_ext

    # SECTION: validation helpers
    def _validate_model_input(self, model_input: Dict) -> None:
        # ? does caller provide the runtime input map?
        if not isinstance(model_input, dict):
            raise TypeError("model_input must be dict")
        if "mole_fraction" not in model_input:
            raise KeyError("mole_fraction is required in model_input")

    def _mole_fraction_array(self, mole_fraction: Dict[str, float]) -> np.ndarray:
        # ? are mole fractions keyed by the configured component names?
        if not isinstance(mole_fraction, dict):
            raise TypeError("mole_fraction must be dict")

        missing = [
            component for component in self.components
            if component not in mole_fraction
        ]
        if missing:
            raise KeyError(f"Missing mole fractions for components: {missing}")

        xi = np.asarray(
            [mole_fraction[component] for component in self.components],
            dtype=float,
        )
        if np.any(xi < 0):
            raise ValueError("mole fractions must be non-negative")
        if not np.isclose(float(np.sum(xi)), 1.0, rtol=0.0, atol=1e-12):
            raise ValueError("mole fractions must sum to 1")

        # NOTE: store the latest composition for excess-Gibbs calls
        self.__xi = xi
        self.__mole_fraction = mole_fraction
        return xi

    def _latest_mole_fraction_array(
        self,
        mole_fraction: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        if mole_fraction is None:
            mole_fraction = self.__mole_fraction
        if mole_fraction is None:
            raise ValueError("mole_fraction is not set")
        return self._mole_fraction_array(mole_fraction)

    def _require_binary(self) -> None:
        # ! these classical polynomial models are implemented for binaries
        if self.comp_num != 2:
            raise ValueError("This model is implemented for binary mixtures only")

    def _matrix_parameter(
        self,
        data: TableMatrixData | np.ndarray | Dict[str, float] | List[List[float]],
        property_name: str,
        symbol_delimiter: Literal["|", "_"] = "|",
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        if isinstance(data, np.ndarray):
            mat = self.to_matrix_ij(data, symbol_delimiter=symbol_delimiter)
            comp = self.to_dict_ij(mat, symbol_delimiter=symbol_delimiter)
        elif isinstance(data, list):
            mat = self.to_matrix_ij(data, symbol_delimiter=symbol_delimiter)
            comp = self.to_dict_ij(mat, symbol_delimiter=symbol_delimiter)
        elif isinstance(data, dict):
            mat = self.to_matrix_ij(data, symbol_delimiter=symbol_delimiter)
            comp = data
        elif isinstance(data, TableMatrixData):
            mat = self.to_matrix_ij(
                data,
                property_name=property_name,
                symbol_delimiter=symbol_delimiter,
            )
            comp = self.to_dict_ij(mat, symbol_delimiter=symbol_delimiter)
        else:
            raise TypeError(
                f"{property_name} must be numpy array, list, dict, or TableMatrixData"
            )

        if mat.shape != (self.comp_num, self.comp_num):
            raise ValueError(
                f"{property_name} shape {mat.shape} does not match component number {self.comp_num}"
            )
        return mat, comp

    def _activity_result(
        self,
        gamma: np.ndarray,
        mole_fraction: np.ndarray,
        model_name: str,
        message: Optional[str],
    ) -> Dict[str, Any]:
        components_str = ", ".join(self.components)
        if message is None:
            message = (
                f"Calculate activity coefficients for {components_str} "
                f"using {model_name} model"
            )

        return {
            "property_name": "activity coefficients",
            "components": self.components,
            "mole_fraction": mole_fraction.tolist(),
            "value": [float(val) for val in gamma],
            "unit": 1,
            "symbol": "AcCo_i",
            "message": message,
        }

    def _excess_result(
        self,
        value: float,
        mole_fraction: np.ndarray,
        message: Optional[str],
        res_format: Literal["str", "json", "dict"] = "dict",
    ) -> Dict[str, Any] | str:
        components_str = ", ".join(self.components)
        if message is None:
            message = f"Excess Gibbs Free Energy for {components_str}"

        res = {
            "property_name": "Excess Molar Gibbs Free Energy (G^E/RT)",
            "components": self.components,
            "mole_fraction": mole_fraction.tolist(),
            "value": float(value),
            "unit": 1,
            "symbol": "ExMoGiFrEn",
            "message": message,
        }

        if res_format == "dict":
            return res
        if res_format in ("json", "str"):
            return json.dumps(res, indent=4)
        raise ValueError("res_format must be 'dict', 'json' or 'str'")

    def _gamma_from_lngamma(self, ln_gamma: np.ndarray) -> np.ndarray:
        return np.asarray([exp(float(val)) for val in ln_gamma], dtype=float)

    def _gibbs_identity(self, xi: np.ndarray, gamma: np.ndarray) -> float:
        return float(np.sum(xi * np.log(gamma)))
