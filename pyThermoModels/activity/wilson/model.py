# import libs
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
from pyThermoDB import TableMatrixData

# local
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from ..classical_base import ClassicalActivityBase


class Wilson(ClassicalActivityBase):
    """
    Wilson activity-coefficient model.

    NOTE: The standard Wilson model cannot represent liquid-liquid phase
    splitting.
    """

    @add_attributes(metadata=ACTIVITY_MODELS["WILSON"])
    def cal(
        self,
        model_input: Dict,
        symbol_delimiter: Literal["|", "_"] = "|",
        message: Optional[str] = None,
        **kwargs
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        try:
            # SECTION: validate and unpack inputs
            self._validate_model_input(model_input)
            xi = self._mole_fraction_array(model_input["mole_fraction"])

            # ? does caller provide Wilson Lambda parameters?
            if "lambda_ij" not in model_input:
                raise KeyError("lambda_ij is required in model_input")

            lambda_ij, lambda_ij_comp = self._matrix_parameter(
                model_input["lambda_ij"],
                property_name="lambda",
                symbol_delimiter=symbol_delimiter,
            )

            # ! Wilson Lambda parameters must be positive for logarithms
            if np.any(lambda_ij <= 0):
                raise ValueError("lambda_ij values must be positive")
            np.fill_diagonal(lambda_ij, 1.0)
            lambda_ij_comp = self.to_dict_ij(
                lambda_ij,
                symbol_delimiter=symbol_delimiter,
            )

            # SECTION: calculate activity coefficients
            ln_gamma = self.CalLnAcCo_V1(xi=xi, lambda_ij=lambda_ij)
            gamma = self._gamma_from_lngamma(ln_gamma)
            AcCo_i_comp = {
                self.components[i]: float(gamma[i])
                for i in range(self.comp_num)
            }

            # SECTION: prepare result
            res = self._activity_result(
                gamma=gamma,
                mole_fraction=xi,
                model_name="Wilson",
                message=message,
            )
            other_values = {
                "AcCo_i_comp": AcCo_i_comp,
                "ln_gamma": ln_gamma,
                "lambda_ij": lambda_ij,
                "lambda_ij_comp": lambda_ij_comp,
                "excess_gibbs_RT": self.CalExcessGibbs_RT(
                    xi=xi,
                    lambda_ij=lambda_ij,
                ),
            }
            return res, other_values
        except Exception as e:
            raise Exception(f"Error in Wilson model cal: {str(e)}")

    def CalLnAcCo_V1(
        self,
        xi: np.ndarray,
        lambda_ij: np.ndarray,
    ) -> np.ndarray:
        try:
            # SECTION: Wilson multicomponent expression
            sums = lambda_ij @ xi
            if np.any(sums <= 0):
                raise ValueError("Wilson logarithm arguments must be positive")

            ln_gamma = np.zeros(self.comp_num)
            for i in range(self.comp_num):
                # NOTE: ln(gamma_i) = 1 - ln(S_i) - sum_j x_j Lambda_ji / S_j
                ln_gamma[i] = 1.0 - np.log(sums[i]) - np.sum(
                    xi * lambda_ij[:, i] / sums
                )
            return ln_gamma
        except Exception as e:
            raise Exception(f"Error in Wilson CalLnAcCo_V1: {str(e)}")

    def CalExcessGibbs_RT(
        self,
        xi: np.ndarray,
        lambda_ij: np.ndarray,
    ) -> float:
        try:
            # SECTION: G^E/RT expression
            sums = lambda_ij @ xi
            if np.any(sums <= 0):
                raise ValueError("Wilson logarithm arguments must be positive")
            return float(-np.sum(xi * np.log(sums)))
        except Exception as e:
            raise Exception(f"Error in Wilson CalExcessGibbs_RT: {str(e)}")

    def excess_gibbs_free_energy(
        self,
        mole_fraction: Optional[Dict[str, float]] = None,
        lambda_ij: Optional[
            TableMatrixData | np.ndarray | Dict[str, float] | List[List[float]]
        ] = None,
        symbol_delimiter: Literal["|", "_"] = "|",
        message: Optional[str] = None,
        res_format: Literal["str", "json", "dict"] = "dict",
    ) -> Dict[str, Any] | str:
        try:
            # SECTION: validate inputs
            xi = self._latest_mole_fraction_array(mole_fraction)
            if lambda_ij is None:
                raise ValueError("lambda_ij is required")
            lambda_ij, _ = self._matrix_parameter(
                lambda_ij,
                property_name="lambda",
                symbol_delimiter=symbol_delimiter,
            )
            if np.any(lambda_ij <= 0):
                raise ValueError("lambda_ij values must be positive")
            np.fill_diagonal(lambda_ij, 1.0)

            # SECTION: calculate excess Gibbs energy
            gE_RT = self.CalExcessGibbs_RT(xi=xi, lambda_ij=lambda_ij)
            return self._excess_result(
                value=gE_RT,
                mole_fraction=xi,
                message=message,
                res_format=res_format,
            )
        except Exception as e:
            raise Exception(f"Error in Wilson excess_gibbs_free_energy: {str(e)}")
