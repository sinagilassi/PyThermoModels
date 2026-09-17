# import libs
from typing import Any, Dict, Literal, Optional, Tuple

import numpy as np

# local
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from ..classical_base import ClassicalActivityBase


class Margules(ClassicalActivityBase):
    """
    Binary two-suffix and three-suffix Margules model.
    """

    @add_attributes(metadata=ACTIVITY_MODELS["MARGULES"])
    def cal(
        self,
        model_input: Dict,
        calculation_mode: Literal["two_suffix", "three_suffix"] = "three_suffix",
        message: Optional[str] = None,
        **kwargs
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        try:
            # SECTION: validate and unpack inputs
            self._require_binary()
            self._validate_model_input(model_input)
            xi = self._mole_fraction_array(model_input["mole_fraction"])
            x1 = float(xi[0])

            # SECTION: calculate model values
            if calculation_mode == "two_suffix":
                # ? does caller provide the symmetric Margules parameter?
                if "A" not in model_input:
                    raise KeyError("A is required for two_suffix Margules")
                params = {"A": float(model_input["A"])}
                gE_RT, dg_dx1 = self._two_suffix(x1=x1, A=params["A"])
            elif calculation_mode == "three_suffix":
                # ? does caller provide ordered A12 and A21 parameters?
                if "A12" not in model_input or "A21" not in model_input:
                    raise KeyError(
                        "A12 and A21 are required for three_suffix Margules"
                    )
                params = {
                    "A12": float(model_input["A12"]),
                    "A21": float(model_input["A21"]),
                }
                gE_RT, dg_dx1 = self._three_suffix(
                    x1=x1,
                    A12=params["A12"],
                    A21=params["A21"],
                )
            else:
                raise ValueError("calculation_mode must be two_suffix or three_suffix")

            ln_gamma = self._binary_lngamma(gE_RT=gE_RT, dg_dx1=dg_dx1, xi=xi)
            gamma = self._gamma_from_lngamma(ln_gamma)
            AcCo_i_comp = {
                self.components[i]: float(gamma[i])
                for i in range(self.comp_num)
            }

            # SECTION: prepare result
            res = self._activity_result(
                gamma=gamma,
                mole_fraction=xi,
                model_name="Margules",
                message=message,
            )
            other_values = {
                "AcCo_i_comp": AcCo_i_comp,
                "ln_gamma": ln_gamma,
                "parameters": params,
                "calculation_mode": calculation_mode,
                "excess_gibbs_RT": gE_RT,
            }
            return res, other_values
        except Exception as e:
            raise Exception(f"Error in Margules model cal: {str(e)}")

    def _two_suffix(self, x1: float, A: float) -> Tuple[float, float]:
        # NOTE: x2 is component 2 in configured component order
        x2 = 1.0 - x1
        gE_RT = A * x1 * x2
        dg_dx1 = A * (1.0 - 2.0 * x1)
        return float(gE_RT), float(dg_dx1)

    def _three_suffix(self, x1: float, A12: float, A21: float) -> Tuple[float, float]:
        # NOTE: convention is G^E/RT = x1*x2*(A12*x1 + A21*x2)
        x2 = 1.0 - x1
        bracket = A12 * x1 + A21 * x2
        gE_RT = x1 * x2 * bracket
        dg_dx1 = (1.0 - 2.0 * x1) * bracket + x1 * x2 * (A12 - A21)
        return float(gE_RT), float(dg_dx1)

    def _binary_lngamma(
        self,
        gE_RT: float,
        dg_dx1: float,
        xi: np.ndarray,
    ) -> np.ndarray:
        # SECTION: binary derivative relation
        ln_gamma = np.zeros(2)
        ln_gamma[0] = gE_RT + xi[1] * dg_dx1
        ln_gamma[1] = gE_RT - xi[0] * dg_dx1
        return ln_gamma

    def excess_gibbs_free_energy(
        self,
        mole_fraction: Optional[Dict[str, float]] = None,
        A: Optional[float] = None,
        A12: Optional[float] = None,
        A21: Optional[float] = None,
        calculation_mode: Literal["two_suffix", "three_suffix"] = "three_suffix",
        message: Optional[str] = None,
        res_format: Literal["str", "json", "dict"] = "dict",
    ) -> Dict[str, Any] | str:
        try:
            # SECTION: validate inputs
            self._require_binary()
            xi = self._latest_mole_fraction_array(mole_fraction)
            x1 = float(xi[0])

            if calculation_mode == "two_suffix":
                if A is None:
                    raise ValueError("A is required")
                gE_RT, _ = self._two_suffix(x1=x1, A=float(A))
            elif calculation_mode == "three_suffix":
                if A12 is None or A21 is None:
                    raise ValueError("A12 and A21 are required")
                gE_RT, _ = self._three_suffix(
                    x1=x1,
                    A12=float(A12),
                    A21=float(A21),
                )
            else:
                raise ValueError("calculation_mode must be two_suffix or three_suffix")

            return self._excess_result(
                value=gE_RT,
                mole_fraction=xi,
                message=message,
                res_format=res_format,
            )
        except Exception as e:
            raise Exception(f"Error in Margules excess_gibbs_free_energy: {str(e)}")
