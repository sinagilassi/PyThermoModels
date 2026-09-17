# import libs
from typing import Any, Dict, Literal, Optional, Tuple

import numpy as np

# local
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from ..classical_base import ClassicalActivityBase


class VanLaar(ClassicalActivityBase):
    """
    Binary van Laar activity-coefficient model.
    """

    @add_attributes(metadata=ACTIVITY_MODELS["VAN_LAAR"])
    def cal(
        self,
        model_input: Dict,
        message: Optional[str] = None,
        **kwargs
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        try:
            # SECTION: validate and unpack inputs
            self._require_binary()
            self._validate_model_input(model_input)
            xi = self._mole_fraction_array(model_input["mole_fraction"])

            # ? does caller provide both ordered van Laar parameters?
            if "a1" not in model_input or "a2" not in model_input:
                raise KeyError("a1 and a2 are required in model_input")
            a1 = float(model_input["a1"])
            a2 = float(model_input["a2"])

            # SECTION: calculate model values
            gE_RT, dg_dx1 = self._gibbs_and_derivative(
                x1=float(xi[0]),
                a1=a1,
                a2=a2,
            )
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
                model_name="van Laar",
                message=message,
            )
            other_values = {
                "AcCo_i_comp": AcCo_i_comp,
                "ln_gamma": ln_gamma,
                "a1": a1,
                "a2": a2,
                "excess_gibbs_RT": gE_RT,
            }
            return res, other_values
        except Exception as e:
            raise Exception(f"Error in van Laar model cal: {str(e)}")

    def _gibbs_and_derivative(
        self,
        x1: float,
        a1: float,
        a2: float,
    ) -> Tuple[float, float]:
        x2 = 1.0 - x1
        denominator = a1 * x1 + a2 * x2
        # ! denominator must remain finite for the reciprocal model
        if np.isclose(denominator, 0.0):
            raise ValueError("van Laar denominator cannot be zero")

        numerator = a1 * a2 * x1 * x2
        gE_RT = numerator / denominator
        dg_dx1 = (
            a1 * a2 * (
                (1.0 - 2.0 * x1) * denominator
                - x1 * x2 * (a1 - a2)
            )
            / denominator**2
        )
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
        a1: Optional[float] = None,
        a2: Optional[float] = None,
        message: Optional[str] = None,
        res_format: Literal["str", "json", "dict"] = "dict",
    ) -> Dict[str, Any] | str:
        try:
            # SECTION: validate inputs
            self._require_binary()
            xi = self._latest_mole_fraction_array(mole_fraction)
            if a1 is None or a2 is None:
                raise ValueError("a1 and a2 are required")

            gE_RT, _ = self._gibbs_and_derivative(
                x1=float(xi[0]),
                a1=float(a1),
                a2=float(a2),
            )
            return self._excess_result(
                value=gE_RT,
                mole_fraction=xi,
                message=message,
                res_format=res_format,
            )
        except Exception as e:
            raise Exception(f"Error in van Laar excess_gibbs_free_energy: {str(e)}")
