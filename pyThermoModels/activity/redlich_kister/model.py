# import libs
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np

# local
from ...plugin import ACTIVITY_MODELS
from ...utils import add_attributes
from ..classical_base import ClassicalActivityBase


class RedlichKister(ClassicalActivityBase):
    """
    Binary Redlich-Kister excess-Gibbs expansion.
    """

    @add_attributes(metadata=ACTIVITY_MODELS["REDLICH_KISTER"])
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

            # ? does caller provide Redlich-Kister expansion coefficients?
            if "a" not in model_input:
                raise KeyError("a is required in model_input")
            a = self._coefficient_array(model_input["a"])

            # SECTION: calculate model values
            gE_RT, dg_dx1 = self._gibbs_and_derivative(
                x1=float(xi[0]),
                a=a,
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
                model_name="Redlich-Kister",
                message=message,
            )
            other_values = {
                "AcCo_i_comp": AcCo_i_comp,
                "ln_gamma": ln_gamma,
                "a": a,
                "excess_gibbs_RT": gE_RT,
            }
            return res, other_values
        except Exception as e:
            raise Exception(f"Error in Redlich-Kister model cal: {str(e)}")

    def _coefficient_array(self, a: List[float] | np.ndarray) -> np.ndarray:
        # ? are coefficients provided as an ordered sequence a[0..n]?
        if not isinstance(a, (list, np.ndarray)):
            raise TypeError("a must be a list or numpy array")
        coeffs = np.asarray(a, dtype=float)
        if coeffs.ndim != 1 or coeffs.size == 0:
            raise ValueError("a must be a non-empty one-dimensional sequence")
        return coeffs

    def _gibbs_and_derivative(
        self,
        x1: float,
        a: np.ndarray,
    ) -> Tuple[float, float]:
        # NOTE: d = x1 - x2 = 2*x1 - 1 follows the NIST/TDE convention
        x2 = 1.0 - x1
        d = x1 - x2
        powers = np.asarray([d**k for k in range(a.size)], dtype=float)
        series = float(np.sum(a * powers))

        # SECTION: derivative of x1*x2*sum_k a_k*(x1-x2)^k
        if a.size == 1:
            dseries_dx1 = 0.0
        else:
            dseries_dx1 = float(np.sum([
                a[k] * k * d**(k - 1) * 2.0
                for k in range(1, a.size)
            ]))

        gE_RT = x1 * x2 * series
        dg_dx1 = (1.0 - 2.0 * x1) * series + x1 * x2 * dseries_dx1
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
        a: Optional[List[float] | np.ndarray] = None,
        message: Optional[str] = None,
        res_format: Literal["str", "json", "dict"] = "dict",
    ) -> Dict[str, Any] | str:
        try:
            # SECTION: validate inputs
            self._require_binary()
            xi = self._latest_mole_fraction_array(mole_fraction)
            if a is None:
                raise ValueError("a is required")

            gE_RT, _ = self._gibbs_and_derivative(
                x1=float(xi[0]),
                a=self._coefficient_array(a),
            )
            return self._excess_result(
                value=gE_RT,
                mole_fraction=xi,
                message=message,
                res_format=res_format,
            )
        except Exception as e:
            raise Exception(
                f"Error in Redlich-Kister excess_gibbs_free_energy: {str(e)}"
            )
