# import libs
from math import log, sqrt
from typing import Dict, Literal, Optional

import numpy as np

# NOTE: supported ENRTL long-range electrostatic formulations
LongRangeModel = Literal["pitzer_debye_huckel", "debye_huckel"]


class ENRTLLongRange:
    """
    Long-range electrostatic contribution for the Electrolyte NRTL (ENRTL) model.

    This class calculates `ln(gamma_i)` for the long-range electrostatic
    contribution using either the extended Debye-Huckel expression or the
    Pitzer-Debye-Huckel expression. Neutral species (charge = 0) always
    contribute zero to this term.
    """

    def cal_ln_gamma_long_range(
        self,
        ionic_strength: float,
        charges: Dict[str, int],
        model: LongRangeModel = "pitzer_debye_huckel",
        ion_size: Optional[Dict[str, float]] = None,
        A: Optional[float] = None,
        B: Optional[float] = None,
        A_phi: Optional[float] = None,
        b: Optional[float] = None,
    ) -> Dict[str, float]:
        """
        Dispatch and calculate `ln(gamma_i^LR)` for the selected long-range model.

        Parameters
        ----------
        ionic_strength : float
            Ionic strength of the true-species composition (basis defined by caller).
        charges : Dict[str, int]
            Mapping of component key to net ionic charge `z_i`.
        model : LongRangeModel
            Long-range electrostatic formulation to use.
        ion_size : Optional[Dict[str, float]]
            Ion-size parameter `a_i` per component, required for `"debye_huckel"`.
        A : Optional[float]
            Debye-Huckel constant `A`, required for `"debye_huckel"`.
        B : Optional[float]
            Debye-Huckel constant `B`, required for `"debye_huckel"`.
        A_phi : Optional[float]
            Pitzer-Debye-Huckel constant, required for `"pitzer_debye_huckel"`.
        b : Optional[float]
            Pitzer-Debye-Huckel closest-approach parameter (default 1.2 if not supplied).

        Returns
        -------
        Dict[str, float]
            Mapping of component key to `ln(gamma_i^LR)`.

        Raises
        ------
        ValueError
            If `ionic_strength` is negative/non-finite or `model` is not supported.
        """
        # SECTION: validate ionic strength
        if ionic_strength < 0 or not np.isfinite(ionic_strength):
            raise ValueError("ionic_strength must be finite and non-negative")

        # SECTION: dispatch to the selected long-range formulation
        # 1. Pitzer-Debye-Huckel
        if model == "pitzer_debye_huckel":
            return self.cal_ln_gamma_pitzer_debye_huckel(
                ionic_strength=ionic_strength,
                charges=charges,
                A_phi=A_phi,
                b=b,
            )

        # 2. extended Debye-Huckel
        if model == "debye_huckel":
            return self.cal_ln_gamma_debye_huckel(
                ionic_strength=ionic_strength,
                charges=charges,
                ion_size=ion_size,
                A=A,
                B=B,
            )

        raise ValueError(f"Unsupported ENRTL long-range model: {model}")

    def cal_ln_gamma_debye_huckel(
        self,
        ionic_strength: float,
        charges: Dict[str, int],
        ion_size: Optional[Dict[str, float]],
        A: Optional[float],
        B: Optional[float],
    ) -> Dict[str, float]:
        """
        Calculate `ln(gamma_i^LR)` using the extended Debye-Huckel expression.

        ```
        ln(gamma_i^LR) = -ln(10) * A * z_i^2 * sqrt(I) / (1 + B * a_i * sqrt(I))
        ```

        Parameters
        ----------
        ionic_strength : float
            Ionic strength `I`.
        charges : Dict[str, int]
            Mapping of component key to net ionic charge `z_i`.
        ion_size : Optional[Dict[str, float]]
            Ion-size parameter `a_i` per charged component.
        A : Optional[float]
            Debye-Huckel constant `A`.
        B : Optional[float]
            Debye-Huckel constant `B`.

        Returns
        -------
        Dict[str, float]
            Mapping of component key to `ln(gamma_i^LR)`. Neutral species and
            zero ionic strength both yield `0.0`.

        Raises
        ------
        ValueError
            If `A`/`B` are missing, `ion_size` is missing for a charged
            component, or `ion_size` is negative/non-finite.
        """
        # NOTE: A and B are required inputs for this formulation (not derived internally)
        if A is None or B is None:
            raise ValueError("Debye-Huckel long-range model requires A and B")

        ion_size = {} if ion_size is None else ion_size
        sqrt_i = sqrt(ionic_strength)
        result = {}

        # SECTION: calculate ln(gamma_i^LR) per component
        for key, charge in charges.items():
            z = int(charge)
            # 1. neutral species or infinite dilution -> no contribution
            if z == 0 or ionic_strength == 0:
                result[key] = 0.0
                continue

            # 2. charged species require an ion-size parameter
            if key not in ion_size:
                raise ValueError(
                    f"Missing ion_size for charged component '{key}'")

            a_i = float(ion_size[key])
            if a_i < 0 or not np.isfinite(a_i):
                raise ValueError(
                    f"ion_size for '{key}' must be finite and non-negative")

            # 3. extended Debye-Huckel expression
            denominator = 1.0 + float(B) * a_i * sqrt_i
            result[key] = -log(10.0) * float(A) * z * z * sqrt_i / denominator

        return result

    def cal_ln_gamma_pitzer_debye_huckel(
        self,
        ionic_strength: float,
        charges: Dict[str, int],
        A_phi: Optional[float],
        b: Optional[float],
    ) -> Dict[str, float]:
        """
        Calculate `ln(gamma_i^LR)` using the Pitzer-Debye-Huckel expression.

        Parameters
        ----------
        ionic_strength : float
            Ionic strength `I`.
        charges : Dict[str, int]
            Mapping of component key to net ionic charge `z_i`.
        A_phi : Optional[float]
            Pitzer-Debye-Huckel constant.
        b : Optional[float]
            Closest-approach parameter; defaults to `1.2` if not supplied.

        Returns
        -------
        Dict[str, float]
            Mapping of component key to `ln(gamma_i^LR)`. All species yield
            `0.0` at zero ionic strength.

        Raises
        ------
        ValueError
            If `A_phi` is missing or `b` is not positive/finite.
        """
        # NOTE: A_phi is a required input for this formulation (not derived internally)
        if A_phi is None:
            raise ValueError(
                "Pitzer-Debye-Huckel long-range model requires A_phi")

        # NOTE: b defaults to the common literature value of 1.2 kg^0.5/mol^0.5
        b = 1.2 if b is None else float(b)
        if b <= 0 or not np.isfinite(b):
            raise ValueError(
                "Pitzer-Debye-Huckel parameter b must be positive")

        # SECTION: infinite-dilution limit -> zero contribution for all species
        if ionic_strength == 0:
            return {key: 0.0 for key in charges}

        # SECTION: Pitzer-Debye-Huckel expression, shared across all species via z_i^2
        sqrt_i = sqrt(ionic_strength)
        term = sqrt_i / (1.0 + b * sqrt_i) + (2.0 / b) * log(1.0 + b * sqrt_i)
        f_gamma = -float(A_phi) * term

        return {
            key: float(int(charge) * int(charge) * f_gamma)
            for key, charge in charges.items()
        }
