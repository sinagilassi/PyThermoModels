from math import log, sqrt
from typing import Dict, Literal, Optional

import numpy as np


LongRangeModel = Literal["pitzer_debye_huckel", "debye_huckel"]


class ENRTLLongRange:
    """Long-range electrostatic contribution for ENRTL."""

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
        if ionic_strength < 0 or not np.isfinite(ionic_strength):
            raise ValueError("ionic_strength must be finite and non-negative")

        if model == "pitzer_debye_huckel":
            return self.cal_ln_gamma_pitzer_debye_huckel(
                ionic_strength=ionic_strength,
                charges=charges,
                A_phi=A_phi,
                b=b,
            )

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
        if A is None or B is None:
            raise ValueError("Debye-Huckel long-range model requires A and B")

        ion_size = {} if ion_size is None else ion_size
        sqrt_i = sqrt(ionic_strength)
        result = {}

        for key, charge in charges.items():
            z = int(charge)
            if z == 0 or ionic_strength == 0:
                result[key] = 0.0
                continue

            if key not in ion_size:
                raise ValueError(f"Missing ion_size for charged component '{key}'")

            a_i = float(ion_size[key])
            if a_i < 0 or not np.isfinite(a_i):
                raise ValueError(f"ion_size for '{key}' must be finite and non-negative")

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
        if A_phi is None:
            raise ValueError("Pitzer-Debye-Huckel long-range model requires A_phi")

        b = 1.2 if b is None else float(b)
        if b <= 0 or not np.isfinite(b):
            raise ValueError("Pitzer-Debye-Huckel parameter b must be positive")

        if ionic_strength == 0:
            return {key: 0.0 for key in charges}

        sqrt_i = sqrt(ionic_strength)
        term = sqrt_i / (1.0 + b * sqrt_i) + (2.0 / b) * log(1.0 + b * sqrt_i)
        f_gamma = -float(A_phi) * term

        return {
            key: float(int(charge) * int(charge) * f_gamma)
            for key, charge in charges.items()
        }
