from typing import Literal

import numpy as np


ActivityBasis = Literal["mole_fraction", "molality", "molarity"]
CompositionRepresentation = Literal["true_species", "apparent_species"]


class ENRTLCore:
    COMPOSITION_ATOL = 1e-12
    CHARGE_BALANCE_ATOL = 1e-12
    R_CONST = 8.31446261815324

    def net_charge(self, composition: np.ndarray, charges: np.ndarray) -> float:
        self._validate_vector_pair(composition, charges)
        return float(np.dot(composition, charges))

    def check_electroneutrality(
        self,
        composition: np.ndarray,
        charges: np.ndarray,
        atol: float = CHARGE_BALANCE_ATOL,
    ) -> bool:
        value = self.net_charge(composition=composition, charges=charges)
        if abs(value) > atol:
            raise ValueError(
                "ENRTL composition is not electrically neutral: "
                f"net charge = {value:.6e}"
            )
        return True

    def ionic_strength(
        self,
        composition: np.ndarray,
        charges: np.ndarray,
        basis: ActivityBasis,
    ) -> float:
        if basis not in ("mole_fraction", "molality", "molarity"):
            raise ValueError(
                "ionic-strength basis must be 'mole_fraction', 'molality', or 'molarity'"
            )

        self._validate_vector_pair(composition, charges)
        value = 0.5 * np.sum(composition * np.square(charges))
        return float(value)

    def combine_ln_gamma(
        self,
        ln_gamma_long_range: np.ndarray,
        ln_gamma_local_composition: np.ndarray,
        ln_gamma_born: np.ndarray | None = None,
    ) -> np.ndarray:
        if ln_gamma_long_range.shape != ln_gamma_local_composition.shape:
            raise ValueError("ln_gamma contribution arrays must have the same shape")

        total = ln_gamma_long_range + ln_gamma_local_composition
        if ln_gamma_born is not None:
            if ln_gamma_born.shape != total.shape:
                raise ValueError("ln_gamma_born must match other contribution shapes")
            total = total + ln_gamma_born

        if not np.all(np.isfinite(total)):
            raise ValueError("combined ln_gamma contains non-finite values")

        return total

    def _validate_vector_pair(self, composition: np.ndarray, charges: np.ndarray) -> None:
        if composition.shape != charges.shape:
            raise ValueError("composition and charges must have the same shape")
        if not np.all(np.isfinite(composition)):
            raise ValueError("composition contains non-finite values")
        if np.any(composition < 0):
            raise ValueError("composition values must be non-negative")
        if not np.all(np.isfinite(charges)):
            raise ValueError("charges contains non-finite values")
