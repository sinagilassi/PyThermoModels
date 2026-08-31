# import libs
from typing import Literal

import numpy as np

# NOTE: composition/reference-state scale used for activity and ionic-strength calculations
ActivityBasis = Literal["mole_fraction", "molality", "molarity"]
# NOTE: distinguishes user-facing feed composition from the true liquid-phase species state
CompositionRepresentation = Literal["true_species", "apparent_species"]


class ENRTLCore:
    """Core electrolyte bookkeeping shared across ENRTL contributions.

    Provides electroneutrality checks, ionic strength, and contribution
    combination utilities that are independent of the specific long-range
    or local-composition formulation used.
    """

    # tolerance for composition sanity checks
    COMPOSITION_ATOL = 1e-12
    # tolerance for the electroneutrality (net charge) check
    CHARGE_BALANCE_ATOL = 1e-12
    # universal gas constant [J/mol/K]
    R_CONST = 8.31446261815324

    def net_charge(self, composition: np.ndarray, charges: np.ndarray) -> float:
        """
        Calculate the net charge `sum_i(z_i * x_i)` of a composition.

        Parameters
        ----------
        composition : np.ndarray
            True-species composition (basis defined by the caller).
        charges : np.ndarray
            Net ionic charge per component, same shape as `composition`.

        Returns
        -------
        float
            Net charge; should be `~0` for an electrically neutral mixture.
        """
        self._validate_vector_pair(composition, charges)
        return float(np.dot(composition, charges))

    def check_electroneutrality(
        self,
        composition: np.ndarray,
        charges: np.ndarray,
        atol: float = CHARGE_BALANCE_ATOL,
    ) -> bool:
        """
        Validate that a true-species composition satisfies electroneutrality.

        Parameters
        ----------
        composition : np.ndarray
            True-species composition (basis defined by the caller).
        charges : np.ndarray
            Net ionic charge per component, same shape as `composition`.
        atol : float
            Absolute tolerance on the net charge.

        Returns
        -------
        bool
            `True` if the composition is electrically neutral within `atol`.

        Raises
        ------
        ValueError
            If the net charge exceeds `atol`.
        """
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
        """
        Calculate the ionic strength `I = 0.5 * sum_i(x_i * z_i^2)`.

        Parameters
        ----------
        composition : np.ndarray
            True-species composition expressed on `basis`.
        charges : np.ndarray
            Net ionic charge per component, same shape as `composition`.
        basis : ActivityBasis
            Composition scale (`"mole_fraction"`, `"molality"`, or
            `"molarity"`); recorded for traceability, not converted here.

        Returns
        -------
        float
            Ionic strength on the supplied `basis`.

        Raises
        ------
        ValueError
            If `basis` is not one of the supported values.
        """
        # NOTE: basis is not converted here; caller must supply composition already on this scale
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
        """
        Sum the `ln(gamma_i)` contributions into the total activity coefficient.

        ```
        ln(gamma_i) = ln(gamma_i^LR) + ln(gamma_i^LC) [+ ln(gamma_i^Born)]
        ```

        Parameters
        ----------
        ln_gamma_long_range : np.ndarray
            Long-range electrostatic contribution.
        ln_gamma_local_composition : np.ndarray
            Local-composition contribution, same shape as `ln_gamma_long_range`.
        ln_gamma_born : Optional[np.ndarray]
            Optional Born/solvation contribution, same shape as the others.

        Returns
        -------
        np.ndarray
            Total `ln(gamma_i)`.

        Raises
        ------
        ValueError
            If contribution shapes disagree or the result is non-finite.
        """
        if ln_gamma_long_range.shape != ln_gamma_local_composition.shape:
            raise ValueError(
                "ln_gamma contribution arrays must have the same shape")

        # NOTE: contributions are summed in log-space for numerical robustness
        total = ln_gamma_long_range + ln_gamma_local_composition
        if ln_gamma_born is not None:
            if ln_gamma_born.shape != total.shape:
                raise ValueError(
                    "ln_gamma_born must match other contribution shapes")
            total = total + ln_gamma_born

        if not np.all(np.isfinite(total)):
            raise ValueError("combined ln_gamma contains non-finite values")

        return total

    def _validate_vector_pair(self, composition: np.ndarray, charges: np.ndarray) -> None:
        """Validate that `composition`/`charges` are finite, non-negative, and shape-matched."""
        if composition.shape != charges.shape:
            raise ValueError(
                "composition and charges must have the same shape")
        if not np.all(np.isfinite(composition)):
            raise ValueError("composition contains non-finite values")
        if np.any(composition < 0):
            raise ValueError("composition values must be non-negative")
        if not np.all(np.isfinite(charges)):
            raise ValueError("charges contains non-finite values")
