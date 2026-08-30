# import libs
from typing import Dict, Literal, Optional, Tuple

import numpy as np

# NOTE: supported ENRTL local-composition literature formulations
ENRTLFormulation = Literal["chen_evans_1986"]


class ENRTLLocalComposition:
    """Chen-Evans ENRTL local-composition contribution boundary.

    The production ENRTL formulation is Chen and Evans (1986). This class
    intentionally does not treat ordinary NRTL as a substitute for ionic
    Chen-Evans local-composition equations.
    """

    def __init__(
        self,
        components: list[str],
        comp_idx: dict[str, int],
        formulation: ENRTLFormulation = "chen_evans_1986",
    ) -> None:
        """
        Initialize the local-composition contribution for a fixed component set.

        Parameters
        ----------
        components : list[str]
            Ordered list of true-species component keys.
        comp_idx : dict[str, int]
            Mapping of component key to its index/position.
        formulation : ENRTLFormulation
            Literature formulation identifier; only `"chen_evans_1986"` is
            supported in version 1.

        Raises
        ------
        ValueError
            If `formulation` is not `"chen_evans_1986"`.
        """
        # NOTE: guard against silently mixing equations from other ENRTL variants
        if formulation != "chen_evans_1986":
            raise ValueError(
                "ENRTL v1 supports only formulation 'chen_evans_1986'."
            )

        self.components = components
        self.comp_idx = comp_idx
        self.comp_num = len(components)
        self.formulation = formulation

    # ! Calculate local-composition contribution (ln(gamma_i^LC))
    def cal_ln_gamma_lc(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
        mode: Literal[
            "chen_evans_1986",
            "neutral_nrtl_limit"
        ] = "chen_evans_1986",
    ) -> np.ndarray:
        """
        Calculate `ln(gamma_i^LC)`, the local-composition contribution.

        Parameters
        ----------
        mole_fraction : np.ndarray
            True-species mole fractions, shape `(comp_num,)`.
        charges : np.ndarray
            Net ionic charge per component, shape `(comp_num,)`.
        tau_ij : np.ndarray
            Binary interaction parameter matrix, shape `(comp_num, comp_num)`.
        alpha_ij : np.ndarray
            Non-randomness parameter matrix, shape `(comp_num, comp_num)`.
        mode : Literal["chen_evans_1986", "neutral_nrtl_limit"]
            Calculation path. `"chen_evans_1986"` requires the full
            electrolyte formulation for any charged species;
            `"neutral_nrtl_limit"` is valid only for all-neutral mixtures.

        Returns
        -------
        np.ndarray
            `ln(gamma_i^LC)`, shape `(comp_num,)`.

        Raises
        ------
        ValueError
            If input shapes/values are invalid or `mode` is unsupported.
        NotImplementedError
            If `mode="chen_evans_1986"` and any species is charged, since the
            ionic Chen-Evans equations are not yet implemented here.
        """
        # SECTION: validate inputs
        self._validate_vectors(mole_fraction, charges)
        self._validate_matrix(tau_ij, "tau_ij")
        self._validate_matrix(alpha_ij, "alpha_ij")

        # SECTION: explicit neutral-limit path (bypasses the Chen-Evans ionic check)
        if mode != "chen_evans_1986":
            if mode == "neutral_nrtl_limit":
                return self.cal_ln_gamma_neutral_nrtl_limit(
                    mole_fraction=mole_fraction,
                    charges=charges,
                    tau_ij=tau_ij,
                    alpha_ij=alpha_ij,
                )
            raise ValueError(
                f"Unsupported ENRTL local-composition mode: {mode}")

        # NOTE: fail loudly instead of silently reusing ordinary NRTL for ions
        if np.any(charges != 0):
            raise NotImplementedError(
                "Chen-Evans 1986 electrolyte local-composition equations require "
                "the full electrolyte-specific parameter structure and equation "
                "implementation. This path intentionally fails instead of using "
                "ordinary NRTL for ionic species."
            )

        # SECTION: all-neutral mixture -> reduces to the ordinary NRTL limit
        return self.cal_ln_gamma_neutral_nrtl_limit(
            mole_fraction=mole_fraction,
            charges=charges,
            tau_ij=tau_ij,
            alpha_ij=alpha_ij,
        )

    # ! Ordinary NRTL local-composition limit for neutral mixtures
    def cal_ln_gamma_neutral_nrtl_limit(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
    ) -> np.ndarray:
        """
        Calculate `ln(gamma_i)` using the ordinary NRTL local-composition equation.

        This is the neutral-mixture limiting case of the ENRTL local-composition
        contribution and is valid only when every component has zero charge.

        ```
        ln(gamma_i) = sum_j(x_j * tau_ji * G_ji) / sum_k(x_k * G_ki)
            + sum_j( x_j * G_ij / sum_k(x_k * G_kj)
                * (tau_ij - sum_k(x_k * tau_kj * G_kj) / sum_k(x_k * G_kj)) )
        ```

        Parameters
        ----------
        mole_fraction : np.ndarray
            True-species mole fractions, shape `(comp_num,)`.
        charges : np.ndarray
            Net ionic charge per component, shape `(comp_num,)`; must all be zero.
        tau_ij : np.ndarray
            Binary interaction parameter matrix, shape `(comp_num, comp_num)`.
        alpha_ij : np.ndarray
            Non-randomness parameter matrix, shape `(comp_num, comp_num)`.

        Returns
        -------
        np.ndarray
            `ln(gamma_i)`, shape `(comp_num,)`.

        Raises
        ------
        ValueError
            If any component has a nonzero charge.
        """
        # NOTE: this equation is only thermodynamically valid for neutral mixtures
        if np.any(charges != 0):
            raise ValueError(
                "neutral_nrtl_limit is valid only when every component has zero charge"
            )

        G_ij = self.cal_G_ij(tau_ij=tau_ij, alpha_ij=alpha_ij)
        ln_gamma = np.zeros(self.comp_num, dtype=float)

        # SECTION: ordinary NRTL summation per component i
        for i in range(self.comp_num):
            denom_i = np.sum(G_ij[:, i] * mole_fraction)
            numer_i = np.sum(tau_ij[:, i] * G_ij[:, i] * mole_fraction)
            first_term = numer_i / denom_i

            second_term = 0.0
            for j in range(self.comp_num):
                denom_j = np.sum(G_ij[:, j] * mole_fraction)
                numer_j = np.sum(mole_fraction * tau_ij[:, j] * G_ij[:, j])
                second_term += (
                    mole_fraction[j] * G_ij[i, j] / denom_j
                ) * (tau_ij[i, j] - numer_j / denom_j)

            ln_gamma[i] = first_term + second_term

        return ln_gamma

    # ! Calculate the NRTL local-composition weighting matrix G_ij
    def cal_G_ij(
        self,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
    ) -> np.ndarray:
        """
        Calculate the NRTL local-composition weighting matrix `G_ij = exp(-alpha_ij * tau_ij)`.

        Parameters
        ----------
        tau_ij : np.ndarray
            Binary interaction parameter matrix, shape `(comp_num, comp_num)`.
        alpha_ij : np.ndarray
            Non-randomness parameter matrix, shape `(comp_num, comp_num)`.

        Returns
        -------
        np.ndarray
            `G_ij` matrix with diagonal forced to `1.0`.
        """
        G_ij = np.exp(-alpha_ij * tau_ij)
        # NOTE: G_ii is not calculated from tau/alpha; it is defined as 1 by convention
        np.fill_diagonal(G_ij, 1.0)
        return G_ij

    # ! Build local-composition diagnostics (G_ij matrix and component-keyed mapping)
    def build_local_composition_diagnostics(
        self,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
        symbol_delimiter: Literal["|", "_"] = "|",
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Build the `G_ij` matrix along with a component-keyed diagnostic mapping.

        Parameters
        ----------
        tau_ij : np.ndarray
            Binary interaction parameter matrix, shape `(comp_num, comp_num)`.
        alpha_ij : np.ndarray
            Non-randomness parameter matrix, shape `(comp_num, comp_num)`.
        symbol_delimiter : Literal["|", "_"]
            Delimiter used to join component-pair keys in the diagnostic dict.

        Returns
        -------
        Tuple[np.ndarray, Dict[str, float]]
            The `G_ij` matrix and a `"comp_i<delimiter>comp_j" -> G_ij` mapping.
        """
        G_ij = self.cal_G_ij(tau_ij=tau_ij, alpha_ij=alpha_ij)
        delimiter = " | " if symbol_delimiter == "|" else "_"
        return G_ij, {
            f"{self.components[i]}{delimiter}{self.components[j]}": float(G_ij[i, j])
            for i in range(self.comp_num)
            for j in range(self.comp_num)
        }

    # ! Validate vectors and matrices for local-composition calculations
    def _validate_vectors(self, mole_fraction: np.ndarray, charges: np.ndarray) -> None:
        """Validate that `mole_fraction`/`charges` have the expected shape and are finite."""
        if mole_fraction.shape != (self.comp_num,):
            raise ValueError(
                f"mole_fraction must have shape ({self.comp_num},)")
        if charges.shape != (self.comp_num,):
            raise ValueError(f"charges must have shape ({self.comp_num},)")
        if not np.all(np.isfinite(mole_fraction)):
            raise ValueError("mole_fraction contains non-finite values")

    # ! Validate a parameter matrix for local-composition calculations
    def _validate_matrix(self, matrix: np.ndarray, name: str) -> None:
        """Validate that a `(comp_num, comp_num)` parameter matrix is well-formed and finite."""
        if matrix.shape != (self.comp_num, self.comp_num):
            raise ValueError(
                f"{name} must have shape ({self.comp_num}, {self.comp_num})")
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} contains non-finite values")
