# import libs
from typing import Dict, Literal, Optional, Tuple

import numpy as np
from .parameter_core import ENRTLInteractionParameters

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
        self._last_diagnostics: Dict[str, object] = {}

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

        parameters = ENRTLInteractionParameters.from_tau_alpha(
            component_keys=self.components,
            charges=charges,
            tau_ij=tau_ij,
            alpha_ij=alpha_ij,
        )

        # SECTION: ionic Chen-Evans local-composition path
        if np.any(charges != 0):
            return self._cal_ln_gamma_chen_evans_1986(
                mole_fraction=mole_fraction,
                charges=charges,
                parameters=parameters,
            )

        # SECTION: all-neutral mixture -> reduces to the ordinary NRTL limit
        return self.cal_ln_gamma_neutral_nrtl_limit(
            mole_fraction=mole_fraction,
            charges=charges,
            tau_ij=tau_ij,
            alpha_ij=alpha_ij,
        )

    def _cal_ln_gamma_chen_evans_1986(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        parameters: ENRTLInteractionParameters,
    ) -> np.ndarray:
        """Calculate ionic Chen-Evans short-range activity coefficients.

        Chen and Evans (1986), short-range/local-composition electrolyte NRTL
        excess-Gibbs formulation. The implementation evaluates the
        local-composition excess Gibbs energy kernel and obtains
        `ln(gamma_i^LC)` from the thermodynamic derivative
        `d(n g^E_LC/RT) / d n_i`.

        Index convention
        ----------------
        `tau[j, i]` is the interaction parameter for local/neighbor species
        `j` around central species `i`; `G[j, i] = exp(-alpha[j, i] tau[j, i])`.

        Composition convention
        ----------------------
        The Chen-Evans charge-weighted local composition is represented as
        `X_i = x_i` for molecular species and `X_i = |z_i| x_i` for ions.
        Like-ion local compositions are excluded through
        `parameters.interaction_mask`.
        """
        cation_idx = np.flatnonzero(charges > 0)
        anion_idx = np.flatnonzero(charges < 0)
        if cation_idx.size == 0 or anion_idx.size == 0:
            raise ValueError(
                "Chen-Evans 1986 ionic local-composition calculations require "
                "at least one cation and one anion."
            )

        g_value, kernel_diagnostics = self._chen_evans_excess_gibbs_kernel(
            mole_fraction=mole_fraction,
            charges=charges,
            parameters=parameters,
        )
        ln_gamma = self._partial_molar_ln_gamma_from_kernel(
            mole_fraction=mole_fraction,
            charges=charges,
            parameters=parameters,
        )
        if not np.all(np.isfinite(ln_gamma)):
            raise ValueError("Chen-Evans ln_gamma_local_composition contains non-finite values")

        self._last_diagnostics = {
            "formulation": self.formulation,
            "gE_local_composition_RT": g_value,
            "component_roles": {
                item.key: item.role for item in parameters.components
            },
            "interaction_type": parameters.interaction_type,
            "interaction_mask": parameters.interaction_mask,
            **kernel_diagnostics,
        }
        return ln_gamma

    def _partial_molar_ln_gamma_from_kernel(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        parameters: ENRTLInteractionParameters,
    ) -> np.ndarray:
        """Differentiate `n g^E_LC/RT` with respect to component mole numbers."""
        n0 = np.asarray(mole_fraction, dtype=float)
        ln_gamma = np.zeros(self.comp_num, dtype=float)

        def total_excess(moles: np.ndarray) -> float:
            total = float(np.sum(moles))
            if total <= 0.0:
                raise ValueError("total mole number must be positive")
            x = moles / total
            value, _ = self._chen_evans_excess_gibbs_kernel(
                mole_fraction=x,
                charges=charges,
                parameters=parameters,
            )
            return total * value

        for i in range(self.comp_num):
            step = 1e-6 * max(1.0, abs(float(n0[i])))
            plus = n0.copy()
            plus[i] += step
            if n0[i] > step:
                minus = n0.copy()
                minus[i] -= step
                ln_gamma[i] = (total_excess(plus) - total_excess(minus)) / (2.0 * step)
            else:
                ln_gamma[i] = (total_excess(plus) - total_excess(n0)) / step

        return ln_gamma

    def _chen_evans_excess_gibbs_kernel(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        parameters: ENRTLInteractionParameters,
    ) -> Tuple[float, Dict[str, object]]:
        """Evaluate the Chen-Evans 1986 local-composition `g^E_LC/RT` kernel.

        The kernel uses three local cell families: molecular centers, cationic
        centers with local anions/molecules, and anionic centers with local
        cations/molecules. This encodes the Chen-Evans local electroneutrality
        and like-ion repulsion postulates at the local-composition level.
        """
        x_eff = self._effective_mole_fraction(mole_fraction, charges)
        G_ij = self.cal_G_ij(tau_ij=parameters.tau, alpha_ij=parameters.alpha)

        molecule_idx = np.flatnonzero(charges == 0)
        cation_idx = np.flatnonzero(charges > 0)
        anion_idx = np.flatnonzero(charges < 0)
        cation_charge_sum = float(np.sum(x_eff[cation_idx]))
        anion_charge_sum = float(np.sum(x_eff[anion_idx]))
        if cation_charge_sum <= 0.0 or anion_charge_sum <= 0.0:
            raise ValueError("charged Chen-Evans systems require positive cation and anion charge fractions")

        value = 0.0
        local_fractions = np.zeros((self.comp_num, self.comp_num), dtype=float)

        for center in molecule_idx:
            neighbors = np.arange(self.comp_num)
            term, fractions = self._local_cell_term(center, neighbors, x_eff, G_ij, parameters)
            value += x_eff[center] * term
            local_fractions[:, center] = fractions

        for center in cation_idx:
            neighbors = np.concatenate((molecule_idx, anion_idx))
            term, fractions = self._local_cell_term(center, neighbors, x_eff, G_ij, parameters)
            value += x_eff[center] * term
            local_fractions[:, center] = fractions

        for center in anion_idx:
            neighbors = np.concatenate((molecule_idx, cation_idx))
            term, fractions = self._local_cell_term(center, neighbors, x_eff, G_ij, parameters)
            value += x_eff[center] * term
            local_fractions[:, center] = fractions

        local_charge_residual = local_fractions.T @ charges
        return float(value), {
            "effective_mole_fraction": x_eff,
            "local_composition_fraction": local_fractions,
            "local_electroneutrality_residual": local_charge_residual,
        }

    def _local_cell_term(
        self,
        center: int,
        neighbors: np.ndarray,
        x_eff: np.ndarray,
        G_ij: np.ndarray,
        parameters: ENRTLInteractionParameters,
    ) -> Tuple[float, np.ndarray]:
        """Return the local NRTL energy term around one Chen-Evans cell center."""
        weights = np.zeros(self.comp_num, dtype=float)
        for neighbor in neighbors:
            if parameters.interaction_mask[neighbor, center]:
                weights[neighbor] = x_eff[neighbor] * G_ij[neighbor, center]

        denominator = float(np.sum(weights))
        if denominator <= 0.0 or not np.isfinite(denominator):
            raise ValueError(
                "Chen-Evans local-composition denominator is zero for "
                f"central species '{self.components[center]}'."
            )

        fractions = weights / denominator
        term = float(np.sum(fractions * parameters.tau[:, center]))
        return term, fractions

    def _effective_mole_fraction(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
    ) -> np.ndarray:
        """Build Chen-Evans charge-weighted local-composition fractions."""
        charge_weight = np.where(charges == 0, 1.0, np.abs(charges).astype(float))
        x_eff = mole_fraction * charge_weight
        if not np.all(np.isfinite(x_eff)) or np.any(x_eff < 0.0):
            raise ValueError("effective mole fractions must be finite and non-negative")
        return x_eff

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

    @property
    def last_diagnostics(self) -> Dict[str, object]:
        """Return diagnostics from the most recent ionic local-composition call."""
        return self._last_diagnostics

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
