from typing import Dict, Literal, Optional, Tuple

import numpy as np


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
        if formulation != "chen_evans_1986":
            raise ValueError(
                "ENRTL v1 supports only formulation 'chen_evans_1986'."
            )

        self.components = components
        self.comp_idx = comp_idx
        self.comp_num = len(components)
        self.formulation = formulation

    def cal_ln_gamma_lc(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
        mode: Literal["chen_evans_1986", "neutral_nrtl_limit"] = "chen_evans_1986",
    ) -> np.ndarray:
        self._validate_vectors(mole_fraction, charges)
        self._validate_matrix(tau_ij, "tau_ij")
        self._validate_matrix(alpha_ij, "alpha_ij")

        if mode != "chen_evans_1986":
            if mode == "neutral_nrtl_limit":
                return self.cal_ln_gamma_neutral_nrtl_limit(
                    mole_fraction=mole_fraction,
                    charges=charges,
                    tau_ij=tau_ij,
                    alpha_ij=alpha_ij,
                )
            raise ValueError(f"Unsupported ENRTL local-composition mode: {mode}")

        if np.any(charges != 0):
            raise NotImplementedError(
                "Chen-Evans 1986 electrolyte local-composition equations require "
                "the full electrolyte-specific parameter structure and equation "
                "implementation. This path intentionally fails instead of using "
                "ordinary NRTL for ionic species."
            )

        return self.cal_ln_gamma_neutral_nrtl_limit(
            mole_fraction=mole_fraction,
            charges=charges,
            tau_ij=tau_ij,
            alpha_ij=alpha_ij,
        )

    def cal_ln_gamma_neutral_nrtl_limit(
        self,
        mole_fraction: np.ndarray,
        charges: np.ndarray,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
    ) -> np.ndarray:
        if np.any(charges != 0):
            raise ValueError(
                "neutral_nrtl_limit is valid only when every component has zero charge"
            )

        G_ij = self.cal_G_ij(tau_ij=tau_ij, alpha_ij=alpha_ij)
        ln_gamma = np.zeros(self.comp_num, dtype=float)

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

    def cal_G_ij(
        self,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
    ) -> np.ndarray:
        G_ij = np.exp(-alpha_ij * tau_ij)
        np.fill_diagonal(G_ij, 1.0)
        return G_ij

    def build_local_composition_diagnostics(
        self,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
        symbol_delimiter: Literal["|", "_"] = "|",
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        G_ij = self.cal_G_ij(tau_ij=tau_ij, alpha_ij=alpha_ij)
        delimiter = " | " if symbol_delimiter == "|" else "_"
        return G_ij, {
            f"{self.components[i]}{delimiter}{self.components[j]}": float(G_ij[i, j])
            for i in range(self.comp_num)
            for j in range(self.comp_num)
        }

    def _validate_vectors(self, mole_fraction: np.ndarray, charges: np.ndarray) -> None:
        if mole_fraction.shape != (self.comp_num,):
            raise ValueError(f"mole_fraction must have shape ({self.comp_num},)")
        if charges.shape != (self.comp_num,):
            raise ValueError(f"charges must have shape ({self.comp_num},)")
        if not np.all(np.isfinite(mole_fraction)):
            raise ValueError("mole_fraction contains non-finite values")

    def _validate_matrix(self, matrix: np.ndarray, name: str) -> None:
        if matrix.shape != (self.comp_num, self.comp_num):
            raise ValueError(f"{name} must have shape ({self.comp_num}, {self.comp_num})")
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} contains non-finite values")
