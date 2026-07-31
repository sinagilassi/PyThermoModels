# SECTION: imports
from math import exp, log, pow
from typing import Dict, Literal, Tuple

import numpy as np

from .local_composition_base import IJData, LocalCompositionBase


class UNIQUACLocalComposition(LocalCompositionBase):
    """UNIQUAC local-composition interaction-parameter calculations."""

    # SECTION: interaction energy
    def cal_dU_ij_M1(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC interaction-energy parameter.

        dU_ij = a_ij + b_ij*T + c_ij*T^2
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij)

            # SECTION: initialize result containers
            dU_ij = np.zeros((self.comp_num, self.comp_num))
            dU_ij_comp = {}

            # SECTION: calculate dU_ij
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        dU_ij[row, col] = 0
                        dU_ij_comp[key] = 0
                        continue

                    # NOTE: UNIQUAC M1 energy correlation
                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) * temperature
                        + self._ij_value(c_ij, "c", i, j, key)
                        * pow(temperature, 2)
                    )
                    dU_ij[row, col] = value
                    dU_ij_comp[key] = value

            return dU_ij, dU_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_dU_ij_M1: {str(e)}")

    # SECTION: tau from energy
    def cal_tau_ij_M1(
        self,
        temperature: float,
        dU_ij: IJData,
        dU_ij_symbol: Literal["dU", "dU_ij"] = "dU",
        R_CONST: float = 8.314,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC tau from interaction energy.

        tau_ij = exp(-dU_ij / (R*T))
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_ij_data(dU_ij, "dU_ij")

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        tau_ij[row, col] = 0
                        tau_ij_comp[key] = 0
                        continue

                    # NOTE: UNIQUAC tau is an exponential weighting factor
                    value = self._ij_value(dU_ij, dU_ij_symbol, i, j, key)
                    value = exp(-value / (R_CONST * temperature))
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M1: {str(e)}")

    # SECTION: ln-tau correlation
    def cal_tau_ij_M2(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        d_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC tau from a four-coefficient ln(tau) correlation.

        ln(tau_ij) = a_ij + b_ij/T + c_ij*ln(T) + d_ij*T
        tau_ij = exp(ln(tau_ij))
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij, d_ij)

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        tau_ij[row, col] = 0
                        tau_ij_comp[key] = 0
                        continue

                    # NOTE: UNIQUAC M2 evaluates ln(tau_ij), then exponentiates
                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key)
                        * log(temperature)
                        + self._ij_value(d_ij, "d", i, j, key) * temperature
                    )
                    value = exp(value)
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M2: {str(e)}")
