# SECTION: imports
from math import log, pow
from typing import Dict, Literal, Tuple

import numpy as np

from .local_composition_base import IJData, LocalCompositionBase


class NRTLLocalComposition(LocalCompositionBase):
    """NRTL local-composition interaction-parameter calculations."""

    # SECTION: interaction energy
    def cal_dg_ij_M1(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL interaction-energy parameter.

        dg_ij = a_ij + b_ij*T + c_ij*T^2
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij)

            # SECTION: initialize result containers
            dg_ij = np.zeros((self.comp_num, self.comp_num))
            dg_ij_comp = {}

            # SECTION: calculate dg_ij
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        dg_ij[row, col] = 0
                        dg_ij_comp[key] = 0
                        continue

                    # NOTE: NRTL M1 energy correlation
                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) * temperature
                        + self._ij_value(c_ij, "c", i, j, key)
                        * pow(temperature, 2)
                    )
                    dg_ij[row, col] = value
                    dg_ij_comp[key] = value

            return dg_ij, dg_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_dg_ij_M1: {str(e)}")

    # SECTION: tau from energy
    def cal_tau_ij_M1(
        self,
        temperature: float,
        dg_ij: IJData,
        dg_ij_symbol: Literal["dg", "dg_ij"] = "dg",
        R_CONST: float = 8.314,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau from interaction energy.

        tau_ij = dg_ij / (R*T)
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_ij_data(dg_ij, "dg_ij")

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

                    # NOTE: NRTL tau is a dimensionless energy difference
                    value = self._ij_value(dg_ij, dg_ij_symbol, i, j, key)
                    value = value / (R_CONST * temperature)
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M1: {str(e)}")

    # SECTION: direct tau correlation
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
        Calculate NRTL tau using a direct four-coefficient correlation.

        tau_ij = a_ij + b_ij/T + c_ij*ln(T) + d_ij*T
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

                    # NOTE: NRTL M2 evaluates tau_ij directly
                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key)
                        * log(temperature)
                        + self._ij_value(d_ij, "d", i, j, key) * temperature
                    )
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M2: {str(e)}")
