# SECTION: imports
from math import log, pow
from typing import Dict, Literal, Tuple

import numpy as np

from ..local_composition_base import IJData, LocalCompositionBase


class UNIQUACLocalComposition(LocalCompositionBase):
    """UNIQUAC local-composition interaction-parameter calculations.

    NOTE: In UNIQUAC, tau_ij is an exponential weighting term. If an energy
    correlation is used, tau_ij = exp(-dU_ij/(R*T)); if an ln(tau_ij)
    correlation is used, tau_ij = exp(ln(tau_ij)).
    """

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
        Calculate UNIQUAC interaction-energy difference matrix.

        SECTION: equation
        dU_ij = a_ij + b_ij*T + c_ij*T^2

        NOTE: dU_ij is typically interpreted as energy per mole.
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij)

            # SECTION: initialize result containers
            dU_ij = np.zeros((self.comp_num, self.comp_num))
            dU_ij_comp = {}

            # SECTION: calculate dU_ij
            temperature_sq = pow(temperature, 2)

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
                        * temperature_sq
                    )
                    self._validate_numeric_result(value, f"dU_ij[{key}]")
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

        SECTION: equation
        tau_ij = exp(-dU_ij / (R*T))

        NOTE: tau_ij must be strictly positive in UNIQUAC, and self
        interactions evaluate to exp(0) = 1.
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_r_const(R_CONST)
            self._validate_ij_data(dU_ij, "dU_ij")

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: UNIQUAC self-interaction tau terms are unity
                    if i == j:
                        tau_ij[row, col] = 1.0
                        tau_ij_comp[key] = 1.0
                        continue

                    # NOTE: UNIQUAC tau is an exponential weighting factor
                    value = self._ij_value(dU_ij, dU_ij_symbol, i, j, key)
                    exponent = -value / (R_CONST * temperature)
                    value = self._exp_checked(exponent, f"tau_ij[{key}]")
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

        SECTION: equation
        ln(tau_ij) = a_ij + b_ij/T + c_ij*ln(T) + d_ij*T
        tau_ij = exp(ln(tau_ij))

        NOTE: the full ln(tau_ij) expression must be evaluated before
        exponentiation. Self interactions evaluate to exp(0) = 1.
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

                    # NOTE: UNIQUAC self-interaction tau terms are unity
                    if i == j:
                        tau_ij[row, col] = 1.0
                        tau_ij_comp[key] = 1.0
                        continue

                    # NOTE: UNIQUAC M2 evaluates ln(tau_ij), then exponentiates
                    ln_tau = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key)
                        * log(temperature)
                        + self._ij_value(d_ij, "d", i, j, key) * temperature
                    )
                    self._validate_numeric_result(ln_tau, f"ln_tau_ij[{key}]")
                    value = self._exp_checked(ln_tau, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M2: {str(e)}")

    # SECTION: tau from (dU/R) correlation
    def cal_tau_ij_M3(
        self,
        temperature: float,
        a_ij: IJData,
        sign: Literal["negative", "positive"] = "negative",
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC tau from an energy-over-R correlation.

        SECTION: equation
        if sign == "negative": tau_ij = exp(-a_ij / T)
        if sign == "positive": tau_ij = exp(+a_ij / T)

        NOTE: use sign='positive' when the source equation follows
        -dU_ij / R = a_ij.
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_ij_data(a_ij, "a_ij")

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            sign_factor = -1.0 if sign == "negative" else 1.0
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: UNIQUAC self-interaction tau terms are unity
                    if i == j:
                        tau_ij[row, col] = 1.0
                        tau_ij_comp[key] = 1.0
                        continue

                    value = self._ij_value(a_ij, "a", i, j, key)
                    exponent = sign_factor * value / temperature
                    value = self._exp_checked(exponent, f"tau_ij[{key}]")
                    self._validate_positive_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M3: {str(e)}")

    # SECTION: ln-tau correlation (A + B/T)
    def cal_tau_ij_M4(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC tau from a two-coefficient ln(tau) correlation.

        SECTION: equation
        ln(tau_ij) = a_ij + b_ij / T
        tau_ij = exp(ln(tau_ij))
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij)

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: UNIQUAC self-interaction tau terms are unity
                    if i == j:
                        tau_ij[row, col] = 1.0
                        tau_ij_comp[key] = 1.0
                        continue

                    ln_tau = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                    )
                    self._validate_numeric_result(ln_tau, f"ln_tau_ij[{key}]")
                    value = self._exp_checked(ln_tau, f"tau_ij[{key}]")
                    self._validate_positive_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M4: {str(e)}")

    # SECTION: ln-tau correlation (A + B/T + C/T^2)
    def cal_tau_ij_M5(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC tau from an inverse-temperature ln(tau) polynomial.

        SECTION: equation
        ln(tau_ij) = a_ij + b_ij / T + c_ij / T^2
        tau_ij = exp(ln(tau_ij))
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij)

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            temperature_sq = pow(temperature, 2)
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: UNIQUAC self-interaction tau terms are unity
                    if i == j:
                        tau_ij[row, col] = 1.0
                        tau_ij_comp[key] = 1.0
                        continue

                    ln_tau = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key) / temperature_sq
                    )
                    self._validate_numeric_result(ln_tau, f"ln_tau_ij[{key}]")
                    value = self._exp_checked(ln_tau, f"tau_ij[{key}]")
                    self._validate_positive_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M5: {str(e)}")

    # SECTION: ln-tau correlation (A + B/T + C*ln(T))
    def cal_tau_ij_M6(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate UNIQUAC tau from a logarithmic ln(tau) correlation.

        SECTION: equation
        ln(tau_ij) = a_ij + b_ij / T + c_ij * ln(T)
        tau_ij = exp(ln(tau_ij))
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij)

            # SECTION: initialize result containers
            tau_ij = np.zeros((self.comp_num, self.comp_num), dtype=float)
            tau_ij_comp = {}

            # SECTION: calculate tau_ij
            log_t = log(temperature)
            for i in range(self.comp_num):
                for j in range(self.comp_num):
                    key = self._pair_key(i, j, symbol_delimiter)
                    row, col = self._component_position(i, j)

                    # NOTE: UNIQUAC self-interaction tau terms are unity
                    if i == j:
                        tau_ij[row, col] = 1.0
                        tau_ij_comp[key] = 1.0
                        continue

                    ln_tau = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key) * log_t
                    )
                    self._validate_numeric_result(ln_tau, f"ln_tau_ij[{key}]")
                    value = self._exp_checked(ln_tau, f"tau_ij[{key}]")
                    self._validate_positive_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M6: {str(e)}")

