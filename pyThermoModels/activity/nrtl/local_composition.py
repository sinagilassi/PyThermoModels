# SECTION: imports
from math import log, pow
from typing import Dict, Literal, Tuple
import numpy as np
from pyThermoDB import TableMatrixData
# locals
from ..local_composition_base import IJData, LocalCompositionBase


class NRTLLocalComposition(LocalCompositionBase):
    """NRTL local-composition interaction-parameter calculations.

    NOTE: In NRTL, tau_ij is dimensionless and defined directly by
    dg_ij / (R*T) or by a declared tau-correlation. The exponential term is
    applied later to compute G_ij, not to redefine tau_ij.
    """

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
        Calculate NRTL interaction-energy difference matrix.

        SECTION: equation
        dg_ij = a_ij + b_ij*T + c_ij*T^2

        NOTE: dg_ij is typically interpreted as energy per mole.
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_same_source_type(a_ij, b_ij, c_ij)

            # SECTION: initialize result containers
            dg_ij = np.zeros((self.comp_num, self.comp_num))
            dg_ij_comp = {}

            # SECTION: calculate dg_ij
            temperature_sq = pow(temperature, 2)

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
                        + self._ij_value(c_ij, "c", i, j, key) * temperature_sq
                    )
                    self._validate_numeric_result(value, f"dg_ij[{key}]")
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
        Calculate NRTL tau from interaction-energy difference.

        SECTION: equation
        tau_ij = dg_ij / (R*T)

        NOTE: tau_ij can be negative in NRTL and remains physically acceptable
        as a dimensionless energy difference.
        """
        try:
            # SECTION: validation
            self._validate_temperature(temperature)
            self._validate_r_const(R_CONST)
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
                    self._validate_numeric_result(value, f"tau_ij[{key}]")
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

        SECTION: equation
        tau_ij = a_ij + b_ij/T + c_ij*ln(T) + d_ij*T

        NOTE: this method computes tau_ij directly (dimensionless).
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
                    self._validate_numeric_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M2: {str(e)}")

    # SECTION: direct tau correlation (A + B/T)
    def cal_tau_ij_M3(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau using a two-coefficient direct correlation.

        SECTION: equation
        tau_ij = a_ij + b_ij / T

        NOTE: this matches common NRTL direct-tau forms where tau_ij is
        dimensionless and evaluated without exponential mapping.
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

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        tau_ij[row, col] = 0
                        tau_ij_comp[key] = 0
                        continue

                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                    )
                    self._validate_numeric_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M3: {str(e)}")

    # SECTION: direct tau correlation (A + B/T + C/T^2)
    def cal_tau_ij_M4(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau using an inverse-temperature polynomial.

        SECTION: equation
        tau_ij = a_ij + b_ij / T + c_ij / T^2

        NOTE: this method corresponds to NRTL direct-tau polynomial forms.
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

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        tau_ij[row, col] = 0
                        tau_ij_comp[key] = 0
                        continue

                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key) / temperature_sq
                    )
                    self._validate_numeric_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M4: {str(e)}")

    # SECTION: direct tau correlation (A + B/T + C*ln(T))
    def cal_tau_ij_M5(
        self,
        temperature: float,
        a_ij: IJData,
        b_ij: IJData,
        c_ij: IJData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau using a logarithmic direct correlation.

        SECTION: equation
        tau_ij = a_ij + b_ij / T + c_ij * ln(T)

        NOTE: this method preserves NRTL semantics where tau_ij is directly
        computed and may be negative.
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

                    # NOTE: self-interaction terms are defined as zero
                    if i == j:
                        tau_ij[row, col] = 0
                        tau_ij_comp[key] = 0
                        continue

                    value = (
                        self._ij_value(a_ij, "a", i, j, key)
                        + self._ij_value(b_ij, "b", i, j, key) / temperature
                        + self._ij_value(c_ij, "c", i, j, key) * log_t
                    )
                    self._validate_numeric_result(value, f"tau_ij[{key}]")
                    tau_ij[row, col] = value
                    tau_ij_comp[key] = value

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M5: {str(e)}")

    # NOTE: universal tau calculation method

    def cal_tau_ij(
        self,
        temperature: float,
        dg_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        d_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        tau_correlation: Literal[
            'M1', 'M2', 'M3', 'M4', 'M5'
        ],
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
    ):
        """
        Calculate interaction parameters `tau_ij` matrix for NRTL model based on the selected correlation method.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin [K].
        dg_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction energy parameter [J/mol] matrix where dg_ij[i][j] between component i and j.
        a_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
        b_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
        c_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
        d_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter d_ij matrix where d_ij[i][j] between component i and j.
        tau_correlation : Literal['M1', 'M2', 'M3', 'M4', 'M5']
            Correlation method to use for calculating tau_ij. Default is 'M1'.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        tau_ij : np.ndarray
            interaction parameters `tau_ij` matrix for NRTL model.
        tau_ij_comp : dict
            Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.

        Notes
        -----
        1. The calculation method for tau_ij depends on the selected correlation method (M1, M2, M3, M4, M5).
        2. All parameters including a_ij, b_ij, c_ij, d_ij must be in the same format (numpy array, dict or TableMatrixData).
        3. The equations for each correlation method are as follows:

            M1: tau_ij = dg_ij / (R * T)
            M2: tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T
            M3: tau_ij = a_ij + b_ij / T
            M4: tau_ij = a_ij + b_ij / T + c_ij / T^2
            M5: tau_ij = a_ij + b_ij / T + c_ij * ln(T)
        """
        try:
            if tau_correlation == 'M1':
                return self.cal_tau_ij_M1(
                    temperature=temperature,
                    dg_ij=dg_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif tau_correlation == 'M2':
                return self.cal_tau_ij_M2(
                    temperature=temperature,
                    a_ij=a_ij,
                    b_ij=b_ij,
                    c_ij=c_ij,
                    d_ij=d_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif tau_correlation == 'M3':
                return self.cal_tau_ij_M3(
                    temperature=temperature,
                    a_ij=a_ij,
                    b_ij=b_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif tau_correlation == 'M4':
                return self.cal_tau_ij_M4(
                    temperature=temperature,
                    a_ij=a_ij,
                    b_ij=b_ij,
                    c_ij=c_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif tau_correlation == 'M5':
                return self.cal_tau_ij_M5(
                    temperature=temperature,
                    a_ij=a_ij,
                    b_ij=b_ij,
                    c_ij=c_ij,
                    symbol_delimiter=symbol_delimiter
                )
            else:
                raise ValueError(
                    f"Unsupported tau_correlation method: {tau_correlation}")
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij: {str(e)}")

