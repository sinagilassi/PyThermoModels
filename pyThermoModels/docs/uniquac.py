# import libs
import numpy as np
import json
import os
from math import pow, exp, log
from typing import List, Dict, Tuple, Any, Literal, Optional, Union
from rich import print
import pyThermoDB
from pyThermoDB import (
    TableMatrixData, TableData, TableEquation, TableMatrixEquation
)


class UNIQUAC:
    """
    The UNIQUAC (`Universal Quasi-Chemical`) model - a thermodynamic framework used to describe the behavior of mixtures,
    particularly in the context of phase equilibria and activity coefficients

    To apply the UNIQUAC model, you'll need the following parameters:

    **Pure Component Parameters**
    - r_i (volume parameter): represents the volume of a molecule in the mixture.
    - q_i (surface area parameter): represents the surface area of a molecule in the mixture.

    **Binary Interaction Parameters**
    - Δu_ij (interaction energy parameter): represents the interaction energy between two molecules [J/mol].
    - τ_ij (binary interaction parameter): represents the interaction energy between two molecules of different components [dimensionless].

    Universal gas constant (R) is defined as 8.314 J/mol/K.

    Z is a constant used in the model, default value is 10.0.
    """

    # universal gas constant [J/mol/K]
    R_CONST = 8.314
    # constant
    Z = 10.0

    def __init__(self,
                 components: List[str],
                 datasource: Dict = {},
                 equationsource: Dict = {}):
        '''
        Initialize the activity model, UNIQUAC (`Universal Quasi-Chemical`) used to calculate activity coefficients in liquid mixtures.

        Parameters
        ----------
        datasource: Dict
            Data source for the model
        equationsource: Dict
            Equation source for the model
        components: List[str]
            List of component names in the mixture

        Raises
        ------
        TypeError
        - If datasource is not a dict
        - If equationsource is not a dict

        Notes
        -----
        The model needs the following parameters:
        - datasource: Data source for the model
        - equationsource: Equation source for the model
        - components: List of component names in the mixture

        The component names define the order of the parameters in the model. The first component in the list is component 1, the second is component 2, and so on.

        Universal gas constant is defined as 8.314 J/mol/K.

        Z is a constant used in the model, default value is 10.0.

        UNIQUAC model parameters are defined as:
        - dU_ij: Interaction energy parameter [J/mol]
        - tau_ij: Interaction parameter [dimensionless]

        The tau_ij equation is defined as:
        - tau_ij = exp(-dU_ij / (R * T))
        - tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T
        - tau_ij = exp(a_ij + b_ij/T)

        The dU_ij equation is defined as:
        - dU_ij = a_ij + b_ij * T + c_ij * T^2
        '''
        # SECTION:
        # Check datasource
        if not isinstance(datasource, dict):
            raise TypeError("datasource must be a dict")

        # Check equationsource
        if not isinstance(equationsource, dict):
            raise TypeError("equationsource must be a dict")

        # Check if components is a list
        if not isinstance(components, list):
            raise TypeError("components must be a list")

        # SECTION: Assign the parameters to instance variables
        self.datasource = datasource
        self.equationsource = equationsource
        self.components = [components.strip() for components in components]

        # SECTION
        # Get the number of components
        self.comp_num = len(components)
        # idx
        self.comp_idx = {components[i]: i for i in range(self.comp_num)}

    def __str__(self):
        model_ = """
        The UNIQUAC (`Universal Quasi-Chemical`) model - a thermodynamic framework used to describe the behavior of mixtures,
        particularly in the context of phase equilibria and activity coefficients

        To apply the UNIQUAC model, you'll need the following parameters:

        **Pure Component Parameters**
        - r_i (volume parameter): represents the volume of a molecule in the mixture.
        - q_i (surface area parameter): represents the surface area of a molecule in the mixture.

        **Binary Interaction Parameters**
        - Δu_ij (interaction energy parameter): represents the interaction energy between two molecules [J/mol].
        - τ_ij (binary interaction parameter): represents the interaction energy between two molecules of different components [dimensionless].

        Universal gas constant (R) is defined as 8.314 J/mol/K.

        Z is a constant used in the model, default value is 10.0.
        """
        return model_

    def to_ij(self, data: TableMatrixData, prop_symbol: str,
              symbol_delimiter: Literal["|", "_"] = "|") -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Convert TableMatrixData to numpy array (mat_ij) and dictionary (dict_ij).

        Parameters
        ----------
        data : TableMatrixData
            Parameter dictionary for data[i][j] between component i and j.
        prop_symbol : str
            Property symbol for the parameter.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        mat_ij : np.ndarray
            Parameter matrix (numpy array).
        dict_ij : dict
            Dictionary of parameters where keys are component pairs and values are their respective values.
        """
        try:
            # NOTE: check dU_ij is TableMatrixData
            if not isinstance(data, TableMatrixData):
                raise TypeError("dict_ij_src must be TableMatrixData")

            # Get the number of components
            comp_num = self.comp_num

            # prop symbol
            prop_symbol = prop_symbol.strip()

            # Initialize
            mat_ij = np.zeros((comp_num, comp_num))
            dict_ij = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # Set the interaction energy parameter matrix
            for i in range(comp_num):
                for j in range(comp_num):
                    # val
                    val = data.ij(
                        f"{prop_symbol}_{self.components[i]}-{self.components[j]}")
                    # to matrix
                    mat_ij[i, j] = val

                    # to dict
                    key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"
                    dict_ij[key_] = val

            # res
            return mat_ij, dict_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_i(self, data: Dict[str, float]):
        """
        Convert data to numpy array with respect to component id.

        Parameters
        ----------
        data : Dict[str, float]
            Parameter dictionary where keys are component names and values are their respective values.

        Returns
        -------
        data_i : np.ndarray
            Parameter array (numpy array) with respect to component id.
        """
        try:
            # NOTE: check
            if not isinstance(data, dict):
                raise TypeError("data must be dict")

            # Get the number of components
            comp_num = self.comp_num

            # Initialize
            data_i = np.zeros(comp_num)

            # Set the interaction energy parameter matrix
            for i in range(comp_num):
                # check if the component id is in the data dictionary
                if self.components[i] not in data:
                    raise KeyError(
                        f"Component {self.components[i]} not found in data dictionary")

                # val
                val = data[self.components[i]]

                # to matrix
                data_i[i] = val

            # res
            return data_i
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_dict_ij(self, data: np.ndarray, symbol_delimiter: Literal["|", "_"] = "|") -> Dict[str, float]:
        """
        Convert to dictionary (dict_ij) according to the component id.

        Parameters
        ----------
        data : np.ndarray
            Parameter matrix (numpy array) with respect to component id.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        dict_ij : Dict[str, float]
            Dictionary of parameters where keys are component pairs and values are their respective values.
        """
        try:
            # NOTE: check
            if not isinstance(data, np.ndarray):
                raise TypeError("data must be numpy array")

            # Get the number of components
            comp_num = self.comp_num

            # Initialize
            dict_ij = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # Set the interaction energy parameter matrix
            for i in range(comp_num):
                for j in range(comp_num):
                    # val
                    val = data[i, j]

                    # to dict
                    key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"
                    dict_ij[key_] = val

            # res
            return dict_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_dict_i(self, data: List[float] | np.ndarray) -> Dict[str, float]:
        """
        Convert to dictionary (dict_i) according to the component id.

        Parameters
        ----------
        data : List[float] | np.ndarray
            Parameter list or numpy 1d array with respect to component id.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        dict_i : Dict[str, float]
            Dictionary of parameters where keys are component pairs and values are their respective values.
        """
        try:
            # NOTE: check
            if not isinstance(data, np.ndarray) and not isinstance(data, list):
                raise TypeError("data must be numpy array or list")

            # Get the number of components
            comp_num = self.comp_num

            # Initialize
            dict_i = {}

            # Set the interaction energy parameter matrix
            for i in range(comp_num):
                # val
                val = data[i]

                # to dict
                key_ = self.components[i].strip()
                dict_i[key_] = val

            # res
            return dict_i
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_matrix_ij(self, data: Dict[str, float] | List[float], symbol_delimiter: Literal["|", "_"] = "|") -> np.ndarray:
        """
        Convert to matrix (mat_ij) according to `the component id`.

        Parameters
        ----------
        data : Dict[str, float] | List[float]
            Dictionary of parameters where keys are component pairs and values are their respective values or list of values according to the component id.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        mat_ij : np.ndarray
            Parameter matrix (numpy array).
        """
        try:
            # NOTE: check
            if not isinstance(data, dict):
                raise TypeError("data must be dict")

            # Get the number of components
            comp_num = self.comp_num

            # Initialize
            mat_ij = np.zeros((comp_num, comp_num))

            # NOTE: check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # Set the interaction energy parameter matrix
            # SECTION:
            if isinstance(data, Dict):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # val
                        val = data[f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"]

                        # find the component id
                        comp_id_i = self.comp_idx[self.components[i]]
                        comp_id_j = self.comp_idx[self.components[j]]

                        # to matrix
                        mat_ij[comp_id_i, comp_id_j] = val
            elif isinstance(data, List):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # val
                        val = data[i][j]

                        # find the component id
                        comp_id_i = self.comp_idx[self.components[i]]
                        comp_id_j = self.comp_idx[self.components[j]]

                        # to matrix
                        mat_ij[comp_id_i, comp_id_j] = val

            # res
            return mat_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def cal_dU_ij_M1(self, temperature: float,
                     a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                     b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                     c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                     symbol_delimiter: Literal["|", "_"] = "|") -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate interaction energy parameter `dU_ij` matrix dependent of temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin [K].
        a_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
        b_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
        c_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        dU_ij : np.ndarray
            Interaction energy parameter `dU_ij` matrix for UNIQUAC model.
        dU_ij_comp : dict
            Dictionary of interaction energy parameters where keys are component pairs and values are their respective dU_ij values.

        Notes
        -----
        1. The interaction energy parameter matrix is calculated using the formula:

            `dU_ij = a_ij + b_ij * T + c_ij * T^2`

            where T is the temperature [K].

        2. All parameters including a_ij, b_ij, c_ij must be in the same format (numpy array, dict or TableMatrixData).
        """
        try:
            # SECTION: check
            if not isinstance(a_ij, np.ndarray) and not isinstance(a_ij, dict) and not isinstance(a_ij, TableMatrixData):
                raise TypeError(
                    "a_ij must be numpy array, dict or TableMatrixData")

            if not isinstance(b_ij, np.ndarray) and not isinstance(b_ij, dict) and not isinstance(b_ij, TableMatrixData):
                raise TypeError(
                    "b_ij must be numpy array, dict or TableMatrixData")

            if not isinstance(c_ij, np.ndarray) and not isinstance(c_ij, dict) and not isinstance(c_ij, TableMatrixData):
                raise TypeError(
                    "c_ij must be numpy array, dict or TableMatrixData")

            # Get the number of components
            comp_num = self.comp_num

            # Initialize dU_ij matrix
            dU_ij = np.zeros((comp_num, comp_num))

            # dU_ij components
            dU_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # SECTION: calculate dU_ij values
            if isinstance(a_ij, np.ndarray) and isinstance(b_ij, np.ndarray) and isinstance(c_ij, np.ndarray):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"

                        # check
                        if i != j:
                            # val
                            val_ = a_ij[i, j] + b_ij[i, j] * \
                                temperature + c_ij[i, j] * pow(temperature, 2)
                            # set
                            dU_ij[i, j] = val_
                            # set by name
                            dU_ij_comp[key_] = val_
                        else:
                            # set
                            dU_ij[i, j] = 0
                            # set by name
                            dU_ij_comp[key_] = 0

            # SECTION: if dU_ij is dict
            elif isinstance(a_ij, dict) and isinstance(b_ij, dict) and isinstance(c_ij, dict):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"

                        # component id
                        comp_id_i = self.comp_idx[self.components[i]]
                        comp_id_j = self.comp_idx[self.components[j]]

                        # check
                        if i != j:
                            # val
                            val_ = a_ij[key_] + b_ij[key_] * \
                                temperature + c_ij[key_] * pow(temperature, 2)
                            # set
                            dU_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            dU_ij_comp[key_] = val_
                        else:
                            # set
                            dU_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            dU_ij_comp[key_] = 0
            # SECTION: if dU_ij is TableMatrixData
            elif isinstance(a_ij, TableMatrixData) and isinstance(b_ij, TableMatrixData) and isinstance(c_ij, TableMatrixData):
                # convert to numpy array and dict
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{self.components[i]}_{self.components[j]}"
                        # dict
                        key_comp = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"

                        # val
                        val_ = a_ij.ij(f"a_{key_}") + b_ij.ij(f"b_{key_}") * \
                            temperature + \
                            c_ij.ij(f"c_{key_}") * pow(temperature, 2)

                        # component id
                        comp_id_i = self.comp_idx[self.components[i]]
                        comp_id_j = self.comp_idx[self.components[j]]

                        # set
                        if i != j:
                            # set
                            dU_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            dU_ij_comp[key_comp] = val_
                        else:
                            # set
                            dU_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            dU_ij_comp[key_comp] = 0
            else:
                raise TypeError(
                    "a_ij, b_ij and c_ij must be numpy array or dict")
            # res
            return dU_ij, dU_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_dU_ij_M1: {str(e)}")

    def cal_tau_ij_M1(self, temperature: float, dU_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                      dU_ij_symbol: Literal['dU', 'dU_ij'] = 'dU', R_CONST: float = 8.314,
                      symbol_delimiter: Literal["|", "_"] = "|") -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate interaction parameters `tau_ij` matrix for UNIQUAC model.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin [K].
        dU_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction energy parameter [J/mol] matrix where dU_ij[i][j] between component i and j.
        dU_ij_symbol : str
            Interaction energy parameter symbol. Default is 'dU'.
        R_CONST : float
            Univeral gas constant [J/mol/K], default R_CONST = 8.314
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        tau_ij : np.ndarray
            interaction parameters `tau_ij` matrix for UNIQUAC model.
        tau_ij_comp : dict
            Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.

        Notes
        -----
        1. The tau_ij matrix is calculated using the formula:

            `tau_ij = exp(-dU_ij / (R * T))`

            where R is the universal gas constant [J/mol/K] and T is the temperature [K].

        2. Interaction energy parameter symbol is `dU` for TableMatrixData as:

        - `dU_{component_i}-{component_j}`
        - `dU | {component_i} | {component_j}`.
        """
        try:
            # check
            if not isinstance(dU_ij, np.ndarray) and not isinstance(dU_ij, dict) and not isinstance(dU_ij, TableMatrixData):
                raise TypeError(
                    "dU_ij must be numpy array, dict or TableMatrixData")

            # Get the number of components
            comp_num = self.comp_num

            # components
            components = self.components

            # Initialize tauij matrix
            tau_ij = np.zeros((comp_num, comp_num))

            # tauij components
            tau_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # Calculate tauij values
            # SECTION: if dU_ij is numpy array
            if isinstance(dU_ij, np.ndarray):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # check
                        if i != j:
                            # val
                            val_ = exp(-1*dU_ij[i, j] /
                                       (R_CONST * temperature))
                            # set
                            tau_ij[i, j] = val_
                            # set by name
                            tau_ij_comp[key_] = val_
                        else:
                            # set
                            tau_ij[i, j] = 0
                            # set by name
                            tau_ij_comp[key_] = 0
            # SECTION: if dU_ij is dict
            elif isinstance(dU_ij, dict):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # component id
                        comp_id_i = self.comp_idx[components[i]]
                        comp_id_j = self.comp_idx[components[j]]

                        # check
                        if i != j:
                            # val
                            val_ = exp(-1*dU_ij[key_] /
                                       (R_CONST * temperature))
                            # set
                            tau_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            tau_ij_comp[key_] = val_
                        else:
                            # set
                            tau_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            tau_ij_comp[key_] = 0
            # SECTION: if dU_ij is TableMatrixData
            elif isinstance(dU_ij, TableMatrixData):
                # convert to numpy array and dict
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{dU_ij_symbol}_{components[i]}_{components[j]}"
                        # dict
                        key_comp = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # val
                        val_ = dU_ij.ij(key_)

                        # component id
                        comp_id_i = self.comp_idx[components[i]]
                        comp_id_j = self.comp_idx[components[j]]

                        # set
                        if i != j:
                            # val
                            val_ = exp(-1*val_ / (R_CONST * temperature))
                            # set
                            tau_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            tau_ij_comp[key_comp] = val_
                        else:
                            # set
                            tau_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            tau_ij_comp[key_comp] = 0
            else:
                raise TypeError("dU_ij must be numpy array or dict")

            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tauij: {str(e)}")

    def cal_tau_ij_M2(self, temperature: float,
                      a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                      b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                      c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                      d_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                      symbol_delimiter: Literal["|", "_"] = "|") -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate interaction parameters `tau_ij` matrix for UNIQUAC model.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin [K].
        a_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter a_ij matrix where a_ij[i][j] between component i and j.
        b_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter b_ij matrix where b_ij[i][j] between component i and j.
        c_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter c_ij matrix where c_ij[i][j] between component i and j.
        d_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction parameter d_ij matrix where d_ij[i][j] between component i and j.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        tau_ij : np.ndarray
            interaction parameters `tau_ij` matrix for UNIQUAC model.
        tau_ij_comp : dict
            Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.

        Notes
        -----
        1. The extended Antoine equation format is used to calculate the interaction parameters using the following formula:

            `tau_ij = a_ij + b_ij / T + c_ij * ln(T) + d_ij * T`

        2. Interaction energy parameter symbol is `X` for TableMatrixData as:

        - `X_{component_i}-{component_j}`
        - `X | {component_i} | {component_j}`.

        2. All parameters including a_ij, b_ij, c_ij, d_ij must be in the same format (numpy array, dict or TableMatrixData).
        """
        try:
            # SECTION: check
            if not isinstance(a_ij, np.ndarray) and not isinstance(a_ij, dict) and not isinstance(a_ij, TableMatrixData):
                raise TypeError(
                    "a_ij must be numpy array, dict or TableMatrixData")

            if not isinstance(b_ij, np.ndarray) and not isinstance(b_ij, dict) and not isinstance(b_ij, TableMatrixData):
                raise TypeError(
                    "b_ij must be numpy array, dict or TableMatrixData")

            if not isinstance(c_ij, np.ndarray) and not isinstance(c_ij, dict) and not isinstance(c_ij, TableMatrixData):
                raise TypeError(
                    "c_ij must be numpy array, dict or TableMatrixData")

            if not isinstance(d_ij, np.ndarray) and not isinstance(d_ij, dict) and not isinstance(d_ij, TableMatrixData):
                raise TypeError(
                    "d_ij must be numpy array, dict or TableMatrixData")

            # Get the number of components
            comp_num = self.comp_num

            # components
            components = self.components

            # Initialize tauij matrix
            tau_ij = np.zeros((comp_num, comp_num))

            # tauij components
            tau_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # SECTION: Calculate tauij values
            if isinstance(a_ij, np.ndarray) and isinstance(b_ij, np.ndarray) and isinstance(c_ij, np.ndarray) and isinstance(d_ij, np.ndarray):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # check
                        if i != j:
                            # val
                            val_ = a_ij[i, j] + b_ij[i, j] / temperature + \
                                c_ij[i, j] * log(temperature) + \
                                d_ij[i, j] * temperature
                            # set
                            tau_ij[i, j] = val_
                            # set by name
                            tau_ij_comp[key_] = val_
                        else:
                            # set
                            tau_ij[i, j] = 0
                            # set by name
                            tau_ij_comp[key_] = 0
            # SECTION: if dU_ij is dict
            elif isinstance(a_ij, dict) and isinstance(b_ij, dict) and isinstance(c_ij, dict) and isinstance(d_ij, dict):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # component id
                        comp_id_i = self.comp_idx[components[i]]
                        comp_id_j = self.comp_idx[components[j]]

                        # check
                        if i != j:
                            # val
                            val_ = a_ij[key_] + b_ij[key_] / temperature + \
                                c_ij[key_] * log(temperature) + \
                                d_ij[key_] * temperature
                            # set
                            tau_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            tau_ij_comp[key_] = val_
                        else:
                            # set
                            tau_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            tau_ij_comp[key_] = 0
            # SECTION: if dU_ij is TableMatrixData
            elif isinstance(a_ij, TableMatrixData) and isinstance(b_ij, TableMatrixData) and isinstance(c_ij, TableMatrixData) and isinstance(d_ij, TableMatrixData):
                # convert to numpy array and dict
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"X_{components[i]}_{components[j]}"
                        # dict
                        key_comp = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # val
                        val_ = a_ij.ij(key_) + b_ij.ij(key_) / temperature + c_ij.ij(
                            key_) * log(temperature) + d_ij.ij(key_) * temperature

                        # component id
                        comp_id_i = self.comp_idx[components[i]]
                        comp_id_j = self.comp_idx[components[j]]

                        # set
                        if i != j:
                            # set
                            tau_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            tau_ij_comp[key_comp] = val_
                        else:
                            # set
                            tau_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            tau_ij_comp[key_comp] = 0
            else:
                raise TypeError(
                    "a_ij, b_ij, c_ij and d_ij must be numpy array or dict")

            # res
            return tau_ij, tau_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_tauij: {str(e)}")

    def cal(self,
            model_input: Dict,
            Z: Optional[float | int] = None,
            calculation_mode: Literal['V1'] = 'V1',
            symbol_delimiter: Literal["|",
                                      "_"] = "|",
            message: Optional[str] = None,
            res_format: Literal['dict', 'str', 'json'] = 'dict') -> Tuple[Dict[str, str | float | Dict], Dict[str, str | float | Dict]] | str:
        """
        Calculate activity coefficients for a multi-component mixture using the UNIQUAC model.

        Parameters
        -----------
        model_input : Dict
            Dictionary of input values where keys are component names and values are their respective values.
                - `mole_fraction`: Dict[str, float]
                    dictionary of mole fractions where keys are component names and values are their respective mole fractions.
                - `tau_ij` : TableMatrixData | np.ndarray | Dict[str, float]
                    Interaction parameters (tau_ij) between component i and j.
                - `r_i` : List[float] | Dict[str, float]
                    relative van der Waals volume of component i
                - `q_i` : List[float] | Dict[str, float]
                    relative surface area of component i
        Z : int | float
            Model constant, Default is 10.
        calculation_mode : Literal['V1', 'V2']
            Mode of calculation. If 'V1', use the first version of the UNIQUAC model. If 'V2', use the second version.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".
        message : Optional[str]
            Message to be displayed. Default is None.
        res_format : Literal['dict', 'str', 'json']
            Format of the result. Default is 'dict'.
        """
        try:
            # SECTION: check
            if not isinstance(model_input, dict):
                raise TypeError("model_input must be dict")

            # SECTION: check model input
            if 'mole_fraction' not in model_input:
                raise KeyError("mole_fraction is required in model_input")
            if 'tau_ij' not in model_input:
                raise KeyError("tau_ij is required in model_input")
            if 'r_i' not in model_input:
                raise KeyError("r_i is required in model_input")
            if 'q_i' not in model_input:
                raise KeyError("q_i is required in model_input")

            # NOTE: get values
            # mole fraction
            mole_fraction = model_input['mole_fraction']
            # tau_ij
            tau_ij = model_input['tau_ij']
            # r_i
            r_i = model_input['r_i']
            # q_i
            q_i = model_input['q_i']

            # SECTION: calculate activity coefficients
            return self.__calculate_activity_coefficients(
                mole_fraction=mole_fraction,
                tau_ij_data=tau_ij,
                r_i_data=r_i,
                q_i_data=q_i,
                Z=Z,
                calculation_mode=calculation_mode,
                symbol_delimiter=symbol_delimiter,
                message=message,
                res_format=res_format)
        except Exception as e:
            raise Exception(f"Error in uniquac model cal: {str(e)}")

    def __calculate_activity_coefficients(self,
                                          mole_fraction: Dict[str, float],
                                          tau_ij_data: TableMatrixData | np.ndarray | Dict[str, float],
                                          r_i_data: List[float] | Dict[str, float],
                                          q_i_data: List[float] | Dict[str, float],
                                          Z: Optional[float | int],
                                          calculation_mode: Literal['V1'],
                                          symbol_delimiter: Literal["|", "_"],
                                          message: Optional[str],
                                          res_format: Literal['dict', 'str', 'json']) -> Tuple[Dict[str, str | float | Dict], Dict[str, str | float | Dict]] | str:
        """
        Calculate activity coefficients for a multi-component mixture using the UNIQUAC model.

        Parameters
        -----------
        mole_fraction : Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        tau_ij_comp : TableMatrixData | np.ndarray | Dict[str, float]
            Interaction parameters (tau_ij) between component i and j.
        r_i_data : List[float] | Dict[str, float]
            relative van der Waals volume of component i
        q_i_data : List[float] | Dict[str, float]
            relative surface area of component i
        Z : int | float
            Model constant, Default is 10.
        calculation_mode : Literal['V1', 'V2']
            Mode of calculation. If 'V1', use the first version of the UNIQUAC model. If 'V2', use the second version.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".
        message : Optional[str]
            Message to be displayed. Default is None.
        res_format : Literal['dict', 'str', 'json']
            Format of the result. Default is 'dict'.

        Returns
        --------
        res : Dict[str, float | Dict]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.

        Notes
        -----
        input_values : Dict[str, float | Dict]
            Dictionary of input values used for the calculation as:
                - tau_ij: np.ndarray
                - tau_ij_comp: Dict[str, float]
                - r_i:  np.ndarray
                - r_i_comp: Dict[str, float]
                - q_i:  np.ndarray
                - q_i_comp: Dict[str, float]
        """
        try:
            # SECTION
            # Get the number of components
            components = self.components
            components_str = ", ".join(components)

            # comp no
            comp_num = self.comp_num

            # check Z
            if Z is None:
                Z = self.Z

            # check message
            if message is None:
                message = f"Calculate activity coefficients for {components_str} using UNIQUAC model"

            # SECTION
            # mole fraction (sorted by component id)
            xi = [mole_fraction[components[i]] for i in range(comp_num)]

            # SECTION
            # set the interaction parameter matrix (tau_ij) for the UNIQUAC model
            if isinstance(tau_ij_data, np.ndarray):
                # set
                tau_ij = tau_ij_data
                # to dict
                tau_ij_comp = self.to_dict_ij(
                    tau_ij_data, symbol_delimiter=symbol_delimiter)
            elif isinstance(tau_ij_data, TableMatrixData):
                # convert to numpy array and dict
                res_ = self.to_ij(data=tau_ij_data)
                # set
                tau_ij = res_[0]
                # to dict
                tau_ij_comp = res_[1]
            elif isinstance(tau_ij_data, dict):
                # convert dict to numpy array
                tau_ij = self.to_matrix_ij(
                    data=tau_ij_data, symbol_delimiter=symbol_delimiter)
                # to dict
                tau_ij_comp = tau_ij_data
            else:
                raise TypeError(
                    "tau_ij_data must be numpy array, dict or TableMatrixData")

            # SECTION
            # set relative van der Waals volume of component i (r_i) for the UNIQUAC model
            if isinstance(r_i_data, np.ndarray):
                # set
                r_i = r_i_data
                # to dict
                r_i_comp = self.to_dict_i(r_i_data)
            elif isinstance(r_i_data, List):
                # set
                r_i = np.array(r_i_data)
                # to dict
                r_i_comp = self.to_dict_i(r_i_data)
            elif isinstance(r_i_data, dict):
                # convert dict to numpy array
                r_i = self.to_i(data=r_i_data)
                # to dict
                r_i_comp = r_i_data
            else:
                raise TypeError("r_i_data must be numpy array, dict or List")

            # SECTION
            # set relative surface area of component i (q_i) for the UNIQUAC model
            if isinstance(q_i_data, np.ndarray):
                # set
                q_i = q_i_data
                # to dict
                q_i_comp = self.to_dict_i(q_i_data)
            elif isinstance(q_i_data, List):
                # set
                q_i = np.array(q_i_data)
                # to dict
                q_i_comp = self.to_dict_i(q_i_data)
            elif isinstance(q_i_data, dict):
                # convert dict to numpy array
                q_i = self.to_i(data=q_i_data)
                # to dict
                q_i_comp = q_i_data
            else:
                raise TypeError("q_i_data must be numpy array, dict or List")

            # SECTION
            # Calculate activity coefficients using the UNIQUAC model
            if calculation_mode == 'V1':
                AcCo_i = self.CalAcCo_V1(
                    xi=xi, tau_ij=tau_ij, r_i=r_i, q_i=q_i, Z=Z)
            else:
                raise ValueError("calculation_mode not supported!")

            # convert to float
            AcCo_i = [float(i) for i in AcCo_i]

            # SECTION
            # init the activity coefficients
            AcCo_i_comp = {components[i]: AcCo_i[i] for i in range(comp_num)}

            # SECTION: prepare result
            # input values
            other_values = {
                "AcCo_i_comp": AcCo_i_comp,
                'tau_ij': tau_ij,
                'tau_ij_comp': tau_ij_comp,
                'r_i': r_i,
                'r_i_comp': r_i_comp,
                'q_i': q_i,
                'q_i_comp': q_i_comp,
                'Z': Z,
                'calculation_mode': calculation_mode,
            }

            # res
            res = {
                'property_name': 'activity coefficients',
                'components': components,
                'mole_fraction': xi,
                'value': AcCo_i,
                'unit': 1,
                'symbol': "AcCo_i",
                'message': message,
            }

            # NOTE: check res_format
            if res_format == 'dict':
                # return as dict
                return res, other_values
            elif res_format == 'json':
                # return as json string
                res = json.dumps(res, indent=4)
                other_values = json.dumps(other_values, indent=4)
                return res, other_values
            elif res_format == 'str':
                # return as string
                res = str(res)
                other_values = str(other_values)
                return res, other_values
            else:
                raise ValueError("res_format not supported!")

        except Exception as e:
            raise Exception(
                f"Error in calculate_activity_coefficients: {str(e)}")

    def CalAcCo_V1(self, xi: List[float], tau_ij: np.ndarray, r_i: np.ndarray, q_i: np.ndarray, Z: int | float) -> List[float]:
        '''
        Calculate activity coefficient (AcCo) using UNIQUAC model.

        Parameters
        -----------
        xi: List[float]
            mole fraction of each component in the mixture
        tau_ij: np.ndarray
            interaction parameters (tau_ij) between component i and j
        r_i: np.ndarray
            relative van der Waals volume of component i
        q_i: np.ndarray
            relative surface area of component i
        Z: int | float
            model constant, default is 10

        Returns
        -------
        AcCoi: List[float]
            activity coefficient for each component

        Notes
        -----
        This function is used to calculate the activity coefficient for each component.

        1. tau_ij: temperature dependent parameters (tau_ij[i,i]=tau_ij[j,j]=1) calculated at temperature T
        '''
        try:
            # component no
            comp_num = self.comp_num
            # check
            if len(xi) != comp_num:
                raise ValueError(
                    f"xi length {len(xi)} does not match component number {comp_num}")

            # SECTION: calculate activity coefficients
            # ∑r[i]x[i]
            Sigma_rx = r_i@xi

            # ∑q[i]x[i]
            Sigma_qx = q_i@xi

            # volume fraction/mole fraction of component i
            Vi = r_i/Sigma_rx

            # surface area/mole fraction of component i
            Fi = q_i/Sigma_qx

            # volume fraction
            phi_i = (r_i*xi)/Sigma_rx

            # surface area fraction
            teta_i = (q_i*xi)/Sigma_qx

            # S
            Si = np.zeros(comp_num)
            for i in range(comp_num):
                Si[i] = np.dot(teta_i, tau_ij[:, i])

            # combinatorial part of the activity coefficient
            gamma_comb_ij = (np.log(phi_i/xi) + 1 - (phi_i/xi) -
                             (Z/2)*q_i*(np.log(phi_i/teta_i)+1-(phi_i/teta_i)))

            # residual part of the activity coefficient
            gamma_res_ij = np.zeros(comp_num)
            for i in range(comp_num):
                gamma_res_ij[i] = q_i[i] * \
                    (1-np.log(Si[i])-np.dot(tau_ij[i, :], (teta_i/Si)))

            # activity coefficient
            AcCoi = np.exp(gamma_comb_ij+gamma_res_ij)
            # to list
            AcCoi = AcCoi.tolist()

            # res
            return AcCoi
        except Exception as e:
            raise Exception(f"Error in CalAcCo_V1: {str(e)}")

    def excess_gibbs_free_energy(self, mole_fraction: Dict[str, float], tau_ij: np.ndarray,
                                 r_i: np.ndarray, q_i: np.ndarray, Z: Optional[int | float] = None,
                                 message: Optional[str] = None, res_format: Literal['str', 'json', 'dict'] = 'dict') -> Dict[str, float | Dict] | str:
        """
        Calculate excess Gibbs energy (G^E/RT) for a multicomponent mixture using the UNIQUAC model.

        Parameters
        -----------
        mole_fraction : dict
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        tau_ij : np.ndarray
            Matrix of tau parameters where tau[i][j] is tau_ij between component i and j.
        r_i : np.ndarray
            Array of relative van der Waals volumes of each component.
        q_i : np.ndarray
            Array of relative surface areas of each component.
        Z : int | float, optional
            Model constant, default is 10.
        message : str, optional
            Message to be printed, default is None.
        res_format : str, optional
            Format of the result, default is 'dict'. Can be 'str' or 'dict'.

        Returns
        --------
        res : dict
            Dictionary containing the excess Gibbs energy and other information.
        """
        try:
            # check
            if not isinstance(mole_fraction, dict):
                raise TypeError("mole_fraction must be a dictionary")

            # check Z
            if Z is None:
                Z = self.Z

            # components
            components = self.components
            components_str = ', '.join(components)

            # comp no
            comp_num = self.comp_num

            # set message
            message = f'Excess Gibbs Free Energy for {components_str}' if message is None else message

            # mole fraction
            xi = [mole_fraction[components[i]] for i in range(len(components))]

            # Normalize mole fractions to ensure they sum to 1
            x = xi / np.sum(xi)
            x = np.array(x)

            # NOTE: check all input
            if len(x) != comp_num:
                raise ValueError(
                    f"mole_fraction length {len(x)} does not match component number {comp_num}")

            if len(r_i) != comp_num:
                raise ValueError(
                    f"r_i length {len(r_i)} does not match component number {comp_num}")

            if len(q_i) != comp_num:
                raise ValueError(
                    f"q_i length {len(q_i)} does not match component number {comp_num}")

            # SECTION
            # excess gibbs free energy
            gE_RT = 0

            # Volume and surface area fractions
            phi = (r_i * x) / np.sum(r_i * x)
            theta = (q_i * x) / np.sum(q_i * x)

            # Combinatorial term
            term1 = np.sum(x * np.log(phi / x))
            term2 = np.sum(q_i * x * np.log(theta / phi))
            gE_C = term1 + (Z / 2) * term2

            # Residual term
            gE_R = 0.0
            for i in range(len(x)):
                # sum_j (theta_j * tau_ji)
                inner_sum = np.sum(theta * tau_ij[:, i])
                gE_R -= q_i[i] * x[i] * np.log(inner_sum)

            gE_RT = gE_C + gE_R

            # SECTION: set result format
            res = {
                "property_name": "Excess Molar Gibbs Free Energy",
                "components": components,
                "mole_fraction": xi,
                "value": float(gE_RT),
                "unit": 1,
                "symbol": "ExMoGiFrEn",
                'message': message
            }

            if res_format == 'dict':
                return res
            elif res_format == 'json' or res_format == 'str':
                return json.dumps(res, indent=4)
            else:
                raise ValueError("res_format must be 'dict', 'json' or 'str'")
        except Exception as e:
            raise Exception(f"Error in excess_gibbs_free_energy: {str(e)}")
