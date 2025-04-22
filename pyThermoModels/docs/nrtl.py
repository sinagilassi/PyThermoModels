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


class NRTL:
    """
    The NRTL (`Non-Random Two-Liquid`) model - a thermodynamic framework used to describe the behavior of mixtures,
    particularly in the context of phase equilibria and activity coefficients.

    The NRTL model relies on several key parameters to describe the interactions between components in a mixture. These parameters are:
    - Δg_ij (interaction energy parameter): represents the interaction energy between two molecules [J/mol].
    - α_ij (non-randomness parameter): represents the non-randomness of the mixture [dimensionless].
    - τ_ij (binary interaction parameter): represents the interaction energy between two molecules of different components [dimensionless].

    Universal gas constant (R) is defined as 8.314 J/mol/K.
    """

    # universal gas constant [J/mol/K]
    R_CONST = 8.314

    def __init__(self,
                 components: List[str],
                 datasource: Dict = {},
                 equationsource: Dict = {}):
        '''
        Initialize the NRTL (`Non-Random Two-Liquid`) model used to calculate activity coefficients in liquid mixtures.

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
        The NRTL model needs the following parameters:
        - datasource: Data source for the model
        - equationsource: Equation source for the model
        - components: List of component names in the mixture

        The component names define the order of the parameters in the model. The first component in the list is component 1, the second is component 2, and so on.
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

    def __repr__(self) -> str:
        return """
        The NRTL (`Non-Random Two-Liquid`) model - a thermodynamic framework used to describe the behavior of mixtures,
        particularly in the context of phase equilibria and activity coefficients.

        The NRTL model relies on several key parameters to describe the interactions between components in a mixture. These parameters are:
        - Δg_ij (interaction energy parameter): represents the interaction energy between two molecules [J/mol].
        - α_ij (non-randomness parameter): represents the non-randomness of the mixture [dimensionless].
        - τ_ij (binary interaction parameter): represents the interaction energy between two molecules of different components [dimensionless].

        Universal gas constant (R) is defined as 8.314 J/mol/K.
        """

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

        Notes
        -----
        1.
        """
        try:
            # NOTE: check dg_ij is TableMatrixData
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

    def to_matrix_ij(self, data: Dict[str, float], symbol_delimiter: Literal["|", "_"] = "|") -> np.ndarray:
        """
        Convert to matrix (mat_ij) according to `the component id`.

        Parameters
        ----------
        data : Dict[str, float]
            Dictionary of parameters where keys are component pairs and values are their respective values.
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
            for i in range(comp_num):
                for j in range(comp_num):
                    # val
                    val = data[f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"]

                    # find the component id
                    comp_id_i = self.comp_idx[self.components[i]]
                    comp_id_j = self.comp_idx[self.components[j]]

                    # to matrix
                    mat_ij[comp_id_i, comp_id_j] = val

            # res
            return mat_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def cal_dg_ij_M1(self, temperature: float,
                     a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                     b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                     c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                     symbol_delimiter: Literal["|", "_"] = "|") -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate interaction energy parameter `dg_ij` matrix dependent of temperature.

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
        dg_ij : np.ndarray
            Interaction energy parameter `dg_ij` matrix for NRTL model.
        dg_ij_comp : dict
            Dictionary of interaction energy parameters where keys are component pairs and values are their respective dg_ij values.

        Notes
        -----
        1. The interaction energy parameter matrix is calculated using the formula:

            `dg_ij = a_ij + b_ij * T + c_ij * T^2`

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

            # Initialize dg_ij matrix
            dg_ij = np.zeros((comp_num, comp_num))

            # dg_ij components
            dg_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # SECTION: calculate dg_ij values
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
                            dg_ij[i, j] = val_
                            # set by name
                            dg_ij_comp[key_] = val_
                        else:
                            # set
                            dg_ij[i, j] = 0
                            # set by name
                            dg_ij_comp[key_] = 0

            # SECTION: if dg_ij is dict
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
                            dg_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            dg_ij_comp[key_] = val_
                        else:
                            # set
                            dg_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            dg_ij_comp[key_] = 0
            # SECTION: if dg_ij is TableMatrixData
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
                            dg_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            dg_ij_comp[key_comp] = val_
                        else:
                            # set
                            dg_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            dg_ij_comp[key_comp] = 0
            else:
                raise TypeError(
                    "a_ij, b_ij and c_ij must be numpy array or dict")
            # res
            return dg_ij, dg_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_dg_ij_M1: {str(e)}")

    def cal_tau_ij_M1(self, temperature: float, dg_ij: np.ndarray | Dict[str, float] | TableMatrixData,
                      dg_ij_symbol: Literal['dg', 'dg_ij'] = 'dg', R_CONST: float = 8.314,
                      symbol_delimiter: Literal["|", "_"] = "|") -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate interaction parameters `tau_ij` matrix for NRTL model.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin [K].
        dg_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction energy parameter [J/mol] matrix where dg_ij[i][j] between component i and j.
        dg_ij_symbol : str
            Interaction energy parameter symbol. Default is 'dg'.
        R_CONST : float
            Univeral gas constant [J/mol/K], default R_CONST = 8.314
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
        1. The tau_ij matrix is calculated using the formula:

            `tau_ij = dg_ij / (R * T)`

            where R is the universal gas constant [J/mol/K] and T is the temperature [K].

        2. Interaction energy parameter symbol is `dg` for TableMatrixData as:

        - `dg_{component_i}-{component_j}`
        - `dg | {component_i} | {component_j}`.
        """
        try:
            # check
            if not isinstance(dg_ij, np.ndarray) and not isinstance(dg_ij, dict) and not isinstance(dg_ij, TableMatrixData):
                raise TypeError(
                    "dg_ij must be numpy array, dict or TableMatrixData")

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
            # SECTION: if dg_ij is numpy array
            if isinstance(dg_ij, np.ndarray):
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # check
                        if i != j:
                            # val
                            val_ = dg_ij[i, j] / (R_CONST * temperature)
                            # set
                            tau_ij[i, j] = val_
                            # set by name
                            tau_ij_comp[key_] = val_
                        else:
                            # set
                            tau_ij[i, j] = 0
                            # set by name
                            tau_ij_comp[key_] = 0
            # SECTION: if dg_ij is dict
            elif isinstance(dg_ij, dict):
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
                            val_ = dg_ij[key_] / (R_CONST * temperature)
                            # set
                            tau_ij[comp_id_i, comp_id_j] = val_
                            # set by name
                            tau_ij_comp[key_] = val_
                        else:
                            # set
                            tau_ij[comp_id_i, comp_id_j] = 0
                            # set by name
                            tau_ij_comp[key_] = 0
            # SECTION: if dg_ij is TableMatrixData
            elif isinstance(dg_ij, TableMatrixData):
                # convert to numpy array and dict
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{dg_ij_symbol}_{components[i]}_{components[j]}"
                        # dict
                        key_comp = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # val
                        val_ = dg_ij.ij(key_)

                        # component id
                        comp_id_i = self.comp_idx[components[i]]
                        comp_id_j = self.comp_idx[components[j]]

                        # set
                        if i != j:
                            # val
                            val_ = val_ / (R_CONST * temperature)
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
                raise TypeError("dg_ij must be numpy array or dict")

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
        Calculate interaction parameters `tau_ij` matrix for NRTL model.

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
            interaction parameters `tau_ij` matrix for NRTL model.
        tau_ij_comp : dict
            Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.

        Notes
        -----
        1. The extended Antoine equation format is used to calculate the interaction parameters using the following formula:

            `tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T`

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
            # SECTION: if dg_ij is dict
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
            # SECTION: if dg_ij is TableMatrixData
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

    def cal_G_ij(self, tau_ij: np.ndarray, alpha_ij: np.ndarray, symbol_delimiter: Literal["|", "_"] = "|"):
        """
        Calculate non-randomness parameters `G_ij` matrix for NRTL model according to `the component id`.

        Parameters
        ----------
        tau_ij : np.ndarray
            Interaction parameters `tau_ij` matrix for NRTL model.
        alpha_ij : np.ndarray
            Non-randomness parameters [dimensionless] matrix where alpha_ij[i][j] between component i and j.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        G_ij : np.ndarray
            Non-randomness parameters `G_ij` matrix for NRTL model.
        G_ij_comp : dict
            Dictionary of non-randomness parameters where keys are component pairs and values are their respective G_ij values.

        Notes
        -----
        The G_ij matrix is calculated using the formula:

        `G_ij = exp(-alpha_ij * tau_ij)`

        where alpha_ij is the non-randomness parameter and tau_ij is the interaction parameter.
        """
        try:
            # check
            if not isinstance(tau_ij, np.ndarray):
                raise TypeError("tau_ij must be numpy array")

            if not isinstance(alpha_ij, np.ndarray):
                raise TypeError("alpha_ij must be numpy array")

            # Get the number of components
            comp_num = self.comp_num

            # components
            components = self.components

            # Initialize Gij matrix
            G_ij = np.ones((comp_num, comp_num))
            # dict
            G_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # Calculate Gij values
            for i in range(comp_num):
                for j in range(comp_num):
                    # key
                    key_ = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                    # check
                    if i != j:
                        # val
                        val_ = exp(-1 * alpha_ij[i, j] * tau_ij[i, j])
                        # set
                        G_ij[i, j] = val_
                        # set by name
                        G_ij_comp[key_] = val_
                    else:
                        # set
                        G_ij[i, j] = 1
                        # set by name
                        G_ij_comp[key_] = 1

            # res
            return G_ij, G_ij_comp
        except Exception as e:
            raise Exception(f"Error in cal_Gij: {str(e)}")

    def cal(self,
            model_input: Dict,
            calculation_mode: Literal['V1', 'V2'] = 'V1',
            symbol_delimiter: Literal["|", "_"] = "|",
            message: Optional[str] = None,
            res_format: Literal['dict', 'str', 'json'] = 'dict',
            **kwargs):
        '''
        Calculate activity coefficients for a multi-component mixture using the NRTL model.

        Parameters
        ----------
        model_input: Dict
            Dictionary of model input values where keys are parameter names and values are their respective values.
                - `mole_fraction`: Dict[str, float]
                    Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
                - `tau_ij`: TableMatrixData | np.ndarray | Dict[str, float]
                    Interaction parameters (tau_ij) between component i and j.
                - `alpha_ij`: TableMatrixData | np.ndarray | Dict[str, float]
                    Non-randomness parameters (alpha_ij) between component i and j.
        calculation_mode: Literal['V1', 'V2']
            Mode of calculation. If 'V1', use the first version of the NRTL model. If 'V2', use the second version.
        symbol_delimiter: Literal["|", "_"]
            Delimiter for the component id. Default is "|".
        message: Optional[str]
            Message to be displayed. Default is None.
        res_format: Literal['dict', 'str', 'json']
            Format of the result. Default is 'dict'.
        **kwargs: Optional
            Additional keyword arguments for the calculation.

        Returns
        -------
        res: Dict[str, float | Dict]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.

        Examples
        --------
        >>> model_input = {
        ...     'mole_fraction': {'A': 0.5, 'B': 0.5},
        ...     'tau_ij': np.array([[0, 1], [1, 0]]),
        ...     'alpha_ij': np.array([[0, 0.5], [0.5, 0]])
        ... }
        >>> calculation_mode = 'V1'
        >>> message = 'Calculating activity coefficients'
        >>> result = cal(model_input, calculation_mode, message)
        >>> print(result)
        '''
        try:
            # SECTION: check
            if not isinstance(model_input, dict):
                raise TypeError("model_input must be dict")

            # SECTION: check keys
            required_keys = [
                'mole_fraction', 'tau_ij', 'alpha_ij']
            for key in required_keys:
                if key not in model_input:
                    raise KeyError(f"{key} is required in model_input")

            # SECTION: get values
            mole_fraction = model_input['mole_fraction']
            tau_ij_data = model_input['tau_ij']
            alpha_ij_data = model_input['alpha_ij']

            # SECTION: calculate activity coefficients
            return self.__calculate_activity_coefficients(
                mole_fraction=mole_fraction,
                tau_ij_data=tau_ij_data,
                alpha_ij_data=alpha_ij_data,
                calculation_mode=calculation_mode,
                symbol_delimiter=symbol_delimiter,
                message=message,
                res_format=res_format
            )
        except Exception as e:
            raise Exception(f"Error in launch_calculation: {str(e)}")

    def __calculate_activity_coefficients(self,
                                          mole_fraction: Dict[str, float],
                                          tau_ij_data: TableMatrixData | np.ndarray | Dict[str, float],
                                          alpha_ij_data: TableMatrixData | np.ndarray | Dict[str, float],
                                          calculation_mode: Literal['V1', 'V2'],
                                          symbol_delimiter: Literal["|", "_"],
                                          message: Optional[str],
                                          res_format: Literal['dict', 'str', 'json']) -> Tuple[Dict[str, str | float | Dict], Dict[str, str | float | Dict]] | str:
        """
        Calculate activity coefficients for a multicomponent mixture using the NRTL model.

        Parameters
        -----------
        mole_fraction : Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        temperature : float
            Temperature in Kelvin.
        tau_ij_comp : TableMatrixData | np.ndarray
            Interaction parameters (tau_ij) between component i and j.
        alpha_ij_comp : TableMatrixData | np.ndarray
            Non-randomness parameters (alpha_ij) between component i and j.
        calculation_mode : Literal['V1', 'V2']
            Mode of calculation. If 'V1', use the first version of the NRTL model. If 'V2', use the second version.
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
                - alpha_ij: np.ndarray
                - alpha_ij_comp: Dict[str, float]
                - G_ij: np.ndarray
                - G_ij_comp: Dict[str, float]
                - calculation_mode: str
        """
        try:
            # SECTION
            # Get the number of components
            components = self.components
            components_str = ", ".join(components)

            # comp no
            comp_num = self.comp_num

            # SECTION
            # mole fraction (sorted by component id)
            xi = [mole_fraction[components[i]] for i in range(comp_num)]

            # check message
            if message is None:
                message = f"Calculate activity coefficients for {components_str} using NRTL model"

            # SECTION
            # set the interaction parameter matrix (tau_ij) for the NRTL model
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
            # set the non-randomness parameter matrix (alpha_ij) for the NRTL model
            if isinstance(alpha_ij_data, np.ndarray):
                # set
                alpha_ij = alpha_ij_data
                # to dict
                alpha_ij_comp = self.to_dict_ij(
                    alpha_ij_data, symbol_delimiter=symbol_delimiter)
            elif isinstance(alpha_ij_data, TableMatrixData):
                # convert to numpy array and dict
                res_ = self.to_ij(data=alpha_ij_data)
                # set
                alpha_ij = res_[0]
                # to dict
                alpha_ij_comp = res_[1]
            elif isinstance(alpha_ij_data, dict):
                # convert dict to numpy array
                alpha_ij = self.to_matrix_ij(
                    data=alpha_ij_data, symbol_delimiter=symbol_delimiter)
                # to dict
                alpha_ij_comp = alpha_ij_data
            else:
                raise TypeError(
                    "alpha_ij_data must be numpy array, dict or TableMatrixData")

            # SECTION
            # set G_ij matrix for NRTL model
            G_ij, G_ij_comp = self.cal_G_ij(
                tau_ij=tau_ij, alpha_ij=alpha_ij, symbol_delimiter=symbol_delimiter)

            # SECTION
            # Calculate activity coefficients using the NRTL model
            if calculation_mode == 'V1':
                AcCo_i = self.CalAcCo_V1(xi=xi, tau_ij=tau_ij, G_ij=G_ij)
            elif calculation_mode == 'V2':
                AcCo_i = self.CalAcCo_V2(xi=xi, tau_ij=tau_ij, G_ij=G_ij)
            else:
                raise ValueError("calculation_mode not supported!")

            # set the activity coefficients float
            AcCo_i = [float(AcCo_i[i]) for i in range(comp_num)]

            # SECTION
            # init the activity coefficients
            AcCo_i_comp = {components[i]: float(
                AcCo_i[i]) for i in range(comp_num)}

            # SECTION: prepare result
            # input values
            other_values = {
                "AcCo_i_comp": AcCo_i_comp,
                'tau_ij': tau_ij,
                'tau_ij_comp': tau_ij_comp,
                'alpha_ij': alpha_ij,
                'alpha_ij_comp': alpha_ij_comp,
                'G_ij': G_ij,
                'G_ij_comp': G_ij_comp,
                'calculation_mode': calculation_mode,
            }

            # res
            res = {
                'property_name': 'activity coefficients',
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

    def CalAcCo_V1(self, xi: List[float], tau_ij: np.ndarray, G_ij: np.ndarray) -> np.ndarray:
        '''
        Calculate activity coefficient (AcCo) using Non-random two-liquid (NRTL) model.

        Parameters
        -----------
        xi: List[float]
            mole fraction of each component in the mixture
        tau_ij: np.ndarray
            interaction parameters (tau_ij) between component i and j
        G_ij: np.ndarray
            non-randomness parameters (G_ij) between component i and j

        Returns
        -------
        AcCoi: np.ndarray
            activity coefficient for each component

        Notes
        -----
        This function is used to calculate the activity coefficient for each component.

        1. tau_ij: temperature dependent parameters (ta[i,i]=ta[j,j]=0) calculated at temperature T

        '''
        try:
            # component no
            comp_num = self.comp_num
            # check
            if len(xi) != comp_num:
                raise ValueError(
                    f"xi length {len(xi)} does not match component number {comp_num}")

            # activity coefficient
            AcCoi = np.zeros(comp_num)

            # activity coefficient
            C0 = np.zeros((comp_num, comp_num))

            for i in range(comp_num):
                _c0 = 0
                for j in range(comp_num):
                    _c0 = tau_ij[j, i]*G_ij[j, i]*xi[j] + _c0

                _c1 = 0
                for k in range(comp_num):
                    _c1 = G_ij[k, i]*xi[k] + _c1

                for j in range(comp_num):
                    _c2 = xi[j]*G_ij[i, j]

                    _c3 = 0
                    for k in range(comp_num):
                        _c3 = G_ij[k, j]*xi[k] + _c3

                    _c4 = 0
                    for n in range(comp_num):
                        _c4 = xi[n]*tau_ij[n, j]*G_ij[n, j] + _c4

                    _c5 = tau_ij[i, j] - (_c4/_c3)

                    # set
                    C0[i, j] = (_c2/_c3)*_c5

                _c6 = (_c0/_c1) + np.sum(C0[i, :])
                AcCoi[i] = exp(_c6)

            # res
            return AcCoi
        except Exception as e:
            raise Exception(f"Error in CalAcCo_V1: {str(e)}")

    def CalAcCo_V2(self, xi: list[float], tau_ij: np.ndarray, G_ij: np.ndarray) -> np.ndarray:
        """
        Calculate activity coefficients for a multicomponent mixture using the NRTL model.

        Parameters:
        -----------
        xi : list[float]
            Mole fractions of each component in the mixture.
        tau_ij : np.ndarray
            Binary interaction parameters (tau_ij) between component i and j.
        G_ij : np.ndarray
            Non-randomness parameters (G_ij) between component i and j.

        Returns:
        --------
        AcCoi : np.ndarray
            activity coefficients for each component, tauij matrix, Gij matrix
        """
        try:
            # component no
            comp_num = self.comp_num
            # check
            if len(xi) != comp_num:
                raise ValueError(
                    f"xi length {len(xi)} does not match component number {comp_num}")

            # activity coefficient
            ln_gamma = np.zeros(comp_num)

            for i in range(comp_num):
                # Calculate the first term: Σj τ_ji*G_ji*x_j / Σk G_ki*x_k
                denom_i = np.sum(G_ij[:, i] * xi)
                numer_i = np.sum(tau_ij[:, i] * G_ij[:, i] * xi)
                first_term = numer_i / denom_i

                # Calculate the second term (the summation)
                second_term = 0
                for j in range(comp_num):
                    denom_j = np.sum(G_ij[:, j] * xi)
                    numer_j = np.sum(xi * tau_ij[:, j] * G_ij[:, j])
                    second_term += (xi[j] * G_ij[i, j] / denom_j) * \
                        (tau_ij[i, j] - numer_j / denom_j)

                # Combine the terms to get ln(gamma_i)
                ln_gamma[i] = first_term + second_term

            # Calculate the activity coefficients (gamma_i)
            AcCoi = np.zeros(comp_num)
            for i in range(comp_num):
                AcCoi[i] = exp(ln_gamma[i])

            # res
            return AcCoi
        except Exception as e:
            raise Exception(f"Error in CalAcCoV2: {str(e)}")

    def excess_gibbs_free_energy(self, mole_fraction: Dict[str, float], G_ij: np.ndarray, tau_ij: np.ndarray,
                                 message: Optional[str] = None, res_format: Literal['str', 'json', 'dict'] = 'dict') -> Dict[str, float | Dict] | str:
        """
        Calculate excess Gibbs energy for a multicomponent mixture using the NRTL model.

        Parameters
        -----------
        mole_fraction : dict
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        G_ij : np.ndarray
            Matrix of G parameters where G[i][j] is G_ij
        tau_ij : np.ndarray
            Matrix of tau parameters where tau[i][j] is tau_ij
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
            # components
            components = self.components
            components_str = ', '.join(components)

            # mole fraction
            xi = [mole_fraction[components[i]] for i in range(len(components))]

            # set message
            message = f'Excess Gibbs Free Energy for {components_str}' if message is None else message

            # Normalize mole fractions to ensure they sum to 1
            x = xi / np.sum(xi)

            n = len(x)  # Number of components
            gE_RT = 0

            for i in range(n):
                # Calculate denominator sum (Σj G_ji*x_j)
                denom = np.sum(G_ij[:, i] * x)

                # Calculate numerator sum (Σj τ_ji*G_ji*x_j)
                numer = np.sum(tau_ij[:, i] * G_ij[:, i] * x)

                # Add to excess Gibbs energy
                gE_RT += xi[i] * numer / denom

            # SECTION: set result format
            res = {
                "property_name": "Excess Molar Gibbs Free Energy",
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
