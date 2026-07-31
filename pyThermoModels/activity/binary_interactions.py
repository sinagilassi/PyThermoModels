# import libs
import logging
import numpy as np
from math import pow, exp, log
from typing import List, Dict, Tuple, Literal, Optional, Any
import pycuc
from pyThermoDB import (
    TableMatrixData,
    TableData,
    TableEquation,
    TableMatrixEquation
)
# locals


class BinaryInteractions:
    def __init__(
        self,
        components: List[str],
        datasource: Dict = {},
        equationsource: Dict = {},
        **kwargs
    ):
        '''
        Initialize shared binary-interaction model state.

        Parameters
        ----------
        datasource: Dict
            Data source for the model.
        equationsource: Dict
            Equation source for the model.
        components: List[str]
            List of component names in the mixture.
        **kwargs: dict
            Additional keyword arguments.
        '''
        # Check datasource
        if not isinstance(datasource, dict):
            raise TypeError("datasource must be a dict")

        # Check equationsource
        if not isinstance(equationsource, dict):
            raise TypeError("equationsource must be a dict")

        # Check if components is a list
        if not isinstance(components, list):
            raise TypeError("components must be a list")

        # Assign the parameters to instance variables
        self.datasource = datasource
        self.equationsource = equationsource

        # components
        self.components = [component.strip() for component in components]

        # Get the number of components
        self.comp_num = len(self.components)

        # component index
        self.comp_idx = {
            self.components[i]: i for i in range(self.comp_num)
        }

    def cal_tau_ij_X(
        self,
        temperature: float,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        d_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
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

            tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T

        2. Interaction energy parameter symbol is `X` for TableMatrixData as:

        - `X_{component_i}-{component_j}`
        - `X | {component_i} | {component_j}`.

        2. All parameters including a_ij, b_ij, c_ij, d_ij must be in the same format (numpy array, dict or TableMatrixData).
        """
        try:
            # SECTION: check
            if (not isinstance(a_ij, np.ndarray) and
                not isinstance(a_ij, dict) and
                    not isinstance(a_ij, TableMatrixData)):
                raise TypeError(
                    "a_ij must be numpy array, dict or TableMatrixData")

            if (not isinstance(b_ij, np.ndarray) and
                not isinstance(b_ij, dict) and
                    not isinstance(b_ij, TableMatrixData)):
                raise TypeError(
                    "b_ij must be numpy array, dict or TableMatrixData")

            if (not isinstance(c_ij, np.ndarray) and
                not isinstance(c_ij, dict) and
                    not isinstance(c_ij, TableMatrixData)):
                raise TypeError(
                    "c_ij must be numpy array, dict or TableMatrixData")

            if (not isinstance(d_ij, np.ndarray) and
                not isinstance(d_ij, dict) and
                    not isinstance(d_ij, TableMatrixData)):
                raise TypeError(
                    "d_ij must be numpy array, dict or TableMatrixData")

            # Get the number of components
            comp_num = self.comp_num

            # components
            components = self.components

            # Initialize tau_ij matrix
            tau_ij = np.zeros((comp_num, comp_num))

            # tau_ij components
            tau_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # SECTION: Calculate tau_ij values
            if (
                isinstance(a_ij, np.ndarray) and
                isinstance(b_ij, np.ndarray) and
                isinstance(c_ij, np.ndarray) and
                    isinstance(d_ij, np.ndarray)
            ):
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
            elif (
                isinstance(a_ij, dict) and
                isinstance(b_ij, dict) and
                isinstance(c_ij, dict) and
                isinstance(d_ij, dict)
            ):
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
            elif (
                isinstance(a_ij, TableMatrixData) and
                isinstance(b_ij, TableMatrixData) and
                isinstance(c_ij, TableMatrixData) and
                isinstance(d_ij, TableMatrixData)
            ):
                # convert to numpy array and dict
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"X_{components[i]}_{components[j]}"
                        # dict
                        key_comp = f"{components[i]}{symbol_delimiter_set}{components[j]}"

                        # TODO: extract val
                        # a_ij
                        a_ij_ = a_ij.ij(key_)
                        if a_ij_ is not None and a_ij_["value"] is not None:
                            a_ij_val = float(a_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for a_ij: {a_ij_} for key: {key_}")

                        # b_ij
                        b_ij_ = b_ij.ij(key_)
                        if b_ij_ is not None and b_ij_["value"] is not None:
                            b_ij_val = float(b_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for b_ij: {b_ij_} for key: {key_}")

                        # c_ij
                        c_ij_ = c_ij.ij(key_)
                        if c_ij_ is not None and c_ij_["value"] is not None:
                            c_ij_val = float(c_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for c_ij: {c_ij_} for key: {key_}")

                        # d_ij
                        d_ij_ = d_ij.ij(key_)
                        if d_ij_ is not None and d_ij_["value"] is not None:
                            d_ij_val = float(d_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for d_ij: {d_ij_} for key: {key_}")

                        # val
                        val_ = a_ij_val + b_ij_val / temperature + c_ij_val * \
                            log(temperature) + d_ij_val * temperature

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
