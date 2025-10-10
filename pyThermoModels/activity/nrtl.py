# import libs
import logging
import numpy as np
import yaml
import json
from math import pow, exp, log
from typing import List, Dict, Tuple, Literal, Optional, Any
import pycuc
from pyThermoDB import (
    TableMatrixData,
    TableData,
    TableEquation,
    TableMatrixEquation
)
from pythermodb_settings.utils import create_mixture_id
# local
from ..utils import add_attributes
from ..plugin import ACTIVITY_MODELS

# NOTE: logger
logger = logging.getLogger(__name__)


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

    # NOTE: variables based on component id
    __tau_ij = None
    __tau_ij_comp = None
    __dg_ij = None
    __dg_ij_comp = None
    __alpha_ij = None
    __alpha_ij_comp = None
    __G_ij = None
    __G_ij_comp = None
    __mole_fraction = None
    __xi = None

    # mixture id
    _mixture_id: Optional[str] = None

    def __init__(
        self,
        components: List[str],
        datasource: Dict = {},
        equationsource: Dict = {},
        **kwargs
    ):
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
        **kwargs: dict
            Additional keyword arguments
            - mixture_id: str, optional
                Mixture ID for the components. If not provided, it will be generated automatically.

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

        # components
        self.components = [components.strip() for components in components]

        # SECTION
        # Get the number of components
        self.comp_num = len(components)
        # idx
        self.comp_idx = {components[i]: i for i in range(self.comp_num)}

        # SECTION: kwargs
        self._mixture_id = kwargs.get('mixture_id', None)

    def __repr__(self) -> str:
        model_ = """
        The NRTL (`Non-Random Two-Liquid`) model - a thermodynamic framework used to describe the behavior of mixtures,
        particularly in the context of phase equilibria and activity coefficients.

        The NRTL model relies on several key parameters to describe the interactions between components in a mixture. These parameters are:
        - Δg_ij (interaction energy parameter): represents the interaction energy between two molecules [J/mol].
        - α_ij (non-randomness parameter): represents the non-randomness of the mixture [dimensionless].
        - τ_ij (binary interaction parameter): represents the interaction energy between two molecules of different components [dimensionless].

        Universal gas constant (R) is defined as 8.314 J/mol/K.
        """
        return model_

    def parse_model_inputs(self, model_inputs: str) -> Dict[str, Any]:
        '''
        Convert model inputs from string to dictionary format.

        Parameters
        ----------
        model_inputs: str
            Model inputs in string format, such as:
            - { mole_fraction: { ethanol: 0.4, butyl-methyl-ether: 0.6 }, temperature: [323.15, 'K'], tau_ij: [[],[]], alpha_ij: [[],[]] }

        Returns
        -------
        model_input_parsed: dict
            Parsed model inputs in dictionary format.
        '''
        try:
            # check if model_inputs is None or model_inputs == 'None'
            if model_inputs is None or model_inputs == 'None':
                raise Exception('Model inputs are not provided!')

            # strip
            model_inputs = model_inputs.strip()
            model_input_parsed = yaml.safe_load(model_inputs)

            return model_input_parsed
        except Exception as e:
            raise Exception("Parsing model inputs failed!, ", e)

    def to_ij(
        self,
        data: TableMatrixData,
        prop_symbol: str,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
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
                    # key
                    key_ = f"{prop_symbol}_{self.components[i]}-{self.components[j]}"
                    # val
                    val = data.ij(key_)

                    # to matrix
                    # ? val["value"] or 0.0 in case of None
                    if val is not None and val["value"] is not None:
                        mat_ij[i, j] = float(val["value"])
                    else:
                        raise ValueError(
                            f"Invalid value for {prop_symbol}: {val} for key: {prop_symbol}_{self.components[i]}-{self.components[j]}")

                    # to dict
                    key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"
                    # to dict
                    # ? val["value"] or 0.0 in case of None
                    if val is not None and val["value"] is not None:
                        dict_ij[key_] = float(val["value"])
                    else:
                        raise ValueError(
                            f"Invalid value for {prop_symbol}: {val} for key: {prop_symbol}_{self.components[i]}-{self.components[j]}")

            # res
            return mat_ij, dict_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_dict_ij(
        self,
        data: np.ndarray,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Dict[str, float]:
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

    def to_matrix_ij(
        self,
        data: Dict[str, float],
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> np.ndarray:
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
                    # key
                    key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"
                    # val
                    val = data[key_]

                    # find the component id
                    comp_id_i = self.comp_idx[self.components[i]]
                    comp_id_j = self.comp_idx[self.components[j]]

                    # to matrix
                    mat_ij[comp_id_i, comp_id_j] = val

            # res
            return mat_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def cal_dg_ij_M1(
        self, temperature: float,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
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
            if (
                isinstance(a_ij, np.ndarray) and
                isinstance(b_ij, np.ndarray) and
                    isinstance(c_ij, np.ndarray)
            ):

                # looping through the components
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
            elif (
                isinstance(a_ij, dict) and
                isinstance(b_ij, dict) and
                isinstance(c_ij, dict)
            ):

                # looping through the components
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
            elif (
                isinstance(a_ij, TableMatrixData) and
                isinstance(b_ij, TableMatrixData) and
                isinstance(c_ij, TableMatrixData)
            ):
                # convert to numpy array and dict
                for i in range(comp_num):
                    for j in range(comp_num):
                        # key
                        key_ = f"{self.components[i]}_{self.components[j]}"
                        # dict
                        key_comp = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"

                        # TODO: extract val
                        # a_ij
                        a_ij_ = a_ij.ij(f"a_{key_}")
                        if a_ij_ is not None and a_ij_["value"] is not None:
                            a_ij_val = float(a_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for a_ij: {a_ij_} for key: {key_}")

                        # b_ij
                        b_ij_ = b_ij.ij(f"b_{key_}")
                        if b_ij_ is not None and b_ij_["value"] is not None:
                            b_ij_val = float(b_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for b_ij: {b_ij_} for key: {key_}")

                        # c_ij
                        c_ij_ = c_ij.ij(f"c_{key_}")
                        if c_ij_ is not None and c_ij_["value"] is not None:
                            c_ij_val = float(c_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for c_ij: {c_ij_} for key: {key_}")

                        # val
                        val_ = a_ij_val + b_ij_val * temperature + \
                            c_ij_val * pow(temperature, 2)

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

    def cal_tau_ij_M1(
        self,
        temperature: float,
        dg_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        dg_ij_symbol: Literal[
            'dg', 'dg_ij'
        ] = 'dg',
        R_CONST: float = 8.314,
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
        dg_ij : np.ndarray | Dict[str, float] | TableMatrixData
            Interaction energy parameter [J/mol] matrix where dg_ij[i][j] between component i and j.
        dg_ij_symbol : str
            Interaction energy parameter symbol. Default is 'dg'.
        R_CONST : float
            Universal gas constant [J/mol/K], default R_CONST = 8.314
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
            if (
                not isinstance(dg_ij, np.ndarray) and
                not isinstance(dg_ij, dict) and
                not isinstance(dg_ij, TableMatrixData)
            ):
                raise TypeError(
                    "dg_ij must be numpy array, dict or TableMatrixData")

            # Get the number of components
            comp_num = self.comp_num

            # components
            components = self.components

            # Initialize tau_ij matrix
            tau_ij = np.zeros((comp_num, comp_num), dtype=float)

            # tau_ij components
            tau_ij_comp = {}

            # check delimiter
            if symbol_delimiter == "|":
                symbol_delimiter_set = " | "
            elif symbol_delimiter == "_":
                symbol_delimiter_set = "_"
            else:
                raise ValueError("symbol_delimiter must be '|' or '_'")

            # Calculate tau_ij values
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
                        dg_ij_ = dg_ij.ij(key_)

                        # val
                        if dg_ij_ is not None and dg_ij_["value"] is not None:
                            val_ = float(dg_ij_["value"])
                        else:
                            raise ValueError(
                                f"Invalid value for dg_ij: {dg_ij_} for key: {key_}")

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

    def cal_tau_ij_M2(
        self, temperature: float,
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

            `tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T`

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

    def cal_G_ij(
        self,
        tau_ij: np.ndarray,
        alpha_ij: np.ndarray,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
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

    @add_attributes(metadata=ACTIVITY_MODELS['NRTL'])
    def cal(
        self,
        model_input: Dict,
        calculation_mode: Literal[
            'V1', 'V2'
        ] = 'V1',
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        message: Optional[str] = None,
        **kwargs
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        '''
        Calculate activity coefficients for a multi-component mixture using the NRTL model.

        Parameters
        ----------
        model_input: Dict
            Dictionary of model input values where keys are parameter names and values are their respective values.
                - `mole_fraction`: Dict[str, float]
                    Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
                - `temperature`: List[str | float], Optional
                    List of temperatures in any units as [300, 'K'], it is automatically converted to Kelvin.
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
        **kwargs: Optional
            Additional keyword arguments for the calculation.

        Returns
        -------
        res: Dict[str, float | Dict]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients as:
                - property_name: str
                    Name of the property calculated.
                - components: List[str]
                    List of component names.
                - mole_fraction: Dict[str, float]
                    Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
                - value: Dict[str, float]
                    Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
                - unit: float
                    Unit of the property calculated.
                - symbol: str
                    Symbol of the property calculated.
                - message: str
                    Message to be displayed.
        other_values: Dict[str, float | Dict]
            Dictionary of other values used for the calculation as:
                - AcCo_i_comp: Dict[str, float]
                - tau_ij: np.ndarray
                - tau_ij_comp: Dict[str, float]
                - alpha_ij: np.ndarray
                - alpha_ij_comp: Dict[str, float]
                - G_ij: np.ndarray
                - G_ij_comp: Dict[str, float]
                - calculation_mode: str

        Notes
        -----
        The activity coefficients are calculated using the NRTL model based on the provided input parameters as:

        - `tau_ij`: np.ndarray
        - `alpha_ij`: np.ndarray

        If the `tau_ij` and `alpha_ij` are not provided, they are generated using the `inputs_generator` method.
        The `inputs_generator` method generates the required input parameters based on the provided temperature and model input values.

        Examples
        --------
        >>> model_input = {
        ...     'mole_fraction': {'A': 0.5, 'B': 0.5},
        ...     'tau_ij': np.array([[0, 1], [1, 0]]),
        ...     'alpha_ij': np.array([[0, 0.5], [0.5, 0]])
        ... }
        >>> calculation_mode = 'V1'
        >>> message = 'Calculating activity coefficients'
        >>> result = activity_nrtl.cal(model_input, calculation_mode, message)
        >>> print(result)

        ```python
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
            'components': components,
            'mole_fraction': xi,
            'value': AcCo_i,
            'unit': 1,
            'symbol': "AcCo_i",
            'message': message,
        }
        ```
        '''
        try:
            # SECTION: check
            if not isinstance(model_input, dict):
                raise TypeError("model_input must be dict")

            # SECTION: check keys
            required_keys = ['tau_ij', 'alpha_ij']

            # ? check mole_fraction
            if 'mole_fraction' not in model_input:
                raise KeyError("mole_fraction is required in model_input")

            # set
            mole_fraction = model_input['mole_fraction']

            # ? checking alpha_ij and tau_ij
            # ! user should provide the required keys
            missed_keys = [
                key for key in required_keys if key not in model_input
            ]

            # check required keys
            if len(missed_keys) > 0:
                # check temperature
                if 'temperature' not in model_input:
                    # error
                    raise KeyError("temperature is required in model_input")

                # check if temperature is list
                if not isinstance(model_input['temperature'], list):
                    # error
                    raise TypeError("temperature must be list")

                # check if temperature is empty
                if len(model_input['temperature']) == 0:
                    # error
                    raise ValueError("temperature list is empty")

                # check format as [300, 'K']
                if not all(isinstance(temp, (int, float, str)) for temp in model_input['temperature']):
                    # error
                    raise TypeError("temperature list must be int or float")

                # call input generator
                inputs_ = self.inputs_generator(
                    temperature=model_input['temperature'],
                    model_input=model_input,
                )

                # looping through the missed keys
                for key in missed_keys:
                    # key value
                    value_ = inputs_[key]

                    # check
                    if value_ is None:
                        # error
                        raise ValueError(
                            f"{key} is required in model_input")

                    # update the model_input
                    model_input[key] = value_

            # SECTION: get values
            tau_ij_data = model_input['tau_ij']
            alpha_ij_data = model_input['alpha_ij']

            # SECTION: calculate activity coefficients
            return self.__calculate_activity_coefficients(
                mole_fraction=mole_fraction,
                tau_ij_data=tau_ij_data,
                alpha_ij_data=alpha_ij_data,
                calculation_mode=calculation_mode,
                symbol_delimiter=symbol_delimiter,
                message=message
            )
        except Exception as e:
            raise Exception(f"Error in launch_calculation: {str(e)}")

    def __calculate_activity_coefficients(
        self,
        mole_fraction: Dict[str, float],
        tau_ij_data: TableMatrixData | np.ndarray | Dict[str, float] | List[List[float | int | str]],
        alpha_ij_data: TableMatrixData | np.ndarray | Dict[str, float] | List[List[float | int | str]],
        calculation_mode: Literal['V1', 'V2'],
        symbol_delimiter: Literal["|", "_"],
        message: Optional[str],
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Calculate activity coefficients for a multi-component mixture using the NRTL model.

        Parameters
        -----------
        mole_fraction : Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
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

        Returns
        --------
        res : Dict[str, float | Dict]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
        other_values : Dict[str, float | Dict]
            Dictionary of other values used for the calculation as:
                - AcCo_i_comp: Dict[str, float]
                - tau_ij: np.ndarray
                - tau_ij_comp: Dict[str, float]
                - alpha_ij: np.ndarray
                - alpha_ij_comp: Dict[str, float]
                - G_ij: np.ndarray
                - G_ij_comp: Dict[str, float]
                - calculation_mode: str

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

            # NOTE: store in class
            self.__xi = xi
            self.__mole_fraction = mole_fraction

            # check message
            if message is None:
                message = f"Calculate activity coefficients for {components_str} using NRTL model"

            # SECTION
            # set the interaction parameter matrix (tau_ij) for the NRTL model
            if isinstance(tau_ij_data, np.ndarray):  # ! numpy array
                # set
                tau_ij = tau_ij_data
                # to dict
                tau_ij_comp = self.to_dict_ij(
                    tau_ij_data,
                    symbol_delimiter=symbol_delimiter
                )
            elif isinstance(tau_ij_data, list):  # ! list
                # convert list to numpy array
                tau_ij = np.array(tau_ij_data)
                # to dict
                tau_ij_comp = self.to_dict_ij(
                    tau_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif isinstance(tau_ij_data, TableMatrixData):  # ! PyThermoDB
                # convert to numpy array and dict
                res_ = self.to_ij(
                    data=tau_ij_data,
                    prop_symbol="tau"
                )
                # set
                tau_ij = res_[0]
                # to dict
                tau_ij_comp = res_[1]
            elif isinstance(tau_ij_data, dict):  # ! dict
                # convert dict to numpy array
                tau_ij = self.to_matrix_ij(
                    data=tau_ij_data,
                    symbol_delimiter=symbol_delimiter
                )
                # to dict
                tau_ij_comp = tau_ij_data
            else:
                raise TypeError(
                    "tau_ij_data must be numpy array, dict or TableMatrixData")

            # NOTE: store in class
            self.__tau_ij = tau_ij
            self.__tau_ij_comp = tau_ij_comp

            # SECTION
            # set the non-randomness parameter matrix (alpha_ij) for the NRTL model
            if isinstance(alpha_ij_data, np.ndarray):  # ! numpy array
                # set
                alpha_ij = alpha_ij_data
                # to dict
                alpha_ij_comp = self.to_dict_ij(
                    alpha_ij_data,
                    symbol_delimiter=symbol_delimiter
                )
            elif isinstance(alpha_ij_data, list):  # ! list
                # convert list to numpy array
                alpha_ij = np.array(alpha_ij_data)
                # to dict
                alpha_ij_comp = self.to_dict_ij(
                    alpha_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif isinstance(alpha_ij_data, TableMatrixData):  # ! PyThermoDB
                # convert to numpy array and dict
                res_ = self.to_ij(
                    data=alpha_ij_data,
                    prop_symbol="alpha")
                # set
                alpha_ij = res_[0]
                # to dict
                alpha_ij_comp = res_[1]
            elif isinstance(alpha_ij_data, dict):  # ! dict
                # convert dict to numpy array
                alpha_ij = self.to_matrix_ij(
                    data=alpha_ij_data,
                    symbol_delimiter=symbol_delimiter
                )
                # to dict
                alpha_ij_comp = alpha_ij_data
            else:
                raise TypeError(
                    "alpha_ij_data must be numpy array, dict or TableMatrixData")

            # NOTE: store in class
            self.__alpha_ij = alpha_ij
            self.__alpha_ij_comp = alpha_ij_comp

            # SECTION
            # set G_ij matrix for NRTL model
            G_ij, G_ij_comp = self.cal_G_ij(
                tau_ij=tau_ij,
                alpha_ij=alpha_ij,
                symbol_delimiter=symbol_delimiter
            )

            # NOTE: store in class
            self.__G_ij = G_ij
            self.__G_ij_comp = G_ij_comp

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
                'components': components,
                'mole_fraction': xi,
                'value': AcCo_i,
                'unit': 1,
                'symbol': "AcCo_i",
                'message': message,
            }

            # res
            return res, other_values

        except Exception as e:
            raise Exception(
                f"Error in calculate_activity_coefficients: {str(e)}")

    def CalAcCo_V1(
        self,
        xi: List[float],
        tau_ij: np.ndarray,
        G_ij: np.ndarray
    ) -> np.ndarray:
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

    def CalAcCo_V2(
        self,
        xi: list[float],
        tau_ij: np.ndarray,
        G_ij: np.ndarray
    ) -> np.ndarray:
        """
        Calculate activity coefficients for a multi-component mixture using the NRTL model.

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

    def excess_gibbs_free_energy(
        self,
        mole_fraction: Optional[
            Dict[str, float]
        ] = None,
        G_ij: Optional[np.ndarray] = None,
        tau_ij: Optional[np.ndarray] = None,
        message: Optional[str] = None,
        res_format: Literal[
            'str', 'json', 'dict'
        ] = 'dict'
    ) -> Dict[str, float | Dict] | str:
        """
        Calculate excess Gibbs energy (G^E/RT) for a multi-component mixture using the NRTL model.

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
            # NOTE: components
            components = self.components
            components_str = ', '.join(components)

            # NOTE: mole fraction
            # check
            if mole_fraction is None:
                mole_fraction = self.__mole_fraction

            # check
            if not isinstance(mole_fraction, dict):
                raise TypeError("mole_fraction must be dict")

            # set
            xi = [mole_fraction[components[i]] for i in range(len(components))]

            # NOTE: G_ij
            # check
            if G_ij is None:
                G_ij = self.__G_ij
            # list
            if isinstance(G_ij, list):
                G_ij = np.array(G_ij)
            # check array
            if not isinstance(G_ij, np.ndarray):
                raise TypeError("G_ij must be numpy array")

            # NOTE: tau_ij
            # check
            if tau_ij is None:
                tau_ij = self.__tau_ij
            # list
            if isinstance(tau_ij, list):
                tau_ij = np.array(tau_ij)
            # check array
            if not isinstance(tau_ij, np.ndarray):
                raise TypeError("tau_ij must be numpy array")

            # NOTE: set message
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
                "property_name": "Excess Molar Gibbs Free Energy (G^E/RT)",
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

    def inputs_generator(
        self,
        temperature: Optional[
            List[float | str]
        ] = None,
        **kwargs
    ):
        '''
        Prepares inputs for the NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: check src
            # check NRTL & nrtl keys in datasource
            if "NRTL" in self.datasource.keys():
                datasource = self.datasource["NRTL"]
            elif "nrtl" in self.datasource.keys():
                datasource = self.datasource["nrtl"]
            elif (
                self._mixture_id is not None and
                self._mixture_id in self.datasource.keys()
            ):
                datasource = self.datasource[self._mixture_id]
            else:
                # log
                logger.warning(
                    "No NRTL or nrtl key found in datasource, using model_input if provided."
                )
                datasource = {}

            # NOTE: check model inputs
            if kwargs.get('model_input') is not None:
                # update the datasource
                datasource.update(kwargs['model_input'])

            # ! set initial values
            a_ij = None
            b_ij = None
            c_ij = None
            d_ij = None
            dg_ij = None
            alpha_ij = None
            tau_ij = None

            # NOTE: check if datasource is a dictionary
            if datasource is not None:
                # check if datasource is a dictionary
                if not isinstance(datasource, dict):
                    raise ValueError(
                        "datasource must be a dictionary."
                    )

                # check if datasource is empty
                if len(datasource) == 0:
                    raise ValueError(
                        "datasource cannot be empty."
                    )

            # NOTE: check temperature
            # init temperature [K]
            T = -1  # invalid temperature

            # >> check
            if temperature is not None:
                # check if temperature is a list
                if not isinstance(temperature, list):
                    raise ValueError(
                        "temperature must be a list of floats or strings."
                    )

                # temperature
                T_value = float(temperature[0])
                T_unit = str(temperature[1])

                # convert temperature to Kelvin
                T = pycuc.convert_from_to(
                    T_value,
                    T_unit,
                    'K'
                )

            # NOTE: method 1
            # ! Δg_ij, interaction energy parameter
            dg_ij_src = datasource.get(
                'dg_ij',
                None
            )
            if dg_ij_src is None:
                dg_ij_src = datasource.get(
                    'dg',
                    None
                )

            # NOTE: method 2
            # ! constants a, b, c, and d
            a_ij_src = datasource.get('a_ij', None)
            if a_ij_src is None:
                a_ij_src = datasource.get('a', None)
            b_ij_src = datasource.get('b_ij', None)
            if b_ij_src is None:
                b_ij_src = datasource.get('b', None)
            c_ij_src = datasource.get('c_ij', None)
            if c_ij_src is None:
                c_ij_src = datasource.get('c', None)
            d_ij_src = datasource.get('d_ij', None)
            if d_ij_src is None:
                d_ij_src = datasource.get('d', None)

            # NOTE: α_ij, non-randomness parameter
            # ! check if alpha_ij is provided
            alpha_ij_src = datasource.get('alpha_ij', None)
            if alpha_ij_src is None:
                alpha_ij_src = datasource.get('alpha', None)

            # NOTE: tau_ij, binary interaction parameter
            # ! check if tau_ij is provided
            tau_ij_src = datasource.get('tau_ij', None)
            if tau_ij_src is None:
                tau_ij_src = datasource.get('tau', None)

            # SECTION: extract data
            # NOTE: check method
            tau_ij_cal_method = 0

            # NOTE: check if dg_ij is provided
            if dg_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if (
                    a_ij_src is None or
                    b_ij_src is None or
                    c_ij_src is None or
                    d_ij_src is None
                ):  # SECTION: check if a_ij, b_ij, c_ij, d_ij are provided
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', self.components)
                elif isinstance(a_ij_src, list):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', self.components)
                elif isinstance(b_ij_src, list):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', self.components)
                elif isinstance(c_ij_src, list):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', self.components)
                elif isinstance(d_ij_src, list):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dg_ij_src is not None:  # SECTION: check if dg_ij is provided
                # ! use dg_ij
                if isinstance(dg_ij_src, TableMatrixData):
                    dg_ij = dg_ij_src.mat('dg', self.components)
                elif isinstance(dg_ij_src, list):
                    dg_ij = np.array(dg_ij_src)
                elif isinstance(dg_ij_src, np.ndarray):
                    dg_ij = dg_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, list of lists, or numpy array.")
                # set method
                tau_ij_cal_method = 1
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, d."
                )

            # SECTION: extract data
            # NOTE: α_ij, non-randomness parameter
            # check
            if alpha_ij_src is not None or alpha_ij_src != 'None':
                if isinstance(alpha_ij_src, TableMatrixData):
                    alpha_ij = alpha_ij_src.mat('alpha', self.components)
                elif isinstance(alpha_ij_src, list):
                    alpha_ij = np.array(alpha_ij_src)
                elif isinstance(alpha_ij_src, np.ndarray):
                    alpha_ij = alpha_ij_src
                else:
                    raise ValueError(
                        "Invalid source for non-randomness parameter (α_ij). Must be TableMatrixData, list of lists, or numpy array.")
            else:
                # set default value
                alpha_ij = None

            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            # check
            if (
                tau_ij_src is None or
                tau_ij_src == 'None'
            ):
                # ! tau_ij is None
                # ? check method
                if tau_ij_cal_method == 1:
                    # Check if dg_ij is None and convert values to float if needed
                    if dg_ij is None:
                        raise ValueError(
                            "dg_ij cannot be None for calculating tau_ij")

                    # If dg_ij is a dictionary with mixed value types, convert all values to float
                    if isinstance(dg_ij, np.ndarray):
                        dg_ij = dg_ij.astype(float)
                    else:
                        raise ValueError(
                            "dg_ij must be a numpy array")

                    # >> calculate
                    tau_ij, tau_ij_comp = self.cal_tau_ij_M1(
                        temperature=T,
                        dg_ij=dg_ij
                    )
                elif tau_ij_cal_method == 2:
                    # check if a_ij, b_ij, c_ij, d_ij are None
                    if (
                        a_ij is None or
                        b_ij is None or
                        c_ij is None or
                        d_ij is None
                    ):
                        raise ValueError(
                            "a_ij, b_ij, c_ij, d_ij cannot be None for calculating tau_ij")

                    # If a_ij, b_ij, c_ij, d_ij are numpy array with mixed value types, convert all values to float
                    if isinstance(a_ij, np.ndarray):
                        a_ij = a_ij.astype(float)
                    else:
                        raise ValueError(
                            "a_ij must be a numpy array")

                    if isinstance(b_ij, np.ndarray):
                        b_ij = b_ij.astype(float)
                    else:
                        raise ValueError(
                            "b_ij must be a numpy array")

                    if isinstance(c_ij, np.ndarray):
                        c_ij = c_ij.astype(float)
                    else:
                        raise ValueError(
                            "c_ij must be a numpy array")

                    if isinstance(d_ij, np.ndarray):
                        d_ij = d_ij.astype(float)
                    else:
                        raise ValueError(
                            "d_ij must be a numpy array")

                    # >> calculate
                    tau_ij, tau_ij_comp = self.cal_tau_ij_M2(
                        temperature=T,
                        a_ij=a_ij,
                        b_ij=b_ij,
                        c_ij=c_ij,
                        d_ij=d_ij
                    )
                else:
                    raise ValueError(
                        "tau_ij_cal_method not supported!"
                    )
            else:
                # ! check if tau_ij is provided
                # check types
                if isinstance(tau_ij_src, TableMatrixData):
                    tau_ij = tau_ij_src.mat('tau', self.components)
                elif isinstance(tau_ij_src, list):
                    tau_ij = np.array(tau_ij_src)
                elif isinstance(tau_ij_src, np.ndarray):
                    tau_ij = tau_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (tau_ij). Must be TableMatrixData, list of lists, or numpy array.")

            # NOTE: nrtl inputs
            inputs = {
                "alpha_ij": alpha_ij,
                "tau_ij": tau_ij,
                "dg_ij": dg_ij,
                "a_ij": a_ij,
                "b_ij": b_ij,
                "c_ij": c_ij,
                "d_ij": d_ij,
            }

            # res
            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to generate NRTL activity inputs: {e}") from e
