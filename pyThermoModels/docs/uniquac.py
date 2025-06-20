# import libs
import numpy as np
import json
import yaml
from math import pow, exp, log
from typing import List, Dict, Tuple, Any, Literal, Optional
import pycuc
from pyThermoDB import (
    TableMatrixData,
)
# local
from ..plugin import ACTIVITY_MODELS
from ..utils import add_attributes


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

    # NOTE: variable based on the component id
    __tau_ij = None
    __tau_ij_comp = None
    __dU_ij = None
    __dU_ij_comp = None
    __r_i = None
    __r_i_comp = None
    __q_i = None
    __q_i_comp = None
    __mole_fraction = None
    __xi = None

    def __init__(
        self,
        components: List[str],
        datasource: Dict = {},
        equationsource: Dict = {}
    ):
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

    def parse_model_inputs(self, model_inputs: str) -> Dict[str, Any]:
        '''
        Convert model inputs from string to dictionary format.

        Parameters
        ----------
        model_inputs: str
            Model inputs in string format, such as:
            - { mole_fraction: { ethanol: 0.4, butyl-methyl-ether: 0.6 }, temperature: [323.15, 'K'], tau_ij: [[],[]], r_i: [[],[]] }

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
        self, data: TableMatrixData,
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

                    # key
                    key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"
                    # val
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

    def to_dict_i(
        self,
        data: List[float] | np.ndarray
    ) -> Dict[str, float]:
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

    def to_matrix_ij(
        self,
        data: Dict[str, float] | List[float],
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> np.ndarray:
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

            # SECTION: Set the interaction energy parameter matrix
            # check if data is dict or list
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

    def cal_dU_ij_M1(
        self, temperature: float,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
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
            if (isinstance(a_ij, np.ndarray) and
                isinstance(b_ij, np.ndarray) and
                    isinstance(c_ij, np.ndarray)):

                # loop over the components
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
            elif (
                isinstance(a_ij, dict) and
                isinstance(b_ij, dict) and
                isinstance(c_ij, dict)
            ):

                # loop over the components
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
            elif (
                isinstance(a_ij, TableMatrixData) and
                isinstance(b_ij, TableMatrixData) and
                isinstance(c_ij, TableMatrixData)
            ):

                # loop over the components
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

    def cal_tau_ij_M1(
        self,
        temperature: float,
        dU_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        dU_ij_symbol:
        Literal[
            'dU', 'dU_ij'
        ] = 'dU',
        R_CONST: float = 8.314,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
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
                        dU_ij_ = dU_ij.ij(key_)

                        # check
                        if dU_ij_ is None or dU_ij_["value"] is None:
                            raise ValueError(
                                f"Invalid value for {dU_ij_symbol}: {dU_ij_} for key: {key_}")

                        val_ = float(dU_ij_["value"])

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
            if (isinstance(a_ij, np.ndarray) and
                isinstance(b_ij, np.ndarray) and
                isinstance(c_ij, np.ndarray) and
                    isinstance(d_ij, np.ndarray)):

                # loop over the components
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
            elif (isinstance(a_ij, dict) and
                  isinstance(b_ij, dict) and
                  isinstance(c_ij, dict) and
                  isinstance(d_ij, dict)):

                # loop over the components
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
            elif (isinstance(a_ij, TableMatrixData) and
                  isinstance(b_ij, TableMatrixData) and
                  isinstance(c_ij, TableMatrixData) and
                  isinstance(d_ij, TableMatrixData)):

                # looping over the components
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

    def __X_ij(
            self,
            ij_data: TableMatrixData | np.ndarray | Dict[str, float] | List[List[float]],
            prop_symbol: str,
            symbol_delimiter: Literal[
                "|", "_"
            ] = "|"):
        """
        Convert interaction parameter data to numpy array and dict.

        Parameters
        ----------
        ij_data : TableMatrixData | np.ndarray | Dict[str, float] | List[List[float]]
            Interaction parameters (tau_ij) between component i and j.
        prop_symbol : str
            Interaction parameter symbol.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        ij_array : np.ndarray
            Interaction parameter matrix (numpy array).
        ij_comp : Dict[str, float]
            Dictionary of interaction parameters where keys are component pairs and values are their respective values.
        """
        try:
            # SECTION
            # set the interaction parameter matrix (such as tau_ij) for the UNIQUAC model
            if isinstance(ij_data, np.ndarray):
                # set
                ij_array = ij_data
                # to dict
                ij_comp = self.to_dict_ij(
                    ij_data, symbol_delimiter=symbol_delimiter)
            elif isinstance(ij_data, TableMatrixData):
                # prop symbol
                prop_symbol = prop_symbol.strip()
                # convert to numpy array and dict
                res_ = self.to_ij(
                    data=ij_data,
                    prop_symbol=prop_symbol)
                # set
                ij_array = res_[0]
                # to dict
                ij_comp = res_[1]
            elif isinstance(ij_data, dict):
                # convert dict to numpy array
                ij_array = self.to_matrix_ij(
                    data=ij_data, symbol_delimiter=symbol_delimiter)
                # to dict
                ij_comp = ij_data
            elif isinstance(ij_data, List):
                # convert list to numpy array
                ij_array = np.array(ij_data)
                # to dict
                ij_comp = self.to_dict_ij(
                    ij_array, symbol_delimiter=symbol_delimiter)
            else:
                raise TypeError(
                    "tau_ij_data must be numpy array, dict or TableMatrixData")

            # res
            return ij_array, ij_comp
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def __X_i(
            self,
            i_data: List[float] | Dict[str, float] | np.ndarray,
            symbol_delimiter: Literal[
                "|", "_"
            ] = "|"):
        '''
        Convert interaction parameter data to numpy array and dict.

        Parameters
        ----------
        i_data : List[float] | Dict[str, float] | np.ndarray
            Interaction parameters (r_i or q_i) for component i.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        i_array : np.ndarray
            Interaction parameter matrix (numpy array).
        i_comp : Dict[str, float]
            Dictionary of interaction parameters where keys are component pairs and values are their respective values.
        '''
        try:
            # SECTION
            # set relative surface area of component i (such as q_i) for the UNIQUAC model
            if isinstance(i_data, np.ndarray):
                # set
                i_array = i_data
                # to dict
                i_comp = self.to_dict_i(i_data)
            elif isinstance(i_data, List):
                # set
                i_array = np.array(i_data)
                # to dict
                i_comp = self.to_dict_i(i_data)
            elif isinstance(i_data, dict):
                # convert dict to numpy array
                i_array = self.to_i(data=i_data)
                # to dict
                i_comp = i_data
            else:
                raise TypeError("q_i_data must be numpy array, dict or List")

            # res
            return i_array, i_comp

        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    @add_attributes(metadata=ACTIVITY_MODELS['UNIQUAC'])
    def cal(
        self,
        model_input: Dict,
        Z: Optional[float | int] = None,
        calculation_mode: Literal['V1'] = 'V1',
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        message: Optional[str] = None,
        **kwargs
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
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
                - `temperature`: List[str | float], Optional
                    List of temperatures in any units as [300, 'K'], it is automatically converted to Kelvin.
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
        **kwargs : Optional
            Additional keyword arguments.

        Returns
        -------
        res: Dict[str, float | Dict]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients as:
                - property_name: str
                    Name of the property. Default is 'activity coefficients'.
                - components: List[str]
                    List of component names.
                - mole_fraction: Dict[str, float]
                    Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
                - value: np.ndarray
                    Activity coefficients (AcCo_i) for each component.
                - unit: int
                    Unit of the property. Default is 1.
                - symbol: str
                    Symbol of the property. Default is 'AcCo_i'.
                - message: str
                    Message to be displayed. Default is None.
        other_values: Dict[str, Any]
            Dictionary of other values used for the calculation as:
                - AcCo_i_comp: Dict[str, float]
                    Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
                - tau_ij: np.ndarray
                    Interaction parameters (tau_ij) between component i and j.
                - tau_ij_comp: Dict[str, float]
                    Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.
                - r_i: np.ndarray
                    relative van der Waals volume of component i
                - r_i_comp: Dict[str, float]
                    Dictionary of relative van der Waals volume of component i where keys are component names and values are their respective values.
                - q_i: np.ndarray
                    relative surface area of component i
                - q_i_comp: Dict[str, float]
                    Dictionary of relative surface area of component i where keys are component names and values are their respective values.
                - calculation_mode: str
                    Mode of calculation. If 'V1', use the first version of the UNIQUAC model. If 'V2', use the second version.

        Examples
        --------
        >>> model_input = {
        ...     'mole_fraction': {'A': 0.5, 'B': 0.5},
        ...     'tau_ij': np.array([[0, 1], [1, 0]]),
        ...     'r_i': r_i,
        ...     'q_i': q_i
        ... }
        >>> calculation_mode = 'V1'
        >>> message = 'Calculating activity coefficients'
        >>> result = activity_uniquac.cal(model_input=model_input)
        >>> print(result)

        ```python
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
        ```
        """
        try:
            # SECTION: check
            if not isinstance(model_input, dict):
                raise TypeError("model_input must be dict")

            # SECTION: check keys
            required_keys = ['tau_ij', 'r_i', 'q_i']

            # ? check mole_fraction
            if 'mole_fraction' not in model_input:
                raise KeyError("mole_fraction is required in model_input")

            # set
            mole_fraction = model_input['mole_fraction']

            # ? checking alpha_ij and tau_ij
            # ! user should provide the required keys
            missed_keys = [
                key for key in required_keys if key not in model_input]

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
                message=message)
        except Exception as e:
            raise Exception(f"Error in uniquac model cal: {str(e)}")

    def __calculate_activity_coefficients(
        self,
        mole_fraction: Dict[str, float],
        tau_ij_data: TableMatrixData | np.ndarray | Dict[str, float] | List[List[float | int | str]],
        r_i_data: List[float] | Dict[str, float] | np.ndarray,
        q_i_data: List[float] | Dict[str, float] | np.ndarray,
        Z: Optional[float | int],
        calculation_mode: Literal['V1'],
        symbol_delimiter: Literal["|", "_"],
        message: Optional[str],
    ) -> Tuple[
        Dict[str, Any], Dict[str, Any]
    ]:
        """
        Calculate activity coefficients for a multi-component mixture using the UNIQUAC model.

        Parameters
        -----------
        mole_fraction : Dict[str, float]
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        tau_ij_comp : TableMatrixData | np.ndarray | Dict[str, float]
            Interaction parameters (tau_ij) between component i and j.
        r_i_data : List[float] | Dict[str, float] | np.ndarray
            relative van der Waals volume of component i
        q_i_data : List[float] | Dict[str, float] | np.ndarray
            relative surface area of component i
        Z : int | float
            Model constant, Default is 10.
        calculation_mode : Literal['V1', 'V2']
            Mode of calculation. If 'V1', use the first version of the UNIQUAC model. If 'V2', use the second version.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".
        message : Optional[str]
            Message to be displayed. Default is None.

        Returns
        --------
        res: Dict[str, float | Dict]
            Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
        other_values: Dict[str, Any]
            Dictionary of other values used for the calculation as:
                - AcCo_i_comp: Dict[str, float]
                    Dictionary of activity coefficients where keys are component names and values are their respective activity coefficients.
                - tau_ij: np.ndarray
                    Interaction parameters (tau_ij) between component i and j.
                - tau_ij_comp: Dict[str, float]
                    Dictionary of interaction parameters where keys are component pairs and values are their respective tau_ij values.
                - r_i: np.ndarray
                    relative van der Waals volume of component i
                - r_i_comp: Dict[str, float]
                    Dictionary of relative van der Waals volume of component i where keys are component names and values are their respective values.
                - q_i: np.ndarray
                    relative surface area of component i
                - q_i_comp: Dict[str, float]
                    Dictionary of relative surface area of component i where keys are component names and values are their respective values.
                - calculation_mode: str
                    Mode of calculation. If 'V1', use the first version of the UNIQUAC model. If 'V2', use the second version.

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

            # NOTE: store class variables
            self.__xi = xi
            self.__mole_fraction = mole_fraction

            # SECTION
            # set the interaction parameter matrix (tau_ij) for the UNIQUAC model
            if isinstance(tau_ij_data, np.ndarray):  # ! numpy array
                # set
                tau_ij = tau_ij_data
                # to dict
                tau_ij_comp = self.to_dict_ij(
                    tau_ij_data, symbol_delimiter=symbol_delimiter)
            elif isinstance(tau_ij_data, TableMatrixData):  # ! TableMatrixData
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
            elif isinstance(tau_ij_data, list):  # ! list
                # convert list to numpy array
                tau_ij = np.array(tau_ij_data)
                # to dict
                tau_ij_comp = self.to_dict_ij(
                    tau_ij, symbol_delimiter=symbol_delimiter)
            else:
                raise TypeError(
                    "tau_ij_data must be numpy array, dict or TableMatrixData")

            # NOTE: store class variables
            self.__tau_ij = tau_ij
            self.__tau_ij_comp = tau_ij_comp

            # SECTION
            # set relative van der Waals volume of component i (r_i) for the UNIQUAC model
            if isinstance(r_i_data, np.ndarray):  # ! numpy array
                # set
                r_i = r_i_data
                # to dict
                r_i_comp = self.to_dict_i(r_i_data)
            elif isinstance(r_i_data, list):  # ! list
                # set
                r_i = np.array(r_i_data)
                # to dict
                r_i_comp = self.to_dict_i(r_i_data)
            elif isinstance(r_i_data, dict):  # ! dict
                # convert dict to numpy array
                r_i = self.to_i(data=r_i_data)
                # to dict
                r_i_comp = r_i_data
            else:
                raise TypeError("r_i_data must be numpy array, dict or List")

            # NOTE: store class variables
            self.__r_i = r_i
            self.__r_i_comp = r_i_comp

            # SECTION
            # set relative surface area of component i (q_i) for the UNIQUAC model
            if isinstance(q_i_data, np.ndarray):  # ! numpy array
                # set
                q_i = q_i_data
                # to dict
                q_i_comp = self.to_dict_i(q_i_data)
            elif isinstance(q_i_data, list):  # ! list
                # set
                q_i = np.array(q_i_data)
                # to dict
                q_i_comp = self.to_dict_i(q_i_data)
            elif isinstance(q_i_data, dict):  # ! dict
                # convert dict to numpy array
                q_i = self.to_i(data=q_i_data)
                # to dict
                q_i_comp = q_i_data
            else:
                raise TypeError("q_i_data must be numpy array, dict or List")

            # NOTE: store class variables
            self.__q_i = q_i
            self.__q_i_comp = q_i_comp

            # SECTION
            # Calculate activity coefficients using the UNIQUAC model
            if calculation_mode == 'V1':
                AcCo_i = self.CalAcCo_V1(
                    xi=xi,
                    tau_ij=tau_ij,
                    r_i=r_i,
                    q_i=q_i,
                    Z=Z)
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

            # res
            return res, other_values

        except Exception as e:
            raise Exception(
                f"Error in calculate_activity_coefficients: {str(e)}")

    def CalAcCo_V1(
        self,
        xi: List[float],
        tau_ij: np.ndarray,
        r_i: np.ndarray,
        q_i: np.ndarray,
        Z: int | float
    ) -> np.ndarray:
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
        AcCoi: np.ndarray
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

            # activity coefficient
            AcCoi = np.zeros(comp_num)

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

            # res
            return AcCoi
        except Exception as e:
            raise Exception(f"Error in CalAcCo_V1: {str(e)}")

    def excess_gibbs_free_energy(
        self,
        mole_fraction: Optional[Dict[str,
                                     float]] = None,
        tau_ij: Optional[np.ndarray] = None,
        r_i:  Optional[np.ndarray] = None,
        q_i:  Optional[np.ndarray] = None,
        Z: Optional[int | float] = None,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        message: Optional[str] = None,
        res_format: Literal[
            'str', 'json', 'dict'
        ] = 'dict'
    ) -> Dict[str, float | Dict] | str:
        """
        Calculate excess Gibbs energy (G^E/RT) for a multi-component mixture using the UNIQUAC model.

        Parameters
        -----------
        mole_fraction : dict | None
            Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
        tau_ij : np.ndarray | None
            Matrix of tau parameters where tau[i][j] is tau_ij between component i and j.
        r_i : np.ndarray | None
            Array of relative van der Waals volumes of each component.
        q_i : np.ndarray | None
            Array of relative surface areas of each component.
        Z : int | float, optional
            Model constant, default is 10.
        symbol_delimiter : Literal["|", "_"], optional
            Delimiter for the component id. Default is "|".
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
            if mole_fraction:
                if not isinstance(mole_fraction, dict):
                    raise TypeError("mole_fraction must be a dictionary")
            else:
                mole_fraction = self.__mole_fraction
                # check
                if mole_fraction is None:
                    raise ValueError("mole_fraction is not set")

            # check Z
            if Z is None:
                Z = self.Z

            # NOTE: components
            components = self.components
            components_str = ', '.join(components)

            # comp no
            comp_num = self.comp_num

            # NOTE: set message
            message = f'Excess Gibbs Free Energy for {components_str}' if message is None else message

            # NOTE: mole fraction
            xi = [mole_fraction[components[i]] for i in range(len(components))]

            # Normalize mole fractions to ensure they sum to 1
            x = xi / np.sum(xi)
            x = np.array(x)

            # NOTE: check all input
            if len(x) != comp_num:
                raise ValueError(
                    f"mole_fraction length {len(x)} does not match component number {comp_num}")

            # NOTE: set tau_ij
            if tau_ij is None:
                tau_ij = self.__tau_ij
                # check
                if tau_ij is None:
                    raise ValueError("tau_ij is not set")

            # check
            if isinstance(tau_ij, List):
                # convert to numpy array
                tau_ij = np.array(tau_ij)
            # check
            if not isinstance(tau_ij, np.ndarray):
                raise TypeError("tau_ij must be numpy array")
            # check
            if tau_ij.shape[0] != comp_num or tau_ij.shape[1] != comp_num:
                raise ValueError(
                    f"tau_ij shape {tau_ij.shape} does not match component number {comp_num}")

            # NOTE: set r_i
            if r_i is None:
                r_i = self.__r_i
                # check
                if r_i is None:
                    raise ValueError("r_i is not set")
            # check
            if len(r_i) != comp_num:
                raise ValueError(
                    f"r_i length {len(r_i)} does not match component number {comp_num}")

            # NOTE: set q_i
            if q_i is None:
                q_i = self.__q_i
                # check
                if q_i is None:
                    raise ValueError("q_i is not set")
            # check
            if len(q_i) != comp_num:
                raise ValueError(
                    f"q_i length {len(q_i)} does not match component number {comp_num}")

            # set array
            r_i = np.array(r_i)
            q_i = np.array(q_i)

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
            **kwargs):
        '''
        Prepares inputs for the UNIQUAC activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.

        Returns
        -------
        inputs : dict
            Dictionary of inputs for the UNIQUAC activity model.
            - tau_ij : np.ndarray
                Interaction parameters (tau_ij) between component i and j.
            - r_i : np.ndarray
                Relative van der Waals volume of component i.
            - q_i : np.ndarray
                Relative surface area of component i.
            a_ij : np.ndarray
                Interaction energy parameter (a_ij) between component i and j.
            b_ij : np.ndarray
                Interaction energy parameter (b_ij) between component i and j.
            c_ij : np.ndarray
                Interaction energy parameter (c_ij) between component i and j.
            d_ij : np.ndarray
                Interaction energy parameter (d_ij) between component i and j.
        '''
        try:
            # SECTION: check src
            # extract activity model inputs
            datasource = self.datasource.get('UNIQUAC', {})

            # NOTE: check model inputs
            if kwargs.get('model_input') is not None:
                # update the datasource
                datasource.update(kwargs['model_input'])

            # ! set initial values
            r_i = None
            q_i = None
            dU_ij = None
            a_ij = None
            b_ij = None
            c_ij = None
            d_ij = None

            # NOTE: check if datasource is a dictionary
            if datasource is not None:
                # check if datasource is a dictionary
                if not isinstance(datasource, dict):
                    raise ValueError(
                        "datasource must be a dictionary.")
                # check if datasource is empty
                if len(datasource) == 0:
                    raise ValueError(
                        "datasource cannot be empty.")

            # NOTE: check temperature
            if temperature is not None:
                # check if temperature is a list
                if not isinstance(temperature, list):
                    raise ValueError(
                        "temperature must be a list of floats or strings.")

                # temperature
                T_value = float(temperature[0])
                T_unit = str(temperature[1])

                # convert temperature to Kelvin
                T = pycuc.convert_from_to(
                    T_value, T_unit, 'K')

            # NOTE: method 1
            # ! Δg_ij, interaction energy parameter
            dU_ij_src = datasource.get('dU_ij', None)
            if dU_ij_src is None:
                dU_ij_src = datasource.get('dU', None)

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

            # NOTE: tau_ij, binary interaction parameter
            tau_ij_src = datasource.get('tau_ij', None)
            if tau_ij_src is None:
                tau_ij_src = datasource.get('tau', None)

            # NOTE: r_i, relative van der Waals volume of component i
            r_i_src = datasource.get('r_i', None)
            if r_i_src is None:
                # set default value
                r_i_src = datasource.get('r', None)

            # check if r_i is a list or numpy array
            if r_i_src is not None:
                if isinstance(r_i_src, list):
                    r_i = np.array(r_i_src)
                elif isinstance(r_i_src, np.ndarray):
                    r_i = r_i_src
                else:
                    raise ValueError(
                        "r_i must be a list or numpy array.")

            # NOTE: q_i, relative van der Waals area of component i
            q_i_src = datasource.get('q_i', None)
            if q_i_src is None:
                # set default value
                q_i_src = datasource.get('q', None)

            # check if q_i is a list or numpy array
            if q_i_src is not None:
                if isinstance(q_i_src, list):
                    q_i = np.array(q_i_src)
                elif isinstance(q_i_src, np.ndarray):
                    q_i = q_i_src
                else:
                    raise ValueError(
                        "q_i must be a list or numpy array.")

            # SECTION: extract data
            # NOTE: check method
            tau_ij_cal_method = 0

            # check if dU_ij, a_ij, b_ij, c_ij, d_ij are provided

            # ! check if dU_ij is None
            if dU_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if (a_ij_src is None or
                    b_ij_src is None or
                    c_ij_src is None or
                        d_ij_src is None):
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
            elif dU_ij_src is not None:
                # ! use dU_ij
                if isinstance(dU_ij_src, TableMatrixData):
                    dU_ij = dU_ij_src.mat('dU', self.components)
                elif isinstance(dU_ij_src, list):
                    dU_ij = np.array(dU_ij_src)
                elif isinstance(dU_ij_src, np.ndarray):
                    dU_ij = dU_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (dU_ij). Must be TableMatrixData, list of lists, or numpy array.")
                # set method
                tau_ij_cal_method = 1
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: calculate tau_ij
            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            # check
            if tau_ij_src is None or tau_ij_src == 'None':
                # ! tau_ij is None
                # ? check method
                if tau_ij_cal_method == 1:
                    # check if dU_ij is None
                    if dU_ij is None:
                        raise ValueError(
                            "dU_ij is not set. Cannot calculate tau_ij.")

                    # convert values to float
                    if isinstance(dU_ij, np.ndarray):
                        dU_ij = dU_ij.astype(float)
                    else:
                        raise ValueError(
                            "dU_ij must be numpy array.")

                    tau_ij, _ = self.cal_tau_ij_M1(
                        temperature=T,
                        dU_ij=dU_ij
                    )
                elif tau_ij_cal_method == 2:
                    # check if a_ij, b_ij, c_ij, d_ij are None
                    if (a_ij is None or
                        b_ij is None or
                        c_ij is None or
                            d_ij is None):
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

                    tau_ij, _ = self.cal_tau_ij_M2(
                        temperature=T,
                        a_ij=a_ij,
                        b_ij=b_ij,
                        c_ij=c_ij,
                        d_ij=d_ij)
                else:
                    raise ValueError(
                        "Invalid tau_ij_cal_method. Must be 1 or 2.")
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
                "r_i": r_i,
                "q_i": q_i,
                "dU_ij": dU_ij,
                "tau_ij": tau_ij,
                "a_ij": a_ij,
                "b_ij": b_ij,
                "c_ij": c_ij,
                "d_ij": d_ij,
            }

            # res
            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to generate UNIQUAC activity inputs: {e}") from e
