# import libs
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union, Literal
from pyThermoDB import (
    TableMatrixData,
    TableData,
    TableEquation,
    TableMatrixEquation
)
from pythermodb_settings.models import Component, ComponentKey
from pythermodb_settings.utils import set_component_id
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)


class ComponentParameterMixin:

    def __init__(
            self,
            components: List[str],
            comp_idx: Dict[str, int]
    ):
        """
        Initialize the ComponentParameterMixin class.

        Parameters
        ----------
        components : List[str]
            List of component names.
        comp_idx : Dict[str, int]
            Dictionary mapping component names to their respective indices.
        """
        self.components = components
        self.comp_idx = comp_idx
        self.comp_num = len(components)

    # SECTION: Transform TableMatrixData to numpy array and dictionary
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
                    key_ = f"{prop_symbol}_{self.components[i]}_{self.components[j]}"
                    # val
                    val = data.ij(key_)

                    # to matrix
                    # ? val["value"] or 0.0 in case of None
                    if val is not None and val["value"] is not None:
                        mat_ij[i, j] = float(val["value"])
                    else:
                        raise ValueError(
                            f"Invalid value for {prop_symbol}: {val} for key: {prop_symbol}_{self.components[i]}_{self.components[j]}")

                    # to dict
                    key_ = f"{self.components[i]}{symbol_delimiter_set}{self.components[j]}"
                    # to dict
                    # ? val["value"] or 0.0 in case of None
                    if val is not None and val["value"] is not None:
                        dict_ij[key_] = float(val["value"])
                    else:
                        raise ValueError(
                            f"Invalid value for {prop_symbol}: {val} for key: {prop_symbol}_{self.components[i]}_{self.components[j]}")

            # res
            return mat_ij, dict_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_i(
        self,
        data: Dict[str, float]
    ):
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
                data_i[i] = float(val)

            # res
            return data_i
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    # SECTION: Transform numpy array to dictionary
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

    # NOTE: wrapper for to_dict_ij with external component sources
    def to_dict_ij_ext(
        self,
        data: np.ndarray,
        components: List[Component],
        component_key: ComponentKey = "Name",
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ):
        try:
            # iterate over components and set component id
            component_list = [set_component_id(
                comp, component_key) for comp in components]

            # Get the number of components
            comp_num = len(component_list)

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
                    key_ = f"{component_list[i]}{symbol_delimiter_set}{component_list[j]}"
                    dict_ij[key_] = val

            # res
            return dict_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    # SECTION: Transform list or numpy array to dictionary
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

    # SECTION: Transform dictionary or list to numpy array
    def to_matrix_ij(
        self,
        data: Union[
            Dict[str, float | int],
            List[List[float | int]],
            TableMatrixData |
            np.ndarray,
        ],
        property_name: Optional[str] = None,
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|"
    ) -> np.ndarray:
        """
        Convert to matrix (mat_ij) according to `the component id`.

        Parameters
        ----------
        data : Dict[str, float | int] | List[List[float | int]] | TableMatrixData | np.ndarray
            Dictionary of parameters where keys are component pairs and values are their respective values or matrix-like list of values according to the component id.
        property_name : Optional[str], default=None
            Property name for the parameter. If provided, it will be used to extract values from TableMatrixData. If not provided, the function will attempt to use the keys in the data dictionary or the indices in the list.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".

        Returns
        -------
        mat_ij : np.ndarray
            Parameter matrix (numpy array).
        """
        try:
            # NOTE: check dict, list, or TableMatrixData
            if (
                not isinstance(data, dict) and
                not isinstance(data, list) and
                not isinstance(data, TableMatrixData) and
                not isinstance(data, np.ndarray)
            ):
                raise TypeError(
                    "data must be dict, list, TableMatrixData, or numpy array"
                )

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
            if isinstance(data, dict):
                # ! dict
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
                        mat_ij[comp_id_i, comp_id_j] = float(val)
            elif isinstance(data, list):
                # ! list
                for i in range(comp_num):
                    for j in range(comp_num):
                        # val
                        val = data[i][j]

                        # find the component id
                        comp_id_i = self.comp_idx[self.components[i]]
                        comp_id_j = self.comp_idx[self.components[j]]

                        # to matrix
                        mat_ij[comp_id_i, comp_id_j] = float(val)
            elif isinstance(data, TableMatrixData):
                # ! TableMatrixData
                if property_name is None:
                    raise ValueError(
                        "property_name must be provided for TableMatrixData")
                mat_ij = np.asarray(
                    data.mat(property_name, self.components),
                    dtype=float
                )
                if mat_ij.shape != (comp_num, comp_num):
                    raise ValueError(
                        f"data shape must be ({comp_num}, {comp_num}) for TableMatrixData"
                    )
            elif isinstance(data, np.ndarray):
                # ! numpy array
                if data.shape != (comp_num, comp_num):
                    raise ValueError(
                        f"data shape must be ({comp_num}, {comp_num}) for numpy array"
                    )

                mat_ij = data
            else:
                raise TypeError("data must be dict, list, or TableMatrixData")

            # res
            return mat_ij
        except Exception as e:
            raise Exception(f"Error in extraction data: {str(e)}")

    def to_matrix_ij_or(
            self,
            data: Union[
                Dict[str, float | int],
                List[List[float | int]],
                TableMatrixData |
                np.ndarray,
                None
            ],
            property_name: Optional[str] = None,
            symbol_delimiter: Literal[
                "|", "_"
            ] = "|"
    ):
        """
        Wrapper for to_matrix_ij with error handling
        """
        try:
            if data is None:
                raise ValueError("data cannot be None")

            return self.to_matrix_ij(
                data=data,
                property_name=property_name,
                symbol_delimiter=symbol_delimiter
            )
        except Exception as e:
            raise Exception(
                f"Failed to convert data to matrix for property '{property_name}': {e}"
            ) from e
