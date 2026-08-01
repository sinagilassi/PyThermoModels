# import libs
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union, Literal, TypeAlias, cast
import pycuc
from pyThermoDB import (
    TableMatrixData,
    TableData,
    TableEquation,
    TableMatrixEquation
)
# locals
from .nrtl_local_composition import NRTLLocalComposition
from .component_parameter_mixin import ComponentParameterMixin


# NOTE: set up logger for this module
logger = logging.getLogger(__name__)


class NRTLParameterBuilder:

    def __init__(
        self,
        components: List[str],
        comp_idx: Dict[str, int],
        datasource: Dict,
        equationsource: Dict,
        **kwargs
    ):
        # set
        self.components = components
        self.comp_idx = comp_idx
        self.comp_num = len(components)

        self.datasource = datasource
        self.equationsource = equationsource

        # SECTION: local composition model
        self.local_composition_model = NRTLLocalComposition(
            components=self.components,
            component_idx=self.comp_idx
        )

        # SECTION: component parameter mixin
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx
        )
        # ! set component access methods
        self.to_ij = self.component_parameter_mixin.to_ij
        self.to_dict_ij = self.component_parameter_mixin.to_dict_ij
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij

    # SECTION: Calculation methods

    def cal_dg_ij_M1(
        self,
        temperature: float,
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

            dg_ij = a_ij + b_ij * T + c_ij * T^2

            where T is the temperature [K].

        2. All parameters including a_ij, b_ij, c_ij must be in the same format (numpy array, dict or TableMatrixData).
        """
        try:
            return self.local_composition_model.cal_dg_ij_M1(
                temperature=temperature,
                a_ij=a_ij,
                b_ij=b_ij,
                c_ij=c_ij,
                symbol_delimiter=symbol_delimiter,
            )
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

        - `dg_{component_i}_{component_j}`
        - `dg | {component_i} | {component_j}`.
        """
        try:
            return self.local_composition_model.cal_tau_ij_M1(
                temperature=temperature,
                dg_ij=dg_ij,
                dg_ij_symbol=dg_ij_symbol,
                R_CONST=R_CONST,
                symbol_delimiter=symbol_delimiter,
            )
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

            tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T

        2. Interaction energy parameter symbol is `X` for TableMatrixData as:

        - `X_{component_i}_{component_j}`
        - `X | {component_i} | {component_j}`.

        2. All parameters including a_ij, b_ij, c_ij, d_ij must be in the same format (numpy array, dict or TableMatrixData).
        """
        try:
            return self.local_composition_model.cal_tau_ij_M2(
                temperature=temperature,
                a_ij=a_ij,
                b_ij=b_ij,
                c_ij=c_ij,
                d_ij=d_ij,
                symbol_delimiter=symbol_delimiter,
            )
        except Exception as e:
            raise Exception(f"Error in cal_tauij: {str(e)}")

    # SECTION: additional direct tau correlations
    def cal_tau_ij_M3(
        self,
        temperature: float,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau_ij using:

        SECTION: equation
        tau_ij = a_ij + b_ij / T
        """
        try:
            return self.local_composition_model.cal_tau_ij_M3(
                temperature=temperature,
                a_ij=a_ij,
                b_ij=b_ij,
                symbol_delimiter=symbol_delimiter,
            )
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M3: {str(e)}")

    def cal_tau_ij_M4(
        self,
        temperature: float,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau_ij using:

        SECTION: equation
        tau_ij = a_ij + b_ij / T + c_ij / T^2
        """
        try:
            return self.local_composition_model.cal_tau_ij_M4(
                temperature=temperature,
                a_ij=a_ij,
                b_ij=b_ij,
                c_ij=c_ij,
                symbol_delimiter=symbol_delimiter,
            )
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M4: {str(e)}")

    def cal_tau_ij_M5(
        self,
        temperature: float,
        a_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        b_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        c_ij: np.ndarray | Dict[str, float] | TableMatrixData,
        symbol_delimiter: Literal["|", "_"] = "|"
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate NRTL tau_ij using:

        SECTION: equation
        tau_ij = a_ij + b_ij / T + c_ij * ln(T)
        """
        try:
            return self.local_composition_model.cal_tau_ij_M5(
                temperature=temperature,
                a_ij=a_ij,
                b_ij=b_ij,
                c_ij=c_ij,
                symbol_delimiter=symbol_delimiter,
            )
        except Exception as e:
            raise Exception(f"Error in cal_tau_ij_M5: {str(e)}")

    # SECTION: inputs generator
    # NOTE: data source generator
    def data_source_generator(
        self,
        mixture_ids: Optional[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Generate a data source dictionary for NRTL activity model.

        Parameters
        ----------
        mixture_ids : Optional[Dict[str, str]], optional
            Dictionary containing mixture identifiers. Default is None.

        Returns
        -------
        Dict[str, Any]
            A dictionary containing the data source for the NRTL activity model.
        """
        try:
            # SECTION: check src
            # check NRTL & nrtl keys in datasource
            if "NRTL" in self.datasource.keys():
                # ! NRTL provided
                datasource = self.datasource["NRTL"]
            elif "nrtl" in self.datasource.keys():
                # ! nrtl provided
                datasource = self.datasource["nrtl"]
            elif mixture_ids is not None:
                # ! mixture_ids provided
                # init datasource
                datasource = {}
                # set datasource by mixture ids
                if 'Name' in mixture_ids.keys():
                    # check not empty
                    if mixture_ids.get('Name', None):
                        key_ = mixture_ids['Name']
                        # check key in datasource
                        if key_ in self.datasource.keys():
                            datasource = self.datasource[key_]
                elif 'Formula' in mixture_ids.keys():
                    # check not empty
                    if mixture_ids.get('Formula', None):
                        key_ = mixture_ids['Formula']
                        # check key in datasource
                        if key_ in self.datasource.keys():
                            datasource = self.datasource[key_]
            else:
                # ! no keys found, use model_input if provided
                # log
                logger.warning(
                    "No NRTL or nrtl key found in datasource, using model_input if provided."
                )
                datasource = {}

            if (
                datasource is not None and
                not isinstance(datasource, dict)
            ):
                raise ValueError(
                    "datasource must be a dictionary."
                )

            # NOTE: set datasource to empty dict if None
            datasource = {} if datasource is None else dict(datasource)

            # NOTE: check model inputs
            # ! when constants are provided in model_input, they override the datasource
            if kwargs.get('model_input') is not None:
                # update the datasource
                datasource.update(kwargs['model_input'])

            # NOTE: final datasource validation
            if datasource is None:
                raise ValueError(
                    "datasource cannot be None."
                )

            if not isinstance(datasource, dict):
                raise ValueError(
                    "datasource must be a dictionary."
                )

            if len(datasource) == 0:
                raise ValueError(
                    "datasource cannot be empty."
                )

            # return
            return datasource
        except Exception as e:
            raise Exception(
                f"Failed to generate NRTL activity data source: {e}"
            ) from e

    # NOTE: extract parameter source
    def extract_parameter_source(
            self,
            datasource: Dict[str, Any]
    ):
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

        return {
            'dg_ij_src': dg_ij_src,
            'a_ij_src': a_ij_src,
            'b_ij_src': b_ij_src,
            'c_ij_src': c_ij_src,
            'd_ij_src': d_ij_src,
            'alpha_ij_src': alpha_ij_src,
            'tau_ij_src': tau_ij_src
        }

    # NOTE: inputs generator

    def inputs_generator(
        self,
        temperature: Optional[
            List[float | str]
        ] = None,
        tau_correlation: Literal[
            'M1', 'M2', 'M3', 'M4', 'M5'
        ] = 'M1',
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        mixture_ids: Optional[Dict[str, str]] = None,
        **kwargs
    ):
        '''
        Prepares inputs for the NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for component-pair dictionary keys. Default is "|".
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: data source
            # generate datasource
            datasource: Dict[str, Any] = self.data_source_generator(
                mixture_ids=mixture_ids,
                **kwargs
            )

            # ! set initial values
            a_ij = None
            b_ij = None
            c_ij = None
            d_ij = None
            dg_ij = None
            alpha_ij = None
            tau_ij = None

            def _has_source(src) -> bool:
                return src is not None and src != 'None'

            def _to_matrix(src, symbol: str, parameter_name: str) -> np.ndarray:
                if isinstance(src, TableMatrixData):
                    return self.to_matrix_ij(
                        src,
                        property_name=symbol,
                        symbol_delimiter=symbol_delimiter
                    )
                if isinstance(src, list):
                    return np.array(src)
                if isinstance(src, np.ndarray):
                    return src
                if isinstance(src, dict):
                    return self.to_matrix_ij(
                        src,
                        symbol_delimiter=symbol_delimiter
                    )
                raise ValueError(
                    f"Invalid source for {parameter_name}. Must be TableMatrixData, dict, list of lists, or numpy array."
                )

            def _require_matrix(
                value: Optional[np.ndarray],
                parameter_name: str
            ) -> np.ndarray:
                if value is None:
                    raise ValueError(f"{parameter_name} cannot be None")
                return value

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

            # SECTION: extract parameter source
            parameter_sources = self.extract_parameter_source(datasource)
            # >> unpack
            dg_ij_src = parameter_sources['dg_ij_src']
            a_ij_src = parameter_sources['a_ij_src']
            b_ij_src = parameter_sources['b_ij_src']
            c_ij_src = parameter_sources['c_ij_src']
            d_ij_src = parameter_sources['d_ij_src']
            alpha_ij_src = parameter_sources['alpha_ij_src']
            tau_ij_src = parameter_sources['tau_ij_src']

            # SECTION: check correlation methods

            allowed_tau_correlations = ('M1', 'M2', 'M3', 'M4', 'M5')
            allowed_dg_correlations = ('M1',)

            if tau_correlation not in allowed_tau_correlations:
                raise ValueError(
                    f"tau_correlation must be one of {allowed_tau_correlations}."
                )

            # SECTION: extract data

            # NOTE: check if dg_ij is provided
            if (
                tau_ij_src is not None
            ):
                # ! use tau_ij
                # set method
                tau_ij_cal_method = 0

            elif (
                tau_correlation == 'M1' and
                _has_source(dg_ij_src)
            ):  # >>> check if dg_ij is provided
                # ! use dg_ij
                # set method
                tau_ij_cal_method = 1

                if isinstance(dg_ij_src, TableMatrixData):
                    dg_ij = dg_ij_src.mat('dg', self.components)
                elif isinstance(dg_ij_src, list):
                    dg_ij = np.array(dg_ij_src)
                elif isinstance(dg_ij_src, np.ndarray):
                    dg_ij = dg_ij_src
                elif isinstance(dg_ij_src, dict):
                    dg_ij = self.to_matrix_ij(
                        dg_ij_src,
                        symbol_delimiter=symbol_delimiter
                    )
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, dict, list of lists, or numpy array."
                    )

            else:
                # ! use constants based on the selected correlation
                tau_ij_cal_method = 2

                tau_required_sources = {
                    'M2': {'a_ij': a_ij_src, 'b_ij': b_ij_src, 'c_ij': c_ij_src, 'd_ij': d_ij_src},
                    'M3': {'a_ij': a_ij_src, 'b_ij': b_ij_src},
                    'M4': {'a_ij': a_ij_src, 'b_ij': b_ij_src, 'c_ij': c_ij_src},
                    'M5': {'a_ij': a_ij_src, 'b_ij': b_ij_src, 'c_ij': c_ij_src},
                }

                required_sources = tau_required_sources[tau_correlation]
                missing_sources = [
                    name for name, src in required_sources.items()
                    if not _has_source(src)
                ]

                if missing_sources:
                    raise ValueError(
                        f"Missing required source(s) for tau_correlation '{tau_correlation}': {', '.join(missing_sources)}."
                    )

                if _has_source(a_ij_src):
                    a_ij = _to_matrix(
                        a_ij_src, 'a', 'interaction parameter (a_ij)')
                if _has_source(b_ij_src):
                    b_ij = _to_matrix(
                        b_ij_src, 'b', 'interaction parameter (b_ij)')
                if _has_source(c_ij_src):
                    c_ij = _to_matrix(
                        c_ij_src, 'c', 'interaction parameter (c_ij)')
                if _has_source(d_ij_src):
                    d_ij = _to_matrix(
                        d_ij_src, 'd', 'interaction parameter (d_ij)')

            # SECTION: extract data
            # NOTE: α_ij, non-randomness parameter
            # check
            if alpha_ij_src is not None:
                if isinstance(alpha_ij_src, TableMatrixData):
                    alpha_ij = alpha_ij_src.mat('alpha', self.components)
                elif isinstance(alpha_ij_src, list):
                    alpha_ij = np.array(alpha_ij_src)
                elif isinstance(alpha_ij_src, np.ndarray):
                    alpha_ij = alpha_ij_src
                elif isinstance(alpha_ij_src, dict):
                    alpha_ij = self.to_matrix_ij(
                        alpha_ij_src,
                        symbol_delimiter=symbol_delimiter
                    )
                else:
                    raise ValueError(
                        "Invalid source for non-randomness parameter (α_ij). Must be TableMatrixData, dict, list of lists, or numpy array."
                    )
            else:
                # set default value
                alpha_ij = None

            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            # check
            if tau_ij_src is None:
                # ! tau_ij is None
                if T <= 0:
                    raise ValueError(
                        "temperature must be provided and greater than 0 K when calculating tau_ij from dg_ij or a/b/c/d."
                    )

                # ? check method
                if tau_ij_cal_method == 1:
                    # >>> method 1
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
                        dg_ij=dg_ij,
                        symbol_delimiter=symbol_delimiter
                    )
                elif tau_ij_cal_method == 2:
                    # >>> method 2
                    required_arrays = {
                        'M2': {'a_ij': a_ij, 'b_ij': b_ij, 'c_ij': c_ij, 'd_ij': d_ij},
                        'M3': {'a_ij': a_ij, 'b_ij': b_ij},
                        'M4': {'a_ij': a_ij, 'b_ij': b_ij, 'c_ij': c_ij},
                        'M5': {'a_ij': a_ij, 'b_ij': b_ij, 'c_ij': c_ij},
                    }[tau_correlation]

                    missing_arrays = [
                        name for name, value in required_arrays.items()
                        if value is None
                    ]

                    if missing_arrays:
                        raise ValueError(
                            f"Missing required coefficient matrix/matrices for tau_correlation '{tau_correlation}': {', '.join(missing_arrays)}."
                        )

                    if a_ij is not None:
                        if isinstance(a_ij, np.ndarray):
                            a_ij = a_ij.astype(float)
                        else:
                            raise ValueError("a_ij must be a numpy array")

                    if b_ij is not None:
                        if isinstance(b_ij, np.ndarray):
                            b_ij = b_ij.astype(float)
                        else:
                            raise ValueError("b_ij must be a numpy array")

                    if c_ij is not None:
                        if isinstance(c_ij, np.ndarray):
                            c_ij = c_ij.astype(float)
                        else:
                            raise ValueError("c_ij must be a numpy array")

                    if tau_correlation == 'M2':
                        a_ij_m2 = _require_matrix(a_ij, "a_ij")
                        b_ij_m2 = _require_matrix(b_ij, "b_ij")
                        c_ij_m2 = _require_matrix(c_ij, "c_ij")
                        d_ij_m2 = _require_matrix(d_ij, "d_ij")
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M2(
                            temperature=T,
                            a_ij=a_ij_m2,
                            b_ij=b_ij_m2,
                            c_ij=c_ij_m2,
                            d_ij=d_ij_m2,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_correlation == 'M3':
                        a_ij_m3 = _require_matrix(a_ij, "a_ij")
                        b_ij_m3 = _require_matrix(b_ij, "b_ij")
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M3(
                            temperature=T,
                            a_ij=a_ij_m3,
                            b_ij=b_ij_m3,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_correlation == 'M4':
                        a_ij_m4 = _require_matrix(a_ij, "a_ij")
                        b_ij_m4 = _require_matrix(b_ij, "b_ij")
                        c_ij_m4 = _require_matrix(c_ij, "c_ij")
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M4(
                            temperature=T,
                            a_ij=a_ij_m4,
                            b_ij=b_ij_m4,
                            c_ij=c_ij_m4,
                            symbol_delimiter=symbol_delimiter
                        )
                    elif tau_correlation == 'M5':
                        a_ij_m5 = _require_matrix(a_ij, "a_ij")
                        b_ij_m5 = _require_matrix(b_ij, "b_ij")
                        c_ij_m5 = _require_matrix(c_ij, "c_ij")
                        tau_ij, tau_ij_comp = self.cal_tau_ij_M5(
                            temperature=T,
                            a_ij=a_ij_m5,
                            b_ij=b_ij_m5,
                            c_ij=c_ij_m5,
                            symbol_delimiter=symbol_delimiter
                        )
                    else:
                        raise ValueError(
                            f"tau_correlation '{tau_correlation}' not supported!"
                        )
                else:
                    raise ValueError(
                        "tau_ij_cal_method not supported!"
                    )
            else:
                # ! check if tau_ij is provided
                # >>> method 0
                # check types
                if isinstance(tau_ij_src, TableMatrixData):
                    tau_ij = tau_ij_src.mat('tau', self.components)
                elif isinstance(tau_ij_src, list):
                    tau_ij = np.array(tau_ij_src)
                elif isinstance(tau_ij_src, np.ndarray):
                    tau_ij = tau_ij_src
                elif isinstance(tau_ij_src, dict):
                    tau_ij = self.to_matrix_ij(
                        tau_ij_src,
                        symbol_delimiter=symbol_delimiter
                    )
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (tau_ij). Must be TableMatrixData, dict, list of lists, or numpy array.")

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
