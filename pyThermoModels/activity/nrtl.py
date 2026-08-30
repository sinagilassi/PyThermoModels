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
)
from pythermodb_settings.utils import create_mixture_id
# local
from ..utils import add_attributes
from ..plugin import ACTIVITY_MODELS
from ..utils.utility import TauCorrelation
from .nrtl_parameter_builder import NRTLParameterBuilder
from .component_parameter_mixin import ComponentParameterMixin

# NOTE: logger
logger = logging.getLogger(__name__)


class NRTL:
    """
    The NRTL (`Non-Random Two-Liquid`) model - a thermodynamic framework used to describe the behavior of mixtures,
    particularly in the context of phase equilibria and activity coefficients.

    The NRTL model relies on several key parameters to describe the interactions between components in a mixture. These parameters are:
    - Î”g_ij (interaction energy parameter): represents the interaction energy between two molecules [J/mol].
    - Î±_ij (non-randomness parameter): represents the non-randomness of the mixture [dimensionless].
    - Ï„_ij (binary interaction parameter): represents the interaction energy between two molecules of different components [dimensionless].

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

    # NOTE: mixture ids
    # ! default ids: Name and Formula
    _mixture_ids: Dict[str, str] = {}

    # mixture id
    _mixture_id: str = ""

    # NOTE: components ids
    _components_ids: Dict[str, List[str]] = {}

    def __init__(
        self,
        components: List[str],
        datasource: Optional[Dict] = None,
        equationsource: Optional[Dict] = None,
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
        datasource = {} if datasource is None else datasource
        equationsource = {} if equationsource is None else equationsource

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
        self._mixture_ids = {}
        self._mixture_id = ""
        self._components_ids = {}

        # NOTE: component configurations
        # components
        self.components = [component.strip() for component in components]

        # Get the number of components
        self.comp_num = len(self.components)

        # component index
        self.comp_idx = {
            self.components[i]: i for i in range(self.comp_num)
        }

        # SECTION: nrtl parameter builder
        self.nrtl_parameter_builder = NRTLParameterBuilder(
            components=self.components,
            comp_idx=self.comp_idx,
            datasource=self.datasource,
            equationsource=self.equationsource,
            **kwargs
        )
        # ! set input generator
        self.inputs_generator = self.nrtl_parameter_builder.inputs_generator
        # ! expose local-composition parameter helpers for compatibility
        self.cal_dg_ij_M1 = self.nrtl_parameter_builder.cal_dg_ij_M1
        self.cal_tau_ij_M1 = self.nrtl_parameter_builder.cal_tau_ij_M1
        self.cal_tau_ij_M2 = self.nrtl_parameter_builder.cal_tau_ij_M2
        self.cal_tau_ij_M3 = self.nrtl_parameter_builder.cal_tau_ij_M3
        self.cal_tau_ij_M4 = self.nrtl_parameter_builder.cal_tau_ij_M4
        self.cal_tau_ij_M5 = self.nrtl_parameter_builder.cal_tau_ij_M5

        # SECTION: component parameter mixin
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx
        )
        # ! set component access methods
        self.to_ij = self.component_parameter_mixin.to_ij
        self.to_dict_ij = self.component_parameter_mixin.to_dict_ij
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij
        self.to_dict_ij_ext = self.component_parameter_mixin.to_dict_ij_ext

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

    @property
    def mixture_ids(self) -> Dict[str, str]:
        '''
        Get the mixture ids.

        Returns
        -------
        mixture_ids: Dict[str, str]
            Dictionary of mixture ids
        '''
        return self._mixture_ids

    @mixture_ids.setter
    def mixture_ids(self, mixture_ids: Dict[str, str]) -> None:
        '''
        Set the mixture ids.

        Parameters
        ----------
        mixture_ids: Dict[str, str]
            Dictionary of mixture ids
        '''
        if not isinstance(mixture_ids, dict):
            raise TypeError("mixture_ids must be a dict")

        # reset
        self._mixture_ids = {}
        # set
        self._mixture_ids = mixture_ids

    @property
    def mixture_id(self) -> str:
        '''
        Get the mixture id.

        Returns
        -------
        mixture_id: str
            Mixture id
        '''
        return self._mixture_id

    @mixture_id.setter
    def mixture_id(self, mixture_id: str) -> None:
        '''
        Set the mixture id.

        Parameters
        ----------
        mixture_id: str
            Mixture id
        '''
        if not isinstance(mixture_id, str):
            raise TypeError("mixture_id must be a str")

        # set
        self._mixture_id = mixture_id

    @property
    def components_ids(self) -> Dict[str, List[str]]:
        '''
        Get the components ids.

        Returns
        -------
        components_ids: Dict[str, List[str]]
            Dictionary of components ids
        '''
        return self._components_ids

    @components_ids.setter
    def components_ids(self, components_ids: Dict[str, List[str]]) -> None:
        '''
        Set the components ids.

        Parameters
        ----------
        components_ids: Dict[str, List[str]]
            Dictionary of components ids
        '''
        if not isinstance(components_ids, dict):
            raise TypeError("components_ids must be a dict")

        # reset
        self._components_ids = {}
        # set
        self._components_ids = components_ids

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

    # SECTION: Calculate Gij matrix
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

    # SECTION: Check & build inputs
    def check_and_build_inputs(
        self,
        model_input: Dict,
        required_keys: List[str] = ['tau_ij', 'alpha_ij'],
        tau_correlation: TauCorrelation = "gibbs_energy",
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        return_all: bool = True,
        **kwargs
    ) -> Dict:
        """
        Check and build the required input parameters for the NRTL model.

        Parameters
        ----------
        model_input : dict
            Dictionary of model input values where keys are parameter names and values are their respective values.
                - `mole_fraction`: Dict[str, float]
                    Dictionary of mole fractions where keys are component names and values are their respective mole fractions.
                - `temperature`: List[str | float], Optional
                    List of temperatures in any units as [300, 'K'], it is automatically converted to Kelvin.
                - `tau_ij`: TableMatrixData | np.ndarray | Dict[str, float]
                    Interaction parameters (tau_ij) between component i and j.
                - `alpha_ij`: TableMatrixData | np.ndarray | Dict[str, float]
                    Non-randomness parameters (alpha_ij) between component i and j.
        tau_correlation : TauCorrelation
            Descriptive tau-correlation name. Default is `gibbs_energy`.
        symbol_delimiter : Literal["|", "_"]
            Delimiter for the component id. Default is "|".
        **kwargs : Optional
            Additional keyword arguments for the calculation.

        Returns
        -------
        model_input : dict
            Updated dictionary of model input values with required parameters added if they were missing.
        """
        try:
            # SECTION: check
            if not isinstance(model_input, dict):
                raise TypeError("model_input must be dict")

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
                # NOTE: check temperature
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

                # NOTE: call input generator
                # ! check and calculate the required keys
                # ! to eventually calculate tau_ij & alpha_ij
                inputs_ = self.inputs_generator(
                    temperature=model_input['temperature'],
                    tau_correlation=tau_correlation,
                    symbol_delimiter=symbol_delimiter,
                    mixture_ids=self.mixture_ids,
                    model_input=model_input,
                    **kwargs
                )

                # looping through the missed keys
                for key in missed_keys:
                    # key value
                    value_ = inputs_[key]

                    # check
                    if value_ is None:
                        # error
                        raise ValueError(
                            f"{key} is required in model_input"
                        )

                    # update the model_input
                    model_input[key] = value_

            # SECTION: get requested values
            res = {
                'mole_fraction': mole_fraction,
            }

            for key in required_keys:
                if key not in model_input:
                    raise ValueError(f"{key} is required in model_input")
                res[key] = model_input[key]

            # check return_all
            if return_all is False:
                # remove mole_fraction
                res.pop('mole_fraction', None)
                # >> res
                return res

            # return all
            return res

        except Exception as e:
            raise Exception(f"Error in check_and_build_inputs: {str(e)}")

    # SECTION: Calculate activity coefficients
    @add_attributes(metadata=ACTIVITY_MODELS['NRTL'])
    def cal(
        self,
        model_input: Dict,
        calculation_mode: Literal[
            'V1', 'V2'
        ] = 'V1',
        tau_correlation: TauCorrelation = "gibbs_energy",
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
        tau_correlation: TauCorrelation
            Descriptive tau-correlation name. Default is `gibbs_energy`.
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

        tau_ij can be calculated using constants or temperature-dependent correlations. The `tau_correlation` parameter specifies the method used to calculate tau_ij.

        - M1: tau_ij = dg_ij / (R * T)
        - M2: tau_ij = a_ij + b_ij / T + c_ij * log(T) + d_ij * T
        - M3: tau_ij = a_ij + b_ij / T
        - M4: tau_ij = a_ij + b_ij / T + c_ij / T^2
        - M5: tau_ij = a_ij + b_ij / T + c_ij * ln(T)

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
            # SECTION: check and build inputs
            inputs_src = self.check_and_build_inputs(
                model_input=model_input,
                tau_correlation=tau_correlation,
                symbol_delimiter=symbol_delimiter,
                return_all=True,
                **kwargs
            )
            # >> unpack
            mole_fraction = inputs_src['mole_fraction']
            tau_ij_data = inputs_src['tau_ij']
            alpha_ij_data = inputs_src['alpha_ij']

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
                tau_ij = np.asarray(tau_ij_data, dtype=float)
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
                alpha_ij = np.asarray(alpha_ij_data, dtype=float)
                # to dict
                alpha_ij_comp = self.to_dict_ij(
                    alpha_ij,
                    symbol_delimiter=symbol_delimiter
                )
            elif isinstance(alpha_ij_data, TableMatrixData):  # ! PyThermoDB
                # convert to numpy array and dict
                res_ = self.to_ij(
                    data=alpha_ij_data,
                    prop_symbol="alpha"
                )
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
            # >>> to np.array
            xi = np.asarray(xi, dtype=float)

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
                gE_RT += x[i] * numer / denom

            # SECTION: set result format
            res = {
                "property_name": "Excess Molar Gibbs Free Energy (G^E/RT)",
                "components": components,
                "mole_fraction": xi.tolist(),
                "mole_fraction_normalized": x.tolist(),
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
