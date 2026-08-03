# import libs
import logging
import numpy as np
import json
import yaml
from typing import (
    List,
    Dict,
    Tuple,
    Any,
    Literal,
    Optional
)
from pyThermoDB import (
    TableMatrixData,
)
# local
from ..plugin import ACTIVITY_MODELS
from ..utils import add_attributes
from ..utils.utility import TauCorrelation
from .uniquac_parameter_builder import UNIQUACParameterBuilder
from .component_parameter_mixin import ComponentParameterMixin

# NOTE: logger
logger = logging.getLogger(__name__)


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
        datasource: Dict = {},
        equationsource: Dict = {},
        **kwargs
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
        **kwargs: dict
            Additional keyword arguments

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
        # SECTION: inputs
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
        self.components = [component.strip() for component in components]

        # SECTION
        # Get the number of components
        self.comp_num = len(self.components)
        # idx
        self.comp_idx = {
            self.components[i]: i for i in range(self.comp_num)
        }

        # SECTION: uniquac parameter builder
        self.uniquac_parameter_builder = UNIQUACParameterBuilder(
            components=self.components,
            comp_idx=self.comp_idx,
            datasource=self.datasource,
            equationsource=self.equationsource,
            **kwargs
        )
        # ! set input generator
        self.inputs_generator = self.uniquac_parameter_builder.inputs_generator

        # SECTION: component parameter mixin
        self.component_parameter_mixin = ComponentParameterMixin(
            components=self.components,
            comp_idx=self.comp_idx
        )
        # ! set component access methods
        self.to_ij = self.component_parameter_mixin.to_ij
        self.to_i = self.component_parameter_mixin.to_i
        self.to_dict_ij = self.component_parameter_mixin.to_dict_ij
        self.to_dict_i = self.component_parameter_mixin.to_dict_i
        self.to_matrix_ij = self.component_parameter_mixin.to_matrix_ij
        self.to_dict_ij_ext = self.component_parameter_mixin.to_dict_ij_ext

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
            - { mole_fraction: { ethanol: 0.4, butyl-methyl-ether: 0.6 },
                temperature: [323.15, 'K'], tau_ij: [[],[]], r_i: [[],[]] }

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
        ij_data : TableMatrixData | np.ndarray | Dict[str,
            float] | List[List[float]]
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

    # SECTION: Check & build inputs
    def check_and_build_inputs(
        self,
        model_input: Dict,
        required_keys: List[str] = ['tau_ij', 'r_i', 'q_i'],
        tau_correlation: TauCorrelation = "extended_temperature",
        symbol_delimiter: Literal[
            "|", "_"
        ] = "|",
        return_all: bool = True,
        **kwargs
    ) -> Dict:
        """
        Check and build the required input parameters for the UNIQUAC model.

        Parameters
        ----------
        model_input : Dict
            Runtime model input. Must include `mole_fraction`; may include
            `tau_ij`, `r_i`, `q_i`, `temperature`, `mixture_ids_dict`, and
            `components_ids_dict`.
        required_keys : List[str]
            Required UNIQUAC keys to validate or generate.
        tau_correlation : TauCorrelation
            Descriptive tau-correlation name used when `tau_ij` must be
            generated. Defaults to `extended_temperature`.
        symbol_delimiter : Literal["|", "_"]
            Delimiter used for component-pair dictionary keys.
        return_all : bool
            If False, omit `mole_fraction` from the returned dictionary.
        **kwargs : dict
            Additional values forwarded to the parameter builder.

        Returns
        -------
        Dict
            Validated and generated UNIQUAC inputs.
        """
        try:
            # SECTION: validate model input
            if not isinstance(model_input, dict):
                raise TypeError("model_input must be dict")

            if 'mole_fraction' not in model_input:
                raise KeyError("mole_fraction is required in model_input")

            mole_fraction = model_input['mole_fraction']

            # SECTION: identify missing required parameters
            missed_keys = [
                key for key in required_keys if key not in model_input
            ]

            if len(missed_keys) > 0:
                # SECTION: validate optional temperature before generation
                temperature = model_input.get('temperature', None)

                if temperature is not None:
                    if not isinstance(temperature, list):
                        raise TypeError("temperature must be list")
                    if len(temperature) == 0:
                        raise ValueError("temperature list is empty")
                    if not all(
                        isinstance(temp, (int, float, str))
                        for temp in temperature
                    ):
                        raise TypeError(
                            "temperature list must be int, float, or str"
                        )

                # SECTION: build missing parameters from datasource/model_input
                inputs_ = self.inputs_generator(
                    temperature=temperature,
                    tau_correlation=tau_correlation,
                    symbol_delimiter=symbol_delimiter,
                    mixture_ids=(
                        self.mixture_ids or
                        model_input.get('mixture_ids_dict', None)
                    ),
                    components_ids=(
                        self.components_ids or
                        model_input.get('components_ids_dict', None)
                    ),
                    model_input=model_input,
                    **kwargs
                )

                # NOTE: generated values are written back for downstream use
                for key in missed_keys:
                    value_ = inputs_[key]
                    if value_ is None:
                        raise ValueError(f"{key} is required in model_input")
                    model_input[key] = value_

            # SECTION: package validated inputs
            res = {
                'mole_fraction': mole_fraction,
                'tau_ij': model_input['tau_ij'],
                'r_i': model_input['r_i'],
                'q_i': model_input['q_i'],
            }

            if return_all is False:
                # NOTE: utility callers can request only parameter matrices
                res.pop('mole_fraction', None)

            return res
        except Exception as e:
            raise Exception(f"Error in check_and_build_inputs: {str(e)}")

    @add_attributes(metadata=ACTIVITY_MODELS['UNIQUAC'])
    def cal(
        self,
        model_input: Dict,
        Z: Optional[float | int] = None,
        calculation_mode: Literal['V1'] = 'V1',
        tau_correlation: TauCorrelation = "extended_temperature",
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
        calculation_mode : Literal['V1']
            Mode of calculation. If 'V1', use the first version of the UNIQUAC model.
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
            inputs_src = self.check_and_build_inputs(
                model_input=model_input,
                tau_correlation=tau_correlation,
                symbol_delimiter=symbol_delimiter,
                return_all=True,
                **kwargs
            )

            mole_fraction = inputs_src['mole_fraction']
            tau_ij = inputs_src['tau_ij']
            r_i = inputs_src['r_i']
            q_i = inputs_src['q_i']

            # check r_i and q_i
            if r_i is None:
                # log
                logger.error("q_i is not provided in model_input")
                raise ValueError("r_i is required in model_input")

            if q_i is None:
                # log
                logger.error("q_i is not provided in model_input")
                raise ValueError("q_i is required in model_input")

            # SECTION: calculate activity coefficients
            return self.__calculate_activity_coefficients(
                mole_fraction=mole_fraction,
                tau_ij_data=tau_ij,
                r_i_data=r_i,
                q_i_data=q_i,
                Z=Z,
                calculation_mode=calculation_mode,
                symbol_delimiter=symbol_delimiter,
                message=message
            )
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
                    tau_ij_data,
                    symbol_delimiter=symbol_delimiter
                )
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
                    Z=Z
                )
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
            gamma_comb_ij = (
                np.log(phi_i/xi) + 1 - (phi_i/xi) - (Z/2)*q_i *
                (np.log(phi_i/teta_i)+1-(phi_i/teta_i))
            )

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
        mole_fraction: Optional[
            Dict[str, float]
        ] = None,
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

