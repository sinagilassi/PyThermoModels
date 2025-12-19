# import libs
import logging
from typing import Dict, List, Any, Optional, Tuple
import pycuc
# local
from .unifac1 import UNIFAC1
from ..utils import sanitize_mole_fractions

# NOTE: logger
logger = logging.getLogger(__name__)


class UNIFAC():
    '''
    UNIFAC (Universal Quasi-Chemical Functional-group Activity Coefficients) activity model class.

    Introduction
    ------------
    The UNIFAC model is a group-contribution method used to estimate activity coefficients in non-ideal liquid mixtures. It divides molecules into functional groups and calculates the activity coefficients based on group interactions, allowing for the prediction of phase equilibria in multi-component systems. This approach is widely used in chemical engineering for the design and analysis of separation processes, such as distillation and extraction, where accurate modeling of non-ideal behavior is essential.

    '''
    # universal gas constant [J/mol/K]
    R_CONST = 8.314
    # constant
    Z = 10.0

    # NOTE: components ids
    _components_ids: Dict[str, List[str]] = {}

    # NOTE: group data
    _group_data: Dict[str, Dict[str, Any]] = {}
    # NOTE: interaction data
    _interaction_data: Dict[str, Dict[str, Any]] = {}

    # SECTION: Initialization model
    _model: Optional[UNIFAC1] = None

    def __init__(
        self,
        components: List[str],
        datasource: Dict = {},
        equationsource: Dict = {},
        **kwargs
    ):
        '''
        Initialize the UNIFAC activity model.

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
        self.components = [components.strip() for components in components]

        # SECTION
        # Get the number of components
        self.comp_num = len(components)
        # idx
        self.comp_idx = {components[i]: i for i in range(self.comp_num)}

    def __str__(self) -> str:
        '''
        String representation of the UNIFAC model.

        Returns
        -------
        str
            String representation of the UNIFAC model.
        '''
        comp_str = ', '.join(self.components)
        return (
            f"UNIFAC Activity Model\n"
            f"----------------------\n"
            f"Components: {comp_str}\n"
            f"Number of components: {self.comp_num}\n"
            f"Data source: {type(self.datasource).__name__}\n"
            f"Equation source: {type(self.equationsource).__name__}\n"
        )

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

    @property
    def group_data(self) -> Dict[str, Dict[str, Any]]:
        '''
        Get the group data.

        Returns
        -------
        group_data: Dict[str, Dict[str, Any]]
            Dictionary of group data
        '''
        return self._group_data

    @group_data.setter
    def group_data(self, group_data: Dict[str, Dict[str, Any]]) -> None:
        '''
        Set the group data.

        Parameters
        ----------
        group_data: Dict[str, Dict[str, Any]]
            Dictionary of group data
        '''
        if not isinstance(group_data, dict):
            raise TypeError("group_data must be a dict")

        # reset
        self._group_data = {}
        # set
        self._group_data = group_data

    @property
    def interaction_data(self) -> Dict[str, Dict[str, Any]]:
        '''
        Get the interaction data.

        Returns
        -------
        interaction_data: Dict[str, Dict[str, Any]]
            Dictionary of interaction data
        '''
        return self._interaction_data

    @interaction_data.setter
    def interaction_data(self, interaction_data: Dict[str, Dict[str, Any]]) -> None:
        '''
        Set the interaction data.

        Parameters
        ----------
        interaction_data: Dict[str, Dict[str, Any]]
            Dictionary of interaction data
        '''
        if not isinstance(interaction_data, dict):
            raise TypeError("interaction_data must be a dict")

        # reset
        self._interaction_data = {}
        # set
        self._interaction_data = interaction_data

    def load_data(
            self,
            group_data: Dict[str, Dict[str, Any]],
            interaction_data: Dict[str, Dict[str, Any]],
            **kwargs: Dict[str, Any]
    ) -> None:
        '''
        Load the UNIFAC model data, including group properties and interaction parameters.

        Parameters
        ----------
        - group_data : Dict[str, Dict[str, Any]]
            Dictionary containing group properties (R, Q, main-group).
        - interaction_data : Dict[str, Dict[str, Any]]
            Dictionary containing interaction parameters between main groups.
        - **kwargs : dict
            Additional keyword arguments for UNIFAC1 model initialization.
            - eps : float, optional
                Small number to avoid division by zero (default is 1e-30).
            - z : float, optional
                Coordination number (default is 10.0).
        '''
        try:
            # SECTION: Config group & interaction data
            self.group_data = group_data
            self.interaction_data = interaction_data

            # SECTION: Initialize UNIFAC1 model
            # NOTE: reset
            self._model = None

            # NOTE: initialize
            self._model = UNIFAC1(
                group_data=self.group_data,
                interaction_data=self.interaction_data,
                **kwargs
            )

        except Exception as e:
            raise Exception(f"Error initializing UNIFAC model: {e}") from e

    def set_component_groups(
            self,
            component_groups: Dict[str, Dict[str, float | int]]
    ):
        '''
        Set the component group contributions for the UNIFAC model.

        Parameters
        ----------
        component_groups : Dict[str, Dict[str, float|int]]
            Dictionary mapping component names to their group contributions.
            Example:
            {
                "Component1": {"Group1": 2, "Group2": 1},
                "Component2": {"Group1": 1, "Group3": 3},
                ...
            }
        '''
        try:
            # NOTE: Check if model is initialized
            if self._model is None:
                raise Exception(
                    "UNIFAC model data not loaded. Call load_data() first.")

            # SECTION: Validate component_groups based on defined components
            # NOTE: init component data
            component_groups_validated: List[Dict[str, float | int]] = []

            # iterate components
            for comp in self.components:
                if comp not in component_groups:
                    raise KeyError(
                        f"Component '{comp}' not found in component_groups."
                    )

                # set
                component_groups_validated.append(component_groups[comp])

            # NOTE: Set component groups in the UNIFAC1 model
            self._model.component_group_data = component_groups_validated

            # SECTION: initialize calculation parameters
            self._model.initialize_calc()

        except Exception as e:
            raise Exception(f"Error setting component groups: {e}") from e

    def cal(
        self,
        model_inputs: Dict[str, Any],
        **kwargs
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        '''
        Calculate activity coefficients using the UNIFAC model.

        Parameters
        ----------
        model_inputs : Dict[str, Any]
            Dictionary containing model input parameters.
            Required keys:
            - mole_fraction : Dict[str, float]
                Mole fractions of each component in the mixture.
            - temperature : List[float, str]
                Temperature as a list containing the value and unit (e.g., [298.15, 'K']).
        **kwargs : dict
            Additional keyword arguments for calculation.
            - eps : float, optional
                Small number to avoid division by zero (default is 1e-30).
            - z : float, optional
                Coordination number (default is 10.0).
            - x_eps : float, optional
                Small value to floor mole fractions to avoid log(0) (default is 1e-30).

        Returns
        -------
        Dict[str, float]
            Dictionary of activity coefficients for each component.

        Notes
        -----
        The model_inputs dictionary must contain the following keys:
        - 'mole_fraction': A dictionary mapping component names to their mole fractions.
        - 'temperature': A list containing the temperature value and its unit.

        Examples
        --------
        ```python
        model_inputs = {
            'mole_fraction': {
                'Component1': 0.5,
                'Component2': 0.5
            },
            'temperature': [298.15, 'K']
        }

        activity_coefficients = unifac_model.calc(model_inputs)
        ```
        '''
        try:
            # SECTION: Prepare model inputs
            # NOTE: mole_fraction_eps
            x_eps = kwargs.get('x_eps', 1e-30)
            # NOTE: eps
            eps = kwargs.get('eps', 1e-30)
            # NOTE: z
            z = kwargs.get('z', 10.0)

            # SECTION: Validate model inputs
            required_keys = ['mole_fraction', 'temperature']

            for key in required_keys:
                if key not in model_inputs:
                    raise KeyError(f"Missing required model input: '{key}'")

            # NOTE: set temperature
            temperature_src = model_inputs['temperature']
            temperature_value = temperature_src[0]
            temperature_unit = temperature_src[1]

            # >> convert temperature to Kelvin if necessary
            temperature_K = pycuc.convert_from_to(
                value=temperature_value,
                from_unit=temperature_unit,
                to_unit='K'
            )

            # NOTE: set mole fractions
            mole_fraction_src = model_inputs['mole_fraction']
            mole_fraction = [
                float(mole_fraction_src[comp]) for comp in self.components
            ]

            # >> sanitize mole fractions
            mole_fraction = sanitize_mole_fractions(
                x=mole_fraction,
                eps=x_eps
            )
            # >> check mole fractions sum to 1
            if not mole_fraction or abs(sum(mole_fraction) - 1.0) > 1e-6:
                raise ValueError(
                    "Mole fractions must sum to 1 after sanitization.")

            # SECTION: Perform calculation using UNIFAC1 model
            # NOTE: Check if model is initialized
            if self._model is None:
                raise Exception(
                    "UNIFAC model data not loaded. Call load_data() first.")

            # NOTE: Set eps and z in the UNIFAC1 model
            self._model.eps = eps
            self._model.z = z

            # NOTE: calculate activity coefficients
            res = self._model.get_activity_coefficients(
                T=temperature_K,
                x=mole_fraction,
                comp=self.components
            )

            # SECTION: Return results
            return res
        except Exception as e:
            raise Exception(
                f"Error calculating activity coefficients: {e}") from e

    def group_data_template(self) -> Dict[str, Dict[str, Any]]:
        '''
        Get a template for UNIFAC group data.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Template dictionary for UNIFAC group data.
        '''
        return {
            "1": {
                'main-group': '1',
                'sub-group': '1',
                'group-id': 'CH3',
                'R': '0.9011',
                'Q': '0.848'
            },
        }

    def interaction_data_template(self) -> Dict[str, Dict[str, Any]]:
        '''
        Get a template for UNIFAC group-group interaction data.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Template dictionary for UNIFAC interaction data.
        '''
        return {
            "1": {
                '': 1.0, '1': 0.0, '2': 86.02, '3': 61.13, '4': 76.5, '5': 986.5, '6': 697.2, '7': 1318.0, '8':
                1333.0, '9': 476.4, '10': 677.0, '11': 232.1, '12': 507.0, '13': 251.5, '14': 391.5, '15': 255.7, '16': 206.6,
                '17': 920.7, '18': 287.8, '19': 597.0, '20': 663.5, '21': 35.93, '22': 53.76, '23': 24.9, '24': 104.3, '25': 11.44,
                '26': 661.5, '27': 543.0, '28': 153.6, '29': 184.4, '30': 354.6, '31': 3025.0, '32': 335.8, '33': 479.5, '34':
                298.9, '35': 526.5, '36': 689.0, '37': -4.189, '38': 125.8, '39': 485.3, '40': -2.859, '41': 387.1, '42': -450.4,
                '43': 252.7, '44': 220.3, '45': -5.869, '46': 390.9, '47': 553.3, '48': 187.0, '49': 216.1, '50': 92.99
            },
        }
