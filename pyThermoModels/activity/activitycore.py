# import packages/modules
from typing import Union, Optional, Dict, List
from math import log
# local
from .nrtl import NRTL
from .uniquac import UNIQUAC


class ActivityCore:
    '''
    ActivityCore class for calculating activity coefficients using different models.
    '''
    # NOTE: activity model
    __nrtl: Optional[NRTL] = None
    __uniquac: Optional[UNIQUAC] = None

    # NOTE: mixture id
    _mixture_id: Optional[str] = None

    def __init__(
        self,
        datasource: dict,
        equationsource: dict,
        components: List[str],
        **kwargs
    ):
        '''
        Initialize the activity core class.

        Parameters
        ----------
        datasource : dict
            Data source object containing thermodynamic data.
        equationsource : dict
            Equation source object containing thermodynamic equations.
        components : list
            List of component names.
        **kwargs : dict
            Additional keyword arguments.
            - mixture_id : str, optional
                Unique identifier for the mixture. Default is None.
        '''
        # NOTE: data/equations
        self.datasource = datasource
        self.equationsource = equationsource
        # components
        self.components = components

        # NOTE: kwargs
        self._mixture_id: Optional[str] = kwargs.get('mixture_id', None)

        # SECTION: init activity models
        # ! nrtl
        self.__nrtl = NRTL(
            components=self.components,
            datasource=self.datasource,
            equationsource=self.equationsource,
            mixture_id=self._mixture_id,
        )
        # ! uniquac
        self.__uniquac = UNIQUAC(
            components=self.components,
            datasource=self.datasource,
            equationsource=self.equationsource,
            mixture_id=self._mixture_id,
        )

    @property
    def mixture_id(self) -> str:
        '''
        Get the unique mixture ID.

        Returns
        -------
        str
            Unique mixture ID.
        '''
        try:
            if self._mixture_id is None:
                return "unknown"
            return self._mixture_id
        except Exception as e:
            raise Exception(f"Error in mixture_id: {e}") from e

    @property
    def nrtl(self):
        '''
        Initialize the NRTL activity model.

        Returns
        -------
        NRTL
            Instance of the NRTL activity model class.
        '''
        try:
            # NOTE: NRTL
            # check if nrtl is None
            if self.__nrtl is None:
                # err
                raise ValueError("NRTL model not initialized.")
            return self.__nrtl
        except Exception as e:
            raise Exception(f"Error in NRTL: {e}") from e

    @property
    def uniquac(self):
        '''
        Initialize the UNIQUAC activity model.

        Returns
        -------
        UNIQUAC
            Instance of the UNIQUAC activity model class.
        '''
        try:
            # NOTE: UNIQUAC
            # check if uniquac is None
            if self.__uniquac is None:
                # err
                raise ValueError("UNIQUAC model not initialized.")
            return self.__uniquac
        except Exception as e:
            raise Exception(f"Error in UNIQUAC: {e}") from e

    def select(self, model_name: str) -> Union[NRTL, UNIQUAC]:
        '''
        Select the activity model based on the model name.

        Parameters
        ----------
        model_name : str
            Name of the activity model to be used.

        Returns
        -------
        NRTL | UNIQUAC
            Instance of the selected activity model class.
        '''
        try:
            if model_name == 'NRTL':
                return NRTL(
                    self.components, self.datasource, self.equationsource)
            elif model_name == 'UNIQUAC':
                return UNIQUAC(
                    self.components, self.datasource, self.equationsource)
            else:
                raise ValueError(f"Model {model_name} not supported.")
        except Exception as e:
            raise Exception(f"Error in activity_cal: {e}") from e

    def general_excess_molar_gibbs_free_energy(
        self,
        mole_fraction: Union[
            Dict[str, float],
            List[float]
        ],
        activity_coefficients: Union[
            Dict[str, float],
            List[float]
        ],
        message: Optional[str] = None,
    ):
        """
        Calculate the general excess molar Gibbs free energy which is not based on any specific activity model.

        (G^E / RT) = x1 * ln(gamma1) + x2 * ln(gamma2) + ... + xn * ln(gamman)

        where:
        - G^E: Excess molar Gibbs free energy
        - R: Universal gas constant
        - T: Temperature in Kelvin
        - x: Mole fraction of the component
        - gamma: Activity coefficient of the component

        Parameters
        ----------
        mole_fraction : dict or list
            Mole fractions of the components in the mixture.
        activity_coefficients : dict or list
            Activity coefficients of the components in the mixture.
        message : str, optional
            Optional message to include in the result. Default is None.

        Returns
        -------
        dict
            A dictionary containing the following keys:
            - 'property_name': Name of the property calculated.
            - 'components': List of component names.
            - 'mole_fraction': Mole fractions of the components.
            - 'value': Calculated excess molar Gibbs free energy (G^E / RT).
            - 'unit': Unit of the calculated value (dimensionless).
            - 'activity_coefficients': Activity coefficients of the components.
            - 'symbol': Symbol for the property.
            - 'message': Optional message.
        """
        try:
            # NOTE: components
            components = self.components
            components_str = ', '.join(components)

            # SECTION: check input types
            # check if mole_fraction and activity_coefficients are of the same type
            if (
                isinstance(mole_fraction, dict) and
                    isinstance(activity_coefficients, dict)
            ):
                # NOTE: convert to list based on component order
                components = list(mole_fraction.keys())
                x_i = [
                    mole_fraction[component] for component in components
                ]
                AcCo_i = [
                    activity_coefficients[component] for component in components
                ]
            elif (
                isinstance(mole_fraction, list) and
                isinstance(activity_coefficients, dict)
            ):
                # NOTE: convert to list based on component order
                components = list(activity_coefficients.keys())
                x_i = [
                    mole_fraction[i] for i in range(len(components))
                ]
                AcCo_i = [
                    activity_coefficients[component] for component in components]
            elif (
                isinstance(mole_fraction, dict) and
                isinstance(activity_coefficients, list)
            ):
                # NOTE: convert to list based on component order
                components = list(mole_fraction.keys())
                AcCo_i = [
                    activity_coefficients[i] for i in range(len(components))
                ]
                x_i = [
                    mole_fraction[component] for component in components
                ]
            elif (
                isinstance(mole_fraction, list) and
                isinstance(activity_coefficients, list)
            ):
                # NOTE: convert to list based on component order
                x_i = mole_fraction
                AcCo_i = activity_coefficients
            else:
                raise TypeError(
                    "mole_fraction and activity_coefficients must be of the same type (dict or list).")

            # SECTION: calculation
            # calculate excess molar Gibbs free energy
            ExMoGiEn = []

            for x, gamma in zip(x_i, AcCo_i):
                if x < 0 or gamma <= 0:
                    raise ValueError(
                        "Mole fraction must be non-negative and activity coefficients must be positive.")

                cal_ = x * log(gamma)
                ExMoGiEn.append(cal_)

            # set
            gE_RT = float(sum(ExMoGiEn))

            # NOTE: set message
            message = f'Excess Gibbs Free Energy for {components_str}' if message is None else message

            # NOTE: result
            res = {
                "property_name": "Excess Molar Gibbs Free Energy (G^E/RT)",
                'components': components,
                'mole_fraction': mole_fraction,
                'value': gE_RT,
                'unit': 'dimensionless',
                'activity_coefficients': activity_coefficients,
                "symbol": "ExMoGiFrEn",
                'message': message
            }
            # return result
            return res
        except Exception as e:
            raise Exception(
                f"Error in general_excess_molar_gibbs_free_energy: {e}") from e
