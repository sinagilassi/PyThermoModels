# import packages/modules
import numpy as np
from math import pow, sqrt, exp, log
from typing import Union, Optional, Any
# local
from ..configs import R_CONST
from .nrtl import NRTL
from .uniquac import UNIQUAC


class ActivityCore:
    '''
    ActivityCore class for calculating activity coefficients using different models.
    '''

    # NOTE: activity model
    __nrtl: Optional[NRTL] = None
    __uniquac: Optional[UNIQUAC] = None

    def __init__(self, datasource, equationsource, components, **kwargs):
        '''
        Initialize the activity core class.

        Parameters
        ----------
        datasource : object
            Data source object containing thermodynamic data.
        equationsource : object
            Equation source object containing thermodynamic equations.
        components : list
            List of component names.
        **kwargs : dict
            Additional keyword arguments.
        '''
        # data/equations
        self.datasource = datasource
        self.equationsource = equationsource
        # components
        self.components = components

        # SECTION: init activity models
        # nrtl
        self.__nrtl = NRTL(
            self.components, self.datasource, self.equationsource)
        # uniquac
        self.__uniquac = UNIQUAC(
            self.components, self.datasource, self.equationsource)

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
