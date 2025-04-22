# import packages/modules
import numpy as np
from math import pow, sqrt, exp, log
# local
from ..configs import R_CONST


class ActivityCore:

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

    def activity_cal(self, mole_fraction: list, params: dict):
        pass
