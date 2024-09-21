# THERMO LIB
# ---------------

# import packages/modules

# local
from .thermo import *
from ..configs import R_CONST


class ThermoLib():

    def __init__(self) -> None:
        pass

    def cal_molar_volume(self, P, T, Z):
        '''
        Calculate molar volume

        Parameters
        ----------
        P : float
            pressure [Pa]
        T : float
            temperature [K]
        Z : float
            compressibility factor

        Returns
        -------
        Vm : float
            molar volume [m^3/mol]
        '''
        try:
            # molar-volume [m^3/mol]
            Vm = Z * ((R_CONST * T) / P)
            return Vm
        except Exception as e:
            raise Exception("Molar-volume calculation failed!, ", e)
