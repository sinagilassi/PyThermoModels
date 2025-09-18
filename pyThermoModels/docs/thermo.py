# THERMODYNAMIC RELATIONS
# ------------------------
# import packages/modules
from math import pow
# internals
from ..configs import Tref, R_CONST


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


def calMolarVolume(P: float, T: float, Z: float, R: float = R_CONST):
    '''
    Calculate molar-volume [m^3/mol]

    Parameters
    ----------
    P: float
        pressure [Pa]
    T: float
        temperature [K]
    Z: float
        compressibility factor [-]
    R: float, optional
        universal gas constant [J/mol.K]. The default is R_CONST.

    Returns
    -------
    Vm: float
        molar volume [m^3/mol]
    '''
    try:
        Vm = Z * ((R * T) / P)
        return Vm
    except Exception as e:
        raise Exception("Error in calMolarVolume: ", e)


def RackettEquation(Vc, Zc, Tr):
    '''
    estimation of saturated liquid volume

    args:
        Vc: critical molar-volume
        Zc: critical compressibility factor
        Tr: reduced temperature
    '''
    _ZcPower = pow(1-Tr, 0.2857)
    return Vc*pow(Zc, _ZcPower)


def ModifiedRackettEquation(T, Pc, Tc, w):
    '''
    estimation of saturated liquid volume

    args:
        T: temperature [K]
        Pc: critical pressure [Pa]
        Tc: critical temperature [K]
        w: acentric factor
    '''
    _c0 = R_CONST*Tc/Pc
    ZRA = 0.2956 - 0.08775*w
    Tr = T/Tc
    _c1 = 1 + pow(1-Tr, 2/7)
    _c2 = _c0*pow(ZRA, _c1)

    return _c2


def SetPhase(state):
    '''
    set phase
    '''
    # define phase
    _phaseSelection = {
        "l": "liquid",
        "g": "gas",
        "s": "solid"
    }

    return _phaseSelection.get(state)
