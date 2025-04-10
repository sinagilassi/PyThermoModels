# import packages/modules

# local
from .configs import __description__, __version__
from .docs import ThermoModelCore, ThermoLib


def init():
    '''
    Initialize app

    Returns
    -------
    ThermoModelCore: object
        ThermoModelCore object for the calculation of fugacity and activity coefficient
    '''
    try:
        # return
        return ThermoModelCore()
    except Exception as e:
        raise Exception("Calculating the Fugacity failed!, ", e)


def thermo_lib():
    '''
    Initialize thermodynamic calculation library

    Parameters
    ----------
    None

    Returns
    -------
    fugacity: list
        fugacity for gas/liquid/solid phase
    '''
    try:
        # manager
        ThermoLibC = ThermoLib()
        # return
        return ThermoLibC
    except Exception as e:
        raise Exception("Initialization failed!, ", e)
