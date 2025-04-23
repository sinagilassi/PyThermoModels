# import packages/modules
from typing import Dict, Optional, Literal, List
# local
from .configs import __description__, __version__
from .docs import ThermoModelCore, ThermoLib


def init():
    '''
    Initialize PyThermoModel app

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


def eos(
        components: List[str],
        model_name: Literal['SRK', 'PR', 'RK', 'vdW'],
        model_source: Optional[Dict] = None,
        **kwargs):
    '''
    Initialize equation of state library for fugacity calculation for single and multi-component systems
    (e.g., SRK, PR, RK)

    Parameters
    ----------
    components: list
        List of component names to be used in the fugacity model, such as ['ethanol', 'butyl-methyl-ether'].
    model_name: str
        Name of the fugacity model to be used (e.g., 'SRK', 'PR', 'RK').
            1. `SRK`: Soave-Redlich-Kwong Equation of State
            2. `PR`: Peng-Robinson Equation of State
            3. `RK`: Redlich-Kwong Equation of State
            4. `vdW`: van der Waals Equation of State
    model_source: dict, optional
        Dictionary containing the source of the fugacity model data.
        If None, default values will be used.
            - `datasource`: dict
                Dictionary containing the data source for the fugacity model.
            - `equationsource`: dict
                Dictionary containing the equation source for the fugacity model.
    **kwargs: dict
        Additional keyword arguments.

    Returns
    -------
    SRK | PR | RK | PCSAFT
        Instance of the selected fugacity model class.
    '''
    try:
        # manager
        FugacityCoreC = ThermoModelCore()
        # return
        return FugacityCoreC.init_fugacity(
            components=components,
            model_name=model_name,
            model_source=model_source,
            **kwargs)
    except Exception as e:
        raise Exception("Initialization failed!, ", e)


def activity(
        components: List[str],
        model_name: Literal['NRTL', 'UNIQUAC'],
        model_source: Optional[Dict] = None,
        **kwargs):
    '''
    Initialize activity calculation library

    Parameters
    ----------
    components: list
        List of component names to be used in the activity model, such as ['ethanol', 'butyl-methyl-ether'].
    model_name: str
        Name of the activity model to be used (e.g., 'NRTL', 'UNIQUAC').
            1. NRTL: Non-Random Two-Liquid Model
            2. UNIQUAC: Universal Quasi-Chemical Model
    model_source: dict, optional
        Dictionary containing the source of the activity model data.
        If None, default values will be used.
            - `datasource`: dict
                Dictionary containing the data source for the activity model.
            - `equationsource`: dict
                Dictionary containing the equation source for the activity model.
    **kwargs: dict
        Additional keyword arguments.

    Returns
    -------
    ActivityCoreC: object
        Instance of the selected activity model class.
    '''
    try:
        # manager
        ActivityCoreC = ThermoModelCore()
        # return
        return ActivityCoreC.init_activity(
            components=components,
            model_name=model_name,
            model_source=model_source,
            **kwargs)
    except Exception as e:
        raise Exception("Initialization failed!, ", e)


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
