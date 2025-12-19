# import packages/modules
from typing import Dict, Optional, Literal, List, Any
# local
from .activity import (
    ActivityCore,
    NRTL,
    UNIQUAC,
    UNIFAC
)
from .docs import ThermoModelCore
from .utils import add_attributes
from .plugin import EQUATION_OF_STATE_MODELS, ACTIVITY_MODELS


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


@add_attributes(metadata=EQUATION_OF_STATE_MODELS)
def eos(**kwargs):
    '''
    Initialize equation of state library for fugacity calculation for single and multi-component systems
    (e.g., SRK, PR, RK)

    Parameters
    ----------
    **kwargs: dict
        Additional keyword arguments.

    Returns
    -------
    eosCoreC: object
        Instance of the selected equation of state class.
    '''
    try:
        # manager
        ThermoModelCore_ = ThermoModelCore()
        # return
        return ThermoModelCore_.init_eos(**kwargs)
    except Exception as e:
        raise Exception("Initialization failed!, ", e)


@add_attributes(metadata=ACTIVITY_MODELS)
def activity(
        components: List[str],
        model_name: Literal['NRTL', 'UNIQUAC', 'UNIFAC'],
        model_source: Optional[Dict[str, Any]] = None,
        **kwargs
) -> ActivityCore:
    '''
    Initialize activity calculation library

    Parameters
    ----------
    components: list
        List of component names to be used in the activity model, such as ['ethanol', 'butyl-methyl-ether'].
    model_name: str
        Name of the activity model to be used (e.g., 'NRTL', 'UNIQUAC', 'UNIFAC').
            1. `NRTL`: Non-Random Two-Liquid Model
            2. `UNIQUAC`: Universal Quasi-Chemical Model
            3. `UNIFAC`: Universal Functional Activity Coefficient Model
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
        ThermoModelCore_ = ThermoModelCore()
        # return
        return ThermoModelCore_.init_activity(
            components=components,
            model_name=model_name,
            model_source=model_source,
            **kwargs
        )
    except Exception as e:
        raise Exception("Initialization failed!, ", e)


@add_attributes(metadata=ACTIVITY_MODELS)
def activities(
        components: List[str],
        model_name: Literal['NRTL', 'UNIQUAC', 'UNIFAC'],
        model_source: Optional[Dict[str, Any]] = None,
        **kwargs
) -> NRTL | UNIQUAC | UNIFAC:
    '''
    Initialize activity calculation library

    Parameters
    ----------
    components: list
        List of component names to be used in the activity model, such as ['ethanol', 'butyl-methyl-ether'].
    model_name: str
        Name of the activity model to be used (e.g., 'NRTL', 'UNIQUAC', 'UNIFAC').
            1. `NRTL`: Non-Random Two-Liquid Model
            2. `UNIQUAC`: Universal Quasi-Chemical Model
            3. `UNIFAC`: Universal Functional Activity Coefficient Model
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
    NRTL | UNIQUAC | UNIFAC: object
        Instance of the selected activity model class.
    '''
    try:
        # init manager
        ThermoModelCore_ = ThermoModelCore()

        # return
        return ThermoModelCore_.init_activities(
            components=components,
            model_name=model_name,
            model_source=model_source,
            **kwargs
        )
    except Exception as e:
        raise Exception("Initialization failed!, ", e)
