# UTILITY FUNCTIONS
# ------------------

# import packages/modules
from typing import Literal, TypeAlias
import numpy as np
# internal
from ..configs import ROUND_FUN_ACCURACY
from ..configs import (
    PENG_ROBINSON,
    SOAVE_REDLICH_KWONG,
    VAN_DER_WAALS,
    REDLICH_KWONG
)


def roundNum(value, ACCURACY=ROUND_FUN_ACCURACY):
    '''
    Round a number, set decimal digit

    Parameters
    ----------
    value : float
        value to round
    ACCURACY : int
        decimal digit

    Returns
    -------
    value : float
        rounded value
    '''
    return np.round(value, ACCURACY)


def removeDuplicatesList(value):
    '''
    Remove duplicates from a list

    Parameters
    ----------
    value : list
        list to remove duplicates

    Returns
    -------
    value : list
        list without duplicates
    '''
    return list(dict.fromkeys(value))


def eos_model_name(
    model_name: str
) -> Literal['PR', 'SRK', 'RK', 'vdW']:
    '''
    Sets eos model name

    Parameters
    ----------
    model_name : str
        name of eos model

    Returns
    -------
    model_name_set : str
        name of eos model
    '''
    try:
        # check
        if model_name is None:
            raise Exception('Empty equation of state name!')
        else:
            # PR
            if (
                model_name == 'PENG_ROBINSON' or
                    model_name == 'PENG-ROBINSON'
            ):
                model_name_set = PENG_ROBINSON

            elif model_name == 'PR':
                model_name_set = PENG_ROBINSON

            # SRK
            elif (
                model_name == 'SOAVE_REDLICH_KWONG' or
                model_name == 'SOAVE-REDLICH-KWONG'
            ):
                model_name_set = SOAVE_REDLICH_KWONG

            elif model_name == 'SRK':
                model_name_set = SOAVE_REDLICH_KWONG

            # VDW
            elif (
                model_name == 'VAN_DER_WAALS' or
                model_name == 'VAN-DEER-WAALS'
            ):
                model_name_set = VAN_DER_WAALS

            elif (
                model_name == 'vdW' or
                model_name == 'VDW'
            ):
                model_name_set = VAN_DER_WAALS

            # RK
            elif (
                model_name == 'REDLICH_KWONG' or
                model_name == 'REDLICH-KWONG'
            ):
                model_name_set = REDLICH_KWONG

            elif model_name == 'RK':
                model_name_set = REDLICH_KWONG

            else:
                raise Exception('Invalid equation of state name!')
        return model_name_set
    except Exception as e:
        raise Exception('Setting eos model failed!, ', e)


# SECTION: Tau Correlation and Method Mapping
TauCorrelation: TypeAlias = Literal[
    "direct_tau",
    "gibbs_energy",
    "extended_temperature",
    "inverse_temperature",
    "inverse_temperature_squared",
    "inverse_log_temperature",
]

TauMethod: TypeAlias = Literal["M0", "M1", "M2", "M3", "M4", "M5"]
UNIQUACTauMethod: TypeAlias = Literal["M0", "M1", "M2", "M4", "M5", "M6"]


def map_tau_correlation_to_method(
    tau_correlation: TauCorrelation,
) -> TauMethod:
    """
    Map a descriptive tau-correlation name to its corresponding method code.

    Parameters
    ----------
    tau_correlation : TauCorrelation
        Descriptive name of the tau correlation.

    Returns
    -------
    TauMethod
        Corresponding method code from M1 to M5.

    Raises
    ------
    ValueError
        If an unsupported correlation name is provided.
    """
    correlation_map: dict[TauCorrelation, TauMethod] = {
        "direct_tau": "M0",
        "gibbs_energy": "M1",
        "extended_temperature": "M2",
        "inverse_temperature": "M3",
        "inverse_temperature_squared": "M4",
        "inverse_log_temperature": "M5",
    }

    try:
        return correlation_map[tau_correlation]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported tau correlation: {tau_correlation!r}. "
            f"Expected one of: {', '.join(correlation_map)}."
        ) from exc


def map_uniquac_tau_correlation_to_method(
    tau_correlation: TauCorrelation,
) -> UNIQUACTauMethod:
    """
    Map a descriptive tau-correlation name to its UNIQUAC method code.

    Parameters
    ----------
    tau_correlation : TauCorrelation
        Descriptive name of the tau correlation.

    Returns
    -------
    UNIQUACTauMethod
        UNIQUAC-specific method code. The coefficient correlations do not map
        one-to-one with NRTL method numbers because UNIQUAC also exposes M3 as
        a direct energy-over-R helper.

    Raises
    ------
    ValueError
        If an unsupported correlation name is provided.
    """
    # SECTION: UNIQUAC-specific public-name to method mapping
    correlation_map: dict[TauCorrelation, UNIQUACTauMethod] = {
        "direct_tau": "M0",
        "gibbs_energy": "M1",
        "extended_temperature": "M2",
        "inverse_temperature": "M4",
        "inverse_temperature_squared": "M5",
        "inverse_log_temperature": "M6",
    }

    try:
        # NOTE: keep raw Mx method names internal to model dispatch
        return correlation_map[tau_correlation]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported UNIQUAC tau correlation: {tau_correlation!r}. "
            f"Expected one of: {', '.join(correlation_map)}."
        ) from exc
