# UTILITY FUNCTIONS
# ------------------

# import packages/modules
import numpy as np
from typing import Literal
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
            if (model_name == 'PENG_ROBINSON' or
                    model_name == 'PENG-ROBINSON'):
                model_name_set = PENG_ROBINSON

            elif model_name == 'PR':
                model_name_set = PENG_ROBINSON

            # SRK
            elif (model_name == 'SOAVE_REDLICH_KWONG' or
                  model_name == 'SOAVE-REDLICH-KWONG'):
                model_name_set = SOAVE_REDLICH_KWONG

            elif model_name == 'SRK':
                model_name_set = SOAVE_REDLICH_KWONG

            # VDW
            elif (model_name == 'VAN_DER_WAALS' or
                  model_name == 'VAN-DEER-WAALS'):
                model_name_set = VAN_DER_WAALS

            elif (model_name == 'vdW' or
                  model_name == 'VDW'):
                model_name_set = VAN_DER_WAALS

            # RK
            elif (model_name == 'REDLICH_KWONG' or
                  model_name == 'REDLICH-KWONG'):
                model_name_set = REDLICH_KWONG

            elif model_name == 'RK':
                model_name_set = REDLICH_KWONG

            else:
                raise Exception('Invalid equation of state name!')
        return model_name_set
    except Exception as e:
        raise Exception('Setting eos model failed!, ', e)
