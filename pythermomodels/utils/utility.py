# UTILITY FUNCTIONS
# ------------------

# import packages/modules
import numpy as np
# internal
from ..configs import ROUND_FUN_ACCURACY


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
